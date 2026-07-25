#!/usr/bin/env python3
"""Incrementally merge two CNVKit reference .cnn files.

CNVKit references contain cohort summaries rather than the per-sample values
used to build them. This script therefore performs a sample-count-weighted
approximation; it cannot reproduce a full rebuild from the original coverage
files exactly.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any


KEY_COLUMNS = ("chromosome", "start", "end", "gene")
MERGED_COLUMNS = ("log2", "depth", "spread")


class ReferenceError(ValueError):
    """Raised when reference files cannot be merged safely."""


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("sample counts must be greater than zero")
    return parsed


def read_reference(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ReferenceError(f"Reference has no header: {path}")
        header = list(reader.fieldnames)
        required = set(KEY_COLUMNS) | set(MERGED_COLUMNS)
        missing = sorted(required.difference(header))
        if missing:
            raise ReferenceError(
                f"Reference {path} is missing required columns: {', '.join(missing)}"
            )
        rows = list(reader)

    if not rows:
        raise ReferenceError(f"Reference has no data rows: {path}")
    return header, rows


def parse_number(row: dict[str, str], column: str, path: Path, line: int) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as exc:
        raise ReferenceError(
            f"Invalid {column!r} value in {path} at data row {line}: "
            f"{row.get(column)!r}"
        ) from exc
    if not math.isfinite(value):
        raise ReferenceError(
            f"Non-finite {column!r} value in {path} at data row {line}: {value}"
        )
    return value


def weighted_mean(old: float, new: float, old_count: int, new_count: int) -> float:
    return (old * old_count + new * new_count) / (old_count + new_count)


def pooled_spread(
    old_mean: float,
    new_mean: float,
    old_spread: float,
    new_spread: float,
    old_count: int,
    new_count: int,
) -> float:
    """Pool within-batch spread and the between-batch log2 displacement."""
    combined_mean = weighted_mean(old_mean, new_mean, old_count, new_count)
    total_count = old_count + new_count
    variance = (
        old_count * (old_spread**2 + (old_mean - combined_mean) ** 2)
        + new_count * (new_spread**2 + (new_mean - combined_mean) ** 2)
    ) / total_count
    return math.sqrt(max(variance, 0.0))


def format_number(value: float) -> str:
    return f"{value:.12g}"


def validate_row_identity(
    old_row: dict[str, str],
    new_row: dict[str, str],
    line: int,
) -> None:
    mismatches = [
        column
        for column in KEY_COLUMNS
        if old_row.get(column) != new_row.get(column)
    ]
    if mismatches:
        old_key = tuple(old_row.get(column, "") for column in KEY_COLUMNS)
        new_key = tuple(new_row.get(column, "") for column in KEY_COLUMNS)
        raise ReferenceError(
            f"Bin mismatch at data row {line} ({', '.join(mismatches)}): "
            f"old={old_key}, new={new_key}"
        )


def merge_references(
    old_path: Path,
    new_path: Path,
    old_count: int,
    new_count: int,
) -> tuple[list[str], list[dict[str, str]]]:
    old_header, old_rows = read_reference(old_path)
    new_header, new_rows = read_reference(new_path)

    if old_header != new_header:
        raise ReferenceError(
            "Reference headers differ. Both references must be produced by the "
            f"same CNVKit setup. old={old_header}, new={new_header}"
        )
    if len(old_rows) != len(new_rows):
        raise ReferenceError(
            f"Reference bin counts differ: old={len(old_rows)}, new={len(new_rows)}"
        )

    merged_rows: list[dict[str, str]] = []
    for line, (old_row, new_row) in enumerate(zip(old_rows, new_rows), start=1):
        validate_row_identity(old_row, new_row, line)

        old_log2 = parse_number(old_row, "log2", old_path, line)
        new_log2 = parse_number(new_row, "log2", new_path, line)
        old_depth = parse_number(old_row, "depth", old_path, line)
        new_depth = parse_number(new_row, "depth", new_path, line)
        old_spread = parse_number(old_row, "spread", old_path, line)
        new_spread = parse_number(new_row, "spread", new_path, line)

        merged = dict(old_row)
        merged["log2"] = format_number(
            weighted_mean(old_log2, new_log2, old_count, new_count)
        )
        merged["depth"] = format_number(
            weighted_mean(old_depth, new_depth, old_count, new_count)
        )
        merged["spread"] = format_number(
            pooled_spread(
                old_log2,
                new_log2,
                old_spread,
                new_spread,
                old_count,
                new_count,
            )
        )
        merged_rows.append(merged)

    return old_header, merged_rows


def atomic_write_tsv(
    output_path: Path,
    header: list[str],
    rows: list[dict[str, str]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        text=True,
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=header,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, output_path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_json(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        text=True,
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, output_path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Merge an existing CNVKit reference with a newly built batch "
            "reference using sample-count weights."
        )
    )
    parser.add_argument("--old-reference", required=True, type=Path)
    parser.add_argument("--new-reference", required=True, type=Path)
    parser.add_argument("--old-sample-count", required=True, type=positive_int)
    parser.add_argument("--new-sample-count", required=True, type=positive_int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        header, merged_rows = merge_references(
            args.old_reference,
            args.new_reference,
            args.old_sample_count,
            args.new_sample_count,
        )
        atomic_write_tsv(args.output, header, merged_rows)
        total_count = args.old_sample_count + args.new_sample_count
        atomic_write_json(
            args.report,
            {
                "approximate": True,
                "bin_count": len(merged_rows),
                "method": "sample_count_weighted_reference_merge",
                "new_reference": str(args.new_reference),
                "new_sample_count": args.new_sample_count,
                "new_weight": args.new_sample_count / total_count,
                "old_reference": str(args.old_reference),
                "old_sample_count": args.old_sample_count,
                "old_weight": args.old_sample_count / total_count,
                "output_reference": str(args.output),
                "total_sample_count": total_count,
                "warning": (
                    "CNVKit reference files contain cohort summaries. This "
                    "incremental result is not identical to rebuilding from "
                    "all original per-sample coverage files."
                ),
            },
        )
    except ReferenceError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
