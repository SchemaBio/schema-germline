"""Reference-bin selection and incremental merge regression tests."""
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from test_cnv_chain import BASH, command
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("update_cnv_reference", ROOT / "scripts/update_cnv_reference.py")
UPDATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(UPDATE)
HEADER = "chromosome\tstart\tend\tgene\tlog2\tdepth\tspread\n"


class ReferenceBinsTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.old = self.root / "old.cnn"
        self.old.write_text(HEADER + (
            "1\t0\t100\tA\t0\t20\t0.2\n"
            "1\t100\t200\tB\t0\t20\t0.2\n"
            "1\t200\t300\tC\t0\t20\t0.2\n"
            "1\t300\t400\tD\t0\t20\t0.2\n"
            "1\t400\t500\tAntitarget\t0\t20\t0.2\n"
            "MT\t0\t100\tMTGENE\t0\t20\t0.2\n"
        ))
        self.bed = self.root / "panel.bed"
        self.selected = self.root / "selected.cnn"
        self.targets = self.root / "targets.bed"
        self.antitargets = self.root / "antitargets.bed"

    def prepare(self, bed=True):
        UPDATE.prepare_reference(self.old, self.bed if bed else None,
                                 self.selected, self.targets, self.antitargets)
        return UPDATE.read_reference(self.selected)[1]

    def test_full_coverage_by_merged_union_preserves_old_coordinates(self):
        self.bed.write_text(
            "track name=panel\n"
            "chr1\t50\t150\nchr1\t0\t50\nchr1\t100\t200\n"
            "chr1\t200\t249\nchr1\t250\t300\n"
            "chr1\t350\t500\nchrM\t0\t100\n"
        )
        rows = self.prepare()
        self.assertEqual([row["gene"] for row in rows], ["A", "B", "Antitarget", "MTGENE"])
        self.assertEqual(self.targets.read_text(), "1\t0\t100\tA\n1\t100\t200\tB\nMT\t0\t100\tMTGENE\n")
        self.assertEqual(self.antitargets.read_text(), "1\t400\t500\tAntitarget\n")

    def test_only_targets_and_compressed_bed(self):
        self.bed = self.root / "panel.bed.gz"
        with gzip.open(self.bed, "wt") as handle:
            handle.write("chr1\t0\t200\n")
        self.assertEqual(len(self.prepare()), 2)
        self.assertEqual(self.antitargets.read_bytes(), b"")

    def test_analysis_reuses_reference_bins_without_rebinning(self):
        self.assertEqual(len(self.prepare(bed=False)), 6)
        self.assertEqual(self.old.read_text(), self.selected.read_text())

    def test_no_complete_bin_and_malformed_bed_fail(self):
        for contents in ["1\t1\t99\n", "1\t-1\t100\n", "1\t100\t100\n", "1\tbad\t100\n", "# empty\n"]:
            with self.subTest(contents=contents):
                self.bed.write_text(contents)
                with self.assertRaises(UPDATE.ReferenceError):
                    self.prepare()
                self.assertFalse(self.selected.exists())

    def test_merge_matches_coordinates_independently_of_row_order(self):
        self.bed.write_text("1\t0\t200\n")
        self.prepare()
        new = self.root / "new.cnn"
        new.write_text(HEADER +
                       "1\t100\t200\tB\t1\t40\t0.4\n"
                       "1\t0\t100\tA\t1\t40\t0.4\n")
        header, rows = UPDATE.merge_references(self.selected, new, 20, 5)
        self.assertEqual(header, HEADER.strip().split("\t"))
        self.assertEqual([row["gene"] for row in rows], ["A", "B"])
        self.assertAlmostEqual(float(rows[0]["log2"]), 0.2)
        self.assertAlmostEqual(float(rows[0]["depth"]), 24)
        self.assertAlmostEqual(float(rows[0]["spread"]), (0.8 * (0.04 + 0.04) + 0.2 * (0.16 + 0.64)) ** 0.5)
        new.write_text(new.read_text().replace("100\t200", "100\t201"))
        with self.assertRaisesRegex(UPDATE.ReferenceError, "bin sets differ"):
            UPDATE.merge_references(self.selected, new, 20, 5)

    def test_duplicate_coordinates_rejected(self):
        with self.old.open("a") as handle:
            handle.write("1\t0\t100\tOTHER\t0\t20\t0.2\n")
        with self.assertRaisesRegex(UPDATE.ReferenceError, "Duplicate"):
            self.prepare(bed=False)


    def test_wdl_preparation_command_uses_selected_bins_and_optional_bed(self):
        self.bed.write_text("1\t0\t200\n")
        for use_bed, expected in [(True, 2), (False, 6)]:
            with self.subTest(use_bed=use_bed):
                body = command("PrepareCNVReference", "tasks/germline.wdl", {
                    "defined(bed)": "true" if use_bed else "false",
                    'default="" bed': self.bed.as_posix() if use_bed else "",
                    "reference": self.old.as_posix(),
                })
                executable = '"' + sys.executable.replace("\\", "/") + '"'
                body = body.replace("python3 /opt/schema-germline/scripts/update_cnv_reference.py",
                                    executable + ' "' + (ROOT / "scripts/update_cnv_reference.py").as_posix() + '"')
                result = subprocess.run([BASH, "-c", body], cwd=self.root, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                rows = UPDATE.read_reference(self.root / "selected.cnvkit.ref.cnn")[1]
                self.assertEqual(len(rows), expected)

    def test_target_only_batch_can_preserve_static_repeat_masker_annotation(self):
        self.bed.write_text("1\t0\t200\n")
        self.prepare()
        with_mask = self.root / "with-mask.cnn"
        header, rows = UPDATE.read_reference(self.selected)
        for row in rows:
            row["rmask"] = "0.1"
        UPDATE.atomic_write_tsv(with_mask, header + ["rmask"], rows)
        _, merged = UPDATE.merge_references(with_mask, self.selected, 20, 3)
        self.assertEqual([row["rmask"] for row in merged], ["0.1", "0.1"])

    def test_cli_prepare_then_merge_produces_filtered_reference_and_report(self):
        self.bed.write_text("1\t0\t200\n")
        script = str(ROOT / "scripts/update_cnv_reference.py")
        result = subprocess.run([sys.executable, script, "--prepare-bins", "--old-reference", str(self.old),
                                 "--bed", str(self.bed), "--output", str(self.selected),
                                 "--target-bed", str(self.targets), "--antitarget-bed", str(self.antitargets)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        output, report = self.root / "merged.cnn", self.root / "report.json"
        result = subprocess.run([sys.executable, script, "--old-reference", str(self.selected),
                                 "--new-reference", str(self.selected), "--old-sample-count", "20",
                                 "--new-sample-count", "3", "--output", str(output), "--report", str(report)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(UPDATE.read_reference(output)[1]), 2)
        self.assertEqual(json.loads(report.read_text())["total_sample_count"], 23)


if __name__ == "__main__":
    unittest.main()
