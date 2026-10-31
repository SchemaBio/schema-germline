version 1.2

task CNVAnno {
    input {
        String prefix
        File cnv_bed
        String assembly
    }

    command <<<
        set -euo pipefail
        # Normalize only the annotator input, then restore the original columns.
        python - "~{cnv_bed}" "~{prefix}.cnvanno.txt" "~{assembly}" <<'PY'
        import os
        from pathlib import Path
        import subprocess
        import sys
        import tempfile


        def table(path):
            lines = Path(path).read_text(encoding="utf-8").splitlines()
            if not lines:
                raise ValueError("empty CNV table")
            header = lines[0].split("\t")
            rows = [line.split("\t") for line in lines[1:] if line.strip()]
            if len(header) < 4 or not rows:
                raise ValueError("CNV table must contain a header and records")
            if any(len(row) != len(header) for row in rows):
                raise ValueError("CNV table column count mismatch")
            return header, rows


        def main():
            source, output, assembly = sys.argv[1:]
            header, original = table(source)
            aliases = {"DEL": "Deletion", "DELETION": "Deletion", "LOSS": "Deletion",
                       "DUP": "Amplification", "DUPLICATION": "Amplification",
                       "AMPLIFICATION": "Amplification", "GAIN": "Amplification", "NORMAL": "Normal"}
            normalized = []
            for number, row in enumerate(original, start=2):
                start, end = int(row[1]), int(row[2])
                if start < 0 or end <= start:
                    raise ValueError("invalid CNV interval on line " + str(number))
                status = aliases.get(row[3].strip().upper())
                if status is None:
                    raise ValueError("unsupported CNV type on line " + str(number))
                normalized.append(row[:3] + [status] + row[4:])
            destination = Path(output)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
                temporary = Path(temporary)
                prepared, annotated = temporary / "input.bed", temporary / "annotated.tsv"
                prepared.write_text("\n".join("\t".join(row) for row in [header] + normalized) + "\n", encoding="utf-8")
                result = subprocess.run([sys.executable, "/app/CNVAnno/cnvanno.py", str(prepared),
                                         "-d", "/app/CNVAnno/data", "-g", assembly, "-f", "tsv", "-o", str(annotated)])
                if result.returncode:
                    raise SystemExit(result.returncode)
                annotated_header, rows = table(annotated)
                if len(rows) != len(original) or len(annotated_header) <= len(header):
                    raise ValueError("CNVAnno changed the record count or omitted annotation columns")
                for number, row in enumerate(rows):
                    if row[:len(header)] != normalized[number]:
                        raise ValueError("CNVAnno changed the input record order or original columns")
                    # Keep the existing downstream DEL/DUP/Normal and ColN contract.
                    row[:len(header)] = original[number]
                restored = temporary / "restored.tsv"
                restored.write_text("\n".join("\t".join(row) for row in [annotated_header] + rows) + "\n", encoding="utf-8")
                os.replace(restored, destination)


        try:
            main()
        except (ValueError, OSError) as error:
            raise SystemExit("CNV annotation contract error: " + str(error))
        PY
    >>>

    output {
        File cnv_anno_result = "~{prefix}.cnvanno.txt"
    }

    runtime {
        docker: "schemabio/cnvanno:v0.0.2"
    }
}
