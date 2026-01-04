#!/usr/bin/env python3
"""
BAM/CRAM format converter using samtools.

Supports conversion between:
- BAM to CRAM
- CRAM to BAM
- Indexed BAM/CRAM generation

Usage:
    python bam_cram_converter.py input.bam --output output.cram --reference hg38.fa
    python bam_cram_converter.py input.cram --output output.bam
    python bam_cram_converter.py input.bam --to cram --reference hg38.fa
    python bam_cram_converter.py input.bam output.cram --reference hg38.fa
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path
from typing import Optional


def check_samtools() -> bool:
    """Check if samtools is available."""
    return shutil.which('samtools') is not None


def get_file_info(file_path: Path) -> dict:
    """Get file information using samtools."""
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    info = {
        'path': str(file_path),
        'name': file_path.name,
        'size': file_path.stat().st_size,
        'format': None,
        'indexed': False
    }

    # Determine format from extension
    suffix = file_path.suffix.lower()
    if suffix == '.bam':
        info['format'] = 'BAM'
    elif suffix == '.cram':
        info['format'] = 'CRAM'

    # Check for index
    index_suffix = '.bai' if info['format'] == 'BAM' else '.crai'
    index_file = file_path.with_suffix(suffix + index_suffix)
    if index_file.exists():
        info['indexed'] = True

    # Get detailed stats with samtools
    try:
        result = subprocess.run(
            ['samtools', 'idxstats', str(file_path)],
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            info['total_reads'] = sum(
                int(line.split('\t')[2])
                for line in lines if line and not line.startswith('@')
            )
            info['chromosomes'] = len([l for l in lines if l and not l.startswith('@')])
    except (subprocess.TimeoutExpired, Exception):
        pass

    return info


def bam_to_cram(
    input_file: Path,
    output_file: Path,
    reference: Path,
    threads: int = 4,
    preserve_level: bool = True,
    extra_args: Optional[list] = None
) -> Path:
    """Convert BAM to CRAM."""
    if not check_samtools():
        raise RuntimeError("samtools not found. Please install samtools first.")

    # Check reference
    if not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")

    # Check for reference index
    if not reference.with_suffix('.fai').exists():
        print(f"Creating reference index for {reference.name}...")
        subprocess.run(
            ['samtools', 'faidx', str(reference)],
            check=True,
            timeout=300
        )

    # Build command
    cmd = [
        'samtools',
        'view',
        '-@', str(threads),
        '-C',
        '--reference', str(reference),
        '-o', str(output_file),
        str(input_file)
    ]

    # Add compression level
    if preserve_level:
        cmd.insert(3, '-l')
        cmd.insert(4, '1')  # Use maximum compression

    if extra_args:
        cmd.extend(extra_args)

    print(f"Converting {input_file.name} -> {output_file.name}...")
    print(f"Reference: {reference.name}")

    subprocess.run(cmd, check=True, timeout=3600)

    # Create index
    index_file = output_file.with_suffix('.cram.crai')
    if not index_file.exists():
        print(f"Creating index for {output_file.name}...")
        subprocess.run(
            ['samtools', 'index', '-@', str(threads), str(output_file)],
            check=True,
            timeout=300
        )

    return output_file


def cram_to_bam(
    input_file: Path,
    output_file: Path,
    reference: Optional[Path] = None,
    threads: int = 4,
    extra_args: Optional[list] = None
) -> Path:
    """Convert CRAM to BAM."""
    if not check_samtools():
        raise RuntimeError("samtools not found. Please install samtools first.")

    # Build command
    cmd = [
        'samtools',
        'view',
        '-@', str(threads),
        '-b',
        '-o', str(output_file),
        str(input_file)
    ]

    if reference:
        cmd.extend(['--reference', str(reference)])

    if extra_args:
        cmd.extend(extra_args)

    print(f"Converting {input_file.name} -> {output_file.name}...")

    subprocess.run(cmd, check=True, timeout=3600)

    # Create index
    index_file = output_file.with_suffix('.bam.bai')
    if not index_file.exists():
        print(f"Creating index for {output_file.name}...")
        subprocess.run(
            ['samtools', 'index', '-@', str(threads), str(output_file)],
            check=True,
            timeout=300
        )

    return output_file


def reindex_file(
    input_file: Path,
    threads: int = 4
) -> Path:
    """Re-index a BAM/CRAM file."""
    if not check_samtools():
        raise RuntimeError("samtools not found. Please install samtools first.")

    suffix = input_file.suffix.lower()
    is_bam = suffix == '.bam'

    index_file = input_file.with_suffix(suffix + ('.bai' if is_bam else '.crai'))

    print(f"Re-indexing {input_file.name}...")

    subprocess.run(
        ['samtools', 'index', '-@', str(threads), str(input_file)],
        check=True,
        timeout=300
    )

    return index_file


def get_format_from_path(path: Path) -> str:
    """Get format from file path."""
    suffix = path.suffix.lower()
    if suffix == '.bam':
        return 'BAM'
    elif suffix == '.cram':
        return 'CRAM'
    return 'UNKNOWN'


def detect_conversion(input_path: Path, output_path: Optional[Path] = None) -> tuple:
    """Detect conversion type from input and output paths."""
    input_format = get_format_from_path(input_path)

    if output_path:
        output_format = get_format_from_path(output_path)
    else:
        # Infer output format from input
        output_format = 'CRAM' if input_format == 'BAM' else 'BAM'

    return input_format, output_format


def batch_convert(
    input_dir: Path,
    output_dir: Path,
    reference: Optional[Path] = None,
    from_format: Optional[str] = None,
    to_format: str = 'CRAM',
    threads: int = 4
) -> list:
    """Batch convert all BAM/CRAM files in a directory."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []

    # Find files
    if from_format:
        pattern = f'*.{from_format.lower()}'
    else:
        pattern = '*.[bc]ram'

    input_files = list(input_dir.glob(pattern))

    if not input_files:
        print(f"No files found in {input_dir} matching {pattern}")
        return results

    for input_file in sorted(input_files):
        input_format = get_format_from_path(input_file)

        # Determine output name and format
        if to_format.upper() == 'CRAM' and input_format == 'BAM':
            output_file = output_dir / input_file.with_suffix('.cram').name
            conv_func = bam_to_cram
            conv_args = {'reference': reference}
        elif to_format.upper() == 'BAM' and input_format == 'CRAM':
            output_file = output_dir / input_file.with_suffix('.bam').name
            conv_func = cram_to_bam
            conv_args = {'reference': reference}
        else:
            print(f"Skipping {input_file.name} - already in {to_format} format")
            continue

        try:
            result = conv_func(input_file, output_file, threads=threads, **conv_args)
            results.append({
                'input': str(input_file),
                'output': str(result),
                'status': 'success'
            })
            print(f"  -> {result.name}")
        except Exception as e:
            results.append({
                'input': str(input_file),
                'output': None,
                'status': 'error',
                'error': str(e)
            })
            print(f"  -> ERROR: {e}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Convert between BAM and CRAM formats using samtools',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # BAM to CRAM
  python bam_cram_converter.py input.bam --output output.cram --reference hg38.fa

  # CRAM to BAM
  python bam_cram_converter.py input.cram --output output.bam

  # Quick conversion (auto-detect format)
  python bam_cram_converter.py input.bam output.cram --reference hg38.fa

  # Batch convert all BAM to CRAM
  python bam_cram_converter.py --batch input_dir/ --to cram --reference hg38.fa

  # Re-index a file
  python bam_cram_converter.py input.bam --reindex
        """
    )

    # Positional arguments (can be used as input/output)
    parser.add_argument('input', nargs='?', type=Path, help='Input BAM/CRAM file')

    # Output
    parser.add_argument('-o', '--output', type=Path, help='Output file path')

    # Format conversion
    parser.add_argument('--to', choices=['bam', 'cram'], dest='to_format',
                        help='Output format (auto-detected if not specified)')

    # Reference
    parser.add_argument('-r', '--reference', type=Path,
                        help='Reference FASTA file (required for CRAM conversion)')

    # Batch mode
    parser.add_argument('--batch', action='store_true',
                        help='Batch mode: convert all files in input directory')
    parser.add_argument('--from', choices=['bam', 'cram'], dest='from_format',
                        help='Input format filter for batch mode')

    # Indexing
    parser.add_argument('--reindex', action='store_true',
                        help='Re-index the input file without conversion')

    # Options
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads (default: 4)')
    parser.add_argument('-c', '--compression', type=int, default=1,
                        choices=range(0, 10), metavar='0-9',
                        help='Compression level for CRAM (0-9, default: 1)')

    # Info
    parser.add_argument('--info', action='store_true',
                        help='Show information about the input file')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show command without executing')

    args = parser.parse_args()

    # Check samtools
    if not check_samtools():
        print("Error: samtools not found in PATH", file=sys.stderr)
        print("Please install samtools: https://github.com/samtools/samtools", file=sys.stderr)
        sys.exit(1)

    # Info mode
    if args.info and args.input:
        try:
            info = get_file_info(args.input)
            print(f"File: {info['name']}")
            print(f"  Format: {info['format']}")
            print(f"  Size: {info['size'] / 1024 / 1024:.2f} MB")
            print(f"  Indexed: {info['indexed']}")
            if 'total_reads' in info:
                print(f"  Total reads: {info['total_reads']:,}")
                print(f"  Chromosomes: {info['chromosomes']}")
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
        sys.exit(0)

    # Reindex mode
    if args.reindex and args.input:
        try:
            reindex_file(args.input, args.threads)
            print(f"Re-indexed: {args.input.name}")
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        sys.exit(0)

    # Validate input
    if not args.input:
        parser.print_help()
        sys.exit(1)

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Batch mode
    if args.batch:
        if not args.input.is_dir():
            print(f"Error: Input must be a directory for batch mode: {args.input}", file=sys.stderr)
            sys.exit(1)

        to_format = args.to_format or 'cram'
        print(f"Batch converting {args.input}/* -> {to_format.upper()}")

        results = batch_convert(
            args.input,
            args.input.parent / f"converted_{to_format}",
            args.reference,
            args.from_format,
            to_format,
            args.threads
        )

        success = sum(1 for r in results if r['status'] == 'success')
        print(f"\nConverted {success}/{len(results)} files")
        sys.exit(0)

    # Determine output path
    output_path = args.output
    if not output_path:
        input_format, output_format = detect_conversion(args.input)
        if args.to_format:
            output_format = args.to_format.upper()
        suffix = f'.{output_format.lower()}'
        output_path = args.input.with_suffix(suffix)

    # Validate CRAM conversion requires reference
    input_format = get_format_from_path(args.input)
    output_format = get_format_from_path(output_path) if output_path.suffix.lower() in ['.bam', '.cram'] else None

    if output_format == 'CRAM' and not args.reference:
        print("Error: Reference file is required for CRAM output", file=sys.stderr)
        print("Use: --reference /path/to/reference.fa", file=sys.stderr)
        sys.exit(1)

    # Build conversion arguments
    extra_args = []
    if args.compression < 9:
        extra_args.extend(['-l', str(args.compression)])

    # Show command if dry-run
    if args.dry_run:
        print(f"Command: samtools view -@ {args.threads} ... {'-C' if output_format == 'CRAM' else '-b'} ...")
        sys.exit(0)

    # Perform conversion
    try:
        if input_format == 'BAM' and output_format == 'CRAM':
            bam_to_cram(
                args.input, output_path, args.reference,
                args.threads, extra_args=extra_args
            )
        elif input_format == 'CRAM' and output_format == 'BAM':
            cram_to_bam(
                args.input, output_path, args.reference,
                args.threads
            )
        else:
            print(f"Error: Cannot convert from {input_format} to {output_format}", file=sys.stderr)
            sys.exit(1)

        print(f"\nOutput: {output_path}")
        print(f"Size: {output_path.stat().st_size / 1024 / 1024:.2f} MB")

    except subprocess.CalledProcessError as e:
        print(f"Error: samtools failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
