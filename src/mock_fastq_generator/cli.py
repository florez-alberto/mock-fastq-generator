from __future__ import annotations

import argparse
import json
import os

from .core import assemble_paired_sequences, assemble_sequences
from .io import read_fasta, write_fastq


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line argument parser.

    Returns:
        Configured ``argparse.ArgumentParser`` instance.
    """
    parser = argparse.ArgumentParser(
        prog="mock-fastq-generator",
        description=(
            "A lightweight Python utility for targeted synthetic FASTQ "
            "generation with parametric quality decay modeling."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Input ---
    parser.add_argument(
        "--template_sequence",
        type=str,
        default=None,
        help="Path to a FASTA file containing the target template sequence.",
    )
    parser.add_argument(
        "--parameters_file",
        type=str,
        default=None,
        help="JSON file with generation parameters (overrides CLI defaults).",
    )

    # Kept for backward compatibility with v1.x scripts; silently ignored.
    parser.add_argument(
        "--with_parameters",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    # --- Generation parameters ---
    parser.add_argument(
        "--upstream_sequence",
        type=str,
        default="GCCGGCCATGGCG",
        help="Upstream adapter / flanking sequence prepended to each read.",
    )
    parser.add_argument(
        "--left_margin",
        type=int,
        default=15,
        help="Number of random bases inserted between the upstream sequence "
        "and the template.",
    )
    parser.add_argument(
        "--total_length",
        type=int,
        default=500,
        help="Total length of each generated read (in bases).",
    )
    parser.add_argument(
        "--number_of_sequences",
        type=int,
        default=10,
        help="Number of reads to generate.",
    )
    parser.add_argument(
        "--center",
        type=int,
        default=150,
        help="Center index of the Gaussian quality peak.",
    )
    parser.add_argument(
        "--min_val",
        type=int,
        default=40,
        help="Minimum Phred quality score (clamped).",
    )
    parser.add_argument(
        "--max_val",
        type=int,
        default=73,
        help="Maximum Phred quality score (clamped).",
    )
    parser.add_argument(
        "--std_dev",
        type=int,
        default=75,
        help="Standard deviation of the Gaussian quality curve.",
    )
    parser.add_argument(
        "--noise_level",
        type=float,
        default=0.1,
        help="Standard deviation of additive Gaussian noise on quality scores.",
    )

    # --- Output mode ---
    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument(
        "--output_file",
        type=str,
        default=None,
        help="Output file path.  For paired-end mode this is used as a "
        "prefix (producing <prefix>_R1.fastq and <prefix>_R2.fastq).",
    )
    output_group.add_argument(
        "--stdout",
        action="store_true",
        default=False,
        help="Stream FASTQ records to standard output instead of a file.",
    )

    # --- Optional flags ---
    parser.add_argument(
        "--gzip",
        action="store_true",
        default=False,
        help="Compress output files with gzip (.fastq.gz).",
    )
    parser.add_argument(
        "--paired-end",
        action="store_true",
        default=False,
        help="Generate paired-end reads (R1 forward + R2 reverse complement).",
    )

    return parser


def _resolve_parameters(parsed_args: argparse.Namespace) -> dict:
    """Return a dict of generation parameters from CLI args or JSON file.

    If ``--parameters_file`` is supplied, values are loaded from the JSON
    file; any keys missing in the file fall back to CLI defaults.
    Otherwise, CLI arguments (which already carry defaults) are used
    directly.

    Args:
        parsed_args: Parsed command-line namespace.

    Returns:
        Dictionary with keys matching ``assemble_sequences`` parameters.
    """
    defaults = {
        "upstream_sequence": parsed_args.upstream_sequence,
        "left_margin": parsed_args.left_margin,
        "total_length": parsed_args.total_length,
        "number_of_sequences": parsed_args.number_of_sequences,
        "center": parsed_args.center,
        "min_val": parsed_args.min_val,
        "max_val": parsed_args.max_val,
        "std_dev": parsed_args.std_dev,
        "noise_level": parsed_args.noise_level,
    }

    if parsed_args.parameters_file is not None:
        with open(parsed_args.parameters_file) as fh:
            file_params = json.load(fh)
        for key in defaults:
            if key in file_params:
                defaults[key] = file_params[key]

    return defaults


def main(args: argparse.Namespace | None = None) -> None:
    """Entry point for the mock-fastq-generator CLI.

    Args:
        args: Pre-parsed namespace for programmatic use.  When *None*,
            arguments are parsed from ``sys.argv``.
    """
    if args is None:
        parser = build_parser()
        parsed_args = parser.parse_args()
    else:
        parsed_args = args

    # --- Resolve template sequence ---
    if parsed_args.template_sequence:
        template_sequence = read_fasta(parsed_args.template_sequence)
    else:
        # Fall back to bundled example shipped with the repository.
        default_fasta = os.path.join(
            os.path.dirname(__file__), "..", "..", "template_example.fasta"
        )
        if os.path.isfile(default_fasta):
            template_sequence = read_fasta(default_fasta)
        else:
            template_sequence = read_fasta("template_example.fasta")

    # --- Resolve generation parameters ---
    params = _resolve_parameters(parsed_args)

    # --- Generate & write ---
    paired_end: bool = getattr(parsed_args, "paired_end", False)

    if paired_end:
        r1_records, r2_records = assemble_paired_sequences(
            template_sequence=template_sequence, **params
        )
        if parsed_args.stdout:
            interleaved: list[list[str]] = []
            for r1, r2 in zip(r1_records, r2_records):
                interleaved.append(r1)
                interleaved.append(r2)
            write_fastq(None, interleaved, use_stdout=True)
        else:
            ext = ".fastq.gz" if parsed_args.gzip else ".fastq"
            out_prefix = parsed_args.output_file
            if out_prefix.endswith(".fastq"):
                out_prefix = out_prefix[: -len(".fastq")]
            elif out_prefix.endswith(".fq"):
                out_prefix = out_prefix[: -len(".fq")]

            write_fastq(f"{out_prefix}_R1{ext}", r1_records, compress=parsed_args.gzip)
            write_fastq(f"{out_prefix}_R2{ext}", r2_records, compress=parsed_args.gzip)
    else:
        assembled_sequences = assemble_sequences(
            template_sequence=template_sequence, **params
        )
        if parsed_args.stdout:
            write_fastq(None, assembled_sequences, use_stdout=True)
        else:
            output_path = parsed_args.output_file
            if parsed_args.gzip and not output_path.endswith(".gz"):
                output_path += ".gz"
            write_fastq(output_path, assembled_sequences, compress=parsed_args.gzip)


if __name__ == "__main__":
    main()
