from __future__ import annotations

import gzip
import sys
from typing import IO


def read_fasta(path: str) -> str:
    """Read a FASTA file and return the concatenated sequence.

    Header lines (starting with ``>``) and blank lines are skipped.
    All sequence lines are uppercased and joined without separators,
    so multi-line FASTA records are handled correctly.

    Args:
        path: Filesystem path to the FASTA file.

    Returns:
        The full template sequence as an uppercase string.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If no sequence lines are found in the file.
    """
    sequence_lines: list[str] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped and not stripped.startswith(">"):
                sequence_lines.append(stripped.upper())
    if not sequence_lines:
        raise ValueError(f"No sequence data found in FASTA file: {path}")
    return "".join(sequence_lines)


def _write_records(handle: IO[str], sequences: list[list[str]]) -> None:
    """Write FASTQ records to an open file handle.

    Args:
        handle: A writable text stream.
        sequences: List of 4-element lists ``[header, seq, "+", qual]``.
    """
    for record in sequences:
        handle.write(record[0] + "\n")
        handle.write(record[1] + "\n")
        handle.write(record[2] + "\n")
        handle.write(record[3] + "\n")


def write_fastq(
    path: str | None,
    sequences: list[list[str]],
    compress: bool = False,
    use_stdout: bool = False,
) -> None:
    """Write FASTQ records to a file, gzip archive, or stdout.

    Exactly one output mode must be selected:

    * **stdout** — ``use_stdout=True`` streams records to
      ``sys.stdout``.
    * **gzip** — ``compress=True`` writes a gzip-compressed file.
    * **plain** — default, writes an uncompressed text file.

    Args:
        path: Output file path.  Required unless *use_stdout* is True.
        sequences: List of FASTQ records (each a 4-element string list).
        compress: If True, write gzip-compressed output.
        use_stdout: If True, write to ``sys.stdout`` instead of a file.

    Raises:
        ValueError: If *path* is ``None`` and *use_stdout* is ``False``.
    """
    if use_stdout:
        _write_records(sys.stdout, sequences)
    elif path is None:
        raise ValueError("path must be provided when use_stdout is False")
    elif compress:
        with gzip.open(path, "wt") as fh:
            _write_records(fh, sequences)
    else:
        with open(path, "w") as fh:
            _write_records(fh, sequences)
