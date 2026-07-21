from __future__ import annotations

__version__ = "2.0.0"

from .core import (
    add_noise_to_sequence,
    assemble_paired_sequences,
    assemble_sequences,
    generate_phred_scores,
    generate_read_header,
    reverse_complement,
)
from .io import read_fasta, write_fastq

__all__ = [
    "__version__",
    "add_noise_to_sequence",
    "assemble_paired_sequences",
    "assemble_sequences",
    "generate_phred_scores",
    "generate_read_header",
    "read_fasta",
    "reverse_complement",
    "write_fastq",
]
