from __future__ import annotations

import random
import string

import numpy as np


def generate_phred_scores(
    num_points: int,
    center: int,
    min_val: int,
    max_val: int,
    std_dev: int,
    noise_level: float = 0.1,
    ascii_base: int = 33,
) -> str:
    """Generate simulated Phred quality scores with a Gaussian decay profile.

    Produces a string of ASCII-encoded Phred quality characters. The ASCII
    codes follow a Gaussian curve centered at *center* with
    standard deviation *std_dev*.  Additive Gaussian noise (controlled by
    *noise_level*) is applied, and values are clamped to
    ``[min_val, max_val]``.

    Args:
        num_points: Number of quality scores to generate (= read length).
        center: Index position of the Gaussian peak.
        min_val: Minimum allowed ASCII code (e.g., 40 for Phred 7).
        max_val: Maximum allowed ASCII code (e.g., 73 for Phred 40).
        std_dev: Standard deviation of the Gaussian curve.
        noise_level: Standard deviation of additive noise.
        ascii_base: ASCII offset for Phred encoding (default 33 = Sanger).

    Returns:
        A string of ASCII characters encoding the quality scores.
    """
    x = np.linspace(0, num_points - 1, num_points)
    gaussian_curve = np.exp(-((x - center) ** 2) / (2 * std_dev**2))

    noise = np.random.normal(0, noise_level, num_points)
    smooth_noisy_curve = gaussian_curve + noise
    smooth_noisy_curve = smooth_noisy_curve * (max_val - min_val) + min_val
    smooth_noisy_curve = np.clip(smooth_noisy_curve, min_val, max_val)
    smooth_noisy_curve = np.round(smooth_noisy_curve).astype(int)

    phred_scores = "".join(chr(value) for value in smooth_noisy_curve)
    return phred_scores


def generate_read_header() -> str:
    """Generate a random FASTQ read header in ``@XXXX:1234`` format.

    Returns:
        A string starting with ``@`` followed by 4 random alphanumeric
        characters, a colon, and 4 random digits.
    """
    instrument = "".join(random.choices(string.ascii_uppercase + string.digits, k=4))
    run_id = "".join(random.choices(string.digits, k=4))
    return f"@{instrument}:{run_id}"


def add_noise_to_sequence(template_sequence: str, noise_level: float = 0.1) -> str:
    """Introduce random substitution noise into a DNA sequence.

    Each base has a probability of *noise_level* of being replaced by a
    uniformly random nucleotide (A, C, G, or T).

    Args:
        template_sequence: Input DNA string (uppercase).
        noise_level: Per-base substitution probability in ``[0, 1]``.

    Returns:
        The (possibly mutated) DNA string.
    """
    bases = ["A", "C", "G", "T"]
    noisy_sequence = [
        random.choice(bases) if random.random() < noise_level else base
        for base in template_sequence
    ]
    return "".join(noisy_sequence)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Handles both upper- and lower-case characters.

    Args:
        sequence: Input DNA string.

    Returns:
        Reverse-complemented DNA string preserving original case.
    """
    trans = str.maketrans("ATCGatcg", "TAGCtagc")
    return sequence.translate(trans)[::-1]


def assemble_sequences(
    template_sequence: str,
    left_margin: int = 110,
    upstream_sequence: str = "GCCGGCCATGGCG",
    total_length: int = 500,
    number_of_sequences: int = 10,
    center: int = 150,
    min_val: int = 40,
    max_val: int = 73,
    std_dev: int = 75,
    noise_level: float = 0.1,
) -> list[list[str]]:
    """Assemble single-end synthetic FASTQ records.

    Each record is a 4-element list:
    ``[header, sequence, "+", quality_string]``.

    The read is built as::

        upstream + random_left_margin + noisy_template + random_right_fill

    so that the total length equals *total_length*.

    Args:
        template_sequence: Target FASTA template (uppercase DNA).
        left_margin: Random bases between upstream and template.
        upstream_sequence: Adapter / flanking sequence prepended.
        total_length: Total read length.
        number_of_sequences: Number of reads to generate.
        center: Center of the Gaussian quality peak.
        min_val: Minimum Phred score.
        max_val: Maximum Phred score.
        std_dev: Std-dev of the Gaussian quality curve.
        noise_level: Per-base substitution probability **and** quality
            noise std-dev.

    Returns:
        A list of FASTQ records, each a 4-element string list.
    """
    bases = ["A", "C", "G", "T"]
    right_margin = (
        total_length - left_margin - len(upstream_sequence) - len(template_sequence)
    )

    assembled_sequences: list[list[str]] = []
    for _ in range(number_of_sequences):
        random_left = "".join(random.choices(bases, k=left_margin))
        random_right = "".join(random.choices(bases, k=right_margin))
        noisy_template = add_noise_to_sequence(template_sequence, noise_level)
        sequence = upstream_sequence + random_left + noisy_template + random_right

        header = generate_read_header()
        seq_len = len(sequence)
        phred_string = generate_phred_scores(
            seq_len, center, min_val, max_val, std_dev, noise_level=noise_level
        )
        assembled_sequences.append([header, sequence, "+", phred_string])

    return assembled_sequences


def assemble_paired_sequences(
    template_sequence: str,
    left_margin: int = 110,
    upstream_sequence: str = "GCCGGCCATGGCG",
    total_length: int = 500,
    number_of_sequences: int = 10,
    center: int = 150,
    min_val: int = 40,
    max_val: int = 73,
    std_dev: int = 75,
    noise_level: float = 0.1,
) -> tuple[list[list[str]], list[list[str]]]:
    """Assemble paired-end synthetic FASTQ records.

    R1 reads are forward-strand constructs identical to those produced by
    :func:`assemble_sequences`.  R2 reads are the reverse complement of
    R1, with independently generated quality scores.

    Headers share the same instrument/run ID, distinguished by ``/1``
    and ``/2`` suffixes.

    Args:
        template_sequence: Target FASTA template (uppercase DNA).
        left_margin: Random bases between upstream and template.
        upstream_sequence: Adapter / flanking sequence prepended.
        total_length: Total read length.
        number_of_sequences: Number of read pairs to generate.
        center: Center of the Gaussian quality peak.
        min_val: Minimum Phred score.
        max_val: Maximum Phred score.
        std_dev: Std-dev of the Gaussian quality curve.
        noise_level: Per-base substitution probability **and** quality
            noise std-dev.

    Returns:
        A tuple ``(r1_records, r2_records)`` where each element is a list
        of 4-element FASTQ records.
    """
    bases = ["A", "C", "G", "T"]
    right_margin = (
        total_length - left_margin - len(upstream_sequence) - len(template_sequence)
    )

    r1_records: list[list[str]] = []
    r2_records: list[list[str]] = []

    for _ in range(number_of_sequences):
        random_left = "".join(random.choices(bases, k=left_margin))
        random_right = "".join(random.choices(bases, k=right_margin))
        noisy_template = add_noise_to_sequence(template_sequence, noise_level)
        sequence = upstream_sequence + random_left + noisy_template + random_right

        header = generate_read_header()
        seq_len = len(sequence)

        # R1 — forward
        phred_r1 = generate_phred_scores(
            seq_len, center, min_val, max_val, std_dev, noise_level=noise_level
        )
        r1_records.append([header + "/1", sequence, "+", phred_r1])

        # R2 — reverse complement
        r2_sequence = reverse_complement(sequence)
        phred_r2 = generate_phred_scores(
            seq_len, center, min_val, max_val, std_dev, noise_level=noise_level
        )
        r2_records.append([header + "/2", r2_sequence, "+", phred_r2])

    return r1_records, r2_records
