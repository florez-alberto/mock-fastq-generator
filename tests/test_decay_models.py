"""Unit tests for exponential, sigmoidal, binning, and homopolymer decay profiles."""

from __future__ import annotations

from mock_fastq_generator.core import assemble_sequences, generate_phred_scores


class TestDecayModels:
    """Verify core decay model calculations and quality adjustments."""

    def test_exponential_decay_model(self) -> None:
        sequence = "ACGTACGTACGT"
        scores_ascii = generate_phred_scores(
            sequence=sequence,
            center=100,
            min_val=40,
            max_val=73,
            std_dev=300,
            noise_level=0.0,
            decay_model="exponential",
            decay_rate=0.001,
        )
        assert len(scores_ascii) == len(sequence)

    def test_sigmoidal_decay_model(self) -> None:
        sequence = "ACGTACGTACGT"
        scores_ascii = generate_phred_scores(
            sequence=sequence,
            center=5,
            min_val=40,
            max_val=73,
            std_dev=300,
            noise_level=0.0,
            decay_model="sigmoidal",
            decay_rate=1.0,
        )
        assert len(scores_ascii) == len(sequence)

    def test_homopolymer_penalty(self) -> None:
        # Sequence with a run of 6 A's (homopolymer > 4)
        sequence = "AAAAAAACGT"
        scores_no_penalty = generate_phred_scores(
            sequence=sequence,
            center=100,
            min_val=40,
            max_val=73,
            std_dev=300,
            noise_level=0.0,
            homopolymer_penalty=False,
        )
        scores_penalty = generate_phred_scores(
            sequence=sequence,
            center=100,
            min_val=40,
            max_val=73,
            std_dev=300,
            noise_level=0.0,
            homopolymer_penalty=True,
        )
        # Position 5 and 6 should be lower in scores_penalty due to penalty (-10)
        ord_no_pen = [ord(c) for c in scores_no_penalty]
        ord_pen = [ord(c) for c in scores_penalty]
        assert ord_pen[5] < ord_no_pen[5]
        assert ord_pen[6] < ord_no_pen[6]

    def test_binned_quality(self) -> None:
        sequence = "ACGTACGTACGTACGT"
        scores_ascii = generate_phred_scores(
            sequence=sequence,
            center=100,
            min_val=33,
            max_val=73,
            std_dev=300,
            noise_level=0.5,
            binned_quality=True,
        )
        # NovaSeq X default 3-state bins in ASCII (33 + [12, 24, 40]) => [45, 57, 73]
        valid_ascii_bins = {45, 57, 73}
        for c in scores_ascii:
            assert ord(c) in valid_ascii_bins

    def test_custom_bins(self) -> None:
        sequence = "ACGTACGTACGTACGT"
        scores_ascii = generate_phred_scores(
            sequence=sequence,
            center=100,
            min_val=33,
            max_val=73,
            std_dev=300,
            noise_level=0.5,
            custom_bins="12,23,37",
        )
        # NovaSeq 6000 custom bins in ASCII (33 + [12, 23, 37]) => [45, 56, 70]
        valid_ascii_bins = {45, 56, 70}
        for c in scores_ascii:
            assert ord(c) in valid_ascii_bins

    def test_assemble_sequences_phred_score_type(self) -> None:
        records = assemble_sequences(
            template_sequence="ACGTACGT",
            total_length=50,
            number_of_sequences=2,
            min_val=7,
            max_val=40,
            score_type="phred",
        )
        assert len(records) == 2
