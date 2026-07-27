"""Benchmark tests evaluating read generation speed and throughput."""

from __future__ import annotations

import time

import pytest

from mock_fastq_generator.core import assemble_sequences

TEMPLATE = "GAAAAAATTAACCAGAGCCTGG"


class TestGenerationBenchmark:
    """Measure read generation throughput for 1,000, 10,000, and 100,000 reads."""

    @pytest.mark.parametrize("num_reads", [1_000, 10_000, 100_000])
    def test_generation_throughput(self, num_reads: int) -> None:
        start_time = time.perf_counter()
        records = assemble_sequences(
            template_sequence=TEMPLATE,
            total_length=150,
            number_of_sequences=num_reads,
            decay_model="gaussian",
            noise_level=0.5,
        )
        elapsed = time.perf_counter() - start_time
        throughput = num_reads / elapsed

        assert len(records) == num_reads
        assert elapsed > 0
        msg = (
            f"\n[BENCHMARK] N={num_reads:,} reads | "
            f"Time: {elapsed:.4f} s | Throughput: {throughput:,.1f} reads/s"
        )
        print(msg)
