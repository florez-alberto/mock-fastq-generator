"""Tests for gzip compression in output files."""

from __future__ import annotations

import gzip
from pathlib import Path

from mock_fastq_generator.io import write_fastq


class TestGzipOutput:
    """Verifies the --gzip functionality works correctly."""

    def test_write_gzip_file(self, tmp_path: Path) -> None:
        """It should write a valid gzip file that can be decompressed."""
        records = [
            ["@READ1:1234", "ATCG", "+", "IIII"],
            ["@READ2:5678", "GCTA", "+", "::::"],
        ]

        out_file = tmp_path / "out.fastq.gz"
        write_fastq(str(out_file), records, compress=True)

        assert out_file.exists()

        # Verify it's actually gzipped by attempting to read with gzip module
        with gzip.open(out_file, "rt") as fh:
            content = fh.read()

        expected = "@READ1:1234\nATCG\n+\nIIII\n@READ2:5678\nGCTA\n+\n::::\n"
        assert content == expected
