# mock-fastq-generator

<a href="https://zenodo.org/doi/10.5281/zenodo.10899656"><img src="https://zenodo.org/badge/756591816.svg" alt="DOI"></a>
![Tests](https://github.com/florez-alberto/mock-fastq-generator/actions/workflows/test.yml/badge.svg)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)

`mock-fastq-generator` is an open-source software utility for generating targeted synthetic FASTQ datasets with controllable Phred quality score profiles. Available both as an installable Python package/CLI and as a standalone client-side web application, it allows bioinformatics developers to benchmark sequence processing tools (trimmers, aligners, denoisers) without relying on heavy empirical reference-genome simulators or exposing proprietary data.

---

## Key Features

- **Parametric Quality Decay Models**: Choose between three distinct decay models:
  - **Gaussian**: Bounded bell-curve quality peak centered at a designated position.
  - **Exponential**: Simulates 3' end quality degradation with a sub-exponential initial plateau.
  - **Sigmoidal**: Sharp quality score drop centered at a specific nucleotide coordinate.
- **Illumina 4-State Quality Binning (`--binned_quality`)**: Snaps Phred scores to discrete bins (`Q2, Q12, Q23, Q37`) emulating modern Illumina NovaSeq platforms.
- **Homopolymer Quality Penalties (`--homopolymer_penalty`)**: Applies context-dependent $-10$ Phred score penalties inside runs of $>4$ identical consecutive nucleotides.
- **Targeted Construct Assembly**: Combines custom FASTA templates with adapter/flanking upstream sequences and randomized nucleotide margins.
- **Paired-End Read Simulation (`--paired-end`)**: Generates matching R1 (forward) and R2 (reverse complement) paired FASTQ files with matching headers (`/1` and `/2`).
- **Inline Gzip Compression (`--gzip`)**: Directly writes `.fastq.gz` archives.
- **Standard Output Streaming (`--stdout`)**: Streams FASTQ records directly to `stdout` for UNIX piping (`|`) into tools like `bwa`, `fastp`, or `fastqc`.
- **JSON Configuration Profiles (`--parameters_file`)**: Load reproducible simulation settings from structured JSON files.
- **Standalone Web Application**: Complete client-side web interface (`public/index.html`) requiring no backend server.

---

![Quality Profiles](./img/smooth_curve_example.png)

*Figure 1: Comparison of synthetic Phred quality score profiles generated across a 500-base construct. (A) Parametric Gaussian decay. (B) Exponential 3' quality dropoff. (C) Sigmoidal drop combined with 4-state NovaSeq binning and homopolymer penalties.*

---

## Installation

### Prerequisites
- Python **3.8+**
- `numpy >= 1.20`

### Package Installation

Install directly from GitHub via `pip`:

```bash
git clone https://github.com/florez-alberto/mock-fastq-generator.git
cd mock-fastq-generator
pip install .
```

For active development (includes `pytest` and `ruff`):

```bash
pip install -e ".[dev]"
```

---

## Usage

### Command-Line Interface (CLI)

#### 1. Basic Single-End Generation
```bash
mock-fastq-generator \
  --template_sequence template_example.fasta \
  --output_file synthetic.fastq
```

#### 2. Advanced Paired-End with NovaSeq Binning & Gzip Compression
```bash
mock-fastq-generator \
  --template_sequence template_example.fasta \
  --output_file sample_data \
  --paired-end \
  --gzip \
  --decay_model sigmoidal \
  --decay_rate 0.05 \
  --center 100 \
  --binned_quality \
  --homopolymer_penalty \
  --number_of_sequences 1000
```
*(Produces `sample_data_R1.fastq.gz` and `sample_data_R2.fastq.gz`)*

#### 3. UNIX Piping (Streaming directly into aligner or trimmer)
```bash
mock-fastq-generator \
  --template_sequence template_example.fasta \
  --number_of_sequences 500 \
  --stdout | bwa mem reference.fasta - > mapped_reads.sam
```

#### 4. JSON Configuration File
Create a `config.json` file:
```json
{
  "upstream_sequence": "GCCGGCCATGGCG",
  "left_margin": 15,
  "total_length": 300,
  "number_of_sequences": 200,
  "decay_model": "exponential",
  "decay_rate": 0.001,
  "min_val": 7,
  "max_val": 40,
  "score_type": "phred"
}
```
And execute:
```bash
mock-fastq-generator --template_sequence template_example.fasta --parameters_file config.json --output_file output.fastq
```

---

### Python API

You can also import `mock-fastq-generator` into Python scripts or automated test suites:

```python
from mock_fastq_generator import assemble_sequences, assemble_paired_sequences
from mock_fastq_generator.io import write_fastq

# Single-end records
records = assemble_sequences(
    template_sequence="GAATTCGCAAGCTT",
    upstream_sequence="GCCGGCCATGGCG",
    left_margin=15,
    total_length=150,
    number_of_sequences=100,
    decay_model="exponential",
    decay_rate=0.0008,
    min_val=7,
    max_val=40,
    score_type="phred",
    binned_quality=True
)

# Write to file
write_fastq("output.fastq.gz", records, compress=True)
```

---

### Standalone Web Application

`mock-fastq-generator` includes a zero-dependency, client-side web application located in `public/index.html`.

- Open `public/index.html` directly in any modern web browser.
- Interactive controls allow configuring all simulation parameters.
- Real-time visualization of positional Phred score decay curves powered by Chart.js.
- Generates and downloads FASTQ files client-side using browser Blob APIs without server computation.

---

## CLI Reference

| Argument | Default | Description |
|---|---|---|
| `--output_file` | `None` | Output file path (or prefix for paired-end mode). Mutually exclusive with `--stdout`. |
| `--stdout` | `False` | Stream FASTQ content directly to stdout. Mutually exclusive with `--output_file`. |
| `--template_sequence` | `None` | Path to FASTA template file. Defaults to bundled `template_example.fasta`. |
| `--parameters_file` | `None` | Path to JSON configuration file overriding defaults. |
| `--paired-end` | `False` | Generate R1 forward and R2 reverse complement read pairs. |
| `--gzip` | `False` | Write gzip-compressed `.fastq.gz` output files. |
| `--upstream_sequence` | `"GCCGGCCATGGCG"` | Upstream adapter/flanking sequence prepended to reads. |
| `--left_margin` | `15` | Number of random bases between upstream adapter and template. |
| `--total_length` | `500` | Total read length in bases. |
| `--number_of_sequences` | `50` | Number of reads (or read pairs) to simulate. |
| `--decay_model` | `"gaussian"` | Quality decay model (`gaussian`, `exponential`, `sigmoidal`). |
| `--decay_rate` | `0.0008` | Steepness parameter for exponential ($1.25$ exponent) or sigmoidal drop. |
| `--center` | `100` | Peak center index (Gaussian) or drop position (Sigmoidal). |
| `--std_dev` | `300` | Standard deviation for Gaussian decay curve. |
| `--min_val` | `40` | Minimum quality score floor (clamped). |
| `--max_val` | `73` | Maximum quality score ceiling. |
| `--score_type` | `"ascii"` | Input format for `min_val`/`max_val` (`"ascii"` offset 33 or `"phred"` absolute). |
| `--noise_level` | `0.5` | Standard deviation of additive Gaussian noise on Phred scores. |
| `--binned_quality` | `False` | Snap quality scores to Illumina 4-state bins (`Q2, Q12, Q23, Q37`). |
| `--homopolymer_penalty` | `False` | Subtract $10$ Phred points within homopolymer runs $>4$ bases. |

---

## Testing

Run the test suite with `pytest`:

```bash
pytest tests/
```

If `fastp` is installed on your system PATH, integration tests will automatically execute fastp quality filtering verification.

---

## License

This project is licensed under the **MIT License**. See `LICENSE` for details.

## References

1. Cock PJ, Fields CJ, Goto N, Heuer ML, Rice PM. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. *Nucleic Acids Research*. 2010;38(6):1767-1771. DOI: [10.1093/nar/gkp1137](https://doi.org/10.1093/nar/gkp1137)