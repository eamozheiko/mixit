# Mixit - Mix your DNA Simulator 2000

A tool for simulating DNA reads with variants from VCF files.

## Installation

```bash
pip install -r requirements.txt
pip install -e .
```

## Usage

```bash
mixit -v input.vcf -r reference.fasta -l 150 -n 10000 -o output.fasta
```

### Command Line Arguments

#### General Options
- `-v, --vcf`: Input VCF file (required)
- `-r, --reference`: Reference FASTA file (required)
- `-l, --length`: Consensus length (default: 150)
- `-n, --number`: Number of consensuses per sample (default: 10000)
- `-o, --output`: Output FASTA file (default: consensus.fasta)

#### Filtering Options
- `--min-qual`: Minimum QUAL score (Phred-scaled, default: 30)
- `--min-dp`: Minimum depth (DP, default: 10)
- `--min-mq`: Minimum mapping quality (MQ, default: 40)
- `--min-af`: Minimum allele frequency (AF, default: 0.0)

## Example

```bash
mixit -v demo/filtered.vcf -r demo/whole_pangenome.fasta -l 1000 -n 10 -o consensus.fasta
```

### Running Tests
```bash
python -m pytest tests/
```

## License

MIT License 