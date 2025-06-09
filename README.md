# Mixit

A tool for simulating DNA reads with variants from VCF files.

## Installation

### Manual Istallation
```bash
git clone https://github.com/eamozheiko/mixit.git
cd mixit
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

Citation
--------

If you use Mixit in your research, the most relevant link to cite is:

* https://github.com/eamozheiko/mixit

### Running Tests
```bash
python -m pytest tests/
```

## License

MIT License 
