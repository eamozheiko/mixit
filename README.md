# Mixit

A bioinformatics tool that generates simulated DNA sequencing reads containing variants from a multi-sample VCF file. The tool produces a specified number of reads with fixed length, randomly distributed across a provided reference genome (FASTA format).

Key Features:
- Processes VCF files containing multiple samples
- Simulates reads incorporating variants (SNVs, indels)
- Works with standard reference genomes (FASTA)
- Configurable read length, number of reads and min allele fraction
- Random positional distribution across genome (without repeats)

## Installation

### Manual Istallation
```bash
git clone https://github.com/eamozheiko/mixit.git
cd mixit
pip install -e .
```

## Usage

```bash
mixit -v demo/filtered.vcf -r demo/whole_pangenome.fasta -l 150 -n 10000 -o consensus.fasta
```

### Command Line Arguments

#### General Options
- `-v, --vcf`: Input VCF file (required)
- `-r, --reference`: Reference FASTA file (required)
- `-l, --length`: Simulated read's length (default: 150)
- `-n, --number`: Nessesary number of reads to simulate (default: 10000)
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
