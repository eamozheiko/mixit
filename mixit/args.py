"""Command line argument parsing module."""

import argparse as ap
import sys

class SilentArgumentParser(ap.ArgumentParser):
    def error(self, message):
        # Print only your custom help message, not usage
        self.print_help(sys.stderr)
        print(f"\nError: {message}\n", file=sys.stderr)
        sys.exit(2)

def prepare_args():
    """Prepare optparser object with grouped arguments."""
    
    description = "Mixit - Mix your DNA Simulator 2000"
    epilog = "Example: %(prog)s -v demo/filtered.vcf -r demo/whole_pangenome.fasta -l 150 -n 10000 -o consensus.fasta"

    argparser = SilentArgumentParser(description=description, epilog=epilog)
    Mixit_VERSION = '1.0'

    # General Options
    general = argparser.add_argument_group("General Options")
    general.add_argument("--version", action="version", version="%(prog)s " + Mixit_VERSION)
    general.add_argument("-v", "--vcf", dest="vcf", type=str, required=True, help="Input VCF file")
    general.add_argument("-r", "--reference", dest="fasta", required=True, help="Reference FASTA file")
    general.add_argument("-l", "--length", dest="length", type=int, default=150, help="Consensus length")
    general.add_argument("-n", "--number", dest="number", type=int, default=10000, help="Number of consensuses per sample")
    general.add_argument("-o", "--output", dest="output", default="consensus.fasta", help="Output FASTA file")

    # Filtering Options
    filters = argparser.add_argument_group("Filtering Options")

    # Quality Filters
    filters.add_argument("--min-qual", type=float, default=30.0, help="Minimum QUAL score (Phred-scaled, default: 30)")
    filters.add_argument("--min-dp", type=int, default=10, help="Minimum depth (DP, default: 10)")
    filters.add_argument("--min-mq", type=float, default=40.0, help="Minimum mapping quality (MQ, default: 40)")
    filters.add_argument("--min-af", type=float, default=0.0, help="Minimum allele frequency (AF, default: 0.0)")
    args = argparser.parse_args()

    return args 