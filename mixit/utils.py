"""Utility functions for the mixit package."""

import random

def get_vcf_samples(filename):
    """Extract sample names from a VCF file header."""
    with open(filename) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                return fields[9:]  # sample names
    return []  # fallback if no header found

def vcf_line_reader(vcf_path):
    with open(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                yield line