"""Utility functions for the mixit package."""

import time
from typing import Dict, Union, Optional

def print_statistics(
    total_reads: int,
    reads_with_variants: int,
    output_path: str,
    start_time: Optional[float] = None,
) -> None:
    """Prints stats."""
    # Calculate derived metrics
    stats = {
        "Total reads generated": total_reads,
        "Reads with applied variants": reads_with_variants,
        "Variant application rate": f"{(reads_with_variants / total_reads) * 100:.2f}%"
    }

    # Header
    print("\n" + "=" * 50)
    print(" SIMULATION STATISTICS ".center(50, "="))
    
    # Statistics rows (aligned)
    max_label_length = max(len(label) for label in stats.keys())
    for label, value in stats.items():
        print(f"• {label.ljust(max_label_length)} : {str(value).rjust(10)}")

    # Footer with execution time and output path
    print("-" * 50)
    if start_time:
        print(f"Execution time: {time.time() - start_time:.2f}s")
    print(f"Output: {output_path}")
    print("=" * 50 + "\n")

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