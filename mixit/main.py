#!/usr/bin/env python3

"""Main entry point for the mixit package."""

import time
from .reads_simulator import ReadSimulator
from .args import prepare_args
from .utils import print_statistics



def main():
    """Main entry point."""
    # Start timing
    start_time = time.time()

    # Parse args
    args = prepare_args()

    # Generate simulated reads
    simulator = ReadSimulator(args)
    simulator.run()

    # Statistics
    print_statistics(
        total_reads=simulator.total_reads,
        reads_with_variants=simulator.reads_with_applied_variants,
        output_path=args.output,
        start_time=start_time,
    )

    print("Finished successfully.")

if __name__ == "__main__":
    main() 