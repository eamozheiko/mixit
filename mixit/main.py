#!/usr/bin/env python3

"""Main entry point for the mixit package."""

import time
from .reads_simulator import ReadSimulator
from .args import prepare_args



def main():
    """Main entry point."""
    # Start timing
    start_time = time.time()

    # Parse args
    args = prepare_args()

    # Generate simulated reads
    simulator = ReadSimulator(args)
    simulator.run()
    print(f"Finished successfully. Please find you output fasta in {args.output}")

    # Output elapsed time
    elapsed_time = time.time() - start_time
    print(f"Total execution time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main() 