#!/usr/bin/env python3

"""Reads simulator module."""

import os
import time
import sys
import random
from .fasta import Fasta
from .vcf import VCFRecord
from .read import Read
from .utils import vcf_line_reader
random.seed(42)

class ReadSimulator:
    def __init__(self, args):
        self.args = args
        self.total_reads = 0
        self.total_simulate_single_read_time = 0
        self.current_vcf_record = None
        self.vcf_reader = None
        self.fasta = None
        self.fasta_iter = None
        self.n_samples = 0
        
        # Simulation data for current contig
        self.ref_contig = None
        self.seq = None
        self.contig_len = 0
        self.reads_per_contig = 0
        self.start_pos = []
        self.end_pos = []
        self.random_samples = []
        
    def initialize_simulation_data(self):
        """Initialize a new contig and simulation data."""
        try:
            self.ref_contig, self.seq = next(self.fasta_iter)
        except StopIteration:
            return False

        self.contig_len = self.fasta.get_length(self.ref_contig)
        self.reads_per_contig = int(self.args.number * self.contig_len / self.fasta.genome_length) + 2

        if self.contig_len < self.reads_per_contig:
            print(f"Warning: contig {self.ref_contig} is too short ({self.contig_len} bp) to sample {self.reads_per_contig} reads. Skipping.")
            return False

        self.start_pos = sorted(random.sample(range(1, self.contig_len), self.reads_per_contig))
        self.end_pos = [s + self.args.length for s in self.start_pos]
        self.random_samples = random.choices(range(0, self.n_samples), k=self.reads_per_contig)

        # Remove last element if end position exceeds contig length
        if self.end_pos[-1] > self.contig_len:
            self.start_pos.pop()
            self.end_pos.pop()
            self.random_samples.pop()
            self.reads_per_contig -= 1
            
            # Return False if no reads can be generated
            if self.reads_per_contig == 0:
                return False

        return True

    def collect_variants_for_read(self, read_idx):
        """Collect VCF records that fall within a read's boundaries."""
        vcf_records = []
        start_pos = self.start_pos[read_idx]
        end_pos = self.end_pos[read_idx]
        
        if self.current_vcf_record.contig != self.ref_contig:
            return vcf_records
            
        # Skip records before read start
        while self.current_vcf_record.pos < start_pos and self.current_vcf_record.contig == self.ref_contig:
            try:
                line = next(self.vcf_reader)
                self.current_vcf_record = VCFRecord(line)
            except StopIteration:
                break
        
        # Collect records within read boundaries
        while self.current_vcf_record.pos < end_pos and self.current_vcf_record.contig == self.ref_contig:
            if self.current_vcf_record.pass_quality_filter(self.args):
                vcf_records.append(self.current_vcf_record)
            try:
                line = next(self.vcf_reader)
                self.current_vcf_record = VCFRecord(line)
            except StopIteration:
                break
                
        return vcf_records

    def simulate_single_read(self, read_idx, out_file):
        """Simulate a single read and write it to the output file."""        
        vcf_records = self.collect_variants_for_read(read_idx)
        read = Read(self.seq, self.start_pos[read_idx], self.end_pos[read_idx], self.ref_contig)
        
        if vcf_records:
            read.apply_variants(vcf_records, self.random_samples[read_idx])
            
        read.write_to_file(out_file)
        self.total_reads += 1

    def run(self):
        """Run the read simulation process."""
        # Setup output file
        if os.path.exists(self.args.output):
            os.remove(self.args.output)

        # Get samples list
        samples = (self.args.vcf)
        self.n_samples = len(samples)

        # Read FASTA
        self.fasta = Fasta(self.args.fasta)
        self.fasta_iter = iter(self.fasta)

        # Init vcf reader
        self.vcf_reader = vcf_line_reader(self.args.vcf)
        
        # Read first VCF record
        try:
            line = next(self.vcf_reader)
            self.current_vcf_record = VCFRecord(line)
        except StopIteration:
            print("Error: VCF file is empty")
            sys.exit(1)

        # Generate reads for each contig
        with open(self.args.output, "w") as out_file:
            for contig in self.fasta.contigs():
                # Initialize simulation data per contig
                if not self.initialize_simulation_data():
                    continue
                
                # Iterate through the reads per contig
                for read_idx in range(self.reads_per_contig):
                    # Exit if we have already met reads number requirements
                    if self.total_reads == self.args.number:
                        print(f"Total reads generated: {self.total_reads}")
                        return
                    
                    # Simulate read
                    self.simulate_single_read(read_idx, out_file)
                    
        print(f"Total reads generated: {self.total_reads}")
