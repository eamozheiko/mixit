#!/usr/bin/env python3

"""Reads simulator module."""

import os
import time
import sys
import random
from .fasta import Fasta
from .vcf import VCFRecord
from .read import Read
from .utils import vcf_line_reader, get_vcf_samples
#random.seed(42)

class ReadSimulator:
    def __init__(self, args):
        self.args = args
        self.total_reads = 0
        self.reads_with_applied_variants = 0
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
    
    def initialize_simulation_data_for_contig(self):
        """Initialize a new contig and simulation data."""
        try:
            self.ref_contig, self.seq = next(self.fasta_iter)
        except StopIteration:
            return False

        # min reads per contig. Ny a che? Tak i zhivem
        min_reads_per_contig = 2

        # Calculate reads per contig
        self.contig_len = self.fasta.get_length(self.ref_contig)
        max_reads_per_contig = int(self.args.number * self.contig_len / self.fasta.genome_length) + min_reads_per_contig
        self.reads_per_contig = min(self.contig_len, max_reads_per_contig)

        # Warning if contig is too short
        if self.contig_len < max_reads_per_contig:
            print(f"Warning: contig {self.ref_contig} is too short ({self.contig_len} bp) to sample {max_reads_per_contig} reads.")

        # Random position in contig
        self.start_pos = sorted(random.sample(range(1, self.contig_len), self.reads_per_contig))
        self.end_pos = [s + self.args.length for s in self.start_pos]
        self.random_samples = random.choices(range(0, self.n_samples), k=self.reads_per_contig)

        return True

    def collect_variants_for_read(self, read_idx):
        """Collect VCF records that fall within a read's boundaries."""
        vcf_records = []
        start_pos = self.start_pos[read_idx]
        end_pos = self.end_pos[read_idx]
        
        if self.current_vcf_record.contig != self.ref_contig:
            return vcf_records
            
        # Read VCF file row by row and Skip records before start of simulated read
        while self.current_vcf_record.pos < start_pos and self.current_vcf_record.contig == self.ref_contig:
            try:
                line = next(self.vcf_reader)
                self.current_vcf_record = VCFRecord(line)
            except StopIteration:
                break
        
        # Continue read VCF file row by row and Collect records within read boundaries
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
            self.reads_with_applied_variants += 1
            
        read.write_to_file(out_file)
        self.total_reads += 1

    def run(self):
        """Run the read simulation process."""
        # Setup output file
        if os.path.exists(self.args.output):
            os.remove(self.args.output)

        # Get samples list
        samples = get_vcf_samples(self.args.vcf)
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
                if not self.initialize_simulation_data_for_contig():
                    continue
                
                # Iterate through the reads per contig
                for read_idx in range(self.reads_per_contig):
                    # Exit if we have already met reads number requirements
                    if self.total_reads == self.args.number:
                        return
                    
                    # Simulate read
                    self.simulate_single_read(read_idx, out_file)
        

