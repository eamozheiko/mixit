"""Read simulation module."""

class Read:
    def __init__(self, full_sequence, start, end, contig):
        self.contig = contig
        self.start = start
        self.end = end
        self.sequence = full_sequence[start:end]  # Extract subsequence

    def apply_variants(self, vcf_records, sample_name):
        """Apply VCF variants to the read sequence."""
        seq = list(self.sequence)
        offset = self.start  # Genomic coordinate of seq[0]
        pos_shift = 0        # Track shift in read index due to insertions/deletions

        for record in vcf_records:
            if not record.has_variant(sample_name):
                continue
            # VCF position is 1-based, convert to 0-based
            in_read_pos = record.pos - 1 - offset + pos_shift

            # Apply the variant (SNP or INDEL)
            seq[in_read_pos:in_read_pos + len(record.ref)] = list(record.alt)

            # Update position shift for downstream variants
            alt_len = len(record.alt)
            ref_len = len(record.ref)
            pos_shift += alt_len - ref_len

        self.sequence = "".join(seq)

    def write_to_file(self, out_file):
        """Write the read in FASTA format to an already open file handle"""
        out_file.write(f">{self.contig}:{self.start}:{self.end}\n{self.sequence}\n")