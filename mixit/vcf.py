"""VCF file handling module."""

class VCFRecord:
    def __init__(self, line):
        fields = line.strip().split("\t")
        self.contig = fields[0]
        self.pos = int(fields[1])
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4].split(",")[0]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = fields[7]
        self.fmt = fields[8]
        self.sample_data = fields[9:] if len(fields) > 9 else []
        self.info_dict = dict(kv.split("=") for kv in self.info.split(";") if "=" in kv)

    def get_gt(self, sample_index):
        """Get genotype string (e.g., '0/1') for a given sample index."""
        return self.sample_data[sample_index].split(":")[0]

    def has_variant(self, sample_index):
        """Check for the presence of a variant in a specific sample."""
        gt = self.get_gt(sample_index)
        if "." in gt or gt == "0/0":
            return False
        return True

    def pass_quality_filter(self, args):
        """Filter VCF record based on quality metrics."""
        try:
            if float(self.qual) < args.min_qual:
                return False

            dp = float(self.info_dict.get("DP", 0))
            if dp < args.min_dp:
                return False

            mq = float(self.info_dict.get("MQ", 0))
            if mq < args.min_mq:
                return False

            if "AF" not in self.info_dict:
                # Compute allele frequency if AF not available
                ac = float(self.info_dict.get("AC", 0))
                an = float(self.info_dict.get("AN", 0))
                af = ac / an if an > 0 else 0.0
            else:
                af = float(self.info_dict.get("AF", 0))

            if af < args.min_af:
                return False

        except ValueError:
            # Handles cases where fields like QUAL or DP are malformed
            return False

        return True 