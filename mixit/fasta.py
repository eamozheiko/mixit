"""FASTA file handling module."""

class Fasta:
    def __init__(self, filename):
        self.sequences = {}
        self._load_fasta(filename)
        self.genome_length = self.get_genome_length()

    def _load_fasta(self, filename):
        with open(filename) as f:
            header = None
            seq_lines = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        self.sequences[header] = ''.join(seq_lines)
                    header = line[1:].split()[0]  # first word of header
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header:
                self.sequences[header] = ''.join(seq_lines)

    def get_sequence(self, contig):
        return self.sequences.get(contig, "")

    def get_length(self, contig):
        return len(self.sequences.get(contig, ""))

    def get_genome_length(self):
        return sum(len(seq) for seq in self.sequences.values())

    def get_total_length(self):
        return sum(len(seq) for seq in self.sequences.values())

    def __iter__(self):
        return iter(self.sequences.items())

    def __getitem__(self, contig):
        return self.sequences[contig]

    def __contains__(self, contig):
        return contig in self.sequences

    def contigs(self):
        return list(self.sequences.keys()) 