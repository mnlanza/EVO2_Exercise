from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# === Input and output filenames ===
input_fasta = "codon83_variants.fasta"  # <-- change to your input file
output_fasta = "reverse_complements_83.fasta"

# === Read input, reverse complement each sequence ===
records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    revcomp_seq = record.seq.reverse_complement()
    revcomp_record = SeqRecord(revcomp_seq, id=record.id + "_rc", description="")
    records.append(revcomp_record)

# === Write to new FASTA ===
SeqIO.write(records, output_fasta, "fasta")

print(f"Wrote {len(records)} reverse complement sequences to {output_fasta}")
