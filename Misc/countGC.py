from Bio import SeqIO

GC = 0 
total_length = 0
for seq_record in SeqIO.parse("lobstah.ests.fa", "fasta"):
    Gs = seq_record.seq.count("G")
    Cs = seq_record.seq.count("C")
    total = len(seq_record)
    GC += Gs
    GC += Cs
    total_length += total
  
print GC
print total_length
print GC/total_length
print float(GC)/float(total_length)