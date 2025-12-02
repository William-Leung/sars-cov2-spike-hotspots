from Bio import Phylo, SeqIO
from Bio.Data import CodonTable
from calculate_sites import construct_dictionary
from scipy.stats import binom

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]

codon_to_amino_acid = codon_table.forward_table

ancestral_tree = Phylo.read("out/ancestral_tree.tree", "newick")
sequences = {}
total_branch_length = 0
s_t = {}
n_t = {}
codon_site_counts = construct_dictionary()
for record in SeqIO.parse("out/ancestral_sequences.fasta", "fasta"):
  sequences[record.id] = str(record.seq)
  sequence_length = len(sequences[record.id])

sequence_site_counts = {}
for i in range(0, sequence_length, 3):
  for seq_id, sequence in sequences.items():
    codon = sequence[i:i+3]
    if codon in codon_site_counts:
      syn, nonsyn = codon_site_counts[codon]
      if seq_id not in sequence_site_counts:
        sequence_site_counts[seq_id] = {}
      sequence_site_counts[seq_id][i] = (syn, nonsyn)
    else:
      sequence_site_counts[seq_id][i] = (0, 0)

for i in range(0, sequence_length, 3):
  s_t[i] = 0
  n_t[i] = 0
  total_branch_length = 0
  for clade in ancestral_tree.find_clades():
    for child in clade.clades:
      branch_length = child.branch_length
      if clade.name not in sequence_site_counts or child.name not in sequence_site_counts:
        print("missing sequence site counts for", clade.name, "or", child.name)
        continue
      s_t[i] += branch_length * (sequence_site_counts[clade.name][i][0] + sequence_site_counts[child.name][i][0]) / 2
      n_t[i] += branch_length * (sequence_site_counts[clade.name][i][1] + sequence_site_counts[child.name][i][1]) / 2
      total_branch_length += branch_length

  s_t[i] /= total_branch_length
  n_t[i] /= total_branch_length

# observed mutation counts
s_c = {}
n_c = {}
for i in range(0, sequence_length, 3):
  s_c[i] = 0
  n_c[i] = 0
  for clade in ancestral_tree.find_clades():
    for child in clade.clades:
      parent_codon = sequences[clade.name][i:i+3]
      child_codon = sequences[child.name][i:i+3]
      if parent_codon not in codon_to_amino_acid or child_codon not in codon_to_amino_acid:
        continue
      if parent_codon == child_codon:
        continue
      if codon_to_amino_acid[parent_codon] == codon_to_amino_acid[child_codon]:
        s_c[i] += 1
      else:
        n_c[i] += 1

significance = 0.05 # 2-sided significance test

with open("out/sg_results.txt", "w") as output_file:
  for i in range(0, sequence_length, 3):
    if s_t[i] + n_t[i] == 0:
      continue
    
    prob_syn = s_t[i] / (s_t[i] + n_t[i])
    p = binom.cdf(s_c[i], s_c[i] + n_c[i], prob_syn)
    if p <= significance / 2:
      output_file.write(f"Positive selection, position {i}-{i+2}: s_c={s_c[i]}, n_c={n_c[i]}, s_t={s_t[i]:.2f}, n_t={n_t[i]:.2f}, p={p:.6f}\n")

    prob_nonsyn = n_t[i] / (s_t[i] + n_t[i])
    p = binom.cdf(n_c[i], s_c[i] + n_c[i], prob_nonsyn)
    if p <= significance / 2:
      output_file.write(f"Negative selection, position {i}-{i+2}: s_c={s_c[i]}, n_c={n_c[i]}, s_t={s_t[i]:.2f}, n_t={n_t[i]:.2f}, p={p:.6f}\n")