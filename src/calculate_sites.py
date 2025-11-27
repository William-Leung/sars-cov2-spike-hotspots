from Bio.Data import CodonTable

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]

codon_to_amino_acid = codon_table.forward_table
stop_codons = set(codon_table.stop_codons)

bases = ["A", "C", "G", "T"]

"""
Given a codon, returns a pair (x, y) where x is the nnumber of synonymous
sites and y is the number of nonsynonymous sites.
"""
def count_sites(codon):
  amino_acid = codon_to_amino_acid[codon]
  syn = 0
  nonsyn = 0
  for i in range(3):
    for base in bases:
      if base != codon[i]:
        new_codon = codon[:i] + base + codon[i+1:]
        if new_codon in stop_codons:
          nonsyn += 1/3
        else:
          if codon_to_amino_acid[new_codon] == amino_acid:
            syn += 1/3
          else:
            nonsyn += 1/3
  return (syn, nonsyn)

"""
Returns a dictionary which maps sequences of 3 bases to the number of synonymous
and non-synonymous sites.
"""
def construct_dictionary():
  output = {}
  for base1 in bases:
    for base2 in bases:
      for base3 in bases:
        codon = base1 + base2 + base3
        if codon not in ["TAA", "TAG", "TGA"]: # exclude termination codons
          output[codon] = count_sites(codon)
  return output
