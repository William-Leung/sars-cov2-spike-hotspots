RAW_SEQUENCES_FILE = "../data/raw_sequences.fasta"
WUHAN_REFERENCE_FILE = "../data/wuhan_reference.fasta"
WUHAN_SPIKE_REFERENCE_FILE = "../data/wuhan_spike_reference.fasta"
ERROR_LOG_FILE = "./out/extraction_errors.txt"
UNALIGNED_DNA_FILE = "./out/spikes_dna.fasta"
UNALIGNED_AA_FILE = "./out/spikes_aa.fasta"
ALIGNED_AA_FILE = "./out/spikes_aa_aligned.fasta"
ALIGNED_CODON_FILE = "./out/spikes_codon_aligned.fasta"
PHYLOGENETIC_TREE_FILE = "./out/phylogenetic_tree.tree"

# Got these numbers from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
# Note that the website uses 1 based indexing but we have to convert it into 0 based indexing for Python.
REFERENCE_SPIKE_START_INDEX = 21562
REFERENCE_SPIKE_END_INDEX = 25384
SPIKE_EXPECTED_LENGTH = 3822
LENGTH_TOLERANCE = 50