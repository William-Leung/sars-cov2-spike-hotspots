RAW_SEQUENCES = "../data/raw_sequences.fasta"
WUHAN_REFERENCE = "../data/wuhan_reference.fasta"
WUHAN_SPIKE_REFERENCE = "../data/wuhan_spike_reference.fasta"
ERROR_LOG = "./out/extraction_errors.txt"
DNA_OUTPUT_FILE = "./out/spikes_dna.fasta"
AA_OUTPUT_FILE = "./out/spikes_aa.fasta"

# Got these numbers from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
# Note that the website uses 1 based indexing but we have to convert it into 0 based indexing for Python.
REFERENCE_SPIKE_START_INDEX = 21562
REFERENCE_SPIKE_END_INDEX = 25384
SPIKE_EXPECTED_LENGTH = 3822
LENGTH_TOLERANCE = 50