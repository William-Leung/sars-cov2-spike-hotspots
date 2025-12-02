from Bio import Entrez, SeqIO

import config

Entrez.email = "wcl53@cornell.edu"
num_sequences = 1000

def extract_spike_from_reference():
    # Extracts just the spike part of the wuhan-hu-1 reference
    record = SeqIO.read(config.WUHAN_REFERENCE_FILE, "fasta")

    spike_sequence = record.seq[config.REFERENCE_SPIKE_START_INDEX:config.REFERENCE_SPIKE_END_INDEX]
    sequence_length = len(spike_sequence)
    assert sequence_length == config.SPIKE_EXPECTED_LENGTH, "Length mismatch with SARS-CoV-2 spike sequence."

    with open(config.WUHAN_SPIKE_REFERENCE_FILE, "w") as output:
        output.write(f">Wuhan_Hu_1_Spike_Gene\n{spike_sequence}")

    print(f"SARS-CoV-2 Spike Sequence saved to {config.WUHAN_SPIKE_REFERENCE_FILE}")

def extract_raw_sequences():
    # Searches for complete SARS-CoV-2 genomes from USA 
    print(f"Fetching {num_sequences} SARS-CoV-2 sequences...")
    search_handle = Entrez.esearch(db="nucleotide", term="SARS-CoV-2[Organism] AND complete genome[Title] AND USA[Geo Location]", retmax=num_sequences)
    record = Entrez.read(search_handle)
    id_list = record["IdList"]

    # query this nucleotide db to retrieve all rows with id in our id_list (~SQL WHERE)
    with Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text") as handle:
        # streams data into our out file
        with open("../data/raw_sequences.fasta", "w") as out:
            out.write(handle.read())
    print("Data acquisition complete.")

if __name__ == "__main__":
    extract_spike_from_reference()
    extract_raw_sequences()

#kc-align --sequences ./data/raw_sequences.fasta --reference ./data/wuhan_reference.fasta --mode biological --gene S --out spike_alignment.fasta