'''
Spike Start/End Source: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
gene            21563..25384
                     /gene="S"
                     /locus_tag="GU280_gp02"
                     /gene_synonym="spike glycoprotein"
                     /db_xref="GeneID:43740568"
'''
from Bio import SeqIO
RAW_SEQUENCES = "../data/raw_sequences.fasta"

def extra_spike_from_reference():
    input_file = "../data/wuhan_reference.fasta"
    output_file = "../data/wuhan_spike_reference.fasta"

    record = SeqIO.read(input_file, "fasta")

    start_index = 21562
    end_index = 25384

    spike_sequence = record.seq[start_index:end_index]
    sequence_length = len(spike_sequence)
    assert sequence_length == 3822, "Length mismatch with SARS-CoV-2 spike sequence."

    with open(output_file, "w") as output:
        output.write(f">Wuhan_Hu_1_Spike_Gene\n{spike_sequence}")

    print(f"SARS-CoV-2 Spike Sequence saved to {output_file}")
    

# def extract_spike(record):
#     # Extract the spike protein part of the sequence
#     # Discard sequences with lengths that are not multiples of 3 and with too many uncalled/unresolved bases
#     full_sequence = str(record.seq)
#     spike_sequence = full_sequence[SPIKE_START:SPIKE_END]

#     if len(spike_sequence) % 3 != 0:
#         return None
    
#     if spike_sequence.count('N') / len(spike_sequence) > 0.01:
#         return None

#     return spike_sequence


# def run_pipeline():
    # num_discarded = 0
    # num_kept = 0
    # # print(f"Extracting spike sequences from {SPIKE_START} to {SPIKE_END}")
    # for record in SeqIO.parse(RAW_SEQUENCES, "fasta"):
    #     dna = extract_spike(record)
    #     if dna:
    #         num_kept += 1
    #     else:
    #         num_discarded += 1
    # print(f"Kept {num_kept} sequences.")
    # print(f"Discarded {num_discarded} sequences.")


if __name__ == "__main__":
    # run_pipeline()
    extra_spike_from_reference()