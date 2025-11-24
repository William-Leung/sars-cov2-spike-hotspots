# cd src && python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord

RAW_SEQUENCES = "../data/raw_sequences.fasta"
WUHAN_REFERENCE = "../data/wuhan_reference.fasta"
WUHAN_SPIKE_REFERENCE = "../data/wuhan_spike_reference.fasta"
ERROR_LOG = "extraction_errors.txt"
DNA_OUTPUT_FILE = "./out/spikes_dna.fasta"
AA_OUTPUT_FILE = "./out/spikes_aa.fasta"

# Got these numbers from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
# Note that the website uses 1 based indexing but we have to convert it into 0 based indexing for Python.
REFERENCE_SPIKE_START_INDEX = 21562
REFERENCE_SPIKE_END_INDEX = 25384
SPIKE_EXPECTED_LENGTH = 3822
LENGTH_TOLERANCE = 50

def extract_spike_from_reference():
    # Extracts just the spike part of the wuhan-hu-1 reference
    record = SeqIO.read(WUHAN_REFERENCE, "fasta")

    spike_sequence = record.seq[REFERENCE_SPIKE_START_INDEX:REFERENCE_SPIKE_END_INDEX]
    sequence_length = len(spike_sequence)
    assert sequence_length == SPIKE_EXPECTED_LENGTH, "Length mismatch with SARS-CoV-2 spike sequence."

    with open(WUHAN_SPIKE_REFERENCE, "w") as output:
        output.write(f">Wuhan_Hu_1_Spike_Gene\n{spike_sequence}")

    print(f"SARS-CoV-2 Spike Sequence saved to {WUHAN_SPIKE_REFERENCE}")

def is_valid_spike(record):
    sequence_length = len(record)

    # Check that the length is within a threshold
    if not (SPIKE_EXPECTED_LENGTH - LENGTH_TOLERANCE <= sequence_length <= SPIKE_EXPECTED_LENGTH + LENGTH_TOLERANCE):
        return False, f"Invalid length {sequence_length}"
    
    # Check that the stop and start codons are correct
    sequence_str = record.seq
    # ATG in a sequence signals that it's time to start translating RNA into protein 
    if not sequence_str.startswith("ATG"): 
        return False, f"No Start Codon (Found {sequence_str[:3]})"
    # TAA,TAG,TGA signal that it's time to stop building the protein and release it
    if sequence_str[-3:] not in ["TAA", "TAG", "TGA"]:
        return False, f"No Stop Codon (Found {sequence_str[-3:]})"

    # Check that there's no premature stops 
    protein = record.seq.translate()
    if protein.count("*") > 1:
        return False, "Premature Stop Codons Detected"
    if protein.count("*") == 1 and not protein.endswith("*"):
        return False, "Stop Codon in middle of sequence"

    return True, " "

def create_aligner():
    # Create local aligner
    """Sets up a local aligner (Smith-Waterman style) optimized for DNA."""
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    return aligner

def extract_spike_from_raw_sequences():
    ref_spike = SeqIO.read(WUHAN_SPIKE_REFERENCE, "fasta")
    aligner = create_aligner()
    valid_count = 0
    total_count = 0

    extracted_dna = []
    extracted_aa = []

    print(f"Processing genomes from {RAW_SEQUENCES}...")

    with open(ERROR_LOG, "w") as err_handle:
        for genome in SeqIO.parse(RAW_SEQUENCES, "fasta"):
            if total_count == 10: # TODO early stop for testing remove this later
                break
            total_count += 1
            if total_count % 50 == 0:
                print(f"Processed {total_count} genomes...")

            # Find the section of the raw sequence that looks most like our spike reference
            alignments = aligner.align(ref_spike.seq, genome.seq)
            if not alignments:
                err_handle.write(f"{genome.id}: No alignment found\n")
                continue
            best_alignment = alignments[0]


            # Locate and cut out that best aligned sequence
            target_ranges = best_alignment.aligned[1]
            target_start = target_ranges[0][0]
            target_end = target_ranges[-1][1]
            extracted_sequence = genome.seq[target_start:target_end]
            new_record = genome[:] 
            new_record.seq = extracted_sequence
            new_record.description = f"{genome.description} | Spike DNA"


            # Validate the extracted sequence
            is_valid, message = is_valid_spike(new_record)
            if is_valid:
                extracted_dna.append(new_record)
                
                aa_seq = new_record.seq.translate()
                aa_record = SeqRecord(
                    aa_seq,
                    id=new_record.id,
                    description="Spike Protein"
                )
                extracted_aa.append(aa_record)
                
                valid_count += 1
            else:
                err_handle.write(f"{genome.id}: {message}\n")

    if extracted_dna:
        SeqIO.write(extracted_dna, DNA_OUTPUT_FILE, "fasta")
        SeqIO.write(extracted_aa, AA_OUTPUT_FILE, "fasta")
        print(f"\nDone! Extracted {valid_count}/{total_count} sequences.")
        print(f"DNA saved to: {DNA_OUTPUT_FILE}")
        print(f"Proteins saved to: {AA_OUTPUT_FILE}")


if __name__ == "__main__":
    extract_spike_from_reference()
    extract_spike_from_raw_sequences()