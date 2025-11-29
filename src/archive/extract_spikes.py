from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
import config

def is_valid_spike(record):
    sequence_length = len(record)

    # Check that the length is within a threshold
    if not (config.SPIKE_EXPECTED_LENGTH - config.LENGTH_TOLERANCE <= sequence_length <= config.SPIKE_EXPECTED_LENGTH + config.LENGTH_TOLERANCE):
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
    ref_spike = SeqIO.read(config.WUHAN_SPIKE_REFERENCE_FILE, "fasta")
    aligner = create_aligner()
    valid_count = 0
    total_count = 0

    extracted_dna = []
    extracted_aa = []

    print(f"Processing genomes from {config.RAW_SEQUENCES_FILE}...")

    with open(config.ERROR_LOG_FILE, "w") as err_handle:
        for genome in SeqIO.parse(config.RAW_SEQUENCES_FILE, "fasta"):
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
        SeqIO.write(extracted_dna, config.UNALIGNED_DNA_FILE, "fasta")
        SeqIO.write(extracted_aa, config.UNALIGNED_AA_FILE, "fasta")
        print(f"\nDone! Extracted {valid_count}/{total_count} sequences.")
        print(f"DNA saved to: {config.UNALIGNED_DNA_FILE}")
        print(f"Proteins saved to: {config.UNALIGNED_AA_FILE}")


if __name__ == "__main__":
    extract_spike_from_raw_sequences()