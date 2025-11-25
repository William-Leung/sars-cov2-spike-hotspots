from Bio import SeqIO
from Bio.Seq import Seq
import config

def thread_dna():
    dna_dict = SeqIO.to_dict(SeqIO.parse(config.UNALIGNED_DNA_FILE, "fasta"))
    codon_aligned_records = []

    print("Threading DNA onto Protein Alignment...")
    
    for protein_record in SeqIO.parse(config.ALIGNED_AA_FILE, "fasta"):
        if protein_record.id not in dna_dict:
            print(f"{protein_record.id} found in protein but not in DNA file.")
            continue
            
        dna_seq = dna_dict[protein_record.id].seq
        aligned_aa = str(protein_record.seq)
        
        aligned_dna = ""
        dna_idx = 0
        
        for aa in aligned_aa:
            if aa == "-":
                # if the protein has a gap, insert a gap codon into DNA
                aligned_dna += "---"
            else:
                # if protein has a letter, then insert the next 3 letters
                codon = dna_seq[dna_idx : dna_idx + 3]
                aligned_dna += str(codon)
                dna_idx += 3
                
        new_record = protein_record[:] # copy ID and description
        new_record.seq = Seq(aligned_dna)
        codon_aligned_records.append(new_record)

    SeqIO.write(codon_aligned_records, config.ALIGNED_CODON_FILE, "fasta")
    print(f"Finished! Saved codon alignment to {config.ALIGNED_CODON_FILE}")

if __name__ == "__main__":
    thread_dna()