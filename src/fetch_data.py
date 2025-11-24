from Bio import Entrez, SeqIO

Entrez.email = "wcl53@cornell.edu"
num_sequences = 1000

def fetch_data():
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
    fetch_data()

#kc-align --sequences ./data/raw_sequences.fasta --reference ./data/wuhan_reference.fasta --mode biological --gene S --out spike_alignment.fasta