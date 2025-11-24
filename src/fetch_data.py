from Bio import Entrez, SeqIO

Entrez.email = "wcl53@cornell.edu"

def fetch_data():
    # Searches for complete SARS-CoV-2 genomes from USA 
    print("Fetching 200 SARS-CoV-2 sequences...")
    search_handle = Entrez.esearch(db="nucleotide", term="SARS-CoV-2[Organism] AND complete genome[Title] AND USA[Geo Location]", retmax=200)
    record = Entrez.read(search_handle)
    id_list = record["IdList"]

    # query this nucleotide db to retrieve all rows with id in our id_list (~SQL WHERE)
    with Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text") as handle:
        # streams data into our out file
        with open("raw_sequences.fasta", "w") as out:
            out.write(handle.read())
    print("Data acquisition complete.")

if __name__ == "__main__":
    fetch_data()