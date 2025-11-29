# SARS-CoV-2 Hotspot Identification

Workflow
cd src && python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt

python3 fetch_data.py \# Fetches raw CoV-Sar-2 sequences and the spike section of the wuhan-hu-1 reference

python3 extract_spikes.py \# Extracts best aligned parts of raw sequences with the spike section

mafft --auto --thread -1 out/spikes_aa.fasta > out/spikes_aa_aligned.fasta \# Aligns the AA spike sequences together

python3 thread_dna.py \# Converts the AA back to DNA because SLAC uses DNA as input

python3 create_tree.py \# Create phylogenetic tree using FastTree

python3 fitch_algorithm.py \# Run Fitch algorithm for ancestral state reconstruction

python3 verify_fitch.py \# (Optional) Verify Fitch algorithm output

You need to install FastTree
# SARS-CoV-2 Hotspot Identification

Workflow
cd src && python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt

python3 fetch_data.py \# Fetches raw CoV-Sar-2 sequences and the spike section of the wuhan-hu-1 reference

python3 align_sequences.py \# Extract the spike sections from the raw sequences, align the amino acid sequences together, and convert them back to DNA.

python3 create_tree.py \# Create phylogenetic tree using FastTree

python3 fitch_algorithm.py \# Run Fitch algorithm for ancestral state reconstruction

python3 verify_fitch.py \# (Optional) Verify Fitch algorithm output

You need to install KC-Align, FastTree