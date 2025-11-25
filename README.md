# SARS-CoV-2 Hotspot Identification

Workflow
cd src && python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt
\# Fetches raw CoV-Sar-2 sequences and the spike section of the wuhan-hu-1 reference
python3 fetch_data.py
\# Extracts best aligned parts of raw sequences with the spike section
python3 extract_spikes.py
\# Aligns the AA spike sequences together
mafft --auto --thread -1 out/spikes_aa.fasta > out/spikes_aa_aligned.fasta
\# Converts the AA back to DNA because SLAC uses DNA as input
python3 thread_dna.py
\# Create phylogenetic tree using FastTree
python3 create_tree.py
\# Run Fitch algorithm for ancestral state reconstruction
python3 fitch_algorithm.py
\# (Optional) Verify Fitch algorithm output
python3 verify_fitch.py

You need to install FastTree
