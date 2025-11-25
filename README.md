# SARS-CoV-2 Hotspot Identification

Workflow
python3 fetch_data.py # Fetches raw CoV-Sar-2 sequences and the spike section of the wuhan-hu-1 reference
python3 extract_spikes.py # Extracts best aligned parts of raw sequences with the spike section
mafft --auto --thread -1 out/spikes_aa.fasta > out/spikes_aa_aligned.fasta # align the AA spike sequences together 