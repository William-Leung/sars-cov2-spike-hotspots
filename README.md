# SARS-CoV-2 Hotspot Identification

### Workflow

First, enter the src directory with `cd src`. The first 5 steps are preparation and pre-processing. However, the necessary processed data is already in the repository (if you don't want to re-run it yourself, you can skip to the last 2 steps with the SG and SLAC methods).

### Install requirements (if necessary)
`python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt`  
You also need to install kc-align and FastTree. kalign is also a dependency for kc-align. Instructions can be found here.  
kc-align: https://github.com/davebx/kc-align  
FastTree: https://anaconda.org/bioconda/fasttree  
kalign: https://github.com/TimoLassmann/kalign  

### Fetch raw CoV-Sar-2 sequences and the spike section of the wuhan-hu-1 reference
`python3 fetch_data.py`

### Extract the spike sections from the raw sequences, align the amino acid sequences together, and convert them back to DNA.
`python3 align_sequences.py`

### Create phylogenetic tree using FastTree
`python3 create_tree.py`

### Run Fitch algorithm for ancestral state reconstruction
`python3 fitch_algorithm.py`

### Run SG Method to detect selective pressure
`python3 sg_method.py`

### Run SLAC to detect selective pressure
`python3 slac_method.py`

