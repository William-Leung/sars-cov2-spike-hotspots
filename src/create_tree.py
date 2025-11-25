import subprocess
import config 

def build_ml_tree():
    print(f"Building Maximum Likelihood Tree from: {config.ALIGNED_CODON_FILE}")

    # Command: FastTree -gtr -nt <input> > <output>
    # -gtr: Generalized Time Reversible model (standard for DNA)
    # -nt:  Nucleotide alignment (since we are using the Codon Aligned DNA)

    cmd = ["FastTree", "-gtr", "-nt", config.ALIGNED_CODON_FILE]
    
    # Run the command and redirect output to the tree file
    with open(config.PHYLOGENETIC_TREE_FILE, "w") as outfile:
        result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        
    if result.returncode == 0:
        print(f"Finished! ML Tree saved to: {config.PHYLOGENETIC_TREE_FILE}")
    else:
        print("Error running FastTree:")
        print(result.stderr)

if __name__ == "__main__":
    build_ml_tree()