import subprocess
import sys

def run_alignment():
    command = [
        "kc-align",
        "--mode", "genome",
        "--reference", "../data/wuhan_reference.fasta",
        "--sequences", "../data/raw_sequences.fasta",
        "--start", "21563",
        "--end", "25384",
    ]

    print("Running KC-Align...")
    try:
        subprocess.run(command, check=True) # raises an error if this fails
        print("Alignment complete! Output saved to out/spikes_aligned.fasta")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("Please install kc-align.")
        sys.exit(1)

if __name__ == "__main__":
    run_alignment()