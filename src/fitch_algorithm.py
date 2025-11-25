from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import config
import os


class FitchNode:
    """Represents a node in the tree with Fitch algorithm state information."""
    
    def __init__(self, name=None, is_leaf=False):
        self.name = name
        self.is_leaf = is_leaf
        self.children = []
        self.parent = None
        self.state_sets = {}
        self.assigned_states = {}
        self.parsimony_scores = {}
        self.sequence = None


def parse_tree_with_sequences(tree_file, sequence_file):
    """
    Parse the phylogenetic tree and align sequences to leaves.
    
    Returns:
        root: Root node of the tree
        leaf_sequences: Dictionary mapping leaf names to sequences
        sequence_length: Length of aligned sequences
    """
    tree = Phylo.read(tree_file, "newick")
    
    leaf_sequences = {}
    sequence_length = None
    
    for record in SeqIO.parse(sequence_file, "fasta"):
        seq_id = record.id.split('.')[0] if '.' in record.id else record.id
        sequence = str(record.seq)
        leaf_sequences[seq_id] = sequence
        
        if sequence_length is None:
            sequence_length = len(sequence)
        elif len(sequence) != sequence_length:
            print(f"Warning: Sequence {seq_id} has length {len(sequence)}, expected {sequence_length}")
    
    def convert_node(bio_node, parent=None):
        """Convert BioPython tree node to FitchNode."""
        is_leaf = bio_node.is_terminal()
        name = bio_node.name if bio_node.name else None
        
        fitch_node = FitchNode(name=name, is_leaf=is_leaf)
        fitch_node.parent = parent
        
        if is_leaf:
            matched = False
            if name:
                if name in leaf_sequences:
                    fitch_node.sequence = leaf_sequences[name]
                    fitch_node.name = name
                    matched = True
                else:
                    name_base = name.split(':')[0].split('.')[0]
                    for seq_id, seq in leaf_sequences.items():
                        seq_base = seq_id.split('.')[0] if '.' in seq_id else seq_id
                        if name_base == seq_base:
                            fitch_node.sequence = seq
                            fitch_node.name = seq_id
                            matched = True
                            break
                    
                    if not matched:
                        for seq_id, seq in leaf_sequences.items():
                            if name_base in seq_id or seq_id.split('.')[0] in name:
                                fitch_node.sequence = seq
                                fitch_node.name = seq_id
                                matched = True
                                break
            
            if not matched:
                print(f"Warning: Could not find sequence for leaf {name}")
        
        for child in bio_node.clades:
            child_node = convert_node(child, fitch_node)
            fitch_node.children.append(child_node)
        
        return fitch_node
    
    root = convert_node(tree.root)
    return root, leaf_sequences, sequence_length


def fitch_bottom_up(node, position):
    """
    Bottom-up pass of Fitch algorithm.
    Computes state sets and parsimony scores for each node.
    
    Returns:
        state_set: Set of possible states at this node
        parsimony_score: Minimum number of mutations required
    """
    if node.is_leaf:
        if node.sequence and position < len(node.sequence):
            state = node.sequence[position]
            if state == '-':
                return set(), 0
            node.state_sets[position] = {state}
            node.parsimony_scores[position] = 0
            return {state}, 0
        else:
            return set(), 0
    
    child_state_sets = []
    child_scores = []
    
    for child in node.children:
        state_set, score = fitch_bottom_up(child, position)
        child_state_sets.append(state_set)
        child_scores.append(score)
    
    intersection = set.intersection(*child_state_sets) if child_state_sets else set()
    
    if intersection:
        node.state_sets[position] = intersection
        node.parsimony_scores[position] = sum(child_scores)
    else:
        union = set.union(*child_state_sets) if child_state_sets else set()
        node.state_sets[position] = union
        node.parsimony_scores[position] = sum(child_scores) + 1
    
    return node.state_sets[position], node.parsimony_scores[position]


def fitch_top_down(node, position, parent_state=None):
    """
    Top-down pass of Fitch algorithm.
    Assigns specific states to each node based on parent state.
    """
    if node.is_leaf:
        if node.sequence and position < len(node.sequence):
            state = node.sequence[position]
            node.assigned_states[position] = state
            return state
        return None
    
    state_set = node.state_sets.get(position, set())
    
    if not state_set:
        node.assigned_states[position] = '-'
        return '-'
    
    if parent_state and parent_state in state_set:
        chosen_state = parent_state
    else:
        non_gap_states = [s for s in state_set if s != '-']
        if non_gap_states:
            chosen_state = non_gap_states[0]
        else:
            chosen_state = list(state_set)[0] if state_set else '-'
    
    node.assigned_states[position] = chosen_state
    
    for child in node.children:
        fitch_top_down(child, position, chosen_state)
    
    return chosen_state


def reconstruct_ancestral_sequences(root, sequence_length):
    """
    Run Fitch algorithm for all positions and reconstruct ancestral sequences.
    """
    for position in range(sequence_length):
        fitch_bottom_up(root, position)
    
    for position in range(sequence_length):
        fitch_top_down(root, position, None)
    
    def build_sequences(node):
        if node.is_leaf:
            if node.sequence:
                node.sequence = node.sequence
        else:
            sequence = ''.join([node.assigned_states.get(i, '-') 
                               for i in range(sequence_length)])
            node.sequence = sequence
        
        for child in node.children:
            build_sequences(child)
    
    build_sequences(root)


def assign_node_names(root):
    """
    Assign unique names to all internal nodes (ancestral nodes) that don't have names.
    This ensures every node can be matched with its sequence.
    """
    node_counter = [0]  # Use list to allow modification in nested function
    
    def assign_names(node):
        if not node.name:
            # Assign a unique name to internal nodes
            node.name = f"Ancestor_{node_counter[0]}"
            node_counter[0] += 1
        
        # Recursively assign to children
        for child in node.children:
            assign_names(child)
    
    assign_names(root)


def build_newick_tree(node):
    """
    Build a Newick format tree string from the FitchNode tree structure.
    Includes all nodes (observed and ancestral) with their names.
    """
    if node.is_leaf:
        name = node.name if node.name else f"Leaf_{id(node)}"
        name = name.replace(":", "_").replace(" ", "_").replace(",", "_")
        return name
    else:
        child_strings = []
        for child in node.children:
            child_strings.append(build_newick_tree(child))
        
        name = node.name if node.name else f"Ancestor_{id(node)}"
        name = name.replace(":", "_").replace(" ", "_").replace(",", "_")
        
        newick = f"({','.join(child_strings)}){name}"
        return newick


def output_results(root, output_dir=None):
    """
    Output SLAC input files (tree + sequences for Suzuki & Gojobori algorithm).
    
    Based on Suzuki & Gojobori (1999) method for detecting positive selection.
    """
    if output_dir is None:
        output_dir = os.path.dirname(config.ALIGNED_CODON_FILE) or "./out"
    
    os.makedirs(output_dir, exist_ok=True)
    assign_node_names(root)
    
    all_records = []
    
    def collect_all_sequences(node):
        """Collect sequences from all nodes (leaves and internal)."""
        if node.sequence:
            name = node.name
            if node.is_leaf:
                record = SeqRecord(Seq(node.sequence), id=name, description="Observed Sequence")
            else:
                record = SeqRecord(Seq(node.sequence), id=name, description="Ancestral Sequence")
            all_records.append(record)
        
        for child in node.children:
            collect_all_sequences(child)
    
    collect_all_sequences(root)
    
    if all_records:
        slac_sequences_file = os.path.join(output_dir, "slac_input_sequences.fasta")
        SeqIO.write(all_records, slac_sequences_file, "fasta")
        
        newick_tree = build_newick_tree(root)
        slac_tree_file = os.path.join(output_dir, "slac_input_tree.tree")
        with open(slac_tree_file, "w") as f:
            f.write(newick_tree + ";\n")
    else:
        print("Warning: No sequences to save")


def run_fitch_algorithm(sequence_file=None, tree_file=None):
    """
    Main function to run the Fitch algorithm.
    Uses codon-aligned sequences (required for SLAC dN/dS calculation).
    
    Args:
        sequence_file: Path to aligned codon sequence file (default: from config)
        tree_file: Path to phylogenetic tree file (default: from config)
    """
    if sequence_file is None:
        sequence_file = config.ALIGNED_CODON_FILE
    
    if tree_file is None:
        tree_file = config.PHYLOGENETIC_TREE_FILE
    
    # Parse tree and sequences
    root, leaf_sequences, sequence_length = parse_tree_with_sequences(tree_file, sequence_file)
    
    if sequence_length is None:
        print("Error: No sequences found or sequences have inconsistent lengths")
        return
    
    reconstruct_ancestral_sequences(root, sequence_length)
    output_results(root)
    
    return root


if __name__ == "__main__":
    run_fitch_algorithm()

