from Bio import SeqIO, Phylo
import config
import os


def verify_tree_structure(tree_file):
    """Verify the tree file is valid and can be parsed."""
    print("1. Verifying Tree Structure")
    
    try:
        tree = Phylo.read(tree_file, "newick")
        
        terminals = tree.get_terminals()
        non_terminals = tree.get_nonterminals()
        
        print(f"✓ Tree file is valid Newick format")
        print(f"  - Terminal nodes (leaves): {len(terminals)}")
        print(f"  - Internal nodes: {len(non_terminals)}")
        print(f"  - Total nodes: {len(terminals) + len(non_terminals)}")
        
        nodes_with_names = sum(1 for node in tree.find_clades() if node.name)
        total_nodes = len(list(tree.find_clades()))
        print(f"  - Nodes with names: {nodes_with_names}/{total_nodes}")
        
        if nodes_with_names < total_nodes:
            print(f"  ⚠ Warning: Some nodes are missing names")
        
        return True, tree
        
    except Exception as e:
        print(f"✗ Error parsing tree file: {e}")
        return False, None


def verify_sequence_tree_matching(tree_file, sequence_file):
    """Verify that all nodes in tree have corresponding sequences."""
    print("\n2. Verifying Sequence-Tree Node Matching")
    
    try:
        tree = Phylo.read(tree_file, "newick")
        
        sequences = {}
        for record in SeqIO.parse(sequence_file, "fasta"):
            sequences[record.id] = str(record.seq)
        
        print(f"✓ Found {len(sequences)} sequences in FASTA file")
        
        tree_node_names = set()
        for node in tree.find_clades():
            if node.name:
                tree_node_names.add(node.name)
        
        print(f"✓ Found {len(tree_node_names)} named nodes in tree")
        
        missing_in_sequences = tree_node_names - set(sequences.keys())
        missing_in_tree = set(sequences.keys()) - tree_node_names
        
        if missing_in_sequences:
            print(f"✗ {len(missing_in_sequences)} tree nodes missing sequences:")
            for name in list(missing_in_sequences)[:10]:
                print(f"    - {name}")
            if len(missing_in_sequences) > 10:
                print(f"    ... and {len(missing_in_sequences) - 10} more")
        
        if missing_in_tree:
            print(f"✗ {len(missing_in_tree)} sequences not found in tree:")
            for name in list(missing_in_tree)[:10]:
                print(f"    - {name}")
            if len(missing_in_tree) > 10:
                print(f"    ... and {len(missing_in_tree) - 10} more")
        
        if not missing_in_sequences and not missing_in_tree:
            print("✓ All tree nodes have matching sequences")
            return True
        
        return False
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def verify_sequence_consistency(sequence_file):
    """Verify sequences are consistent (same length, valid characters)."""
    print("\n3. Verifying Sequence Consistency")
    
    try:
        sequences = {}
        lengths = []
        
        for record in SeqIO.parse(sequence_file, "fasta"):
            seq = str(record.seq)
            sequences[record.id] = seq
            lengths.append(len(seq))
        
        if not sequences:
            print("✗ No sequences found")
            return False
        
        unique_lengths = set(lengths)
        if len(unique_lengths) == 1:
            print(f"✓ All sequences have consistent length: {lengths[0]}")
        else:
            print(f"✗ Sequences have inconsistent lengths:")
            for length in sorted(unique_lengths):
                count = lengths.count(length)
                print(f"    - Length {length}: {count} sequences")
            return False
        
        valid_chars = set('ATCGN-')
        invalid_seqs = []
        
        for name, seq in sequences.items():
            seq_chars = set(seq.upper())
            invalid = seq_chars - valid_chars
            if invalid:
                invalid_seqs.append((name, invalid))
        
        if invalid_seqs:
            print(f"⚠ {len(invalid_seqs)} sequences contain invalid characters:")
            for name, chars in invalid_seqs[:5]:
                print(f"    - {name}: {chars}")
        else:
            print("✓ All sequences contain valid characters (ATCGN-)")
        
        return True
        
    except Exception as e:
        print(f"✗ Error: {e}")
        return False


def verify_parsimony_properties(tree_file, sequence_file):
    """Verify basic parsimony properties of ancestral sequences."""
    print("\n4. Verifying Parsimony Properties")
    
    try:
        tree = Phylo.read(tree_file, "newick")
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(sequence_file, "fasta")}
        
        seq_length = len(list(sequences.values())[0])
        issues = []
        
        def check_node(node, position):
            """Check if node's state is consistent with children."""
            if node.name and node.name in sequences:
                node_seq = sequences[node.name]
                if position < len(node_seq):
                    node_state = node_seq[position]
                    
                    child_states = []
                    for child in node.clades:
                        if child.name and child.name in sequences:
                            child_seq = sequences[child.name]
                            if position < len(child_seq):
                                child_states.append(child_seq[position])
                    
                    if child_states:
                        unique_child_states = set(child_states) - {'-'}
                        if unique_child_states and node_state not in {'-', ''}:
                            if node_state not in unique_child_states and len(unique_child_states) == 1:
                                issues.append((node.name, position, node_state, child_states))
        
        sample_positions = [0, seq_length // 4, seq_length // 2, 3 * seq_length // 4, seq_length - 1]
        sample_positions = [p for p in sample_positions if p < seq_length]
        
        for pos in sample_positions[:5]:
            for node in tree.find_clades():
                if not node.is_terminal():
                    check_node(node, pos)
        
        if issues:
            print(f"⚠ Found {len(issues)} potential parsimony issues (sampling {len(sample_positions)} positions)")
        else:
            print(f"✓ Sampled positions show reasonable parsimony properties")
        
        return True
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all verification checks."""
    print("\nFitch Algorithm Output Verification\n")
    
    output_dir = os.path.dirname(config.ALIGNED_CODON_FILE) or "./out"
    tree_file = os.path.join(output_dir, "slac_input_tree.tree")
    sequence_file = os.path.join(output_dir, "slac_input_sequences.fasta")
    
    if not os.path.exists(tree_file):
        print(f"✗ Tree file not found: {tree_file}")
        print("  Run: python3 fitch_algorithm.py")
        return
    
    if not os.path.exists(sequence_file):
        print(f"✗ Sequence file not found: {sequence_file}")
        print("  Run: python3 fitch_algorithm.py")
        return
    
    results = []
    
    valid, tree = verify_tree_structure(tree_file)
    results.append(("Tree Structure", valid))
    
    if valid:
        match_ok = verify_sequence_tree_matching(tree_file, sequence_file)
        results.append(("Sequence-Tree Matching", match_ok))
    
    seq_ok = verify_sequence_consistency(sequence_file)
    results.append(("Sequence Consistency", seq_ok))
    
    if valid and seq_ok:
        parsimony_ok = verify_parsimony_properties(tree_file, sequence_file)
        results.append(("Parsimony Properties", parsimony_ok))
    
    print("\nVerification Summary")
    
    for check_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {check_name}")
    
    all_passed = all(result[1] for result in results)
    
    if all_passed:
        print("\n✓ All critical checks passed!")
    else:
        print("\n⚠ Some checks failed - review output above")


if __name__ == "__main__":
    main()

