import numpy as np
from Bio import Phylo, SeqIO
from Bio.Data import CodonTable
from scipy.stats import chi2
from scipy.optimize import minimize
from calculate_sites import construct_dictionary

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
codon_to_amino_acid = codon_table.forward_table

def load_data():
    ancestral_tree = Phylo.read("out/ancestral_tree.tree", "newick")
    sequences = {}
    for record in SeqIO.parse("out/ancestral_sequences.fasta", "fasta"):
        sequences[record.id] = str(record.seq)
    return ancestral_tree, sequences

def get_codon_site_counts(sequences):
    codon_site_counts = construct_dictionary()
    sequence_length = len(next(iter(sequences.values())))
    sequence_site_counts = {}
    
    for i in range(0, sequence_length, 3):
        for seq_id, sequence in sequences.items():
            codon = sequence[i:i+3]
            if seq_id not in sequence_site_counts:
                sequence_site_counts[seq_id] = {}
            if codon in codon_site_counts:
                syn, nonsyn = codon_site_counts[codon]
                sequence_site_counts[seq_id][i] = (syn, nonsyn)
            else:
                sequence_site_counts[seq_id][i] = (0, 0)
    
    return sequence_site_counts, sequence_length

def calculate_expected_sites(tree, sequence_site_counts):
    sequence_length = max(max(sites.keys()) for sites in sequence_site_counts.values()) + 3
    s_t_total = {}
    n_t_total = {}
    s_t_avg = {}
    n_t_avg = {}
    
    for i in range(0, sequence_length, 3):
        s_t_total[i] = 0
        n_t_total[i] = 0
        total_branch_length = 0
        
        for clade in tree.find_clades():
            for child in clade.clades:
                branch_length = child.branch_length or 0.0
                if clade.name not in sequence_site_counts or child.name not in sequence_site_counts:
                    continue
                if i not in sequence_site_counts[clade.name] or i not in sequence_site_counts[child.name]:
                    continue
                
                parent_syn, parent_nonsyn = sequence_site_counts[clade.name][i]
                child_syn, child_nonsyn = sequence_site_counts[child.name][i]
                
                avg_syn = (parent_syn + child_syn) / 2
                avg_nonsyn = (parent_nonsyn + child_nonsyn) / 2
                
                s_t_total[i] += branch_length * avg_syn
                n_t_total[i] += branch_length * avg_nonsyn
                total_branch_length += branch_length
        
        if total_branch_length > 0:
            s_t_avg[i] = s_t_total[i] / total_branch_length
            n_t_avg[i] = n_t_total[i] / total_branch_length
        else:
            s_t_avg[i] = 0
            n_t_avg[i] = 0
    
    return s_t_total, n_t_total, s_t_avg, n_t_avg

def count_observed_substitutions(tree, sequences):
    sequence_length = len(next(iter(sequences.values())))
    s_c = {}
    n_c = {}
    
    for i in range(0, sequence_length, 3):
        s_c[i] = 0
        n_c[i] = 0
        
        for clade in tree.find_clades():
            for child in clade.clades:
                if clade.name not in sequences or child.name not in sequences:
                    continue
                
                parent_codon = sequences[clade.name][i:i+3]
                child_codon = sequences[child.name][i:i+3]
                
                if parent_codon not in codon_to_amino_acid or child_codon not in codon_to_amino_acid:
                    continue
                if parent_codon == child_codon:
                    continue
                
                if codon_to_amino_acid[parent_codon] == codon_to_amino_acid[child_codon]:
                    s_c[i] += 1
                else:
                    n_c[i] += 1
    
    return s_c, n_c

def calculate_total_branch_length(tree):
    total = 0.0
    for clade in tree.find_clades():
        for child in clade.clades:
            total += child.branch_length or 0.0
    return total

def estimate_ml_rates(s_c, n_c, s_t_total, n_t_total, codon_position):
    obs_s = s_c[codon_position]
    obs_n = n_c[codon_position]
    exp_s = s_t_total[codon_position]
    exp_n = n_t_total[codon_position]
    
    if exp_s + exp_n == 0:
        return None, None, None, None
    
    if obs_s + obs_n == 0:
        return 0.0, 0.0, 1.0, 0.0
    
    def neg_log_likelihood_alt(params):
        dS, dN = params
        if dS <= 0 or dN <= 0:
            return 1e10
        
        lambda_s = exp_s * dS
        lambda_n = exp_n * dN
        
        if lambda_s <= 0 or lambda_n <= 0:
            return 1e10
        
        ll = obs_s * np.log(lambda_s) - lambda_s + obs_n * np.log(lambda_n) - lambda_n
        return -ll
    
    def neg_log_likelihood_null(params):
        d = params[0]
        if d <= 0:
            return 1e10
        
        lambda_s = exp_s * d
        lambda_n = exp_n * d
        
        if lambda_s <= 0 or lambda_n <= 0:
            return 1e10
        
        ll = obs_s * np.log(lambda_s) - lambda_s + obs_n * np.log(lambda_n) - lambda_n
        return -ll
    
    initial_d = (obs_s + obs_n) / (exp_s + exp_n) if (exp_s + exp_n) > 0 else 0.01
    
    result_alt = minimize(
        neg_log_likelihood_alt,
        x0=[initial_d, initial_d],
        bounds=[(1e-8, 100), (1e-8, 100)],
        method="L-BFGS-B"
    )
    
    result_null = minimize(
        neg_log_likelihood_null,
        x0=[initial_d],
        bounds=[(1e-8, 100)],
        method="L-BFGS-B"
    )
    
    if not result_alt.success or not result_null.success:
        return None, None, None, None
    
    ll_alt = -result_alt.fun
    ll_null = -result_null.fun
    
    lrt_statistic = 2 * (ll_alt - ll_null)
    
    if lrt_statistic < 0:
        lrt_statistic = 0
    
    p_value = 1 - chi2.cdf(lrt_statistic, df=1)
    
    dS_ml, dN_ml = result_alt.x
    
    return dS_ml, dN_ml, p_value, lrt_statistic

def run_slac_analysis():
    ancestral_tree, sequences = load_data()
    sequence_site_counts, sequence_length = get_codon_site_counts(sequences)
    s_t_total, n_t_total, s_t_avg, n_t_avg = calculate_expected_sites(ancestral_tree, sequence_site_counts)
    s_c, n_c = count_observed_substitutions(ancestral_tree, sequences)
    
    significance = 0.05
    positive_selection_sites = []
    negative_selection_sites = []
    
    for i in range(0, sequence_length, 3):
        if s_t_total[i] + n_t_total[i] == 0:
            continue
        
        dS, dN, p_value, lrt = estimate_ml_rates(s_c, n_c, s_t_total, n_t_total, i)
        
        if dS is None or dN is None:
            continue
        
        dN_dS = dN / dS if dS > 0 else float("inf")
        
        if dN > dS and p_value < significance:
            positive_selection_sites.append({
                "position": f"{i}-{i+2}",
                "dN": dN,
                "dS": dS,
                "dN_dS": dN_dS,
                "s_c": s_c[i],
                "n_c": n_c[i],
                "s_t": s_t_avg[i],
                "n_t": n_t_avg[i],
                "p_value": p_value,
                "lrt": lrt
            })
        elif dS > dN and p_value < significance:
            negative_selection_sites.append({
                "position": f"{i}-{i+2}",
                "dN": dN,
                "dS": dS,
                "dN_dS": dN_dS,
                "s_c": s_c[i],
                "n_c": n_c[i],
                "s_t": s_t_avg[i],
                "n_t": n_t_avg[i],
                "p_value": p_value,
                "lrt": lrt
            })
    
    print(f"\nSLAC Analysis Results (significance level: {significance})")
    print("=" * 80)
    
    if positive_selection_sites:
        print(f"\nPositive Selection Sites ({len(positive_selection_sites)} found):")
        print("-" * 80)
        for site in positive_selection_sites:
            print(f"Position {site['position']}: dN={site['dN']:.4f}, dS={site['dS']:.4f}, "
                  f"dN/dS={site['dN_dS']:.4f}, p={site['p_value']:.6f}, LRT={site['lrt']:.4f}")
            print(f"  Observed: S={site['s_c']}, N={site['n_c']} | Expected sites: S={site['s_t']:.2f}, N={site['n_t']:.2f}")
    else:
        print("\nNo sites under positive selection detected.")
    
    if negative_selection_sites:
        print(f"\nNegative Selection Sites ({len(negative_selection_sites)} found):")
        print("-" * 80)
        for site in negative_selection_sites:
            print(f"Position {site['position']}: dN={site['dN']:.4f}, dS={site['dS']:.4f}, "
                  f"dN/dS={site['dN_dS']:.4f}, p={site['p_value']:.6f}, LRT={site['lrt']:.4f}")
            print(f"  Observed: S={site['s_c']}, N={site['n_c']} | Expected sites: S={site['s_t']:.2f}, N={site['n_t']:.2f}")
    else:
        print("\nNo sites under negative selection detected.")
    
    return positive_selection_sites, negative_selection_sites

if __name__ == "__main__":
    run_slac_analysis()
