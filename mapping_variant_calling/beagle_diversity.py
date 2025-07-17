#!/usr/bin/env python3

import numpy as np
import gzip
import sys

def calculate_diversity_from_beagle(beagle_file, output_file):
    """
    Calculate genetic diversity metrics from ANGSD beagle file
    """
    print(f"Reading beagle file: {beagle_file}")
    
    # Read the beagle file
    with gzip.open(beagle_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        
        # Extract sample names (every 3rd column starting from index 3)
        samples = []
        for i in range(3, len(header), 3):
            sample_name = header[i].replace('_0', '')
            samples.append(sample_name)
        
        n_samples = len(samples)
        print(f"Found {n_samples} samples")
        
        # Initialize variables for calculations
        pi_sum = 0  # Nucleotide diversity
        theta_w_sum = 0  # Watterson's theta
        sites_processed = 0
        segregating_sites = 0
        
        for line_num, line in enumerate(f):
            if line_num % 50000 == 0:
                print(f"Processing site {line_num}")
            
            parts = line.strip().split('\t')
            
            # Extract genotype probabilities for all samples
            allele_freqs = []
            for i in range(3, len(parts), 3):
                # Get probabilities for this sample (0/0, 0/1, 1/1)
                p00 = float(parts[i])
                p01 = float(parts[i+1]) 
                p11 = float(parts[i+2])
                
                # Calculate expected allele frequency for this individual
                af = p01 * 0.5 + p11  # 0.5 * heterozygote + 1.0 * homozygote alt
                allele_freqs.append(af)
            
            # Population allele frequency
            p = np.mean(allele_freqs)
            
            # Only use polymorphic sites with reasonable frequency
            if 0.01 <= p <= 0.99:
                # Nucleotide diversity (π) - average pairwise differences
                pi_site = 2 * p * (1 - p)
                pi_sum += pi_site
                
                # Count as segregating site for Watterson's theta
                segregating_sites += 1
                
                sites_processed += 1
    
    print(f"Processed {sites_processed} polymorphic sites")
    print(f"Segregating sites: {segregating_sites}")
    
    # Calculate final statistics
    pi = pi_sum / sites_processed if sites_processed > 0 else 0
    
    # Watterson's theta calculation
    # θw = S / Σ(1/i) where S = segregating sites, sum from i=1 to n-1
    harmonic_number = sum(1/i for i in range(1, n_samples))
    theta_w = segregating_sites / harmonic_number / sites_processed if sites_processed > 0 else 0
    
    # Tajima's D
    tajima_d = (pi - theta_w) / np.sqrt(calculate_tajima_variance(n_samples, segregating_sites, sites_processed)) if sites_processed > 0 else 0
    
    # Estimate effective population size
    # Ne = θ / (4 * μ), using typical avian mutation rate
    mu = 2.3e-9  # per bp per generation
    generation_time = 10  # years (typical for gulls)
    
    Ne_pi = pi / (4 * mu)
    Ne_theta = theta_w / (4 * mu)
    
    # Save comprehensive results
    with open(f"{output_file}_diversity.txt", 'w') as f:
        f.write("LAVA GULL GENETIC DIVERSITY ANALYSIS\n")
        f.write("=" * 60 + "\n")
        f.write(f"Sample size: {n_samples} individuals\n")
        f.write(f"Sites analyzed: {sites_processed:,}\n")
        f.write(f"Segregating sites: {segregating_sites:,}\n")
        f.write(f"Proportion polymorphic: {segregating_sites/sites_processed:.4f}\n\n")
        
        f.write("DIVERSITY METRICS:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Nucleotide diversity (π): {pi:.6f}\n")
        f.write(f"Watterson's theta (θw): {theta_w:.6f}\n")
        f.write(f"Tajima's D: {tajima_d:.4f}\n\n")
        
        f.write("EFFECTIVE POPULATION SIZE ESTIMATES:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Ne (from π): {Ne_pi:,.0f}\n")
        f.write(f"Ne (from θw): {Ne_theta:,.0f}\n")
        f.write(f"Mean Ne: {(Ne_pi + Ne_theta)/2:,.0f}\n\n")
        
        f.write("INTERPRETATION:\n")
        f.write("-" * 20 + "\n")
        if tajima_d > 0.5:
            f.write("Tajima's D > 0: Suggests population structure or balancing selection\n")
        elif tajima_d < -0.5:
            f.write("Tajima's D < 0: Suggests population expansion or directional selection\n")
        else:
            f.write("Tajima's D ≈ 0: Suggests neutral evolution at mutation-drift equilibrium\n")
        
        f.write(f"\nGenetic diversity level: ")
        if pi < 0.001:
            f.write("Low (π < 0.001)\n")
        elif pi < 0.005:
            f.write("Moderate (0.001 ≤ π < 0.005)\n")
        else:
            f.write("High (π ≥ 0.005)\n")
        
        f.write(f"\nEffective population size: ")
        mean_ne = (Ne_pi + Ne_theta) / 2
        if mean_ne < 1000:
            f.write("Very small (Ne < 1,000) - Conservation concern\n")
        elif mean_ne < 10000:
            f.write("Small (1,000 ≤ Ne < 10,000) - Monitor closely\n")
        elif mean_ne < 50000:
            f.write("Moderate (10,000 ≤ Ne < 50,000) - Stable population\n")
        else:
            f.write("Large (Ne ≥ 50,000) - Healthy population\n")
        
        f.write(f"\nASSUMPTIONS:\n")
        f.write(f"- Mutation rate: {mu} per bp per generation\n")
        f.write(f"- Generation time: {generation_time} years\n")
        f.write(f"- Neutral evolution model\n")
    
    print(f"\nRESULTS SUMMARY:")
    print(f"Nucleotide diversity (π): {pi:.6f}")
    print(f"Watterson's theta (θw): {theta_w:.6f}")
    print(f"Tajima's D: {tajima_d:.4f}")
    print(f"Effective population size: {(Ne_pi + Ne_theta)/2:,.0f}")
    print(f"\nDetailed results saved to: {output_file}_diversity.txt")

def calculate_tajima_variance(n, S, L):
    """Calculate variance for Tajima's D statistic"""
    a1 = sum(1/i for i in range(1, n))
    a2 = sum(1/(i**2) for i in range(1, n))
    
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    
    c1 = b1 - 1/a1
    c2 = b2 - (n + 2)/(a1 * n) + a2/(a1**2)
    
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    
    return (e1 * S + e2 * S * (S - 1)) / L

if __name__ == "__main__":
    beagle_file = "/home/askehel/Sequencing_Combined/angsd_variant_calling/lava_gulls_illumina_only.beagle.gz"
    output_file = "lava_gulls"
    
    calculate_diversity_from_beagle(beagle_file, output_file)
