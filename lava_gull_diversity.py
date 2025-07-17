#!/usr/bin/env python3
"""
Simple Genetic Diversity Calculator for Galapagos Lava Gulls
"""

import numpy as np
import pandas as pd

def calculate_diversity_from_vcf(vcf_file):
    """
    Calculate basic diversity metrics from a VCF file
    """
    # Read VCF file (skip header lines starting with #)
    variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                variants.append(line.strip().split('\t'))
    
    # Convert to DataFrame
    df = pd.DataFrame(variants)
    genotype_columns = df.iloc[:, 9:]  # Genotype data starts at column 9
    
    n_samples = len(genotype_columns.columns)
    n_variants = len(df)
    
    print(f"Loaded {n_samples} samples with {n_variants} variants")
    
    # Calculate basic metrics
    het_counts = 0
    total_genotypes = 0
    polymorphic_sites = 0
    
    for _, row in genotype_columns.iterrows():
        site_genotypes = row.values
        
        # Count heterozygotes (genotypes like 0/1, 1/0)
        site_het = 0
        site_total = 0
        alleles = set()
        
        for gt in site_genotypes:
            if '/' in gt or '|' in gt:
                # Split genotype
                allele1, allele2 = gt.replace('|', '/').split('/')[:2]
                if allele1 != '.' and allele2 != '.':
                    alleles.add(allele1)
                    alleles.add(allele2)
                    site_total += 1
                    if allele1 != allele2:
                        site_het += 1
        
        het_counts += site_het
        total_genotypes += site_total
        
        # Check if site is polymorphic (more than 1 allele)
        if len(alleles) > 1:
            polymorphic_sites += 1
    
    # Calculate final metrics
    observed_heterozygosity = het_counts / total_genotypes if total_genotypes > 0 else 0
    proportion_polymorphic = polymorphic_sites / n_variants if n_variants > 0 else 0
    
    return {
        'n_samples': n_samples,
        'n_variants': n_variants,
        'observed_heterozygosity': observed_heterozygosity,
        'polymorphic_sites': polymorphic_sites,
        'proportion_polymorphic': proportion_polymorphic
    }

def calculate_diversity_from_matrix(genotype_matrix, sample_names=None):
    """
    Calculate diversity from a simple genotype matrix
    genotype_matrix: 2D array where rows=variants, cols=samples
    Values: 0=homozygous ref, 1=heterozygous, 2=homozygous alt
    """
    n_variants, n_samples = genotype_matrix.shape
    
    # Count heterozygotes (value = 1)
    het_count = np.sum(genotype_matrix == 1)
    total_genotypes = n_variants * n_samples
    
    # Count polymorphic sites (sites with more than one genotype type)
    polymorphic_sites = 0
    for i in range(n_variants):
        unique_genotypes = np.unique(genotype_matrix[i, :])
        if len(unique_genotypes) > 1:
            polymorphic_sites += 1
    
    observed_heterozygosity = het_count / total_genotypes
    proportion_polymorphic = polymorphic_sites / n_variants
    
    return {
        'n_samples': n_samples,
        'n_variants': n_variants,
        'observed_heterozygosity': observed_heterozygosity,
        'polymorphic_sites': polymorphic_sites,
        'proportion_polymorphic': proportion_polymorphic
    }

def print_diversity_report(results):
    """
    Print a simple diversity report
    """
    print("\n" + "="*50)
    print("LAVA GULL GENETIC DIVERSITY REPORT")
    print("="*50)
    print(f"Number of samples: {results['n_samples']}")
    print(f"Number of variants: {results['n_variants']}")
    print(f"Observed heterozygosity: {results['observed_heterozygosity']:.4f}")
    print(f"Polymorphic sites: {results['polymorphic_sites']}")
    print(f"Proportion polymorphic: {results['proportion_polymorphic']:.4f}")
    print("="*50)

# Simple usage examples
def analyze_vcf(vcf_file):
    """
    One-line analysis of VCF file
    """
    results = calculate_diversity_from_vcf(vcf_file)
    print_diversity_report(results)
    return results

def create_test_data():
    """
    Create simple test data
    """
    # Simulate 35 samples, 1000 variants
    # Low diversity typical of endangered species
    np.random.seed(42)
    
    # Most sites are monomorphic (0), some heterozygous (1), few homozygous alt (2)
    genotype_matrix = np.random.choice([0, 1, 2], 
                                     size=(1000, 35), 
                                     p=[0.95, 0.04, 0.01])  # Low diversity
    
    return genotype_matrix

# Main usage
if __name__ == "__main__":
    print("Simple Lava Gull Genetic Diversity Analysis")
    
    # Test with simulated data
    print("\nTesting with simulated data:")
    test_data = create_test_data()
    results = calculate_diversity_from_matrix(test_data)
    print_diversity_report(results)
    
    print("\nTo use with your VCF file:")
    print("results = analyze_vcf('your_file.vcf')")
