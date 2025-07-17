#!/bin/bash

import pysam
import numpy as np

def determine_sex(bam_file, ref_genome, sex_determination_system="XY"):
    # Calculate coverage across genome
    autosomal_coverage = calculate_autosomal_coverage(bam_file, ref_genome)
    x_coverage = calculate_x_coverage(bam_file, ref_genome)
    y_or_w_coverage = calculate_y_or_w_coverage(bam_file, ref_genome)
    
    # Calculate ratios
    x_to_autosomal_ratio = x_coverage / autosomal_coverage
    
    # Determine sex based on system and ratios
    if sex_determination_system == "XY":
        if x_to_autosomal_ratio < 0.7:  # Threshold for XY (male)
            return "male"
        else:
            return "female"
    elif sex_determination_system == "ZW":
        # Inverse logic for ZW systems
        if x_to_autosomal_ratio > 0.7:  # Threshold for ZZ (male)
            return "male"
        else:
            return "female"

