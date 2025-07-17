#!/bin/bash

echo "Determining sex for all 35 LVGU samples..."
echo "Sample Z_Coverage Auto_Coverage Ratio Predicted_Sex" > sex_determination_results.txt

# Get all sample names from the VCF
bcftools query -l ~/Sequencing_Round2/variant_calling/variants.vcf.gz | while read sample; do
    
    bam_file="/home/askehel/Sequencing_Round2/complete_mapping/${sample}.bam"
    
    # Get Z chromosome coverage (OZ118750.1)
    z_cov=$(samtools coverage $bam_file | awk '$1=="OZ118750.1" {printf "%.2f", $7}')
    
    # Get average autosomal coverage (chromosomes 1-5)
    auto_cov=$(samtools coverage $bam_file | awk '
        $1~/OZ118746.1|OZ118747.1|OZ118748.1|OZ118749.1|OZ118751.1/ {
            sum+=$7; count++
        } 
        END {printf "%.2f", sum/count}')
    
    # Calculate ratio
    ratio=$(echo "$z_cov $auto_cov" | awk '{printf "%.2f", $1/$2}')
    
    # Predict sex (Males ZZ should have ratio ~1.0, Females ZW should have ratio ~0.5)
    sex=$(echo "$ratio" | awk '{if($1 > 0.75) print "Male"; else print "Female"}')
    
    echo "$sample ${z_cov}X ${auto_cov}X $ratio $sex" >> sex_determination_results.txt
    echo "Processed $sample: $sex (ratio: $ratio)"
    
done

echo ""
echo "Sex determination complete! Results saved to: sex_determination_results.txt"
echo ""
echo "Summary:"
cat sex_determination_results.txt | tail -n +2 | awk '{print $4}' | sort | uniq -c
