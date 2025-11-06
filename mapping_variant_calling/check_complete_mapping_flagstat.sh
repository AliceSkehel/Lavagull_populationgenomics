#!/bin/bash

echo "Running flagstat on all complete_mapping BAMs..."
echo "Sample,Total_Reads,Mapped,Mapped_%,Duplicates,Properly_Paired,Properly_Paired_%"

for bam in /home/askehel/Sequencing_Combined/complete_mapping/LVGU_{1,2,3,4,5,6,7,9,10,11,12,13,14,15,18,19,20,23,24,25,30,34,37,39,40,42,46,47,48,52,54,56}.bam; do
    sample=$(basename $bam .bam)
    
    # Get flagstat output
    stats=$(samtools flagstat $bam)
    
    total=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped=$(echo "$stats" | grep "mapped (" | head -1 | awk '{print $1}')
    mapped_pct=$(echo "$stats" | grep "mapped (" | head -1 | awk '{print $5}' | tr -d '(%:')
    dups=$(echo "$stats" | grep "duplicates" | head -1 | awk '{print $1}')
    paired=$(echo "$stats" | grep "properly paired" | awk '{print $1}')
    paired_pct=$(echo "$stats" | grep "properly paired" | awk '{print $5}' | tr -d '(%:')
    
    echo "$sample,$total,$mapped,$mapped_pct,$dups,$paired,$paired_pct"
done
