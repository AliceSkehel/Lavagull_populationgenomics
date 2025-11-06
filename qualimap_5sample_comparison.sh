#!/bin/bash
# QualiMap on 5 samples from each mapping approach

OUT_DIR=~/Sequencing_Combined/qualimap_5sample_comparison
mkdir -p $OUT_DIR/complete_mapping
mkdir -p $OUT_DIR/paleomix

# 5 complete_mapping samples (old minimap2)
echo "Running QualiMap on 5 complete_mapping samples..."
for sample in LVGU_1 LVGU_10 LVGU_20 LVGU_30 LVGU_40; do
    echo "[$(date)] QualiMap: $sample (complete_mapping)"
    qualimap bamqc \
      -bam /home/askehel/Sequencing_Combined/complete_mapping/${sample}.bam \
      -outdir $OUT_DIR/complete_mapping/$sample \
      -nt 8 \
      --java-mem-size=8G
done

# 5 PALEOMIX samples (new BWA)
echo "Running QualiMap on 5 PALEOMIX samples..."
for sample in LVGU_1 LVGU_10 LVGU_20 LVGU_30 LVGU_40; do
    echo "[$(date)] QualiMap: $sample (PALEOMIX)"
    qualimap bamqc \
      -bam /home/askehel/Sequencing_Combined/modern_remap_paleomix_Oct25/${sample}.YellowLeggedGull_ONT.bam \
      -outdir $OUT_DIR/paleomix/$sample \
      -nt 8 \
      --java-mem-size=8G
done

echo "[$(date)] QualiMap comparison complete!"
echo "Results: $OUT_DIR"
