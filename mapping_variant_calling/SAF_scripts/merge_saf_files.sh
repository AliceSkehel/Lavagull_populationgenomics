> WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis

cd $WORK_DIR

echo "Merging all SAF files..."
echo "Start time: $(date)"

realSFS cat lava_gulls_saf_group1.saf.idx \
            lava_gulls_saf_group2.saf.idx \
            lava_gulls_saf_group3.saf.idx \
            lava_gulls_saf_group4a.saf.idx \
            lava_gulls_saf_group4b.saf.idx \
            lava_gulls_saf_group4c.saf.idx \
            -outnames lava_gulls_all_merged

echo "Merged SAF file created!"
echo "End time: $(date)"
