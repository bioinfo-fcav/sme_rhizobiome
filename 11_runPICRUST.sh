#!/bin/bash --login

in_dir=$1

mkdir -p $in_dir/results

echo "Running PICRUSt2 in: $in_dir"
echo ""

#conda activate picrust2

picrust2_pipeline.py \
        -s $in_dir/dna-sequences.fasta \
        -i $in_dir/feature-table.tsv \
        -o $in_dir/picrust2_out_pipeline \
        -p 12 \
        --in_traits COG,EC,KO 2> /dev/null

#conda deactivate

echo "Adding descriptions ..."

add_descriptions.py -i $in_dir/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o $in_dir/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz 2> /dev/null
add_descriptions.py -i $in_dir/picrust2_out_pipeline/COG_metagenome_out/pred_metagenome_unstrat.tsv.gz -m COG -o $in_dir/picrust2_out_pipeline/COG_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 2> /dev/null
add_descriptions.py -i $in_dir/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o $in_dir/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 2> /dev/null
add_descriptions.py -i $in_dir/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO -o $in_dir/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 2> /dev/null

zcat $in_dir/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz > $in_dir/results/PATHWAYS_pred.tsv
zcat $in_dir/picrust2_out_pipeline/COG_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz > $in_dir/results/COG_pred.tsv
zcat $in_dir/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz > $in_dir/results/EC_pred.tsv
zcat $in_dir/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz > $in_dir/results/KO_pred.tsv

echo "Done!"
