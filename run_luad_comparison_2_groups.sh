#!/bin/bash

cat ./submitter_run_luad.txt | while read -r -a myArray
do
name="${myArray[0]}"
comparison_a="${myArray[1]}"
group1a="${myArray[2]}"
group2a="${myArray[3]}"
comparison_b="${myArray[4]}"
group1b="${myArray[5]}"
group2b="${myArray[6]}"

Rscript luad-v2.0_comparison_pipeline.R \
--up ./outlier_analysis_results/luad2_phos_outlier_up.tsv \
--down ./outlier_analysis_results/luad2_phos_outlier_down.tsv \
--meta ./outlier_analysis_results/luad2_phos_meta.tsv \
-o luad2_phos_${name} \
--comparison_a ${comparison_a} \
--group_1a ${group1a} \
--group_2a ${group2a} \
--comparison_a ${comparison_b} \
--group_1a ${group1b} \
--group_2a ${group2b} \
--group_comp two \
--normal_tumor both \
--tumor_column Type \
--tumor_group Tumor

Rscript luad-v2.0_comparison_pipeline.R \
--up ./outlier_analysis_results/luad2_prot_outlier_up.tsv \
--down ./outlier_analysis_results/luad2_prot_outlier_down.tsv \
--meta ./outlier_analysis_results/luad2_prot_meta.tsv \
-o luad2_prot_${name} \
--comparison_a ${comparison_a} \
--group_1a ${group1a} \
--group_2a ${group2a} \
--comparison_a ${comparison_b} \
--group_1a ${group1b} \
--group_2a ${group2b} \
--group_comp two \
--normal_tumor both \
--tumor_column Type \
--tumor_group Tumor

Rscript luad-v2.0_comparison_pipeline.R \
--up ./outlier_analysis_results/luad2_rna_log2_outlier_up.tsv \
--down ./outlier_analysis_results/luad2_rna_log2_outlier_down.tsv \
--meta ./outlier_analysis_results/luad2_rna_log2_meta.tsv \
-o luad2_phos_test_${name} \
--comparison_a ${comparison_a} \
--group_1a ${group1a} \
--group_2a ${group2a} \
--comparison_a ${comparison_b} \
--group_1a ${group1b} \
--group_2a ${group2b} \
--group_comp two \
--normal_tumor both \
--tumor_column Type \
--tumor_group Tumor

Rscript luad-v2.0_comparison_pipeline.R \
--up ./outlier_analysis_results/luad2_phos_outlier_up.tsv \
--down ./outlier_analysis_results/luad2_phos_outlier_down.tsv \
--meta ./outlier_analysis_results/luad2_phos_meta.tsv \
-o luad2_phos_test_${name} \
--comparison_a ${comparison_a} \
--group_1a ${group1a} \
--group_2a ${group2a} \
--comparison_a ${comparison_b} \
--group_1a ${group1b} \
--group_2a ${group2b} \
--group_comp two \
--normal_tumor both \
--tumor_column Type \
--tumor_group Tumor ; done
