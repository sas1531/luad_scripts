CPTAC Project Outlier Analysis and Comparison Pipeline
Examples: LUAD


#### Outlier Analysis:

Identifies aberrant expression (outliers) in genomic, phosphoproteomic, and proteomic data. 

Input file: 
- GCT file (includes meta data)
OR
- TSV file (does not include meta data) * note that row_gene and column_gene will be 1 or 2 depending on format

Output files:
- Meta data (tsv)
- Feature data (tsv)
- List of sample names (txt)
- Distribution plot (png)
- Up-regulated outlier file (tsv)
- Down-regulated outlier file (tsv)

Flags:
-f/--file		Input file path
--skip			If there are extra lines after the header, skip this number
--row_gene		Row number where features begin (after meta data) 
--column_gene		Column number where features begin (after meta data) 
--meta_column 		Column number where gene IDs or symbols are located 
-o/-out			Output file prefix for all file outputs
--aes			Color aesthetic of distribution plot
-u/--upper_x_limit	Upper x limit for distribution plot
-l/--lower_x_limit	Lower x limit for distribution plot
--tag_phospho		Tag for phospho data: write 'phospho' if the data is 				
			phosphoproteomic and 'not_phospho' if it isn't	
-b/--base_log		If data is in log10, specify here to change to log2 options: 
			c(not_log10, log10)

Example Command Line:
Rscript outlier_analysis.R \
-f ./data/proteome-ratio-norm-NArm.gct \
--skip 2 \
--row_gene 54 \
--column_gene 18 \
--meta_column 17 \
-o proteome \
--tag_phospho not_phospho



#### Cohort Comparison:

Compares two cohorts within samples and identifies significantly enriched genes in one (first) cohort against the other (second) cohort. Comparison can include one cohort (ex. Male vs Female) or multiple cohorts (ex. Western Male vs. Western Female).

Input files (all from the output of outlier analysis): 
- Meta data (tsv) * rows must be meta categories, columns are samples
- Up-regulated outlier file (tsv)
- Down-regulated outlier file (tsv)
- Gene list (txt) * optional 

Output files:
- Up-regulated significant genes & q-values (txt)
- Up-regulated sohort comparison heatmap (png)
- Down-regulated significant genes & q-values (txt)
- Down-regulated cohort comparison heatmap (png)

Flags:
--up			up-regulated data file path from outlier analysis
--down			down-regulated data file path from outlier analysis
--meta			meta data file path from outlier analysis 
-o/-out			Output file prefix for all file outputs
--comparison_a		Column name that identifies comparison A (e.g. Region, Gender, Type)
--comparison_b		Column name that identifies comparison B (e.g. Region, Gender, Type)
--comparison_c		Column name that identifies comparison C (e.g. Region, Gender, Type)
--group_1a		Group one of comparison A
--group_1b		Group one of comparison B
--group_1c		Group one of comparison C
--group_1a		Group two of comparison A
--group_1b		Group two of comparison B
--group_1c		Group two of comparison C
--group_comp		Number of groups: comparing one group to one, input 'one' or comparing
              		two groups combined, input 'two', or comparing three groups combined, input 
			'three'
--normal_tumor		If your samples include both tumor and normal samples AND you want to 
			isolate tumor samples insert 'both'. If you want to compare normal and tumor 
			samples insert 'both_normal'. If you only have tumor samples insert 'tumor_only'
--tumor_column		Name of column that indicates whether or not the sample is a tumor sample 
			or normal sample (only required if there are normal and tumor samples
--tumor_group		Name of tumor group in tumor column (only required if there are normal 
			and tumor samples
--gene_list		Path to gene list: only required if you want to isolate specific genes 
			(txt file format)

Example command line:

Rscript comparison_analysis.R \
--up ./outlier_analysis_results/luad2_prot_outlier_up.tsv \
--down ./outlier_analysis_results/luad2_prot_outlier_down.tsv \
--meta ./outlier_analysis_results/luad2_prot_meta.tsv \
-o phospho \
--comparison_a Region.of.Origin \
--group_1a Western \
--group_2a Western \
--comparison_b Gender \
--group_1b Male \
--group_2b Female \
--group_comp two \
--normal_tumor both \
--tumor_column Type \
--tumor_group Tumor \
--gene_list ./gene_drug_interactions_unique_list.txt

