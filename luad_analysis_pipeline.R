#!/Applications/anaconda3/pkgs/r-base-3.4.3-h9b62496_0/bin/Rscript

library(tidyverse)
library(dplyr)
library(ggplot2)
library(optparse, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(gridExtra)

options(warn=1)

# Set up command line
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--row_gene"), type="integer", default=1, 
              help="index number(s) of row where gene values begins", metavar="character"),
  make_option(c("-c", "--column_gene"), type="integer", default=1, 
              help="index number(s) of column where gene values begins", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="gene_out.tsv", 
              help="output file name for geneID/samples dataframe", metavar="character"),
  make_option(c("-s", "--sample_names"), type="character", default="names_out.txt", 
              help="output file name of list of samples", metavar="character"),
  make_option(c("-p", "--plus_outliers"), type="character", default="outlier_plus.tsv", 
              help="output file name of plus outliers", metavar="character"),
  make_option(c("-m", "--minus_outliers"), type="character", default="outlier_minus.tsv", 
              help="output file name of minus outliers", metavar="character"),
  make_option(c("-b", "--both_outliers"), type="character", default="outlier_both.tsv", 
              help="output file name of plus", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# test script:
# Rscript luad_analysis_pipeline_test.R -f luad-v1.0-phosphodfeome-ratio-norm-NArm.gct -r 27 -c 23 -o luad_test.tsv -s names_test.txt -p outliers_plus_test.tsv -m outliers_minus_test.tsv -b outliers_both_test.tsv

### Tidy Data

# Isolate only the gene ID and sample columns
# Export this as a new dataframe 
df <- read.table(opt$f, skip = 2, header = TRUE, sep = "\t", fill = TRUE)
df_gene <- df[c(opt$r:nrow(df)), c(opt$c:ncol(df))]
as.tibble(df_gene)
write.table(df_gene, file=opt$o, quote=FALSE, sep='\t', col.names = NA)

# Get a list of sample names for the outlier analysis
df_names <- colnames(df_gene[3:ncol(df_gene)])
write.table(df_names, opt$s, sep=',',row.names=F, col.names = F)


### Visualize Distribution 





### Perform Outlier Analysis 

# Create dataframe without gene ids and one with only gene ids
df_stats <- df_gene[,2:ncol(df_gene)]
df_genesymbol <- df_gene[,1:1]
df_genesymbol <- as.data.frame(df_genesymbol)
as.tibble(df_genesymbol)
as.tibble(df_stats)

# Change all values to numeric
df_stats[] <- lapply(df_stats, function(x) as.numeric(as.character(x)))

# Calculate the row median and iqr, make sure they are numeric
df_stats$row_medians <- apply(df_stats, 1 , median, na.rm = TRUE)
df_stats$row_medians <- as.numeric(as.character((df_stats$row_medians)))
df_stats$row_iqr <- apply(df_stats, 1 , IQR, na.rm = TRUE)
df_stats$row_iqr <- as.numeric(as.character((df_stats$row_iqr)))

# Calculate the outliers above and below 1.5*iqr
df_stats$row_med_plus = df_stats$row_medians + (1.5*df_stats$row_iqr)
df_stats$row_med_minus = df_stats$row_medians - (1.5*df_stats$row_iqr)

# Change all values calculated values to numeric
df_stats$row_med_plus <- as.numeric(as.character((df_stats$row_med_plus)))
df_stats$row_med_minus <- as.numeric(as.character((df_stats$row_med_minus)))

# Replace the outliers above 1.5*iqr with 1, all else replace with zero
df_outlier_plus <- apply(df_stats, 2, function(x) ifelse(x > df_stats$row_med_plus, 1, 0))

# Replace the outliers below 1.5*iqr with 1, all else replace with zero
df_outlier_minus <- apply(df_stats, 2, function(x) ifelse(x < df_stats$row_med_minus, 1, 0))

# Replace the outliers above 1.5*iqr with 1 and replace the outliers below 1.5*iqr with -1, all else replace with zero
df_outlier_both <- apply(df_stats, 2, function(x) ifelse(x > df_stats$row_med_plus, 1, ifelse(x < df_stats$row_med_minus, -1, 0)))

# Remove the calculation columns and add a count column to the plus/minus dataframes
df_outlier_plus <- df_outlier_plus[,-c((ncol(df_outlier_plus)-3):ncol(df_outlier_plus))]
df_outlier_plus <- cbind(df_outlier_plus, count = rowSums(df_outlier_plus))
df_outlier_minus <- df_outlier_minus[,-c((ncol(df_outlier_plus)-3):ncol(df_outlier_minus))]
df_outlier_minus <- cbind(df_outlier_minus, count = rowSums(df_outlier_minus))

# Bind dataframe with original gene IDs
df_outlier_gene_plus <- cbind(df_genesymbol, df_outlier_plus)
df_outlier_gene_minus <- cbind(df_genesymbol, df_outlier_minus)
df_outlier_gene_both <- cbind(df_genesymbol, df_outlier_both)

# Aggregate genes by gene symbol to remove duplicates
# When aggregating both plus and minus take the mean of the values rather than the sum
df_outlier_plus_agg <- aggregate(. ~ df_genesymbol, data = df_outlier_gene_plus, sum)
df_outlier_minus_agg <- aggregate(. ~ df_genesymbol, data = df_outlier_gene_minus, sum)

# should I aggregate both?
df_outlier_both_agg <- aggregate(. ~ df_genesymbol, data = df_outlier_gene_both, sum)

# Write out outlier dataframes
write.table(df_outlier_plus_agg, opt$p, quote=FALSE, sep='\t', col.names = NA)
write.table(df_outlier_minus_agg, opt$m, quote=FALSE, sep='\t', col.names = NA)
write.table(df_outlier_both_agg, opt$b, quote=FALSE, sep='\t', col.names = NA)



