#!/Applications/anaconda3/pkgs/r-base-3.4.3-h9b62496_0/bin/Rscript

library(tidyverse)
library(dplyr)
library(ggplot2)
library(optparse, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(gridExtra, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

options(warn=1)

# Set up command line
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--row_gene"), type="integer", default=1, 
              help="index number(s) of row where gene values begins", metavar="character"),
  make_option(c("-c", "--column_gene"), type="integer", default=1, 
              help="index number(s) of column where gene values begins", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="output file name for all outputs, do not include file extension", metavar="character"),
  make_option(c("-a", "--aes"), type="character", default="#B2182B", 
              help="aesthetic of distribution plot: color", metavar="character"),
  make_option(c("-u", "--upper_x_lime"), type="integer", default="20", 
              help="upper x limit for distribution plot", metavar="character"),
  make_option(c("-l", "--lower_x_lim"), type="integer", default="-20", 
              help="lower x limit for distribution plot", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# test script:
# Rscript luad_analysis_pipeline_test.R -f luad-v1.0-phosphoproteome-ratio-norm-NArm.gct -r 27 -c 23 -o xxxxxxxx


### Establish names 
out_dataframe <- paste(opt$o, ".tsv", sep="")
out_names <- paste(opt$o,"_names.txt", sep="")
out_distribution <- paste(opt$o, "_distribution.tiff", sep = "")
out_outlier_plus <- paste(opt$o,"_outlier_plus.tsv", sep="")
out_outlier_minus <- paste(opt$o,"_outlier_minus.tsv", sep="")
out_outlier_both <- paste(opt$o,"_outlier_both.tsv", sep="")
ggplot_title <- paste(opt$o, "Distribution", sep=" ")

### Tidy Data

# Isolate only the gene ID and sample columns
# Export this as a new dataframe 
df <- read.table(opt$f, skip = 2, header = TRUE, sep = "\t", fill = TRUE)
df_gene <- df[c(opt$r:nrow(df)), c(opt$c:ncol(df))]
as.tibble(df_gene)
write.table(df_gene, file=out_dataframe, quote=FALSE, sep='\t', col.names = NA)

# Get a list of sample names for the outlier analysis
df_names <- colnames(df_gene[3:ncol(df_gene)])
write.table(df_names, out_names, sep=',',row.names=F, col.names = F)


### Visualize Distribution 

# Gather only the values in the dataframe, use the sample name as a key
df_long <- gather(df_gene[3:ncol(df_gene)], na.rm = TRUE)
df_long$value <- as.numeric(as.character(df_long$value))

# Create plot
# Options given to optimize 
df_hist <- ggplot(data = df_long) +
  geom_histogram(mapping = aes(x = value), fill = opt$a, binwidth = 0.01) + 
  xlim(opt$l,opt$u) + ylab("Count") + xlab('Normalized Value') + 
  ggtitle(ggplot_title) + theme(plot.title = element_text(hjust = 0.5))

# Write out distribution plot
tiff(out_distribution, units="in", width=5, height=5, res=300)
df_hist
dev.off()


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
write.table(df_outlier_plus_agg, out_outlier_plus, quote=FALSE, sep='\t', col.names = NA)
write.table(df_outlier_minus_agg, out_outlier_minus, quote=FALSE, sep='\t', col.names = NA)
write.table(df_outlier_both_agg, out_outlier_both, quote=FALSE, sep='\t', col.names = NA)



