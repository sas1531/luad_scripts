#!/Applications/anaconda3/pkgs/r-base-3.4.3-h9b62496_0/bin/Rscript

library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(optparse, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(gridExtra, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(data.table)
library(ComplexHeatmap, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(circlize, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")


options(warn=1)

# Set up command line
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input file name", metavar="character"),
  make_option(c("--skip"), type="character", default=2, 
              help="How many lines to skip, not including the header", metavar="integer"),
  make_option(c("--row_gene"), type="integer", default=1, 
              help="Row number where the genes and values begin", metavar="character"),
  make_option(c("--column_gene"), type="integer", default=1, 
              help="Column number where the genes and values begin", metavar="character"),
  make_option(c("--meta_column"), type="integer", default=1, 
              help="Column number where gene IDs are located for labelling", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="Output file prefix for all outlier outputs, do not include file extension", 
              metavar="character"),
  make_option(c("--aes"), type="character", default="#B2182B", 
              help="aesthetic of distribution plot: color", metavar="character"),
  make_option(c("-u", "--upper_x_lim"), type="integer", default="20", 
              help="Upper x limit for distribution plot", metavar="character"),
  make_option(c("-l", "--lower_x_lim"), type="integer", default="-20", 
              help="Lower x limit for distribution plot", metavar="character"),
  make_option(c("--tag_phospho"), type="character", default="", 
              help="Tag for phospho data: write 'phospho' if the data is phosphoproteomic
              and 'not_phospho' if it isn't", metavar="character"),
  make_option(c("-b", "--base_log"), type="character", default="not_log10", 
              help="if data is in 'log10' specify here and it will be transformed 
              into log2: options c(not_log10, log10)", metavar="character"),
  make_option(c("--prot"), type="character", default="no", 
              help="Yes if the input data is proteomic and includes isoforms with the 
              same gene name", metavar="character")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

### Establish names 
out_dataframe <- paste(opt$o, ".tsv", sep="")
out_names <- paste(opt$o,"_names.txt", sep="")
out_meta_data <- paste(opt$o,"_meta.tsv", sep="")
out_distribution <- paste(opt$o, "_distribution.png", sep = "")
out_outlier_plus <- paste(opt$o,"_outlier_up.tsv", sep="")
out_outlier_minus <- paste(opt$o,"_outlier_down.tsv", sep="")
ggplot_title <- paste(opt$o, "Distribution", sep=" ")

### Functions

### Log transformation function
# If log is log10, change to log2
log10_transform <- function(df_gene, log){
  if (log == 'log10'){
    GeneSymbol <- df_gene[, 2]
    df_gene_log10 <- lapply(df_gene[,c(3:ncol(df_gene))], 
                            function(x) as.numeric(as.character(x)))
    df_gene_log <- lapply(df_gene_log10, function(x) (10**(x)))
    df_gene_log2 <- lapply(df_gene_log,function(x) (log2(x)))
    df_gene_log2_final <- cbind(as.data.frame(GeneSymbol), as.data.frame(df_gene_log2))
    return(df_gene_log2_final)
  } else {
    return(df_gene)
  }
}

### Phosphosite aggregation function
# If data is phosphoproteomic, aggregate and sum outlier values 
agg_phospho <- function(df_outlier_agg, df_not_outlier_agg, tag){
  if (tag == 'phospho'){
    df_final_outlier_agg <- aggregate(. ~ GeneSymbol.out, data = df_outlier_agg, 
                                      count_outliers, na.action = na.pass)
    df_final_not_outlier_agg <- aggregate(. ~ GeneSymbol, data = df_not_outlier_agg, 
                                          count_not_outliers, na.action = na.pass)
    df_final_outlier <- cbind(df_final_not_outlier_agg, df_final_outlier_agg)
    df_final_outlier <- dplyr::select(df_final_outlier, -GeneSymbol.out)
    df_final_outlier$GeneSymbol <- as.factor(as.character(df_final_outlier$GeneSymbol))
    return(df_final_outlier)
  } else 
    df_not_outlier_agg <- dplyr::select(df_not_outlier_agg, -GeneSymbol)
    df_outlier_agg <- dplyr::select(df_outlier_agg, -GeneSymbol.out)
    df_not_outlier_agg <- lapply(df_not_outlier_agg, function(x) ifelse(x == 0, 1, 0))
    df_final_outlier <- cbind(df_genesymbol, df_not_outlier_agg)
    df_final_outlier <- cbind(df_final_outlier, df_outlier_agg)
    colnames(df_final_outlier)[1] <- "GeneSymbol"
    return(df_final_outlier)
}

### Tidy Data

# Isolate only the gene ID and sample columns
# Export this as a new dataframe 
df <- read.table(opt$f, skip = opt$skip, header = TRUE, sep = "\t", fill = TRUE)
df <- as.data.frame(df)
colnames(df)[opt$meta_column] <- "GeneSymbol"
df_gene <- df[c(opt$row_gene:nrow(df)), c(1, opt$meta_column, opt$column_gene:ncol(df))]

# If the data is proteomic and has isoforms, combine the genesymbol with id (isoform)
if (opt$prot == "yes"){
  df_gene$GeneSymbol <- paste(df_gene$GeneSymbol, df_gene$id, sep=" - ")
}


# Change dataframe values from log10 to log2 if needed
df_gene <- log10_transform(df_gene, opt$b)

# Get a list of sample names for the outlier analysis
df_names <- colnames(df_gene[3:ncol(df_gene)])

# Isolate meta data
meta_df <- df[c(1:((opt$row_gene)-1)), c(1, opt$meta_column, opt$column_gene:ncol(df))]

### Visualize Distribution 

# Gather only the values in the dataframe, use the sample name as a key
df_long <- gather(df_gene[3:ncol(df_gene)], na.rm = TRUE)
df_long$value <- as.numeric(as.character(df_long$value))

# Create plot
# Options given to optimize 
df_hist <- ggplot(data = df_long) +
  geom_histogram(mapping = aes(x = value), fill = opt$aes, binwidth = 0.01) + 
  xlim(opt$l,opt$u) + ylab("Count") + xlab('Normalized Value') + 
  ggtitle(ggplot_title) + theme(plot.title = element_text(hjust = 0.5))


### Perform Outlier Analysis 

# Create dataframe without gene ids and one with only gene ids
df_stats <- df_gene[,3:ncol(df_gene)]
df_genesymbol <- dplyr::select(df_gene, GeneSymbol)

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
#df_outlier_both <- apply(df_stats, 2, function(x) ifelse(x > df_stats$row_med_plus, 1, ifelse(x < df_stats$row_med_minus, -1, 0)))

# Remove the calculation columns and add a count column to the plus/minus dataframes
df_outlier_plus <- df_outlier_plus[,-c(((ncol(df_outlier_plus))-3):ncol(df_outlier_plus))]
df_outlier_plus <- cbind(df_outlier_plus, count = rowSums(df_outlier_plus))
df_outlier_minus <- df_outlier_minus[,-c(((ncol(df_outlier_minus))-3):ncol(df_outlier_minus))]
df_outlier_minus <- cbind(df_outlier_minus, count = rowSums(df_outlier_minus))

# Bind dataframe with original gene IDs
df_outlier_gene_plus <- cbind(df_genesymbol, df_outlier_plus)
df_outlier_gene_minus <- cbind(df_genesymbol, df_outlier_minus)

# Aggregate phosphosites 
count_not_outliers <-  function(x) {length(which(x==0))}
count_outliers <- function(x) {length(which(x==1))}

# Plus
df_final_outlier_plus <- dplyr::select(df_outlier_gene_plus, -count)
colnames(df_final_outlier_plus) <- paste(colnames(df_final_outlier_plus), "outlier", sep = "_")
colnames(df_final_outlier_plus)[1] <- "GeneSymbol.out"
df_final_not_outlier_plus <- dplyr::select(df_outlier_gene_plus, -count)
colnames(df_final_not_outlier_plus) <- paste(colnames(df_final_not_outlier_plus), "not_outlier", sep = "_")
colnames(df_final_not_outlier_plus)[1] <- "GeneSymbol"

# Minus
df_final_outlier_minus <- dplyr::select(df_outlier_gene_minus, -count)
colnames(df_final_outlier_minus) <- paste(colnames(df_final_outlier_minus), "outlier", sep = "_")
colnames(df_final_outlier_minus)[1] <- "GeneSymbol.out"
df_final_not_outlier_minus <- dplyr::select(df_outlier_gene_minus, -count)
colnames(df_final_not_outlier_minus) <- paste(colnames(df_final_not_outlier_minus), "not_outlier", sep = "_")
colnames(df_final_not_outlier_minus)[1] <- "GeneSymbol"

# Run aggregation function
df_final_plus <- agg_phospho(df_final_outlier_plus, df_final_not_outlier_plus, opt$tag_phospho)
df_final_minus <- agg_phospho(df_final_outlier_minus, df_final_not_outlier_minus, opt$tag_phospho)

# Export tables and figures
# Write out gene and sample dataframe
write.table(df_gene, file=out_dataframe, quote=FALSE, sep='\t', col.names = NA)
# Write out meta data
write.table(meta_df, file=out_meta_data, quote=FALSE, sep='\t', col.names = NA)
# Write out sample names
write.table(df_names, out_names, sep=',', row.names=F, col.names = F)
# Write out outlier dataframes
write.table(df_final_plus, file = out_outlier_plus, quote=FALSE, sep='\t', col.names = NA)
write.table(df_final_minus, file = out_outlier_minus, quote=FALSE, sep='\t', col.names = NA)
# Write out distribution plot
png(out_distribution, units="in", width=5, height=5, res=300)
df_hist
dev.off()

print("Outlier Analyis Complete")
