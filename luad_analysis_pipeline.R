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
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--row_gene"), type="integer", default=1, 
              help="index number of row where gene values begins", metavar="character"),
  make_option(c("--column_gene"), type="integer", default=1, 
              help="index number of column where sample values begins", metavar="character"),
  make_option(c("-m", "--meta_column"), type="integer", default=1, 
              help="index number of column where gene IDs are located", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="output file name for all outlier outputs, do not include file extension", metavar="character"),
  make_option(c("-a", "--aes"), type="character", default="#B2182B", 
              help="aesthetic of distribution plot: color", metavar="character"),
  make_option(c("-u", "--upper_x_lime"), type="integer", default="20", 
              help="upper x limit for distribution plot", metavar="character"),
  make_option(c("-l", "--lower_x_lim"), type="integer", default="-20", 
              help="lower x limit for distribution plot", metavar="character"),
  make_option(c("--tag_phospho"), type="character", default="", 
              help="tag for phospho data: write 'phospho' if the data is phosproteome 
              and 'not_phospho' if it isn't", metavar="character"),
  make_option(c("-b", "--base_log"), type="character", default="not_log10", 
              help="if data is in 'log10' specify here and it will be transformed 
              into log2: options c(not_log10, log10)", metavar="character"),
  make_option(c("--comparison"), type="character", default="Type", 
              help="column name that identifies comparison (e.g. Region, Gender, Type)", metavar="character"),
  make_option(c("--group_1"), type="character", default="Tumor", 
              help="Group one of comparison", metavar="character"),
  make_option(c("--group_2"), type="character", default="Normal", 
              help="Group 2 of comparison", metavar="character"),
  make_option(c("--normal_tumor"), type="character", default="both", 
              help="If your samples include both tumor and normal samples tag 'both'", metavar="character"),
  make_option(c("--tumor_column"), type="character", default="NULL", 
              help="Name of tumor column (only required if there are normal and tumor samples", metavar="character"),
  make_option(c("--tumor_group"), type="character", default="NULL", 
              help="Name of tumor group (only reguired if there are normal and tumor samples", metavar="character"),
  make_option(c("--prot"), type="character", default="yes", 
              help="yes if the input data is proteomic and includes isoforms with the same gene name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


# test script:
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-phosphoproteome-ratio-norm-NArm.gct -r 27 --column_gene 24 -m 23 -o luad_phos_tumor --tag_phospho phospho --comparison Region.of.Origin --group_1 asian --group_2 western --normal_tumor both --tumor_column Type --tumor_group Tumor 
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-proteome-ratio-norm-NArm.gct -r 27 --column_gene 18 -m 17 -o luad_prot_tumor --tag_phospho not_phospho --comparison Region.of.Origin --group_1 asian --group_2 western --normal_tumor both --tumor_column Type --tumor_group Tumor --prot yes
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-rnaseq-linear-gene-fpkm-uq-log10-NArm-row-norm.gct -r 22 --column_gene 5 -m 4 -o luad_rna_tumor_log2 --tag_phospho not_phospho -b log10 --comparison Region.of.Origin --group_1 asian --group_2 western --normal_tumor both --tumor_column Type --tumor_group Tumor
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-rnaseq-linear-gene-fpkm-uq-log10-NArm-row-norm.gct -r 22 --column_gene 5 -m 4 -o luad_rna_tumor_log10 --tag_phospho not_phospho --comparison Region.of.Origin --group_1 asian --group_2 western --normal_tumor both --tumor_column Type --tumor_group Tumor
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-wgs-cnv-somatic.luad.all.hg38-v2.0-20190121.gct -r 22 --column_gene 6 -m 2 -o luad_cnv_tumor --tag_phospho not_phospho --comparison Region.of.Origin --group_1 asian --group_2 western --normal_tumor both --tumor_column Type --tumor_group Tumor

### Establish names 
out_dataframe <- paste(opt$o, ".tsv", sep="")
out_names <- paste(opt$o,"_names.txt", sep="")
out_meta_data <- paste(opt$o,"_meta.tsv", sep="")
out_distribution <- paste(opt$o, "_distribution.tiff", sep = "")
out_outlier_plus <- paste(opt$o,"_outlier_up.tsv", sep="")
out_outlier_minus <- paste(opt$o,"_outlier_down.tsv", sep="")
out_outlier_both <- paste(opt$o,"_outlier_both.tsv", sep="")
ggplot_title <- paste(opt$o, "Distribution", sep=" ")
plus_sig_output <- paste(opt$o, "_up_", opt$group_1, "_", opt$group_2, "_sig_outliers.txt", sep="")
minus_sig_output <- paste(opt$o, "_down_", opt$group_1, "_", opt$group_2, "_sig_outliers.txt", sep="")
plus_heatmap <- paste(opt$o,"_up_heatmap_sig_", opt$group_1, "_", opt$group_2, ".tiff", sep="")
minus_heatmap <- paste(opt$o,"_down_heatmap_sig_", opt$group_1, "_", opt$group_2, ".tiff", sep="")


### Tidy Data

# Isolate only the gene ID and sample columns
# Export this as a new dataframe 
df <- read.table(opt$f, skip = 2, header = TRUE, sep = "\t", fill = TRUE)
df <- as.data.frame(df)
colnames(df)[opt$m] <- "GeneSymbol"
df_gene <- df[c(opt$r:nrow(df)), c(1, opt$m, opt$column_gene:ncol(df))]

# If the data is proteomic and has isoforms, combine the genesymbol with id (isoform)
if (opt$prot == "yes"){
  df_gene$GeneSymbol <- paste(df_gene$GeneSymbol, df_gene$id, sep=" - ")
}


# Change dataframe values from log10 to log2 if needed
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
df_gene <- log10_transform(df_gene, opt$b)

# Write out gene and sample dataframe
write.table(df_gene, file=out_dataframe, quote=FALSE, sep='\t', col.names = NA)

# Get a list of sample names for the outlier analysis
df_names <- colnames(df_gene[3:ncol(df_gene)])
write.table(df_names, out_names, sep=',',row.names=F, col.names = F)

# Isolate meta data
meta_df <- df[c(1:((opt$r)-1)), c(1, opt$m, opt$column_gene:ncol(df))]
write.table(meta_df, file=out_meta_data, quote=FALSE, sep='\t', col.names = NA)

### Modify meta data for downstream filtering

# Change meta data factors to characters for word manipulation
meta_df[] <- lapply(meta_df, function(x) as.character(x))
# Transpose 
meta_df_t <- dplyr::select(meta_df,-GeneSymbol)
rownames(meta_df_t) <- meta_df$id
meta_df_t$id <- NULL
meta_df_t <- as.data.frame(t(meta_df_t))
meta_df_t <- setDT(meta_df_t, keep.rownames = TRUE)[]

# Make outliers and not outliers columns
# Add extra column stating whether or not the column is an outlier or not
meta_df_outlier <- dplyr::select(meta_df,-GeneSymbol)
colnames(meta_df_outlier) <- paste(colnames(meta_df_outlier), "outlier", sep = "_")
colnames(meta_df_outlier)[1] <- c("id.out")
outlier_list <- list(c("yes"))
outlier_list <- rep(outlier_list, ncol(meta_df_outlier))
meta_df_outlier[nrow(meta_df_outlier)+1,] <- outlier_list
meta_df_outlier$id.out <- as.character(meta_df_outlier$id.out)
meta_df_outlier$id.out[nrow(meta_df_outlier)] <- "outlier"
meta_df_not_outlier <- dplyr::select(meta_df, -GeneSymbol)
colnames(meta_df_not_outlier) <- paste(colnames(meta_df_not_outlier), "not_outlier", sep = "_")
colnames(meta_df_not_outlier)[1] <- c("id")
not_outlier_list <- list(c("no"))
not_outlier_list <- rep(not_outlier_list, ncol(meta_df_not_outlier))
meta_df_not_outlier[nrow(meta_df_not_outlier)+1,] <- not_outlier_list
meta_df_not_outlier$id <- as.character(meta_df_not_outlier$id)
meta_df_not_outlier$id[nrow(meta_df_not_outlier)] <- "outlier"
meta_df_final <- cbind(meta_df_not_outlier, meta_df_outlier)
meta_df_final <- dplyr::select(meta_df_final, -id.out)

# Transpose for selection, keep rownames
meta_df_final_t <- meta_df_final
rownames(meta_df_final_t) <- meta_df_final_t$id
meta_df_final_t$id <- NULL
meta_df_final_t <- as.data.frame(t(meta_df_final_t))
meta_df_final_t <- setDT(meta_df_final_t, keep.rownames = TRUE)[]

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
#df_outlier_both <- df_outlier_both[,-c(((ncol(df_outlier_both))-3):ncol(df_outlier_both))]
#df_outlier_both <- cbind(df_outlier_both, count = rowSums(df_outlier_both))

# Bind dataframe with original gene IDs
df_outlier_gene_plus <- cbind(df_genesymbol, df_outlier_plus)
df_outlier_gene_minus <- cbind(df_genesymbol, df_outlier_minus)
#df_outlier_gene_both <- cbind(df_genesymbol, df_outlier_both)

# Aggregate dfdfites 
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

agg_phospho <- function(df_outlier_agg, df_not_outlier_agg, tag){
  if (tag == 'phospho'){
    # Plus
    df_final_outlier_agg <- aggregate(. ~ GeneSymbol.out, data = df_outlier_agg, 
                                       count_outliers, na.action = na.pass)
    df_final_not_outlier_agg <- aggregate(. ~ GeneSymbol, data = df_not_outlier_agg, 
                                           count_not_outliers, na.action = na.pass)
    df_final_outlier <- cbind(df_final_not_outlier_agg, df_final_outlier_agg)
    df_final_outlier <- dplyr::select(df_final_outlier, -GeneSymbol.out)
    df_final_outlier$GeneSymbol <- as.factor(as.character(df_final_outlier$GeneSymbol))
    return(df_final_outlier)
  } else 
    # Replace df_not_outlier_agg
    df_not_outlier_agg <- dplyr::select(df_not_outlier_agg, -GeneSymbol)
    df_outlier_agg <- dplyr::select(df_outlier_agg, -GeneSymbol.out)
    df_not_outlier_agg <- lapply(df_not_outlier_agg, function(x) ifelse(x == 0, 1, 0))
    df_final_outlier <- cbind(df_genesymbol, df_not_outlier_agg)
    df_final_outlier <- cbind(df_final_outlier, df_outlier_agg)
    colnames(df_final_outlier)[1] <- "GeneSymbol"
    return(df_final_outlier)
}

# Run aggregation function
df_final_plus <- agg_phospho(df_final_outlier_plus, df_final_not_outlier_plus, opt$tag_phospho)
df_final_minus <- agg_phospho(df_final_outlier_minus, df_final_not_outlier_minus, opt$tag_phospho)

# Write out outlier dataframes
write.table(df_final_plus, file = out_outlier_plus, quote=FALSE, sep='\t', col.names = NA)
write.table(df_final_minus, file = out_outlier_minus, quote=FALSE, sep='\t', col.names = NA)
#write.table(df_outlier_both_agg, out_outlier_both, quote=FALSE, sep='\t', col.names = NA)


### Group Comparison

# Isolate Group Samples

comparison <- grep(opt$comparison, colnames(meta_df_t))
colnames(meta_df_t)[comparison] <- c('comp')
sample <- grep(opt$tumor_column, colnames(meta_df_t))
colnames(meta_df_t)[sample] <- c('samp')
comparison <- grep(opt$comparison, colnames(meta_df_final_t))
colnames(meta_df_final_t)[comparison] <- c('comp')
sample <- grep(opt$tumor_column, colnames(meta_df_final_t))
colnames(meta_df_final_t)[sample] <- c('samp')


if (opt$normal_tumor == 'both'){
  all_col <- meta_df_t[((meta_df_t$comp == opt$group_1 | meta_df_t$comp == opt$group_2) &
                         (meta_df_t$samp == opt$tumor_group)), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_1 & meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_1 & meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_2 & meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_2 & meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
  } else {
    all_col <- meta_df_t[(meta_df_t$comp == opt$group_1 | meta_df_t$comp == opt$group_2), ]
    all_col <- as.vector(all_col[,1])
    outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_1 & meta_df_final_t$outlier == 'yes'), ]
    outlier_group1_col <- (outlier_group1_col[,1])
    not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_1 & meta_df_final_t$outlier == 'no'), ]
    not_outlier_group1_col <- (not_outlier_group1_col[,1])
    outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_2 & meta_df_final_t$outlier == 'yes'), ]
    outlier_group2_col <- (outlier_group2_col[,1])
    not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp == opt$group_2 & meta_df_final_t$outlier == 'no'), ]
    not_outlier_group2_col <- (not_outlier_group2_col[,1])
}

all_col <- as.vector(all_col[[1]])
outlier_group1_col <- as.vector(outlier_group1_col[[1]])
not_outlier_group1_col <- as.vector(not_outlier_group1_col[[1]])
outlier_group2_col <- as.vector(outlier_group2_col[[1]])
not_outlier_group2_col <- as.vector(not_outlier_group2_col[[1]])


# Plus comparison pipeline 

# sum all four groups separately
# put these values in 4 new columns
# take the values in each of these columns to create a contingecy table
# run fischer on this table 
# record the p-value from this table in a new column

# Only include targeted samples
df_final_plus_sample <- df_final_plus[, c("GeneSymbol", (outlier_group1_col), 
                                          (outlier_group2_col), 
                                          (not_outlier_group1_col),
                                          (not_outlier_group2_col))]

# Outlier group1 
outlier_group1_col_number <- which(colnames(df_final_plus) %in% (outlier_group1_col))
outlier_group_1_df <- df_final_plus[outlier_group1_col_number]
outlier_group_1_df_final <- data.frame(outlier_group1 = apply(outlier_group_1_df, 1, sum, na.rm = TRUE))
outlier_group_1_df_final$outlier_group1 <- as.numeric(outlier_group_1_df_final$outlier_group1)
df_final_plus <- cbind(df_final_plus_sample, outlier_group_1_df_final)
# Not Outlier group1 
not_outlier_group1_col_number <- which(colnames(df_final_plus) %in% not_outlier_group1_col)
not_outlier_group_1_df <- df_final_plus[not_outlier_group1_col_number]
not_outlier_group_1_df_final <- data.frame(not_outlier_group1 = apply(not_outlier_group_1_df, 1, sum, na.rm = TRUE))
not_outlier_group_1_df_final$not_outlier_group1 <- as.numeric(not_outlier_group_1_df_final$not_outlier_group1)
df_final_plus <- cbind(df_final_plus, not_outlier_group_1_df_final)
# Outlier group2
outlier_group2_col_number <- which(colnames(df_final_plus) %in% outlier_group2_col)
outlier_group_2_df <- df_final_plus[outlier_group2_col_number]
outlier_group_2_df_final <- data.frame(outlier_group2 = apply(outlier_group_2_df, 1, sum, na.rm = TRUE))
outlier_group_2_df_final$outlier_group2 <- as.numeric(outlier_group_2_df_final$outlier_group2)
df_final_plus <- cbind(df_final_plus, outlier_group_2_df_final)
# Not Outlier group2
not_outlier_group2_col_number <- which(colnames(df_final_plus) %in% not_outlier_group2_col)
not_outlier_group_2_df <- df_final_plus[not_outlier_group2_col_number]
not_outlier_group_2_df_final <- data.frame(not_outlier_group2 = apply(not_outlier_group_2_df, 1, sum, na.rm = TRUE))
not_outlier_group_2_df_final$not_outlier_group2 <- as.numeric(not_outlier_group_2_df_final$not_outlier_group2)
df_final_plus <- cbind(df_final_plus, not_outlier_group_2_df_final)
df_final_plus <- as.data.frame(df_final_plus)

# Add outlier count column to dataframe - this is the count of number of outlier samples in group 1
outlier_count <- df_final_plus[, c(outlier_group1_col)]
df_final_plus$count_outliers <- rowSums((outlier_count[, 2:(ncol(outlier_count))])!=0, na.rm = TRUE)
df_final_plus$count_outliers <- as.numeric(df_final_plus$count_outliers)

# Create outlier fractions (outlier divided by total)
#### Filter by fraction
# Filter so that you only have the genes that are upregulated in group
df_names <- as.vector(c(outlier_group1_col, outlier_group2_col))
df_names <- gsub('.{8}$', '', df_names)
combine_samples <- sapply(df_names, function(xx) rowSums(df_final_plus[,grep(xx, names(df_final_plus)), drop=FALSE], na.rm = TRUE))
combine_samples <- as.data.frame(combine_samples)
outlier_samples_only <- df_final_plus[, c(1, 2:((ncol(combine_samples))+1))]
fraction_outlier <- cbind(outlier_samples_only[1], round((outlier_samples_only[-1]/combine_samples), 6))
df_plus_fraction_outlier <- cbind(fraction_outlier, df_final_plus[,(ncol(df_final_plus)-4):ncol(df_final_plus)])

### Create group 1 fraction mean and group 2 fraction mean columns 
# Group 1
fraction_group1_col_number <- which(colnames(df_plus_fraction_outlier) %in% outlier_group1_col)
fraction_mean_group1 <- df_plus_fraction_outlier[fraction_group1_col_number]
fraction_mean_group1 <- data.frame(group1_fraction_mean = apply(fraction_mean_group1, 1, mean, na.rm=TRUE))
fraction_mean_group1[is.na(fraction_mean_group1)] <- 0
df_plus_fraction_outlier <- cbind(df_plus_fraction_outlier, fraction_mean_group1)

# Group 2
fraction_group2_col_number <- which(colnames(df_plus_fraction_outlier) %in% outlier_group2_col)
fraction_mean_group2 <- df_plus_fraction_outlier[fraction_group2_col_number]
fraction_mean_group2 <- data.frame(group2_fraction_mean = apply(fraction_mean_group2, 1, mean, na.rm=TRUE))
fraction_mean_group2[is.na(fraction_mean_group2)] <- 0
df_plus_fraction_outlier <- cbind(df_plus_fraction_outlier, fraction_mean_group2)

# create contingency table 
# If number of outlier in group 1 is greater than 0.3 * number of outliers

fish_exact <- function(y){
  min_num_outliers <- as.numeric(length(outlier_group1_col)*0.3)
  num_outlier_samples <- as.numeric(y['count_outliers'])
  if (num_outlier_samples > min_num_outliers){
    table <- (matrix(as.numeric(c(y['outlier_group1'], y['outlier_group2'],
                                  y['not_outlier_group1'], y['not_outlier_group2'])),
                     nrow = 2))
    p <- fisher.test(table)$p.value
    p
  } else {
    p <- NA
    p
  }
}


# Create a p-value column
fisher_df <- data.frame(p_value = apply(df_plus_fraction_outlier, 1, fish_exact))
df_plus_fraction_outlier <- cbind(df_plus_fraction_outlier, fisher_df)

# Create a column with the corrected p-values 
p_adjust_df <- data.frame(lapply(df_plus_fraction_outlier['p_value'], p.adjust, method = "BH"))
colnames(p_adjust_df)[1] <- c("q_value")
df_plus_fraction_outlier <- cbind(df_plus_fraction_outlier, p_adjust_df)

# Get significant outliers
# If p_value_corrected is greater than FDR cutoff, insert yes (significant outlier)
df_plus_fraction_outlier$sig_outlier_q <- NA
df_plus_fraction_outlier$sig_outlier_q[df_plus_fraction_outlier$q_value <= 0.05] <- "yes"
df_plus_fraction_outlier$sig_outlier_q[df_plus_fraction_outlier$q_value > 0.05] <- "no"

# If fraction group 1 mean is greater than group 2 fraction mean, insert yes (significant outlier)
df_plus_fraction_outlier$sig_outlier_mean[(df_plus_fraction_outlier$group1_fraction_mean) > (df_plus_fraction_outlier$group2_fraction_mean)] <- "yes"
df_plus_fraction_outlier$sig_outlier_mean[df_plus_fraction_outlier$group1_fraction_mean <= df_plus_fraction_outlier$group2_fraction_mean ] <- "no"

## Length of list of signfiicant outliers will be the number of outliers = print this
df_plus_fraction_outlier <- filter(df_plus_fraction_outlier,
                                   sig_outlier_q == 'yes' & sig_outlier_mean == 'yes')
df_plus_fraction_outlier <- df_plus_fraction_outlier[order(df_plus_fraction_outlier$q_value, decreasing = FALSE), ]

# Export list of genes and FDR 
sig_outlier_genes_q_value <- df_plus_fraction_outlier[c(1, ncol(df_plus_fraction_outlier)-2)]
write.table(sig_outlier_genes_q_value, file = plus_sig_output, sep='\t',row.names=F, col.names = T, quote=FALSE)

### Make heat map of

# Define heatmap dataframe
# Make gene symbol row.names 
# Isolate values
# Make values into matrix
plus_fraction_heat_1 <- df_plus_fraction_outlier
plus_fraction_heat <- plus_fraction_heat_1[,-1]
rownames(plus_fraction_heat) <- plus_fraction_heat_1[,1]
plus_fraction_heat <- dplyr::select(plus_fraction_heat, -outlier_group1, -not_outlier_group1, -outlier_group2,
                               -not_outlier_group2, -count_outliers, -group1_fraction_mean, -group2_fraction_mean,
                               -p_value, -q_value, -sig_outlier_q, -sig_outlier_mean)
is.na(plus_fraction_heat) <- sapply(plus_fraction_heat, is.infinite)
plus_fraction_heat[is.na(plus_fraction_heat)]<-0

plus_fraction_heat[] <- lapply(plus_fraction_heat, function(x) as.numeric(x))
plus_fraction_heat <- as.matrix(plus_fraction_heat)

## Make annotation bar (pull from meta data)
if (opt$normal_tumor == 'both'){
  plus_heat_filter <- meta_df_t[(meta_df_t$samp == opt$tumor_group), ]
  } else {
    plus_heat_filter <- meta_df_t
}


#comparison <- grep(opt$comparison, colnames(meta_df_t))
#colnames(meta_df_t)[comparison] <- c('comp')

plus_heat_annotation <- as.data.frame(plus_heat_filter$comp)
colnames(plus_heat_annotation)[1] <- c('comp')
plus_heat_annotation <- as.data.frame(plus_heat_annotation[with(plus_heat_annotation, order(comp)),])
colnames(plus_heat_annotation)[1] <- c('comp')

heat_color <- c("red4", "gray85")
names(heat_color) <- levels(plus_heat_annotation$comp)
plus_heat_annotation <- HeatmapAnnotation(df = data.frame(Heat = plus_heat_annotation$comp),
                                          col = list(Heat = heat_color))

if (nrow(plus_fraction_heat) < 25){
  heat_font <- 6
  } else {
  heat_font <- 0
}

frac_heat_1 <- Heatmap(plus_fraction_heat,
                     top_annotation = plus_heat_annotation,
                     show_column_names = F, 
                     show_row_names = T,
                     cluster_rows = F,
                     cluster_columns = F,
                     row_names_gp = gpar(fontsize = heat_font),
                     top_annotation_height = unit(1, "cm"),
                     heatmap_legend_param = list(legend_height = unit(3, "in"),
                                                 title = NULL),
                     col = colorRamp2(c(0, 0.5, 1), c("white", "red2", "red4")),
                     row_title = "Fraction Outliers",
                     row_title_gp = gpar(fontsize = 10))

h1 <- draw(frac_heat_1, heatmap_legend_side = "left")
tiff(plus_heatmap, units="in", width=6, height=5, res=600)
h1
dev.off()


# Minus comparison pipeline 

# sum all four groups separately
# put these values in 4 new columns
# take the values in each of these columns to create a contingecy table
# run fischer on this table 
# record the p-value from this table in a new column

# Only include targeted samples
df_final_minus_sample <- df_final_minus[, c("GeneSymbol", (outlier_group1_col),
                                            (outlier_group2_col),
                                            (not_outlier_group1_col),
                                            (not_outlier_group2_col))]

# Outlier group1 
outlier_group1_col_number <- which(colnames(df_final_minus) %in% (outlier_group1_col))
outlier_group_1_df <- df_final_minus[outlier_group1_col_number]
outlier_group_1_df_final <- data.frame(outlier_group1 = apply(outlier_group_1_df, 1, sum, na.rm = TRUE))
outlier_group_1_df_final$outlier_group1 <- as.numeric(outlier_group_1_df_final$outlier_group1)
df_final_minus <- cbind(df_final_minus_sample, outlier_group_1_df_final)
# Not Outlier group1 
not_outlier_group1_col_number <- which(colnames(df_final_minus) %in% not_outlier_group1_col)
not_outlier_group_1_df <- df_final_minus[not_outlier_group1_col_number]
not_outlier_group_1_df_final <- data.frame(not_outlier_group1 = apply(not_outlier_group_1_df, 1, sum, na.rm = TRUE))
not_outlier_group_1_df_final$not_outlier_group1 <- as.numeric(not_outlier_group_1_df_final$not_outlier_group1)
df_final_minus <- cbind(df_final_minus, not_outlier_group_1_df_final)
# Outlier group2
outlier_group2_col_number <- which(colnames(df_final_minus) %in% outlier_group2_col)
outlier_group_2_df <- df_final_minus[outlier_group2_col_number]
outlier_group_2_df_final <- data.frame(outlier_group2 = apply(outlier_group_2_df, 1, sum, na.rm = TRUE))
outlier_group_2_df_final$outlier_group2 <- as.numeric(outlier_group_2_df_final$outlier_group2)
df_final_minus <- cbind(df_final_minus, outlier_group_2_df_final)
# Not Outlier group2
not_outlier_group2_col_number <- which(colnames(df_final_minus) %in% not_outlier_group2_col)
not_outlier_group_2_df <- df_final_minus[not_outlier_group2_col_number]
not_outlier_group_2_df_final <- data.frame(not_outlier_group2 = apply(not_outlier_group_2_df, 1, sum, na.rm = TRUE))
not_outlier_group_2_df_final$not_outlier_group2 <- as.numeric(not_outlier_group_2_df_final$not_outlier_group2)
df_final_minus <- cbind(df_final_minus, not_outlier_group_2_df_final)
df_final_minus <- as.data.frame(df_final_minus)

# Add outlier count column to dataframe - this is the count of number of outlier samples in group 1
outlier_count <- df_final_minus[, c(outlier_group1_col)]
df_final_minus$count_outliers <- rowSums((outlier_count[, 2:(ncol(outlier_count))])!=0, na.rm = TRUE)
df_final_minus$count_outliers <- as.numeric(df_final_minus$count_outliers)

# Create outlier fractions (outlier divided by total)
#### Filter by fraction
# Filter so that you only have the genes that are upregulated in group
df_names <- as.vector(c(outlier_group1_col, outlier_group2_col))
df_names <- gsub('.{8}$', '', df_names)
combine_samples <- sapply(df_names, function(xx) rowSums(df_final_minus[,grep(xx, names(df_final_minus)), drop=FALSE], na.rm = TRUE))
combine_samples <- as.data.frame(combine_samples)
outlier_samples_only <- df_final_minus[, c(1, 2:((ncol(combine_samples))+1))]
fraction_outlier <- cbind(outlier_samples_only[1], round((outlier_samples_only[-1]/combine_samples), 6))
df_minus_fraction_outlier <- cbind(fraction_outlier, df_final_minus[,(ncol(df_final_minus)-4):ncol(df_final_minus)])

### Create group 1 fraction mean and group 2 fraction mean columns 
# Group 1
fraction_group1_col_number <- which(colnames(df_minus_fraction_outlier) %in% outlier_group1_col)
fraction_mean_group1 <- df_minus_fraction_outlier[fraction_group1_col_number]
fraction_mean_group1 <- data.frame(group1_fraction_mean = apply(fraction_mean_group1, 1, mean, na.rm=TRUE))
df_minus_fraction_outlier <- cbind(df_minus_fraction_outlier, fraction_mean_group1)

# Group 2
fraction_group2_col_number <- which(colnames(df_minus_fraction_outlier) %in% outlier_group2_col)
fraction_mean_group2 <- df_minus_fraction_outlier[fraction_group2_col_number]
fraction_mean_group2 <- data.frame(group2_fraction_mean = apply(fraction_mean_group2, 1, mean, na.rm=TRUE))
df_minus_fraction_outlier <- cbind(df_minus_fraction_outlier, fraction_mean_group2)

# create contingency table 
# If number of outlier in group 1 is greater than 0.3 * number of outliers

fish_exact <- function(y){
  min_num_outliers <- as.numeric(length(outlier_group1_col)*0.3)
  num_outlier_samples <- as.numeric(y['count_outliers'])
  if (num_outlier_samples > min_num_outliers){
    table <- (matrix(as.numeric(c(y['outlier_group1'], y['outlier_group2'],
                                  y['not_outlier_group1'], y['not_outlier_group2'])),
                     nrow = 2))
    p <- fisher.test(table)$p.value
    p
  } else {
    p <- NA
    p
  }
}


# Create a p-value column
fisher_df <- data.frame(p_value = apply(df_minus_fraction_outlier, 1, fish_exact))
df_minus_fraction_outlier <- cbind(df_minus_fraction_outlier, fisher_df)

# Create a column with the corrected p-values 
p_adjust_df <- data.frame(lapply(df_minus_fraction_outlier['p_value'], p.adjust, method = "BH"))
colnames(p_adjust_df)[1] <- c("q_value")
df_minus_fraction_outlier <- cbind(df_minus_fraction_outlier, p_adjust_df)

# Get significant outliers
# If p_value_corrected is greater than FDR cutoff, insert yes (significant outlier)
df_minus_fraction_outlier$sig_outlier_q <- NA
df_minus_fraction_outlier$sig_outlier_q[df_minus_fraction_outlier$q_value <= 0.05] <- "yes"
df_minus_fraction_outlier$sig_outlier_q[df_minus_fraction_outlier$q_value > 0.05] <- "no"

# If fraction group 1 mean is greater than group 2 fraction mean, insert yes (significant outlier)
df_minus_fraction_outlier$sig_outlier_mean[(df_minus_fraction_outlier$group1_fraction_mean) > (df_minus_fraction_outlier$group2_fraction_mean)] <- "yes"
df_minus_fraction_outlier$sig_outlier_mean[df_minus_fraction_outlier$group1_fraction_mean <= df_minus_fraction_outlier$group2_fraction_mean ] <- "no"

## Length of list of signfiicant outliers will be the number of outliers = print this
df_minus_fraction_outlier <- filter(df_minus_fraction_outlier,
                                   sig_outlier_q == 'yes' & sig_outlier_mean == 'yes')
df_minus_fraction_outlier <- df_minus_fraction_outlier[order(df_minus_fraction_outlier$q_value, decreasing = FALSE), ]

# Export list of genes and FDR 
sig_outlier_genes_q_value <- df_minus_fraction_outlier[c(1, ncol(df_minus_fraction_outlier)-2)]
write.table(sig_outlier_genes_q_value, file = minus_sig_output, sep='\t',row.names=F, col.names = T, quote=FALSE)

### Make heat map of

# Define heatmap dataframe
# Make gene symbol row.names 
# Isolate values
# Make values into matrix
minus_fraction_heat_1 <- df_minus_fraction_outlier
minus_fraction_heat <- minus_fraction_heat_1[,-1]
rownames(minus_fraction_heat) <- minus_fraction_heat_1[,1]
minus_fraction_heat <- dplyr::select(minus_fraction_heat, -outlier_group1, -not_outlier_group1, -outlier_group2,
                                    -not_outlier_group2, -count_outliers, -group1_fraction_mean, -group2_fraction_mean,
                                    -p_value, -q_value, -sig_outlier_q, -sig_outlier_mean)
is.na(minus_fraction_heat) <- sapply(minus_fraction_heat, is.infinite)
minus_fraction_heat[is.na(minus_fraction_heat)]<-0
minus_fraction_heat[] <- lapply(minus_fraction_heat, function(x) as.numeric(x))
minus_fraction_heat <- as.matrix(minus_fraction_heat)

## Make annotation bar (pull from meta data)
if (opt$normal_tumor == 'both'){
  minus_heat_filter <- meta_df_t[(meta_df_t$samp == opt$tumor_group), ]
} else {
  minus_heat_filter <- meta_df_t
}


#comparison <- grep(opt$comparison, colnames(meta_df_t))
#colnames(meta_df_t)[comparison] <- c('comp')

minus_heat_annotation <- as.data.frame(minus_heat_filter$comp)
colnames(minus_heat_annotation)[1] <- c('comp')
minus_heat_annotation <- as.data.frame(minus_heat_annotation[with(minus_heat_annotation, order(comp)),])
colnames(minus_heat_annotation)[1] <- c('comp')

heat_color <- c("darkblue", "gray85")
names(heat_color) <- levels(minus_heat_annotation$comp)
minus_heat_annotation <- HeatmapAnnotation(df = data.frame(Heat = minus_heat_annotation$comp),
                                          col = list(Heat = heat_color))

if (nrow(minus_fraction_heat) < 25){
  heat_font <- 6
} else {
  heat_font <- 0
}

frac_heat_2 <- Heatmap(minus_fraction_heat,
                       top_annotation = minus_heat_annotation,
                       show_column_names = F, 
                       show_row_names = T,
                       cluster_rows = F,
                       cluster_columns = F,
                       row_names_gp = gpar(fontsize = heat_font),
                       top_annotation_height = unit(1, "cm"),
                       heatmap_legend_param = list(legend_height = unit(3, "in"),
                                                   title = NULL),
                       col = colorRamp2(c(0, 0.5, 1), c("white", "mediumblue", "darkblue")),
                       row_title = "Fraction Outliers",
                       row_title_gp = gpar(fontsize = 10))

h2 <- draw(frac_heat_2, heatmap_legend_side = "left")
tiff(minus_heatmap, units="in", width=6, height=5, res=600)
h2
dev.off()



