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
  make_option(c("--up"), type="character", default=NULL, 
              help="up-regulated dataset file path from outlier analysis", metavar="character"),
  make_option(c("--down"), type="character", default=NULL, 
              help="down-regulated dataset file path from outlier analysis", metavar="character"),
  make_option(c("--meta"), type="character", default=NULL, 
              help="meta dataset file path from outlier analysis", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="output file prefix for all comparison outputs, do not include file extension", 
              metavar="character"),
  make_option(c("--sig_cutoff"), type="character", default='0.05', 
              help="false discovery rate cutoff as decimal (ex. 0.05)", 
              metavar="character"),
  make_option(c("--comparison_a"), type="character", default="Type", 
              help="Column name that identifies comparison A (e.g. Region, Gender, Type)", 
              metavar="character"),
  make_option(c("--comparison_b"), type="character", default=NULL, 
              help="Column name that identifies comparison B (e.g. Region, Gender, Type)", 
              metavar="character"),
  make_option(c("--comparison_c"), type="character", default=NULL, 
              help="column name that identifies comparison C (e.g. Region, Gender, Type)", 
              metavar="character"),
  make_option(c("--group_1a"), type="character", default=NULL, 
              help="Group one of comparison A", metavar="character"),
  make_option(c("--group_1b"), type="character", default=NULL, 
              help="Group 1 of comparison B", metavar="character"),
  make_option(c("--group_1c"), type="character", default=NULL, 
              help="Group 1 of comparison C", metavar="character"),
  make_option(c("--group_2a"), type="character", default=NULL, 
              help="Group two of comparison A", metavar="character"),
  make_option(c("--group_2b"), type="character", default=NULL, 
              help="Group two of comparison B", metavar="character"),
  make_option(c("--group_2c"), type="character", default=NULL, 
              help="Group two of comparison C", metavar="character"),
  make_option(c("--group_comp"), type="character", default="one", 
              help="Number of groups: comparing one group to one, input 'one' or comparing 
              two groups combined, input 'two', or comparing three groups combined, input 
              'three'", metavar="character"),
  make_option(c("--normal_tumor"), type="character", default="both", 
              help="If your samples include both tumor and normal samples AND you want to 
              isolate tumor samples insert 'both'. If you want to compare normal and tumor 
              samples insert 'both_normal'. If you only have tumor samples insert 'tumor_only"
              , metavar="character"),
  make_option(c("--tumor_column"), type="character", default="Type", 
              help="Name of column that indicates whether or not the sample is a tumor sample 
              or normal sample (only required if there are normal and tumor samples", 
              metavar="character"),
  make_option(c("--tumor_group"), type="character", default="Tumor", 
              help="Name of tumor group in tumor column (only reguired if there are normal 
              and tumor samples", metavar="character"),
  make_option(c("--prot"), type="character", default="no", 
              help="'yes' if the input data is proteomic and includes isoforms with the 
              same gene name", metavar="character"),
  make_option(c("--gene_list"), type="character", default="no", 
              help="Path to gene list: only required if you want to isolate specific genes. 
              Text file format", metavar="character")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$up)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

### Establish names 
plus_sig_output <- paste(opt$o, "_up_", opt$group_1a, "_", opt$group_1b, "_",
                         opt$group_2a, "_", opt$group_2b, "_", opt$group_3a, "_", 
                         opt$group_3b, "_sig_outliers.txt", sep="")
minus_sig_output <- paste(opt$o, "_down_", opt$group_1a, "_", opt$group_1b, "_", 
                          opt$group_2a, "_", opt$group_2b, "_", opt$group_3a, "_", 
                          opt$group_3b,  "_sig_outliers.txt", sep="")
plus_heatmap <- paste(opt$o,"_up_heatmap_sig_", opt$group_1a, "_", opt$group_1b, "_",
                      opt$group_2a, "_", opt$group_2b, "_", opt$group_3a, "_", 
                      opt$group_3b, ".png", sep="")
minus_heatmap <- paste(opt$o,"_down_heatmap_sig_", opt$group_1a, "_", opt$group_1b, 
                       "_", opt$group_2a, "_", opt$group_2b, "_", opt$group_3a, 
                       "_", opt$group_3b,  ".png", sep="")

### Functions

### Fisher Exact function
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



### Import Outlier Analysis Dataframes 
df_final_plus <- read.table(opt$up, sep = "\t", header = TRUE)
df_final_plus <- dplyr::select(df_final_plus, -X)
df_final_minus <- read.table(opt$down, sep = "\t", header = TRUE)
df_final_minus <- dplyr::select(df_final_minus, -X)

### Modify meta data for downstream filtering

# Import meta dataframe
meta_df <- read.table(opt$meta, sep = "\t", header = TRUE)
meta_df <- dplyr::select(meta_df, -X)

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


### Group Comparison

# Isolate Group Samples
comparison_1 <- grep(opt$comparison_a, colnames(meta_df_t))
colnames(meta_df_t)[comparison_1] <- c('comp_a')
sample_1 <- grep(opt$tumor_column, colnames(meta_df_t))
colnames(meta_df_t)[sample_1] <- c('samp')

comparison_2 <- grep(opt$comparison_a, colnames(meta_df_final_t))
colnames(meta_df_final_t)[comparison_2] <- c('comp_a')
sample_2 <- grep(opt$tumor_column, colnames(meta_df_final_t))
colnames(meta_df_final_t)[sample_2] <- c('samp')

if (opt$group_comp == 'two'){
  comparison_3 <- grep(opt$comparison_b, colnames(meta_df_t))
  colnames(meta_df_t)[comparison_3] <- c('comp_b')
  comparison_4 <- grep(opt$comparison_b, colnames(meta_df_final_t))
  colnames(meta_df_final_t)[comparison_4] <- c('comp_b')
} else if (opt$group_comp == 'three'){
  comparison_3 <- grep(opt$comparison_b, colnames(meta_df_t))
  colnames(meta_df_t)[comparison_3] <- c('comp_b')
  comparison_4 <- grep(opt$comparison_b, colnames(meta_df_final_t))
  colnames(meta_df_final_t)[comparison_4] <- c('comp_b')
  comparison_5 <- grep(opt$comparison_c, colnames(meta_df_t))
  colnames(meta_df_t)[comparison_5] <- c('comp_c')
  comparison_6 <- grep(opt$comparison_c, colnames(meta_df_final_t))
  colnames(meta_df_final_t)[comparison_6] <- c('comp_c')
}

if (opt$normal_tumor == 'both' & opt$group_comp == 'one'){
  all_col <- meta_df_t[(((meta_df_t$comp_a == opt$group_1a & meta_df_t$comp_a == opt$group_2a)) &
                          (meta_df_t$samp == opt$tumor_group)), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & 
                                           meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a &
                                               meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a &
                                           meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a &
                                               meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
} else if (opt$normal_tumor == 'both' & opt$group_comp == 'two'){
  all_col <- meta_df_t[(((meta_df_t$comp_a == opt$group_1a & meta_df_t$comp_b == opt$group_1b) | 
                           (meta_df_t$comp_a == opt$group_2a & meta_df_t$comp_b == opt$group_2b)) &
                          (meta_df_t$samp == opt$tumor_group)), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & 
                                           meta_df_final_t$comp_b == opt$group_1b &
                                           meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a &
                                               meta_df_final_t$comp_b == opt$group_1b &
                                               meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a &
                                           meta_df_final_t$comp_b == opt$group_2b &
                                           meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a &
                                               meta_df_final_t$comp_b == opt$group_2b &
                                               meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
} else if (opt$normal_tumor == 'both' & opt$group_comp == 'three'){
  all_col <- meta_df_t[(((meta_df_t$comp_a == opt$group_1a & meta_df_t$comp_b == opt$group_1b & meta_df_t$comp_c == opt$group_1c) | 
                           (meta_df_t$comp_a == opt$group_2a & meta_df_t$comp_b == opt$group_2b & meta_df_t$comp_c == opt$group_2c)) &
                          (meta_df_t$samp == opt$tumor_group)), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & 
                                           meta_df_final_t$comp_b == opt$group_1b &
                                           meta_df_final_t$comp_c == opt$group_1c &
                                           meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a &
                                               meta_df_final_t$comp_b == opt$group_1b &
                                               meta_df_final_t$comp_c == opt$group_1c &
                                               meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a &
                                           meta_df_final_t$comp_b == opt$group_2b &
                                           meta_df_final_t$comp_c == opt$group_2c &
                                           meta_df_final_t$outlier == 'yes' &
                                           meta_df_final_t$samp == opt$tumor_group), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a &
                                               meta_df_final_t$comp_b == opt$group_2b &
                                               meta_df_final_t$comp_c == opt$group_2c &
                                               meta_df_final_t$outlier == 'no' &
                                               meta_df_final_t$samp == opt$tumor_group), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
} else if ((opt$normal_tumor == 'both_normal' | opt$normal_tumor == 'tumor_only') & opt$group_comp == 'one'){
  all_col <- meta_df_t[(meta_df_t$comp_a == opt$group_1a | meta_df_t$comp_a == opt$group_2a), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & meta_df_final_t$outlier == 'yes'), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & meta_df_final_t$outlier == 'no'), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a & meta_df_final_t$outlier == 'yes'), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a & meta_df_final_t$outlier == 'no'), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
} else if ((opt$normal_tumor == 'both_normal' | opt$normal_tumor == 'tumor_only') & opt$group_comp == 'two'){
  all_col <- meta_df_t[((meta_df_t$comp_a == opt$group_1a & meta_df_t$comp_b == opt$group_1b) | 
                          (meta_df_t$comp_a == opt$group_2a & meta_df_t$comp_b == opt$group_2b)), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & 
                                           meta_df_final_t$comp_b == opt$group_1b &
                                           meta_df_final_t$outlier == 'yes'), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a &
                                               meta_df_final_t$comp_b == opt$group_1b &
                                               meta_df_final_t$outlier == 'no'), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a & 
                                           meta_df_final_t$comp_b == opt$group_2b &
                                           meta_df_final_t$outlier == 'yes'), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a & 
                                               meta_df_final_t$comp_b == opt$group_2b &
                                               meta_df_final_t$outlier == 'no'), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
} else if ((opt$normal_tumor == 'both_normal' | opt$normal_tumor == 'tumor_only') & opt$group_comp == 'three'){
  all_col <- meta_df_t[((meta_df_t$comp_a == opt$group_1a & meta_df_t$comp_b == opt$group_1b & meta_df_final_t$comp_c == opt$group_1c) | 
                          (meta_df_t$comp_a == opt$group_2a & meta_df_t$comp_b == opt$group_2b & meta_df_final_t$comp_c == opt$group_2c)), ]
  all_col <- as.vector(all_col[,1])
  outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a & 
                                           meta_df_final_t$comp_b == opt$group_1b &
                                           meta_df_final_t$comp_c == opt$group_1c &
                                           meta_df_final_t$outlier == 'yes'), ]
  outlier_group1_col <- (outlier_group1_col[,1])
  not_outlier_group1_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_1a &
                                               meta_df_final_t$comp_b == opt$group_1b &
                                               meta_df_final_t$comp_c == opt$group_1c &
                                               meta_df_final_t$outlier == 'no'), ]
  not_outlier_group1_col <- (not_outlier_group1_col[,1])
  outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a & 
                                           meta_df_final_t$comp_b == opt$group_2b &
                                           meta_df_final_t$comp_c == opt$group_2c &
                                           meta_df_final_t$outlier == 'yes'), ]
  outlier_group2_col <- (outlier_group2_col[,1])
  not_outlier_group2_col <- meta_df_final_t[(meta_df_final_t$comp_a == opt$group_2a & 
                                               meta_df_final_t$comp_b == opt$group_2b &
                                               meta_df_final_t$comp_c == opt$group_2c &
                                               meta_df_final_t$outlier == 'no'), ]
  not_outlier_group2_col <- (not_outlier_group2_col[,1])
}

all_col <- as.vector(all_col[[1]])
outlier_group1_col <- as.vector(outlier_group1_col[[1]])
not_outlier_group1_col <- as.vector(not_outlier_group1_col[[1]])
outlier_group2_col <- as.vector(outlier_group2_col[[1]])
not_outlier_group2_col <- as.vector(not_outlier_group2_col[[1]])

# Create list of dataframes for comparison: upregulated and downregulated
dataframes <- list(df_final_plus, df_final_minus)

for (df_final in dataframes){
  
  # Exit if there are no samples in any group
  if (length(outlier_group1_col) == 0 | length(not_outlier_group1_col) == 0) { 
    print("Not enough observations in Group 1: length = 0")
    next
  } else if (length(outlier_group2_col) == 0 | length(not_outlier_group2_col) == 0) { 
    print("Not enough observations in Group 2: length = 0")
    next
  }
  
  # Set Parameters for plus/minus output 
  if (identical(df_final, df_final_plus) == TRUE){
    color1 <-  "red4"
    color2 <-  "red2"
    heatmap_name <- plus_heatmap
    sig_output <- plus_sig_output
  } else {
    color1 <- "darkblue" 
    color2 <- "mediumblue" 
    heatmap_name <- minus_heatmap
    sig_output <- minus_sig_output
  }
  
  # Plus comparison pipeline 
  # sum all four groups separately
  # put these values in 4 new columns
  # take the values in each of these columns to create a contingecy table
  # run fischer on this table 
  # record the p-value from this table in a new column
  
  # Only include targeted samples
  df_final_sample <- df_final[, c("GeneSymbol", 
                                  (outlier_group1_col),
                                  (outlier_group2_col),
                                  (not_outlier_group1_col),
                                  (not_outlier_group2_col))]
  
  # Outlier group1 
  outlier_group1_col_number <- which(colnames(df_final) %in% (outlier_group1_col))
  outlier_group_1_df <- df_final[outlier_group1_col_number]
  outlier_group_1_df_final <- data.frame(outlier_group1 = apply(outlier_group_1_df, 1, sum, na.rm = TRUE))
  outlier_group_1_df_final$outlier_group1 <- as.numeric(outlier_group_1_df_final$outlier_group1)
  df_final <- cbind(df_final_sample, outlier_group_1_df_final)
  # Not Outlier group1 
  not_outlier_group1_col_number <- which(colnames(df_final) %in% not_outlier_group1_col)
  not_outlier_group_1_df <- df_final[not_outlier_group1_col_number]
  not_outlier_group_1_df_final <- data.frame(not_outlier_group1 = apply(not_outlier_group_1_df, 1, sum, na.rm = TRUE))
  not_outlier_group_1_df_final$not_outlier_group1 <- as.numeric(not_outlier_group_1_df_final$not_outlier_group1)
  df_final <- cbind(df_final, not_outlier_group_1_df_final)
  # Outlier group2
  outlier_group2_col_number <- which(colnames(df_final) %in% outlier_group2_col)
  outlier_group_2_df <- df_final[outlier_group2_col_number]
  outlier_group_2_df_final <- data.frame(outlier_group2 = apply(outlier_group_2_df, 1, sum, na.rm = TRUE))
  outlier_group_2_df_final$outlier_group2 <- as.numeric(outlier_group_2_df_final$outlier_group2)
  df_final <- cbind(df_final, outlier_group_2_df_final)
  # Not Outlier group2
  not_outlier_group2_col_number <- which(colnames(df_final) %in% not_outlier_group2_col)
  not_outlier_group_2_df <- df_final[not_outlier_group2_col_number]
  not_outlier_group_2_df_final <- data.frame(not_outlier_group2 = apply(not_outlier_group_2_df, 1, sum, na.rm = TRUE))
  not_outlier_group_2_df_final$not_outlier_group2 <- as.numeric(not_outlier_group_2_df_final$not_outlier_group2)
  df_final <- cbind(df_final, not_outlier_group_2_df_final)
  df_final <- as.data.frame(df_final)
  
  # Add outlier count column to dataframe - this is the count of number of outlier samples in group 1
  outlier_count <- df_final[, c(outlier_group1_col)]
  
  # Exit if there are no outliers in the outlier count list
  if (ncol(outlier_count) == 0) {
    print("There are no significant values for this comparison.")
    next
  }
  
  df_final$count_outliers <- rowSums((outlier_count[, 2:(ncol(outlier_count))])!=0, na.rm = TRUE)
  df_final$count_outliers <- as.numeric(df_final$count_outliers)
  
  #print(head(df_final))
  
  ### Create outlier fractions (outlier divided by total)
  # Filter so that you have the combined samples (summed) for outlier and not_outlier (use grep)
  # Divide the outlier genes by the combined samples to create outlier fractions
  # Isolate genes that are upregulated in group

  df_names <- as.vector(c(outlier_group1_col))
  df_names <- gsub('.{8}$', '', df_names)
  collapse_1 <- cbind(outlier_group_1_df, not_outlier_group_1_df)
  combine_samples_1 <- sapply(df_names, function(xx) rowSums(collapse_1[,grep(xx, names(collapse_1)), drop=FALSE], na.rm = TRUE))
  combine_samples_1 <- as.data.frame(combine_samples_1)
  df_names <- as.vector(c(outlier_group2_col))
  df_names <- gsub('.{8}$', '', df_names)
  collapse_2 <- cbind(outlier_group_2_df, not_outlier_group_2_df)
  combine_samples_2 <- sapply(df_names, function(xx) rowSums(collapse_2[,grep(xx, names(collapse_2)), drop=FALSE], na.rm = TRUE))
  combine_samples_2 <- as.data.frame(combine_samples_2)
  combine_samples <- cbind(combine_samples_1, combine_samples_2)
  
  #df_names <- as.vector(c(outlier_group1_col, outlier_group2_col))
  #df_names <- gsub('.{8}$', '', df_names)
  #combine_samples <- sapply(df_names, function(xx) rowSums(df_final[,grep(xx, names(df_final)), drop=FALSE], na.rm = TRUE))

  ###
  combine_samples <- as.data.frame(combine_samples)
  outlier_samples_only <- df_final[, c(1, 2:((ncol(combine_samples))+1))]
  fraction_outlier <- cbind(outlier_samples_only[1], round((outlier_samples_only[-1]/combine_samples), 6))
  df_fraction_outlier <- cbind(fraction_outlier, df_final[,(ncol(df_final)-4):ncol(df_final)])
  
  ### Create group 1 fraction mean and group 2 fraction mean columns 
  # Group 1
  fraction_group1_col_number <- which(colnames(df_fraction_outlier) %in% outlier_group1_col)
  fraction_mean_group1 <- df_fraction_outlier[fraction_group1_col_number]
  fraction_mean_group1 <- data.frame(group1_fraction_mean = apply(fraction_mean_group1, 1, mean, na.rm=TRUE))
  fraction_mean_group1[is.na(fraction_mean_group1)] <- 0
  df_fraction_outlier <- cbind(df_fraction_outlier, fraction_mean_group1)
  
  # Group 2
  fraction_group2_col_number <- which(colnames(df_fraction_outlier) %in% outlier_group2_col)
  fraction_mean_group2 <- df_fraction_outlier[fraction_group2_col_number]
  fraction_mean_group2 <- data.frame(group2_fraction_mean = apply(fraction_mean_group2, 1, mean, na.rm=TRUE))
  fraction_mean_group2[is.na(fraction_mean_group2)] <- 0
  df_fraction_outlier <- cbind(df_fraction_outlier, fraction_mean_group2)
  
  # Run fischer exact function
  # Create a p-value column
  fisher_df <- data.frame(p_value = apply(df_fraction_outlier, 1, fish_exact))
  df_fraction_outlier <- cbind(df_fraction_outlier, fisher_df)
  
  # Create a column with the corrected p-values 
  p_adjust_df <- data.frame(lapply(df_fraction_outlier['p_value'], p.adjust, method = "BH"))
  colnames(p_adjust_df)[1] <- c("q_value")
  df_fraction_outlier <- cbind(df_fraction_outlier, p_adjust_df)
  
  # Get significant outliers
  # If p_value_corrected is greater than FDR cutoff, insert yes (significant outlier)
  fdr <- opt$sig_cutoff
  df_fraction_outlier$sig_outlier_q <- NA
  df_fraction_outlier$sig_outlier_q[df_fraction_outlier$q_value <= 0.01] <- "yes"
  df_fraction_outlier$sig_outlier_q[df_fraction_outlier$q_value > 0.01] <- "no"
  
  # If fraction group 1 mean is greater than group 2 fraction mean, insert yes (significant outlier)
  df_fraction_outlier$sig_outlier_mean[(df_fraction_outlier$group1_fraction_mean) > (df_fraction_outlier$group2_fraction_mean)] <- "yes"
  df_fraction_outlier$sig_outlier_mean[df_fraction_outlier$group1_fraction_mean <= df_fraction_outlier$group2_fraction_mean ] <- "no"
  
  ## Length of list of signfiicant outliers will be the number of outliers = print this
  df_fraction_outlier <- filter(df_fraction_outlier,
                                sig_outlier_q == 'yes' & sig_outlier_mean == 'yes')
  df_fraction_outlier <- df_fraction_outlier[order(df_fraction_outlier$q_value, decreasing = FALSE), ]
  
  # Create list of genes and FDR
  if (opt$gene_list == "no"){
    sig_outlier_genes_q_value <- df_fraction_outlier[c(1, ncol(df_fraction_outlier)-2)]
  } else {
    genes <- scan(opt$gene_list, character(), quote = "")
    sig_outlier_genes_q_value <- df_fraction_outlier[df_fraction_outlier$GeneSymbol %in% genes,]
    sig_outlier_genes_q_value <- sig_outlier_genes_q_value[c(1, ncol(df_fraction_outlier)-2)]
  }
  
  if (nrow(sig_outlier_genes_q_value) == 0) {
    print("There are no significant values for this comparison.")
    next
  }
  
  ### Make heat map
  # Define heatmap dataframe
  # Make gene symbol row.names 
  # Isolate values
  # Make values into matrix
  
  
  #fraction_heat_1 <- df_fraction_outlier
  
  if (opt$gene_list == "no"){
    fraction_heat_1 <- df_fraction_outlier
  } else {
    genes <- scan(opt$gene_list, character(), quote = "")
    fraction_heat_1 <- df_fraction_outlier[df_fraction_outlier$GeneSymbol %in% genes,]
  }
  
  fraction_heat <- fraction_heat_1[,-1]
  rownames(fraction_heat) <- fraction_heat_1[,1]
  fraction_heat <- dplyr::select(fraction_heat, -outlier_group1, -not_outlier_group1, -outlier_group2,
                                 -not_outlier_group2, -count_outliers, -group1_fraction_mean, -group2_fraction_mean,
                                 -p_value, -q_value, -sig_outlier_q, -sig_outlier_mean)
  is.na(fraction_heat) <- sapply(fraction_heat, is.infinite)
  fraction_heat[is.na(fraction_heat)]<-0
  fraction_heat[] <- lapply(fraction_heat, function(x) as.numeric(x))
  fraction_heat <- as.matrix(fraction_heat)
  
  ## Make annotation bar (pull from meta data)
  if (opt$normal_tumor == 'both'){
    meta_df_t$comp[meta_df_t$comp_a == ""] <- NA
    meta_df_t$comp[meta_df_t$comp_b == ""] <- NA
    meta_df_t$comp[meta_df_t$comp_c == ""] <- NA
    heat_filter <- meta_df_t[(meta_df_t$samp == opt$tumor_group), ]
  } else if (opt$normal_tumor == 'both_normal' | opt$normal_tumor == 'tumor_only') {
    meta_df_t$comp[meta_df_t$comp_a == ""] <- NA
    meta_df_t$comp[meta_df_t$comp_b == ""] <- NA
    meta_df_t$comp[meta_df_t$comp_c == ""] <- NA
    heat_filter <- meta_df_t
  }
  
  
  #comparison <- grep(opt$comparison, colnames(meta_df_t))
  #colnames(meta_df_t)[comparison] <- c('comp')
  
  if (opt$group_comp == "one"){
    heat_annotation <- heat_filter
    heat_annotation <- heat_annotation[((heat_annotation$comp_a == opt$group_1a) |
                                          (heat_annotation$comp_a == opt$group_2a)), ]
    heat_annotation <- as.data.frame(heat_annotation)
    heat_annotation[, 'comp_a'] <- as.factor(heat_annotation[, 'comp_a'])
    heat_annotation <- as.data.frame(heat_annotation$comp_a)
    colnames(heat_annotation)[1] <- c('comp')
    levels(heat_annotation)[levels(heat_annotation)==opt$group_1a] <- opt$group_1a
    levels(heat_annotation)[levels(heat_annotation)==opt$group_2a] <- opt$group_2a
    ordering <- c(opt$group_1a, opt$group_2a)
    heat_annotation <- as.data.frame(heat_annotation[order(match(heat_annotation$comp, ordering)), ])
    colnames(heat_annotation)[1] <- c('comp')
    heat_color <- c(color1, "gray85")
    heat_annotation$comp <- relevel(heat_annotation$comp, opt$group_1a)
    heat_annotation[, 'comp'] <- factor(heat_annotation[, 'comp'])
    names(heat_color) <- levels(heat_annotation$comp)
    heat_annotation_final <- HeatmapAnnotation(df = data.frame(Comparison = heat_annotation$comp),
                                               col = list(Comparison = heat_color),
                                               annotation_legend_param = list(title_gp = gpar(fontsize = 12),
                                                                              labels_gp = gpar(fontsize = 10)))
  } else if (opt$group_comp == "two") {
    heat_annotation <- heat_filter
    heat_annotation <- heat_annotation[((heat_annotation$comp_a == opt$group_1a & 
                                           heat_annotation$comp_b == opt$group_1b) |
                                          (heat_annotation$comp_a == opt$group_2a &
                                             heat_annotation$comp_b == opt$group_2b)), ]
    heat_annotation$comp <- paste(heat_annotation$comp_a, heat_annotation$comp_b, sep = "_")
    heat_annotation <- as.data.frame(heat_annotation)
    heat_annotation[, 'comp'] <- as.factor(heat_annotation[, 'comp'])
    heat_annotation <- as.data.frame(heat_annotation$comp)
    colnames(heat_annotation)[1] <- c('comp')
    heat_annotation$comp <- as.factor(as.character(heat_annotation$comp))
    comp_a <- paste(opt$group_1a, opt$group_1b, sep = "_")
    comp_b <- paste(opt$group_2a, opt$group_2b, sep = "_")
    levels(heat_annotation)[levels(heat_annotation)=="comp_a"] <- comp_a
    levels(heat_annotation)[levels(heat_annotation)=="comp_b"] <- comp_b
    ordering <- c(comp_a, comp_b)
    heat_annotation <- as.data.frame(heat_annotation[order(match(heat_annotation$comp, ordering)), ])
    colnames(heat_annotation)[1] <- c('comp')
    heat_color <- c(color1, "gray85")
    heat_annotation[, 'comp'] <- factor(heat_annotation[, 'comp'])
    names(heat_color) <- levels(heat_annotation$comp)
    heat_annotation_final <- HeatmapAnnotation(df = data.frame(Comparison = heat_annotation$comp), 
                                               col = list(Comparison = heat_color),
                                               annotation_legend_param = list(title_gp = gpar(fontsize = 12),
                                                                              labels_gp = gpar(fontsize = 10)))
  } else if (opt$group_comp == "three") {
    heat_annotation <- heat_filter
    heat_annotation <- heat_annotation[((heat_annotation$comp_a == opt$group_1a & 
                                           heat_annotation$comp_b == opt$group_1b &
                                           heat_annotation$comp_c == opt$group_1c) |
                                          (heat_annotation$comp_a == opt$group_2a &
                                             heat_annotation$comp_b == opt$group_2b &
                                             heat_annotation$comp_c == opt$group_2c)), ]
    heat_annotation$comp <- paste(heat_annotation$comp_a, heat_annotation$comp_b, sep = "_")
    heat_annotation <- as.data.frame(heat_annotation)
    heat_annotation[, 'comp'] <- as.factor(heat_annotation[, 'comp'])
    heat_annotation <- as.data.frame(heat_annotation$comp)
    colnames(heat_annotation)[1] <- c('comp')
    heat_annotation$comp <- as.factor(as.character(heat_annotation$comp))
    comp_a <- paste(opt$group_1a, opt$group_1b, opt$group_1c, sep = "_")
    comp_b <- paste(opt$group_2a, opt$group_2b, opt$group_2c, sep = "_")
    levels(heat_annotation)[levels(heat_annotation)=="comp_a"] <- comp_a
    levels(heat_annotation)[levels(heat_annotation)=="comp_b"] <- comp_b
    ordering <- c(comp_a, comp_b)
    heat_annotation <- as.data.frame(heat_annotation[order(match(heat_annotation$comp, ordering)), ])
    colnames(heat_annotation)[1] <- c('comp')
    heat_color <- c(color1, "gray85")
    heat_annotation$comp <- relevel(heat_annotation$comp, comp_a)
    heat_annotation[, 'comp'] <- factor(heat_annotation[, 'comp'])
    names(heat_color) <- levels(heat_annotation$comp)
    heat_annotation_final <- HeatmapAnnotation(df = data.frame(Comparison = heat_annotation$comp), 
                                               col = list(Comparison = heat_color),
                                               annotation_legend_param = list(title_gp = gpar(fontsize = 12),
                                                                              labels_gp = gpar(fontsize = 10)))
  }
  
  # Set gene printing options
  if (nrow(fraction_heat) < 50){
    heat_font <- 12
  } else {
    heat_font <- 0
  }
  
  frac_heat_1 <- Heatmap(fraction_heat,
                         top_annotation = heat_annotation_final,
                         show_column_names = F, 
                         show_row_names = T,
                         cluster_rows = F,
                         cluster_columns = F,
                         row_names_gp = gpar(fontsize = heat_font),
                         top_annotation_height = unit(1, "cm"),
                         heatmap_legend_param = list(legend_height = unit(3, "in"),
                                                     title = NULL),
                         col = colorRamp2(c(0, 0.5, 1), c("white", color2, color1)),
                         row_title = "Fraction Outliers",
                         row_title_gp = gpar(fontsize = 12))
  
  h1 <- draw(frac_heat_1, heatmap_legend_side = "left")
  
  # Write out sig genes and q values
  write.table(sig_outlier_genes_q_value, file = sig_output, sep='\t',row.names=F, col.names = T, quote=FALSE)
  # Write out heatmap
  png(heatmap_name, units="in", width=6, height=5, res=600)
  print(h1)
  dev.off()
  
  print("Comparison Analysis Complete")
}