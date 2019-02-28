#!/Applications/anaconda3/pkgs/r-base-3.4.3-h9b62496_0/bin/Rscript

library(tidyverse)
library(plyr)
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
              help="index number of row where gene values begins", metavar="character"),
  make_option(c("-c", "--column_gene"), type="integer", default=1, 
              help="index number of column where gene values begins", metavar="character"),
  make_option(c("-m", "--meta_column"), type="integer", default=1, 
              help="index number of column where gene IDs are located", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="output file name for all outputs, do not include file extension", metavar="character"),
  make_option(c("-a", "--aes"), type="character", default="#B2182B", 
              help="aesthetic of distribution plot: color", metavar="character"),
  make_option(c("-u", "--upper_x_lime"), type="integer", default="20", 
              help="upper x limit for distribution plot", metavar="character"),
  make_option(c("-l", "--lower_x_lim"), type="integer", default="-20", 
              help="lower x limit for distribution plot", metavar="character"),
  make_option(c("-t", "--tag_phospho"), type="character", default="", 
              help="tag for phospho data: write 'phospho' if the data is phosphoproteome 
              and 'not_phospho' if it isn't", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


### Required Functions
# moveme ("HanjoStudy/quotidieR")
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}


# test script:
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-phosphoproteome-ratio-norm-NArm.gct -r 27 -c 24 -m 23 -o xxxxxxxx -t phospho
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-proteome-ratio-norm-NArm.gct -r 27 -c 18 -m 17 -o xxxxxxxx -t not_phospho
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-rnaseq-linear-gene-fpkm-uq-log10-NArm-row-norm.gct -r 22 -c 6 -m 4 -o xxxxxxxx -t not_phospho
# Rscript luad_analysis_pipeline_test.R -f ./luad_original_data/luad-v1.0-wgs-cnv-somatic.luad.all.hg38-v2.0-20190121.gct -r 22 -c 6 -m 2 -o xxxxxxxx -t not_phospho

### Establish names 
out_dataframe <- paste(opt$o, ".tsv", sep="")
out_names <- paste(opt$o,"_names.txt", sep="")
out_meta_data <- paste(opt$o,"_meta.tsv", sep="")
out_distribution <- paste(opt$o, "_distribution.tiff", sep = "")
out_outlier_plus <- paste(opt$o,"_outlier_plus.tsv", sep="")
out_outlier_minus <- paste(opt$o,"_outlier_minus.tsv", sep="")
out_outlier_both <- paste(opt$o,"_outlier_both.tsv", sep="")
ggplot_title <- paste(opt$o, "Distribution", sep=" ")

### Tidy Data

# Isolate only the gene ID and sample columns
# Export this as a new dataframe 
df <- read.table(opt$f, skip = 2, header = TRUE, sep = "\t", fill = TRUE)
colnames(df)[opt$m] <- "GeneSymbol"
df_gene <- df[c(opt$r:nrow(df)), c(1, opt$m, opt$c:ncol(df))]
write.table(df_gene, file=out_dataframe, quote=FALSE, sep='\t', col.names = NA)

# Get a list of sample names for the outlier analysis
df_names <- colnames(df_gene[3:ncol(df_gene)])
write.table(df_names, out_names, sep=',',row.names=F, col.names = F)

# Isolate meta data
df_meta <- df[c(1:((opt$r)-1)), c(1, opt$m, opt$c:ncol(df))]
write.table(df_meta, file=out_meta_data, quote=FALSE, sep='\t', col.names = NA)

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
df_stats <- df_gene[,4:ncol(df_gene)]
df_genesymbol <- dplyr::select(df_gene, GeneSymbol)

########



####### CHANGE GENESYMBOL IN COMMAND INPUT


#######

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
count_not_outliers <-  function(x){ length(which(x==0))}
count_outliers <- function(x){ length(which(x==1))}

# Plus
df_final_outlier_plus <- dplyr::select(df_outlier_gene_plus, -count)
colnames(df_final_outlier_plus) <- paste(colnames(df_final_outlier_plus), "not.outlier", sep = ".")
colnames(df_final_outlier_plus)[1] <- "GeneSymbol"
df_final_not_outlier_plus <- dplyr::select(df_outlier_gene_plus, -count)
colnames(df_final_not_outlier_plus) <- paste(colnames(df_final_not_outlier_plus), "outlier", sep = ".")
colnames(df_final_not_outlier_plus)[1] <- "GeneSymbol.not"

# Minus
df_final_outlier_minus <- dplyr::select(df_outlier_gene_minus, -count)
colnames(df_final_outlier_minus) <- paste(colnames(df_final_outlier_minus), "not.outlier", sep = ".")
colnames(df_final_outlier_minus)[1] <- "GeneSymbol"
df_final_not_outlier_minus <- dplyr::select(df_outlier_gene_minus, -count)
colnames(df_final_not_outlier_minus) <- paste(colnames(df_final_not_outlier_minus), "outlier", sep = ".")
colnames(df_final_not_outlier_minus)[1] <- "GeneSymbol.not"

agg_phospho <- function(df_outlier_agg, df_not_outlier_agg, tag){
  if (tag == 'phospho'){
    # Plus
    df_final_outlier_agg <- aggregate(. ~ GeneSymbol, data = df_outlier_agg, 
                                       count_outliers, na.action = na.pass)
    df_final_not_outlier_agg <- aggregate(. ~ GeneSymbol.not, data = df_not_outlier_agg, 
                                           count_not_outliers, na.action = na.pass)
    df_final_outlier <- cbind(df_final_outlier_agg, df_final_not_outlier_agg)
    df_final_outlier <- dplyr::select(df_final_outlier, -GeneSymbol.not)
    df_final_outlier$GeneSymbol <- as.factor(as.character(df_final_outlier$GeneSymbol))
    return(df_final_outlier)
  } else 
    df_final_outlier <- cbind(df_outlier_agg, df_not_outlier_agg)
    df_final_outlier <- dplyr::select(df_final_outlier, -GeneSymbol.not)
    return(df_final_outlier)
}


write.table(df_outlier_gene_plus, file = 'yyyy.tsv', quote=FALSE, sep='\t', col.names = NA)


# Run aggregation function
df_final_plus <- agg_phospho(df_final_outlier_plus ,df_final_not_outlier_plus, opt$t)
df_final_minus <- agg_phospho(df_final_outlier_minus ,df_final_not_outlier_minus, opt$t)

# Write out outlier dataframes
write.table(df_final_plus, file = out_outlier_plus, quote=FALSE, sep='\t', col.names = NA)
write.table(df_final_minus, file = out_outlier_minus, quote=FALSE, sep='\t', col.names = NA)
#write.table(df_outlier_both_agg, out_outlier_both, quote=FALSE, sep='\t', col.names = NA)


### Group Comparison

# Add Count column filled with NAs to the end of the meta dataframe
#df_meta$Count <- NA

# Combine outlier dataframes with meta data
#meta_outlier_plus <- rbind.fill(df_meta, df_outlier_plus_agg)
#meta_outlier_minus <- rbind.fill(df_meta, df_outlier_minus_agg)
#meta_outlier_both <- rbind.fill(df_meta, df_outlier_both_agg)

# Call group dataframes



