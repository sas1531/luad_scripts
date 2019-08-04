

### Shaleigh Smith
# Final figure script


# Load libraries
library(tidyverse)
library(dplyr)
library(rDGIdb)
library(plyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyr)
library(RCircos)
library(ggstance)
library(org.Hs.eg.db)
library(biomaRt)


# Modified data only includes:
##### up-regulated comparisons
##### enriched genes for the mutated (1) population
##### updated names (i.e. cna_egfr, phos_alk, etc.)


# Read in data
data_10_kinase <- read.table("./outliers_qvalues_mutations_frac_filter_0.1.csv", header = TRUE, sep = ",")
colnames(data_10_kinase)[1] <- "gene"

# Fill NA with zero and remove decimals from ensembl gene ids
data_10_kinase[is.na(data_10_kinase)] <- 0
data_10_kinase$gene <- as.character(data_10_kinase$gene)
data_10_kinase$gene <- ifelse(nchar(data_10_kinase$gene) > 15, substring(data_10_kinase$gene, 1, 15), data_10_kinase$gene)

# Check against original data frame
test <- filter(data_10_kinase, gene == "EGFR")
test


# Read in kinase list
kinase <- read.table("./kinase_list.txt", header = FALSE, sep = "\t")
colnames(kinase)[1] <- "gene"
kinase$gene <- as.character(kinase$gene)

# Extract ensembl gene ids for kinases
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
kinase_ensembl <- getBM(attributes='ensembl_gene_id', 
                        filters = 'hgnc_symbol', 
                        values = kinase$gene, 
                        mart = ensembl)
colnames(kinase_ensembl)[1] <- "gene"
kinase_list <- kinase$gene
kinase_ensembl_list <- kinase_ensembl$gene

#kinase_master <- rbind(kinase, kinase_ensembl)
#kinase_master$gene
#master_gene_list <- kinase_master$gene

# Filter data frame for kinase list (hugo and ensembl ids)
data_10_kinase_filter <- data_10_kinase[data_10_kinase$gene %in% kinase_list,]
data_10_kinase_ensembl_filter <- data_10_kinase[data_10_kinase$gene %in% kinase_ensembl_list,]


# Match hugo genes with ensemble genes
gene_2 <- getBM(attributes='hgnc_symbol', 
                filters = 'ensembl_gene_id', 
                values = data_10_kinase_ensembl_filter$gene, 
                mart = ensembl)

# Add hugo genes to ensemble list then remove ensemble column
data_10_kinase_ensembl_filter$gene_2 <- gene_2$hgnc_symbol
data_10_kinase_ensembl_filter <- data_10_kinase_ensembl_filter[, 2:ncol(data_10_kinase_ensembl_filter)]
ensembl_col <- colnames(data_10_kinase_ensembl_filter)
data_10_kinase_ensembl_filter <- dplyr::select(data_10_kinase_ensembl_filter, gene_2, ensembl_col)
colnames(data_10_kinase_ensembl_filter)[1] <- "gene"

# Create master kinase filtered dataframe with hugo gene symbols 
data_10_kinase_master <- rbind(data_10_kinase_ensembl_filter, data_10_kinase_filter)
data_10_kinase_master <- aggregate(. ~ gene, data = data_10_kinase_master, sum)
data_10_kinase_master

### NEG RNAI CRISPR FILTER
#data_10_kinase_master <- data_10_kinase_master[data_10_kinase_master$gene %in% rnai_crispr_gene_list,]
#data_10_kinase_master

# Remove rows that only contain zeros
data_10_kinase_master <- data_10_kinase_master[ rowSums(data_10_kinase_master[,2:ncol(data_10_kinase_master)])!=0, ]
data_10_kinase_master


# Change order and add empty columns
data_10_kinase_master <- dplyr::select(data_10_kinase_master, gene,
                                       phos_egfr, prot_egfr, rna_egfr, cna_egfr,
                                       phos_kras, prot_kras, rna_kras, cna_kras,
                                       phos_tp53, prot_tp53, rna_tp53, cna_tp53,
                                       phos_stk11, prot_stk11, rna_stk11, cna_stk11,
                                       phos_keap1, prot_keap1, rna_keap1, cna_keap1,
                                       phos_alk, prot_alk, rna_alk, cna_alk)

# Extract gene list and make it into the data frame row names
gene_list <- as.data.frame(data_10_kinase_master$gene)
row.names(data_10_kinase_master) <- data_10_kinase_master$gene
data_10_kinase_master <- dplyr::select(data_10_kinase_master, -gene)

# Filter for FDR cutoff of 0.01
data_10_kinase_0.01 <- data_10_kinase_master
data_10_kinase_0.01[data_10_kinase_0.01 > 0.01] <- 0
data_10_kinase_0.01 <- data_10_kinase_0.01[rowSums(data_10_kinase_0.01)!=0, ]
data_10_kinase_0.01

# Remove rows where phos doesn't have a q value 
# every gene in heatmap will have at least one value in a phos column
data_10_kinase_0.01 <- data_10_kinase_0.01[rowSums(data_10_kinase_0.01[,c(1,5,9,13,17,21)])!=0, ]
data_10_kinase_0.01


# Query DGIdb for kinases

# Filterd
druggable_kinase <- rownames(data_10_kinase_0.01)
length(druggable_kinase) #63

# Query 
result_drug_kinase <- queryDGIdb(druggable_kinase)


## Result summary
# Score is the total number of source databases listing the interaction (including PubMed articles)
summary_drug_kinase <- resultSummary(result_drug_kinase) 
summary_drug_kinase_score <- dplyr::select(summary_drug_kinase, Gene, Score)

# Change to numeric and aggregate for total score
summary_drug_kinase_score$Score <- as.numeric(summary_drug_kinase_score$Score)
summary_drug_kinase_score <- aggregate(. ~ Gene, data = summary_drug_kinase_score, sum)

# Get log score
summary_drug_kinase_score$Log_Score <- log(summary_drug_kinase_score$Score)

# Order rows to match heatmap
target <- as.vector(rownames(data_10_kinase_0.01))
summary_drug_kinase_score <- summary_drug_kinase_score[match(target, summary_drug_kinase_score$Gene),]

# Label row names for heatmap
row.names(summary_drug_kinase_score) <- target

# Replace NA with 0 for plotting
summary_drug_kinase_score$Score[is.na(summary_drug_kinase_score$Score) == TRUE] <- 0
summary_drug_kinase_score$Log_Score[is.na(summary_drug_kinase_score$Log_Score) == TRUE] <- 0

# Confirm filtering
summary_drug_kinase_score


# Import manual interactions then aggregate and separate into 1/0
interactions <- read.delim("./druggable/DGIDB_gene_drug_interactions_20190417_modified_dup_removed.txt", sep = "\t")
interactions <- interactions[interactions$gene_name %in% druggable_kinase ,]
interactions
interactions <- filter(interactions, interaction_types == "inhibitor")
interactions

# Select gene name and interaction type only
interactions <- dplyr::select(interactions, gene_name, interaction_types)

# Aggregate to a list
interactions <- aggregate(interaction_types ~ gene_name, interactions, paste, collapse = ",")

# Order rows to match heatmap
target_int <- as.vector(druggable_kinase)
interactions <- interactions[match(target_int, interactions$gene_name),]
interactions$gene_name_2 <- as.character(row.names(summary_drug_kinase_score))

# Change to Yes and No
interactions$interaction_types <- as.factor(interactions$interaction_types)
interactions$interaction_types <- ifelse(is.na(interactions$interaction_types) == TRUE, "False", "True")

# Select
interactions <- dplyr::select(interactions, gene_name_2, interaction_types)
colnames(interactions)[] <- c("gene_name", "interaction_types")
interactions



# FDA approved drugs 
fda <- read.table("./druggable/druggable_FDA_approved_2019-04-17_list.txt")
fda$FDA_approved <- "True"
colnames(fda)[1] <- c("gene")

# Intersect with druggable genes
fda_int <- as.data.frame(interactions$gene_name)
colnames(fda_int) <- c("gene")
fda_int <- left_join(fda_int, fda, by = "gene")

# Change NA to 
fda_int$FDA_approved[is.na(fda_int$FDA_approved) == TRUE] <- "False"
fda_int


# Create CRISPR and RNAi panels

# RNAi
rnai <- read.table("./lung_mean_RNAi.txt", sep = "\t", fill = TRUE, header = TRUE)
rnai

# Intersect with druggable genes
int_rnai <- as.data.frame(interactions$gene_name)
colnames(int_rnai) <- c("gene")
int_rnai <- left_join(int_rnai, rnai, by = "gene")
colnames(int_rnai)[2] <- "RNAi"

# Creat final data frame 
rnai_df <- as.data.frame(int_rnai$RNAi)
row.names(rnai_df) <- int_rnai$gene
colnames(rnai_df)[1] <- "RNAi"
rnai_df

# CRISPR
crispr <- read.table("./lung_mean_CRISPR.txt", sep = "\t", fill = TRUE, header = TRUE)
crispr

# Intersect with druggable genes
int_crispr <- as.data.frame(interactions$gene_name)
colnames(int_crispr) <- c("gene")
int_crispr <- left_join(int_crispr, crispr, by = "gene")
colnames(int_crispr)[2] <- "CRISPR"

# Create final data frame
crispr_df <- as.data.frame(int_crispr$CRISPR)
row.names(crispr_df) <- int_crispr$gene
colnames(crispr_df)[1] <- "CRISPR"
crispr_df

# Combine CRISPR and RNAi
int_rnai_crispr <- inner_join(int_rnai, int_crispr, by = "gene")
int_rnai_crispr

# Remove values greater than or equal to zero
# We only want to visualize the neg values
int_rnai_crispr_neg <- int_rnai_crispr
int_rnai_crispr_neg <- subset(int_rnai_crispr_neg, int_rnai_crispr_neg[, 2] < 0)
int_rnai_crispr_neg <- subset(int_rnai_crispr_neg, int_rnai_crispr_neg[, 3] < 0)

# Create filter list for data frame 
rnai_crispr_gene_list <- int_rnai_crispr_neg$gene
rnai_crispr_gene_list

# Final filtered rna and cripr data frame
rnai_final <- dplyr::select(int_rnai_crispr_neg, gene, RNAi)
crispr_final <- dplyr::select(int_rnai_crispr_neg, gene, CRISPR)


# Apply neg RNAI CRISPR filter
data_10_kinase_0.01 <- data_10_kinase_0.01[row.names(data_10_kinase_0.01) %in% rnai_crispr_gene_list,]
data_10_kinase_0.01

summary_drug_kinase_score <- summary_drug_kinase_score[row.names(summary_drug_kinase_score) %in% rnai_crispr_gene_list,]
summary_drug_kinase_score

interactions <- interactions[interactions$gene %in% rnai_crispr_gene_list,]
interactions

fda_int <- fda_int[fda_int$gene %in% rnai_crispr_gene_list,]
fda_int


# Create function that removes rows with x amount columns containing values greater than zero 
# This will be used to fill out the data frame and shorten it
zeroes = function(row) {sum(row !=0 ) >= 1}

# Apply zeros function
data_10_kinase_0.01 <-  data_10_kinase_0.01[apply(data_10_kinase_0.01, 1, zeroes),]
data_10_kinase_0.01

# Manually cluster the dataframe phos to cna 
data_10_kinase_0.01[data_10_kinase_0.01 == 0] <- NA

data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$cna_alk, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$rna_alk, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$prot_alk, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$phos_alk, decreasing = FALSE),]

data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$cna_keap1, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$rna_keap1, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$prot_keap1, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$phos_keap1, decreasing = FALSE),]

data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$cna_stk11, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$rna_stk11, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$prot_stk11, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$phos_stk11, decreasing = FALSE),]

data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$cna_tp53, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$rna_tp53, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$prot_tp53, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$phos_tp53, decreasing = FALSE),]

data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$cna_kras, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$rna_kras, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$prot_kras, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$phos_kras, decreasing = FALSE),]

data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$cna_egfr, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$rna_egfr, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$prot_egfr, decreasing = FALSE),]
data_10_kinase_0.01 <- data_10_kinase_0.01[order(data_10_kinase_0.01$phos_egfr, decreasing = FALSE),]


# Replace na with zero for the heatmap
#data_10_kinase_0.01 <- log10(data_10_kinase_0.01)
data_10_kinase_0.01[is.na(data_10_kinase_0.01)] <- 0
data_10_kinase_0.01


### Create final figure 

# Main heatmap annotation
data_drug_kinase_annotation <- HeatmapAnnotation(Comparison = c(rep('EGFR', 4), rep('KRAS', 4), 
                                                                rep('TP53', 4), rep('STK11', 4), 
                                                                rep('KEAP1', 4), rep('EML4.ALK', 4)), 
                                                 Data = c("Phosphoprotein", "Protein", "RNA", "CNV",
                                                          "Phosphoprotein", "Protein", "RNA", "CNV",
                                                          "Phosphoprotein", "Protein", "RNA", "CNV",
                                                          "Phosphoprotein", "Protein", "RNA", "CNV",
                                                          "Phosphoprotein", "Protein", "RNA", "CNV",
                                                          "Phosphoprotein", "Protein", "RNA", "CNV"),
                                                 col = list(Comparison = c(EGFR = "#EF7C12",
                                                                           KRAS = "#FFB90F",
                                                                           TP53 = "#F7E690",
                                                                           STK11 = "#EFC7E6",
                                                                           KEAP1 = "#CB87B4",
                                                                           EML4.ALK = "#B25D91"),
                                                            Data = c(Phosphoprotein = "#27408B",
                                                                     Protein = "#0076C0" ,
                                                                     RNA = "#54BCD1",
                                                                     CNV = "#AFDFEF")),
                                                 annotation_legend_param = list(legend_height = unit(1, "cm"),
                                                                                legend_width = unit(0.20, "cm"),
                                                                                title_gp = gpar(fontsize = 11),
                                                                                labels_gp = gpar(fontsize = 9)),
                                                 name = "data_heat_anno")

# Heatmap q values
data_heat_drug_kinase <- Heatmap(data_10_kinase_0.01,
                                 name = "data_heat",
                                 top_annotation = data_drug_kinase_annotation,
                                 show_column_names = F, 
                                 show_row_names = T,
                                 cluster_rows = F,
                                 cluster_columns = F,
                                 row_names_gp = gpar(fontsize = 9),
                                 top_annotation_height = unit(1, "cm"),
                                 heatmap_legend_param = list(legend_height = unit(1, "cm"),
                                                             legend_width = unit(0.20, "cm"),
                                                             title = "Q Value",
                                                             title_gp = gpar(fontsize = 11),
                                                             labels_gp = gpar(fontsize = 9)),
                                 col = colorRamp2(c(0, 1.0e-50, 0.0001, 0.01), 
                                                  c("gray90", "red3", "orangered2", "oldlace")),
                                 row_title_gp = gpar(fontsize = 10),
                                 rect_gp = gpar(col = "white", lwd = 0.75),
                                 width = unit(13, "cm"))

# Create heatmap legend
data_heat_drug_kinase_fig <- draw(data_heat_drug_kinase, heatmap_legend_side = "right")

# Barplot
row_anno_kinase <- rowAnnotation(barplot = row_anno_barplot(summary_drug_kinase_score$Log_Score,
                                                            axis_direction = "reverse",
                                                            baseline = 0,
                                                            border = TRUE, 
                                                            axis = TRUE, 
                                                            axis_side = "top",
                                                            gp = gpar(fill = "black")), 
                                 width = unit(2, "cm"))


# Colors for interactions (druggable) panel
colors_3 <-  structure(c("darkorchid4", "gray90"), 
                       names = c("True", "False"))

# Create interactions heatmap
interactions_drug_kinase_heat <- Heatmap(interactions$interaction_types,
                                         show_column_names = F, 
                                         show_row_names = F,
                                         cluster_rows = F,
                                         cluster_columns = F,
                                         row_names_gp = gpar(fontsize = 0),
                                         column_names_gp = gpar(fontsize = 2),
                                         top_annotation_height = unit(1, "cm"),
                                         show_heatmap_legend = TRUE,
                                         heatmap_legend_param = list(legend_height = unit(1, "cm"),
                                                                     legend_width = unit(0.20, "cm"),
                                                                     title = "Druggable",
                                                                     title_gp = gpar(fontsize = 11),
                                                                     labels_gp = gpar(fontsize = 9)),
                                         col = colors_3,
                                         row_title_gp = gpar(fontsize = 10),
                                         rect_gp = gpar(col = "white", lwd = 0.75),
                                         name = "data_interactions",
                                         width = unit(0.5, "cm"))

# Colors for FDA panel
colors_2 <-  structure(c("#0099D5", "gray90"), 
                       names = c("True", "False"))


# Create FDA heatmap
fda_drug_kinase_heat <- Heatmap(fda_int$FDA_approved,
                                show_column_names = F, 
                                show_row_names = F,
                                cluster_rows = F,
                                cluster_columns = F,
                                row_names_gp = gpar(fontsize = 0),
                                column_names_gp = gpar(fontsize = 2),
                                top_annotation_height = unit(1, "cm"),
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = list(legend_height = unit(1, "cm"),
                                                            legend_width = unit(0.20, "cm"),
                                                            title = "FDA Approved",
                                                            title_gp = gpar(fontsize = 11),
                                                            labels_gp = gpar(fontsize = 9)),
                                col = colors_2,
                                row_title_gp = gpar(fontsize = 10),
                                rect_gp = gpar(col = "white", lwd = 0.75),
                                name = "data_fda",
                                width = unit(0.5, "cm"))

# Create color ramps for crispr and rnai
sequence <- seq(from = -0.5, to = 0, by = 0.1)
col_seq <- c("#276419", "#4D9221", "#7FBC41", "#B8E186", "#E6F5D0", "#F7F7F7")

# Create crispr heatmap
crispr_drug_kinase_heat <- Heatmap(crispr_final$CRISPR,
                                   show_column_names = F,
                                   column_names_side = "top",
                                   show_row_names = F,
                                   cluster_rows = F,
                                   cluster_columns = F,
                                   row_names_gp = gpar(fontsize = 0),
                                   column_names_gp = gpar(fontsize = 2),
                                   top_annotation_height = unit(1, "cm"),
                                   show_heatmap_legend = TRUE,
                                   heatmap_legend_param = list(legend_height = unit(1, "cm"),
                                                               legend_width = unit(0.20, "cm"),
                                                               title = "CRISPR",
                                                               title_gp = gpar(fontsize = 11),
                                                               labels_gp = gpar(fontsize = 9)),
                                   col = colorRamp2(c(sequence),
                                                    col_seq),
                                   row_title_gp = gpar(fontsize = 10),
                                   rect_gp = gpar(col = "white", lwd = 0.75),
                                   name = "data_crispr",
                                   width = unit(0.5, "cm"),
                                   na_col = "gray90")

# Create crispr heatmap
rnai_drug_kinase_heat <- Heatmap(rnai_final$RNAi,
                                 show_column_names = F,
                                 column_names_side = "top",
                                 show_row_names = F,
                                 cluster_rows = F,
                                 cluster_columns = F,
                                 row_names_gp = gpar(fontsize = 0),
                                 column_names_gp = gpar(fontsize = 2),
                                 top_annotation_height = unit(1, "cm"),
                                 show_heatmap_legend = TRUE,
                                 heatmap_legend_param = list(legend_height = unit(1, "cm"),
                                                             legend_width = unit(0.20, "cm"),
                                                             title = "RNAi",
                                                             title_gp = gpar(fontsize = 11),
                                                             labels_gp = gpar(fontsize = 9)),
                                 col = colorRamp2(c(sequence),
                                                  col_seq),
                                 row_title_gp = gpar(fontsize = 10),
                                 rect_gp = gpar(col = "white", lwd = 0.75),
                                 name = "data_rnai",
                                 width = unit(0.5, "cm"),
                                 na_col = "gray90")


# create score legend
lgd <- Legend(labels = "Log10(# of References)", 
              legend_gp = gpar(fill = "black"),
              title_gp = gpar(fontsize = 11),
              labels_gp = gpar(fontsize = 9))

# create NA legend
na_lgd <- Legend(labels = "NA", 
                 legend_gp = gpar(fill = "gray90"),
                 title_gp = gpar(fontsize = 11),
                 labels_gp = gpar(fontsize = 9))



# Create pdf and add boxes around heatmap panels 
pdf("updated_1x_neg_phos_kinase_heatmap_0.01_fdr_0.1_filter.pdf", width=12, height=6)

final_druggable <- draw(row_anno_kinase + 
                          rnai_drug_kinase_heat +
                          crispr_drug_kinase_heat +
                          fda_drug_kinase_heat +
                          interactions_drug_kinase_heat + 
                          data_heat_drug_kinase_fig,
                        gap = unit(c(2, 2, 2, 2, 2, 0), "mm"), annotation_legend_list = list(lgd, na_lgd))

decorate_heatmap_body("data_heat", {
  grid.rect(0, 0, width = 4/24, height = 1.08, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_heat", {
  grid.rect(4/24, 0, width = 4/24, height = 1.08, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_heat", {
  grid.rect(8/24, 0, width = 4/24, height = 1.08, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_heat", {
  grid.rect(12/24, 0, width = 4/24, height = 1.08, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_heat", {
  grid.rect(16/24, 0, width = 4/24, height = 1.08, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_heat", {
  grid.rect(20/24, 0, width = 4/24, height = 1.08, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_interactions", {
  grid.rect(0, 0, width = 1, height = 1, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_fda", {
  grid.rect(0, 0, width = 1, height = 1, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_crispr", {
  grid.rect(0, 0, width = 1, height = 1, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

decorate_heatmap_body("data_rnai", {
  grid.rect(0, 0, width = 1, height = 1, just = c("left", "bottom"), 
            gp = gpar(lwd = 0.75, col = "black", fill = "transparent"))
})

dev.off()





