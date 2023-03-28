## template script for DILA Shiny App

## load libraries
library("Seurat")
library("dplyr")
library("tibble")
library("ggplot2")
library("patchwork")
library("plotly")
library("ggpubr")
# library("tidyverse")

## file path to wherever the RDS object is
file_path <- "~/Documents/tmp/202204_lodi/test_data/"

## Alternative to AddModuleScore that calculates the signature without creating
## a new Seurat object. This is just a simple average, so probably only for 
## prototyping purposes...
CalculateSignatureScore = function (
  object, 
  features, 
  assay = "RNA",
  slot = "data"
) {
  
  # Get the data from the Seurat object
  data = GetAssayData(object,
                      slot = slot,
                      assay = assay)[features, ]
  
  # If using log-normalized data then perform average in non-log space
  if (slot == "data")
    data = expm1(data)
  
  # Calculate average per cell
  data = data %>% Matrix::colMeans(data)
  
  # Return named vector with the signature score
  return(data)
}

## read RDS Seurat object
# setwd("OneDrive - KU Leuven/sc-dev/PanCancer2_app/")
rds = readRDS(file = paste0(file_path, "pancancer_1000.rds"))
# rds <- SetIdent(object = rds, value = "CellType")
rds@meta.data %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename(ID = names(.)[1]) -> meta_data
sort(unique(meta_data$CellType_lev6))
# sort(unique(meta_data$TumorType))

## gene list
gene_list = c("CD4", "CD8A", "CD3D", "CD3E", "CD2")

read.csv(file = paste0(file_path, "test_genes.csv")) %>%
          dplyr::rename(ID = names(.)[1]) %>% as.list() -> a1

## calculate scores
as.data.frame(CalculateSignatureScore(rds, gene_list)) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(ID = names(.)[1],
                sig_score = names(.)[2]) %>%
  dplyr::left_join(meta_data, by = "ID") %>%
  dplyr::rename(cell_id = ID,
                sample_id = `orig.ident`) %>%
  dplyr::mutate(sample_id = as.factor(sample_id),
                CellType_lev5 = as.factor(CellType_lev5)) -> df1

# df = s@meta.data
# df$SampleID = factor(df$SampleID)
# df$CellType_lev5 = factor(df$CellType_lev5)
# 
# ## summarize scores by cell type
df1 %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(mean_sig_score = mean(sig_score)) %>%
  dplyr::arrange(desc(mean_sig_score)) -> df2

## plot scores by cell type
df2 %>%
  dplyr::mutate(CellType = factor(CellType, levels = unique(CellType))) %>%
  ggplot(aes(x = mean_sig_score, y = CellType, fill = CellType)) +
  geom_bar(stat = "identity") +
  coord_flip(expand = TRUE) +
  theme_classic() +
  ylab("Cell Types") +
  xlab("Mean Expression Score") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d() +
  ggtitle("Mean Expression Score by Cell Type") -> p1
ggplotly(p1)

## correlation of signature score vs percent cell type
# df1 %>%
#   dplyr::group_by(orig.ident) %>% 
#   dplyr::summarise(mean_score = mean(scores)) %>%
#   dplyr::arrange(desc(mean_score)) %>% 
#   dplyr::left_join((df1 %>% dplyr::select())) -> df3

## Sam plots

target_cell = "Fibroblast" # Name of cell type of interest
score = "nCount_RNA" # This should be the name of the signature score
ncol = 3

# Some columns should be factors
# df = s@meta.data
# df$SampleID = factor(df$SampleID)
# df$CellType_lev5 = factor(df$CellType_lev5)

# Calculate the frequencies
df1 %>%
  dplyr::group_by(sample_id, TumorType3, CellType_lev5, .drop = FALSE) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n/sum(n)) %>%
  dplyr::filter(CellType_lev5 == target_cell) -> res1

# Calculate the average score per sample and merge with the other dataframe
df1 %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise_at(vars(sig_score), mean) %>%
  merge(res1, by = "sample_id") %>%
  dplyr::filter(CellType_lev5 == target_cell) %>%
  as.data.frame() -> res2

# Create the separate plots
# (ggplot(data = res2, aes(x = freq, y = sig_score)) +
#         geom_point(size = 0.5, aes(color = TumorType3)) +
#         stat_smooth(method = "lm", se = TRUE,
#                     formula = y ~ poly(x, 1, raw = TRUE),
#                     aes(color = TumorType3, fill = TumorType3)) +
#         ggpubr::stat_cor(method = "spearman", label.y = 1) +
#         facet_wrap(paste0("~", "TumorType3"), ncol = ncol, scales = "fixed") +
#         theme_linedraw() +
#         theme(panel.grid = element_blank(),
#               strip.text.x = element_text(size = 16),
#               axis.title = element_text(size = 20),
#               axis.title.y = element_blank()) +
#         labs(color = "Tumor type", y = score,
#              x = paste0(target_cell, " proportion"))) -> p_tumortypes

# Create the main plot of all cancer types
ggplot(data = res2, aes(x = freq, y = sig_score)) +
  geom_point(size = 1, aes(color = TumorType3)) +
  stat_smooth(method = "lm", se = TRUE, fill = "gray", color = "black",
              formula = y ~ poly(x, 1, raw = TRUE)) +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  # labs(y = score, x = paste0(target_cell, " proportion")) -> p_all
  labs(y = "Signature score", x = paste0(target_cell, " proportion")) -> p_all
ggplotly(p_all)

# Combine the two plots into a single figure panel
# p_total <- ((p_all | p_tumortypes) &
#          plot_annotation(tag_levels = "A")) +
#   plot_layout(heights = c(2, 5)) &
#   theme(legend.position = "none")
