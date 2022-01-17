library(dplyr)
library(readr)
source("src/clone_expansion_plots.R")
library(ggplot2)
library(Seurat)
library(SeuratDisk)

sc_data <- read_csv(file = 'data/scTCR_data_merge.csv') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv')
obj <- LoadH5Seurat('data/seurat_results.h5Seurat')

top_clone_id <- get_top_expansion_id(clone_exp, 10)

data <- sc_data %>% filter(clone_id %in% top_clone_id)
#obj <- ScaleData(obj)

obj_exp <- subset(obj, subset = unique_index %in% data$unique_index)
obj_exp <- ScaleData(obj_exp)

FeaturePlot(obj_exp, features = c('CD3E', 'CD4', 'CD8A', 'ITGAE'))

tmp <- obj_exp[[]] %>% left_join(data, by = 'unique_index')
obj_exp$clone_id <- tmp$clone_id
Idents(obj_exp) <- 'clone_id'
all_features <- FindAllMarkers(obj_exp)
top_features <- all_features %>% 
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC)
DimPlot(obj_exp, group.by = 'clone_id')

tmp <- GetAssayData(obj_exp, slot = 'scale.data')
d <- hclust(dist(tmp[unique(top_features$gene),], method = 'euclidean'), method = 'average')
DoHeatmap(obj_exp, group.by = 'clone_id', features = c(d$labels[d$order]))
ggsave('figures/top_clone_phenotype.pdf', width = 15, height = 10, units = 'in')


