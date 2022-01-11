library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(sctransform)

proj_id <- c('P18953_2001', 'P23359_1001', 'P23359_1002', 'P23359_1003')

# load
source('src/load_BD_scTCR.R')
glob_path <- '/Users/tan/OneDrive - KI.SE/TCR_processed_data/single cell/*/*/*'
raw_tcr <- lapply(proj_id, BD_load_VDJ, dir_path = glob_path) %>% do.call(what = rbind)
sample_tag <- lapply(proj_id, BD_load_sample_tag, dir_path = glob_path) %>% do.call(what = rbind)
raw_gene <- lapply(proj_id, BD_load_gene_exp, dir_path = glob_path, norm_method = 'DBEC') %>% do.call(what = rbind)

# process by seurat
raw_gene_merge <- raw_gene %>% 
  left_join(sample_tag, by = 'unique_index') %>% 
  left_join(raw_tcr %>% select(unique_index, TCR_Paired_Chains, at_least_one_chain, is_gdT, proj_id), by = 'unique_index') %>%
  mutate(Sample_Name = case_when(
    is.na(Sample_Name) ~ proj_id,
    TRUE ~ Sample_Name
  ))
data_gene <- raw_gene_merge %>% 
  filter(!Sample_Name %in% c('Multiplet', 'Undetermined'))
counts_seurat <- data_gene %>% 
  select(-c('Sample_Tag', 'at_least_one_chain', 'proj_id',
            'Sample_Name', 'TCR_Paired_Chains', 'unique_index', 'is_gdT')) %>%
  t()
colnames(counts_seurat) <- data_gene$unique_index

meta_seurat <- data_gene %>%
  select('Sample_Name', 'TCR_Paired_Chains', 'at_least_one_chain', 'is_gdT', 'unique_index') %>%
  as.data.frame()
rownames(meta_seurat) <- data_gene$unique_index

obj <- CreateSeuratObject(counts = counts_seurat,
                          meta.data = meta_seurat)
obj <- SCTransform(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, verbose = FALSE)
g1 <- DimPlot(obj, label = TRUE) + NoLegend()
g2 <- DimPlot(obj, group.by='TCR_Paired_Chains')
g3 <- DimPlot(obj, group.by='at_least_one_chain')
g4 <- DimPlot(obj, group.by='is_gdT')
(g1 + g2) / (g3 + g4)
ggsave(filename = 'umap_all_with_healthy.pdf')
FeaturePlot(obj, features = c('CD3E', 'CD4', 'CD8A'))
ggsave(filename = 'CD4_CD8_expression_with_healthy.pdf')

# manual CD4T/CD8T discrimination
df <- obj[[]] %>%
  as_tibble() %>%
  mutate(cell_type = case_when(
    seurat_clusters %in% c('1', '2', '3', '7', '8', '9', '14', '15') ~ 'CD4T',
    seurat_clusters %in% c('4', '5', '13', '11', '6', '10') ~ 'CD8T',
    is_gdT ~ 'gdT',
    TRUE ~ 'others'
  )) %>%
  cbind(obj@reductions$umap[[]])
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point()
ggsave(filename = 'CD4_CD8_types_with_healthy.pdf')

df %>% 
  select(unique_index, cell_type, UMAP_1, UMAP_2, seurat_clusters) %>%
  write_csv(file = 'data/cell_types_with_healthy.csv')


