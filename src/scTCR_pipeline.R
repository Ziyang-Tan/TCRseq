library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
# 'P18953_2001', 
proj_id <- c('P18953_2001', 'P23359_1001', 'P23359_1002', 'P23359_1003')
sample_list = list(ISAC99=c('ISAC99_1','ISAC99_2','ISAC99_3','ISAC99_4','ISAC99_5','ISAC99_6'), 
                   ISAC35=c('ISAC35_1','ISAC35_2','ISAC35_3','ISAC35_4','ISAC35_5','ISAC35_6','ISAC35_7'))

# load data
source('src/load_BD_scTCR.R')
glob_path <- '/Users/tan/OneDrive - KI.SE/TCR_processed_data/single cell/*/*/*'
raw_tcr <- lapply(proj_id, BD_load_VDJ, dir_path = glob_path) %>% do.call(what = rbind)
sample_tag <- lapply(proj_id, BD_load_sample_tag, dir_path = glob_path) %>% do.call(what = rbind)

cell_type <- read_csv('data/cell_types.csv')

raw_tcr_merge <- left_join(raw_tcr, sample_tag, by = 'unique_index') %>%
  left_join(cell_type, by = 'unique_index')

# summarise 
df <- raw_tcr_merge %>% filter(!Sample_Tag %in% c('Multiplet', 'Undetermined'))
table_summary <- df %>% group_by(proj_id) %>% summarise(demultiplexed = n()) %>% left_join(
  df %>% filter(at_least_one_chain) %>% group_by(proj_id) %>% summarise(at_least_one_CDR3 = n())
) %>% left_join(
  df %>% filter(TCR_Paired_Chains) %>% group_by(proj_id) %>% summarise(paired = n())
)

data <- raw_tcr_merge %>%
  filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
  filter(!is.na(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant)) %>%
  filter(!Sample_Tag %in% c('Multiplet', 'Undetermined')) %>%
  mutate(CDR3_concat = paste0(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, '_', 
                              TCR_Beta_Delta_CDR3_Nucleotide_Dominant),
         CDR3aa_concat = paste0(TCR_Alpha_Gamma_CDR3_Translation_Dominant, '_', 
                                TCR_Beta_Delta_CDR3_Translation_Dominant)) %>%
  mutate(clone_id = as.character(as.numeric(as.factor(CDR3_concat))))
clone_id_map <- data %>% select(CDR3_concat, clone_id) %>% unique()
clone_exp <- data %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')

# donut chart 
source("scTCR and targeted mRNA (BD)/clone_expansion_plots.R")
patient<- 'ISAC35'
for (sub_name in c('CD4T', 'CD8T', 'gdT')){
  fig_dir <- file.path('figures', patient)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  clone_exp_sub <- data %>%
    filter(cell_type == sub_name) %>%
    group_by(CDR3_concat, Sample_Name) %>%
    summarise(clone_count = n()) %>%
    ungroup() %>%
    inner_join(clone_id_map, by='CDR3_concat')
  g_list1 <- lapply(sample_list[[patient]],function(x){clone_expansion_donut(x, clone_exp_sub)})
  ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
    ggexport(filename = file.path(fig_dir, paste0(patient, '_clone_expansion_', sub_name, '.pdf')), 
             width = 10, height = 20)
  g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_gdT'))
  ggsave(plot = g, filename = file.path(fig_dir, paste0(patient, '_top_clone_changes_', sub_name, '.pdf')))
}












# alluvium
#g_list2 <- lapply(c('ISAC99', 'ISAC35'), function(x){clone_expansion_alluvium(x, clone_exp)})
#ggarrange(plotlist = g_list2, ncol = 1, nrow = 2) %>%
#  ggexport(filename = 'shared_clones_ratio.pdf', width = 10, height = 20)
# expansion at different visit
data <- data %>%
  left_join(clone_exp, by = c('CDR3_concat', 'Sample_Name')) %>%
  group_by(Sample_Name) %>%
  mutate(clone_ratio = clone_count/n())

df_exp <- data %>%
  select(clone_id, clone_count, clone_ratio, Sample_Name) %>%
  filter(clone_count > 1) %>%
  mutate(clone_ratio_bin = case_when(
    clone_ratio < 0.003 ~ '<0.3%',
    clone_ratio >= 0.003 & clone_ratio < 0.005 ~ '0.3-0.5%',
    clone_ratio >= 0.005 & clone_ratio < 0.025 ~ '0.5%-2.5%',
    clone_ratio >= 0.025 ~ '>=2.5%'
  )) %>%
  mutate(clone_ratio_bin = factor(clone_ratio_bin, levels = c('<0.3%', '0.3-0.5%', '0.5%-2.5%', '>=2.5%'))) %>%
  unique()

ggplot(df_exp, aes(fill=clone_ratio_bin, x=Sample_Name)) + 
  geom_bar(position="fill", stat="count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# typing
tmp = data %>%
  filter(TCR_Paired_Chains)%>%
  select(TCR_Alpha_Gamma_V_gene_Dominant, TCR_Beta_Delta_V_gene_Dominant, Sample_Name) %>% 
  filter(!is.na(TCR_Alpha_Gamma_V_gene_Dominant)) %>% 
  filter(!is.na(TCR_Beta_Delta_V_gene_Dominant)) %>% 
  mutate(type_a=str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}")) %>% 
  mutate(type_b=str_extract(TCR_Beta_Delta_V_gene_Dominant, "^.{4}")) %>% 
  mutate(type = paste0(type_a, '_', type_b))%>%
  group_by(type) %>% 
  count()

# dat <- data_tcr %>% 
#   select(Cell_Index, 
#          TCR_Alpha_Gamma_V_gene_Dominant, 
#          TCR_Beta_Delta_V_gene_Dominant, 
#          TCR_Paired_Chains) %>%
#   left_join(data_gene, by = 'Cell_Index') %>%
#   filter(TCR_Paired_Chains) %>%
#   mutate(type=case_when(
#     str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}") == 'TRGV' ~ 'gdT',
#     str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}") == 'TRAV' ~ 'abT',
#     TRUE ~ 'unknown'
#     ))
# 
# ggplot(dat, aes(x = CD4, y = CD8A)) +
#   geom_point()
# 
# dat %>% filter(CD8A >= 9) %>% dim()
# dat %>% filter(CD4 >= 4) %>% dim()
# dat %>% filter(CD8A<9 & CD4<4) %>% dim()


# export for vdjview
# raw <- read_csv(path, skip = 6) %>%  # skip the head info
#   filter(TCR_Paired_Chains) %>%
#   filter(!(grepl('TRD', TCR_Beta_Delta_V_gene_Dominant) | 
#              grepl('TRG', TCR_Alpha_Gamma_V_gene_Dominant))) # remove gdT
# TRA <- raw %>% 
#   select(Cell_Index, 
#          TCR_Alpha_Gamma_V_gene_Dominant, 
#          TCR_Alpha_Gamma_J_gene_Dominant, 
#          TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, 
#          TCR_Alpha_Gamma_CDR3_Translation_Dominant) %>% 
#   rename(CellID = Cell_Index, 
#          v.segment = TCR_Alpha_Gamma_V_gene_Dominant, 
#          j.segment = TCR_Alpha_Gamma_J_gene_Dominant, 
#          cdr3nt = TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, 
#          cdr3aa = TCR_Alpha_Gamma_CDR3_Translation_Dominant) %>%
#   tibble::add_column(d.segment=NA, 
#                      "Isotype/Constant"=NA, 
#                      "Membrane/Secreted"=NA, 
#                      "Expression"=NA)
# TRB <- raw %>% 
#   select(Cell_Index, 
#          TCR_Beta_Delta_V_gene_Dominant, 
#          TCR_Beta_Delta_D_gene_Dominant, 
#          TCR_Beta_Delta_J_gene_Dominant, 
#          TCR_Beta_Delta_CDR3_Nucleotide_Dominant, 
#          TCR_Beta_Delta_CDR3_Translation_Dominant) %>% 
#   rename(CellID = Cell_Index, 
#          v.segment = TCR_Beta_Delta_V_gene_Dominant, 
#          d.segment = TCR_Beta_Delta_D_gene_Dominant,
#          j.segment = TCR_Beta_Delta_J_gene_Dominant, 
#          cdr3nt = TCR_Beta_Delta_CDR3_Nucleotide_Dominant, 
#          cdr3aa = TCR_Beta_Delta_CDR3_Translation_Dominant) %>%
#   tibble::add_column("Isotype/Constant"=NA, 
#                      "Membrane/Secreted"=NA, 
#                      "Expression"=NA)
# 
# write_csv(TRA, '/Users/tan/TCRseq/VDJView/data/TRA.csv')
# write_csv(TRB, '/Users/tan/TCRseq/VDJView/data/TRB.csv')