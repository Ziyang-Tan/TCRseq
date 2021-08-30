library(dplyr)
library(readr)
library(ggplot2)

path_tcr <- file.path('scTCR and targeted mRNA (BD)', 'data', 'P18953_2001_VDJ_perCell.csv')
path_gene <- file.path('scTCR and targeted mRNA (BD)', 'data', 'P18953_2001_DBEC_MolsPerCell.csv')

data_tcr <- read_csv(path_tcr, skip = 6)
data_gene <- read_csv(path_gene, skip = 6)

dat <- data_tcr %>% 
  select(Cell_Index, 
         TCR_Alpha_Gamma_V_gene_Dominant, 
         TCR_Beta_Delta_V_gene_Dominant, 
         TCR_Paired_Chains) %>%
  left_join(data_gene, by = 'Cell_Index') %>%
  filter(TCR_Paired_Chains) %>%
  mutate(type=case_when(
    str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}") == 'TRGV' ~ 'gdT',
    str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}") == 'TRAV' ~ 'abT',
    TRUE ~ 'unknown'
    ))

ggplot(dat, aes(x = CD4, y = CD8A)) +
  geom_point()

dat %>% filter(CD8A >= 9) %>% dim()
dat %>% filter(CD4 >= 4) %>% dim()
dat %>% filter(CD8A<9 & CD4<4) %>% dim()


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

