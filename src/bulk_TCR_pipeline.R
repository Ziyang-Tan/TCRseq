library(dplyr)
library(readr)
library(tibble)
source("src/clone_expansion_plots.R")

sc_data <- read_csv(file = 'data/scTCR_data_merge.csv')
clone_exp <- read_csv(file = 'data/clone_expansion.csv')

proj_id <- 'P23556'
dir_path <- '/Users/tan/OneDrive - KI.SE/TCR_processed_data/bulk'

sample_info <- read_tsv(Sys.glob(paste0(dir_path, '/', proj_id, '/*_sample_info.txt')), show_col_types = FALSE) %>%
  rename(sample_id = `NGI ID`, Sample_Name = `User ID`) %>%
  mutate(Sample_Name = gsub('-', '_', Sample_Name))

paths <- Sys.glob(paste0(dir_path, '/*/*/*', proj_id, '*.TRB.txt'))
data <- lapply(paths, function(x){
  sample_id = strsplit(x, split='/')[[1]][8]
  read_tsv(x, show_col_types = FALSE) %>%
    add_column(sample_id = sample_id)
}) %>% do.call(what = rbind) %>%
  left_join(sample_info, by='sample_id')



top_clone_id <- get_top_expansion_id(clone_exp, 10)
top_clone_seq <- clone_exp %>% 
  filter(clone_id %in% top_clone_id) %>% 
  select(CDR3_concat, clone_id) %>% 
  unique() %>%
  left_join(sc_data %>% 
              filter(clone_id %in% top_clone_id) %>% 
              select(CDR3aa_concat, clone_id) %>% 
              unique(), by='clone_id') %>%
  tidyr::separate(CDR3_concat, into = c('nt_p1', 'nt_p2')) %>%
  tidyr::separate(CDR3aa_concat, into = c('aa_p1', 'aa_p2'))

data %>% filter(grepl('ASSGGLGQVELRNTIY', aaSeqCDR3))

#set_sc = unique(c(top_clone_seq$aa_p1, top_clone_seq$aa_p2))
#set_bulk = unique(data$aaSeqCDR3)
#intersect(set_sc, set_bulk)







