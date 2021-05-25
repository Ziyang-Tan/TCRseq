library(bigrquery)
sql <- "SELECT * FROM `cradle-259115.HLA.typing`"

tb <- bq_project_query('cradle-259115', sql)
HLA <- bq_table_download(tb)

library(dplyr)
library(tibble)
library(readr)
HLA0201 <- HLA %>% 
  filter(grepl('^02:01', A1) | grepl('^02:01',A2))
#write_csv(HLA0201, 'HLA_0201.csv')

HLA_frequency <- function(HLA_dat, allele){
  # output a frequency table showing that in how many subjects a HLA type appears in 
  # at least one of the alleles
  # HLA_dat is HLA typing table from bigQuery
  # allele is a vector of allele names
  nameList <- HLA_dat %>% select(one_of(allele)) %>% unlist(use.names = F) %>% unique()
  table <- lapply(nameList, function(x){
    tmp <- HLA_dat %>% filter(get(allele[1]) == x | get(allele[2]) == x) %>% count()
    return(tmp)
  }) %>% bind_rows() %>% add_column(HLA = nameList) %>% arrange(desc(n))
  return(table)
}

DRB <- HLA_frequency(HLA0201, c('DRB11', 'DRB12'))
DPB <- HLA_frequency(HLA0201, c('DPB11', 'DPB12'))
DQB <- HLA_frequency(HLA0201, c('DQB11', 'DQB12'))

# select the samples with A02:01:01 B07:02:01 and DRB15:01:01
HLASelect <- HLA %>% 
  filter(grepl('^02:01', A1) | grepl('^02:01',A2)) %>%
  filter(grepl('^15:01', DRB11) | grepl('^15:01', DRB12)) %>%
  filter(grepl('^07:02', B1) | grepl('^07:02', B2))
write_csv(HLASelect, 'HLA_select.csv')



