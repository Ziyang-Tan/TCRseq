library(bigrquery)
sql <- "SELECT * FROM `cradle-259115.HLA.typing`"

tb <- bq_project_query('cradle-259115', sql)
HLA <- bq_table_download(tb)

library(dplyr)
library(readr)
tmp <- HLA %>% 
  select(SampleID, A1, A2) %>% 
  filter(grepl('^02:01', A1) | grepl('^02:01',A2))
write_csv(tmp, 'HLA_0201.csv')
