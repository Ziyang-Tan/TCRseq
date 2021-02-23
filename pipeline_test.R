library(immunarch)

path = 'data/P18853_1001_result.clonotypes.TRB.txt'
immdata = repLoad(path)





### 
#
raw =mm.fastread('/Users/tan/TCRseq/test_pipeline_with_public_dataset/data/2017/Mtb single cell sequenced TCR clones.xlsx')
dat <- raw %>% 
  filter(!is.na(CDR3beta)) %>% 
  select(CDR3beta, Vbeta, Jbeta, CDR3alpha, BetaReads, Donor, Stim) %>% 
  mutate(subject=paste0(Donor, ':', 'Stim')) %>% 
  select(CDR3beta, Vbeta, Jbeta, CDR3alpha, subject, BetaReads)

write_delim(dat, 'Test_2017_cdr3_file.txt', delim = '\t')
