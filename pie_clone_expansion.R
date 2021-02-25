library(readr)
library(ggrepel)
library(ggpubr)
library(stringr)

path <- file.path('scTCR and targeted mRNA (BD)', 'data', 'P18953_2001_VDJ_perCell.csv') 
raw <- read_csv(path, skip = 6) # skip the head info

# summary of type
summaryType <- bind_rows(raw %>%
                           filter(TCR_Paired_Chains) %>%
                           select(TCR_Alpha_Gamma_V_gene_Dominant) %>% 
                           filter(!is.na(TCR_Alpha_Gamma_V_gene_Dominant)) %>% 
                           distinct() %>% 
                           mutate(type=str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}")) %>% 
                           group_by(type) %>% 
                           count(),
                         raw %>% 
                           filter(TCR_Paired_Chains) %>%
                           select(TCR_Beta_Delta_V_gene_Dominant) %>% 
                           filter(!is.na(TCR_Beta_Delta_V_gene_Dominant)) %>% 
                           distinct() %>% 
                           mutate(type=str_extract(TCR_Beta_Delta_V_gene_Dominant, "^.{4}")) %>% 
                           group_by(type) %>% 
                           count())

t1 <- ggtexttable(summaryType,rows = NULL, 
                  theme = ttheme("mOrange"))

dat <- raw %>% 
  select(Cell_Index, 
         TCR_Alpha_Gamma_V_gene_Dominant,
         TCR_Beta_Delta_V_gene_Dominant,
         TCR_Paired_Chains)
# pair vs not pair
p1 <- ggplot(dat %>% group_by(TCR_Paired_Chains) %>% count() %>% ungroup() %>%
               arrange(desc(TCR_Paired_Chains)) %>%              
               mutate(prop = n / count(dat)[[1]] *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),           # compute position
             aes(fill=TCR_Paired_Chains,x="", y=prop)) +
  geom_bar(stat = "identity",width = 1, color='white') +
  coord_polar("y",start=0) +
  theme_void() +
  geom_text(aes(y = ypos, label = n), color = "white") +
  labs(title='overall single cells')

# unique pairs
uniquePair <- dat %>% filter(TCR_Paired_Chains) %>% 
  mutate(unique_pair=paste0(TCR_Alpha_Gamma_V_gene_Dominant, 
                            ' and ', 
                            TCR_Beta_Delta_V_gene_Dominant)) %>% 
  group_by(unique_pair) %>% 
  count() %>%
  arrange(desc(n))

p2 <- ggplot(uniquePair, aes(x=n)) + geom_histogram() +
  labs(title='unique pair distribution')

p3 <- ggplot(uniquePair %>% 
               ungroup() %>% 
               filter(n>4) %>% 
               bind_rows(uniquePair %>% 
                           group_by(n) %>% 
                           count() %>% 
                           ungroup() %>% 
                           filter(n<=4) %>% 
                           mutate(unique_pair = paste0(n, ' clone(s)'))%>% 
                           select(unique_pair, nn) %>% 
                           rename(n=nn)) %>%
               arrange(desc(n)) %>%
               mutate(pos=cumsum(n)-0.5*n) %>%  # calculate position
               mutate(label=case_when(n>10 ~ paste0(unique_pair,'(n=',n,')'))), # show labels with n>=10
             aes(x="", y=n, fill=reorder(unique_pair,n))) +
  geom_bar(stat = "identity", width=.01, color='white') +
  geom_text_repel(aes(label = label, y = pos),
                  hjust = "left", 
                  fontface = "bold", 
                  size = 3, 
                  nudge_x = .3, 
                  direction = "y",
                  force = 0.3) +
  coord_polar("y")+
  theme_void() + 
  theme(legend.position = "none")

ggarrange(ggarrange(p1,p2, nrow = 2), ggarrange(p3, t1, nrow=2, heights = c(1,0.5)), nrow=1, ncol=2) %>%
  ggexport(filename='P18953_2001_pie_clone_exansion.pdf')

