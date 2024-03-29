library(ggplot2)
library(ggalluvial)
library(dplyr)
library(ggrepel)
library(wesanderson)

clone_expansion_donut <- function(sample_name, clone_exp){
  #sample_name = 'ISAC99_1' # for debug purpose
  clone_exp <- clone_exp %>%
    mutate(clone_id = case_when(
      clone_count == 1 ~ 'unique',
      clone_count >1 ~ clone_id
    )) %>%
    arrange(desc(clone_count))
  df <- clone_exp %>% 
    filter(Sample_Name == sample_name) %>%
    arrange(clone_count)
  df <- rbind(tibble(CDR3_concat='unique', 
                     Sample_Name = sample_name, 
                     clone_count = dim(df%>%filter(clone_id=='unique'))[1], 
                     clone_id = 'unique'),
              df %>% filter(clone_id != 'unique')) %>%
    mutate(label = factor(case_when(
      clone_id == 'unique' ~ 'Unique',
      clone_id != 'unique' & clone_count == 2 ~ '2',
      clone_id != 'unique' & clone_count == 3 ~ '3',
      clone_id != 'unique' & clone_count == 4 ~ '4',
      clone_id != 'unique' & clone_count == 5 ~ '5',
      clone_id != 'unique' & clone_count > 5 & clone_count <=9 ~ '6-9',
      clone_id != 'unique' & clone_count > 9 & clone_count <= 19 ~ '10-19',
      clone_id != 'unique' & clone_count > 19 ~ '>=20'
    ), levels=c('Unique', '2', '3', '4', '5', '6-9', '10-19', '>=20')),
    fraction = clone_count/sum(clone_count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    )
  breaks <- df %>% group_by(label) %>% summarise(label_frac=sum(fraction)) %>% 
    mutate(max = cumsum(label_frac), 
           min = c(0,head(max,n=-1)),
           pos = (max+min)/2
    )
  # make a named vector for color
  cols = c('#FFFFFF', wes_palette("Rushmore1", nlevels(breaks$label)-1, type = "continuous"))
  names(cols) = levels(breaks$label)
  g <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.7, fill=label)) +
    geom_rect(color = 1, size=0.2) +
    annotate("text", x=2, y=0, label = sum(df$clone_count), size = 8) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    scale_y_continuous(breaks = breaks$pos, labels = breaks$label) +
    scale_fill_manual(values = cols) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          axis.text = element_text(size = 12), 
          legend.position = "none",
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = sample_name)
  return(g)
}

get_top_expansion_id <- function(clone_exp, top_number=10){
  top_clone_id <- clone_exp %>% 
    group_by(clone_id) %>% 
    summarise(max_count = max(clone_count)) %>% 
    slice_max(order_by = max_count, n = top_number) %>% 
    select(clone_id) %>% 
    unlist(use.names = F)
  return(top_clone_id)
}

clone_expansion_alluvium<- function(individual_name, clone_exp, top_number=10){
  clone_exp <- clone_exp  %>%
    tidyr::separate(Sample_Name, into = c('individual', 'time_point')) %>%
    filter(individual == individual_name)
  top_clone_id <- get_top_expansion_id(clone_exp, top_number)
  df <- clone_exp %>%
    group_by(time_point) %>%
    mutate(clone_ratio = clone_count/sum(clone_count)) %>%
    filter(clone_id %in% top_clone_id) 
  g <- ggplot(df, aes(x = time_point, stratum = clone_id, alluvium = clone_id, y= clone_ratio,fill = clone_id)) +
    geom_alluvium() +
    geom_stratum(size = 0.1) +
    #theme(legend.position = "none")+
    labs(title = individual_name)
  return(g)
}


