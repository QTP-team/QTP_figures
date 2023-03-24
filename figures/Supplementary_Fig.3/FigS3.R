# Supplementary Fig. 3. Taxonomic statistics of 19,251 SGBs.
library(tidyverse)
library(readxl)
library(RColorBrewer)

taxon_comp <- function(infile, outfile, top_num){
  rank_order <-  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  #infile <- "0.data/Tibetan.ar.count.tsv"
  dset <- read_tsv(infile) 
  dset <- dset %>%
    #filter(!is.na(Taxon) & Rank != "Species") %>%
    mutate(Percentage = Counts/dset[1,3][[1]]*100,
           Rank = factor(Rank, levels = rank_order)) %>%
    arrange(Rank, desc(Percentage))
  
  # label other taxa
  ranks_selected = c("Phylum", "Class", "Order", "Family", "Genus")
  Top = top_num
  index = 1
  res_lst <-  list()
  
  for (r in ranks_selected){
    #r = "Genus"
    dset.rank <- dset %>% filter(Rank == r)
    dset.rank.res <- dset.rank %>% filter(Taxon != "Unclassified") %>%
      .[1:Top,] %>%
      mutate(Color = brewer.pal(Top,"Set2")) %>%
      arrange(Percentage)
    
    dset.rank.other <- dset.rank %>% filter(Taxon != "Unclassified") %>%
      anti_join(dset.rank.res, by = "Taxon") 
    
    dset.rank.unclassified <- dset.rank %>% filter(Taxon == "Unclassified") %>%
      mutate(Color = "white")
    dset.rank.unclassified$Taxon <- paste0(r, "__unclassified")
    
    dset.rank.res <- dset.rank.unclassified %>% 
      add_row(Rank = r, Taxon = paste0(r,"__others"), Counts = sum(dset.rank.other$Counts), Percentage = sum(dset.rank.other$Percentage), Color = "grey") %>%
      rbind(dset.rank.res)
    
    res_lst[[index]] <- dset.rank.res
    index = index + 1
  }
  
  plot.dset <- dplyr::bind_rows(res_lst) %>%
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           Rank = factor(Rank, levels=c("Phylum", "Class", "Order", "Family", "Genus")))
  
  # plot stacked plot taxa counts
  SGB_taxon <- print(ggplot(plot.dset, aes(x=Rank, y=Percentage, fill=Taxon)) 
                     + geom_bar(stat="identity", colour="black", alpha=0.5, size=0.2)
                     #+ geom_bar(stat="identity", alpha=0.7, size=0.2)
                     + theme_classic()
                     + ylab("Proportion (%)")
                     + scale_fill_manual(values=as.vector(plot.dset$Color))
                     + scale_y_continuous(expand = c(0, 0), limits = c(0,101))
                     + theme(axis.title.y = element_text(size=6))
                     + theme(axis.text.y = element_text(size=6))
                     + theme(axis.title.x = element_blank())
                     #+ theme(panel.grid =element_blank(), panel.border = element_blank())
                     + theme(legend.position="none")
                     + theme(axis.text.x = element_text(size=6, hjust=0.5, vjust=0.5)))
  ggsave(outfile, SGB_taxon, width = 6, height = 8)
}

taxon_comp("0.data/Tibetan.bac.count.new.tsv", "1.result/FigS3A.Tibetan_bac_taxa.new.pdf", 5)
taxon_comp("0.data/Tibetan.ar.count.new.tsv", "1.result/FigS3B.Tibetan_ar_taxa.new.pdf", 3)
