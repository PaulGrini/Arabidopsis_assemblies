#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("usage: col2syn.R <file.gff> <file.collinearity>", call.=FALSE)
}

require(tidyverse)

file.gff <- args[1]
file.collinearity <- args[2]
file.synteny <- str_replace(file.collinearity, "collinearity", "synteny")
  

collinearity <- read_tsv(file.collinearity, 
                         comment = "#",
                         col_names = c("link_id", "gene_a", "gene_b", "e_value")) %>%
  mutate(link_id = str_remove_all(link_id, " |:")) %>% 
  separate(link_id, c("block","gene"), sep = "-", remove = FALSE, convert = TRUE)  

gff <- read_tsv(file.gff,
                col_names = c("chr", "gene", "start", "end")) %>% 
  mutate(chr_n = str_extract(chr, "\\d+") %>% as.integer())


synteny <- collinearity %>% 
  group_by(block) %>% 
  summarise(n_genes = n(),
            first_gene_a = first(gene_a),
            first_gene_b = first(gene_b),
            last_gene_a = last(gene_a),
            last_gene_b = last(gene_b),
            #max_e = max(e_value),
            #median_e = median(e_value)
  ) %>% 
  left_join(select(gff, c("chr_n","start", "gene")), by = c("first_gene_a" = "gene")) %>% 
  left_join(select(gff, c("end", "gene")), by = c("last_gene_a" = "gene")) %>% 
  left_join(select(gff, c("chr_n","start", "gene")), by = c("first_gene_b" = "gene")) %>% 
  left_join(select(gff, c("end", "gene")), by = c("last_gene_b" = "gene")) %>% 
  transmute(Species_1 = chr_n.x,
            Start_1 = start.x,
            End_1 = end.x,
            Species_2 = chr_n.y,
            Start_2 = start.y,
            End_2 = end.y,
            fill = "cccccc")

write_tsv(synteny, path = file.synteny)