library(phyloseq)
library(tidyverse)
library(magrittr)
library(readxl)
library(Biostrings)
library(ape)
library(micInt)
ss_physeq <- readRDS(file = 'selection_switch/reactor_swapped.rds')
# We do not want to keep the original phylogentic tree as we do not know how
# it was made
ss_tree <- read.tree(file = 'selection_switch/ss_sina_nj_tree.nwk')
# ss_tree$edge.length <- rep(1,ss_tree$edge %>% nrow())
phy_tree(ss_physeq) <- ss_tree
phyloseq::plot_tree(ss_physeq,color = 'Phylum')
ss_physeq %<>% {prune_taxa(taxa = taxa_sums(.) != 0,x =.)}
ss_sam_data <- sample_data(ss_physeq)
ss_sam_data$From.experiment <- NULL
sample_data(ss_physeq) <- ss_sam_data
ss_ref_seq <- Biostrings::readDNAStringSet(filepath = 'selection_switch/reference_sequences.fac')
ss_ref_seq %<>% extract(names(.) %in% taxa_names(ss_physeq))
ss_physeq@refseq <- ss_ref_seq
# Writes the reference sequences to a file for construction of a new phylogenetic tree
writeXStringSet(ss_ref_seq,filepath = 'selection_switch/ss_filtered_sequences.fasta')
writeXStringSet(ss_ref_seq %>% RNAStringSet(),
                filepath = 'selection_switch/ss_filtered_rna_sequences.fasta')
# Saving the finished phyloseq object prior to further analysis
saveRDS(object = ss_physeq,file = 'selection_switch/ss_physeq.rds',compress = 'xz')
