
#----Set up---------------------------------------------------------------------
library(tidyverse)
library(rvest)
library(bears)
library(here)
source(here("R", "utils.R"))

#----Create a data frame of included species and their info---------------------
# Only species with >50% of complete BUSCOs will be included
species_info <- readr::read_tsv(
    "https://raw.githubusercontent.com/almeidasilvaf/pg_genomedata/master/data/genome_info_table.tsv",
    show_col_types = FALSE
) %>%
    filter(BUSCO > 50)

rnaseq_samples <- unlist(lapply(species_info$Species, count_available_samples))
species_info$RNAseq_samples <- rnaseq_samples
species_info <- species_info %>%
    filter(RNAseq_samples > 0) %>%
    arrange(Family, -RNAseq_samples) %>%
    select(Family, Species, BUSCO, RNAseq_samples)

readr::write_tsv(species_info, file = here("data", "rnaseq_samples_count.tsv"))
