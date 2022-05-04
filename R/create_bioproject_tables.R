
#----Setup----------------------------------------------------------------------
library(here)
library(tidyverse)
source(here("R", "utils.R"))
species_info <- read_tsv(
    here("data", "rnaseq_samples_count.tsv"), show_col_types = FALSE
)

# When was the last update?
date_file <- here("data", "date.rda")
date <- "2000/01/01"
if(file.exists(date_file)) {
    load(date_file)
}
pdats <- paste0(date, ":2025/01/01")

#----Create bioproject table for Brassicaceae-----------------------------------
brassicaceae_projects <- NULL
if(file.exists(here("data", "brassicaceae_projects.rda"))) {
    load(here("data", "brassicaceae_projects.rda"))
}
species_df_brassicaceae <- species_info %>%
    filter(Family == "Brassicaceae")

retmax <- 50000
brassicaceae_old <- brassicaceae_projects

## Delete later
brassicaceae_new1 <- get_bp_table(species_df_brassicaceae$Species[1], pdats)
brassicaceae_new2 <- get_bp_table(species_df_brassicaceae$Species[2:5], pdats)
brassicaceae_new3 <- get_bp_table(species_df_brassicaceae$Species[6:20], pdats)
brassicaceae_new4 <- get_bp_table(species_df_brassicaceae$Species[21:36], pdats)

brassicaceae_new <- rbind(
    brassicaceae_new1, brassicaceae_new2,
    brassicaceae_new3, brassicaceae_new4
)
##################

brassicaceae_new <- get_bp_table(species_df_brassicaceae$Species, pdats)

brassicaceae_projects <- rbind(
    brassicaceae_old, brassicaceae_new
)

#----Create bioproject table for Poaceae-----------------------------------
poaceae_projects <- NULL
if(file.exists(here("data", "poaceae_projects.rda"))) {
    load(here("data", "poaceae_projects.rda"))
}
species_df_poaceae <- species_info %>%
    filter(Family == "Poaceae")

retmax <- 50000
poaceae_old <- poaceae_projects

## Delete later
poaceae_new1 <- get_bp_table(species_df_poaceae$Species[1], pdats)
poaceae_new2 <- get_bp_table(species_df_poaceae$Species[2], pdats)
poaceae_new3 <- get_bp_table(species_df_poaceae$Species[3], pdats)
poaceae_new4 <- get_bp_table(species_df_poaceae$Species[4], pdats)
poaceae_new5 <- get_bp_table(species_df_poaceae$Species[5], pdats)
poaceae_new6 <- get_bp_table(species_df_poaceae$Species[6], pdats)
poaceae_new7 <- get_bp_table(species_df_poaceae$Species[7], pdats)
poaceae_new8 <- get_bp_table(species_df_poaceae$Species[8], pdats)
poaceae_new9 <- get_bp_table(species_df_poaceae$Species[9], pdats)
poaceae_new10 <- get_bp_table(species_df_poaceae$Species[10], pdats)
poaceae_new11 <- get_bp_table(species_df_poaceae$Species[11:48], pdats)

poaceae_new <- rbind(
    poaceae_new1, poaceae_new2, poaceae_new3, poaceae_new4,
    poaceae_new5, poaceae_new6, poaceae_new7, poaceae_new8,
    poaceae_new9, poaceae_new10, poaceae_new11
)
###############

poaceae_new <- get_bp_table(species_df_poaceae$Species, pdats)

poaceae_projects <- rbind(
    poaceae_old, poaceae_new
)

#----Create bioproject table for Fabaceae---------------------------------------
fabaceae_projects <- NULL
if(file.exists(here("data", "fabaceae_projects.rda"))) {
    load(here("data", "fabaceae_projects.rda"))
}
species_df_fabaceae <- species_info %>%
    filter(Family == "Fabaceae")

retmax <- 30000
fabaceae_old <- fabaceae_projects

# Delete later
fabaceae_new1 <- get_bp_table(species_df_fabaceae$Species[1], pdats)
fabaceae_new2 <- get_bp_table(species_df_fabaceae$Species[2:4], pdats)
fabaceae_new3 <- get_bp_table(species_df_fabaceae$Species[5:30], pdats)

fabaceae_new <- rbind(
    fabaceae_new1, fabaceae_new2, fabaceae_new3
)

####################

fabaceae_new <- get_bp_table(species_df_fabaceae$Species, pdats)

fabaceae_projects <- rbind(
    fabaceae_old, fabaceae_new
)

#----Create bioproject table for all other families-----------------------------
otherfam_projects <- NULL
if(file.exists(here("data", "otherfam_projects.rda"))) {
    load(here("data", "otherfam_projects.rda"))
}
species_df_otherfam <- species_info %>%
    filter(!Family %in% c("Poaceae", "Fabaceae", "Brassicaceae")) %>%
    arrange(-RNAseq_samples)

retmax <- 30000
otherfam_old <- otherfam_projects

# Delete later
otherfam_new1 <- get_bp_table(species_df_otherfam$Species[1], pdats)
otherfam_new2 <- get_bp_table(species_df_otherfam$Species[2], pdats)
otherfam_new3 <- get_bp_table(species_df_otherfam$Species[3], pdats)
otherfam_new4 <- get_bp_table(species_df_otherfam$Species[4], pdats)
otherfam_new5 <- get_bp_table(species_df_otherfam$Species[5], pdats)
otherfam_new6 <- get_bp_table(species_df_otherfam$Species[6], pdats)
otherfam_new7 <- get_bp_table(species_df_otherfam$Species[7], pdats)
otherfam_new8 <- get_bp_table(species_df_otherfam$Species[8], pdats)
otherfam_new9 <- get_bp_table(species_df_otherfam$Species[9], pdats)
otherfam_new10 <- get_bp_table(species_df_otherfam$Species[10], pdats)
otherfam_new11 <- get_bp_table(species_df_otherfam$Species[11], pdats)
otherfam_new12 <- get_bp_table(species_df_otherfam$Species[12], pdats)
otherfam_new13 <- get_bp_table(species_df_otherfam$Species[13:48], pdats)

otherfam_new <- rbind(
    otherfam_new1, otherfam_new2, otherfam_new3, otherfam_new4,
    otherfam_new5, otherfam_new6, otherfam_new7, otherfam_new8,
    otherfam_new9, otherfam_new10, otherfam_new11, otherfam_new12,
    otherfam_new3
)


######################

otherfam_new <- get_bp_table(species_df_otherfam$Species, pdats)

otherfam_projects <- rbind(
    otherfam_old, otherfam_new
)

#----Save updated tables--------------------------------------------------------
brassicaceae_projects <- dplyr::distinct(brassicaceae_projects, .keep_all = TRUE)
fabaceae_projects <- dplyr::distinct(fabaceae_projects, .keep_all = TRUE)
poaceae_projects <- dplyr::distinct(poaceae_projects, .keep_all = TRUE)
otherfam_projects <- dplyr::distinct(otherfam_projects, .keep_all = TRUE)

save(brassicaceae_projects,
     file = here("data", "brassicaceae_projects.rda"),
     compress = "xz")

save(poaceae_projects,
     file = here("data", "poaceae_projects.rda"),
     compress = "xz")

save(fabaceae_projects,
     file = here("data", "fabaceae_projects.rda"),
     compress = "xz")

save(otherfam_projects,
     file = here("data", "otherfam_projects.rda"),
     compress = "xz")

# Finally, update date
date <- gsub("-", "/", Sys.Date())
save(date,
     file = here::here("data", "date.rda"),
     compress = "xz")





