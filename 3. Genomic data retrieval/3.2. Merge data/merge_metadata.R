#!/usr/bin/env Rscript

# data of 3754 individuals from Wohns et al., 2022 database
all <- read.csv("ind_list.csv")

setwd("Master-Thesis-MBC/2. Selection of human populations and genomic samples/Environmental divergence index/selected_samples.csv")
# data of 330 individuals selected for this study
samples <- read.csv("selected_samples.csv", sep = "\t")
# merge both data 
merged <- merge(all, samples, by.x = "FINAL")

setwd("Master-Thesis-MBC/3. Genomic data retrieval/3.2. Merge data/")
# write a table with the merged data (287 individuals)
write.csv(merged, "merged.csv")
