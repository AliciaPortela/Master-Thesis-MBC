#!/usr/bin/env Rscript

all <- read.csv("ind_list.csv")

setwd("Master-Thesis-MBC/2. Selection of human populations and genomic samples/Environmental divergence index/selected_samples.csv")
samples <- read.csv("selected_samples.csv", sep = "\t")
merged <- merge(all, samples, by.x = "FINAL")

setwd("Master-Thesis-MBC/3. Genomic data retrieval/3.2. Merge data/")
write.csv(merged, "merged.csv")
