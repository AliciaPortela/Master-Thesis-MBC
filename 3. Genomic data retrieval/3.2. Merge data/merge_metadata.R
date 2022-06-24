#!/usr/bin/env Rscript

all <- read.csv("ind_list.csv")

setwd("2. Selection of human populations and genomic samples/Environmental divergence index/selected_samples.csv")
samples <- read.csv("selected_samples.csv", sep = "\t")

setwd("2. Selection of human populations and genomic samples/Environmental divergence index/selected_samples.csv")
merged <- merge(all, samples, by.x = "FINAL")

write.csv(merged, "merged.csv")
