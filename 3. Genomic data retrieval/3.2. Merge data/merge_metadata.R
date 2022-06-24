#!/usr/bin/env Rscript

all <- read.csv("ind_list.csv")

setwd("Master-s-Thesis-/Outlier_test/Output_GRoSS")
samples <- read.csv("selected_samples.csv", sep = "\t")

merged <- merge(all, samples, by.x = "FINAL")

write.csv(merged, "merged.csv")
