# read the table 
Frecuencies <- read.csv("allele_freq_database_IS_LD")

# keep SNPs whose allele freq is > 0.05 in all populations  
Frecuencies <- Frecuencies[(Frecuencies$CLM > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$PUR > 0.05),]

Frecuencies <- Frecuencies[(Frecuencias$FRN > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$SAR > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$ADY > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$RUS > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$ESN > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$MZB > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$HZP > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$KLS > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$MAN > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$MSL > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$IBS > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$BAS > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$PEL > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$PTH > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$CEU > 0.05),]

Frecuencies <- Frecuencies[(Frecuencies$MXL > 0.05),]

# write new table only with SNPs whose allele frequency is 0.05<allele freq<0.95
write.csv(Frecuencies, "Path/allele_freq_database_IS_LD_MAF")
