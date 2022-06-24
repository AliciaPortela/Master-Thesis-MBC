# charge functions 
source("Master-s-Thesis-/GEAM/Pre-GEAM_filter_sig._SNPs_and_add_environmental_variables/1.Filtering_sig_SNPs_add_env_var_functions.R")

###################################
## 1. FILTERING SIGNIFICANT SNPs ##
###################################

# open GRoSS results 
setwd("Master-s-Thesis-/Outlier_test/Output_GRoSS")

h1 <- read.table("human_chr12_1.tsv", header=T, sep="\t")
h2 <- read.table("human_chr12_2.tsv", header=T, sep="\t")

# min p-value for each genetic position (SNP)
# create an empty vector 
mins <- numeric() 
for(i in 1:nrow(h1)){
  # save the min p-value for a given position (SNP)
  mins[i] <- min(h1[i,38:71]) 
print(i)
}

mins2 <- numeric()
for(i in 1:nrow(h2)){
  mins2[i] <- min(h2[i,38:71])
  print(i)
}

# add min p-value to the previous tables only with columns one ("CHR"=nº chromosome) and two ("START"=SNP position)
h1f <- as.data.frame(cbind(h1[,1], h1[,2], mins))
h2f <- as.data.frame(cbind(h2[,1], h2[,2], mins2))
colnames(h1f) <- colnames(h2f) <- c("CHR", "START", "min_Pvalue")

# plot histogram of p-values for Scree plot
x11()
hist(h1f$min_Pvalue, breaks=5000, xlim=c(0,0.01))
hist(h2f$min_Pvalue, breaks=5000, xlim=c(0,0.01))

# number of SNPs with p-value < 10E-4
sum(h1f$min_Pvalue<0.00001)
sum(h2f$min_Pvalue<0.00001)

# number of SNPs with p-value < 10E-5
sum(h1f$min_Pvalue<0.000001)
sum(h2f$min_Pvalue<0.000001)

# P-value=10E-5 was selected as cut-off for giving the lowest number of variables (significant SNPs)

# save only significant SNPs (with P-value < 10E-5
h1s <- h1[which(mins < 0.00001), ]
h2s <- h2[which(mins2 < 0.00001), ]

# filter allele frequency files by significant SNPs ()
# enter the directory where allele frequency files (.afreq) are
# (each file corresponds to a population --> 18 in total)
setwd("Master-s-Thesis-/GEAM/Pre-GEAM_filter_sig._SNPs_and_add_environmental_variables/Allele_freq_all_SNPs/")

# create output empty lists. Each element of the lists will correspond to the filtering  
# of allele frequency tables of each population (fr) by significant SNPs (those in h1s and h2s)
frh1_list <- list(); frh2_list <- list()

for(i in 1:length(list.files(pattern="afreq"))){

  fr <- read.table(list.files(pattern="afreq")[[i]])
  
  # Guardo los resultados como elementos de las listas-output (=una población)
  frh1_list[[i]] <- filtering(fr, h1s)
  frh2_list[[i]] <- filtering(fr, h2s)
  
  print(i)
}

# convert previous results from list to table
# first column = genetic position (SNP) (the same in all the elements of the lists)
frh1_final <- frh1_list[[1]][,1]; frh2_final <- frh2_list[[1]][, 1]
for(i in 1:length(frh1_list)){
  
  # paste the second column (allelic frequencies) of each element of the list (each population) 
  # to the first column (genetic position)
  frh1_final <- as.data.frame(cbind(frh1_final, frh1_list[[i]][, 2]))
  colnames(frh1_final)[i + 1] <- colnames(frh1_list[[i]])[2]
  
  frh2_final <- as.data.frame(cbind(frh2_final, frh2_list[[i]][, 2]))
  colnames(frh2_final)[i + 1] <- colnames(frh2_list[[i]])[2] 
}

##################################
## 2. ADDING ENVIRONMENTAL DATA ##
##################################

# to join them with environmental data, in which population names appear in first column, 
# transpose genetic data tables and 
# establish population names in the first column 
rownames(frh1_final) <- frh1_final[,1]; frh1_final <- frh1_final[, -1] # genetic positions (SNPs) as rownames to make the transposition
gen1 <- as.data.frame(t(frh1_final)) # transpose 
gen1 <- as.data.frame(cbind("pop" = rownames(gen1), gen1)) # as first column population names 
rownames(gen1) <- NULL # remove rownames 

rownames(frh2_final) <- frh2_final[,1]; frh2_final <- frh2_final[, -1] # genetic positions (SNPs) as rownames to make the transposition
gen2 <- as.data.frame(t(frh2_final)) # transpose
gen2 <- as.data.frame(cbind("pop" = rownames(gen2), gen2)) # as first column population names 
rownames(gen2) <- NULL # remove rownames

# open environmental data 
setwd("Master-s-Thesis-/GEAM/Pre-GEAM_filter_sig._SNPs_and_add_environmental_variables/Environmental data/")
# gen.env = environment of populations in 1000genomes database 
load("3.environment_1000genomes_populations.RData")
# sim.env = environment of populations in Simons database
load("3.environment_simons_populations.RData") 

# change 3 letter codes by complete population names in genomic 
# tables to be able to join to the environmental ones by merge function 
corr <- read.table("correspondence_names_populations.txt", header=T, sep="\t")

gen1 <- merge(gen1, corr, by.x="pop", by.y="Code")
gen1 <- as.data.frame(cbind(gen1$pop.y, gen1[, 2 : 259])) # column with complete pop names the first and join with genetic position data 
colnames(gen1)[1] <- "pop" # rename the column with complete names 

gen2 <- merge(gen2, corr, by.x="pop", by.y="Code")
gen2 <- as.data.frame(cbind(gen2$pop.y, gen2[, 2 : 273])) # column with complete pop names the first and join with genetic position data 
colnames(gen2)[1] <- "pop" # rename the column with complete names

# check if all populations are in both databases 
all(gen1$pop %in% gen.env$pop)
all(gen1$pop %in% sim.env$pop)

# not all populations are in Simons database
# use 1000genomes database
gen1 <- merge(gen1, gen.env, by="pop")
gen2 <- merge(gen2, gen.env, by="pop")

# save complete databases with significant SNPs and environmental data 
write.csv(gen1, "1.gen1.csv", row.names=FALSE)
write.csv(gen2, "1.gen2.csv", row.names=FALSE)

# there are too much variables to include in a PCA 
# make 1000 random subsamples of 5% of the total of SNPs
(5 * (ncol(gen1)-6)) / 100 # 5% of the gen1 SNPs are aprox. 13 SNPs by table
(5 * (ncol(gen2)-6)) / 100 # 5% of the gen 2 SNPs are aprox. 14 SNPs by table

# create output empty lists
gs1 <- list(); gs2 <- list()
for(i in 1 : 1000){
  
  # remove from tables all that do not correspond to genetic positions 
  # and make the subsample 
  gs1[[i]] <- gen1[, -c(1, 260 : 264)][, sample(1 : ncol(gen1[, -c(1, 260 : 264)]), 13)] 
  gs2[[i]] <- gen2[, -c(1, 274 : 278)][, sample(1 : ncol(gen2[, -c(1, 274 : 278)]), 14)]
  
  # add again environmental variables 
  gs1[[i]] <- as.data.frame(cbind(gen1[, 262:264], gs1[[i]]))
  gs2[[i]] <- as.data.frame(cbind(gen2[, 276:278], gs2[[i]]))
  
  # add X to the name of genetic position (easier work with them as a character
  # and not as an integer)
  colnames(gs1[[i]])[-c(1 : 3)] <- paste("X", colnames(gs1[[i]])[-c(1 : 3)], sep="")
  colnames(gs2[[i]])[-c(1 : 3)] <- paste("X", colnames(gs2[[i]][-c(1 : 3)]), sep="")
  
  # save each subsample as an element of the output lists
  gs1[[i]] <- as.data.frame(apply(gs1[[i]], 2, as.numeric))
  gs2[[i]] <- as.data.frame(apply(gs2[[i]], 2, as.numeric))
  
}

# save random tables (input of GEAM)
setwd("Master-s-Thesis-/GEAM/Pre-GEAM_filter_sig._SNPs_and_add_environmental_variables/Databases_sig_SNPs_env_data/")
save(gs1, file = "1000_random_SNPs_gen1.RData")
save(gs2, file = "1000_random_SNPs_gen2.RData")
