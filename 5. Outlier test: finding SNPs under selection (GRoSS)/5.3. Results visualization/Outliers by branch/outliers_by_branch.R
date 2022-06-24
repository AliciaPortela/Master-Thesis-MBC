# import libraries 
library(dplyr)
library(ggplot2)


# set working directory
getwd()
setwd("")
dir()


#### Phylogeny 1 ####

# read the table 
tabletsv1 <- read.csv("human_chr12_11KG_CLM_PUR_FRN_SAR_ADY_RUS_ESN_MZB_HZP_KLS_MAN_MSL_IBS_BAS_PEL_PTH_CEU_MXL.tsv", sep = '\t', header = TRUE)

# check the table dimensions 
head(tabletsv1)
count(tabletsv1)
nrow(tabletsv1)
ncol(tabletsv1)

# find out the number of columns of interest 
which(colnames(tabletsv1)=="Pval_a_r" )
which(colnames(tabletsv1)=="Pval_MZB_p")

# create an empty vector to include the outlier counts by branch
outliers_branch_1 <- c()
outliers_branch_1
typeof(outliers_branch_1)

# iterate over each branch column 
range <- 38:ncol(tabletsv1)
range
length(range)

for (i in range)
{
  # append counts of outliers by branch to the empty vector
  outliers_branch_1[i-37] <- sum(tabletsv1[,i] < 0.00001)
}

# check vector content and length
outliers_branch_1
length(outliers_branch_1)

# find out the number of columns of interest 
which(colnames(tabletsv1)=="a_r" )
which(colnames(tabletsv1)=="MZB_p")

# create a vector with the names of branches 
names_branches_1 <- colnames(tabletsv1)[4:37]
names_branches_1
as.vector(names_branches_1)
class(names_branches_1)

#### create dataframe 1 with data of phylogeny 1 ####
dataframe1 <- data.frame(names_branches_1, outliers_branch_1)
dataframe1

# change the name of columns
colnames(dataframe1)[1] <- "Branches"
colnames(dataframe1)[2] <- "Number_outliers"
dataframe1

# plot data
ggplot(dataframe1, aes(x=Branches, y=Number_outliers)) + geom_bar(stat ="identity", fill="purple", colour="black") + geom_text(aes(label=Number_outliers), vjust=-0.2, position = position_dodge(.9), size=3) + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text.x = element_text(angle = 90, vjust = 0.5)) + ggtitle("Number of SNPs under positive selection by branch-Phylogeny 1") + theme(plot.title = element_text(hjust = 0.5))



#### Phylogeny 2 ####

# read the table 
tabletsv2 <- read.csv("human_chr12_21KG_CLM_PUR_FRN_SAR_ADY_RUS_ESN_MZB_HZP_KLS_MAN_MSL_IBS_BAS_PEL_PTH_CEU_MXL.tsv", sep = '\t', header = TRUE)

# check the table dimensions 
head(tabletsv2)
count(tabletsv2)
nrow(tabletsv2)
ncol(tabletsv2)

# find out the number of columns of interest 
which(colnames(tabletsv2)=="Pval_u_r" )
which(colnames(tabletsv2)=="Pval_KLS_h")

# create an empty vector to include the outlier counts by branch
outliers_branch_2 <- c()
outliers_branch_2
typeof(outliers_branch_2)

# iterate over each branch column 
range <- 38:ncol(tabletsv1)
range
length(range)

for (i in range)
{
  # append counts of outliers by branch to the empty vector
  outliers_branch_2[i-37] <- sum(tabletsv2[,i] < 0.00001)
}

# check vector content and length
outliers_branch_2
length(outliers_branch_2)

# find out the number of columns of interest 
which(colnames(tabletsv2)=="u_r" )
which(colnames(tabletsv2)=="KLS_h")

# create a vector with the names of branches 
names_branches_2 <- colnames(tabletsv2)[4:37]
names_branches_2
as.vector(names_branches_2)
class(names_branches_2)

#### create dataframe 1 with data of phylogeny 1 ####
dataframe2 <- data.frame(names_branches_2, outliers_branch_2)
dataframe2

# change the name of columns
colnames(dataframe2)[1] <- "Branches"
colnames(dataframe2)[2] <- "Number_outliers"
dataframe2

# plot data
ggplot(dataframe2, aes(x=Branches, y=Number_outliers)) + geom_bar(stat ="identity", fill="green", colour="black") + geom_text(aes(label=Number_outliers), vjust=-0.2, position = position_dodge(.9), size=3) + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text.x = element_text(angle = 90, vjust = 0.5)) + ggtitle("Number of SNPs under positive selection by branch-Phylogeny 2") + theme(plot.title = element_text(hjust = 0.5))
