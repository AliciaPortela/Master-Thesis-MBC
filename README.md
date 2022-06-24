# Master-Thesis-MBC
Final Master's Thesis of the Master in Computational Biology (UPM)

# Procedures 
1. Environmental data retrieval 

The environmental data retrieval of populations belonging to 1000 Genomes and Simons Genome Diversity Project (SGDP) genomics databases was carried out prior to the realization of this Master’s Thesis and remains outside the limits that it covers. However, it is necessary to give a brief explanation of how this process was accomplished, which will help to understand the rest of the procedures. 
Environmental data for three environmental variables were downloaded: 1) radiation (rad) from Worldclim; 2) temperature and humidity (th) from ecoClimate and 3) elevation (ele) from ETOPO1 Global Relief Model. To reduce the number of correlated variables in rad and th, a Principal Component Analysis (PCA) was carried out. In the case of ele, original variables were used. With data of each of the three environmental variables, three environmental raster images were built showing environmental gradients across the whole world. A raster image is a two-dimensional grid of square pixels. In this case, each pixel was assigned an environmental value. The interest in using raster for data representation lies in its greater efficiency to manage environmental data.

2. Selection of human populations and genomic samples 

Genomic human samples were selected according to the population of origin to which they belonged. The first step was to search for those populations, based on four criteria: 
	
  - Environmental divergence among populations
  
This criterion aims to maximize the environmental variance and minimize the geographic distance between the populations (Rellstab et al., 2015). With this, it is possible to correct the effect of isolation by distance (IBD) between them, so that if genome variations are found among them, they will be due to adaptive processes with greater probability.  To achieve that, the environmental divergence index was calculated for every pair of populations and every environmental variable as follows:
Environmental divergence index = (Environmental distance between populations)/(Geographic distance between populations)
Then, a ranking with the maximum punctuations of the environmental divergence index of pairs of populations was established, one for each environmental variable. 

  - Population size
  
From those rankings, only pairs of populations with a minimum number of genotyped individuals were selected. This criterion is based on the results thrown by Graph-aware Retrieval of Selective Sweeps (GRoSS) (Refoyo-Martínez et al., 2019), the tool that is used in this study to find loci under selection. For low sample sizes (N = 4), the authors found high rates of type I and II errors, whereas it does not happen for greater sample sizes. In order to ensure the minimization of errors, the threshold was established at N = 12. For that reason, any pair of populations that did not meet this criterion was excluded.

   - Availability of their genomic data in genomic databases of interest
  
At first, the genomic data was intended to be extracted from different independent databases, to later unify them. The genomic databases from which it is easier to extract genomic data in VCF format are 1000 Genomes and Human Genome Diversity Project (HGDP). HGDP is a genomic database similar to SGDP, with the addition that it contains a greater number of samples by population. For that reason, only those pairs of populations whose individual’s genomic data were included in the 1000 Genomes and HGDP were selected (rankings3). However, a recent article (Wohns et al., 2022) was published while this study was being developed. In it, genomic data from 1000 Genomes, HGDP and SGDP are unified in the same database. This fact facilitated the process of genomic data extraction (see the section ‘Genomic data retrieval’). Due to the change in the methodological strategy, the application of this criterion was no longer necessary. Even so, it was maintained, since the selection of the populations was carried out prior to the publication of the aforementioned article.

   - One pair of populations by great continent 
  
Continents are the geographic accidents that contribute the most to the genetic structure. In order to have a representative sample of humankind, only four pairs of populations (each belonging to a great continent: Africa, Europe, Asia and America) were selected for each environmental variable. The exception was rad, for which five pairs were selected. It was considered interesting to include two populations from different continents. 

After applying all these four criteria, 18 populations were included in the study. From each population, 20 individuals were randomly selected, provided that the population size allowed it. A total of 330 samples were selected, although not all of them were finally included in the analysis (see the section ‘Genomic data retrieval’).

