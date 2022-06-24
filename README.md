# Master-Thesis-MBC
Final Master's Thesis of the Master in Computational Biology (UPM)

# Procedures 
1. Environmental data retrieval 

The environmental data retrieval of populations belonging to 1000 Genomes and Simons Genome Diversity Project (SGDP) genomics databases was carried out prior to the realization of this Master’s Thesis and remains outside the limits that it covers. However, it is necessary to give a brief explanation of how this process was accomplished, which will help to understand the rest of the procedures. 
Environmental data for three environmental variables were downloaded: 1) radiation (rad) from Worldclim; 2) temperature and humidity (th) from ecoClimate and 3) elevation (ele) from ETOPO1 Global Relief Model. To reduce the number of correlated variables in rad and th, a Principal Component Analysis (PCA) was carried out. In the case of ele, original variables were used. With data of each of the three environmental variables (envdata1000Gen, envdataSGDP), three environmental raster images were built showing environmental gradients across the whole world (Fig.2). A raster image is a two-dimensional grid of square pixels. In this case, each pixel was assigned an environmental value. The interest in using raster for data representation lies in its greater efficiency to manage environmental data.
