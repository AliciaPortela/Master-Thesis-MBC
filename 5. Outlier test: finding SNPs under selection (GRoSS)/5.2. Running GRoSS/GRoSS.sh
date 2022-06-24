#!/bin/bash

output="Master-Thesis-MBC/5. Outlier test: finding SNPs under selection (GRoSS)/5.2. Running GRoSS/Outputs"


# for Treemix topology
Rscript GRoSS.R -e human_chr12.gross -r human_chr12_1.graph -o $output/"human_chr12_1.tsv"

# for our topology
Rscript GRoSS.R -e human_chr12.gross -r human_chr12_2.graph -o $output/"human_chr12_2.tsv"
