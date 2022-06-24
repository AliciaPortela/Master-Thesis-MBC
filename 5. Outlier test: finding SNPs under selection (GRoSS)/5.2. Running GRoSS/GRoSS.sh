#!/bin/bash

output="Master-s-Thesis-/Outlier_test/Output_GRoSS/"


# for Treemix topology
Rscript GRoSS.R -e human_chr12.gross -r human_chr12_1.graph -o $output/"human_chr12_1.tsv"

# for our topology
Rscript GRoSS.R -e human_chr12.gross -r human_chr12_2.graph -o $output/"human_chr12_2.tsv"
