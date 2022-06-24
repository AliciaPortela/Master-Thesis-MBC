#!/bin/bash

output="Master-s-Thesis-/Outlier_test/Run_GRoSS"

./glactools/glactools vcfm2acf --onlyGT --fai human_g1k_v37.fasta.fai hgdp_tgp_sgdp_chr12_p.dated_287.vcf | ./glactools/glactools meld -f panel.txt - | ./glactools/glactools acf2gross --noroot - | gzip > $output/"human_chr12.gross.gz"
