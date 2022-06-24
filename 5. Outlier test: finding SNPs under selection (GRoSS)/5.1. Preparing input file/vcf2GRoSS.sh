#!/bin/bash

out_file="Master-Thesis-MBC/5. Outlier test: finding SNPs under selection (GRoSS)/5.2. Running GRoSS/Inputs/"

./glactools/glactools vcfm2acf --onlyGT --fai human_g1k_v37.fasta.fai hgdp_tgp_sgdp_chr12_p.dated_287.vcf | ./glactools/glactools meld -f panel.txt - | ./glactools/glactools acf2gross --noroot - | gzip > $out_file/"human_chr12.gross.gz"
