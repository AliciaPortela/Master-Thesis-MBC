#!/bin/bash

./plink2/plink2 --allow-extra-chr --freq --loop-cats POPULATION --pheno popfile2.txt --set-missing-var-ids @:# --vcf hgdp_tgp_sgdp_chr12_p_invsites_LD_pruned.vcf --out allele-freq_chr12 --threads 15
