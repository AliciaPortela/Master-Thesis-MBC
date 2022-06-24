#!/bin/bash

#filter invariant sites 
bcftools view --min-ac=1 --no-update hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz > hgdp_tgp_sgdp_chr12_p_invsites.vcf

#convert vcf into ped and map files 
vcftools --vcf hgdp_tgp_sgdp_chr12_p_invsites.vcf --plink --out hgdp_tgp_sgdp_chr12_p_invsites

#filter LD sites 
plink1.9 --file hgdp_tgp_sgdp_chr12_p_invsites --indep-pairwise 50 5 0.5 --recode --out hgdp_tgp_sgdp_chr12_p_invsites_LD

#new prune file
plink1.9 --file hgdp_tgp_sgdp_chr12_p_invsites_LD --extract hgdp_tgp_sgdp_chr12_p_invsites_LD.prune.in --make-bed --out hgdp_tgp_sgdp_chr12_p_invsites_LD_pruned

#convert to vcf
plink1.9 --bfile hgdp_tgp_sgdp_chr12_p_invsites_LD_pruned --recode vcf --out hgdp_tgp_sgdp_chr12_p_invsites_LD_pruned




