#!/bin/bash

# filtering the VCF file by 287 individuals (merged individuals)
vcftools --vcf hgdp_tgp_sgdp_chr12_p.dated.vcf --keep samples.txt --recode --stdout | bgzip -c > hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz

# hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz is the final VCF file 
