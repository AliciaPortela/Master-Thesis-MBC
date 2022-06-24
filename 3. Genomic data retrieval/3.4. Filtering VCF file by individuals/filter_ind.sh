#!/bin/bash

vcftools --vcf hgdp_tgp_sgdp_chr12_p.dated.vcf --keep samples.txt --recode --stdout | bgzip -c > hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz
# This is our final vcf
