#!/bin/bash

# hgdp_tgp_sgdp_chr12_p.dated.trees downloaded from Zenodo
# convert it to VCF format 
python3 -m tskit vcf hgdp_tgp_sgdp_chr12_p.dated.trees > hgdp_tgp_sgdp_chr12_p.dated.vcf
