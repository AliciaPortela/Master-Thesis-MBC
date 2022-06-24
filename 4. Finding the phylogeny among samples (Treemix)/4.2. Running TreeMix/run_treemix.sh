#!/bin/bash

# from an allele-counts file, this command outputs an unrooted tree
~/treemix-1.13/src/treemix -i input_treemix_chr12_0.05.prn.gz -bootstrap -k 500 -o chr12_0.05
