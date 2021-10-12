#!/bin/bash
# This script generates an alignment, builds a tree.

# Data file is newfasta2019-10-25.fas which contains 2012 pol sequences corresponding to 2012/2125 genotyped samples (the remaining host genotypes are not associated with a pol sequence)
# The sequences are not aligned (some are longer than others) even though gap characters have been introduced --> what does this mean??

# Generate alignment from fasta file
### PUT IN RIGHT ALIGNMENT!
./muscle3.8.31_i86linux64 \
-in sequences.fasta \
-out output/pathogen.fasta \
-fasta \
-maxiters 3

# Build a ML tree with the GTR+F+R4 model and 100 bootstrap replicates
iqtree \
-s output/pathogen.fasta \
-m GTR+F+R4 \
#-b 100 \
-pre output/pathogen \
-nt AUTO

# Hodcroft stripped sequences of codons in positions associated with drug resistance mutations (but no sig. differences to when these left in apparently)
# 38 reference pol sequences subtypes A-K from Los Alamos HIV database used as outgroup
# Used RAxML to build the tree, RAxML default model is hard to tell from quick look at documentation
# Used 100 bootstrap replicates



