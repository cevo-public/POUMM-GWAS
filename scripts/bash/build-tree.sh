#!/bin/bash
# This script generates an alignment, builds a tree.
# Run with: bash /scripts/bash/build-tree.sh

# Generate alignment from fasta file
./muscle3.8.31_i86linux64 \
-in sequences.fasta \
-out output/pathogen.fasta \
-fasta \
-maxiters 3

# Build a ML tree with the GTR+F+R4 model and (TODO) 100 bootstrap replicates
iqtree \
-s output/pathogen.fasta \
-m GTR+F+R4 \
-redo \
-pre output/pathogen \
-nt AUTO

# Hodcroft stripped sequences of codons in positions associated with drug resistance mutations (but no sig. differences to when these left in apparently)
# 38 reference pol sequences subtypes A-K from Los Alamos HIV database used as outgroup
# Used RAxML to build the tree, RAxML default model is hard to tell from quick look at documentation
# Used 100 bootstrap replicates



