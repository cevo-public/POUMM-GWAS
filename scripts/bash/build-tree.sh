#!/bin/bash
# This script generates an alignment, builds a tree.
# Run with: bash /scripts/bash/build-tree.sh

# Generate alignment from fasta file
./muscle3.8.31_i86linux64 \
-in sequences.fasta \
-out output/pathogen.fasta \
-fasta \
-maxiters 3

# Trim the alignment, first removing wrapping
sed -e 's/\(^>.*$\)/#\1#/' output/pathogen.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > output/pathogen_oneline.fasta
cut -c 1-1505 output/pathogen_oneline.fasta > output/pathogen_trimmed.fasta

# Build a ML tree with the GTR+F+R4 model and (TODO) 100 bootstrap replicates
iqtree \
-s output/pathogen_trimmed.fasta \
-m GTR+F+R4 \
-redo \
-pre output/pathogen \
-nt AUTO

# Hodcroft stripped sequences of codons in positions associated with drug resistance mutations (but no sig. differences to when these left in apparently)
# 38 reference pol sequences subtypes A-K from Los Alamos HIV database used as outgroup
# Used RAxML to build the tree, RAxML default model is hard to tell from quick look at documentation
# Used 100 bootstrap replicates



