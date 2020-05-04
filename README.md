# Reference-Genome-Assembly
Guide for creation of reference genomes in Linux with dDocent and the Stacks package.

A reference genome (aka reference assembly) serves as a representative sample of a set of genes for an idealized organism of a species and are typically used as guides for developing new genomes. Their creation often comprises assembly of reference contigs, which are sets of overlapping DNA segments that form a consensus region of DNA.

This guide describes the process of using genetic sequence data to plot reference contig data, generate a fasta file with complete restriction enzyme associated DNA (RAD) fragments, and output the optimal interval of cutoff values to capture maximum diversity and accuracy.

**Requirements**
* Linux OS (I used Windows Subsystem for Linux in Windows 10)
* dDocent (a bash wrapper to QC, assemble, map, and call SNPs from almost any kind of RAD sequencing)
  * with conda installed in your Linux terminal:
```
# add the bioconda channel
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# install dDocent and activate dDocent environment
conda install ddocent
conda create -n ddocent_env ddocent
source activate ddocent_env
```

