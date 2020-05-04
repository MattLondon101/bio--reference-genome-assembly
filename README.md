# Reference-Genome-Assembly
Guide for creation of reference genomes in Linux with dDocent and the Stacks package.

A reference genome (aka reference assembly) serves as a representative sample of a set of genes for an idealized organism of a species and are typically used as guides for developing new genomes. Their creation often comprises assembly of reference contigs, which are sets of overlapping DNA segments that form a consensus region of DNA.

This guide describes the process of using genetic sequence data to plot reference contig data, generate a fasta file with complete restriction enzyme associated DNA (RAD) fragments, and output the optimal interval of cutoff values to capture maximum diversity and accuracy.

**Requirements**
* Linux OS (I used Windows Subsystem for Linux in Windows 10)
* dDocent (a bash wrapper to QC, assemble, map, and call SNPs from almost any kind of RAD sequencing)
  * with conda installed in your Linux terminal:
```
# create a working directory
mkdir ddocentdir
cd ddocentdir
# add the bioconda channel
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# install dDocent and activate dDocent environment
conda install ddocent
conda create -n ddocent_env ddocent
source activate ddocent_env
```
**Test Dataset**
The SimRAD dataset contains DNA sequences from double digest restriction-site associated DNA (ddRAD). 
Download the dataset, unzip it, and view contents
```
curl -L -o data.zip https://www.dropbox.com/s/t09xjuudev4de72/data.zip?dl=0
unzip data.zip
ll
```
Create list of just barcodes
```
cut -f2 SimRAD.barcodes > barcodes
head barcodes
```
**Install Stacks package with process_radtags**

[Download Stacks from] (http://catchenlab.life.illinois.edu/stacks/)
