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

Download Stacks from http://catchenlab.life.illinois.edu/stacks/ to ddocentdir
Untar and install Stacks
```
tar xfvz stacks-2.xx.tar.gz
cd stacks-2.xx.tar.gz
./configure --prefix=/home/joe/ddocentdir
sudo make install
```
Run process_radtags
```
process_radtags -1 SimRAD_R1.fastq.gz -2 SimRAD_R2.fastq.gz -b barcodes -e ecoRI --renz_2 mspI -r -i gzfastq
```
The option -e specifies the 5' restriction site and --renze_2 specifes the 3' restriction site. -i states the format of the input sequences. The -r option tells the program to fix cut sites and barcodes that have up to 1-2 mutations in them.

When the program is completed, ddocentdir should have several files that look like: sample_AAGAGG.1.fq.gz, sample_AAGAGG.2.fq.gz, sample_AAGAGG.rem.1.fq.gz, and sample_AAGAGG.rem.2.fq.gz

The .rem..fq.gz files would normally have files that fail process_radtags (bad barcode, ambitious cut sites), but we have simulated data and none of those bad reads. We can delete.
```
rm *rem*
```
The individual files are currently only names by barcode sequence. They can be renamed in an easier convention with the following bash script. It reads the names into an array and all barcodes into second array, gets length of both arrays, then iterates the task of renaming the samples.
```
curl -L -O https://github.com/jpuritz/dDocent/raw/master/Rename_for_dDocent.sh
bash Rename_for_dDocent.sh SimRAD.barcodes
ls *.fq.gz
```
There should now be 40 individually labeled .F.fq.gz and 40 .R.fq.gz. Twenty from PopA and Twenty from PopB. 

**Assemble RAD data**
Create a set of unique
