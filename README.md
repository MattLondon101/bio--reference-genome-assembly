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

Create a set of unique reads with counts for each individual
```
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
```
The first four lines set shell variables for various bits of AWK and perl code, to make parallelization with GNU-parallel easier. The first line after the variables, creates a set of forward reads for each individual by using mawk (a faster, c++ version of awk) to sort through the fastq file and strip away the quality scores. The second line does the same for the PE reads. Lastly, the final line concatentates the forward and PE reads together (with 10 Ns between them) and then find the unique reads within that individual and counts the occurences (coverage).

Sequences with very small levels of coverage within an individual are likely to be sequencing errors. So for assembly we can eliminate reads with low copy numbers to remove non-informative data. 

Sum up the numbers within individual coverage levels of unique reads in our data set. Use mawk to query the first column and select data above a certain copy number (from 2-20) and print to a file.
```
cat *.uniq.seqs > uniq.seqs
for i in {2..20};
do 
echo $i >> pfile
done
cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
```
Plot to terminal with gnuplot
```
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [2:20] 
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.data' with lines notitle
pause -1
EOF
```
![Coverage](https://github.com/MattLondon101/Reference-Genome-Assembly/blob/master/Images/coverage_simrad1.png)

Choose a cutoff value that captures as much of the diversity of the data as possible while simultaneously eliminating sequences that are likely errors. Try 4.
```
cat pfile | parallel --no-notice mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while (($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
wc -l uniqCperindv
```
The data has been reduced down to 7598 sequences. Yet, it can be reduced further. Restrict data by the number of different individuals a sequence appears within.
```
for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
```
Plot the data
```
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences present in more than X Individuals"
set xlabel "Number of Individuals"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.peri.data' with lines notitle
pause -1
EOF
```
![Individuals](https://github.com/MattLondon101/Reference-Genome-Assembly/blob/master/Images/individuals_simrad1.png)

Another cutoff value is chosen to capture maximum diversity in the data, while eliminating sequences with little value on the population scale. Try 4.
```
mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
wc -l uniq.k.4.c.4.seqs
```
The data has been reduced to 3840 sequences!

Convert sequences back to fasta format. This reads the totaluniqseq file line by line and adds a sequence header of >Contig X.
```
cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta
```
At this point, dDocent also checks for reads that have a substantial amount of Illumina adapter in them.
Our data is simulated and does not contain adapter, so we'll skip that step for the time being.

**Assemble Reference Contigs**  
Extract forward reads


