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
Extract forward reads with the sed command to replace the 10N separator into a tab character and then use the cut function to split the files into forward reads.
```
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
```
Cluster all forward reads by 80% similarity. 
```
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1
```
Convert output of CD-hit to match output of first phase of rainbow.
```
mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster
```
The output should follow a simple text format of:
```
Read_ID	Cluster_ID	Forward_Read	Reverse_Read
```
Get the exact number of clusters
```
cut -f2 rcluster | uniq | wc -l 
```
In this case, the number of clusters is 1000.

**Split clusters into smaller clusters that represent significant variants**
The first clustering steps found RAD loci, now split the loci into alleles.
```
rainbow div -i rcluster -o rbdiv.out
```
The output of the div process is similar to the previous output with the exception that the second column is now the new divided cluster_ID (this value is numbered sequentially) and there was a column added to the end of the file that holds the original first cluster ID The parameter -f can be set to control what is the minimum frequency of an allele necessary to divide it into its own cluster Since this is from multiple individuals, we want to lower this from the default of 0.2.
```
rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
```
A parameter of interest to add here is the -r parameter, which is the minimum number of reads to assemble. The default is 5 which works well if assembling reads from a single individual. However, we are assembling a reduced data set, so there may only be one copy of a locus. Therefore, it's more appropriate to use a cutoff of 2.
```
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
```
The rbasm output lists optimal and suboptimal contigs. Previous versions of dDocent used rainbow's included perl scripts to retrieve optimal contigs. However, as of version 2.0, dDocent uses customized AWK code to extract optimal contigs for RAD sequencing.
```
cat rbasm.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
if (NR == 1) e=$2;
else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
else if ($1 ~/C/) clus=$2;
else if ($1 ~/L/) len=$2;
else if ($1 ~/S/) seq=$2;
else if ($1 ~/N/) freq=$2;
else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
}' > rainbow.fasta
```
First, the script looks at all the contigs assembled for a cluster. If any of the contigs contain forward and PE reads, then that contig is output as optimal. If no overlap contigs exists (the usual for most RAD data sets), then the contig with the most assembled reads PE (most common) is output with the forward read contig with a 10 N spacer. If two contigs have equal number of reads, the longer contig is output.

Align and cluster data by sequence similarity with cd-hit.
```
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c 0.9
```
The -M and -T flags instruct the program on memory usage (-M) and number of threads (-T). Setting the value to 0 uses all available. The real parameter of significan is the -c parameter which sets the percentage of sequence similarity to group contigs by. The above code uses 90%. Try using 95%, 85%, 80%, and 99%. Since this is simulated data, we know the real number of contigs, 1000. By choosing an cutoffs of 4 and 4, we are able to get the real number of contigs, no matter what the similarty cutoff.

In this example, it's easy to know the correct number of reference contigs, but with real data this is less obvious. As you just demonstrated, varying the uniq sequence copy cutoff and the final clustering similarity have the the largest effect on the number of final contigs. You could go back and retype all the steps from above to explore the data, but scripting makes this easier. Place remake_reference.sh (in this repository) in ddocentdir, or whatever your working directory is.

Remake reference by calling the script along with a new cutoff value and similarity.
```
bash remake_reference.sh 4 4 0.90 PE 2
```
This command will remake the reference with a cutoff of 20 copies of a unique sequence to use for assembly and a final clustering value of 90%. It will output the number of reference sequences and create a new, indexed reference with the given parameters. The output from the code above should be "1000" Experiment with some different values on your own.
What you choose for a final number of contigs will be something of a judgement call. However, we could try to heuristically search the parameter space to find an optimal value.

Use this repository's ReferenceOpt.sh to automate this process.
This script uses different loops to assemble references from an interval of cutoff values and c values from 0.8-0.98. 

Run ReferenceOpt.sh. It take as a while to run (5-20 mins).
```
bash ReferenceOpt.sh 4 8 4 8 PE 16
```
![RefContigs](https://github.com/MattLondon101/Reference-Genome-Assembly/blob/master/Images/refContigs_simrad1.png)

You can see that the most common number of contigs across all iteration is 1000, but also that the top three occuring and the average are all within 1% of the true value Again, this is simulated data and with real data, the number of exact reference contigs is unknown and you will ultimately have to make a judgement call.

Examine the reference
```
bash remake_reference.sh 4 4 0.90 PE 2
head reference.fasta
```
```
>dDocent_Contig_1
NAATTCCTCCGACATTGTCGGCTTTAAATAGCTCATAACTTGAGCCCAGGTAAAGACTTTAGTATACTCGCACCTTCCGCTTATCCCCCGGCCGCNNNNNNNNNNATTCAACCGCGGGACCTGAACTAACATAGCGTTGTGTATACCATCCGAGGTAACCTTATAACTCTCTGCCATTCGGACAGGTAACACGGCATATCGTCCGN
>dDocent_Contig_2
NAATTCAGAATGGTCATACAGGGCGGTAGAATGGAATCCTGAATCGAATGGCGGTTGCATTGAGAACCTGGTACCAGATAGGATCTGGATTAAATNNNNNNNNNNGTCGGGTACTAATTATCTATTGGGTCCAAACCCTCCGCCCCGTTTACTGCCCACCCGGCATGCAGTCATGAGAATTCCAAGGAACTAAGATAAGAGACCGN
>dDocent_Contig_3
NAATTCGGGCTCCTTGGAGAGATTCTTTCAATTATGCCCCCTACGTGGGAAACAGGGTCGGAAGTGGTCGGCTGAGAATTACTCGAAAGCCGCTCNNNNNNNNNNCCACCAGCATGATAGGACTTCAAGCTTGCCGTTTGTTGGGAGGACCGGTCGCTACGGAGCTGACGCTATCTCCCGCATCGGACCTCGTGGACAAAAACCGN
>dDocent_Contig_4
NAATTCAAAAGTCGCCCATAGGTACGTGATGAATTAGGTCAAGCGGGGACGTCGCATAGATGCGTGACGTCTGGAGCATGATGTTGTTTCTAACCNNNNNNNNNNAATCACTCGGTCAACGTGGTCCGTGCTCTGCAACGAAAAAAACTTCGCATGTGAACGATGATGCCTATAGGTGCGACCGCCGTCAGAGGCCCGTTGACCGN
>dDocent_Contig_5
NAATTCATACGGATATGATACTTCGTCTGGCAGGGTGGCTAGCGAGTTTAAGGATTCTTGGATAAAGGTAGGTAAAATTCTCGAGATTCTGATCTNNNNNNNNNNTAGAGGTGCTGGCGGGGCCTAGACGTGTTTCTACGCTTACTGATCAAATTAGCTAGCTTAGGTTCCTATAGTCTACGCTGGATTGTCCTTAGATGCACCGN
```
You can now see that we have complete RAD fragments starting with our EcoRI restriction site (AATT), followed by R1, then a filler of 10Ns, and then R2 ending with the mspI restriction site (CCG). The start and end of the sequence are buffered with a single N

We can use simple shell commands to query this data. Find out how many lines in the file (this is double the number of sequences)
```
wc -l reference.fasta
```
Find out how many sequences there are directly by counting lines that only start with the header character ">"
```
mawk '/>/' reference.fasta | wc -l
```
We can test that all sequences follow the expected format.
```
mawk '/^NAATT.*N*.*CCGN$/' reference.fasta | wc -l
grep '^NAATT.*N*.*CCGN$' reference.fasta | wc -l
```
1000 it is!

**Assemble references**
Assemble references across cutoff values, then map 20 random samples and evaluate mappings to the reference, along with number of contigs and coverage.

Run this repository's RefMapOpt.sh. It take as a while to run (5-20 mins).
```
RefMapOpt.sh 4 8 4 8 0.9 64 PE
```
This will loop across cutoffs of 4-8 using a similarity of 90% for clustering, parellized across 64 processors, using PE assembly technique.

The output is stored in a file called mapping.results
```
cat mapping.results
```
The output contains the average coverage per contig, the average coverage per contig not counting zero coverage contigs, the number of contigs, the mean number of contigs mapped, the two cutoff values used, the sum of all mapped reads, the sum of all properly mapped reads, the mean number of mapped reads, the mean number of properly mapped reads, and the number of reads that are mapped to mismatching contigs. Here, we are looking to values that maximize properly mapped reads, the mean number of contigs mapped, and the coverage. In this example, it's easy. Values 4,7 produce the highes number of properly mapped reads, coverage, and contigs.
Real data will involve a judgement call.




