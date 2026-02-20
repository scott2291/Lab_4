## Login in OSC, either by using SSH in your terminal, a command line instance in your web browser, or using VSCode from OSC
## Check where you are and make sure you create your working directory in the right place

```shell
cd /fs/scratch/PAS3260/User/
mkdir Lab_4
```

## Copy files to your working directory, for example:

```shell
cp /fs/scratch/PAS3260/Jonathan/* /fs/scratch/PAS3260/mia/Lab_4
```

## Go to your working directory and execute:
```shell
cd Lab_4
ls -l
ls -Llh
ls -lh
```
### What do you see?

For `ls -l` the output is:

```shell
drwxr-xr-x 2 scott2291 PAS2880 4096 Feb 19 19:18 data
-rw-r--r-- 1 scott2291 PAS2880 3319 Feb 19 19:58 Lab_4.md
-rwxr-xr-x 1 scott2291 PAS2880  351 Feb  4 16:47 Lab5.md5
drwxr-xr-x 2 scott2291 PAS2880 4096 Feb 19 19:18 results
```
For `ls -Llh` the output is:

```shell
total 2.0K
drwxr-xr-x 2 scott2291 PAS2880 4.0K Feb 19 19:18 data
-rw-r--r-- 1 scott2291 PAS2880 3.6K Feb 19 19:59 Lab_4.md
-rwxr-xr-x 1 scott2291 PAS2880  351 Feb  4 16:47 Lab5.md5
drwxr-xr-x 2 scott2291 PAS2880 4.0K Feb 19 19:18 results
```
For `ls -lh` the output is:

```shell
total 2.0K
drwxr-xr-x 2 scott2291 PAS2880 4.0K Feb 19 19:18 data
-rw-r--r-- 1 scott2291 PAS2880 3.8K Feb 19 20:00 Lab_4.md
-rwxr-xr-x 1 scott2291 PAS2880  351 Feb  4 16:47 Lab5.md5
drwxr-xr-x 2 scott2291 PAS2880 4.0K Feb 19 19:18 results
```

## There is a md5 file, what are md5 files for? Explain

md5sum are used to check the integrity of the data files in the 

## Reorganizing files

Before I work with the md5 file, I want to reorganize my data by moving it into its own directory alongside the md5 file. I also want the text file I pipe the logging info into to be in a `results` directory.

```shell
mkdir data results
mv ERR3638927_1.fastq.gz data
mv ref_Run1CB3334.fasta data
mv Rice.gff3 data
mv Run* data
mv Lab5.md5 data

ls -lh data
```
The output of these commands are:

```shell
total 2.7G
-rw-r--r-- 1 scott2291 PAS2880  49M Feb 18 16:51 ERR3638927_1.fastq.gz
-rwxr-xr-x 1 scott2291 PAS2880  351 Feb  4 16:47 Lab5.md5
-rw-r--r-- 1 scott2291 PAS2880  279 Feb 18 16:21 ref_Run1CB3334.fasta
-rw-r--r-- 1 scott2291 PAS2880  54M Feb  4 16:26 Rice.gff3
-rw-r--r-- 1 scott2291 PAS2880 156M Feb  4 16:26 Run1CB3334.bam
-rw-r--r-- 1 scott2291 PAS2880 1.1G Feb  4 16:26 Run1CB3334.fastq
-rw-r--r-- 1 scott2291 PAS2880 134M Feb  4 16:26 Run1CB3334.fastq.gz
-rw-r--r-- 1 scott2291 PAS2880 1.2G Feb  4 16:26 Run1CB3334.sam
-rw-r--r-- 1 scott2291 PAS2880 420K Feb  4 16:26 Run1CB3334.vcf
```

## Let's work with this file:
```shell
cd data
md5sum -c Lab5.md5 > results/Checking.txt
grep "FAIL" Checking.txt
## Parenthesis on integrity checking
touch myDNA.fasta
nano myDNA.fasta
md5sum myDNA.fasta
# Inspect the files
```
The output of this command is:

```shell
68b329da9893e34099c7d8ad5cb9c940  myDNA.fasta
```
The md5sum command generated a unique cryptographic key to identify this file. I want to see what would happen if I change this file by adding text into it

```shell
nano myDNA.fasta
md5sum myDNA.fasta
```

The output to this command is:

```shell
e59ff97941044f85df5297e1c302d260  myDNA.fasta
```

Now there is a new cryptographic key for this file which shows that the file has been altered in some way.

## Now, let's explore the files a little:
```shell
head ref_Run1CB3334.fasta
head Rice.gff3
head Run1CB3334.fastq
head Run1CB3334.vcf
head Run1CB3334.fastq.gz
head Run1CB3334.bam
```
### Do you see something weird?

The gzipped fastq file is machine readable instead of human readable, so the output is not legible.

## For `Run1CB3334.fastq.gz` try:
```shell
zcat Run1CB3334.fastq.gz | head
```
This command made the output readable, and similar to the fastq file format I have previously seen.

## Let's check the fastq file in a little more detail
```shell
head -100 Run1CB3334.fastq | less
```
I see a lot of capital letters in the quality score line of each read, which indicates the reads are not of high quality.

### How would you count the number of reads in the fastq file?

```shell
wc -l Run1CB3334.fastq
```

The output to this command was `7560972 Run1CB3334.fastq` which indicated to me that there are approximately 7.5 million reads within this fastq file and ~30 million lines within the file.

# Is it? Try:
```shell
zcat Run1CB3334.fastq.gz | echo "$((`wc -l` / 4))"
```
The output to this command was `1890243` which indicates to me that there are ~8 million lines within this file.

Why the difference?

>I think there is a difference between the number of lines within these files due to the nature of compressing files. When looking at the sizes of these files, the compressed file is signifincantly smaller. Therefore, the number of lines in the machine readable, compressed version of the fastq file may not hold much meaning compared to the unzipped human-readable version.


### Let's download some data from NCBI
```
module load sratoolkit/3.0.2
```
# What is sratoolkit?

It is a module that can be loaded into the Pitzer cluster of OSC that is capable of storing raw sequence data and moving it to your files.
https://www.osc.edu/resources/available_software/software_list/sra_toolkit

```shell
fastq-dump --gzip --split-files --readids --origfmt ERR3638927
```
(a faster option is fasterq-dump, how can you look for information about this command?)

The website I linked above provides information for each command within the `sratoolkit`. Within this link, there is a link to the github for this tool which provides even more information. When I don't know how to use a tool and the name of the tool and `--help` does not work, I can google information about the tool.

https://github.com/ncbi/sra-tools/wiki

```shell
zcat ERR3638927_1.fastq.gz | head
```
## Let's download a SAM file
```shell
sam-dump --output-file ERR3638927.sam ERR3638927
```
To manipulate the SAM file we will need samtools:
```shell
module load samtools/1.21
```
Before working with `samtools` I will try to learn more about it using the `--help` option

```shell

```

Then, we can make a sorted BAM file:

```shell
samtools view -bS ERR3638927.sam | samtools sort -o ERR3638927_sorted.bam
```
Index it:
```shell
samtools index ERR3638927_sorted.bam
```
Make it a FASTA file:
```shell
samtools fasta ERR3638927_sorted.bam > ERR3638927.fasta
head ERR3638927.fasta
```
And index it:
```shell
samtools faidx ERR3638927.fasta
```
## Finally, let's check the VCF file
```shell
head -100 Run1CB3334.vcf | less
module load vcftools/0.1.16
vcftools --vcf Run1CB3334.vcf
vcftools --vcf Run1CB3334.vcf --maf 0.05
vcftools --vcf Run1CB3334.vcf --maf 0.05 --recode --out Run1CB3334_Filtered_maf05
vcftools --vcf Run1CB3334_Filtered_maf05.recode.vcf
```
## Can we check the VCF file using a container as in Lab 3?
```shell
module unload vcftools/0.1.16
vcftools --help
apptainer pull vcftools.sif docker://quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_9
apptainer exec vcftools.sif vcftools --help
apptainer exec vcftools.sif vcftools --vcf Run1CB3334.vcf --maf 0.05 --recode --out Run1CB3334_Filtered_maf05_from_container
apptainer exec vcftools.sif vcftools --vcf Run1CB3334_Filtered_maf05_from_container.recode.vcf
```
## How could we use a container of samtools?
```shell
module unload samtools/1.10
samtools --help
apptainer pull samtools.sif docker://
apptainer exec samtools.sif samtools --help
```
Followign the same structure as VCFtools, define the commands to sort and index fasta and BAM files using the samtools container