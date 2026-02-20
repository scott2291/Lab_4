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
cd /fs/scratch/PAS3260/mia/Lab_4
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
cd data
sam-dump --output-file ERR3638927.sam ERR3638927
```
To manipulate the SAM file we will need samtools:
```shell
cd ..
module load samtools/1.21
```
Before working with `samtools` I will try to learn more about it using the `--help` option

```shell

Program: samtools (Tools for alignments in the SAM format)
Version: 1.21 (using htslib 1.21)

Usage:   samtools <command> [options]

Commands:
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     fqidx          index/extract FASTQ
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate information
     reheader       replace BAM header
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags
     markdup        mark duplicates
     ampliconclip   clip oligos from the end of reads

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     consensus      produce a consensus Pileup/FASTA/FASTQ
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA
     import         Converts FASTA or FASTQ files to SAM/BAM/CRAM
     reference      Generates a reference from aligned data
     reset          Reverts aligner changes in reads

  -- Statistics
     bedcov         read depth per BED region
     coverage       alignment depth and percent coverage
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     cram-size      list CRAM Content-ID and Data-Series sizes
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)
     ampliconstats  generate amplicon specific stats

  -- Viewing
     flags          explain BAM flags
     head           header viewer
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM
     samples        list the samples in a set of SAM/BAM/CRAM files

  -- Misc
     help [cmd]     display this help message or help for [cmd]
     version        detailed version information
```

Then, we can make a sorted BAM file:

```shell
samtools view -bS data/ERR3638927.sam | samtools sort -o data/ERR3638927_sorted.bam
```
Command breakdown using help info:

`view`: used to convert between BAM, SAM, and CRAM files

`-bS`: According to the help command for `view`, the `b` indicates the ouput file will be a BAM file, and the `S` portion is ignored becasuse the input format is autodetected.

Then, the sam file is inputed and piped into the following commands

`sort`: used to sort alignment files. After looking further into the command, it sorts the alignment file by the leftmost coordinate when used in default settings. https://www.htslib.org/doc/samtools-sort.html

`-o`: writes to an output file rather than standard format

The file is then named what is specified after the `-o` option.

Index it:
```shell
samtools index data/ERR3638927_sorted.bam
```

Command breakdown:

`index`: :creates an index of the alignment for fast random access

Make it a FASTA file

```shell
samtools fasta data/ERR3638927_sorted.bam > data/ERR3638927.fasta
head data/ERR3638927.fasta
```
Command breakdown:

`fasta`: converts the BAM/SAM file to a fasta type file

And index it:
```shell
samtools faidx data/ERR3638927.fasta
```

Command Breakdown:

`faidx`: creates an index of a fasta type file

## Finally, let's check the VCF file

Before working with `vcftools` I want to learn more about the command using the `man vcftools`

```shell
head -100 data/Run1CB3334.vcf | less
module load vcftools/0.1.16
vcftools --vcf Run1CB3334.vcf
```
Command breakdown: Based on the output, it appears it completed some type of filtering on the vcf file, but the file we inputed did not have any data that was filtered out on the Individual or Site level. There was also several warnings that came up stating that there was not enough parts specified in the INFO or FORMAT entry for each variant

`--vcf`: processes the vcf file (look more into this later?)

```shell
vcftools --vcf data/Run1CB3334.vcf --maf 0.05
```

Command breakdown:  I got the same warnings in the output, but this time the filtering kept all the Individual but only 5 out of 40 Sites

`maf 0.05`: includes only sites with minor allele frequency greater than or equal to 0.05

```shell
vcftools --vcf data/Run1CB3334.vcf --maf 0.05 --recode --out results/Run1CB3334_Filtered_maf05
mv Run1CB3334_Filtered_maf05.log results/
```
Command breakdown: I got the exact same output to my terminal as the previous command, but there are now two additional files in my directory.

`--recode`: creates a new vcf file based on the input and filtering options specified

`--out`: defines the output files' basename

The output files include a filtered vcf file only containing sites with a minor allele frequenct greater than or equal to 0.05, and a log file that contains the information that was printed to the screen.

```shell
vcftools --vcf results/Run1CB3334_Filtered_maf05.recode.vcf
```
The output of this command resulted in all the Individuals and Sites remaining in the file after filtering the specfied file in default settings.

## Can we check the VCF file using a container as in Lab 3?

```shell
module unload vcftools/0.1.16
vcftools --help
apptainer pull vcftools.sif docker://quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_9
apptainer exec vcftools.sif vcftools --help
apptainer exec vcftools.sif vcftools --vcf data/Run1CB3334.vcf --maf 0.05 --recode --out results/Run1CB3334_Filtered_maf05_from_container
apptainer exec vcftools.sif vcftools --vcf results/Run1CB3334_Filtered_maf05_from_container.recode.vcf
mv Run1CB3334_Filtered_maf05_from_container.log results/
```
I got the exact same results as when I used the built in `vcftools` module.

## How could we use a container of samtools?
```shell
module unload samtools/1.10
samtools --help
apptainer pull samtools.sif docker://
apptainer exec samtools.sif samtools --help
```
Followign the same structure as VCFtools, define the commands to sort and index fasta and BAM files using the samtools container