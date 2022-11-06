# RNA-seq_pipeline_MSc_assignment1
An RNA-seq pipeline written with bash scripts for an MSc assignment

## Overview:

### 1. FASTQC quality check:
All paired-end RNAseq data in the fastq directory are processed by `fastqc` producing result zip files and extracted `fastqc` result directories for all sample files as an output. The resulting summary.txt files in these directories are processed to display the total numbers of FAILS and WARNINGS to the user and output a summary table file for all samples. The user then has the option to view the filenames and `fastqc` modules of all FAILS and WARNINGS. Based on these results the user then has the option to continue to the alignment and quantification steps of the pipeline or exit the pipeline if the quality of the raw sequence data requires further evaluation.

### 2. Alignment and quantification:
Taking the T. congolense genome sequence provided in fasta format, `bowtie2-build` produce an efficient data structure or ‘genome index’ for searching. `bowtie2` then takes each of the paired-end fastq files and aligns them against this index to produce alignments in .sam format – these are piped directly into `samtools view` which converts to compressed binary .bam format to reduce disk space. These bam files are then used as inputs for `samtools sort` to produce .bam files sorted by the left most alignment, which are subsequently used as inputs for `samtools index` to produce an index file for fast random access of the alignments. `bedtools coverage` then takes each of the sorted .bam files as inputs in addition to the T. congolense .bed file gene region information and counts the number of aligned reads that align to different gene regions, producing tab-delimited counts data files.

### 3. Mean counts data:
The 100k.fqfiles file with sample information is then sorted by all other experimental variables before sorting by replicate number to ensure that replicates are listed consecutively in triplets. The sample names of this sorted information are then extracted to
combine read counts for replicates into mean counts for each group of samples producing tab-delimited mean counts data files.

### 4. Fold change data:
All mean counts files are then compared to all other files individually to generate fold change data for all possible group-wise comparisons. This is done by determining the ratio of the mean counts for each gene between the files. If a division by zero is encountered, ‘NA’ is printed – it is not possible to determine the fold change between these values. Fold changes of gene counts are sorted such that they are listed in decreasing order. The output is stored in tab-delimited fold change data files.
