# CONTENTS
Code for finding and analysing tandem repeats in DNA long reads obtained by Naonopore ONT.

The code is organized as : 
* PIPELINE : suite of bash shells
* STRops : an R package for downstream computation

## PREQUISITES
The Shell commands use the following tools : 
* bedtools
* bgzip
* 
* minimap2
* samtools
* seqtk
* zcat

## PIPELINE
* stepA.sh : obtain GENE coordinates from T2T genome 
```
stepA.sh GENE repeat_left repeat_right
```
Modify in file : 
* $REFS directory : contains reference files, including  
* $GENOME_GFF : a gff3 T2T genome file named.
* $GENOME_FA : a fastq file for the T2T genome.

Looks up GENE for chromosome and position in $GENOME.gff3 file. 
Extracts flanking sequences flanking the repeat: to the left of repeat_left and to the right of repeat_right. 
Sequences of various sizes are extracted from $GENOME.fa. 
Sequence files are stored in $PIPELINEDIR/GENE

* stepM.sh : submit batch jobs to map initial fastq to T2T
```
stepM.sh HERE GO DEP
```
Write .job files to perform mapping for each file found in $HERE/BAM-ORIG (unaligned BAM file) or $HERE/FASTQ (fastq file). 
If GO is not null, submit .job to slurm. 
If DEP is not null, add DEP as a dependency in sbatch

* stepO.job : extract reads overlaping GENE
```
stepO.job HERE GENE
```
Extract reads overlaping GENE position from bam files in $HERE/BAM

* stepR.job : extract reads mapping before AND after repeat region
```
stepR.job HERE GENE
```
Align each read in BAM file with repeat flanking sequences taken at increasing distances 
Select result providing the largest number of reads
* stepC.job : run charONT PIPELINE
```
stepC HERE GENE MODEL
```
MODEL refers to chemistry fro medaka
* stepH.job : compute haplotype from reads
```
stepH.job HERE GENE
```
compute haplotypes from reads in $HERE/OVERLAPS/overlap_GENE.bam 
* stepS.job : extract repeated regions from reads
```
stepS.job HERE GENE 
```
extract repeated sequences from reads, computing summary measures.
* stepM.job : compute methylation
```
stepM.job HERE GENE MOD
```
extract methylation information from BAM files in HERE/METHYL for GENE and modification MOD
## STRops
A R package processing the ouptut files of the PIPELINE to draw waterfall and methylation plots

```
codeGeneAnalysis.R
codePlotMethylation.R
```
contains the code in the R folder of the package.

