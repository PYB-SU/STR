# CONTENTS
Code for finding and analysing tandem repeats in DNA long reads obtained by Naonopore ONT.

The code is organized as : 
* PIPELINE : suite of bash shells
* STRops : an R package for downstream computation

## PIPELINE
* stepA.sh : obtain GENE coordinates from T2T genome
```
stepA.sh GENE coord_left coord_right
```
Lookup GENE in T2T genome to find coordinates. Create sequence files centered around coord_left-coord_right
with various flanking sizes.
* stepM.sh : submit batch jobs to map fastq to T2T
```
stepM.sh HERE GO DEP
```
Write .job files to perform mapping for each fastsq file found in $HERE/FASTQ. 
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
