# CONTENTS
Code for finding and analysing tandem repeats in DNA long reads obtained by Naonopore ONT.

The code is organized as : 
* PIPELINE : suite of bash shells
* STRops : an R package for downstream computation

## PREQUISITES
A gff3 T2T genome file named  chm13.draft_v2.0.gene_annotation.gff3. Its location can be changed by changing REFS.
A fastq file for the T2T genome in REFS directory.

The Shell commands use the following tools : 


## PIPELINE
* stepA.sh : obtain GENE coordinates from T2T genome
```
stepA.sh GENE repeat_left repeat_right
```
GENE chromosome and position is looked for in $GENOME.gff3 file. 
Then flanking sequences ending in repeat_left and starting in repeat_right with various sizes are extracted from $GENOME.fa. 
Sequence files are stored in $PIPELINEDIR/GENE

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

```
codeGeneAnalysis.R
codePlotMethylation.R
```
contains the code in the R folder of the package.

