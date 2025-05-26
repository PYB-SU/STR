# CONTENTS
Code for finding and analysing tandem repeats in DNA long reads obtained by Naonopore ONT.

The code is organized as : 
* PIPELINE : suite of bash shells
* STRops : an R package for downstream computation

## PIPELINE
* stepA.sh : obtain GENE coordinates from T2T genome
```
code
```
* stepM.sh : map fastq to T2T
* stepO.job : extract reads overlaping GENE
* stepR.job : extract reads mapping before AND after repeat region
* stepC.job : run charONT PIPELINE  
* stepH.job : compute haplotype from reads
* stepS.job : extract repeated regions from reads
* stepM.job : compute methylation

## STRops
An R package taking the ouptut files of the PIPELINE t draw waterfall and methylation
