# Identifying repeat counts 

## 1) pipelineSTR.sh
    - sh pipeline_STR.sh DIR GENE WHATTODO BONITO_MODEL MEDAKA_MODEL MODIF

DIR base directory containing a FAST5 subdir

GENE of interest. GENE information must be available in directory STR (using get_primer_and_flanking.sh)

WHATTODO is a character chain with:
- (B)onito for basecalling
- (P)orechop+zip+minimap2
- (O)verlap for extracting bam
- (R)eads extraction
- (C)harONT
- (S)equence from reads
- (M)ethylation from fast5

BONITO_MODEL and MEDAKA_MODEL is :
- HAC10 - uses dna_r10.4_e8.1_hac@v3.4 (bonito) and r1041_e82_400bps_hac (medaka)
- FAST10 - uses dna_r10.4_e8.1_fast@v3.4 (bonito) and r1041_e82_400bps_hac (medaka)
- SUP10 - uses dna_r10.4_e8.1_sup@v3.4 (bonito) and  r1041_e82_400bps_sup_g615 (medaka)
- HAC - uses dna_r9.4.1_e8.1_hac@v3.3 (bonito) and r941_e81_hac_g514 (medaka)
- FAST - uses dna_r9.4.1_e8.1_fast@v3.4 (bonito) and r941_e81_fast_g514 (medaka)
- SUP - uses dna_r9.4.1_e8.1_sup@v3.3 (bonito) and r941_e81_sup_g514 (medaka)

MODIF specifies the type of methylation :
-  5m : 5mc(bonito) , 5mCG (dorado)
-  5m+5hmc : 5mc_5hmc (bonito) , 5mCG_5hmCG (dorado)
-  6mA : 6mA (dorado)


example: 

    - sh pipelineSTR.sh 2022_GEN_DM1_pilote01 DMPK BORCSM HAC10 HAC10 5mCG_5hmCG

## 2) graphical representation

### geneAnalysis.R
takes files from the (S) step above and makes a summary of identified repeats. Makes a plot of repeats found in reads.

process.one.gene(DIR,root_dir, gene, query=NULL, ncommon=NULL, caller="guppy", QUALITY=20,IDENTITY=0.9, reverse=FALSE, complement=FALSE, first.codon=NULL,repeat.codon=NULL)
- DIR : base directory containing caller subdirectories
- CALLER : caller subdirectory (guppy, dorado_hac, bonito_hac)
- query: list of sequences to highlight, these will be matched in order and colored differently
- ncommon : frequent codons will be looked for in the sequences and add to query
- QUALITY : minimum average quality of reads to keep
- IDENTITY : % of identity in the left+right flanking sequences
- reverse : should sequences be reversed
- complement : should sequences be complemented
- first.codon : first codon to register sequences 
- repeat.codon : repeat codon that will be looked for to find stop in sequences

### plotMethylation.R
takes files from the (M) step and plots methylation % and heatmap along genomic information using Gviz

plot.methyl.fancy <- function(DIR,CALLER,gene,methyl, root="C:/Users/boelle/nextcloud_SU/paper/",
                              from=NULL, to=NULL, type="s", reads=c("short","long"),coverage=FALSE, repeat.rows=NULL,focus=NULL) {
- DIR : base directory containing caller subdirectories
- CALLER : caller subdirectory (guppy, dorado_hac, bonito_hac)
- gene : gene symbol to look for in genomic information
- root : root directory
- from, to: coordinates to limit the graph - superseded if (focus not null)
- methyl is type of methylation, found in files
- reads : name of file with subgroup data
- coverage : add coverage plot
- repeat.rows : rowname of repeat to show in repeat track
- focus = number of bases to show on each side of the repeat location in focussed plot
