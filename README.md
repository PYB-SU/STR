# Identifying repeat counts 

## 1) pipelineSTR.sh
sh pipeline_STR.sh DIR GENE WHATTODO BONITO_MODEL MEDAKA_MODEL MODIF

DIR base directory containing a FAST5 subdir

GENE of interest. GENE information must be available in directory STR (using get_primer_and_flanking.sh)

WHATTODO is a character chain with:
-        (B)onito for basecalling
-        (P)orechop+zip+minimap2
-        (O)verlap for extracting bam
-        (R)eads extraction
-        (C)harONT
-        (S)equence from reads
-        (M)ethylation from fast5

BONITO_MODEL and MEDAKA_MODEL is :
-        HAC10 - uses dna_r10.4_e8.1_hac@v3.4 (bonito) and r1041_e82_400bps_hac (medaka)
-        FAST10 - uses dna_r10.4_e8.1_fast@v3.4 (bonito) and r1041_e82_400bps_hac (medaka)
-        SUP10 - uses dna_r10.4_e8.1_sup@v3.4 (bonito) and  r1041_e82_400bps_sup_g615 (medaka)
-        HAC - uses dna_r9.4.1_e8.1_hac@v3.3 (bonito) and r941_e81_hac_g514 (medaka)
-        FAST - uses dna_r9.4.1_e8.1_fast@v3.4 (bonito) and r941_e81_fast_g514 (medaka)
-        SUP - uses dna_r9.4.1_e8.1_sup@v3.3 (bonito) and r941_e81_sup_g514 (medaka)

MODIF specifies the type of methylation :
-  5m : 5mc(bonito) , 5mCG (dorado)
-  5m+5hmc : 5mc_5hmc (bonito) , 5mCG_5hmCG (dorado)
-  6mA : 6mA (dorado)
