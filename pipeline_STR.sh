#! /bin/bash

# submits all treatment chain
# basecalling with bonito from fast5 files

if [ $1 = "-h" ]
then
    echo "sh pipeline_STR.sh DIR GENE WHATTODO BONITO_MODEL MEDAKA_MODEL"
    echo "DIR must contain FAST5 subdir"
    echo "GENE must be already extracted to STR"
    echo "WHATTODO is :"
    echo -e "\t(B)onito for basecalling"
    echo -e "\t(P)orechop+zip+minimap"
    echo -e "\t(O)verlap for extracting bam"
    echo -e "\t(R)eads extraction"
    echo -e "\t(C)harONT"
    echo -e "\t(S)equence from reads"
    echo -e "\t(M)ethylation from fast5"
    echo "MODEL is : "
    echo -e "\tHAC10 - uses dna_r10.4_e8.1_hac@v3.4 (bonito) and r1041_e82_400bps_hac (medaka)"
    echo -e "\tFAST10 - uses dna_r10.4_e8.1_fast@v3.4 (bonito) and r1041_e82_400bps_hac (medaka)"
    echo -e "\tSUP10 - uses dna_r10.4_e8.1_sup@v3.4 (bonito) and  r1041_e82_400bps_sup_g615 (medaka)"
    echo -e "\tHAC - uses dna_r9.4.1_e8.1_hac@v3.3 (bonito) and r941_e81_hac_g514 (medaka)"
    echo -e "\tFAST - uses dna_r9.4.1_e8.1_fast@v3.4 (bonito) and r941_e81_fast_g514 (medaka)"
    echo -e "\tSUP - uses dna_r9.4.1_e8.1_sup@v3.3 (bonito) and r941_e81_sup_g514 (medaka)"
    echo -e "\tMIN - uses dna_r9.4.1_e8_fast@v3.4 (bonito) and r941_min_fast_g507 (medaka)"
    echo -e "\tTELO - uses dna_r9.4.1_telomere (bonito) and r941_e81_sup_g514 (medaka)"
    exit
fi


# 

HERE=$1 

if [ "x"$HERE = "x" ]
then
    echo "provide ROOT directory"
    exit
fi

cd $HERE
echo "cd to $HERE"

GENE=$2
echo "GENE is $GENE"

if [ "x"$2 = "x" ]
then
    echo "provide GENE"
    exit
fi

WHATTODO=$3
if [ "x"$3 = "x" ]
then
    echo "provide WHATTODO (BPORCSM)"
    exit
else
    echo "doing $WHATTODO"
fi

do_BONITO=$(echo $WHATTODO | grep "B" | wc -l)
do_MINIMAP=$(echo $WHATTODO | grep "P" | wc -l)
do_OVERLAP=$(echo $WHATTODO | grep "O" | wc -l)
do_READS=$(echo $WHATTODO | grep "R" | wc -l)
do_CHARONT=$(echo $WHATTODO | grep "C" | wc -l)
do_SEQUENCES=$(echo $WHATTODO | grep "S" | wc -l)
do_METHYLATION=$(echo $WHATTODO | grep "M" | wc -l)

PIPELINEDIR=/scratch/ONT/STR

BONITO_MODEL=$4
if [ "x"$BONITO_MODEL = "x" ]
then
    echo "no BONITO MODEL provided."
    echo "choices are HAC10|FAST10|SUP10|HAC81|FAST81|SUP81|HAC8|FAST8|SUP8"
    exit
elif [ $BONITO_MODEL = "HAC10" ]
then
    echo "using BONITO_MODEL HAC10 dna_r10.4_e8.1_hac@v3.4"
elif [ $BONITO_MODEL = "FAST10" ]
then
    echo "using BONITO_MODEL FAST10 dna_r10.4_e8.1_fast@v3.4"
elif [ $BONITO_MODEL = "SUP10" ]
then
    echo "using BONITO_MODEL SUP10 dna_r10.4_e8.1_sup@v3.4"
elif [ $BONITO_MODEL = "HAC81" ]
then
    echo "using BONITO_MODEL HAC81 dna_r9.4.1_e8.1_hac@v3.3"
elif [ $BONITO_MODEL = "FAST81" ]
then
    echo "using BONITO_MODEL FAST81 dna_r9.4.1_e8.1_fast@v3.4"
elif [ $BONITO_MODEL = "SUP81" ]
then
    echo "using BONITO_MODEL SUP81 dna_r9.4.1_e8.1_sup@v3.3"
elif [ $BONITO_MODEL = "HAC8" ]
then
    echo "using BONITO_MODEL HAC8 dna_r9.4.1_e8_hac@v3.3"
elif [ $BONITO_MODEL = "FAST8" ]
then
    echo "using BONITO_MODEL FAST8 dna_r9.4.1_e8_fast@v3.4"
elif [ $BONITO_MODEL = "SUP8" ]
then
    echo "using BONITO_MODEL SUP8 dna_r9.4.1_e8_sup@v3.3"
elif [ $BONITO_MODEL = "TELO" ]
then
    echo "using BONITO_MODEL TELO dna_r9.4.1_telomere"
else
    echo "unknown BONITO_MODEL $BONITO_MODEL"
    echo "change in call_zip_trim.sh"
    exit
fi

MEDAKA_MODEL=$5
if [ "x"$MEDAKA_MODEL = "x" ]
then
    echo "no MEDAKA_MODEL provided."
    echo "choices are HAC10|HAC|FAST|SUP|MIN|PROM_HAC|PROM_SUP"
    exit
elif [ $MEDAKA_MODEL = "HAC10" ]
then
    echo "using MEDAKA_MODEL r1041_e81_hac_g514"
elif [ $MEDAKA_MODEL = "HAC" ]
then
    echo "using MEDAKA_MODEL r941_e81_hac_g514"
elif [ $MEDAKA_MODEL = "FAST" ]
then
    echo "using MEDAKA_MODEL r941_e81_fast_g514"
elif [ $MEDAKA_MODEL = "SUP" ]
then
    echo "using MEDAKA_MODEL r941_e81_sup_g514"
elif [ $MEDAKA_MODEL = "MIN" ]
then
    echo "using MEDAKA_MODEL r941_min_hac_g507"
elif [ $MEDAKA_MODEL = "PROM_HAC" ]
then
    echo "using MEDAKA_MODEL r941_prom_hac_g507"
elif [ $MEDAKA_MODEL = "PROM_SUP" ]
then
    echo "using MEDAKA_MODEL r941_prom_sup_g507"
else
    echo "unknown MODEL $MODEL"
    echo "change in config_CharONT.R"
    exit
fi

MOD_NAME=$6
if [ "x"$MOD_NAME = "x" ]
then
    if [ $do_METHYLATION = "1" ]
    then
	echo "provide MOD_NAME for METHYLATION step"
	echo "accepted values are 5mCG, 5mCG_5hmCG, 6mA"
	exit
    fi
elif [ $MOD_NAME = "5mCG" ]
then
    echo "doing 5mCG methylation"
elif [ $MOD_NAME = "5mCG_5hmCG" ]
then
    echo "doing 5mCG_5hmCG methylation"
elif [ $MOD_NAME = "6mA" ]
then
    echo "doing 6mA methylation"
else
    echo "unknown model $MOD_NAME for methylation"
fi





if [ -d /scratch/ONT/STR/$GENE ]
then
    echo "information on $GENE seems to exist in /scratch/ONT/STR/$GENE"
else
    echo "information on $GENE not found in /scratch/ONT/STR"
    echo "you must run :"
    echo "/scratch/ONT/STR/get_primer_and_flanking_for_gene.sh $GENE REPEAT_POS"
fi

echo "pipeline started on $(date)" > $HERE/pipeline.log

LAST5=$(echo $HERE | sed -n 's/.*pilote\([^/]*\).*/\1/p')

if [[ $do_BONITO == "1" ]]
then
    echo "making BONITO step"
    echo  "bash $PIPELINEDIR/call_trim_zip.sh $HERE $MODEL"
    if [ ! -d $HERE/JOBS ]
    then
	mkdir $HERE/JOBS
    fi
    #BONITO_PID=$(bash $PIPELINEDIR/call.sh $HERE $BONITO_MODEL)
    BONITO_PID=$(bash $PIPELINEDIR/call_trim_zip_map.sh $HERE $BONITO_MODEL)
    DEP_BONITO="--dependency=afterok$BONITO_PID"
    echo "BONITO step submitted ($BONITO_PID)" >>$HERE/pipeline.log
else
    DEP_BONITO=""
fi

if [[ $do_MINIMAP == "1" ]]
then
    echo "making MINIMAP step"
    echo  "bash $PIPELINEDIR/map.sh $HERE $MODEL $DEP_BONITO"
    if [ ! -d $HERE/JOBS ]
    then
	mkdir $HERE/JOBS
    fi
    MINIMAP_PID=$(bash $PIPELINEDIR/map.sh $HERE $DEP_BONITO)
    DEP_MINIMAP="--dependency=afterok$MINIMAP_PID"
    echo "MINIMAP step submitted ($MINIMAP_PID)" >>$HERE/pipeline.log
else
    DEP_MINIMAP=""
fi

if [[ $do_OVERLAP == "1" ]]
then
    echo "making OVERLAP step"
    echo "sbatch -J o$LAST5$GENE --parsable $DEP_MINIMAP $PIPELINEDIR/extract_overlap_GENE.job $HERE $GENE"
    OVERLAP_PID=$(sbatch -J o$LAST5$GENE --parsable $DEP_MINIMAP $PIPELINEDIR/extract_overlap_GENE.job $HERE $GENE)
    DEP_OVERLAP="--dependency=afterok:$OVERLAP_PID"
    echo "OVERLAP step submitted ($OVERLAP_PID)" >>$HERE/pipeline.log
else
    DEP_OVERLAP=""
fi


# here we aim to find the pos and cutoff that maximise 
if [[ $do_READS == "1" ]]
then
    echo "making READS step"
    echo "sbatch  -J r$LAST5$GENE --parsable $DEP_OVERLAP $PIPELINEDIR/find_best_overlaping_repeat_primer.job $HERE $GENE"
    READS_PID=$(sbatch  -J r$LAST5$GENE --parsable $DEP_OVERLAP $PIPELINEDIR/find_best_overlaping_primer.job $HERE $GENE)
    DEP_READS="--dependency=afterok:$READS_PID"
    echo "READS step submitted ($READS_PID)" >>$HERE/pipeline.log
else
    DEP_READS=""
fi

if [[ $do_CHARONT == "1" ]]
then
    echo "making CHARONT step"
    echo "sbatch  -J c$LAST5$GENE --parsable $DEP_READS  $PIPELINEDIR/charONT.job $HERE $GENE $MODEL"
    CHARONT_PID=$(sbatch  -J c$LAST5$GENE --parsable  $DEP_READS $PIPELINEDIR/charONT.job $HERE $GENE $MEDAKA_MODEL)
    DEP_CHARONT="--dependency=afterok:$CHARONT_PID"
    echo "CHARONT step submitted ($CHARONT_PID)" >>$HERE/pipeline.log
else
    DEP_CHARONT=""
fi

if [[ $do_SEQUENCES == "1" ]]
then
    echo "making SEQUENCES step"
    echo "sbatch  -J s$LAST5$GENE --parsable $DEP_CHARONT  $PIPELINEDIR/repeat_sequence_from_each_read.job $HERE $GENE 3"
    SEQUENCES_PID=$(sbatch  -J s$LAST5$GENE --parsable $DEP_CHARONT $PIPELINEDIR/repeat_sequence_from_each_read.job $HERE $GENE 3)
    echo "SEQUENCES step submitted ($SEQUENCES_PID)" >>$HERE/pipeline.log
fi

if [[ $do_METHYLATION == "1" ]]
then
    echo "making METHYLATION step"
    echo "sbatch  -J m$LAST5$GENE --parsable  $DEP_CHARONT $PIPELINEDIR/extract_reads_for_methylation_NEW.job $HERE $GENE $MOD_NAME"
    METHYLATION_PID=$(sbatch  -J m$LAST5$GENE --parsable  $DEP_CHARONT $PIPELINEDIR/extract_reads_for_methylation_modkit.job $HERE $GENE $MOD_NAME)
    echo "METHYLATION step submitted ($METHYLATION_PID)" >>$HERE/pipeline.log
fi




