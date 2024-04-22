#! /bin/bash

THREADS=32
REFS=/softs/references/Homo_sapiens.T2T/

HERE=$1

if [ ! -d JOBS ]
then
    mkdir $HERE/$JOBS
fi
   

LOG=$HERE/call_trim_zip_map.log
echo "writing to log">$LOG

if [ $HERE = "-h" ] || [ "x"$HERE = "x" ]
then
    echo -e "call_with_bonito.job HERE MODEL STEPS"
    echo -e "INPUT : HERE is root directory with subdirectory FAST5"
    echo -e "OUTPUT : BAM FASTQ in separate subdirectories"
    echo -e  "\t\tuse fast5 files in HERE/FAST5"
    echo -e  "\t\trun bonito with each using MODEL"
    echo -e  "\t\toutputs tmp*fastq.gz files to HERE/FASTQ"    
    echo -e  "\t\tif MODEL is default -> dna_r9.4.1_e8.1_fast@v3.4"
    echo -e  "\t\tuse tmp*fastq.gz files in HERE/FASTQ"
    echo -e  "\t\tuse porechop with each"
    echo -e "\t\toutputs final fastq.gz to HERE/FASTQ"
    echo -e  "submits one job to slurm for each file"
    echo -e  "models are : "
    echo -e  "\t1) FAST10 dna_r10.4_e8.1_fast@v3.4"
    echo -e  "\t2) HAC10 dna_r10.4_e8.1_hac@v3.4"
    echo -e  "\t3) SUP10 dna_r10.4_e8.1_sup@v3.4"
    echo -e  "\t4) FAST81 dna_r9.4.1_e8.1_fast@v3.4"
    echo -e  "\t5) HAC81 dna_r9.4.1_e8.1_hac@v3.3"
    echo -e  "\t6) SUP81 dna_r9.4.1_e8.1_sup@v3.3"
    echo -e  "\t7) FAST8 dna_r9.4.1_e8_fast@v3.4"
    echo -e  "\t8) HAC8 dna_r9.4.1_e8_hac@v3.3"
    echo -e  "\t9) SUP8 dna_r9.4.1_e8_sup@v3.3"
    echo -e  "\t10) SUP8 dna_r9.4.1_telomere"    
    exit
fi

SAMPLEDIR=$(dirname $1)

if [ -d $HERE/FAST5 ]
then
    echo "using fast5 files in $HERE/FAST5">>$LOG
else
    echo "$HERE/FAST5 directory does not exist">>$LOG
    echo "should be linked to $SAMPLEDIR/FAST5_single">>$LOG
    exit
fi

if [ -d $HERE/FASTQ ]
then
    echo "removing old $HERE/FASTQ">>$LOG
    rm -rf $HERE/FASTQ
fi
mkdir $HERE/FASTQ

if [ -d $HERE/BAM ]
then
    echo "removing old $HERE/BAM">>$LOG
    rm -rf $HERE/BAM
fi
mkdir $HERE/BAM

echo "working in $HERE">>$LOG

echo "source files in $HERE/FAST5" >> $LOG

MOD=$2

if [ "x"$MOD = "x" ]
then
    echo "bonito MODEL must be provided"
    exit
elif [ $MOD = 1 ] | [ $MOD = "FAST10" ] ; then
    MODEL=dna_r10.4_e8.1_fast@v3.4
elif [ $MOD = 2 ] | [ $MOD = "HAC10" ] ; then
    MODEL=dna_r10.4_e8.1_hac@v3.4
elif [ $MOD = 3 ] | [ $MOD = "SUP10" ] ; then
    MODEL=dna_r10.4_e8.1_sup@v3.4
elif [ $MOD = 4 ] | [ $MOD = "FAST81" ] ; then
    MODEL=dna_r9.4.1_e8.1_fast@v3.4
elif [ $MOD = 5 ] | [ $MOD = "HAC81" ]; then
    MODEL=dna_r9.4.1_e8.1_hac@v3.3
elif [ $MOD = 6 ] | [ $MOD = "SUP81" ]; then
    MODEL=dna_r9.4.1_e8.1_sup@v3.3
elif [ $MOD = 7 ] | [ $MOD = "FAST8" ]; then
    MODEL=dna_r9.4.1_e8_fast@v3.4
elif [ $MOD = 8 ] | [ $MOD = "HAC8" ]; then
    MODEL=dna_r9.4.1_e8_hac@v3.3
elif [ $MOD = 9 ] | [ $MOD = "SUP8" ]; then
    MODEL=dna_r9.4.1_e8_sup@v3.3
elif [ $MOD = 10 ] | [ $MOD = "TELO" ]; then
    MODEL=dna_r9.4.1_telomere
fi

echo "using MODEL : $MODEL" >> $LOG
EXT=".fast5"

cd $HERE
FAST5FILES=$(ls $SAMPLEDIR/FAST5/*.fast5)

# we test whether there is a FAST5_single directory
if [ ! -d ../FAST5_single ]
then
    mkdir ../FAST5_single
fi

NBFAST5FILES=$(echo $FAST5FILES | wc -l)

echo "found $NBFAST5FILES files with $EXT extension in $SAMPLENAME/FAST5)" >>  $LOG
cd $HERE

afterPID=""
for FILEEXT in $FAST5FILES
do
    FILE=$(basename -s $EXT $FILEEXT)
    echo "working on $FILE">> $LOG
    if [ ! -d $HERE/FAST5/$FILE ]
    then
	echo "cannot find $HERE/FAST5/$FILE"
	echo "you should run split_fast5.sh first"
	exit
    fi
    
    for SUBDIR in {0..3}
    do
	JOBFILE=$HERE/JOBS/pipeline${FILE}_$SUBDIR.job
	echo "#! /bin/bash" > $JOBFILE
	echo "#SBATCH --mem 250G" >> $JOBFILE
	echo "#SBATCH -n 32" >> $JOBFILE
	echo "#SBATCH --time 0:40:00" >> $JOBFILE
	echo "module load bonito" >>$JOBFILE
	echo "module load htslib" >>$JOBFILE
	echo "bonito basecaller --use_openvino --device cpu $MODEL $HERE/FAST5/${FILE}/$SUBDIR | bgzip -@ $THREADS > $HERE/FASTQ/tmp${FILE}_$SUBDIR.fastq.gz" >> $JOBFILE
	echo "module unload bonito" >>$JOBFILE
	echo "module unload htslib" >>$JOBFILE
	echo "cd $HERE/FASTQ" >>$JOBFILE
	echo "module load porechop" >>$JOBFILE
	echo "porechop-runner.py -t $THREADS -i $HERE/FASTQ/tmp${FILE}_$SUBDIR.fastq.gz -o $HERE/FASTQ/${FILE}_$SUBDIR.fastq.gz">>$JOBFILE
#    echo "rm -f $HERE/FASTQ/tmp$FILE.fastq.gz" >>$JOBFILE 
	echo "module unload porechop" >>$JOBFILE
	pID=$(sbatch --parsable $JOBFILE)
	echo "FILE $FILE in process $pID" >> $LOG
	afterPID=$afterPID":"$pID
    done
done
echo $afterPID



