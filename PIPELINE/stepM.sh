#! /bin/bash

THREADS=4
REFS=/softs/cinbios/references/Homo_sapiens.T2T/

HERE=$1
GO=$2
DEP=$3

echo "going to $HERE"
cd $HERE
if [ ! -d JOBS ]
then
    mkdir $HERE/JOBS
fi
   
LOG=$HERE/map.log
echo "writing to log">$LOG

if [ $HERE = "-h" ] || [ "x"$HERE = "x" ]
then
    echo -e "call_with_bonito.job HERE MODEL STEPS"
    echo -e "INPUT : HERE is root directory with subdirectory FASTQ"
    echo -e "OUTPUT : BAM in separate subdirectories"
    echo -e  "\t\tuses *fastq.gz files in HERE/FASTQ"    
    echo -e  "\t\tuse minimap2 to map files"
    echo -e "\t\toutputs bam files and index to HERE/BAM"
    echo -e  "submits one job to slurm for each file"
    exit
fi

SAMPLEDIR=$(dirname $1)
echo "extracted $SAMPLEDIR from $HERE">>$LOG

if [ -d $HERE/FASTQ ]
then
    echo "using fastq files in $HERE/FASTQ">>$LOG
else
    echo "$HERE/FASTQ directory does not exist">>$LOG
    echo "should call fast5 files first">>$LOG
    exit
fi

if [ -d $HERE/BAM ]
then
    echo "removing old $HERE/BAM">>$LOG
    rm -rf $HERE/BAM
fi
mkdir $HERE/BAM

echo "working in $HERE">>$LOG

echo "source files in $HERE/FASTQ" >> $LOG

cd $HERE

FASTQFILES=$(ls $HERE/FASTQ/*.fastq*)

NBFASTQFILES=$(ls -1 $HERE/FASTQ/*.fastq* | wc -l)

FIRSTFILE=$(ls -1 $HERE/FASTQ/*.fastq* | head -n1)
EXT="${FIRSTFILE#*.}"

echo "found $NBFASTQFILES files with $EXT extension in $HERE/FASTQ)" >>  $LOG
cd $HERE

afterPID=""
for FILEEXT in $FASTQFILES
do
    FILE=$(basename -s $EXT $FILEEXT)
    echo "working on $FILE">> $LOG
    JOBFILE=$HERE/JOBS/pipelineM${FILE}job
    echo "putting this in $JOBFILE" >> $LOG
    echo "#! /bin/bash" > $JOBFILE
    echo "#SBATCH --mem 30G" >> $JOBFILE
    echo "#SBATCH --nodes 1" >> $JOBFILE
    echo "#SBATCH --tasks 1" >> $JOBFILE
    echo "#SBATCH --cpus-per-task 4" >> $JOBFILE
    echo "#SBATCH --time 0:40:00" >> $JOBFILE
    echo "cd $HERE/FASTQ" >>$JOBFILE
    echo "module load minimap2">>$JOBFILE
    echo "module load samtools">>$JOBFILE
    echo "cd $HERE/BAM">>$JOBFILE
    echo "minimap2 -a -x map-ont -t $THREADS $REFS/chm13v2.0.mmi $HERE/FASTQ/${FILE}$EXT | samtools sort -@ $THREADS -o $HERE/BAM/${FILE}bam" >>$JOBFILE
    echo "samtools index  $HERE/BAM/${FILE}bam" >> $JOBFILE
    echo "module unload minimap2">>$JOBFILE
    echo "module unload samtools">>$JOBFILE
    if [ "x"$GO != "x" ]
    then
	pID=$(sbatch --parsable $DEP $JOBFILE)
	echo "FILE $FILE in process $pID" >> $LOG
	afterPID=$afterPID":"$pID
    else
	echo "not run" >> $LOG
    fi
done
echo $afterPID



