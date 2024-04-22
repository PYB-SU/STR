#! /bin/bash

THREADS=8
REFS=/softs/references/Homo_sapiens.T2T/

HERE=$1

if [ ! -d JOBS ]
then
    mkdir $HERE/$JOBS
fi
   
LOG=$HERE/map.log
echo "writing to log">$LOG

if [ $HERE = "-h" ] || [ "x"$HERE = "x" ]
then
    echo -e "call_with_bonito.job HERE MODEL STEPS"
    echo -e "INPUT : HERE is root directory with subdirectory FASTQ"
    echo -e "OUTPUT : BAM FASTQ in separate subdirectories"
    echo -e  "\t\tuses tmp*fastq.gz files in HERE/FASTQ"    
    echo -e  "\t\tuse porechop with each"
    echo -e "\t\toutputs final fastq.gz to HERE/FASTQ"
    echo -e  "\t\tuse minimap2 to map files"
    echo -e "\t\toutputs bam files and index to HERE/BAM"
    echo -e  "submits one job to slurm for each file"
    exit
fi

SAMPLEDIR=$1

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

cd $HERE/FASTQ/
FASTQFILES=$(ls *.fastq.gz)

NBFASTQFILES=$(echo $FASTQFILES | wc -w)

EXT=".fastq.gz"

echo "found $NBFASTQFILES files with $EXT extension in ${SAMPLENAME}/FASTQ)" >>  $LOG
cd $HERE

afterPID=""
for FILEEXT in $FASTQFILES
do
    FILE=$(basename -s $EXT $FILEEXT)
    echo "working on $FILE">> $LOG
    
    JOBFILE=$HERE/JOBS/pipelineM${FILE}.job
    echo "#! /bin/bash" > $JOBFILE
    echo "#SBATCH --mem 60G" >> $JOBFILE
    echo "#SBATCH -n 8" >> $JOBFILE
    echo "#SBATCH --time 0:40:00" >> $JOBFILE
    echo "cd $HERE/FASTQ" >>$JOBFILE
    echo "module load minimap2">>$JOBFILE
    echo "module load samtools">>$JOBFILE
    echo "cd $HERE/BAM">>$JOBFILE
    echo "minimap2 -a -x map-ont -t $THREADS $REFS/chm13v2.0.mmi $HERE/FASTQ/${FILEEXT} | samtools sort -@ $THREADS -o $HERE/BAM/${FILE}.bam" >>$JOBFILE
    echo "samtools index  $HERE/BAM/${FILE}.bam" >> $JOBFILE
    echo "module unload minimap2">>$JOBFILE
    echo "module unload samtools">>$JOBFILE
    pID=$(sbatch --parsable $JOBFILE)
    echo "FILE $FILE in process $pID" >> $LOG
    afterPID=$afterPID":"$pID
done
echo $afterPID



