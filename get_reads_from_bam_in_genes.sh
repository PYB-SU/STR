#! /bin/bash

ROOTDIR=$1
# here must find FAST5_single directory

cd $ROOTDIR

# where to store individual fast5 files
mkdir $ROOTDIR/extract_FAST5

ROOTDIR_NO_BAR=$(echo $ROOTDIR | tr '/' '_')

PATH=/softs/apps/fast5_research/4.1.2/bin:/softs/apps/samtools/1.16/bin:$PATH
export PATH

# find BAM files indexing reads in arorescence - -L to follow symlinks
# there may be several BAM files
BAMFILE=$(ls $ROOTDIR/guppy/BAM/*.bam)
echo "found $BAMFILE"
module load samtools


# check version of reference geneome for BAM files
# take first BAM FILE
FIRST_BAMFILE=$(echo $BAMFILE | cut -d " " -f 1)
COORDS_CHR10=$(samtools view -H $FIRST_BAMFILE | grep chr10 | awk '{print $3}')

echo "looking in $FIRST_BAMFILE"
echo "coords CHR10 : $COORDS_CHR10"

# T2T
if [ $COORDS_CHR10 ==  "LN:134758134" ]
then
    # coord of genes to look for
    coordCNBP="chr3:131912733-131928817"
    coordHTT="chr4:3040281-3085772"
    coordDMPK="chr19:48597244-48610042"
    coordGIPC1="chr19:14604393-14622766"
    coordLRP12="chr8:105616729-105716466"
    coordNOTCH2NLC="chr1:148519507-148596912"
    coordRILPL1="chr12:123468958-123532572"
elif [ $COORDS_CHR10 == "LN:133797422" ]
then
    #hg38
    coordCNBP="chr3:129167827-129183896"
    coordHTT="chr4:3074681-3243957"
    coordDMPK="chr19:45769717-45782490"
    coordGIPC1="chr19:14477762-14496127"
    coordLRP12="chr8:104489236-104589258"
    coordNOTCH2NLC="chr1:149390621-149471833"
    coordRILPL1="chr12:123470054-123533719"

else
    echo "check version of reference genome in BAM files $BAMFILE"
fi

# coord of genomes are OK
READSFILE=$(basename $ROOTDIR)
echo "read_id" >$ROOTDIR/$READSFILE
mkdir $ROOTDIR/tmpFAST5list
rm $ROOTDIR/tmpFAST5list/*
for BF in $BAMFILE
do
    BF_NO_DIR=$(basename $BF .bam)
    echo "working with $BF"
    samtools index $BF
    for GENE in GIPC1 LRP12 NOTCH2NLC CNBP DMPK HTT RILPL1
    do
	echo "doing $GENE"
	coordGENE=coord$GENE
	echo "samtools view $BF ${!coordGENE}"
	samtools view $BF ${!coordGENE} |awk -v ROOTDIR=$ROOTDIR -v GENE=$GENE '{ if ($2 == 0 || $2 == 16) { print $1}}' | tr -d '@' >> $ROOTDIR/tmpFAST5list/$BF_NO_DIR
    done
done

cat tmpFAST5list/* | sort | uniq >> $ROOTDIR/$READSFILE

cd $ROOTDIR/extract_FAST5
echo "removing files in $ROOTDIR/extract_FAST5"
rm -f *.fast5
cd $ROOTDIR



module load fast5_research

fast5_subset -i $ROOTDIR/FAST5 -s $ROOTDIR/extract_FAST5 -l $ROOTDIR/$READSFILE -f ${READSFILE}_

