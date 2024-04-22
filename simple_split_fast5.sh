#! /bin/bash

ROOTDIR=$1

cd $ROOTDIR
mkdir FAST5_single

cd FAST5
FAST5FILES=$(ls *.fast5)

cd $ROOTDIR

echo "converting multi fast5 file to single fast5"
module load fast5_research

for FI in $FAST5FILES
do
    FILENOEXT=$(basename -s .fast5 $FI)
    multi_to_single_fast5 -i $ROOTDIR/FAST5/$FI -s $ROOTDIR/FAST5_single/$FILENOEXT
done

module unload fast5_research
