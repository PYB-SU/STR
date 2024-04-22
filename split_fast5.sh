#! /bin/bash

ROOTDIR=$1

cd $ROOTDIR
mkdir FAST5_single

cd FAST5
FAST5FILES=$(ls *.fast5)

cd $ROOTDIR

echo "converting multi fast5 file to single fast5" >>$LOG
module load fast5_research

for FI in $FAST5FILES
do
    FILENOEXT=$(basename -s .fast5 $FI)
    multi_to_single_fast5 -i $ROOTDIR/FAST5/$FI -s $ROOTDIR/FAST5_single/$FILENOEXT
done

module unload fast5_research

cd $ROOTDIR

for DIR in $(ls FAST5_single)
do
    cd $ROOTDIR/FAST5_single/$DIR
    mkdir $ROOTDIR/FAST5_single/$DIR/1
    mkdir $ROOTDIR/FAST5_single/$DIR/2
    mkdir $ROOTDIR/FAST5_single/$DIR/3
    
    FAST5FILES=$(ls $ROOTDIR/FAST5_single/$DIR/0/*.fast5)
    count=0
    for FI in $FAST5FILES
    do
	NEWDIR=$(( $count % 4 ))
	if [ $NEWDIR -gt 0 ]
	then
	    mv $FI $ROOTDIR/FAST5_single/$DIR/$NEWDIR
	fi
	count=$(( $count + 1 ))
    done
done
