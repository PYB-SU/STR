#!/bin/bash

GENE=$1

if [ $GENE = "-h" ]
then
    echo "get_primer_and_flanking_for_gene.sh GENE REPEATLEFT REPEATRIGHT"
    echo "extract sequence from reference fasta file"
    echo "REPEATLEFT < REPEATRIGHT are coordinates of repeat in T2T"
    echo "primer sequences are flanking 100bp at +/-400 and +/-500 from (REPEATLEFT+REPEATPOS)/2"
    echo "flanking sequence is 20/25/35 bp flanking REPEATLEFT/REPEATRIGHT"
    
    echo "output in /scratch/ONT/STR/GENE"
    
    echo "file primer-[left|right]-400.fasta file primer-[left|right]-500.fasta"
    echo "file left-[left|right]-400.fasta file primer-[left|right]-500.fasta"
    exit
fi

REPEATleft=$2
REPEATright=$3

if [ $REPEATleft = "" ]
then
    "please provide repeat region coords in T2T"
    exit
fi

CHR_POS=`grep "gene_name=$GENE" /softs/references/Homo_sapiens.T2T/chm13.draft_v2.0.gene_annotation.gff3 | head -n 1 | awk '{print $1" "$4" "$5}'`
CHR=`echo $CHR_POS | awk '{print $1}'`
POSLEFT=`echo $CHR_POS | awk '{print $2}'`
POSRIGHT=`echo $CHR_POS | awk '{print $3}'`



echo "found $GENE on $CHR:$POSLEFT-$POSRIGHT"

if [ $REPEATleft -lt $POSLEFT ]
then
    echo "REPEATleft outside GENE coords"
fi
if [ $REPEATright -gt $POSRIGHT ]
then
    echo "REPEATright outside GENE coords"
fi

if [ -d /scratch/ONT/STR/$GENE ]
then
    echo "writing to /scratch/ONT/STR/$GENE"
else
    mkdir /scratch/ONT/STR/$GENE
    echo "writing to /scratch/ONT/STR/$GENE"
fi

cd /scratch/ONT/STR/$GENE

echo "$GENE $CHR $POSLEFT $POSRIGHT $REPEATleft $REPEATright" > $GENE.specs

module load bedtools
module load seqtk

for LEN in $(seq 50 50 1000)
do
    for DIR in left right
    do
	echo "primer +/-$LEN $DIR"
	echo $CHR $REPEATleft $REPEATright $LEN $DIR | awk '{mid=int($2+$3)/2; if ($5 == "left") {print $1"\t"sprintf("%d",mid-100-$4)"\t"sprintf("%d",mid-$4)} else if ($5=="right") {print $1"\t"sprintf("%d",mid+$4)"\t"sprintf("%d",mid+$4+100)}}' > /scratch/ONT/STR/$GENE/primers-$DIR-$LEN.bed
	echo "running bedtools getfasta -fi /softs/references/Homo_sapiens.T2T/chm13v2.0.fa -bed /scratch/ONT/STR/$GENE/primers-$DIR-$LEN.bed > /scratch/ONT/STR/$GENE/primers-$DIR-$LEN.fa"
	bedtools getfasta -fi /softs/references/Homo_sapiens.T2T/chm13v2.0.fa -bed /scratch/ONT/STR/$GENE/primers-$DIR-$LEN.bed | sed -e 's/\([a-z]\)/\U\1/g' > /scratch/ONT/STR/$GENE/primers-$DIR-$LEN.fa
    done
done

for LEN in 30 40 50 60 80 100
do
    for DIR in left right
    do
	echo "flank +/-$LEN $DIR"
	if [ $DIR = "left" ]
	then
	    REPEAT=$REPEATleft
	else
	    REPEAT=$REPEATright
	fi	
	echo $CHR $REPEAT $LEN $DIR | awk '{ if ($4 == "left") {print $1"\t"$2-$3"\t"$2} else if ($4=="right") {print $1"\t"$2"\t"$2+$3}}' > /scratch/ONT/STR/$GENE/flank-$DIR-$LEN.bed
	echo "running bedtools getfasta -fi /softs/references/Homo_sapiens.T2T/chm13v2.0.fa -bed /scratch/ONT/STR/$GENE/flank-$DIR-$LEN.bed > /scratch/ONT/STR/$GENE/flank-$DIR-$LEN.fa"
	bedtools getfasta -fi /softs/references/Homo_sapiens.T2T/chm13v2.0.fa -bed /scratch/ONT/STR/$GENE/flank-$DIR-$LEN.bed | sed -e 's/\([a-z]\)/\U\1/g' > /scratch/ONT/STR/$GENE/flank-$DIR-$LEN.fa
	seqtk seq -r  /scratch/ONT/STR/$GENE/flank-$DIR-$LEN.fa >  /scratch/ONT/STR/$GENE/flank-$DIR-${LEN}_rev.fa
    done
done


