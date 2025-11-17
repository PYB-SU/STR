#!/bin/bash

GENE=$1

if [ "x"$GENE = "x" ]
then
    echo "provide GENE"
    exit 1
fi

PIPELINEDIR=/home/boelle/STR/PIPELINE/
REFS=/softs/cinbios/references/Homo_sapiens.T2T
GENOME_GFF=chm13.draft_v2.0.gene_annotation.gff3
GENOME_FA=chm13.v2.fa
if [ $GENE = "-h" ]
then
    echo "get_primer_and_flanking_for_gene.sh GENE REPEATLEFT REPEATRIGHT"
    echo "extract sequence from reference fasta file"
    echo "REPEATLEFT < REPEATRIGHT are coordinates of repeat in T2T"
    echo "primer sequences are flanking 100bp at +/-400 and +/-500 from (REPEATLEFT+REPEATPOS)/2"
    echo "flanking sequence is 20/25/35 bp flanking REPEATLEFT/REPEATRIGHT"
    
    echo "output in $PIPELINEDIR/GENE"
    
    echo "file primer-[left|right]-400.fasta file primer-[left|right]-500.fasta"
    echo "file left-[left|right]-400.fasta file primer-[left|right]-500.fasta"
    exit
fi

REPEATleft=$2
REPEATright=$3



if [ "x"$REPEATright = "x" ]
then
    "please provide repeat region coords in T2T"
    exit 2
fi

CHR_POS=$(grep "gene_name=$GENE" $REFS/$GENOME_GFF | head -n 1 | awk '{print $1" "$4" "$5}')
CHR=$(echo $CHR_POS | awk '{print $1}')
POSLEFT=$(echo $CHR_POS | awk '{print $2}')
POSRIGHT=$(echo $CHR_POS | awk '{print $3}')

echo "found $GENE on $CHR:$POSLEFT-$POSRIGHT"

if [ $REPEATleft -lt $POSLEFT ]
then
    echo "REPEATleft outside GENE coords"
fi
if [ $REPEATright -gt $POSRIGHT ]
then
    echo "REPEATright outside GENE coords"
fi

if [ -d $PIPELINEDIR/$GENE ]
then
    echo "writing to $PIPELINEDIR/$GENE"
else
    mkdir $PIPELINEDIR/$GENE
    echo "writing to /scratch/ONT/STR/$GENE"
fi

cd $PIPELINEDIR/$GENE

echo "$GENE $CHR $POSLEFT $POSRIGHT $REPEATleft $REPEATright" > $GENE.specs

module load bedtools
module load seqtk

for LEN in $(seq 50 50 1000)
do
    for DIR in left right
    do
	echo "primer +/-$LEN $DIR"
	echo $CHR $REPEATleft $REPEATright $LEN $DIR | awk '{mid=int($2+$3)/2; if ($5 == "left") {print $1"\t"sprintf("%d",mid-100-$4)"\t"sprintf("%d",mid-$4)} else if ($5=="right") {print $1"\t"sprintf("%d",mid+$4)"\t"sprintf("%d",mid+$4+100)}}' > $PIPELINEDIR/$GENE/primers-$DIR-$LEN.bed
	echo "running bedtools getfasta -fi $REFS/$GENOME_FA -bed $PIPELINEDIR/$GENE/primers-$DIR-$LEN.bed > $PIPELINEDIR/$GENE/primers-$DIR-$LEN.fa"
	bedtools getfasta -fi $REFS/$GENOME_FA -bed $PIPELINEDIR/$GENE/primers-$DIR-$LEN.bed | sed -e 's/\([a-z]\)/\U\1/g' > $PIPELINEDIR/$GENE/primers-$DIR-$LEN.fa
	seqtk seq -r  $PIPELINEDIR/$GENE/primers-$DIR-$LEN.fa >  $PIPELINEDIR/$GENE/primers-$DIR-$LEN-rev.fa
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
	echo $CHR $REPEAT $LEN $DIR | awk '{ if ($4 == "left") {print $1"\t"$2-$3"\t"$2} else if ($4=="right") {print $1"\t"$2"\t"$2+$3}}'
	echo $CHR $REPEAT $LEN $DIR | awk '{ if ($4 == "left") {print $1"\t"$2-$3"\t"$2} else if ($4=="right") {print $1"\t"$2"\t"$2+$3}}' > $PIPELINEDIR/$GENE/flank-$DIR-$LEN.bed
	cat $PIPELINEDIR/$GENE/flank-$DIR-$LEN.bed
	echo "running bedtools getfasta -fi $REFS/$GENOME_FA -bed $PIPELINEDIR/$GENE/flank-$DIR-$LEN.bed > $PIPELINEDIR/$GENE/flank-$DIR-$LEN.fa"
	bedtools getfasta -fi $REFS/$GENOME_FA -bed $PIPELINEDIR/$GENE/flank-$DIR-$LEN.bed | sed -e 's/\([a-z]\)/\U\1/g' > $PIPELINEDIR/$GENE/flank-$DIR-$LEN.fa
	seqtk seq -r  $PIPELINEDIR/$GENE/flank-$DIR-$LEN.fa >  $PIPELINEDIR/$GENE/flank-$DIR-${LEN}_rev.fa
    done
done


module unload bedtools
module unload seqtk
