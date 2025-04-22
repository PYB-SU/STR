
PIPELINEDIR=/home/boelle/STR/PIPELINE
ROOTDIR=/store/cinbios/ONT

FILE=$1
if [ "x"$1 = "x" ]
then
    echo "provide FILE"
    exit 1
fi

GENE=$2
if [ "x"$2 = "x" ]
then
    echo "provide GENE"
    exit 2
fi

STEPS=$3
if [ "x"$3 = "x" ]
then
    echo "provide STEPS"
    exit 3
fi

CALLER=$4
if [ "x"$4 = "x" ]
then
    echo "provide caller"
    exit 4
fi

GO=$5
if [ "x"$5 = "x" ]
then
    echo "files will not be submitted"
else
    echo "files will be submitted"
fi

while IFS="" read -r DIR
do
    sh ${PIPELINEDIR}/pipeline_STR.sh ${ROOTDIR}/$DIR/$CALLER $GENE $STEPS HAC10 HAC10 all $GO
done < $FILE

    
