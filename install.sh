#! /bin/bash
echo -------- Installing ReALLEN --------
CDIR="$0"
echo $CDIR
if [ -L $CDIR ] ; then
    CDIR=$(readlink $CDIR)
fi
echo $CDIR
CDIR=$(cd $(dirname $CDIR); pwd)

echo "DIR=${CDIR}"
echo "Checking scripts:"
LIST="AlignDDP.rb calc_coverage.rb DDP.rb filtRead.rb mergeSAM.rb SAMReader.rb splitSoftClip.rb CSVReader.rb FastqReader.rb filtMapped.rb prepDDP.rb SeqTools.rb bp2table.rb filtBP.rb filtNonSpecific.rb GenomeAnnotation.rb splitRead.rb"
for i in $LIST
do
    echo "$i"
    if [ ! -e $CDIR/src/$i ]
    then
	echo "File $CDIR/src/$i is missing"
	exit
    fi
done
echo Checking Reference:
LIST="bwa-hg19.amb bwa-hg19.ann bwa-hg19.bwt bwa-hg19.pac bwa-hg19.sa hg19.fa hg19.fa.fai"
for i in $LIST
do
    echo "$i"
    if [ ! -e $CDIR/reference/$i ]
    then
	echo "File $CDIR/reference/$i is missing"
	exit
    fi
done
echo Checking Resource:
LIST="cytoband.csv fusion_cosmic.bed hg19_exon_region.csv hg19_gene_region.csv"
for i in $LIST
do
    echo "$i"
    if [ ! -e $CDIR/resource/$i ]
    then
	echo "File $CDIR/resource/$i is missing"
	exit
    fi
done
echo Checking BWA
if [ -x $(which bwa) ]
then
    echo OK
else
    echo "BWA is not found"
    echo "Please specify the name in reallen.env" 
fi
echo Checking Samtools
if [ -x $(which samtools) ]
then
    echo OK
else
    echo "Samtools is not found"
    echo "Please specify the name in reallen.env" 
fi
