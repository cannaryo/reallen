#! /bin/bash

function usage()
{
    echo "usage:ReALLEN [-SF] [-b <in.bed>] <in.bam>"
    echo "  -S          input file is SAM format instead of BAM format"
    echo "  -F          fast mode"
    echo "  -b FILE     set region file"
    exit 0
}

# Check option
while getopts FSb:h option
do
    case ${option} in
	b)
	    ARG1=${OPTARG}
	    ;;
	S)
	    ARG2=TRUE
	    ;;
	F)
	    ARG3=TRUE
	    ;;
	h)
	    usage
	    ;;
	*)
	    usage
	    ;;
    esac
done
shift $((OPTIND -1))

if [ $# -le 0 ]
then
    usage
fi

# Set variables

REALLENROOT="/home/clc/test/reallen"
REALLENDIR="$REALLENROOT/src"
REFHG19="$REALLENROOT/reference/hg19.fa"
REFHG19I="$REALLENROOT/reference/hg19.fa.fai"
BWAHG19="$REALLENROOT/reference/bwa-hg19/bwa-hg19"
ANNOTATION="$REALLENROOT/resource"
TMPDIR=`pwd`/tmp
SAMTLS=samtools

if [ $ARG1 ]
then
    BEDOPT="-b $ARG1"
fi

if [ ! -e $TMPDIR ]
then
    mkdir $TMPDIR
fi

if [ $ARG2 ]
then
    rt=`basename $1 .sam`
    FIRSTDATA=${rt}.sam
else
    rt=`basename $1 .bam`
    $SAMTLS view -h ${rt}.bam > $TMPDIR/${rt}.sam 
    FIRSTDATA=$TMPDIR/${rt}.sam
fi

if [ $ARG3 ]
then
    FILTOPT="--indel 10 -s 20"
else
    FILTOPT="--indel 8 -s 8 -u"
fi

# Run ReALLEN & BWA

$REALLENDIR/filtRead.rb $FILTOPT $FIRSTDATA -o $TMPDIR/${rt}.filt.sam
$REALLENDIR/splitSoftClip.rb $TMPDIR/${rt}.filt.sam -e $TMPDIR/${rt}.um.sam -o $TMPDIR/${rt}.sc.fq
echo
echo "performing BWA for ${rt}.sc.fq"
bwa mem -t 4 -O 3 -E 1 -T 20 -a $BWAHG19  $TMPDIR/${rt}.sc.fq > $TMPDIR/${rt}.sc_re.sam
$REALLENDIR/filtNonSpecific.rb -s 10 -f 0.7 $TMPDIR/${rt}.sc_re.sam -o $TMPDIR/${rt}.sc_ns.sam
$REALLENDIR/filtMapped.rb -l 500 -r $TMPDIR/${rt}.sc_ns.sam -o $TMPDIR/${rt}.sc_map.sam

$REALLENDIR/splitRead.rb -l 100 -s 45 $TMPDIR/${rt}.um.sam -o $TMPDIR/${rt}.um.fq
echo
echo "performing BWA for ${rt}.um.fq"
bwa mem -t 4 -O 3 -E 1 -T 35 -a $BWAHG19 $TMPDIR/${rt}.um.fq > $TMPDIR/${rt}.um_re.sam
$REALLENDIR/filtNonSpecific.rb -s 15 -f 0.6 $TMPDIR/${rt}.um_re.sam -o $TMPDIR/${rt}.um_ns.sam
$REALLENDIR/filtMapped.rb -l 500 -r $TMPDIR/${rt}.um_ns.sam -o $TMPDIR/${rt}.um_map.sam
$REALLENDIR/prepDDP.rb $TMPDIR/${rt}.um_map.sam -o $TMPDIR/${rt}.um_ddp.pos
$REALLENDIR/DDP.rb $REFHG19 -i $REFHG19I $TMPDIR/${rt}.um.sam $TMPDIR/${rt}.um_ddp.pos -o $TMPDIR/${rt}.um_ddp.sam

$REALLENDIR/mergeSAM.rb $TMPDIR/${rt}.sc_map.sam $TMPDIR/${rt}.um_ddp.sam -o $TMPDIR/${rt}.merge.sam
$REALLENDIR/filtBP.rb -r 30 $TMPDIR/${rt}.merge.sam -b $TMPDIR/${rt}.bp.bed -o $TMPDIR/${rt}.bp.sam -p ${rt}.bp

if [ -e ${rt}.bam ]
then
    echo
    echo calculate original coverage from ${rt}.bam
    $REALLENDIR/calc_coverage.sh $TMPDIR/${rt}.bp.bed ${rt}.bam $TMPDIR/${rt}.cov.csv
    COVOPT="-c $TMPDIR/${rt}.cov.csv"
else
    echo
    echo "original bam file cannot be specified"
    echo "skip coverage calculation"
    echo "do not use -S option or prepare ${rt}.bam for calculating coverage"
fi

$REALLENDIR/bp2table.rb -d $COVOPT \
    -g $ANNOTATION/hg19_gene_region.csv -e $ANNOTATION/hg19_exon_region.csv -i $ANNOTATION/cytoband.csv \
    -o ${rt}.csv ${rt}.bp $BEDOPT
