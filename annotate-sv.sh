#! /bin/bash

function usage()
{
    echo "usage: annotate-sv.sh [-r <in.bed>] [-b <in.bam>] <in.bp>"
    echo "  -b FILE     specify original BAM file (default: root name of bp file)"
    echo "  -r FILE     set region BED file"
    exit 0
}

# Check option
while getopts b:r:h option
do
    case ${option} in
	r)
	    ARG1=${OPTARG}
	    ;;
	b)
	    ARG2=${OPTARG}
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
SDIR="$0"
if [ -L $SDIR ] ; then
    SDIR=$(readlink $SDIR)
fi
REALLENROOT=$(cd $(dirname $SDIR); pwd)

REALLENDIR="$REALLENROOT/src"
ANNOTATION="$REALLENROOT/resource"
TMPDIR=`pwd`/temporary_files_reallen
SAMTLS=samtools

if [ $ARG1 ]
then
    BEDOPT="-b $ARG1"
fi

rt=`basename $1 .bp`

if [ $ARG2 ]
then
    COVBAM=$ARG2
else
    COVBAM=${rt}.bam
fi
    
if [ ! -e $TMPDIR ]
then
    mkdir $TMPDIR
fi

if [ -e $COVBAM ]
then
    if [ ! -e $TMPDIR/${rt}.bp.bed ]
    then
	$REALLENDIR/bp2bed.rb -s 1 -o $TMPDIR/${rt}.bp.bed $1
    fi
    echo
    echo calculate original coverage from $COVBAM
    $REALLENDIR/calc_coverage.rb $TMPDIR/${rt}.bp.bed $COVBAM -o $TMPDIR/${rt}.cov.csv --samtools $SAMTLS
    COVOPT="-c $TMPDIR/${rt}.cov.csv"
else
    echo
    echo "cannot find original bam file: $COVBAM"
    echo "skip coverage calculation"
fi

$REALLENDIR/bp2table.rb -d $COVOPT --annotation $ANNOTATION -o ${rt}.csv ${rt}.bp $BEDOPT
