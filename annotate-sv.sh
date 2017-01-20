#! /bin/bash

function usage()
{
    echo "usage: annotate-sv.sh [-b <in.bed>] [-d <in.bam>] <in.bp>"
    echo "  -B          output by BLAST-like format"
    echo "  -d FILE     specify original BAM file (default: <root>.bam)"
    echo "  -c FILE     specify config file [reallen.config]"
    echo "  -b FILE     set region BED file"
    echo "  -s          apply stringent coverage filter"
    echo "  -m          just mark filter-records (do not remove)"
    echo "  -o option   directly give other options (string) to bp2table"
    exit 0
}

# Check option
while getopts Bd:c:b:o:smh option
do
    case ${option} in
	B)
	    BLASTOPT="--blast"
	    ;;
	b)
	    ARG1=${OPTARG}
	    ;;
	d)
	    ARG2=${OPTARG}
	    ;;
	c)
	    ARG3=${OPTARG}
	    ;;
	s)
	    FILTOPT="--cov-filter 4 --cov-rate-filter 0.02"
	    ;;
	o)
	    RESTOPT=${OPTARG}
	    ;;
	m)
	    KEEPOPT="-m"
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
SCRIPTPATH=$(cd $(dirname $SDIR); pwd)

if [ -z $ARG3 ] ; then
    CONFIGFILE="$SCRIPTPATH/reallen.config"
else
    CONFIGFILE=$ARG3
fi
source $CONFIGFILE

if [ $ARG1 ]; then
    BEDOPT="-b $ARG1"
fi

rt=`basename $1 .bp`

if [ $ARG2 ]; then
    COVBAM=$ARG2
else
    COVBAM=${rt}.bam
fi
    
if [ ! -e $TMPDIR ]
then
    mkdir $TMPDIR
fi

if [ -e $TMPDIR/${rt}.cov.csv ]
then
    echo "$TMPDIR/${rt}.cov.csv was found"
    echo "use existing coverage bed file"
    COVOPT="-c $TMPDIR/${rt}.cov.csv"
else
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
fi

echo "$REALLENDIR/bp2table.rb -d -t 100000 $COVOPT --annotation $ANNOTATION -o ${rt}.csv ${rt}.bp $BEDOPT $FILTOPT $KEEPOPT $BLASTOPT $RESTOPT"
$REALLENDIR/bp2table.rb -d -t 100000 $COVOPT --annotation $ANNOTATION -o ${rt}.csv ${rt}.bp $BEDOPT $FILTOPT $KEEPOPT $BLASTOPT $RESTOPT
