#! /bin/bash

function usage()
{
    echo "usage: reallen-fast [-p <NP>][-s <size>][-c <FILE>] <in.bam>"
    echo "  -p NUMBER   Number of processor [1]"
    echo "  -c FILE     specify config file [reallen.config]"
    echo "  -s SIZE     Minimum SV size [40]"
    echo "  -h          Show help message and exit"
    exit 0
}

# Check option
while getopts p:c:s:h option
do
    case ${option} in
	p)
	    ARG1=${OPTARG}
	    ;;
	c)
	    ARG2=${OPTARG}
	    ;;
	s)
	    ARG3=${OPTARG}
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

if [ -z $ARG2 ] ; then
    CONFIGFILE="$SCRIPTPATH/reallen.config"
else
    CONFIGFILE=$ARG2
fi
source $CONFIGFILE

if [ ! -e $TMPDIR ]
then
    mkdir $TMPDIR
fi

rt=`basename $1 .bam`

if [ $ARG1 ]
then
    PROCOPT="-t ${ARG4}"
fi

if [ $ARG3 ]
then
    SVOPT="-l $ARG3"
else
    SVOPT="-l 40"
fi

# Run ReALLEN & BWA
$BAMFILTER -u -p -b 8 -s 8 ${rt}.bam --softclip $TMPDIR/${rt}.sc.fq -o $TMPDIR/${rt}.um.sam --fixed $TMPDIR/${rt}.um.fq -f 45 -l 100 -k
echo
echo "performing BWA for ${rt}.sc.fq"
$BWA mem ${PROCOPT} -O3 -E1 -T20 -a $BWAHG19 $TMPDIR/${rt}.sc.fq > $TMPDIR/${rt}.sc_re.sam
$REALLENDIR/filtNonSpecific.rb -s 10 -f 0.7 $TMPDIR/${rt}.sc_re.sam -o $TMPDIR/${rt}.sc_ns.sam
$REALLENDIR/pairTest.rb $SVOPT -r $TMPDIR/${rt}.sc_ns.sam -o $TMPDIR/${rt}.sc_map.sam

echo
echo "performing BWA for ${rt}.um.fq"
$BWA mem ${PROCOPT} -O3 -E1 -T35 -a $BWAHG19 $TMPDIR/${rt}.um.fq > $TMPDIR/${rt}.um_re.sam
$REALLENDIR/filtNonSpecific.rb -s 15 -f 0.6 $TMPDIR/${rt}.um_re.sam -o $TMPDIR/${rt}.um_ns.sam
$REALLENDIR/extractLength.rb $TMPDIR/${rt}.um.sam -o  $TMPDIR/${rt}.um.len
$REALLENDIR/pairTest.rb $SVOPT -f $TMPDIR/${rt}.um.len -r $TMPDIR/${rt}.um_ns.sam -o $TMPDIR/${rt}.um_map.sam
$REALLENDIR/prepDDP.rb $TMPDIR/${rt}.um_map.sam -o $TMPDIR/${rt}.um_ddp.pos
$REALLENDIR/DDP.rb $REFHG19 -i $REFHG19I -f 0.7 $TMPDIR/${rt}.um.sam $TMPDIR/${rt}.um_ddp.pos -o $TMPDIR/${rt}.um_ddp.sam

$REALLENDIR/mergeSAM.rb $TMPDIR/${rt}.sc_map.sam $TMPDIR/${rt}.um_ddp.sam -o $TMPDIR/${rt}.merge.sam
$REALLENDIR/filtBP.rb -r 30 $TMPDIR/${rt}.merge.sam -b $TMPDIR/${rt}.bp.bed -o $TMPDIR/${rt}.bp.sam -p ${rt}.bp
