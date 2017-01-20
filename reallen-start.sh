#! /bin/bash

function usage()
{
    echo "usage: reallen-start.sh [-SF][-p <NP>][-s <size>][-c <FILE>] <in.bam>"
    echo "  -S          input file is SAM format instead of BAM format"
    echo "  -F          rough filtering option (a little fast)"
    echo "  -p NUMBER   Number of processor [1]"
    echo "  -c FILE     specify config file [reallen.config]"
    echo "  -s SIZE     Minimum SV size [20]"
    echo "  -h          Show help message and exit"
    exit 0
}

# Check option
while getopts FSp:c:s:h option
do
    case ${option} in
	S)
	    ARG2=TRUE
	    ;;
	F)
	    ARG3=TRUE
	    ;;
	p)
	    ARG4=${OPTARG}
	    ;;
	c)
	    ARG5=${OPTARG}
	    ;;
	s)
	    ARG6=${OPTARG}
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

if [ -z $ARG5 ] ; then
    CONFIGFILE="$SCRIPTPATH/reallen.config"
else
    CONFIGFILE=$ARG5
fi
source $CONFIGFILE

if [ ! -e $TMPDIR ] ; then
    mkdir $TMPDIR
fi

if [ $ARG2 ] ; then
    rt=`basename $1 .sam`
    FIRSTDATA=${rt}.sam
else
    rt=`basename $1 .bam`
    $SAMTLS view -h ${rt}.bam > $TMPDIR/${rt}.sam 
    FIRSTDATA=$TMPDIR/${rt}.sam
fi

if [ $ARG3 ] ; then
    FILTOPT="--indel 10 -s 20"
else
    FILTOPT="--indel 8 -s 8 -u"
fi

if [ $ARG4 ] ; then
    PROCOPT="-t ${ARG4} "
fi

if [ $ARG6 ]
then
    SVOPT="-l $ARG6"
else
    SVOPT="-l 20"
fi


# Run ReALLEN & BWA
$REALLENDIR/filtRead.rb $FILTOPT $FIRSTDATA -o $TMPDIR/${rt}.filt.sam
$REALLENDIR/splitSoftClip.rb $TMPDIR/${rt}.filt.sam -e $TMPDIR/${rt}.um.sam -o $TMPDIR/${rt}.sc.fq
echo
echo "performing BWA for ${rt}.sc.fq"
$BWA mem ${PROCOPT} -O3 -E1 -T20 -a $BWAHG19  $TMPDIR/${rt}.sc.fq > $TMPDIR/${rt}.sc_re.sam
$REALLENDIR/filtNonSpecific.rb -s 10 -f 0.7 $TMPDIR/${rt}.sc_re.sam -o $TMPDIR/${rt}.sc_ns.sam
$REALLENDIR/pairTest.rb $SVOPT -r $TMPDIR/${rt}.sc_ns.sam -o $TMPDIR/${rt}.sc_map.sam

$REALLENDIR/splitRead.rb -l 100 -s 45 $TMPDIR/${rt}.um.sam -o $TMPDIR/${rt}.um.fq
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
