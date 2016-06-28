#! /bin/bash

function usage()
{
    echo "usage: reallen-fast [-p <NP>] <in.bam>"
    echo "  -p NUMBER   Number of processor"
    exit 0
}

# Check option
while getopts p:h option
do
    case ${option} in
	p)
	    ARG1=${OPTARG}
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
BWAHG19="$REALLENROOT/reference/bwa-hg19"
FASTFILT="/home/clc/test/bamfilter/bin/bamfilter"
TMPDIR=`pwd`/temporary_files_reallen
SAMTLS=samtools
BWA=bwa

if [ ! -e $TMPDIR ]
then
    mkdir $TMPDIR
fi

rt=`basename $1 .bam`

if [ $ARG1 ]
then
    PROCOPT="-t ${ARG4}"
fi

# Run ReALLEN & BWA
$FASTFILT -u -p -b 8 -s 8 ${rt}.bam --softclip $TMPDIR/${rt}.sc.fq -o $TMPDIR/${rt}.um.sam --fixed $TMPDIR/${rt}.um.fq -f 45 -l 100 -k
echo
echo "performing BWA for ${rt}.sc.fq"
$BWA mem ${PROCOPT} -O3 -E1 -T20 -a $BWAHG19 $TMPDIR/${rt}.sc.fq > $TMPDIR/${rt}.sc_re.sam
$REALLENDIR/filtNonSpecific.rb -s 10 -f 0.7 $TMPDIR/${rt}.sc_re.sam -o $TMPDIR/${rt}.sc_ns.sam
$REALLENDIR/pairTest.rb -l 40 -r $TMPDIR/${rt}.sc_ns.sam -o $TMPDIR/${rt}.sc_map.sam

echo
echo "performing BWA for ${rt}.um.fq"
$BWA mem ${PROCOPT} -O3 -E1 -T35 -a $BWAHG19 $TMPDIR/${rt}.um.fq > $TMPDIR/${rt}.um_re.sam
$REALLENDIR/filtNonSpecific.rb -s 15 -f 0.6 $TMPDIR/${rt}.um_re.sam -o $TMPDIR/${rt}.um_ns.sam
$REALLENDIR/extractLength.rb $TMPDIR/${rt}.um.sam -o  $TMPDIR/${rt}.um.len
$REALLENDIR/pairTest.rb -l 40 -f $TMPDIR/${rt}.um.len -r $TMPDIR/${rt}.um_ns.sam -o $TMPDIR/${rt}.um_map.sam
$REALLENDIR/prepDDP.rb $TMPDIR/${rt}.um_map.sam -o $TMPDIR/${rt}.um_ddp.pos
$REALLENDIR/DDP.rb $REFHG19 -i $REFHG19I $TMPDIR/${rt}.um.sam $TMPDIR/${rt}.um_ddp.pos -o $TMPDIR/${rt}.um_ddp.sam

$REALLENDIR/mergeSAM.rb $TMPDIR/${rt}.sc_map.sam $TMPDIR/${rt}.um_ddp.sam -o $TMPDIR/${rt}.merge.sam
$REALLENDIR/filtBP.rb -r 30 $TMPDIR/${rt}.merge.sam -b $TMPDIR/${rt}.bp.bed -o $TMPDIR/${rt}.bp.sam -p ${rt}.bp
