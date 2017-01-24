#! /bin/bash

echo "Thank you for downloading ReALLEN :-)"
echo "This script checks dependent binaries and resources"
echo

CDIR="$0"
if [ -L $CDIR ] ; then
    CDIR=$(readlink $CDIR)
fi
CDIR=$(cd $(dirname $CDIR); pwd)

echo "Install Directory is ${CDIR}"

echo
echo "Searching for BWA:"
if [ -x "$(which bwa)" ]
then
    echo "OK"
else
    echo "BWA is not found!"
    echo "Please confirm the path to BWA in reallen.config"
    echo
fi
echo
echo "Searching for Samtools:"
if [ -x "$(which samtools)" ]
then
    echo "OK"
else
    echo "Samtools is not found!"
    echo "Please confirm the path to Samtools in reallen.config"
    echo
fi

REMOTEURL="mars.kmc.gr.jp/~kanno/download/reallen"
echo
echo "Remote file path: $REMOTEURL/hg19/reference/"
echo "Start downloading reference and resource files from remote host (Y/n)?"
read ANSWER
case $ANSWER in
    "" | "Y" | "y" | "yes" | "Yes" | "YES" )
	if [ ! -e "reference" ] ; then
	    mkdir reference
	fi
	cd reference
	wget $REMOTEURL/hg19/reference/hg19.fa
	wget $REMOTEURL/hg19/reference/hg19.fa.fai
	wget $REMOTEURL/hg19/reference/bwa-hg19.amb
	wget $REMOTEURL/hg19/reference/bwa-hg19.ann
	wget $REMOTEURL/hg19/reference/bwa-hg19.bwt
	wget $REMOTEURL/hg19/reference/bwa-hg19.pac
	wget $REMOTEURL/hg19/reference/bwa-hg19.sa
	cd ..
	if [ ! -e "resource" ] ; then
	    mkdir resource
	fi
	cd resource
	wget $REMOTEURL/hg19/resource/chr_name.txt
	wget $REMOTEURL/hg19/resource/fusion_cosmic.bed
	wget $REMOTEURL/hg19/resource/cytoband.csv
	wget $REMOTEURL/hg19/resource/hg19_gene_region.csv
	wget $REMOTEURL/hg19/resource/hg19_exon_region.csv
	cd ..
	if [ ! -e "testdata" ] ; then
	    mkdir testdata
	fi
	cd testdata
	wget $REMOTEURL/testdata/sample.bam
	wget $REMOTEURL/testdata/sample.bam.bai
	cd ..	
	;;
    * ) echo "Skip downloadging";;
esac

echo
echo "Setup finished."
