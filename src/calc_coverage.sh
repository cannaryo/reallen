#! /bin/bash

if [ $# -le 2 ]
    then
    echo "Usage: calc_coverage.sh input1[bed] input2[BAM] output[csv]"
    echo "Input file1 is estimated as output of filtBP.rb" 
    echo "Input file2 is original BAM file" 
    exit 0
fi

if [ -e coverage.csv.tmp ]
    then
    echo "Cannot write coverage.csv.tmp"
    echo "File exist"
    exit 0
fi

samtools bedcov $1 $2 > coverage_tmp.bed
/home/clc/test/GS-SATv2/calc_coverage_sub.rb coverage_tmp.bed > $3
rm coverage_tmp.bed
