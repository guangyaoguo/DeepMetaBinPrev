#!/bin/bash
set -xe

contigs=`pwd`/$1
reads1=`pwd`/$2
reads2=`pwd`/$3
out=$4


bwa index $contigs
bwa mem $contigs $reads1 $reads2 -t 100 | samtools sort -@ 40 -o $out.bam
samtools index -@ 50 $out.bam
