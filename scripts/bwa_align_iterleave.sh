#!/bin/bash
set -xe

contigs=`pwd`/$1
reads=`pwd`/$2

if [ ! -f $contigs.amb ]
then
    bwa index $contigs
fi
bwa mem $contigs $reads -t 100 -p | samtools sort -@ 40 -n -o contigs.map.sorted.bam