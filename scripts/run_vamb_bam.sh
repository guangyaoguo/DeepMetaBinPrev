#!/bin/bash
set -xe
bam=`pwd`/$1
contigs=`pwd`/$2
output=$3

mkdir $output
cd $output
cp $contigs contigs.fastas
samtools sort $bam -@ 40 -n -o align.bam
/home/comp/zmzhang/software/bin/time -v vamb --outdir vamb_out --fasta contigs.fasta --bamfiles align.bam -m 1000
