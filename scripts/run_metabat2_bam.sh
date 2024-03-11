#!/bin/bash
set -xe
bam=`pwd`/$1
contigs=`pwd`/$2
output=$3

mkdir $output
cd $output
cp $contigs contigs.fasta

jgi_summarize_bam_contig_depths --outputDepth depth.txt $bam
/home/comp/zmzhang/software/bin/time -v metabat2 -i contigs.fasta -a depth.txt -o bin -m 1500