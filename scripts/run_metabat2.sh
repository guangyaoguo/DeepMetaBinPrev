#!/bin/bash
set -xe
reads1=`pwd`/$1
reads2=`pwd`/$2
contigs=`pwd`/$3
output=$4

mkdir $output
cd $output
cp $contigs contigs.fasta

bwa index contigs.fasta
bwa mem -t 100 contigs.fasta $reads1 $reads2 | samtools sort -@ 100 -o align.bam
samtools index -@ 100 align.bam

jgi_summarize_bam_contig_depths --outputDepth depth.txt align.bam
/home/comp/zmzhang/software/bin/time -v metabat2 -i contigs.fasta -a depth.txt -o bin -m 1500