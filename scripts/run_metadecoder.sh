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
bwa mem contigs.fasta $reads1 $reads2 -t 100 | samtools sort -@ 40 -o align.sam
metadecoder coverage -s align.sam -o METADECODER.COVERAGE
metadecoder seed --threads 50 -f contigs.fasta -o METADECODER.SEED
/home/comp/zmzhang/software/bin/time -v metadecoder cluster -f contigs.fasta -c METADECODER.COVERAGE -s METADECODER.SEED -o METADECODER --min_sequence_length 2000
