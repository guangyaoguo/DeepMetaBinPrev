#!/bin/bash
set -xe
reads1=`pwd`/$1
reads2=`pwd`/$2
contigs=`pwd`/$3
output=$4

mkdir $output
cd $output
cp $contigs contigs.fasta

run_MaxBin.pl -contig contigs.fasta -out maxbin -reads $reads1 -reads $reads2 -min_contig_length 1000