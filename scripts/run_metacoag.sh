#!/bin/bash
set -xe
reads1=`pwd`/$1
reads2=`pwd`/$2
contigs=`pwd`/$3
output=$4
graph=`pwd`/$5
path=`pwd`/$6

export PATH=/home/comp/zmzhang/mingxing/MetaCoAG:$PATH
mkdir $output
cd $output
cp $contigs contigs.fasta
coverm contig -1 $reads1 -2 $reads2 -r contigs.fasta -o abundance.tsv -t 8
sed -i '1d' abundance.tsv
/home/comp/zmzhang/software/bin/time -v metacoag --assembler spades --graph $graph --contigs contigs.fasta --paths $path --abundance abundance.tsv --output output_folder