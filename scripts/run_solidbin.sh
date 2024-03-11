#!/bin/bash
set -xe
reads1=`pwd`/$1
reads2=`pwd`/$2
contigs=`pwd`/$3
output=$4

mkdir $output
cd $output
cp $contigs contigs.fasta

mkdir reads_dir
ln -s $reads1 $reads2 reads_dir/
bash /home/comp/zmzhang/mingxing/SolidBin/scripts/gen_cov.sh `pwd`
bash /home/comp/zmzhang/mingxing/SolidBin/scripts/run.sh contigs.fasta 1000 4
/home/comp/zmzhang/software/bin/time -v python /home/comp/zmzhang/mingxing/SolidBin/SolidBin.py --contig_file contigs_1000.fa --coverage_profiles coverage_f1000.tsv --composition_profiles kmer_4_f1000.csv --output result.tsv --log log.txt > time.log
