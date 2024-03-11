#!/bin/bash
set -xe
reads1=`pwd`/$1
reads2=`pwd`/$2
contigs=`pwd`/$3
output=$4

if [ ! -d $output  ]; then
	mkdir $output
fi
cd $output
if [ ! -f contigs.fasta ]; then
	cp $contigs contigs.fasta
fi
if [ ! -f contigs.fasta.amb ]; then
	bwa index contigs.fasta
fi
if [ ! -f align.bam ]; then
	bwa mem contigs.fasta $reads1 $reads2 -t 100 | samtools sort -@ 40 -n -o align.bam
fi
if [ ! -d vamb_out/bins ]; then
	/home/comp/zmzhang/software/bin/time -v vamb --outdir vamb_out --fasta contigs.fasta --bamfiles align.bam --minfasta 200000
	/home/comp/zmzhang/software/bin/time -v vamb --outdir vamb_out_2.5k --fasta contigs.fasta --bamfiles align.bam --minfasta 200000 -m 2500
fi
checkm lineage_wf vamb_out/bins -x fna -f vamb_out.checkm.out vamb_out.checkm -t 100
checkm lineage_wf vamb_out_2.5k/bins -x fna -f vamb_out_2.5k.checkm.out vamb_out_2.5k.checkm -t 100
