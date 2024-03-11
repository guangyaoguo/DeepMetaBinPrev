#!/bin/bash
contigs=`pwd`/$1
bins=`pwd`/$2
output=$3
graph=`pwd`/$4
path=`pwd`/$5

export PATH=/home/comp/zmzhang/code:$PATH

mkdir $output
cd $output

python /home/comp/zmzhang/software/GraphBin2/support/prepResult.py --binned $bins --output ./
graphbin2 --assembler spades --contigs --graph $graph --paths $path --binned initial_contig_bins.csv --output graphbin_output
cd graphbin_output
get_binning_results.py $contigs graphbin2_output.csv bins
cd bins
checkm lineage_wf -t 100 -x fasta --tab_table -f checkm.tsv ./ ./