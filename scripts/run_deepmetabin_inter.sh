#!/bin/bash
set -xe
# reads=`pwd`/$1
# contigs=`pwd`/$3

reads=$1
contigs=$3
output=$4

if [ ! -d $output ]
then
    mkdir $output
    cd $output
else
	cd $output
fi

cp $contigs contigs.fasta

if [ ! -f "contigs.map.sorted.bam" ]
then
    bwa index contigs.fasta
    bwa mem contigs.fasta $reads -t 100 | samtools sort -@ 40 -n -o contigs.map.sorted.bam
fi

source ~/.bashrc; conda deactivate; conda activate deepmetabin
python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/preprocessing.py --outdir ./ --fasta contigs.fasta --bamfiles contigs.map.sorted.bam -m 1000
#python /datahome/datasets/ericteam/zmzhang/csmxrao/DeepMetaBin/mingxing/deepmetabin/run.py datamodule.zarr_dataset_path=data.zarr datamodule.output=deepmetabin_out model.contignames_path=contignames.npz model.contig_path=contigs.fasta
python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/train.py -data data.zarr --contignames_path contignames.npz --contig_path contigs.fasta --output deepmetabin_out

cd deepmetabin_out/results/pre_bins && checkm lineage_wf -t 100 -x fasta --tab_table -f checkm.tsv ./ ./
cd ../../../
python /datahome/datasets/ericteam/csgyguo/deepmetabin/secondary_clustering.py --primary_out deepmetabin_out --contigname_path contignames.npz --output_path deepmetabin_out/results --binned_length 1000 --contig_path contigs.fasta
