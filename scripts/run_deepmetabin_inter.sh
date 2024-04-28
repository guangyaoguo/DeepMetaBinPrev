#!/bin/bash
set -xe
# reads=`pwd`/$1
# contigs=`pwd`/$3

#cmd
#nohup scripts/run_deepmetabin_inter.sh /datahome/datasets/ericteam/zmzhang/Benchmarking/reads/CAMI/medium/RM2_S001__insert_270.fq.gz /datahome/datasets/ericteam/zmzhang/Benchmarking/reads/CAMI/medium/metaspades_S001_done/contigs.fasta deepmetabin_out

reads=$1
contigs=$2
output=$3

if [ ! -d $output ]
then
    mkdir $output
    cd $output
else
	cd $output
fi

source ~/.bashrc; conda deactivate; conda activate deepmetabin
python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/get_multiviews.py --n_views 6 --contig_file $contigs --out_augdata_path /datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out --contig_len 1000

# tips to use get_multiviews.py
# python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/get_multiviews.py --n_views 6 --contig_file /datahome/datasets/ericteam/zmzhang/Benchmarking/reads/CAMI/medium/metaspades_S002_done/contigs.fasta --out_augdata_path /datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out --contig_len 1000

cp $contigs contigs.fasta

if [ ! -f "contigs.map.sorted.bam" ]
then
    bwa index contigs.fasta
    bwa mem contigs.fasta $reads -t 100 | samtools sort -@ 40 -n -o contigs.map.sorted.bam
fi

python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/preprocessing.py --outdir ./ --fasta /datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out/contigs.fasta --bamfiles /datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out/contigs.map.sorted.bam -m 1000
#python /datahome/datasets/ericteam/zmzhang/csmxrao/DeepMetaBin/mingxing/deepmetabin/run.py datamodule.zarr_dataset_path=data.zarr datamodule.output=deepmetabin_out model.contignames_path=contignames.npz model.contig_path=contigs.fasta
python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/train.py -data data.zarr --contignames_path contignames.npz --contig_path /datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out/contigs.fasta --output deepmetabin_out

cd deepmetabin_out/results/pre_bins && checkm lineage_wf -t 100 -x fasta --tab_table -f checkm.tsv ./ ./
cd ../../../
python /datahome/datasets/ericteam/csgyguo/DeepMetaBin/secondary_clustering.py --primary_out deepmetabin_out --contigname_path contignames.npz --output_path deepmetabin_out/results --binned_length 1000 --contig_path /datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out/contigs.fasta