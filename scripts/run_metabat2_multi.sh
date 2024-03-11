#!/bin/bash

set -xe

trap clean EXIT SIGTERM
clean(){
    # 全部结束后解绑文件描述符并删除管道
    exec 4<&-
    exec 4>&-
    rm -f mylist
    kill -9 -$$
}

cd $(dirname $0)
bamfiles=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/binning_results/Airways/metabat2/bamfiles
samfiles=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/binning_results/Airways/metabat2/samfiles
reads_dir=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/reads/Airways/short_read
contig_dir=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/spades/Airways
combined_contig_fasta=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/binning_results/Airways/vamb/catalogue.fna.gz
[ ! -d $bamfiles ] && mkdir $bamfiles
[ ! -d $samfiles ] && mkdir $samfiles
reads_names=`ls $reads_dir`

source /home/comp/zmzhang/software/anaconda3/bin/activate SemiBin
thread_num=10

# make FIFO
mkfifo mylist
# bind file descriptor 4 to mylist
exec 4<>mylist
# write the thread number to the fifo
for ((i=0; i < $thread_num; i++)); do
    echo $i >&4
done

# for read_name in $reads_names; do
#     tmp=(${read_name//_/ })
#     cp $contig_dir/$read_name/contigs.fasta S${tmp[-1]}.fasta
# done

# fasta=`ls S*.fasta`
# SemiBin concatenate_fasta --input-fasta $fasta --output .

# bowtie2-build \
#     -f concatenated.fa concatenated.fa

# for read_name in $reads_names; do
#     read p_idx <&4
#     {
#         bowtie2 -q --fr \
#             -x concatenated.fa \
#             --interleaved $reads_dir/$read_name/reads/anonymous_reads.fq.gz \
#             -S $samfiles/$read_name.sam \
#             -p 150

#         samtools view -h -b -S $samfiles/$read_name.sam -o $bamfiles/$read_name.bam -@ 150

#         samtools view -b -F 4 $bamfiles/$read_name.bam -o $bamfiles/$read_name.mapped.bam -@ 150

#         samtools sort \
#             -m 1000000000 $bamfiles/$read_name.mapped.bam \
#             -o $bamfiles/$read_name.mapped.sorted.bam -@ 150
#         echo p_idx>&4
#     } &
# done
# wait

source /home/comp/zmzhang/software/anaconda3/bin/activate base
# bams=`find $bamfiles/*.sorted.bam`
# jgi_summarize_bam_contig_depths --outputDepth depth.txt $bams
python metabat2_split.py --abdfile depth.txt --out sample_dir
mkdir -p all_outputs
sample_names=`ls sample_dir`
for sample_name in $sample_names; do
    read p_idx <&4
    {
        mv $sample_name.fasta sample_dir/$sample_name
        cd sample_dir/$sample_name
        metabat2 -i $sample_name.fasta -a depth.txt -o bin -m 1500
        rm $sample_name.fasta
        for file in `ls bin.*.fa`; do mv $file $sample_name-$file; done
        cp *bin.*.fa $(dirname $bamfiles)/all_outputs
        echo p_idx>&4
    } &
done
wait
