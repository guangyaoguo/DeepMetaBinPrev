#!/bin/bash

set -ex

trap clean EXIT SIGTERM
clean(){
    # 全部结束后解绑文件描述符并删除管道
    exec 4<&-
    exec 4>&-
    rm -f mylist
    kill -9 -$$
}

cd $(dirname $0)
reads_dir=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/reads/Airways/short_read
contig_dir=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/spades/Airways
bams_dir=/datahome/datasets/ericteam/csmxrao/DeepMetaBin/CAMI2/binning_results/vamb1/bamfiles
[ ! -d $bams_dir ] && mkdir $bams_dir
reads_names=`ls $reads_dir`




thread_num=10


main() {
    # 创建一个管道
    mkfifo mylist
    # 给管道绑定文件描述符4
    exec 4<>mylist
    # 事先往管道中写入数据，要开启几个子进程就写入多少条数据
    for ((i=0; i < $thread_num; i++)); do
        echo $i >&4
    done


    contigs=`ls $contig_dir/*/contigs.fasta`
    concatenate.py $(dirname $bams_dir)/catalogue.fna.gz $contigs -m 1000
    bwa index $(dirname $bams_dir)/catalogue.fna.gz
    for read_name in $reads_names; do
        read p_idx <&4
        # 这里的 & 会开启一个子进程执行
        {
            # [ ! -d $bams_dir/read_name ] && mkdir $bams_dir/read_name
            call $reads_dir/$read_name/reads/anonymous_reads.fq.gz $(dirname $bams_dir)/catalogue.fna.gz $bams_dir/$read_name.bam
            echo p_idx>&4
        } &
    done
    # 使用 wait 命令阻塞当前进程，直到所有子进程全部执行完
    wait
    /home/comp/zmzhang/software/bin/time -v vamb --outdir $(dirname $bams_dir)/vamb_out --fasta $(dirname $bams_dir)/catalogue.fna.gz --bamfiles $bams_dir/*.bam --minfasta 200000 -m 1000 -o C
}

call() {
    bwa mem -p $2 $1 -t 100 | samtools sort -@ 40 -n -o $3
}

main
