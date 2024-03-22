#!/bin/bash

# 自动获取当前目录的路径
current_dir=$(pwd)

# 创建一个大文件夹来存放所有软件的 checkm2 结果
mkdir -p "$current_dir/all_checkm2_results"

# 获取可用线程数量
# total_threads=$(nproc)*0.9
threads=240
# 定义软件名称和对应的bin路径
low/high
software_bins=(
    "deepmetabin:deepmetabin/postprocess/post_bins/*.fasta"
    "metabat2:metabat2/*.fa"
    "metadecoder:metadecoder/*.fasta"
    "vamb_1000:vamb_1000/vamb_out/bins/*.fasta"
    "metacoag:metacoag/output_folder/bins/*.fasta"
    "maxbin:maxbin/*.fasta"
)
# medium
# software_bins=(
#     "deepmetabin_001:deepmetabin_001/postprocess/post_bins/*.fasta"
#     "deepmetabin_002:deepmetabin_002/postprocess/post_bins/*.fasta"
#     "metabat2_001:metabat2_001/*.fa"
#     "metabat2_002:metabat2_002/*.fa"
#     "metadecoder_001:metadecoder_001/*.fasta"
#     "metadecoder_002:metadecoder_002/*.fasta"
#     "vamb_1000_S001:vamb_1000_S001/vamb_out/bins/*.fasta"
#     "vamb_1000_S002:vamb_1000_S002/vamb_out/bins/*.fasta"
#     "metacoag_001:metacoag_001/output_folder/bins/*.fasta"
#     "metacoag_002:metacoag_002/output_folder/bins/*.fasta"
#     "maxbin_S001:maxbin_S001/*.fasta"
#     "maxbin_S002:maxbin_S002/*.fasta"
# )
# 循环处理每个软件的bin路径
for software_bin in "${software_bins[@]}"; do
    software_name=$(echo $software_bin | cut -d':' -f1)
    bin_path="$current_dir/$(echo $software_bin | cut -d':' -f2)"
    fasta_files=$(ls $bin_path)
    # 执行 checkm2 命令来预测
    checkm2 predict --threads $threads --input $fasta_files --output-directory "$current_dir/all_checkm2_results/$software_name"
done

# 创建结果表格
echo -e "Software\tBin_Count_Completeness_GT_90_Contamination_LT_5" > "$current_dir/all_checkm2_results/comparison_table.tsv"

# 循环处理每个软件的结果
for software_bin in "$current_dir/all_checkm2_results"/*; do

    if [[ $software_bin == $current_dir'/all_checkm2_results/comparison_table.tsv' ]]; then
        continue
    fi
    software_name=$(basename $software_bin)
    bin_count=$(awk -F'\t' 'NR>1 && $2 > 90 && $3 < 5 {count++} END {print count}' "$software_bin/quality_report.tsv")
    echo -e "$software_name\t$bin_count" >> "$current_dir/all_checkm2_results/comparison_table.tsv"
done
