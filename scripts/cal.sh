#!/bin/bash

# 检查是否提供了文件路径参数
if [ $# -ne 1 ]; then
    echo "请提供表格文件的路径作为参数"
    exit 1
fi

# 初始化计数器
count=0

# 读取表格中的数据并逐行处理
while IFS=$'\t' read -r name completeness contamination _ _ _ _ _ _ _ _ _ _ _
do
    # 检查 completeness 大于 90% 且 contamination 小于 5%
    if [ "$(echo "$completeness > 90" | bc -l)" -eq 1 ] && [ "$(echo "$contamination < 5" | bc -l)" -eq 1 ]; then
        count=$((count+1))
    fi
done < "$1"

# 输出符合条件的数量
echo "Completeness 大于 90% 且 Contamination 小于 5% 的数量为: $count"