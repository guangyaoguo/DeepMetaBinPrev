import sys
import pandas as pd

# 检查是否有足够的命令行参数
if len(sys.argv) < 2:
    print("Usage: python script.py path_to_your_tsv_file.tsv")
    sys.exit(1)

# 读取TSV文件路径
file_path = sys.argv[1]

# 读取TSV文件
df = pd.read_csv(file_path, sep='\t')

# 筛选出完整性大于90%且污染度小于5%的行
filtered_df = df[(df['Completeness'] > 90) & (df['Contamination'] < 5)]

# 打印筛选后的DataFrame
print(filtered_df)

# 可选：将筛选后的结果保存到一个新的TSV文件中
# filtered_df.to_csv('filtered_clusters.tsv', sep='\t', index=False)
