import sys
import re

# 获取命令行参数
csv_file = sys.argv[1]
old_fasta_file = sys.argv[2]
new_fasta_file = sys.argv[3]

# 读取CSV文件中的内容
with open(csv_file, 'r') as file:
    csv_lines = file.readlines()

# 读取旧的FASTA文件
with open(old_fasta_file, 'r') as file:
    old_fasta_lines = file.readlines()

# 创建一个新的FASTA文件
new_fasta_lines = []

# 遍历CSV文件中的每一行
i = 0
for line in old_fasta_lines:
    if not line.startswith('>NODE_'):
        new_fasta_lines.append(line)
    else:
        if i == len(csv_lines):
            new_fasta_lines.append(line.strip() + '_kraken:taxid|0\n')
            continue

        match_fa = re.search(r'>NODE_(\d+)_', line)
        match_csv = re.search(r'>NODE_(\d+)_', csv_lines[i])

        if match_fa.group(1) == match_csv.group(1):
            new_fasta_lines.append(csv_lines[i].replace(' ', '_'))
            i += 1
        else:
            new_fasta_lines.append(line.strip() + '_kraken:taxid|0\n')
                
# 将新的FASTA行写入新的FASTA文件
with open(new_fasta_file, 'w') as file:
    file.writelines(new_fasta_lines)