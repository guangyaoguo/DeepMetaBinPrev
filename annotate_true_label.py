import sys
import re

# cmd
# python annotate_true_label.py /datahome/datasets/ericteam/zmzhang/csmxrao/DeepMetaBin/CAMI1/medium/true_labels_S001/labels.csv /datahome/datasets/ericteam/zmzhang/csmxrao/DeepMetaBin/CAMI1/medium/metaspades_S001/contigs.fasta /datahome/datasets/ericteam/zmzhang/Benchmarking/reads/CAMI/medium/metaspades_S001_done/contigs_t.fasta

# 获取命令行参数
csv_file = sys.argv[1]
old_fasta_file = sys.argv[2]
new_fasta_file = sys.argv[3]

def sort_nodes(filename):
    node_list = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            match = re.search(r'NODE_(\d+)', line)
            if match:
                node_number = int(match.group(1))
                node_list.append((node_number, line.strip()))
    
    sorted_nodes = [node[1] for node in sorted(node_list, key=lambda x: x[0])]
    return sorted_nodes

sorted_nodes = sort_nodes(csv_file)
print(sorted_nodes[:5])

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
        if i == len(sorted_nodes):
            new_fasta_lines.append(line.strip() + '_taxid|000000000' + '\n')
            continue

        match_fa = re.search(r'>NODE_(\d+)_', line)
        match_csv = re.search(r'NODE_(\d+)_', sorted_nodes[i])

        if match_fa.group(1) == match_csv.group(1):
            new_fasta_lines.append('>' + sorted_nodes[i].replace(',', '_taxid|') + '\n')
            i += 1
            
        else:
            new_fasta_lines.append(line.strip() + '_taxid|000000000' + '\n')

# 将新的FASTA行写入新的FASTA文件
with open(new_fasta_file, 'w') as file:
    file.writelines(new_fasta_lines)