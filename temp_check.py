import numpy as np

# 文件路径
file_path = '/datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out/deepmetabin_out/results/latent.npy'

# 加载.npy文件
data = np.load(file_path)

# 打印数据
print(data.shape)