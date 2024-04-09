import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py <first_table_path> <second_table_path>")
    sys.exit(1)

first_table_path = sys.argv[1]
second_table_path = sys.argv[2]

df1 = pd.read_csv(first_table_path, sep='\t', index_col=0)

df2 = pd.read_csv(second_table_path, sep='\t', index_col=0)


max_values1 = df1.max(axis=1)

max_values2 = df2.max(axis=1)

merged_max_values = pd.concat([max_values1, max_values2], axis=1)
merged_max_values.columns = ['prev_max', 'new_max']

merged_max_values.to_csv('/datahome/datasets/ericteam/csgyguo/DeepMetaBin/deepmetabin_out/compare_gf.csv')