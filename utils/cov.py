import pandas as pd

def coverage_score(path):
    df = pd.read_csv(path, sep='\t', header=None, names=['contig', 'cluster'])

    df['length'] = df['contig'].str.extract(r'length_(\d+)_').astype(int)
    df['taxid'] = df['contig'].str.extract(r'taxid\|(.*)')
    df = df[df['taxid'] != '000000000']

    # df['taxid'] = df['contig'].str.extract(r'kraken:taxid\|(\d+)')
    # df = df[df['taxid'] != '0']

    df2 = df.copy()
    max_length_per_taxid_cluster = df.groupby(['cluster', 'taxid'])['length'].sum().reset_index()
    total_length_per_taxid = df.groupby('taxid')['length'].sum().reset_index()
    total_length_per_cluster = df.groupby('cluster')['length'].sum().reset_index()
    max_length_per_taxid_cluster_refined = max_length_per_taxid_cluster.groupby('cluster')['length'].idxmax()
    max_length_per_taxid_cluster_filted = max_length_per_taxid_cluster.loc[max_length_per_taxid_cluster_refined]
    merged_df_precision = pd.merge(max_length_per_taxid_cluster_filted, total_length_per_cluster, on='cluster', suffixes=('_major', '_total'))
    merged_df_precision['precision'] = merged_df_precision['length_major']/ merged_df_precision['length_total']
    precision = merged_df_precision['precision'].mean()

    merged_df = pd.merge(max_length_per_taxid_cluster, total_length_per_taxid, on='taxid', suffixes=('_max_cluster', '_total'))

    merged_df['length_ratio'] = merged_df['length_max_cluster'] / merged_df['length_total']

    idx = merged_df.groupby('cluster')['length_max_cluster'].idxmax()
    filtered_df = merged_df.loc[idx]
    
    recall = filtered_df['length_ratio'].mean()



    return [precision, recall, 2*recall*precision/(recall+precision)]


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_csv>")
        sys.exit(1)

    path = sys.argv[1]
    
    result = coverage_score(path)

    print('Precision, Recall, F1 Score:')
    print(result)