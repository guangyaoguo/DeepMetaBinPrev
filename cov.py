import pandas as pd

def coverage_score(path):
    df = pd.read_csv(path, sep='\t', header=None, names=['contig', 'cluster'])

    df['length'] = df['contig'].str.extract(r'length_(\d+)_').astype(int)
    df['taxid'] = df['contig'].str.extract(r'kraken:taxid\|(\d+)')
    df = df[df['taxid'] != '0']

    max_length_per_taxid_cluster = df.groupby(['cluster', 'taxid'])['length'].sum().reset_index()
    total_length_per_taxid = df.groupby('taxid')['length'].sum().reset_index()
    

    merged_df = pd.merge(max_length_per_taxid_cluster, total_length_per_taxid, on='taxid', suffixes=('_max_cluster', '_total'))

    merged_df['length_ratio'] = merged_df['length_max_cluster'] / merged_df['length_total']

    idx = merged_df.groupby('cluster')['length_max_cluster'].idxmax()
    filtered_df = merged_df.loc[idx]
    
    cluster_avg_ratios = filtered_df['length_ratio'].mean()

    return filtered_df

if __name__ == '__main__':
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_csv>")
        sys.exit(1)

    path = sys.argv[1]
    result = coverage_score(path)
    print(result)