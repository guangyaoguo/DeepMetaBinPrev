from leidenalg import RBERVertexPartition, Optimiser
import scanpy as sc
from igraph import Graph
import hnswlib
import numpy as np
import pandas as pd
import os
import logging
import multiprocessing
from sklearn.preprocessing import normalize

logger = logging.getLogger('Leiden')
logger.setLevel(logging.INFO)

def hnsw_index(features,num_threads, ef: int = 100, M: int = 16,
                   space: str = 'l2', save_index_file: bool = False):
    """
    Fit an HNSW index with the given features using the HNSWlib library; Convenience function to create HNSW graph.

    :param logger: The logger object for logging messages.
    :param features: A list of lists containing the embeddings.
    :param ef: The ef parameter to tune the HNSW algorithm (default: 100).
    :param M: The M parameter to tune the HNSW algorithm (default: 16).
    :param space: The space in which the index operates (default: 'l2').
    :param save_index_file: The path to save the HNSW index file (optional).

    :return: The HNSW index created using the given features.

    This function fits an HNSW index to the provided features, allowing efficient similarity search in high-dimensional spaces.
    """
        
    dim = features.shape[1]
    num_elements = len(features)
    label_id = np.arange(num_elements)

    p = hnswlib.Index(space=space, dim=dim)
    p.init_index(max_elements=num_elements, ef_construction=ef, M=M)
    int_labels = p.add_items(features, label_id, num_threads=num_threads)
    p.set_ef(ef)
    return p

def leiden_clustering_scanpy(anndata, resolution=1.0, random_state = 0, n_iterations = -1):
    anndata_neibour = sc.pp.neighbors(anndata, n_neighbors=5, n_pcs=None, use_rep=None, knn=True, method='gauss', transformer=None, metric='l2', random_state=0, key_added=None, copy=True)
    anndata_clustered = sc.tl.leiden(anndata_neibour, resolution=resolution, restrict_to=None, random_state=random_state, key_added='leiden', adjacency=None, directed=None, use_weights=True, n_iterations=n_iterations, partition_type=None, neighbors_key=None, obsp=None, copy=True, flavor='leidenalg')
    return anndata_clustered
    
def leiden_clustering_alg(namelist, ann_neighbor_indices, ann_distances, length_weight, max_edges, latent, bandwidth: float=0.1, lmode='l2', initial_list=None, is_membership_fixed=None, resolution_parameter=1.0, partgraph_ratio=50):
    sources = np.repeat(np.arange(len(latent)), max_edges)
    targets_indices = ann_neighbor_indices[:, 1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:, 1:]
    wei = wei.flatten()
    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    # if lmode == 'l1':
    #     wei = np.sqrt(wei)
    #     wei = np.exp(-wei / bandwidth)

    # if lmode == 'l2':
    #     wei = np.exp(-wei / bandwidth)

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]

    vcount = len(latent)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    res = RBERVertexPartition(g, weights=wei, initial_membership=None, resolution_parameter=resolution_parameter, node_sizes=length_weight)

    optimiser = Optimiser()
    optimiser.optimise_partition(res, is_membership_fixed=is_membership_fixed, n_iterations=-1)

    part = list(res)

    # contig_labels_dict = {}
    # numnode = 0
    # for ci in range(len(part)):
    #     numnode = numnode + len(part[ci])
    #     for id in part[ci]:
    #         contig_labels_dict[namelist[id]] = 'group' + str(ci)
    
    # return contig_labels_dict

    contig_labels_list = []  # 存储标签的列表
    for ci in range(len(part)):
        for id in part[ci]:
            contig_labels_list.append(str(ci))

    return contig_labels_list


def get_length_weight(contignames):
    length_list = np.array([int(np.char.find(name, 'length_') + 7) for name in contignames])
    return length_list

def cluster(latent, contignames, threads, max_edges = 100, prefix=None):

    p = hnsw_index(latent, threads, ef=max_edges * 10)

    ann_neighbor_indices, ann_distances = p.knn_query(latent, max_edges + 1, num_threads=threads)
    length_weight = get_length_weight(contignames)

    norm_embeddings = normalize(latent)
    
    return leiden_clustering_alg(contignames, ann_neighbor_indices, ann_distances, length_weight, max_edges, norm_embeddings,
                                                                        bandwidth = 0.1, lmode = 'l2', partgraph_ratio = 50)