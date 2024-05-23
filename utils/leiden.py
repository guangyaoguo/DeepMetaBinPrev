import leidenalg
from leidenalg import RBERVertexPartition, Optimiser
import scanpy as sc
from igraph import Graph
import hnswlib
import numpy as np
import pandas as pd
import functools
import time
import os
import logging
import multiprocessing
from sklearn.preprocessing import normalize
import scipy.sparse as sp
from sklearn.cluster._kmeans import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms, MiniBatchKMeans
from typing import List, Optional, Union
from utils.utils import gen_seed
logger = logging.getLogger('Leiden')
logger.setLevel(logging.INFO)

def fit_hnsw_index(logger, features,num_threads, ef: int = 100, M: int = 16,
                   space: str = 'l2', save_index_file: bool = False) -> hnswlib.Index:
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

    time_start = time.time()
    num_elements = len(features)
    labels_index = np.arange(num_elements)
    EMBEDDING_SIZE = len(features[0])

    # Declaring index
    # possible space options are l2, cosine or ip
    p = hnswlib.Index(space=space, dim=EMBEDDING_SIZE)

    # Initing index - the maximum number of elements should be known
    p.init_index(max_elements=num_elements, ef_construction=ef, M=M)

    # Element insertion
    int_labels = p.add_items(features, labels_index, num_threads=num_threads)

    # Controlling the recall by setting ef
    # ef should always be > k
    p.set_ef(ef)

    # If you want to save the graph to a file
    if save_index_file:
        p.save_index(save_index_file)
    time_end = time.time()
    logger.info('Time cost:\t' +str(time_end - time_start) + "s")
    return p

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

def gen_seed_idx(seedURL: str, contig_id_list: List[str]) -> List[int]:
    """
    Generate a list of indices corresponding to seed contig IDs from a given URL.

    :param seedURL: The URL or path to the file containing seed contig names.
    :param contig_id_list: List of all contig IDs to match with the seed contig names.
    :return: List[int]
    """
    seed_list = []
    with open(seedURL) as f:
        for line in f:
            if line.rstrip('\n') in contig_id_list:
                seed_list.append(line.rstrip('\n'))
    name_map = dict(zip(contig_id_list, range(len(contig_id_list))))
    seed_idx = [name_map[seed_name] for seed_name in seed_list]
    return seed_idx

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

def run_leiden(output_file: str, namelist: List[str],
               ann_neighbor_indices: np.ndarray, ann_distances: np.ndarray,
               length_weight: List[float], max_edges: int, norm_embeddings: np.ndarray,
               bandwidth: float = 0.1, lmode: str = 'l2', initial_list: Optional[List[Union[int, None]]] = None,
               is_membership_fixed: Optional[bool] = None, resolution_parameter: float = 1.0,
               partgraph_ratio: int = 50):
    """
    Run Leiden community detection algorithm and save the results.

    :param output_file: The path to the output file.
    :param namelist: A list of contig names.
    :param ann_neighbor_indices: Array of ANN neighbor indices.
    :param ann_distances: Array of ANN distances.
    :param length_weight: List of length weights.
    :param max_edges: Maximum number of edges.
    :param norm_embeddings: Array of normalized embeddings.
    :param bandwidth: Bandwidth parameter (default: 0.1).
    :param lmode: Distance mode ('l1' or 'l2', default: 'l2').
    :param initial_list: Initial membership list (default: None).
    :param is_membership_fixed: Whether membership is fixed (default: None).
    :param resolution_parameter: Resolution parameter (default: 1.0).
    :param partgraph_ratio: Partition graph ratio (default: 50).

    :return: None
    """

    sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    targets_indices = ann_neighbor_indices[:,1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:,1:]
    wei = wei.flatten()

    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    if lmode == 'l1':
        wei = np.sqrt(wei)
        wei = np.exp(-wei / bandwidth)

    if lmode == 'l2':
        wei = np.exp(-wei / bandwidth)

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    res = leidenalg.RBERVertexPartition(g,
                                        weights=wei, initial_membership = initial_list,
                                        resolution_parameter = resolution_parameter,node_sizes=length_weight)

    optimiser = leidenalg.Optimiser()
    optimiser.optimise_partition(res, is_membership_fixed=is_membership_fixed,n_iterations=-1)

    part = list(res)


    contig_labels_dict ={}
    # dict of communities
    numnode = 0
    rang = []
    for ci in range(len(part)):
        rang.append(ci)
        numnode = numnode+len(part[ci])
        for id in part[ci]:
            contig_labels_dict[namelist[id]] = 'group'+str(ci)

    logger.info(output_file)
    f = open(output_file, 'w')
    for contigIdx in range(len(contig_labels_dict)):
        f.write(namelist[contigIdx] + "\t" + str(contig_labels_dict[namelist[contigIdx]]) + "\n")
    f.close()

def leiden_tuning(namelist, contig_path, latent, output_path):
    import multiprocessing

    # Create a logger instance
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Create a formatter
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    console_handler.setFormatter(formatter)

    # Add the handler to the logger
    logger.addHandler(console_handler)

    # Log an information message
    logger.info('Logger has been set up for logging information.')

    num_workers = 240
    norm_embeddings = normalize(latent)
    contig_length_threshold = 1000
    marker_name = "bacar_marker"
    quarter="2quarter"

    seed_file = gen_seed(logger, contig_path, num_workers, contig_length_threshold,
             marker_name, quarter)
    ##### #########  hnswlib_method
    parameter_list = [1, 5,10,30,50,70, 90, 110]
    bandwidth_list = [0.05, 0.1,0.15, 0.2,0.3]
    partgraph_ratio_list =[50,100,80]
    max_edges_list = [100]
    for max_edges in max_edges_list:
        p = fit_hnsw_index(logger, norm_embeddings, num_workers, ef=max_edges * 10)
        seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
        initial_list = list(np.arange(len(namelist)))
        is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

        time_start = time.time()
        ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges+1, num_threads=num_workers)
        #ann_distances is l2 distance's square
        time_end = time.time()
        logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")

        with multiprocessing.Pool(num_workers) as multiprocess:
            for partgraph_ratio in partgraph_ratio_list:
                for bandwidth in bandwidth_list:
                    for para in parameter_list:
                        output_file = output_path + 'Leiden_bandwidth_' + str(
                            bandwidth) + '_res_maxedges' + str(max_edges) + 'respara_'+str(para)+'_partgraph_ratio_'+str(partgraph_ratio)+'.tsv'

                        if not (os.path.exists(output_file)):
                            multiprocess.apply_async(run_leiden, (output_file, namelist, ann_neighbor_indices, ann_distances, length_weight, max_edges, norm_embeddings,
                                                                        bandwidth, 'l2', initial_list,is_membership_fixed,
                                                                        para, partgraph_ratio))

            multiprocess.close()
            multiprocess.join()
        logger.info('multiprocess Done')
