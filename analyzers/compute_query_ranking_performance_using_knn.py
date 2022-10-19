import sys
sys.path.append("../remote_homologs")

from sklearn.neighbors import NearestNeighbors
import numpy as np
import collections
from sklearn import metrics
import utils as Utils

def compute_f1_at_k(positive_label, num_of_positive_examples, ranked_neighbor_indices, K=1):
    # the ground-truth labels of K best neighbors
    ones = np.repeat([1], num_of_positive_examples)[:K] # first K ranked neighbors should be relevent as 1
    y_true = np.repeat([0], K)
    y_true[:ones.shape[0]] = ones # others should be non-relevent as 0
    # print(y_true)
    
    ranked_neighbor_indices_upto_k = ranked_neighbor_indices[:K]
    y_pred = []
    for label in sfam_labels[ranked_neighbor_indices_upto_k]:
        if positive_label==label:
            y_pred.append(1) # if the ranked neighbors have the same label as query label
        else:
            y_pred.append(0)
    y_pred = np.array(y_pred)
    # print(y_pred)

    f1 = metrics.f1_score(y_true, y_pred)
    # print(f1)
    return f1
        
        
# the ranking should put the homologs at the first than the later
def compute_query_ranking_performance_metric(seqs_rep, sfam_labels):
    sfam_freq_dict = collections.Counter(sfam_labels)

    num_of_queries = seqs_rep.shape[0]
    neigh = NearestNeighbors(metric="cosine")

    f1_scores = {"at1": [], "at5": [], "at10": []}
    for i in range(num_of_queries):
        a_query, sfam = seqs_rep[i], sfam_labels[i]
        
        # the dataset and labels of which against the query will be run
        X = np.concatenate((seqs_rep[:i, :], seqs_rep[i+1:, :])) 
        Y = sfam_labels[:i], sfam_labels[i+1:]
        # print(sfam)#a_query, , X)
        
        # computing the neighbors and distances of the query
        neigh.fit(X)
        k_neighbors = neigh.kneighbors([a_query], 10, return_distance=True) # it returns fixed 10 neighbors
        distances, ranked_neighbor_indices = k_neighbors[0][0], k_neighbors[1][0]
        # print(distances, ranked_neighbor_indices) # the idxs are sorted by the distance
        
        # compute f1 score for each query at K
        # K=1 #n_neighbors 1, 5, 10
        f1at1 = compute_f1_at_k(sfam, sfam_freq_dict[sfam], ranked_neighbor_indices, K=1)
        f1at5 = compute_f1_at_k(sfam, sfam_freq_dict[sfam], ranked_neighbor_indices, K=5)
        f1at10 = compute_f1_at_k(sfam, sfam_freq_dict[sfam], ranked_neighbor_indices, K=10)
        print(f"query no: {i}: f1at1={f1at1}, f1at5={f1at5}, f1at10={f1at10}")
        f1_scores["at1"].append(f1at1)
        f1_scores["at5"].append(f1at5)
        f1_scores["at10"].append(f1at10)
        
        # if i==1000: break
    print(np.mean(f1_scores["at1"]), np.mean(f1_scores["at5"]), np.mean(f1_scores["at10"]))
    return f1_scores
    
# seqs_rep = np.array([[0., 0., 0.], [0., .5, 0.], [1., 1., .5], [0., 0., 0.1], [0., 0.1, 0.]])
# sfam_labels = np.array(["0", "1", "1","0", "0"])
# compute_query_ranking_performance_metric(seqs_rep, sfam_labels)


def accumulate_seq_reps_and_sf_labels(seq_identity_th):
    remote_homolog_dict = Utils.load_pickle(f"data/generated/remote_homolog_dicts/at_{seq_identity_th}_sfam.pkl") #{sf1: {p1, p2}, sf2: {p3, p4}}
    saved_seq_embedding_path = "data/generated/seqs_aalevel_embeddings_by_esm/"

    X, Y = [], []
    for i, (sfam_id, prot_ids_set) in enumerate(remote_homolog_dict.items()):
        # print(i, sfam_id, len(prot_ids_set))
        for prot_id in prot_ids_set:
            Y.append(sfam_id)
            
            seq_rep = Utils.load_pickle(f"{saved_seq_embedding_path}/{prot_id}.pkl") # expects amino-acid level numpy array
            seq_rep = np.mean(seq_rep, axis=0)
            X.append(seq_rep)
            # break
        
        # break

    X, Y = np.array(X), np.array(Y)
    print(X.shape, Y.shape)
    # Utils.save_as_pickle(X, f"data/generated/esm_seqs_rep_with_sfam_at_seq_identity/{seq_identity_th}_seqs_rep.pkl")
    # Utils.save_as_pickle(Y, f"data/generated/esm_seqs_rep_with_sfam_at_seq_identity/{seq_identity_th}_sfam.pkl")
    return X, Y

seq_identity_th=10
seqs_rep, sfam_labels = accumulate_seq_reps_and_sf_labels(seq_identity_th)
f1_scores = compute_query_ranking_performance_metric(seqs_rep, sfam_labels)
Utils.save_as_pickle(f1_scores, f"outputs/remote_homology_f1_scores_for_each_query_at_k_neighbors_using_esm/seq_identity_th_{seq_identity_th}.pkl")