import sys
sys.path.append("../remote_homologs")

from sklearn.neighbors import NearestNeighbors
import numpy as np
import collections
from sklearn import metrics
import utils as Utils


seqs_rep = np.array([[0., 0., 0.], [0., .5, 0.], [1., 1., .5], [0., 0., 0.1], [0., 0.1, 0.]])
sfam_labels = np.array(["0", "1", "1","0", "0"])

sfam_freq_dict = collections.Counter(sfam_labels)

num_of_queries = seqs_rep.shape[0]
neigh = NearestNeighbors(metric="cosine")

f1_scores = []
for i in range(num_of_queries):
    a_query, sfam = seqs_rep[i], sfam_labels[i]
    
    # the dataset and labels of which against the query will be run
    X = np.concatenate((seqs_rep[:i, :], seqs_rep[i+1:, :])) 
    Y = sfam_labels[:i], sfam_labels[i+1:]
    # print(sfam)#a_query, , X)
    
    # computing the neighbors and distances of the query
    neigh.fit(X)
    k_neighbors = neigh.kneighbors([a_query], 10, return_distance=True) # it returns fixed 10 neighbors
    distances, idxs = k_neighbors[0][0], k_neighbors[1][0]
    # print(distances, idxs) # the idxs are sorted by the distance
    
    
    # compute f1 score for each query at K
    K=1 #n_neighbors 1, 5, 10
    
    # the ground-truth labels of K best neighbors
    ones = np.repeat([1], sfam_freq_dict[sfam])[:K]
    y_true = np.repeat([0], K)
    y_true[:ones.shape[0]] = ones
    # print(y_true)
    
    y_pred = []
    for label in sfam_labels[idxs]:
        if sfam==label:
            y_pred.append(1) # if the ranked neighbors have the same label as query label
        else:
            y_pred.append(0)
    y_pred = np.array(y_pred)
    print(y_pred)

    f1 = metrics.f1_score(y_true, y_pred)
    print(f1)
    f1_scores.append(f1)

print(np.mean(f1_scores), np.std(f1_scores))


seq_identity_th = 10
remote_homolog_dict = Utils.load_pickle(f"data/generated/remote_homolog_dicts/at_{seq_identity_th}_sfam.pkl") #{sf1: {p1, p2}, sf2: {p3, p4}}
saved_seq_embedding_path = "data/generated/seqs_aalevel_embeddings_by_esm/"

X, Y = [], []
for i, (sfam_id, prot_ids_set) in enumerate(remote_homolog_dict.items()):
    print(i, sfam_id, len(prot_ids_set))
    for prot_id in prot_ids_set:
        Y.append(sfam_id)
        
        seq_rep = Utils.load_pickle(f"{saved_seq_embedding_path}/{prot_id}.pkl") # expects amino-acid level numpy array
        seq_rep = np.mean(seq_rep, axis=0)
        X.append(seq_rep)
        # break
    
    # break

X, Y = np.array(X), np.array(Y)
print(X.shape, Y.shape)
Utils.save_as_pickle(X, f"data/generated/esm_seqs_rep_with_sfam/{seq_identity_th}_seqs_rep.pkl")
Utils.save_as_pickle(Y, f"data/generated/esm_seqs_rep_with_sfam/{seq_identity_th}_sfam.pkl")