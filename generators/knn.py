import sys
sys.path.append("../remote_homologs")

from sklearn.neighbors import NearestNeighbors
import numpy as np
import collections
from sklearn import metrics
import utils as Utils

# K=1 #n_neighbors
# seqs_rep = np.array([[0., 0., 0.], [0., .5, 0.], [1., 1., .5], [0., 0., 0.1], [0., 0.1, 0.]])
# sfam_labels = np.array(["0", "1", "1","0", "0"])

# sfam_freq_dict = collections.Counter(sfam_labels)

# num_of_queries = seqs_rep.shape[0]
# neigh = NearestNeighbors(metric="cosine")

# f1_scores = []
# for i in range(num_of_queries):
#     a_query, sfam = seqs_rep[i], sfam_labels[i]
    
#     # the dataset and labels of which against the query will be run
#     X = np.concatenate((seqs_rep[:i, :], seqs_rep[i+1:, :])) 
#     Y = sfam_labels[:i], sfam_labels[i+1:]
#     # print(sfam)#a_query, , X)
    
#     # computing the neighbors and distances of the query
#     neigh.fit(X)
#     k_neighbors = neigh.kneighbors([a_query], K, return_distance=True)
#     distances, idxs = k_neighbors[0][0], k_neighbors[1][0]
#     # print(distances, idxs)
    
#     # the ground-truth labels of K neighbors
#     ones = np.repeat([1], sfam_freq_dict[sfam])[:K]
#     y_true = np.repeat([0], K)
#     y_true[:ones.shape[0]] = ones
#     # print(y_true)
    
#     y_pred = []
#     for label in sfam_labels[idxs]:
#         if sfam==label:
#             y_pred.append(1)
#         else:
#             y_pred.append(0)
#     y_pred = np.array(y_pred)
#     print(y_pred)

#     f1 = metrics.f1_score(y_true, y_pred)
#     print(f1)
#     f1_scores.append(f1)

# print(np.mean(f1_scores), np.std(f1_scores))

import esm
import torch

device = "cpu"#"cuda" if torch.cuda.is_available() else "cpu" # "cpu"#
esm1b, alphabet = esm.pretrained.esm1_t12_85M_UR50S() # model
esm1b = esm1b.to(device)
esm1b_batch_converter = alphabet.get_batch_converter()

def get_seq_rep(id, seq):
    uniprotid, batch_strs, seq_tokens = esm1b_batch_converter([(id, seq)])

    with torch.no_grad():
        results = esm1b(seq_tokens, repr_layers=[12], return_contacts=False)
    seq_rep = results["representations"][12] #1, max_seq_len, esmb_embed_dim
    seq_rep.squeeze_(0)
    seq_rep = seq_rep[1 : len(seq) + 1].mean(0) # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    return seq_rep.detach().cpu().numpy()


from dict_computation import get_fasta_seq_dict
fasta_seq_dict = get_fasta_seq_dict()

seq_identity_th = 10
remote_homolog_dict = Utils.load_pickle(f"data/generated/remote_homolog_dicts/at_{seq_identity_th}_sfam.pkl") #{sf1: {p1, p2}, sf2: {p3, p4}}
X, Y = [], []
for i, (sfam_id, prot_ids_set) in enumerate(remote_homolog_dict.items()):
    print(i, sfam_id, len(prot_ids_set))
    for prot_id in prot_ids_set:
        Y.append(sfam_id)
        
        seq = fasta_seq_dict[prot_id]["seq"]
        seq_rep = get_seq_rep(prot_id, seq)
        X.append(seq_rep)
        # break
    
    # break

X, Y = np.array(X), np.array(Y)
print(X.shape, Y.shape)
Utils.save_as_pickle(X, f"data/generated/esm_seqs_rep_with_sfam/{seq_identity_th}_seqs_rep.pkl")
Utils.save_as_pickle(Y, f"data/generated/esm_seqs_rep_with_sfam/{seq_identity_th}_sfam.pkl")