

#### Arbitrary
* Downloaded files description can be found here: http://scop.mrc-lmb.cam.ac.uk/download

# Workflow
* To compute fasta seq dictionary and corresponding embeddings: `python generators/precompute_seq_dict_and_embeddings.py`
* To generate the remote homologs using exclusion of seq-identity-based clusters and family tag: `python generators/precompute_remote_homologs_dataset.py`


# Analysis
* To compute the query ranking performance (f1) at K (1, 5, 10) neighbors for each query: `python analyzers/compute_query_ranking_performance_using_knn.py`
