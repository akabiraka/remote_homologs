

#### Remote homology analysis
* Downloaded the family (fa) and superfamily (sf) representative SCOP sequences from here: http://scop.mrc-lmb.cam.ac.uk/download
    File: scop_fa_represeq_lib_latest.fa scop family representative domain sequences
        `wget http://scop.mrc-lmb.cam.ac.uk/files/scop_fa_represeq_lib_latest.fa -P data/downloaded/`
    File: scop_sf_represeq_lib_latest.fa scop superfamily representative domain sequences
        `wget http://scop.mrc-lmb.cam.ac.uk/files/scop_sf_represeq_lib_latest.fa -P data/downloaded/`
    File: scop-cla-latest.txt
        `wget http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt -P data/downloaded/`
    File: scop-des-latest.txt
        `wget http://scop.mrc-lmb.cam.ac.uk/files/scop-des-latest.txt -P data/downloaded/`

* CD-HIT clustering
    Install web-server can be found in the following two links: 
        https://sites.google.com/view/cd-hit/web-server?authuser=0 and https://github.com/weizhongli/cdhit-web-server
        sudo snap install docker
        docker run -d -h cdhit-server --name cdhit-server -p 8088:80 weizhongli1987/cdhit-server:latest (this requires super-user permission)
        http://localhost:8088/cdhit-web-server on the browser

    running CD-HIT:
        Sequence identity cut-off: .9 - .1
        And other parameters are per default, such as 
            -G: use global sequence identity: Yes
            -b bandwidth of alignment: 20
            -l length of sequence to skip: 10
        Then we downloaded the cluster file that is sorted by cluster size.
        wget http://localhost:8088/cdhit-web-server/output/1681921330/1681921330.fas.1.clstr.sorted -O data/cdhit_clusters_at_th_family/at_90_seq_identity.txt
        ... and so on (link will be different).


# Workflow
* To compute fasta seq dictionary and corresponding embeddings: `python generators/precompute_seq_dict_and_embeddings.py`
* To generate the remote homologs using exclusion of seq-identity-based clusters and family tag: `python generators/precompute_remote_homologs_dataset.py`


# Analysis
* To compute the query ranking performance (f1) at K (1, 5, 10) neighbors for each query: `python analyzers/compute_query_ranking_performance_using_knn.py`
