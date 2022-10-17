import sys
sys.path.append("../remote_homologs")

import os.path as osp
import utils as Utils

def get_cluster_dict_from_cdhitclusters(th=90):
    # if the cluster dict is already computed for usage, return that
    clus_dict_path = f"data/tmp/scop_sf_seq_cluster_dict_at_{th}_seq_identity.pkl"
    if osp.exists(clus_dict_path):
        cluster_dict = Utils.load_pickle(clus_dict_path)
        print(f"    num of clusters: {len(cluster_dict)}") #18575 at 40% identity
        return cluster_dict
    
    #else
    print("    Computing cluster dictionary.")
    
    f = open(f"data/generated/cdhit_clusters/scop_sf_seq_cluster_at_{th}_seq_identity.txt", mode="r")
    # out = open(f"data/generated/cdhit_clusters/scop_sf_seq_cluster_dict_at_{th}_seq_identity.txt", mode="w")
    flag=False
    cluster_dict = dict()
    cluster_no = 0
    for i, line in enumerate(f.readlines()):
        # print(i, line)
        
        if ">Cluster" not in line:
            line_items = line.split(" ")
            cluster.append(str(line_items[1][1:-3]))
            # print(cluster)

        elif ">Cluster" in line:
            if flag: 
                # out.writelines(" ".join(cluster)) #if want to see the clusters, comment out this line
                # out.write("\n")
                
                cluster_dict[f"c{cluster_no}"] = cluster
                cluster_no += 1
            
            cluster = []
            flag = True
            # print(" ".join(cluster))
            
    print(f"    num of clusters: {len(cluster_dict)}") #18575 at 40% identity
    Utils.save_as_pickle(cluster_dict, clus_dict_path)
    # print(cluster_dict)
    return cluster_dict

# cluster_dict = get_cluster_dict_from_cdhitclusters(90)
# print(cluster_dict["c0"]) # printing first cluster