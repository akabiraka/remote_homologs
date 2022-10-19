import sys
sys.path.append("../remote_homologs")

import os.path as osp
import random
import utils as Utils
from dict_computation import get_fam_and_sfam_dict
from get_cluster_dict_from_cdhitclusters import get_cluster_dict_from_cdhitclusters


def exclude_SF_members(sfam_domain_id_set, domain_id_set, level="Family"):
    x = sfam_domain_id_set.intersection(domain_id_set)
    if len(x) > 1:                  
        print(f"    excluding {len(x)-1} SF members that belong to same {level}...")                    
        sfam_domain_id_set = sfam_domain_id_set - x
        sfam_domain_id_set = sfam_domain_id_set | set(random.sample(list(x), 1))
        # if level=="Cluster": # debug: cluster level redundency check
        #     print("cluster")
    return sfam_domain_id_set


def compute_remote_homologs_for_each_superfamily(seq_identity_th=90):
    out = f"data/generated/remote_homolog_dicts/at_{seq_identity_th}_sfam.pkl"
    if osp.exists(out):
        return Utils.load_pickle(out)
    
    cluster_dict = get_cluster_dict_from_cdhitclusters(seq_identity_th)
    fam_dict, sfam_dict = get_fam_and_sfam_dict()
    n_sfams = len(sfam_dict)

    new_sfam_dict = dict()
    for i, (sf, sfam_domain_id_set) in enumerate(sfam_dict.items()):
        n_mem_b4_exlusion = len(sfam_domain_id_set)

        # for fam, fam_domain_id_set in fam_dict.items():
        #     sfam_domain_id_set = exclude_SF_members(sfam_domain_id_set, fam_domain_id_set, level="Family")

        for clus, clus_domain_id_set in cluster_dict.items():
            sfam_domain_id_set = exclude_SF_members(sfam_domain_id_set, clus_domain_id_set, level="Cluster")
            
        
        print(f"{i}/{n_sfams}: {n_mem_b4_exlusion} and {len(sfam_domain_id_set)} SFam members b4 and after exclusion")
        print()

        if len(sfam_domain_id_set) >= 3: # excluding superfamilies that has no remote homologs.
            new_sfam_dict[sf] = sfam_domain_id_set
        # if i==1000: break

    new_sfam_dict = dict(sorted(new_sfam_dict.items(), key=lambda item: len(item[1]), reverse=True)) # sort the dict by the number of proteins in it by descending order
    Utils.save_as_pickle(new_sfam_dict, out)
    # print(new_sfam_dict)
    return new_sfam_dict

new_sfam_dict = compute_remote_homologs_for_each_superfamily(10) # [10, 20, 30, 40, 50, 60, 70, 80, 90]

count = 0
for sf, sfam_domain_id_set in new_sfam_dict.items():
    print(len(sfam_domain_id_set))
    count += len(sfam_domain_id_set)
print(f"number of prots: {count}")
print(f"num of superfamilies: {len(new_sfam_dict)}")