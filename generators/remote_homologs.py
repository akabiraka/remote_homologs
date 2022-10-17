import sys
sys.path.append("../remote_homologs")

import random
import utils as Utils
from dict_computation import get_fam_and_sfam_dict, get_cluster_dict


def exclude_SF_members(sfam_domain_id_set, domain_id_set, level="Family"):
    x = sfam_domain_id_set.intersection(domain_id_set)
    if len(x) > 1:                  
        print(f"    excluding {len(x)-1} SF members that belong to same {level}...")                    
        sfam_domain_id_set = sfam_domain_id_set - x
        sfam_domain_id_set = sfam_domain_id_set | set(random.sample(list(x), 1))
        if level=="Cluster": # debug: cluster level redundency check
            print("cluster")
    return sfam_domain_id_set


seq_identity_th = 35 # [20, 25, 35, 40, 95]
cluster_dict = get_cluster_dict(seq_identity_th)
fam_dict, sfam_dict = get_fam_and_sfam_dict()
n_sfams = len(sfam_dict)

new_sfam_dict = dict()
for i, (sf, sfam_domain_id_set) in enumerate(sfam_dict.items()):
    n_mem_b4_exlusion = len(sfam_domain_id_set)

    for fam, fam_domain_id_set in fam_dict.items():
        sfam_domain_id_set = exclude_SF_members(sfam_domain_id_set, fam_domain_id_set, level="Family")

    for clus, clus_domain_id_set in cluster_dict.items():
        sfam_domain_id_set = exclude_SF_members(sfam_domain_id_set, clus_domain_id_set, level="Cluster")
        
    
    print(f"{i}/{n_sfams}: {n_mem_b4_exlusion} and {len(sfam_domain_id_set)} SF members b4 and after exclusion")
    print()

    if len(sfam_domain_id_set) > 1: # excluding superfamilies that has no remote homologs.
        new_sfam_dict[sf] = sfam_domain_id_set
    # if i==1000: break

print(len(new_sfam_dict))
Utils.save_as_pickle(new_sfam_dict, f"data/tmp/remote_homolog_dict_at_{seq_identity_th}_sfam.pkl")
# print(new_sfam_dict)