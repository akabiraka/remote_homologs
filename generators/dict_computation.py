import sys
sys.path.append("../remote_homologs")

import pandas as pd
import random
import os.path as osp
from Bio import SeqIO

import utils as Utils

def update_dict(x:dict, key, value):
    if key not in x:
        x[key] = set()
        x[key].add(value)
    else:
        x[key].add(value)
    return x


def get_fam_and_sfam_dict():
    sfam_dict_path, fam_dict_path = "data/tmp/sfam_dict.pkl", "data/tmp/fam_dict.pkl"
    if osp.exists(sfam_dict_path) and osp.exists(fam_dict_path): 
        sfam_dict = Utils.load_pickle(sfam_dict_path)
        fam_dict = Utils.load_pickle(fam_dict_path)
        print(f"    num of families: {len(fam_dict)}") # num of FA: 5936
        print(f"    num of superfamilies: {len(sfam_dict)}") # num of SF: 2816
        return fam_dict, sfam_dict

    print("Computing Family and Superfamily dictionary.")
    col_names = "FA-DOMID FA-PDBID FA-PDBREG FA-UNIID FA-UNIREG SF-DOMID SF-PDBID SF-PDBREG SF-UNIID SF-UNIREG SCOPCLA".split(" ")
    df = pd.read_csv("data/downloads/scop-cla-latest.txt", names=col_names, delim_whitespace=True, skiprows=6, header=None)
    # print(df.columns)
    print(df.head())

    fam_dict, sfam_dict = dict(), dict()
    for i, row in df.iterrows():
        x = row["SCOPCLA"].split(",") 
        TP, CL, CF, SF, FA = x[0][3:], x[1][3:], x[2][3:], x[3][3:], x[4][3:]
        # print(TP, CL, CF, SF, FA)

        fam_dict = update_dict(fam_dict, FA, str(row["SF-DOMID"]))
        sfam_dict = update_dict(sfam_dict, SF, str(row["SF-DOMID"]))

        # if i==1000: break
    print(f"    num of families: {len(fam_dict)}") # num of FA: 5936
    print(f"    num of superfamilies: {len(sfam_dict)}") # num of SF: 2816
    # print(fam_dict)

    Utils.save_as_pickle(fam_dict, fam_dict_path)
    Utils.save_as_pickle(sfam_dict, sfam_dict_path)
    return fam_dict, sfam_dict


def get_cluster_dict(seq_identity_th):
    sfam_clus_file_path = f"data/downloads/clusters/scop_sf_represeq_lib_latest.{seq_identity_th}clust"
    f = open(sfam_clus_file_path, mode="r")

    clus_dict_path = f"data/tmp/cluster_dict_at_seq_identity_th_{seq_identity_th}.pkl"
    if osp.exists(clus_dict_path):
        cluster_dict = Utils.load_pickle(clus_dict_path)
        print(f"    num of clusters: {len(cluster_dict)}") #18575 at 40% identity
        return cluster_dict

    print("Computing cluster dictionary.")
    cluster_dict = dict()
    for i, line in enumerate(f.readlines()):
        cluster_dict[f"c{i}"] = set(line.rstrip().split(" "))
        # break

    print(f"    num of clusters: {len(cluster_dict)}") #18575 at 40% identity
    Utils.save_as_pickle(cluster_dict, clus_dict_path)
    # print(cluster_dict)
    return cluster_dict
# sample usage
# get_cluster_dict(95)

def get_fasta_seq_dict():
    fasta_seq_dict_path = f"data/tmp/fasta_seq_dict.pkl"
    if osp.exists(fasta_seq_dict_path):
        fasta_seq_dict = Utils.load_pickle(fasta_seq_dict_path)
        print(f"    num of seq: {len(fasta_seq_dict)}")
        return fasta_seq_dict


    inp_file_path = f"data/downloads/scop_sf_represeq_lib_latest.fa"

    fasta_seq_dict = {}
    for record in SeqIO.parse(inp_file_path, "fasta"):
        # print(record)
        fasta_seq_dict[record.id] = {"description": record.description, "seq": str(record.seq)}
        # break
        

    Utils.save_as_pickle(fasta_seq_dict, fasta_seq_dict_path)
    print(f"    num of seq: {len(fasta_seq_dict)}")
    return fasta_seq_dict


def get_SCOP_node_label_dict():
    df_cls_des = pd.read_csv("data/downloads/scop-des-latest.txt", skiprows=6)
    scop_node_label_dict = {}
    for i in range(len(df_cls_des)):
        key = str(df_cls_des.loc[i].values[0].split(" ")[0])
        label = " ".join(df_cls_des.loc[i].values[0].split(" ")[1:])
        scop_node_label_dict[key] = label
        # if i==5: break
    print(f"    num of SCOP nodes: {len(scop_node_label_dict)}")
    return scop_node_label_dict