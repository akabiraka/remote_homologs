
import sys
sys.path.append("../remote_homologs")

import os.path as osp
from Bio import SeqIO

import utils as Utils
from classes.Embedding import SequenceEmbeddingLayer



def get_fasta_seq_dict(inp_file_path="data/downloaded/scop_sf_represeq_lib_latest.fa", fasta_seq_dict_path="data/generated/fasta_seq_dict.pkl"):
    
    if osp.exists(fasta_seq_dict_path):
        fasta_seq_dict = Utils.load_pickle(fasta_seq_dict_path)
        print(f"    num of seq: {len(fasta_seq_dict)}")
        return fasta_seq_dict

    fasta_seq_dict = {}
    for record in SeqIO.parse(inp_file_path, "fasta"):
        # print(record)
        fasta_seq_dict[record.id] = {"description": record.description, "seq": str(record.seq)}
        # break
        
    Utils.save_as_pickle(fasta_seq_dict, fasta_seq_dict_path)
    print(f"    num of seq: {len(fasta_seq_dict)}")
    return fasta_seq_dict

    
def compute_all_seqs_embedding(fasta_seq_dict, seq_embed_func, out_path):
    num_of_seq = len(fasta_seq_dict)
    for i, (key, item) in enumerate(fasta_seq_dict.items()):
        # print(key, item["seq"])    
        seq_rep = seq_embed_func(key, item["seq"], requires_grad=False).detach().cpu().numpy() # this is on amino acid level
        Utils.save_as_pickle(seq_rep, f"{out_path}/{key}.pkl")
        # seq_rep = Utils.load_pickle(f"{out_path}/{key}.pkl") # tesing if saving is ok
        print(f"{i}|{num_of_seq}, {seq_rep.shape}")
        # if i==10: break
        
model_name = "esm2_t33_650M_UR50D" #esm2_t6_8M_UR50D, esm2_t33_650M_UR50D, esm1_t12_85M_UR50S()
fasta_seq_dict = get_fasta_seq_dict()
seq_embed_func = SequenceEmbeddingLayer(embed_format="esm")  
compute_all_seqs_embedding(fasta_seq_dict, seq_embed_func, out_path=f"data/generated/seqs_aalevel_embeddings_by_{model_name}/")