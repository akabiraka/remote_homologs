import sys
sys.path.append("../remote_homologs")

import time
import torch
import esm
import numpy as np
import utils as Utils

seq_identity_th = 35 #[20, 25, 35, 40, 95]
remote_homolog_dict = Utils.load_pickle(f"data/tmp/remote_homolog_dict_at_{seq_identity_th}_sfam.pkl")
fasta_seq_dict = Utils.load_pickle(f"data/tmp/fasta_seq_dict.pkl")

# print(fasta_seq_dict['8037501']["seq"])

device = "cpu"#"cuda" if torch.cuda.is_available() else "cpu" # "cpu"#
esm1b, alphabet = esm.pretrained.esm1_t12_85M_UR50S() # model
esm1b = esm1b.to(device)
esm1b_batch_converter = alphabet.get_batch_converter()

def get_seq_rep(idwithseq):
    uniprotid, batch_strs, seq_tokens = esm1b_batch_converter(idwithseq)

    with torch.no_grad():
        results = esm1b(seq_tokens, repr_layers=[12], return_contacts=False)
    seq_rep = results["representations"][12] #1, max_seq_len, esmb_embed_dim
    seq_rep.squeeze_(0)
    seq_rep = seq_rep[1 : len(seq) + 1].mean(0) # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    return seq_rep.detach().cpu().numpy()

start = time.time()
at_least_n_proteins_per_sfam = [1, 5, 10, 20, 30]
for n in at_least_n_proteins_per_sfam:
    data = []
    for i, (sfam, remote_homologs) in enumerate(remote_homolog_dict.items()):
        if len(remote_homologs) < n: continue
        print(i, sfam, remote_homologs)
        
        for remhom in remote_homologs:
            seq = fasta_seq_dict[remhom]["seq"]
            idwithseq = [(remhom, seq)]
            seq_rep = get_seq_rep(idwithseq) # shape: 768

            data.append(np.array([remhom, sfam, seq_rep], dtype=object))
        #     break
        # if i>50: break
    end = time.time()
    print(f"time required: {(end-start)/60}") # in minutes
    print(len(data))

    Utils.save_as_pickle(np.array(data), f"data/tmp/remote_homologs/remote_homologs_at_th_{seq_identity_th}_n_{n}.pkl")