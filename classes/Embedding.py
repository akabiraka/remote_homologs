import sys
sys.path.append("../remote_homologs")

import torch
import torch.nn as nn
import numpy as np


class SequenceIntEmbeddingLayer(nn.Module):
    def __init__(self) -> None:
        super(SequenceIntEmbeddingLayer, self).__init__()
        self.AA = "ARNDCQEGHILKMFPSTWYV"
        self.AA_dict = {aa:i+1 for i, aa in enumerate(self.AA)}
        self.vocab_size = len(self.AA)+1 # 20 amino acid + 1 padding index (0)

    def truncate_and_padd(self, x:torch.Tensor, max_len):
        x = x[:max_len]
        temp = torch.zeros(max_len)
        temp[:x.shape[0]] = x
        return temp.to(torch.long)

    def forward(self, seq:str, max_len=None):
        #int amino-acid encoding [1, 20], 0 as padding index
        int_enc = torch.tensor([self.AA_dict[aa] for aa in seq])
        if max_len is not None:
            int_enc = self.truncate_and_padd(int_enc, max_len)
        return int_enc #shape: seq_len or max_len

# sample usage
# seq_int_embed_layer = SequenceIntEmbeddingLayer(5)
# print(seq_int_embed_layer("ARN"))#, 10))



class SequenceLearnableEmbeddingLayer(nn.Module):
    def __init__(self, embed_dim=5) -> None:
        super(SequenceLearnableEmbeddingLayer, self).__init__()
        self.seq_int_embed_layer = SequenceIntEmbeddingLayer()
        self.embed = nn.Embedding(self.seq_int_embed_layer.vocab_size, embed_dim, padding_idx=0)
        self.embed_dim = embed_dim

    def forward(self, seq, requires_grad=True):
        seq_rep = self.seq_int_embed_layer(seq)
        if requires_grad:
            seq_rep = self.embed(seq_rep) * np.sqrt(self.embed_dim)
        else: 
            with torch.no_grad():
                seq_rep = self.embed(seq_rep) * np.sqrt(self.embed_dim)
            
        return seq_rep

# sample usage
# seq = "ARN"
# model = SequenceLearnableEmbeddingLayer()
# seq_rep = model(seq, requires_grad=False)
# print(seq_rep)



class SequenceOnehotEmbeddingLayer(nn.Module):
    def __init__(self) -> None:
        super(SequenceOnehotEmbeddingLayer, self).__init__()
        self.seq_int_embed_layer = SequenceIntEmbeddingLayer()

    def forward(self, seq, requires_grad=True):
        seq_rep = self.seq_int_embed_layer(seq)
        seq_rep = torch.nn.functional.one_hot(seq_rep, self.seq_int_embed_layer.vocab_size)
        seq_rep = seq_rep.to(dtype=torch.float32)
        seq_rep.requires_grad_(requires_grad)
        
        return seq_rep 
        
# sample usage
# seq = "ARN"
# model = SequenceOnehotEmbeddingLayer()
# seq_rep = model(seq, requires_grad=False)
# print(seq_rep)

import esm
class SequenceESMEmbeddingLayer(nn.Module):
    def __init__(self) -> None:
        super(SequenceESMEmbeddingLayer, self).__init__()
        self.esm_model, alphabet = esm.pretrained.esm2_t33_650M_UR50D() #esm2_t6_8M_UR50D() #esm1_t12_85M_UR50S()
        self.last_layer_num = 33 #12, 6 
        self.esm_batch_converter = alphabet.get_batch_converter()

    def forward(self, id, seq, requires_grad=False):
        uniprotid, batch_strs, seq_tokens = self.esm_batch_converter([(id, seq)])

        if requires_grad:
            results = self.esm_model(seq_tokens, repr_layers=[self.last_layer_num], return_contacts=False)
        else:
            with torch.no_grad():
                results = self.esm_model(seq_tokens, repr_layers=[self.last_layer_num], return_contacts=False)
        
        token_reps = results["representations"][self.last_layer_num] # tokens are amino-acid, mask or special tokens
        # print(token_reps)

        seq_rep = token_reps[0, 1:len(seq)+1]
        return seq_rep

# esm_layer = SequenceESMEmbeddingLayer()
# seq_rep = esm_layer("id1", "MKTV", requires_grad=False)
# print(seq_rep.shape) #[4, 768]


class SequenceEmbeddingLayer(nn.Module):
    def __init__(self, embed_format, embed_dim=None):
        super(SequenceEmbeddingLayer, self).__init__()
        if embed_format == "onehot":
            self.seq_embed_layer = SequenceOnehotEmbeddingLayer()
        elif embed_format == "learnable":
            if embed_dim==None:
                raise BaseException("embed_dim' of 'int' type is expected.")
            self.seq_embed_layer =  SequenceLearnableEmbeddingLayer(embed_dim)
        elif embed_format == "esm":
            self.seq_embed_layer =  SequenceESMEmbeddingLayer()
        else:
            raise NotImplementedError(f"'{embed_format}'" + "not in 'onehot, learnable, esm'")

    def forward(self, id, seq, requires_grad=None):
        # i1, s1, true
        
        if isinstance(self.seq_embed_layer, SequenceESMEmbeddingLayer):
            seq_rep = self.seq_embed_layer(id, seq, requires_grad)
        else:
            seq_rep = self.seq_embed_layer(seq, requires_grad)
        
        return seq_rep

# test cases
# seqs = [("id1", "MKTV"),
#         ("id2", "KA")]

# seq_embed_layer = SequenceEmbeddingLayer(embed_format="onehot")  
# seq_reps = seq_embed_layer(seqs, max_len=3, requires_grad=True)        
# print(seq_reps)

# seq_embed_layer = SequenceEmbeddingLayer(embed_format="learnable", embed_dim=10)  
# seq_reps = seq_embed_layer(seqs, max_len=3, requires_grad=True)        
# print(seq_reps)

# seq_embed_layer = SequenceEmbeddingLayer(embed_format="esm")  
# seq_reps = seq_embed_layer(seqs, max_len=3, requires_grad=True)        
# print(seq_reps)