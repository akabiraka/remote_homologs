import sys
sys.path.append("../remote_homologs")

import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import utils as Utils
from generators.dict_computation import get_SCOP_node_label_dict


def plot_embedding(X, Y, th):
    scop_node_label_dict = get_SCOP_node_label_dict()
    _, ax = plt.subplots()
    Y_set = set(Y)
    print(f"    num of superfamilies: {len(Y_set)}")
    print(f"    num of datapoints: {len(Y)}")
    for i, y in enumerate(Y_set):
        embedding = X[Y==y]
        # print(embedding.shape)
        # colors = np.repeat(color_names[i], len(embedding))
        ax.scatter(embedding[:, 0], embedding[:, 1], label=scop_node_label_dict[y], alpha=0.425)
    # ax.legend()
    # ax.legend(bbox_to_anchor=(.5, 1.6), loc='lower center', ncol=3)
    plt.show()
    plt.axis('off')
    # plt.savefig(f"outputs/images/remote_homologs_at_th_{th}.png", dpi=300, format="png", bbox_inches='tight', pad_inches=0.0)

for seq_identity_th in  [20]: # seq_identity_th #[20, 25, 35, 40, 95]
    data = Utils.load_pickle(f"data/generated/remote_homolog_dicts/at_{seq_identity_th}_sfam.pkl")
    X = np.stack(data[:, 2])
    Y = data[:, 1]

    tsne= TSNE(
            n_components=2,
            init="pca",
            learning_rate="auto",
            n_iter=5000,
            n_iter_without_progress=300,
            random_state=0,
            # perplexity=max(2, n*.25),
            verbose=1,
        )

    X = tsne.fit_transform(X, Y)
    plot_embedding(X, Y, seq_identity_th)
    # print(seq_identity_th)
    
    # print(f"    KL: {tsne.kl_divergence_}\n")