from openTSNE.callbacks import ErrorLogger
import umap
import pyVectorizer
import pandas as pd
from sklearn.decomposition import PCA
import metagenomics_processor
import numpy as np
# from sklearn.manifold import TSNE
from openTSNE import TSNE

def process(input_path, k,pca_dim, dim, iter, lr, perp, output_path):
    vecs = np.array(pyVectorizer.vectorize_file(input_path, k)).astype(np.float32)
    metagenomics_processor.normalize_over_axis1(vecs)
    # pca = PCA(n_components=pca_dim)
    # pca_data = pca.fit_transform(vecs)
    # print(pca_data)

    # tsne = TSNE(
    #         n_components=dim,
    #         n_iter=max(iter,250),
    #         learning_rate=lr,
    #         metric="euclidean",
    #         callbacks=ErrorLogger(),
    #         n_jobs=8,
    #         perplexity=perp,
    #         random_state=1131,
    #     )

    # r = tsne.fit_transform(vecs)
    # r = tsne.fit(vecs)

    reducer = umap.UMAP()
    r = reducer.fit_transform(vecs)

    print(r)

    dataset = pd.DataFrame({'a0': r[:, 0], 'a1': r[:, 1]}) if (r.shape[1] == 2) else pd.DataFrame(
        {'a0': r[:, 0], 'a1': r[:, 1], 'a2': r[:, 2]})


    if output_path: dataset.to_csv(output_path)

