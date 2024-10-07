from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder

def ohe_seqs(seqs):
    ohe = OneHotEncoder(handle_unknown='ignore')
    seqs_ohe = ohe.fit_transform(seqs)

    return seqs_ohe

def pca_ohe_seqs(seqs_ohe):
    pca = PCA(n_components=2, svd_solver='arpack')
    pca.fit(seqs_ohe)
    return pca.transform(seqs_ohe)

def plot_pca(pca, pca_title, alpha=1):
    fig, ax = plt.subplots()
    ax.set_title(pca_title)
    ax.scatter(pca[:, 0], pca[:, 1], s=1, alpha=alpha)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    return fig, ax