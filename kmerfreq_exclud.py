from Bio import SeqIO
import os
import math
import argparse
import numpy as np
import itertools
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy.spatial.distance import mahalanobis
from numpy.linalg import inv
import numpy as np
from scipy.stats import chi2
from multiprocessing import Pool

# 生成所有4-mer

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer", type=int, default=4, help="kmer")
    parser.add_argument("-t", "--threads", type=int, default=4, help="threads")
    return parser.parse_args()

def kmer_count(seq):
    vec = np.zeros(len(kmers))
    seq = str(seq).upper()
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        if kmer in kmer_index:
            vec[kmer_index[kmer]] += 1
    if vec.sum() > 0:
        vec = vec / vec.sum()
    return vec

def process_genome(filepath):
    names = []
    X = []

    # ====== 读序列 ======
    for record in SeqIO.parse(filepath, "fasta"):
        if len(record.seq) > 2000:
            names.append(record.id)
            X.append(kmer_count(record.seq))

    # 如果contig太少，跳过
    if len(X) < 5:
        return fasta_path, [], []

    X = np.array(X)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)

    mean = np.mean(X_pca, axis=0)
    cov_matrix = np.cov(X_pca.T)

    try:
        inv_cov = inv(cov_matrix)
    except:
        return filepath, [], []
        # Mahalanobis distance
    dist = []
    for x in X_pca:
        d = mahalanobis(x, mean, inv_cov)
        dist.append(d)

    dist = np.array(dist)
    #threshold = np.mean(dist) + 2 * np.std(dist)
    df = X_pca.shape[1]
    threshold = chi2.ppf(0.95, df=df)
    keep_idx = dist <= threshold
    remove_idx = dist > threshold
    good_contigs = list(np.array(names)[keep_idx])
    bad_contigs = list(np.array(names)[remove_idx])
    
    return filepath, good_contigs, bad_contigs,X_pca, keep_idx, remove_idx

if __name__ == "__main__":

    args = parse_args()

    k = args.kmer
    n_jobs = args.threads

    global kmers, kmer_index

    kmers = [''.join(p) for p in itertools.product('ACGT', repeat=k)]
    kmer_index = {k:i for i,k in enumerate(kmers)}


    input_dir = 'path/to/genome'
    
    fafiles = []
    for file in os.listdir(input_dir):
        if file.endswith((".fa", ".fasta", ".fna")):
            filepath = os.path.join(input_dir, file)
            fafiles.append(filepath)

    with Pool(n_jobs) as pool:
        results = pool.map(process_genome, fafiles)

    #with open(f"kmerfreq_total_{k}mer", "w") as fout:
        #for genome, good, bad in results:
            #fout.write(f"{os.path.basename(genome).split('.')[0]}\t{','.join(bad)}\n")
    
    batch_size = 25
    n = len(results)
    cols = 5
    rows = 5
    for i in range(0, n, batch_size):
        subset = results[i:i+batch_size]
        fig, axes = plt.subplots(rows, cols, figsize=(cols*3, rows*3))
        axes = axes.flatten()

        for j, (genome, good_contigs, bad_contigs, X_pca, keep_idx, remove_idx) in enumerate(subset):
            ax = axes[j]

            # good contigs
            ax.scatter(
                X_pca[keep_idx, 0],
                X_pca[keep_idx, 1],
                s=5,
                alpha=0.7
            )

            # bad contigs
            ax.scatter(
                X_pca[remove_idx, 0],
                X_pca[remove_idx, 1],
                s=20,
                marker="x",
                alpha=0.9
            )

            ax.set_title(os.path.basename(genome).split('.')[0], fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
        
        for k in range(len(subset), len(axes)):
            fig.delaxes(axes[k])

        plt.tight_layout()
        
        outname = f"{k}mer_pca_part_{i//batch_size}.png"
        plt.savefig(outname, dpi=300)

        plt.close()  # 很重要，防止内存爆
