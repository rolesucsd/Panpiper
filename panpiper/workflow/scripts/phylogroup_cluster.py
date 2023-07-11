import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import permanova
from scipy.cluster.hierarchy import linkage, fcluster, cut_tree
from skbio import DistanceMatrix
import random
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
from sklearn.manifold import TSNE
import os
import umap

def cluster_samples(mash_file, output_folder):
    matplotlib.use('Agg')

    # Load data
    mash = pd.read_csv(mash_file, sep="\t", index_col=0)
    # Update row and column names
    mash.index = [os.path.splitext(os.path.basename(row))[0] for row in mash.index]
    mash.columns = [os.path.splitext(os.path.basename(col))[0] for col in mash.columns]

    # Convert to distance matrix
    dist = pdist(mash, metric='euclidean')

    # Cluster using hierarchical clustering
    linkage_matrix = linkage(dist, method='ward')
    plt.title('Hierarchical clustering')

    # Perform clustering with the best number of clusters
    max_d = 0.23 * np.max(linkage_matrix[:, 2])
    clusters = cut_tree(linkage_matrix, height=max_d)

    cluster_labels = [str(i + 1) for i in clusters.flatten()]

    # Plot clusters
    cluster_labels = [str(label).strip('[]') for label in clusters]  # Remove "[" and "]" characters
    phylogroup = pd.DataFrame({'Phylogroup': cluster_labels, 'Sample': mash.index})
    phylogroup.to_csv(f"{output_folder}/phylogroups.txt", sep="\t", index=False)
    
    phylogroup_wide = phylogroup.copy()

    # Convert the 'Phylogroup' column into dummy variables
    phylogroup_dummies = pd.get_dummies(phylogroup_wide['Phylogroup'])

    # Concatenate the 'Sample' column with the dummy variables
    phylogroup_wide = pd.concat([phylogroup_wide['Sample'], phylogroup_dummies], axis=1)

    # Set 'Sample' column as the DataFrame index
    phylogroup_wide.set_index('Sample', inplace=True)

    # Replace the dummy variable values (0 and 1) with the actual phylogroup names
    phylogroup_wide = phylogroup_wide.apply(lambda x: x.astype(str).str.replace('0', '').str.replace('1', x.name))
    
    # Replace all numbers with 1
    phylogroup_wide = phylogroup_wide.replace(np.nan, 0).replace('[0-9]+', 1, regex=True)

    # Fill empty spaces with 0s
    phylogroup_wide = phylogroup_wide.mask(phylogroup_wide != 1, other=0)
    phylogroup_wide = phylogroup_wide.fillna(0)
    
    phylogroup_wide.to_csv(f"{output_folder}/phylogroup_wide.txt", sep="\t", index=True)

    return(phylogroup)
