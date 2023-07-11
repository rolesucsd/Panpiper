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


    # Load the color information from the phylogroups
    # Generate a random color for each unique phylogroup
    unique_phylogroups = phylogroup['Phylogroup'].unique()
    random.shuffle(unique_phylogroups)
    phylogroup_colors = {group: '#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for group in unique_phylogroups}
   # phylogroup_colors = {'1': "#89B4AE", '0':"#B3A3CA"}

    # Map the phylogroup colors to each sample in the dataframe
    row_colors = phylogroup.set_index('Sample')['Phylogroup'].map(phylogroup_colors)


    # Plot clusters using PCA
    pca = PCA(n_components=2)
    pca.fit(mash)
    pca_scores = pca.transform(mash)
    pca_df = pd.DataFrame(pca_scores, columns=['PC1', 'PC2'])
    pca_df['Cluster'] = cluster_labels
    plt.figure()
    sns.scatterplot(x='PC1', y='PC2', hue='Cluster', palette=phylogroup_colors, data=pca_df, legend=False)  # Add row_colors
    plt.title('Clusters of samples by PCA')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)')  # Add explained variance ratio to x-axis label
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)')  # Add explained variance ratio to y-axis label

    # Move the legend to the right of the graph
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

    # Save the plot to a PNG file
    plt.savefig(f"{output_folder}/mash_PCA.png", dpi=300, bbox_inches='tight')


    # Perform MDS on the distance matrix
    mds = MDS(n_components=2, dissimilarity='precomputed')
    embedding = mds.fit_transform(mash)

    # Create a DataFrame with the MDS coordinates and cluster labels
    mds_df = pd.DataFrame(embedding, columns=['MDS1', 'MDS2'])
    mds_df['Cluster'] = cluster_labels

    # Plot the clusters using MDS
    plt.figure()
    sns.scatterplot(x='MDS1', y='MDS2', hue='Cluster', palette=phylogroup_colors, data=mds_df, legend=False)
    plt.title('Clusters of samples by MDS')
    plt.xlabel('MDS1')
    plt.ylabel('MDS2')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.savefig(f"{output_folder}/mash_MDS.png", dpi=300, bbox_inches='tight')

    # Perform t-SNE on the distance matrix
    perplexity_value = min(30, len(mash) - 1)  # Choose a maximum perplexity value based on the number of samples
    tsne = TSNE(n_components=2, perplexity=perplexity_value)
    embedding = tsne.fit_transform(mash)

    # Create a DataFrame with the t-SNE coordinates and cluster labels
    tsne_df = pd.DataFrame(embedding, columns=['TSNE1', 'TSNE2'])
    tsne_df['Cluster'] = cluster_labels

    # Plot the clusters using t-SNE
    plt.figure()
    sns.scatterplot(x='TSNE1', y='TSNE2', hue='Cluster', palette=phylogroup_colors, data=tsne_df, legend=False)
    plt.title('Clusters of samples by t-SNE')
    plt.xlabel('TSNE1')
    plt.ylabel('TSNE2')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.savefig(f"{output_folder}/mash_tSNE.png", dpi=300, bbox_inches='tight')

    # Plot heatmap of clusters with color annotations and inverted colormap
    clustermap = sns.clustermap(np.log1p(mash), cmap='binary_r', row_cluster=True,
                            row_colors=row_colors, row_linkage=linkage_matrix, col_linkage=linkage_matrix)

    # Adjust the size of the legend
    cbar = clustermap.cax.get_position()
    clustermap.cax.set_position([cbar.x0 + 0, cbar.y0 + 0.1, cbar.width * 0.5, cbar.height * 0.5])  # Make the legend 5 times smaller

    # Remove tick marks from the axis
    clustermap.ax_heatmap.set_xticks([])
    clustermap.ax_heatmap.set_yticks([])

    # Remove axis names
    clustermap.ax_heatmap.set_xlabel('')
    clustermap.ax_heatmap.set_ylabel('')
    # Save the heatmap to a file
    plt.savefig(f"{output_folder}/mash_heatmap.png", dpi=300, bbox_inches='tight')
    
    # Perform UMAP on the distance matrix
    reducer = umap.UMAP(n_components=2)
    embedding = reducer.fit_transform(mash)

    # Create a DataFrame with the UMAP coordinates and cluster labels
    umap_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
    umap_df['Cluster'] = cluster_labels

    # Plot the clusters using UMAP
    plt.figure()
    sns.scatterplot(x='UMAP1', y='UMAP2', hue='Cluster', palette=phylogroup_colors, data=umap_df, legend=False)
    plt.title('Clusters of samples by UMAP')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.savefig(f"{output_folder}/mash_UMAP.png", dpi=300, bbox_inches='tight')

    return(phylogroup)