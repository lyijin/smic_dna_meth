#!/usr/bin/env python3

"""
> plot_2d_pca.py <

Plot PCA of all methylation levels from all samples.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def assign_color(label):
    if label[:2] == '26':
        if label[-1] in ['A', 'B', 'C']:
            return '#abcee3'
        else:
            return '#1f78b4'
    elif label[:2] == '29':
        return '#33a02c'
    else:
        if label[-1] in ['A', 'B', 'C']:
            return '#fb9a99'
        else:
            return '#e31a1c'

# read in methylation %s: [0-100]
data = pd.read_table('../../all.filt.pct.tsv', usecols=range(2,17))

# PCA stuff
pca = PCA(n_components=3)
pca_array = pca.fit_transform(data.transpose())
pca_evr = [round(x * 100, 1) for x in pca.explained_variance_ratio_]
    
# plot PCA
fig, ax_array = plt.subplots(1, 3, figsize=(13, 4))
point_labels = data.columns.tolist()

labels = [x[:2] for x in point_labels]
markers = {'26': 'o', '29': '^', '32': 's'}

pc1_values = [x[0] for x in pca_array]
offset = 0.02 * (max(pc1_values) - min(pc1_values))

for m, ax in enumerate(ax_array):
    x_pc = [0, 0, 1]
    y_pc = [1, 2, 2]
    for l in sorted(set(labels)):
        x = []
        y = []
        for n in range(len(pca_array)):
            if labels[n] == l:
                x.append(pca_array[n][x_pc[m]])
                y.append(pca_array[n][y_pc[m]])
                ax.text(pca_array[n][x_pc[m]] + offset, pca_array[n][y_pc[m]],
                        point_labels[n], ha='left', va='center')
        
        ax.scatter(x, y, c=assign_color(l), marker=markers[l], label=l)
    
    ax.set_xlabel('PC{} ({}%)'.format(x_pc[m] + 1, pca_evr[x_pc[m]]))
    ax.set_ylabel('PC{} ({}%)'.format(y_pc[m] + 1, pca_evr[y_pc[m]]))

plt.tight_layout()
#plt.show()

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'all_strains.2d_pca.pdf'
fig.savefig(output_filename, bbox_inches='tight')
