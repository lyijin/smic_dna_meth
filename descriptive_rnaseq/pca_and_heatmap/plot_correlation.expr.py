#!/usr/bin/env python3

"""
> plot_correlation.expr.py <

Plots a correlation heatmap using seaborn.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# read in expression TPMs: [0, inf)
data = pd.read_table('../../rnaseq/normalised_abundances.all.tsv',
                     index_col=0)

# remove rows that contain a single 0 (i.e., gene must be expressed in all 
# replicates)
data = data.loc[~(data==0).any(axis=1)]

# log10-transform TPMs
data = data.applymap(lambda x: np.log10(x + 1))

# do kendall correlation
data_corr = data.corr(method='kendall')

mask = np.zeros((len(data_corr), len(data_corr)), 'int8')
np.fill_diagonal(mask, 1)

# plot the heatmap!
fig, ax = plt.subplots(figsize=(4, 4))
sns.set(font_scale=1.5)

cg = sns.clustermap(data_corr, method='ward',
                    vmin=0.80, vmax=0.92,
                    mask=mask,
                    cmap='YlGnBu', linewidth=0.5,
                    cbar_kws={'ticks': [0.80, 0.84, 0.88, 0.92]})

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'corr_plot.expr.pdf'
fig.savefig(output_filename, bbox_inches='tight')
