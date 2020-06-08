#!/usr/bin/env python3

"""
> plot_correlation.meth.py <

Plots a correlation heatmap using seaborn.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# read in methylation %s: [0-100]
data = pd.read_table('../../all.filt.pct.tsv', usecols=range(2,17))
data_corr = data.corr(method='kendall')

mask = np.zeros((len(data_corr), len(data_corr)), 'int8')
np.fill_diagonal(mask, 1)

# plot the heatmap!
fig, ax = plt.subplots(figsize=(6, 6))
sns.set(font_scale=1.3)

cg = sns.clustermap(data_corr, method='ward',
                    vmin=0.45, vmax=0.75, mask=mask,
                    cmap='YlGnBu', annot=False, linewidth=0.5,
                    cbar_kws={'ticks': np.linspace(0.45, 0.75, 4)})

# override the ytick rotation
plt.setp(cg.ax_heatmap.get_yticklabels(), rotation=0)

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'corr_plot.meth.pdf'
fig.savefig(output_filename, bbox_inches='tight')
