#!/usr/bin/env python3

"""
> plot_expr_vs_meth_status.py <

Plot a histogram to show that methylated genes tend to be higher expressed 
than unmethylated ones!
"""
import csv
import math
import statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import scipy.stats

# get list of methylated genes
meth_genes = []
for line in open('smic_1pos.txt'):
    meth_genes.append(line.strip().split(' ')[1])

# read in expression TPMs: [0, inf)
data = pd.read_table('../../rnaseq/normalised_abundances.all.tsv',
                     index_col=0)

# remove rows that contain a single 0 (i.e., gene must be expressed in all 
# replicates)
data = data.loc[~(data==0).any(axis=1)]

# log10-transform TPMs
data = data.applymap(lambda x: np.log10(x + 1))

# calculate mean
data['mean'] = data.mean(axis=1)

# split dataset into methylated vs. unmethylated genes
meth_data = data[data.index.isin(meth_genes)]
unmeth_data = data[~data.index.isin(meth_genes)]

# do t-test to show that distributions means are significantly different
t_stat, p_val = scipy.stats.ttest_ind(meth_data['mean'], unmeth_data['mean'])

print ('geomean of meth genes: ', meth_data['mean'].mean())
print ('geomean of unmeth genes: ', unmeth_data['mean'].mean())
print ('t-test p value: ', p_val)

# seaborn
#sns.set_style('white')
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(3, 3))

sns.distplot(meth_data['mean'],
             hist=False,
             color='#fc9272',
             label='Methylated genes')
sns.distplot(unmeth_data['mean'],
             hist=False,
             color='#fee0d2',
             label='Unmethylated genes')

ax.set_xlim(0, 3)
ax.set_ylim(0, 1.2)
sns.despine(offset=10, trim=True)
ax.set_xlabel('log10 (mean expression + 1)')
ax.set_ylabel('Relative frequency')

plt.legend(loc=1, ncol=1)

plt.tight_layout()

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'expr_vs_meth_status.pdf'
fig.savefig(output_filename, bbox_inches='tight')
