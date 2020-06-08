#!/usr/bin/env python3

"""
> plot_coverage_histo.py <

Uses seaborn to plot a generic coverage histogram of the high-quality methylated
positions common across in all post_datasets.
"""
import csv
import glob
import gzip
import re

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

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

# read post-subsampling data
subsamp_covs = sorted(glob.glob('../../2??.cov') + glob.glob('../../3??.cov'))
post_data = {}
valid_pos = []
for i in subsamp_covs:
    j = re.search('(\w{3})\.cov', i)[1]
    post_data[j] = []
    tsv_reader = csv.reader(open(i), delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        
        # append position to valid_pos
        scaf = row[0]
        pos = row[1]
        valid_pos.append((scaf, pos))
        
        cov = int(row[4]) + int(row[5])
        post_data[j].append(cov)
    
    post_data[j] = np.array(post_data[j])
valid_pos = list(set(valid_pos))

# read pre-subsampling data
presamp_file = gzip.open('../../filter_meth_pos_CX/merged_unfilt_covs/' + \
                         'compiled_coverage.pre-filt.meth_unmeth.tsv.gz', 'rt')
tsv_reader = csv.reader(presamp_file, delimiter='\t')
conds = next(tsv_reader)[2::2]
conds = [x.replace('.cov', '') for x in conds]
pre_data = {x:[] for x in conds}
for row in tsv_reader:
    if not row: continue
    
    # only process positions that are within valid_pos
    scaf = row[0]
    pos = row[1]
    if (scaf, pos) not in valid_pos: continue
    
    # compute, then store coverage in pre_data
    meth_cov = row[2::2]
    unmeth_cov = row[3::2]
    total_cov = [int(x) + int(y) for x, y in zip(meth_cov, unmeth_cov)]
    
    for n in range(len(conds)):
        pre_data[conds[n]].append(total_cov[n])

pre_data = {x:np.array(pre_data[x]) for x in pre_data}

# contrast pre-subsampled and post-subsampled plots OF THE SAME POSITIONS
sns.set_style('ticks')

fig, ax = plt.subplots(1, 2, figsize=(7, 3), sharey=True)
for c in conds:
    sns.distplot(pre_data[c][pre_data[c] <= 250],
                 hist=False,
                 kde=True, 
                 kde_kws={'alpha': 0.6},
                 color=assign_color(c), label=c,
                 ax=ax[0])
    sns.distplot(post_data[c][post_data[c] <= 250],
                 hist=False,
                 kde=True, 
                 kde_kws={'alpha': 0.6},
                 color=assign_color(c), label=c,
                 ax=ax[1])

ax[0].legend(loc=1, ncol=2)
ax[1].legend(loc=1, ncol=2)

ax[0].set_xlim(0, 250)
ax[1].set_xlim(0, 250)

sns.despine(offset=5, trim=True)
ax[0].set_xlabel('Coverage')
ax[1].set_xlabel('Coverage')
ax[0].set_ylabel('Density')

fig = plt.gcf()
    
# without bbox_inches, the saved figure has truncated axes.
fig.savefig('coverage_histo.pdf', bbox_inches='tight')
