#!/usr/bin/env python3

"""
> plot_meth_levels.sequence_context.py <

Uses seaborn to check whether methylation levels of CpG vs. CHG vs. CHH
positions are different from each other.
"""
import csv
import glob
import re

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def assign_color(label):
    if label[:3] == 'cpg':
        return '#fed98e'
    elif label[:3] == 'chg':
        return '#fe9929'
    elif label[:3] == 'chh':
        return '#d95f0e'

# read data
input_files = sorted(glob.glob('../../all.*.c20.*.cov'))
data = {}
for i in input_files:
    j = re.search('\.c20\.(\w{3})', i).group(1)
    data[j] = []
    tsv_reader = csv.reader(open(i), delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        
        meth_pct = float(row[3])
        data[j].append(meth_pct)
    
    data[j] = np.array(data[j])

# seaborn
sns.set_style('white')
sns.set_style('ticks')

fig, ax = plt.subplots(figsize=(4, 2))
for d in sorted(data):
    sns.distplot(data[d],
                 hist=False,
                 kde=True, 
                 hist_kws={'alpha': 0.3},
                 color=assign_color(d), label=d)

plt.legend(loc=1, ncol=1)

ax.set_xlim(0, 100)

sns.despine(offset=5, trim=True)
ax.set_xlabel('Methylation level (%)')
ax.set_ylabel('Density')

fig = plt.gcf()
    
# without bbox_inches, the saved figure has truncated axes
fig.savefig('meth_levels.sequence_context.pdf', bbox_inches='tight')
