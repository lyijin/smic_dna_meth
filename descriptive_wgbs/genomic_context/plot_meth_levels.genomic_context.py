#!/usr/bin/env python3

"""
> plot_meth_levels.genomic_context.py <

Uses seaborn to plot a histogram of the methylation to contrast methylation
across genomic contexts (exon/intron/intergenic).
"""
import csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def define_genomic_feature(descriptor):
    """
    From descriptor, parse and return whether position is 'intergenic',
    'exon' or 'intron'.
    """
    if 'Exon_' in descriptor:
        return 'exon'
    elif 'Intron_' in descriptor:
        return 'intron'
    else:
        return 'intergenic'

data = pd.read_table('../../all.bona_fide_meth_pos.annot.c20.cov',
                     header=None,
                     usecols=[3, 6, 10],
                     names=['meth_pct', 'gene', 'ei'])

data['temp'] = data['gene'] + data['ei']
data['temp'] = data['temp'].fillna(value='intergenic')
data['genomic_feature'] = data['temp'].apply(define_genomic_feature)

data = data.drop('temp', axis=1)
data = data.drop('gene', axis=1)
data = data.drop('ei', axis=1)

# seaborn
sns.set_style('white')
sns.set_style('ticks')

fig, ax = plt.subplots(figsize=(5, 5))
color = {'exon': '#2166ac', 'intron': '#67a9cf', 'intergenic': '#f4a582'}
for f in ['exon', 'intron', 'intergenic']:
    sns.distplot(data[data['genomic_feature'] == f]['meth_pct'],
                 hist=False,
                 kde=True,
                 kde_kws={'bw': 2},
                 color=color[f],
                 label=f.title())

plt.legend(loc=1, ncol=1)

plt.xlim(0, 100)

sns.despine(offset=5, trim=True)
ax.set_xlabel('Methylation level (%)')
ax.set_ylabel('Density')

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'meth_levels.genomic_context.pdf'
fig.savefig(output_filename, bbox_inches='tight')
