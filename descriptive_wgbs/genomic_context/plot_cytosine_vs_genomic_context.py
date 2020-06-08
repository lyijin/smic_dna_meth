#!/usr/bin/env python3

"""
> plot_cytosine_vs_genomic_context.py <

Uses seaborn to plot a histogram of the methylation to contrast methylation
across genomic contexts (exon/intron/intergenic) and across methylation 
contexts (CpG/CHG/CHH).
"""
import collections
import csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
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

context_tally = pd.DataFrame()
for cxx in ['CpG', 'CHG', 'CHH']:
    data = pd.read_table(
        f'../../all.bona_fide_meth_pos.annot.c20.{cxx.lower()}.cov',
        header=None,
        usecols=[3, 6, 10],
        names=['meth_pct', 'gene', 'ei'])
    
    data['temp'] = data['gene'] + data['ei']
    data['temp'] = data['temp'].fillna(value='intergenic')
    data['genomic_feature'] = data['temp'].apply(define_genomic_feature)
    
    data = data.drop('temp', axis=1)
    data = data.drop('gene', axis=1)
    data = data.drop('ei', axis=1)
    
    # add counter to context_tally DataFrame
    context_tally[cxx] = pd.Series(collections.Counter(data['genomic_feature']))

# do a bit of magic to plot stacked bar charts. exon is at the bottom, intron
# in the middle (i.e. its height is exon + intron), and intergenic at the top 
# of the stack (i.e. its height is exon + intron + intergenic)
context_tally = context_tally.transpose()
context_tally['intron_stack'] = context_tally['intron'] + context_tally['exon']
context_tally['intergenic_stack'] = context_tally['intergenic'] + \
                                    context_tally['intron_stack']

# seaborn
sns.set_style('whitegrid')
f, ax = plt.subplots(figsize=(4, 4))

sns.barplot(x=context_tally.index, y=context_tally['intergenic_stack'],
            label='Intergenic', color='#f4a582')
sns.barplot(x=context_tally.index, y=context_tally['intron_stack'],
            label='Intron', color='#67a9cf')
sns.barplot(x=context_tally.index, y=context_tally['exon'],
            label='Exon', color='#2166ac')

plt.legend(loc=9, ncol=1)

ax.set_xlabel('Sequence context')
ax.set_ylabel('Number of methylated cytosines')

# rescale ticks from xx,000 to xxk
ticks = [f'{int(x)}k' if x else '0' for x in ax.get_yticks() / 1000]
ax.set_yticklabels(ticks)

sns.despine(left=True, bottom=True, offset=10, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'cytosine_vs_genomic_context.pdf'
fig.savefig(output_filename, bbox_inches='tight')
