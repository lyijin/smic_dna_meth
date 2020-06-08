#!/usr/bin/env python3

"""
> plot_dotplot.spurious.meth_vs_unmeth.py <

Takes in output file from compile_overall_coverage.py, and plots the relative
expression of exon 1 vs. subsequent exons, up to NUM_EXONS.

Trends of unmethylated genes and densely methylated genes (# meth pos >= 10)
are plotted, to see whether methylation helps suppress spurious transcription.
"""
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns

NUM_EXONS = 6

# read in median methylation levels
median_meth = pd.read_table('../../descriptive_wgbs/bias_density_medians/all.bdm.c.tsv')

# get lists of meth genes
unmeth_genes = median_meth[median_meth['meth_pos'] == 0]['gene'].tolist()
meth_genes = median_meth[median_meth['meth_pos'] > 0]['gene'].tolist()
densely_meth_genes = median_meth[median_meth['meth_pos'] >= 10]['gene'].tolist()

# read in mean expression coverages
mean_cov = pd.read_table('mean_cov.all.compiled.tsv',
                         usecols=range(NUM_EXONS+1),
                         header=None, 
                         names=['gene'] + \
                               [str(x) for x in range(1, NUM_EXONS+1)])

# colours for plot: create empty dict
colours = {}

unmeth_mean_cov = mean_cov[mean_cov['gene'].isin(unmeth_genes)]
n = len(unmeth_mean_cov)
unmeth_mean_cov = unmeth_mean_cov.drop(['gene'], axis=1)
um_label = 'Unmethylated (n={:,})'.format(n)
unmeth_mean_cov = pd.melt(unmeth_mean_cov, value_name='relative_expr', var_name='exon')
unmeth_mean_cov['meth_status'] = um_label
colours[um_label] = '#fee0d2'

meth_mean_cov = mean_cov[mean_cov['gene'].isin(meth_genes)]
n = len(meth_mean_cov)
meth_mean_cov = meth_mean_cov.drop(['gene'], axis=1)
m_label = 'Methylated (n={:,})'.format(n)
meth_mean_cov = pd.melt(meth_mean_cov, value_name='relative_expr', var_name='exon')
meth_mean_cov['meth_status'] = m_label
colours[m_label] = '#fc9272'

densely_meth_mean_cov = mean_cov[mean_cov['gene'].isin(densely_meth_genes)]
n = len(densely_meth_mean_cov)
densely_meth_mean_cov = densely_meth_mean_cov.drop(['gene'], axis=1)
dm_label = 'Densely methylated (n={:,})'.format(n)
densely_meth_mean_cov = pd.melt(densely_meth_mean_cov, value_name='relative_expr', var_name='exon')
densely_meth_mean_cov['meth_status'] = dm_label
colours[dm_label] = '#de2d26'

mean_cov = pd.concat(
    [unmeth_mean_cov, meth_mean_cov, densely_meth_mean_cov], ignore_index=True)

# seaborn
sns.set_style('white')
sns.set_style('ticks')
fig, ax= plt.subplots(figsize=(6, 4))

sns.pointplot(x='exon', y='relative_expr', hue='meth_status', label='',
              data=mean_cov, palette=colours,
              ci=68, n_boot=10000, dodge=True)

ax.set_xlabel('Exon')
ax.set_ylabel('ln (relative expression vs. exon 1)')
ax.legend(loc='lower right')
sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'spurious.meth_vs_unmeth.pdf'
fig.savefig(output_filename, bbox_inches='tight')

# perform statistical tests (t-tests)
for n in range(1, NUM_EXONS+1):
    um = unmeth_mean_cov[unmeth_mean_cov['exon'] == str(n)]['relative_expr']
    m = meth_mean_cov[meth_mean_cov['exon'] == str(n)]['relative_expr']
    dm = densely_meth_mean_cov[densely_meth_mean_cov['exon'] == str(n)]['relative_expr']
    
    print ('Exon', n)
    print ('t test p for um vs. m: ', scipy.stats.ttest_ind(um, m))
    print ('t test p for um vs. dm: ', scipy.stats.ttest_ind(um, dm))
    print ()
