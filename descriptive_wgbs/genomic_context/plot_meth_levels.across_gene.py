#!/usr/bin/env python3

"""
> plot_meth_levels.across_gene.py <

Plots methylation levels across the gene, from -4 kb upstream, across a
normalised gene (following mean gene length of 12,898 bp in S. microadriaticum),
and ending with a 4 kb downstream region.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

SMIC_GENE_LENGTH = 12898
START = -4000
END = SMIC_GENE_LENGTH + 4000

# read the cov files
context = ['cpg', 'chg', 'chh']
all_data = {}
combined = {}
for c in context:
    all_data[c] = \
        pd.read_table(f'../../all.bona_fide_meth_pos.annot.c20.{c}.cov',
                      header=None,
                      names=['scaf', 'start_pos', 'end_pos', 'meth_pct', 
                             'meth', 'unmeth', 'gene', 'gene_relpos',
                             'gene_5p', 'gene_3p', 'ei', 'ei_relpos', 
                             'ei_5p', 'ei_3p'])

    data = all_data[c]
    
    # remove lines with 'no_info'
    data = data[data['gene'] != 'no_info'][data['ei'] != 'no_info']
    
    # grab genic context
    tmp = data[data['gene'].str.match('SmicGene')]
    genic = pd.DataFrame()
    genic['meth_pct'] = tmp['meth_pct']
    genic['recalibrated_x'] = tmp['gene_relpos'] * SMIC_GENE_LENGTH
    
    # grab upstream region
    upstreams = []
    tmp = data[data['ei'].str.match('upstream_crick\|')]
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = -tmp['gene_5p'] - 1
    upstreams.append(tmp2)
    
    tmp = data[data['ei_relpos'].str.match('downstream_watson\|')]
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = tmp['gene_3p']
    upstreams.append(tmp2)
    
    upstreams = pd.concat(upstreams)
    upstreams = upstreams[upstreams['recalibrated_x'] >= START]
    
    # grab downstream region
    downstreams = []
    tmp = data[data['ei'].str.match('upstream_watson\|')]
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = SMIC_GENE_LENGTH + tmp['gene_5p'] + 1
    downstreams.append(tmp2)
    
    tmp = data[data['ei_relpos'].str.match('downstream_crick\|')]
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = SMIC_GENE_LENGTH - tmp['gene_3p']
    downstreams.append(tmp2)
    
    downstreams = pd.concat(downstreams)
    downstreams = downstreams[downstreams['recalibrated_x'] <= END]
    
    # combine everything together
    combined[c] = pd.concat([upstreams, genic, downstreams])

    # testing: randomly pick 2000 rows to generate approximate plot
    #combined[c] = combined[c].sample(2000)

# plot the heatmap!
sns.set_style('white')
sns.set_style('ticks')

fig, ax = plt.subplots(figsize=(4, 2))
colour = {'cpg': '#fed98e', 'chg': '#fe9929', 'chh': '#d95f0e'}

for c in context:
    sns.regplot(x='recalibrated_x', y='meth_pct', data=combined[c],
                scatter=False, lowess=True,
                label=c, color=colour[c])

plt.xlabel('Genomic context')
plt.ylabel('Methylation level (%)')

plt.xlim(START, END)
#plt.ylim(0, 100)
xtick_labels = ['4 kb upstream', 'normalised gene', '4 kb downstream']
xt = [-2000, 0.5 * SMIC_GENE_LENGTH, END - 2000]
plt.xticks(xt, xtick_labels, ha='center')

plt.legend(loc='upper right')

# mark genic region out with solid lines
plt.axvline(x=0, lw=1, color='k')
plt.axvline(x=SMIC_GENE_LENGTH, lw=1, color='k')

sns.despine(offset=10, bottom=True, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'meth_levels.across_gene.pdf'
fig.savefig(output_filename, bbox_inches='tight')
