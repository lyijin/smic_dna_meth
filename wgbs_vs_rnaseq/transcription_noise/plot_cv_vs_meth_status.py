#!/usr/bin/env python3

"""
> plot_cv_vs_meth_status.py <

Does methylation have an effect on transcriptional noise (measured as inverse
coeffficient of variation here)?
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

# read data
data = pd.read_table('../../rnaseq/normalised_abundances.all.tsv', index_col=0)

# drop rows containing a single 0
data = data.loc[(data!=0).all(axis=1)]

# apply log10 (tpm+1)
data = data.applymap(lambda x: np.log10(x + 1))

# calculate mean and sample stdev
data2 = pd.DataFrame()
data2['mean'] = data.mean(axis=1)
data2['stdev'] = data.apply(lambda x: np.std(x, ddof=1), axis=1)

# cv = sample stdev / mean; inv_cv = mean / sample stdev
data['inv_cv'] = data2['mean']/data2['stdev']

# pull in expression means
data['mean'] = data2['mean']

# figure out which genes are methylated and which aren't
bdm_data = pd.read_table('../../descriptive_wgbs/bias_density_medians/all.bdm.c.tsv')
unmeth_genes = bdm_data[bdm_data['meth_pos'] == 0]['gene'].tolist()
meth_genes = bdm_data[bdm_data['meth_pos'] > 0]['gene'].tolist()

# partition data into meth_ and unmeth_data
unmeth_data = data[data.index.isin(unmeth_genes)]
meth_data = data[data.index.isin(meth_genes)]

# model as inverse power graph: y = k/x
# note: there's no +c in the equation, as a gene with 0 expression has,
#       by definition, 0 noise.
x_meth = meth_data['mean']
y_meth = meth_data['inv_cv']
x_meth_mean = np.mean(x_meth)
# X = sm.add_constant(X)
model = sm.OLS(y_meth, x_meth)
results = model.fit()
rsq_meth = results.rsquared
rmaxp_meth = '{:.0e}'.format(max([x for x in results.pvalues]))
if rmaxp_meth == '0e+00': rmaxp_meth = '1e-300'
k_meth = results.params[0]

# do the same for unmethylated genes
x_unmeth = unmeth_data['mean']
y_unmeth = unmeth_data['inv_cv']
x_unmeth_mean = np.mean(x_unmeth)
# X = sm.add_constant(X)
model = sm.OLS(y_unmeth, x_unmeth)
results = model.fit()
rsq_unmeth = results.rsquared
rmaxp_unmeth = '{:.0e}'.format(max([x for x in results.pvalues]))
if rmaxp_unmeth == '0e+00': rmaxp_unmeth = '1e-300'
k_unmeth = results.params[0]

# plot using vanilla matplotlib. mpl v2 is pretty!
# note: seaborn's regplot plots using y = mx + c, can't find way to suppress
#       the +c portion.
sns.set_style('white')
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(6, 4))

# plot the dots
su = plt.scatter(unmeth_data['mean'], unmeth_data['inv_cv'], s=5,
                 c='#fee0d2', marker='o', alpha=0.15)
sm = plt.scatter(meth_data['mean'], meth_data['inv_cv'], s=5,
                 c='#fc9272', marker='^', alpha=0.15)

# plot the lines
plt.plot([0,3], [0,3*k_unmeth], linewidth=2, color='#fee0d2')
plt.plot([0,3], [0,3*k_meth], linewidth=2, color='#fc9272')

# plot the text
ax.text(2.5, 2.5 * k_meth + 5.5,
        '$r^2 =${:.3f}\n$p <${}'.format(rsq_meth, rmaxp_meth),
        fontsize=10)
ax.text(2.5, 2.5 * k_unmeth - 5.5, 
        '$r^2 =${:.3f}\n$p <${}'.format(rsq_unmeth, rmaxp_unmeth),
        fontsize=10)

# prettifiers
ax.set_xlim(0, 3)
ax.set_ylim(0, 60)
sns.despine(offset=10, trim=True)
ax.set_xlabel('$\log_{10}$ (mean expression + 1)')
ax.set_ylabel('Coefficient of variation$^{-1}$')

plt.legend((sm, su),
           (f'Methylated genes (n={len(meth_data)})',
            f'Unmethylated genes (n={len(unmeth_data)})'), 
           loc='upper left', ncol=1, scatterpoints=3, frameon=True)
legend = ax.get_legend()
legend.legendHandles[0].set_alpha(1)
legend.legendHandles[1].set_alpha(1)

plt.tight_layout()

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'cv_vs_meth_status.pdf'
fig.savefig(output_filename, bbox_inches='tight')
