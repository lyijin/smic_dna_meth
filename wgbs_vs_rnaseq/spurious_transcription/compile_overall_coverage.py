#!/usr/bin/env python3

"""
> compile_overall_coverage.py <

Given the per-position coverages from the desired input files, assign equal
weightage to each file, and compile overall trends from all datasets.

Implements several filters:
I.   Gene has to have >= NUM_EXONS exon
II.  Genes must be expressed in all tested exons.
III. Gene has to have mean RPKM > 0.5 overall
IV.  Gene has to pass all criteria above in all files.
"""
import csv

import numpy as np
import pandas as pd

NUM_EXONS = 6

exon_covs = {}
norm_cov = {}

# read input files
data = pd.read_table(open('all.ei_norm.merged.tsv'))

# drop intron lines
data = data[data['feature'] == 'Exon']

# drop genes with fewer than NUM_EXONS (criteria I)
data['num_exons'] = data['feature_no'] - data['feature_no_rev'] - 1
data = data[data['num_exons'] >= NUM_EXONS]

# check whether gene has sufficient overall coverage (criteria III)
# passing_genes = []
# for gene in set(data['gene']):
    # if all(data[data['gene'] == gene].mean().filter(regex='26|29|32') > 0.5):
        # passing_genes.append(gene)

# data = data[data['gene'].isin(passing_genes)]

# check whether genes are expressed in all tested exons (criteria II)
data = data[data['feature_no'] <= NUM_EXONS]
passing_genes = []
for gene in set(data['gene']):
    if data[data['gene'] == gene].filter(regex='26|29|32').min().min() > 0:
        passing_genes.append(gene)

data = data[data['gene'].isin(passing_genes)]

data = data.sort_values(by=['gene', 'feature_no'])
# create sub_data containing expression values for exon 1
sub_data = data[data['feature_no'] == 1]

# hack: repeat this dataframe to have NUM_EXONS line per gene, so that
# normalisation can be done using data[samp] / sub_data[samp]
sub_data = sub_data.append([sub_data] * (NUM_EXONS - 1), ignore_index=True)
sub_data = sub_data.sort_values(by='gene')

for samp in ['26D', '26E', '26F', '29D', '29E', '29F', '32D', '32E', '32F']:
    # calculate ln (FC of first exon expression)
    data[f'{samp}_norm'] = np.log(np.array(data[samp]) / np.array(sub_data[samp]))

# calculate mean ln FC
data['mean_norm'] = data.filter(regex='_norm').mean(axis=1)

# retain only important columns for subsequent printing to stdout
data = data[['gene', 'feature_no', 'mean_norm']]
data = data.pivot_table(values='mean_norm', index='gene', columns='feature_no')
data = data.applymap(lambda x: round(x, 4))

data.to_csv('mean_cov.all.compiled.tsv', sep='\t', header=False)
