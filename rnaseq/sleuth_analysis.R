# > sleuth_analysis.R <
#
# Carries out the sleuth bit of the differential expression analysis. Performs
# two pairwise analyses: 26Cv29C and 26Cv32C.
#
# Contains bits of code adapted from
# https://gist.github.com/jaquol/03f41f57dc6b0eacef101e9920f24d78
#
# The bit that performs LRT and WT, then merges the output tables, was a
# suggestion from Guoxin Cui.

library(sleuth)

# hardcode folders that contain kallisto results
base_dir <- '~/kaust/symb/dna_meth/rnaseq/kallisto'
sample_id <- dir(file.path(base_dir, 'results'))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, 'results', id))

# read experimental setup
s2c <- read.table(file.path(base_dir, 'expt_setup.tsv'), 
                  header=TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path=kal_dirs)

# check that s2c is set up correctly
print(s2c)

# start by printing TPMs normalised across all replicates: this table is
# used in plotting PCAs.

# load data
so_all <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE)

# print normalised TPMs
tpms <- kallisto_table(so_all)[, c('target_id', 'sample', 'tpm')]
tpms <- reshape2::dcast(tpms, target_id~sample, value.var='tpm')
write.table(tpms, 'normalised_abundances.all.tsv', sep='\t',
            quote=FALSE, row.names=FALSE)

# create a loop to fork analysis for 26Cv29C and 26Cv32C
for (temp in c('29C', '32C')) {
    # subselect control and temperature of interest
    s2c_subset <- dplyr::filter(s2c, condition %in% c('26C', temp))
    
    # load data
    so <- sleuth_prep(s2c_subset, extra_bootstrap_summary=TRUE)
    
    # print normalised TPMs
    tpms <- kallisto_table(so)[, c('target_id', 'sample', 'tpm')]
    tpms <- reshape2::dcast(tpms, target_id~sample, value.var='tpm')
    write.table(tpms, paste0('normalised_abundances.26Cv', temp, '.tsv'), 
                sep='\t', quote=FALSE, row.names=FALSE)
    
    # estimate parameters for the response error measurement model that is
    # dependent on 'condition'
    so <- sleuth_fit(so, ~condition, 'full')
    
    # run another model where gene expression is independent of any factors
    # (i.e. null expectations)
    so <- sleuth_fit(so, ~1, 'reduced')
    
    # calculate likelihood ratios that gene expression is dependent
    # on 'condition'
    so <- sleuth_lrt(so, 'reduced', 'full')
    
    # the current drawback with sleuth_lrt is that it does not provide the
    # direction of differential expression (over/underexpression). while LRT
    # is technically better than wald's test, wald's test is run here to
    # calculate beta (approx equals ln FC)
    so <- sleuth_wt(so, paste0('condition', temp))
    
    # print DEG table
    res_lrt <- sleuth_results(so, 'reduced:full', test_type='lrt')
    res_wt <- sleuth_results(so, paste0('condition', temp))
    res <- merge(res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], 
                 res_lrt, on='target_id', sort=TRUE)
    write.table(res, paste0('sleuth_results.26Cv', temp, '.tsv'),
                sep="\t", quote=FALSE, row.names=FALSE)
}