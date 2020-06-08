==================================
Filtering for methylated positions
==================================

This is gonna be a bit confusing, apologies in advance.

Underlying observations
-----------------------
We noticed early on that low mapping rates + considerable non-CpG methylation means that

1. not many positions were commonly methylated across all analysed samples.
2. coverage has an oversized effect on methylation levels.

For (1), we attempted to re-sequence the WGBS libraries at greater depth, but the mapping and duplication rates were high. We tried to troubleshoot this (see Supplementary Discussion), but the low numbers do not seem to be errors on our part.

On the other hand, (2) can be countered by subsampling everything down to the least-covered library (similar to how microbial studies normalise between libraries of different depths).

Deviations from typical pipeline
--------------------------------
So far, the methylation pipeline is

1. do the usual Bismark-y things: ``bismark``, ``deduplicate_bismark``, ``bismark_methylation_extractor``. This produces the files in ``merged_unfilt_covs/`` (not provided, not used in any analytical scripts).

2. feed all ``*.cov`` files into the script ``subsample_bismark_covs.py``. This produces the files in ``subsampled_covs/`` (also not provided, as not used in downstream analyses).

3. at this point, the filtering forks to produce (files in both forks are provided this folder, as they are used in downstream analyses)

   a. well-covered positions that are present in all 15 samples. Files have ``.filt`` somewhere in their filenames, and contain 6,748 lines/positions. Unannotated files are in ``subsampled_filt_covs``, while annotated ones are in ``./``, as ``*.subsamp.filt.annot.cov``. Line counts::
   
       lyijin@mirrordin:/mnt/c/Yi/github/smic_dna_meth/filter_meth_pos_CX/subsampled_filt_covs$ for a in *.gz; do printf "${a}:    "; zcat ${a} | wc -l; done
       26A.subsamp.filt.cov.gz:    6748
       26B.subsamp.filt.cov.gz:    6748
       26C.subsamp.filt.cov.gz:    6748
       26D.subsamp.filt.cov.gz:    6748
       26E.subsamp.filt.cov.gz:    6748
       26F.subsamp.filt.cov.gz:    6748
       29D.subsamp.filt.cov.gz:    6748
       29E.subsamp.filt.cov.gz:    6748
       29F.subsamp.filt.cov.gz:    6748
       32A.subsamp.filt.cov.gz:    6748
       32B.subsamp.filt.cov.gz:    6748
       32C.subsamp.filt.cov.gz:    6748
       32D.subsamp.filt.cov.gz:    6748
       32E.subsamp.filt.cov.gz:    6748
       32F.subsamp.filt.cov.gz:    6748
       
       lyijin@mirrordin:/mnt/c/Yi/github/smic_dna_meth/filter_meth_pos_CX$ for a in *subsamp.filt.annot.cov.gz; do printf "${a}:    "; zcat ${a} | wc -l; done
       26A.subsamp.filt.annot.cov.gz:    6748
       26B.subsamp.filt.annot.cov.gz:    6748
       26C.subsamp.filt.annot.cov.gz:    6748
       26D.subsamp.filt.annot.cov.gz:    6748
       26E.subsamp.filt.annot.cov.gz:    6748
       26F.subsamp.filt.annot.cov.gz:    6748
       29D.subsamp.filt.annot.cov.gz:    6748
       29E.subsamp.filt.annot.cov.gz:    6748
       29F.subsamp.filt.annot.cov.gz:    6748
       32A.subsamp.filt.annot.cov.gz:    6748
       32B.subsamp.filt.annot.cov.gz:    6748
       32C.subsamp.filt.annot.cov.gz:    6748
       32D.subsamp.filt.annot.cov.gz:    6748
       32E.subsamp.filt.annot.cov.gz:    6748
       32F.subsamp.filt.annot.cov.gz:    6748
      
      These files are primarily used in analysis that cannot tolerate NAs: PCAs, correlation analysis.
      
   b. statistically bona-fide positions for each sample. Methylation in that position is most likely "real", with binomial *p* < 0.05. Files are at ``subsampled_bona_fide/``. They are further combined to create an "all" file, and a coverage filter of 20 picks out positions that are better covered (this is why the "all" files have ".c20" in them). The "all" file is further split by C context, using ``split_bismark_cov_by_context.py``. Line counts::
   
       lyijin@mirrordin:/mnt/c/Yi/github/smic_dna_meth/filter_meth_pos_CX/subsampled_bona_fide$ for a in *.gz; do printf "${a}:    "; zcat ${a} | wc -l; done
       26A.subsamp.bona_fide.cov.gz:    1297783
       26B.subsamp.bona_fide.cov.gz:    1531685
       26C.subsamp.bona_fide.cov.gz:    1556262
       26D.subsamp.bona_fide.cov.gz:    1271749
       26E.subsamp.bona_fide.cov.gz:    1414944
       26F.subsamp.bona_fide.cov.gz:    1245523
       29D.subsamp.bona_fide.cov.gz:    1140256
       29E.subsamp.bona_fide.cov.gz:    979973
       29F.subsamp.bona_fide.cov.gz:    1300358
       32A.subsamp.bona_fide.cov.gz:    1236713
       32B.subsamp.bona_fide.cov.gz:    1428229
       32C.subsamp.bona_fide.cov.gz:    1173380
       32D.subsamp.bona_fide.cov.gz:    936243
       32E.subsamp.bona_fide.cov.gz:    834840
       32F.subsamp.bona_fide.cov.gz:    822445
       
       lyijin@mirrordin:/mnt/c/Yi/github/smic_dna_meth/filter_meth_pos_CX$ for a in all.bona_fide*.cov.gz; do printf "${a}:    "; zcat ${a} | wc -l; done
       all.bona_fide_meth_pos.annot.c20.chg.cov.gz:    86666
       all.bona_fide_meth_pos.annot.c20.chh.cov.gz:    213296
       all.bona_fide_meth_pos.annot.c20.cov.gz:    478717
       all.bona_fide_meth_pos.annot.c20.cpg.cov.gz:    178755
     
      These files are used in all other analyses e.g. methylation of repeats, methylation vs. expression, methylation vs. spurious transcription, methylation vs. transcriptional noise...

Original file sizes
-------------------
This is why all files couldn't be uploaded::

    lyijin@mirrordin:/mnt/c/Yi/KAUST_postdoc/lithium/symb/dna_meth/filter_meth_pos_CX$ du -d 1
    5.3G    ./merged_unfilt_covs
    411M    ./subsampled_bona_fide   (*)
    1.7G    ./subsampled_covs
    948K    ./subsampled_filt_covs   (*)
    8.4G    .
    8.4G    total

Lines with asterisks are provided in this subfolder.
