==================================================
Scripts to analyse differentially methylated genes
==================================================

Scripts in all subfolders can be executed in any order.

Brief description of subfolder contents
---------------------------------------
1. ``coverage_plots/`` visualises the effect of the subsampling technique, aimed to reduce effect of coverage on methylation %. Produces Supplementary Fig. 2.

2. ``expt_setup/`` has the artisanal illustration of the experiment--Supplementary Fig. 1.

3. ``genomic_context/`` illustrates where the methylated positions lie within the genome. Are they in intergenic regions? Genic regions? For genic--introns, or exons? etc. Fig. 1 is produced by code here.

4. ``overlap_rep_elements/`` checks whether methylated positions are more commonly found in repeat elements. While the analysis relies on Python scripts, the results were unified in Excel (!), plotted in Excel (!!!) and prettified in Illustrator. Makes Figs. 2a, 2b, 2d.

5. ``pca_and_heatmap/`` does what it says on the tin. Produces Fig. 3.
