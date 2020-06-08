=============================
Raw data used in this project
=============================

Most files are too large to be uploaded, but here's a list of files used::

  lyijin@mirrordin:/mnt/c/Yi/KAUST_postdoc/lithium/symb/dna_meth/raw_data$ du -d 1
  25G     ./bis_genomes
  99M     ./repeatmasker
  25G     .
  25G     total
  
  lyijin@mirrordin:/mnt/c/Yi/KAUST_postdoc/lithium/symb/dna_meth/raw_data$ ls -l
  total 400M
  drwxrwxrwx 1 lyijin lyijin  512 May 26  2019 bis_genomes/
  drwxrwxrwx 1 lyijin lyijin  512 Apr 26 21:37 repeatmasker/
  -rwxrwxrwx 1 lyijin lyijin  44M Jan 18  2017 Smic.genome.annotation.gff3.gz*
  -rwxrwxrwx 1 lyijin lyijin 509K Oct 12  2016 Smic.genome.scaffold.final.fa.fai*
  -rwxrwxrwx 1 lyijin lyijin 199M Sep  7  2014 Smic.genome.scaffold.final.fa.gz*
  -rwxrwxrwx 1 lyijin lyijin 157M Sep  7  2014 Smic.genome.scaffold.final.masked.fa.gz*

``Smic.*.gz`` files can be downloaded from http://smic.reefgenomics.org/download/

``bis_genomes/`` contains the genomes that WGBS were mapped against. See ``../mapping_against_non_symb`` for a list of non-Smic genomes that we mapped against, in order to understand why mapping rates were low.

``repeatmasker/`` contains repeat annotations for the **S. microadriaticum** genome (provided here; courtesy Xin Wang).
