===========================================
Mapping of WGBS reads against other genomes
===========================================

In Supplementary Discussion, we mused about the low mapping rates when the WGBS reads were mapped against the *S. microadriaticum* reads. Perhaps the reads were contaminated by other unknown organisms? Are we actually dealing with *S. microadriaticum* here? Would be facepalm-worthy if we weren't...

These subfolders contain the mapping logs of 26A vs. multiple genomes.

Bacterial, to check for contamination:

1. **acine**: gi|613736623|gb|CP007549.1| Acinetobacter baumannii AC12, complete genome

2. **ecoli**: Escherichia_coli_K-12_complete_genome

3. **salmo**: gi|478427060|ref|NZ_KB731355.1| Salmonella enterica subsp. enterica serovar Dublin str. UC16 genomic scaffold Scaffold138, whole genome shotgun sequence

4. **strep**: gi|633274572|gb|CP007628.1| Streptococcus sp. VT 162, complete genome

Technical, to check for human error:

1. **phix**: gi|9626372|ref|NC_001422.1| Enterobacteria phage phiX174 sensu lato, complete genome (did we accidentally flood the whole run with PhiX?)

2. **smic_genome_masked**: masked *S. microadriaticum* genome (maybe we're getting massive amounts of reads mapping to repeat elements?)

3. **smin**: *S. minutum* (now known as *Breviolum minutum*) genome (maybe there was a confusion in *Symbiodinium* strain?)

The reads were mapped in paired-end mode and single-end mode, in case the mapper (bowtie) had technical errors... yes we were grasping at straws.

Results
-------

Too lazy to create a formatted table, here's the grep output of the logs:: 

  lyijin@mirrordin:/mnt/c/Yi/github/smic_dna_meth/mapping_against_non_symb$ grep 'Mapping efficiency' */*.txt
  acine_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:      1.0%
  acine_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      1.2%
  acine_map/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      1.1%
  ecoli_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:      0.0%
  ecoli_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      0.1%
  ecoli_map/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      0.0%
  phix_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:       1.0%
  phix_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:       1.1%
  phix_map/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:       1.1%
  salmo_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:      1.0%
  salmo_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      1.1%
  salmo_map/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      1.1%
  smic_genome_masked/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:     17.9%
  smic_genome_masked/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:     26.7%
  smic_genome_masked/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:     25.9%
  smin_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:       3.0%
  smin_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:       4.1%
  smin_map/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:       3.6%
  strep_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:      1.0%
  strep_map/26-Alpha_S1_L001_R1_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      1.4%
  strep_map/26-Alpha_S1_L001_R2_001.trim.fastq_bismark_bt2_SE_report.txt:Mapping efficiency:      1.1%


Observations
------------

1. PE vs. SE: not much difference. SE of course maps slightly better
2. Bacterial genomes: low mapping % throughout
3. smic vs. smin: we're most likely dealing with *S. microadriaticum* here
4. smic vs. smic_masked: similar %s are obtained from mapping against unmasked genome, i.e. repeat elements not a big deal here

Conclusions
-----------

We're most likely doing things correctly (phew). Jury's still out on why mapping rates are low.
