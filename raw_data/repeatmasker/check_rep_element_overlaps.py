#!/usr/bin/env python3

"""
> check_rep_element_overlaps.py <

How often do repeat elements overlap (i.e. same genomic loci has two different
repeat categories assigned to it)?

Answer: 145,853 of 2,039,087 = 7.2%
"""
import gzip
import re

counter = 0
with gzip.open('smic.fa.out.gz', 'rt') as f:
    # skip first three lines
    f.readline()
    f.readline()
    f.readline()
    prev_line = re.split(r'\s+', f.readline().strip())
    
    for line in f:
        line = re.split(r'\s+', line.strip())
        
        prev_scaf = prev_line[4]
        prev_start = prev_line[5]
        prev_end = prev_line[6]
        prev_element = prev_line[10]
        
        curr_scaf = line[4]
        curr_start = line[5]
        curr_end = line[6]
        curr_element = line[10]
        
        prev_element = prev_element.split('/')[0]
        curr_element = curr_element.split('/')[0]
        
        if int(curr_start) < int(prev_end) and prev_scaf == curr_scaf and \
                prev_element != curr_element:
            counter += 1
            print (f'--{counter}--')
            print ('\t'.join(prev_line))
            print ('\t'.join(line))
        
        prev_line = line
