#!/usr/bin/env python

# methseg-parser 
# Song Qiang <qiang.song@usc.edu> 2010
# this program parses the output from methdiff, selects domains satisfying
# certain criteria. It uses stdin and stdout.

import sys, string

min_posterier_score = 15				# almost equivalent to domains containg more than 15 CpGs 

for line in open(sys.stdin):
	(chrom, start, end, readName, score, strand) = line.split()
	score = float(score)
	if score > min_posterier_score:
		print line

	
