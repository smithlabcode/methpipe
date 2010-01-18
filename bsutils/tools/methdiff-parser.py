#!/usr/bin/env python

# methdiff-parser 
# Song Qiang <qiang.song@usc.edu> 2010
# this program parse the output from methdiff, outputs a bedGraph file

import sys, string

output_file = open(sys.argv[2], "w")

for line in open(sys.argv[1]):
	(chrom, start, end, readName, score, strand) = line.split()
	score = str( float(score) - 0.5)
	output_file.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, score) )

output_file.close()
	
