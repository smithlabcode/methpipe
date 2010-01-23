#!/usr/bin/env python

# methcounts-parser 
# Song Qiang <qiang.song@usc.edu> 2010
# this program parse the output from methcounts, outputs two files:
# a BedGraph file containing reads coverage for each CpG and 
# another file containing the methylation probability.


import sys, string

coverage_output_file = open(sys.argv[2], "w")
meth_prob_output_file = open(sys.argv[3], "w")

for line in open(sys.argv[1]):
	(chrom, start, end, readName, score, strand) = line.split()
	coverage = int(readName.split(":")[1])
	coverage_output_file.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, coverage) )
	meth_prob_output_file.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, score) )

coverage_output_file.close()
meth_prob_output_file.close()
	
