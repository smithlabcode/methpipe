#!/usr/bin/env python

# methseg-parser 
# Song Qiang <qiang.song@usc.edu> 2010
# this program parses the output from methdiff, selects domains satisfying
# certain criteria. It uses stdin and stdout.

import sys, string
from optparse import OptionParser

opt_parser = OptionParser()
opt_parser.add_option("-c", "--min_score", dest = "min_posterier_score", default = 20, type="float",\
					  help="Minimum score for a selected domain")
opt_parser.add_option("-b", "--bed_file", help="output bed file corresponding to original BED format",\
					  default = False, action="store_true", dest="output_bed")
(options, args) = opt_parser.parse_args()

# min_posterier_score = 20				# almost equivalent to domains containg more than 15 CpGs 

for line in sys.stdin:
	(chrom, start, end, readName, score, strand) = line.split()
	score = float(score)
	if score > options.min_posterier_score:
		if options.output_bed:
			sys.stdout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, readName)) 
		else:
			sys.stdout.write(line)

	
