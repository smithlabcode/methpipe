#! /usr/bin/python

# intersection.py
# Song Qiang <qiang.song@usc.edu> 2010
# this program produce the intersection of genomic regions defined in
# two BED files. It assumes the regions the file are sorted and non-overlapping.
# it outputs the result to stdout

import os, sys
from optparse import OptionParser

# parse options
opt_parser = OptionParser()
opt_parser.add_option("-o", "--output", dest = "output_file_name",  default = "", \
					  help = "File name of output file name", metavar = "FILE")
(options, args) = opt_parser.parse_args()

class BEDEntry:
	def __init__(self, bed_line = ""):
		if bed_line:
			bed_line = bed_line.strip()
			parts = bed_line.split("\t")
			self.chrom = parts[0]
			self.start = int(parts[1])
			self.end = int(parts[2])
			self.other = "\t".join(parts[3:])
			self.bed_line = bed_line
		else:
			self.chrom = ""
			self.start = 0
			self.end = 0
			self.other = ""
			self.bed_line = bed_line

def intersect(bed_entry_1, bed_entry_2):
	if not bed_entry_1.chrom == bed_entry_2.chrom:
		return BEDEntry()
	chrom = bed_entry_1.chrom
	start = max(bed_entry_1.start, bed_entry_2.start)
	end = min(bed_entry_1.end, bed_entry_2.end)
	if start >= end:
		return BEDEntry()
	else:
		return BEDEntry("%s\t%d\t%d\t%s\t%s" % (chrom, \
												start, \
												end, \
												bed_entry_1.bed_line, \
												bed_entry_2.bed_line))
		
# load all files into memory: stupid but easy	
bed_entry_list_1 = list()
for line in open(args[0]):
	bed_entry_list_1.append( BEDEntry(line) )

bed_entry_list_2 = list()
for line in open(args[1]):
	bed_entry_list_2.append( BEDEntry(line) )

# do the intersection
len_1  = len(bed_entry_list_1)
len_2  = len(bed_entry_list_2)
i1 = 0
i2 = 0
while(i1 < len_1 and i2 < len_2):
	entry_1 = bed_entry_list_1[i1]
	entry_2 = bed_entry_list_2[i2]
	inter_bed_entry = intersect(entry_1, \
								entry_2)
	if inter_bed_entry.chrom:
		sys.stdout.write(inter_bed_entry.bed_line + "\n")
	if ( entry_1.chrom < entry_2.chrom ) \
	   or ( entry_1.chrom == entry_2.chrom and entry_1.end < entry_2.end):
		i1 += 1
	else:
		i2 += 1
		
