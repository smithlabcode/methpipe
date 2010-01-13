#!/usr/bin/env python

# jackpot-del
# Jin Park <jinhpark@usc.edu> 2010
# This program remove duplicated reads mapped to the same location. It there are multiple
# reads mapped to the same location, it choose the one with the fewest errors.
# It assumes the input bed file is sorted

import sys, string, random

prevName = '' # empty string
same = [] # empty list
val = 10000000 # Just a large value

## Iterate over the lines in a file (does not read all in at once)
for i in open(sys.argv[1]):
    (chrom, start, end, readName, score, strand) = i.split()
    # Unique identifyer for mapping location: (end is redundant as all reads have same length)
    name = chrom + start
    score = float(score)
    if prevName != '' and prevName != name: #--when new name started
        #  among reads at same location (and with same score) randomly choose one
        print random.choice(same),
        val = 10000000
    #  among reads mapping to the same location, select for lower score
    if score < val: #--initially comming to here, and also when new score is LT old one with same name
        same = [i]  #--list with only one element - current one
        val = score #--update value
    ###else: same.append(i) #--this is BUG!! --append also bigger scored ones
    elif score == val: #--if same score with same name
        same.append(i)  #--append current one to the list, for random choice later
    #--otherwise (bigger score case), disregard
    prevName = name
    
if prevName != '':
    print random.choice(same),
