#!/usr/bin/env python

# bed2html 
# Song Qiang <qiang.song@usc.edu> 2010
# this program convert a BED file to html file
# it is used to display DMR files and / or hypo-methylated regions 
# in html format

import sys, os

from optparse import OptionParser

opt_parser = OptionParser()
opt_parser.add_option("-p", "--proportion", dest = "proportion", default = 1 - 0.618, type="float",\
                          help="Proportion of width occupied by region")
opt_parser.add_option("-o", "--output", help="Output file name", dest = "output_file_name", default="", metavar = "FILE")
(options, args) = opt_parser.parse_args()

output_file_name = options.output_file_name
proportion = options.proportion


if args:
    input_file_name = args[0]
    if not output_file_name:
        if len(args) >= 2:
            output_file_name = args[1]
        else:
            output_file_name = os.path.splitext(args[0])[0] + ".html"
    input_file = open(input_file_name)
    output_file = open(output_file_name, "w")
else:
    input_file = sys.stdin
    output_file = sys.stdout

#  html header and footer
html_header = """
<html>
<head>
<title> </title>
<link rel="stylesheet" href="css/table.css" type="text/css">
</head>
<body>
Last updated: <!--#echo var="LAST_MODIFIED"\ -->
<br>
<table class="dmrtable" cellspacing="0">
"""
html_footer = """
</table>
</body>
</html>
"""

output_file.write(html_header)

odd_row_start = """\t<tr class="odd"><td>%d</td>"""
even_row_start = """\t<tr><td>%d</td>"""
location_td = """<td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%(CHROM)s%%3A%(START)d-%(END)d">%(CHROM)s:%(START)d-%(END)d</a></td>"""

index = 1

for line in input_file:
    if index % 2 == 0:
        row = even_row_start % index
    else:
        row = odd_row_start % index
	    
    parts = line.split()
	# build link to Genome Browser
    chrom = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    extension = int((end - start) / proportion * (1 - proportion) / 2)
    start -= extension
    end += extension
    row += location_td % {"CHROM":chrom, "START":start, "END":end}
    
		# add other parts in bed file
    for part in parts[3:]:
        row += "<td>%s</td>" % part
        
		# ending td		
    row += "</tr>\n"
		
    output_file.write(row)
    index += 1
    

output_file.write(html_footer)

if not output_file == sys.stdout:
    output_file.close()

if not input_file == sys.stdin:    
    input_file.close()	
