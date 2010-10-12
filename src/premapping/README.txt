1. read-quality-prof.cpp
This program take a fastq file as input and output the base composition and 
quality scores for each column

2. quality_prof.R
This R script define a function that takes the output of read-quality_prof.cpp as input 
and draw the figure of base composition

3. trim-adapter.cpp 
This program expects a fastq file as input and trim the adapter sequence from the 3' end
of reads if there is.

4. visireads.cpp
This program takes a fastq file as input and output a BED file displaying Cs in 
the sequences
