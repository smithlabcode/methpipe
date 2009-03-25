The programs in this directory are for processing data from a
bisulfite capture experiment. They require that the RMAP source be
present in order to compile, and that an RMAP environment variable be
set so the location of the RMAP common source is known at compile
time.

cpgcaller
=========
OUTPUT: The calls have the following meaning in the names of each BED
line. Fields are separated by colons, and they start with the name of
the island/cpg:
1) unmeth reads
2) meth reads
3) confident comparison (enough reads map so that 90% CI would be in
   0.25 window at a meth freq of 0.5, for possibly some other CpG)
4) confident call (90% CI in 0.25 window at meth freq for this CpG)
The score column of the BED line indicates the category of the CpG
according to whether it is called unmethylated (0), methylated (1),
partially methylated (and confident; 2) or no confident call (3).

cgicaller
=========
OUTPUT: The name is also divided as for cpgcaller with fields
delimited by colons (':'):
1) Expected freq of meth
2) CI bound
3) other CI bound
4) confident comparison (90% of CpGs covered in the island)
5) confident call (90% CI within 0.25 window)
The score column has the exact same meaning as for cpgcaller, above.
