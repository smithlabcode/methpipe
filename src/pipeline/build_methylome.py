#!/usr/bin/env python

import sys, os.path
from optparse import OptionParser

SCRIPT_HEADER = """#PBS -S /bin/sh
#PBS -e %(stderr)s
#PBS -o %(stdout)s
#PBS -l mem=%(memSize)dG
#PBS -l pmem=%(coreMemSize)dG
#PBS -l nodes=1:ppn=%(cores)d
#PBS -l walltime=48:00:00
#PBS -q cmb

# Set the path just in case...
export PATH=${HOME}/bin:${PATH}

# And this one is no joke...
export LC_ALL=C

"""

#### Need to make changes:
####
#### (1) The --batch-size=nmerge must be large enough so that no temporary file is used when merging
#### (2) The -T or --temporary-directory=tmpdir should be used instread of changing TMPDIR
#### (3) Use a pipe to avoid creating an input file for the duplicate-remover program
####

MAIN_SCRIPT_TEMPLATE = """
# Merge all the reads files, which are assumed to have been sorted based on their first ends
sort -T %(tmpDirName)s --batch-size=1000 -S %(memSize)dG -m -k 1,1 -k 2,2g -k 3,3g -k 6,6 \\
%(inFileNames)s | \\
%(binDir)s/duplicate-remover -stdin \\
-S %(outFileName)s.u1_dup_stats | \\
%(binDir)s/reorder -stdin | \\
%(binDir)s/duplicate-remover -stdin -B \\
-o %(outFileName)s.mr \\
-S %(outFileName)s.u2_dup_stats
"""

SCRIPT_FOOTER = """
# Now merge all the library output files
sort -T %(tmpDirName)s --batch-size=1000 -m -S %(memSize)dG -k 1,1 -k 3,3g -k 2,2g -k 6,6 \\
-o %(outFileName)s.mr \\
%(toMerge)s
""" 

def parseLibraryFiles(fileOfFiles):
    return ' '.join([os.path.abspath(i.strip()) 
                     for i in open(fileOfFiles) if len(i.strip()) > 0])

def main(argv):                         

    cores = 4

    parser = OptionParser()
    parser.add_option("-o", "--outdir", action="store", type="string",
                      dest="outDir", help="output directory",
                      metavar="<dir>")
    parser.add_option("-L", "--libs", action="store", type="string",
                      dest="libraryFiles", help="comma separated library files",
                      metavar="<files>")
    parser.add_option("-T", "--tmp", action="store", type="string",
                      dest="tmpDir", help="temporary files dir", metavar="<dir>")
    parser.add_option("-B", "--bin", action="store", type="string",
                      dest="binDir", help="directory for binaries", metavar="<dir>")
    parser.add_option("-N", "--name", action="store", type="string",
                      dest="methylomeName", help="name of methylome", 
                      metavar="<string>")
    parser.add_option("-S", "--size", action="store", type="int",
                      dest="memorySize", help="memory size (in GB)", 
                      metavar="<int>")
    parser.add_option("-A", "--run", action="store_true",
                      dest="actuallyRun", help="actually run the script")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 13:
        parser.print_help()
        sys.exit(1)

    opt.binDir = os.path.abspath(opt.binDir)
    if not os.access(opt.binDir, os.R_OK):
        print >>sys.stderr, 'ERROR: bin dir'
        sys.exit(1)

    opt.tmpDir = os.path.abspath(opt.tmpDir)
    if not os.access(opt.tmpDir, os.W_OK):
        print >>sys.stderr, 'ERROR: tmp dir'
        sys.exit(1)

    opt.outDir = os.path.abspath(opt.outDir)
    if not os.access(opt.outDir, os.W_OK):
        os.makedirs(opt.outDir)

    scriptFileName = os.path.join(opt.outDir, 'scr.qsub')
    stderr = os.path.join(opt.outDir, opt.methylomeName + '.err')
    stdout = os.path.join(opt.outDir, opt.methylomeName + '.out')

    scriptText = SCRIPT_HEADER % { 'tmpDirName' : opt.tmpDir,
                                   'stderr' : stderr,
                                   'stdout' : stdout,
                                   'memSize' : opt.memorySize,
                                   'coreMemSize' : opt.memorySize/cores,
                                   'cores' : cores }

    counter = 1
    mergedFiles = []
    for i in opt.libraryFiles.split(','):
        scriptText += MAIN_SCRIPT_TEMPLATE % \
            { 'tmpDirName' : opt.tmpDir,
              'binDir' : opt.binDir,
              'tmpDirName' : opt.tmpDir, 
              'inFileNames' : parseLibraryFiles(i),
              'outFileName' : os.path.join(opt.tmpDir, opt.methylomeName + str(counter)),
              'memSize' : opt.memorySize }
        mergedFiles.append(os.path.join(opt.tmpDir, opt.methylomeName + str(counter)))
        counter += 1

    scriptText += SCRIPT_FOOTER % \
        { 'tmpDirName' : opt.tmpDir,
          'memSize' : opt.memorySize,
          'outFileName' : os.path.join(opt.outDir, opt.methylomeName),
          'toMerge' : ' '.join([i + '.mr' for i in mergedFiles]) }

    f = open(scriptFileName, "w")
    f.write(scriptText)
    f.close()

    if opt.actuallyRun:
        os.system("qsub -q cmb %s" % scriptFileName)
        
if __name__ == "__main__":
    main(sys.argv)
