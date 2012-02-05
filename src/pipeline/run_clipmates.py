#!/usr/bin/env python

import sys, os.path, string
from datetime import date
from optparse import OptionParser

SCRIPT_HEADER = """#PBS -S /bin/sh
#PBS -e %(stderr)s
#PBS -o %(stdout)s
#PBS -l vmem=10100M
#PBS -l mem=10100M
#PBS -l pmem=10100M
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -q cmb

export LC_ALL=C

# Set the temporary dir
export TMPDIR=%(tmpDirName)s
"""

MAIN_SCRIPT_TEMPLATE = """

# Sort the files according to name
sort -S 10000M -k 4 \\
-o %(endOne)s.tmp \\
%(endOne)s

sort -S 13000M -k 4 \\
-o %(endTwo)s.tmp \\
%(endTwo)s

# Clip the ends of the mates
%(binDir)s/clipmates -s %(suffLen)d -L 1000 \\
-S %(statsOut)s \\
-o %(mergedOut)s.tmp \\
-T %(endOne)s.tmp \\
-A %(endTwo)s.tmp

# Remove the temporary files
rm %(endOne)s.tmp %(endTwo)s.tmp

# Now sort the ouptut by location
sort -S 10000M -k 1,1 -k 2,2g -k 6,6 -k 3,3g \\
-o %(mergedOut)s \\
%(mergedOut)s.tmp

# Remove the temporary merged file
rm %(mergedOut)s.tmp

"""

jobCount = 0

def getMergedName(endOneFile, endTwoFile):
    return os.path.splitext(string.replace(endOneFile, '_1_seq', '_seq'))[0]

def writeScript(noSuff, binDir, workDir, endOneFile, endTwoFile):
    global jobCount
    suffLen = 1
    if noSuff:
        suffLen = 0
    scrBase = "clipmates.%s.%s.%d" % \
        (date.today(), os.getpid(), jobCount)
    scriptFileName = scrBase + '.qsub'
    scriptFileName = os.path.join(workDir, scriptFileName)
    f = open(scriptFileName, "w")
    f.write(SCRIPT_HEADER % \
                { "tmpDirName" : workDir,
                  "stderr" : os.path.join(workDir, scrBase + '.ER'),
                  "stdout" : os.path.join(workDir, scrBase + '.OU') })

    mergedName = getMergedName(endOneFile, endTwoFile)
    
    statsOut = os.path.join(workDir, mergedName + '.clipstats')
    mergedOut = os.path.join(workDir, mergedName + '.mr')
    
    f.write(MAIN_SCRIPT_TEMPLATE % \
                { "binDir" : binDir,
                  "suffLen" : suffLen,
                  "statsOut" : statsOut,
                  "mergedOut" : mergedOut,
                  "endOne" : endOneFile,
                  "endTwo" : endTwoFile })
    f.close()
    jobCount += 1
    return scriptFileName

def msimatchOnly1vs2(a, b):
    a = os.path.basename(a)
    if len(a) != len(b): return False
    if b.find('_2_seq') < 0: return False
    dist = 0
    for i in range(len(a)):
        if a[i] == '1' and b[i] == '2':
            dist += 1
        elif a[i] != b[i]:
            dist += 2
    return (dist == 1)

def maybeFirstEndMappedFile(i):
    return i.find('_1_seq') > 0

def findSecondEndFile(firstEndFile):
    for i in os.listdir(os.path.dirname(firstEndFile)):
        if i.find('_2_seq') > 0 and msimatchOnly1vs2(firstEndFile, i):
            return i
    return ''

def findMatePairFiles(dirName):
    pairs = []
    files = [os.path.join(dirName, i) for i in os.listdir(dirName)]
    for i in files:
        if maybeFirstEndMappedFile(i):
            secondEndFile = os.path.join(dirName, findSecondEndFile(i))
            if len(secondEndFile) > 0:
                pairs.append((i, secondEndFile))
    return pairs

def main(argv):                         

    parser = OptionParser()
    parser.add_option("-W", "--workdir", action="store", type="string",
                      dest="workDir", help="input/output/work directory",
                      metavar="<dir>")
    parser.add_option("-B", "--bin", action="store", type="string",
                      dest="binDir", help="directory for binaries", metavar="<dir>")
    parser.add_option("-N", "--no-suff", action="store_true",
                      dest="noSuff", help="no read name suffix", metavar="<bool>")
    parser.add_option("-A", "--run", action="store_true",
                      dest="actuallyRun", help="actually run the script")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 5:
        parser.print_help()
        sys.exit(1)

    opt.workDir = os.path.abspath(opt.workDir)
    if not os.access(opt.workDir, os.W_OK):
        print >>sys.stderr, 'ERROR: work dir'
        sys.exit(1)

    for i in findMatePairFiles(opt.workDir):
        scriptFileName = writeScript(opt.noSuff,
                                     opt.binDir, opt.workDir, i[0], i[1])
        if opt.actuallyRun:
            os.system("qsub %s" % scriptFileName)

if __name__ == "__main__":
    main(sys.argv)
