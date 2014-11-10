#!/usr/bin/env python
# bigWig_to_methcounts.py: a tool to conver MethBase tracks to
# methcounts format.
# 
# Copyright (C) 2014 University of Southern California and
#                          Meng Zhou
# 
# Authors: Meng Zhou
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""This script is used for converting tracks of MethBase in bigWig (.bw)
format to methcounts format.
"""

import sys, os
import subprocess, tempfile
from optparse import OptionParser

def parse_line(line):
  field = line.split()
  chr = field[0]
  start = field[1]
  value = field[3]

  return ([chr, start], value)

def write_line(fh, chr, start, meth_level, coverage):
  outline = "\t".join((chr, start, "+", "CpG", \
    str(meth_level), str(coverage))) + "\n"
  fh.write(outline)

def is_exe(file):
  return os.path.isfile(file) and os.access(file, os.X_OK)

def which(program):
  """Do the same thing as linux "which" command.
  """
  for path in os.environ["PATH"].split(os.pathsep):
    path = path.strip('"')
    exe_file = os.path.join(path, program)
    if is_exe(exe_file):
      return exe_file

  return None

def opt_validation(parser, opt):
  if not opt.meth or not opt.read:
    parser.print_help()
    sys.exit(0)
  if not opt.bwtool:
    if which("bigWigToBedGraph"):
      opt.bwtool = which("bigWigToBedGraph")
    else:
      sys.stderr.write("Cannot locate bigWigToBedGraph. Please specify path.\n")
      sys.exit(1)
  else:
    opt.bwtool = os.path.abspath(opt.bwtool)
  if not is_exe(opt.bwtool):
    sys.stderr.write(\
      "%s is not a proper executable file. Please check your path!\n"%opt.bwtool)
    sys.exit(1)

def main():
  usage = "Usage: %prog -m <methylation track> -r <coverage track>" + \
      " -o <output methcounts> [-p <path to bigWigToBedGraph>]"
  parser = OptionParser(usage=usage)
  parser.add_option("-m", "--methylation", action="store", type="string",
    dest="meth", help="MethBase methylation track file.", \
    metavar="<TRACK.meth.bw>")
  parser.add_option("-r", "--read-coverage", action="store", type="string",
    dest="read", help="MethBase read coverage track file.", \
    metavar="<TRACK.read.bw>")
  parser.add_option("-p", "--path", action="store", type="string", \
    dest="bwtool", \
    help="Path to bigWigToBedGraph executable file. " + \
    "Leave blank if you already have it in environment path.", metavar="<PATH>")
  parser.add_option("-o", "--output", action="store", type="string", \
    dest="output", \
    help="Output methcounts file.")
  (opt, args) = parser.parse_args(sys.argv)
  opt_validation(parser, opt)

  # conver bw files
  methtmp = tempfile.NamedTemporaryFile()
  readtmp = tempfile.NamedTemporaryFile()
  convert_args = [opt.bwtool, opt.meth, methtmp.name]
  try:
    subprocess.check_call(convert_args)
  except subprocess.CalledProcessError:
    sys.stderr.write(
      "An error occured in converting track file %s\n"%opt.meth)
    sys.exit(1)
  convert_args = [opt.bwtool, opt.read, readtmp.name]
  try:
    subprocess.check_call(convert_args)
  except subprocess.CalledProcessError:
    sys.stderr.write(
      "An error occured in converting track file %s\n"%opt.read)
    sys.exit(1)

  # combine converted files
  outputfh = open(opt.output, 'w')
  meth_line = methtmp.readline()
  read_line = readtmp.readline()
  while read_line and meth_line:
    (meth_coordinate, meth_value) = parse_line(meth_line)
    (read_coordinate, read_value) = parse_line(read_line)
    order = cmp(meth_coordinate, read_coordinate)
    if order == 0:
      write_line(outputfh, read_coordinate[0], read_coordinate[1], \
        meth_value, read_value)
      meth_line = methtmp.readline()
      read_line = readtmp.readline()
    elif order == 1:
      # site missing in methylation track
      write_line(outputfh, read_coordinate[0], read_coordinate[1], \
        0, 0)
      read_line = readtmp.readline()
    else:
      # site missing in read track
      write_line(outputfh, meth_coordinate[0], meth_coordinate[1], \
        0, 0)
      meth_line = methtmp.readline()

  methtmp.close()
  readtmp.close()
  outputfh.close()

if __name__ == '__main__':
  main()
