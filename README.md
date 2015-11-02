The MethPipe software package is a computational pipeline for
analyzing bisulfite sequencing data (WGBS and RRBS). MethPipe
provides tools for mapping bisulfite sequencing read and estimating
methylation levels at individual cytosine sites. Additionally,
MethPipe also includes tools for identifying higher-level methylation
features, such as hypo-methylated regions (HMR), partially methylated
domains (PMD), hyper-methylated regions (HyperMR), and allele-specific
methylated regions (AMR).

Building and Installing 
=======================

You may download the latest stable release from:http://smithlabresearch.org/software/methpipe/  
This software package has been designed to operate in a UNIX-like environment.
It has been tested on MacOS X Snow Leopard and Linux. 

* Step 0
  
  This software package requires a functioning installation of the GNU 
  Scientific Library (GSL). If you don't already have this installed, you 
  will need to download and install it from http://www.gnu.org/software/gsl/

  If gsl is not installed in the default path, 
  ```
  export CPATH=/path_to_my_gsl/include
  export LIBRARY_PATH=/path_to_my_gsl/lib
  ```
  will add search paths for compiling and linking. 

* Step 1
  
  To build the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution:
  
      > make all

* Step 2
  
  To install the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution:
  
      > make install
  
  This will place the binaries in the bin directory under the package root.
  They can be used directly from there without any additional steps. You can
  add that directory to your PATH environment variable to avoid having to 
  specify their full paths, or you can copy the binaries to another directory
  of your choice in your PATH 

For advanced users who are interested in the newest features, you may obtain the 
latest source code by cloning the MethPipe repository:

    > git clone --recursive https://github.com/smithlabcode/methpipe.git

After you clone the latest source code, follow the above steps for installation.

Usage
=====

Read methpipe-manual.pdf in the docs directory.

Contacts and bug reports
========================

Andrew D. Smith
andrewds@usc.edu

Ben Decato
decato@usc.edu

MethPipe and MethBase Users' Mailinglist
methpipe@googlegroups.com
http://groups.google.com/group/methpipe?hl=en

Copyright and License Information
=================================

Copyright (C) 2013-2015
University of Southern California,
Andrew D. Smith
  
Current Authors:  Andrew D. Smith, Ben Decato, Meng Zhou, Liz Ji, Jenny Qu, Egor Dolzhenko
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

SamTools
Copyright (c) 2008-2009 Genome Research Ltd.
SamTools software is free software distributed under the MIT License.
Refer to the COPYING file in src/samtoos/ for further information.
