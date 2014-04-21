
*********************************
Copyright and License Information
*******************************************************************************
Copyright (C) 2013
University of Southern California,
Andrew D. Smith, Qiang Song, Fang Fang, Ben Decato and Meng Zhou
  
Authors:  Andrew D. Smith, Qiang Song, Fang Fang, Ben Decato and Meng Zhou
  
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


***********************
Building and Installing 
*******************************************************************************
This software package has been designed to operate in a UNIX-like environment.
It has been tested on MacOS X Snow Leopard and Linux. 

Step 0
------
  This software package requires a functioning installation of the GNU 
  Scientific Library (GSL). If you don't already have this installed, you 
  will need to download and install it from http://www.gnu.org/software/gsl/

Step 1
------
  To build the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution  
  
  > make all 
  
Step 2
------
  To install the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution
  
  > make install
  
  This will place the binaries in the bin directory under the package root.
  They can be used directly from there without any additional steps. You can
  add that directory to your PATH environment variable to avoid having to 
  specify their full paths, or you can copy the binaries to another directory
  of your choice in your PATH 
  
  
*****
Usage
*******************************************************************************

Read methpipe-manual.pdf in the docs directory.

************************
Contacts and bug reports
*******************************************************************************
Andrew D. Smith
andrewds@usc.edu

Qiang Song
qiang.song@usc.edu

MethPipe and MethBase Users' Mailinglist
methpipe@googlegroups.com
http://groups.google.com/group/methpipe?hl=en





