The MethPipe software package is a computational pipeline for
analyzing bisulfite sequencing data (WGBS and RRBS). MethPipe provides
tools methylation-specific technical evaluation of sequencing data,
and for estimating methylation levels at individual cytosines.
Additionally, MethPipe includes tools for identifying higher-level
methylation features, such as hypo-methylated regions (HMR), partially
methylated domains (PMD), hyper-methylated regions (HyperMR), and
allele-specific methylated regions (AMR).

Major release 4.0.0
===================

This release may be unstable. Much has changed in 4.0.0 compared with
earlier releases. The main new functionality is the capacity for
reading most (large) input files in gzip format. Installing version
4.0.0 may be less convenient. It requires building and installing the
`smithlab_cpp` library first. The instructions for installing earlier
releases are below.

## Installing release 4.0.0

### Required libraries

* The GNU Scientific Library: this has always been required. It can be
  installed using `apt` on Linux, using `brew` on macOS, or from
  source available [here](http://www.gnu.org/software/gsl).
* The Zlib compression library. Most likely you already have this
  installed on your system. If not, it can be installed using `apt`
  on Linux through the package `zlib1g-dev`. On macOS, Zlib can be
  installed with `brew`.
* The `smithlab_cpp` library. We suggest you install a release, and do
  not attempt to simply clone the source and build the library unless
  you know what you are doing. Currently the best way to install the
  `smithlab_cpp` library is using `./configure && make && make
  install` but please see the instructions
  [here](https://github.com/smithlabcode/smithlab_cpp)
* Optional: The HTSlib library, which can be installed through `brew`
  on macOS, through `apt` on Linux, or from source downloadable
  [here](https://github.com/samtools/htslib). This is only required
  for using the `to-mr` tool, but it you plan to build `methpipe` with
  HTSlib support, then make sure you also build `smithlab_cpp` with
  HTSlib support.

### Configuration

1. Dowload methpipe-4.0.0.tar.gz [here](https://github.com/smithlabcode/methpipe).
2. Unpack the archive:
```
$ tar -zxvf methpipe-4.0.0.tar.gz
```
3. Move into the methpipe directory and create a build directory:
```
$ cd methpipe-4.0.0
$ mkdir build
$ cd build
```
4. Run the configuration script:
```
$ ../configure
```
If you do not want to install the methpipe system-wide, or if you do
not have admin privileges, specify a prefix directory:
```
$ ../configure --prefix=/some/reasonable/place
```
If you are specifying a particular location for installing `methpipe`
then you likely did the same thing when installing `smithlab_cpp`. If
so, you must configure as follows:
```
$ ../configure --prefix=/some/reasonable/place \
    CPPFLAGS='-I /path/to/smithlab_cpp/headers' \
    LDFLAGS='-L/path/to/smithlab_cpp/lib'
```
Notice that the argment to `LDFLAGS` does not have a space after the
`-L`. This next argument is not required, but will make the code run
slightly faster:
```
$ ../configure CXXFLAGS='-O3 -Wall'
```
Finally, if you want to build with HTSlib support (for the `to-mr`
program) then you need to specify the following:
```
$ ../configure --enable-hts
```
And if you installed HTSlib yourself in some non-standard directory,
you must specify the location like this:
```
$ ../configure --enable-hts CPPFLAGS='-I /path/to/htslib/headers' \
    LDFLAGS='-L/path/to/htslib/lib'
```
If you need to specify locations for both HTSlib and `smithlab_cpp`,
then they must be in the same sets of quotes, separated by spaces.
Hopefully this process will be easier soon.

### Building and installing the tools

If you are still in the `build` directory, run `make` to compile the
tools, and then `make install` to install them. If your HTSlib is not
installed system-wide, then you might need to udpate your library
path:
```
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/htslib/lib
```

Previous releases of Methpipe
=============================

Building and Installing
-----------------------

You may download the latest stable release from
http://smithlabresearch.org/software/methpipe/ This software package
has been designed to operate in a UNIX-like environment.  It has been
tested on MacOS X Snow Leopard and Linux.

* Step 0

  This software package requires a functioning installation of the GNU
  Scientific Library (GSL). If you don't already have this installed,
  you will need to download and install it from
  http://www.gnu.org/software/gsl/

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

For advanced users who are interested in the newest features, you may
obtain the latest source code by cloning the MethPipe repository:
```
> git clone --recursive https://github.com/smithlabcode/methpipe.git
```
After you clone the latest source code, follow the above steps for
installation.

Usage
=====

Read methpipe-manual.pdf in the docs directory.

Contacts and bug reports
========================

Andrew D. Smith
andrewds@usc.edu

Ben Decato
decato@usc.edu

Meng Zhou
mengzhou@usc.edu

MethPipe and MethBase Users' Mailinglist
methpipe@googlegroups.com
http://groups.google.com/group/methpipe

Copyright and License Information
=================================

Copyright (C) 2018-2020
University of Southern California,
Andrew D. Smith

Current Authors: Andrew D. Smith, Ben Decato, Meng Zhou, Liz Ji, Jenny Qu, Egor Dolzhenko

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
