The MethPipe software package is a computational pipeline for
analyzing bisulfite sequencing data (WGBS and RRBS). MethPipe provides
tools methylation-specific technical evaluation of sequencing data,
and for estimating methylation levels at individual cytosines.
Additionally, MethPipe includes tools for identifying higher-level
methylation features, such as hypo-methylated regions (HMR), partially
methylated domains (PMD), hyper-methylated regions (HyperMR), and
allele-specific methylated regions (AMR).

Release 4.1.0
===================

This release should be more stable than 4.0.0 and should be easier to
build and install.

## Installing release 4.1.0

### Required libraries

* The GNU Scientific Library: this has always been required. It can be
  installed using `apt` on Linux, using `brew` on macOS, or from
  source available [here](http://www.gnu.org/software/gsl).
* The Zlib compression library. Most likely you already have this
  installed on your system. If not, it can be installed using `apt`
  on Linux through the package `zlib1g-dev`. On macOS, Zlib can be
  installed with `brew`.
* Optional: The HTSlib library, which can be installed through `brew`
  on macOS, through `apt` on Linux, or from source downloadable
  [here](https://github.com/samtools/htslib). This is only required
  for using the `to-mr` tool.

### Configuration

1. Download methpipe-4.1.0.tar.gz [here](https://github.com/smithlabcode/methpipe).
2. Unpack the archive:
```
$ tar -zxvf methpipe-4.1.0.tar.gz
```
3. Move into the methpipe directory and create a build directory:
```
$ cd methpipe-4.1.0
$ mkdir build && cd build
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

### Building and installing the tools

If you are still in the `build` directory, run `make` to compile the
tools, and then `make install` to install them. If your HTSlib is not
installed system-wide, then you might need to udpate your library
path:
```
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/htslib/lib
```

### Building and installing from source

We strongly recommend using methpipe through the latest stable release
under the releases section on GitHub. However, developers who wish to
work on the latest commits, which are potentially unstable, can
compile the cloned repository using the `Makefile` available in the
repository. To compile without programs that require HTSlib, simply
run:
```
make
```
If HTSlib is available and users which to compile all programs,
including `to-mr`, run:
```
make HAVE_HTSLIB=1
```

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
