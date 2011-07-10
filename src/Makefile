#    Makefile from methpipe software
#
#    Copyright (C) 2010 University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef METHPIPE_ROOT
$(error Must define METHPIPE_ROOT variable)
endif

all_subdirs=analysis cytosines mapper premapping common postmapping smithlab_cpp 
lib_subdirs=common smithlab_cpp
app_subdirs=analysis cytosines mapper premapping postmapping 

## COMMENT: The smithlab_cpp directory is only appearing in the
## "clean" target below because the compilation of object files inside
## that directory is achieved through dependency structure -- the
## reason we do it this way is because much in smithlab_cpp is not
## needed for methpipe.

all:
	@for i in $(lib_subdirs); do \
		make -C $${i} METHPIPE_ROOT=$(METHPIPE_ROOT) OPT=2; \
	done
	@for i in $(app_subdirs); do \
		make -C $${i} METHPIPE_ROOT=$(METHPIPE_ROOT)  OPT=2; \
	done;

install:
	@for i in $(lib_subdirs); do \
		make -C $${i} METHPIPE_ROOT=$(METHPIPE_ROOT) OPT=2; \
	done
	@for i in $(app_subdirs); do \
		make -C $${i} METHPIPE_ROOT=$(METHPIPE_ROOT)  OPT=2 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} METHPIPE_ROOT=$(METHPIPE_ROOT) test; \
	done;

clean:
	@for i in $(all_subdirs); do \
		make -C $${i} METHPIPE_ROOT=$(METHPIPE_ROOT) clean; \
	done;
