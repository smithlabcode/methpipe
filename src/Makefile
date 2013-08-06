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

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(METHPIPE_ROOT)/src/smithlab_cpp
endif

all_subdirs=smithlab_cpp common postmapping analysis
lib_subdirs=smithlab_cpp common
app_subdirs=postmapping analysis utils

## COMMENT: The smithlab_cpp directory is only appearing in the
## "clean" target below because the compilation of object files inside
## that directory is achieved through dependency structure -- the
## reason we do it this way is because much in smithlab_cpp is not
## needed for methpipe.

all:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(METHPIPE_ROOT) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(METHPIPE_ROOT) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(METHPIPE_ROOT) test; \
	done;

clean:
	@for i in $(all_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(METHPIPE_ROOT) clean; \
	done;
