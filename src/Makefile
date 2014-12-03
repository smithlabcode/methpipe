#    Makefile from methpipe software
#
#    Copyright (C) 2010-2014 University of Southern California and
#                            Andrew D. Smith
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

all_subdirs=common utils analysis amrfinder mlml samtools radmeth
lib_subdirs=common
app_subdirs=analysis utils amrfinder mlml radmeth

# check if a global copy of smithlab_cpp cannot be found and try to
# use a copy that is in the current directory
ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/smithlab_cpp
ifeq ("$(wildcard $(SMITHLAB_CPP))","")
$(error SMITHLAB_CPP not set and smithlab_cpp not found)
endif
all_subdirs += $(SMITHLAB_CPP)
lib_subdirs += $(SMITHLAB_CPP)
endif

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
