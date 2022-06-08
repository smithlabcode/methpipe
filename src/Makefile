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

# check if recursive clone was done
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/smithlab_cpp

ifeq (,$(wildcard $(SMITHLAB_CPP)/Makefile))
$(error src/smithlab_cpp does not have a Makefile. \
	      Did you use --recursive when running git clone?)
endif

all_subdirs=common utils analysis rseg amrfinder mlml radmeth
lib_subdirs=common
app_subdirs=analysis rseg utils amrfinder mlml radmeth

all_subdirs += $(SMITHLAB_CPP)
lib_subdirs += $(SMITHLAB_CPP)

all:
	make -C $(SMITHLAB_CPP) HAVE_HTSLIB=1 all
	@for i in $(app_subdirs); do \
		make -C $${i} \
			SMITHLAB_CPP=$(SMITHLAB_CPP) \
			SRC_ROOT=$(METHPIPE_ROOT); \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} \
			SMITHLAB_CPP=$(SMITHLAB_CPP) \
			SRC_ROOT=$(METHPIPE_ROOT) install; \
	done;

clean:
	make -C $(SMITHLAB_CPP) clean
	@for i in $(all_subdirs); do \
		make -C $${i} \
			SMITHLAB_CPP=$(SMITHLAB_CPP) \
			SRC_ROOT=$(METHPIPE_ROOT) clean; \
	done;
