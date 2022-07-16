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

ifndef SRC_ROOT
SRC_ROOT=..
endif

CXX = g++
CXXFLAGS = -Wall -std=c++11
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

all_subdirs=common utils analysis amrfinder mlml radmeth
lib_subdirs=common
app_subdirs=analysis utils amrfinder mlml radmeth

all_subdirs += $(SMITHLAB_CPP)
lib_subdirs += $(SMITHLAB_CPP)
LDLIBS = -lgsl -lgslcblas -lz -lhts
LDFLAGS = $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

COMMON_DIR = common
INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))
override CPPFLAGS += $(INCLUDEARGS)

PROGS = methpipe
OBJS = amrfinder/allelicmeth.o amrfinder/amrfinder.o amrfinder/amrtester.o \
			 analysis/bsrate.o analysis/hmr.o analysis/hmr_rep.o analysis/hypermr.o \
			 analysis/levels.o analysis/methcounts.o analysis/methentropy.o \
			 analysis/methstates.o analysis/multimethstat.o analysis/pmd.o \
			 analysis/roimethstat.o mlml/mlml.o radmeth/dmr.o radmeth/methdiff.o \
			 radmeth/radmeth-adjust.o radmeth/radmeth-merge.o radmeth/radmeth.o \
			 utils/clean-hairpins.o utils/duplicate-remover.o utils/fast-liftover.o \
			 utils/format_reads.o utils/guessprotocol.o utils/lc_approx.o \
			 utils/lift-filter.o utils/merge-bsrate.o utils/merge-methcounts.o \
			 utils/symmetric-cpgs.o utils/selectsites.o

COMMON_OBJS = $(addprefix $(COMMON_DIR)/, \
              BetaBin.o bsutils.o Distro.o EmissionDistribution.o Epiread.o \
							EpireadStats.o LevelsCounter.o MethpipeSite.o numerical_utils.o \
							Smoothing.o ThreeStateHMM.o TwoStateHMM.o TwoStateHMM_PMD.o)


all: $(PROGS)
	make -C $(SMITHLAB_CPP) HAVE_HTSLIB=1 all
	@for i in $(app_subdirs); do \
		make -C $${i} \
			SMITHLAB_CPP=$(SMITHLAB_CPP) \
			SRC_ROOT=$(METHPIPE_ROOT); \
	done;

$(PROGS): methpipe.cpp $(OBJS) $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $<  $(OBJS) $(COMMON_OBJS) $(LDFLAGS) $(CPPFLAGS) $(LDLIBS)

install:
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

clean:
	make -C $(SMITHLAB_CPP) clean
	@for i in $(all_subdirs); do \
		make -C $${i} \
			SMITHLAB_CPP=$(SMITHLAB_CPP) \
			SRC_ROOT=$(METHPIPE_ROOT) clean; \
	done;
	@-rm -f $(PROGS) *.o *.so *.a *~
