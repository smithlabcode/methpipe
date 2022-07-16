#  Copyright (C) 2014 University of Southern California
#                     and Andrew D. Smith and Benjamin E. Decato
#
#  Authors: Andrew D. Smith and Benjamin E. Decato
#
#  This is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301 USA

PROGS = pmd.o methcounts.o bsrate.o hmr.o hypermr.o levels.o roimethstat.o \
	methstates.o methentropy.o hmr_rep.o multimethstat.o

CXX = g++
CXXFLAGS = -Wall -std=c++11
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifndef SRC_ROOT
SRC_ROOT=../..
endif

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = ../common
INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))
override CPPFLAGS += $(INCLUDEARGS)

LDLIBS = -lgsl -lgslcblas -lz -lhts

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

levels.o hmr_rep.o hmr.o methcounts.o roimethstat.o hypermr.o pmd.o: \
	$(addprefix $(COMMON_DIR)/, MethpipeSite.o)

hmr.o hmr_rep.o: $(addprefix $(COMMON_DIR)/, TwoStateHMM.o)

pmd.o: $(addprefix $(COMMON_DIR)/, bsutils.o \
	TwoStateHMM_PMD.o EmissionDistribution.o)

hypermr.o: $(addprefix $(COMMON_DIR)/, ThreeStateHMM.o Smoothing.o \
	Distro.o BetaBin.o numerical_utils.o)

levels.o: $(addprefix $(COMMON_DIR)/, LevelsCounter.o)

methcounts.o: $(addprefix $(COMMON_DIR)/, MethpipeSite.o)

roimethstat.o levels.o: $(addprefix $(COMMON_DIR)/, bsutils.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS) $(LDLIBS) $(LDFLAGS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
