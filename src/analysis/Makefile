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

PROGS = pmd methcounts bsrate hmr hypermr levels roimethstat \
	methstates methentropy hmr_rep multimethstat

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

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

levels hmr_rep hmr methcounts roimethstat hypermr pmd: \
	$(addprefix $(COMMON_DIR)/, MethpipeSite.o)

methcounts_sam: \
	$(addprefix $(COMMON_DIR)/, MethpipeSite.o)

hmr hmr_rep: $(addprefix $(COMMON_DIR)/, TwoStateHMM.o)

pmd: $(addprefix $(COMMON_DIR)/, bsutils.o \
	TwoStateHMM_PMD.o EmissionDistribution.o)

hypermr: $(addprefix $(COMMON_DIR)/, ThreeStateHMM.o Smoothing.o \
	Distro.o BetaBin.o numerical_utils.o)

levels: $(addprefix $(COMMON_DIR)/, LevelsCounter.o)

methcounts: $(addprefix $(COMMON_DIR)/, MethpipeSite.o)

roimethstat levels: $(addprefix $(COMMON_DIR)/, bsutils.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS) $(LDLIBS) $(LDFLAGS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
