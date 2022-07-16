#  Copyright (C) 2011 University of Southern California
#                     and Andrew D. Smith
#
#  Authors: Andrew D. Smith
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
#

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

PROGS = lc_approx.o fast-liftover.o lift-filter.o \
	merge-bsrate.o merge-methcounts.o \
	symmetric-cpgs.o clean-hairpins.o selectsites.o guessprotocol.o \
	format_reads.o duplicate-remover.o

COMMON_DIR = ../common
INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))
override CPPFLAGS += $(INCLUDEARGS)

LDLIBS = -lgsl -lgslcblas -lz -lhts

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

merge-methcounts.o symmetric-cpgs.o selectsites.o lift-filter.o \
fast-liftover.o guessprotocol.o: \
	$(addprefix $(COMMON_DIR)/, MethpipeSite.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS) $(LDLIBS) $(LDFLAGS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
