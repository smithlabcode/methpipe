#    Makefile from common code directory of RSEG source tree
#
#    Copyright (C) 2011 University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Qiang Song and Andrew D. Smith
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
$(error SMITHLAB_CPP variable undefined)
endif

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS =  deadzones rseg rseg-diff

INCLUDEDIRS = $(SMITHLAB_CPP)  $(SAMTOOL_LIB)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lz -lm -lgsl -lgslcblas -lhts

CXX = g++
CC = gcc
CCFLAGS = -Wall -std=c++11
CXXFLAGS = -O3 -Wall
DEBUGFLAGS = -g -lefence -L$(HOME)/lib

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

rseg: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o) \
	$(addprefix $(COMMON_LIB)/, numerical_utils.o rseg_utils.o \
	Distro.o SplitDistro.o Smoothing.o BAMFile.o \
	SortGenomicRegion.o LoadReadsByRegion.o SelectBinSize.o ReadCounts.o \
	TwoStateScaleHMM.o TwoStateScaleResolveMixture.o \
	EvaluateBoundaries.o)

rseg-diff: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o) \
	$(addprefix $(COMMON_LIB)/, numerical_utils.o rseg_utils.o \
	Distro.o SplitDistro.o Smoothing.o BAMFile.o \
	SortGenomicRegion.o LoadReadsByRegion.o SelectBinSize.o ReadCounts.o \
	TwoStateScaleSplitHMM.o TwoStateScaleSplitResolveMixture.o \
	ThreeStateScaleSplitHMM.o ThreeStateScaleSplitResolveMixture.o \
	EvaluateBoundaries.o)

deadzones: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS) $(LIBS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *~

.PHONY: clean
