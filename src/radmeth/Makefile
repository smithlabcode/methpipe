#    Copyright (C) 2013 University of Southern California and
#                       Egor Dolzhenko
#                       Andrew D Smith
#
#    Authors: Andrew D. Smith and Egor Dolzhenko
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

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

PROGS = radmeth methdiff dmr

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50 -std=c++11
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = ../common
INCLUDEDIRS =  $(SMITHLAB_CPP) $(COMMON_DIR)

INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lgsl -lgslcblas

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, \
	smithlab_os.o smithlab_utils.o GenomicRegion.o OptionParser.o) \
	$(addprefix $(COMMON_DIR)/, MethpipeFiles.o)

OBJS= regression.o combine_pvals.o merge.o

radmeth: radmeth.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
