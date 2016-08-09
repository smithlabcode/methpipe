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
#

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

PROGS = allelicmeth amrfinder amrtester

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

amrfinder: $(addprefix $(COMMON_DIR)/, EpireadStats.o Epiread.o) \
	$(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

amrtester: $(addprefix $(COMMON_DIR)/, EpireadStats.o Epiread.o) \
		$(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

allelicmeth:   $(addprefix $(COMMON_DIR)/, Epiread.o) \
    $(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
