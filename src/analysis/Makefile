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

ifndef SRC_ROOT
$(error Must define SRC_ROOT variable)
endif

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

PROGS = methcounts bsrate hmr amr dmr methdiff roimethstat

COMMON_DIR = $(SRC_ROOT)/src/common
INCLUDEDIRS =  $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lgsl -lgslcblas 

CXX = g++
CFLAGS = -Wall -fmessage-length=50
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

amr: $(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o \
	FileIterator.o MappedRead.o) 

bsrate:	$(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o \
	FileIterator.o MappedRead.o QualityScore.o)

dmr: $(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o)

hmr: $(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o) \
	$(addprefix $(COMMON_DIR)/, TwoStateHMM.o)

methcounts: $(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o \
	FileIterator.o MappedRead.o QualityScore.o)		

methdiff: $(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o)

roimethstat: $(addprefix $(SMITHLAB_CPP)/, \
	GenomicRegion.o OptionParser.o rmap_os.o rmap_utils.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
