#    Copyright (C) 2010 University of Southern California and
#                       Andrew D. Smith
#                       Tianyin Zhou
#
#    Authors: Andrew D. Smith and Tianyin Zhou
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

PROGS =	methcounts bsrate hmr methdiff dmr amr \
	sortreads mergelanes checkoverlap

LIBS = -lgsl -lgslcblas
CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g
SMITHLAB_CPP = $(METHPIPE)/src/smithlab_cpp/
COMMON = $(METHPIPE)/src/common/
TEST_DIR = $(METHPIPE)/test

# VPATH = $(SMITHLAB_CPP):$(COMMON)

ifndef METHPIPE
$(error must define METHPIPE variable)
endif

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< -I$(COMMON) -I$(SMITHLAB_CPP)

$(PROGS): \
	$(addprefix $(SMITHLAB_CPP), GenomicRegion.o rmap_os.o \
	rmap_utils.o OptionParser.o)

mergelanes: \
	$(addprefix $(SMITHLAB_CPP), RNG.o)

methcounts bsrate: \
	$(addprefix $(SMITHLAB_CPP), QualityScore.o)

checkoverlap: \
	$(addprefix $(SMITHLAB_CPP), MappedRead.o)

methcounts methdiff dmr amr bsrate: \
	$(addprefix $(SMITHLAB_CPP), FileIterator.o MappedRead.o) 

hmr : $(addprefix $(COMMON), TwoStateHMM.o)

install: all
	@mkdir -p $(METHPIPE)/bin
	@install -m 755 $(PROGS) $(METHPIPE)/bin

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -I$(SMITHLAB_CPP) -I$(COMMON) $(LIBS)

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
