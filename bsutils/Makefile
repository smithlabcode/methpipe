#    Copyright (C) 2008 University of Southern California and
#                       Andrew D. Smith
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

PROGS = cpgcaller cpgcaller2 \
	cgicaller methcounts \
	bsrate    bsrate2    \
	coverage  coverage2  \
	sortreads cpgsmoother \
	methdiff  mergelanes \
	dmrfinder  dmr\
	allelic

LIBS = -lgsl -lgslcblas
LIBDIR = $(RMAP)/lib
INCLUDEDIR = $(RMAP)/include
CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g
COMMON_DIR = $(RMAP)/src/common/
TEST_DIR = $(RMAP)/test

ifndef RMAP
$(error must define RMAP variable)
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
	$(CXX) $(CFLAGS) -c -o $@ $< -I$(COMMON_DIR)

cpgcaller cpgcaller2 cgicaller cgicaller2 \
bsrate bsrate2 coverage coverage2 methdiff \
sortreads cpgsmoother methcounts mergelanes \
dmrfinder dmr: \
	$(addprefix $(COMMON_DIR), GenomicRegion.o rmap_os.o \
	rmap_utils.o OptionParser.o)

mergelanes cgicaller: $(addprefix $(COMMON_DIR), RNG.o)

methcounts bsrate2 coverage2 dmrfinder: FileIterator.o

allelic : \
	$(addprefix $(COMMON_DIR), GenomicRegion.o rmap_os.o \
	rmap_utils.o OptionParser.o) FileIterator.o

install: all
	@mkdir -p $(RMAP)/bin
	@install -m 755 $(PROGS) $(RMAP)/bin

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -I$(COMMON_DIR) -L$(LIBDIR) $(LIBS)

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
