#    Makefile from rmap software
#
#    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
#                       University of Southern California and
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

ifndef SRC_ROOT
SRC_ROOT=../..
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=../smithlab_cpp
endif

PROGS = rmapbs rmapbs-pe

CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS = -g
COMMON_DIR = ../common
TEST_DIR = $(SRC_ROOT)/test

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
	$(CXX) $(CFLAGS) -c -o $@ $< -I$(COMMON_DIR) -I$(SMITHLAB_CPP)

rmapbs rmapbs-pe : $(addprefix $(COMMON_DIR)/, SeedMaker.o FastRead.o load_reads.o \
	clip_adaptor_from_reads.o) \
	$(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o smithlab_utils.o \
	OptionParser.o QualityScore.o)

rmapbs-pe: $(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -I$(COMMON_DIR) -I$(SMITHLAB_CPP)

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
