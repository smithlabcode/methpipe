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

PROGS = radmeth.o radmeth-adjust.o radmeth-merge.o methdiff.o dmr.o

CXX = g++
CXXFLAGS = -Wall -std=c++11
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = ../common
INCLUDEDIRS =  $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))
override CPPFLAGS += $(INCLUDEARGS)

LDLIBS = -lgsl -lgslcblas -lz

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

methdiff.o: $(addprefix $(COMMON_DIR)/, MethpipeSite.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS) $(LDLIBS) $(LDFLAGS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
