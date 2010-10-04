#    Makefile from smithlab_cpp C++ code library
#
#    Copyright (C) 2010 University of Southern California and
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

HEADERS = $(shell ls *.hpp)
SOURCES = $(shell ls *.cpp)
OBJECTS = $(subst .cpp,.o,$(SOURCES))

BINDIR = $(ROOT)/bin
LIBDIR = $(ROOT)/lib
INCLUDEDIR = $(ROOT)/include

STATIC_LIB = libmethpipe.a
DYNAMIC_LIB = libmethpipe.so

CXX = g++
CXXFLAGS = -Wall -fPIC -fmessage-length=50 
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif


install: $(DYNAMIC_LIB) $(STATIC_LIB)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(INCLUDEDIR)
	@install -m 755 $(DYNAMIC_LIB) $(LIBDIR)
	@install -m 644 $(STATIC_LIB) $(LIBDIR)
	@install -m 644 $(HEADERS) $(INCLUDEDIR)


$(DYNAMIC_LIB): $(OBJECTS)
	$(CXX) -shared $(CXXFLAGS) -o $@ $^

$(STATIC_LIB): $(OBJECTS)
	ar cru $@ $^

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(INCLUDEDIR)

.PHONY: clean
clean: 
	@-rm -f *.o *.so *.a *~

