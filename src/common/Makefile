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

ifndef METHPIPE_ROOT
$(error Must define METHPIPE_ROOT variable)
endif

HEADERS = $(wildcard *.hpp)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(subst .cpp,.o,$(SOURCES))

SMITHLAB_DIR = ../smithlab_cpp/
COMMON_DIR = ../common/
INCLUDEDIRS =  $(SMITHLAB_DIR) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50 
OPTFLAGS = -O2
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

all: $(OBJECTS)

install: all

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDEARGS) -c  -o $@ $<

.PHONY: clean
clean: 
	@-rm -f *.o *.so *.a *~

