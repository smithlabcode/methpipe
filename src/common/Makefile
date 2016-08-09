#    Copyright (C) 2011 University of Southern California and
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

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50 -std=c++11
OPTFLAGS = -O2
DEBUGFLAGS = -g

INCLUDEDIRS =  $(SMITHLAB_CPP)/
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

clean: 
	@-rm -f *.o *~
.PHONY: clean
