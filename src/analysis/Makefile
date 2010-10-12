# 
#  Copyright (C) 2008 Cold Spring Harbor Laboratory
#  Authors: Andrew D. Smith
# 
#  This file is part of rmap.
# 
#  rmap is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  rmap is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with rmap; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
### ==========================

ifndef ROOT
$(error Must define ROOT variable)
endif

SOURCES = $(shell find . -name "*.cpp")
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS =  $(patsubst %.cpp,%,$(SOURCES))

STATLIBS = $(shell find $(ROOT) -name "*.a")
INCLUDEDIRS = $(shell find $(ROOT) -name "*.hpp"|xargs -n 1 dirname|sort|uniq)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -Wl,-Bdynamic  -lgsl -lgslcblas 

CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
CXXFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g -lefence -L$(HOME)/lib

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(PROGS) $(ROOT)/bin

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCLUDEARGS) $(STATLIBS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean

