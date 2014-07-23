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

CXXFLAGS += -g -Wall -Wextra -pthread

ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

LIBS += -lgsl -lgslcblas

INCLUDE_DIRS = $(SMITHLAB_CPP)
INCLUDE_ARGS = $(addprefix -I,$(INCLUDE_DIRS))

SMITH_OBJ = $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
            smithlab_utils.o OptionParser.o)

OBJ = design.o regression.o loglikratio_test.o table_row.o pipeline.o \
      gsl_fitter.o $(SMITH_OBJ)

wand : wand.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDE_ARGS) $(LIBS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE_ARGS)

clean:
	rm -f wand *.o
