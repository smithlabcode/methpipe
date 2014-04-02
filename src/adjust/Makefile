
ifndef SMITHLAB_CPP
SMITHLAB_CPP=../../smithlab_cpp/
endif

BIN_DIR = ../../bin

CXXFLAGS += -g -Wall -Wextra -pthread

ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif


LIBS += -lgsl -lgslcblas

INCLUDEDIRS = $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

SMITH_OBJ = $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
							smithlab_utils.o OptionParser.o)
							
OBJ = locus.o bin_for_distance.o distance_correlation.o proximal_loci.o \
			stoufferliptak.o combine_pvals.o fdr.o pvallocus.o $(SMITH_OBJ)

all: $(BIN_DIR)/adjust $(BIN_DIR)/dmrs

$(BIN_DIR)/adjust : adjust.cpp $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

$(BIN_DIR)/dmrs : dmrs.cpp merge.o locus.o $(SMITH_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

clean:
	rm -f $(BIN_DIR)/adjust $(BIN_DIR)/dmrs *.o
