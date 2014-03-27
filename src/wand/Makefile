
ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

BIN_DIR = ../../bin

CXXFLAGS += -g -Wall -Wextra -pthread 

LIBS += -lgsl -lgslcblas

INCLUDE_DIRS = $(SMITHLAB_CPP)
INCLUDE_ARGS = $(addprefix -I,$(INCLUDE_DIRS))

SMITH_OBJ = $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
            smithlab_utils.o OptionParser.o)
							
OBJ = design.o regression.o loglikratio_test.o table_row.o pipeline.o \
      gsl_fitter.o $(SMITH_OBJ)

$(BIN_DIR)/wand : wand.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDE_ARGS) $(LIBS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE_ARGS)
	
clean:
	rm -f $(BIN_DIR)/wand *.o
