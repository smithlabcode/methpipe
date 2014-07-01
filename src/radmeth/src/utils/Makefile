
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

all: $(BIN_DIR)/make_table $(BIN_DIR)/split_table

$(BIN_DIR)/make_table : make_table.cpp merge.o cpg.o $(SMITH_OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

$(BIN_DIR)/split_table : split_table.cpp table_splitter.o $(SMITH_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

clean:
	rm -f $(BIN_DIR)/make_table $(BIN_DIR)/split_table *.o