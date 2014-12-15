program_NAME := blazar++
program_CXX_SRCS := $(wildcard src/*.cpp)
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_CXX_OBJS)

CXX := g++ -g

program_INCLUDE_DIRS += include
LDFLAGS := $(shell pkg-config --libs gsl)

program_INCLUDE_DIRS += $(SCFGP_DIR)/include
program_LIBRARY_DIRS += $(SCFGP_DIR)/lib
program_INCLUDE_DIRS += $(BAZINGA_DIR)/include
program_LIBRARY_DIRS += $(BAZINGA_DIR)/lib
program_LIBRARIES += scfgp bazinga

ROOT=NO
ifeq ($(ROOT),YES)	
	program_INCLUDE_DIRS += $(shell root-config --incdir)
	CPPFLAGS += -DUSE_ROOT
	LDFLAGS += $(shell root-config --libs)
endif

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir)) 
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

.PHONY: all clean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(program_OBJS) -o $(program_NAME)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	@- $(RM) scripts/createJob/*pyc
