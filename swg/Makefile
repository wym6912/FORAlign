CPP=g++
###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../lib

###############################################################################
# Modules
###############################################################################
MODULES = parallel_swg
SRCS=$(addsuffix .cpp, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.cpp=.o))

CPP_FLAGS= -O3 -pthread
GPPVERSIONGTEQ8 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 8)
ifeq ($(GPPVERSIONGTEQ8), 1)
	CPP_FLAGS += -std=c++2a
else
	CPP_FLAGS += -std=c++17
endif

###############################################################################
# Rules
###############################################################################

all: $(FOLDER_BUILD)/libswg.a

$(FOLDER_BUILD)/libswg.a: $(FOLDER_BUILD)/parallel_swg.o
	$(AR) rc $@ $(FOLDER_BUILD)/parallel_swg.o
	rm $(FOLDER_BUILD)/parallel_swg.o 

# General building rule
$(FOLDER_BUILD)/%.o : %.cpp
	$(CPP) $(CPP_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@
