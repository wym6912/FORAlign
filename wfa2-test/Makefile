###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BIN=../bin
FOLDER_BUILD=build
FOLDER_BUILD_CPP=build/cpp
FOLDER_LIB=lib
FOLDER_TESTS=tests
FOLDER_ROOT=.
FOLDER_B1=b1
FOLDER_B2=b2
FOLDER_B3=b3
UNAME=$(shell uname)

CC:=$(CC)
CPP:=$(CXX)

CC_FLAGS=-Wall -g -fPIE

AR=ar
AR_FLAGS=-rsc

BUILD_WFA_PARALLEL=1
STATIC_TEST=0

ifeq ($(STATIC_TEST), 0)
	PTHREAD_FLAGS=-pthread -lm -fopenmp
else
	PTHREAD_FLAGS=-static -Wl,--whole-archive -lc -lstdc++ -ldl -lpthread -lm -fopenmp -Wl,--no-whole-archive
endif

###############################################################################
# Configuration rules
###############################################################################
LIB_WFA=$(FOLDER_LIB)/libwfa.a
LIB_WFA_CPP=$(FOLDER_LIB)/libwfacpp.a
SUBDIRS=alignment \
		benchmark \
		edit \
		gap_affine \
		gap_affine2p \
		gap_linear \
		indel \
        system \
        utils \
        wavefront \
		../foralign \
		../swg

APPS+=align_benchmark multiswg_benchmark error_measure
OBJS+=$(FOLDER_BUILD)/*.o


all: CC_FLAGS+=-O3 -march=native #-flto -ffat-lto-objects
all: build

debug: build

ASAN_OPT=-fsanitize=address -fsanitize=undefined -fsanitize=shift -fsanitize=alignment
ASAN_OPT+=-fsanitize=signed-integer-overflow -fsanitize=bool -fsanitize=enum
ASAN_OPT+=-fsanitize=pointer-compare -fsanitize=pointer-overflow -fsanitize=builtin

# ASAN: ASAN_OPTIONS=detect_leaks=1:symbolize=1 LSAN_OPTIONS=verbosity=2:log_threads=1
asan: CC_FLAGS+=$(ASAN_OPT) -fno-omit-frame-pointer -fno-common
asan: build

###############################################################################
# Build rules
###############################################################################
build: setup
build: $(SUBDIRS) 
build: lib_wfa 
build: $(APPS)

setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_BUILD_CPP) $(FOLDER_LIB) $(FOLDER_B1) $(FOLDER_B2) $(FOLDER_B3)
	
lib_wfa: $(SUBDIRS)
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null

clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB) $(FOLDER_B1) $(FOLDER_B2) $(FOLDER_B3) 2> /dev/null

$(FOLDER_B1)/align_benchmark.o: */*.c */*.h align_benchmark.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c align_benchmark.c -o $@ $(LD_FLAGS) -L$(FOLDER_LIB) -fopenmp -I$(FOLDER_ROOT) -I../foralign -lm

$(FOLDER_B1)/align_benchmark_params.o: */*.c */*.h align_benchmark.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c align_benchmark_params.c -o $@ $(LD_FLAGS) -L$(FOLDER_LIB) -fopenmp -I$(FOLDER_ROOT) -lm

align_benchmark: $(FOLDER_B1)/align_benchmark.o $(FOLDER_B1)/align_benchmark_params.o $(LIB_WFA)
	$(CPP) -o $(FOLDER_BIN)/$@ $(FOLDER_BUILD)/*.o $(FOLDER_B1)/*.o -L../lib -lforalign -L./lib -lswg $(PTHREAD_FLAGS)

$(FOLDER_B2)/multiswg_benchmark.o: */*.c */*.h align_benchmark.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c multiswg_benchmark.c -o $@ $(LD_FLAGS) -L$(FOLDER_LIB) -fopenmp -I$(FOLDER_ROOT) -I../foralign -I../parallel_swg -lm

$(FOLDER_B2)/multiswg_benchmark_params.o: */*.c */*.h align_benchmark.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c multiswg_benchmark_params.c -o $@ $(LD_FLAGS) -L$(FOLDER_LIB) -fopenmp -I$(FOLDER_ROOT) -lm

multiswg_benchmark: $(FOLDER_B2)/multiswg_benchmark.o $(FOLDER_B2)/multiswg_benchmark_params.o $(LIB_WFA)
	$(CPP) -o $(FOLDER_BIN)/$@ $(FOLDER_BUILD)/*.o $(FOLDER_B2)/*.o -L../lib -lswg $(PTHREAD_FLAGS)

$(FOLDER_B3)/error_measure.o: */*.c */*.h error_measure.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c error_measure.c -o $@ $(LD_FLAGS) -L$(FOLDER_LIB) -fopenmp -I$(FOLDER_ROOT) -I../foralign -lm

$(FOLDER_B3)/error_measure_params.o: */*.c */*.h error_measure_params.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c error_measure_params.c -o $@ $(LD_FLAGS) -L$(FOLDER_LIB) -fopenmp -I$(FOLDER_ROOT) -lm

error_measure: $(FOLDER_B3)/error_measure.o $(FOLDER_B3)/error_measure_params.o $(LIB_WFA)
	$(CPP) -o $(FOLDER_BIN)/$@ $(FOLDER_BUILD)/*.o $(FOLDER_B3)/*.o -L../lib -lforalign -L./lib -lswg $(PTHREAD_FLAGS)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	

.PHONY: $(SUBDIRS)

