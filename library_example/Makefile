CC = gcc
CPP = g++

STATIC_LINK = 1

ifeq ($(STATIC_LINK), 1)
	STATIC_ARG = -static
else
	STATIC_ARG = 
endif
WFLAGS= -Wall -Wextra  -DNDEBUG -Wno-ignored-qualifiers
#WFLAGS=-Wall -Werror -Wno-unused-function
OPT_FLAGS= -g -O3 -ffast-math -funroll-loops -m64 -fomit-frame-pointer -D_FILE_OFFSET_BITS=64 -I../include/boost/include -I../include/threadpool/include -I../wfa-test/
CFLAGS = $(WFLAGS) $(OPT_FLAGS) $(STATIC_ARG)
CPPFLAGS = ${CFLAGS} -std=c++17 $(OPT_FLAGS)
LDFLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
SYSTEM?=$(shell uname -s)
DIR = ../bin

PREFIX = /usr/local
BINDIR = $(PREFIX)/bin
HAS_64 = 1

THREADS = 1

#valgrind
#CFLAGS += -g -O0
##


TARGETS = $(DIR)/foralign-cpp $(DIR)/foralign-c $(DIR)/swg

FORALIGNFLIB = ../lib/libforalign.a
SWGLIB     = ../lib/libswg.a
##

TARGETS_ALL = $(TARGETS)

DEBUG = 0
M64 = 0
DNA = 0
TERMINATOR = 1
DEFINES = -DDEBUG=$(DEBUG)
DEFINES += -DTERMINATOR=$(TERMINATOR)
CFLAGS += $(DEFINES)
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0
endif
##

all:${TARGETS_ALL}

$(DIR)/foralign-cpp: foralign.cpp $(FORALIGNFLIB)
	mkdir -p $(DIR)
	$(CPP) -o $@ foralign.cpp $(CPPFLAGS) $(LDFLAGS) $(OPT_FLAGS) $(FORALIGNFLIB)

$(DIR)/foralign-c: foralign.c $(FORALIGNFLIB)
	mkdir -p $(DIR)
	$(CC) -o $@ foralign.c $(CFLAGS) $(LDFLAGS) $(OPT_FLAGS) $(FORALIGNFLIB) -lstdc++

$(DIR)/swg: swg.cpp $(SWGLIB) $(FORALIGNFLIB)
	mkdir -p $(DIR)
	$(CPP) -o $@ swg.cpp $(CPPFLAGS) $(LDFLAGS) $(OPT_FLAGS) $(SWGLIB) $(FORALIGNFLIB)

$(FORALIGNFLIB):
	make -C ../foralign

$(SWGLIB):
	make -C ../swg

clean:
	${RM} -rf $(DIR)

install: all
	mkdir -p $(BINDIR)
	cp $(TARGETS) $(SUBMAFFT_TARGETS) $(BINDIR)
	@echo Install to $(BINDIR) done.

###

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(OPT_FLAGS) -c -o $@ $<

%.o:%.c
	 $(CC) $(CFLAGS) -c -o $@ $<