CC?=gcc
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
LOCATION = $(shell pwd)
MAKE_PLACE = $(LOCATION)/$(firstword $(MAKEFILE_LIST))
DIR = $(shell dirname $(MAKE_PLACE))

PREFIX = /usr/local
BINDIR = $(PREFIX)/bin
HAS_64 = 1

THREADS = 1

#valgrind
#CFLAGS += -g -O0
##


LIBOBJ_MAIN = hirschberg_single.o dp.o hirschberg.o mtset.o mtxutl.o seqio.o russians.o mtdp.o

TARGETS = libforalign.a
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

INPUT = input.txt
DOCS = 3

##

all:${TARGETS_ALL}

libforalign.a: $(LIBOBJ_MAIN)
	$(AR) rcs $@ $(LIBOBJ_MAIN)
	mkdir -p ../lib
	cp $@ ../lib

clean:
	${RM} *.o ${TARGETS} *.a ../lib/*.a

install: all
	mkdir -p $(BINDIR)
	cp $(TARGETS) $(SUBMAFFT_TARGETS) $(BINDIR)
	@echo Install to $(BINDIR) done.

###

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(OPT_FLAGS) -c -o $@ $<

%.o:%.c
	 $(CC) $(CFLAGS) -c -o $@ $<

##