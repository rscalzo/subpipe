# RS 2013/02/21:  Modified this makefile to compile on maipenrai (64-bit).
# Removed "-m32" flag, and changed libcfitsio.a to ours (in $SUBPIPEHOME/ext).
LIBDIR  := ${SUBPIPEHOME}/ext/lib
INCDIRS := . ${SUBPIPEHOME}/ext/include

CLIBS := -L$(LIBDIR) -Wl,-rpath=$(LIBDIR) -lcfitsio -lm

# You can use either gcc or g++ (the later is faster and gives better warnings)
#CC = gcc
CC = g++
CFLAGS = -Wall -pedantic -O3

# This block added for SDP
ifndef SUBPIPEHOME
  platform = $(shell uname -m)
  include ../Makefile.common.$(platform)
  INCDIRS := . ${CFITSIO}/include
  CLIBS := $(CLIBS) -L$(CFITSIO) -L$(CFITSIO)/lib -lcfitsio
endif

SRC := fits64to32_cal.c biascorr.c

OBJ := $(SRC:%.c=%.o)
INCLUDES := $(INCDIRS:%=-I%)

fits64to32_cal: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(CLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $^ $(INCLUDES)

biascorr.c: biascorr.h

.PHONY: clean

clean :
	-rm -f fits64to32_cal $(OBJ)
