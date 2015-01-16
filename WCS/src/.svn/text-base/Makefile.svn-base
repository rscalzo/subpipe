# RS 2011/10/28:  Adapted from Brian's original makefile.

LIBS = -L$(EXTROOT)/lib -lm -lcfitsio -lwcs
INC = -I$(EXTROOT)/include -I../include
CC = gcc
CFLAGS = -g -ansi $(INC)
CIR_OBJS = cir_memory.o cir_wcssubs.o cir_stamp.o

BINARIES = brian_fitwcs brian_xi_eta finalmatch starmatch_new

all:  $(BINARIES)

brian_xi_eta: brian_xi_eta.o $(CIR_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

brian_fitwcs: brian_fitwcs.o invert_c.o minimrq.o select.o nrutil.o $(CIR_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

finalmatch:  finalmatch.o nrutil.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

starmatch_new:  starmatch_new.o indexx.o hunt.o nrutil.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

install:  $(BINARIES)
	cp $(BINARIES) ../bin

clean:
	rm -f $(BINARIES) *.o *~
