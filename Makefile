ROBOTDIR=$(HOME)/robolab
IDIR1=$(ROBOTDIR)/include
IDIR2=$(ROBOTDIR)/shared/barrett/include
CC=gcc
LIBS=-lm -lnlopt
CFLAGS=-I$(IDIR1) -I$(IDIR2)
#DEPS=$(IDIR)/SL.h

all:
	$(CC) $(CFLAGS) src/optimpoly.c $(LIBS) -o test 

clean:
	rm -rf *.o test
