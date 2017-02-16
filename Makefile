HEADER=$(HOME)/optimpoly/include
ROBOTDIR=$(HOME)/robolab
IDIR=$(ROBOTDIR)/shared/barrett/math #math folder for kinematics computations
CC=gcc
LIBS=-lm -lnlopt -lpthread
CFLAGS= -g -DUNIX -I$(HEADER) -I$(IDIR) #g is necessary for debugging

all:
	$(CC) $(CFLAGS) src/optimpoly.c \
					src/table_tennis.c \
					src/kinematics.c \
					src/utils.c \
					src/example.c \
	      $(LIBS) -o test 

clean:
	rm -rf *.o test

.PHONY: all test clean