ROBOTDIR=$(HOME)/robolab
IDIR1=$(ROBOTDIR)/include
IDIR2=$(ROBOTDIR)/shared/barrett/include
IDIR3=$(ROBOTDIR)/shared/include
IDIR4=$(ROBOTDIR)/shared/barrett/math #math folder for kinematics computations
CC=gcc
LIBS=-lm -lnlopt
CFLAGS= -DUNIX -I$(IDIR1) -I$(IDIR2) -I$(IDIR3) -I$(IDIR4)
#DEPS=$(IDIR)/SL.h

all:
	$(CC) $(CFLAGS) src/optimpoly.c \
					src/example.c \
	                $(ROBOTDIR)/shared/utilities/src/utility.c \
	      $(LIBS) -o test 

clean:
	rm -rf *.o test
