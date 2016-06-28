HEADER=$(HOME)/optimpoly/include
ROBOTDIR=$(HOME)/robolab
IDIR1=$(ROBOTDIR)/include
IDIR2=$(ROBOTDIR)/shared/barrett/include
IDIR3=$(ROBOTDIR)/shared/include
IDIR4=$(ROBOTDIR)/shared/barrett/math #math folder for kinematics computations
IDIR5=$(ROBOTDIR)/barrett/include
CC=gcc
LIBS=-lm -lnlopt -lpthread
CFLAGS= -g -DUNIX -I$(HEADER) -I$(IDIR1) -I$(IDIR2) -I$(IDIR3) -I$(IDIR4) -I$(IDIR5) #g is necessary for debugging

all:
	$(CC) $(CFLAGS) src/optimpoly.c \
					src/table_tennis.c \
					src/kinematics.c \
					src/extra.c \
					src/example.c \
	                $(ROBOTDIR)/shared/utilities/src/utility.c \
	                $(ROBOTDIR)/shared/utilities/src/ludcmp.c \
	                $(ROBOTDIR)/shared/utilities/src/lubksb.c \
	                $(ROBOTDIR)/shared/utilities/src/svdcmp.c \
	                $(ROBOTDIR)/shared/utilities/src/pythag.c \
	      $(LIBS) -o test 

clean:
	rm -rf *.o test
