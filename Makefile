HEADER=$(HOME)/optimpoly/include
ROBOTDIR=$(HOME)/robolab
CC=g++
LIBS=-lm -lnlopt -lpthread
CFLAGS= -g -Wall -DUNIX -O3 -I$(HEADER) #g is necessary for debugging

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