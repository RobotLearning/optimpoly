DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/carma/include
HEADER2=$(DIR)/optim/include
CC=g++
LIBS=-larmadillo -lm
INSTALLFLAGS=-fPIC -g -I$(HEADER1) -shared -pthread -std=c++11#-O2
TESTFLAGS=-g --std=c++11 #-O2 -pthread
OPTIMFLAGS=-g -Wall -O3 -I$(HEADER2)

opt:
	$(CC) $(OPTIMFLAGS) optim/src/optimpoly.c \
					optim/src/table_tennis.c \
					optim/src/kinematics.c \
					optim/src/utils.c \
	      -lm -lnlopt -o optim.o

test-opt:
	

test-carma:
	$(CC) $(TESTFLAGS) carma/test/table_tennis.cpp \
					 carma/test/optim.cpp \
		             carma/test/kinematics.cpp \
	                 carma/test/kalman.cpp \
	                  -o unit_tests.o \
	                   $(LIBS) -I$(HEADER1) -I/usr/local/include \
	                   /usr/local/lib/libboost_unit_test_framework.a ./libcarma.so -lnlopt                               

install:
	$(CC) $(INSTALLFLAGS) carma/src/player.cpp \
					      carma/src/optim.cpp \
					      carma/src/kinematics.cpp \
					      carma/src/kalman.cpp \
	                      carma/src/extkalman.cpp \
	                      carma/src/table_tennis.cpp $(LIBS) -o libcarma.so
	#gcc -c -I$(HEADER) src/carma.c -o carma.o
	#g++ -Wall -c -I$(HEADER) src/carma.cpp $(LIBS) -o carma.o
	#ar rcs libcarma.a carma.o

clean:
	rm -rf *.a *.o *.so

.PHONY: all test clean
