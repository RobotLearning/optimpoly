DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/carma/include
HEADER2=$(DIR)/optim/include
CC=g++
LIBS=-larmadillo -lm
INSTALLFLAGS=-fPIC -g -I$(HEADER1) -shared -pthread -std=c++11#-O2
TESTFLAGS=-g --std=c++11 #-O2 -pthread
OPTIMFLAGS=-fPIC -g -Wall -O3 -shared -I$(HEADER2)

carma:
	$(CC) $(INSTALLFLAGS) carma/src/player.cpp \
					      carma/src/kalman.cpp \
	                      carma/src/extkalman.cpp \
	                      carma/src/table_tennis.cpp $(LIBS) -o libcarma.so
	#gcc -c -I$(HEADER) src/carma.c -o carma.o
	#g++ -Wall -c -I$(HEADER) src/carma.cpp $(LIBS) -o carma.o
	#ar rcs libcarma.a carma.o

optim:
	$(CC) $(OPTIMFLAGS) optim/src/optimpoly.c \
					optim/src/kinematics.c \
					optim/src/utils.c \
	      -lm -o liboptim.so

test-optim:
	$(CC) $(TESTFLAGS) optim/test/optim.cpp \
	                 -o opt_unit_test.o -lm -larmadillo \
	                   -I$(HEADER1) -I$(HEADER2) -I/usr/local/include \
	                   /usr/local/lib/libboost_unit_test_framework.a ./liboptim.so -lnlopt ./libcarma.so

test-carma:
	$(CC) $(TESTFLAGS) carma/test/table_tennis.cpp \
					 carma/test/optim.cpp \
		             carma/test/kinematics.cpp \
	                 carma/test/kalman.cpp \
	                  -o carma_unit_tests.o \
	                   $(LIBS) -I$(HEADER1) -I/usr/local/include \
	                   /usr/local/lib/libboost_unit_test_framework.a ./libcarma.so -lnlopt	                  

clean:
	rm -rf *.a *.o *.so

.PHONY: all test clean carma optim
