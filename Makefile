DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/include/player
HEADER2=$(DIR)/include/optim
LIBDIR=$(DIR)/lib
CC=g++
LIBS=-larmadillo -lm
INSTALLFLAGS=-fPIC -g -Wall -I$(HEADER1) -I$(HEADER2) -shared -pthread -std=c++11 -O0
TESTFLAGS=-g --std=c++11 -pthread
OPTIMFLAGS=-fPIC -g -Wall -shared -I$(HEADER2) -O3 -std=c++11

# for compiling everything 
all: install interface lookup kinematics 

# for compiling only necessary stuff to play table tennis (in test mode)
install: player filter tabletennis optim

tabletennis:
	$(CC) $(INSTALLFLAGS) src/player/table_tennis.cpp $(LIBS) -lboost_program_options -o $(LIBDIR)/libtennis.so

kinematics:
	$(CC) $(INSTALLFLAGS) src/player/kinematics.cpp $(LIBS) -o $(LIBDIR)/libkin.so

lookup:
	$(CC) $(INSTALLFLAGS) src/player/lookup.cpp $(LIBS) -o $(LIBDIR)/liblookup.so
	
filter:
	$(CC) $(INSTALLFLAGS) src/player/kalman.cpp src/player/extkalman.cpp $(LIBS) -o $(LIBDIR)/libfilter.so		  
						  
player:
	$(CC) $(INSTALLFLAGS) src/player/player.cpp $(LIBS) -o $(LIBDIR)/libplayer.so

interface:
	$(CC) $(INSTALLFLAGS) src/player/sl_interface.cpp $(LIBS) -lboost_program_options -o $(LIBDIR)/libinterface.so

optim:
	$(CC) $(OPTIMFLAGS) src/optim/optimpoly.cpp \
						src/optim/lazyoptim.cpp \
						src/optim/invkin.cpp \
					    src/optim/kinematics.c \
					    src/optim/utils.c \
	                    -lm -o $(LIBDIR)/liboptim.so

test:
	$(CC) $(TESTFLAGS) test/test_optim.cpp -o unit_tests.o \
	                   $(LIBS) /usr/local/lib/libboost_unit_test_framework.a -I$(HEADER1) -I$(HEADER2) \
	                   $(LIBDIR)/liblookup.so $(LIBDIR)/libfilter.so $(LIBDIR)/libplayer.so \
	                   $(LIBDIR)/libtennis.so $(LIBDIR)/libkin.so $(LIBDIR)/liboptim.so -lnlopt
					    
clean:
	rm -rf *.a *.o lib/*.so

.PHONY: all test clean player optim

