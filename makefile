DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/include/player
HEADER2=$(DIR)/include/optim
CC=g++
LIBS=-larmadillo -lm
INSTALLFLAGS=-fPIC -g -Wall -I$(HEADER1) -I$(HEADER2) -shared -pthread -std=c++11 -O0
TESTFLAGS=-g --std=c++11 -pthread
OPTIMFLAGS=-fPIC -g -Wall -shared -I$(HEADER2) -O3

# for compiling everything 
all: install interface lookup kinematics

# for compiling only necessary stuff to play table tennis (in test mode)
install: player filter tabletennis optim

tabletennis:
	$(CC) $(INSTALLFLAGS) src/player/table_tennis.cpp $(LIBS) \
	-lboost_program_options -o bin/libtennis.so

kinematics:
	$(CC) $(INSTALLFLAGS) src/player/kinematics.cpp $(LIBS) -o bin/libkin.so

lookup:
	$(CC) $(INSTALLFLAGS) src/player/lookup.cpp $(LIBS) -o bin/liblookup.so
	
filter:
	$(CC) $(INSTALLFLAGS) src/player/kalman.cpp src/player/extkalman.cpp $(LIBS) -o bin/libfilter.so		  
						  
player:
	$(CC) $(INSTALLFLAGS) src/player/player.cpp $(LIBS) -o bin/libplayer.so

interface:
	$(CC) $(INSTALLFLAGS) src/player/sl_interface.cpp $(LIBS) -lboost_program_options -o bin/libinterface.so

optim:
	$(CC) $(OPTIMFLAGS) src/optim/optimpoly.cpp \
						src/optim/lazyoptim.cpp \
						src/optim/invkin.cpp \
					    src/optim/kinematics.c \
					    src/optim/utils.c \
	                    -lm -o bin/liboptim.so

test:
	$(CC) $(TESTFLAGS) test/test_table_tennis.cpp \
	                  -o unit_tests.o -lm -larmadillo \
	                   $(LIBS) -I$(HEADER1) -I$(HEADER2) -I/usr/local/include \
	                   /usr/local/lib/libboost_unit_test_framework.a \
	                   bin/liblookup.so bin/libfilter.so bin/libplayer.so \
	                   bin/libtennis.so bin/libkin.so bin/liboptim.so -lnlopt
						#test/test_kinematics.cpp \
					    #test/test_optim.cpp \
					    #test/test_table_tennis.cpp \
					    #test/test_kalman.cpp
					    
clean:
	rm -rf *.a *.o *.so

.PHONY: all test clean player optim
