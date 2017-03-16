DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/player/include
HEADER2=$(DIR)/optim/include
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
	$(CC) $(INSTALLFLAGS) player/src/table_tennis.cpp $(LIBS) -o libtennis.so

kinematics:
	$(CC) $(INSTALLFLAGS) player/src/kinematics.cpp $(LIBS) -o libkin.so

lookup:
	$(CC) $(INSTALLFLAGS) player/src/lookup.cpp $(LIBS) -o liblookup.so
	
filter:
	$(CC) $(INSTALLFLAGS) player/src/kalman.cpp player/src/extkalman.cpp $(LIBS) -o libfilter.so				  
						  
player:
	$(CC) $(INSTALLFLAGS) player/src/player.cpp $(LIBS) -o libplayer.so

interface:
	$(CC) $(INSTALLFLAGS) player/src/carma.cpp $(LIBS) -o libcarma.so

optim:
	$(CC) $(OPTIMFLAGS) optim/src/optimpoly.cpp \
						optim/src/lazyoptim.cpp \
						optim/src/invkin.cpp \
					    optim/src/kinematics.cpp \
					    optim/src/utils.c \
	                    -lm -o liboptim.so

test:
	$(CC) $(TESTFLAGS) player/test/kalman.cpp \
	                  -o unit_tests.o -lm -larmadillo \
	                   $(LIBS) -I$(HEADER1) -I$(HEADER2) -I/usr/local/include \
	                   /usr/local/lib/libboost_unit_test_framework.a \
	                   ./liblookup.so ./libfilter.so ./libplayer.so ./libtennis.so ./libkin.so ./liboptim.so -lnlopt
						#optim/test/kinematics.cpp \
					    #optim/test/optim.cpp \
					    #player/test/table_tennis.cpp \
					    
clean:
	rm -rf *.a *.o *.so

.PHONY: all test clean player optim
