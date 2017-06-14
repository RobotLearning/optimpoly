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

##### ALL SHARED LIBRARIES FOR POLYOPTIM
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
						src/optim/vhp.cpp \
					    src/optim/kinematics.c \
					    src/optim/utils.c \
	                    -lm -o $(LIBDIR)/liboptim.so
	                    
##### ALL TESTS ARE INCLUDED HERE
test:
	$(CC) $(TESTFLAGS) test/test_kinematics.cpp -o unit_tests.o \
	                   $(LIBS) /usr/local/lib/libboost_unit_test_framework.a -I$(HEADER1) -I$(HEADER2) \
	                   $(LIBDIR)/liblookup.so $(LIBDIR)/libplayer.so $(LIBDIR)/libfilter.so \
	                   $(LIBDIR)/libtennis.so $(LIBDIR)/libkin.so $(LIBDIR)/liboptim.so -lnlopt
	                   
test-ipopt: # example for testing IPOPT optimization library
	$(CC) $(TESTFLAGS) src/optim/test_ipopt.cpp -o ipopt_ex -I/is/ei/okoc/Downloads/CoinIpopt/build/include/coin \
						-L/is/ei/okoc/Downloads/CoinIpopt/build/lib -lipopt -lcoinmetis -lcoinmumps

test-autodiff: # example for testing automatic differentiation using ADOL-C library
	$(CC) $(TESTFLAGS) test/test_autodiff.cpp -o autodiff.o -I/is/ei/okoc/adolc_base/include \
						-L/is/ei/okoc/adolc_base/lib64 -ladolc
test-nlopt: # example for testing NLOPT + autodiff
	$(CC) $(TESTLFAGS) test/test_nlopt.cpp -o nlopt_autodiff.o -lnlopt						
					    
clean:
	rm -rf *.a *.o lib/*.so

.PHONY: all test clean player optim

