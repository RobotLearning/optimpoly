DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/include/player
HEADER2=$(DIR)/include/optim
LIBDIR=$(DIR)/lib
CC=g++
LIBS=-larmadillo -lm
INSTALLFLAGS=-fPIC -g -Wall -I$(HEADER1) -I$(HEADER2) -shared -pthread -std=c++11 -O0
TESTFLAGS=-g --std=c++11 -pthread
OPTIMFLAGS=-fPIC -g -Wall -shared -I$(HEADER2) -O3 -std=c++11

# Definitions for ADOL-C and IPOPT example
ADOLCPATH = $(HOME)/adolc_base/include
ADOLCLIBDIR = $(HOME)/adolc_base/lib64
IPOPTLIBS = -Wl,--rpath -Wl,$(ADOLCLIBDIR) -L$(ADOLCLIBDIR) -ladolc -L$(HOME)/Downloads/Ipopt-3.12.8/lib -lipopt
IPOPTFLAGS = -I$(ADOLCPATH) -I$(HOME)/Downloads/Ipopt-3.12.8/include/coin -I$(HOME)/Downloads/Ipopt-3.12.8/include


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
	$(CC) $(TESTFLAGS) test/test_ipopt.cpp -o ipopt_ex -I$(HOME)/Downloads/Ipopt-3.12.8/include/coin \
						-L/home/okoc/Downloads/Ipopt-3.12.8/lib -lipopt 

test-autodiff-ipopt: # example for testing IPOPT optimization library + AUTODIFF 	                   
	$(CC) $(TESTFLAGS) test/test_auto_ipopt.cpp -o ipopt_autodiff.o $(IPOPTFLAGS) $(IPOPTLIBS)				

test-nlopt: # example for testing NLOPT + autodiff
	$(CC) $(TESTFLAGS) test/test_nlopt.cpp -o nlopt_autodiff.o -lnlopt -I$(HOME)/adolc_base/include \
						-L$(HOME)/adolc_base/lib64 -ladolc				
					    
clean:
	rm -rf *.a *.o lib/*.so

.PHONY: all test clean player optim

