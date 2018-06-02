CC=g++
TESTFLAGS=-g --std=c++11 -pthread
OPTFLAGS=-O3 -DNDEBUG -std=c++11 -pthread
DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/include/player
HEADER2=$(DIR)/include/optim
# Definitions for ADOL-C and IPOPT example
IPOPTBUILD = $(HOME)/Downloads/CoinIpopt/build
ADOLCPATH = $(HOME)/adolc_base/include
ADOLCLIBDIR = $(HOME)/adolc_base/lib64
IPOPTLIBS = -Wl,--rpath -Wl,$(ADOLCLIBDIR) -L$(ADOLCLIBDIR) -ladolc -L$(IPOPTBUILD)/lib -lipopt
IPOPTFLAGS = -I$(ADOLCPATH) -I$(IPOPTBUILD)/include/coin -I$(IPOPTBUILD)/include

test-ipopt: # example for testing IPOPT optimization library 
	$(CC) $(TESTFLAGS) ex_ipopt.cpp -o ipopt_ex -I$(IPOPTBUILD)/include/coin \
						-L$(IPOPTBUILD)/lib -lipopt

test-auto-ipopt: # example for testing IPOPT optimization library + AUTODIFF 	                   
	$(CC) $(TESTFLAGS) ex_auto_ipopt.cpp -o ipopt_autodiff $(IPOPTFLAGS) $(IPOPTLIBS)				

test-auto-nlopt: # example for testing NLOPT + autodiff
	$(CC) $(TESTFLAGS) ex_auto_nlopt.cpp -o nlopt_autodiff -lnlopt -I$(ADOLCPATH) \
						-L$(ADOLCLIBDIR) -ladolc	

test-est-init: # example for testing NLOPT + ball initial state estimation
	$(CC) $(TESTFLAGS) ex_estimate_prior.cpp -o est_prior -I$(HEADER1) -I$(HEADER2) \
	                  -lm -larmadillo -lnlopt ../lib/libplayer.so				

test-bvp: # example for testing NLOPT + Boundary Value Problem (BVP)
	$(CC) $(OPTFLAGS) ex_bvp.cpp -o solve_bvp -I$(HEADER1) -I$(HEADER2) \
	                  -lm -larmadillo ../lib/libplayer.so -lnlopt
					    
test-zmq: # example for testing ZMQ interface to vision data
	$(CC) $(OPTFLAGS) ex_zmq.cpp -o ex_zmq -I. -lzmq -lzmqpp

					    
clean:
	rm -rf *.a *.o 