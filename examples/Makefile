CC=g++
TESTFLAGS=-g --std=c++11 -pthread

# Definitions for ADOL-C and IPOPT example
IPOPTBUILD = $(HOME)/Downloads/Ipopt-3.12.8			#CoinIpopt/build
ADOLCPATH = $(HOME)/adolc_base/include
ADOLCLIBDIR = $(HOME)/adolc_base/lib64
IPOPTLIBS = -Wl,--rpath -Wl,$(ADOLCLIBDIR) -L$(ADOLCLIBDIR) -ladolc -L$(IPOPTBUILD)/lib -lipopt
IPOPTFLAGS = -I$(ADOLCPATH) -I$(IPOPTBUILD)/include/coin -I$(IPOPTBUILD)/include

ipopt: # example for testing IPOPT optimization library 
	$(CC) $(TESTFLAGS) ex_ipopt.cpp -o ipopt_ex -I$(IPOPTBUILD)/include/coin \
						-L$(IPOPTBUILD)/lib -lipopt

auto-ipopt: # example for testing IPOPT optimization library + AUTODIFF 	                   
	$(CC) $(TESTFLAGS) ex_auto_ipopt.cpp -o ipopt_autodiff $(IPOPTFLAGS) $(IPOPTLIBS)				

auto-nlopt: # example for testing NLOPT + autodiff
	$(CC) $(TESTFLAGS) ex_auto_nlopt.cpp -o nlopt_autodiff -lnlopt -I$(ADOLCPATH) \
						-L$(ADOLCLIBDIR) -ladolc			
					    
clean:
	rm -rf *.a *.o 