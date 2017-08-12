HOST=$(shell hostname)
ifeq ($(HOST),sill) # new machine in the new MPI building
	BOOSTL=$(HOME)/install/lib
else
	BOOSTL=/usr/local/lib
endif

DIR=$(HOME)/polyoptim
HEADER1=$(DIR)/include/player
HEADER2=$(DIR)/include/optim
LIBDIR=$(DIR)/lib
CC=g++
LIBS=-larmadillo -lm -lboost_program_options
FLAGS=-I$(HEADER1) -I$(HEADER2) -pthread -std=c++11
EXTRA_WARNINGS=-Werror -Wextra -Weffc++ -pedantic -pedantic-errors \
-Waggregate-return -Wcast-align -Wcast-qual -Wconversion \
-Wdisabled-optimization -Wfloat-equal\
-Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k \
-Wimport  -Winline -Winvalid-pch -Wlong-long \
-Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn \
-Wpacked  -Wpadded -Wpointer-arith \
-Wredundant-decls -Wshadow -Wstack-protector \
-Wstrict-aliasing=2 -Wswitch-default \
-Wswitch-enum -Wunreachable-code -Wunused -Wundef \
-Wvariadic-macros -Wwrite-strings
RELEASE=-O3 -DNDEBUG
DEBUG=-DDEBUG -g -Wall -pedantic -Wextra #$(EXTRA_WARNINGS)
SHARED_OBJECT = $(LIBDIR)/libplayer.so
PLAYER_DIR = $(DIR)/src/player
OPTIM_DIR = $(DIR)/src/optim
OBJ_PLAYER_DIR = $(DIR)/obj/player
OBJ_OPTIM_DIR = $(DIR)/obj/optim
SRC_PLAYER = $(wildcard $(PLAYER_DIR)/*.cpp)
SRC_OPTIM = $(wildcard $(OPTIM_DIR)/*.cpp $(OPTIM_DIR)/*.c))
OBJS_PLAYER = $(addprefix $(OBJ_PLAYER_DIR)/,$(basename $(notdir $(SRC_PLAYER))))
OBJS_OPTIM = $(addprefix $(OBJ_OPTIM_DIR)/,$(basename $(notdir $(SRC_OPTIM))))
OBJS = $(addsuffix .o,$(OBJS_PLAYER) $(OBJS_OPTIM))
#$(info $$OBJS_PLAYER is [${OBJS}])

# for compiling everything, release and debug modes
release: FLAGS += $(RELEASE)
release: all
debug: FLAGS += $(DEBUG)
debug: all

all: $(SHARED_OBJECT)

$(SHARED_OBJECT) : $(OBJS)
	$(CC) -shared $(FLAGS) -o $@ $^ $(LIBS)

# ADD HEADER PREREQUISITES!
$(OBJ_PLAYER_DIR)/%.o : $(PLAYER_DIR)/%.cpp
	$(CC) -c -fPIC $(FLAGS) -o $@ $<

$(OBJ_OPTIM_DIR)/%.o : $(OPTIM_DIR)/%.cpp
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJ_OPTIM_DIR)/%.o : $(OPTIM_DIR)/%.c
	$(CC) -c -fPIC $(FLAGS) -o $@ $<

	                    
##### ALL TESTS ARE INCLUDED HERE
test:
	$(CC) $(FLAGS) test/test_table_tennis.cpp -o unit_tests.o \
	               $(SHARED_OBJECT) $(LIBS) $(BOOSTL)/libboost_unit_test_framework.a -lnlopt
					    
clean:
	rm -rf obj/player/*.o obj/optim/*.o lib/*.so

.PHONY: all release debug test