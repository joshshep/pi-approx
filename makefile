# created by Joshua Shepherd

CC = mpic++ #nvcc
CPPFLAGS = -Wall -lm -O3 #-DDEBUG

#CFLAGS=-std=c99 -Wall -O0 -m32# -DVERBOSE# -Iinclude

#CFLAGS += #-Wno-unused-variable
#VPATH = src include
#CPPFLAGS=-std=c++11 -O3 -Wall
#LDFLAGS=
SRC=parallel_prefix.cpp
#INCLUDE=
OBJ = $(SRC:.cpp=.o)
TARGET=parallel_prefix
SUB_SCRIPT = sub.bash

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(TARGET) $(OBJ)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<
	
run: all
	sbatch -N8 $(SUB_SCRIPT)

watch_run: all
	sbatch -N8 $(SUB_SCRIPT)
	./monitor.bash

watch_run_all: all
	sbatch -N1 sub1.bash
	sbatch -N2 sub2.bash
	sbatch -N4 sub4.bash
	sbatch -N8 sub8.bash
	sbatch -N8 sub16.bash
	sbatch -N8 sub32.bash
	sbatch -N8 sub64.bash
	./monitor.bash
#concise:

#disk:

clean:
	rm -rf $(OBJ) $(TARGET) *.out
