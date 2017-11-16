# created by Joshua Shepherd

CC = g++ #mpic++ #nvcc
CPPFLAGS = -Wall -fopenmp -O3 -lpthread -std=c++11 #-DVERBOSE #-DDEBUG

#CFLAGS=-std=c99 -Wall -O0 -m32# -DVERBOSE# -Iinclude

#CFLAGS += #-Wno-unused-variable
#VPATH = src include
#CPPFLAGS=-std=c++11 -O3 -Wall
#LDFLAGS=
SRC = pi_approx.cpp
#INCLUDE=
OBJ = $(SRC:.cpp=.o)
TARGET = pi_approx

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(TARGET) $(OBJ)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<
	
test: all
	./$(TARGET) 1024 16

clean:
	rm -rf $(OBJ) $(TARGET)
