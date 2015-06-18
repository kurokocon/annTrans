CC=g++
CFLAGS=-c -Wall -O2 -std=c++11

all: anno_transfer
        $(CC) anno_transfer.o -o anno_transfer
anno_transfer: anno_transfer.cpp
        $(CC) anno_transfer.cpp $(CFLAGS)
clean:
        rm -f anno_transfer anno_transfer.o

