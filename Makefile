
CC = g++
FLAGS = -std=c++11 -stdlib=libc++ -Wall -Wextra -pedantic
CC_OPT = -O3 -ffast-math -fno-common
BIN = sw
SRC = main.cpp matrix.cpp wavelet.cpp

compile:
	$(CC) $(FLAGS) $(CC_OPT) -o $(BIN) $(SRC)

test: compile
	time ./$(BIN)
