
CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -pedantic -ffast-math -fno-common
BIN = sw
SRC = main.cpp matrix.cpp wavelet.cpp

ifeq ($(OS),Windows_NT)
	# TODO
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CFLAGS += -O3
	endif
	ifeq ($(UNAME_S),Darwin)
		CFLAGS += -stdlib=libc++ -Ofast
	endif
endif

compile:
	$(CC) $(CFLAGS) -o $(BIN) main.cpp matrix.cpp wavelet.cpp

clean:
	rm $(BIN)

test: compile
	./$(BIN) 64 64 4
