
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -ffast-math -fno-common
OBJS = main.o wavelet.o matrix.o
TARGET = sw

ifeq ($(OS),Windows_NT)
	# TODO
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CXXFLAGS += -O3
	endif
	ifeq ($(UNAME_S),Darwin)
		CXXFLAGS += -stdlib=libc++ -Ofast
	endif
endif

# Default rule
all: compile
	./$(TARGET) 64 64 4

# Dependency rules
main.o : wavelet.h matrix.h vector.h
wavelet.o : wavelet.h matrix.h
matrix.o : vector.h

# Pattern rule to create an object file from a cpp file
# $@ expands to the target
# $^ expands to the dependencies
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

compile: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $^

clean:
	rm $(TARGET) $(OBJS)
