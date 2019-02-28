# Makefile for particular-solutions

SHELL = /bin/sh

INCS = -I/usr/include/eigen3/ -I$(CURDIR)/src/
LIBS = -L$(CURDIR) -larb -lflint -lmpfr -lgmp

CC ?= gcc
CXX ?= g++

override CFLAGS := $(CFLAGS) -Wall -O3

AT = @

BUILD_DIRS = src

export

SOURCES = $(wildcard $(patsubst %, %/*.cpp, $(BUILD_DIRS)))

HEADERS = $(patsubst %, %/*.h, $(BUILD_DIRS))

OBJS = $(patsubst src/%.cpp, build/%.o, $(SOURCES))

EXMP_SOURCES = $(wildcard examples/*.cpp)
EXMPS = $(patsubst %.cpp, build/%, $(EXMP_SOURCES))

TEST_SOURCES = $(wildcard tests/*.cpp)
TESTS = $(patsubst %.cpp, build/%, $(TEST_SOURCES))

.SECONDARY: $(OBJS)

all: examples

clean:
	rm -rf build

examples: $(EXMPS)

tests: $(TESTS)

check: $(TESTS)
	$(AT)$(foreach prog, $(TESTS), $(prog) || exit $$?;)

build:
	mkdir -p build

build/%.o: src/%.cpp $(HEADERS) | build
	$(CXX) $(CFLAGS) $(INCS) -c $< -o $@

build/tests/%: tests/%.cpp $(OBJS) | build/tests
	$(CXX) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

build/tests:
	mkdir -p build/tests

build/examples/%: examples/%.cpp $(OBJS) | build/examples
	$(CXX) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

build/examples:
	mkdir -p build/examples
