# Makefile for particular-solutions

SHELL = /bin/sh

INCS = -I$(HOME)/Sources/eigen/ -I$(CURDIR)/src/
LIBS = -L$(CURDIR) -larb -lflint -lmpfr -lgmp -lm -lstdc++

CC ?= gcc
CXX ?= g++

override CFLAGS := $(CFLAGS) -Wall -O3

PIC_FLAG = -fPIC

AT = @

BUILD_DIRS = src

export

SOURCES = $(wildcard $(patsubst %, %/*.c, $(BUILD_DIRS)))

SOURCESXX = $(wildcard $(patsubst %, %/*.cpp, $(BUILD_DIRS)))

HEADERS = $(patsubst %, %/*.h, $(BUILD_DIRS)) $(patsubst %, %/*.hpp, $(BUILD_DIRS))

OBJS = $(patsubst src/%.c, build/%.o, $(SOURCES)) $(patsubst src/%.cpp, build/%.o, $(SOURCESXX))

EXMP_SOURCES = $(wildcard examples/*.c)
EXMPS = $(patsubst %.c, build/%, $(EXMP_SOURCES))

TEST_SOURCES = $(wildcard tests/*.c)
TESTS = $(patsubst %.c, build/%, $(TEST_SOURCES))

PS_LIB = build/particular_solution.so

.SECONDARY: $(OBJS)

all: examples

clean:
	rm -rf build

examples: $(EXMPS)

tests: $(TESTS)

shared: $(PS_LIB)

check: $(TESTS)
	$(AT)$(foreach prog, $(TESTS), $(prog) || exit $$?;)

build:
	mkdir -p build

build/sigma_eigen.o: src/sigma_eigen.cpp $(HEADERS) | build
	$(CXX) $(PIC_FLAG) $(CFLAGS) $(INCS) -c $< -o $@

build/%.o: src/%.c $(HEADERS) | build
	$(CC) $(PIC_FLAG) $(CFLAGS) $(INCS) -c $< -o $@

build/examples/%: examples/%.c $(OBJS) | build/examples
	$(CC) $(PIC_FLAG) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

build/examples:
	mkdir -p build/examples

build/tests/%: tests/%.c $(OBJS) | build/tests
	$(CC) $(PIC_FLAG) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

build/tests:
	mkdir -p build/tests

$(PS_LIB): $(OBJS) | build
	$(CC) -shared $(OBJS) $(INCS) -o $@ $(LIBS)
