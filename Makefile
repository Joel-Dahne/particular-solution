CXX=g++
CFLAGS=-Wall -O3

INCS=-I/home/urathai/Sources/eigen-git-mirror/ -I/home/urathai/Sources/eigen-git-mirror/unsupported/test/mpreal/
LIBS=-L$(CURDIR) -larb -lflint -lmpfr -lgmp

LINK_TARGETS = build/particular-solution
OBJS = build/geom.o build/generate-matrix.o build/sigma.o build/enclose.o build/particular-solution.o
REBUILDABLES = $(OBJS) $(LINK_TARGETS)

.SECONDARY: $(OBJS)

all: $(LINK_TARGETS)

clean:
	rm -f $(REBUILDABLES)

build/%.o: src/%.cpp | build
	$(CXX) $(CFLAGS) $(INCS) -c $^ -o $@

build/particular-solution: $(OBJS) | build
	$(CXX) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

build:
	mkdir -p build
