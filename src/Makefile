.PHONY: all clean

all: do
	./do test/folded.gro test/test.xtc prefix

do: do.cc
	g++ do.cc -I ../deps/xtc/install/ -o do -L ../deps/xtc/install -lccxtc_x86_64

do.cc:

clean:
	-rm do
