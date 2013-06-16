.PHONY: all clean

all: dihedralPCA	
	./dihedralPCA prefix test/rama.xvg

dihedralPCA: dihedralPCA.cc
	g++ dihedralPCA.cc -I ../deps/ezOptionParser/ -o dihedralPCA

dihedralPCA.cc:

clean:
	-rm dihedralPCA
