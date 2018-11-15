RDynState.so: DynState

DynState: matalloc.o RDynState.o
	g++ -shared  -o RDynState.so matalloc.o RDynState.o -lgomp

RDynState.o: RDynState.cpp
	g++ -I${R_INCLUDE_DIR} -fpic  -g -O2 -c -fopenmp RDynState.cpp -o RDynState.o

matalloc.o: matalloc.cpp
	g++ -I${R_INCLUDE_DIR} -fpic  -g -O2 -c -fopenmp matalloc.cpp -o matalloc.o
