
objects = TNT.o TNT_gold.o

all: $(objects)
	nvcc -arch=sm_20 -lmpi -L/uaopt/openmpi/1.10.2/lib -I /uaopt/openmpi/1.10.2/include $(objects) -o mpi_implementation

TNT_gold.o: TNT_gold.cpp
	g++ -std=c++0x -c TNT_gold.cpp
TNT.o: TNT.cu
	nvcc -lmpi -c -L/uaopt/openmpi/1.10.2/lib -I /uaopt/openmpi/1.10.2/include -I /uaopt/cuda/7.0.28/samples/common/inc -I. TNT.cu

clean:
	rm -f *.o 
