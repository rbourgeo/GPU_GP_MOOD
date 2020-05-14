C = nvcc
NVCCFLAGS = -arch=sm_60
CFLAGS = -std=c++11 -rdc=true

all: transport

transport: transport.cu
	$(C) $(NVCCFLAGS) $(CFLAGS) -o transport.exe 2D_bij.cu bc.cu set_equal.cu init.cu FE.cu transport.cu
clean:
	rm -f transport.exe *.dat *.o
