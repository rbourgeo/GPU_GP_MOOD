C = nvcc
NVCCFLAGS = -arch=sm_60 --compiler-options '-Wuninitialized,-Wall,-Wextra'
CFLAGS = -std=c++11 -rdc=true -lquadmath

all: transport

transport: transport.cu
	$(C) $(NVCCFLAGS) $(CFLAGS) -o transport.exe bc.cu Update.cu init.cu set_equal.cu recons_flux_computation.cu transport.cu
clean:
	rm -f transport.exe *.dat *.o
