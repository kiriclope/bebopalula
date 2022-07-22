SRC = main.cu matrix_utils.cu globals.h

a.out : main.cu cuda_utils.h utils.h cuda_globals.h globals.h
	nvcc -arch=sm_35 -std=c++11 --compiler-options -mcmodel=medium -O3 main.cu 

debug : main.cu
	nvcc -arch=sm_35 -g -G -lineinfo $<

clean:
	-rm *.out

.PHONY: clean
