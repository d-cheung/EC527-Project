default:
	nvcc -I/ad/eng/support/software/linux/all/x86_64/cuda/cuda_sdk/C/common/inc -I /ad/eng/support/software/linux/all/x86_64/cuda/cuda/include/ -L/ad/eng/support/software/linux/all/x86_64/cuda/cuda/lib64/ *.cu -o surf -lcuda -lcudart -lm -lrt

gdb:
	nvcc -I/ad/eng/support/software/linux/all/x86_64/cuda/cuda_sdk/C/common/inc -I /ad/eng/support/software/linux/all/x86_64/cuda/cuda/include/ -L/ad/eng/support/software/linux/all/x86_64/cuda/cuda/lib64/ -g *.cu -o surf -lcuda -lcudart -lm -lrt

clean:
	rm -f surf *~
