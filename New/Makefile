Cpp=g++
## 1.	This is for SSE2 support (on GCC).
CcOpt=-msse2 -O3 -DAE_OS=AE_LINUX -DAE_CPU=AE_INTEL -DAE_HAVE_STDINT
## 2.	This is for SSE2 and multi-core support (on GCC).
##	It does not seem to provide any additional benefit.
## CcOpt=-mtune=core2 -msse2 -O3 -DAE_OS=AE_LINUX -DAE_CPU=AE_INTEL
## 3.	Without -O3, nothing worthwhile happens with GCC (on Linux)
##	and the CPU-specific initializations are no longer needed either.
#CcOpt=-DAE_OS=AE_LINUX
Libs=-lstdc++ -lm -lpthread
X=
#X=.exe
O=.o
#O=.obj

all: Idw0$X Idw1$X
Idw0$X: Idw0.cpp
	${Cpp} ${CcOpt} $^ -o $@
Idw1$X: Idw1.cpp
	${Cpp} ${CcOpt} $^ -o $@
test: Idw0$X Idw1$X
	echo "Idw0 Demo"
	./Idw0$X
	echo "Idw1 Demo"
	./Idw1$X
clean:
	rm -f Idw0$X Idw1$X
