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
#ModA=Ap KernelsAvx2 KernelsFma KernelsSse2
ModA=Ap
Mods=${ModA} AlgLibInternal AlgLibMisc LinAlg Solvers Optimization Integration Interpolation SpecialFunctions DataAnalysis Statistics DiffEquations FastTransforms
Objs=$(Mods:%=%$O)
ModX=${ModA} AlgLibInternal AlgLibMisc LinAlg Solvers Optimization Integration Interpolation SpecialFunctions DataAnalysis Statistics
ObjX=$(ModX:%=%$O)
ModY=${ModA} AlgLibInternal AlgLibMisc LinAlg
ObjY=$(ModY:%=%$O)
ModZ=${ModA} AlgLibInternal AlgLibMisc LinAlg Solvers Optimization Integration Interpolation SpecialFunctions
SrcZ=$(ModZ:%=%.cpp)

all: test
Ap$O: Ap.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
AlgLibInternal$O: AlgLibInternal.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
AlgLibMisc$O: AlgLibMisc.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
DataAnalysis$O: DataAnalysis.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
DiffEquations$O: DiffEquations.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
FastTransforms$O: FastTransforms.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
Integration$O: Integration.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
Interpolation$O: Interpolation.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
#KernelsSse2$O: KernelsSse2.cpp
#	${Cpp} ${CcOpt} -c $^ -o $@
#KernelsAvx2$O: KernelsAvx2.cpp
#	${Cpp} ${CcOpt} -c $^ -o $@
#KernelsFma$O: KernelsFma.cpp
#	${Cpp} ${CcOpt} -c $^ -o $@
LinAlg$O: LinAlg.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
Optimization$O: Optimization.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
Solvers$O: Solvers.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
SpecialFunctions$O: SpecialFunctions.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
Statistics$O: Statistics.cpp
	${Cpp} ${CcOpt} -c $^ -o $@
TestC$O: TestC.cpp
	${Cpp} ${CcOpt} -c $^ -DAE_DEBUG4POSIX -o $@
TestI$O: TestI.cpp
	${Cpp} ${CcOpt} -c $^ -DAE_USE_ALLOC_COUNTER
TestX$O: TestX.cpp
	${Cpp} ${CcOpt} -c $^ -DAE_DEBUG4POSIX -DAE_USE_ALLOC_COUNTER -o $@
TestY$O: TestY.cpp
	${Cpp} ${CcOpt} -c $^ -o $@

TestC$X: ${Objs} TestC$O
	${Cpp} ${CcOpt} $^ ${Libs} -o $@
TestI$X: ${Objs} TestI$O
	${Cpp} ${CcOpt} $^ ${Libs} -o $@
TestX$X: ${ObjX} TestX$O
	${Cpp} ${CcOpt} $^ ${Libs} -o $@
TestY$X: ${ObjY} TestY$O
	${Cpp} ${CcOpt} $^ ${Libs} -o $@
TestZ$X: ${SrcZ} TestZ.cpp
	${Cpp} ${CcOpt} $^ ${Libs} -DAE_DEBUG4POSIX -DAE_USE_ALLOC_COUNTER -DAE_NO_EXCEPTIONS -DAE_THREADING=NonTH -o $@
test:	TestI$X TestY$X TestX$X TestC$X TestZ$X
	echo "TestI: API Interface"
	./TestI$X
	echo "TestY: C++ Wrapper Test"
	./TestY$X
	echo "TestZ: Exception-Free Wrapper Test"
	./TestZ$X
	echo "TestX: Speed and Consistency"
	./TestX$X
	echo "TestC: Internal Routines"
	./TestC$X
clean:
	rm -f Test{I,X,Y,Z,C}$O ${Objs}
clobber: clean
	rm -f Test{I,X,Y,Z,C}$X

## Source - Header dependencies:
Ap.cpp: Ap.h
AlgLibInternal.cpp: AlgLibInternal.h
AlgLibMisc.cpp: AlgLibMisc.h
DataAnalysis.cpp: DataAnalysis.h
DiffEquations.cpp: DiffEquations.h
FastTransforms.cpp: FastTransforms.h
Integration.cpp: Integration.h
Interpolation.cpp: Interpolation.h
#KernelsAvx2.cpp: KerrnelsAvx2.h
#KernelsFma.cpp: KernelsFma.h
#KernelsSse2.cpp: KernelsSse2.h
LinAlg.cpp: LinAlg.h
Optimization.cpp: Optimization.h
Solvers.cpp: Solvers.h
SpecialFunctions.cpp: SpecialFunctions.h
Statistics.cpp: Statistics.h
TestC.cpp TestI.cpp:	DataAnalysis.h DiffEquations.h FastTransforms.h Interpolation.h
TestX.cpp:	DataAnalysis.h Interpolation.h
TestY.cpp:	LinAlg.h
TestZ.cpp:	Interpolation.h

## Header - Header dependencies:
## KernelsAvx2, KernelsFma, KernelsSse2 -> Ap
## Optimization -> Solvers -> LinAlg -> AlgLibMisc -> AlgLibInternal -> Ap
## DiffEquations, FastTransforms -> AlgLibInternal
## SpecialFunctions -> AlgLibMisc
## Statistics, Integration -> SpecialFunctions, LinAlg
## DataAnalysis -> Statistics
## Interpolation -> Integration
## DataAnalysis, Interpolation -> Optimization
## TestC, TestI -> DiffEquations, FastTransforms, DataAnalysis, Interpolation
## TestX -> DataAnalysis, Interpolation
## TestY -> LinAlg
## TestZ -> Interpolation
AlgLibInternal.h:	Ap.h
#KernelsAvx2.h KernelsFma.h KernelsSse2.h: Ap.h
AlgLibMisc.h DiffEquations.h FastTransforms.h:	AlgLibInternal.h
LinAlg.h SpecialFunctions.h:	AlgLibMisc.h
Integration.h Statistics.h Solvers.h:	LinAlg.h
Integration.h Statistics.h:	SpecialFunctions.h
Optimization.h:		Solvers.h
DataAnalysis.h Interpolation.h:	Optimization.h
DataAnalysis.h:		Statistics.h
Interpolation.h:	Integration.h
