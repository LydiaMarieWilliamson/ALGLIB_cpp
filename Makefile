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
Mods=Ap AlgLibInternal AlgLibMisc DataAnalysis DiffEquations FastTransforms Integration Interpolation LinAlg Optimization Solvers SpecialFunctions Statistics
Objs=$(Mods:%=%$O)
Srcs=$(Mods:%=%.cpp)
ModsF=Ap AlgLibInternal AlgLibMisc DataAnalysis LinAlg Optimization Solvers SpecialFunctions Statistics
ObjF=$(ModsF:%=%.o)
ModsY=Ap AlgLibInternal AlgLibMisc LinAlg
ObjY=$(ModsY:%=%.o)

all: test
Ap$O: Ap.cpp
	${Cpp} ${CcOpt} -c Ap.cpp -o Ap$O
AlgLibInternal$O: AlgLibInternal.cpp
	${Cpp} ${CcOpt} -c AlgLibInternal.cpp -o AlgLibInternal$O
AlgLibMisc$O: AlgLibMisc.cpp
	${Cpp} ${CcOpt} -c AlgLibMisc.cpp -o AlgLibMisc$O
DataAnalysis$O: DataAnalysis.cpp
	${Cpp} ${CcOpt} -c DataAnalysis.cpp -o DataAnalysis$O
DiffEquations$O: DiffEquations.cpp
	${Cpp} ${CcOpt} -c DiffEquations.cpp -o DiffEquations$O
FastTransforms$O: FastTransforms.cpp
	${Cpp} ${CcOpt} -c FastTransforms.cpp -o FastTransforms$O
Integration$O: Integration.cpp
	${Cpp} ${CcOpt} -c Integration.cpp -o Integration$O
Interpolation$O: Interpolation.cpp
	${Cpp} ${CcOpt} -c Interpolation.cpp -o Interpolation$O
LinAlg$O: LinAlg.cpp
	${Cpp} ${CcOpt} -c LinAlg.cpp -o LinAlg$O
Optimization$O: Optimization.cpp
	${Cpp} ${CcOpt} -c Optimization.cpp -o Optimization$O
Solvers$O: Solvers.cpp
	${Cpp} ${CcOpt} -c Solvers.cpp -o Solvers$O
SpecialFunctions$O: SpecialFunctions.cpp
	${Cpp} ${CcOpt} -c SpecialFunctions.cpp -o SpecialFunctions$O
Statistics$O: Statistics.cpp
	${Cpp} ${CcOpt} -c Statistics.cpp -o Statistics$O
TestC$O: TestC.cpp
	${Cpp} ${CcOpt} -c TestC.cpp -DAE_DEBUG4POSIX -o TestC$O
TestI$O: TestI.cpp
	${Cpp} ${CcOpt} -c TestI.cpp -DAE_USE_ALLOC_COUNTER
TestX$O: TestX.cpp
	${Cpp} ${CcOpt} -c TestX.cpp -DAE_DEBUG4POSIX -DAE_USE_ALLOC_COUNTER -o TestX$O
TestY$O: TestY.cpp
	${Cpp} ${CcOpt} -c TestY.cpp -o TestY$O

TestC$X: ${Objs} TestC$O
	${Cpp} ${CcOpt} ${Objs} TestC$O ${Libs} -o TestC$X
TestI$X: ${Objs} TestI$O
	${Cpp} ${CcOpt} ${Objs} TestI$O ${Libs} -o TestI$X
TestX$X: ${Objs} TestX$O
	${Cpp} ${CcOpt} ${Objs} TestX$O ${Libs} -o TestX$X
TestY$X: ${ObjY} TestY.o
	${Cpp} ${CcOpt} ${ObjY} TestY$O ${Libs} -o TestY$X
TestZ$X: $(Srcs) TestZ.cpp
	${Cpp} ${CcOpt} ${Srcs} TestZ.cpp ${Libs} -O3 -DAE_OS=AE_LINUX -DAE_DEBUG4POSIX -DAE_USE_ALLOC_COUNTER -DAE_NO_EXCEPTIONS -DAE_THREADING=AE_SERIAL_UNSAFE -o TestZ$X
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
	rm -f Test{I,X,Y,Z,C}$X Test{I,X,Y,Z,C}$O ${Objs}

## Source - Header dependencies
Ap.cpp: Ap.h
AlgLibInternal.cpp: AlgLibInternal.h
AlgLibMisc.cpp: AlgLibMisc.h
DataAnalysis.cpp: DataAnalysis.h
DiffEquations.cpp: DiffEquations.h
FastTransforms.cpp: FastTransforms.h
Integration.cpp: Integration.h
Interpolation.cpp: Interpolation.h
LinAlg.cpp: LinAlg.h
Optimization.cpp: Optimization.h
Solvers.cpp: Solvers.h
SpecialFunctions.cpp: SpecialFunctions.h
Statistics.cpp: Statistics.h

## Header - Header dependencies
## Optimization -> Solvers -> LinAlg -> AlgLibMisc -> AlgLibInternal -> Ap
## DiffEquations, FastTransforms -> AlgLibInternal
## SpecialFunctions -> AlgLibMisc
## Statistics, Interpolation -> SpecialFunctions, LinAlg
## DataAnalysis, Interpolation -> Optimization
## TestC, TestI, TestX, TestZ -> DiffEquations, FastTransforms, DataAnalysis, Interpolation
## TestY -> LinAlg
AlgLibInternal.h:	Ap.h
AlgLibMisc.h DiffEquations.h FastTransforms.h:	AlgLibInternal.h
LinAlg.h SpecialFunctions.h:	AlgLibMisc.h
Integration.h Statistics.h Solvers.h:	LinAlg.h
Integration.h Statistics.h:	SpecialFunctions.h
Optimization.h:		Solvers.h
DataAnalysis.h Interpolation.h:	Optimization.h
DataAnalysis.h:		Statistics.h
Interpolation.h:	Integration.h
TestC.cpp TestI.cpp TestX.cpp TestZ.cpp:	DataAnalysis.h DiffEquations.h FastTransforms.h Interpolation.h
TestY.cpp:	LinAlg.h
