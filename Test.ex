echo "TestI: API Interface"
TestI: API Interface
./TestI
CPUID: sse2 avx2 fma
OS: POSIX
C++ tests. Please wait...
Allocation counter activated...
0/151
50/151
100/151
150/151
151/151
Allocation counter checked... OK
echo "TestY: C++ Wrapper Test"
TestY: C++ Wrapper Test
./TestY
echo "TestZ: Exception-Free Wrapper Test"
TestZ: Exception-Free Wrapper Test
./TestZ
Test exception-free error handling:
Allocation counter activated...
* default flag value          OK
* 1D arrays                   OK
* 2D arrays                   OK
* ALGLIB objects              OK
* ALGLIB functions            OK
Allocation counter checked... OK
echo "TestX: Speed and Consistency"
TestX: Speed and Consistency
./TestX
System:
* cores count                  1
Allocation counter activated...
Basic functions:
* 1D arrays                   OK
* 2D arrays                   OK
* CSV support                 OK
* Serialization (kd-tree)     OK
* Serialization (RBF)         OK
* Progress/termination (RBF)  OK
* Exceptions in constructors  OK
SMP settings vs GEMM speedup:
* test skipped (no SMP)       ??
Issues:
* issue 505                   OK
* issue 478                   OK
* issue 528                   OK
* issue 591                   OK
* issue 594                   OK
* issue 764                   OK
* issue 813                   OK
* issue 824                   OK
Performance:
* RGEMM-SEQ-16   (MFLOPS)   7519
* RGEMM-MTN-16              1.0x
* RGEMM-SEQ-32   (MFLOPS)  12049
* RGEMM-MTN-32              1.0x
* RGEMM-SEQ-64   (MFLOPS)  11909
* RGEMM-MTN-64              1.0x
* RGEMM-SEQ-1024 (MFLOPS)   9587
* RGEMM-MTN-1024            1.0x
Allocation counter checked... OK
echo "TestC: Internal Routines"
TestC: Internal Routines
./TestC
SEED: 1635369021
COMPILER: GCC
HARDWARE: 64-bit
HARDWARE: little-endian
CPU:   Intel
CORES: 1 (serial version)
LIBS:  
CPUID: sse2 avx2 fma
OS: POSIX
TESTING MODE: single core
ablasf                           OK
hqrnd                            OK
ablas                            OK
hblas                            OK
creflections                     OK
sblas                            OK
ortfac                           OK
matgen                           OK
tsort                            OK
sparse                           OK
blas                             OK
evd                              OK
trfac                            OK
polynomialsolver                 OK
bdsvd                            OK
svd                              OK
trlinsolve                       OK
safesolve                        OK
rcond                            OK
xblas                            OK
directdensesolvers               OK
directsparsesolvers              OK
fbls                             OK
iterativesparse                  OK
lincg                            OK
normestimator                    OK
linlsqr                          OK
linmin                           OK
nleq                             OK
matinv                           OK
optserv                          OK
minlbfgs                         OK
cqmodels                         OK
snnls                            OK
sactivesets                      OK
minbleic                         OK
minqp                            OK
minlm                            OK
mincg                            OK
minlp                            OK
minnlc                           OK
minns                            OK
minbc                            OK
nearestneighbor                  OK
odesolver                        OK
inverseupdate                    OK
schur                            OK
spdgevd                          OK
gammafunc                        OK
gq                               OK
gkq                              OK
autogk                           OK
normaldistr                      OK
basestat                         OK
wsr                              OK
mannwhitneyu                     OK
stest                            OK
studentttests                    OK
ratint                           OK
idw                              OK
polint                           OK
spline1d                         OK
lsfit                            OK
fitsphere                        OK
parametric                       OK
spline2d                         OK
spline3d                         OK
rbf                              OK
fft                              OK
fht                              OK
conv                             OK
corr                             OK
chebyshev                        OK
hermite                          OK
legendre                         OK
laguerre                         OK
pca                              OK
bdss                             OK
mlpbase                          OK
mlpe                             OK
clustering                       OK
dforest                          OK
linreg                           OK
filters                          OK
ssa                              OK
lda                              OK
mcpd                             OK
knn                              OK
mlptrain                         OK
alglibbasics                     OK
Done in 83 seconds
