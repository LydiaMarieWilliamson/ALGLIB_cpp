cp ../AlgLibInternal.cpp ../AlgLibInternal.h ../AlgLibMisc.cpp ../AlgLibMisc.h ../Ap.cpp ../Ap.h ../LinAlg.cpp ../LinAlg.h .
echo "Test 1"
g++ -I. -o Demo.out *.cpp && ./Demo.out && rm Demo.out
echo "Test 2"
g++ -I. -o Demo.out -O3 *.cpp && ./Demo.out && rm Demo.out
echo "Test 3"
g++ -I. -o Demo.out -lpthread -O3 -DAE_OS=AE_LINUX -DAE_CPU=AE_INTEL -msse2 *.cpp && ./Demo.out && rm Demo.out
rm AlgLibInternal.cpp AlgLibInternal.h AlgLibMisc.cpp AlgLibMisc.h Ap.cpp Ap.h LinAlg.cpp LinAlg.h
