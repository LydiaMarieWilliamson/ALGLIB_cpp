#include "LinAlg.h"

using namespace alglib;

int main() {
   real_2d_array a("[[1]]");
   spdmatrixcholesky(a, 1, true);
   return 0;
}
