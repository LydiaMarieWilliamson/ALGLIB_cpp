#include <sys/time.h>
#include "LinAlg.h"

using namespace alglib;

double counter() {
   struct timeval now;
   alglib_impl::ae_int64_t r, v;
   gettimeofday(&now, NULL);
   v = now.tv_sec;
   r = v*1000;
   v = now.tv_usec/1000;
   r += v;
   return 0.001*r;
}

int main() {
   real_2d_array a, b, c;
   int n = 2000;
   int i, j;
   double timeneeded, flops;
// Initialize the arrays.
   a.setlength(n, n);
   b.setlength(n, n);
   c.setlength(n, n);
   for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
      a[i][j] = randomreal() - 0.5;
      b[i][j] = randomreal() - 0.5;
      c[i][j] = 0.0;
   }
// Set global threading settings (applied to all ALGLIB/C++ functions);
// the default is to perform serial computations, unless parallel execution is activated.
// Parallel execution tries to utilize all cores; this behavior can be changed with setnworkers().
   setglobalthreading(parallel);
// Perform matrix-matrix product.
   flops = 2*pow(n, 3);
   timeneeded = counter();
   rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 1, 0.0, c, 0, 0);
   timeneeded = counter() - timeneeded;
// Evaluate the performance.
   printf("Performance is %.1f GFLOPS\n", 1.0E-9*flops/timeneeded);
   return 0;
}
