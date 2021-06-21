// ALGLIB++
// Based on ALGLIB: Copyright (c) Sergey Bochkanov (ALGLIB project).
// Revisions Copyright (c) Lydia Marie Williamson, Mark Hopkins Consulting
// Source License:
//	This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
//	as published by the Free Software Foundation (www.fsf.org);
//	either version 2 of the License, or (at your option) any later version.
//
//	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//	See the GNU General Public License for more details.
//
//	A copy of the GNU General Public License is available at http://www.fsf.org/licensing/licenses
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
   setglobalthreading(ParTH);
// Perform matrix-matrix product.
   flops = 2*pow(n, 3);
   timeneeded = counter();
   rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 1, 0.0, c, 0, 0);
   timeneeded = counter() - timeneeded;
// Evaluate the performance.
   printf("Performance is %.1f GFLOPS\n", 1.0E-9*flops/timeneeded);
   return 0;
}
