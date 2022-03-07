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
#define InAlgLib
#include "FastTransforms.h"

// === FFT Package ===
// Depends on: (AlgLibInternal) FTBASE
namespace alglib_impl {
// 1-dimensional complex FFT.
//
// Array size N may be arbitrary number (composite or prime).  Composite  N's
// are handled with cache-oblivious variation of  a  Cooley-Tukey  algorithm.
// Small prime-factors are transformed using hard coded  codelets (similar to
// FFTW codelets, but without low-level  optimization),  large  prime-factors
// are handled with Bluestein's algorithm.
//
// Fastests transforms are for smooth N's (prime factors are 2, 3,  5  only),
// most fast for powers of 2. When N have prime factors  larger  than  these,
// but orders of magnitude smaller than N, computations will be about 4 times
// slower than for nearby highly composite N's. When N itself is prime, speed
// will be 6 times lower.
//
// Algorithm has O(N*logN) complexity for any N (composite or prime).
//
// Inputs:
//     A   -   array[0..N-1] - complex function to be transformed
//     N   -   problem size
//
// Outputs:
//     A   -   DFT of a input array, array[0..N-1]
//             A_out[j] = SUM(A_in[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)
// ALGLIB: Copyright 29.05.2009 by Sergey Bochkanov
// API: void fftc1d(complex_1d_array &a, const ae_int_t n);
// API: void fftc1d(complex_1d_array &a);
void fftc1d(CVector *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewObj(fasttransformplan, plan);
   NewVector(buf, 0, DT_REAL);
   ae_assert(n > 0, "FFTC1D: incorrect N!");
   ae_assert(a->cnt >= n, "FFTC1D: Length(A)<N!");
   ae_assert(isfinitecvector(a, n), "FFTC1D: A contains infinite or NAN values!");
// Special case: N=1, FFT is just identity transform.
// After this block we assume that N is strictly greater than 1.
   if (n == 1) {
      ae_frame_leave();
      return;
   }
// convert input array to the more convinient format
   ae_vector_set_length(&buf, 2 * n);
   for (i = 0; i < n; i++) {
      buf.xR[2 * i] = a->xC[i].x;
      buf.xR[2 * i + 1] = a->xC[i].y;
   }
// Generate plan and execute it.
//
// Plan is a combination of a successive factorizations of N and
// precomputed data. It is much like a FFTW plan, but is not stored
// between subroutine calls and is much simpler.
   ftcomplexfftplan(n, 1, &plan);
   ftapplyplan(&plan, &buf, 0, 1);
// result
   for (i = 0; i < n; i++) {
      a->xC[i] = complex_from_d(buf.xR[2 * i], buf.xR[2 * i + 1]);
   }
   ae_frame_leave();
}

// 1-dimensional complex inverse FFT.
//
// Array size N may be arbitrary number (composite or prime).  Algorithm  has
// O(N*logN) complexity for any N (composite or prime).
//
// See FFTC1D() description for more information about algorithm performance.
//
// Inputs:
//     A   -   array[0..N-1] - complex array to be transformed
//     N   -   problem size
//
// Outputs:
//     A   -   inverse DFT of a input array, array[0..N-1]
//             A_out[j] = SUM(A_in[k]/N*exp(+2*pi*sqrt(-1)*j*k/N), k = 0..N-1)
// ALGLIB: Copyright 29.05.2009 by Sergey Bochkanov
// API: void fftc1dinv(complex_1d_array &a, const ae_int_t n);
// API: void fftc1dinv(complex_1d_array &a);
void fftc1dinv(CVector *a, ae_int_t n) {
   ae_int_t i;
   ae_assert(n > 0, "FFTC1DInv: incorrect N!");
   ae_assert(a->cnt >= n, "FFTC1DInv: Length(A)<N!");
   ae_assert(isfinitecvector(a, n), "FFTC1DInv: A contains infinite or NAN values!");
// Inverse DFT can be expressed in terms of the DFT as
//
//     invfft(x) = fft(x')'/N
//
// here x' means conj(x).
   for (i = 0; i < n; i++) {
      a->xC[i].y = -a->xC[i].y;
   }
   fftc1d(a, n);
   for (i = 0; i < n; i++) {
      a->xC[i] = complex_from_d(a->xC[i].x / n, -a->xC[i].y / n);
   }
}

// 1-dimensional real FFT.
//
// Algorithm has O(N*logN) complexity for any N (composite or prime).
//
// Inputs:
//     A   -   array[0..N-1] - real function to be transformed
//     N   -   problem size
//
// Outputs:
//     F   -   DFT of a input array, array[0..N-1]
//             F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)
//
// NOTE:
//     F[] satisfies symmetry property F[k] = conj(F[N-k]),  so just one half
// of  array  is  usually needed. But for convinience subroutine returns full
// complex array (with frequencies above N/2), so its result may be  used  by
// other FFT-related subroutines.
// ALGLIB: Copyright 01.06.2009 by Sergey Bochkanov
// API: void fftr1d(const real_1d_array &a, const ae_int_t n, complex_1d_array &f);
// API: void fftr1d(const real_1d_array &a, complex_1d_array &f);
void fftr1d(RVector *a, ae_int_t n, CVector *f) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_int_t idx;
   complex hn;
   complex hmnc;
   complex v;
   ae_frame_make(&_frame_block);
   SetVector(f);
   NewVector(buf, 0, DT_REAL);
   NewObj(fasttransformplan, plan);
   ae_assert(n > 0, "FFTR1D: incorrect N!");
   ae_assert(a->cnt >= n, "FFTR1D: Length(A)<N!");
   ae_assert(isfinitevector(a, n), "FFTR1D: A contains infinite or NAN values!");
// Special cases:
// * N=1, FFT is just identity transform.
// * N=2, FFT is simple too
//
// After this block we assume that N is strictly greater than 2
   if (n == 1) {
      ae_vector_set_length(f, 1);
      f->xC[0] = complex_from_d(a->xR[0]);
      ae_frame_leave();
      return;
   }
   if (n == 2) {
      ae_vector_set_length(f, 2);
      f->xC[0] = complex_from_d(a->xR[0] + a->xR[1]);
      f->xC[1] = complex_from_d(a->xR[0] - a->xR[1]);
      ae_frame_leave();
      return;
   }
// Choose between odd-size and even-size FFTs
   if (n % 2 == 0) {
   // even-size real FFT, use reduction to the complex task
      n2 = n / 2;
      ae_vector_set_length(&buf, n);
      ae_v_move(buf.xR, 1, a->xR, 1, n);
      ftcomplexfftplan(n2, 1, &plan);
      ftapplyplan(&plan, &buf, 0, 1);
      ae_vector_set_length(f, n);
      for (i = 0; i <= n2; i++) {
         idx = 2 * (i % n2);
         hn = complex_from_d(buf.xR[idx], buf.xR[idx + 1]);
         idx = 2 * ((n2 - i) % n2);
         hmnc = complex_from_d(buf.xR[idx], -buf.xR[idx + 1]);
         v = complex_from_d(-sin(-2 * pi * i / n), cos(-2 * pi * i / n));
         f->xC[i] = ae_c_sub(ae_c_add(hn, hmnc), ae_c_mul(v, ae_c_sub(hn, hmnc)));
         f->xC[i].x *= 0.5;
         f->xC[i].y *= 0.5;
      }
      for (i = n2 + 1; i < n; i++) {
         f->xC[i] = conj(f->xC[n - i]);
      }
   } else {
   // use complex FFT
      ae_vector_set_length(f, n);
      for (i = 0; i < n; i++) {
         f->xC[i] = complex_from_d(a->xR[i]);
      }
      fftc1d(f, n);
   }
   ae_frame_leave();
}

// 1-dimensional real inverse FFT.
//
// Algorithm has O(N*logN) complexity for any N (composite or prime).
//
// Inputs:
//     F   -   array[0..floor(N/2)] - frequencies from forward real FFT
//     N   -   problem size
//
// Outputs:
//     A   -   inverse DFT of a input array, array[0..N-1]
//
// NOTE:
//     F[] should satisfy symmetry property F[k] = conj(F[N-k]), so just  one
// half of frequencies array is needed - elements from 0 to floor(N/2).  F[0]
// is ALWAYS real. If N is even F[floor(N/2)] is real too. If N is odd,  then
// F[floor(N/2)] has no special properties.
//
// Relying on properties noted above, FFTR1DInv subroutine uses only elements
// from 0th to floor(N/2)-th. It ignores imaginary part of F[0],  and in case
// N is even it ignores imaginary part of F[floor(N/2)] too.
//
// When you call this function using full arguments list - "FFTR1DInv(F,N,A)"
// - you can pass either either frequencies array with N elements or  reduced
// array with roughly N/2 elements - subroutine will  successfully  transform
// both.
//
// If you call this function using reduced arguments list -  "FFTR1DInv(F,A)"
// - you must pass FULL array with N elements (although higher  N/2 are still
// not used) because array size is used to automatically determine FFT length
// ALGLIB: Copyright 01.06.2009 by Sergey Bochkanov
// API: void fftr1dinv(const complex_1d_array &f, const ae_int_t n, real_1d_array &a);
// API: void fftr1dinv(const complex_1d_array &f, real_1d_array &a);
void fftr1dinv(CVector *f, ae_int_t n, RVector *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   SetVector(a);
   NewVector(h, 0, DT_REAL);
   NewVector(fh, 0, DT_COMPLEX);
   ae_assert(n > 0, "FFTR1DInv: incorrect N!");
   ae_int_t Nq = n / 2, Nr = n % 2;
   ae_assert(f->cnt >= Nq + 1, "FFTR1DInv: Length(F)<Floor(N/2)+1!");
   ae_assert(isfinite(f->xC[0].x), "FFTR1DInv: F contains infinite or NAN values!");
   for (i = 1; i < Nq; i++) {
      ae_assert(isfinite(f->xC[i].x) && isfinite(f->xC[i].y), "FFTR1DInv: F contains infinite or NAN values!");
   }
   ae_assert(isfinite(f->xC[Nq].x), "FFTR1DInv: F contains infinite or NAN values!");
   if (Nr) {
      ae_assert(isfinite(f->xC[Nq].y), "FFTR1DInv: F contains infinite or NAN values!");
   }
// Special case: N=1, FFT is just identity transform.
// After this block we assume that N is strictly greater than 1.
   if (n == 1) {
      ae_vector_set_length(a, 1);
      a->xR[0] = f->xC[0].x;
      ae_frame_leave();
      return;
   }
// inverse real FFT is reduced to the inverse real FHT,
// which is reduced to the forward real FHT,
// which is reduced to the forward real FFT.
//
// Don't worry, it is really compact and efficient reduction :)
   ae_vector_set_length(&h, n);
   ae_vector_set_length(a, n);
   h.xR[0] = f->xC[0].x;
   for (i = 1; i < Nq; i++) {
      h.xR[i] = f->xC[i].x - f->xC[i].y;
      h.xR[n - i] = f->xC[i].x + f->xC[i].y;
   }
   if (!Nr) {
      h.xR[Nq] = f->xC[Nq].x;
   } else {
      h.xR[Nq] = f->xC[Nq].x - f->xC[Nq].y;
      h.xR[Nq + 1] = f->xC[Nq].x + f->xC[Nq].y;
   }
   fftr1d(&h, n, &fh);
   for (i = 0; i < n; i++) {
      a->xR[i] = (fh.xC[i].x - fh.xC[i].y) / n;
   }
   ae_frame_leave();
}

// Internal subroutine. Never call it directly!
// ALGLIB: Copyright 01.06.2009 by Sergey Bochkanov
void fftr1dinternaleven(RVector *a, ae_int_t n, RVector *buf, fasttransformplan *plan) {
   double x;
   double y;
   ae_int_t i;
   ae_int_t n2;
   ae_int_t idx;
   complex hn;
   complex hmnc;
   complex v;
   ae_assert(n > 0 && n % 2 == 0, "FFTR1DEvenInplace: incorrect N!");
// Special cases:
// * N=2
//
// After this block we assume that N is strictly greater than 2
   if (n == 2) {
      x = a->xR[0] + a->xR[1];
      y = a->xR[0] - a->xR[1];
      a->xR[0] = x;
      a->xR[1] = y;
      return;
   }
// even-size real FFT, use reduction to the complex task
   n2 = n / 2;
   ae_v_move(buf->xR, 1, a->xR, 1, n);
   ftapplyplan(plan, buf, 0, 1);
   a->xR[0] = buf->xR[0] + buf->xR[1];
   for (i = 1; i < n2; i++) {
      idx = 2 * (i % n2);
      hn = complex_from_d(buf->xR[idx], buf->xR[idx + 1]);
      idx = 2 * (n2 - i);
      hmnc = complex_from_d(buf->xR[idx], -buf->xR[idx + 1]);
      v = complex_from_d(-sin(-2 * pi * i / n), cos(-2 * pi * i / n));
      v = ae_c_sub(ae_c_add(hn, hmnc), ae_c_mul(v, ae_c_sub(hn, hmnc)));
      a->xR[2 * i] = 0.5 * v.x;
      a->xR[2 * i + 1] = 0.5 * v.y;
   }
   a->xR[1] = buf->xR[0] - buf->xR[1];
}

// Internal subroutine. Never call it directly!
// ALGLIB: Copyright 01.06.2009 by Sergey Bochkanov
void fftr1dinvinternaleven(RVector *a, ae_int_t n, RVector *buf, fasttransformplan *plan) {
   double x;
   double y;
   double t;
   ae_int_t i;
   ae_int_t n2;
   ae_assert(n > 0 && n % 2 == 0, "FFTR1DInvInternalEven: incorrect N!");
// Special cases:
// * N=2
//
// After this block we assume that N is strictly greater than 2
   if (n == 2) {
      x = 0.5 * (a->xR[0] + a->xR[1]);
      y = 0.5 * (a->xR[0] - a->xR[1]);
      a->xR[0] = x;
      a->xR[1] = y;
      return;
   }
// inverse real FFT is reduced to the inverse real FHT,
// which is reduced to the forward real FHT,
// which is reduced to the forward real FFT.
//
// Don't worry, it is really compact and efficient reduction :)
   n2 = n / 2;
   buf->xR[0] = a->xR[0];
   for (i = 1; i < n2; i++) {
      x = a->xR[2 * i];
      y = a->xR[2 * i + 1];
      buf->xR[i] = x - y;
      buf->xR[n - i] = x + y;
   }
   buf->xR[n2] = a->xR[1];
   fftr1dinternaleven(buf, n, a, plan);
   a->xR[0] = buf->xR[0] / n;
   t = 1.0 / n;
   for (i = 1; i < n2; i++) {
      x = buf->xR[2 * i];
      y = buf->xR[2 * i + 1];
      a->xR[i] = t * (x - y);
      a->xR[n - i] = t * (x + y);
   }
   a->xR[n2] = buf->xR[1] / n;
}
} // end of namespace alglib_impl

namespace alglib {
void fftc1d(complex_1d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftc1d(ConstT(ae_vector, a), n);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void fftc1d(complex_1d_array &a) {
   ae_int_t n = a.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftc1d(ConstT(ae_vector, a), n);
   alglib_impl::ae_state_clear();
}
#endif

void fftc1dinv(complex_1d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftc1dinv(ConstT(ae_vector, a), n);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void fftc1dinv(complex_1d_array &a) {
   ae_int_t n = a.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftc1dinv(ConstT(ae_vector, a), n);
   alglib_impl::ae_state_clear();
}
#endif

void fftr1d(const real_1d_array &a, const ae_int_t n, complex_1d_array &f) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftr1d(ConstT(ae_vector, a), n, ConstT(ae_vector, f));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void fftr1d(const real_1d_array &a, complex_1d_array &f) {
   ae_int_t n = a.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftr1d(ConstT(ae_vector, a), n, ConstT(ae_vector, f));
   alglib_impl::ae_state_clear();
}
#endif

void fftr1dinv(const complex_1d_array &f, const ae_int_t n, real_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftr1dinv(ConstT(ae_vector, f), n, ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void fftr1dinv(const complex_1d_array &f, real_1d_array &a) {
   ae_int_t n = f.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fftr1dinv(ConstT(ae_vector, f), n, ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}
#endif
} // end of namespace alglib

// === FHT Package ===
// Depends on: FFT
namespace alglib_impl {
// 1-dimensional Fast Hartley Transform.
//
// Algorithm has O(N*logN) complexity for any N (composite or prime).
//
// Inputs:
//     A   -   array[0..N-1] - real function to be transformed
//     N   -   problem size
//
// Outputs:
//     A   -   FHT of a input array, array[0..N-1],
//             A_out[k] = sum(A_in[j]*(cos(2*pi*j*k/N)+sin(2*pi*j*k/N)), j=0..N-1)
// ALGLIB: Copyright 04.06.2009 by Sergey Bochkanov
// API: void fhtr1d(real_1d_array &a, const ae_int_t n);
void fhtr1d(RVector *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(fa, 0, DT_COMPLEX);
   ae_assert(n > 0, "FHTR1D: incorrect N!");
// Special case: N=1, FHT is just identity transform.
// After this block we assume that N is strictly greater than 1.
   if (n == 1) {
      ae_frame_leave();
      return;
   }
// Reduce FHt to real FFT
   fftr1d(a, n, &fa);
   for (i = 0; i < n; i++) {
      a->xR[i] = fa.xC[i].x - fa.xC[i].y;
   }
   ae_frame_leave();
}

// 1-dimensional inverse FHT.
//
// Algorithm has O(N*logN) complexity for any N (composite or prime).
//
// Inputs:
//     A   -   array[0..N-1] - complex array to be transformed
//     N   -   problem size
//
// Outputs:
//     A   -   inverse FHT of a input array, array[0..N-1]
// ALGLIB: Copyright 29.05.2009 by Sergey Bochkanov
// API: void fhtr1dinv(real_1d_array &a, const ae_int_t n);
void fhtr1dinv(RVector *a, ae_int_t n) {
   ae_int_t i;
   ae_assert(n > 0, "FHTR1DInv: incorrect N!");
// Special case: N=1, iFHT is just identity transform.
// After this block we assume that N is strictly greater than 1.
   if (n == 1) {
      return;
   }
// Inverse FHT can be expressed in terms of the FHT as
//
//     invfht(x) = fht(x)/N
   fhtr1d(a, n);
   for (i = 0; i < n; i++) {
      a->xR[i] /= n;
   }
}
} // end of namespace alglib_impl

namespace alglib {
void fhtr1d(real_1d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fhtr1d(ConstT(ae_vector, a), n);
   alglib_impl::ae_state_clear();
}

void fhtr1dinv(real_1d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fhtr1dinv(ConstT(ae_vector, a), n);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === CONV Package ===
// Depends on: FFT
namespace alglib_impl {
// 1-dimensional complex convolution.
//
// Extended subroutine which allows to choose convolution algorithm.
// Intended for internal use, ALGLIB users should call ConvC1D()/ConvC1DCircular().
//
// Inputs:
//     A   -   array[0..M-1] - complex function to be transformed
//     M   -   problem size
//     B   -   array[0..N-1] - complex function to be transformed
//     N   -   problem size, N <= M
//     Alg -   algorithm type:
//             *-2     auto-select Q for overlap-add
//             *-1     auto-select algorithm and parameters
//             * 0     straightforward formula for small N's
//             * 1     general FFT-based code
//             * 2     overlap-add with length Q
//     Q   -   length for overlap-add
//
// Outputs:
//     R   -   convolution: A*B. array[0..N+M-1].
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
void convc1dx(CVector *a, ae_int_t m, CVector *b, ae_int_t n, bool circular, ae_int_t alg, ae_int_t q, CVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t p;
   ae_int_t ptotal;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t j1;
   ae_int_t j2;
   complex v;
   double ax;
   double ay;
   double bx;
   double by;
   double t;
   double tx;
   double ty;
   double flopcand;
   double flopbest;
   ae_int_t algbest;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(bbuf, 0, DT_COMPLEX);
   NewObj(fasttransformplan, plan);
   NewVector(buf, 0, DT_REAL);
   NewVector(buf2, 0, DT_REAL);
   ae_assert(n > 0 && m > 0, "ConvC1DX: incorrect N or M!");
   ae_assert(n <= m, "ConvC1DX: N<M assumption is false!");
// Auto-select
   if (alg == -1 || alg == -2) {
   // Initial candidate: straightforward implementation.
   //
   // If we want to use auto-fitted overlap-add,
   // flop count is initialized by large real number - to force
   // another algorithm selection
      algbest = 0;
      if (alg == -1) {
         flopbest = (double)(2 * m * n);
      } else {
         flopbest = maxrealnumber;
      }
   // Another candidate - generic FFT code
      if (alg == -1) {
         if (circular && ftbaseissmooth(m)) {
         // special code for circular convolution of a sequence with a smooth length
            flopcand = 3 * ftbasegetflopestimate(m) + 6 * m;
            if (flopcand < flopbest) {
               algbest = 1;
               flopbest = flopcand;
            }
         } else {
         // general cyclic/non-cyclic convolution
            p = ftbasefindsmooth(m + n - 1);
            flopcand = 3 * ftbasegetflopestimate(p) + 6 * p;
            if (flopcand < flopbest) {
               algbest = 1;
               flopbest = flopcand;
            }
         }
      }
   // Another candidate - overlap-add
      q = 1;
      ptotal = 1;
      while (ptotal < n) {
         ptotal *= 2;
      }
      while (ptotal < m + n) {
         p = ptotal - n + 1;
         flopcand = iceil((double)m / p) * (2 * ftbasegetflopestimate(ptotal) + 8 * ptotal);
         if (flopcand < flopbest) {
            flopbest = flopcand;
            algbest = 2;
            q = p;
         }
         ptotal *= 2;
      }
      alg = algbest;
      convc1dx(a, m, b, n, circular, alg, q, r);
      ae_frame_leave();
      return;
   }
// straightforward formula for
// circular and non-circular convolutions.
//
// Very simple code, no further comments needed.
   if (alg == 0) {
   // Special case: N=1
      if (n == 1) {
         ae_vector_set_length(r, m);
         v = b->xC[0];
         ae_v_cmovec(r->xC, 1, a->xC, 1, "N", m, v);
         ae_frame_leave();
         return;
      }
   // use straightforward formula
      if (circular) {
      // circular convolution
         ae_vector_set_length(r, m);
         v = b->xC[0];
         ae_v_cmovec(r->xC, 1, a->xC, 1, "N", m, v);
         for (i = 1; i < n; i++) {
            v = b->xC[i];
            i1 = 0;
            i2 = i - 1;
            j1 = m - i;
            j2 = m - 1;
            ae_v_caddc(&r->xC[i1], 1, &a->xC[j1], 1, "N", i2 - i1 + 1, v);
            i1 = i;
            i2 = m - 1;
            j1 = 0;
            j2 = m - i - 1;
            ae_v_caddc(&r->xC[i1], 1, &a->xC[j1], 1, "N", i2 - i1 + 1, v);
         }
      } else {
      // non-circular convolution
         ae_vector_set_length(r, m + n - 1);
         for (i = 0; i < m + n - 1; i++) {
            r->xC[i] = complex_from_i(0);
         }
         for (i = 0; i < n; i++) {
            v = b->xC[i];
            ae_v_caddc(&r->xC[i], 1, a->xC, 1, "N", m, v);
         }
      }
      ae_frame_leave();
      return;
   }
// general FFT-based code for
// circular and non-circular convolutions.
//
// First, if convolution is circular, we test whether M is smooth or not.
// If it is smooth, we just use M-length FFT to calculate convolution.
// If it is not, we calculate non-circular convolution and wrap it arount.
//
// IF convolution is non-circular, we use zero-padding + FFT.
   if (alg == 1) {
      if (circular && ftbaseissmooth(m)) {
      // special code for circular convolution with smooth M
         ftcomplexfftplan(m, 1, &plan);
         ae_vector_set_length(&buf, 2 * m);
         for (i = 0; i < m; i++) {
            buf.xR[2 * i] = a->xC[i].x;
            buf.xR[2 * i + 1] = a->xC[i].y;
         }
         ae_vector_set_length(&buf2, 2 * m);
         for (i = 0; i < n; i++) {
            buf2.xR[2 * i] = b->xC[i].x;
            buf2.xR[2 * i + 1] = b->xC[i].y;
         }
         for (i = n; i < m; i++) {
            buf2.xR[2 * i] = 0.0;
            buf2.xR[2 * i + 1] = 0.0;
         }
         ftapplyplan(&plan, &buf, 0, 1);
         ftapplyplan(&plan, &buf2, 0, 1);
         for (i = 0; i < m; i++) {
            ax = buf.xR[2 * i];
            ay = buf.xR[2 * i + 1];
            bx = buf2.xR[2 * i];
            by = buf2.xR[2 * i + 1];
            tx = ax * bx - ay * by;
            ty = ax * by + ay * bx;
            buf.xR[2 * i] = tx;
            buf.xR[2 * i + 1] = -ty;
         }
         ftapplyplan(&plan, &buf, 0, 1);
         t = 1.0 / m;
         ae_vector_set_length(r, m);
         for (i = 0; i < m; i++) {
            r->xC[i] = complex_from_d(t * buf.xR[2 * i], -t * buf.xR[2 * i + 1]);
         }
      } else {
      // M is non-smooth, general code (circular/non-circular):
      // * first part is the same for circular and non-circular
      //   convolutions. zero padding, FFTs, inverse FFTs
      // * second part differs:
      //   * for non-circular convolution we just copy array
      //   * for circular convolution we add array tail to its head
         p = ftbasefindsmooth(m + n - 1);
         ftcomplexfftplan(p, 1, &plan);
         ae_vector_set_length(&buf, 2 * p);
         for (i = 0; i < m; i++) {
            buf.xR[2 * i] = a->xC[i].x;
            buf.xR[2 * i + 1] = a->xC[i].y;
         }
         for (i = m; i < p; i++) {
            buf.xR[2 * i] = 0.0;
            buf.xR[2 * i + 1] = 0.0;
         }
         ae_vector_set_length(&buf2, 2 * p);
         for (i = 0; i < n; i++) {
            buf2.xR[2 * i] = b->xC[i].x;
            buf2.xR[2 * i + 1] = b->xC[i].y;
         }
         for (i = n; i < p; i++) {
            buf2.xR[2 * i] = 0.0;
            buf2.xR[2 * i + 1] = 0.0;
         }
         ftapplyplan(&plan, &buf, 0, 1);
         ftapplyplan(&plan, &buf2, 0, 1);
         for (i = 0; i < p; i++) {
            ax = buf.xR[2 * i];
            ay = buf.xR[2 * i + 1];
            bx = buf2.xR[2 * i];
            by = buf2.xR[2 * i + 1];
            tx = ax * bx - ay * by;
            ty = ax * by + ay * bx;
            buf.xR[2 * i] = tx;
            buf.xR[2 * i + 1] = -ty;
         }
         ftapplyplan(&plan, &buf, 0, 1);
         t = 1.0 / p;
         if (circular) {
         // circular, add tail to head
            ae_vector_set_length(r, m);
            for (i = 0; i < m; i++) {
               r->xC[i] = complex_from_d(t * buf.xR[2 * i], -t * buf.xR[2 * i + 1]);
            }
            for (i = m; i < m + n - 1; i++) {
               r->xC[i - m].x += t * buf.xR[2 * i];
               r->xC[i - m].y -= t * buf.xR[2 * i + 1];
            }
         } else {
         // non-circular, just copy
            ae_vector_set_length(r, m + n - 1);
            for (i = 0; i < m + n - 1; i++) {
               r->xC[i] = complex_from_d(t * buf.xR[2 * i], -t * buf.xR[2 * i + 1]);
            }
         }
      }
      ae_frame_leave();
      return;
   }
// overlap-add method for
// circular and non-circular convolutions.
//
// First part of code (separate FFTs of input blocks) is the same
// for all types of convolution. Second part (overlapping outputs)
// differs for different types of convolution. We just copy output
// when convolution is non-circular. We wrap it around, if it is
// circular.
   if (alg == 2) {
      ae_vector_set_length(&buf, 2 * (q + n - 1));
   // prepare R
      if (circular) {
         ae_vector_set_length(r, m);
         for (i = 0; i < m; i++) {
            r->xC[i] = complex_from_i(0);
         }
      } else {
         ae_vector_set_length(r, m + n - 1);
         for (i = 0; i < m + n - 1; i++) {
            r->xC[i] = complex_from_i(0);
         }
      }
   // pre-calculated FFT(B)
      ae_vector_set_length(&bbuf, q + n - 1);
      ae_v_cmove(bbuf.xC, 1, b->xC, 1, "N", n);
      for (j = n; j < q + n - 1; j++) {
         bbuf.xC[j] = complex_from_i(0);
      }
      fftc1d(&bbuf, q + n - 1);
   // prepare FFT plan for chunks of A
      ftcomplexfftplan(q + n - 1, 1, &plan);
   // main overlap-add cycle
      i = 0;
      while (i < m) {
         p = imin2(q, m - i);
         for (j = 0; j < p; j++) {
            buf.xR[2 * j] = a->xC[i + j].x;
            buf.xR[2 * j + 1] = a->xC[i + j].y;
         }
         for (j = p; j < q + n - 1; j++) {
            buf.xR[2 * j] = 0.0;
            buf.xR[2 * j + 1] = 0.0;
         }
         ftapplyplan(&plan, &buf, 0, 1);
         for (j = 0; j < q + n - 1; j++) {
            ax = buf.xR[2 * j];
            ay = buf.xR[2 * j + 1];
            bx = bbuf.xC[j].x;
            by = bbuf.xC[j].y;
            tx = ax * bx - ay * by;
            ty = ax * by + ay * bx;
            buf.xR[2 * j] = tx;
            buf.xR[2 * j + 1] = -ty;
         }
         ftapplyplan(&plan, &buf, 0, 1);
         t = 1.0 / (q + n - 1);
         if (circular) {
            j1 = imin2(i + p + n - 2, m - 1) - i;
            j2 = j1 + 1;
         } else {
            j1 = p + n - 2;
            j2 = j1 + 1;
         }
         for (j = 0; j <= j1; j++) {
            r->xC[i + j].x += buf.xR[2 * j] * t;
            r->xC[i + j].y -= buf.xR[2 * j + 1] * t;
         }
         for (j = j2; j < p + n - 1; j++) {
            r->xC[j - j2].x += buf.xR[2 * j] * t;
            r->xC[j - j2].y -= buf.xR[2 * j + 1] * t;
         }
         i += p;
      }
      ae_frame_leave();
      return;
   }
   ae_frame_leave();
}

// 1-dimensional complex convolution.
//
// For given A/B returns conv(A,B) (non-circular). Subroutine can automatically
// choose between three implementations: straightforward O(M*N)  formula  for
// very small N (or M), overlap-add algorithm for  cases  where  max(M,N)  is
// significantly larger than min(M,N), but O(M*N) algorithm is too slow,  and
// general FFT-based formula for cases where two previois algorithms are  too
// slow.
//
// Algorithm has max(M,N)*log(max(M,N)) complexity for any M/N.
//
// Inputs:
//     A   -   array[0..M-1] - complex function to be transformed
//     M   -   problem size
//     B   -   array[0..N-1] - complex function to be transformed
//     N   -   problem size
//
// Outputs:
//     R   -   convolution: A*B. array[0..N+M-2].
//
// NOTE:
//     It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
// functions have non-zero values at negative T's, you  can  still  use  this
// subroutine - just shift its result correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convc1d(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
void convc1d(CVector *a, ae_int_t m, CVector *b, ae_int_t n, CVector *r) {
   SetVector(r);
   ae_assert(n > 0 && m > 0, "ConvC1D: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer that B.
   if (m < n) {
      convc1d(b, n, a, m, r);
      return;
   }
   convc1dx(a, m, b, n, false, -1, 0, r);
}

// 1-dimensional complex non-circular deconvolution (inverse of ConvC1D()).
//
// Algorithm has M*log(M)) complexity for any M (composite or prime).
//
// Inputs:
//     A   -   array[0..M-1] - convolved signal, A = conv(R, B)
//     M   -   convolved signal length
//     B   -   array[0..N-1] - response
//     N   -   response length, N <= M
//
// Outputs:
//     R   -   deconvolved signal. array[0..M-N].
//
// NOTE:
//     deconvolution is unstable process and may result in division  by  zero
// (if your response function is degenerate, i.e. has zero Fourier coefficient).
//
// NOTE:
//     It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
// functions have non-zero values at negative T's, you  can  still  use  this
// subroutine - just shift its result correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convc1dinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
void convc1dinv(CVector *a, ae_int_t m, CVector *b, ae_int_t n, CVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t p;
   complex c1;
   complex c2;
   complex c3;
   double t;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(buf, 0, DT_REAL);
   NewVector(buf2, 0, DT_REAL);
   NewObj(fasttransformplan, plan);
   ae_assert(n > 0 && m > 0 && n <= m, "ConvC1DInv: incorrect N or M!");
   p = ftbasefindsmooth(m);
   ftcomplexfftplan(p, 1, &plan);
   ae_vector_set_length(&buf, 2 * p);
   for (i = 0; i < m; i++) {
      buf.xR[2 * i] = a->xC[i].x;
      buf.xR[2 * i + 1] = a->xC[i].y;
   }
   for (i = m; i < p; i++) {
      buf.xR[2 * i] = 0.0;
      buf.xR[2 * i + 1] = 0.0;
   }
   ae_vector_set_length(&buf2, 2 * p);
   for (i = 0; i < n; i++) {
      buf2.xR[2 * i] = b->xC[i].x;
      buf2.xR[2 * i + 1] = b->xC[i].y;
   }
   for (i = n; i < p; i++) {
      buf2.xR[2 * i] = 0.0;
      buf2.xR[2 * i + 1] = 0.0;
   }
   ftapplyplan(&plan, &buf, 0, 1);
   ftapplyplan(&plan, &buf2, 0, 1);
   for (i = 0; i < p; i++) {
      c1 = complex_from_d(buf.xR[2 * i], buf.xR[2 * i + 1]);
      c2 = complex_from_d(buf2.xR[2 * i], buf2.xR[2 * i + 1]);
      c3 = ae_c_div(c1, c2);
      buf.xR[2 * i] = c3.x;
      buf.xR[2 * i + 1] = -c3.y;
   }
   ftapplyplan(&plan, &buf, 0, 1);
   t = 1.0 / p;
   ae_vector_set_length(r, m - n + 1);
   for (i = 0; i <= m - n; i++) {
      r->xC[i] = complex_from_d(t * buf.xR[2 * i], -t * buf.xR[2 * i + 1]);
   }
   ae_frame_leave();
}

// 1-dimensional circular complex convolution.
//
// For given S/R returns conv(S,R) (circular). Algorithm has linearithmic
// complexity for any M/N.
//
// IMPORTANT:  normal convolution is commutative,  i.e.   it  is symmetric  -
// conv(A,B)=conv(B,A).  Cyclic convolution IS NOT.  One function - S - is  a
// signal,  periodic function, and another - R - is a response,  non-periodic
// function with limited length.
//
// Inputs:
//     S   -   array[0..M-1] - complex periodic signal
//     M   -   problem size
//     B   -   array[0..N-1] - complex non-periodic response
//     N   -   problem size
//
// Outputs:
//     R   -   convolution: A*B. array[0..M-1].
//
// NOTE:
//     It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
// negative T's, you can still use this subroutine - just  shift  its  result
// correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convc1dcircular(const complex_1d_array &s, const ae_int_t m, const complex_1d_array &r, const ae_int_t n, complex_1d_array &c);
void convc1dcircular(CVector *s, ae_int_t m, CVector *r, ae_int_t n, CVector *c) {
   ae_frame _frame_block;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t j2;
   ae_frame_make(&_frame_block);
   SetVector(c);
   NewVector(buf, 0, DT_COMPLEX);
   ae_assert(n > 0 && m > 0, "ConvC1DCircular: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer (at least - not shorter) that B.
   if (m < n) {
      ae_vector_set_length(&buf, m);
      for (i1 = 0; i1 < m; i1++) {
         buf.xC[i1] = complex_from_i(0);
      }
      i1 = 0;
      while (i1 < n) {
         i2 = imin2(i1 + m - 1, n - 1);
         j2 = i2 - i1;
         ae_v_cadd(buf.xC, 1, &r->xC[i1], 1, "N", j2 + 1);
         i1 += m;
      }
      convc1dcircular(s, m, &buf, m, c);
      ae_frame_leave();
      return;
   }
   convc1dx(s, m, r, n, true, -1, 0, c);
   ae_frame_leave();
}

// 1-dimensional circular complex deconvolution (inverse of ConvC1DCircular()).
//
// Algorithm has M*log(M)) complexity for any M (composite or prime).
//
// Inputs:
//     A   -   array[0..M-1] - convolved periodic signal, A = conv(R, B)
//     M   -   convolved signal length
//     B   -   array[0..N-1] - non-periodic response
//     N   -   response length
//
// Outputs:
//     R   -   deconvolved signal. array[0..M-1].
//
// NOTE:
//     deconvolution is unstable process and may result in division  by  zero
// (if your response function is degenerate, i.e. has zero Fourier coefficient).
//
// NOTE:
//     It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
// negative T's, you can still use this subroutine - just  shift  its  result
// correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convc1dcircularinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
void convc1dcircularinv(CVector *a, ae_int_t m, CVector *b, ae_int_t n, CVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t j2;
   complex c1;
   complex c2;
   complex c3;
   double t;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(buf, 0, DT_REAL);
   NewVector(buf2, 0, DT_REAL);
   NewVector(cbuf, 0, DT_COMPLEX);
   NewObj(fasttransformplan, plan);
   ae_assert(n > 0 && m > 0, "ConvC1DCircularInv: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer (at least - not shorter) that B.
   if (m < n) {
      ae_vector_set_length(&cbuf, m);
      for (i = 0; i < m; i++) {
         cbuf.xC[i] = complex_from_i(0);
      }
      i1 = 0;
      while (i1 < n) {
         i2 = imin2(i1 + m - 1, n - 1);
         j2 = i2 - i1;
         ae_v_cadd(cbuf.xC, 1, &b->xC[i1], 1, "N", j2 + 1);
         i1 += m;
      }
      convc1dcircularinv(a, m, &cbuf, m, r);
      ae_frame_leave();
      return;
   }
// Task is normalized
   ftcomplexfftplan(m, 1, &plan);
   ae_vector_set_length(&buf, 2 * m);
   for (i = 0; i < m; i++) {
      buf.xR[2 * i] = a->xC[i].x;
      buf.xR[2 * i + 1] = a->xC[i].y;
   }
   ae_vector_set_length(&buf2, 2 * m);
   for (i = 0; i < n; i++) {
      buf2.xR[2 * i] = b->xC[i].x;
      buf2.xR[2 * i + 1] = b->xC[i].y;
   }
   for (i = n; i < m; i++) {
      buf2.xR[2 * i] = 0.0;
      buf2.xR[2 * i + 1] = 0.0;
   }
   ftapplyplan(&plan, &buf, 0, 1);
   ftapplyplan(&plan, &buf2, 0, 1);
   for (i = 0; i < m; i++) {
      c1 = complex_from_d(buf.xR[2 * i], buf.xR[2 * i + 1]);
      c2 = complex_from_d(buf2.xR[2 * i], buf2.xR[2 * i + 1]);
      c3 = ae_c_div(c1, c2);
      buf.xR[2 * i] = c3.x;
      buf.xR[2 * i + 1] = -c3.y;
   }
   ftapplyplan(&plan, &buf, 0, 1);
   t = 1.0 / m;
   ae_vector_set_length(r, m);
   for (i = 0; i < m; i++) {
      r->xC[i] = complex_from_d(t * buf.xR[2 * i], -t * buf.xR[2 * i + 1]);
   }
   ae_frame_leave();
}

// 1-dimensional real convolution.
//
// Extended subroutine which allows to choose convolution algorithm.
// Intended for internal use, ALGLIB users should call ConvR1D().
//
// Inputs:
//     A   -   array[0..M-1] - complex function to be transformed
//     M   -   problem size
//     B   -   array[0..N-1] - complex function to be transformed
//     N   -   problem size, N <= M
//     Alg -   algorithm type:
//             *-2     auto-select Q for overlap-add
//             *-1     auto-select algorithm and parameters
//             * 0     straightforward formula for small N's
//             * 1     general FFT-based code
//             * 2     overlap-add with length Q
//     Q   -   length for overlap-add
//
// Outputs:
//     R   -   convolution: A*B. array[0..N+M-1].
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
void convr1dx(RVector *a, ae_int_t m, RVector *b, ae_int_t n, bool circular, ae_int_t alg, ae_int_t q, RVector *r) {
   ae_frame _frame_block;
   double v;
   ae_int_t i;
   ae_int_t j;
   ae_int_t p;
   ae_int_t ptotal;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t j1;
   ae_int_t j2;
   double ax;
   double ay;
   double bx;
   double by;
   double tx;
   double ty;
   double flopcand;
   double flopbest;
   ae_int_t algbest;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewObj(fasttransformplan, plan);
   NewVector(buf, 0, DT_REAL);
   NewVector(buf2, 0, DT_REAL);
   NewVector(buf3, 0, DT_REAL);
   ae_assert(n > 0 && m > 0, "ConvC1DX: incorrect N or M!");
   ae_assert(n <= m, "ConvC1DX: N<M assumption is false!");
// handle special cases
   if (imin2(m, n) <= 2) {
      alg = 0;
   }
// Auto-select
   if (alg < 0) {
   // Initial candidate: straightforward implementation.
   //
   // If we want to use auto-fitted overlap-add,
   // flop count is initialized by large real number - to force
   // another algorithm selection
      algbest = 0;
      if (alg == -1) {
         flopbest = 0.15 * m * n;
      } else {
         flopbest = maxrealnumber;
      }
   // Another candidate - generic FFT code
      if (alg == -1) {
         if (circular && ftbaseissmooth(m) && m % 2 == 0) {
         // special code for circular convolution of a sequence with a smooth length
            flopcand = 3 * ftbasegetflopestimate(m / 2) + (double)(6 * m) / 2.0;
            if (flopcand < flopbest) {
               algbest = 1;
               flopbest = flopcand;
            }
         } else {
         // general cyclic/non-cyclic convolution
            p = ftbasefindsmootheven(m + n - 1);
            flopcand = 3 * ftbasegetflopestimate(p / 2) + (double)(6 * p) / 2.0;
            if (flopcand < flopbest) {
               algbest = 1;
               flopbest = flopcand;
            }
         }
      }
   // Another candidate - overlap-add
      q = 1;
      ptotal = 1;
      while (ptotal < n) {
         ptotal *= 2;
      }
      while (ptotal < m + n) {
         p = ptotal - n + 1;
         flopcand = iceil((double)m / p) * (2 * ftbasegetflopestimate(ptotal / 2) + 1 * (ptotal / 2));
         if (flopcand < flopbest) {
            flopbest = flopcand;
            algbest = 2;
            q = p;
         }
         ptotal *= 2;
      }
      alg = algbest;
      convr1dx(a, m, b, n, circular, alg, q, r);
      ae_frame_leave();
      return;
   }
// straightforward formula for
// circular and non-circular convolutions.
//
// Very simple code, no further comments needed.
   if (alg == 0) {
   // Special case: N=1
      if (n == 1) {
         ae_vector_set_length(r, m);
         v = b->xR[0];
         ae_v_moved(r->xR, 1, a->xR, 1, m, v);
         ae_frame_leave();
         return;
      }
   // use straightforward formula
      if (circular) {
      // circular convolution
         ae_vector_set_length(r, m);
         v = b->xR[0];
         ae_v_moved(r->xR, 1, a->xR, 1, m, v);
         for (i = 1; i < n; i++) {
            v = b->xR[i];
            i1 = 0;
            i2 = i - 1;
            j1 = m - i;
            j2 = m - 1;
            ae_v_addd(&r->xR[i1], 1, &a->xR[j1], 1, i2 - i1 + 1, v);
            i1 = i;
            i2 = m - 1;
            j1 = 0;
            j2 = m - i - 1;
            ae_v_addd(&r->xR[i1], 1, &a->xR[j1], 1, i2 - i1 + 1, v);
         }
      } else {
      // non-circular convolution
         ae_vector_set_length(r, m + n - 1);
         for (i = 0; i < m + n - 1; i++) {
            r->xR[i] = 0.0;
         }
         for (i = 0; i < n; i++) {
            v = b->xR[i];
            ae_v_addd(&r->xR[i], 1, a->xR, 1, m, v);
         }
      }
      ae_frame_leave();
      return;
   }
// general FFT-based code for
// circular and non-circular convolutions.
//
// First, if convolution is circular, we test whether M is smooth or not.
// If it is smooth, we just use M-length FFT to calculate convolution.
// If it is not, we calculate non-circular convolution and wrap it arount.
//
// If convolution is non-circular, we use zero-padding + FFT.
//
// We assume that M+N-1>2 - we should call small case code otherwise
   if (alg == 1) {
      ae_assert(m + n - 1 > 2, "ConvR1DX: internal error!");
      if (circular && ftbaseissmooth(m) && m % 2 == 0) {
      // special code for circular convolution with smooth even M
         ae_vector_set_length(&buf, m);
         ae_v_move(buf.xR, 1, a->xR, 1, m);
         ae_vector_set_length(&buf2, m);
         ae_v_move(buf2.xR, 1, b->xR, 1, n);
         for (i = n; i < m; i++) {
            buf2.xR[i] = 0.0;
         }
         ae_vector_set_length(&buf3, m);
         ftcomplexfftplan(m / 2, 1, &plan);
         fftr1dinternaleven(&buf, m, &buf3, &plan);
         fftr1dinternaleven(&buf2, m, &buf3, &plan);
         buf.xR[0] *= buf2.xR[0];
         buf.xR[1] *= buf2.xR[1];
         for (i = 1; i < m / 2; i++) {
            ax = buf.xR[2 * i];
            ay = buf.xR[2 * i + 1];
            bx = buf2.xR[2 * i];
            by = buf2.xR[2 * i + 1];
            tx = ax * bx - ay * by;
            ty = ax * by + ay * bx;
            buf.xR[2 * i] = tx;
            buf.xR[2 * i + 1] = ty;
         }
         fftr1dinvinternaleven(&buf, m, &buf3, &plan);
         ae_vector_set_length(r, m);
         ae_v_move(r->xR, 1, buf.xR, 1, m);
      } else {
      // M is non-smooth or non-even, general code (circular/non-circular):
      // * first part is the same for circular and non-circular
      //   convolutions. zero padding, FFTs, inverse FFTs
      // * second part differs:
      //   * for non-circular convolution we just copy array
      //   * for circular convolution we add array tail to its head
         p = ftbasefindsmootheven(m + n - 1);
         ae_vector_set_length(&buf, p);
         ae_v_move(buf.xR, 1, a->xR, 1, m);
         for (i = m; i < p; i++) {
            buf.xR[i] = 0.0;
         }
         ae_vector_set_length(&buf2, p);
         ae_v_move(buf2.xR, 1, b->xR, 1, n);
         for (i = n; i < p; i++) {
            buf2.xR[i] = 0.0;
         }
         ae_vector_set_length(&buf3, p);
         ftcomplexfftplan(p / 2, 1, &plan);
         fftr1dinternaleven(&buf, p, &buf3, &plan);
         fftr1dinternaleven(&buf2, p, &buf3, &plan);
         buf.xR[0] *= buf2.xR[0];
         buf.xR[1] *= buf2.xR[1];
         for (i = 1; i < p / 2; i++) {
            ax = buf.xR[2 * i];
            ay = buf.xR[2 * i + 1];
            bx = buf2.xR[2 * i];
            by = buf2.xR[2 * i + 1];
            tx = ax * bx - ay * by;
            ty = ax * by + ay * bx;
            buf.xR[2 * i] = tx;
            buf.xR[2 * i + 1] = ty;
         }
         fftr1dinvinternaleven(&buf, p, &buf3, &plan);
         if (circular) {
         // circular, add tail to head
            ae_vector_set_length(r, m);
            ae_v_move(r->xR, 1, buf.xR, 1, m);
            if (n >= 2) {
               ae_v_add(r->xR, 1, &buf.xR[m], 1, n - 1);
            }
         } else {
         // non-circular, just copy
            ae_vector_set_length(r, m + n - 1);
            ae_v_move(r->xR, 1, buf.xR, 1, m + n - 1);
         }
      }
      ae_frame_leave();
      return;
   }
// overlap-add method
   if (alg == 2) {
      ae_assert((q + n - 1) % 2 == 0, "ConvR1DX: internal error!");
      ae_vector_set_length(&buf, q + n - 1);
      ae_vector_set_length(&buf2, q + n - 1);
      ae_vector_set_length(&buf3, q + n - 1);
      ftcomplexfftplan((q + n - 1) / 2, 1, &plan);
   // prepare R
      if (circular) {
         ae_vector_set_length(r, m);
         for (i = 0; i < m; i++) {
            r->xR[i] = 0.0;
         }
      } else {
         ae_vector_set_length(r, m + n - 1);
         for (i = 0; i < m + n - 1; i++) {
            r->xR[i] = 0.0;
         }
      }
   // pre-calculated FFT(B)
      ae_v_move(buf2.xR, 1, b->xR, 1, n);
      for (j = n; j < q + n - 1; j++) {
         buf2.xR[j] = 0.0;
      }
      fftr1dinternaleven(&buf2, q + n - 1, &buf3, &plan);
   // main overlap-add cycle
      i = 0;
      while (i < m) {
         p = imin2(q, m - i);
         ae_v_move(buf.xR, 1, &a->xR[i], 1, p);
         for (j = p; j < q + n - 1; j++) {
            buf.xR[j] = 0.0;
         }
         fftr1dinternaleven(&buf, q + n - 1, &buf3, &plan);
         buf.xR[0] *= buf2.xR[0];
         buf.xR[1] *= buf2.xR[1];
         for (j = 1; j < (q + n - 1) / 2; j++) {
            ax = buf.xR[2 * j];
            ay = buf.xR[2 * j + 1];
            bx = buf2.xR[2 * j];
            by = buf2.xR[2 * j + 1];
            tx = ax * bx - ay * by;
            ty = ax * by + ay * bx;
            buf.xR[2 * j] = tx;
            buf.xR[2 * j + 1] = ty;
         }
         fftr1dinvinternaleven(&buf, q + n - 1, &buf3, &plan);
         if (circular) {
            j1 = imin2(i + p + n - 2, m - 1) - i;
            j2 = j1 + 1;
         } else {
            j1 = p + n - 2;
            j2 = j1 + 1;
         }
         ae_v_add(&r->xR[i], 1, buf.xR, 1, j1 + 1);
         if (p + n - 2 >= j2) {
            ae_v_add(r->xR, 1, &buf.xR[j2], 1, p + n - 1 - j2);
         }
         i += p;
      }
      ae_frame_leave();
      return;
   }
   ae_frame_leave();
}

// 1-dimensional real convolution.
//
// Analogous to ConvC1D(), see ConvC1D() comments for more details.
//
// Inputs:
//     A   -   array[0..M-1] - real function to be transformed
//     M   -   problem size
//     B   -   array[0..N-1] - real function to be transformed
//     N   -   problem size
//
// Outputs:
//     R   -   convolution: A*B. array[0..N+M-2].
//
// NOTE:
//     It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
// functions have non-zero values at negative T's, you  can  still  use  this
// subroutine - just shift its result correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convr1d(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
void convr1d(RVector *a, ae_int_t m, RVector *b, ae_int_t n, RVector *r) {
   SetVector(r);
   ae_assert(n > 0 && m > 0, "ConvR1D: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer that B.
   if (m < n) {
      convr1d(b, n, a, m, r);
      return;
   }
   convr1dx(a, m, b, n, false, -1, 0, r);
}

// 1-dimensional real deconvolution (inverse of ConvC1D()).
//
// Algorithm has M*log(M)) complexity for any M (composite or prime).
//
// Inputs:
//     A   -   array[0..M-1] - convolved signal, A = conv(R, B)
//     M   -   convolved signal length
//     B   -   array[0..N-1] - response
//     N   -   response length, N <= M
//
// Outputs:
//     R   -   deconvolved signal. array[0..M-N].
//
// NOTE:
//     deconvolution is unstable process and may result in division  by  zero
// (if your response function is degenerate, i.e. has zero Fourier coefficient).
//
// NOTE:
//     It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
// functions have non-zero values at negative T's, you  can  still  use  this
// subroutine - just shift its result correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convr1dinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
void convr1dinv(RVector *a, ae_int_t m, RVector *b, ae_int_t n, RVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t p;
   complex c1;
   complex c2;
   complex c3;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(buf, 0, DT_REAL);
   NewVector(buf2, 0, DT_REAL);
   NewVector(buf3, 0, DT_REAL);
   NewObj(fasttransformplan, plan);
   ae_assert(n > 0 && m > 0 && n <= m, "ConvR1DInv: incorrect N or M!");
   p = ftbasefindsmootheven(m);
   ae_vector_set_length(&buf, p);
   ae_v_move(buf.xR, 1, a->xR, 1, m);
   for (i = m; i < p; i++) {
      buf.xR[i] = 0.0;
   }
   ae_vector_set_length(&buf2, p);
   ae_v_move(buf2.xR, 1, b->xR, 1, n);
   for (i = n; i < p; i++) {
      buf2.xR[i] = 0.0;
   }
   ae_vector_set_length(&buf3, p);
   ftcomplexfftplan(p / 2, 1, &plan);
   fftr1dinternaleven(&buf, p, &buf3, &plan);
   fftr1dinternaleven(&buf2, p, &buf3, &plan);
   buf.xR[0] /= buf2.xR[0];
   buf.xR[1] /= buf2.xR[1];
   for (i = 1; i < p / 2; i++) {
      c1 = complex_from_d(buf.xR[2 * i], buf.xR[2 * i + 1]);
      c2 = complex_from_d(buf2.xR[2 * i], buf2.xR[2 * i + 1]);
      c3 = ae_c_div(c1, c2);
      buf.xR[2 * i] = c3.x;
      buf.xR[2 * i + 1] = c3.y;
   }
   fftr1dinvinternaleven(&buf, p, &buf3, &plan);
   ae_vector_set_length(r, m - n + 1);
   ae_v_move(r->xR, 1, buf.xR, 1, m - n + 1);
   ae_frame_leave();
}

// 1-dimensional circular real convolution.
//
// Analogous to ConvC1DCircular(), see ConvC1DCircular() comments for more details.
//
// Inputs:
//     S   -   array[0..M-1] - real signal
//     M   -   problem size
//     B   -   array[0..N-1] - real response
//     N   -   problem size
//
// Outputs:
//     R   -   convolution: A*B. array[0..M-1].
//
// NOTE:
//     It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
// negative T's, you can still use this subroutine - just  shift  its  result
// correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convr1dcircular(const real_1d_array &s, const ae_int_t m, const real_1d_array &r, const ae_int_t n, real_1d_array &c);
void convr1dcircular(RVector *s, ae_int_t m, RVector *r, ae_int_t n, RVector *c) {
   ae_frame _frame_block;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t j2;
   ae_frame_make(&_frame_block);
   SetVector(c);
   NewVector(buf, 0, DT_REAL);
   ae_assert(n > 0 && m > 0, "ConvC1DCircular: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer (at least - not shorter) that B.
   if (m < n) {
      ae_vector_set_length(&buf, m);
      for (i1 = 0; i1 < m; i1++) {
         buf.xR[i1] = 0.0;
      }
      i1 = 0;
      while (i1 < n) {
         i2 = imin2(i1 + m - 1, n - 1);
         j2 = i2 - i1;
         ae_v_add(buf.xR, 1, &r->xR[i1], 1, j2 + 1);
         i1 += m;
      }
      convr1dcircular(s, m, &buf, m, c);
      ae_frame_leave();
      return;
   }
// reduce to usual convolution
   convr1dx(s, m, r, n, true, -1, 0, c);
   ae_frame_leave();
}

// 1-dimensional complex deconvolution (inverse of ConvC1D()).
//
// Algorithm has M*log(M)) complexity for any M (composite or prime).
//
// Inputs:
//     A   -   array[0..M-1] - convolved signal, A = conv(R, B)
//     M   -   convolved signal length
//     B   -   array[0..N-1] - response
//     N   -   response length
//
// Outputs:
//     R   -   deconvolved signal. array[0..M-N].
//
// NOTE:
//     deconvolution is unstable process and may result in division  by  zero
// (if your response function is degenerate, i.e. has zero Fourier coefficient).
//
// NOTE:
//     It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
// negative T's, you can still use this subroutine - just  shift  its  result
// correspondingly.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void convr1dcircularinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
void convr1dcircularinv(RVector *a, ae_int_t m, RVector *b, ae_int_t n, RVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t j2;
   complex c1;
   complex c2;
   complex c3;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(buf, 0, DT_REAL);
   NewVector(buf2, 0, DT_REAL);
   NewVector(buf3, 0, DT_REAL);
   NewVector(cbuf, 0, DT_COMPLEX);
   NewVector(cbuf2, 0, DT_COMPLEX);
   NewObj(fasttransformplan, plan);
   ae_assert(n > 0 && m > 0, "ConvR1DCircularInv: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer (at least - not shorter) that B.
   if (m < n) {
      ae_vector_set_length(&buf, m);
      for (i = 0; i < m; i++) {
         buf.xR[i] = 0.0;
      }
      i1 = 0;
      while (i1 < n) {
         i2 = imin2(i1 + m - 1, n - 1);
         j2 = i2 - i1;
         ae_v_add(buf.xR, 1, &b->xR[i1], 1, j2 + 1);
         i1 += m;
      }
      convr1dcircularinv(a, m, &buf, m, r);
      ae_frame_leave();
      return;
   }
// Task is normalized
   ae_int_t Mq = m /2, Mr = m % 2;
   if (Mr == 0) {
   // size is even, use fast even-size FFT
      ae_vector_set_length(&buf, m);
      ae_v_move(buf.xR, 1, a->xR, 1, m);
      ae_vector_set_length(&buf2, m);
      ae_v_move(buf2.xR, 1, b->xR, 1, n);
      for (i = n; i < m; i++) {
         buf2.xR[i] = 0.0;
      }
      ae_vector_set_length(&buf3, m);
      ftcomplexfftplan(Mq, 1, &plan);
      fftr1dinternaleven(&buf, m, &buf3, &plan);
      fftr1dinternaleven(&buf2, m, &buf3, &plan);
      buf.xR[0] /= buf2.xR[0];
      buf.xR[1] /= buf2.xR[1];
      for (i = 1; i < Mq; i++) {
         c1 = complex_from_d(buf.xR[2 * i], buf.xR[2 * i + 1]);
         c2 = complex_from_d(buf2.xR[2 * i], buf2.xR[2 * i + 1]);
         c3 = ae_c_div(c1, c2);
         buf.xR[2 * i] = c3.x;
         buf.xR[2 * i + 1] = c3.y;
      }
      fftr1dinvinternaleven(&buf, m, &buf3, &plan);
      ae_vector_set_length(r, m);
      ae_v_move(r->xR, 1, buf.xR, 1, m);
   } else {
   // odd-size, use general real FFT
      fftr1d(a, m, &cbuf);
      ae_vector_set_length(&buf2, m);
      ae_v_move(buf2.xR, 1, b->xR, 1, n);
      for (i = n; i < m; i++) {
         buf2.xR[i] = 0.0;
      }
      fftr1d(&buf2, m, &cbuf2);
      for (i = 0; i <= Mq; i++) {
         cbuf.xC[i] = ae_c_div(cbuf.xC[i], cbuf2.xC[i]);
      }
      fftr1dinv(&cbuf, m, r);
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void convc1d(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convc1d(ConstT(ae_vector, a), m, ConstT(ae_vector, b), n, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void convc1dinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convc1dinv(ConstT(ae_vector, a), m, ConstT(ae_vector, b), n, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void convc1dcircular(const complex_1d_array &s, const ae_int_t m, const complex_1d_array &r, const ae_int_t n, complex_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convc1dcircular(ConstT(ae_vector, s), m, ConstT(ae_vector, r), n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}

void convc1dcircularinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convc1dcircularinv(ConstT(ae_vector, a), m, ConstT(ae_vector, b), n, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void convr1d(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convr1d(ConstT(ae_vector, a), m, ConstT(ae_vector, b), n, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void convr1dinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convr1dinv(ConstT(ae_vector, a), m, ConstT(ae_vector, b), n, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void convr1dcircular(const real_1d_array &s, const ae_int_t m, const real_1d_array &r, const ae_int_t n, real_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convr1dcircular(ConstT(ae_vector, s), m, ConstT(ae_vector, r), n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}

void convr1dcircularinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::convr1dcircularinv(ConstT(ae_vector, a), m, ConstT(ae_vector, b), n, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === CORR Package ===
// Depends on: CONV
namespace alglib_impl {
// 1-dimensional complex cross-correlation.
//
// For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).
//
// Correlation is calculated using reduction to  convolution.  Algorithm with
// max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
// about performance).
//
// IMPORTANT:
//     for  historical reasons subroutine accepts its parameters in  reversed
//     order: CorrC1D(Signal, Pattern) = Pattern x Signal (using  traditional
//     definition of cross-correlation, denoting cross-correlation as "x").
//
// Inputs:
//     Signal  -   array[0..N-1] - complex function to be transformed,
//                 signal containing pattern
//     N       -   problem size
//     Pattern -   array[0..M-1] - complex function to be transformed,
//                 pattern to search withing signal
//     M       -   problem size
//
// Outputs:
//     R       -   cross-correlation, array[0..N+M-2]:
//                 * positive lags are stored in R[0..N-1],
//                   R[i] = sum(conj(pattern[j])*signal[i+j]
//                 * negative lags are stored in R[N..N+M-2],
//                   R[N+M-1-i] = sum(conj(pattern[j])*signal[-i+j]
//
// NOTE:
//     It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
// on [-K..M-1],  you can still use this subroutine, just shift result by K.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void corrc1d(const complex_1d_array &signal, const ae_int_t n, const complex_1d_array &pattern, const ae_int_t m, complex_1d_array &r);
void corrc1d(CVector *signal, ae_int_t n, CVector *pattern, ae_int_t m, CVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(p, 0, DT_COMPLEX);
   NewVector(b, 0, DT_COMPLEX);
   ae_assert(n > 0 && m > 0, "CorrC1D: incorrect N or M!");
   ae_vector_set_length(&p, m);
   for (i = 0; i < m; i++) {
      p.xC[m - 1 - i] = conj(pattern->xC[i]);
   }
   convc1d(&p, m, signal, n, &b);
   ae_vector_set_length(r, m + n - 1);
   ae_v_cmove(r->xC, 1, &b.xC[m - 1], 1, "N", n);
   if (m + n - 2 >= n) {
      ae_v_cmove(&r->xC[n], 1, b.xC, 1, "N", m - 1);
   }
   ae_frame_leave();
}

// 1-dimensional circular complex cross-correlation.
//
// For given Pattern/Signal returns corr(Pattern,Signal) (circular).
// Algorithm has linearithmic complexity for any M/N.
//
// IMPORTANT:
//     for  historical reasons subroutine accepts its parameters in  reversed
//     order:   CorrC1DCircular(Signal, Pattern) = Pattern x Signal    (using
//     traditional definition of cross-correlation, denoting cross-correlation
//     as "x").
//
// Inputs:
//     Signal  -   array[0..N-1] - complex function to be transformed,
//                 periodic signal containing pattern
//     N       -   problem size
//     Pattern -   array[0..M-1] - complex function to be transformed,
//                 non-periodic pattern to search withing signal
//     M       -   problem size
//
// Outputs:
//     R   -   convolution: A*B. array[0..M-1].
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void corrc1dcircular(const complex_1d_array &signal, const ae_int_t m, const complex_1d_array &pattern, const ae_int_t n, complex_1d_array &c);
void corrc1dcircular(CVector *signal, ae_int_t m, CVector *pattern, ae_int_t n, CVector *c) {
   ae_frame _frame_block;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t i;
   ae_int_t j2;
   ae_frame_make(&_frame_block);
   SetVector(c);
   NewVector(p, 0, DT_COMPLEX);
   NewVector(b, 0, DT_COMPLEX);
   ae_assert(n > 0 && m > 0, "ConvC1DCircular: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer (at least - not shorter) that B.
   if (m < n) {
      ae_vector_set_length(&b, m);
      for (i1 = 0; i1 < m; i1++) {
         b.xC[i1] = complex_from_i(0);
      }
      i1 = 0;
      while (i1 < n) {
         i2 = imin2(i1 + m - 1, n - 1);
         j2 = i2 - i1;
         ae_v_cadd(b.xC, 1, &pattern->xC[i1], 1, "N", j2 + 1);
         i1 += m;
      }
      corrc1dcircular(signal, m, &b, m, c);
      ae_frame_leave();
      return;
   }
// Task is normalized
   ae_vector_set_length(&p, n);
   for (i = 0; i < n; i++) {
      p.xC[n - 1 - i] = conj(pattern->xC[i]);
   }
   convc1dcircular(signal, m, &p, n, &b);
   ae_vector_set_length(c, m);
   ae_v_cmove(c->xC, 1, &b.xC[n - 1], 1, "N", m - n + 1);
   if (1 < n) {
      ae_v_cmove(&c->xC[m - n + 1], 1, b.xC, 1, "N", n - 1);
   }
   ae_frame_leave();
}

// 1-dimensional real cross-correlation.
//
// For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).
//
// Correlation is calculated using reduction to  convolution.  Algorithm with
// max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
// about performance).
//
// IMPORTANT:
//     for  historical reasons subroutine accepts its parameters in  reversed
//     order: CorrR1D(Signal, Pattern) = Pattern x Signal (using  traditional
//     definition of cross-correlation, denoting cross-correlation as "x").
//
// Inputs:
//     Signal  -   array[0..N-1] - real function to be transformed,
//                 signal containing pattern
//     N       -   problem size
//     Pattern -   array[0..M-1] - real function to be transformed,
//                 pattern to search withing signal
//     M       -   problem size
//
// Outputs:
//     R       -   cross-correlation, array[0..N+M-2]:
//                 * positive lags are stored in R[0..N-1],
//                   R[i] = sum(pattern[j]*signal[i+j]
//                 * negative lags are stored in R[N..N+M-2],
//                   R[N+M-1-i] = sum(pattern[j]*signal[-i+j]
//
// NOTE:
//     It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
// on [-K..M-1],  you can still use this subroutine, just shift result by K.
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void corrr1d(const real_1d_array &signal, const ae_int_t n, const real_1d_array &pattern, const ae_int_t m, real_1d_array &r);
void corrr1d(RVector *signal, ae_int_t n, RVector *pattern, ae_int_t m, RVector *r) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   SetVector(r);
   NewVector(p, 0, DT_REAL);
   NewVector(b, 0, DT_REAL);
   ae_assert(n > 0 && m > 0, "CorrR1D: incorrect N or M!");
   ae_vector_set_length(&p, m);
   for (i = 0; i < m; i++) {
      p.xR[m - 1 - i] = pattern->xR[i];
   }
   convr1d(&p, m, signal, n, &b);
   ae_vector_set_length(r, m + n - 1);
   ae_v_move(r->xR, 1, &b.xR[m - 1], 1, n);
   if (m + n - 2 >= n) {
      ae_v_move(&r->xR[n], 1, b.xR, 1, m - 1);
   }
   ae_frame_leave();
}

// 1-dimensional circular real cross-correlation.
//
// For given Pattern/Signal returns corr(Pattern,Signal) (circular).
// Algorithm has linearithmic complexity for any M/N.
//
// IMPORTANT:
//     for  historical reasons subroutine accepts its parameters in  reversed
//     order:   CorrR1DCircular(Signal, Pattern) = Pattern x Signal    (using
//     traditional definition of cross-correlation, denoting cross-correlation
//     as "x").
//
// Inputs:
//     Signal  -   array[0..N-1] - real function to be transformed,
//                 periodic signal containing pattern
//     N       -   problem size
//     Pattern -   array[0..M-1] - real function to be transformed,
//                 non-periodic pattern to search withing signal
//     M       -   problem size
//
// Outputs:
//     R   -   convolution: A*B. array[0..M-1].
// ALGLIB: Copyright 21.07.2009 by Sergey Bochkanov
// API: void corrr1dcircular(const real_1d_array &signal, const ae_int_t m, const real_1d_array &pattern, const ae_int_t n, real_1d_array &c);
void corrr1dcircular(RVector *signal, ae_int_t m, RVector *pattern, ae_int_t n, RVector *c) {
   ae_frame _frame_block;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t i;
   ae_int_t j2;
   ae_frame_make(&_frame_block);
   SetVector(c);
   NewVector(p, 0, DT_REAL);
   NewVector(b, 0, DT_REAL);
   ae_assert(n > 0 && m > 0, "ConvC1DCircular: incorrect N or M!");
// normalize task: make M >= N,
// so A will be longer (at least - not shorter) that B.
   if (m < n) {
      ae_vector_set_length(&b, m);
      for (i1 = 0; i1 < m; i1++) {
         b.xR[i1] = 0.0;
      }
      i1 = 0;
      while (i1 < n) {
         i2 = imin2(i1 + m - 1, n - 1);
         j2 = i2 - i1;
         ae_v_add(b.xR, 1, &pattern->xR[i1], 1, j2 + 1);
         i1 += m;
      }
      corrr1dcircular(signal, m, &b, m, c);
      ae_frame_leave();
      return;
   }
// Task is normalized
   ae_vector_set_length(&p, n);
   for (i = 0; i < n; i++) {
      p.xR[n - 1 - i] = pattern->xR[i];
   }
   convr1dcircular(signal, m, &p, n, &b);
   ae_vector_set_length(c, m);
   ae_v_move(c->xR, 1, &b.xR[n - 1], 1, m - n + 1);
   if (1 < n) {
      ae_v_move(&c->xR[m - n + 1], 1, b.xR, 1, n - 1);
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void corrc1d(const complex_1d_array &signal, const ae_int_t n, const complex_1d_array &pattern, const ae_int_t m, complex_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::corrc1d(ConstT(ae_vector, signal), n, ConstT(ae_vector, pattern), m, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void corrc1dcircular(const complex_1d_array &signal, const ae_int_t m, const complex_1d_array &pattern, const ae_int_t n, complex_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::corrc1dcircular(ConstT(ae_vector, signal), m, ConstT(ae_vector, pattern), n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}

void corrr1d(const real_1d_array &signal, const ae_int_t n, const real_1d_array &pattern, const ae_int_t m, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::corrr1d(ConstT(ae_vector, signal), n, ConstT(ae_vector, pattern), m, ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void corrr1dcircular(const real_1d_array &signal, const ae_int_t m, const real_1d_array &pattern, const ae_int_t n, real_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::corrr1dcircular(ConstT(ae_vector, signal), m, ConstT(ae_vector, pattern), n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib
