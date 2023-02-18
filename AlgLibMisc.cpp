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
#include "AlgLibMisc.h"

// === HQRND Package ===
// Depends on: (AlgLibInternal) ABLASF
namespace alglib_impl {
static const ae_int_t hqrnd_hqrndm1 = 1 + 2 * 3 * 7 * 631 * 81031; // == 2^31 - 85
static const ae_int_t hqrnd_hqrndm2 = 1 + 2 * 19 * 31 * 1019 * 1789; // == 2^31 - 249
static const ae_int_t hqrnd_hqrndmax = hqrnd_hqrndm1 - 2;
static const ae_int_t hqrnd_hqrndmagic = 1634357784;

// HQRNDState initialization with seed values
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: void hqrndseed(const ae_int_t s1, const ae_int_t s2, hqrndstate &state);
void hqrndseed(ae_int_t s1, ae_int_t s2, hqrndstate *state) {
   SetObj(hqrndstate, state);
// Protection against negative seeds:
//
//     SEED = -(SEED+1)
//
// We can't use just "-SEED" because there exists an integer number  N
// such that -N < 0 and +N < 0. (This number is equal to 0x800...000).
// The need  to handle such  a seed correctly  forces us to  use a bit
// more complicated formula.
   if (s1 < 0) {
      s1 = -(s1 + 1);
   }
   if (s2 < 0) {
      s2 = -(s2 + 1);
   }
   state->s1 = s1 % (hqrnd_hqrndm1 - 1) + 1;
   state->s2 = s2 % (hqrnd_hqrndm2 - 1) + 1;
   state->magicv = hqrnd_hqrndmagic;
}

// HQRNDState  initialization  with  random  values  which come from standard
// RNG.
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: void hqrndrandomize(hqrndstate &state);
void hqrndrandomize(hqrndstate *state) {
   ae_int_t s0;
   ae_int_t s1;
   SetObj(hqrndstate, state);
   s0 = randominteger(hqrnd_hqrndm1);
   s1 = randominteger(hqrnd_hqrndm2);
   hqrndseed(s0, s1, state);
}

// This function returns random integer in [0,HQRNDMax]
//
// P.L'Ecuyer, ‟Efficient and portable combined random number generators”, CACM 31(6), June 1988, 742-751.
static ae_int_t hqrnd_hqrndintegerbase(hqrndstate *state) {
   const ae_int_t a1 = 40014, q1 = hqrnd_hqrndm1 / a1, r1 = hqrnd_hqrndm1 % a1;
   const ae_int_t a2 = 40692, q2 = hqrnd_hqrndm2 / a2, r2 = hqrnd_hqrndm2 % a2;
// This requires:
// ∙    (hqrnd_hqrndm1 - 1)/2 and (hqrnd_hqrndm2 - 1)/2 be relatively prime,
// ∙    a1^2 < hqrnd_hqrndm1 and a2^2 < hqrnd_hqrndm2,
// ∙    0 <= state->s1 < hqrnd_hqrndm1, 0 <= state->s2 < hqrnd_hqrndm2.
// Under these conditions, the following is equivalent to:
// ∙    state->s1 = a1 * state->s1 % hqrnd_hqrndm1, state->s2 = a2 * state->s2 % hqrnd_hqrndm2;
// and 0 <= state->s1 < hqrnd_hqrndm1, 0 <= state->s2 < hqrnd_hqrndm2 continues to hold.
// The period of the generator is:
//      (hqrnd_hqrndm1 - 1)(hqrnd_hqrndm2 - 1)/2 == 2 * 3 * 7 * 19 * 31 * 631 * 1019 * 1789 * 81031 == 2^61 - 360777242114
   ae_int_t k;
   ae_int_t result;
   ae_assert(state->magicv == hqrnd_hqrndmagic, "hqrnd_hqrndintegerbase: state is not correctly initialized!");
   k = state->s1 / q1;
   state->s1 = a1 * (state->s1 - k * q1) - k * r1;
   if (state->s1 < 0) {
      state->s1 += hqrnd_hqrndm1;
   }
   k = state->s2 / q2;
   state->s2 = a2 * (state->s2 - k * q2) - k * r2;
   if (state->s2 < 0) {
      state->s2 += hqrnd_hqrndm2;
   }
// Result
   result = state->s1 - state->s2;
   if (result < 1) {
      result += hqrnd_hqrndmax;
   } else {
      result--;
   }
   return result;
}

// This function generates random integer number in [0, N)
//
// 1. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()
// 2. N can be any positive number except for very large numbers:
//    * close to 2^31 on 32-bit systems
//    * close to 2^62 on 64-bit systems
//    An exception will be generated if N is too large.
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: ae_int_t hqrnduniformi(const hqrndstate &state, const ae_int_t n);
ae_int_t hqrnduniformi(hqrndstate *state, ae_int_t n) {
   ae_int_t maxcnt;
   ae_int_t mx;
   ae_int_t a;
   ae_int_t b;
   ae_int_t result;
   ae_assert(n > 0, "hqrnduniformi: N <= 0!");
   maxcnt = hqrnd_hqrndmax + 1;
   if (n > maxcnt) {
   // We have two options here:
   // a) N is exactly divisible by MaxCnt
   // b) N is not divisible by MaxCnt
   // In both cases we reduce problem on interval spanning [0,N)
   // to several subproblems on intervals spanning [0,MaxCnt).
   // [0,N) range is divided into ceil(N/MaxCnt) bins,
   // each of them having length equal to MaxCnt.
      if (n % maxcnt == 0) {
      // We generate:
      // * random bin number B
      // * random offset within bin A
      // Both random numbers are generated by recursively
      // calling hqrnduniformi().
      //
      // Result is equal to A+MaxCnt*B.
         ae_assert(n / maxcnt <= maxcnt, "hqrnduniformi: N is too large");
         a = hqrnduniformi(state, maxcnt);
         b = hqrnduniformi(state, n / maxcnt);
         result = a + maxcnt * b;
      } else {
      // We generate:
      // * random bin number B in [0, ceil(N/MaxCnt)-1]
      // * random offset within bin A
      // * if both of what is below is true
      //   1) bin number B is that of the last bin
      //   2) A >= N mod MaxCnt
      //   then we repeat generation of A/B.
      //   This stage is essential in order to avoid bias in the result.
      // * otherwise, we return A*MaxCnt+N
         ae_assert(n / maxcnt + 1 <= maxcnt, "hqrnduniformi: N is too large");
         result = -1;
         do {
            a = hqrnduniformi(state, maxcnt);
            b = hqrnduniformi(state, n / maxcnt + 1);
            if (b == n / maxcnt && a >= n % maxcnt) {
               continue;
            }
            result = a + maxcnt * b;
         } while (result < 0);
      }
   } else {
   // Code below is a bit complicated because we can not simply
   // return "hqrnd_hqrndintegerbase() mod N" - it will be skewed for
   // large N's in [0.1*hqrnd_hqrndmax...hqrnd_hqrndmax].
      mx = maxcnt - maxcnt % n;
      do {
         result = hqrnd_hqrndintegerbase(state);
      } while (result >= mx);
      result %= n;
   }
   return result;
}

// This function generates random real number in (0,1),
// not including interval boundaries
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: double hqrnduniformr(const hqrndstate &state);
double hqrnduniformr(hqrndstate *state) {
   double result;
   result = (hqrnd_hqrndintegerbase(state) + 1.0) / (hqrnd_hqrndmax + 2);
   return result;
}

// This function generates random real number in (-1,+1),
// not including interval boundaries
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: double hqrndmiduniformr(const hqrndstate &state);
double hqrndmiduniformr(hqrndstate *state) {
   double result;
   result = (2.0 * hqrnd_hqrndintegerbase(state) - hqrnd_hqrndmax) / (hqrnd_hqrndmax + 2);
   return result;
}

// Random number generator: normal numbers
//
// This function generates two independent random numbers from normal
// distribution. Its performance is equal to that of HQRNDNormal()
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: void hqrndnormal2(const hqrndstate &state, double &x1, double &x2);
void hqrndnormal2(hqrndstate *state, double *x1, double *x2) {
   double u;
   double v;
   double s;
   *x1 = 0;
   *x2 = 0;
   while (true) {
      u = hqrndmiduniformr(state);
      v = hqrndmiduniformr(state);
      s = sqr(u) + sqr(v);
      if (s > 0.0 && s < 1.0) {
      // two sqrt's instead of one to
      // avoid overflow when S is too small
         s = sqrt(-2 * log(s)) / sqrt(s);
         *x1 = u * s;
         *x2 = v * s;
         return;
      }
   }
}

// Random number generator: normal numbers
//
// This function generates one random number from normal distribution.
// Its performance is equal to that of HQRNDNormal2()
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: double hqrndnormal(const hqrndstate &state);
double hqrndnormal(hqrndstate *state) {
   double v1;
   double v2;
   double result;
   hqrndnormal2(state, &v1, &v2);
   result = v1;
   return result;
}

// Random number generator: vector with random entries (normal distribution)
//
// This function generates N random numbers from normal distribution.
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: void hqrndnormalv(const hqrndstate &state, const ae_int_t n, real_1d_array &x);
void hqrndnormalv(hqrndstate *state, ae_int_t n, RVector *x) {
   ae_int_t i;
   ae_int_t n2;
   double v1;
   double v2;
   SetVector(x);
   n2 = n / 2;
   allocv(n, x);
   for (i = 0; i < n2; i++) {
      hqrndnormal2(state, &v1, &v2);
      x->xR[2 * i] = v1;
      x->xR[2 * i + 1] = v2;
   }
   if (n % 2 != 0) {
      hqrndnormal2(state, &v1, &v2);
      x->xR[n - 1] = v1;
   }
}

// Random number generator: matrix with random entries (normal distribution)
//
// This function generates MxN random matrix.
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: void hqrndnormalm(const hqrndstate &state, const ae_int_t m, const ae_int_t n, real_2d_array &x);
void hqrndnormalm(hqrndstate *state, ae_int_t m, ae_int_t n, RMatrix *x) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t n2;
   double v1;
   double v2;
   SetMatrix(x);
   n2 = n / 2;
   ae_matrix_set_length(x, m, n);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n2; j++) {
         hqrndnormal2(state, &v1, &v2);
         x->xyR[i][2 * j] = v1;
         x->xyR[i][2 * j + 1] = v2;
      }
      if (n % 2 != 0) {
         hqrndnormal2(state, &v1, &v2);
         x->xyR[i][n - 1] = v1;
      }
   }
}

// Random number generator: random X and Y such that X^2+Y^2 == 1
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
// API: void hqrndunit2(const hqrndstate &state, double &x, double &y);
void hqrndunit2(hqrndstate *state, double *x, double *y) {
   double v;
   double mx;
   double mn;
   *x = 0;
   *y = 0;
   do {
      hqrndnormal2(state, x, y);
   } while (!(*x != 0.0 || *y != 0.0));
   mx = rmax2(fabs(*x), fabs(*y));
   mn = rmin2(fabs(*x), fabs(*y));
   v = mx * sqrt(1 + sqr(mn / mx));
   *x /= v;
   *y /= v;
}

// Random number generator: exponential distribution
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 11.08.2007 by Sergey Bochkanov
// API: double hqrndexponential(const hqrndstate &state, const double lambdav);
double hqrndexponential(hqrndstate *state, double lambdav) {
   double result;
   ae_assert(lambdav > 0.0, "hqrndexponential: LambdaV <= 0!");
   result = -log(hqrnduniformr(state)) / lambdav;
   return result;
}

// This function generates  random number from discrete distribution given by
// finite sample X.
//
// Inputs:
//     State   -   high quality random number generator, must be
//                 initialized with HQRNDRandomize() or HQRNDSeed().
//         X   -   finite sample
//         N   -   number of elements to use, N >= 1
//
// Result:
//     this function returns one of the X[i] for random i = 0..N-1
// ALGLIB: Copyright 08.11.2011 by Sergey Bochkanov
// API: double hqrnddiscrete(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
double hqrnddiscrete(hqrndstate *state, RVector *x, ae_int_t n) {
   double result;
   ae_assert(n > 0, "hqrnddiscrete: N <= 0");
   ae_assert(n <= x->cnt, "hqrnddiscrete: Length(X) < N");
   result = x->xR[hqrnduniformi(state, n)];
   return result;
}

// This function generates random number from continuous  distribution  given
// by finite sample X.
//
// Inputs:
//     State   -   high quality random number generator, must be
//                 initialized with HQRNDRandomize() or HQRNDSeed().
//         X   -   finite sample, array[N] (can be larger, in this  case only
//                 leading N elements are used). THIS ARRAY MUST BE SORTED BY
//                 ASCENDING.
//         N   -   number of elements to use, N >= 1
//
// Result:
//     this function returns random number from continuous distribution which
//     tries to approximate X as mush as possible. min(X) <= Result <= max(X).
// ALGLIB: Copyright 08.11.2011 by Sergey Bochkanov
// API: double hqrndcontinuous(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
double hqrndcontinuous(hqrndstate *state, RVector *x, ae_int_t n) {
   double mx;
   double mn;
   ae_int_t i;
   double result;
   ae_assert(n > 0, "hqrndcontinuous: N <= 0");
   ae_assert(n <= x->cnt, "hqrndcontinuous: Length(X) < N");
   if (n == 1) {
      result = x->xR[0];
      return result;
   }
   i = hqrnduniformi(state, n - 1);
   mn = x->xR[i];
   mx = x->xR[i + 1];
   ae_assert(mx >= mn, "hqrndcontinuous: X is not sorted by ascending");
   if (mx != mn) {
      result = (mx - mn) * hqrnduniformr(state) + mn;
   } else {
      result = mn;
   }
   return result;
}

void hqrndstate_init(void *_p, bool make_automatic) {
}

void hqrndstate_copy(void *_dst, void *_src, bool make_automatic) {
   hqrndstate *dst = (hqrndstate *)_dst;
   hqrndstate *src = (hqrndstate *)_src;
   dst->s1 = src->s1;
   dst->s2 = src->s2;
   dst->magicv = src->magicv;
}

void hqrndstate_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
// Portable high quality random number generator state.
// Initialized with HQRNDRandomize() or HQRNDSeed().
//
// Fields:
//     S1, S2      -   seed values
//     MagicV      -   'magic' value used to determine whether State structure
//                     was correctly initialized.
DefClass(hqrndstate, )

void hqrndseed(const ae_int_t s1, const ae_int_t s2, hqrndstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hqrndseed(s1, s2, ConstT(hqrndstate, state));
   alglib_impl::ae_state_clear();
}

void hqrndrandomize(hqrndstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hqrndrandomize(ConstT(hqrndstate, state));
   alglib_impl::ae_state_clear();
}

ae_int_t hqrnduniformi(const hqrndstate &state, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::hqrnduniformi(ConstT(hqrndstate, state), n);
   alglib_impl::ae_state_clear();
   return Z;
}

double hqrnduniformr(const hqrndstate &state) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hqrnduniformr(ConstT(hqrndstate, state));
   alglib_impl::ae_state_clear();
   return D;
}

double hqrndmiduniformr(const hqrndstate &state) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hqrndmiduniformr(ConstT(hqrndstate, state));
   alglib_impl::ae_state_clear();
   return D;
}

void hqrndnormal2(const hqrndstate &state, double &x1, double &x2) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hqrndnormal2(ConstT(hqrndstate, state), &x1, &x2);
   alglib_impl::ae_state_clear();
}

double hqrndnormal(const hqrndstate &state) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hqrndnormal(ConstT(hqrndstate, state));
   alglib_impl::ae_state_clear();
   return D;
}

void hqrndnormalv(const hqrndstate &state, const ae_int_t n, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hqrndnormalv(ConstT(hqrndstate, state), n, ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void hqrndnormalm(const hqrndstate &state, const ae_int_t m, const ae_int_t n, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hqrndnormalm(ConstT(hqrndstate, state), m, n, ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void hqrndunit2(const hqrndstate &state, double &x, double &y) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hqrndunit2(ConstT(hqrndstate, state), &x, &y);
   alglib_impl::ae_state_clear();
}

double hqrndexponential(const hqrndstate &state, const double lambdav) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hqrndexponential(ConstT(hqrndstate, state), lambdav);
   alglib_impl::ae_state_clear();
   return D;
}

double hqrnddiscrete(const hqrndstate &state, const real_1d_array &x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hqrnddiscrete(ConstT(hqrndstate, state), ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}

double hqrndcontinuous(const hqrndstate &state, const real_1d_array &x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hqrndcontinuous(ConstT(hqrndstate, state), ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === XDEBUG Package ===
namespace alglib_impl {
// Debug functions to test the ALGLIB interface generator: not meant for use in production code.

// Creates and returns XDebugRecord1 structure:
// * integer and complex fields of Rec1 are set to 1 and 1+i correspondingly
// * array field of Rec1 is set to [2,3]
// ALGLIB: Copyright 27.05.2014 by Sergey Bochkanov
// API: void xdebuginitrecord1(xdebugrecord1 &rec1);
void xdebuginitrecord1(xdebugrecord1 *rec1) {
   SetObj(xdebugrecord1, rec1);
   rec1->i = 1;
   rec1->c = complex_from_d(1.0, 1.0);
   ae_vector_set_length(&rec1->a, 2);
   rec1->a.xR[0] = 2.0;
   rec1->a.xR[1] = 3.0;
}

// Counts number of True values in the boolean 1D array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: ae_int_t xdebugb1count(const boolean_1d_array &a);
ae_int_t xdebugb1count(BVector *a) {
   ae_int_t i;
   ae_int_t result;
   result = 0;
   for (i = 0; i < a->cnt; i++) {
      if (a->xB[i]) {
         result++;
      }
   }
   return result;
}

// Replace all values in array by NOT(a[i]).
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugb1not(const boolean_1d_array &a);
void xdebugb1not(BVector *a) {
   ae_int_t i;
   for (i = 0; i < a->cnt; i++) {
      a->xB[i] = !a->xB[i];
   }
}

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugb1appendcopy(boolean_1d_array &a);
void xdebugb1appendcopy(BVector *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(b, 0, DT_BOOL);
   ae_vector_set_length(&b, a->cnt);
   for (i = 0; i < b.cnt; i++) {
      b.xB[i] = a->xB[i];
   }
   ae_vector_set_length(a, 2 * b.cnt);
   for (i = 0; i < a->cnt; i++) {
      a->xB[i] = b.xB[i % b.cnt];
   }
   ae_frame_leave();
}

// Generate N-element array with even-numbered elements set to True.
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugb1outeven(const ae_int_t n, boolean_1d_array &a);
void xdebugb1outeven(ae_int_t n, BVector *a) {
   ae_int_t i;
   SetVector(a);
   ae_vector_set_length(a, n);
   for (i = 0; i < a->cnt; i++) {
      a->xB[i] = i % 2 == 0;
   }
}

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: ae_int_t xdebugi1sum(const integer_1d_array &a);
ae_int_t xdebugi1sum(ZVector *a) {
   ae_int_t i;
   ae_int_t result;
   result = 0;
   for (i = 0; i < a->cnt; i++) {
      result += a->xZ[i];
   }
   return result;
}

// Replace all values in array by -A[I]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugi1neg(const integer_1d_array &a);
void xdebugi1neg(ZVector *a) {
   ae_int_t i;
   for (i = 0; i < a->cnt; i++) {
      a->xZ[i] = -a->xZ[i];
   }
}

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugi1appendcopy(integer_1d_array &a);
void xdebugi1appendcopy(ZVector *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(b, 0, DT_INT);
   ae_vector_set_length(&b, a->cnt);
   for (i = 0; i < b.cnt; i++) {
      b.xZ[i] = a->xZ[i];
   }
   ae_vector_set_length(a, 2 * b.cnt);
   for (i = 0; i < a->cnt; i++) {
      a->xZ[i] = b.xZ[i % b.cnt];
   }
   ae_frame_leave();
}

// Generate N-element array with even-numbered A[I] set to I, and odd-numbered
// ones set to 0.
//
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugi1outeven(const ae_int_t n, integer_1d_array &a);
void xdebugi1outeven(ae_int_t n, ZVector *a) {
   ae_int_t i;
   SetVector(a);
   ae_vector_set_length(a, n);
   for (i = 0; i < a->cnt; i++) {
      if (i % 2 == 0) {
         a->xZ[i] = i;
      } else {
         a->xZ[i] = 0;
      }
   }
}

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: double xdebugr1sum(const real_1d_array &a);
double xdebugr1sum(RVector *a) {
   ae_int_t i;
   double result;
   result = 0.0;
   for (i = 0; i < a->cnt; i++) {
      result += a->xR[i];
   }
   return result;
}

// Replace all values in array by -A[I]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugr1neg(const real_1d_array &a);
void xdebugr1neg(RVector *a) {
   ae_int_t i;
   for (i = 0; i < a->cnt; i++) {
      a->xR[i] = -a->xR[i];
   }
}

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugr1appendcopy(real_1d_array &a);
void xdebugr1appendcopy(RVector *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(b, 0, DT_REAL);
   ae_vector_set_length(&b, a->cnt);
   for (i = 0; i < b.cnt; i++) {
      b.xR[i] = a->xR[i];
   }
   ae_vector_set_length(a, 2 * b.cnt);
   for (i = 0; i < a->cnt; i++) {
      a->xR[i] = b.xR[i % b.cnt];
   }
   ae_frame_leave();
}

// Generate N-element array with even-numbered A[I] set to I*0.25,
// and odd-numbered ones are set to 0.
//
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugr1outeven(const ae_int_t n, real_1d_array &a);
void xdebugr1outeven(ae_int_t n, RVector *a) {
   ae_int_t i;
   SetVector(a);
   ae_vector_set_length(a, n);
   for (i = 0; i < a->cnt; i++) {
      if (i % 2 == 0) {
         a->xR[i] = i * 0.25;
      } else {
         a->xR[i] = 0.0;
      }
   }
}

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: complex xdebugc1sum(const complex_1d_array &a);
complex xdebugc1sum(CVector *a) {
   ae_int_t i;
   complex result;
   result = complex_from_i(0);
   for (i = 0; i < a->cnt; i++) {
      result = ae_c_add(result, a->xC[i]);
   }
   return result;
}

// Replace all values in array by -A[I]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugc1neg(const complex_1d_array &a);
void xdebugc1neg(CVector *a) {
   ae_int_t i;
   for (i = 0; i < a->cnt; i++) {
      a->xC[i] = ae_c_neg(a->xC[i]);
   }
}

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugc1appendcopy(complex_1d_array &a);
void xdebugc1appendcopy(CVector *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(b, 0, DT_COMPLEX);
   ae_vector_set_length(&b, a->cnt);
   for (i = 0; i < b.cnt; i++) {
      b.xC[i] = a->xC[i];
   }
   ae_vector_set_length(a, 2 * b.cnt);
   for (i = 0; i < a->cnt; i++) {
      a->xC[i] = b.xC[i % b.cnt];
   }
   ae_frame_leave();
}

// Generate N-element array with even-numbered A[K] set to (x,y) = (K*0.25, K*0.125)
// and odd-numbered ones are set to 0.
//
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugc1outeven(const ae_int_t n, complex_1d_array &a);
void xdebugc1outeven(ae_int_t n, CVector *a) {
   ae_int_t i;
   SetVector(a);
   ae_vector_set_length(a, n);
   for (i = 0; i < a->cnt; i++) {
      if (i % 2 == 0) {
         a->xC[i] = complex_from_d(i * 0.250, i * 0.125);
      } else {
         a->xC[i] = complex_from_i(0);
      }
   }
}

// Counts number of True values in the boolean 2D array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: ae_int_t xdebugb2count(const boolean_2d_array &a);
ae_int_t xdebugb2count(BMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t result;
   result = 0;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         if (a->xyB[i][j]) {
            result++;
         }
      }
   }
   return result;
}

// Replace all values in array by NOT(a[i]).
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugb2not(const boolean_2d_array &a);
void xdebugb2not(BMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyB[i][j] = !a->xyB[i][j];
      }
   }
}

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugb2transpose(boolean_2d_array &a);
void xdebugb2transpose(BMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   NewMatrix(b, 0, 0, DT_BOOL);
   ae_matrix_set_length(&b, a->rows, a->cols);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         b.xyB[i][j] = a->xyB[i][j];
      }
   }
   ae_matrix_set_length(a, b.cols, b.rows);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         a->xyB[j][i] = b.xyB[i][j];
      }
   }
   ae_frame_leave();
}

// Generate MxN matrix with elements set to "sin(3*I+5*J) > 0"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugb2outsin(const ae_int_t m, const ae_int_t n, boolean_2d_array &a);
void xdebugb2outsin(ae_int_t m, ae_int_t n, BMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(a);
   ae_matrix_set_length(a, m, n);
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyB[i][j] = sin(3.0 * i + 5.0 * j) > 0.0;
      }
   }
}

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: ae_int_t xdebugi2sum(const integer_2d_array &a);
ae_int_t xdebugi2sum(ZMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t result;
   result = 0;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         result += a->xyZ[i][j];
      }
   }
   return result;
}

// Replace all values in array by -a[i,j]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugi2neg(const integer_2d_array &a);
void xdebugi2neg(ZMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyZ[i][j] = -a->xyZ[i][j];
      }
   }
}

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugi2transpose(integer_2d_array &a);
void xdebugi2transpose(ZMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   NewMatrix(b, 0, 0, DT_INT);
   ae_matrix_set_length(&b, a->rows, a->cols);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         b.xyZ[i][j] = a->xyZ[i][j];
      }
   }
   ae_matrix_set_length(a, b.cols, b.rows);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         a->xyZ[j][i] = b.xyZ[i][j];
      }
   }
   ae_frame_leave();
}

// Generate MxN matrix with elements set to "sign(sin(3*I+5*J))"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugi2outsin(const ae_int_t m, const ae_int_t n, integer_2d_array &a);
void xdebugi2outsin(ae_int_t m, ae_int_t n, ZMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(a);
   ae_matrix_set_length(a, m, n);
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyZ[i][j] = sign(sin(3.0 * i + 5.0 * j));
      }
   }
}

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: double xdebugr2sum(const real_2d_array &a);
double xdebugr2sum(RMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   double result;
   result = 0.0;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         result += a->xyR[i][j];
      }
   }
   return result;
}

// Replace all values in array by -a[i,j]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugr2neg(const real_2d_array &a);
void xdebugr2neg(RMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyR[i][j] = -a->xyR[i][j];
      }
   }
}

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugr2transpose(real_2d_array &a);
void xdebugr2transpose(RMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   NewMatrix(b, 0, 0, DT_REAL);
   ae_matrix_set_length(&b, a->rows, a->cols);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         b.xyR[i][j] = a->xyR[i][j];
      }
   }
   ae_matrix_set_length(a, b.cols, b.rows);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         a->xyR[j][i] = b.xyR[i][j];
      }
   }
   ae_frame_leave();
}

// Generate MxN matrix with elements set to "sin(3*I+5*J)"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugr2outsin(const ae_int_t m, const ae_int_t n, real_2d_array &a);
void xdebugr2outsin(ae_int_t m, ae_int_t n, RMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(a);
   ae_matrix_set_length(a, m, n);
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyR[i][j] = sin(3.0 * i + 5.0 * j);
      }
   }
}

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: complex xdebugc2sum(const complex_2d_array &a);
complex xdebugc2sum(CMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   complex result;
   result = complex_from_i(0);
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         result = ae_c_add(result, a->xyC[i][j]);
      }
   }
   return result;
}

// Replace all values in array by -a[i,j]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugc2neg(const complex_2d_array &a);
void xdebugc2neg(CMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyC[i][j] = ae_c_neg(a->xyC[i][j]);
      }
   }
}

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugc2transpose(complex_2d_array &a);
void xdebugc2transpose(CMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   NewMatrix(b, 0, 0, DT_COMPLEX);
   ae_matrix_set_length(&b, a->rows, a->cols);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         b.xyC[i][j] = a->xyC[i][j];
      }
   }
   ae_matrix_set_length(a, b.cols, b.rows);
   for (i = 0; i < b.rows; i++) {
      for (j = 0; j < b.cols; j++) {
         a->xyC[j][i] = b.xyC[i][j];
      }
   }
   ae_frame_leave();
}

// Generate MxN matrix with elements set to "sin(3*I+5*J),cos(3*I+5*J)"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: void xdebugc2outsincos(const ae_int_t m, const ae_int_t n, complex_2d_array &a);
void xdebugc2outsincos(ae_int_t m, ae_int_t n, CMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(a);
   ae_matrix_set_length(a, m, n);
   for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
         a->xyC[i][j] = complex_from_d(sin(3.0 * i + 5.0 * j), cos(3.0 * i + 5.0 * j));
      }
   }
}

// Returns sum of a[i,j]*(1+b[i,j]) such that c[i,j] is True
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
// API: double xdebugmaskedbiasedproductsum(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const real_2d_array &b, const boolean_2d_array &c);
double xdebugmaskedbiasedproductsum(ae_int_t m, ae_int_t n, RMatrix *a, RMatrix *b, BMatrix *c) {
   ae_int_t i;
   ae_int_t j;
   double result;
   ae_assert(m >= a->rows, "Assertion failed");
   ae_assert(m >= b->rows, "Assertion failed");
   ae_assert(m >= c->rows, "Assertion failed");
   ae_assert(n >= a->cols, "Assertion failed");
   ae_assert(n >= b->cols, "Assertion failed");
   ae_assert(n >= c->cols, "Assertion failed");
   result = 0.0;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (c->xyB[i][j]) {
            result += a->xyR[i][j] * (1 + b->xyR[i][j]);
         }
      }
   }
   return result;
}

void xdebugrecord1_init(void *_p, bool make_automatic) {
   xdebugrecord1 *p = (xdebugrecord1 *)_p;
   ae_vector_init(&p->a, 0, DT_REAL, make_automatic);
}

void xdebugrecord1_copy(void *_dst, void *_src, bool make_automatic) {
   xdebugrecord1 *dst = (xdebugrecord1 *)_dst;
   xdebugrecord1 *src = (xdebugrecord1 *)_src;
   dst->i = src->i;
   dst->c = src->c;
   ae_vector_copy(&dst->a, &src->a, make_automatic);
}

void xdebugrecord1_free(void *_p, bool make_automatic) {
   xdebugrecord1 *p = (xdebugrecord1 *)_p;
   ae_vector_free(&p->a, make_automatic);
}
} // end of namespace alglib_impl

namespace alglib {
// This is a debug class intended for testing ALGLIB interface generator.
// Never use it in any real life project.
// ALGLIB: Copyright 20.07.2021 by Sergey Bochkanov
DefClass(xdebugrecord1, DecVal(i) DecComplex(c) DecVar(a))

void xdebuginitrecord1(xdebugrecord1 &rec1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebuginitrecord1(ConstT(xdebugrecord1, rec1));
   alglib_impl::ae_state_clear();
}

ae_int_t xdebugb1count(const boolean_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::xdebugb1count(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
   return Z;
}

void xdebugb1not(const boolean_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugb1not(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugb1appendcopy(boolean_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugb1appendcopy(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugb1outeven(const ae_int_t n, boolean_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugb1outeven(n, ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

ae_int_t xdebugi1sum(const integer_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::xdebugi1sum(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
   return Z;
}

void xdebugi1neg(const integer_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugi1neg(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugi1appendcopy(integer_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugi1appendcopy(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugi1outeven(const ae_int_t n, integer_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugi1outeven(n, ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

double xdebugr1sum(const real_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::xdebugr1sum(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
   return D;
}

void xdebugr1neg(const real_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugr1neg(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugr1appendcopy(real_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugr1appendcopy(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugr1outeven(const ae_int_t n, real_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugr1outeven(n, ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

complex xdebugc1sum(const complex_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(complex(0.0))
   alglib_impl::complex C = alglib_impl::xdebugc1sum(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
   return ComplexOf(C);
}

void xdebugc1neg(const complex_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugc1neg(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugc1appendcopy(complex_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugc1appendcopy(ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

void xdebugc1outeven(const ae_int_t n, complex_1d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugc1outeven(n, ConstT(ae_vector, a));
   alglib_impl::ae_state_clear();
}

ae_int_t xdebugb2count(const boolean_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::xdebugb2count(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
   return Z;
}

void xdebugb2not(const boolean_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugb2not(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugb2transpose(boolean_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugb2transpose(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugb2outsin(const ae_int_t m, const ae_int_t n, boolean_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugb2outsin(m, n, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

ae_int_t xdebugi2sum(const integer_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::xdebugi2sum(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
   return Z;
}

void xdebugi2neg(const integer_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugi2neg(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugi2transpose(integer_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugi2transpose(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugi2outsin(const ae_int_t m, const ae_int_t n, integer_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugi2outsin(m, n, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

double xdebugr2sum(const real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::xdebugr2sum(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
   return D;
}

void xdebugr2neg(const real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugr2neg(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugr2transpose(real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugr2transpose(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugr2outsin(const ae_int_t m, const ae_int_t n, real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugr2outsin(m, n, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

complex xdebugc2sum(const complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch(complex(0.0))
   alglib_impl::complex C = alglib_impl::xdebugc2sum(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
   return ComplexOf(C);
}

void xdebugc2neg(const complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugc2neg(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugc2transpose(complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugc2transpose(ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void xdebugc2outsincos(const ae_int_t m, const ae_int_t n, complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::xdebugc2outsincos(m, n, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

double xdebugmaskedbiasedproductsum(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const real_2d_array &b, const boolean_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::xdebugmaskedbiasedproductsum(m, n, ConstT(ae_matrix, a), ConstT(ae_matrix, b), ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === NEARESTNEIGHBOR Package ===
// Depends on: (AlgLibInternal) SCODES, TSORT
namespace alglib_impl {
static const ae_int_t nearestneighbor_splitnodesize = 6;
static const ae_int_t nearestneighbor_kdtreefirstversion = 0;

// This function allocates all dataset-independend array  fields  of  KDTree,
// i.e.  such  array  fields  that  their dimensions do not depend on dataset
// size.
//
// This function do not sets KDT.NX or KDT.NY - it just allocates arrays
// ALGLIB: Copyright 14.03.2011 by Sergey Bochkanov
static void nearestneighbor_kdtreeallocdatasetindependent(kdtree *kdt, ae_int_t nx, ae_int_t ny) {
   ae_assert(kdt->n > 0, "KDTreeAllocDatasetIndependent: internal error");
   ae_vector_set_length(&kdt->boxmin, nx);
   ae_vector_set_length(&kdt->boxmax, nx);
}

// This function allocates all dataset-dependent array fields of KDTree, i.e.
// such array fields that their dimensions depend on dataset size.
//
// This function do not sets KDT.N, KDT.NX or KDT.NY -
// it just allocates arrays.
// ALGLIB: Copyright 14.03.2011 by Sergey Bochkanov
static void nearestneighbor_kdtreeallocdatasetdependent(kdtree *kdt, ae_int_t n, ae_int_t nx, ae_int_t ny) {
   ae_assert(n > 0, "KDTreeAllocDatasetDependent: internal error");
   ae_matrix_set_length(&kdt->xy, n, 2 * nx + ny);
   ae_vector_set_length(&kdt->tags, n);
   ae_vector_set_length(&kdt->nodes, nearestneighbor_splitnodesize * 2 * n);
   ae_vector_set_length(&kdt->splits, 2 * n);
}

// Rearranges nodes [I1,I2) using partition in D-th dimension with S as threshold.
// Returns split position I3: [I1,I3) and [I3,I2) are created as result.
//
// This subroutine doesn't create tree structures, just rearranges nodes.
static ae_int_t nearestneighbor_kdtreesplit(kdtree *kdt, ae_int_t i1, ae_int_t i2, ae_int_t d, double s) {
   ae_int_t i;
   ae_int_t ileft;
   ae_int_t iright;
   ae_assert(kdt->n > 0, "KDTreeSplit: internal error");
// split XY/Tags in two parts:
// * [ILeft,IRight] is non-processed part of XY/Tags
//
// After cycle is done, we have Ileft == IRight. We deal with
// this element separately.
//
// After this, [I1,ILeft) contains left part, and [ILeft,I2)
// contains right part.
   ileft = i1;
   iright = i2 - 1;
   while (ileft < iright) {
      if (kdt->xy.xyR[ileft][d] <= s) {
      // XY[ILeft] is on its place.
      // Advance ILeft.
         ileft++;
      } else {
      // XY[ILeft,..] must be at IRight.
      // Swap and advance IRight.
         for (i = 0; i < 2 * kdt->nx + kdt->ny; i++) {
            swapr(&kdt->xy.xyR[ileft][i], &kdt->xy.xyR[iright][i]);
         }
         swapi(&kdt->tags.xZ[ileft], &kdt->tags.xZ[iright]);
         iright--;
      }
   }
   if (kdt->xy.xyR[ileft][d] <= s) {
      ileft++;
   } else {
      iright--;
   }
   return ileft;
}

// Recursive kd-tree generation subroutine.
//
// Parameters:
//     KDT         tree
//     NodesOffs   unused part of Nodes[] which must be filled by tree
//     SplitsOffs  unused part of Splits[]
//     I1, I2      points from [I1,I2) are processed
//
// NodesOffs[] and SplitsOffs[] must be large enough.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
static void nearestneighbor_kdtreegeneratetreerec(kdtree *kdt, ae_int_t *nodesoffs, ae_int_t *splitsoffs, ae_int_t i1, ae_int_t i2, ae_int_t maxleafsize) {
   ae_int_t n;
   ae_int_t nx;
   ae_int_t ny;
   ae_int_t i;
   ae_int_t oldoffs;
   ae_int_t i3;
   ae_int_t cntless;
   ae_int_t cntgreater;
   double minv;
   double maxv;
   ae_int_t minidx;
   ae_int_t maxidx;
   ae_int_t d;
   double ds;
   double s;
   double v;
   double v0;
   double v1;
   ae_assert(kdt->n > 0, "KDTreeGenerateTreeRec: internal error");
   ae_assert(i2 > i1, "KDTreeGenerateTreeRec: internal error");
// Generate leaf if needed
   if (i2 - i1 <= maxleafsize) {
      kdt->nodes.xZ[*nodesoffs] = i2 - i1;
      kdt->nodes.xZ[*nodesoffs + 1] = i1;
      *nodesoffs += 2;
      return;
   }
// Load values for easier access
   nx = kdt->nx;
   ny = kdt->ny;
// Select dimension to split:
// * D is a dimension number
// In case bounding box has zero size, we enforce creation of the leaf node.
   d = 0;
   ds = kdt->innerbuf.curboxmax.xR[0] - kdt->innerbuf.curboxmin.xR[0];
   for (i = 1; i < nx; i++) {
      v = kdt->innerbuf.curboxmax.xR[i] - kdt->innerbuf.curboxmin.xR[i];
      if (v > ds) {
         ds = v;
         d = i;
      }
   }
   if (ds == 0.0) {
      kdt->nodes.xZ[*nodesoffs] = i2 - i1;
      kdt->nodes.xZ[*nodesoffs + 1] = i1;
      *nodesoffs += 2;
      return;
   }
// Select split position S using sliding midpoint rule,
// rearrange points into [I1,I3) and [I3,I2).
//
// In case all points has same value of D-th component
// (MinV == MaxV) we enforce D-th dimension of bounding
// box to become exactly zero and repeat tree construction.
   s = kdt->innerbuf.curboxmin.xR[d] + 0.5 * ds;
   ae_v_move(kdt->innerbuf.buf.xR, 1, &kdt->xy.xyR[i1][d], kdt->xy.stride, i2 - i1);
   n = i2 - i1;
   cntless = 0;
   cntgreater = 0;
   minv = kdt->innerbuf.buf.xR[0];
   maxv = kdt->innerbuf.buf.xR[0];
   minidx = i1;
   maxidx = i1;
   for (i = 0; i < n; i++) {
      v = kdt->innerbuf.buf.xR[i];
      if (v < minv) {
         minv = v;
         minidx = i1 + i;
      }
      if (v > maxv) {
         maxv = v;
         maxidx = i1 + i;
      }
      if (v < s) {
         cntless++;
      }
      if (v > s) {
         cntgreater++;
      }
   }
   if (minv == maxv) {
   // In case all points has same value of D-th component
   // (MinV == MaxV) we enforce D-th dimension of bounding
   // box to become exactly zero and repeat tree construction.
      v0 = kdt->innerbuf.curboxmin.xR[d];
      v1 = kdt->innerbuf.curboxmax.xR[d];
      kdt->innerbuf.curboxmin.xR[d] = minv;
      kdt->innerbuf.curboxmax.xR[d] = maxv;
      nearestneighbor_kdtreegeneratetreerec(kdt, nodesoffs, splitsoffs, i1, i2, maxleafsize);
      kdt->innerbuf.curboxmin.xR[d] = v0;
      kdt->innerbuf.curboxmax.xR[d] = v1;
      return;
   }
   if (cntless > 0 && cntgreater > 0) {
   // normal midpoint split
      i3 = nearestneighbor_kdtreesplit(kdt, i1, i2, d, s);
   } else {
   // sliding midpoint
      if (cntless == 0) {
      // 1. move split to MinV,
      // 2. place one point to the left bin (move to I1),
      //    others - to the right bin
         s = minv;
         if (minidx != i1) {
            for (i = 0; i < 2 * nx + ny; i++) {
               swapr(&kdt->xy.xyR[minidx][i], &kdt->xy.xyR[i1][i]);
            }
            swapi(&kdt->tags.xZ[minidx], &kdt->tags.xZ[i1]);
         }
         i3 = i1 + 1;
      } else {
      // 1. move split to MaxV,
      // 2. place one point to the right bin (move to I2-1),
      //    others - to the left bin
         s = maxv;
         if (maxidx != i2 - 1) {
            for (i = 0; i < 2 * nx + ny; i++) {
               swapr(&kdt->xy.xyR[maxidx][i], &kdt->xy.xyR[i2 - 1][i]);
            }
            swapi(&kdt->tags.xZ[maxidx], &kdt->tags.xZ[i2 - 1]);
         }
         i3 = i2 - 1;
      }
   }
// Generate 'split' node
   kdt->nodes.xZ[*nodesoffs] = 0;
   kdt->nodes.xZ[*nodesoffs + 1] = d;
   kdt->nodes.xZ[*nodesoffs + 2] = *splitsoffs;
   kdt->splits.xR[*splitsoffs] = s;
   oldoffs = *nodesoffs;
   *nodesoffs += nearestneighbor_splitnodesize;
   ++*splitsoffs;
// Recursive generation:
// * update CurBox
// * call subroutine
// * restore CurBox
   kdt->nodes.xZ[oldoffs + 3] = *nodesoffs;
   v = kdt->innerbuf.curboxmax.xR[d];
   kdt->innerbuf.curboxmax.xR[d] = s;
   nearestneighbor_kdtreegeneratetreerec(kdt, nodesoffs, splitsoffs, i1, i3, maxleafsize);
   kdt->innerbuf.curboxmax.xR[d] = v;
   kdt->nodes.xZ[oldoffs + 4] = *nodesoffs;
   v = kdt->innerbuf.curboxmin.xR[d];
   kdt->innerbuf.curboxmin.xR[d] = s;
   nearestneighbor_kdtreegeneratetreerec(kdt, nodesoffs, splitsoffs, i3, i2, maxleafsize);
   kdt->innerbuf.curboxmin.xR[d] = v;
// Zero-fill unused portions of the node (avoid false warnings by Valgrind
// about attempt to serialize uninitialized values)
   ae_assert(nearestneighbor_splitnodesize == 6, "KDTreeGenerateTreeRec: node size has unexpectedly changed");
   kdt->nodes.xZ[oldoffs + 5] = 0;
}

// This function creates buffer  structure  which  can  be  used  to  perform
// parallel KD-tree requests.
//
// KD-tree subpackage provides two sets of request functions - ones which use
// internal buffer of KD-tree object  (these  functions  are  single-threaded
// because they use same buffer, which can not shared between  threads),  and
// ones which use external buffer.
//
// This function is used to initialize external buffer.
//
// Inputs:
//     KDT         -   KD-tree which is associated with newly created buffer
//
// Outputs:
//     Buf         -   external buffer.
//
// IMPORTANT: KD-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use buffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
// API: void kdtreecreaterequestbuffer(const kdtree &kdt, kdtreerequestbuffer &buf);
void kdtreecreaterequestbuffer(kdtree *kdt, kdtreerequestbuffer *buf) {
   SetObj(kdtreerequestbuffer, buf);
   ae_vector_set_length(&buf->x, kdt->nx);
   ae_vector_set_length(&buf->boxmin, kdt->nx);
   ae_vector_set_length(&buf->boxmax, kdt->nx);
   ae_vector_set_length(&buf->idx, kdt->n);
   ae_vector_set_length(&buf->r, kdt->n);
   ae_vector_set_length(&buf->buf, imax2(kdt->n, kdt->nx));
   ae_vector_set_length(&buf->curboxmin, kdt->nx);
   ae_vector_set_length(&buf->curboxmax, kdt->nx);
   buf->kcur = 0;
}

// KD-tree creation
//
// This  subroutine  creates  KD-tree  from set of X-values, integer tags and
// optional Y-values
//
// Inputs:
//     XY      -   dataset, array[0..N-1,0..NX+NY-1].
//                 one row corresponds to one point.
//                 first NX columns contain X-values, next NY (NY may be zero)
//                 columns may contain associated Y-values
//     Tags    -   tags, array[0..N-1], contains integer tags associated
//                 with points.
//     N       -   number of points, N >= 0
//     NX      -   space dimension, NX >= 1.
//     NY      -   number of optional Y-values, NY >= 0.
//     NormType-   norm type:
//                 * 0 denotes infinity-norm
//                 * 1 denotes 1-norm
//                 * 2 denotes 2-norm (Euclidean norm)
//
// Outputs:
//     KDT     -   KD-tree
//
// NOTES
//
// 1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
//    requirements.
// 2. Although KD-trees may be used with any combination of N  and  NX,  they
//    are more efficient than brute-force search only when N >> 4^NX. So they
//    are most useful in low-dimensional tasks (NX == 2, NX == 3). NX == 1  is another
//    inefficient case, because  simple  binary  search  (without  additional
//    structures) is much more efficient in such tasks than KD-trees.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
// API: void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuildtagged(RMatrix *xy, ZVector *tags, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree *kdt) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t nodesoffs;
   ae_int_t splitsoffs;
   SetObj(kdtree, kdt);
   ae_assert(n >= 0, "KDTreeBuildTagged: N < 0");
   ae_assert(nx >= 1, "KDTreeBuildTagged: NX < 1");
   ae_assert(ny >= 0, "KDTreeBuildTagged: NY < 0");
   ae_assert(normtype >= 0 && normtype <= 2, "KDTreeBuildTagged: incorrect NormType");
   ae_assert(xy->rows >= n, "KDTreeBuildTagged: rows(X) < N");
   ae_assert(xy->cols >= nx + ny || n == 0, "KDTreeBuildTagged: cols(X) < NX+NY");
   ae_assert(apservisfinitematrix(xy, n, nx + ny), "KDTreeBuildTagged: XY contains infinite or NaN values");
// initialize
   kdt->n = n;
   kdt->nx = nx;
   kdt->ny = ny;
   kdt->normtype = normtype;
   kdt->innerbuf.kcur = 0;
// N == 0 => quick exit
   if (n == 0) {
      return;
   }
// Allocate
   nearestneighbor_kdtreeallocdatasetindependent(kdt, nx, ny);
   nearestneighbor_kdtreeallocdatasetdependent(kdt, n, nx, ny);
   kdtreecreaterequestbuffer(kdt, &kdt->innerbuf);
// Initial fill
   for (i = 0; i < n; i++) {
      ae_v_move(kdt->xy.xyR[i], 1, xy->xyR[i], 1, nx);
      ae_v_move(&kdt->xy.xyR[i][nx], 1, xy->xyR[i], 1, nx + ny);
      kdt->tags.xZ[i] = tags->xZ[i];
   }
// Determine bounding box
   ae_v_move(kdt->boxmin.xR, 1, kdt->xy.xyR[0], 1, nx);
   ae_v_move(kdt->boxmax.xR, 1, kdt->xy.xyR[0], 1, nx);
   for (i = 1; i < n; i++) {
      for (j = 0; j < nx; j++) {
         kdt->boxmin.xR[j] = rmin2(kdt->boxmin.xR[j], kdt->xy.xyR[i][j]);
         kdt->boxmax.xR[j] = rmax2(kdt->boxmax.xR[j], kdt->xy.xyR[i][j]);
      }
   }
// Generate tree
   nodesoffs = 0;
   splitsoffs = 0;
   ae_v_move(kdt->innerbuf.curboxmin.xR, 1, kdt->boxmin.xR, 1, nx);
   ae_v_move(kdt->innerbuf.curboxmax.xR, 1, kdt->boxmax.xR, 1, nx);
   nearestneighbor_kdtreegeneratetreerec(kdt, &nodesoffs, &splitsoffs, 0, n, 8);
   ivectorresize(&kdt->nodes, nodesoffs);
   rvectorresize(&kdt->splits, splitsoffs);
}

// KD-tree creation
//
// This subroutine creates KD-tree from set of X-values and optional Y-values
//
// Inputs:
//     XY      -   dataset, array[0..N-1,0..NX+NY-1].
//                 one row corresponds to one point.
//                 first NX columns contain X-values, next NY (NY may be zero)
//                 columns may contain associated Y-values
//     N       -   number of points, N >= 0.
//     NX      -   space dimension, NX >= 1.
//     NY      -   number of optional Y-values, NY >= 0.
//     NormType-   norm type:
//                 * 0 denotes infinity-norm
//                 * 1 denotes 1-norm
//                 * 2 denotes 2-norm (Euclidean norm)
//
// Outputs:
//     KDT     -   KD-tree
//
// NOTES
//
// 1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
//    requirements.
// 2. Although KD-trees may be used with any combination of N  and  NX,  they
//    are more efficient than brute-force search only when N >> 4^NX. So they
//    are most useful in low-dimensional tasks (NX == 2, NX == 3). NX == 1  is another
//    inefficient case, because  simple  binary  search  (without  additional
//    structures) is much more efficient in such tasks than KD-trees.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreebuild(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
// API: void kdtreebuild(const real_2d_array &xy, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuild(RMatrix *xy, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree *kdt) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   SetObj(kdtree, kdt);
   NewVector(tags, 0, DT_INT);
   ae_assert(n >= 0, "KDTreeBuild: N < 0");
   ae_assert(nx >= 1, "KDTreeBuild: NX < 1");
   ae_assert(ny >= 0, "KDTreeBuild: NY < 0");
   ae_assert(normtype >= 0 && normtype <= 2, "KDTreeBuild: incorrect NormType");
   ae_assert(xy->rows >= n, "KDTreeBuild: rows(X) < N");
   ae_assert(xy->cols >= nx + ny || n == 0, "KDTreeBuild: cols(X) < NX+NY");
   ae_assert(apservisfinitematrix(xy, n, nx + ny), "KDTreeBuild: XY contains infinite or NaN values");
   if (n > 0) {
      ae_vector_set_length(&tags, n);
      for (i = 0; i < n; i++) {
         tags.xZ[i] = 0;
      }
   }
   kdtreebuildtagged(xy, &tags, n, nx, ny, normtype, kdt);
   ae_frame_leave();
}

// This  function   checks  consistency  of  request  buffer  structure  with
// dimensions of kd-tree object.
// ALGLIB: Copyright 02.04.2016 by Sergey Bochkanov
static void nearestneighbor_checkrequestbufferconsistency(kdtree *kdt, kdtreerequestbuffer *buf) {
   ae_assert(buf->x.cnt >= kdt->nx, "KDTree: dimensions of kdtreerequestbuffer are inconsistent with kdtree structure");
   ae_assert(buf->idx.cnt >= kdt->n, "KDTree: dimensions of kdtreerequestbuffer are inconsistent with kdtree structure");
   ae_assert(buf->r.cnt >= kdt->n, "KDTree: dimensions of kdtreerequestbuffer are inconsistent with kdtree structure");
   ae_assert(buf->buf.cnt >= imax2(kdt->n, kdt->nx), "KDTree: dimensions of kdtreerequestbuffer are inconsistent with kdtree structure");
   ae_assert(buf->curboxmin.cnt >= kdt->nx, "KDTree: dimensions of kdtreerequestbuffer are inconsistent with kdtree structure");
   ae_assert(buf->curboxmax.cnt >= kdt->nx, "KDTree: dimensions of kdtreerequestbuffer are inconsistent with kdtree structure");
}

// Copies X[] to Buf.X[]
// Loads distance from X[] to bounding box.
// Initializes Buf.CurBox[].
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
static void nearestneighbor_kdtreeinitbox(kdtree *kdt, RVector *x, kdtreerequestbuffer *buf) {
   ae_int_t i;
   double vx;
   double vmin;
   double vmax;
   ae_assert(kdt->n > 0, "KDTreeInitBox: internal error");
// calculate distance from point to current bounding box
   buf->curdist = 0.0;
   if (kdt->normtype == 0) {
      for (i = 0; i < kdt->nx; i++) {
         vx = x->xR[i];
         vmin = kdt->boxmin.xR[i];
         vmax = kdt->boxmax.xR[i];
         buf->x.xR[i] = vx;
         buf->curboxmin.xR[i] = vmin;
         buf->curboxmax.xR[i] = vmax;
         if (vx < vmin) {
            buf->curdist = rmax2(buf->curdist, vmin - vx);
         } else {
            if (vx > vmax) {
               buf->curdist = rmax2(buf->curdist, vx - vmax);
            }
         }
      }
   }
   if (kdt->normtype == 1) {
      for (i = 0; i < kdt->nx; i++) {
         vx = x->xR[i];
         vmin = kdt->boxmin.xR[i];
         vmax = kdt->boxmax.xR[i];
         buf->x.xR[i] = vx;
         buf->curboxmin.xR[i] = vmin;
         buf->curboxmax.xR[i] = vmax;
         if (vx < vmin) {
            buf->curdist += vmin - vx;
         } else {
            if (vx > vmax) {
               buf->curdist += vx - vmax;
            }
         }
      }
   }
   if (kdt->normtype == 2) {
      for (i = 0; i < kdt->nx; i++) {
         vx = x->xR[i];
         vmin = kdt->boxmin.xR[i];
         vmax = kdt->boxmax.xR[i];
         buf->x.xR[i] = vx;
         buf->curboxmin.xR[i] = vmin;
         buf->curboxmax.xR[i] = vmax;
         if (vx < vmin) {
            buf->curdist += sqr(vmin - vx);
         } else {
            if (vx > vmax) {
               buf->curdist += sqr(vx - vmax);
            }
         }
      }
   }
}

// Recursive subroutine for NN queries.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
static void nearestneighbor_kdtreequerynnrec(kdtree *kdt, kdtreerequestbuffer *buf, ae_int_t offs) {
   double ptdist;
   ae_int_t i;
   ae_int_t j;
   ae_int_t nx;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t d;
   double s;
   double v;
   double t1;
   ae_int_t childbestoffs;
   ae_int_t childworstoffs;
   ae_int_t childoffs;
   double prevdist;
   bool todive;
   bool bestisleft;
   bool updatemin;
   ae_assert(kdt->n > 0, "KDTreeQueryNNRec: internal error");
// Leaf node.
// Process points.
   if (kdt->nodes.xZ[offs] > 0) {
      i1 = kdt->nodes.xZ[offs + 1];
      i2 = i1 + kdt->nodes.xZ[offs];
      for (i = i1; i < i2; i++) {
      // Calculate distance
         ptdist = 0.0;
         nx = kdt->nx;
         if (kdt->normtype == 0) {
            for (j = 0; j < nx; j++) {
               ptdist = rmax2(ptdist, fabs(kdt->xy.xyR[i][j] - buf->x.xR[j]));
            }
         }
         if (kdt->normtype == 1) {
            for (j = 0; j < nx; j++) {
               ptdist += fabs(kdt->xy.xyR[i][j] - buf->x.xR[j]);
            }
         }
         if (kdt->normtype == 2) {
            for (j = 0; j < nx; j++) {
               ptdist += sqr(kdt->xy.xyR[i][j] - buf->x.xR[j]);
            }
         }
      // Skip points with zero distance if self-matches are turned off
         if (ptdist == 0 && !buf->selfmatch) {
            continue;
         }
      // We CAN'T process point if R-criterion isn't satisfied,
      // i.e. (RNeeded != 0) AND (PtDist > R).
         if (buf->rneeded == 0 || ptdist <= buf->rneeded) {
         // R-criterion is satisfied, we must either:
         // * replace worst point, if (KNeeded != 0) AND (KCur == KNeeded)
         //   (or skip, if worst point is better)
         // * add point without replacement otherwise
            if (buf->kcur < buf->kneeded || buf->kneeded == 0) {
            // add current point to heap without replacement
               tagheappushi(&buf->r, &buf->idx, &buf->kcur, ptdist, i);
            } else {
            // New points are added or not, depending on their distance.
            // If added, they replace element at the top of the heap
               if (ptdist < buf->r.xR[0]) {
                  if (buf->kneeded == 1) {
                     buf->idx.xZ[0] = i;
                     buf->r.xR[0] = ptdist;
                  } else {
                     tagheapreplacetopi(&buf->r, &buf->idx, buf->kneeded, ptdist, i);
                  }
               }
            }
         }
      }
      return;
   }
// Simple split
   if (kdt->nodes.xZ[offs] == 0) {
   // Load:
   // * D  dimension to split
   // * S  split position
      d = kdt->nodes.xZ[offs + 1];
      s = kdt->splits.xR[kdt->nodes.xZ[offs + 2]];
   // Calculate:
   // * ChildBestOffs      child box with best chances
   // * ChildWorstOffs     child box with worst chances
      if (buf->x.xR[d] <= s) {
         childbestoffs = kdt->nodes.xZ[offs + 3];
         childworstoffs = kdt->nodes.xZ[offs + 4];
         bestisleft = true;
      } else {
         childbestoffs = kdt->nodes.xZ[offs + 4];
         childworstoffs = kdt->nodes.xZ[offs + 3];
         bestisleft = false;
      }
   // Navigate through childs
      for (i = 0; i <= 1; i++) {
      // Select child to process:
      // * ChildOffs      current child offset in Nodes[]
      // * UpdateMin      whether minimum or maximum value
      //                  of bounding box is changed on update
         if (i == 0) {
            childoffs = childbestoffs;
            updatemin = !bestisleft;
         } else {
            updatemin = bestisleft;
            childoffs = childworstoffs;
         }
      // Update bounding box and current distance
         if (updatemin) {
            prevdist = buf->curdist;
            t1 = buf->x.xR[d];
            v = buf->curboxmin.xR[d];
            if (t1 <= s) {
               if (kdt->normtype == 0) {
                  buf->curdist = rmax2(buf->curdist, s - t1);
               }
               if (kdt->normtype == 1) {
                  buf->curdist -= rmax2(v, t1) - s;
               }
               if (kdt->normtype == 2) {
                  buf->curdist -= sqr(rmax2(v - t1, 0.0)) - sqr(s - t1);
               }
            }
            buf->curboxmin.xR[d] = s;
         } else {
            prevdist = buf->curdist;
            t1 = buf->x.xR[d];
            v = buf->curboxmax.xR[d];
            if (t1 >= s) {
               if (kdt->normtype == 0) {
                  buf->curdist = rmax2(buf->curdist, t1 - s);
               }
               if (kdt->normtype == 1) {
                  buf->curdist -= s - rmin2(v, t1);
               }
               if (kdt->normtype == 2) {
                  buf->curdist -= sqr(rmax2(t1 - v, 0.0)) - sqr(t1 - s);
               }
            }
            buf->curboxmax.xR[d] = s;
         }
      // Decide: to dive into cell or not to dive
         if (buf->rneeded != 0 && buf->curdist > buf->rneeded) {
            todive = false;
         } else {
            if (buf->kcur < buf->kneeded || buf->kneeded == 0) {
            // KCur < KNeeded (i.e. not all points are found)
               todive = true;
            } else {
            // KCur == KNeeded, decide to dive or not to dive
            // using point position relative to bounding box.
               todive = buf->curdist <= buf->r.xR[0] * buf->approxf;
            }
         }
         if (todive) {
            nearestneighbor_kdtreequerynnrec(kdt, buf, childoffs);
         }
      // Restore bounding box and distance
         if (updatemin) {
            buf->curboxmin.xR[d] = v;
         } else {
            buf->curboxmax.xR[d] = v;
         }
         buf->curdist = prevdist;
      }
      return;
   }
}

// K-NN query: approximate K nearest neighbors, using thread-local buffer.
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//     Eps         -   approximation factor, Eps >= 0. eps-approximate  nearest
//                     neighbor  is  a  neighbor  whose distance from X is at
//                     most (1+eps) times distance of true nearest neighbor.
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// NOTES
//     significant performance gain may be achieved only when Eps  is  is  on
//     the order of magnitude of 1 or larger.
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
// API: ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
// API: ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const double eps);
ae_int_t kdtreetsqueryaknn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, ae_int_t k, bool selfmatch, double eps) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t result;
   ae_assert(k > 0, "KDTreeTsQueryAKNN: incorrect K!");
   ae_assert(eps >= 0.0, "KDTreeTsQueryAKNN: incorrect Eps!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeTsQueryAKNN: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeTsQueryAKNN: X contains infinite or NaN values!");
// Handle special case: KDT.N == 0
   if (kdt->n == 0) {
      buf->kcur = 0;
      result = 0;
      return result;
   }
// Check consistency of request buffer
   nearestneighbor_checkrequestbufferconsistency(kdt, buf);
// Prepare parameters
   k = imin2(k, kdt->n);
   buf->kneeded = k;
   buf->rneeded = 0.0;
   buf->selfmatch = selfmatch;
   if (kdt->normtype == 2) {
      buf->approxf = 1.0 / sqr(1 + eps);
   } else {
      buf->approxf = 1.0 / (1 + eps);
   }
   buf->kcur = 0;
// calculate distance from point to current bounding box
   nearestneighbor_kdtreeinitbox(kdt, x, buf);
// call recursive search
// results are returned as heap
   nearestneighbor_kdtreequerynnrec(kdt, buf, 0);
// pop from heap to generate ordered representation
//
// last element is non pop'ed because it is already in
// its place
   result = buf->kcur;
   j = buf->kcur;
   for (i = buf->kcur; i >= 2; i--) {
      tagheappopi(&buf->r, &buf->idx, &j);
   }
   return result;
}

// K-NN query: approximate K nearest neighbors
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: KDTreeTsQueryAKNN() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//     Eps         -   approximation factor, Eps >= 0. eps-approximate  nearest
//                     neighbor  is  a  neighbor  whose distance from X is at
//                     most (1+eps) times distance of true nearest neighbor.
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// NOTES
//     significant performance gain may be achieved only when Eps  is  is  on
//     the order of magnitude of 1 or larger.
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// these results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
// API: ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const double eps);
ae_int_t kdtreequeryaknn(kdtree *kdt, RVector *x, ae_int_t k, bool selfmatch, double eps) {
   ae_int_t result;
   result = kdtreetsqueryaknn(kdt, &kdt->innerbuf, x, k, selfmatch, eps);
   return result;
}

// K-NN query: K nearest neighbors, using external thread-local buffer.
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   kd-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
// API: ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
// API: ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k);
ae_int_t kdtreetsqueryknn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, ae_int_t k, bool selfmatch) {
   ae_int_t result;
   ae_assert(k >= 1, "KDTreeTsQueryKNN: K < 1!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeTsQueryKNN: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeTsQueryKNN: X contains infinite or NaN values!");
   result = kdtreetsqueryaknn(kdt, buf, x, k, selfmatch, 0.0);
   return result;
}

// K-NN query: K nearest neighbors
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: KDTreeTsQueryKNN() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// these results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
// API: ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k);
ae_int_t kdtreequeryknn(kdtree *kdt, RVector *x, ae_int_t k, bool selfmatch) {
   ae_int_t result;
   ae_assert(k >= 1, "KDTreeQueryKNN: K < 1!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeQueryKNN: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeQueryKNN: X contains infinite or NaN values!");
   result = kdtreetsqueryaknn(kdt, &kdt->innerbuf, x, k, selfmatch, 0.0);
   return result;
}

// R-NN query: all points within  R-sphere  centered  at  X,  using  external
// thread-local buffer, sorted by distance between point and X (by ascending)
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// NOTE: it is also possible to perform undordered queries performed by means
//       of kdtreequeryrnnu() and kdtreetsqueryrnnu() functions. Such queries
//       are faster because we do not have to use heap structure for sorting.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
static ae_int_t nearestneighbor_tsqueryrnn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch, bool orderedbydist) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t result;
// Handle special case: KDT.N == 0
   if (kdt->n == 0) {
      buf->kcur = 0;
      result = 0;
      return result;
   }
// Check consistency of request buffer
   nearestneighbor_checkrequestbufferconsistency(kdt, buf);
// Prepare parameters
   buf->kneeded = 0;
   if (kdt->normtype != 2) {
      buf->rneeded = r;
   } else {
      buf->rneeded = sqr(r);
   }
   buf->selfmatch = selfmatch;
   buf->approxf = 1.0;
   buf->kcur = 0;
// calculate distance from point to current bounding box
   nearestneighbor_kdtreeinitbox(kdt, x, buf);
// call recursive search
// results are returned as heap
   nearestneighbor_kdtreequerynnrec(kdt, buf, 0);
   result = buf->kcur;
// pop from heap to generate ordered representation
//
// last element is not pop'ed because it is already in
// its place
   if (orderedbydist) {
      j = buf->kcur;
      for (i = buf->kcur; i >= 2; i--) {
         tagheappopi(&buf->r, &buf->idx, &j);
      }
   }
   return result;
}

// R-NN query: all points within  R-sphere  centered  at  X,  using  external
// thread-local buffer, sorted by distance between point and X (by ascending)
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// NOTE: it is also possible to perform undordered queries performed by means
//       of kdtreequeryrnnu() and kdtreetsqueryrnnu() functions. Such queries
//       are faster because we do not have to use heap structure for sorting.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
// API: ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch);
// API: ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r);
ae_int_t kdtreetsqueryrnn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch) {
   ae_int_t result;
   ae_assert(isfinite(r) && r > 0.0, "KDTreeTsQueryRNN: incorrect R!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeTsQueryRNN: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeTsQueryRNN: X contains infinite or NaN values!");
   result = nearestneighbor_tsqueryrnn(kdt, buf, x, r, selfmatch, true);
   return result;
}

// R-NN query: all points within R-sphere centered at X, ordered by  distance
// between point and X (by ascending).
//
// NOTE: it is also possible to perform undordered queries performed by means
//       of kdtreequeryrnnu() and kdtreetsqueryrnnu() functions. Such queries
//       are faster because we do not have to use heap structure for sorting.
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: kdtreetsqueryrnn() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// actual results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
// API: ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r);
ae_int_t kdtreequeryrnn(kdtree *kdt, RVector *x, double r, bool selfmatch) {
   ae_int_t result;
   ae_assert(r > 0.0, "KDTreeQueryRNN: incorrect R!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeQueryRNN: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeQueryRNN: X contains infinite or NaN values!");
   result = kdtreetsqueryrnn(kdt, &kdt->innerbuf, x, r, selfmatch);
   return result;
}

// R-NN query: all points within  R-sphere  centered  at  X,  using  external
// thread-local buffer, no ordering by distance as undicated  by  "U"  suffix
// (faster than ordered query, for large queries - significantly faster).
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// As indicated by "U" suffix, this function returns unordered results.
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
// API: ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch);
// API: ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r);
ae_int_t kdtreetsqueryrnnu(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch) {
   ae_int_t result;
   ae_assert(isfinite(r) && r > 0.0, "KDTreeTsQueryRNNU: incorrect R!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeTsQueryRNNU: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeTsQueryRNNU: X contains infinite or NaN values!");
   result = nearestneighbor_tsqueryrnn(kdt, buf, x, r, selfmatch, false);
   return result;
}

// R-NN query: all points within R-sphere  centered  at  X,  no  ordering  by
// distance as undicated by "U" suffix (faster than ordered query, for  large
// queries - significantly faster).
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: kdtreetsqueryrnn() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// actual results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
//
// As indicated by "U" suffix, this function returns unordered results.
// ALGLIB: Copyright 01.11.2018 by Sergey Bochkanov
// API: ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
// API: ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r);
ae_int_t kdtreequeryrnnu(kdtree *kdt, RVector *x, double r, bool selfmatch) {
   ae_int_t result;
   ae_assert(r > 0.0, "KDTreeQueryRNNU: incorrect R!");
   ae_assert(x->cnt >= kdt->nx, "KDTreeQueryRNNU: Length(X) < NX!");
   ae_assert(isfinitevector(x, kdt->nx), "KDTreeQueryRNNU: X contains infinite or NaN values!");
   result = kdtreetsqueryrnnu(kdt, &kdt->innerbuf, x, r, selfmatch);
   return result;
}

// Recursive subroutine for box queries.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
static void nearestneighbor_kdtreequeryboxrec(kdtree *kdt, kdtreerequestbuffer *buf, ae_int_t offs) {
   bool inbox;
   ae_int_t nx;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t i;
   ae_int_t j;
   ae_int_t d;
   double s;
   double v;
   ae_assert(kdt->n > 0, "KDTreeQueryBoxRec: internal error");
   nx = kdt->nx;
// Check that intersection of query box with bounding box is non-empty.
// This check is performed once for Offs == 0 (tree root).
   if (offs == 0) {
      for (j = 0; j < nx; j++) {
         if (buf->boxmin.xR[j] > buf->curboxmax.xR[j]) {
            return;
         }
         if (buf->boxmax.xR[j] < buf->curboxmin.xR[j]) {
            return;
         }
      }
   }
// Leaf node.
// Process points.
   if (kdt->nodes.xZ[offs] > 0) {
      i1 = kdt->nodes.xZ[offs + 1];
      i2 = i1 + kdt->nodes.xZ[offs];
      for (i = i1; i < i2; i++) {
      // Check whether point is in box or not
         inbox = true;
         for (j = 0; j < nx; j++) {
            inbox = inbox && kdt->xy.xyR[i][j] >= buf->boxmin.xR[j];
            inbox = inbox && kdt->xy.xyR[i][j] <= buf->boxmax.xR[j];
         }
         if (!inbox) {
            continue;
         }
      // Add point to unordered list
         buf->r.xR[buf->kcur] = 0.0;
         buf->idx.xZ[buf->kcur] = i;
         buf->kcur++;
      }
      return;
   }
// Simple split
   if (kdt->nodes.xZ[offs] == 0) {
   // Load:
   // * D  dimension to split
   // * S  split position
      d = kdt->nodes.xZ[offs + 1];
      s = kdt->splits.xR[kdt->nodes.xZ[offs + 2]];
   // Check lower split (S is upper bound of new bounding box)
      if (s >= buf->boxmin.xR[d]) {
         v = buf->curboxmax.xR[d];
         buf->curboxmax.xR[d] = s;
         nearestneighbor_kdtreequeryboxrec(kdt, buf, kdt->nodes.xZ[offs + 3]);
         buf->curboxmax.xR[d] = v;
      }
   // Check upper split (S is lower bound of new bounding box)
      if (s <= buf->boxmax.xR[d]) {
         v = buf->curboxmin.xR[d];
         buf->curboxmin.xR[d] = s;
         nearestneighbor_kdtreequeryboxrec(kdt, buf, kdt->nodes.xZ[offs + 4]);
         buf->curboxmin.xR[d] = v;
      }
      return;
   }
}

// Box query: all points within user-specified box, using thread-local buffer.
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     BoxMin      -   lower bounds, array[0..NX-1].
//     BoxMax      -   upper bounds, array[0..NX-1].
//
// Result:
//     number of actual neighbors found (in [0,N]).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "ts" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() returns zeros for this query
//
// NOTE: this particular query returns unordered results, because there is no
//       meaningful way of  ordering  points.  Furthermore,  no 'distance' is
//       associated with points - it is either INSIDE  or OUTSIDE (so request
//       for distances will return zeros).
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 14.05.2016 by Sergey Bochkanov
// API: ae_int_t kdtreetsquerybox(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &boxmin, const real_1d_array &boxmax);
ae_int_t kdtreetsquerybox(kdtree *kdt, kdtreerequestbuffer *buf, RVector *boxmin, RVector *boxmax) {
   ae_int_t j;
   ae_int_t result;
   ae_assert(boxmin->cnt >= kdt->nx, "KDTreeTsQueryBox: Length(BoxMin) < NX!");
   ae_assert(boxmax->cnt >= kdt->nx, "KDTreeTsQueryBox: Length(BoxMax) < NX!");
   ae_assert(isfinitevector(boxmin, kdt->nx), "KDTreeTsQueryBox: BoxMin contains infinite or NaN values!");
   ae_assert(isfinitevector(boxmax, kdt->nx), "KDTreeTsQueryBox: BoxMax contains infinite or NaN values!");
// Check consistency of request buffer
   nearestneighbor_checkrequestbufferconsistency(kdt, buf);
// Quick exit for degenerate boxes
   for (j = 0; j < kdt->nx; j++) {
      if (boxmin->xR[j] > boxmax->xR[j]) {
         buf->kcur = 0;
         result = 0;
         return result;
      }
   }
// Prepare parameters
   for (j = 0; j < kdt->nx; j++) {
      buf->boxmin.xR[j] = boxmin->xR[j];
      buf->boxmax.xR[j] = boxmax->xR[j];
      buf->curboxmin.xR[j] = boxmin->xR[j];
      buf->curboxmax.xR[j] = boxmax->xR[j];
   }
   buf->kcur = 0;
// call recursive search
   nearestneighbor_kdtreequeryboxrec(kdt, buf, 0);
   result = buf->kcur;
   return result;
}

// Box query: all points within user-specified box.
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: KDTreeTsQueryBox() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     BoxMin      -   lower bounds, array[0..NX-1].
//     BoxMax      -   upper bounds, array[0..NX-1].
//
// Result:
//     number of actual neighbors found (in [0,N]).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// these results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() returns zeros for this request
//
// NOTE: this particular query returns unordered results, because there is no
//       meaningful way of  ordering  points.  Furthermore,  no 'distance' is
//       associated with points - it is either INSIDE  or OUTSIDE (so request
//       for distances will return zeros).
// ALGLIB: Copyright 14.05.2016 by Sergey Bochkanov
// API: ae_int_t kdtreequerybox(const kdtree &kdt, const real_1d_array &boxmin, const real_1d_array &boxmax);
ae_int_t kdtreequerybox(kdtree *kdt, RVector *boxmin, RVector *boxmax) {
   ae_int_t result;
   result = kdtreetsquerybox(kdt, &kdt->innerbuf, boxmin, boxmax);
   return result;
}

// X-values from last query associated with kdtreerequestbuffer object.
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     X       -   possibly preallocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     X       -   rows are filled with X-values
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreetsqueryresultsx(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &x);
void kdtreetsqueryresultsx(kdtree *kdt, kdtreerequestbuffer *buf, RMatrix *x) {
   ae_int_t i;
   ae_int_t k;
   if (buf->kcur == 0) {
      return;
   }
   matrixsetlengthatleast(x, buf->kcur, kdt->nx);
   k = buf->kcur;
   for (i = 0; i < k; i++) {
      ae_v_move(x->xyR[i], 1, &kdt->xy.xyR[buf->idx.xZ[i]][kdt->nx], 1, kdt->nx);
   }
}

// X-values from last query.
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultsx().
//
// Inputs:
//     KDT     -   KD-tree
//     X       -   possibly preallocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     X       -   rows are filled with X-values
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultsx(const kdtree &kdt, real_2d_array &x);
void kdtreequeryresultsx(kdtree *kdt, RMatrix *x) {
   kdtreetsqueryresultsx(kdt, &kdt->innerbuf, x);
}

// X- and Y-values from last query associated with kdtreerequestbuffer object.
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     XY      -   possibly preallocated buffer. If XY is too small to store
//                 result, it is resized. If size(XY) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     XY      -   rows are filled with points: first NX columns with
//                 X-values, next NY columns - with Y-values.
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreetsqueryresultsxy(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &xy);
void kdtreetsqueryresultsxy(kdtree *kdt, kdtreerequestbuffer *buf, RMatrix *xy) {
   ae_int_t i;
   ae_int_t k;
   if (buf->kcur == 0) {
      return;
   }
   matrixsetlengthatleast(xy, buf->kcur, kdt->nx + kdt->ny);
   k = buf->kcur;
   for (i = 0; i < k; i++) {
      ae_v_move(xy->xyR[i], 1, &kdt->xy.xyR[buf->idx.xZ[i]][kdt->nx], 1, kdt->nx + kdt->ny);
   }
}

// X- and Y-values from last query
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultsxy().
//
// Inputs:
//     KDT     -   KD-tree
//     XY      -   possibly preallocated buffer. If XY is too small to store
//                 result, it is resized. If size(XY) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     XY      -   rows are filled with points: first NX columns with
//                 X-values, next NY columns - with Y-values.
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultsxy(const kdtree &kdt, real_2d_array &xy);
void kdtreequeryresultsxy(kdtree *kdt, RMatrix *xy) {
   kdtreetsqueryresultsxy(kdt, &kdt->innerbuf, xy);
}

// Tags from last query associated with kdtreerequestbuffer object.
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - KDTreeTsqueryresultstags().
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     Tags    -   possibly preallocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     Tags    -   filled with tags associated with points,
//                 or, when no tags were supplied, with zeros
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreetsqueryresultstags(const kdtree &kdt, const kdtreerequestbuffer &buf, integer_1d_array &tags);
void kdtreetsqueryresultstags(kdtree *kdt, kdtreerequestbuffer *buf, ZVector *tags) {
   ae_int_t i;
   ae_int_t k;
   if (buf->kcur == 0) {
      return;
   }
   vectorsetlengthatleast(tags, buf->kcur);
   k = buf->kcur;
   for (i = 0; i < k; i++) {
      tags->xZ[i] = kdt->tags.xZ[buf->idx.xZ[i]];
   }
}

// Tags from last query
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultstags().
//
// Inputs:
//     KDT     -   KD-tree
//     Tags    -   possibly preallocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     Tags    -   filled with tags associated with points,
//                 or, when no tags were supplied, with zeros
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultstags(const kdtree &kdt, integer_1d_array &tags);
void kdtreequeryresultstags(kdtree *kdt, ZVector *tags) {
   kdtreetsqueryresultstags(kdt, &kdt->innerbuf, tags);
}

// Distances from last query associated with kdtreerequestbuffer object.
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - KDTreeTsqueryresultsdistances().
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     R       -   possibly preallocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     R       -   filled with distances (in corresponding norm)
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreetsqueryresultsdistances(const kdtree &kdt, const kdtreerequestbuffer &buf, real_1d_array &r);
void kdtreetsqueryresultsdistances(kdtree *kdt, kdtreerequestbuffer *buf, RVector *r) {
   ae_int_t i;
   ae_int_t k;
   if (buf->kcur == 0) {
      return;
   }
   vectorsetlengthatleast(r, buf->kcur);
   k = buf->kcur;
// unload norms
//
// Abs() call is used to handle cases with negative norms
// (generated during KFN requests)
   if (kdt->normtype == 0) {
      for (i = 0; i < k; i++) {
         r->xR[i] = fabs(buf->r.xR[i]);
      }
   }
   if (kdt->normtype == 1) {
      for (i = 0; i < k; i++) {
         r->xR[i] = fabs(buf->r.xR[i]);
      }
   }
   if (kdt->normtype == 2) {
      for (i = 0; i < k; i++) {
         r->xR[i] = sqrt(fabs(buf->r.xR[i]));
      }
   }
}

// Distances from last query
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultsdistances().
//
// Inputs:
//     KDT     -   KD-tree
//     R       -   possibly preallocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     R       -   filled with distances (in corresponding norm)
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultsdistances(const kdtree &kdt, real_1d_array &r);
void kdtreequeryresultsdistances(kdtree *kdt, RVector *r) {
   kdtreetsqueryresultsdistances(kdt, &kdt->innerbuf, r);
}

// X-values from last query; 'interactive' variant for languages like  Python
// which   support    constructs   like  "X = KDTreeQueryResultsXI(KDT)"  and
// interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultsxi(const kdtree &kdt, real_2d_array &x);
void kdtreequeryresultsxi(kdtree *kdt, RMatrix *x) {
   SetMatrix(x);
   kdtreequeryresultsx(kdt, x);
}

// XY-values from last query; 'interactive' variant for languages like Python
// which   support    constructs   like "XY = KDTreeQueryResultsXYI(KDT)" and
// interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultsxyi(const kdtree &kdt, real_2d_array &xy);
void kdtreequeryresultsxyi(kdtree *kdt, RMatrix *xy) {
   SetMatrix(xy);
   kdtreequeryresultsxy(kdt, xy);
}

// Tags  from  last  query;  'interactive' variant for languages like  Python
// which  support  constructs  like "Tags = KDTreeQueryResultsTagsI(KDT)" and
// interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultstagsi(const kdtree &kdt, integer_1d_array &tags);
void kdtreequeryresultstagsi(kdtree *kdt, ZVector *tags) {
   SetVector(tags);
   kdtreequeryresultstags(kdt, tags);
}

// Distances from last query; 'interactive' variant for languages like Python
// which  support  constructs   like  "R = KDTreeQueryResultsDistancesI(KDT)"
// and interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
// API: void kdtreequeryresultsdistancesi(const kdtree &kdt, real_1d_array &r);
void kdtreequeryresultsdistancesi(kdtree *kdt, RVector *r) {
   SetVector(r);
   kdtreequeryresultsdistances(kdt, r);
}

// It is informational function which returns bounding box for entire dataset.
// This function is not visible to ALGLIB users, only ALGLIB itself  may  use
// it.
//
// This function assumes that output buffers are preallocated by caller.
// ALGLIB: Copyright 20.06.2016 by Sergey Bochkanov
void kdtreeexplorebox(kdtree *kdt, RVector *boxmin, RVector *boxmax) {
   ae_int_t i;
   vectorsetlengthatleast(boxmin, kdt->nx);
   vectorsetlengthatleast(boxmax, kdt->nx);
   for (i = 0; i < kdt->nx; i++) {
      boxmin->xR[i] = kdt->boxmin.xR[i];
      boxmax->xR[i] = kdt->boxmax.xR[i];
   }
}

// It is informational function which allows to get  information  about  node
// type. Node index is given by integer value, with 0  corresponding  to root
// node and other node indexes obtained via exploration.
//
// You should not expect that serialization/unserialization will retain  node
// indexes. You should keep in  mind  that  future  versions  of  ALGLIB  may
// introduce new node types.
//
// Output Value:
//     NodeType    -   node type:
//                     * 0 corresponds to leaf node, which can be explored by
//                       kdtreeexploreleaf() function
//                     * 1 corresponds to split node, which can  be  explored
//                       by kdtreeexploresplit() function
// ALGLIB: Copyright 20.06.2016 by Sergey Bochkanov
void kdtreeexplorenodetype(kdtree *kdt, ae_int_t node, ae_int_t *nodetype) {
   *nodetype = 0;
   ae_assert(node >= 0, "KDTreeExploreNodeType: incorrect node");
   ae_assert(node < kdt->nodes.cnt, "KDTreeExploreNodeType: incorrect node");
   if (kdt->nodes.xZ[node] > 0) {
   // Leaf node
      *nodetype = 0;
      return;
   }
   if (kdt->nodes.xZ[node] == 0) {
   // Split node
      *nodetype = 1;
      return;
   }
   ae_assert(false, "KDTreeExploreNodeType: integrity check failure");
}

// It is informational function which allows to get  information  about  leaf
// node. Node index is given by integer value, with 0  corresponding  to root
// node and other node indexes obtained via exploration.
//
// You should not expect that serialization/unserialization will retain  node
// indexes. You should keep in  mind  that  future  versions  of  ALGLIB  may
// introduce new node types.
//
// Output Values:
//     XT      -   output buffer is reallocated (if too small) and filled by
//                 XY values
//     K       -   number of rows in XY
// ALGLIB: Copyright 20.06.2016 by Sergey Bochkanov
void kdtreeexploreleaf(kdtree *kdt, ae_int_t node, RMatrix *xy, ae_int_t *k) {
   ae_int_t offs;
   ae_int_t i;
   ae_int_t j;
   *k = 0;
   ae_assert(node >= 0, "KDTreeExploreLeaf: incorrect node index");
   ae_assert(node + 1 < kdt->nodes.cnt, "KDTreeExploreLeaf: incorrect node index");
   ae_assert(kdt->nodes.xZ[node] > 0, "KDTreeExploreLeaf: incorrect node index");
   *k = kdt->nodes.xZ[node];
   offs = kdt->nodes.xZ[node + 1];
   ae_assert(offs >= 0, "KDTreeExploreLeaf: integrity error");
   ae_assert(offs + *k - 1 < kdt->xy.rows, "KDTreeExploreLeaf: integrity error");
   matrixsetlengthatleast(xy, *k, kdt->nx + kdt->ny);
   for (i = 0; i < *k; i++) {
      for (j = 0; j < kdt->nx + kdt->ny; j++) {
         xy->xyR[i][j] = kdt->xy.xyR[offs + i][kdt->nx + j];
      }
   }
}

// It is informational function which allows to get  information  about split
// node. Node index is given by integer value, with 0  corresponding  to root
// node and other node indexes obtained via exploration.
//
// You should not expect that serialization/unserialization will retain  node
// indexes. You should keep in  mind  that  future  versions  of  ALGLIB  may
// introduce new node types.
//
// Output Values:
//	*d = Nodes[idx+1]:		the dimension to split
//	*s = Nodes[idx+2]:		the offset of splitting point in Splits[]
//	*nodele = Nodes[idx+3]:		the position of left child in Nodes[]
//	*nodege = Nodes[idx+4]:		the position of right child in Nodes[]
// ALGLIB: Copyright 20.06.2016 by Sergey Bochkanov
void kdtreeexploresplit(kdtree *kdt, ae_int_t node, ae_int_t *d, double *s, ae_int_t *nodele, ae_int_t *nodege) {
   *d = 0;
   *s = 0;
   *nodele = 0;
   *nodege = 0;
   ae_assert(node >= 0, "KDTreeExploreSplit: incorrect node index");
   ae_assert(node + 4 < kdt->nodes.cnt, "KDTreeExploreSplit: incorrect node index");
   ae_assert(kdt->nodes.xZ[node] == 0, "KDTreeExploreSplit: incorrect node index");
   *d = kdt->nodes.xZ[node + 1];
   *s = kdt->splits.xR[kdt->nodes.xZ[node + 2]];
   *nodele = kdt->nodes.xZ[node + 3];
   *nodege = kdt->nodes.xZ[node + 4];
   ae_assert(*d >= 0, "KDTreeExploreSplit: integrity failure");
   ae_assert(*d < kdt->nx, "KDTreeExploreSplit: integrity failure");
   ae_assert(isfinite(*s), "KDTreeExploreSplit: integrity failure");
   ae_assert(*nodele >= 0, "KDTreeExploreSplit: integrity failure");
   ae_assert(*nodele < kdt->nodes.cnt, "KDTreeExploreSplit: integrity failure");
   ae_assert(*nodege >= 0, "KDTreeExploreSplit: integrity failure");
   ae_assert(*nodege < kdt->nodes.cnt, "KDTreeExploreSplit: integrity failure");
}

// Serializer: allocation
// ALGLIB: Copyright 14.03.2011 by Sergey Bochkanov
void kdtreealloc(ae_serializer *s, kdtree *tree) {
// Header
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
// Data
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
   allocrealmatrix(s, &tree->xy, -1, -1);
   allocintegerarray(s, &tree->tags, -1);
   allocrealarray(s, &tree->boxmin, -1);
   allocrealarray(s, &tree->boxmax, -1);
   allocintegerarray(s, &tree->nodes, -1);
   allocrealarray(s, &tree->splits, -1);
}

// Serializer: serialization
// These functions serialize a data structure to a C++ string or stream.
// * serialization can be freely moved across 32-bit and 64-bit systems,
//   and different byte orders. For example, you can serialize a string
//   on a SPARC and unserialize it on an x86.
// * ALGLIB++ serialization is compatible with serialization in ALGLIB,
//   in both directions.
// Important properties of s_out:
// * it contains alphanumeric characters, dots, underscores, minus signs
// * these symbols are grouped into words, which are separated by spaces
//   and Windows-style (CR+LF) newlines
// ALGLIB: Copyright 14.03.2011 by Sergey Bochkanov
// API: void kdtreeserialize(kdtree &obj, std::string &s_out);
// API: void kdtreeserialize(kdtree &obj, std::ostream &s_out);
void kdtreeserialize(ae_serializer *s, kdtree *tree) {
// Header
   ae_serializer_serialize_int(s, getkdtreeserializationcode());
   ae_serializer_serialize_int(s, nearestneighbor_kdtreefirstversion);
// Data
   ae_serializer_serialize_int(s, tree->n);
   ae_serializer_serialize_int(s, tree->nx);
   ae_serializer_serialize_int(s, tree->ny);
   ae_serializer_serialize_int(s, tree->normtype);
   serializerealmatrix(s, &tree->xy, -1, -1);
   serializeintegerarray(s, &tree->tags, -1);
   serializerealarray(s, &tree->boxmin, -1);
   serializerealarray(s, &tree->boxmax, -1);
   serializeintegerarray(s, &tree->nodes, -1);
   serializerealarray(s, &tree->splits, -1);
}

// Serializer: unserialization
// These functions unserialize a data structure from a C++ string or stream.
// Important properties of s_in:
// * any combination of spaces, tabs, Windows or Unix stype newlines can
//   be used as separators, so as to allow flexible reformatting of the
//   stream or string from text or XML files.
// * But you should not insert separators into the middle of the "words"
//   nor you should change case of letters.
// ALGLIB: Copyright 14.03.2011 by Sergey Bochkanov
// API: void kdtreeunserialize(const std::string &s_in, kdtree &obj);
// API: void kdtreeunserialize(const std::istream &s_in, kdtree &obj);
void kdtreeunserialize(ae_serializer *s, kdtree *tree) {
   ae_int_t i0;
   ae_int_t i1;
   SetObj(kdtree, tree);
// check correctness of header
   i0 = ae_serializer_unserialize_int(s);
   ae_assert(i0 == getkdtreeserializationcode(), "kdtreeunserialize: stream header corrupted");
   i1 = ae_serializer_unserialize_int(s);
   ae_assert(i1 == nearestneighbor_kdtreefirstversion, "kdtreeunserialize: stream header corrupted");
// Unserialize data
   tree->n = ae_serializer_unserialize_int(s);
   tree->nx = ae_serializer_unserialize_int(s);
   tree->ny = ae_serializer_unserialize_int(s);
   tree->normtype = ae_serializer_unserialize_int(s);
   unserializerealmatrix(s, &tree->xy);
   unserializeintegerarray(s, &tree->tags);
   unserializerealarray(s, &tree->boxmin);
   unserializerealarray(s, &tree->boxmax);
   unserializeintegerarray(s, &tree->nodes);
   unserializerealarray(s, &tree->splits);
   kdtreecreaterequestbuffer(tree, &tree->innerbuf);
}

void kdtreerequestbuffer_init(void *_p, bool make_automatic) {
   kdtreerequestbuffer *p = (kdtreerequestbuffer *)_p;
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->boxmin, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->boxmax, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->idx, 0, DT_INT, make_automatic);
   ae_vector_init(&p->r, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->buf, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->curboxmin, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->curboxmax, 0, DT_REAL, make_automatic);
}

void kdtreerequestbuffer_copy(void *_dst, void *_src, bool make_automatic) {
   kdtreerequestbuffer *dst = (kdtreerequestbuffer *)_dst;
   kdtreerequestbuffer *src = (kdtreerequestbuffer *)_src;
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   ae_vector_copy(&dst->boxmin, &src->boxmin, make_automatic);
   ae_vector_copy(&dst->boxmax, &src->boxmax, make_automatic);
   dst->kneeded = src->kneeded;
   dst->rneeded = src->rneeded;
   dst->selfmatch = src->selfmatch;
   dst->approxf = src->approxf;
   dst->kcur = src->kcur;
   ae_vector_copy(&dst->idx, &src->idx, make_automatic);
   ae_vector_copy(&dst->r, &src->r, make_automatic);
   ae_vector_copy(&dst->buf, &src->buf, make_automatic);
   ae_vector_copy(&dst->curboxmin, &src->curboxmin, make_automatic);
   ae_vector_copy(&dst->curboxmax, &src->curboxmax, make_automatic);
   dst->curdist = src->curdist;
}

void kdtreerequestbuffer_free(void *_p, bool make_automatic) {
   kdtreerequestbuffer *p = (kdtreerequestbuffer *)_p;
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->boxmin, make_automatic);
   ae_vector_free(&p->boxmax, make_automatic);
   ae_vector_free(&p->idx, make_automatic);
   ae_vector_free(&p->r, make_automatic);
   ae_vector_free(&p->buf, make_automatic);
   ae_vector_free(&p->curboxmin, make_automatic);
   ae_vector_free(&p->curboxmax, make_automatic);
}

void kdtree_init(void *_p, bool make_automatic) {
   kdtree *p = (kdtree *)_p;
   ae_matrix_init(&p->xy, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tags, 0, DT_INT, make_automatic);
   ae_vector_init(&p->boxmin, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->boxmax, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->nodes, 0, DT_INT, make_automatic);
   ae_vector_init(&p->splits, 0, DT_REAL, make_automatic);
   kdtreerequestbuffer_init(&p->innerbuf, make_automatic);
}

void kdtree_copy(void *_dst, void *_src, bool make_automatic) {
   kdtree *dst = (kdtree *)_dst;
   kdtree *src = (kdtree *)_src;
   dst->n = src->n;
   dst->nx = src->nx;
   dst->ny = src->ny;
   dst->normtype = src->normtype;
   ae_matrix_copy(&dst->xy, &src->xy, make_automatic);
   ae_vector_copy(&dst->tags, &src->tags, make_automatic);
   ae_vector_copy(&dst->boxmin, &src->boxmin, make_automatic);
   ae_vector_copy(&dst->boxmax, &src->boxmax, make_automatic);
   ae_vector_copy(&dst->nodes, &src->nodes, make_automatic);
   ae_vector_copy(&dst->splits, &src->splits, make_automatic);
   kdtreerequestbuffer_copy(&dst->innerbuf, &src->innerbuf, make_automatic);
   dst->debugcounter = src->debugcounter;
}

void kdtree_free(void *_p, bool make_automatic) {
   kdtree *p = (kdtree *)_p;
   ae_matrix_free(&p->xy, make_automatic);
   ae_vector_free(&p->tags, make_automatic);
   ae_vector_free(&p->boxmin, make_automatic);
   ae_vector_free(&p->boxmax, make_automatic);
   ae_vector_free(&p->nodes, make_automatic);
   ae_vector_free(&p->splits, make_automatic);
   kdtreerequestbuffer_free(&p->innerbuf, make_automatic);
}
} // end of namespace alglib_impl

namespace alglib {
// Buffer object which is used to perform nearest neighbor  requests  in  the
// multithreaded mode (multiple threads working with same KD-tree object).
//
// This object should be created with KDTreeCreateRequestBuffer().
DefClass(kdtreerequestbuffer, )

// KD-tree object.
DefClass(kdtree, )

void kdtreeserialize(kdtree &obj, std::string &s_out) {
   alglib_impl::ae_state_init();
   TryCatch()
   NewSerializer(serializer);
   alglib_impl::ae_serializer_alloc_start(&serializer);
   alglib_impl::kdtreealloc(&serializer, obj.c_ptr());
   ae_int_t ssize = alglib_impl::ae_serializer_get_alloc_size(&serializer);
   s_out.clear();
   s_out.reserve((size_t)(ssize + 1));
   alglib_impl::ae_serializer_sstart_str(&serializer, &s_out);
   alglib_impl::kdtreeserialize(&serializer, obj.c_ptr());
   alglib_impl::ae_serializer_stop(&serializer);
   alglib_impl::ae_assert(s_out.length() <= (size_t)ssize, "kdtreeserialize: serialization integrity error");
   alglib_impl::ae_state_clear();
}
void kdtreeserialize(kdtree &obj, std::ostream &s_out) {
   alglib_impl::ae_state_init();
   TryCatch()
   NewSerializer(serializer);
   alglib_impl::ae_serializer_alloc_start(&serializer);
   alglib_impl::kdtreealloc(&serializer, obj.c_ptr());
   alglib_impl::ae_serializer_get_alloc_size(&serializer); // not actually needed, but we have to ask
   alglib_impl::ae_serializer_sstart_stream(&serializer, &s_out);
   alglib_impl::kdtreeserialize(&serializer, obj.c_ptr());
   alglib_impl::ae_serializer_stop(&serializer);
   alglib_impl::ae_state_clear();
}

void kdtreeunserialize(const std::string &s_in, kdtree &obj) {
   alglib_impl::ae_state_init();
   TryCatch()
   NewSerializer(serializer);
   alglib_impl::ae_serializer_ustart_str(&serializer, &s_in);
   alglib_impl::kdtreeunserialize(&serializer, obj.c_ptr());
   alglib_impl::ae_serializer_stop(&serializer);
   alglib_impl::ae_state_clear();
}
void kdtreeunserialize(const std::istream &s_in, kdtree &obj) {
   alglib_impl::ae_state_init();
   TryCatch()
   NewSerializer(serializer);
   alglib_impl::ae_serializer_ustart_stream(&serializer, &s_in);
   alglib_impl::kdtreeunserialize(&serializer, obj.c_ptr());
   alglib_impl::ae_serializer_stop(&serializer);
   alglib_impl::ae_state_clear();
}

void kdtreecreaterequestbuffer(const kdtree &kdt, kdtreerequestbuffer &buf) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreecreaterequestbuffer(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf));
   alglib_impl::ae_state_clear();
}

void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreebuildtagged(ConstT(ae_matrix, xy), ConstT(ae_vector, tags), n, nx, ny, normtype, ConstT(kdtree, kdt));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt) {
   if (xy.rows() != tags.length()) ThrowError("Error while calling 'kdtreebuildtagged': looks like one of arguments has wrong size");
   ae_int_t n = xy.rows();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreebuildtagged(ConstT(ae_matrix, xy), ConstT(ae_vector, tags), n, nx, ny, normtype, ConstT(kdtree, kdt));
   alglib_impl::ae_state_clear();
}
#endif

void kdtreebuild(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreebuild(ConstT(ae_matrix, xy), n, nx, ny, normtype, ConstT(kdtree, kdt));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void kdtreebuild(const real_2d_array &xy, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt) {
   ae_int_t n = xy.rows();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreebuild(ConstT(ae_matrix, xy), n, nx, ny, normtype, ConstT(kdtree, kdt));
   alglib_impl::ae_state_clear();
}
#endif

ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryaknn(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), k, selfmatch, eps);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const double eps) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryaknn(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), k, selfmatch, eps);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryaknn(ConstT(kdtree, kdt), ConstT(ae_vector, x), k, selfmatch, eps);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const double eps) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryaknn(ConstT(kdtree, kdt), ConstT(ae_vector, x), k, selfmatch, eps);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryknn(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), k, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryknn(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), k, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryknn(ConstT(kdtree, kdt), ConstT(ae_vector, x), k, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryknn(ConstT(kdtree, kdt), ConstT(ae_vector, x), k, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryrnn(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryrnn(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryrnn(ConstT(kdtree, kdt), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryrnn(ConstT(kdtree, kdt), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryrnnu(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsqueryrnnu(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryrnnu(ConstT(kdtree, kdt), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#if !defined AE_NO_EXCEPTIONS
ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r) {
   bool selfmatch = true;
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequeryrnnu(ConstT(kdtree, kdt), ConstT(ae_vector, x), r, selfmatch);
   alglib_impl::ae_state_clear();
   return Z;
}
#endif

ae_int_t kdtreetsquerybox(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &boxmin, const real_1d_array &boxmax) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreetsquerybox(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, boxmin), ConstT(ae_vector, boxmax));
   alglib_impl::ae_state_clear();
   return Z;
}

ae_int_t kdtreequerybox(const kdtree &kdt, const real_1d_array &boxmin, const real_1d_array &boxmax) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::kdtreequerybox(ConstT(kdtree, kdt), ConstT(ae_vector, boxmin), ConstT(ae_vector, boxmax));
   alglib_impl::ae_state_clear();
   return Z;
}

void kdtreetsqueryresultsx(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreetsqueryresultsx(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultsx(const kdtree &kdt, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultsx(ConstT(kdtree, kdt), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void kdtreetsqueryresultsxy(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &xy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreetsqueryresultsxy(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_matrix, xy));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultsxy(const kdtree &kdt, real_2d_array &xy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultsxy(ConstT(kdtree, kdt), ConstT(ae_matrix, xy));
   alglib_impl::ae_state_clear();
}

void kdtreetsqueryresultstags(const kdtree &kdt, const kdtreerequestbuffer &buf, integer_1d_array &tags) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreetsqueryresultstags(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, tags));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultstags(const kdtree &kdt, integer_1d_array &tags) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultstags(ConstT(kdtree, kdt), ConstT(ae_vector, tags));
   alglib_impl::ae_state_clear();
}

void kdtreetsqueryresultsdistances(const kdtree &kdt, const kdtreerequestbuffer &buf, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreetsqueryresultsdistances(ConstT(kdtree, kdt), ConstT(kdtreerequestbuffer, buf), ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultsdistances(const kdtree &kdt, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultsdistances(ConstT(kdtree, kdt), ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultsxi(const kdtree &kdt, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultsxi(ConstT(kdtree, kdt), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultsxyi(const kdtree &kdt, real_2d_array &xy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultsxyi(ConstT(kdtree, kdt), ConstT(ae_matrix, xy));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultstagsi(const kdtree &kdt, integer_1d_array &tags) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultstagsi(ConstT(kdtree, kdt), ConstT(ae_vector, tags));
   alglib_impl::ae_state_clear();
}

void kdtreequeryresultsdistancesi(const kdtree &kdt, real_1d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::kdtreequeryresultsdistancesi(ConstT(kdtree, kdt), ConstT(ae_vector, r));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib
