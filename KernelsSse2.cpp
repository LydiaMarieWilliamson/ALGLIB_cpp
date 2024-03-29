#define InAlgLib
// Must be defined before we include kernel header.
#define _ALGLIB_IMPL_DEFINES
#define _ALGLIB_INTEGRITY_CHECKS_ONCE
#include "KernelsSse2.h"

namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_SSE2_INTRINSICS
double sse2_rdotv(const ae_int_t n, const Real *x, const Real *y) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const ae_int_t unrollLen = (sse2len >> 3) << 3;
   const __m128d *__restrict pX = (const __m128d *)x;
   const __m128d *__restrict pY = (const __m128d *)y;
   __m128d ans;
   if (unrollLen >= 8) {
      __m128d unroll0 = _mm_mul_pd(pX[0], pY[0]);
      __m128d unroll1 = _mm_mul_pd(pX[1], pY[1]);
      __m128d unroll2 = _mm_mul_pd(pX[2], pY[2]);
      __m128d unroll3 = _mm_mul_pd(pX[3], pY[3]);
      __m128d unroll4 = _mm_mul_pd(pX[4], pY[4]);
      __m128d unroll5 = _mm_mul_pd(pX[5], pY[5]);
      __m128d unroll6 = _mm_mul_pd(pX[6], pY[6]);
      __m128d unroll7 = _mm_mul_pd(pX[7], pY[7]);
      for (i = 8; i < unrollLen; i += 8) {
         unroll0 = _mm_add_pd(_mm_mul_pd(pX[i], pY[i]), unroll0);
         unroll1 = _mm_add_pd(_mm_mul_pd(pX[i + 1], pY[i + 1]), unroll1);
         unroll2 = _mm_add_pd(_mm_mul_pd(pX[i + 2], pY[i + 2]), unroll2);
         unroll3 = _mm_add_pd(_mm_mul_pd(pX[i + 3], pY[i + 3]), unroll3);
         unroll4 = _mm_add_pd(_mm_mul_pd(pX[i + 4], pY[i + 4]), unroll4);
         unroll5 = _mm_add_pd(_mm_mul_pd(pX[i + 5], pY[i + 5]), unroll5);
         unroll6 = _mm_add_pd(_mm_mul_pd(pX[i + 6], pY[i + 6]), unroll6);
         unroll7 = _mm_add_pd(_mm_mul_pd(pX[i + 7], pY[i + 7]), unroll7);
      }
      switch (sse2len - unrollLen) {
         case 7:
            unroll6 = _mm_add_pd(_mm_mul_pd(pX[i + 6], pY[i + 6]), unroll6);
         case 6:
            unroll5 = _mm_add_pd(_mm_mul_pd(pX[i + 5], pY[i + 5]), unroll5);
         case 5:
            unroll4 = _mm_add_pd(_mm_mul_pd(pX[i + 4], pY[i + 4]), unroll4);
         case 4:
            unroll3 = _mm_add_pd(_mm_mul_pd(pX[i + 3], pY[i + 3]), unroll3);
         case 3:
            unroll2 = _mm_add_pd(_mm_mul_pd(pX[i + 2], pY[i + 2]), unroll2);
         case 2:
            unroll1 = _mm_add_pd(_mm_mul_pd(pX[i + 1], pY[i + 1]), unroll1);
         case 1:
            unroll0 = _mm_add_pd(_mm_mul_pd(pX[i], pY[i]), unroll0);
      }
      ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(unroll0, unroll1), _mm_add_pd(unroll2, unroll3)), _mm_add_pd(_mm_add_pd(unroll4, unroll5), _mm_add_pd(unroll6, unroll7)));
   } else {
      switch (sse2len) {
         case 0:
            if (n == 0) {
               return 0;
            } else {
               return x[0] * y[0];
            }
         case 1:
            ans = _mm_mul_pd(pX[0], pY[0]);
         break;
         case 2:
            ans = _mm_add_pd(_mm_mul_pd(pX[0], pY[0]), _mm_mul_pd(pX[1], pY[1]));
         break;
         case 3:
            ans = _mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pY[0]), _mm_mul_pd(pX[1], pY[1])), _mm_mul_pd(pX[2], pY[2]));
         break;
         case 4:
            ans = _mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pY[0]), _mm_mul_pd(pX[1], pY[1])), _mm_add_pd(_mm_mul_pd(pX[2], pY[2]), _mm_mul_pd(pX[3], pY[3])));
         break;
         case 5:
            ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pY[0]), _mm_mul_pd(pX[1], pY[1])), _mm_add_pd(_mm_mul_pd(pX[2], pY[2]), _mm_mul_pd(pX[3], pY[3]))), _mm_mul_pd(pX[4], pY[4]));
         break;
         case 6:
            ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pY[0]), _mm_mul_pd(pX[1], pY[1])), _mm_add_pd(_mm_mul_pd(pX[2], pY[2]), _mm_mul_pd(pX[3], pY[3]))), _mm_add_pd(_mm_mul_pd(pX[4], pY[4]), _mm_mul_pd(pX[5], pY[5])));
         break;
         case 7:
            ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pY[0]), _mm_mul_pd(pX[1], pY[1])), _mm_add_pd(_mm_mul_pd(pX[2], pY[2]), _mm_mul_pd(pX[3], pY[3]))), _mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[4], pY[4]), _mm_mul_pd(pX[5], pY[5])), _mm_mul_pd(pX[6], pY[6])));
         break;
      }
   }
   const double *pComps = (const double *)&ans;
   double scalar = pComps[0] + pComps[1];
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      return scalar + x[tail] * y[tail];
   } else {
      return scalar;
   }
}

double sse2_rdotv2(const ae_int_t n, const Real *x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const ae_int_t unrollLen = (sse2len >> 3) << 3;
   const __m128d *__restrict pX = (const __m128d *)x;
   __m128d ans;
   if (unrollLen >= 8) {
      __m128d unroll0 = _mm_mul_pd(pX[0], pX[0]);
      __m128d unroll1 = _mm_mul_pd(pX[1], pX[1]);
      __m128d unroll2 = _mm_mul_pd(pX[2], pX[2]);
      __m128d unroll3 = _mm_mul_pd(pX[3], pX[3]);
      __m128d unroll4 = _mm_mul_pd(pX[4], pX[4]);
      __m128d unroll5 = _mm_mul_pd(pX[5], pX[5]);
      __m128d unroll6 = _mm_mul_pd(pX[6], pX[6]);
      __m128d unroll7 = _mm_mul_pd(pX[7], pX[7]);
      for (i = 8; i < unrollLen; i += 8) {
         unroll0 = _mm_add_pd(_mm_mul_pd(pX[i], pX[i]), unroll0);
         unroll1 = _mm_add_pd(_mm_mul_pd(pX[i + 1], pX[i + 1]), unroll1);
         unroll2 = _mm_add_pd(_mm_mul_pd(pX[i + 2], pX[i + 2]), unroll2);
         unroll3 = _mm_add_pd(_mm_mul_pd(pX[i + 3], pX[i + 3]), unroll3);
         unroll4 = _mm_add_pd(_mm_mul_pd(pX[i + 4], pX[i + 4]), unroll4);
         unroll5 = _mm_add_pd(_mm_mul_pd(pX[i + 5], pX[i + 5]), unroll5);
         unroll6 = _mm_add_pd(_mm_mul_pd(pX[i + 6], pX[i + 6]), unroll6);
         unroll7 = _mm_add_pd(_mm_mul_pd(pX[i + 7], pX[i + 7]), unroll7);
      }
      switch (sse2len - unrollLen) {
         case 7:
            unroll6 = _mm_add_pd(_mm_mul_pd(pX[i + 6], pX[i + 6]), unroll6);
         case 6:
            unroll5 = _mm_add_pd(_mm_mul_pd(pX[i + 5], pX[i + 5]), unroll5);
         case 5:
            unroll4 = _mm_add_pd(_mm_mul_pd(pX[i + 4], pX[i + 4]), unroll4);
         case 4:
            unroll3 = _mm_add_pd(_mm_mul_pd(pX[i + 3], pX[i + 3]), unroll3);
         case 3:
            unroll2 = _mm_add_pd(_mm_mul_pd(pX[i + 2], pX[i + 2]), unroll2);
         case 2:
            unroll1 = _mm_add_pd(_mm_mul_pd(pX[i + 1], pX[i + 1]), unroll1);
         case 1:
            unroll0 = _mm_add_pd(_mm_mul_pd(pX[i], pX[i]), unroll0);
      }
      ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(unroll0, unroll1), _mm_add_pd(unroll2, unroll3)), _mm_add_pd(_mm_add_pd(unroll4, unroll5), _mm_add_pd(unroll6, unroll7)));
   } else {
      switch (sse2len) {
         case 0:
            if (n == 0) {
               return 0;
            } else {
               return x[0] * x[0];
            }
         case 1:
            ans = _mm_mul_pd(pX[0], pX[0]);
         break;
         case 2:
            ans = _mm_add_pd(_mm_mul_pd(pX[0], pX[0]), _mm_mul_pd(pX[1], pX[1]));
         break;
         case 3:
            ans = _mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pX[0]), _mm_mul_pd(pX[1], pX[1])), _mm_mul_pd(pX[2], pX[2]));
         break;
         case 4:
            ans = _mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pX[0]), _mm_mul_pd(pX[1], pX[1])), _mm_add_pd(_mm_mul_pd(pX[2], pX[2]), _mm_mul_pd(pX[3], pX[3])));
         break;
         case 5:
            ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pX[0]), _mm_mul_pd(pX[1], pX[1])), _mm_add_pd(_mm_mul_pd(pX[2], pX[2]), _mm_mul_pd(pX[3], pX[3]))), _mm_mul_pd(pX[4], pX[4]));
         break;
         case 6:
            ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pX[0]), _mm_mul_pd(pX[1], pX[1])), _mm_add_pd(_mm_mul_pd(pX[2], pX[2]), _mm_mul_pd(pX[3], pX[3]))), _mm_add_pd(_mm_mul_pd(pX[4], pX[4]), _mm_mul_pd(pX[5], pX[5])));
         break;
         case 7:
            ans = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[0], pX[0]), _mm_mul_pd(pX[1], pX[1])), _mm_add_pd(_mm_mul_pd(pX[2], pX[2]), _mm_mul_pd(pX[3], pX[3]))), _mm_add_pd(_mm_add_pd(_mm_mul_pd(pX[4], pX[4]), _mm_mul_pd(pX[5], pX[5])), _mm_mul_pd(pX[6], pX[6])));
         break;
      }
   }
   const double *pComps = (const double *)&ans;
   double scalar = pComps[0] + pComps[1];
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      return scalar + x[tail] * x[tail];
   } else {
      return scalar;
   }
}

void sse2_rcopyv(const ae_int_t n, const Real *__restrict x, Real *__restrict y) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const ae_int_t tail = sse2len << 1;
   const __m128d *__restrict pSrc = (const __m128d *)x;
   __m128d *__restrict pDest = (__m128d *)y;
   for (i = 0; i < sse2len; i++)
      pDest[i] = pSrc[i];
   if (n - tail)
      *(double *)(pDest + i) = *(const double *)(pSrc + i);
}

void sse2_rcopymulv(const ae_int_t n, const double v, const Real *__restrict x, Real *__restrict y) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)x;
   __m128d *__restrict pDest = (__m128d *)y;
   const __m128d sse2v = _mm_set1_pd(v);
   const ae_int_t tail = sse2len << 1;
   for (i = 0; i < sse2len; i++) {
      pDest[i] = _mm_mul_pd(sse2v, pSrc[i]);
   }
   if (n - tail) {
      *(double *)(pDest + i) = v * *(const double *)(pSrc + i);
   }
}

void sse2_icopyv(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y) {
   const ae_int_t tail = (n * sizeof(ae_int_t)) & 15;
   const ae_int_t even = (n * sizeof(ae_int_t)) - tail;
   __m128i *__restrict pDest = (__m128i *)y;
   const __m128i *__restrict pSrc = (const __m128i *)x;
   const ae_int_t nVec = even >> 4;
   ae_int_t i;
   for (i = 0; i < nVec; i++) {
      pDest[i] = pSrc[i];
   }
   i = even / sizeof(ae_int_t);
   if (tail & 8) {
      *(ae_int64_t *)(y + i) = *(const ae_int64_t *)(x + i);
      i += 8 / sizeof(ae_int_t);
   }
   if (tail & 4) {
      *(ae_int32_t *)(y + i) = *(const ae_int32_t *)(x + i);
   }
}

void sse2_bcopyv(const ae_int_t n, const bool *__restrict x, bool *__restrict y) {
   const ae_int_t tail = n & 15;
   const ae_int_t even = n - tail;
   __m128i *__restrict pDest = (__m128i *)y;
   const __m128i *__restrict pSrc = (const __m128i *)x;
   const ae_int_t nVec = even >> 4;
   ae_int_t i;
   for (i = 0; i < nVec; i++) {
      pDest[i] = pSrc[i];
   }
   i = even;
   if (tail & 8) {
      *(ae_int64_t *)(y + i) = *(const ae_int64_t *)(x + i);
      i += 8;
   }
   if (tail & 4) {
      *(ae_int32_t *)(y + i) = *(const ae_int32_t *)(x + i);
      i += 4;
   }
   if (tail & 2) {
      *(y + i) = *(x + i);
      *(y + i + 1) = *(x + i + 1);
      i += 2;
   }
   if (tail & 1) {
      *(y + i) = *(x + i);
   }
}

void sse2_rsetv(const ae_int_t n, const double v, Real *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   __m128d *__restrict pDest = (__m128d *)x;
   const __m128d sse2v = _mm_set1_pd(v);
   for (i = 0; i < sse2len; i++) {
      pDest[i] = sse2v;
   }
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      *(double *)(pDest + i) = v;
   }
}

void sse2_rsetvx(const ae_int_t n, const double v, Real *__restrict x) {
   if (n <= 4) {
      ae_int_t j;
      for (j = 0; j < n; j++)
         x[j] = v;
      return;
   }
   if (((ptrdiff_t)x & 15) == 0) {
      sse2_rsetv(n, v, x);
      return;
   }
   x[0] = v;
   sse2_rsetv(n - 1, v, x + 1);
}

void sse2_isetv(const ae_int_t n, const ae_int_t v, ae_int_t *__restrict x) {
   const ae_int_t tail = (n * sizeof(ae_int_t)) & 15;
   const ae_int_t even = (n * sizeof(ae_int_t)) - tail;
   __m128i *__restrict pDest = (__m128i *)x;
   const ae_int_t v2[2] = { v, v };
   const __m128i sse2v = ((sizeof(v) == 4) ? _mm_set1_epi32((ae_int32_t)v) : _mm_loadu_si128((const __m128i *)v2));
   const ae_int_t nVec = even >> 4;
   ae_int_t i;
   for (i = 0; i < nVec; i++) {
      pDest[i] = sse2v;
   }
   memmove(pDest + i, &sse2v, tail);
}

void sse2_bsetv(const ae_int_t n, const bool v, bool *__restrict x) {
   const ae_int_t tail = n & 15;
   const ae_int_t even = n - tail;
   __m128i *__restrict pDest = (__m128i *)x;
   const __m128i sse2v = _mm_set1_epi8(v);
   const ae_int_t nVec = even >> 4;
   ae_int_t i;
   for (i = 0; i < nVec; i++) {
      pDest[i] = sse2v;
   }
// _mm_storel_epi64() has a too high latency and too low throughput on the recent (Skylake+) processors
   memset(x + even, v, tail);
}

void sse2_rmulv(const ae_int_t n, const double v, double *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   __m128d *__restrict pDest = (__m128d *)x;
   const __m128d sse2v = _mm_set1_pd(v);
   for (i = 0; i < sse2len; i++) {
      pDest[i] = _mm_mul_pd(sse2v, pDest[i]);
   }
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      *(double *)(pDest + i) = v * *(const double *)(pDest + i);
   }
}

void sse2_rmulvx(const ae_int_t n, const double v, double *__restrict x) {
   if (n <= 4) {
      ae_int_t i;
      for (i = 0; i < n; i++)
         x[i] *= v;
      return;
   }
   if (((ptrdiff_t)x & 15) == 0) {
      sse2_rmulv(n, v, x);
      return;
   }
   x[0] *= v;
   sse2_rmulv(n - 1, v, x + 1);
}

void sse2_raddv(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)y;
   __m128d *__restrict pDest = (__m128d *)x;
   const __m128d sse2alpha = _mm_set1_pd(alpha);
   for (i = 0; i < sse2len; i++) {
      pDest[i] = _mm_add_pd(_mm_mul_pd(sse2alpha, pSrc[i]), pDest[i]);
   }
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      *(double *)(pDest + i) = alpha * *(const double *)(pSrc + i) + *(const double *)(pDest + i);
   }
}

static void sse2_raddvx_xaligned(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x) {
   ae_int_t i;
   const ae_int_t vecLen = (n >> 1) << 1;
   const __m128d sse2alpha = _mm_set1_pd(alpha);
   __m128d *__restrict pDest = (__m128d *)x;
   for (i = 0; i < vecLen; i += 2) {
      const ae_int_t iDest = i >> 1;
      pDest[iDest] = _mm_add_pd(_mm_mul_pd(sse2alpha, _mm_loadu_pd(y + i)), pDest[iDest]);
   }
   if (n - vecLen)
      x[i] += alpha * y[i];
}

void sse2_raddvx(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x) {
   if (n <= 4) {
      ae_int_t i;
      for (i = 0; i < n; i++)
         x[i] += alpha * y[i];
      return;
   }
   if (((ptrdiff_t)x & 15) == 0) {
      sse2_raddvx_xaligned(n, alpha, y, x);
      return;
   }
   x[0] += alpha * y[0];
   sse2_raddvx_xaligned(n - 1, alpha, y + 1, x + 1);
}

void sse2_rmergemulv(const ae_int_t n, const Real *__restrict y, Real *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)y;
   __m128d *__restrict pDest = (__m128d *)x;
   for (i = 0; i < sse2len; i++) {
      pDest[i] = _mm_mul_pd(pSrc[i], pDest[i]);
   }
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      *(double *)(pDest + i) = *(const double *)(pSrc + i) * *(const double *)(pDest + i);
   }
}

void sse2_rmergemaxv(const ae_int_t n, const Real *__restrict y, Real *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)y;
   __m128d *__restrict pDest = (__m128d *)x;
   for (i = 0; i < sse2len; i++) {
      pDest[i] = _mm_max_pd(pSrc[i], pDest[i]);
   }
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      *(double *)(pDest + i) = rmax2(*(const double *)(pSrc + i), *(const double *)(pDest + i));
   }
}

void sse2_rmergeminv(const ae_int_t n, const Real *__restrict y, Real *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)y;
   __m128d *__restrict pDest = (__m128d *)x;
   for (i = 0; i < sse2len; i++) {
      pDest[i] = _mm_min_pd(pSrc[i], pDest[i]);
   }
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      *(double *)(pDest + i) = rmin2(*(const double *)(pSrc + i), *(const double *)(pDest + i));
   }
}

double sse2_rmaxv(ae_int_t n, const Real *__restrict x) {
   ae_int_t i;
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)x;
   if (n <= 4) {
      double result;
      if (n == 0)
         return 0.0;
      result = x[0];
      for (i = 1; i < n; i++) {
         double v = x[i];
         if (v > result)
            result = v;
      }
      return result;
   }
   __m128d curMax = pSrc[0];
   for (i = 1; i < sse2len; i++) {
      curMax = _mm_max_pd(curMax, pSrc[i]);
   }
   const double *pComps = (const double *)&curMax;
   const double dMax = (pComps[0] > pComps[1]) ? pComps[0] : pComps[1];
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      const double candidate = *(const double *)(pSrc + i);
      return (candidate > dMax) ? candidate : dMax;
   } else {
      return dMax;
   }
}

double sse2_rmaxabsv(ae_int_t n, const Real *__restrict x) {
   const __m128d signMask = _mm_set1_pd(-0.); // -0. = 1 << 63
   const ae_int_t sse2len = n >> 1;
   const __m128d *__restrict pSrc = (const __m128d *)x;
   if (n <= 4) {
      double result;
      ae_int_t i;
      result = 0;
      for (i = 0; i < n; i++) {
         double v = fabs(x[i]);
         if (v > result)
            result = v;
      }
      return result;
   }
   __m128d curMax = _mm_andnot_pd(signMask, pSrc[0]); // abs
   ae_int_t i;
   for (i = 1; i < sse2len; i++)
      curMax = _mm_max_pd(curMax, _mm_andnot_pd(signMask, pSrc[i])); // abs
   const double *pComps = (const double *)&curMax;
   const double dMax = (pComps[0] > pComps[1]) ? pComps[0] : pComps[1];
   const ae_int_t tail = sse2len << 1;
   if (n - tail) {
      const double candidate = fabs(*(const double *)(pSrc + i));
      return (candidate > dMax) ? candidate : dMax;
   } else {
      return dMax;
   }
}

static void sse2_rcopyvx_xaligned(const ae_int_t n, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   const ae_int_t vecLen = (n >> 1) << 1;
   const __m128d *__restrict pSrc = (const __m128d *)x;
   for (i = 0; i < vecLen; i += 2) {
      const ae_int_t iSrc = i >> 1;
      _mm_storeu_pd(y + i, pSrc[iSrc]);
   }
   if (n - vecLen) {
      y[i] = x[i];
   }
}

void sse2_rcopyvx(const ae_int_t n, const double *__restrict x, double *__restrict y) {
   if (((ptrdiff_t)x & 15) == 0) {
      sse2_rcopyvx_xaligned(n, x, y);
      return;
   }
   y[0] = x[0];
   sse2_rcopyvx_xaligned(n - 1, x + 1, y + 1);
}

static void sse2_icopyvx_xaligned(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y) {
   const ae_int_t tail = (n * sizeof(ae_int_t)) & 15;
   const ae_int_t even = (n * sizeof(ae_int_t)) - tail;
   const __m128i *__restrict pSrc = (const __m128i *)x;
   const ae_int_t nVec = even >> 4;
   const ae_int_t shift_by = 2 - sizeof(ae_int_t) / 8;
   ae_int_t i;
   for (i = 0; i < nVec; i++) {
      const ae_int_t j = i << shift_by;
      _mm_storeu_si128((__m128i *)(y + j), pSrc[i]);
   }
   i = even / sizeof(ae_int_t);
   if (tail & 8) {
      *(ae_int64_t *)(y + i) = *(const ae_int64_t *)(x + i);
      i += 8 / sizeof(ae_int_t);
   }
   if (tail & 4) {
      *(ae_int32_t *)(y + i) = *(const ae_int32_t *)(x + i);
   }
}

void sse2_icopyvx(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y) {
   const ptrdiff_t unal = (ptrdiff_t)x & 15;
   if (n <= 8) {
      ae_int_t j;
      for (j = 0; j < n; j++)
         y[j] = x[j];
      return;
   }
   if (unal == 0) {
      sse2_icopyvx_xaligned(n, x, y);
      return;
   }
   const ae_int_t offset = 16 - unal;
   memmove(y, x, offset);
   const ae_int_t nDone = offset / sizeof(ae_int_t);
   sse2_icopyvx_xaligned(n - nDone, x + nDone, y + nDone);
}
#endif // ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_SSE2_INTRINSICS
} // end of namespace alglib_impl
