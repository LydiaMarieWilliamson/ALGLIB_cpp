#define InAlgLib
// Must be defined before we include kernel header.
#define _ALGLIB_IMPL_DEFINES
#define _ALGLIB_INTEGRITY_CHECKS_ONCE
#include "KernelsFma.h"

namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_FMA_INTRINSICS
double fma_rdotv(const ae_int_t n, const Real *__restrict x, const Real *__restrict y) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const ae_int_t fmaLen = (avx2len >> 2) << 2;
   const __m256d *__restrict pX = (const __m256d *)x;
   const __m256d *__restrict pY = (const __m256d *)y;
   __m256d ans;
   if (fmaLen >= 4) {
      __m256d fmaUnroll[4];
      fmaUnroll[0] = _mm256_mul_pd(pX[0], pY[0]);
      fmaUnroll[1] = _mm256_mul_pd(pX[1], pY[1]);
      fmaUnroll[2] = _mm256_mul_pd(pX[2], pY[2]);
      fmaUnroll[3] = _mm256_mul_pd(pX[3], pY[3]);
      for (i = 4; i < fmaLen; i += 4) {
         fmaUnroll[0] = _mm256_fmadd_pd(pX[i], pY[i], fmaUnroll[0]);
         fmaUnroll[1] = _mm256_fmadd_pd(pX[i + 1], pY[i + 1], fmaUnroll[1]);
         fmaUnroll[2] = _mm256_fmadd_pd(pX[i + 2], pY[i + 2], fmaUnroll[2]);
         fmaUnroll[3] = _mm256_fmadd_pd(pX[i + 3], pY[i + 3], fmaUnroll[3]);
      }
      switch (avx2len - fmaLen) {
         case 3:
            fmaUnroll[2] = _mm256_fmadd_pd(pX[i + 2], pY[i + 2], fmaUnroll[2]);
         case 2:
            fmaUnroll[1] = _mm256_fmadd_pd(pX[i + 1], pY[i + 1], fmaUnroll[1]);
         case 1:
            fmaUnroll[0] = _mm256_fmadd_pd(pX[i], pY[i], fmaUnroll[0]);
      }
      ans = _mm256_add_pd(_mm256_add_pd(fmaUnroll[0], fmaUnroll[1]), _mm256_add_pd(fmaUnroll[2], fmaUnroll[3]));
   } else {
      ans = _mm256_setzero_pd();
      switch (avx2len) {
         case 3:
            ans = _mm256_mul_pd(pX[2], pY[2]);
         case 2:
            ans = _mm256_fmadd_pd(pX[1], pY[1], ans);
         case 1:
            ans = _mm256_fmadd_pd(pX[0], pY[0], ans);
      }
   }
   const __m128d s = _mm_add_pd(_mm256_extractf128_pd(ans, 0), _mm256_extractf128_pd(ans, 1));
   const double *pComps = (const double *)&s;
   double dot = pComps[0] + pComps[1];
   const ae_int_t tail = avx2len << 2;
   switch (n - tail) {
      case 1:
         dot += x[tail] * y[tail];
      break;
      case 2:
         dot += x[tail] * y[tail];
         dot += x[tail + 1] * y[tail + 1];
      break;
      case 3:
         dot += x[tail] * y[tail];
         dot += x[tail + 1] * y[tail + 1];
         dot += x[tail + 2] * y[tail + 2];
      break;
   }
   return dot;
}

double fma_rdotv2(const ae_int_t n, const Real *__restrict x) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const ae_int_t fmaLen = (avx2len >> 2) << 2;
   const __m256d *__restrict pX = (const __m256d *)x;
   __m256d ans;
   if (fmaLen >= 4) {
   //TODO: this can be further unrolled to 8 because AVX(2) provides 16 registers
      __m256d fmaUnroll[4];
      fmaUnroll[0] = _mm256_mul_pd(pX[0], pX[0]);
      fmaUnroll[1] = _mm256_mul_pd(pX[1], pX[1]);
      fmaUnroll[2] = _mm256_mul_pd(pX[2], pX[2]);
      fmaUnroll[3] = _mm256_mul_pd(pX[3], pX[3]);
      for (i = 4; i < fmaLen; i += 4) {
         fmaUnroll[0] = _mm256_fmadd_pd(pX[i], pX[i], fmaUnroll[0]);
         fmaUnroll[1] = _mm256_fmadd_pd(pX[i + 1], pX[i + 1], fmaUnroll[1]);
         fmaUnroll[2] = _mm256_fmadd_pd(pX[i + 2], pX[i + 2], fmaUnroll[2]);
         fmaUnroll[3] = _mm256_fmadd_pd(pX[i + 3], pX[i + 3], fmaUnroll[3]);
      }
      switch (avx2len - fmaLen) {
         case 3:
            fmaUnroll[2] = _mm256_fmadd_pd(pX[i + 2], pX[i + 2], fmaUnroll[2]);
         case 2:
            fmaUnroll[1] = _mm256_fmadd_pd(pX[i + 1], pX[i + 1], fmaUnroll[1]);
         case 1:
            fmaUnroll[0] = _mm256_fmadd_pd(pX[i], pX[i], fmaUnroll[0]);
      }
      ans = _mm256_add_pd(_mm256_add_pd(fmaUnroll[0], fmaUnroll[1]), _mm256_add_pd(fmaUnroll[2], fmaUnroll[3]));
   } else {
      ans = _mm256_setzero_pd();
      switch (avx2len) {
         case 3:
            ans = _mm256_mul_pd(pX[2], pX[2]);
         case 2:
            ans = _mm256_fmadd_pd(pX[1], pX[1], ans);
         case 1:
            ans = _mm256_fmadd_pd(pX[0], pX[0], ans);
      }
   }
   const __m128d s = _mm_add_pd(_mm256_extractf128_pd(ans, 0), _mm256_extractf128_pd(ans, 1));
   const double *pComps = (const double *)&s;
   double dot = pComps[0] + pComps[1];
   const ae_int_t tail = avx2len << 2;
   switch (n - tail) {
      case 1:
         dot += x[tail] * x[tail];
      break;
      case 2:
         dot += x[tail] * x[tail];
         dot += x[tail + 1] * x[tail + 1];
      break;
      case 3:
         dot += x[tail] * x[tail];
         dot += x[tail + 1] * x[tail + 1];
         dot += x[tail + 2] * x[tail + 2];
      break;
   }
   return dot;
}

void fma_raddv(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const __m256d *__restrict pSrc = (const __m256d *)y;
   __m256d *__restrict pDest = (__m256d *)x;
   const __m256d avx2alpha = _mm256_set1_pd(alpha);
   for (i = 0; i < avx2len; i++) {
      pDest[i] = _mm256_fmadd_pd(avx2alpha, pSrc[i], pDest[i]);
   }
   const ae_int_t tail = avx2len << 2;
   switch (n - tail) {
      case 1:
         *(double *)(pDest + i) += alpha * *(const double *)(pSrc + i);
      break;
      case 2:
         *(__m128d *)(pDest + i) = _mm_fmadd_pd(_mm256_extractf128_pd(avx2alpha, 0), *(const __m128d *)(pSrc + i), *(const __m128d *)(pDest + i));
      break;
      case 3:
         *(__m128d *)(pDest + i) = _mm_fmadd_pd(_mm256_extractf128_pd(avx2alpha, 0), *(const __m128d *)(pSrc + i), *(const __m128d *)(pDest + i));
         x[tail + 2] += alpha * y[tail + 2];
      break;
   }
}

static void raddvx_fma3avx_xaligned(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x) {
   ae_int_t i;
   const ae_int_t vecLen = (n >> 2) << 2;
   const __m256d avx2alpha = _mm256_set1_pd(alpha);
   __m256d *__restrict pDest = (__m256d *)x;
   for (i = 0; i < vecLen; i += 4) {
      const ae_int_t iDest = i >> 2;
      pDest[iDest] = _mm256_fmadd_pd(avx2alpha, _mm256_loadu_pd(y + i), pDest[iDest]);
   }
   switch (n - vecLen) {
      case 1:
         x[i] += alpha * y[i];
      break;
      case 2: {
         const ae_int_t iDest = i >> 2;
         *(__m128d *)(pDest + iDest) = _mm_fmadd_pd(_mm256_extractf128_pd(avx2alpha, 0), _mm_loadu_pd(y + i), *(const __m128d *)(pDest + iDest));
      }
      break;
      case 3: {
         const ae_int_t iDest = i >> 2;
         *(__m128d *)(pDest + iDest) = _mm_fmadd_pd(_mm256_extractf128_pd(avx2alpha, 0), _mm_loadu_pd(y + i), *(const __m128d *)(pDest + iDest));
         x[i + 2] += alpha * y[i + 2];
      }
      break;
   }
}

void fma_raddvx(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x) {
   const ptrdiff_t unal = (ptrdiff_t)x & 31;
   if (n <= 4) {
      ae_int_t i;
      for (i = 0; i < n; i++)
         x[i] += alpha * y[i];
      return;
   }
   switch (unal) {
      case 0:
         raddvx_fma3avx_xaligned(n, alpha, y, x);
      return;
      case 8:
         x[2] += alpha * y[2];
      case 16:
         x[1] += alpha * y[1];
      case 24: {
         x[0] += alpha * y[0];
         const ptrdiff_t nDone = 4 - (unal >> 3);
         raddvx_fma3avx_xaligned(n - nDone, alpha, y + nDone, x + nDone);
      }
      return;
   }
}

void fma_rmuladdv(const ae_int_t n, const double *__restrict y, const double *__restrict z, double *__restrict x) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const __m256d *__restrict pSrc0 = (const __m256d *)y;
   const __m256d *__restrict pSrc1 = (const __m256d *)z;
   __m256d *__restrict pDest = (__m256d *)x;
   for (i = 0; i < avx2len; i++)
      pDest[i] = _mm256_fmadd_pd(pSrc0[i], pSrc1[i], pDest[i]);
   const ae_int_t tail = avx2len << 2;
   for (i = tail; i < n; i++)
      x[i] += y[i] * z[i];
}

void fma_rnegmuladdv(const ae_int_t n, const double *__restrict y, const double *__restrict z, double *__restrict x) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const __m256d *__restrict pSrc0 = (const __m256d *)y;
   const __m256d *__restrict pSrc1 = (const __m256d *)z;
   __m256d *__restrict pDest = (__m256d *)x;
   for (i = 0; i < avx2len; i++)
      pDest[i] = _mm256_fnmadd_pd(pSrc0[i], pSrc1[i], pDest[i]);
   const ae_int_t tail = avx2len << 2;
   for (i = tail; i < n; i++)
      x[i] -= y[i] * z[i];
}

void fma_rcopymuladdv(const ae_int_t n, const double *__restrict y, const double *__restrict z, const double *__restrict x, double *__restrict r) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const __m256d *__restrict pSrc0 = (const __m256d *)y;
   const __m256d *__restrict pSrc1 = (const __m256d *)z;
   const __m256d *__restrict pSrc2 = (const __m256d *)x;
   __m256d *__restrict pDest = (__m256d *)r;
   for (i = 0; i < avx2len; i++)
      pDest[i] = _mm256_fmadd_pd(pSrc0[i], pSrc1[i], pSrc2[i]);
   const ae_int_t tail = avx2len << 2;
   for (i = tail; i < n; i++)
      r[i] = x[i] + y[i] * z[i];
}

void fma_rcopynegmuladdv(const ae_int_t n, const double *__restrict y, const double *__restrict z, const double *__restrict x, double *__restrict r) {
   ae_int_t i;
   const ae_int_t avx2len = n >> 2;
   const __m256d *__restrict pSrc0 = (const __m256d *)y;
   const __m256d *__restrict pSrc1 = (const __m256d *)z;
   const __m256d *__restrict pSrc2 = (const __m256d *)x;
   __m256d *__restrict pDest = (__m256d *)r;
   for (i = 0; i < avx2len; i++)
      pDest[i] = _mm256_fnmadd_pd(pSrc0[i], pSrc1[i], pSrc2[i]);
   const ae_int_t tail = avx2len << 2;
   for (i = tail; i < n; i++)
      r[i] = x[i] - y[i] * z[i];
}

void fma_rgemv_straight(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   ae_int_t j;
   const __m256d *__restrict pX = (const __m256d *)x;
   const ae_int_t nVec = n >> 2;
   const ae_int_t nUnroll = nVec >> 3;
   for (i = 0; i < m; i++) {
      const __m256d *__restrict pRow = (const __m256d *)a->xyR[i];
      __m256d sum;
      if (nUnroll >= 1) {
         __m256d u0 = _mm256_mul_pd(pRow[0], pX[0]);
         __m256d u1 = _mm256_mul_pd(pRow[1], pX[1]);
         __m256d u2 = _mm256_mul_pd(pRow[2], pX[2]);
         __m256d u3 = _mm256_mul_pd(pRow[3], pX[3]);
         __m256d u4 = _mm256_mul_pd(pRow[4], pX[4]);
         __m256d u5 = _mm256_mul_pd(pRow[5], pX[5]);
         __m256d u6 = _mm256_mul_pd(pRow[6], pX[6]);
         __m256d u7 = _mm256_mul_pd(pRow[7], pX[7]);
         for (j = 1; j < nUnroll; j++) {
            const ae_int_t at = j << 3;
            u0 = _mm256_fmadd_pd(pRow[at], pX[at], u0);
            u1 = _mm256_fmadd_pd(pRow[at + 1], pX[at + 1], u1);
            u2 = _mm256_fmadd_pd(pRow[at + 2], pX[at + 2], u2);
            u3 = _mm256_fmadd_pd(pRow[at + 3], pX[at + 3], u3);
            u4 = _mm256_fmadd_pd(pRow[at + 4], pX[at + 4], u4);
            u5 = _mm256_fmadd_pd(pRow[at + 5], pX[at + 5], u5);
            u6 = _mm256_fmadd_pd(pRow[at + 6], pX[at + 6], u6);
            u7 = _mm256_fmadd_pd(pRow[at + 7], pX[at + 7], u7);
         }
         const ae_int_t at = j << 3;
         switch (nVec - at) {
            case 7:
               u6 = _mm256_fmadd_pd(pX[at + 6], pRow[at + 6], u6);
            case 6:
               u5 = _mm256_fmadd_pd(pX[at + 5], pRow[at + 5], u5);
            case 5:
               u4 = _mm256_fmadd_pd(pX[at + 4], pRow[at + 4], u4);
            case 4:
               u3 = _mm256_fmadd_pd(pX[at + 3], pRow[at + 3], u3);
            case 3:
               u2 = _mm256_fmadd_pd(pX[at + 2], pRow[at + 2], u2);
            case 2:
               u1 = _mm256_fmadd_pd(pX[at + 1], pRow[at + 1], u1);
            case 1:
               u0 = _mm256_fmadd_pd(pX[at], pRow[at], u0);
         }
         sum = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(u0, u1), _mm256_add_pd(u2, u3)), _mm256_add_pd(_mm256_add_pd(u4, u5), _mm256_add_pd(u6, u7)));
      } else {
         switch (nVec) {
            case 0:
               sum = _mm256_setzero_pd();
            break;
            case 1:
               sum = _mm256_mul_pd(pX[0], pRow[0]);
            break;
            case 2:
               sum = _mm256_fmadd_pd(pX[0], pRow[0], _mm256_mul_pd(pX[1], pRow[1]));
            break;
            case 3:
               sum = _mm256_fmadd_pd(pX[2], pRow[2], _mm256_fmadd_pd(pX[0], pRow[0], _mm256_mul_pd(pX[1], pRow[1])));
            break;
            case 4:
               sum = _mm256_add_pd(_mm256_fmadd_pd(pX[0], pRow[0], _mm256_mul_pd(pX[1], pRow[1])), _mm256_fmadd_pd(pX[2], pRow[2], _mm256_mul_pd(pX[3], pRow[3])));
            break;
            case 5:
               sum = _mm256_fmadd_pd(pX[4], pRow[4], _mm256_add_pd(_mm256_fmadd_pd(pX[0], pRow[0], _mm256_mul_pd(pX[1], pRow[1])), _mm256_fmadd_pd(pX[2], pRow[2], _mm256_mul_pd(pX[3], pRow[3]))));
            break;
            case 6:
               sum = _mm256_add_pd(_mm256_add_pd(_mm256_fmadd_pd(pX[0], pRow[0], _mm256_mul_pd(pX[1], pRow[1])), _mm256_fmadd_pd(pX[2], pRow[2], _mm256_mul_pd(pX[3], pRow[3]))), _mm256_fmadd_pd(pX[4], pRow[4], _mm256_mul_pd(pX[5], pRow[5])));
            break;
            case 7:
               sum = _mm256_add_pd(_mm256_add_pd(_mm256_fmadd_pd(pX[0], pRow[0], _mm256_mul_pd(pX[1], pRow[1])), _mm256_fmadd_pd(pX[2], pRow[2], _mm256_mul_pd(pX[3], pRow[3]))), _mm256_fmadd_pd(pX[6], pRow[6], _mm256_fmadd_pd(pX[4], pRow[4], _mm256_mul_pd(pX[5], pRow[5]))));
            break;
         }
      }
      const __m128d t = _mm_add_pd(_mm256_extractf128_pd(sum, 0), _mm256_extractf128_pd(sum, 1));
      const double *pComps = (const double *)&t;
      double ans = pComps[0] + pComps[1];
      const ae_int_t tail = nVec << 2;
      for (j = tail; j < n; j++) {
         ans += a->xyR[i][j] * x[j];
      }
      y[i] += alpha * ans;
   }
}

void fma_rgemv_transposed(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   ae_int_t j;
   __m256d *__restrict pY = (__m256d *)y;
   const ae_int_t nVec = m >> 2;
   for (i = 0; i < n; i++) {
      const __m256d *__restrict pRow = (const __m256d *)a->xyR[i];
      const double v = alpha * x[i];
      const __m256d vV = _mm256_set1_pd(v);
      for (j = 0; j < nVec; j++) {
         pY[j] = _mm256_fmadd_pd(vV, pRow[j], pY[j]);
      }
      const ae_int_t tail = nVec << 2;
      for (j = tail; j < m; j++) {
         y[j] += v * a->xyR[i][j];
      }
   }
}

static void fma_rgemvx_straight_xaligned(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   ae_int_t j;
   const __m256d *__restrict pX = (const __m256d *)x;
   const ae_int_t nVec = n >> 2;
   const ae_int_t nUnroll = nVec >> 3;
   for (i = 0; i < m; i++) {
      const __m256d *__restrict pRow = (const __m256d *)(a->xyR[i + ia] + ja);
      __m256d sum;
      if (nUnroll >= 1) {
         __m256d u0 = _mm256_mul_pd(ULOAD256PD(pRow[0]), pX[0]);
         __m256d u1 = _mm256_mul_pd(ULOAD256PD(pRow[1]), pX[1]);
         __m256d u2 = _mm256_mul_pd(ULOAD256PD(pRow[2]), pX[2]);
         __m256d u3 = _mm256_mul_pd(ULOAD256PD(pRow[3]), pX[3]);
         __m256d u4 = _mm256_mul_pd(ULOAD256PD(pRow[4]), pX[4]);
         __m256d u5 = _mm256_mul_pd(ULOAD256PD(pRow[5]), pX[5]);
         __m256d u6 = _mm256_mul_pd(ULOAD256PD(pRow[6]), pX[6]);
         __m256d u7 = _mm256_mul_pd(ULOAD256PD(pRow[7]), pX[7]);
         for (j = 1; j < nUnroll; j++) {
            const ae_int_t at = j << 3;
            u0 = _mm256_fmadd_pd(ULOAD256PD(pRow[at]), pX[at], u0);
            u1 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 1]), pX[at + 1], u1);
            u2 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 2]), pX[at + 2], u2);
            u3 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 3]), pX[at + 3], u3);
            u4 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 4]), pX[at + 4], u4);
            u5 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 5]), pX[at + 5], u5);
            u6 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 6]), pX[at + 6], u6);
            u7 = _mm256_fmadd_pd(ULOAD256PD(pRow[at + 7]), pX[at + 7], u7);
         }
         const ae_int_t at = j << 3;
         switch (nVec - at) {
            case 7:
               u6 = _mm256_fmadd_pd(pX[at + 6], ULOAD256PD(pRow[at + 6]), u6);
            case 6:
               u5 = _mm256_fmadd_pd(pX[at + 5], ULOAD256PD(pRow[at + 5]), u5);
            case 5:
               u4 = _mm256_fmadd_pd(pX[at + 4], ULOAD256PD(pRow[at + 4]), u4);
            case 4:
               u3 = _mm256_fmadd_pd(pX[at + 3], ULOAD256PD(pRow[at + 3]), u3);
            case 3:
               u2 = _mm256_fmadd_pd(pX[at + 2], ULOAD256PD(pRow[at + 2]), u2);
            case 2:
               u1 = _mm256_fmadd_pd(pX[at + 1], ULOAD256PD(pRow[at + 1]), u1);
            case 1:
               u0 = _mm256_fmadd_pd(pX[at], ULOAD256PD(pRow[at]), u0);
         }
         sum = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(u0, u1), _mm256_add_pd(u2, u3)), _mm256_add_pd(_mm256_add_pd(u4, u5), _mm256_add_pd(u6, u7)));
      } else {
         switch (nVec) {
            case 0:
               sum = _mm256_setzero_pd();
            break;
            case 1:
               sum = _mm256_mul_pd(pX[0], ULOAD256PD(pRow[0]));
            break;
            case 2:
               sum = _mm256_fmadd_pd(pX[0], ULOAD256PD(pRow[0]), _mm256_mul_pd(pX[1], ULOAD256PD(pRow[1])));
            break;
            case 3:
               sum = _mm256_fmadd_pd(pX[2], ULOAD256PD(pRow[2]), _mm256_fmadd_pd(pX[0], ULOAD256PD(pRow[0]), _mm256_mul_pd(pX[1], ULOAD256PD(pRow[1]))));
            break;
            case 4:
               sum = _mm256_add_pd(_mm256_fmadd_pd(pX[0], ULOAD256PD(pRow[0]), _mm256_mul_pd(pX[1], ULOAD256PD(pRow[1]))), _mm256_fmadd_pd(pX[2], ULOAD256PD(pRow[2]), _mm256_mul_pd(pX[3], ULOAD256PD(pRow[3]))));
            break;
            case 5:
               sum = _mm256_fmadd_pd(pX[4], ULOAD256PD(pRow[4]), _mm256_add_pd(_mm256_fmadd_pd(pX[0], ULOAD256PD(pRow[0]), _mm256_mul_pd(pX[1], ULOAD256PD(pRow[1]))), _mm256_fmadd_pd(pX[2], ULOAD256PD(pRow[2]), _mm256_mul_pd(pX[3], ULOAD256PD(pRow[3])))));
            break;
            case 6:
               sum = _mm256_add_pd(_mm256_add_pd(_mm256_fmadd_pd(pX[0], ULOAD256PD(pRow[0]), _mm256_mul_pd(pX[1], ULOAD256PD(pRow[1]))), _mm256_fmadd_pd(pX[2], ULOAD256PD(pRow[2]), _mm256_mul_pd(pX[3], ULOAD256PD(pRow[3])))), _mm256_fmadd_pd(pX[4], ULOAD256PD(pRow[4]), _mm256_mul_pd(pX[5], ULOAD256PD(pRow[5]))));
            break;
            case 7:
               sum = _mm256_add_pd(_mm256_add_pd(_mm256_fmadd_pd(pX[0], ULOAD256PD(pRow[0]), _mm256_mul_pd(pX[1], ULOAD256PD(pRow[1]))), _mm256_fmadd_pd(pX[2], ULOAD256PD(pRow[2]), _mm256_mul_pd(pX[3], ULOAD256PD(pRow[3])))), _mm256_fmadd_pd(pX[6], ULOAD256PD(pRow[6]), _mm256_fmadd_pd(pX[4], ULOAD256PD(pRow[4]), _mm256_mul_pd(pX[5], ULOAD256PD(pRow[5])))));
            break;
         }
      }
      const __m128d t = _mm_add_pd(_mm256_extractf128_pd(sum, 0), _mm256_extractf128_pd(sum, 1));
      const double *pComps = (const double *)&t;
      double ans = pComps[0] + pComps[1];
      const ae_int_t tail = nVec << 2;
      for (j = tail; j < n; j++) {
         ans += a->xyR[i + ia][j + ja] * x[j];
      }
      y[i] += alpha * ans;
   }
}

void fma_rgemvx_straight(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   ae_int_t j;
   if (n <= 3) {
      for (i = 0; i < m; i++) {
         const double *p_a = a->xyR[ia + i] + ja;
         double v = 0.0;
         for (j = 0; j < n; j++) {
            v += p_a[j] * x[j];
         }
         y[i] += alpha * v;
      }
      return;
   }
   const ptrdiff_t unal = (ptrdiff_t)x & 31;
   if (unal == 0) {
      fma_rgemvx_straight_xaligned(m, n, alpha, a, ia, ja, x, y);
      return;
   }
   const ptrdiff_t shift = 4 - (unal >> 3);
   for (i = 0; i < m; i++) {
      const double *p_a = a->xyR[ia + i] + ja;
      double v = 0.0;
      for (j = 0; j < shift; j++) {
         v += p_a[j] * x[j];
      }
      y[i] += alpha * v;
   }
   fma_rgemvx_straight_xaligned(m, n - shift, alpha, a, ia, ja + shift, x + shift, y);
}

static void fma_rgemvx_transposed_yaligned(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   ae_int_t j;
   __m256d *__restrict pY = (__m256d *)y;
   const ae_int_t nVec = m >> 2;
   for (i = 0; i < n; i++) {
      const __m256d *__restrict pRow = (const __m256d *)(a->xyR[i + ia] + ja);
      const double v = alpha * x[i];
      const __m256d vV = _mm256_set1_pd(v);
      for (j = 0; j < nVec; j++) {
         pY[j] = _mm256_fmadd_pd(vV, ULOAD256PD(pRow[j]), pY[j]);
      }
      const ae_int_t tail = nVec << 2;
      for (j = tail; j < m; j++) {
         y[j] += v * a->xyR[i + ia][j + ja];
      }
   }
}

void fma_rgemvx_transposed(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y) {
   ae_int_t i;
   ae_int_t j;
   if (m <= 3) {
      for (i = 0; i < n; i++) {
         const double *p_a = a->xyR[ia + i] + ja;
         const double v = alpha * x[i];
         for (j = 0; j < m; j++) {
            y[j] += v * p_a[j];
         }
      }
      return;
   }
   const ptrdiff_t unal = (ptrdiff_t)y & 31;
   if (unal == 0) {
      fma_rgemvx_transposed_yaligned(m, n, alpha, a, ia, ja, x, y);
      return;
   }
   const ptrdiff_t shift = 4 - (unal >> 3);
   for (i = 0; i < n; i++) {
      const double *p_a = a->xyR[ia + i] + ja;
      const double v = alpha * x[i];
      for (j = 0; j < shift; j++) {
         y[j] += v * p_a[j];
      }
   }
   fma_rgemvx_transposed_yaligned(m - shift, n, alpha, a, ia, ja + shift, x, y + shift);
}

// Computes  product   A*transpose(B)  of two MICRO_SIZE*ROUND_LENGTH rowwise
// 'horizontal' matrices, stored with stride == block_size, and writes it to the
// row-wise matrix C.
//
// ROUND_LENGTH is expected to be properly SIMD-rounded length,  as  returned
// by avx2_ablasf_packblkh().
//
// Present version of the function supports only MICRO_SIZE == 2,  the  behavior
// is undefined for other micro sizes.
//
// Requires AVX2, does NOT check its presense.
// ALGLIB Routine: Copyright 19.07.2021 by Sergey Bochkanov
void fma_ablasf_dotblkh(const double *src_a, const double *src_b, ae_int_t round_length, ae_int_t block_size, ae_int_t micro_size, double *dst, ae_int_t dst_stride) {
   ae_int_t z;
   __m256d r00 = _mm256_setzero_pd(), r01 = _mm256_setzero_pd(), r10 = _mm256_setzero_pd(), r11 = _mm256_setzero_pd();
   if (round_length & 0x7) {
   // round_length is multiple of 4, but not multiple of 8
      for (z = 0; z < round_length; z += 4, src_a += 4, src_b += 4) {
         __m256d a0 = _mm256_load_pd(src_a);
         __m256d a1 = _mm256_load_pd(src_a + block_size);
         __m256d b0 = _mm256_load_pd(src_b);
         __m256d b1 = _mm256_load_pd(src_b + block_size);
         r00 = _mm256_fmadd_pd(a0, b0, r00);
         r01 = _mm256_fmadd_pd(a0, b1, r01);
         r10 = _mm256_fmadd_pd(a1, b0, r10);
         r11 = _mm256_fmadd_pd(a1, b1, r11);
      }
   } else {
   // round_length is multiple of 8
      for (z = 0; z < round_length; z += 8, src_a += 8, src_b += 8) {
         __m256d a0 = _mm256_load_pd(src_a);
         __m256d a1 = _mm256_load_pd(src_a + block_size);
         __m256d b0 = _mm256_load_pd(src_b);
         __m256d b1 = _mm256_load_pd(src_b + block_size);
         __m256d c0 = _mm256_load_pd(src_a + 4);
         __m256d c1 = _mm256_load_pd(src_a + block_size + 4);
         __m256d d0 = _mm256_load_pd(src_b + 4);
         __m256d d1 = _mm256_load_pd(src_b + block_size + 4);
         r00 = _mm256_fmadd_pd(c0, d0, _mm256_fmadd_pd(a0, b0, r00));
         r01 = _mm256_fmadd_pd(c0, d1, _mm256_fmadd_pd(a0, b1, r01));
         r10 = _mm256_fmadd_pd(c1, d0, _mm256_fmadd_pd(a1, b0, r10));
         r11 = _mm256_fmadd_pd(c1, d1, _mm256_fmadd_pd(a1, b1, r11));
      }
   }
   __m256d sum0 = _mm256_hadd_pd(r00, r01);
   __m256d sum1 = _mm256_hadd_pd(r10, r11);
   _mm_store_pd(dst, _mm_add_pd(_mm256_castpd256_pd128(sum0), _mm256_extractf128_pd(sum0, 1)));
   _mm_store_pd(dst + dst_stride, _mm_add_pd(_mm256_castpd256_pd128(sum1), _mm256_extractf128_pd(sum1, 1)));
}

// Solving linear system: propagating computed supernode.
//
// Propagates computed supernode to the rest of the RHS  using  SIMD-friendly
// RHS storage format.
// ALGLIB Routine: Copyright 08.09.2021 by Sergey Bochkanov
void fma_spchol_propagatefwd(RVector *x, ae_int_t cols0, ae_int_t blocksize, ZVector *superrowidx, ae_int_t rbase, ae_int_t offdiagsize, RVector *rowstorage, ae_int_t offss, ae_int_t sstride, RVector *simdbuf, ae_int_t simdwidth) {
   ae_int_t k;
   ae_assert(simdwidth == 4, "SPCHOL: unexpected stride in propagatefwd()");
   if (sstride == 4) {
   // blocksize is either 3 or 4
      if (CurCPU & CPU_FMA) {
         __m256d x_simd;
         double *p_mat_row = rowstorage->xR + offss + blocksize * 4;
         double *p_simd_buf = simdbuf->xR;
         ae_int_t *p_rowidx = superrowidx->xZ + rbase;
         if (blocksize == 3) {
            double xx[4];
            xx[0] = x->xR[cols0];
            xx[1] = x->xR[cols0 + 1];
            xx[2] = x->xR[cols0 + 2];
            xx[3] = 0;
            x_simd = _mm256_loadu_pd(xx);
         } else
            x_simd = _mm256_loadu_pd(x->xR + cols0);
         for (k = 0; k < offdiagsize; k++) {
            double *p_buf_row = p_simd_buf + p_rowidx[k] * 4;
            _mm256_store_pd(p_buf_row, _mm256_fnmadd_pd(_mm256_load_pd(p_mat_row), x_simd, _mm256_load_pd(p_buf_row)));
            p_mat_row += 4;
         }
         return;
      }
   }
   if (blocksize == 2 && sstride == 2) {
   // blocksize is 2, stride is 2
      if (CurCPU & CPU_FMA) {
         __m128d x_simd = _mm_loadu_pd(x->xR + cols0);
         double *p_mat_row = rowstorage->xR + offss + 2 * 2;
         double *p_simd_buf = simdbuf->xR;
         ae_int_t *p_rowidx = superrowidx->xZ + rbase;
         for (k = 0; k < offdiagsize; k++) {
            double *p_buf_row = p_simd_buf + p_rowidx[k] * 4;
            _mm_store_pd(p_buf_row, _mm_fnmadd_pd(_mm_load_pd(p_mat_row), x_simd, _mm_load_pd(p_buf_row)));
            p_mat_row += 2;
         }
         return;
      }
   }
   ae_assert(false, "fma_spchol_propagatefwd: unsupported input");
}

bool fma_spchol_updatekernelabc4(double *rowstorage, ae_int_t offss, ae_int_t twidth, ae_int_t offsu, ae_int_t uheight, ae_int_t urank, ae_int_t urowstride, ae_int_t uwidth, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase) {
   ae_int_t k;
   ae_int_t targetrow;
   ae_int_t targetcol;
// Filter out unsupported combinations (ones that are too sparse for the non-SIMD code).
   if (twidth < 3 || twidth > 4) {
      return false;
   }
   if (uwidth < 1 || uwidth > 4) {
      return false;
   }
   if (urank > 4) {
      return false;
   }
// Shift input arrays to the beginning of the working area.
// Prepare SIMD masks
   __m256i v_rankmask = _mm256_cmpgt_epi64(_mm256_set_epi64x(urank, urank, urank, urank), _mm256_set_epi64x(3, 2, 1, 0));
   double *update_storage = rowstorage + offsu;
   double *target_storage = rowstorage + offss;
   superrowidx += urbase;
// Load head of the update matrix
   __m256d v_d0123 = _mm256_maskload_pd(diagd + offsd, v_rankmask);
   __m256d u_0_0123 = _mm256_setzero_pd();
   __m256d u_1_0123 = _mm256_setzero_pd();
   __m256d u_2_0123 = _mm256_setzero_pd();
   __m256d u_3_0123 = _mm256_setzero_pd();
   for (k = 0; k < uwidth; k++) {
      targetcol = raw2smap[superrowidx[k]];
      if (targetcol == 0)
         u_0_0123 = _mm256_mul_pd(v_d0123, _mm256_maskload_pd(update_storage + k * urowstride, v_rankmask));
      if (targetcol == 1)
         u_1_0123 = _mm256_mul_pd(v_d0123, _mm256_maskload_pd(update_storage + k * urowstride, v_rankmask));
      if (targetcol == 2)
         u_2_0123 = _mm256_mul_pd(v_d0123, _mm256_maskload_pd(update_storage + k * urowstride, v_rankmask));
      if (targetcol == 3)
         u_3_0123 = _mm256_mul_pd(v_d0123, _mm256_maskload_pd(update_storage + k * urowstride, v_rankmask));
   }
// Transpose head
   __m256d u01_lo = _mm256_unpacklo_pd(u_0_0123, u_1_0123);
   __m256d u01_hi = _mm256_unpackhi_pd(u_0_0123, u_1_0123);
   __m256d u23_lo = _mm256_unpacklo_pd(u_2_0123, u_3_0123);
   __m256d u23_hi = _mm256_unpackhi_pd(u_2_0123, u_3_0123);
   __m256d u_0123_0 = _mm256_permute2f128_pd(u01_lo, u23_lo, 0x20);
   __m256d u_0123_1 = _mm256_permute2f128_pd(u01_hi, u23_hi, 0x20);
   __m256d u_0123_2 = _mm256_permute2f128_pd(u23_lo, u01_lo, 0x13);
   __m256d u_0123_3 = _mm256_permute2f128_pd(u23_hi, u01_hi, 0x13);
// Run update
   if (urank == 1) {
      for (k = 0; k < uheight; k++) {
         targetrow = raw2smap[superrowidx[k]] * 4;
         double *update_row = rowstorage + offsu + k * urowstride;
         _mm256_store_pd(target_storage + targetrow, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row), u_0123_0, _mm256_load_pd(target_storage + targetrow)));
      }
   }
   if (urank == 2) {
      for (k = 0; k < uheight; k++) {
         targetrow = raw2smap[superrowidx[k]] * 4;
         double *update_row = rowstorage + offsu + k * urowstride;
         _mm256_store_pd(target_storage + targetrow, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row + 1), u_0123_1, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row), u_0123_0, _mm256_load_pd(target_storage + targetrow))));
      }
   }
   if (urank == 3) {
      for (k = 0; k < uheight; k++) {
         targetrow = raw2smap[superrowidx[k]] * 4;
         double *update_row = rowstorage + offsu + k * urowstride;
         _mm256_store_pd(target_storage + targetrow, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row + 2), u_0123_2, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row + 1), u_0123_1, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row), u_0123_0, _mm256_load_pd(target_storage + targetrow)))));
      }
   }
   if (urank == 4) {
      for (k = 0; k < uheight; k++) {
         targetrow = raw2smap[superrowidx[k]] * 4;
         double *update_row = rowstorage + offsu + k * urowstride;
         _mm256_store_pd(target_storage + targetrow, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row + 3), u_0123_3, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row + 2), u_0123_2, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row + 1), u_0123_1, _mm256_fnmadd_pd(_mm256_broadcast_sd(update_row), u_0123_0, _mm256_load_pd(target_storage + targetrow))))));
      }
   }
   return true;
}

bool fma_spchol_updatekernel4444(double *rowstorage, ae_int_t offss, ae_int_t sheight, ae_int_t offsu, ae_int_t uheight, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase) {
   ae_int_t k;
   ae_int_t targetrow;
   ae_int_t offsk;
   __m256d v_negd_u0, v_negd_u1, v_negd_u2, v_negd_u3, v_negd;
   __m256d v_w0, v_w1, v_w2, v_w3, u01_lo, u01_hi, u23_lo, u23_hi;
// Compute W = -D*transpose(U[0:3]).
   v_negd = _mm256_mul_pd(_mm256_loadu_pd(diagd + offsd), _mm256_set1_pd(-1.0));
   v_negd_u0 = _mm256_mul_pd(_mm256_load_pd(rowstorage + offsu), v_negd);
   v_negd_u1 = _mm256_mul_pd(_mm256_load_pd(rowstorage + offsu + 4), v_negd);
   v_negd_u2 = _mm256_mul_pd(_mm256_load_pd(rowstorage + offsu + 8), v_negd);
   v_negd_u3 = _mm256_mul_pd(_mm256_load_pd(rowstorage + offsu + 12), v_negd);
   u01_lo = _mm256_unpacklo_pd(v_negd_u0, v_negd_u1);
   u01_hi = _mm256_unpackhi_pd(v_negd_u0, v_negd_u1);
   u23_lo = _mm256_unpacklo_pd(v_negd_u2, v_negd_u3);
   u23_hi = _mm256_unpackhi_pd(v_negd_u2, v_negd_u3);
   v_w0 = _mm256_permute2f128_pd(u01_lo, u23_lo, 0x20);
   v_w1 = _mm256_permute2f128_pd(u01_hi, u23_hi, 0x20);
   v_w2 = _mm256_permute2f128_pd(u23_lo, u01_lo, 0x13);
   v_w3 = _mm256_permute2f128_pd(u23_hi, u01_hi, 0x13);
//
// Compute update S = S + row_scatter(U*W).
//
   if (sheight == uheight) {
   // No row scatter, the most efficient code
      for (k = 0; k < uheight; k++) {
         targetrow = offss + k * 4;
         offsk = offsu + k * 4;
         _mm256_store_pd(rowstorage + targetrow, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk + 3), v_w3, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk + 2), v_w2, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk + 1), v_w1, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk), v_w0, _mm256_load_pd(rowstorage + targetrow))))));
      }
   } else {
   // Row scatter is performed, less efficient code using double mapping to determine target row index
      for (k = 0; k < uheight; k++) {
         targetrow = offss + raw2smap[superrowidx[urbase + k]] * 4;
         offsk = offsu + k * 4;
         _mm256_store_pd(rowstorage + targetrow, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk + 3), v_w3, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk + 2), v_w2, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk + 1), v_w1, _mm256_fmadd_pd(_mm256_broadcast_sd(rowstorage + offsk), v_w0, _mm256_load_pd(rowstorage + targetrow))))));
      }
   }
   return true;
}
#endif // ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_AVX2_INTRINSICS
} // end of namespace alglib_impl
