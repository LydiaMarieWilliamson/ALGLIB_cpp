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
#include "Statistics.h"

// === BASESTAT Package ===
// Depends on: (AlgLibInternal) BASICSTATOPS
// Depends on: (LinAlg) ABLAS
namespace alglib_impl {
// Calculation of the distribution moments: mean, variance, skewness, kurtosis.
//
// Inputs:
//     X       -   sample
//     N       -   N >= 0, sample size:
//                 * if given, only leading N elements of X are processed
//                 * if not given, automatically determined from size of X
//
// Outputs:
//     Mean    -   mean.
//     Variance-   variance.
//     Skewness-   skewness (if variance != 0; zero otherwise).
//     Kurtosis-   kurtosis (if variance != 0; zero otherwise).
//
// NOTE: variance is calculated by dividing sum of squares by N-1, not N.
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis);
// API: void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis);
void samplemoments(RVector *x, ae_int_t n, double *mean, double *variance, double *skewness, double *kurtosis) {
   ae_int_t i;
   double v;
   double v1;
   double v2;
   double stddev;
   *mean = 0;
   *variance = 0;
   *skewness = 0;
   *kurtosis = 0;
   ae_assert(n >= 0, "SampleMoments: N<0");
   ae_assert(x->cnt >= n, "SampleMoments: Length(X)<N!");
   ae_assert(isfinitevector(x, n), "SampleMoments: X is not finite vector");
// Init, special case 'N=0'
   *mean = 0.0;
   *variance = 0.0;
   *skewness = 0.0;
   *kurtosis = 0.0;
   stddev = 0.0;
   if (n <= 0) {
      return;
   }
// Mean
   for (i = 0; i < n; i++) {
      *mean += x->xR[i];
   }
   *mean /= n;
// Variance (using corrected two-pass algorithm)
   if (n != 1) {
      v1 = 0.0;
      for (i = 0; i < n; i++) {
         v1 += sqr(x->xR[i] - (*mean));
      }
      v2 = 0.0;
      for (i = 0; i < n; i++) {
         v2 += x->xR[i] - (*mean);
      }
      v2 = sqr(v2) / n;
      *variance = (v1 - v2) / (n - 1);
      if (*variance < 0.0) {
         *variance = 0.0;
      }
      stddev = sqrt(*variance);
   }
// Skewness and kurtosis
   if (stddev != 0.0) {
      for (i = 0; i < n; i++) {
         v = (x->xR[i] - (*mean)) / stddev;
         v2 = sqr(v);
         *skewness += v2 * v;
         *kurtosis += sqr(v2);
      }
      *skewness /= n;
      *kurtosis = *kurtosis / n - 3;
   }
}

// Calculation of the mean.
//
// Inputs:
//     X       -   sample
//     N       -   N >= 0, sample size:
//                 * if given, only leading N elements of X are processed
//                 * if not given, automatically determined from size of X
//
// NOTE:
//
// This function return result  which calculated by 'SampleMoments' function
// and stored at 'Mean' variable.
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: double samplemean(const real_1d_array &x, const ae_int_t n);
// API: double samplemean(const real_1d_array &x);
double samplemean(RVector *x, ae_int_t n) {
   double mean;
   double tmp0;
   double tmp1;
   double tmp2;
   double result;
   samplemoments(x, n, &mean, &tmp0, &tmp1, &tmp2);
   result = mean;
   return result;
}

// Calculation of the variance.
//
// Inputs:
//     X       -   sample
//     N       -   N >= 0, sample size:
//                 * if given, only leading N elements of X are processed
//                 * if not given, automatically determined from size of X
//
// NOTE:
//
// This function return result  which calculated by 'SampleMoments' function
// and stored at 'Variance' variable.
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: double samplevariance(const real_1d_array &x, const ae_int_t n);
// API: double samplevariance(const real_1d_array &x);
double samplevariance(RVector *x, ae_int_t n) {
   double variance;
   double tmp0;
   double tmp1;
   double tmp2;
   double result;
   samplemoments(x, n, &tmp0, &variance, &tmp1, &tmp2);
   result = variance;
   return result;
}

// Calculation of the skewness.
//
// Inputs:
//     X       -   sample
//     N       -   N >= 0, sample size:
//                 * if given, only leading N elements of X are processed
//                 * if not given, automatically determined from size of X
//
// NOTE:
//
// This function return result  which calculated by 'SampleMoments' function
// and stored at 'Skewness' variable.
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: double sampleskewness(const real_1d_array &x, const ae_int_t n);
// API: double sampleskewness(const real_1d_array &x);
double sampleskewness(RVector *x, ae_int_t n) {
   double skewness;
   double tmp0;
   double tmp1;
   double tmp2;
   double result;
   samplemoments(x, n, &tmp0, &tmp1, &skewness, &tmp2);
   result = skewness;
   return result;
}

// Calculation of the kurtosis.
//
// Inputs:
//     X       -   sample
//     N       -   N >= 0, sample size:
//                 * if given, only leading N elements of X are processed
//                 * if not given, automatically determined from size of X
//
// NOTE:
//
// This function return result  which calculated by 'SampleMoments' function
// and stored at 'Kurtosis' variable.
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: double samplekurtosis(const real_1d_array &x, const ae_int_t n);
// API: double samplekurtosis(const real_1d_array &x);
double samplekurtosis(RVector *x, ae_int_t n) {
   double kurtosis;
   double tmp0;
   double tmp1;
   double tmp2;
   double result;
   samplemoments(x, n, &tmp0, &tmp1, &tmp2, &kurtosis);
   result = kurtosis;
   return result;
}

// ADev
//
// Inputs:
//     X   -   sample
//     N   -   N >= 0, sample size:
//             * if given, only leading N elements of X are processed
//             * if not given, automatically determined from size of X
//
// Outputs:
//     ADev-   ADev
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev);
// API: void sampleadev(const real_1d_array &x, double &adev);
void sampleadev(RVector *x, ae_int_t n, double *adev) {
   ae_int_t i;
   double mean;
   *adev = 0;
   ae_assert(n >= 0, "SampleADev: N<0");
   ae_assert(x->cnt >= n, "SampleADev: Length(X)<N!");
   ae_assert(isfinitevector(x, n), "SampleADev: X is not finite vector");
// Init, handle N=0
   mean = 0.0;
   *adev = 0.0;
   if (n <= 0) {
      return;
   }
// Mean
   for (i = 0; i < n; i++) {
      mean += x->xR[i];
   }
   mean /= n;
// ADev
   for (i = 0; i < n; i++) {
      *adev += fabs(x->xR[i] - mean);
   }
   *adev /= n;
}

// Median calculation.
//
// Inputs:
//     X   -   sample (array indexes: [0..N-1])
//     N   -   N >= 0, sample size:
//             * if given, only leading N elements of X are processed
//             * if not given, automatically determined from size of X
//
// Outputs:
//     Median
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
// API: void samplemedian(const real_1d_array &x, const ae_int_t n, double &median);
// API: void samplemedian(const real_1d_array &x, double &median);
void samplemedian(RVector *x, ae_int_t n, double *median) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t ir;
   ae_int_t j;
   ae_int_t l;
   ae_int_t midp;
   ae_int_t k;
   double a;
   ae_frame_make(&_frame_block);
   DupVector(x);
   *median = 0;
   ae_assert(n >= 0, "SampleMedian: N<0");
   ae_assert(x->cnt >= n, "SampleMedian: Length(X)<N!");
   ae_assert(isfinitevector(x, n), "SampleMedian: X is not finite vector");
// Some degenerate cases
   *median = 0.0;
   if (n <= 0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      *median = x->xR[0];
      ae_frame_leave();
      return;
   }
   if (n == 2) {
      *median = 0.5 * (x->xR[0] + x->xR[1]);
      ae_frame_leave();
      return;
   }
// Common case, N >= 3.
// Choose X[(N-1)/2]
   l = 0;
   ir = n - 1;
   k = (n - 1) / 2;
   while (true) {
      if (ir <= l + 1) {
      // 1 or 2 elements in partition
         if (ir == l + 1 && x->xR[ir] < x->xR[l]) {
            swapr(&x->xR[l], &x->xR[ir]);
         }
         break;
      } else {
         midp = (l + ir) / 2;
         swapr(&x->xR[midp], &x->xR[l + 1]);
         if (x->xR[l] > x->xR[ir]) {
            swapr(&x->xR[l], &x->xR[ir]);
         }
         if (x->xR[l + 1] > x->xR[ir]) {
            swapr(&x->xR[l + 1], &x->xR[ir]);
         }
         if (x->xR[l] > x->xR[l + 1]) {
            swapr(&x->xR[l], &x->xR[l + 1]);
         }
         i = l + 1;
         j = ir;
         a = x->xR[l + 1];
         while (true) {
            do {
               i++;
            } while (x->xR[i] < a);
            do {
               j--;
            } while (x->xR[j] > a);
            if (j < i) {
               break;
            }
            swapr(&x->xR[i], &x->xR[j]);
         }
         x->xR[l + 1] = x->xR[j];
         x->xR[j] = a;
         if (j >= k) {
            ir = j - 1;
         }
         if (j <= k) {
            l = i;
         }
      }
   }
// If N is odd, return result
   if (n % 2 == 1) {
      *median = x->xR[k];
      ae_frame_leave();
      return;
   }
   a = x->xR[n - 1];
   for (i = k + 1; i < n; i++) {
      if (x->xR[i] < a) {
         a = x->xR[i];
      }
   }
   *median = 0.5 * (x->xR[k] + a);
   ae_frame_leave();
}

// Percentile calculation.
//
// Inputs:
//     X   -   sample (array indexes: [0..N-1])
//     N   -   N >= 0, sample size:
//             * if given, only leading N elements of X are processed
//             * if not given, automatically determined from size of X
//     P   -   percentile (0 <= P <= 1)
//
// Outputs:
//     V   -   percentile
// ALGLIB: Copyright 01.03.2008 by Sergey Bochkanov
// API: void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v);
// API: void samplepercentile(const real_1d_array &x, const double p, double &v);
void samplepercentile(RVector *x, ae_int_t n, double p, double *v) {
   ae_frame _frame_block;
   ae_int_t i1;
   double t;
   ae_frame_make(&_frame_block);
   DupVector(x);
   *v = 0;
   NewVector(rbuf, 0, DT_REAL);
   ae_assert(n >= 0, "SamplePercentile: N<0");
   ae_assert(x->cnt >= n, "SamplePercentile: Length(X)<N!");
   ae_assert(isfinitevector(x, n), "SamplePercentile: X is not finite vector");
   ae_assert(isfinite(p), "SamplePercentile: incorrect P!");
   ae_assert(p >= 0.0 && p <= 1.0, "SamplePercentile: incorrect P!");
   tagsortfast(x, &rbuf, n);
   if (p == 0.0) {
      *v = x->xR[0];
      ae_frame_leave();
      return;
   }
   if (p == 1.0) {
      *v = x->xR[n - 1];
      ae_frame_leave();
      return;
   }
   t = p * (n - 1);
   i1 = FloorZ(t);
   t -= FloorZ(t);
   *v = x->xR[i1] * (1 - t) + x->xR[i1 + 1] * t;
   ae_frame_leave();
}

// Basecase code for RankData(), performs actual work on subset of data using
// temporary buffer passed as parameter.
//
// Inputs:
//     XY      -   array[NPoints,NFeatures], dataset
//     I0      -   index of first row to process
//     I1      -   index of past-the-last row to process;
//                 this function processes half-interval [I0,I1).
//     NFeatures-  number of features
//     IsCentered- whether ranks are centered or not:
//                 * True      -   ranks are centered in such way that  their
//                                 within-row sum is zero
//                 * False     -   ranks are not centered
//     Buf0    -   temporary buffers, may be empty (this function automatically
//                 allocates/reuses buffers).
//     Buf1    -   temporary buffers, may be empty (this function automatically
//                 allocates/reuses buffers).
//
// Outputs:
//     XY      -   data in [I0,I1) are replaced by their within-row ranks;
//                 ranking starts from 0, ends at NFeatures-1
// ALGLIB: Copyright 18.04.2013 by Sergey Bochkanov
static void basestat_rankdatabasecase(RMatrix *xy, ae_int_t i0, ae_int_t i1, ae_int_t nfeatures, bool iscentered, apbuffers *buf0, apbuffers *buf1) {
   ae_int_t i;
   ae_assert(i1 >= i0, "RankDataBasecase: internal error");
   if (buf1->ra0.cnt < nfeatures) {
      ae_vector_set_length(&buf1->ra0, nfeatures);
   }
   for (i = i0; i < i1; i++) {
      ae_v_move(buf1->ra0.xR, 1, xy->xyR[i], 1, nfeatures);
      rankx(&buf1->ra0, nfeatures, iscentered, buf0);
      ae_v_move(xy->xyR[i], 1, buf1->ra0.xR, 1, nfeatures);
   }
}

// Recurrent code for RankData(), splits problem into  subproblems  or  calls
// basecase code (depending on problem complexity).
//
// Inputs:
//     XY      -   array[NPoints,NFeatures], dataset
//     I0      -   index of first row to process
//     I1      -   index of past-the-last row to process;
//                 this function processes half-interval [I0,I1).
//     NFeatures-  number of features
//     IsCentered- whether ranks are centered or not:
//                 * True      -   ranks are centered in such way that  their
//                                 within-row sum is zero
//                 * False     -   ranks are not centered
//     Pool    -   shared pool which holds temporary buffers
//                 (APBuffers structure)
//     BasecaseCost-minimum cost of the problem which will be split
//
// Outputs:
//     XY      -   data in [I0,I1) are replaced by their within-row ranks;
//                 ranking starts from 0, ends at NFeatures-1
// ALGLIB: Copyright 18.04.2013 by Sergey Bochkanov
static void basestat_rankdatarec(RMatrix *xy, ae_int_t i0, ae_int_t i1, ae_int_t nfeatures, bool iscentered, ae_shared_pool *pool, ae_int_t basecasecost) {
   ae_frame _frame_block;
   double problemcost;
   ae_int_t im;
   ae_frame_make(&_frame_block);
   RefObj(apbuffers, buf0);
   RefObj(apbuffers, buf1);
   ae_assert(i1 >= i0, "RankDataRec: internal error");
// Try to activate parallelism
// Parallelism was tried if: i1 - i0 >= 4 && (double)(i1 - i0) * nfeatures * logbase2((double)nfeatures) >= smpactivationlevel()
// Recursively split problem, if it is too large
   problemcost = (double)(i1 - i0) * nfeatures * logbase2((double)nfeatures);
   if (i1 - i0 >= 2 && problemcost > spawnlevel()) {
      im = (i1 + i0) / 2;
      basestat_rankdatarec(xy, i0, im, nfeatures, iscentered, pool, basecasecost);
      basestat_rankdatarec(xy, im, i1, nfeatures, iscentered, pool, basecasecost);
      ae_frame_leave();
      return;
   }
// Retrieve buffers from pool, call serial code, return buffers to pool
   ae_shared_pool_retrieve(pool, &_buf0);
   ae_shared_pool_retrieve(pool, &_buf1);
   basestat_rankdatabasecase(xy, i0, i1, nfeatures, iscentered, buf0, buf1);
   ae_shared_pool_recycle(pool, &_buf0);
   ae_shared_pool_recycle(pool, &_buf1);
   ae_frame_leave();
}

// This function replaces data in XY by their ranks:
// * XY is processed row-by-row
// * rows are processed separately
// * tied data are correctly handled (tied ranks are calculated)
// * ranking starts from 0, ends at NFeatures-1
// * sum of within-row values is equal to (NFeatures-1)*NFeatures/2
//
// Inputs:
//     XY      -   array[NPoints,NFeatures], dataset
//     NPoints -   number of points
//     NFeatures-  number of features
//
// Outputs:
//     XY      -   data are replaced by their within-row ranks;
//                 ranking starts from 0, ends at NFeatures-1
// ALGLIB: Copyright 18.04.2013 by Sergey Bochkanov
// API: void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
// API: void rankdata(const real_2d_array &xy);
void rankdata(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures) {
   ae_frame _frame_block;
   ae_int_t basecasecost;
   ae_frame_make(&_frame_block);
   NewObj(apbuffers, buf0);
   NewObj(apbuffers, buf1);
   NewObj(ae_shared_pool, pool);
   ae_assert(npoints >= 0, "RankData: NPoints < 0");
   ae_assert(nfeatures >= 1, "RankData: NFeatures<1");
   ae_assert(xy->rows >= npoints, "RankData: Rows(XY)<NPoints");
   ae_assert(xy->cols >= nfeatures || npoints == 0, "RankData: Cols(XY)<NFeatures");
   ae_assert(apservisfinitematrix(xy, npoints, nfeatures), "RankData: XY contains infinite/NAN elements");
// Basecase cost is a maximum cost of basecase problems.
// Problems harded than that cost will be split.
//
// Problem cost is assumed to be NPoints*NFeatures*log2(NFeatures),
// which is proportional, but NOT equal to number of FLOPs required
// to solve problem.
//
// Try to use serial code for basecase problems, no SMP functionality, no shared pools.
   basecasecost = 10000;
   if ((double)npoints * nfeatures * logbase2((double)nfeatures) < (double)basecasecost) {
      basestat_rankdatabasecase(xy, 0, npoints, nfeatures, false, &buf0, &buf1);
      ae_frame_leave();
      return;
   }
// Parallel code
   ae_shared_pool_set_seed(&pool, &buf0, sizeof(buf0), apbuffers_init, apbuffers_copy, apbuffers_free);
   basestat_rankdatarec(xy, 0, npoints, nfeatures, false, &pool, basecasecost);
   ae_frame_leave();
}

// This function replaces data in XY by their CENTERED ranks:
// * XY is processed row-by-row
// * rows are processed separately
// * tied data are correctly handled (tied ranks are calculated)
// * centered ranks are just usual ranks, but centered in such way  that  sum
//   of within-row values is equal to 0.0.
// * centering is performed by subtracting mean from each row, i.e. it changes
//   mean value, but does NOT change higher moments
//
// Inputs:
//     XY      -   array[NPoints,NFeatures], dataset
//     NPoints -   number of points
//     NFeatures-  number of features
//
// Outputs:
//     XY      -   data are replaced by their within-row ranks;
//                 ranking starts from 0, ends at NFeatures-1
// ALGLIB: Copyright 18.04.2013 by Sergey Bochkanov
// API: void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
// API: void rankdatacentered(const real_2d_array &xy);
void rankdatacentered(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures) {
   ae_frame _frame_block;
   ae_int_t basecasecost;
   ae_frame_make(&_frame_block);
   NewObj(apbuffers, buf0);
   NewObj(apbuffers, buf1);
   NewObj(ae_shared_pool, pool);
   ae_assert(npoints >= 0, "RankData: NPoints < 0");
   ae_assert(nfeatures >= 1, "RankData: NFeatures<1");
   ae_assert(xy->rows >= npoints, "RankData: Rows(XY)<NPoints");
   ae_assert(xy->cols >= nfeatures || npoints == 0, "RankData: Cols(XY)<NFeatures");
   ae_assert(apservisfinitematrix(xy, npoints, nfeatures), "RankData: XY contains infinite/NAN elements");
// Basecase cost is a maximum cost of basecase problems.
// Problems harded than that cost will be split.
//
// Problem cost is assumed to be NPoints*NFeatures*log2(NFeatures),
// which is proportional, but NOT equal to number of FLOPs required
// to solve problem.
//
// Try to use serial code, no SMP functionality, no shared pools.
   basecasecost = 10000;
   if ((double)npoints * nfeatures * logbase2((double)nfeatures) < (double)basecasecost) {
      basestat_rankdatabasecase(xy, 0, npoints, nfeatures, true, &buf0, &buf1);
      ae_frame_leave();
      return;
   }
// Parallel code
   ae_shared_pool_set_seed(&pool, &buf0, sizeof(buf0), apbuffers_init, apbuffers_copy, apbuffers_free);
   basestat_rankdatarec(xy, 0, npoints, nfeatures, true, &pool, basecasecost);
   ae_frame_leave();
}

// 2-sample covariance
//
// Inputs:
//     X       -   sample 1 (array indexes: [0..N-1])
//     Y       -   sample 2 (array indexes: [0..N-1])
//     N       -   N >= 0, sample size:
//                 * if given, only N leading elements of X/Y are processed
//                 * if not given, automatically determined from input sizes
//
// Result:
//     covariance (zero for N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
// API: double cov2(const real_1d_array &x, const real_1d_array &y);
double cov2(RVector *x, RVector *y, ae_int_t n) {
   ae_int_t i;
   double xmean;
   double ymean;
   double v;
   double x0;
   double y0;
   double s;
   bool samex;
   bool samey;
   double result;
   ae_assert(n >= 0, "Cov2: N<0");
   ae_assert(x->cnt >= n, "Cov2: Length(X)<N!");
   ae_assert(y->cnt >= n, "Cov2: Length(Y)<N!");
   ae_assert(isfinitevector(x, n), "Cov2: X is not finite vector");
   ae_assert(isfinitevector(y, n), "Cov2: Y is not finite vector");
// Special case
   if (n <= 1) {
      result = 0.0;
      return result;
   }
// Calculate mean.
//
// Additionally we calculate SameX and SameY -
// flag variables which are set to True when
// all X[] (or Y[]) contain exactly same value.
//
// If at least one of them is True, we return zero
// (othwerwise we risk to get nonzero covariation
// because of roundoff).
   xmean = 0.0;
   ymean = 0.0;
   samex = true;
   samey = true;
   x0 = x->xR[0];
   y0 = y->xR[0];
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      s = x->xR[i];
      samex = samex && s == x0;
      xmean += s * v;
      s = y->xR[i];
      samey = samey && s == y0;
      ymean += s * v;
   }
   if (samex || samey) {
      result = 0.0;
      return result;
   }
// covariance
   v = 1.0 / (n - 1);
   result = 0.0;
   for (i = 0; i < n; i++) {
      result += v * (x->xR[i] - xmean) * (y->xR[i] - ymean);
   }
   return result;
}

// Pearson product-moment correlation coefficient
//
// Inputs:
//     X       -   sample 1 (array indexes: [0..N-1])
//     Y       -   sample 2 (array indexes: [0..N-1])
//     N       -   N >= 0, sample size:
//                 * if given, only N leading elements of X/Y are processed
//                 * if not given, automatically determined from input sizes
//
// Result:
//     Pearson product-moment correlation coefficient
//     (zero for N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
// API: double pearsoncorr2(const real_1d_array &x, const real_1d_array &y);
double pearsoncorr2(RVector *x, RVector *y, ae_int_t n) {
   ae_int_t i;
   double xmean;
   double ymean;
   double v;
   double x0;
   double y0;
   double s;
   bool samex;
   bool samey;
   double xv;
   double yv;
   double t1;
   double t2;
   double result;
   ae_assert(n >= 0, "PearsonCorr2: N<0");
   ae_assert(x->cnt >= n, "PearsonCorr2: Length(X)<N!");
   ae_assert(y->cnt >= n, "PearsonCorr2: Length(Y)<N!");
   ae_assert(isfinitevector(x, n), "PearsonCorr2: X is not finite vector");
   ae_assert(isfinitevector(y, n), "PearsonCorr2: Y is not finite vector");
// Special case
   if (n <= 1) {
      result = 0.0;
      return result;
   }
// Calculate mean.
//
// Additionally we calculate SameX and SameY -
// flag variables which are set to True when
// all X[] (or Y[]) contain exactly same value.
//
// If at least one of them is True, we return zero
// (othwerwise we risk to get nonzero correlation
// because of roundoff).
   xmean = 0.0;
   ymean = 0.0;
   samex = true;
   samey = true;
   x0 = x->xR[0];
   y0 = y->xR[0];
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      s = x->xR[i];
      samex = samex && s == x0;
      xmean += s * v;
      s = y->xR[i];
      samey = samey && s == y0;
      ymean += s * v;
   }
   if (samex || samey) {
      result = 0.0;
      return result;
   }
// numerator and denominator
   s = 0.0;
   xv = 0.0;
   yv = 0.0;
   for (i = 0; i < n; i++) {
      t1 = x->xR[i] - xmean;
      t2 = y->xR[i] - ymean;
      xv += sqr(t1);
      yv += sqr(t2);
      s += t1 * t2;
   }
   if (xv == 0.0 || yv == 0.0) {
      result = 0.0;
   } else {
      result = s / (sqrt(xv) * sqrt(yv));
   }
   return result;
}

// Spearman's rank correlation coefficient
//
// Inputs:
//     X       -   sample 1 (array indexes: [0..N-1])
//     Y       -   sample 2 (array indexes: [0..N-1])
//     N       -   N >= 0, sample size:
//                 * if given, only N leading elements of X/Y are processed
//                 * if not given, automatically determined from input sizes
//
// Result:
//     Spearman's rank correlation coefficient
//     (zero for N=0 or N=1)
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
// API: double spearmancorr2(const real_1d_array &x, const real_1d_array &y);
double spearmancorr2(RVector *x, RVector *y, ae_int_t n) {
   ae_frame _frame_block;
   double result;
   ae_frame_make(&_frame_block);
   DupVector(x);
   DupVector(y);
   NewObj(apbuffers, buf);
   ae_assert(n >= 0, "SpearmanCorr2: N<0");
   ae_assert(x->cnt >= n, "SpearmanCorr2: Length(X)<N!");
   ae_assert(y->cnt >= n, "SpearmanCorr2: Length(Y)<N!");
   ae_assert(isfinitevector(x, n), "SpearmanCorr2: X is not finite vector");
   ae_assert(isfinitevector(y, n), "SpearmanCorr2: Y is not finite vector");
// Special case
   if (n <= 1) {
      result = 0.0;
      ae_frame_leave();
      return result;
   }
   rankx(x, n, false, &buf);
   rankx(y, n, false, &buf);
   result = pearsoncorr2(x, y, n);
   ae_frame_leave();
   return result;
}

// Obsolete function, we recommend to use PearsonCorr2().
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double pearsoncorrelation(RVector *x, RVector *y, ae_int_t n) {
   double result;
   result = pearsoncorr2(x, y, n);
   return result;
}

// Obsolete function, we recommend to use SpearmanCorr2().
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double spearmanrankcorrelation(RVector *x, RVector *y, ae_int_t n) {
   double result;
   result = spearmancorr2(x, y, n);
   return result;
}

// Covariance matrix
//
// Inputs:
//     X   -   array[N,M], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     N   -   N >= 0, number of observations:
//             * if given, only leading N rows of X are used
//             * if not given, automatically determined from input size
//     M   -   M > 0, number of variables:
//             * if given, only leading M columns of X are used
//             * if not given, automatically determined from input size
//
// Outputs:
//     C   -   array[M,M], covariance matrix (zero if N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
// API: void covm(const real_2d_array &x, real_2d_array &c);
void covm(RMatrix *x, ae_int_t n, ae_int_t m, RMatrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_frame_make(&_frame_block);
   DupMatrix(x);
   SetMatrix(c);
   NewVector(t, 0, DT_REAL);
   NewVector(x0, 0, DT_REAL);
   NewVector(same, 0, DT_BOOL);
   ae_assert(n >= 0, "CovM: N<0");
   ae_assert(m >= 1, "CovM: M < 1");
   ae_assert(x->rows >= n, "CovM: Rows(X)<N!");
   ae_assert(x->cols >= m || n == 0, "CovM: Cols(X)<M!");
   ae_assert(apservisfinitematrix(x, n, m), "CovM: X contains infinite/NAN elements");
// N <= 1, return zero
   if (n <= 1) {
      ae_matrix_set_length(c, m, m);
      for (i = 0; i < m; i++) {
         for (j = 0; j < m; j++) {
            c->xyR[i][j] = 0.0;
         }
      }
      ae_frame_leave();
      return;
   }
// Calculate means,
// check for constant columns
   ae_vector_set_length(&t, m);
   ae_vector_set_length(&x0, m);
   ae_vector_set_length(&same, m);
   ae_matrix_set_length(c, m, m);
   for (i = 0; i < m; i++) {
      t.xR[i] = 0.0;
      same.xB[i] = true;
   }
   ae_v_move(x0.xR, 1, x->xyR[0], 1, m);
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      ae_v_addd(t.xR, 1, x->xyR[i], 1, m, v);
      for (j = 0; j < m; j++) {
         same.xB[j] = same.xB[j] && x->xyR[i][j] == x0.xR[j];
      }
   }
// * center variables;
// * if we have constant columns, these columns are
//   artificially zeroed (they must be zero in exact arithmetics,
//   but unfortunately floating point ops are not exact).
// * calculate upper half of symmetric covariance matrix
   for (i = 0; i < n; i++) {
      ae_v_sub(x->xyR[i], 1, t.xR, 1, m);
      for (j = 0; j < m; j++) {
         if (same.xB[j]) {
            x->xyR[i][j] = 0.0;
         }
      }
   }
   rmatrixsyrk(m, n, 1.0 / (n - 1), x, 0, 0, 1, 0.0, c, 0, 0, true);
   rmatrixenforcesymmetricity(c, m, true);
   ae_frame_leave();
}

// Pearson product-moment correlation matrix
//
// Inputs:
//     X   -   array[N,M], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     N   -   N >= 0, number of observations:
//             * if given, only leading N rows of X are used
//             * if not given, automatically determined from input size
//     M   -   M > 0, number of variables:
//             * if given, only leading M columns of X are used
//             * if not given, automatically determined from input size
//
// Outputs:
//     C   -   array[M,M], correlation matrix (zero if N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
// API: void pearsoncorrm(const real_2d_array &x, real_2d_array &c);
void pearsoncorrm(RMatrix *x, ae_int_t n, ae_int_t m, RMatrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_frame_make(&_frame_block);
   SetMatrix(c);
   NewVector(t, 0, DT_REAL);
   ae_assert(n >= 0, "PearsonCorrM: N<0");
   ae_assert(m >= 1, "PearsonCorrM: M < 1");
   ae_assert(x->rows >= n, "PearsonCorrM: Rows(X)<N!");
   ae_assert(x->cols >= m || n == 0, "PearsonCorrM: Cols(X)<M!");
   ae_assert(apservisfinitematrix(x, n, m), "PearsonCorrM: X contains infinite/NAN elements");
   ae_vector_set_length(&t, m);
   covm(x, n, m, c);
   for (i = 0; i < m; i++) {
      if (c->xyR[i][i] > 0.0) {
         t.xR[i] = 1 / sqrt(c->xyR[i][i]);
      } else {
         t.xR[i] = 0.0;
      }
   }
   for (i = 0; i < m; i++) {
      v = t.xR[i];
      for (j = 0; j < m; j++) {
         c->xyR[i][j] *= v * t.xR[j];
      }
   }
   ae_frame_leave();
}

// Spearman's rank correlation matrix
//
// Inputs:
//     X   -   array[N,M], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     N   -   N >= 0, number of observations:
//             * if given, only leading N rows of X are used
//             * if not given, automatically determined from input size
//     M   -   M > 0, number of variables:
//             * if given, only leading M columns of X are used
//             * if not given, automatically determined from input size
//
// Outputs:
//     C   -   array[M,M], correlation matrix (zero if N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
// API: void spearmancorrm(const real_2d_array &x, real_2d_array &c);
void spearmancorrm(RMatrix *x, ae_int_t n, ae_int_t m, RMatrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double vv;
   double x0;
   bool b;
   ae_frame_make(&_frame_block);
   SetMatrix(c);
   NewObj(apbuffers, buf);
   NewMatrix(xc, 0, 0, DT_REAL);
   NewVector(t, 0, DT_REAL);
   ae_assert(n >= 0, "SpearmanCorrM: N<0");
   ae_assert(m >= 1, "SpearmanCorrM: M < 1");
   ae_assert(x->rows >= n, "SpearmanCorrM: Rows(X)<N!");
   ae_assert(x->cols >= m || n == 0, "SpearmanCorrM: Cols(X)<M!");
   ae_assert(apservisfinitematrix(x, n, m), "SpearmanCorrM: X contains infinite/NAN elements");
// N <= 1, return zero
   if (n <= 1) {
      ae_matrix_set_length(c, m, m);
      for (i = 0; i < m; i++) {
         for (j = 0; j < m; j++) {
            c->xyR[i][j] = 0.0;
         }
      }
      ae_frame_leave();
      return;
   }
// Allocate
   ae_vector_set_length(&t, imax2(n, m));
   ae_matrix_set_length(c, m, m);
// Replace data with ranks
   ae_matrix_set_length(&xc, m, n);
   rmatrixtranspose(n, m, x, 0, 0, &xc, 0, 0);
   rankdata(&xc, m, n);
// 1. Calculate means, check for constant columns
// 2. Center variables, constant  columns are
//   artificialy zeroed (they must be zero in exact arithmetics,
//   but unfortunately floating point is not exact).
   for (i = 0; i < m; i++) {
   // Calculate:
   // * V - mean value of I-th variable
   // * B - True in case all variable values are same
      v = 0.0;
      b = true;
      x0 = xc.xyR[i][0];
      for (j = 0; j < n; j++) {
         vv = xc.xyR[i][j];
         v += vv;
         b = b && vv == x0;
      }
      v /= n;
   // Center/zero I-th variable
      if (b) {
      // Zero
         for (j = 0; j < n; j++) {
            xc.xyR[i][j] = 0.0;
         }
      } else {
      // Center
         for (j = 0; j < n; j++) {
            xc.xyR[i][j] -= v;
         }
      }
   }
// Calculate upper half of symmetric covariance matrix
   rmatrixsyrk(m, n, 1.0 / (n - 1), &xc, 0, 0, 0, 0.0, c, 0, 0, true);
// Calculate Pearson coefficients (upper triangle)
   for (i = 0; i < m; i++) {
      if (c->xyR[i][i] > 0.0) {
         t.xR[i] = 1 / sqrt(c->xyR[i][i]);
      } else {
         t.xR[i] = 0.0;
      }
   }
   for (i = 0; i < m; i++) {
      v = t.xR[i];
      for (j = i; j < m; j++) {
         c->xyR[i][j] *= v * t.xR[j];
      }
   }
// force symmetricity
   rmatrixenforcesymmetricity(c, m, true);
   ae_frame_leave();
}

// Cross-covariance matrix
//
// Inputs:
//     X   -   array[N,M1], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     Y   -   array[N,M2], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     N   -   N >= 0, number of observations:
//             * if given, only leading N rows of X/Y are used
//             * if not given, automatically determined from input sizes
//     M1  -   M1 > 0, number of variables in X:
//             * if given, only leading M1 columns of X are used
//             * if not given, automatically determined from input size
//     M2  -   M2 > 0, number of variables in Y:
//             * if given, only leading M1 columns of X are used
//             * if not given, automatically determined from input size
//
// Outputs:
//     C   -   array[M1,M2], cross-covariance matrix (zero if N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
// API: void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
void covm2(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_frame_make(&_frame_block);
   DupMatrix(x);
   DupMatrix(y);
   SetMatrix(c);
   NewVector(t, 0, DT_REAL);
   NewVector(x0, 0, DT_REAL);
   NewVector(y0, 0, DT_REAL);
   NewVector(samex, 0, DT_BOOL);
   NewVector(samey, 0, DT_BOOL);
   ae_assert(n >= 0, "CovM2: N<0");
   ae_assert(m1 >= 1, "CovM2: M1<1");
   ae_assert(m2 >= 1, "CovM2: M2<1");
   ae_assert(x->rows >= n, "CovM2: Rows(X)<N!");
   ae_assert(x->cols >= m1 || n == 0, "CovM2: Cols(X)<M1!");
   ae_assert(apservisfinitematrix(x, n, m1), "CovM2: X contains infinite/NAN elements");
   ae_assert(y->rows >= n, "CovM2: Rows(Y)<N!");
   ae_assert(y->cols >= m2 || n == 0, "CovM2: Cols(Y)<M2!");
   ae_assert(apservisfinitematrix(y, n, m2), "CovM2: X contains infinite/NAN elements");
// N <= 1, return zero
   if (n <= 1) {
      ae_matrix_set_length(c, m1, m2);
      for (i = 0; i < m1; i++) {
         for (j = 0; j < m2; j++) {
            c->xyR[i][j] = 0.0;
         }
      }
      ae_frame_leave();
      return;
   }
// Allocate
   ae_vector_set_length(&t, imax2(m1, m2));
   ae_vector_set_length(&x0, m1);
   ae_vector_set_length(&y0, m2);
   ae_vector_set_length(&samex, m1);
   ae_vector_set_length(&samey, m2);
   ae_matrix_set_length(c, m1, m2);
// * calculate means of X
// * center X
// * if we have constant columns, these columns are
//   artificially zeroed (they must be zero in exact arithmetics,
//   but unfortunately floating point ops are not exact).
   for (i = 0; i < m1; i++) {
      t.xR[i] = 0.0;
      samex.xB[i] = true;
   }
   ae_v_move(x0.xR, 1, x->xyR[0], 1, m1);
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      ae_v_addd(t.xR, 1, x->xyR[i], 1, m1, v);
      for (j = 0; j < m1; j++) {
         samex.xB[j] = samex.xB[j] && x->xyR[i][j] == x0.xR[j];
      }
   }
   for (i = 0; i < n; i++) {
      ae_v_sub(x->xyR[i], 1, t.xR, 1, m1);
      for (j = 0; j < m1; j++) {
         if (samex.xB[j]) {
            x->xyR[i][j] = 0.0;
         }
      }
   }
// Repeat same steps for Y
   for (i = 0; i < m2; i++) {
      t.xR[i] = 0.0;
      samey.xB[i] = true;
   }
   ae_v_move(y0.xR, 1, y->xyR[0], 1, m2);
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      ae_v_addd(t.xR, 1, y->xyR[i], 1, m2, v);
      for (j = 0; j < m2; j++) {
         samey.xB[j] = samey.xB[j] && y->xyR[i][j] == y0.xR[j];
      }
   }
   for (i = 0; i < n; i++) {
      ae_v_sub(y->xyR[i], 1, t.xR, 1, m2);
      for (j = 0; j < m2; j++) {
         if (samey.xB[j]) {
            y->xyR[i][j] = 0.0;
         }
      }
   }
// calculate cross-covariance matrix
   rmatrixgemm(m1, m2, n, 1.0 / (n - 1), x, 0, 0, 1, y, 0, 0, 0, 0.0, c, 0, 0);
   ae_frame_leave();
}

// Pearson product-moment cross-correlation matrix
//
// Inputs:
//     X   -   array[N,M1], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     Y   -   array[N,M2], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     N   -   N >= 0, number of observations:
//             * if given, only leading N rows of X/Y are used
//             * if not given, automatically determined from input sizes
//     M1  -   M1 > 0, number of variables in X:
//             * if given, only leading M1 columns of X are used
//             * if not given, automatically determined from input size
//     M2  -   M2 > 0, number of variables in Y:
//             * if given, only leading M1 columns of X are used
//             * if not given, automatically determined from input size
//
// Outputs:
//     C   -   array[M1,M2], cross-correlation matrix (zero if N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
// API: void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
void pearsoncorrm2(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_frame_make(&_frame_block);
   DupMatrix(x);
   DupMatrix(y);
   SetMatrix(c);
   NewVector(t, 0, DT_REAL);
   NewVector(x0, 0, DT_REAL);
   NewVector(y0, 0, DT_REAL);
   NewVector(sx, 0, DT_REAL);
   NewVector(sy, 0, DT_REAL);
   NewVector(samex, 0, DT_BOOL);
   NewVector(samey, 0, DT_BOOL);
   ae_assert(n >= 0, "PearsonCorrM2: N<0");
   ae_assert(m1 >= 1, "PearsonCorrM2: M1<1");
   ae_assert(m2 >= 1, "PearsonCorrM2: M2<1");
   ae_assert(x->rows >= n, "PearsonCorrM2: Rows(X)<N!");
   ae_assert(x->cols >= m1 || n == 0, "PearsonCorrM2: Cols(X)<M1!");
   ae_assert(apservisfinitematrix(x, n, m1), "PearsonCorrM2: X contains infinite/NAN elements");
   ae_assert(y->rows >= n, "PearsonCorrM2: Rows(Y)<N!");
   ae_assert(y->cols >= m2 || n == 0, "PearsonCorrM2: Cols(Y)<M2!");
   ae_assert(apservisfinitematrix(y, n, m2), "PearsonCorrM2: X contains infinite/NAN elements");
// N <= 1, return zero
   if (n <= 1) {
      ae_matrix_set_length(c, m1, m2);
      for (i = 0; i < m1; i++) {
         for (j = 0; j < m2; j++) {
            c->xyR[i][j] = 0.0;
         }
      }
      ae_frame_leave();
      return;
   }
// Allocate
   ae_vector_set_length(&t, imax2(m1, m2));
   ae_vector_set_length(&x0, m1);
   ae_vector_set_length(&y0, m2);
   ae_vector_set_length(&sx, m1);
   ae_vector_set_length(&sy, m2);
   ae_vector_set_length(&samex, m1);
   ae_vector_set_length(&samey, m2);
   ae_matrix_set_length(c, m1, m2);
// * calculate means of X
// * center X
// * if we have constant columns, these columns are
//   artificially zeroed (they must be zero in exact arithmetics,
//   but unfortunately floating point ops are not exact).
// * calculate column variances
   for (i = 0; i < m1; i++) {
      t.xR[i] = 0.0;
      samex.xB[i] = true;
      sx.xR[i] = 0.0;
   }
   ae_v_move(x0.xR, 1, x->xyR[0], 1, m1);
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      ae_v_addd(t.xR, 1, x->xyR[i], 1, m1, v);
      for (j = 0; j < m1; j++) {
         samex.xB[j] = samex.xB[j] && x->xyR[i][j] == x0.xR[j];
      }
   }
   for (i = 0; i < n; i++) {
      ae_v_sub(x->xyR[i], 1, t.xR, 1, m1);
      for (j = 0; j < m1; j++) {
         if (samex.xB[j]) {
            x->xyR[i][j] = 0.0;
         }
         sx.xR[j] += x->xyR[i][j] * x->xyR[i][j];
      }
   }
   for (j = 0; j < m1; j++) {
      sx.xR[j] = sqrt(sx.xR[j] / (n - 1));
   }
// Repeat same steps for Y
   for (i = 0; i < m2; i++) {
      t.xR[i] = 0.0;
      samey.xB[i] = true;
      sy.xR[i] = 0.0;
   }
   ae_v_move(y0.xR, 1, y->xyR[0], 1, m2);
   v = 1.0 / n;
   for (i = 0; i < n; i++) {
      ae_v_addd(t.xR, 1, y->xyR[i], 1, m2, v);
      for (j = 0; j < m2; j++) {
         samey.xB[j] = samey.xB[j] && y->xyR[i][j] == y0.xR[j];
      }
   }
   for (i = 0; i < n; i++) {
      ae_v_sub(y->xyR[i], 1, t.xR, 1, m2);
      for (j = 0; j < m2; j++) {
         if (samey.xB[j]) {
            y->xyR[i][j] = 0.0;
         }
         sy.xR[j] += y->xyR[i][j] * y->xyR[i][j];
      }
   }
   for (j = 0; j < m2; j++) {
      sy.xR[j] = sqrt(sy.xR[j] / (n - 1));
   }
// calculate cross-covariance matrix
   rmatrixgemm(m1, m2, n, 1.0 / (n - 1), x, 0, 0, 1, y, 0, 0, 0, 0.0, c, 0, 0);
// Divide by standard deviations
   for (i = 0; i < m1; i++) {
      if (sx.xR[i] != 0.0) {
         sx.xR[i] = 1 / sx.xR[i];
      } else {
         sx.xR[i] = 0.0;
      }
   }
   for (i = 0; i < m2; i++) {
      if (sy.xR[i] != 0.0) {
         sy.xR[i] = 1 / sy.xR[i];
      } else {
         sy.xR[i] = 0.0;
      }
   }
   for (i = 0; i < m1; i++) {
      v = sx.xR[i];
      for (j = 0; j < m2; j++) {
         c->xyR[i][j] *= v * sy.xR[j];
      }
   }
   ae_frame_leave();
}

// Spearman's rank cross-correlation matrix
//
// Inputs:
//     X   -   array[N,M1], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     Y   -   array[N,M2], sample matrix:
//             * J-th column corresponds to J-th variable
//             * I-th row corresponds to I-th observation
//     N   -   N >= 0, number of observations:
//             * if given, only leading N rows of X/Y are used
//             * if not given, automatically determined from input sizes
//     M1  -   M1 > 0, number of variables in X:
//             * if given, only leading M1 columns of X are used
//             * if not given, automatically determined from input size
//     M2  -   M2 > 0, number of variables in Y:
//             * if given, only leading M1 columns of X are used
//             * if not given, automatically determined from input size
//
// Outputs:
//     C   -   array[M1,M2], cross-correlation matrix (zero if N=0 or N=1)
// ALGLIB: Copyright 28.10.2010 by Sergey Bochkanov
// API: void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
// API: void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
void spearmancorrm2(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double v2;
   double vv;
   bool b;
   double x0;
   double y0;
   ae_frame_make(&_frame_block);
   SetMatrix(c);
   NewVector(t, 0, DT_REAL);
   NewVector(sx, 0, DT_REAL);
   NewVector(sy, 0, DT_REAL);
   NewMatrix(xc, 0, 0, DT_REAL);
   NewMatrix(yc, 0, 0, DT_REAL);
   NewObj(apbuffers, buf);
   ae_assert(n >= 0, "SpearmanCorrM2: N<0");
   ae_assert(m1 >= 1, "SpearmanCorrM2: M1<1");
   ae_assert(m2 >= 1, "SpearmanCorrM2: M2<1");
   ae_assert(x->rows >= n, "SpearmanCorrM2: Rows(X)<N!");
   ae_assert(x->cols >= m1 || n == 0, "SpearmanCorrM2: Cols(X)<M1!");
   ae_assert(apservisfinitematrix(x, n, m1), "SpearmanCorrM2: X contains infinite/NAN elements");
   ae_assert(y->rows >= n, "SpearmanCorrM2: Rows(Y)<N!");
   ae_assert(y->cols >= m2 || n == 0, "SpearmanCorrM2: Cols(Y)<M2!");
   ae_assert(apservisfinitematrix(y, n, m2), "SpearmanCorrM2: X contains infinite/NAN elements");
// N <= 1, return zero
   if (n <= 1) {
      ae_matrix_set_length(c, m1, m2);
      for (i = 0; i < m1; i++) {
         for (j = 0; j < m2; j++) {
            c->xyR[i][j] = 0.0;
         }
      }
      ae_frame_leave();
      return;
   }
// Allocate
   ae_vector_set_length(&t, imax2(imax2(m1, m2), n));
   ae_vector_set_length(&sx, m1);
   ae_vector_set_length(&sy, m2);
   ae_matrix_set_length(c, m1, m2);
// Replace data with ranks
   ae_matrix_set_length(&xc, m1, n);
   ae_matrix_set_length(&yc, m2, n);
   rmatrixtranspose(n, m1, x, 0, 0, &xc, 0, 0);
   rmatrixtranspose(n, m2, y, 0, 0, &yc, 0, 0);
   rankdata(&xc, m1, n);
   rankdata(&yc, m2, n);
// 1. Calculate means, variances, check for constant columns
// 2. Center variables, constant  columns are
//   artificialy zeroed (they must be zero in exact arithmetics,
//   but unfortunately floating point is not exact).
//
// Description of variables:
// * V - mean value of I-th variable
// * V2- variance
// * VV-temporary
// * B - True in case all variable values are same
   for (i = 0; i < m1; i++) {
      v = 0.0;
      v2 = 0.0;
      b = true;
      x0 = xc.xyR[i][0];
      for (j = 0; j < n; j++) {
         vv = xc.xyR[i][j];
         v += vv;
         b = b && vv == x0;
      }
      v /= n;
      if (b) {
         for (j = 0; j < n; j++) {
            xc.xyR[i][j] = 0.0;
         }
      } else {
         for (j = 0; j < n; j++) {
            vv = xc.xyR[i][j];
            xc.xyR[i][j] = vv - v;
            v2 += (vv - v) * (vv - v);
         }
      }
      sx.xR[i] = sqrt(v2 / (n - 1));
   }
   for (i = 0; i < m2; i++) {
      v = 0.0;
      v2 = 0.0;
      b = true;
      y0 = yc.xyR[i][0];
      for (j = 0; j < n; j++) {
         vv = yc.xyR[i][j];
         v += vv;
         b = b && vv == y0;
      }
      v /= n;
      if (b) {
         for (j = 0; j < n; j++) {
            yc.xyR[i][j] = 0.0;
         }
      } else {
         for (j = 0; j < n; j++) {
            vv = yc.xyR[i][j];
            yc.xyR[i][j] = vv - v;
            v2 += (vv - v) * (vv - v);
         }
      }
      sy.xR[i] = sqrt(v2 / (n - 1));
   }
// calculate cross-covariance matrix
   rmatrixgemm(m1, m2, n, 1.0 / (n - 1), &xc, 0, 0, 0, &yc, 0, 0, 1, 0.0, c, 0, 0);
// Divide by standard deviations
   for (i = 0; i < m1; i++) {
      if (sx.xR[i] != 0.0) {
         sx.xR[i] = 1 / sx.xR[i];
      } else {
         sx.xR[i] = 0.0;
      }
   }
   for (i = 0; i < m2; i++) {
      if (sy.xR[i] != 0.0) {
         sy.xR[i] = 1 / sy.xR[i];
      } else {
         sy.xR[i] = 0.0;
      }
   }
   for (i = 0; i < m1; i++) {
      v = sx.xR[i];
      for (j = 0; j < m2; j++) {
         c->xyR[i][j] *= v * sy.xR[j];
      }
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::samplemoments(ConstT(ae_vector, x), n, &mean, &variance, &skewness, &kurtosis);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::samplemoments(ConstT(ae_vector, x), n, &mean, &variance, &skewness, &kurtosis);
   alglib_impl::ae_state_clear();
}
#endif

double samplemean(const real_1d_array &x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::samplemean(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double samplemean(const real_1d_array &x) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::samplemean(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double samplevariance(const real_1d_array &x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::samplevariance(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double samplevariance(const real_1d_array &x) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::samplevariance(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double sampleskewness(const real_1d_array &x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::sampleskewness(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double sampleskewness(const real_1d_array &x) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::sampleskewness(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double samplekurtosis(const real_1d_array &x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::samplekurtosis(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double samplekurtosis(const real_1d_array &x) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::samplekurtosis(ConstT(ae_vector, x), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sampleadev(ConstT(ae_vector, x), n, &adev);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void sampleadev(const real_1d_array &x, double &adev) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sampleadev(ConstT(ae_vector, x), n, &adev);
   alglib_impl::ae_state_clear();
}
#endif

void samplemedian(const real_1d_array &x, const ae_int_t n, double &median) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::samplemedian(ConstT(ae_vector, x), n, &median);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void samplemedian(const real_1d_array &x, double &median) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::samplemedian(ConstT(ae_vector, x), n, &median);
   alglib_impl::ae_state_clear();
}
#endif

void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::samplepercentile(ConstT(ae_vector, x), n, p, &v);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void samplepercentile(const real_1d_array &x, const double p, double &v) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::samplepercentile(ConstT(ae_vector, x), n, p, &v);
   alglib_impl::ae_state_clear();
}
#endif

void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rankdata(ConstT(ae_matrix, xy), npoints, nfeatures);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void rankdata(const real_2d_array &xy) {
   ae_int_t npoints = xy.rows();
   ae_int_t nfeatures = xy.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rankdata(ConstT(ae_matrix, xy), npoints, nfeatures);
   alglib_impl::ae_state_clear();
}
#endif

void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rankdatacentered(ConstT(ae_matrix, xy), npoints, nfeatures);
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void rankdatacentered(const real_2d_array &xy) {
   ae_int_t npoints = xy.rows();
   ae_int_t nfeatures = xy.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rankdatacentered(ConstT(ae_matrix, xy), npoints, nfeatures);
   alglib_impl::ae_state_clear();
}
#endif

double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cov2(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double cov2(const real_1d_array &x, const real_1d_array &y) {
   if (x.length() != y.length()) ThrowError("Error while calling 'cov2': looks like one of arguments has wrong size");
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cov2(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::pearsoncorr2(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y) {
   if (x.length() != y.length()) ThrowError("Error while calling 'pearsoncorr2': looks like one of arguments has wrong size");
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::pearsoncorr2(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spearmancorr2(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double spearmancorr2(const real_1d_array &x, const real_1d_array &y) {
   if (x.length() != y.length()) ThrowError("Error while calling 'spearmancorr2': looks like one of arguments has wrong size");
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spearmancorr2(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::pearsoncorrelation(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}

double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spearmanrankcorrelation(ConstT(ae_vector, x), ConstT(ae_vector, y), n);
   alglib_impl::ae_state_clear();
   return D;
}

void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::covm(ConstT(ae_matrix, x), n, m, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void covm(const real_2d_array &x, real_2d_array &c) {
   ae_int_t n = x.rows();
   ae_int_t m = x.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::covm(ConstT(ae_matrix, x), n, m, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#endif

void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::pearsoncorrm(ConstT(ae_matrix, x), n, m, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void pearsoncorrm(const real_2d_array &x, real_2d_array &c) {
   ae_int_t n = x.rows();
   ae_int_t m = x.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::pearsoncorrm(ConstT(ae_matrix, x), n, m, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#endif

void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spearmancorrm(ConstT(ae_matrix, x), n, m, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void spearmancorrm(const real_2d_array &x, real_2d_array &c) {
   ae_int_t n = x.rows();
   ae_int_t m = x.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spearmancorrm(ConstT(ae_matrix, x), n, m, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#endif

void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::covm2(ConstT(ae_matrix, x), ConstT(ae_matrix, y), n, m1, m2, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c) {
   if (x.rows() != y.rows()) ThrowError("Error while calling 'covm2': looks like one of arguments has wrong size");
   ae_int_t n = x.rows();
   ae_int_t m1 = x.cols();
   ae_int_t m2 = y.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::covm2(ConstT(ae_matrix, x), ConstT(ae_matrix, y), n, m1, m2, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#endif

void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::pearsoncorrm2(ConstT(ae_matrix, x), ConstT(ae_matrix, y), n, m1, m2, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c) {
   if (x.rows() != y.rows()) ThrowError("Error while calling 'pearsoncorrm2': looks like one of arguments has wrong size");
   ae_int_t n = x.rows();
   ae_int_t m1 = x.cols();
   ae_int_t m2 = y.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::pearsoncorrm2(ConstT(ae_matrix, x), ConstT(ae_matrix, y), n, m1, m2, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#endif

void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spearmancorrm2(ConstT(ae_matrix, x), ConstT(ae_matrix, y), n, m1, m2, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c) {
   if (x.rows() != y.rows()) ThrowError("Error while calling 'spearmancorrm2': looks like one of arguments has wrong size");
   ae_int_t n = x.rows();
   ae_int_t m1 = x.cols();
   ae_int_t m2 = y.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spearmancorrm2(ConstT(ae_matrix, x), ConstT(ae_matrix, y), n, m1, m2, ConstT(ae_matrix, c));
   alglib_impl::ae_state_clear();
}
#endif
} // end of namespace alglib

// === CORRELATIONTESTS Package ===
// Depends on: (SpecialFunctions) STUDENTTDISTR
// Depends on: BASESTAT
namespace alglib_impl {
// Pearson's correlation coefficient significance test
//
// This test checks hypotheses about whether X  and  Y  are  samples  of  two
// continuous  distributions  having  zero  correlation  or   whether   their
// correlation is non-zero.
//
// The following tests are performed:
//     * two-tailed test (null hypothesis - X and Y have zero correlation)
//     * left-tailed test (null hypothesis - the correlation  coefficient  is
//       greater than or equal to 0)
//     * right-tailed test (null hypothesis - the correlation coefficient  is
//       less than or equal to 0).
//
// Requirements:
//     * the number of elements in each sample is not less than 5
//     * normality of distributions of X and Y.
//
// Inputs:
//     R   -   Pearson's correlation coefficient for X and Y
//     N   -   number of elements in samples, N >= 5.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);
void pearsoncorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail) {
   double t;
   double p;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
// Some special cases
   if (r >= 1.0) {
      *bothtails = 0.0;
      *lefttail = 1.0;
      *righttail = 0.0;
      return;
   }
   if (r <= -1.0) {
      *bothtails = 0.0;
      *lefttail = 0.0;
      *righttail = 1.0;
      return;
   }
   if (n < 5) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// General case
   t = r * sqrt((n - 2) / (1 - sqr(r)));
   p = studenttdistribution(n - 2, t);
   *bothtails = 2 * rmin2(p, 1 - p);
   *lefttail = p;
   *righttail = 1 - p;
}

// Tail(T,N), accepts T < 0
static double correlationtests_spearmantail(double t, ae_int_t n) {
   switch (n) {
      case 5:
         return
            t <= -7.584e-01 ? (
               t <= -1.704e+00 ? (
                  t <= -3.580e+00 ? 8.304e-03 : t <= -2.322e+00 ? 4.163e-02: 6.641e-02
               ) : t <= -1.303e+00 ? 1.164e-01 : t <= -1.003e+00 ? 1.748e-01 : 2.249e-01
            ) : t <= -1.759e-01 ? (
               t <= -5.468e-01 ? 2.581e-01 : t <= -3.555e-01 ? 3.413e-01 : 3.911e-01
            ) : t <= -1.741e-03 ? 4.747e-01 : t <= -0.000e+00 ? 5.248e-01 : studenttdistribution(3, t);
      case 6:
         return
            t <= -2.045e+00 ? (
               t <= -3.834e+00 ? (
                  t <= -5.663e+00 ? 1.366e-03: 8.350e-03
               ) : t <= -2.968e+00 ? 1.668e-02 : t <= -2.430e+00 ? 2.921e-02 : 5.144e-02
            ) : t <= -1.295e+00 ? (
               t <= -1.747e+00 ? 6.797e-02 : t <= -1.502e+00 ? 8.752e-02 : 1.210e-01
            ) : t <= -1.113e+00 ? 1.487e-01 : t <= -1.001e+00 ? 1.780e-01 : studenttdistribution(4, t);
      case 7:
         return
            t <= -2.068e+00 ? (
               t <= -3.728e+00 ? (
                  t <= -5.620e+00 ? (
                     t <= -8.159e+00 ? 2.081e-04 : 1.393e-03
                  ) : t <= -4.445e+00 ? 3.398e-03 : 6.187e-03
               ) : t <= -2.844e+00 ? (
                  t <= -3.226e+00 ? 1.200e-02: 1.712e-02
               ) : t <= -2.539e+00 ? 2.408e-02 : t <= -2.285e+00 ? 3.320e-02 : 4.406e-02
            ) : t <= -1.420e+00 ? (
               t <= -1.710e+00 ? (
                  t <= -1.879e+00 ? 5.478e-02 : 6.946e-02
               ) : t <= -1.559e+00 ? 8.331e-02 : 1.001e-01
            ) : t <= -1.173e+00 ? (
               t <= -1.292e+00 ? 1.180e-01: 1.335e-01
            ) : t <= -1.062e+00 ? 1.513e-01 : t <= -1.001e+00 ? 1.770e-01 : studenttdistribution(5, t);
      case 8:
         return
            t <= -3.381e+00 ? (
               t <= -5.213e+00 ? (
                  t <= -7.685e+00 ? (
                     t <= -1.103e+01 ? 2.194e-05 : 2.008e-04
                  ) : t <= -6.143e+00 ? 5.686e-04 : 1.138e-03
               ) : t <= -4.081e+00 ? (
                  t <= -4.567e+00 ? 2.310e-03 : 3.634e-03
               ) : t <= -3.697e+00 ? 5.369e-03 : 7.708e-03
            ) : t <= -2.502e+00 ? (
               t <= -2.884e+00 ? (
                  t <= -3.114e+00 ? 1.087e-02 : 1.397e-02
               ) : t <= -2.682e+00 ? 1.838e-02 : 2.288e-02
            ) : t <= -2.192e+00 ? (
               t <= -2.340e+00 ? 2.883e-02 : 3.469e-02
            ) : t <= -2.057e+00 ? 4.144e-02 : t <= -2.001e+00 ? 4.804e-02 : studenttdistribution(6, t);
     case 9:
         return
            t <= -3.336e+00 ? (
               t <= -4.991e+00 ? (
                  t <= -6.890e+00 ? (
                     t <= -9.989e+00 ? 2.306e-05 : t <= -8.069e+00 ? 8.167e-05 : 1.744e-04
                  ) : t <= -6.077e+00 ? 3.625e-04 : t <= -5.469e+00 ? 6.450e-04 : 1.001e-03
               ) : t <= -3.991e+00 ? (
                  t <= -4.600e+00 ? 1.514e-03 : t <= -4.272e+00 ? 2.213e-03 : 2.990e-03
               ) : t <= -3.746e+00 ? 4.101e-03 : t <= -3.530e+00 ? 5.355e-03 : 6.887e-03
            ) : t <= -2.477e+00 ? (
               t <= -2.855e+00 ? (
                  t <= -3.161e+00 ? 8.598e-03 : t <= -3.002e+00 ? 1.065e-02 : 1.268e-02
               ) : t <= -2.720e+00 ? 1.552e-02 : t <= -2.595e+00 ? 1.836e-02 : 2.158e-02
            ) : t <= -2.166e+00 ? (
               t <= -2.368e+00 ? 2.512e-02 : t <= -2.264e+00 ? 2.942e-02 : 3.325e-02
            ) : t <= -2.073e+00 ? 3.800e-02 : t <= -2.001e+00 ? 4.285e-02 : studenttdistribution(7, t);
      default: return studenttdistribution(n - 2, t);
   }
}

// Spearman's rank correlation coefficient significance test
//
// This test checks hypotheses about whether X  and  Y  are  samples  of  two
// continuous  distributions  having  zero  correlation  or   whether   their
// correlation is non-zero.
//
// The following tests are performed:
//     * two-tailed test (null hypothesis - X and Y have zero correlation)
//     * left-tailed test (null hypothesis - the correlation  coefficient  is
//       greater than or equal to 0)
//     * right-tailed test (null hypothesis - the correlation coefficient  is
//       less than or equal to 0).
//
// Requirements:
//     * the number of elements in each sample is not less than 5.
//
// The test is non-parametric and doesn't require distributions X and Y to be
// normal.
//
// Inputs:
//     R   -   Spearman's rank correlation coefficient for X and Y
//     N   -   number of elements in samples, N >= 5.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);
void spearmanrankcorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail) {
   double t;
   double p;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
// Special case
   if (n < 5) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// General case
   if (r >= 1.0) {
      t = 1.0E10;
   } else {
      if (r <= -1.0) {
         t = -1.0E10;
      } else {
         t = r * sqrt((n - 2) / (1 - sqr(r)));
      }
   }
   if (t < 0.0) {
      p = correlationtests_spearmantail(t, n);
      *bothtails = 2 * p;
      *lefttail = p;
      *righttail = 1 - p;
   } else {
      p = correlationtests_spearmantail(-t, n);
      *bothtails = 2 * p;
      *lefttail = 1 - p;
      *righttail = p;
   }
}
} // end of namespace alglib_impl

namespace alglib {
void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::pearsoncorrelationsignificance(r, n, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}

void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spearmanrankcorrelationsignificance(r, n, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === JARQUEBERA Package ===
namespace alglib_impl {
static void jarquebera_jbcheb(double x, double c, double *tj, double *tj1, double *r) {
   *r += c * (*tj);
   double t = 2.0 * x * (*tj1) - (*tj);
   *tj = *tj1, *tj1 = t;
}

static double jarquebera_jbtbl5(double s) {
   if (s <= 0.4000) {
      double x = 2 * (s - 0.000000) / 0.400000 - 1, x2 = x * x;
      const double c0 = -1.097885e-20, c1 = -2.854501e-20, c2 = -1.756616e-20;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 1.1000) {
      double x = 2 * (s - 0.400000) / 0.700000 - 1, x2 = x * x;
      const double c0 = -1.324545e+00, c1 = -1.075941e+00, c2 = -9.772272e-01, c3 = +3.175686e-01;
      const double c4 = -1.576162e-01, c5 = +1.126861e-01, c6 = -3.434425e-02, c7 = -2.790359e-01;
      const double c8 = +2.809178e-02, c9 = -5.479704e-01, ca = +3.717040e-02, cb = -5.294170e-01;
      const double cc = +2.880632e-02, cd = -3.023344e-01, ce = +1.601531e-02, cf = -7.920403e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8, tb = x2 * ta - t9;
      double tc = x2 * tb - ta, td = x2 * tc - tb, te = x2 * td - tc, tf = x2 * te - td;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta + cb * tb + cc * tc + cd * td + ce * te + cf * tf;
      return result > 0.0 ? 0.0 : result;
   } else return -5.188419e+02 * (s - 1.100000e+00) - 4.767297e+00;
}

static double jarquebera_jbtbl6(double s) {
   if (s <= 0.2500) {
      double x = 2 * (s - 0.000000) / 0.250000 - 1, x2 = x * x;
      const double c0 = -2.274707e-04, c1 = -5.700471e-04, c2 = -3.425764e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 1.3000) {
      double x = 2 * (s - 0.250000) / 1.050000 - 1, x2 = x * x;
      const double c0 = -1.339000e+00, c1 = -2.011104e+00, c2 = -8.168177e-01, c3 = -1.085666e-01;
      const double c4 = +7.738606e-02, c5 = +7.022876e-02, c6 = +3.462402e-02, c7 = +6.908270e-03;
      const double c8 = -8.230772e-03, c9 = -1.006996e-02, ca = -5.410222e-03, cb = -2.893768e-03;
      const double cc = +8.114564e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8, tb = x2 * ta - t9;
      double tc = x2 * tb - ta;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta + cb * tb + cc * tc;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 1.8500) {
      double x = 2 * (s - 1.300000) / 0.550000 - 1, x2 = x * x;
      const double c0 = -6.794311e+00, c1 = -3.578700e+00, c2 = -1.394664e+00, c3 = -7.928290e-01;
      const double c4 = -4.813273e-01, c5 = -3.076063e-01, c6 = -1.835380e-01, c7 = -1.013013e-01;
      const double c8 = -5.058903e-02, c9 = -1.856915e-02, ca = -6.710887e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -1.770029e+02 * (s - 1.850000e+00) - 1.371015e+01;
}

static double jarquebera_jbtbl7(double s) {
   if (s <= 1.4000) {
      double x = 2 * (s - 0.000000) / 1.400000 - 1, x2 = x * x;
      const double c0 = -1.093681e+00, c1 = -1.695911e+00, c2 = -7.473192e-01, c3 = -1.203236e-01;
      const double c4 = +6.590379e-02, c5 = +6.291876e-02, c6 = +3.132007e-02, c7 = +9.411147e-03;
      const double c8 = -1.180067e-03, c9 = -3.487610e-03, ca = -2.436561e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 3.0000) {
      double x = 2 * (s - 1.400000) / 1.600000 - 1, x2 = x * x;
      const double c0 = -5.947854e+00, c1 = -2.772675e+00, c2 = -4.707912e-01, c3 = -1.691171e-01;
      const double c4 = -4.132795e-02, c5 = -1.481310e-02, c6 = +2.867536e-03, c7 = +8.772327e-04;
      const double c8 = +5.033387e-03, c9 = -1.378277e-03, ca = -2.497964e-03, cb = -3.636814e-03;
      const double cc = -9.581640e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8, tb = x2 * ta - t9;
      double tc = x2 * tb - ta;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta + cb * tb + cc * tc;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 3.2000) {
      double x = 2 * (s - 3.000000) / 0.200000 - 1, x2 = x * x;
      const double c0 = -7.511008e+00, c1 = -8.140472e-01, c2 = +1.682053e+00, c3 = -2.568561e-02;
      const double c4 = -1.933930e+00, c5 = -8.140472e-01, c6 = -3.895025e+00, c7 = -8.140472e-01;
      const double c8 = -1.933930e+00, c9 = -2.568561e-02, ca = +1.682053e+00;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -1.824116e+03 * (s - 3.200000e+00) - 1.440330e+01;
}

static double jarquebera_jbtbl8(double s) {
   if (s <= 1.3000) {
      double x = 2 * (s - 0.000000) / 1.300000 - 1, x2 = x * x;
      const double c0 = -7.199015e-01, c1 = -1.095921e+00, c2 = -4.736828e-01, c3 = -1.047438e-01;
      const double c4 = -2.484320e-03, c5 = +7.937923e-03, c6 = +4.810470e-03, c7 = +2.139780e-03;
      const double c8 = +6.708443e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 2.0000) {
      double x = 2 * (s - 1.300000) / 0.700000 - 1, x2 = x * x;
      const double c0 = -3.378966e+00, c1 = -7.802461e-01, c2 = +1.547593e-01, c3 = -6.241042e-02;
      const double c4 = +1.203274e-02, c5 = +5.201990e-03, c6 = -5.125597e-03, c7 = +1.584426e-03;
      const double c8 = +2.546069e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 5.0000) {
      double x = 2 * (s - 2.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -6.828366e+00, c1 = -3.137533e+00, c2 = -5.016671e-01, c3 = -1.745637e-01;
      const double c4 = -5.189801e-02, c5 = -1.621610e-02, c6 = -6.741122e-03, c7 = -4.516368e-03;
      const double c8 = +3.552085e-04, c9 = +2.787029e-03, ca = +5.359774e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -5.087028e+00 * (s - 5.000000e+00) - 1.071300e+01;
}

static double jarquebera_jbtbl9(double s) {
   if (s <= 1.3000) {
      double x = 2 * (s - 0.000000) / 1.300000 - 1, x2 = x * x;
      const double c0 = -6.279320e-01, c1 = -9.277151e-01, c2 = -3.669339e-01, c3 = -7.086149e-02;
      const double c4 = -1.333816e-03, c5 = +3.871249e-03, c6 = +2.007048e-03, c7 = +7.482245e-04;
      const double c8 = +2.355615e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 2.0000) {
      double x = 2 * (s - 1.300000) / 0.700000 - 1, x2 = x * x;
      const double c0 = -2.981430e+00, c1 = -7.972248e-01, c2 = +1.747737e-01, c3 = -3.808530e-02;
      const double c4 = -7.888305e-03, c5 = +9.001302e-03, c6 = -1.378767e-03, c7 = -1.108510e-03;
      const double c8 = +5.915372e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 7.0000) {
      double x = 2 * (s - 2.000000) / 5.000000 - 1, x2 = x * x;
      const double c0 = -6.387463e+00, c1 = -2.845231e+00, c2 = -1.809956e-01, c3 = -7.543461e-02;
      const double c4 = -4.880397e-03, c5 = -1.160074e-02, c6 = -7.356527e-03, c7 = -4.394428e-03;
      const double c8 = +9.619892e-04, c9 = -2.758763e-04, ca = +4.790977e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -2.020952e+00 * (s - 7.000000e+00) - 9.516623e+00;
}

static double jarquebera_jbtbl10(double s) {
   if (s <= 1.2000) {
      double x = 2 * (s - 0.000000) / 1.200000 - 1, x2 = x * x;
      const double c0 = -4.590993e-01, c1 = -6.562730e-01, c2 = -2.353934e-01, c3 = -4.069933e-02;
      const double c4 = -1.849151e-03, c5 = +8.931406e-04, c6 = +3.636295e-04, c7 = +1.178340e-05;
      const double c8 = -8.917749e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 2.0000) {
      double x = 2 * (s - 1.200000) / 0.800000 - 1, x2 = x * x;
      const double c0 = -2.537658e+00, c1 = -9.962401e-01, c2 = +1.838715e-01, c3 = +1.055792e-02;
      const double c4 = -2.580316e-02, c5 = +1.781701e-03, c6 = +3.770362e-03, c7 = -4.838983e-04;
      const double c8 = -6.999052e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 7.0000) {
      double x = 2 * (s - 2.000000) / 5.000000 - 1, x2 = x * x;
      const double c0 = -5.337524e+00, c1 = -1.877029e+00, c2 = +4.734650e-02, c3 = -4.249254e-02;
      const double c4 = +3.320250e-03, c5 = -6.432266e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5;
      return result > 0.0 ? 0.0 : result;
   } else return -8.711035e-01 * (s - 7.000000e+00) - 7.212811e+00;
}

static double jarquebera_jbtbl11(double s) {
   if (s <= 1.2000) {
      double x = 2 * (s - 0.000000) / 1.200000 - 1, x2 = x * x;
      const double c0 = -4.339517e-01, c1 = -6.051558e-01, c2 = -2.000992e-01, c3 = -3.022547e-02;
      const double c4 = -9.808401e-04, c5 = +5.592870e-04, c6 = +3.575081e-04, c7 = +2.086173e-04;
      const double c8 = +6.089011e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 2.2500) {
      double x = 2 * (s - 1.200000) / 1.050000 - 1, x2 = x * x;
      const double c0 = -2.523221e+00, c1 = -1.068388e+00, c2 = +2.179661e-01, c3 = -1.555524e-03;
      const double c4 = -3.238964e-02, c5 = +7.364320e-03, c6 = +4.895771e-03, c7 = -1.762774e-03;
      const double c8 = -8.201340e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 8.0000) {
      double x = 2 * (s - 2.250000) / 5.750000 - 1, x2 = x * x;
      const double c0 = -5.212179e+00, c1 = -1.684579e+00, c2 = +8.299519e-02, c3 = -3.606261e-02;
      const double c4 = +7.310869e-03, c5 = -3.320115e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5;
      return result > 0.0 ? 0.0 : result;
   } else return -5.715445e-01 * (s - 8.000000e+00) - 6.845834e+00;
}

static double jarquebera_jbtbl12(double s) {
   if (s <= 1.0000) {
      double x = 2 * (s - 0.000000) / 1.000000 - 1, x2 = x * x;
      const double c0 = -2.736742e-01, c1 = -3.657836e-01, c2 = -1.047209e-01, c3 = -1.319599e-02;
      const double c4 = -5.545631e-04, c5 = +9.280445e-05, c6 = +2.815679e-05, c7 = -2.213519e-05;
      const double c8 = +1.256838e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 3.0000) {
      double x = 2 * (s - 1.000000) / 2.000000 - 1, x2 = x * x;
      const double c0 = -2.573947e+00, c1 = -1.515287e+00, c2 = +3.611880e-01, c3 = -3.271311e-02;
      const double c4 = -6.495815e-02, c5 = +4.141186e-02, c6 = +7.180886e-04, c7 = -1.388211e-02;
      const double c8 = +4.890761e-03, c9 = +3.233175e-03, ca = -2.946156e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 12.0000) {
      double x = 2 * (s - 3.000000) / 9.000000 - 1, x2 = x * x;
      const double c0 = -5.947819e+00, c1 = -2.034157e+00, c2 = +6.878986e-02, c3 = -4.078603e-02;
      const double c4 = +6.990977e-03, c5 = -2.866215e-03, c6 = +3.897866e-03, c7 = +2.512252e-03;
      const double c8 = +2.073743e-03, c9 = +3.022621e-03, ca = +1.501343e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -2.877243e-01 * (s - 1.200000e+01) - 7.936839e+00;
}

static double jarquebera_jbtbl13(double s) {
   if (s <= 1.0000) {
      double x = 2 * (s - 0.000000) / 1.000000 - 1, x2 = x * x;
      const double c0 = -2.713276e-01, c1 = -3.557541e-01, c2 = -9.459092e-02, c3 = -1.044145e-02;
      const double c4 = -2.546132e-04, c5 = +1.002374e-04, c6 = +2.349456e-05, c7 = -7.025669e-05;
      const double c8 = -1.590242e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 3.0000) {
      double x = 2 * (s - 1.000000) / 2.000000 - 1, x2 = x * x;
      const double c0 = -2.454383e+00, c1 = -1.467539e+00, c2 = +3.270774e-01, c3 = -8.075763e-03;
      const double c4 = -6.611647e-02, c5 = +2.990785e-02, c6 = +8.109212e-03, c7 = -1.135031e-02;
      const double c8 = +5.915919e-04, c9 = +3.522390e-03, ca = -1.144701e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 13.0000) {
      double x = 2 * (s - 3.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -5.736127e+00, c1 = -1.920809e+00, c2 = +1.175858e-01, c3 = -4.002049e-02;
      const double c4 = +1.158966e-02, c5 = -3.157781e-03, c6 = +2.762172e-03, c7 = +5.780347e-04;
      const double c8 = -1.193310e-03, c9 = -2.442421e-05, ca = +2.547756e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -2.799944e-01 * (s - 1.300000e+01) - 7.566269e+00;
}

static double jarquebera_jbtbl14(double s) {
   if (s <= 1.0000) {
      double x = 2 * (s - 0.000000) / 1.000000 - 1, x2 = x * x;
      const double c0 = -2.698527e-01, c1 = -3.479081e-01, c2 = -8.640733e-02, c3 = -8.466899e-03;
      const double c4 = -1.469485e-04, c5 = +2.150009e-05, c6 = +1.965975e-05, c7 = -4.710210e-05;
      const double c8 = -1.327808e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 3.0000) {
      double x = 2 * (s - 1.000000) / 2.000000 - 1, x2 = x * x;
      const double c0 = -2.350359e+00, c1 = -1.421365e+00, c2 = +2.960468e-01, c3 = +1.149167e-02;
      const double c4 = -6.361109e-02, c5 = +1.976022e-02, c6 = +1.082700e-02, c7 = -8.563328e-03;
      const double c8 = -1.453123e-03, c9 = +2.917559e-03, ca = -1.151067e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 3.000000) / 12.000000 - 1, x2 = x * x;
      const double c0 = -5.746892e+00, c1 = -2.010441e+00, c2 = +1.566146e-01, c3 = -5.129690e-02;
      const double c4 = +1.929724e-02, c5 = -2.524227e-03, c6 = +3.192933e-03, c7 = -4.254730e-04;
      const double c8 = +1.620685e-03, c9 = +7.289618e-04, ca = -2.112350e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -2.590621e-01 * (s - 1.500000e+01) - 7.632238e+00;
}

static double jarquebera_jbtbl15(double s) {
   if (s <= 2.0000) {
      double x = 2 * (s - 0.000000) / 2.000000 - 1, x2 = x * x;
      const double c0 = -1.043660e+00, c1 = -1.361653e+00, c2 = -3.009497e-01, c3 = +4.951784e-02;
      const double c4 = +4.377903e-02, c5 = +1.003253e-02, c6 = -1.271309e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 5.0000) {
      double x = 2 * (s - 2.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -3.582778e+00, c1 = -8.349578e-01, c2 = +9.476514e-02, c3 = -2.717385e-02;
      const double c4 = +1.222591e-02, c5 = -6.635124e-03, c6 = +2.815993e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 17.0000) {
      double x = 2 * (s - 5.000000) / 12.000000 - 1, x2 = x * x;
      const double c0 = -6.115476e+00, c1 = -1.655936e+00, c2 = +8.404310e-02, c3 = -2.663794e-02;
      const double c4 = +8.868618e-03, c5 = +1.381447e-03, c6 = +9.444801e-04, c7 = -1.581503e-04;
      const double c8 = -9.468696e-04, c9 = +1.728509e-03, ca = +1.206470e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -1.927937e-01 * (s - 1.700000e+01) - 7.700983e+00;
}

static double jarquebera_jbtbl16(double s) {
   if (s <= 2.0000) {
      double x = 2 * (s - 0.000000) / 2.000000 - 1, x2 = x * x;
      const double c0 = -1.002570e+00, c1 = -1.298141e+00, c2 = -2.832803e-01, c3 = +3.877026e-02;
      const double c4 = +3.539436e-02, c5 = +8.439658e-03, c6 = -4.756911e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 5.0000) {
      double x = 2 * (s - 2.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -3.486198e+00, c1 = -8.242944e-01, c2 = +1.020002e-01, c3 = -3.130531e-02;
      const double c4 = +1.512373e-02, c5 = -8.054876e-03, c6 = +3.556839e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 20.0000) {
      double x = 2 * (s - 5.000000) / 15.000000 - 1, x2 = x * x;
      const double c0 = -6.241608e+00, c1 = -1.832655e+00, c2 = +1.340545e-01, c3 = -3.361143e-02;
      const double c4 = +1.283219e-02, c5 = +3.484549e-03, c6 = +1.805968e-03, c7 = -2.057243e-03;
      const double c8 = -1.454439e-03, c9 = -2.177513e-03, ca = -1.819209e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -2.391580e-01 * (s - 2.000000e+01) - 7.963205e+00;
}

static double jarquebera_jbtbl17(double s) {
   if (s <= 3.0000) {
      double x = 2 * (s - 0.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -1.566973e+00, c1 = -1.810330e+00, c2 = -4.840039e-02, c3 = +2.337294e-01;
      const double c4 = -5.383549e-04, c5 = -5.556515e-02, c6 = -8.656965e-03, c7 = +1.404569e-02;
      const double c8 = +6.447867e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 6.0000) {
      double x = 2 * (s - 3.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -3.905684e+00, c1 = -6.222920e-01, c2 = +4.146667e-02, c3 = -4.809176e-03;
      const double c4 = +1.057028e-03, c5 = -1.211838e-04, c6 = -4.099683e-04, c7 = +1.161105e-04;
      const double c8 = +2.225465e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 24.0000) {
      double x = 2 * (s - 6.000000) / 18.000000 - 1, x2 = x * x;
      const double c0 = -6.594282e+00, c1 = -1.917838e+00, c2 = +1.455980e-01, c3 = -2.999589e-02;
      const double c4 = +5.604263e-03, c5 = -3.484445e-03, c6 = -1.819937e-03, c7 = -2.930390e-03;
      const double c8 = +2.771761e-04, c9 = -6.232581e-04, ca = -7.029083e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -2.127771e-01 * (s - 2.400000e+01) - 8.400197e+00;
}

static double jarquebera_jbtbl18(double s) {
   if (s <= 3.0000) {
      double x = 2 * (s - 0.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -1.526802e+00, c1 = -1.762373e+00, c2 = -5.598890e-02, c3 = +2.189437e-01;
      const double c4 = +5.971721e-03, c5 = -4.823067e-02, c6 = -1.064501e-02, c7 = +1.014932e-02;
      const double c8 = +5.953513e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 6.0000) {
      double x = 2 * (s - 3.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -3.818669e+00, c1 = -6.070918e-01, c2 = +4.277196e-02, c3 = -4.879817e-03;
      const double c4 = +6.887357e-04, c5 = +1.638451e-05, c6 = +1.502800e-04, c7 = -3.165796e-05;
      const double c8 = +5.034960e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 20.0000) {
      double x = 2 * (s - 6.000000) / 14.000000 - 1, x2 = x * x;
      const double c0 = -6.010656e+00, c1 = -1.496296e+00, c2 = +1.002227e-01, c3 = -2.338250e-02;
      const double c4 = +4.137036e-03, c5 = -2.586202e-03, c6 = -9.736384e-04, c7 = +1.332251e-03;
      const double c8 = +1.877982e-03, c9 = -1.160963e-05, ca = -2.547247e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -1.684623e-01 * (s - 2.000000e+01) - 7.428883e+00;
}

static double jarquebera_jbtbl19(double s) {
   if (s <= 3.0000) {
      double x = 2 * (s - 0.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -1.490213e+00, c1 = -1.719633e+00, c2 = -6.459123e-02, c3 = +2.034878e-01;
      const double c4 = +1.113868e-02, c5 = -4.030922e-02, c6 = -1.054022e-02, c7 = +7.525623e-03;
      const double c8 = +5.277360e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 6.0000) {
      double x = 2 * (s - 3.000000) / 3.000000 - 1, x2 = x * x;
      const double c0 = -3.744750e+00, c1 = -5.977749e-01, c2 = +4.223716e-02, c3 = -5.363889e-03;
      const double c4 = +5.711774e-04, c5 = -5.557257e-04, c6 = +4.254794e-04, c7 = +9.034207e-05;
      const double c8 = +5.498107e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 20.0000) {
      double x = 2 * (s - 6.000000) / 14.000000 - 1, x2 = x * x;
      const double c0 = -5.872768e+00, c1 = -1.430689e+00, c2 = +1.136575e-01, c3 = -1.726627e-02;
      const double c4 = +3.421110e-03, c5 = -1.581510e-03, c6 = -5.559520e-04, c7 = -6.838208e-04;
      const double c8 = +8.428839e-04, c9 = -7.170682e-04, ca = -6.006647e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else return -1.539373e-01 * (s - 2.000000e+01) - 7.206941e+00;
}

static double jarquebera_jbtbl20(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.854794e+00, c1 = -1.948947e+00, c2 = +1.632184e-01, c3 = +2.139397e-01;
      const double c4 = -1.006237e-01, c5 = -3.810031e-02, c6 = +3.573620e-02, c7 = +9.951242e-03;
      const double c8 = -1.274092e-02, c9 = -3.464196e-03, ca = +4.882139e-03, cb = +1.575144e-03;
      const double cc = -1.822804e-03, cd = -7.061348e-04, ce = +5.908404e-04, cf = +1.978353e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8, tb = x2 * ta - t9;
      double tc = x2 * tb - ta, td = x2 * tc - tb, te = x2 * td - tc, tf = x2 * te - td;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta + cb * tb + cc * tc + cd * td + ce * te + cf * tf;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -5.030989e+00, c1 = -1.327151e+00, c2 = +1.346404e-01, c3 = -2.840051e-02;
      const double c4 = +7.578551e-03, c5 = -9.813886e-04, c6 = +5.905973e-05, c7 = -5.358489e-04;
      const double c8 = -3.450795e-04, c9 = -6.941157e-04, ca = -7.432418e-04, cb = -2.070537e-04;
      const double cc = +9.375654e-04, cd = +5.367378e-04, ce = +9.890859e-04, cf = +6.679782e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8, tb = x2 * ta - t9;
      double tc = x2 * tb - ta, td = x2 * tc - tb, te = x2 * td - tc, tf = x2 * te - td;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta + cb * tb + cc * tc + cd * td + ce * te + cf * tf;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -7.015854e+00, c1 = -7.487737e-01, c2 = +2.244254e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.318007e-01 * (s - 2.500000e+01) - 7.742185e+00;
}

static double jarquebera_jbtbl30(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.630822e+00, c1 = -1.724298e+00, c2 = +7.872756e-02, c3 = +1.658268e-01;
      const double c4 = -3.573597e-02, c5 = -2.994157e-02, c6 = +5.994825e-03, c7 = +7.394303e-03;
      const double c8 = -5.785029e-04, c9 = -1.990264e-03, ca = -1.037838e-04, cb = +6.755546e-04;
      const double cc = +1.774473e-04, cd = -2.821395e-04, ce = -1.392603e-04, cf = +1.353313e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8, tb = x2 * ta - t9;
      double tc = x2 * tb - ta, td = x2 * tc - tb, te = x2 * td - tc, tf = x2 * te - td;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta + cb * tb + cc * tc + cd * td + ce * te + cf * tf;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.539322e+00, c1 = -1.197018e+00, c2 = +1.396848e-01, c3 = -2.804293e-02;
      const double c4 = +6.867928e-03, c5 = -2.768758e-03, c6 = +5.211792e-04, c7 = +4.925799e-04;
      const double c8 = +5.046235e-04, c9 = -9.536469e-05, ca = -6.489642e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -6.263462e+00, c1 = -6.177316e-01, c2 = +2.590637e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.028212e-01 * (s - 2.500000e+01) - 6.855288e+00;
}

static double jarquebera_jbtbl50(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.436279e+00, c1 = -1.519711e+00, c2 = +1.148699e-02, c3 = +1.001204e-01;
      const double c4 = -3.207620e-03, c5 = -1.034778e-02, c6 = -1.220322e-03, c7 = +1.033260e-03;
      const double c8 = +2.588280e-04, c9 = -1.851653e-04, ca = -1.287733e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.234645e+00, c1 = -1.189127e+00, c2 = +1.429738e-01, c3 = -3.058822e-02;
      const double c4 = +9.086776e-03, c5 = -1.445783e-03, c6 = +1.311671e-03, c7 = -7.261298e-04;
      const double c8 = +6.496987e-04, c9 = +2.605249e-04, ca = +8.162282e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -5.921095e+00, c1 = -5.888603e-01, c2 = +3.080113e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -9.313116e-02 * (s - 2.500000e+01) - 6.479154e+00;
}

static double jarquebera_jbtbl65(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.360024e+00, c1 = -1.434631e+00, c2 = -6.514580e-03, c3 = +7.332038e-02;
      const double c4 = +1.158197e-03, c5 = -5.121233e-03, c6 = -1.051056e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.148601e+00, c1 = -1.214233e+00, c2 = +1.487977e-01, c3 = -3.424720e-02;
      const double c4 = +1.116715e-02, c5 = -4.043152e-03, c6 = +1.718149e-03, c7 = -1.313701e-03;
      const double c8 = +3.097305e-04, c9 = +2.181031e-04, ca = +1.256975e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4, t7 = x2 * t6 - t5;
      double t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -5.858951e+00, c1 = -5.895179e-01, c2 = +2.933237e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -9.443768e-02 * (s - 2.500000e+01) - 6.419137e+00;
}

static double jarquebera_jbtbl100(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.257021e+00, c1 = -1.313418e+00, c2 = -1.628931e-02, c3 = +4.264287e-02;
      const double c4 = +1.518487e-03, c5 = -1.499826e-03, c6 = -4.836044e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.056508e+00, c1 = -1.279690e+00, c2 = +1.665746e-01, c3 = -4.290012e-02;
      const double c4 = +1.487632e-02, c5 = -5.704465e-03, c6 = +2.211669e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -5.866099e+00, c1 = -6.399767e-01, c2 = +2.498208e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.080097e-01 * (s - 2.500000e+01) - 6.481094e+00;
}

static double jarquebera_jbtbl130(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.207999e+00, c1 = -1.253864e+00, c2 = -1.618032e-02, c3 = +3.112729e-02;
      const double c4 = +1.210546e-03, c5 = -4.732602e-04, c6 = -2.410527e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.026324e+00, c1 = -1.331990e+00, c2 = +1.779129e-01, c3 = -4.674749e-02;
      const double c4 = +1.669077e-02, c5 = -5.679136e-03, c6 = +8.833221e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -5.893951e+00, c1 = -6.475304e-01, c2 = +3.116734e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.045722e-01 * (s - 2.500000e+01) - 6.510314e+00;
}

static double jarquebera_jbtbl200(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.146155e+00, c1 = -1.177398e+00, c2 = -1.297970e-02, c3 = +1.869745e-02;
      const double c4 = +1.717288e-04, c5 = -1.982108e-04, c6 = +6.427636e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.034235e+00, c1 = -1.455006e+00, c2 = +1.942996e-01, c3 = -4.973795e-02;
      const double c4 = +1.418812e-02, c5 = -3.156778e-03, c6 = +4.896705e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -6.086071e+00, c1 = -7.152176e-01, c2 = +3.725393e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.132404e-01 * (s - 2.500000e+01) - 6.764034e+00;
}

static double jarquebera_jbtbl301(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.104290e+00, c1 = -1.125800e+00, c2 = -9.595847e-03, c3 = +1.219666e-02;
      const double c4 = +1.502210e-04, c5 = -6.414543e-05, c6 = +6.754115e-05;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.065955e+00, c1 = -1.582060e+00, c2 = +2.004472e-01, c3 = -4.709092e-02;
      const double c4 = +1.105779e-02, c5 = +1.197391e-03, c6 = -8.386780e-04;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3, t6 = x2 * t5 - t4;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -6.311384e+00, c1 = -7.918763e-01, c2 = +3.626584e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.293626e-01 * (s - 2.500000e+01) - 7.066995e+00;
}

static double jarquebera_jbtbl501(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.067426e+00, c1 = -1.079765e+00, c2 = -5.463005e-03, c3 = +6.875659e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.127574e+00, c1 = -1.740694e+00, c2 = +2.044502e-01, c3 = -3.746714e-02;
      const double c4 = +3.810594e-04, c5 = +1.197111e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -6.628194e+00, c1 = -8.846221e-01, c2 = +4.386405e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.418332e-01 * (s - 2.500000e+01) - 7.468952e+00;
}

static double jarquebera_jbtbl701(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.050999e+00, c1 = -1.059769e+00, c2 = -3.922680e-03, c3 = +4.847054e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.192182e+00, c1 = -1.860007e+00, c2 = +1.963942e-01, c3 = -2.838711e-02;
      const double c4 = -2.893112e-04, c5 = +2.159788e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -6.917851e+00, c1 = -9.817020e-01, c2 = +5.383727e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -1.532706e-01 * (s - 2.500000e+01) - 7.845715e+00;
}

static double jarquebera_jbtbl1401(double s) {
   if (s <= 4.0000) {
      double x = 2 * (s - 0.000000) / 4.000000 - 1, x2 = x * x;
      const double c0 = -1.026266e+00, c1 = -1.030061e+00, c2 = -1.259222e-03, c3 = +2.536254e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 15.0000) {
      double x = 2 * (s - 4.000000) / 11.000000 - 1, x2 = x * x;
      const double c0 = -4.329849e+00, c1 = -2.095443e+00, c2 = +1.759363e-01, c3 = -7.751359e-03;
      const double c4 = -6.124368e-03, c5 = -1.793114e-03;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1;
      double t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
      double result = c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5;
      return result > 0.0 ? 0.0 : result;
   } else if (s <= 25.0000) {
      double x = 2 * (s - 15.000000) / 10.000000 - 1, x2 = x * x;
      const double c0 = -7.544330e+00, c1 = -1.225382e+00, c2 = +5.392349e-02;
      double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0;
      double result = c0 * t0 + c1 * t1 + c2 * t2;
      return result > 0.0 ? 0.0 : result;
   } else return -2.019375e-01 * (s - 2.500000e+01) - 8.715788e+00;
}

static double jarquebera_jarqueberaapprox(ae_int_t n, double s) {
   ae_frame _frame_block;
   double x = s;
   ae_frame_make(&_frame_block);
   NewVector(vx, 0, DT_REAL);
   NewVector(vy, 0, DT_REAL);
   NewMatrix(ctbl, 0, 0, DT_REAL);
   double result;
// N = 5, 6, ..., 20, 30, 50, 65, 100, 130, 200, 301, 501, 701, 1401 are tabulated.
// In-between values up to 1400 are interpolated using an interpolating polynomial of the second degree.
// Anything under 5 or beyond 1400 is tabulated by an asymptotic expansion.
   if (n <= 20) switch (n) {
      default: ae_frame_leave(); return 1.0;
      case 5: result = jarquebera_jbtbl5(x); break;
      case 6: result = jarquebera_jbtbl6(x); break;
      case 7: result = jarquebera_jbtbl7(x); break;
      case 8: result = jarquebera_jbtbl8(x); break;
      case 9: result = jarquebera_jbtbl9(x); break;
      case 10: result = jarquebera_jbtbl10(x); break;
      case 11: result = jarquebera_jbtbl11(x); break;
      case 12: result = jarquebera_jbtbl12(x); break;
      case 13: result = jarquebera_jbtbl13(x); break;
      case 14: result = jarquebera_jbtbl14(x); break;
      case 15: result = jarquebera_jbtbl15(x); break;
      case 16: result = jarquebera_jbtbl16(x); break;
      case 17: result = jarquebera_jbtbl17(x); break;
      case 18: result = jarquebera_jbtbl18(x); break;
      case 19: result = jarquebera_jbtbl19(x); break;
      case 20: result = jarquebera_jbtbl20(x); break;
   } else if (n > 20 && n <= 50) {
      double t1 = -1.0 / 20.0, t2 = -1.0 / 30.0, t3 = -1.0 / 50.0, t = -1.0 / n;
      double f1 = jarquebera_jbtbl20(x), f2 = jarquebera_jbtbl30(x), f3 = jarquebera_jbtbl50(x);
      double f12 = ((t - t2) * f1 + (t1 - t) * f2) / (t1 - t2), f23 = ((t - t3) * f2 + (t2 - t) * f3) / (t2 - t3);
      result = ((t - t3) * f12 + (t1 - t) * f23) / (t1 - t3);
      if (result > 0.0) result = 0.0;
   } else if (n > 50 && n <= 100) {
      double t1 = -1.0 / 50.0, t2 = -1.0 / 65.0, t3 = -1.0 / 100.0, t = -1.0 / n;
      double f1 = jarquebera_jbtbl50(x), f2 = jarquebera_jbtbl65(x), f3 = jarquebera_jbtbl100(x);
      double f12 = ((t - t2) * f1 + (t1 - t) * f2) / (t1 - t2), f23 = ((t - t3) * f2 + (t2 - t) * f3) / (t2 - t3);
      result = ((t - t3) * f12 + (t1 - t) * f23) / (t1 - t3);
      if (result > 0.0) result = 0.0;
   } else if (n > 100 && n <= 200) {
      double t1 = -1.0 / 100.0, t2 = -1.0 / 130.0, t3 = -1.0 / 200.0, t = -1.0 / n;
      double f1 = jarquebera_jbtbl100(x), f2 = jarquebera_jbtbl130(x), f3 = jarquebera_jbtbl200(x);
      double f12 = ((t - t2) * f1 + (t1 - t) * f2) / (t1 - t2), f23 = ((t - t3) * f2 + (t2 - t) * f3) / (t2 - t3);
      result = ((t - t3) * f12 + (t1 - t) * f23) / (t1 - t3);
      if (result > 0.0) result = 0.0;
   } else if (n > 200 && n <= 501) {
      double t1 = -1.0 / 200.0, t2 = -1.0 / 301.0, t3 = -1.0 / 501.0, t = -1.0 / n;
      double f1 = jarquebera_jbtbl200(x), f2 = jarquebera_jbtbl301(x), f3 = jarquebera_jbtbl501(x);
      double f12 = ((t - t2) * f1 + (t1 - t) * f2) / (t1 - t2), f23 = ((t - t3) * f2 + (t2 - t) * f3) / (t2 - t3);
      result = ((t - t3) * f12 + (t1 - t) * f23) / (t1 - t3);
      if (result > 0.0) result = 0.0;
   } else if (n > 501 && n <= 1401) {
      double t1 = -1.0 / 501.0, t2 = -1.0 / 701.0, t3 = -1.0 / 1401.0, t = -1.0 / n;
      double f1 = jarquebera_jbtbl501(x), f2 = jarquebera_jbtbl701(x), f3 = jarquebera_jbtbl1401(x);
      double f12 = ((t - t2) * f1 + (t1 - t) * f2) / (t1 - t2), f23 = ((t - t3) * f2 + (t2 - t) * f3) / (t2 - t3);
      result = ((t - t3) * f12 + (t1 - t) * f23) / (t1 - t3);
      if (result > 0.0) result = 0.0;
   } else {
      result = -0.5 * x + (jarquebera_jbtbl1401(x) + 0.5 * x) * sqrt(1401.0 / n);
      if (result > 0.0) result = 0.0;
   }
   result = exp(result);
   ae_frame_leave();
   return result;
}

// Jarque-Bera test
//
// This test checks hypotheses about the fact that a  given  sample  X  is  a
// sample of normal random variable.
//
// Requirements:
//     * the number of elements in the sample is not less than 5.
//
// Inputs:
//     X   -   sample. Array whose index goes from 0 to N-1.
//     N   -   size of the sample. N >= 5
//
// Outputs:
//     P           -   p-value for the test
//
// Accuracy of the approximation used (5 <= N <= 1951):
//
// p-value          relative error (5 <= N <= 1951)
// [1, 0.1]            < 1%
// [0.1, 0.01]         < 2%
// [0.01, 0.001]       < 6%
// [0.001, 0]          wasn't measured
//
// For N > 1951 accuracy wasn't measured but it shouldn't be sharply  different
// from table values.
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p);
void jarqueberatest(RVector *x, ae_int_t n, double *p) {
   *p = 0;
// N is too small.
   if (n < 5) {
      *p = 1.0;
      return;
   }
// N is large enough.
   double mean = 0.0;
   double variance = 0.0;
   double skewness = 0.0;
   double kurtosis = 0.0;
   double stddev = 0.0;
   ae_assert(n > 1, "Assertion failed");
// Mean.
   for (ae_int_t i = 0; i < n; i++) mean += x->xR[i];
   mean /= n;
// Variance (using corrected two-pass algorithm).
   if (n != 1) {
      double v1 = 0.0;
      for (ae_int_t i = 0; i < n; i++) v1 += sqr(x->xR[i] - mean);
      double v2 = 0.0;
      for (ae_int_t i = 0; i < n; i++) v2 += x->xR[i] - mean;
      v2 = sqr(v2) / n;
      variance = (v1 - v2) / (n - 1);
      if (variance < 0.0) variance = 0.0;
      stddev = sqrt(variance);
   }
// Skewness and kurtosis.
   if (stddev != 0.0) {
      for (ae_int_t i = 0; i < n; i++) {
         double v = (x->xR[i] - mean) / stddev;
         double v2 = v * v;
         skewness += v2 * v;
         kurtosis += v2 * v2;
      }
      skewness /= n;
      kurtosis = kurtosis / n - 3;
   }
// Statistic.
   *p = jarquebera_jarqueberaapprox(n, (double)n / 6.0 *(sqr(skewness) + sqr(kurtosis) / 4.0));
}
} // end of namespace alglib_impl

namespace alglib {
void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::jarqueberatest(ConstT(ae_vector, x), n, &p);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === VARIANCETESTS Package ===
// Depends on: (SpecialFunctions) FDISTR, CHISQUAREDISTR
namespace alglib_impl {
// Two-sample F-test
//
// This test checks three hypotheses about dispersions of the given  samples.
// The following tests are performed:
//     * two-tailed test (null hypothesis - the dispersions are equal)
//     * left-tailed test (null hypothesis  -  the  dispersion  of  the first
//       sample is greater than or equal to  the  dispersion  of  the  second
//       sample).
//     * right-tailed test (null hypothesis - the  dispersion  of  the  first
//       sample is less than or equal to the dispersion of the second sample)
//
// The test is based on the following assumptions:
//     * the given samples have normal distributions
//     * the samples are independent.
//
// Inputs:
//     X   -   sample 1. Array whose index goes from 0 to N-1.
//     N   -   sample size.
//     Y   -   sample 2. Array whose index goes from 0 to M-1.
//     M   -   sample size.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
// ALGLIB: Copyright 19.09.2006 by Sergey Bochkanov
// API: void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
void ftest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail) {
   ae_int_t i;
   double xmean;
   double ymean;
   double xvar;
   double yvar;
   ae_int_t df1;
   ae_int_t df2;
   double stat;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   if (n <= 2 || m <= 2) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Mean
   xmean = 0.0;
   for (i = 0; i < n; i++) {
      xmean += x->xR[i];
   }
   xmean /= n;
   ymean = 0.0;
   for (i = 0; i < m; i++) {
      ymean += y->xR[i];
   }
   ymean /= m;
// Variance (using corrected two-pass algorithm)
   xvar = 0.0;
   for (i = 0; i < n; i++) {
      xvar += sqr(x->xR[i] - xmean);
   }
   xvar /= n - 1;
   yvar = 0.0;
   for (i = 0; i < m; i++) {
      yvar += sqr(y->xR[i] - ymean);
   }
   yvar /= m - 1;
   if (xvar == 0.0 || yvar == 0.0) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Statistic
   df1 = n - 1;
   df2 = m - 1;
   stat = rmin2(xvar / yvar, yvar / xvar);
   *bothtails = 1 - (fdistribution(df1, df2, 1 / stat) - fdistribution(df1, df2, stat));
   *lefttail = fdistribution(df1, df2, xvar / yvar);
   *righttail = 1 - (*lefttail);
}

// One-sample chi-square test
//
// This test checks three hypotheses about the dispersion of the given sample
// The following tests are performed:
//     * two-tailed test (null hypothesis - the dispersion equals  the  given
//       number)
//     * left-tailed test (null hypothesis - the dispersion is  greater  than
//       or equal to the given number)
//     * right-tailed test (null hypothesis  -  dispersion is  less  than  or
//       equal to the given number).
//
// Test is based on the following assumptions:
//     * the given sample has a normal distribution.
//
// Inputs:
//     X           -   sample 1. Array whose index goes from 0 to N-1.
//     N           -   size of the sample.
//     Variance    -   dispersion value to compare with.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
// ALGLIB: Copyright 19.09.2006 by Sergey Bochkanov
// API: void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail);
void onesamplevariancetest(RVector *x, ae_int_t n, double variance, double *bothtails, double *lefttail, double *righttail) {
   ae_int_t i;
   double xmean;
   double xvar;
   double s;
   double stat;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   if (n <= 1) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Mean
   xmean = 0.0;
   for (i = 0; i < n; i++) {
      xmean += x->xR[i];
   }
   xmean /= n;
// Variance
   xvar = 0.0;
   for (i = 0; i < n; i++) {
      xvar += sqr(x->xR[i] - xmean);
   }
   xvar /= n - 1;
   if (xvar == 0.0) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Statistic
   stat = (n - 1) * xvar / variance;
   s = chisquaredistribution((double)(n - 1), stat);
   *bothtails = 2 * rmin2(s, 1 - s);
   *lefttail = s;
   *righttail = 1 - (*lefttail);
}
} // end of namespace alglib_impl

namespace alglib {
void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ftest(ConstT(ae_vector, x), n, ConstT(ae_vector, y), m, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}

void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::onesamplevariancetest(ConstT(ae_vector, x), n, variance, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === WSR Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
// Tail(S, 5)
static double wsr_w5(double s) {
   static double tab[] = {
      -3.466e+00, -2.773e+00, -2.367e+00, -1.856e+00, -1.520e+00, -1.163e+00, -9.008e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-3.708099e+00 * s + 7.500000e+00);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 6)
static double wsr_w6(double s) {
   static double tab[] = {
      -4.159e+00, -3.466e+00, -3.060e+00, -2.549e+00, -2.213e+00, -1.856e+00, -1.520e+00, -1.269e+00,
      -1.068e+00, -8.630e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-4.769696e+00 * s + 1.050000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 7)
static double wsr_w7(double s) {
   static double tab[] = {
      -4.852e+00, -4.159e+00, -3.753e+00, -3.243e+00, -2.906e+00, -2.549e+00, -2.213e+00, -1.908e+00,
      -1.674e+00, -1.451e+00, -1.241e+00, -1.068e+00, -9.008e-01, -7.577e-01, -6.325e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-5.916080e+00 * s + 1.400000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 8)
static double wsr_w8(double s) {
   static double tab[] = {
      -5.545e+00, -4.852e+00, -4.447e+00, -3.936e+00, -3.599e+00, -3.243e+00, -2.906e+00, -2.601e+00,
      -2.326e+00, -2.079e+00, -1.856e+00, -1.653e+00, -1.468e+00, -1.297e+00, -1.138e+00, -9.913e-01,
      -8.630e-01, -7.494e-01, -6.399e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-7.141428e+00 * s + 1.800000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 9)
static double wsr_w9(double s) {
   static double tab[] = {
      -6.238e+00, -5.545e+00, -5.140e+00, -4.629e+00, -4.292e+00, -3.936e+00, -3.599e+00, -3.294e+00,
      -3.019e+00, -2.742e+00, -2.501e+00, -2.287e+00, -2.079e+00, -1.895e+00, -1.717e+00, -1.547e+00,
      -1.394e+00, -1.255e+00, -1.120e+00, -1.002e+00, -8.912e-01, -7.873e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-8.440972e+00 * s + 2.250000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 10)
static double wsr_w10(double s) {
   static double tab[] = {
      -6.931e+00, -6.238e+00, -5.833e+00, -5.322e+00, -4.986e+00, -4.629e+00, -4.292e+00, -3.987e+00,
      -3.713e+00, -3.435e+00, -3.170e+00, -2.942e+00, -2.727e+00, -2.525e+00, -2.336e+00, -2.152e+00,
      -1.983e+00, -1.826e+00, -1.674e+00, -1.533e+00, -1.402e+00, -1.279e+00, -1.163e+00, -1.057e+00,
      -9.551e-01, -8.607e-01, -7.745e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-9.810708e+00 * s + 2.750000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 11)
static double wsr_w11(double s) {
   static double tab[] = {
      -7.625e+00, -6.931e+00, -6.526e+00, -6.015e+00, -5.679e+00, -5.322e+00, -4.986e+00, -4.680e+00,
      -4.406e+00, -4.128e+00, -3.863e+00, -3.617e+00, -3.391e+00, -3.182e+00, -2.980e+00, -2.788e+00,
      -2.607e+00, -2.437e+00, -2.273e+00, -2.119e+00, -1.972e+00, -1.832e+00, -1.700e+00, -1.577e+00,
      -1.459e+00, -1.348e+00, -1.243e+00, -1.143e+00, -1.050e+00, -9.615e-01, -8.782e-01, -8.002e-01,
      -7.279e-01, -6.595e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-1.124722e+01 * s + 3.300000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 12)
static double wsr_w12(double s) {
   static double tab[] = {
      -8.318e+00, -7.625e+00, -7.219e+00, -6.708e+00, -6.372e+00, -6.015e+00, -5.679e+00, -5.373e+00,
      -5.099e+00, -4.821e+00, -4.557e+00, -4.310e+00, -4.069e+00, -3.852e+00, -3.645e+00, -3.443e+00,
      -3.255e+00, -3.076e+00, -2.902e+00, -2.738e+00, -2.581e+00, -2.429e+00, -2.285e+00, -2.148e+00,
      -2.017e+00, -1.893e+00, -1.774e+00, -1.660e+00, -1.552e+00, -1.449e+00, -1.350e+00, -1.256e+00,
      -1.168e+00, -1.083e+00, -1.003e+00, -9.276e-01, -8.556e-01, -7.878e-01, -7.239e-01, -6.633e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-1.274755e+01 * s + 3.900000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 13)
static double wsr_w13(double s) {
   static double tab[] = {
      -9.011e+00, -8.318e+00, -7.912e+00, -7.401e+00, -7.065e+00, -6.708e+00, -6.372e+00, -6.066e+00,
      -5.792e+00, -5.514e+00, -5.250e+00, -5.004e+00, -4.762e+00, -4.534e+00, -4.320e+00, -4.113e+00,
      -3.917e+00, -3.733e+00, -3.551e+00, -3.380e+00, -3.215e+00, -3.055e+00, -2.902e+00, -2.755e+00,
      -2.614e+00, -2.479e+00, -2.349e+00, -2.224e+00, -2.104e+00, -1.990e+00, -1.879e+00, -1.773e+00,
      -1.672e+00, -1.574e+00, -1.481e+00, -1.392e+00, -1.306e+00, -1.224e+00, -1.146e+00, -1.072e+00,
      -1.001e+00, -9.328e-01, -8.683e-01, -8.068e-01, -7.486e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-1.430909e+01 * s + 4.550000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 14)
static double wsr_w14(double s) {
   static double tab[] = {
      -9.704e+00, -9.011e+00, -8.605e+00, -8.095e+00, -7.758e+00, -7.401e+00, -7.065e+00, -6.760e+00,
      -6.485e+00, -6.208e+00, -5.943e+00, -5.697e+00, -5.456e+00, -5.227e+00, -5.004e+00, -4.791e+00,
      -4.592e+00, -4.401e+00, -4.215e+00, -4.038e+00, -3.866e+00, -3.700e+00, -3.541e+00, -3.387e+00,
      -3.238e+00, -3.095e+00, -2.956e+00, -2.823e+00, -2.694e+00, -2.570e+00, -2.450e+00, -2.334e+00,
      -2.223e+00, -2.115e+00, -2.011e+00, -1.911e+00, -1.815e+00, -1.722e+00, -1.632e+00, -1.546e+00,
      -1.463e+00, -1.383e+00, -1.306e+00, -1.233e+00, -1.162e+00, -1.094e+00, -1.029e+00, -9.664e-01,
      -9.067e-01, -8.495e-01, -7.950e-01, -7.428e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-1.592953e+01 * s + 5.250000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 15)
static double wsr_w15(double s) {
   static double tab[] = {
      -1.040e+01, -9.704e+00, -9.299e+00, -8.788e+00, -8.451e+00, -8.095e+00, -7.758e+00, -7.453e+00,
      -7.178e+00, -6.901e+00, -6.636e+00, -6.390e+00, -6.149e+00, -5.920e+00, -5.697e+00, -5.477e+00,
      -5.273e+00, -5.079e+00, -4.888e+00, -4.707e+00, -4.531e+00, -4.359e+00, -4.195e+00, -4.036e+00,
      -3.881e+00, -3.732e+00, -3.587e+00, -3.446e+00, -3.310e+00, -3.179e+00, -3.051e+00, -2.928e+00,
      -2.809e+00, -2.693e+00, -2.581e+00, -2.472e+00, -2.366e+00, -2.265e+00, -2.166e+00, -2.070e+00,
      -1.977e+00, -1.888e+00, -1.801e+00, -1.717e+00, -1.636e+00, -1.558e+00, -1.482e+00, -1.409e+00,
      -1.339e+00, -1.270e+00, -1.205e+00, -1.142e+00, -1.081e+00, -1.022e+00, -9.656e-01, -9.114e-01,
      -8.593e-01, -8.093e-01, -7.613e-01, -7.154e-01, -6.714e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-1.760682e+01 * s + 6.000000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 16)
static double wsr_w16(double s) {
   static double tab[] = {
      -1.109e+01, -1.040e+01, -9.992e+00, -9.481e+00, -9.144e+00, -8.788e+00, -8.451e+00, -8.146e+00,
      -7.871e+00, -7.594e+00, -7.329e+00, -7.083e+00, -6.842e+00, -6.613e+00, -6.390e+00, -6.170e+00,
      -5.960e+00, -5.762e+00, -5.569e+00, -5.383e+00, -5.204e+00, -5.029e+00, -4.860e+00, -4.697e+00,
      -4.538e+00, -4.383e+00, -4.234e+00, -4.088e+00, -3.947e+00, -3.810e+00, -3.676e+00, -3.546e+00,
      -3.420e+00, -3.298e+00, -3.179e+00, -3.064e+00, -2.951e+00, -2.842e+00, -2.735e+00, -2.632e+00,
      -2.532e+00, -2.434e+00, -2.339e+00, -2.247e+00, -2.158e+00, -2.071e+00, -1.986e+00, -1.904e+00,
      -1.825e+00, -1.748e+00, -1.673e+00, -1.600e+00, -1.530e+00, -1.462e+00, -1.395e+00, -1.331e+00,
      -1.270e+00, -1.210e+00, -1.152e+00, -1.096e+00, -1.042e+00, -9.895e-01, -9.391e-01, -8.905e-01,
      -8.437e-01, -7.986e-01, -7.551e-01, -7.134e-01, -6.733e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-1.933908e+01 * s + 6.800000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 17)
static double wsr_w17(double s) {
   static double tab[] = {
      -1.178e+01, -1.109e+01, -1.068e+01, -1.017e+01, -9.838e+00, -9.481e+00, -9.144e+00, -8.839e+00,
      -8.565e+00, -8.287e+00, -8.022e+00, -7.776e+00, -7.535e+00, -7.306e+00, -7.083e+00, -6.864e+00,
      -6.654e+00, -6.451e+00, -6.254e+00, -6.066e+00, -5.884e+00, -5.706e+00, -5.534e+00, -5.367e+00,
      -5.204e+00, -5.047e+00, -4.893e+00, -4.743e+00, -4.597e+00, -4.456e+00, -4.317e+00, -4.183e+00,
      -4.052e+00, -3.924e+00, -3.799e+00, -3.678e+00, -3.560e+00, -3.445e+00, -3.332e+00, -3.223e+00,
      -3.116e+00, -3.012e+00, -2.911e+00, -2.812e+00, -2.715e+00, -2.621e+00, -2.530e+00, -2.440e+00,
      -2.353e+00, -2.269e+00, -2.186e+00, -2.106e+00, -2.028e+00, -1.951e+00, -1.877e+00, -1.805e+00,
      -1.735e+00, -1.666e+00, -1.600e+00, -1.536e+00, -1.473e+00, -1.412e+00, -1.353e+00, -1.295e+00,
      -1.240e+00, -1.185e+00, -1.133e+00, -1.082e+00, -1.033e+00, -9.853e-01, -9.392e-01, -8.946e-01,
      -8.514e-01, -8.097e-01, -7.695e-01, -7.306e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-2.112463e+01 * s + 7.650000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 18)
static double wsr_w18(double s) {
   static double tab[] = {
      -1.248e+01, -1.178e+01, -1.138e+01, -1.087e+01, -1.053e+01, -1.017e+01, -9.838e+00, -9.532e+00,
      -9.258e+00, -8.980e+00, -8.715e+00, -8.469e+00, -8.228e+00, -7.999e+00, -7.776e+00, -7.557e+00,
      -7.347e+00, -7.144e+00, -6.943e+00, -6.753e+00, -6.569e+00, -6.388e+00, -6.213e+00, -6.044e+00,
      -5.878e+00, -5.717e+00, -5.561e+00, -5.408e+00, -5.258e+00, -5.113e+00, -4.971e+00, -4.833e+00,
      -4.698e+00, -4.565e+00, -4.437e+00, -4.311e+00, -4.188e+00, -4.068e+00, -3.950e+00, -3.836e+00,
      -3.724e+00, -3.615e+00, -3.508e+00, -3.403e+00, -3.301e+00, -3.201e+00, -3.104e+00, -3.008e+00,
      -2.915e+00, -2.824e+00, -2.735e+00, -2.648e+00, -2.564e+00, -2.481e+00, -2.400e+00, -2.321e+00,
      -2.244e+00, -2.168e+00, -2.095e+00, -2.023e+00, -1.953e+00, -1.885e+00, -1.818e+00, -1.753e+00,
      -1.690e+00, -1.628e+00, -1.568e+00, -1.509e+00, -1.452e+00, -1.396e+00, -1.342e+00, -1.289e+00,
      -1.238e+00, -1.188e+00, -1.140e+00, -1.093e+00, -1.047e+00, -1.003e+00, -9.597e-01, -9.179e-01,
      -8.774e-01, -8.381e-01, -8.001e-01, -7.633e-01, -7.276e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-2.296193e+01 * s + 8.550000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 19)
static double wsr_w19(double s) {
   static double tab[] = {
      -1.317e+01, -1.248e+01, -1.207e+01, -1.156e+01, -1.122e+01, -1.087e+01, -1.053e+01, -1.023e+01,
      -9.951e+00, -9.673e+00, -9.409e+00, -9.162e+00, -8.921e+00, -8.692e+00, -8.469e+00, -8.250e+00,
      -8.040e+00, -7.837e+00, -7.636e+00, -7.443e+00, -7.256e+00, -7.074e+00, -6.897e+00, -6.726e+00,
      -6.558e+00, -6.394e+00, -6.235e+00, -6.080e+00, -5.928e+00, -5.780e+00, -5.634e+00, -5.493e+00,
      -5.355e+00, -5.219e+00, -5.086e+00, -4.957e+00, -4.830e+00, -4.706e+00, -4.585e+00, -4.466e+00,
      -4.350e+00, -4.236e+00, -4.125e+00, -4.015e+00, -3.909e+00, -3.804e+00, -3.702e+00, -3.601e+00,
      -3.503e+00, -3.407e+00, -3.313e+00, -3.220e+00, -3.130e+00, -3.042e+00, -2.955e+00, -2.870e+00,
      -2.787e+00, -2.706e+00, -2.626e+00, -2.548e+00, -2.472e+00, -2.398e+00, -2.325e+00, -2.253e+00,
      -2.183e+00, -2.115e+00, -2.048e+00, -1.983e+00, -1.919e+00, -1.857e+00, -1.796e+00, -1.736e+00,
      -1.678e+00, -1.621e+00, -1.565e+00, -1.511e+00, -1.458e+00, -1.407e+00, -1.356e+00, -1.307e+00,
      -1.259e+00, -1.213e+00, -1.167e+00, -1.123e+00, -1.080e+00, -1.038e+00, -9.968e-01, -9.571e-01,
      -9.185e-01, -8.809e-01, -8.445e-01, -8.090e-01, -7.747e-01, -7.413e-01, -7.089e-01, -6.776e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-2.484955e+01 * s + 9.500000e+01);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 20)
static double wsr_w20(double s) {
   static double tab[] = {
      -1.386e+01, -1.317e+01, -1.276e+01, -1.225e+01, -1.192e+01, -1.156e+01, -1.122e+01, -1.092e+01,
      -1.064e+01, -1.037e+01, -1.010e+01, -9.856e+00, -9.614e+00, -9.386e+00, -9.162e+00, -8.943e+00,
      -8.733e+00, -8.530e+00, -8.330e+00, -8.136e+00, -7.947e+00, -7.763e+00, -7.584e+00, -7.411e+00,
      -7.242e+00, -7.076e+00, -6.915e+00, -6.757e+00, -6.603e+00, -6.453e+00, -6.305e+00, -6.161e+00,
      -6.020e+00, -5.882e+00, -5.746e+00, -5.614e+00, -5.484e+00, -5.356e+00, -5.232e+00, -5.109e+00,
      -4.990e+00, -4.872e+00, -4.757e+00, -4.644e+00, -4.534e+00, -4.425e+00, -4.318e+00, -4.214e+00,
      -4.111e+00, -4.011e+00, -3.912e+00, -3.815e+00, -3.721e+00, -3.627e+00, -3.536e+00, -3.446e+00,
      -3.358e+00, -3.272e+00, -3.187e+00, -3.104e+00, -3.023e+00, -2.943e+00, -2.865e+00, -2.788e+00,
      -2.713e+00, -2.639e+00, -2.566e+00, -2.495e+00, -2.426e+00, -2.357e+00, -2.290e+00, -2.225e+00,
      -2.161e+00, -2.098e+00, -2.036e+00, -1.976e+00, -1.916e+00, -1.859e+00, -1.802e+00, -1.746e+00,
      -1.692e+00, -1.639e+00, -1.587e+00, -1.536e+00, -1.486e+00, -1.438e+00, -1.390e+00, -1.344e+00,
      -1.299e+00, -1.254e+00, -1.211e+00, -1.169e+00, -1.128e+00, -1.087e+00, -1.048e+00, -1.010e+00,
      -9.726e-01, -9.363e-01, -9.010e-01, -8.665e-01, -8.330e-01, -8.004e-01, -7.686e-01, -7.378e-01,
      -7.078e-01, -6.787e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-2.678619e+01 * s + 1.050000e+02);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 21)
static double wsr_w21(double s) {
   static double tab[] = {
      -1.456e+01, -1.386e+01, -1.346e+01, -1.295e+01, -1.261e+01, -1.225e+01, -1.192e+01, -1.161e+01,
      -1.134e+01, -1.106e+01, -1.079e+01, -1.055e+01, -1.031e+01, -1.008e+01, -9.856e+00, -9.636e+00,
      -9.426e+00, -9.223e+00, -9.023e+00, -8.829e+00, -8.640e+00, -8.454e+00, -8.274e+00, -8.099e+00,
      -7.928e+00, -7.762e+00, -7.599e+00, -7.439e+00, -7.283e+00, -7.131e+00, -6.981e+00, -6.835e+00,
      -6.692e+00, -6.551e+00, -6.413e+00, -6.278e+00, -6.146e+00, -6.016e+00, -5.889e+00, -5.763e+00,
      -5.641e+00, -5.520e+00, -5.402e+00, -5.286e+00, -5.172e+00, -5.060e+00, -4.950e+00, -4.842e+00,
      -4.736e+00, -4.632e+00, -4.530e+00, -4.429e+00, -4.330e+00, -4.233e+00, -4.138e+00, -4.044e+00,
      -3.952e+00, -3.861e+00, -3.772e+00, -3.685e+00, -3.599e+00, -3.515e+00, -3.432e+00, -3.350e+00,
      -3.270e+00, -3.192e+00, -3.115e+00, -3.039e+00, -2.964e+00, -2.891e+00, -2.819e+00, -2.748e+00,
      -2.679e+00, -2.611e+00, -2.544e+00, -2.478e+00, -2.414e+00, -2.350e+00, -2.288e+00, -2.227e+00,
      -2.167e+00, -2.108e+00, -2.051e+00, -1.994e+00, -1.939e+00, -1.884e+00, -1.831e+00, -1.779e+00,
      -1.728e+00, -1.677e+00, -1.628e+00, -1.580e+00, -1.533e+00, -1.486e+00, -1.441e+00, -1.397e+00,
      -1.353e+00, -1.311e+00, -1.269e+00, -1.229e+00, -1.189e+00, -1.150e+00, -1.112e+00, -1.075e+00,
      -1.039e+00, -1.003e+00, -9.685e-01, -9.348e-01, -9.018e-01, -8.697e-01, -8.383e-01, -8.077e-01,
      -7.779e-01, -7.489e-01, -7.207e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-2.877064e+01 * s + 1.155000e+02);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 22)
static double wsr_w22(double s) {
   static double tab[] = {
      -1.525e+01, -1.456e+01, -1.415e+01, -1.364e+01, -1.330e+01, -1.295e+01, -1.261e+01, -1.230e+01,
      -1.203e+01, -1.175e+01, -1.149e+01, -1.124e+01, -1.100e+01, -1.077e+01, -1.055e+01, -1.033e+01,
      -1.012e+01, -9.917e+00, -9.716e+00, -9.522e+00, -9.333e+00, -9.147e+00, -8.965e+00, -8.789e+00,
      -8.617e+00, -8.449e+00, -8.285e+00, -8.124e+00, -7.966e+00, -7.813e+00, -7.661e+00, -7.513e+00,
      -7.368e+00, -7.226e+00, -7.086e+00, -6.949e+00, -6.815e+00, -6.683e+00, -6.553e+00, -6.426e+00,
      -6.300e+00, -6.177e+00, -6.057e+00, -5.938e+00, -5.821e+00, -5.706e+00, -5.594e+00, -5.483e+00,
      -5.374e+00, -5.266e+00, -5.161e+00, -5.057e+00, -4.955e+00, -4.855e+00, -4.756e+00, -4.659e+00,
      -4.563e+00, -4.469e+00, -4.376e+00, -4.285e+00, -4.196e+00, -4.108e+00, -4.021e+00, -3.935e+00,
      -3.851e+00, -3.769e+00, -3.687e+00, -3.607e+00, -3.529e+00, -3.451e+00, -3.375e+00, -3.300e+00,
      -3.226e+00, -3.153e+00, -3.082e+00, -3.012e+00, -2.943e+00, -2.875e+00, -2.808e+00, -2.742e+00,
      -2.677e+00, -2.614e+00, -2.551e+00, -2.490e+00, -2.429e+00, -2.370e+00, -2.312e+00, -2.254e+00,
      -2.198e+00, -2.142e+00, -2.088e+00, -2.034e+00, -1.982e+00, -1.930e+00, -1.880e+00, -1.830e+00,
      -1.781e+00, -1.733e+00, -1.686e+00, -1.640e+00, -1.594e+00, -1.550e+00, -1.506e+00, -1.464e+00,
      -1.422e+00, -1.381e+00, -1.340e+00, -1.301e+00, -1.262e+00, -1.224e+00, -1.187e+00, -1.151e+00,
      -1.115e+00, -1.080e+00, -1.046e+00, -1.013e+00, -9.805e-01, -9.486e-01, -9.175e-01, -8.871e-01,
      -8.573e-01, -8.283e-01, -7.999e-01, -7.722e-01, -7.452e-01, -7.189e-01, -6.931e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-3.080179e+01 * s + 1.265000e+02);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 23)
static double wsr_w23(double s) {
   static double tab[] = {
      -1.594e+01, -1.525e+01, -1.484e+01, -1.433e+01, -1.400e+01, -1.364e+01, -1.330e+01, -1.300e+01,
      -1.272e+01, -1.245e+01, -1.218e+01, -1.194e+01, -1.169e+01, -1.147e+01, -1.124e+01, -1.102e+01,
      -1.081e+01, -1.061e+01, -1.041e+01, -1.022e+01, -1.003e+01, -9.840e+00, -9.658e+00, -9.481e+00,
      -9.308e+00, -9.139e+00, -8.974e+00, -8.811e+00, -8.653e+00, -8.498e+00, -8.345e+00, -8.196e+00,
      -8.049e+00, -7.905e+00, -7.764e+00, -7.625e+00, -7.489e+00, -7.355e+00, -7.224e+00, -7.094e+00,
      -6.967e+00, -6.842e+00, -6.719e+00, -6.598e+00, -6.479e+00, -6.362e+00, -6.247e+00, -6.133e+00,
      -6.022e+00, -5.912e+00, -5.804e+00, -5.697e+00, -5.592e+00, -5.489e+00, -5.388e+00, -5.287e+00,
      -5.189e+00, -5.092e+00, -4.996e+00, -4.902e+00, -4.809e+00, -4.718e+00, -4.628e+00, -4.539e+00,
      -4.451e+00, -4.365e+00, -4.280e+00, -4.197e+00, -4.114e+00, -4.033e+00, -3.953e+00, -3.875e+00,
      -3.797e+00, -3.721e+00, -3.645e+00, -3.571e+00, -3.498e+00, -3.426e+00, -3.355e+00, -3.285e+00,
      -3.216e+00, -3.149e+00, -3.082e+00, -3.016e+00, -2.951e+00, -2.888e+00, -2.825e+00, -2.763e+00,
      -2.702e+00, -2.642e+00, -2.583e+00, -2.525e+00, -2.468e+00, -2.412e+00, -2.356e+00, -2.302e+00,
      -2.248e+00, -2.195e+00, -2.144e+00, -2.093e+00, -2.042e+00, -1.993e+00, -1.944e+00, -1.897e+00,
      -1.850e+00, -1.804e+00, -1.758e+00, -1.714e+00, -1.670e+00, -1.627e+00, -1.585e+00, -1.543e+00,
      -1.502e+00, -1.462e+00, -1.423e+00, -1.384e+00, -1.347e+00, -1.309e+00, -1.273e+00, -1.237e+00,
      -1.202e+00, -1.168e+00, -1.134e+00, -1.101e+00, -1.069e+00, -1.037e+00, -1.006e+00, -9.755e-01,
      -9.457e-01, -9.166e-01, -8.880e-01, -8.601e-01, -8.328e-01, -8.061e-01, -7.800e-01, -7.544e-01,
      -7.295e-01, -7.051e-01, -6.813e-01
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-3.287856e+01 * s + 1.380000e+02);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 24)
static double wsr_w24(double s) {
   static double tab[] = {
      -1.664e+01, -1.594e+01, -1.554e+01, -1.503e+01, -1.469e+01, -1.433e+01, -1.400e+01, -1.369e+01,
      -1.342e+01, -1.314e+01, -1.287e+01, -1.263e+01, -1.239e+01, -1.216e+01, -1.194e+01, -1.172e+01,
      -1.151e+01, -1.130e+01, -1.110e+01, -1.091e+01, -1.072e+01, -1.053e+01, -1.035e+01, -1.017e+01,
      -1.000e+01, -9.830e+00, -9.664e+00, -9.501e+00, -9.341e+00, -9.185e+00, -9.031e+00, -8.881e+00,
      -8.733e+00, -8.588e+00, -8.445e+00, -8.305e+00, -8.167e+00, -8.032e+00, -7.899e+00, -7.768e+00,
      -7.639e+00, -7.512e+00, -7.387e+00, -7.265e+00, -7.144e+00, -7.025e+00, -6.907e+00, -6.792e+00,
      -6.678e+00, -6.566e+00, -6.456e+00, -6.347e+00, -6.240e+00, -6.134e+00, -6.030e+00, -5.927e+00,
      -5.826e+00, -5.727e+00, -5.628e+00, -5.531e+00, -5.436e+00, -5.342e+00, -5.249e+00, -5.157e+00,
      -5.067e+00, -4.978e+00, -4.890e+00, -4.803e+00, -4.718e+00, -4.633e+00, -4.550e+00, -4.468e+00,
      -4.387e+00, -4.308e+00, -4.229e+00, -4.151e+00, -4.075e+00, -4.000e+00, -3.925e+00, -3.852e+00,
      -3.779e+00, -3.708e+00, -3.637e+00, -3.568e+00, -3.499e+00, -3.432e+00, -3.365e+00, -3.300e+00,
      -3.235e+00, -3.171e+00, -3.108e+00, -3.046e+00, -2.985e+00, -2.925e+00, -2.865e+00, -2.806e+00,
      -2.749e+00, -2.692e+00, -2.636e+00, -2.580e+00, -2.526e+00, -2.472e+00, -2.419e+00, -2.367e+00,
      -2.316e+00, -2.265e+00, -2.215e+00, -2.166e+00, -2.118e+00, -2.070e+00, -2.024e+00, -1.978e+00,
      -1.932e+00, -1.888e+00, -1.844e+00, -1.800e+00, -1.758e+00, -1.716e+00, -1.675e+00, -1.634e+00,
      -1.594e+00, -1.555e+00, -1.517e+00, -1.479e+00, -1.441e+00, -1.405e+00, -1.369e+00, -1.334e+00,
      -1.299e+00, -1.265e+00, -1.231e+00, -1.198e+00, -1.166e+00, -1.135e+00, -1.103e+00, -1.073e+00,
      -1.043e+00, -1.014e+00, -9.849e-01, -9.567e-01, -9.291e-01, -9.020e-01, -8.754e-01, -8.494e-01,
      -8.239e-01, -7.990e-01, -7.746e-01, -7.507e-01, -7.273e-01, -7.044e-01, -6.820e-01,
   };
   const size_t n = sizeof tab/sizeof tab[0] - 1;
   ae_int_t w = RoundZ(-3.500000e+01 * s + 1.500000e+02);
   return tab[w < 0 ? 0 : w > n ? n : w];
}

// Tail(S, 25).
static double wsr_w25(double s) {
   const double c0 = -5.150509e+00, c1 = -5.695528e+00, c2 = -1.437637e+00, c3 = -2.611906e-01;
   const double c4 = -7.625722e-02, c5 = -2.579892e-02, c6 = -1.086876e-02, c7 = -2.906543e-03;
   const double c8 = -2.354881e-03, c9 = +1.007195e-04, ca = -8.437327e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 26).
static double wsr_w26(double s) {
   const double c0 = -5.117622e+00, c1 = -5.635159e+00, c2 = -1.395167e+00, c3 = -2.382823e-01;
   const double c4 = -6.531987e-02, c5 = -2.060112e-02, c6 = -8.203697e-03, c7 = -1.516523e-03;
   const double c8 = -1.431364e-03, c9 = +6.384553e-04, ca = -3.238369e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 27).
static double wsr_w27(double s) {
   const double c0 = -5.089731e+00, c1 = -5.584248e+00, c2 = -1.359966e+00, c3 = -2.203696e-01;
   const double c4 = -5.753344e-02, c5 = -1.761891e-02, c6 = -7.096897e-03, c7 = -1.419108e-03;
   const double c8 = -1.581214e-03, c9 = +3.033766e-04, ca = -5.901441e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 28).
static double wsr_w28(double s) {
   const double c0 = -5.065046e+00, c1 = -5.539163e+00, c2 = -1.328939e+00, c3 = -2.046376e-01;
   const double c4 = -5.061515e-02, c5 = -1.469271e-02, c6 = -5.711578e-03, c7 = -8.389153e-04;
   const double c8 = -1.250575e-03, c9 = +4.047245e-04, ca = -5.128555e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 29).
static double wsr_w29(double s) {
   const double c0 = -5.043413e+00, c1 = -5.499756e+00, c2 = -1.302137e+00, c3 = -1.915129e-01;
   const double c4 = -4.516329e-02, c5 = -1.260064e-02, c6 = -4.817269e-03, c7 = -5.478130e-04;
   const double c8 = -1.111668e-03, c9 = +4.093451e-04, ca = -5.135860e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 30).
static double wsr_w30(double s) {
   const double c0 = -5.024071e+00, c1 = -5.464515e+00, c2 = -1.278342e+00, c3 = -1.800030e-01;
   const double c4 = -4.046294e-02, c5 = -1.076162e-02, c6 = -3.968677e-03, c7 = -1.911679e-04;
   const double c8 = -8.619185e-04, c9 = +5.125362e-04, ca = -3.984370e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 40).
static double wsr_w40(double s) {
   const double c0 = -4.904809e+00, c1 = -5.248327e+00, c2 = -1.136698e+00, c3 = -1.170982e-01;
   const double c4 = -1.824427e-02, c5 = -3.888648e-03, c6 = -1.344929e-03, c7 = +2.790407e-04;
   const double c8 = -4.619858e-04, c9 = +3.359121e-04, ca = -2.883026e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 60).
static double wsr_w60(double s) {
   const double c0 = -4.809656e+00, c1 = -5.077191e+00, c2 = -1.029402e+00, c3 = -7.507931e-02;
   const double c4 = -6.506226e-03, c5 = -1.391278e-03, c6 = -4.263635e-04, c7 = +2.302271e-04;
   const double c8 = -2.384348e-04, c9 = +1.865587e-04, ca = -1.622355e-04;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 120).
static double wsr_w120(double s) {
   const double c0 = -4.729426e+00, c1 = -4.934426e+00, c2 = -9.433231e-01, c3 = -4.492504e-02;
   const double c4 = +1.673948e-05, c5 = -6.077014e-04, c6 = -7.215768e-05, c7 = +9.086734e-05;
   const double c8 = -8.447980e-05, c9 = +6.705028e-05, ca = -5.828507e-05;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S, 200).
static double wsr_w200(double s) {
   const double c0 = -4.700240e+00, c1 = -4.883080e+00, c2 = -9.132168e-01, c3 = -3.512684e-02;
   const double c4 = +1.726342e-03, c5 = -5.189796e-04, c6 = -1.628659e-06, c7 = +4.261786e-05;
   const double c8 = -4.002498e-05, c9 = +3.146287e-05, ca = -2.727576e-05;
   double x = rmin2(s / 2.0 - 1.0, 1.0), x2 = 2.0 * x;
   double t0 = 1.0, t1 = x, t2 = x2 * t1 - t0, t3 = x2 * t2 - t1, t4 = x2 * t3 - t2, t5 = x2 * t4 - t3;
   double t6 = x2 * t5 - t4, t7 = x2 * t6 - t5, t8 = x2 * t7 - t6, t9 = x2 * t8 - t7, ta = x2 * t9 - t8;
   return c0 * t0 + c1 * t1 + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5 + c6 * t6 + c7 * t7 + c8 * t8 + c9 * t9 + ca * ta;
}

// Tail(S,N), S >= 0
static double wsr_wsigma(double s, ae_int_t n) {
   switch (n) {
      case 5: return wsr_w5(s);
      case 6: return wsr_w6(s);
      case 7: return wsr_w7(s);
      case 8: return wsr_w8(s);
      case 9: return wsr_w9(s);
      case 10: return wsr_w10(s);
      case 11: return wsr_w11(s);
      case 12: return wsr_w12(s);
      case 13: return wsr_w13(s);
      case 14: return wsr_w14(s);
      case 15: return wsr_w15(s);
      case 16: return wsr_w16(s);
      case 17: return wsr_w17(s);
      case 18: return wsr_w18(s);
      case 19: return wsr_w19(s);
      case 20: return wsr_w20(s);
      case 21: return wsr_w21(s);
      case 22: return wsr_w22(s);
      case 23: return wsr_w23(s);
      case 24: return wsr_w24(s);
      case 25: return wsr_w25(s);
      case 26: return wsr_w26(s);
      case 27: return wsr_w27(s);
      case 28: return wsr_w28(s);
      case 29: return wsr_w29(s);
      case 30: return wsr_w30(s);
      default:
         if (n > 30) {
            double x = 1.0 / n;
            const double x0 = 1.0 / 30, x1 = 1.0 / 40, x2 = 1.0 / 60, x3 = 1.0 / 120, x4 = 1.0 / 200;
            double f0 = wsr_w30(s), f1 = wsr_w40(s), f2 = wsr_w60(s), f3 = wsr_w120(s), f4 = wsr_w200(s);
            f1 = ((x - x0) * f1 - (x - x1) * f0) / (x1 - x0);
            f2 = ((x - x0) * f2 - (x - x2) * f0) / (x2 - x0);
            f3 = ((x - x0) * f3 - (x - x3) * f0) / (x3 - x0);
            f4 = ((x - x0) * f4 - (x - x4) * f0) / (x4 - x0);
            f2 = ((x - x1) * f2 - (x - x2) * f1) / (x2 - x1);
            f3 = ((x - x1) * f3 - (x - x3) * f1) / (x3 - x1);
            f4 = ((x - x1) * f4 - (x - x4) * f1) / (x4 - x1);
            f3 = ((x - x2) * f3 - (x - x3) * f2) / (x3 - x2);
            f4 = ((x - x2) * f4 - (x - x4) * f2) / (x4 - x2);
            f4 = ((x - x3) * f4 - (x - x4) * f3) / (x4 - x3);
            return f4;
         }
      return 0.0;
   }
}

// Wilcoxon signed-rank test
//
// This test checks three hypotheses about the median  of  the  given sample.
// The following tests are performed:
//     * two-tailed test (null hypothesis - the median is equal to the  given
//       value)
//     * left-tailed test (null hypothesis - the median is  greater  than  or
//       equal to the given value)
//     * right-tailed test (null hypothesis  -  the  median  is  less than or
//       equal to the given value)
//
// Requirements:
//     * the scale of measurement should be ordinal, interval or  ratio (i.e.
//       the test could not be applied to nominal variables).
//     * the distribution should be continuous and symmetric relative to  its
//       median.
//     * number of distinct values in the X array should be greater than 4
//
// The test is non-parametric and doesn't require distribution X to be normal
//
// Inputs:
//     X       -   sample. Array whose index goes from 0 to N-1.
//     N       -   size of the sample.
//     Median  -   assumed median value.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
//
// To calculate p-values, special approximation is used. This method lets  us
// calculate p-values with two decimal places in interval [0.0001, 1].
//
// "Two decimal places" does not sound very impressive, but in  practice  the
// relative error of less than 1% is enough to make a decision.
//
// There is no approximation outside the [0.0001, 1] interval. Therefore,  if
// the significance level outlies this interval, the test returns 0.0001.
// ALGLIB: Copyright 08.09.2006 by Sergey Bochkanov
// API: void wilcoxonsignedranktest(const real_1d_array &x, const ae_int_t n, const double e, double &bothtails, double &lefttail, double &righttail);
void wilcoxonsignedranktest(RVector *x, ae_int_t n, double e, double *bothtails, double *lefttail, double *righttail) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t t;
   ae_int_t ns;
   double w;
   double p;
   double mp;
   double s;
   double sigma;
   double mu;
   ae_frame_make(&_frame_block);
   DupVector(x);
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   NewVector(r, 0, DT_REAL);
   NewVector(c, 0, DT_INT);
// Prepare
   if (n < 5) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      ae_frame_leave();
      return;
   }
   ns = 0;
   for (i = 0; i < n; i++) {
      if (x->xR[i] == e) {
         continue;
      }
      x->xR[ns] = x->xR[i];
      ns++;
   }
   if (ns < 5) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      ae_frame_leave();
      return;
   }
   ae_vector_set_length(&r, ns);
   ae_vector_set_length(&c, ns);
   for (i = 0; i < ns; i++) {
      r.xR[i] = fabs(x->xR[i] - e);
      c.xZ[i] = i;
   }
// sort {R, C}
   if (ns != 1) {
      i = 2;
      do {
         t = i;
         while (t != 1) {
            k = t / 2;
            if (r.xR[k - 1] >= r.xR[t - 1]) {
               t = 1;
            } else {
               swapr(&r.xR[k - 1], &r.xR[t - 1]);
               swapi(&c.xZ[k - 1], &c.xZ[t - 1]);
               t = k;
            }
         }
         i++;
      } while (i <= ns);
      i = ns - 1;
      do {
         swapr(&r.xR[i], &r.xR[0]);
         swapi(&c.xZ[i], &c.xZ[0]);
         t = 1;
         while (t != 0) {
            k = 2 * t;
            if (k > i) {
               t = 0;
            } else {
               if (k < i) {
                  if (r.xR[k] > r.xR[k - 1]) {
                     k++;
                  }
               }
               if (r.xR[t - 1] >= r.xR[k - 1]) {
                  t = 0;
               } else {
                  swapr(&r.xR[k - 1], &r.xR[t - 1]);
                  swapi(&c.xZ[k - 1], &c.xZ[t - 1]);
                  t = k;
               }
            }
         }
         i--;
      } while (i >= 1);
   }
// compute tied ranks
   i = 0;
   while (i < ns) {
      j = i + 1;
      while (j < ns) {
         if (r.xR[j] != r.xR[i]) {
            break;
         }
         j++;
      }
      for (k = i; k < j; k++) {
         r.xR[k] = 1 + (double)(i + j - 1) / 2.0;
      }
      i = j;
   }
// Compute W+
   w = 0.0;
   for (i = 0; i < ns; i++) {
      if (x->xR[c.xZ[i]] > e) {
         w += r.xR[i];
      }
   }
// Result
   mu = (double)ns * (ns + 1) / 4;
   sigma = sqrt(mu * (2 * ns + 1) / 6);
   s = (w - mu) / sigma;
   if (s <= 0.0) {
      p = exp(wsr_wsigma(-(w - mu) / sigma, ns));
      mp = 1 - exp(wsr_wsigma(-(w - 1 - mu) / sigma, ns));
   } else {
      mp = exp(wsr_wsigma((w - mu) / sigma, ns));
      p = 1 - exp(wsr_wsigma((w + 1 - mu) / sigma, ns));
   }
   *lefttail = rmax2(p, 1.0E-4);
   *righttail = rmax2(mp, 1.0E-4);
   *bothtails = 2 * rmin2(*lefttail, *righttail);
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void wilcoxonsignedranktest(const real_1d_array &x, const ae_int_t n, const double e, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::wilcoxonsignedranktest(ConstT(ae_vector, x), n, e, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === MANNWHITNEYU Package ===
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
// Sequential Chebyshev interpolation.
static void mannwhitneyu_ucheb(double x, double c, double *tj, double *tj1, double *r) {
   double t;
   *r += c * (*tj);
   t = 2 * x * (*tj1) - (*tj);
   *tj = *tj1;
   *tj1 = t;
}

// Tail(S, 5, 5)
static double mannwhitneyu_utbln5n5(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 2.611165e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -2.596264e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.412086e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.858542e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.614282e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.372686e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.524731e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.435331e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.284665e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.184141e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.298360e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 7.447272e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.938769e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.276205e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.138481e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.684625e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.558104e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 6)
static double mannwhitneyu_utbln5n6(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 2.738613e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -2.810459e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.684429e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.712858e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.009324e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.644391e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.034173e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.953498e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.279293e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.563485e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.971952e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.506309e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.541406e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.283205e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.016347e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.221626e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.286752e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 7)
static double mannwhitneyu_utbln5n7(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 2.841993e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -2.994677e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.923264e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.506190e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.054280e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.794587e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.726290e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.534180e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.517845e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.904428e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.882443e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.482988e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.114875e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.515082e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.996056e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.293581e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.349444e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 8)
static double mannwhitneyu_utbln5n8(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 2.927700e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.155727e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.135078e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.247203e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.309697e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.993725e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.567219e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.383704e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.002188e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.487322e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.443899e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.688270e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.600339e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.874948e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.811593e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.072353e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.659457e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 9)
static double mannwhitneyu_utbln5n9(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.000000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.298162e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.325016e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.939852e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.563029e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.222652e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.195200e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.445665e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.204792e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.775217e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.527781e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.221948e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.242968e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.607959e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.771285e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.694026e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.481190e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 10)
static double mannwhitneyu_utbln5n10(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.061862e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.425360e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.496710e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.587658e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.812005e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.427637e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.515702e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.406867e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.796295e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.237591e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.654249e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.181165e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.011665e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.417927e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.534880e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.791255e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.871512e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 11)
static double mannwhitneyu_utbln5n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.115427e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.539959e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.652998e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.196503e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.054363e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.618848e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.109411e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.786668e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.215648e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.484220e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.935991e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.396191e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.894177e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.206979e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.519055e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.210326e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.189679e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 12)
static double mannwhitneyu_utbln5n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.162278e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.644007e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.796173e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.771177e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.290043e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.794686e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.702110e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.185959e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.416259e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.592056e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.201530e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.754365e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.978945e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.012032e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.304579e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.100378e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.728269e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 13)
static double mannwhitneyu_utbln5n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.203616e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.739120e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.928117e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.031605e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.519403e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.962648e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.292183e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.809293e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.465156e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.456278e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.446055e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.109490e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.218256e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.941479e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.058603e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.824402e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.830947e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 14)
static double mannwhitneyu_utbln5n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.240370e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.826559e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.050370e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.083408e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.743164e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.012030e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.884686e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.059656e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.327521e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.134026e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.584201e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.440618e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.524133e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.990007e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.887334e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.534977e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.705395e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 15)
static double mannwhitneyu_utbln5n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.851572e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.082033e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.095983e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.814595e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.073148e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.420213e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.517175e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.344180e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.371393e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.711443e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.228569e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.683483e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.267112e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.156044e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.131316e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.301023e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 16)
static double mannwhitneyu_utbln5n16(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.852210e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.077482e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.091186e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.797282e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.084994e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.667054e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.843909e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.456732e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.039830e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.723508e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.940608e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.478285e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.649144e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.237703e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.707410e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.874293e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 17)
static double mannwhitneyu_utbln5n17(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.851752e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.071259e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.084700e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.758898e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.073846e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.684838e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.964936e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.782442e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.956362e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.984727e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.196936e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.558262e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.690746e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.364855e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.401006e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.546748e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 18)
static double mannwhitneyu_utbln5n18(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.850840e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.064799e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.077651e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.712659e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.049217e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.571333e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.929809e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.752044e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.949464e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.896101e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.614460e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.384357e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.489113e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.445725e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.945636e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.424653e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 19)
static double mannwhitneyu_utbln5n19(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.850027e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.059159e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.071106e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.669960e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.022780e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.442555e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.851335e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.433865e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.514465e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.332989e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.606099e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.341945e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.402164e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.039761e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.512831e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.284427e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 20)
static double mannwhitneyu_utbln5n20(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.849651e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.054729e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.065747e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.636243e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.003234e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.372789e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.831551e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.763090e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.830626e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.122384e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.108328e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.557983e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.945666e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.965696e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.493236e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.162591e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 21)
static double mannwhitneyu_utbln5n21(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.849649e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.051155e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.061430e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.608869e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.902788e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.346562e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.874709e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.682887e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.026206e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.534551e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.990575e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.713334e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.737011e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.304571e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.133110e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.123457e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 22)
static double mannwhitneyu_utbln5n22(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.849598e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.047605e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.057264e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.579513e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.749602e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.275137e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.881768e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.177374e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.981056e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.696290e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.886803e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.085378e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.675242e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.426367e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.039613e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.662378e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 23)
static double mannwhitneyu_utbln5n23(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.849269e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.043761e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.052735e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.544683e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.517503e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.112082e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.782070e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.549483e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.747329e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.694263e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.147141e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.526209e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.039173e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.235615e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.656546e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.014423e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 24)
static double mannwhitneyu_utbln5n24(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.848925e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.040178e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.048355e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.510198e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.261134e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.915864e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.627423e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.307345e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.732992e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.869652e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.494176e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.047533e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.178439e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.424171e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.829195e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.840810e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 25)
static double mannwhitneyu_utbln5n25(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.848937e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.037512e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.044866e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.483269e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.063682e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.767778e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.508540e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.332756e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.881511e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.124041e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.368456e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.930499e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.779630e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.029528e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.658678e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.289695e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 26)
static double mannwhitneyu_utbln5n26(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.849416e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.035915e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.042493e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.466021e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.956432e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.698914e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.465689e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.035254e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.674614e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.492734e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.014021e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.944953e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.255750e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.075841e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.989330e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.134862e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 27)
static double mannwhitneyu_utbln5n27(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.850070e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.034815e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.040650e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.453117e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.886426e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.661702e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.452346e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.002476e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.720126e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.001400e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.729826e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.740640e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.206333e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.366093e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.193471e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.804091e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 28)
static double mannwhitneyu_utbln5n28(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.850668e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.033786e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.038853e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.440281e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.806020e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.612883e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.420436e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.787982e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.535230e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.263121e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.849609e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.863967e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.391610e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.720294e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.952273e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.901413e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 29)
static double mannwhitneyu_utbln5n29(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.851217e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.032834e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.037113e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.427762e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.719146e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.557172e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.375498e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.452033e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.187516e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.916936e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.065533e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.067301e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.615824e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.432244e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.417795e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.710038e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 30)
static double mannwhitneyu_utbln5n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.851845e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.032148e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.035679e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.417758e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.655330e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.522132e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.352106e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.326911e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.064969e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.813321e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.683881e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.813346e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.627085e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.832107e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.519336e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.888530e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 5, 100)
static double mannwhitneyu_utbln5n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.250000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.877940e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.039324e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.022243e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.305825e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.960119e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.112000e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.138868e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.418164e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.174520e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.489617e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.878301e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.302233e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.054113e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.458862e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.186591e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.623412e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 6)
static double mannwhitneyu_utbln6n6(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 2.882307e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.054075e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.998804e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.681518e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.067578e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.709435e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.952661e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.641700e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.304572e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.336275e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.770385e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.401891e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.246148e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.442663e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.502866e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.105855e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.739371e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 7)
static double mannwhitneyu_utbln6n7(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.000000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.265287e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.274613e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.582352e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.334293e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.915502e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.108091e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.546701e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.298827e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.891501e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.313717e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.989501e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.914594e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.062372e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.158841e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.596443e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.185662e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 8)
static double mannwhitneyu_utbln6n8(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.098387e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.450954e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.520462e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.420299e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.604853e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.165840e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.008756e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.723402e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.843521e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.883405e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.720980e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.301709e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.948034e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.776243e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.623736e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.742068e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.796927e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 9)
static double mannwhitneyu_utbln6n9(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.181981e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.616113e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.741650e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.204487e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.873068e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.446794e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.632286e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.266481e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.280067e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.780687e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.480242e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.592200e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.581019e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.264231e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.347174e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.167535e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.092185e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 10)
static double mannwhitneyu_utbln6n10(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.253957e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.764382e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.942366e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.939896e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.137812e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.720270e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.281070e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.901060e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.824937e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.802812e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.258132e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.233536e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.085530e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.212151e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.001329e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.226048e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.035298e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 11)
static double mannwhitneyu_utbln6n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.316625e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.898597e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.125710e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.063297e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.396852e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.990126e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.927977e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.726500e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.858745e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.654590e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.217736e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.989770e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.768493e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.924364e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.140215e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.647914e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.924802e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 12)
static double mannwhitneyu_utbln6n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.371709e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.020941e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.294250e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.128842e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.650389e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.248611e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.578510e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.162852e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.746982e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.454209e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.128042e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.936650e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.530794e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.665192e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.994144e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.662249e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.368541e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 13)
static double mannwhitneyu_utbln6n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.420526e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.133167e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.450016e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.191088e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.898220e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.050249e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.226901e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.471113e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.007470e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.049420e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.059074e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.881249e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.452780e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.441805e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.787493e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.483957e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.481590e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 14)
static double mannwhitneyu_utbln6n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.450000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.201268e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.542568e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.226965e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.046029e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.136657e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.786757e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.843748e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.588022e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.253029e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.667188e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.788330e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.474545e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.540494e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.951188e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.863323e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.220904e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 15)
static double mannwhitneyu_utbln6n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.450000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.195689e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.526567e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.213617e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.975035e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.118480e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.859142e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.083312e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.298720e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.766708e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.026356e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.093113e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.135168e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.136376e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.190870e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.435972e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.413129e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 30)
static double mannwhitneyu_utbln6n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.450000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.166269e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.427399e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.118239e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.360847e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.745885e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.025041e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.187179e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.432089e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.408451e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.388774e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.795560e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.304136e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.258516e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.180236e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.388679e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.836027e-06, &tj, &tj1, &result);
   return result;
}

// Tail(S, 6, 100)
static double mannwhitneyu_utbln6n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.450000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.181350e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.417919e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.094201e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.195883e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.818937e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.514202e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.125047e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.022148e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.284181e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.157766e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.023752e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.127985e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.221690e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.516179e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.501398e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.380220e-06, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 7)
static double mannwhitneyu_utbln7n7(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.130495e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.501264e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.584790e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.577311e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.617002e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.145186e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.023462e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.408251e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.626515e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.072492e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.722926e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.095445e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.842602e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.751427e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.008927e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.892431e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.772386e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 8)
static double mannwhitneyu_utbln7n8(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.240370e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.709965e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.862154e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.504541e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.900195e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.439995e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.678028e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.485540e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.437047e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.440092e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.114227e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.516569e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.829457e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.787550e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.761866e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.991911e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.533481e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 9)
static double mannwhitneyu_utbln7n9(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.334314e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.896550e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.112671e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.037277e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.181695e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.765190e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.360116e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.695960e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.780578e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.963843e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.616148e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.852104e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.390744e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.014041e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.888101e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.467474e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.004611e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 10)
static double mannwhitneyu_utbln7n10(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.415650e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.064844e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.340749e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.118888e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.459730e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.097781e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.057688e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.097406e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.209262e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.065641e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.196677e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.313994e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.827157e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.822284e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.389090e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.340850e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.395172e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 11)
static double mannwhitneyu_utbln7n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.486817e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.217795e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.549783e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.195905e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.733093e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.428447e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.760093e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.431676e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.717152e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.032199e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.832423e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.905979e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.302799e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.464371e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.456211e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.736244e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.140712e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 12)
static double mannwhitneyu_utbln7n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.500000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.235822e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.564100e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.190813e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.686546e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.395083e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.967359e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.747096e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.304144e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.903198e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.134906e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.175035e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.266224e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.892931e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.604706e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.070459e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.427010e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 13)
static double mannwhitneyu_utbln7n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.500000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.222204e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.532300e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.164642e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.523768e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.531984e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.467857e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.483804e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.524136e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.077740e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.745218e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.602085e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.828831e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.994070e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.873879e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.341937e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.706444e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 14)
static double mannwhitneyu_utbln7n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.500000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.211763e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.507542e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.143640e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.395755e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.808020e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.044259e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.182308e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.057325e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.724255e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.303900e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.113148e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.102514e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.559442e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.634986e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.776476e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.054489e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 15)
static double mannwhitneyu_utbln7n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.500000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.204898e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.489960e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.129172e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.316741e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.506107e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.983676e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.258013e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.262515e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.984156e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.912108e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.974023e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.056195e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.090842e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.232620e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.816339e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.020421e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 30)
static double mannwhitneyu_utbln7n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.500000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.176536e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.398705e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.045481e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.821982e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.962304e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.698132e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.062667e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.282353e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.014836e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.035683e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.004137e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.801453e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.920705e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.518735e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.821501e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.801008e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 7, 100)
static double mannwhitneyu_utbln7n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.500000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.188337e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.386949e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.022834e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.686517e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.323516e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.399392e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.644333e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.617044e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.031396e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.792066e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.675457e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.673416e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.258552e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.174214e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.073644e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.349958e-06, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 8)
static double mannwhitneyu_utbln8n8(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.360672e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -3.940217e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.168913e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.051485e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.195325e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.775196e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.385506e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.244902e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.525632e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.771275e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.332874e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.079599e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.882551e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.407944e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.769844e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.062433e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.872535e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 9)
static double mannwhitneyu_utbln8n9(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.464102e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.147004e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.446939e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.146155e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.488561e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.144561e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.116917e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.205667e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.515661e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.618616e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.599011e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.457324e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.482917e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.488267e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.469823e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.957591e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.058326e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 10)
static double mannwhitneyu_utbln8n10(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.554093e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.334282e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.700860e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.235253e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.778489e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.527324e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.862885e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.589781e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.507355e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.717526e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.215726e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.848696e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.918854e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.219614e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.753761e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.573688e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.602177e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 11)
static double mannwhitneyu_utbln8n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.421882e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.812457e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.266153e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.849344e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.971527e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.258944e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.944820e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.894685e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.031836e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.514330e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.351660e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.206748e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.492600e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.005338e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.780099e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.673599e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 12)
static double mannwhitneyu_utbln8n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.398211e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.762214e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.226296e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.603837e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.643223e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.502438e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.544574e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.647734e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.442259e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.011484e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.384758e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.998259e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.659985e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.331046e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.638478e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.056785e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 13)
static double mannwhitneyu_utbln8n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.380670e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.724511e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.195851e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.420511e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.609928e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.893999e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.115919e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.291410e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.339664e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.801548e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.534710e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.793250e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.806718e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.384624e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.120582e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.936453e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 14)
static double mannwhitneyu_utbln8n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.368494e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.697171e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.174440e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.300621e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.087393e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.685826e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.085254e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.525658e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.966647e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.453388e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.826066e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.501958e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.336297e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.251972e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.118456e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.415959e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 15)
static double mannwhitneyu_utbln8n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.358397e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.674485e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.155941e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.195780e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.544830e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.426183e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.309902e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.650956e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.068874e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.538544e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.192525e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.073905e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.079673e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.423572e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.579647e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.765904e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 30)
static double mannwhitneyu_utbln8n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.318823e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.567159e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.064864e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.688413e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.153712e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.309389e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.226861e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.523815e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.780987e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.166866e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.922431e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.466397e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.690036e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.008185e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.271903e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.534751e-06, &tj, &tj1, &result);
   return result;
}

// Tail(S, 8, 100)
static double mannwhitneyu_utbln8n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.600000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.324531e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.547071e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.038129e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.541549e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.525605e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.044992e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.085713e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.017871e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.459226e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.092064e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.024349e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 7.366347e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.385637e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.321722e-08, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.439286e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.058079e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 9)
static double mannwhitneyu_utbln9n9(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.576237e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.372857e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.750859e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.248233e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.792868e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.559372e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.894941e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.643256e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.091370e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.285034e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.112997e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.806229e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.150741e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.509825e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.891051e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.485013e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.343653e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 10)
static double mannwhitneyu_utbln9n10(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.516726e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.939333e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.305046e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.935326e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.029141e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.420592e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.053140e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.065930e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.523581e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.544888e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.813741e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.510631e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.536057e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.833815e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.189692e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.615050e-03, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 11)
static double mannwhitneyu_utbln9n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.481308e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.867483e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.249072e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.591790e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.400128e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.341992e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.463680e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.487211e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.671196e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.343472e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.544146e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.802335e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.117084e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.217443e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.858766e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.193687e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 12)
static double mannwhitneyu_utbln9n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.456776e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.817037e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.209788e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.362108e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.171356e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.661557e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.026141e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.361908e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.093885e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.298389e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.663603e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.768522e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.579015e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.868677e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.440652e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.523037e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 13)
static double mannwhitneyu_utbln9n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.438840e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.779308e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.180614e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.196489e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.346621e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.234857e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.796211e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.575715e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.525647e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.964651e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.275235e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.299124e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.397416e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.295781e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.237619e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 7.269692e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 14)
static double mannwhitneyu_utbln9n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.425981e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.751545e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.159543e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.086570e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.917446e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.120112e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.175519e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.515473e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.727772e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.070629e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.677569e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.876953e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.233502e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.508182e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.120389e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.847212e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 15)
static double mannwhitneyu_utbln9n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.414952e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.727612e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.140634e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.981231e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.382635e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.853575e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.571051e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.567625e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.214197e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.448700e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.712669e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.015050e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.438610e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.301363e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.309386e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.164772e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 30)
static double mannwhitneyu_utbln9n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.370720e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.615712e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.050023e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.504775e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.318265e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.646826e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.741492e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.735360e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.966911e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.100738e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.348991e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.527687e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.917286e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.397466e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.360175e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.892252e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 9, 100)
static double mannwhitneyu_utbln9n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.372506e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.590966e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.021758e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.359849e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.755519e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.533166e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.936659e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.634913e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.730053e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.791845e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.030682e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.228663e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.631175e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.636749e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.404599e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.789872e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 10)
static double mannwhitneyu_utbln10n10(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.468831e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.844398e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.231728e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.486073e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.781321e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.971425e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.215371e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.828451e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.419872e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.430165e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.740363e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.049211e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.269371e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.211393e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.232314e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.016081e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 11)
static double mannwhitneyu_utbln10n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.437998e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.782296e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.184732e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.219585e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.457012e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.296008e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.481501e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.527940e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.953426e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.563840e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.574403e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.535775e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.338037e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.002654e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.852676e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.318132e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 12)
static double mannwhitneyu_utbln10n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.416082e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.737458e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.150952e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.036884e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.609030e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.908684e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.439666e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.162647e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.451601e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.148757e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.803981e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.731621e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.346903e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.013151e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.956148e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.438381e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 13)
static double mannwhitneyu_utbln10n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.399480e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.702863e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.124829e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.897428e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.979802e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.634368e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.180461e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.484926e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.864376e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.186576e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.886925e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.836828e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.074756e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.209547e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.883266e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.380143e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 14)
static double mannwhitneyu_utbln10n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.386924e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.676124e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.104740e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.793826e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.558886e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.492462e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.052903e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.917782e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.878696e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.576046e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.764551e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.288778e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.757658e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.299101e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.265197e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.384503e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 15)
static double mannwhitneyu_utbln10n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.376846e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.654247e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.088083e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.705945e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.169677e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.317213e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.264836e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.548024e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.633910e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.505621e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.658588e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.320254e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.175277e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.122317e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.675688e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.661363e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 30)
static double mannwhitneyu_utbln10n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.333977e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.548099e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.004444e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.291014e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.523674e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.828211e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.716917e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.894256e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.433371e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.522675e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.764192e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.140235e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.629230e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.541895e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.944946e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.726360e-06, &tj, &tj1, &result);
   return result;
}

// Tail(S, 10, 100)
static double mannwhitneyu_utbln10n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.650000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.334008e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.522316e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.769627e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.158110e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.053650e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.242235e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.173571e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.033661e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.824732e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.084420e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.610036e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.728155e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.217130e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.340966e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.001235e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.694052e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 11)
static double mannwhitneyu_utbln11n11(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.519760e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.880694e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.200698e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.174092e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.072304e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.054773e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.506613e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.813942e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.223644e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.417416e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.499166e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.194332e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 7.369096e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.968590e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.630532e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.061000e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 12)
static double mannwhitneyu_utbln11n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.495790e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.832622e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.165420e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.987306e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.265621e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.723537e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.347406e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.353464e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.613369e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.102522e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.237709e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.665652e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.626903e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.167518e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.564455e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.047320e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 13)
static double mannwhitneyu_utbln11n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.477880e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.796242e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.138769e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.851739e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.722104e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.548304e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.176683e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.817895e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.842451e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.935870e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.421777e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.238831e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.867026e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.458255e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.306259e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.961487e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 14)
static double mannwhitneyu_utbln11n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.463683e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.766969e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.117082e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.739574e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.238865e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.350306e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.425871e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.640172e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.660633e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.879883e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.349658e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.271795e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.304544e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.024201e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.816867e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.596787e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 15)
static double mannwhitneyu_utbln11n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.452526e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.743570e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.099705e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.650612e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.858285e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.187036e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.689241e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.294360e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.072623e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.278008e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.322382e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.131558e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.305669e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.825627e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.332689e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.120973e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 30)
static double mannwhitneyu_utbln11n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.402621e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.627440e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.011333e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.224126e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.232856e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.859347e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.377381e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.756709e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.033230e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.875472e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.608399e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.102943e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.740693e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.343139e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.196878e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.658062e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 11, 100)
static double mannwhitneyu_utbln11n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.398795e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.596486e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.814761e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.085187e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.766529e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.379425e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.986351e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.214705e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.360075e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.260869e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.033307e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.727087e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.393883e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.242989e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.111928e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.898823e-09, &tj, &tj1, &result);
   return result;
}

// Tail(S, 12, 12)
static double mannwhitneyu_utbln12n12(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.472616e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.786627e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.132099e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.817523e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.570179e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.479511e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.799492e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.565350e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.530139e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.380132e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.242761e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.576269e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.018771e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.933911e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.002799e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.022048e-06, &tj, &tj1, &result);
   return result;
}

// Tail(S, 12, 13)
static double mannwhitneyu_utbln12n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.454800e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.750794e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.105988e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.684754e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.011826e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.262579e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.044492e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.478741e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.322165e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.621104e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.068753e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.468396e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.056235e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.327375e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.914877e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.784191e-04, &tj, &tj1, &result);
   return result;
}

// Tail(S, 12, 14)
static double mannwhitneyu_utbln12n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.440910e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.722404e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.085254e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.579439e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.563738e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.066730e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.129346e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.014531e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.129679e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.000909e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.996174e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.377924e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.936304e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.051098e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.025820e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 8.730585e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 12, 15)
static double mannwhitneyu_utbln12n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.430123e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.700008e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.068971e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.499725e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.250897e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.473145e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.680008e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.483350e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.766992e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.891081e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.015140e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.977756e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.707414e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.114786e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.238865e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.381445e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 12, 30)
static double mannwhitneyu_utbln12n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.380023e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.585782e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.838583e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.103394e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.834015e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.635212e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.948212e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.574169e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.747980e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.833672e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.722433e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.181038e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.206473e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.716003e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.476434e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.217700e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 12, 100)
static double mannwhitneyu_utbln12n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.700000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.374567e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.553481e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.541334e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.701907e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.414757e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.404103e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.234388e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.453762e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.311060e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.317501e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.713888e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.309583e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.019804e-08, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.224829e-09, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.349019e-08, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.893302e-08, &tj, &tj1, &result);
   return result;
}

// Tail(S, 13, 13)
static double mannwhitneyu_utbln13n13(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.541046e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.859047e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.130164e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.689719e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.950693e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.231455e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.976550e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.538455e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.245603e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.142647e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.831434e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.032483e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.488405e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.156927e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.949279e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.532700e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 13, 14)
static double mannwhitneyu_utbln13n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.525655e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.828341e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.108110e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.579552e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.488307e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.032328e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.988741e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.766394e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.388950e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.338179e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.133440e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.023518e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.110570e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.202332e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.056132e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.536323e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 13, 15)
static double mannwhitneyu_utbln13n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.513585e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.803952e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.090686e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.495310e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.160314e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.073124e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.480313e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.478239e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.140914e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.311541e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.677105e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.115464e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.578563e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.044604e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.888939e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 2.395644e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 13, 30)
static double mannwhitneyu_utbln13n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.455999e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.678434e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.995491e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.078100e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.705220e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.258739e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.671526e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.185458e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.507764e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.411446e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.044355e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.285765e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.345282e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.066940e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.962037e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.723644e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 13, 100)
static double mannwhitneyu_utbln13n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.446787e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.640804e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.671552e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.364990e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.274444e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.047440e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.161439e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.171729e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.562171e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.359762e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.275494e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.747635e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.700292e-08, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.565559e-09, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 5.005396e-09, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 3.335794e-09, &tj, &tj1, &result);
   return result;
}

// Tail(S, 14, 14)
static double mannwhitneyu_utbln14n14(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.510624e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.798584e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.087107e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.478532e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.098050e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.855986e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.409083e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.299536e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.176177e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.479417e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.812761e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -5.225872e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 4.516521e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 6.730551e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 9.237563e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.611820e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 14, 15)
static double mannwhitneyu_utbln14n15(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.498681e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.774668e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.070267e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.399348e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.807239e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.845763e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.071773e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.261698e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.011695e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.305946e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.879295e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.999439e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.904438e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.944986e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.373908e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.140794e-05, &tj, &tj1, &result);
   return result;
}

// Tail(S, 14, 30)
static double mannwhitneyu_utbln14n30(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.440378e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.649587e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.807829e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.989753e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.463646e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.586580e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -6.745917e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.635398e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.923172e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.446699e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.613892e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.214073e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.651683e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.272777e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.464988e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.109803e-07, &tj, &tj1, &result);
   return result;
}

// Tail(S, 14, 100)
static double mannwhitneyu_utbln14n100(double s) {
   double x;
   double tj;
   double tj1;
   double result;
   result = 0.0;
   x = rmin2(2 * (s - 0.000000e+00) / 3.750000e+00 - 1, 1.0);
   tj = 1.0;
   tj1 = x;
   mannwhitneyu_ucheb(x, -4.429701e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -4.610577e+00, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -9.482675e-01, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.605550e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.062151e-02, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.525154e-03, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.835983e-04, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -8.411440e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.744901e-05, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.318850e-06, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.692100e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -1.536270e-07, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -3.705888e-08, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -7.999599e-09, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, -2.908395e-09, &tj, &tj1, &result);
   mannwhitneyu_ucheb(x, 1.546923e-09, &tj, &tj1, &result);
   return result;
}

// Three-point polynomial interpolation.
static double mannwhitneyu_uninterpolate(double p1, double p2, double p3, ae_int_t n) {
   double t1;
   double t2;
   double t3;
   double t;
   double p12;
   double p23;
   double result;
   t1 = 1.0 / 15.0;
   t2 = 1.0 / 30.0;
   t3 = 1.0 / 100.0;
   t = 1.0 / n;
   p12 = ((t - t2) * p1 + (t1 - t) * p2) / (t1 - t2);
   p23 = ((t - t3) * p2 + (t2 - t) * p3) / (t2 - t3);
   result = ((t - t3) * p12 + (t1 - t) * p23) / (t1 - t3);
   return result;
}

// Tail(0, N1, N2)
static double mannwhitneyu_usigma000(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-6.76984e-01, -6.83700e-01, -6.89873e-01, n2);
   p2 = mannwhitneyu_uninterpolate(-6.83700e-01, -6.87311e-01, -6.90957e-01, n2);
   p3 = mannwhitneyu_uninterpolate(-6.89873e-01, -6.90957e-01, -6.92175e-01, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(0.75, N1, N2)
static double mannwhitneyu_usigma075(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-1.44500e+00, -1.45906e+00, -1.47063e+00, n2);
   p2 = mannwhitneyu_uninterpolate(-1.45906e+00, -1.46856e+00, -1.47644e+00, n2);
   p3 = mannwhitneyu_uninterpolate(-1.47063e+00, -1.47644e+00, -1.48100e+00, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(1.5, N1, N2)
static double mannwhitneyu_usigma150(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-2.65380e+00, -2.67352e+00, -2.69011e+00, n2);
   p2 = mannwhitneyu_uninterpolate(-2.67352e+00, -2.68591e+00, -2.69659e+00, n2);
   p3 = mannwhitneyu_uninterpolate(-2.69011e+00, -2.69659e+00, -2.70192e+00, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(2.25, N1, N2)
static double mannwhitneyu_usigma225(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-4.41465e+00, -4.42260e+00, -4.43702e+00, n2);
   p2 = mannwhitneyu_uninterpolate(-4.42260e+00, -4.41639e+00, -4.41928e+00, n2);
   p3 = mannwhitneyu_uninterpolate(-4.43702e+00, -4.41928e+00, -4.41030e+00, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(3.0, N1, N2)
static double mannwhitneyu_usigma300(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-6.89839e+00, -6.83477e+00, -6.82340e+00, n2);
   p2 = mannwhitneyu_uninterpolate(-6.83477e+00, -6.74559e+00, -6.71117e+00, n2);
   p3 = mannwhitneyu_uninterpolate(-6.82340e+00, -6.71117e+00, -6.64929e+00, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(3.33, N1, N2)
static double mannwhitneyu_usigma333(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-8.31272e+00, -8.17096e+00, -8.13125e+00, n2);
   p2 = mannwhitneyu_uninterpolate(-8.17096e+00, -8.00156e+00, -7.93245e+00, n2);
   p3 = mannwhitneyu_uninterpolate(-8.13125e+00, -7.93245e+00, -7.82502e+00, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(3.66, N1, N2)
static double mannwhitneyu_usigma367(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-9.98837e+00, -9.70844e+00, -9.62087e+00, n2);
   p2 = mannwhitneyu_uninterpolate(-9.70844e+00, -9.41156e+00, -9.28998e+00, n2);
   p3 = mannwhitneyu_uninterpolate(-9.62087e+00, -9.28998e+00, -9.11686e+00, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(4.0, N1, N2)
static double mannwhitneyu_usigma400(ae_int_t n1, ae_int_t n2) {
   double p1;
   double p2;
   double p3;
   double result;
   p1 = mannwhitneyu_uninterpolate(-1.20250e+01, -1.14911e+01, -1.13231e+01, n2);
   p2 = mannwhitneyu_uninterpolate(-1.14911e+01, -1.09927e+01, -1.07937e+01, n2);
   p3 = mannwhitneyu_uninterpolate(-1.13231e+01, -1.07937e+01, -1.05285e+01, n2);
   result = mannwhitneyu_uninterpolate(p1, p2, p3, n1);
   return result;
}

// Tail(S, N1, N2)
static double mannwhitneyu_usigma(double s, ae_int_t n1, ae_int_t n2) {
   double f0;
   double f1;
   double f2;
   double f3;
   double f4;
   double s0;
   double s1;
   double s2;
   double s3;
   double s4;
   double result;
   result = 0.0;
// N1=5, N2 = 5, 6, 7, ...
   if (imin2(n1, n2) == 5) {
      if (imax2(n1, n2) == 5) {
         result = mannwhitneyu_utbln5n5(s);
      }
      if (imax2(n1, n2) == 6) {
         result = mannwhitneyu_utbln5n6(s);
      }
      if (imax2(n1, n2) == 7) {
         result = mannwhitneyu_utbln5n7(s);
      }
      if (imax2(n1, n2) == 8) {
         result = mannwhitneyu_utbln5n8(s);
      }
      if (imax2(n1, n2) == 9) {
         result = mannwhitneyu_utbln5n9(s);
      }
      if (imax2(n1, n2) == 10) {
         result = mannwhitneyu_utbln5n10(s);
      }
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln5n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln5n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln5n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln5n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln5n15(s);
      }
      if (imax2(n1, n2) == 16) {
         result = mannwhitneyu_utbln5n16(s);
      }
      if (imax2(n1, n2) == 17) {
         result = mannwhitneyu_utbln5n17(s);
      }
      if (imax2(n1, n2) == 18) {
         result = mannwhitneyu_utbln5n18(s);
      }
      if (imax2(n1, n2) == 19) {
         result = mannwhitneyu_utbln5n19(s);
      }
      if (imax2(n1, n2) == 20) {
         result = mannwhitneyu_utbln5n20(s);
      }
      if (imax2(n1, n2) == 21) {
         result = mannwhitneyu_utbln5n21(s);
      }
      if (imax2(n1, n2) == 22) {
         result = mannwhitneyu_utbln5n22(s);
      }
      if (imax2(n1, n2) == 23) {
         result = mannwhitneyu_utbln5n23(s);
      }
      if (imax2(n1, n2) == 24) {
         result = mannwhitneyu_utbln5n24(s);
      }
      if (imax2(n1, n2) == 25) {
         result = mannwhitneyu_utbln5n25(s);
      }
      if (imax2(n1, n2) == 26) {
         result = mannwhitneyu_utbln5n26(s);
      }
      if (imax2(n1, n2) == 27) {
         result = mannwhitneyu_utbln5n27(s);
      }
      if (imax2(n1, n2) == 28) {
         result = mannwhitneyu_utbln5n28(s);
      }
      if (imax2(n1, n2) == 29) {
         result = mannwhitneyu_utbln5n29(s);
      }
      if (imax2(n1, n2) > 29) {
         f0 = mannwhitneyu_utbln5n15(s);
         f1 = mannwhitneyu_utbln5n30(s);
         f2 = mannwhitneyu_utbln5n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=6, N2 = 6, 7, 8, ...
   if (imin2(n1, n2) == 6) {
      if (imax2(n1, n2) == 6) {
         result = mannwhitneyu_utbln6n6(s);
      }
      if (imax2(n1, n2) == 7) {
         result = mannwhitneyu_utbln6n7(s);
      }
      if (imax2(n1, n2) == 8) {
         result = mannwhitneyu_utbln6n8(s);
      }
      if (imax2(n1, n2) == 9) {
         result = mannwhitneyu_utbln6n9(s);
      }
      if (imax2(n1, n2) == 10) {
         result = mannwhitneyu_utbln6n10(s);
      }
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln6n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln6n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln6n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln6n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln6n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln6n15(s);
         f1 = mannwhitneyu_utbln6n30(s);
         f2 = mannwhitneyu_utbln6n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=7, N2 = 7, 8, ...
   if (imin2(n1, n2) == 7) {
      if (imax2(n1, n2) == 7) {
         result = mannwhitneyu_utbln7n7(s);
      }
      if (imax2(n1, n2) == 8) {
         result = mannwhitneyu_utbln7n8(s);
      }
      if (imax2(n1, n2) == 9) {
         result = mannwhitneyu_utbln7n9(s);
      }
      if (imax2(n1, n2) == 10) {
         result = mannwhitneyu_utbln7n10(s);
      }
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln7n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln7n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln7n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln7n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln7n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln7n15(s);
         f1 = mannwhitneyu_utbln7n30(s);
         f2 = mannwhitneyu_utbln7n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=8, N2 = 8, 9, 10, ...
   if (imin2(n1, n2) == 8) {
      if (imax2(n1, n2) == 8) {
         result = mannwhitneyu_utbln8n8(s);
      }
      if (imax2(n1, n2) == 9) {
         result = mannwhitneyu_utbln8n9(s);
      }
      if (imax2(n1, n2) == 10) {
         result = mannwhitneyu_utbln8n10(s);
      }
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln8n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln8n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln8n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln8n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln8n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln8n15(s);
         f1 = mannwhitneyu_utbln8n30(s);
         f2 = mannwhitneyu_utbln8n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=9, N2 = 9, 10, ...
   if (imin2(n1, n2) == 9) {
      if (imax2(n1, n2) == 9) {
         result = mannwhitneyu_utbln9n9(s);
      }
      if (imax2(n1, n2) == 10) {
         result = mannwhitneyu_utbln9n10(s);
      }
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln9n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln9n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln9n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln9n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln9n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln9n15(s);
         f1 = mannwhitneyu_utbln9n30(s);
         f2 = mannwhitneyu_utbln9n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=10, N2 = 10, 11, ...
   if (imin2(n1, n2) == 10) {
      if (imax2(n1, n2) == 10) {
         result = mannwhitneyu_utbln10n10(s);
      }
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln10n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln10n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln10n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln10n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln10n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln10n15(s);
         f1 = mannwhitneyu_utbln10n30(s);
         f2 = mannwhitneyu_utbln10n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=11, N2 = 11, 12, ...
   if (imin2(n1, n2) == 11) {
      if (imax2(n1, n2) == 11) {
         result = mannwhitneyu_utbln11n11(s);
      }
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln11n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln11n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln11n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln11n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln11n15(s);
         f1 = mannwhitneyu_utbln11n30(s);
         f2 = mannwhitneyu_utbln11n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=12, N2 = 12, 13, ...
   if (imin2(n1, n2) == 12) {
      if (imax2(n1, n2) == 12) {
         result = mannwhitneyu_utbln12n12(s);
      }
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln12n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln12n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln12n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln12n15(s);
         f1 = mannwhitneyu_utbln12n30(s);
         f2 = mannwhitneyu_utbln12n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=13, N2 = 13, 14, ...
   if (imin2(n1, n2) == 13) {
      if (imax2(n1, n2) == 13) {
         result = mannwhitneyu_utbln13n13(s);
      }
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln13n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln13n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln13n15(s);
         f1 = mannwhitneyu_utbln13n30(s);
         f2 = mannwhitneyu_utbln13n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1=14, N2 = 14, 15, ...
   if (imin2(n1, n2) == 14) {
      if (imax2(n1, n2) == 14) {
         result = mannwhitneyu_utbln14n14(s);
      }
      if (imax2(n1, n2) == 15) {
         result = mannwhitneyu_utbln14n15(s);
      }
      if (imax2(n1, n2) > 15) {
         f0 = mannwhitneyu_utbln14n15(s);
         f1 = mannwhitneyu_utbln14n30(s);
         f2 = mannwhitneyu_utbln14n100(s);
         result = mannwhitneyu_uninterpolate(f0, f1, f2, imax2(n1, n2));
      }
      return result;
   }
// N1 >= 15, N2 >= 15
   if (s > 4.0) {
      s = 4.0;
   }
   if (s < 3.0) {
      s0 = 0.000000e+00;
      f0 = mannwhitneyu_usigma000(n1, n2);
      s1 = 7.500000e-01;
      f1 = mannwhitneyu_usigma075(n1, n2);
      s2 = 1.500000e+00;
      f2 = mannwhitneyu_usigma150(n1, n2);
      s3 = 2.250000e+00;
      f3 = mannwhitneyu_usigma225(n1, n2);
      s4 = 3.000000e+00;
      f4 = mannwhitneyu_usigma300(n1, n2);
      f1 = ((s - s0) * f1 - (s - s1) * f0) / (s1 - s0);
      f2 = ((s - s0) * f2 - (s - s2) * f0) / (s2 - s0);
      f3 = ((s - s0) * f3 - (s - s3) * f0) / (s3 - s0);
      f4 = ((s - s0) * f4 - (s - s4) * f0) / (s4 - s0);
      f2 = ((s - s1) * f2 - (s - s2) * f1) / (s2 - s1);
      f3 = ((s - s1) * f3 - (s - s3) * f1) / (s3 - s1);
      f4 = ((s - s1) * f4 - (s - s4) * f1) / (s4 - s1);
      f3 = ((s - s2) * f3 - (s - s3) * f2) / (s3 - s2);
      f4 = ((s - s2) * f4 - (s - s4) * f2) / (s4 - s2);
      f4 = ((s - s3) * f4 - (s - s4) * f3) / (s4 - s3);
      result = f4;
   } else {
      s0 = 3.000000e+00;
      f0 = mannwhitneyu_usigma300(n1, n2);
      s1 = 3.333333e+00;
      f1 = mannwhitneyu_usigma333(n1, n2);
      s2 = 3.666667e+00;
      f2 = mannwhitneyu_usigma367(n1, n2);
      s3 = 4.000000e+00;
      f3 = mannwhitneyu_usigma400(n1, n2);
      f1 = ((s - s0) * f1 - (s - s1) * f0) / (s1 - s0);
      f2 = ((s - s0) * f2 - (s - s2) * f0) / (s2 - s0);
      f3 = ((s - s0) * f3 - (s - s3) * f0) / (s3 - s0);
      f2 = ((s - s1) * f2 - (s - s2) * f1) / (s2 - s1);
      f3 = ((s - s1) * f3 - (s - s3) * f1) / (s3 - s1);
      f3 = ((s - s2) * f3 - (s - s3) * f2) / (s3 - s2);
      result = f3;
   }
   return result;
}

// Mann-Whitney U-test
//
// This test checks hypotheses about whether X  and  Y  are  samples  of  two
// continuous distributions of the same shape  and  same  median  or  whether
// their medians are different.
//
// The following tests are performed:
//     * two-tailed test (null hypothesis - the medians are equal)
//     * left-tailed test (null hypothesis - the median of the  first  sample
//       is greater than or equal to the median of the second sample)
//     * right-tailed test (null hypothesis - the median of the first  sample
//       is less than or equal to the median of the second sample).
//
// Requirements:
//     * the samples are independent
//     * X and Y are continuous distributions (or discrete distributions well-
//       approximating continuous distributions)
//     * distributions of X and Y have the  same  shape.  The  only  possible
//       difference is their position (i.e. the value of the median)
//     * the number of elements in each sample is not less than 5
//     * the scale of measurement should be ordinal, interval or ratio  (i.e.
//       the test could not be applied to nominal variables).
//
// The test is non-parametric and doesn't require distributions to be normal.
//
// Inputs:
//     X   -   sample 1. Array whose index goes from 0 to N-1.
//     N   -   size of the sample. N >= 5
//     Y   -   sample 2. Array whose index goes from 0 to M-1.
//     M   -   size of the sample. M >= 5
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
//
// To calculate p-values, special approximation is used. This method lets  us
// calculate p-values with satisfactory  accuracy  in  interval  [0.0001, 1].
// There is no approximation outside the [0.0001, 1] interval. Therefore,  if
// the significance level outlies this interval, the test returns 0.0001.
//
// Relative precision of approximation of p-value:
//
// N          M          Max.err.   Rms.err.
// 5..10      N..10      1.4e-02    6.0e-04
// 5..10      N..100     2.2e-02    5.3e-06
// 10..15     N..15      1.0e-02    3.2e-04
// 10..15     N..100     1.0e-02    2.2e-05
// 15..100    N..100     6.1e-03    2.7e-06
//
// For N, M > 100 accuracy checks weren't put into  practice,  but  taking  into
// account characteristics of asymptotic approximation used, precision should
// not be sharply different from the values for interval [5, 100].
//
// NOTE: P-value approximation was  optimized  for  0.0001 <= p <= 0.2500.  Thus,
//       P's outside of this interval are enforced to these bounds. Say,  you
//       may quite often get P equal to exactly 0.25 or 0.0001.
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
// API: void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
void mannwhitneyutest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t t;
   double tmp;
   ae_int_t ns;
   double u;
   double p;
   double mp;
   double s;
   double sigma;
   double mu;
   ae_int_t tiecount;
   ae_frame_make(&_frame_block);
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   NewVector(r, 0, DT_REAL);
   NewVector(c, 0, DT_INT);
   NewVector(tiesize, 0, DT_INT);
// Prepare
   if (n <= 4 || m <= 4) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      ae_frame_leave();
      return;
   }
   ns = n + m;
   ae_vector_set_length(&r, ns);
   ae_vector_set_length(&c, ns);
   for (i = 0; i < n; i++) {
      r.xR[i] = x->xR[i];
      c.xZ[i] = 0;
   }
   for (i = 0; i < m; i++) {
      r.xR[n + i] = y->xR[i];
      c.xZ[n + i] = 1;
   }
// sort {R, C}
   if (ns != 1) {
      i = 2;
      do {
         t = i;
         while (t != 1) {
            k = t / 2;
            if (r.xR[k - 1] >= r.xR[t - 1]) {
               t = 1;
            } else {
               swapr(&r.xR[k - 1], &r.xR[t - 1]);
               swapi(&c.xZ[k - 1], &c.xZ[t - 1]);
               t = k;
            }
         }
         i++;
      } while (i <= ns);
      i = ns - 1;
      do {
         swapr(&r.xR[i], &r.xR[0]);
         swapi(&c.xZ[i], &c.xZ[0]);
         t = 1;
         while (t != 0) {
            k = 2 * t;
            if (k > i) {
               t = 0;
            } else {
               if (k < i) {
                  if (r.xR[k] > r.xR[k - 1]) {
                     k++;
                  }
               }
               if (r.xR[t - 1] >= r.xR[k - 1]) {
                  t = 0;
               } else {
                  swapr(&r.xR[k - 1], &r.xR[t - 1]);
                  swapi(&c.xZ[k - 1], &c.xZ[t - 1]);
                  t = k;
               }
            }
         }
         i--;
      } while (i >= 1);
   }
// compute tied ranks
   i = 0;
   tiecount = 0;
   ae_vector_set_length(&tiesize, ns);
   while (i < ns) {
      j = i + 1;
      while (j < ns) {
         if (r.xR[j] != r.xR[i]) {
            break;
         }
         j++;
      }
      for (k = i; k < j; k++) {
         r.xR[k] = 1 + (double)(i + j - 1) / 2.0;
      }
      tiesize.xZ[tiecount] = j - i;
      tiecount++;
      i = j;
   }
// Compute U
   u = 0.0;
   for (i = 0; i < ns; i++) {
      if (c.xZ[i] == 0) {
         u += r.xR[i];
      }
   }
   u = (double)n * m + 0.5 * n * (n + 1) - u;
// Result
   mu = (double)n * m / 2;
   tmp = ns * (sqr((double)ns) - 1) / 12;
   for (i = 0; i < tiecount; i++) {
      tmp -= tiesize.xZ[i] * (sqr((double)tiesize.xZ[i]) - 1) / 12;
   }
   sigma = sqrt((double)n * m / ns / (ns - 1) * tmp);
   s = (u - mu) / sigma;
   if (s <= 0.0) {
      p = exp(mannwhitneyu_usigma(-(u - mu) / sigma, n, m));
      mp = 1 - exp(mannwhitneyu_usigma(-(u - 1 - mu) / sigma, n, m));
   } else {
      mp = exp(mannwhitneyu_usigma((u - mu) / sigma, n, m));
      p = 1 - exp(mannwhitneyu_usigma((u + 1 - mu) / sigma, n, m));
   }
   *lefttail = rboundval(rmax2(mp, 1.0E-4), 0.0001, 0.2500);
   *righttail = rboundval(rmax2(p, 1.0E-4), 0.0001, 0.2500);
   *bothtails = 2 * rmin2(*lefttail, *righttail);
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::mannwhitneyutest(ConstT(ae_vector, x), n, ConstT(ae_vector, y), m, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === STEST Package ===
// Depends on: (SpecialFunctions) BINOMIALDISTR
namespace alglib_impl {
// Sign test
//
// This test checks three hypotheses about the median of  the  given  sample.
// The following tests are performed:
//     * two-tailed test (null hypothesis - the median is equal to the  given
//       value)
//     * left-tailed test (null hypothesis - the median is  greater  than  or
//       equal to the given value)
//     * right-tailed test (null hypothesis - the  median  is  less  than  or
//       equal to the given value)
//
// Requirements:
//     * the scale of measurement should be ordinal, interval or ratio  (i.e.
//       the test could not be applied to nominal variables).
//
// The test is non-parametric and doesn't require distribution X to be normal
//
// Inputs:
//     X       -   sample. Array whose index goes from 0 to N-1.
//     N       -   size of the sample.
//     Median  -   assumed median value.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
//
// While   calculating   p-values   high-precision   binomial    distribution
// approximation is used, so significance levels have about 15 exact digits.
// ALGLIB: Copyright 08.09.2006 by Sergey Bochkanov
// API: void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail);
void onesamplesigntest(RVector *x, ae_int_t n, double median, double *bothtails, double *lefttail, double *righttail) {
   ae_int_t i;
   ae_int_t gtcnt;
   ae_int_t necnt;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   if (n <= 1) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Calculate:
// GTCnt - count of x[i]>Median
// NECnt - count of x[i] != Median
   gtcnt = 0;
   necnt = 0;
   for (i = 0; i < n; i++) {
      if (x->xR[i] > median) {
         gtcnt++;
      }
      if (x->xR[i] != median) {
         necnt++;
      }
   }
   if (necnt == 0) {
   // all x[i] are equal to Median.
   // So we can conclude that Median is a true median :)
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
   *bothtails = rmin2(2 * binomialdistribution(imin2(gtcnt, necnt - gtcnt), necnt, 0.5), 1.0);
   *lefttail = binomialdistribution(gtcnt, necnt, 0.5);
   *righttail = binomialcdistribution(gtcnt - 1, necnt, 0.5);
}
} // end of namespace alglib_impl

namespace alglib {
void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::onesamplesigntest(ConstT(ae_vector, x), n, median, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === STUDENTTTESTS Package ===
// Depends on: (SpecialFunctions) STUDENTTDISTR
namespace alglib_impl {
// One-sample t-test
//
// This test checks three hypotheses about the mean of the given sample.  The
// following tests are performed:
//     * two-tailed test (null hypothesis - the mean is equal  to  the  given
//       value)
//     * left-tailed test (null hypothesis - the  mean  is  greater  than  or
//       equal to the given value)
//     * right-tailed test (null hypothesis - the mean is less than or  equal
//       to the given value).
//
// The test is based on the assumption that  a  given  sample  has  a  normal
// distribution and  an  unknown  dispersion.  If  the  distribution  sharply
// differs from normal, the test will work incorrectly.
//
// Inputs:
//     X       -   sample. Array whose index goes from 0 to N-1.
//     N       -   size of sample, N >= 0
//     Mean    -   assumed value of the mean.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
//
// NOTE: this function correctly handles degenerate cases:
//       * when N=0, all p-values are set to 1.0
//       * when variance of X[] is exactly zero, p-values are set
//         to 1.0 or 0.0, depending on difference between sample mean and
//         value of mean being tested.
// ALGLIB: Copyright 08.09.2006 by Sergey Bochkanov
// API: void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail);
void studentttest1(RVector *x, ae_int_t n, double mean, double *bothtails, double *lefttail, double *righttail) {
   ae_int_t i;
   double xmean;
   double x0;
   double v;
   bool samex;
   double xvariance;
   double xstddev;
   double v1;
   double v2;
   double stat;
   double s;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   if (n <= 0) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Mean
   xmean = 0.0;
   x0 = x->xR[0];
   samex = true;
   for (i = 0; i < n; i++) {
      v = x->xR[i];
      xmean += v;
      samex = samex && v == x0;
   }
   if (samex) {
      xmean = x0;
   } else {
      xmean /= n;
   }
// Variance (using corrected two-pass algorithm)
   xvariance = 0.0;
   xstddev = 0.0;
   if (n != 1 && !samex) {
      v1 = 0.0;
      for (i = 0; i < n; i++) {
         v1 += sqr(x->xR[i] - xmean);
      }
      v2 = 0.0;
      for (i = 0; i < n; i++) {
         v2 += x->xR[i] - xmean;
      }
      v2 = sqr(v2) / n;
      xvariance = (v1 - v2) / (n - 1);
      if (xvariance < 0.0) {
         xvariance = 0.0;
      }
      xstddev = sqrt(xvariance);
   }
   if (xstddev == 0.0) {
      if (xmean == mean) {
         *bothtails = 1.0;
      } else {
         *bothtails = 0.0;
      }
      if (xmean >= mean) {
         *lefttail = 1.0;
      } else {
         *lefttail = 0.0;
      }
      if (xmean <= mean) {
         *righttail = 1.0;
      } else {
         *righttail = 0.0;
      }
      return;
   }
// Statistic
   stat = (xmean - mean) / (xstddev / sqrt((double)n));
   s = studenttdistribution(n - 1, stat);
   *bothtails = 2 * rmin2(s, 1 - s);
   *lefttail = s;
   *righttail = 1 - s;
}

// Two-sample pooled test
//
// This test checks three hypotheses about the mean of the given samples. The
// following tests are performed:
//     * two-tailed test (null hypothesis - the means are equal)
//     * left-tailed test (null hypothesis - the mean of the first sample  is
//       greater than or equal to the mean of the second sample)
//     * right-tailed test (null hypothesis - the mean of the first sample is
//       less than or equal to the mean of the second sample).
//
// Test is based on the following assumptions:
//     * given samples have normal distributions
//     * dispersions are equal
//     * samples are independent.
//
// Inputs:
//     X       -   sample 1. Array whose index goes from 0 to N-1.
//     N       -   size of sample.
//     Y       -   sample 2. Array whose index goes from 0 to M-1.
//     M       -   size of sample.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
//
// NOTE: this function correctly handles degenerate cases:
//       * when N=0 or M=0, all p-values are set to 1.0
//       * when both samples has exactly zero variance, p-values are set
//         to 1.0 or 0.0, depending on difference between means.
// ALGLIB: Copyright 18.09.2006 by Sergey Bochkanov
// API: void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
void studentttest2(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail) {
   ae_int_t i;
   bool samex;
   bool samey;
   double x0;
   double y0;
   double xmean;
   double ymean;
   double v;
   double stat;
   double s;
   double p;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   if (n <= 0 || m <= 0) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Mean
   xmean = 0.0;
   x0 = x->xR[0];
   samex = true;
   for (i = 0; i < n; i++) {
      v = x->xR[i];
      xmean += v;
      samex = samex && v == x0;
   }
   if (samex) {
      xmean = x0;
   } else {
      xmean /= n;
   }
   ymean = 0.0;
   y0 = y->xR[0];
   samey = true;
   for (i = 0; i < m; i++) {
      v = y->xR[i];
      ymean += v;
      samey = samey && v == y0;
   }
   if (samey) {
      ymean = y0;
   } else {
      ymean /= m;
   }
// S
   s = 0.0;
   if (n + m > 2) {
      for (i = 0; i < n; i++) {
         s += sqr(x->xR[i] - xmean);
      }
      for (i = 0; i < m; i++) {
         s += sqr(y->xR[i] - ymean);
      }
      s = sqrt(s * (1.0 / n + 1.0 / m) / (n + m - 2));
   }
   if (s == 0.0) {
      if (xmean == ymean) {
         *bothtails = 1.0;
      } else {
         *bothtails = 0.0;
      }
      if (xmean >= ymean) {
         *lefttail = 1.0;
      } else {
         *lefttail = 0.0;
      }
      if (xmean <= ymean) {
         *righttail = 1.0;
      } else {
         *righttail = 0.0;
      }
      return;
   }
// Statistic
   stat = (xmean - ymean) / s;
   p = studenttdistribution(n + m - 2, stat);
   *bothtails = 2 * rmin2(p, 1 - p);
   *lefttail = p;
   *righttail = 1 - p;
}

// Two-sample unpooled test
//
// This test checks three hypotheses about the mean of the given samples. The
// following tests are performed:
//     * two-tailed test (null hypothesis - the means are equal)
//     * left-tailed test (null hypothesis - the mean of the first sample  is
//       greater than or equal to the mean of the second sample)
//     * right-tailed test (null hypothesis - the mean of the first sample is
//       less than or equal to the mean of the second sample).
//
// Test is based on the following assumptions:
//     * given samples have normal distributions
//     * samples are independent.
// Equality of variances is NOT required.
//
// Inputs:
//     X - sample 1. Array whose index goes from 0 to N-1.
//     N - size of the sample.
//     Y - sample 2. Array whose index goes from 0 to M-1.
//     M - size of the sample.
//
// Outputs:
//     BothTails   -   p-value for two-tailed test.
//                     If BothTails is less than the given significance level
//                     the null hypothesis is rejected.
//     LeftTail    -   p-value for left-tailed test.
//                     If LeftTail is less than the given significance level,
//                     the null hypothesis is rejected.
//     RightTail   -   p-value for right-tailed test.
//                     If RightTail is less than the given significance level
//                     the null hypothesis is rejected.
//
// NOTE: this function correctly handles degenerate cases:
//       * when N=0 or M=0, all p-values are set to 1.0
//       * when both samples has zero variance, p-values are set
//         to 1.0 or 0.0, depending on difference between means.
//       * when only one sample has zero variance, test reduces to 1-sample
//         version.
// ALGLIB: Copyright 18.09.2006 by Sergey Bochkanov
// API: void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
void unequalvariancettest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail) {
   ae_int_t i;
   bool samex;
   bool samey;
   double x0;
   double y0;
   double xmean;
   double ymean;
   double xvar;
   double yvar;
   double v;
   double df;
   double p;
   double stat;
   double c;
   *bothtails = 0;
   *lefttail = 0;
   *righttail = 0;
   if (n <= 0 || m <= 0) {
      *bothtails = 1.0;
      *lefttail = 1.0;
      *righttail = 1.0;
      return;
   }
// Mean
   xmean = 0.0;
   x0 = x->xR[0];
   samex = true;
   for (i = 0; i < n; i++) {
      v = x->xR[i];
      xmean += v;
      samex = samex && v == x0;
   }
   if (samex) {
      xmean = x0;
   } else {
      xmean /= n;
   }
   ymean = 0.0;
   y0 = y->xR[0];
   samey = true;
   for (i = 0; i < m; i++) {
      v = y->xR[i];
      ymean += v;
      samey = samey && v == y0;
   }
   if (samey) {
      ymean = y0;
   } else {
      ymean /= m;
   }
// Variance (using corrected two-pass algorithm)
   xvar = 0.0;
   if (n >= 2 && !samex) {
      for (i = 0; i < n; i++) {
         xvar += sqr(x->xR[i] - xmean);
      }
      xvar /= n - 1;
   }
   yvar = 0.0;
   if (m >= 2 && !samey) {
      for (i = 0; i < m; i++) {
         yvar += sqr(y->xR[i] - ymean);
      }
      yvar /= m - 1;
   }
// Handle different special cases
// (one or both variances are zero).
   if (xvar == 0.0 && yvar == 0.0) {
      if (xmean == ymean) {
         *bothtails = 1.0;
      } else {
         *bothtails = 0.0;
      }
      if (xmean >= ymean) {
         *lefttail = 1.0;
      } else {
         *lefttail = 0.0;
      }
      if (xmean <= ymean) {
         *righttail = 1.0;
      } else {
         *righttail = 0.0;
      }
      return;
   }
   if (xvar == 0.0) {
   // X is constant, unpooled 2-sample test reduces to 1-sample test.
   //
   // NOTE: right-tail and left-tail must be passed to 1-sample
   //       t-test in reverse order because we reverse order of
   //       of samples.
      studentttest1(y, m, xmean, bothtails, righttail, lefttail);
      return;
   }
   if (yvar == 0.0) {
   // Y is constant, unpooled 2-sample test reduces to 1-sample test.
      studentttest1(x, n, ymean, bothtails, lefttail, righttail);
      return;
   }
// Statistic
   stat = (xmean - ymean) / sqrt(xvar / n + yvar / m);
   c = xvar / n / (xvar / n + yvar / m);
   df = (double)(n - 1) * (m - 1) / ((m - 1) * sqr(c) + (n - 1) * sqr(1 - c));
   if (stat > 0.0) {
      p = 1 - 0.5 * incompletebeta(df / 2, 0.5, df / (df + sqr(stat)));
   } else {
      p = 0.5 * incompletebeta(df / 2, 0.5, df / (df + sqr(stat)));
   }
   *bothtails = 2 * rmin2(p, 1 - p);
   *lefttail = p;
   *righttail = 1 - p;
}
} // end of namespace alglib_impl

namespace alglib {
void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::studentttest1(ConstT(ae_vector, x), n, mean, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}

void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::studentttest2(ConstT(ae_vector, x), n, ConstT(ae_vector, y), m, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}

void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::unequalvariancettest(ConstT(ae_vector, x), n, ConstT(ae_vector, y), m, &bothtails, &lefttail, &righttail);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib
