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
#ifndef OnceOnlyStatistics_h
#define OnceOnlyStatistics_h

#include "LinAlg.h"
#include "SpecialFunctions.h"

// === BASESTAT Package ===
// Depends on: (AlgLibInternal) BASICSTATOPS
// Depends on: (LinAlg) ABLAS
namespace alglib_impl {
void samplemoments(RVector x, ae_int_t n, double *mean, double *variance, double *skewness, double *kurtosis);
double samplemean(RVector x, ae_int_t n);
double samplevariance(RVector x, ae_int_t n);
double sampleskewness(RVector x, ae_int_t n);
double samplekurtosis(RVector x, ae_int_t n);
void sampleadev(RVector x, ae_int_t n, double *adev);
void samplemedian(RVector x, ae_int_t n, double *median);
void samplepercentile(RVector x, ae_int_t n, double p, double *v);
double cov2(RVector x, RVector y, ae_int_t n);
double pearsoncorr2(RVector x, RVector y, ae_int_t n);
double spearmancorr2(RVector x, RVector y, ae_int_t n);
void covm(RMatrix x, ae_int_t n, ae_int_t m, RMatrix c);
void pearsoncorrm(RMatrix x, ae_int_t n, ae_int_t m, RMatrix c);
void spearmancorrm(RMatrix x, ae_int_t n, ae_int_t m, RMatrix c);
void covm2(RMatrix x, RMatrix y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix c);
void pearsoncorrm2(RMatrix x, RMatrix y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix c);
void spearmancorrm2(RMatrix x, RMatrix y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix c);
void rankdata(RMatrix xy, ae_int_t npoints, ae_int_t nfeatures);
void rankdatacentered(RMatrix xy, ae_int_t npoints, ae_int_t nfeatures);
double pearsoncorrelation(RVector x, RVector y, ae_int_t n);
double spearmanrankcorrelation(RVector x, RVector y, ae_int_t n);
} // end of namespace alglib_impl

namespace alglib {
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
void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis);
void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis);

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
//
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
double samplemean(const real_1d_array &x, const ae_int_t n);
double samplemean(const real_1d_array &x);

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
//
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
double samplevariance(const real_1d_array &x, const ae_int_t n);
double samplevariance(const real_1d_array &x);

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
//
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
double sampleskewness(const real_1d_array &x, const ae_int_t n);
double sampleskewness(const real_1d_array &x);

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
//
// ALGLIB: Copyright 06.09.2006 by Sergey Bochkanov
double samplekurtosis(const real_1d_array &x, const ae_int_t n);
double samplekurtosis(const real_1d_array &x);

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
void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev);
void sampleadev(const real_1d_array &x, double &adev);

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
void samplemedian(const real_1d_array &x, const ae_int_t n, double &median);
void samplemedian(const real_1d_array &x, double &median);

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
void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v);
void samplepercentile(const real_1d_array &x, const double p, double &v);

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
double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double cov2(const real_1d_array &x, const real_1d_array &y);

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
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y);

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
double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double spearmancorr2(const real_1d_array &x, const real_1d_array &y);

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
void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void covm(const real_2d_array &x, real_2d_array &c);

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
void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void pearsoncorrm(const real_2d_array &x, real_2d_array &c);

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
void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void spearmancorrm(const real_2d_array &x, real_2d_array &c);

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
void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);

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
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);

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
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);

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
void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
void rankdata(real_2d_array &xy);

// This function replaces data in XY by their CENTERED ranks:
// * XY is processed row-by-row
// * rows are processed separately
// * tied data are correctly handled (tied ranks are calculated)
// * centered ranks are just usual ranks, but centered in such way  that  sum
//   of within-row values is equal to 0.0.
// * centering is performed by subtracting mean from each row, i.e it changes
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
void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
void rankdatacentered(real_2d_array &xy);

// Obsolete function, we recommend to use PearsonCorr2().
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);

// Obsolete function, we recommend to use SpearmanCorr2().
// ALGLIB: Copyright 09.04.2007 by Sergey Bochkanov
double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
} // end of namespace alglib

// === WSR Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
void wilcoxonsignedranktest(RVector x, ae_int_t n, double e, double *bothtails, double *lefttail, double *righttail);
} // end of namespace alglib_impl

namespace alglib {
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
void wilcoxonsignedranktest(const real_1d_array &x, const ae_int_t n, const double e, double &bothtails, double &lefttail, double &righttail);
} // end of namespace alglib

// === STEST Package ===
// Depends on: (SpecialFunctions) BINOMIALDISTR
namespace alglib_impl {
void onesamplesigntest(RVector x, ae_int_t n, double median, double *bothtails, double *lefttail, double *righttail);
} // end of namespace alglib_impl

namespace alglib {
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
void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail);
} // end of namespace alglib

// === CORRELATIONTESTS Package ===
// Depends on: (SpecialFunctions) STUDENTTDISTR
// Depends on: BASESTAT
namespace alglib_impl {
void pearsoncorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail);
void spearmanrankcorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail);
} // end of namespace alglib_impl

namespace alglib {
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
void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);

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
void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);
} // end of namespace alglib

// === STUDENTTTESTS Package ===
// Depends on: (SpecialFunctions) STUDENTTDISTR
namespace alglib_impl {
void studentttest1(RVector x, ae_int_t n, double mean, double *bothtails, double *lefttail, double *righttail);
void studentttest2(RVector x, ae_int_t n, RVector y, ae_int_t m, double *bothtails, double *lefttail, double *righttail);
void unequalvariancettest(RVector x, ae_int_t n, RVector y, ae_int_t m, double *bothtails, double *lefttail, double *righttail);
} // end of namespace alglib_impl

namespace alglib {
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
//
// ALGLIB: Copyright 08.09.2006 by Sergey Bochkanov
void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail);

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
void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);

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
void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
} // end of namespace alglib

// === MANNWHITNEYU Package ===
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
void mannwhitneyutest(RVector x, ae_int_t n, RVector y, ae_int_t m, double *bothtails, double *lefttail, double *righttail);
} // end of namespace alglib_impl

namespace alglib {
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
void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
} // end of namespace alglib

// === JARQUEBERA Package ===
namespace alglib_impl {
void jarqueberatest(RVector x, ae_int_t n, double *p);
} // end of namespace alglib_impl

namespace alglib {
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
void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p);
} // end of namespace alglib

// === VARIANCETESTS Package ===
// Depends on: (SpecialFunctions) FDISTR
// Depends on: (SpecialFunctions) CHISQUAREDISTR
namespace alglib_impl {
void ftest(RVector x, ae_int_t n, RVector y, ae_int_t m, double *bothtails, double *lefttail, double *righttail);
void onesamplevariancetest(RVector x, ae_int_t n, double variance, double *bothtails, double *lefttail, double *righttail);
} // end of namespace alglib_impl

namespace alglib {
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
void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);

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
void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail);
} // end of namespace alglib

#endif // OnceOnly
