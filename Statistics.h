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
void samplemoments(RVector *x, ae_int_t n, double *mean, double *variance, double *skewness, double *kurtosis, ae_state *_state);
double samplemean(RVector *x, ae_int_t n, ae_state *_state);
double samplevariance(RVector *x, ae_int_t n, ae_state *_state);
double sampleskewness(RVector *x, ae_int_t n, ae_state *_state);
double samplekurtosis(RVector *x, ae_int_t n, ae_state *_state);
void sampleadev(RVector *x, ae_int_t n, double *adev, ae_state *_state);
void samplemedian(RVector *x, ae_int_t n, double *median, ae_state *_state);
void samplepercentile(RVector *x, ae_int_t n, double p, double *v, ae_state *_state);
double cov2(RVector *x, RVector *y, ae_int_t n, ae_state *_state);
double pearsoncorr2(RVector *x, RVector *y, ae_int_t n, ae_state *_state);
double spearmancorr2(RVector *x, RVector *y, ae_int_t n, ae_state *_state);
void covm(RMatrix *x, ae_int_t n, ae_int_t m, RMatrix *c, ae_state *_state);
void pearsoncorrm(RMatrix *x, ae_int_t n, ae_int_t m, RMatrix *c, ae_state *_state);
void spearmancorrm(RMatrix *x, ae_int_t n, ae_int_t m, RMatrix *c, ae_state *_state);
void covm2(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix *c, ae_state *_state);
void pearsoncorrm2(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix *c, ae_state *_state);
void spearmancorrm2(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t m1, ae_int_t m2, RMatrix *c, ae_state *_state);
void rankdata(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
bool _trypexec_rankdata(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
void rankdatacentered(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
bool _trypexec_rankdatacentered(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
double pearsoncorrelation(RVector *x, RVector *y, ae_int_t n, ae_state *_state);
double spearmanrankcorrelation(RVector *x, RVector *y, ae_int_t n, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis, const xparams _xparams = xdefault);
void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis, const xparams _xparams = xdefault);
double samplemean(const real_1d_array &x, const ae_int_t n, const xparams _xparams = xdefault);
double samplemean(const real_1d_array &x, const xparams _xparams = xdefault);
double samplevariance(const real_1d_array &x, const ae_int_t n, const xparams _xparams = xdefault);
double samplevariance(const real_1d_array &x, const xparams _xparams = xdefault);
double sampleskewness(const real_1d_array &x, const ae_int_t n, const xparams _xparams = xdefault);
double sampleskewness(const real_1d_array &x, const xparams _xparams = xdefault);
double samplekurtosis(const real_1d_array &x, const ae_int_t n, const xparams _xparams = xdefault);
double samplekurtosis(const real_1d_array &x, const xparams _xparams = xdefault);
void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev, const xparams _xparams = xdefault);
void sampleadev(const real_1d_array &x, double &adev, const xparams _xparams = xdefault);
void samplemedian(const real_1d_array &x, const ae_int_t n, double &median, const xparams _xparams = xdefault);
void samplemedian(const real_1d_array &x, double &median, const xparams _xparams = xdefault);
void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v, const xparams _xparams = xdefault);
void samplepercentile(const real_1d_array &x, const double p, double &v, const xparams _xparams = xdefault);
double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = xdefault);
double cov2(const real_1d_array &x, const real_1d_array &y, const xparams _xparams = xdefault);
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = xdefault);
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const xparams _xparams = xdefault);
double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = xdefault);
double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const xparams _xparams = xdefault);
void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c, const xparams _xparams = xdefault);
void covm(const real_2d_array &x, real_2d_array &c, const xparams _xparams = xdefault);
void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c, const xparams _xparams = xdefault);
void pearsoncorrm(const real_2d_array &x, real_2d_array &c, const xparams _xparams = xdefault);
void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c, const xparams _xparams = xdefault);
void spearmancorrm(const real_2d_array &x, real_2d_array &c, const xparams _xparams = xdefault);
void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c, const xparams _xparams = xdefault);
void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c, const xparams _xparams = xdefault);
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c, const xparams _xparams = xdefault);
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c, const xparams _xparams = xdefault);
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c, const xparams _xparams = xdefault);
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c, const xparams _xparams = xdefault);
void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const xparams _xparams = xdefault);
void rankdata(real_2d_array &xy, const xparams _xparams = xdefault);
void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const xparams _xparams = xdefault);
void rankdatacentered(real_2d_array &xy, const xparams _xparams = xdefault);
double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = xdefault);
double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = xdefault);
} // end of namespace alglib

// === CORRELATIONTESTS Package ===
// Depends on: (SpecialFunctions) STUDENTTDISTR
// Depends on: BASESTAT
namespace alglib_impl {
void pearsoncorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void spearmanrankcorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
} // end of namespace alglib

// === JARQUEBERA Package ===
namespace alglib_impl {
void jarqueberatest(RVector *x, ae_int_t n, double *p, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p, const xparams _xparams = xdefault);
} // end of namespace alglib

// === VARIANCETESTS Package ===
// Depends on: (SpecialFunctions) FDISTR, CHISQUAREDISTR
namespace alglib_impl {
void ftest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void onesamplevariancetest(RVector *x, ae_int_t n, double variance, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
} // end of namespace alglib

// === WSR Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
void wilcoxonsignedranktest(RVector *x, ae_int_t n, double e, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void wilcoxonsignedranktest(const real_1d_array &x, const ae_int_t n, const double e, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MANNWHITNEYU Package ===
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
void mannwhitneyutest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
} // end of namespace alglib

// === STEST Package ===
// Depends on: (SpecialFunctions) BINOMIALDISTR
namespace alglib_impl {
void onesamplesigntest(RVector *x, ae_int_t n, double median, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
} // end of namespace alglib

// === STUDENTTTESTS Package ===
// Depends on: (SpecialFunctions) STUDENTTDISTR
namespace alglib_impl {
void studentttest1(RVector *x, ae_int_t n, double mean, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void studentttest2(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void unequalvariancettest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = xdefault);
} // end of namespace alglib

#endif // OnceOnly
