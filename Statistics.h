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

// Declarations for the computational core: datatypes.
namespace alglib_impl {
// === BASESTAT Package ===

// === CORRELATIONTESTS Package ===

// === JARQUEBERA Package ===

// === VARIANCETESTS Package ===

// === WSR Package ===

// === MANNWHITNEYU Package ===

// === STEST Package ===

// === STUDENTTTESTS Package ===
} // end of namespace alglib_impl

// Declarations for the C++ interface.
namespace alglib {
// === BASESTAT Package ===

// === CORRELATIONTESTS Package ===

// === JARQUEBERA Package ===

// === VARIANCETESTS Package ===

// === WSR Package ===

// === MANNWHITNEYU Package ===

// === STEST Package ===

// === STUDENTTTESTS Package ===

// === BASESTAT Package ===
void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis, const xparams _xparams = alglib::xdefault);
void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis, const xparams _xparams = alglib::xdefault);
double samplemean(const real_1d_array &x, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double samplemean(const real_1d_array &x, const xparams _xparams = alglib::xdefault);
double samplevariance(const real_1d_array &x, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double samplevariance(const real_1d_array &x, const xparams _xparams = alglib::xdefault);
double sampleskewness(const real_1d_array &x, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double sampleskewness(const real_1d_array &x, const xparams _xparams = alglib::xdefault);
double samplekurtosis(const real_1d_array &x, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double samplekurtosis(const real_1d_array &x, const xparams _xparams = alglib::xdefault);
void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev, const xparams _xparams = alglib::xdefault);
void sampleadev(const real_1d_array &x, double &adev, const xparams _xparams = alglib::xdefault);
void samplemedian(const real_1d_array &x, const ae_int_t n, double &median, const xparams _xparams = alglib::xdefault);
void samplemedian(const real_1d_array &x, double &median, const xparams _xparams = alglib::xdefault);
void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v, const xparams _xparams = alglib::xdefault);
void samplepercentile(const real_1d_array &x, const double p, double &v, const xparams _xparams = alglib::xdefault);
double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double cov2(const real_1d_array &x, const real_1d_array &y, const xparams _xparams = alglib::xdefault);
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const xparams _xparams = alglib::xdefault);
double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const xparams _xparams = alglib::xdefault);
void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void covm(const real_2d_array &x, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void pearsoncorrm(const real_2d_array &x, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void spearmancorrm(const real_2d_array &x, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c, const xparams _xparams = alglib::xdefault);
void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const xparams _xparams = alglib::xdefault);
void rankdata(real_2d_array &xy, const xparams _xparams = alglib::xdefault);
void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const xparams _xparams = alglib::xdefault);
void rankdatacentered(real_2d_array &xy, const xparams _xparams = alglib::xdefault);
double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = alglib::xdefault);
double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const xparams _xparams = alglib::xdefault);

// === CORRELATIONTESTS Package ===
void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);
void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);

// === JARQUEBERA Package ===
void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p, const xparams _xparams = alglib::xdefault);

// === VARIANCETESTS Package ===
void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);
void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);

// === WSR Package ===
void wilcoxonsignedranktest(const real_1d_array &x, const ae_int_t n, const double e, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);

// === MANNWHITNEYU Package ===
void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);

// === STEST Package ===
void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);

// === STUDENTTTESTS Package ===
void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);
void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);
void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail, const xparams _xparams = alglib::xdefault);
} // end of namespace alglib

// Declarations for the computational core: functions.
namespace alglib_impl {
// === BASESTAT Package ===
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

// === CORRELATIONTESTS Package ===
void pearsoncorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void spearmanrankcorrelationsignificance(double r, ae_int_t n, double *bothtails, double *lefttail, double *righttail, ae_state *_state);

// === JARQUEBERA Package ===
void jarqueberatest(RVector *x, ae_int_t n, double *p, ae_state *_state);

// === VARIANCETESTS Package ===
void ftest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void onesamplevariancetest(RVector *x, ae_int_t n, double variance, double *bothtails, double *lefttail, double *righttail, ae_state *_state);

// === WSR Package ===
void wilcoxonsignedranktest(RVector *x, ae_int_t n, double e, double *bothtails, double *lefttail, double *righttail, ae_state *_state);

// === MANNWHITNEYU Package ===
void mannwhitneyutest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);

// === STEST Package ===
void onesamplesigntest(RVector *x, ae_int_t n, double median, double *bothtails, double *lefttail, double *righttail, ae_state *_state);

// === STUDENTTTESTS Package ===
void studentttest1(RVector *x, ae_int_t n, double mean, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void studentttest2(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
void unequalvariancettest(RVector *x, ae_int_t n, RVector *y, ae_int_t m, double *bothtails, double *lefttail, double *righttail, ae_state *_state);
} // end of namespace alglib_impl

#endif // OnceOnly
