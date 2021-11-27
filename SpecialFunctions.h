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
#ifndef OnceOnlySpecialFunctions_h
#define OnceOnlySpecialFunctions_h

#include "AlgLibMisc.h"

// === GAMMAFUNC Package ===
namespace alglib_impl {
double gammafunction(double x, ae_state *_state);
double lngamma(double x, double *sgngam, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double gammafunction(const double x, const xparams _xparams = xdefault);
double lngamma(const double x, double &sgngam, const xparams _xparams = xdefault);
} // end of namespace alglib

// === NORMALDISTR Package ===
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
double errorfunction(double x, ae_state *_state);
double errorfunctionc(double x, ae_state *_state);
double normaldistribution(double x, ae_state *_state);
double normalpdf(double x, ae_state *_state);
double normalcdf(double x, ae_state *_state);
double inverf(double e, ae_state *_state);
double invnormaldistribution(double y0, ae_state *_state);
double invnormalcdf(double y0, ae_state *_state);
double bivariatenormalpdf(double x, double y, double rho, ae_state *_state);
double bivariatenormalcdf(double x, double y, double rho, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double errorfunction(const double x, const xparams _xparams = xdefault);
double errorfunctionc(const double x, const xparams _xparams = xdefault);
double normaldistribution(const double x, const xparams _xparams = xdefault);
double normalpdf(const double x, const xparams _xparams = xdefault);
double normalcdf(const double x, const xparams _xparams = xdefault);
double inverf(const double e, const xparams _xparams = xdefault);
double invnormaldistribution(const double y0, const xparams _xparams = xdefault);
double invnormalcdf(const double y0, const xparams _xparams = xdefault);
double bivariatenormalpdf(const double x, const double y, const double rho, const xparams _xparams = xdefault);
double bivariatenormalcdf(const double x, const double y, const double rho, const xparams _xparams = xdefault);
} // end of namespace alglib

// === IBETAF Package ===
// Depends on: GAMMAFUNC, NORMALDISTR
namespace alglib_impl {
double incompletebeta(double a, double b, double x, ae_state *_state);
double invincompletebeta(double a, double b, double y, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double incompletebeta(const double a, const double b, const double x, const xparams _xparams = xdefault);
double invincompletebeta(const double a, const double b, const double y, const xparams _xparams = xdefault);
} // end of namespace alglib

// === STUDENTTDISTR Package ===
// Depends on: IBETAF
namespace alglib_impl {
double studenttdistribution(ae_int_t k, double t, ae_state *_state);
double invstudenttdistribution(ae_int_t k, double p, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double studenttdistribution(const ae_int_t k, const double t, const xparams _xparams = xdefault);
double invstudenttdistribution(const ae_int_t k, const double p, const xparams _xparams = xdefault);
} // end of namespace alglib

// === FDISTR Package ===
// Depends on: IBETAF
namespace alglib_impl {
double fdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
double fcdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
double invfdistribution(ae_int_t a, ae_int_t b, double y, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double fdistribution(const ae_int_t a, const ae_int_t b, const double x, const xparams _xparams = xdefault);
double fcdistribution(const ae_int_t a, const ae_int_t b, const double x, const xparams _xparams = xdefault);
double invfdistribution(const ae_int_t a, const ae_int_t b, const double y, const xparams _xparams = xdefault);
} // end of namespace alglib

// === IGAMMAF Package ===
// Depends on: GAMMAFUNC, NORMALDISTR
namespace alglib_impl {
double incompletegamma(double a, double x, ae_state *_state);
double incompletegammac(double a, double x, ae_state *_state);
double invincompletegammac(double a, double y0, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double incompletegamma(const double a, const double x, const xparams _xparams = xdefault);
double incompletegammac(const double a, const double x, const xparams _xparams = xdefault);
double invincompletegammac(const double a, const double y0, const xparams _xparams = xdefault);
} // end of namespace alglib

// === CHISQUAREDISTR Package ===
// Depends on: IGAMMAF
namespace alglib_impl {
double chisquaredistribution(double v, double x, ae_state *_state);
double chisquarecdistribution(double v, double x, ae_state *_state);
double invchisquaredistribution(double v, double y, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double chisquaredistribution(const double v, const double x, const xparams _xparams = xdefault);
double chisquarecdistribution(const double v, const double x, const xparams _xparams = xdefault);
double invchisquaredistribution(const double v, const double y, const xparams _xparams = xdefault);
} // end of namespace alglib

// === BINOMIALDISTR Package ===
// Depends on: (AlgLibInternal) NEARUNITYUNIT
// Depends on: IBETAF
namespace alglib_impl {
double binomialdistribution(ae_int_t k, ae_int_t n, double p, ae_state *_state);
double binomialcdistribution(ae_int_t k, ae_int_t n, double p, ae_state *_state);
double invbinomialdistribution(ae_int_t k, ae_int_t n, double y, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p, const xparams _xparams = xdefault);
double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p, const xparams _xparams = xdefault);
double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y, const xparams _xparams = xdefault);
} // end of namespace alglib

// === EXPINTEGRALS Package ===
namespace alglib_impl {
double exponentialintegralei(double x, ae_state *_state);
double exponentialintegralen(double x, ae_int_t n, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double exponentialintegralei(const double x, const xparams _xparams = xdefault);
double exponentialintegralen(const double x, const ae_int_t n, const xparams _xparams = xdefault);
} // end of namespace alglib

// === JACOBIANELLIPTIC Package ===
namespace alglib_impl {
void jacobianellipticfunctions(double u, double m, double *sn, double *cn, double *dn, double *ph, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph, const xparams _xparams = xdefault);
} // end of namespace alglib

// === TRIGINTEGRALS Package ===
namespace alglib_impl {
void sinecosineintegrals(double x, double *si, double *ci, ae_state *_state);
void hyperbolicsinecosineintegrals(double x, double *shi, double *chi, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void sinecosineintegrals(const double x, double &si, double &ci, const xparams _xparams = xdefault);
void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi, const xparams _xparams = xdefault);
} // end of namespace alglib

// === CHEBYSHEV Package ===
namespace alglib_impl {
double chebyshevcalculate(ae_int_t r, ae_int_t n, double x, ae_state *_state);
double chebyshevsum(RVector *c, ae_int_t r, ae_int_t n, double x, ae_state *_state);
void chebyshevcoefficients(ae_int_t n, RVector *c, ae_state *_state);
void fromchebyshev(RVector *a, ae_int_t n, RVector *b, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x, const xparams _xparams = xdefault);
double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x, const xparams _xparams = xdefault);
void chebyshevcoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = xdefault);
void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b, const xparams _xparams = xdefault);
} // end of namespace alglib

// === POISSONDISTR Package ===
// Depends on: IGAMMAF
namespace alglib_impl {
double poissondistribution(ae_int_t k, double m, ae_state *_state);
double poissoncdistribution(ae_int_t k, double m, ae_state *_state);
double invpoissondistribution(ae_int_t k, double y, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double poissondistribution(const ae_int_t k, const double m, const xparams _xparams = xdefault);
double poissoncdistribution(const ae_int_t k, const double m, const xparams _xparams = xdefault);
double invpoissondistribution(const ae_int_t k, const double y, const xparams _xparams = xdefault);
} // end of namespace alglib

// === BETAF Package ===
// Depends on: GAMMAFUNC
namespace alglib_impl {
double beta(double a, double b, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double beta(const double a, const double b, const xparams _xparams = xdefault);
} // end of namespace alglib

// === FRESNEL Package ===
namespace alglib_impl {
void fresnelintegral(double x, double *c, double *s, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void fresnelintegral(const double x, double &c, double &s, const xparams _xparams = xdefault);
} // end of namespace alglib

// === PSIF Package ===
namespace alglib_impl {
double psi(double x, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double psi(const double x, const xparams _xparams = xdefault);
} // end of namespace alglib

// === AIRYF Package ===
namespace alglib_impl {
void airy(double x, double *ai, double *aip, double *bi, double *bip, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void airy(const double x, double &ai, double &aip, double &bi, double &bip, const xparams _xparams = xdefault);
} // end of namespace alglib

// === DAWSON Package ===
namespace alglib_impl {
double dawsonintegral(double x, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double dawsonintegral(const double x, const xparams _xparams = xdefault);
} // end of namespace alglib

// === HERMITE Package ===
namespace alglib_impl {
double hermitecalculate(ae_int_t n, double x, ae_state *_state);
double hermitesum(RVector *c, ae_int_t n, double x, ae_state *_state);
void hermitecoefficients(ae_int_t n, RVector *c, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double hermitecalculate(const ae_int_t n, const double x, const xparams _xparams = xdefault);
double hermitesum(const real_1d_array &c, const ae_int_t n, const double x, const xparams _xparams = xdefault);
void hermitecoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LEGENDRE Package ===
namespace alglib_impl {
double legendrecalculate(ae_int_t n, double x, ae_state *_state);
double legendresum(RVector *c, ae_int_t n, double x, ae_state *_state);
void legendrecoefficients(ae_int_t n, RVector *c, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double legendrecalculate(const ae_int_t n, const double x, const xparams _xparams = xdefault);
double legendresum(const real_1d_array &c, const ae_int_t n, const double x, const xparams _xparams = xdefault);
void legendrecoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = xdefault);
} // end of namespace alglib

// === BESSEL Package ===
namespace alglib_impl {
double besselj0(double x, ae_state *_state);
double besselj1(double x, ae_state *_state);
double besseljn(ae_int_t n, double x, ae_state *_state);
double bessely0(double x, ae_state *_state);
double bessely1(double x, ae_state *_state);
double besselyn(ae_int_t n, double x, ae_state *_state);
double besseli0(double x, ae_state *_state);
double besseli1(double x, ae_state *_state);
double besselk0(double x, ae_state *_state);
double besselk1(double x, ae_state *_state);
double besselkn(ae_int_t nn, double x, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double besselj0(const double x, const xparams _xparams = xdefault);
double besselj1(const double x, const xparams _xparams = xdefault);
double besseljn(const ae_int_t n, const double x, const xparams _xparams = xdefault);
double bessely0(const double x, const xparams _xparams = xdefault);
double bessely1(const double x, const xparams _xparams = xdefault);
double besselyn(const ae_int_t n, const double x, const xparams _xparams = xdefault);
double besseli0(const double x, const xparams _xparams = xdefault);
double besseli1(const double x, const xparams _xparams = xdefault);
double besselk0(const double x, const xparams _xparams = xdefault);
double besselk1(const double x, const xparams _xparams = xdefault);
double besselkn(const ae_int_t nn, const double x, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LAGUERRE Package ===
namespace alglib_impl {
double laguerrecalculate(ae_int_t n, double x, ae_state *_state);
double laguerresum(RVector *c, ae_int_t n, double x, ae_state *_state);
void laguerrecoefficients(ae_int_t n, RVector *c, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double laguerrecalculate(const ae_int_t n, const double x, const xparams _xparams = xdefault);
double laguerresum(const real_1d_array &c, const ae_int_t n, const double x, const xparams _xparams = xdefault);
void laguerrecoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = xdefault);
} // end of namespace alglib

// === ELLIPTIC Package ===
namespace alglib_impl {
double ellipticintegralk(double m, ae_state *_state);
double ellipticintegralkhighprecision(double m1, ae_state *_state);
double incompleteellipticintegralk(double phi, double m, ae_state *_state);
double ellipticintegrale(double m, ae_state *_state);
double incompleteellipticintegrale(double phi, double m, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double ellipticintegralk(const double m, const xparams _xparams = xdefault);
double ellipticintegralkhighprecision(const double m1, const xparams _xparams = xdefault);
double incompleteellipticintegralk(const double phi, const double m, const xparams _xparams = xdefault);
double ellipticintegrale(const double m, const xparams _xparams = xdefault);
double incompleteellipticintegrale(const double phi, const double m, const xparams _xparams = xdefault);
} // end of namespace alglib

#endif // OnceOnly
