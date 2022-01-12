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
double gammafunction(double x);
double lngamma(double x, double *sgngam);
} // end of namespace alglib_impl

namespace alglib {
double gammafunction(const double x);
double lngamma(const double x, double &sgngam);
} // end of namespace alglib

// === NORMALDISTR Package ===
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
double errorfunction(double x);
double errorfunctionc(double x);
double normaldistribution(double x);
double normalpdf(double x);
double normalcdf(double x);
double invnormalcdf(double y0);
double invnormaldistribution(double y0);
double inverf(double e);
double bivariatenormalpdf(double x, double y, double rho);
double bivariatenormalcdf(double x, double y, double rho);
} // end of namespace alglib_impl

namespace alglib {
double errorfunction(const double x);
double errorfunctionc(const double x);
double normaldistribution(const double x);
double normalpdf(const double x);
double normalcdf(const double x);
double invnormalcdf(const double y0);
double invnormaldistribution(const double y0);
double inverf(const double e);
double bivariatenormalpdf(const double x, const double y, const double rho);
double bivariatenormalcdf(const double x, const double y, const double rho);
} // end of namespace alglib

// === IBETAF Package ===
// Depends on: GAMMAFUNC, NORMALDISTR
namespace alglib_impl {
double incompletebeta(double a, double b, double x);
double invincompletebeta(double a, double b, double y);
} // end of namespace alglib_impl

namespace alglib {
double incompletebeta(const double a, const double b, const double x);
double invincompletebeta(const double a, const double b, const double y);
} // end of namespace alglib

// === STUDENTTDISTR Package ===
// Depends on: IBETAF
namespace alglib_impl {
double studenttdistribution(ae_int_t k, double t);
double invstudenttdistribution(ae_int_t k, double p);
} // end of namespace alglib_impl

namespace alglib {
double studenttdistribution(const ae_int_t k, const double t);
double invstudenttdistribution(const ae_int_t k, const double p);
} // end of namespace alglib

// === FDISTR Package ===
// Depends on: IBETAF
namespace alglib_impl {
double fdistribution(ae_int_t a, ae_int_t b, double x);
double fcdistribution(ae_int_t a, ae_int_t b, double x);
double invfdistribution(ae_int_t a, ae_int_t b, double y);
} // end of namespace alglib_impl

namespace alglib {
double fdistribution(const ae_int_t a, const ae_int_t b, const double x);
double fcdistribution(const ae_int_t a, const ae_int_t b, const double x);
double invfdistribution(const ae_int_t a, const ae_int_t b, const double y);
} // end of namespace alglib

// === IGAMMAF Package ===
// Depends on: GAMMAFUNC, NORMALDISTR
namespace alglib_impl {
double incompletegamma(double a, double x);
double incompletegammac(double a, double x);
double invincompletegammac(double a, double y0);
} // end of namespace alglib_impl

namespace alglib {
double incompletegamma(const double a, const double x);
double incompletegammac(const double a, const double x);
double invincompletegammac(const double a, const double y0);
} // end of namespace alglib

// === CHISQUAREDISTR Package ===
// Depends on: IGAMMAF
namespace alglib_impl {
double chisquaredistribution(double v, double x);
double chisquarecdistribution(double v, double x);
double invchisquaredistribution(double v, double y);
} // end of namespace alglib_impl

namespace alglib {
double chisquaredistribution(const double v, const double x);
double chisquarecdistribution(const double v, const double x);
double invchisquaredistribution(const double v, const double y);
} // end of namespace alglib

// === BINOMIALDISTR Package ===
// Depends on: (AlgLibInternal) NEARUNITYUNIT
// Depends on: IBETAF
namespace alglib_impl {
double binomialdistribution(ae_int_t k, ae_int_t n, double p);
double binomialcdistribution(ae_int_t k, ae_int_t n, double p);
double invbinomialdistribution(ae_int_t k, ae_int_t n, double y);
} // end of namespace alglib_impl

namespace alglib {
double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p);
double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p);
double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y);
} // end of namespace alglib

// === EXPINTEGRALS Package ===
namespace alglib_impl {
double exponentialintegralei(double x);
double exponentialintegralen(double x, ae_int_t n);
} // end of namespace alglib_impl

namespace alglib {
double exponentialintegralei(const double x);
double exponentialintegralen(const double x, const ae_int_t n);
} // end of namespace alglib

// === JACOBIANELLIPTIC Package ===
namespace alglib_impl {
void jacobianellipticfunctions(double u, double m, double *sn, double *cn, double *dn, double *ph);
} // end of namespace alglib_impl

namespace alglib {
void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph);
} // end of namespace alglib

// === TRIGINTEGRALS Package ===
namespace alglib_impl {
void sinecosineintegrals(double x, double *si, double *ci);
void hyperbolicsinecosineintegrals(double x, double *shi, double *chi);
} // end of namespace alglib_impl

namespace alglib {
void sinecosineintegrals(const double x, double &si, double &ci);
void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi);
} // end of namespace alglib

// === CHEBYSHEV Package ===
namespace alglib_impl {
double chebyshevcalculate(ae_int_t r, ae_int_t n, double x);
double chebyshevsum(RVector *c, ae_int_t r, ae_int_t n, double x);
void chebyshevcoefficients(ae_int_t n, RVector *c);
void fromchebyshev(RVector *a, ae_int_t n, RVector *b);
} // end of namespace alglib_impl

namespace alglib {
double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x);
double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x);
void chebyshevcoefficients(const ae_int_t n, real_1d_array &c);
void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b);
} // end of namespace alglib

// === POISSONDISTR Package ===
// Depends on: IGAMMAF
namespace alglib_impl {
double poissondistribution(ae_int_t k, double m);
double poissoncdistribution(ae_int_t k, double m);
double invpoissondistribution(ae_int_t k, double y);
} // end of namespace alglib_impl

namespace alglib {
double poissondistribution(const ae_int_t k, const double m);
double poissoncdistribution(const ae_int_t k, const double m);
double invpoissondistribution(const ae_int_t k, const double y);
} // end of namespace alglib

// === BETAF Package ===
// Depends on: GAMMAFUNC
namespace alglib_impl {
double beta(double a, double b);
} // end of namespace alglib_impl

namespace alglib {
double beta(const double a, const double b);
} // end of namespace alglib

// === FRESNEL Package ===
namespace alglib_impl {
void fresnelintegral(double x, double *c, double *s);
} // end of namespace alglib_impl

namespace alglib {
void fresnelintegral(const double x, double &c, double &s);
} // end of namespace alglib

// === PSIF Package ===
namespace alglib_impl {
double psi(double x);
} // end of namespace alglib_impl

namespace alglib {
double psi(const double x);
} // end of namespace alglib

// === AIRYF Package ===
namespace alglib_impl {
void airy(double x, double *ai, double *aip, double *bi, double *bip);
} // end of namespace alglib_impl

namespace alglib {
void airy(const double x, double &ai, double &aip, double &bi, double &bip);
} // end of namespace alglib

// === DAWSON Package ===
namespace alglib_impl {
double dawsonintegral(double x);
} // end of namespace alglib_impl

namespace alglib {
double dawsonintegral(const double x);
} // end of namespace alglib

// === HERMITE Package ===
namespace alglib_impl {
double hermitecalculate(ae_int_t n, double x);
double hermitesum(RVector *c, ae_int_t n, double x);
void hermitecoefficients(ae_int_t n, RVector *c);
} // end of namespace alglib_impl

namespace alglib {
double hermitecalculate(const ae_int_t n, const double x);
double hermitesum(const real_1d_array &c, const ae_int_t n, const double x);
void hermitecoefficients(const ae_int_t n, real_1d_array &c);
} // end of namespace alglib

// === LEGENDRE Package ===
namespace alglib_impl {
double legendrecalculate(ae_int_t n, double x);
double legendresum(RVector *c, ae_int_t n, double x);
void legendrecoefficients(ae_int_t n, RVector *c);
} // end of namespace alglib_impl

namespace alglib {
double legendrecalculate(const ae_int_t n, const double x);
double legendresum(const real_1d_array &c, const ae_int_t n, const double x);
void legendrecoefficients(const ae_int_t n, real_1d_array &c);
} // end of namespace alglib

// === BESSEL Package ===
namespace alglib_impl {
double besselj0(double x);
double besselj1(double x);
double besseljn(ae_int_t n, double x);
double bessely0(double x);
double bessely1(double x);
double besselyn(ae_int_t n, double x);
double besseli0(double x);
double besseli1(double x);
double besselk0(double x);
double besselk1(double x);
double besselkn(ae_int_t nn, double x);
} // end of namespace alglib_impl

namespace alglib {
double besselj0(const double x);
double besselj1(const double x);
double besseljn(const ae_int_t n, const double x);
double bessely0(const double x);
double bessely1(const double x);
double besselyn(const ae_int_t n, const double x);
double besseli0(const double x);
double besseli1(const double x);
double besselk0(const double x);
double besselk1(const double x);
double besselkn(const ae_int_t nn, const double x);
} // end of namespace alglib

// === LAGUERRE Package ===
namespace alglib_impl {
double laguerrecalculate(ae_int_t n, double x);
double laguerresum(RVector *c, ae_int_t n, double x);
void laguerrecoefficients(ae_int_t n, RVector *c);
} // end of namespace alglib_impl

namespace alglib {
double laguerrecalculate(const ae_int_t n, const double x);
double laguerresum(const real_1d_array &c, const ae_int_t n, const double x);
void laguerrecoefficients(const ae_int_t n, real_1d_array &c);
} // end of namespace alglib

// === ELLIPTIC Package ===
namespace alglib_impl {
double ellipticintegralkhighprecision(double m1);
double ellipticintegralk(double m);
double incompleteellipticintegralk(double phi, double m);
double ellipticintegrale(double m);
double incompleteellipticintegrale(double phi, double m);
} // end of namespace alglib_impl

namespace alglib {
double ellipticintegralkhighprecision(const double m1);
double ellipticintegralk(const double m);
double incompleteellipticintegralk(const double phi, const double m);
double ellipticintegrale(const double m);
double incompleteellipticintegrale(const double phi, const double m);
} // end of namespace alglib

#endif // OnceOnly
