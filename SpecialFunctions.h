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

// Declarations for the computational core: datatypes.
namespace alglib_impl {
// === GAMMAFUNC Package ===

// === NORMALDISTR Package ===

// === IBETAF Package ===

// === STUDENTTDISTR Package ===

// === FDISTR Package ===

// === IGAMMAF Package ===

// === CHISQUAREDISTR Package ===

// === BINOMIALDISTR Package ===

// === EXPINTEGRALS Package ===

// === JACOBIANELLIPTIC Package ===

// === TRIGINTEGRALS Package ===

// === CHEBYSHEV Package ===

// === POISSONDISTR Package ===

// === BETAF Package ===

// === FRESNEL Package ===

// === PSIF Package ===

// === AIRYF Package ===

// === DAWSON Package ===

// === HERMITE Package ===

// === LEGENDRE Package ===

// === BESSEL Package ===

// === LAGUERRE Package ===

// === ELLIPTIC Package ===
} // end of namespace alglib_impl

// Declarations for the C++ interface.
namespace alglib {
// === GAMMAFUNC Package ===

// === NORMALDISTR Package ===

// === IBETAF Package ===

// === STUDENTTDISTR Package ===

// === FDISTR Package ===

// === IGAMMAF Package ===

// === CHISQUAREDISTR Package ===

// === BINOMIALDISTR Package ===

// === EXPINTEGRALS Package ===

// === JACOBIANELLIPTIC Package ===

// === TRIGINTEGRALS Package ===

// === CHEBYSHEV Package ===

// === POISSONDISTR Package ===

// === BETAF Package ===

// === FRESNEL Package ===

// === PSIF Package ===

// === AIRYF Package ===

// === DAWSON Package ===

// === HERMITE Package ===

// === LEGENDRE Package ===

// === BESSEL Package ===

// === LAGUERRE Package ===

// === ELLIPTIC Package ===

// === GAMMAFUNC Package ===
double gammafunction(const double x, const xparams _xparams = alglib::xdefault);
double lngamma(const double x, double &sgngam, const xparams _xparams = alglib::xdefault);

// === NORMALDISTR Package ===
double errorfunction(const double x, const xparams _xparams = alglib::xdefault);
double errorfunctionc(const double x, const xparams _xparams = alglib::xdefault);
double normaldistribution(const double x, const xparams _xparams = alglib::xdefault);
double normalpdf(const double x, const xparams _xparams = alglib::xdefault);
double normalcdf(const double x, const xparams _xparams = alglib::xdefault);
double inverf(const double e, const xparams _xparams = alglib::xdefault);
double invnormaldistribution(const double y0, const xparams _xparams = alglib::xdefault);
double invnormalcdf(const double y0, const xparams _xparams = alglib::xdefault);
double bivariatenormalpdf(const double x, const double y, const double rho, const xparams _xparams = alglib::xdefault);
double bivariatenormalcdf(const double x, const double y, const double rho, const xparams _xparams = alglib::xdefault);

// === IBETAF Package ===
double incompletebeta(const double a, const double b, const double x, const xparams _xparams = alglib::xdefault);
double invincompletebeta(const double a, const double b, const double y, const xparams _xparams = alglib::xdefault);

// === STUDENTTDISTR Package ===
double studenttdistribution(const ae_int_t k, const double t, const xparams _xparams = alglib::xdefault);
double invstudenttdistribution(const ae_int_t k, const double p, const xparams _xparams = alglib::xdefault);

// === FDISTR Package ===
double fdistribution(const ae_int_t a, const ae_int_t b, const double x, const xparams _xparams = alglib::xdefault);
double fcdistribution(const ae_int_t a, const ae_int_t b, const double x, const xparams _xparams = alglib::xdefault);
double invfdistribution(const ae_int_t a, const ae_int_t b, const double y, const xparams _xparams = alglib::xdefault);

// === IGAMMAF Package ===
double incompletegamma(const double a, const double x, const xparams _xparams = alglib::xdefault);
double incompletegammac(const double a, const double x, const xparams _xparams = alglib::xdefault);
double invincompletegammac(const double a, const double y0, const xparams _xparams = alglib::xdefault);

// === CHISQUAREDISTR Package ===
double chisquaredistribution(const double v, const double x, const xparams _xparams = alglib::xdefault);
double chisquarecdistribution(const double v, const double x, const xparams _xparams = alglib::xdefault);
double invchisquaredistribution(const double v, const double y, const xparams _xparams = alglib::xdefault);

// === BINOMIALDISTR Package ===
double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p, const xparams _xparams = alglib::xdefault);
double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p, const xparams _xparams = alglib::xdefault);
double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y, const xparams _xparams = alglib::xdefault);

// === EXPINTEGRALS Package ===
double exponentialintegralei(const double x, const xparams _xparams = alglib::xdefault);
double exponentialintegralen(const double x, const ae_int_t n, const xparams _xparams = alglib::xdefault);

// === JACOBIANELLIPTIC Package ===
void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph, const xparams _xparams = alglib::xdefault);

// === TRIGINTEGRALS Package ===
void sinecosineintegrals(const double x, double &si, double &ci, const xparams _xparams = alglib::xdefault);
void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi, const xparams _xparams = alglib::xdefault);

// === CHEBYSHEV Package ===
double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
void chebyshevcoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = alglib::xdefault);
void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b, const xparams _xparams = alglib::xdefault);

// === POISSONDISTR Package ===
double poissondistribution(const ae_int_t k, const double m, const xparams _xparams = alglib::xdefault);
double poissoncdistribution(const ae_int_t k, const double m, const xparams _xparams = alglib::xdefault);
double invpoissondistribution(const ae_int_t k, const double y, const xparams _xparams = alglib::xdefault);

// === BETAF Package ===
double beta(const double a, const double b, const xparams _xparams = alglib::xdefault);

// === FRESNEL Package ===
void fresnelintegral(const double x, double &c, double &s, const xparams _xparams = alglib::xdefault);

// === PSIF Package ===
double psi(const double x, const xparams _xparams = alglib::xdefault);

// === AIRYF Package ===
void airy(const double x, double &ai, double &aip, double &bi, double &bip, const xparams _xparams = alglib::xdefault);

// === DAWSON Package ===
double dawsonintegral(const double x, const xparams _xparams = alglib::xdefault);

// === HERMITE Package ===
double hermitecalculate(const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
double hermitesum(const real_1d_array &c, const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
void hermitecoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = alglib::xdefault);

// === LEGENDRE Package ===
double legendrecalculate(const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
double legendresum(const real_1d_array &c, const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
void legendrecoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = alglib::xdefault);

// === BESSEL Package ===
double besselj0(const double x, const xparams _xparams = alglib::xdefault);
double besselj1(const double x, const xparams _xparams = alglib::xdefault);
double besseljn(const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
double bessely0(const double x, const xparams _xparams = alglib::xdefault);
double bessely1(const double x, const xparams _xparams = alglib::xdefault);
double besselyn(const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
double besseli0(const double x, const xparams _xparams = alglib::xdefault);
double besseli1(const double x, const xparams _xparams = alglib::xdefault);
double besselk0(const double x, const xparams _xparams = alglib::xdefault);
double besselk1(const double x, const xparams _xparams = alglib::xdefault);
double besselkn(const ae_int_t nn, const double x, const xparams _xparams = alglib::xdefault);

// === LAGUERRE Package ===
double laguerrecalculate(const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
double laguerresum(const real_1d_array &c, const ae_int_t n, const double x, const xparams _xparams = alglib::xdefault);
void laguerrecoefficients(const ae_int_t n, real_1d_array &c, const xparams _xparams = alglib::xdefault);

// === ELLIPTIC Package ===
double ellipticintegralk(const double m, const xparams _xparams = alglib::xdefault);
double ellipticintegralkhighprecision(const double m1, const xparams _xparams = alglib::xdefault);
double incompleteellipticintegralk(const double phi, const double m, const xparams _xparams = alglib::xdefault);
double ellipticintegrale(const double m, const xparams _xparams = alglib::xdefault);
double incompleteellipticintegrale(const double phi, const double m, const xparams _xparams = alglib::xdefault);
} // end of namespace alglib

// Declarations for the computational core: functions.
namespace alglib_impl {
// === GAMMAFUNC Package ===
double gammafunction(double x, ae_state *_state);
double lngamma(double x, double *sgngam, ae_state *_state);

// === NORMALDISTR Package ===
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

// === IBETAF Package ===
double incompletebeta(double a, double b, double x, ae_state *_state);
double invincompletebeta(double a, double b, double y, ae_state *_state);

// === STUDENTTDISTR Package ===
double studenttdistribution(ae_int_t k, double t, ae_state *_state);
double invstudenttdistribution(ae_int_t k, double p, ae_state *_state);

// === FDISTR Package ===
double fdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
double fcdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
double invfdistribution(ae_int_t a, ae_int_t b, double y, ae_state *_state);

// === IGAMMAF Package ===
double incompletegamma(double a, double x, ae_state *_state);
double incompletegammac(double a, double x, ae_state *_state);
double invincompletegammac(double a, double y0, ae_state *_state);

// === CHISQUAREDISTR Package ===
double chisquaredistribution(double v, double x, ae_state *_state);
double chisquarecdistribution(double v, double x, ae_state *_state);
double invchisquaredistribution(double v, double y, ae_state *_state);

// === BINOMIALDISTR Package ===
double binomialdistribution(ae_int_t k, ae_int_t n, double p, ae_state *_state);
double binomialcdistribution(ae_int_t k, ae_int_t n, double p, ae_state *_state);
double invbinomialdistribution(ae_int_t k, ae_int_t n, double y, ae_state *_state);

// === EXPINTEGRALS Package ===
double exponentialintegralei(double x, ae_state *_state);
double exponentialintegralen(double x, ae_int_t n, ae_state *_state);

// === JACOBIANELLIPTIC Package ===
void jacobianellipticfunctions(double u, double m, double *sn, double *cn, double *dn, double *ph, ae_state *_state);

// === TRIGINTEGRALS Package ===
void sinecosineintegrals(double x, double *si, double *ci, ae_state *_state);
void hyperbolicsinecosineintegrals(double x, double *shi, double *chi, ae_state *_state);

// === CHEBYSHEV Package ===
double chebyshevcalculate(ae_int_t r, ae_int_t n, double x, ae_state *_state);
double chebyshevsum(RVector *c, ae_int_t r, ae_int_t n, double x, ae_state *_state);
void chebyshevcoefficients(ae_int_t n, RVector *c, ae_state *_state);
void fromchebyshev(RVector *a, ae_int_t n, RVector *b, ae_state *_state);

// === POISSONDISTR Package ===
double poissondistribution(ae_int_t k, double m, ae_state *_state);
double poissoncdistribution(ae_int_t k, double m, ae_state *_state);
double invpoissondistribution(ae_int_t k, double y, ae_state *_state);

// === BETAF Package ===
double beta(double a, double b, ae_state *_state);

// === FRESNEL Package ===
void fresnelintegral(double x, double *c, double *s, ae_state *_state);

// === PSIF Package ===
double psi(double x, ae_state *_state);

// === AIRYF Package ===
void airy(double x, double *ai, double *aip, double *bi, double *bip, ae_state *_state);

// === DAWSON Package ===
double dawsonintegral(double x, ae_state *_state);

// === HERMITE Package ===
double hermitecalculate(ae_int_t n, double x, ae_state *_state);
double hermitesum(RVector *c, ae_int_t n, double x, ae_state *_state);
void hermitecoefficients(ae_int_t n, RVector *c, ae_state *_state);

// === LEGENDRE Package ===
double legendrecalculate(ae_int_t n, double x, ae_state *_state);
double legendresum(RVector *c, ae_int_t n, double x, ae_state *_state);
void legendrecoefficients(ae_int_t n, RVector *c, ae_state *_state);

// === BESSEL Package ===
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

// === LAGUERRE Package ===
double laguerrecalculate(ae_int_t n, double x, ae_state *_state);
double laguerresum(RVector *c, ae_int_t n, double x, ae_state *_state);
void laguerrecoefficients(ae_int_t n, RVector *c, ae_state *_state);

// === ELLIPTIC Package ===
double ellipticintegralk(double m, ae_state *_state);
double ellipticintegralkhighprecision(double m1, ae_state *_state);
double incompleteellipticintegralk(double phi, double m, ae_state *_state);
double ellipticintegrale(double m, ae_state *_state);
double incompleteellipticintegrale(double phi, double m, ae_state *_state);
} // end of namespace alglib_impl

#endif // OnceOnly
