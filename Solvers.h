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
#ifndef OnceOnlySolvers_h
#define OnceOnlySolvers_h

#include "LinAlg.h"

// === POLYNOMIALSOLVER Package ===
// Depends on: (LinAlg) EVD, TRFAC
namespace alglib_impl {
struct polynomialsolverreport {
   double maxerr;
};
void polynomialsolverreport_init(void *_p, ae_state *_state, bool make_automatic);
void polynomialsolverreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void polynomialsolverreport_free(void *_p, bool make_automatic);

void polynomialsolve(RVector *a, ae_int_t n, CVector *x, polynomialsolverreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(polynomialsolverreport, double &maxerr;);

void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === DIRECTDENSESOLVERS Package ===
// Depends on: (AlgLibInternal) XBLAS
// Depends on: (LinAlg) SVD, RCOND
namespace alglib_impl {
struct densesolverreport {
   double r1;
   double rinf;
};
void densesolverreport_init(void *_p, ae_state *_state, bool make_automatic);
void densesolverreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void densesolverreport_free(void *_p, bool make_automatic);

struct densesolverlsreport {
   double r2;
   ae_matrix cx;
   ae_int_t n;
   ae_int_t k;
};
void densesolverlsreport_init(void *_p, ae_state *_state, bool make_automatic);
void densesolverlsreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void densesolverlsreport_free(void *_p, bool make_automatic);

void rmatrixsolvem(RMatrix *a, ae_int_t n, RMatrix *b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, RMatrix *x, ae_state *_state);
void rmatrixsolvemfast(RMatrix *a, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void rmatrixsolve(RMatrix *a, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x, ae_state *_state);
void rmatrixsolvefast(RMatrix *a, ae_int_t n, RVector *b, ae_int_t *info, ae_state *_state);
void rmatrixlusolvem(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x, ae_state *_state);
void rmatrixlusolvemfast(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void rmatrixlusolve(RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x, ae_state *_state);
void rmatrixlusolvefast(RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, ae_state *_state);
void rmatrixmixedsolvem(RMatrix *a, RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x, ae_state *_state);
void rmatrixmixedsolve(RMatrix *a, RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x, ae_state *_state);
void cmatrixsolvem(CMatrix *a, ae_int_t n, CMatrix *b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, CMatrix *x, ae_state *_state);
void cmatrixsolvemfast(CMatrix *a, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void cmatrixsolve(CMatrix *a, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x, ae_state *_state);
void cmatrixsolvefast(CMatrix *a, ae_int_t n, CVector *b, ae_int_t *info, ae_state *_state);
void cmatrixlusolvem(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x, ae_state *_state);
void cmatrixlusolvemfast(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void cmatrixlusolve(CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x, ae_state *_state);
void cmatrixlusolvefast(CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, ae_state *_state);
void cmatrixmixedsolvem(CMatrix *a, CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x, ae_state *_state);
void cmatrixmixedsolve(CMatrix *a, CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x, ae_state *_state);
void spdmatrixsolvem(RMatrix *a, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x, ae_state *_state);
void spdmatrixsolvemfast(RMatrix *a, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void spdmatrixsolve(RMatrix *a, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x, ae_state *_state);
void spdmatrixsolvefast(RMatrix *a, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, ae_state *_state);
void spdmatrixcholeskysolvem(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x, ae_state *_state);
void spdmatrixcholeskysolvemfast(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void spdmatrixcholeskysolve(RMatrix *cha, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x, ae_state *_state);
void spdmatrixcholeskysolvefast(RMatrix *cha, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, ae_state *_state);
void hpdmatrixsolvem(CMatrix *a, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x, ae_state *_state);
void hpdmatrixsolvemfast(CMatrix *a, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void hpdmatrixsolve(CMatrix *a, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x, ae_state *_state);
void hpdmatrixsolvefast(CMatrix *a, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, ae_state *_state);
void hpdmatrixcholeskysolvem(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x, ae_state *_state);
void hpdmatrixcholeskysolvemfast(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, ae_state *_state);
void hpdmatrixcholeskysolve(CMatrix *cha, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x, ae_state *_state);
void hpdmatrixcholeskysolvefast(CMatrix *cha, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, ae_state *_state);
void rmatrixsolvels(RMatrix *a, ae_int_t nrows, ae_int_t ncols, RVector *b, double threshold, ae_int_t *info, densesolverlsreport *rep, RVector *x, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(densesolverreport, double &r1; double &rinf;);
DecClass(densesolverlsreport, double &r2; real_2d_array cx; ae_int_t &n; ae_int_t &k;);

void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = xdefault);
void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = xdefault);
void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = xdefault);
void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = xdefault);
void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = xdefault);
void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = xdefault);
void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = xdefault);
void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = xdefault);
void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = xdefault);
void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = xdefault);
void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = xdefault);
void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = xdefault);
void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = xdefault);
void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = xdefault);
void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = xdefault);
void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = xdefault);
void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = xdefault);
void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = xdefault);
void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = xdefault);
void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = xdefault);
void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = xdefault);
void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = xdefault);
void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x, const xparams _xparams = xdefault);
} // end of namespace alglib

// === DIRECTSPARSESOLVERS Package ===
// Depends on: (LinAlg) TRFAC
namespace alglib_impl {
struct sparsesolverreport {
   ae_int_t terminationtype;
   ae_int_t nmv;
   ae_int_t iterationscount;
   double r2;
};
void sparsesolverreport_init(void *_p, ae_state *_state, bool make_automatic);
void sparsesolverreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sparsesolverreport_free(void *_p, bool make_automatic);

void initsparsesolverreport(sparsesolverreport *rep, ae_state *_state);
void sparsespdsolvesks(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsespdsolve(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsespdcholeskysolve(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsesolve(sparsematrix *a, RVector *b, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparselusolve(sparsematrix *a, ZVector *p, ZVector *q, RVector *b, RVector *x, sparsesolverreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(sparsesolverreport, ae_int_t &terminationtype; ae_int_t &nmv; ae_int_t &iterationscount; double &r2;);

void sparsespdsolvesks(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsespdsolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsespdcholeskysolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsesolve(const sparsematrix &a, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparselusolve(const sparsematrix &a, const integer_1d_array &p, const integer_1d_array &q, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === ITERATIVESPARSE Package ===
// Depends on: (LinAlg) FBLS
// Depends on: DIRECTSPARSESOLVERS
namespace alglib_impl {
struct sparsesolverstate {
   ae_int_t n;
   ae_vector x0;
   double epsf;
   ae_int_t maxits;
   ae_int_t algotype;
   ae_int_t gmresk;
   bool xrep;
   bool running;
   bool userterminationneeded;
   ae_vector b;
   ae_vector xf;
   ae_int_t repiterationscount;
   ae_int_t repnmv;
   ae_int_t repterminationtype;
   double repr2;
   ae_int_t requesttype;
   ae_vector x;
   ae_vector ax;
   double reply1;
   ae_vector wrkb;
   sparsematrix convbuf;
   fblsgmresstate gmressolver;
   rcommstate rstate;
};
void sparsesolverstate_init(void *_p, ae_state *_state, bool make_automatic);
void sparsesolverstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sparsesolverstate_free(void *_p, bool make_automatic);

void sparsesolversetalgogmres(sparsesolverstate *state, ae_int_t k, ae_state *_state);
void sparsesolversetstartingpoint(sparsesolverstate *state, RVector *x, ae_state *_state);
void sparsesolversetcond(sparsesolverstate *state, double epsf, ae_int_t maxits, ae_state *_state);
void sparsesolversetxrep(sparsesolverstate *state, bool needxrep, ae_state *_state);
void sparsesolvercreate(ae_int_t n, sparsesolverstate *state, ae_state *_state);
void sparsesolveroocstart(sparsesolverstate *state, RVector *b, ae_state *_state);
bool sparsesolverooccontinue(sparsesolverstate *state, ae_state *_state);
void sparsesolveroocgetrequestinfo(sparsesolverstate *state, ae_int_t *requesttype, ae_state *_state);
void sparsesolveroocgetrequestdata(sparsesolverstate *state, RVector *x, ae_state *_state);
void sparsesolveroocgetrequestdata1(sparsesolverstate *state, double *v, ae_state *_state);
void sparsesolveroocsendresult(sparsesolverstate *state, RVector *ax, ae_state *_state);
void sparsesolveroocstop(sparsesolverstate *state, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsesolverresults(sparsesolverstate *state, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsesolversolvesymmetric(sparsesolverstate *state, sparsematrix *a, bool isupper, RVector *b, ae_state *_state);
void sparsesolvesymmetricgmres(sparsematrix *a, bool isupper, RVector *b, ae_int_t k, double epsf, ae_int_t maxits, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsesolversolve(sparsesolverstate *state, sparsematrix *a, RVector *b, ae_state *_state);
void sparsesolvegmres(sparsematrix *a, RVector *b, ae_int_t k, double epsf, ae_int_t maxits, RVector *x, sparsesolverreport *rep, ae_state *_state);
void sparsesolverrequesttermination(sparsesolverstate *state, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(sparsesolverstate, );

void sparsesolversetalgogmres(const sparsesolverstate &state, const ae_int_t k, const xparams _xparams = xdefault);
void sparsesolversetstartingpoint(const sparsesolverstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void sparsesolversetcond(const sparsesolverstate &state, const double epsf, const ae_int_t maxits, const xparams _xparams = xdefault);
void sparsesolversetxrep(const sparsesolverstate &state, const bool needxrep, const xparams _xparams = xdefault);
void sparsesolvercreate(const ae_int_t n, sparsesolverstate &state, const xparams _xparams = xdefault);
void sparsesolveroocstart(const sparsesolverstate &state, const real_1d_array &b, const xparams _xparams = xdefault);
bool sparsesolverooccontinue(const sparsesolverstate &state, const xparams _xparams = xdefault);
void sparsesolveroocgetrequestinfo(const sparsesolverstate &state, ae_int_t &requesttype, const xparams _xparams = xdefault);
void sparsesolveroocgetrequestdata(const sparsesolverstate &state, real_1d_array &x, const xparams _xparams = xdefault);
void sparsesolveroocgetrequestdata1(const sparsesolverstate &state, double &v, const xparams _xparams = xdefault);
void sparsesolveroocsendresult(const sparsesolverstate &state, const real_1d_array &ax, const xparams _xparams = xdefault);
void sparsesolveroocstop(const sparsesolverstate &state, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsesolverresults(const sparsesolverstate &state, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsesolversolvesymmetric(const sparsesolverstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b, const xparams _xparams = xdefault);
void sparsesolvesymmetricgmres(const sparsematrix &a, const bool isupper, const real_1d_array &b, const ae_int_t k, const double epsf, const ae_int_t maxits, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsesolversolve(const sparsesolverstate &state, const sparsematrix &a, const real_1d_array &b, const xparams _xparams = xdefault);
void sparsesolvegmres(const sparsematrix &a, const real_1d_array &b, const ae_int_t k, const double epsf, const ae_int_t maxits, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = xdefault);
void sparsesolverrequesttermination(const sparsesolverstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LINCG Package ===
// Depends on: (LinAlg) MATGEN, SPARSE
namespace alglib_impl {
struct lincgstate {
   ae_vector rx;
   ae_vector b;
   ae_int_t n;
   ae_int_t prectype;
   ae_vector cx;
   ae_vector cr;
   ae_vector cz;
   ae_vector p;
   ae_vector r;
   ae_vector z;
   double alpha;
   double beta;
   double r2;
   double meritfunction;
   ae_vector x;
   ae_vector mv;
   ae_vector pv;
   double vmv;
   ae_vector startx;
   double epsf;
   ae_int_t maxits;
   ae_int_t itsbeforerestart;
   ae_int_t itsbeforerupdate;
   bool xrep;
   bool xupdated;
   bool needmv;
   bool needmtv;
   bool needmv2;
   bool needvmv;
   bool needprec;
   ae_int_t repiterationscount;
   ae_int_t repnmv;
   ae_int_t repterminationtype;
   bool running;
   ae_vector tmpd;
   rcommstate rstate;
};
void lincgstate_init(void *_p, ae_state *_state, bool make_automatic);
void lincgstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void lincgstate_free(void *_p, bool make_automatic);

struct lincgreport {
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
   double r2;
};
void lincgreport_init(void *_p, ae_state *_state, bool make_automatic);
void lincgreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void lincgreport_free(void *_p, bool make_automatic);

void lincgcreate(ae_int_t n, lincgstate *state, ae_state *_state);
void lincgsetstartingpoint(lincgstate *state, RVector *x, ae_state *_state);
void lincgsetb(lincgstate *state, RVector *b, ae_state *_state);
void lincgsetprecunit(lincgstate *state, ae_state *_state);
void lincgsetprecdiag(lincgstate *state, ae_state *_state);
void lincgsetcond(lincgstate *state, double epsf, ae_int_t maxits, ae_state *_state);
bool lincgiteration(lincgstate *state, ae_state *_state);
void lincgrestart(lincgstate *state, ae_state *_state);
void lincgsolvesparse(lincgstate *state, sparsematrix *a, bool isupper, RVector *b, ae_state *_state);
void lincgresults(lincgstate *state, RVector *x, lincgreport *rep, ae_state *_state);
void lincgsetrestartfreq(lincgstate *state, ae_int_t srf, ae_state *_state);
void lincgsetrupdatefreq(lincgstate *state, ae_int_t freq, ae_state *_state);
void lincgsetxrep(lincgstate *state, bool needxrep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(lincgstate, );
DecClass(lincgreport, ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype; double &r2;);

void lincgcreate(const ae_int_t n, lincgstate &state, const xparams _xparams = xdefault);
void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void lincgsetprecunit(const lincgstate &state, const xparams _xparams = xdefault);
void lincgsetprecdiag(const lincgstate &state, const xparams _xparams = xdefault);
void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits, const xparams _xparams = xdefault);
void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b, const xparams _xparams = xdefault);
void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep, const xparams _xparams = xdefault);
void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf, const xparams _xparams = xdefault);
void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq, const xparams _xparams = xdefault);
void lincgsetxrep(const lincgstate &state, const bool needxrep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LINLSQR Package ===
// Depends on: (LinAlg) SVD, NORMESTIMATOR
namespace alglib_impl {
struct linlsqrstate {
   normestimatorstate nes;
   ae_vector rx;
   ae_vector b;
   ae_int_t n;
   ae_int_t m;
   ae_int_t prectype;
   ae_vector ui;
   ae_vector uip1;
   ae_vector vi;
   ae_vector vip1;
   ae_vector omegai;
   ae_vector omegaip1;
   double alphai;
   double alphaip1;
   double betai;
   double betaip1;
   double phibari;
   double phibarip1;
   double phii;
   double rhobari;
   double rhobarip1;
   double rhoi;
   double ci;
   double si;
   double theta;
   double lambdai;
   ae_vector d;
   double anorm;
   double bnorm2;
   double dnorm;
   double r2;
   ae_vector x;
   ae_vector mv;
   ae_vector mtv;
   double epsa;
   double epsb;
   double epsc;
   ae_int_t maxits;
   bool xrep;
   bool xupdated;
   bool needmv;
   bool needmtv;
   bool needmv2;
   bool needvmv;
   bool needprec;
   ae_int_t repiterationscount;
   ae_int_t repnmv;
   ae_int_t repterminationtype;
   bool running;
   bool userterminationneeded;
   ae_vector tmpd;
   ae_vector tmpx;
   rcommstate rstate;
};
void linlsqrstate_init(void *_p, ae_state *_state, bool make_automatic);
void linlsqrstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void linlsqrstate_free(void *_p, bool make_automatic);

struct linlsqrreport {
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
};
void linlsqrreport_init(void *_p, ae_state *_state, bool make_automatic);
void linlsqrreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void linlsqrreport_free(void *_p, bool make_automatic);

void linlsqrcreatebuf(ae_int_t m, ae_int_t n, linlsqrstate *state, ae_state *_state);
void linlsqrcreate(ae_int_t m, ae_int_t n, linlsqrstate *state, ae_state *_state);
void linlsqrsetb(linlsqrstate *state, RVector *b, ae_state *_state);
void linlsqrsetprecunit(linlsqrstate *state, ae_state *_state);
void linlsqrsetprecdiag(linlsqrstate *state, ae_state *_state);
void linlsqrsetlambdai(linlsqrstate *state, double lambdai, ae_state *_state);
bool linlsqriteration(linlsqrstate *state, ae_state *_state);
void linlsqrrestart(linlsqrstate *state, ae_state *_state);
void linlsqrsolvesparse(linlsqrstate *state, sparsematrix *a, RVector *b, ae_state *_state);
void linlsqrsetcond(linlsqrstate *state, double epsa, double epsb, ae_int_t maxits, ae_state *_state);
void linlsqrresults(linlsqrstate *state, RVector *x, linlsqrreport *rep, ae_state *_state);
void linlsqrsetxrep(linlsqrstate *state, bool needxrep, ae_state *_state);
ae_int_t linlsqrpeekiterationscount(linlsqrstate *s, ae_state *_state);
void linlsqrrequesttermination(linlsqrstate *state, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(linlsqrstate, );
DecClass(linlsqrreport, ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype;);

void linlsqrcreatebuf(const ae_int_t m, const ae_int_t n, const linlsqrstate &state, const xparams _xparams = xdefault);
void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state, const xparams _xparams = xdefault);
void linlsqrsetprecunit(const linlsqrstate &state, const xparams _xparams = xdefault);
void linlsqrsetprecdiag(const linlsqrstate &state, const xparams _xparams = xdefault);
void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai, const xparams _xparams = xdefault);
void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b, const xparams _xparams = xdefault);
void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits, const xparams _xparams = xdefault);
void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep, const xparams _xparams = xdefault);
void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep, const xparams _xparams = xdefault);
ae_int_t linlsqrpeekiterationscount(const linlsqrstate &s, const xparams _xparams = xdefault);
void linlsqrrequesttermination(const linlsqrstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === NLEQ Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: (LinAlg) FBLS
namespace alglib_impl {
struct nleqstate {
   ae_int_t n;
   ae_int_t m;
   double epsf;
   ae_int_t maxits;
   bool xrep;
   double stpmax;
   ae_vector x;
   double f;
   ae_vector fi;
   ae_matrix j;
   bool needf;
   bool needfij;
   bool xupdated;
   rcommstate rstate;
   ae_int_t repiterationscount;
   ae_int_t repnfunc;
   ae_int_t repnjac;
   ae_int_t repterminationtype;
   ae_vector xbase;
   double fbase;
   double fprev;
   ae_vector candstep;
   ae_vector rightpart;
   ae_vector cgbuf;
};
void nleqstate_init(void *_p, ae_state *_state, bool make_automatic);
void nleqstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void nleqstate_free(void *_p, bool make_automatic);

struct nleqreport {
   ae_int_t iterationscount;
   ae_int_t nfunc;
   ae_int_t njac;
   ae_int_t terminationtype;
};
void nleqreport_init(void *_p, ae_state *_state, bool make_automatic);
void nleqreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void nleqreport_free(void *_p, bool make_automatic);

void nleqsetcond(nleqstate *state, double epsf, ae_int_t maxits, ae_state *_state);
void nleqsetxrep(nleqstate *state, bool needxrep, ae_state *_state);
void nleqsetstpmax(nleqstate *state, double stpmax, ae_state *_state);
void nleqrestartfrom(nleqstate *state, RVector *x, ae_state *_state);
void nleqcreatelm(ae_int_t n, ae_int_t m, RVector *x, nleqstate *state, ae_state *_state);
bool nleqiteration(nleqstate *state, ae_state *_state);
void nleqresultsbuf(nleqstate *state, RVector *x, nleqreport *rep, ae_state *_state);
void nleqresults(nleqstate *state, RVector *x, nleqreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(nleqstate, bool &needf; bool &needfij; bool &xupdated; double &f; real_1d_array fi; real_2d_array j; real_1d_array x;);
DecClass(nleqreport, ae_int_t &iterationscount; ae_int_t &nfunc; ae_int_t &njac; ae_int_t &terminationtype;);

void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits, const xparams _xparams = xdefault);
void nleqsetxrep(const nleqstate &state, const bool needxrep, const xparams _xparams = xdefault);
void nleqsetstpmax(const nleqstate &state, const double stpmax, const xparams _xparams = xdefault);
void nleqrestartfrom(const nleqstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state, const xparams _xparams = xdefault);
void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state, const xparams _xparams = xdefault);
bool nleqiteration(const nleqstate &state, const xparams _xparams = xdefault);
void nleqsolve(nleqstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep, const xparams _xparams = xdefault);
void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

#endif // OnceOnly
