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
void polynomialsolverreport_init(void *_p, bool make_automatic);
void polynomialsolverreport_copy(void *_dst, void *_src, bool make_automatic);
void polynomialsolverreport_free(void *_p, bool make_automatic);

void polynomialsolve(RVector *a, ae_int_t n, CVector *x, polynomialsolverreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(polynomialsolverreport, double &maxerr;);

void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep);
} // end of namespace alglib

// === DIRECTDENSESOLVERS Package ===
// Depends on: (AlgLibInternal) XBLAS
// Depends on: (LinAlg) SVD, RCOND
namespace alglib_impl {
struct densesolverreport {
   double r1;
   double rinf;
};
void densesolverreport_init(void *_p, bool make_automatic);
void densesolverreport_copy(void *_dst, void *_src, bool make_automatic);
void densesolverreport_free(void *_p, bool make_automatic);

struct densesolverlsreport {
   double r2;
   ae_matrix cx;
   ae_int_t n;
   ae_int_t k;
};
void densesolverlsreport_init(void *_p, bool make_automatic);
void densesolverlsreport_copy(void *_dst, void *_src, bool make_automatic);
void densesolverlsreport_free(void *_p, bool make_automatic);

void rmatrixsolvem(RMatrix *a, ae_int_t n, RMatrix *b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, RMatrix *x);
void rmatrixsolvemfast(RMatrix *a, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info);
void rmatrixsolve(RMatrix *a, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x);
void rmatrixsolvefast(RMatrix *a, ae_int_t n, RVector *b, ae_int_t *info);
void rmatrixlusolvem(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x);
void rmatrixlusolvemfast(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info);
void rmatrixlusolve(RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x);
void rmatrixlusolvefast(RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info);
void rmatrixmixedsolvem(RMatrix *a, RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x);
void rmatrixmixedsolve(RMatrix *a, RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x);
void cmatrixsolvem(CMatrix *a, ae_int_t n, CMatrix *b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, CMatrix *x);
void cmatrixsolvemfast(CMatrix *a, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info);
void cmatrixsolve(CMatrix *a, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x);
void cmatrixsolvefast(CMatrix *a, ae_int_t n, CVector *b, ae_int_t *info);
void cmatrixlusolvem(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x);
void cmatrixlusolvemfast(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info);
void cmatrixlusolve(CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x);
void cmatrixlusolvefast(CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info);
void cmatrixmixedsolvem(CMatrix *a, CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x);
void cmatrixmixedsolve(CMatrix *a, CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x);
void spdmatrixsolvem(RMatrix *a, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x);
void spdmatrixsolvemfast(RMatrix *a, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info);
void spdmatrixsolve(RMatrix *a, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x);
void spdmatrixsolvefast(RMatrix *a, ae_int_t n, bool isupper, RVector *b, ae_int_t *info);
void spdmatrixcholeskysolvem(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x);
void spdmatrixcholeskysolvemfast(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info);
void spdmatrixcholeskysolve(RMatrix *cha, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x);
void spdmatrixcholeskysolvefast(RMatrix *cha, ae_int_t n, bool isupper, RVector *b, ae_int_t *info);
void hpdmatrixsolvem(CMatrix *a, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x);
void hpdmatrixsolvemfast(CMatrix *a, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info);
void hpdmatrixsolve(CMatrix *a, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x);
void hpdmatrixsolvefast(CMatrix *a, ae_int_t n, bool isupper, CVector *b, ae_int_t *info);
void hpdmatrixcholeskysolvem(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x);
void hpdmatrixcholeskysolvemfast(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info);
void hpdmatrixcholeskysolve(CMatrix *cha, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x);
void hpdmatrixcholeskysolvefast(CMatrix *cha, ae_int_t n, bool isupper, CVector *b, ae_int_t *info);
void rmatrixsolvels(RMatrix *a, ae_int_t nrows, ae_int_t ncols, RVector *b, double threshold, ae_int_t *info, densesolverlsreport *rep, RVector *x);
} // end of namespace alglib_impl

namespace alglib {
DecClass(densesolverreport, double &r1; double &rinf;);
DecClass(densesolverlsreport, double &r2; real_2d_array cx; ae_int_t &n; ae_int_t &k;);

void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);
} // end of namespace alglib

// === DIRECTSPARSESOLVERS Package ===
// Depends on: (LinAlg) TRFAC
namespace alglib_impl {
struct sparsesolverreport {
   ae_int_t terminationtype;
};
void sparsesolverreport_init(void *_p, bool make_automatic);
void sparsesolverreport_copy(void *_dst, void *_src, bool make_automatic);
void sparsesolverreport_free(void *_p, bool make_automatic);

void sparsespdsolvesks(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep);
void sparsespdsolve(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep);
void sparsespdcholeskysolve(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep);
void sparsesolve(sparsematrix *a, RVector *b, RVector *x, sparsesolverreport *rep);
void sparselusolve(sparsematrix *a, ZVector *p, ZVector *q, RVector *b, RVector *x, sparsesolverreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(sparsesolverreport, ae_int_t &terminationtype;);

void sparsespdsolvesks(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsespdsolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsespdcholeskysolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsesolve(const sparsematrix &a, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparselusolve(const sparsematrix &a, const integer_1d_array &p, const integer_1d_array &q, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
} // end of namespace alglib

// === ITERATIVESPARSE Package ===
// Depends on: (LinAlg) FBLS
// Depends on: DIRECTSPARSESOLVERS

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
   ae_int_t PQ;
};
void lincgstate_init(void *_p, bool make_automatic);
void lincgstate_copy(void *_dst, void *_src, bool make_automatic);
void lincgstate_free(void *_p, bool make_automatic);

struct lincgreport {
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
   double r2;
};
void lincgreport_init(void *_p, bool make_automatic);
void lincgreport_copy(void *_dst, void *_src, bool make_automatic);
void lincgreport_free(void *_p, bool make_automatic);

void lincgcreate(ae_int_t n, lincgstate *state);
void lincgsetstartingpoint(lincgstate *state, RVector *x);
void lincgsetb(lincgstate *state, RVector *b);
void lincgsetprecunit(lincgstate *state);
void lincgsetprecdiag(lincgstate *state);
void lincgsetcond(lincgstate *state, double epsf, ae_int_t maxits);
bool lincgiteration(lincgstate *state);
void lincgrestart(lincgstate *state);
void lincgsolvesparse(lincgstate *state, sparsematrix *a, bool isupper, RVector *b);
void lincgresults(lincgstate *state, RVector *x, lincgreport *rep);
void lincgsetrestartfreq(lincgstate *state, ae_int_t srf);
void lincgsetrupdatefreq(lincgstate *state, ae_int_t freq);
void lincgsetxrep(lincgstate *state, bool needxrep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(lincgstate, );
DecClass(lincgreport, ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype; double &r2;);

void lincgcreate(const ae_int_t n, lincgstate &state);
void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x);
void lincgsetb(const lincgstate &state, const real_1d_array &b);
void lincgsetprecunit(const lincgstate &state);
void lincgsetprecdiag(const lincgstate &state);
void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits);
void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b);
void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep);
void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf);
void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq);
void lincgsetxrep(const lincgstate &state, const bool needxrep);
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
   ae_int_t PQ;
};
void linlsqrstate_init(void *_p, bool make_automatic);
void linlsqrstate_copy(void *_dst, void *_src, bool make_automatic);
void linlsqrstate_free(void *_p, bool make_automatic);

struct linlsqrreport {
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
};
void linlsqrreport_init(void *_p, bool make_automatic);
void linlsqrreport_copy(void *_dst, void *_src, bool make_automatic);
void linlsqrreport_free(void *_p, bool make_automatic);

void linlsqrcreatebuf(ae_int_t m, ae_int_t n, linlsqrstate *state);
void linlsqrcreate(ae_int_t m, ae_int_t n, linlsqrstate *state);
void linlsqrsetb(linlsqrstate *state, RVector *b);
void linlsqrsetprecunit(linlsqrstate *state);
void linlsqrsetprecdiag(linlsqrstate *state);
void linlsqrsetlambdai(linlsqrstate *state, double lambdai);
bool linlsqriteration(linlsqrstate *state);
void linlsqrrestart(linlsqrstate *state);
void linlsqrsolvesparse(linlsqrstate *state, sparsematrix *a, RVector *b);
void linlsqrsetcond(linlsqrstate *state, double epsa, double epsb, ae_int_t maxits);
void linlsqrresults(linlsqrstate *state, RVector *x, linlsqrreport *rep);
void linlsqrsetxrep(linlsqrstate *state, bool needxrep);
ae_int_t linlsqrpeekiterationscount(linlsqrstate *s);
void linlsqrrequesttermination(linlsqrstate *state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(linlsqrstate, );
DecClass(linlsqrreport, ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype;);

void linlsqrcreatebuf(const ae_int_t m, const ae_int_t n, const linlsqrstate &state);
void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state);
void linlsqrsetprecunit(const linlsqrstate &state);
void linlsqrsetprecdiag(const linlsqrstate &state);
void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai);
void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b);
void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits);
void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep);
void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep);
ae_int_t linlsqrpeekiterationscount(const linlsqrstate &s);
void linlsqrrequesttermination(const linlsqrstate &state);
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
   ae_int_t PQ;
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
void nleqstate_init(void *_p, bool make_automatic);
void nleqstate_copy(void *_dst, void *_src, bool make_automatic);
void nleqstate_free(void *_p, bool make_automatic);

struct nleqreport {
   ae_int_t iterationscount;
   ae_int_t nfunc;
   ae_int_t njac;
   ae_int_t terminationtype;
};
void nleqreport_init(void *_p, bool make_automatic);
void nleqreport_copy(void *_dst, void *_src, bool make_automatic);
void nleqreport_free(void *_p, bool make_automatic);

void nleqsetcond(nleqstate *state, double epsf, ae_int_t maxits);
void nleqsetxrep(nleqstate *state, bool needxrep);
void nleqsetstpmax(nleqstate *state, double stpmax);
void nleqrestartfrom(nleqstate *state, RVector *x);
void nleqcreatelm(ae_int_t n, ae_int_t m, RVector *x, nleqstate *state);
bool nleqiteration(nleqstate *state);
void nleqresultsbuf(nleqstate *state, RVector *x, nleqreport *rep);
void nleqresults(nleqstate *state, RVector *x, nleqreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(nleqstate, bool &needf; bool &needfij; bool &xupdated; double &f; real_1d_array fi; real_2d_array j; real_1d_array x;);
DecClass(nleqreport, ae_int_t &iterationscount; ae_int_t &nfunc; ae_int_t &njac; ae_int_t &terminationtype;);

void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits);
void nleqsetxrep(const nleqstate &state, const bool needxrep);
void nleqsetstpmax(const nleqstate &state, const double stpmax);
void nleqrestartfrom(const nleqstate &state, const real_1d_array &x);
void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state);
void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state);
bool nleqiteration(const nleqstate &state);
void nleqsolve(nleqstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL);
void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep);
void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep);
} // end of namespace alglib

#endif // OnceOnly
