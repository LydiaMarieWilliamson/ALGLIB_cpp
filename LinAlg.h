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
#ifndef OnceOnlyLinAlg_h
#define OnceOnlyLinAlg_h

#include "AlgLibMisc.h"

// === ABLAS Package ===
// Depends on: (AlgLibInternal) APSERV, ABLASF, ABLASMKL
namespace alglib_impl {
ae_int_t ablasblocksize(RMatrix *a, ae_state *_state);
ae_int_t ablascomplexblocksize(CMatrix *a, ae_state *_state);
ae_int_t ablasmicroblocksize(ae_state *_state);
void ablassplitlength(RMatrix *a, ae_int_t n, ae_int_t *n1, ae_int_t *n2, ae_state *_state);
void ablascomplexsplitlength(CMatrix *a, ae_int_t n, ae_int_t *n1, ae_int_t *n2, ae_state *_state);
ae_int_t gemmparallelsize(ae_state *_state);
void generatereflection(RVector *x, ae_int_t n, double *tau, ae_state *_state);
void rmatrixger(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double alpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state);
void rmatrixgemv(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
void applyreflectionfromtheleft(RMatrix *c, double tau, RVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *work, ae_state *_state);
void applyreflectionfromtheright(RMatrix *c, double tau, RVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *work, ae_state *_state);
void rmatrixtranspose(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_state *_state);
void cmatrixtranspose(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_state *_state);
void rmatrixenforcesymmetricity(RMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
void rvectorcopy(ae_int_t n, RVector *a, ae_int_t ia, RVector *b, ae_int_t ib, ae_state *_state);
void rmatrixcopy(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_state *_state);
void cmatrixcopy(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_state *_state);
void rmatrixgencopy(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, double beta, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_state *_state);
void rmatrixrank1(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state);
void cmatrixrank1(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv, ae_state *_state);
void rmatrixmv(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, RVector *y, ae_int_t iy, ae_state *_state);
void cmatrixmv(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CVector *x, ae_int_t ix, CVector *y, ae_int_t iy, ae_state *_state);
void rmatrixsymv(ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
double rmatrixsyvmv(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, RVector *tmp, ae_state *_state);
void rmatrixtrsv(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix, ae_state *_state);
void rmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void cmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void rmatrixrighttrsm(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
void cmatrixrighttrsm(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
void rmatrixlefttrsm(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
void cmatrixlefttrsm(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
void rmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
void cmatrixherk(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
void cmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
void rowwisegramschmidt(RMatrix *q, ae_int_t m, ae_int_t n, RVector *x, RVector *qx, bool needqx, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void rmatrixger(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const double alpha, const real_1d_array &u, const ae_int_t iu, const real_1d_array &v, const ae_int_t iv, const xparams _xparams = NonTH);
void rmatrixgemv(const ae_int_t m, const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy, const xparams _xparams = NonTH);
void rmatrixtranspose(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const xparams _xparams = NonTH);
void cmatrixtranspose(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const xparams _xparams = NonTH);
void rmatrixenforcesymmetricity(const real_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
void rvectorcopy(const ae_int_t n, const real_1d_array &a, const ae_int_t ia, const real_1d_array &b, const ae_int_t ib, const xparams _xparams = NonTH);
void rmatrixcopy(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const xparams _xparams = NonTH);
void cmatrixcopy(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const xparams _xparams = NonTH);
void rmatrixgencopy(const ae_int_t m, const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const double beta, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const xparams _xparams = NonTH);
void rmatrixrank1(const ae_int_t m, const ae_int_t n, real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_1d_array &u, const ae_int_t iu, real_1d_array &v, const ae_int_t iv, const xparams _xparams = NonTH);
void cmatrixrank1(const ae_int_t m, const ae_int_t n, complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_1d_array &u, const ae_int_t iu, complex_1d_array &v, const ae_int_t iv, const xparams _xparams = NonTH);
void rmatrixmv(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, const real_1d_array &y, const ae_int_t iy, const xparams _xparams = NonTH);
void cmatrixmv(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const complex_1d_array &x, const ae_int_t ix, complex_1d_array &y, const ae_int_t iy, const xparams _xparams = NonTH);
void rmatrixsymv(const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy, const xparams _xparams = NonTH);
double rmatrixsyvmv(const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const real_1d_array &x, const ae_int_t ix, const real_1d_array &tmp, const xparams _xparams = NonTH);
void rmatrixtrsv(const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, const ae_int_t ix, const xparams _xparams = NonTH);
void rmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc, const xparams _xparams = NonTH);
void cmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const complex alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const complex beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const xparams _xparams = NonTH);
void rmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2, const xparams _xparams = NonTH);
void cmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2, const xparams _xparams = NonTH);
void rmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2, const xparams _xparams = NonTH);
void cmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2, const xparams _xparams = NonTH);
void rmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper, const xparams _xparams = NonTH);
void cmatrixherk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper, const xparams _xparams = NonTH);
void cmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper, const xparams _xparams = NonTH);
} // end of namespace alglib

// === ORTFAC Package ===
// Depends on: (AlgLibInternal) HBLAS, CREFLECTIONS, SBLAS
// Depends on: (AlgLibMisc) HQRND
// Depends on: ABLAS
namespace alglib_impl {
void rmatrixqrbasecase(RMatrix *a, ae_int_t m, ae_int_t n, RVector *work, RVector *t, RVector *tau, ae_state *_state);
void rmatrixlqbasecase(RMatrix *a, ae_int_t m, ae_int_t n, RVector *work, RVector *t, RVector *tau, ae_state *_state);
void rmatrixqr(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau, ae_state *_state);
void cmatrixqr(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau, ae_state *_state);
void rmatrixlq(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau, ae_state *_state);
void cmatrixlq(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau, ae_state *_state);
void rmatrixqrunpackq(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau, ae_int_t qcolumns, RMatrix *q, ae_state *_state);
void cmatrixqrunpackq(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau, ae_int_t qcolumns, CMatrix *q, ae_state *_state);
void rmatrixqrunpackr(RMatrix *a, ae_int_t m, ae_int_t n, RMatrix *r, ae_state *_state);
void cmatrixqrunpackr(CMatrix *a, ae_int_t m, ae_int_t n, CMatrix *r, ae_state *_state);
void rmatrixlqunpackq(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau, ae_int_t qrows, RMatrix *q, ae_state *_state);
void cmatrixlqunpackq(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau, ae_int_t qrows, CMatrix *q, ae_state *_state);
void rmatrixlqunpackl(RMatrix *a, ae_int_t m, ae_int_t n, RMatrix *l, ae_state *_state);
void cmatrixlqunpackl(CMatrix *a, ae_int_t m, ae_int_t n, CMatrix *l, ae_state *_state);
void rmatrixbd(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tauq, RVector *taup, ae_state *_state);
void rmatrixbdmultiplybyq(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool fromtheright, bool dotranspose, ae_state *_state);
void rmatrixbdunpackq(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, ae_int_t qcolumns, RMatrix *q, ae_state *_state);
void rmatrixbdmultiplybyp(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *taup, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool fromtheright, bool dotranspose, ae_state *_state);
void rmatrixbdunpackpt(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *taup, ae_int_t ptrows, RMatrix *pt, ae_state *_state);
void rmatrixbdunpackdiagonals(RMatrix *b, ae_int_t m, ae_int_t n, bool *isupper, RVector *d, RVector *e, ae_state *_state);
void rmatrixhessenberg(RMatrix *a, ae_int_t n, RVector *tau, ae_state *_state);
void rmatrixhessenbergunpackq(RMatrix *a, ae_int_t n, RVector *tau, RMatrix *q, ae_state *_state);
void rmatrixhessenbergunpackh(RMatrix *a, ae_int_t n, RMatrix *h, ae_state *_state);
void smatrixtd(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RVector *d, RVector *e, ae_state *_state);
void hmatrixtd(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, RVector *d, RVector *e, ae_state *_state);
void smatrixtdunpackq(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RMatrix *q, ae_state *_state);
void hmatrixtdunpackq(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, CMatrix *q, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void rmatrixqr(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau, const xparams _xparams = NonTH);
void cmatrixqr(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau, const xparams _xparams = NonTH);
void rmatrixlq(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau, const xparams _xparams = NonTH);
void cmatrixlq(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau, const xparams _xparams = NonTH);
void rmatrixqrunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qcolumns, real_2d_array &q, const xparams _xparams = NonTH);
void cmatrixqrunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qcolumns, complex_2d_array &q, const xparams _xparams = NonTH);
void rmatrixqrunpackr(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &r, const xparams _xparams = NonTH);
void cmatrixqrunpackr(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &r, const xparams _xparams = NonTH);
void rmatrixlqunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qrows, real_2d_array &q, const xparams _xparams = NonTH);
void cmatrixlqunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qrows, complex_2d_array &q, const xparams _xparams = NonTH);
void rmatrixlqunpackl(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &l, const xparams _xparams = NonTH);
void cmatrixlqunpackl(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &l, const xparams _xparams = NonTH);
void rmatrixbd(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tauq, real_1d_array &taup, const xparams _xparams = NonTH);
void rmatrixbdmultiplybyq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose, const xparams _xparams = NonTH);
void rmatrixbdunpackq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, const ae_int_t qcolumns, real_2d_array &q, const xparams _xparams = NonTH);
void rmatrixbdmultiplybyp(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose, const xparams _xparams = NonTH);
void rmatrixbdunpackpt(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, const ae_int_t ptrows, real_2d_array &pt, const xparams _xparams = NonTH);
void rmatrixbdunpackdiagonals(const real_2d_array &b, const ae_int_t m, const ae_int_t n, bool &isupper, real_1d_array &d, real_1d_array &e, const xparams _xparams = NonTH);
void rmatrixhessenberg(real_2d_array &a, const ae_int_t n, real_1d_array &tau, const xparams _xparams = NonTH);
void rmatrixhessenbergunpackq(const real_2d_array &a, const ae_int_t n, const real_1d_array &tau, real_2d_array &q, const xparams _xparams = NonTH);
void rmatrixhessenbergunpackh(const real_2d_array &a, const ae_int_t n, real_2d_array &h, const xparams _xparams = NonTH);
void smatrixtd(real_2d_array &a, const ae_int_t n, const bool isupper, real_1d_array &tau, real_1d_array &d, real_1d_array &e, const xparams _xparams = NonTH);
void hmatrixtd(complex_2d_array &a, const ae_int_t n, const bool isupper, complex_1d_array &tau, real_1d_array &d, real_1d_array &e, const xparams _xparams = NonTH);
void smatrixtdunpackq(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &tau, real_2d_array &q, const xparams _xparams = NonTH);
void hmatrixtdunpackq(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &tau, complex_2d_array &q, const xparams _xparams = NonTH);
} // end of namespace alglib

// === MATGEN Package ===
// Depends on: (AlgLibInternal) CREFLECTIONS
// Depends on: (AlgLibMisc) HQRND
// Depends on: ABLAS
namespace alglib_impl {
void rmatrixrndorthogonalfromtheright(RMatrix *a, ae_int_t m, ae_int_t n, ae_state *_state);
void cmatrixrndorthogonalfromtheright(CMatrix *a, ae_int_t m, ae_int_t n, ae_state *_state);
void rmatrixrndorthogonalfromtheleft(RMatrix *a, ae_int_t m, ae_int_t n, ae_state *_state);
void cmatrixrndorthogonalfromtheleft(CMatrix *a, ae_int_t m, ae_int_t n, ae_state *_state);
void rmatrixrndorthogonal(ae_int_t n, RMatrix *a, ae_state *_state);
void cmatrixrndorthogonal(ae_int_t n, CMatrix *a, ae_state *_state);
void rmatrixrndcond(ae_int_t n, double c, RMatrix *a, ae_state *_state);
void cmatrixrndcond(ae_int_t n, double c, CMatrix *a, ae_state *_state);
void smatrixrndmultiply(RMatrix *a, ae_int_t n, ae_state *_state);
void hmatrixrndmultiply(CMatrix *a, ae_int_t n, ae_state *_state);
void smatrixrndcond(ae_int_t n, double c, RMatrix *a, ae_state *_state);
void hmatrixrndcond(ae_int_t n, double c, CMatrix *a, ae_state *_state);
void spdmatrixrndcond(ae_int_t n, double c, RMatrix *a, ae_state *_state);
void hpdmatrixrndcond(ae_int_t n, double c, CMatrix *a, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void rmatrixrndorthogonalfromtheright(real_2d_array &a, const ae_int_t m, const ae_int_t n, const xparams _xparams = NonTH);
void cmatrixrndorthogonalfromtheright(complex_2d_array &a, const ae_int_t m, const ae_int_t n, const xparams _xparams = NonTH);
void rmatrixrndorthogonalfromtheleft(real_2d_array &a, const ae_int_t m, const ae_int_t n, const xparams _xparams = NonTH);
void cmatrixrndorthogonalfromtheleft(complex_2d_array &a, const ae_int_t m, const ae_int_t n, const xparams _xparams = NonTH);
void rmatrixrndorthogonal(const ae_int_t n, real_2d_array &a, const xparams _xparams = NonTH);
void cmatrixrndorthogonal(const ae_int_t n, complex_2d_array &a, const xparams _xparams = NonTH);
void rmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a, const xparams _xparams = NonTH);
void cmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a, const xparams _xparams = NonTH);
void smatrixrndmultiply(real_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
void hmatrixrndmultiply(complex_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
void smatrixrndcond(const ae_int_t n, const double c, real_2d_array &a, const xparams _xparams = NonTH);
void hmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a, const xparams _xparams = NonTH);
void spdmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a, const xparams _xparams = NonTH);
void hpdmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a, const xparams _xparams = NonTH);
} // end of namespace alglib

// === SPARSE Package ===
// Depends on: (AlgLibInternal) ABLASMKL, SCODES, TSORT
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
struct sparsematrix {
   ae_vector vals;
   ae_vector idx;
   ae_vector ridx;
   ae_vector didx;
   ae_vector uidx;
   ae_int_t matrixtype;
   ae_int_t m;
   ae_int_t n;
   ae_int_t nfree;
   ae_int_t ninitialized;
   ae_int_t tablesize;
};
void sparsematrix_init(void *_p, ae_state *_state, bool make_automatic);
void sparsematrix_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sparsematrix_free(void *_p, bool make_automatic);
void sparsealloc(ae_serializer *s, sparsematrix *a, ae_state *_state);
void sparseserialize(ae_serializer *s, sparsematrix *a, ae_state *_state);
void sparseunserialize(ae_serializer *s, sparsematrix *a, ae_state *_state);

struct sparsebuffers {
   ae_vector d;
   ae_vector u;
   sparsematrix s;
};
void sparsebuffers_init(void *_p, ae_state *_state, bool make_automatic);
void sparsebuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sparsebuffers_free(void *_p, bool make_automatic);

void sparseinitduidx(sparsematrix *s, ae_state *_state);
void sparsecreatebuf(ae_int_t m, ae_int_t n, ae_int_t k, sparsematrix *s, ae_state *_state);
void sparsecreate(ae_int_t m, ae_int_t n, ae_int_t k, sparsematrix *s, ae_state *_state);
void sparsecreatecrsbuf(ae_int_t m, ae_int_t n, ZVector *ner, sparsematrix *s, ae_state *_state);
void sparsecreatecrs(ae_int_t m, ae_int_t n, ZVector *ner, sparsematrix *s, ae_state *_state);
void sparsecreatesksbuf(ae_int_t m, ae_int_t n, ZVector *d, ZVector *u, sparsematrix *s, ae_state *_state);
void sparsecreatesks(ae_int_t m, ae_int_t n, ZVector *d, ZVector *u, sparsematrix *s, ae_state *_state);
void sparsecreatesksbandbuf(ae_int_t m, ae_int_t n, ae_int_t bw, sparsematrix *s, ae_state *_state);
void sparsecreatesksband(ae_int_t m, ae_int_t n, ae_int_t bw, sparsematrix *s, ae_state *_state);
void sparsecopybuf(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopy(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparseswap(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
bool sparserewriteexisting(sparsematrix *s, ae_int_t i, ae_int_t j, double v, ae_state *_state);
void sparseresizematrix(sparsematrix *s, ae_state *_state);
void sparseset(sparsematrix *s, ae_int_t i, ae_int_t j, double v, ae_state *_state);
void sparseadd(sparsematrix *s, ae_int_t i, ae_int_t j, double v, ae_state *_state);
double sparseget(sparsematrix *s, ae_int_t i, ae_int_t j, ae_state *_state);
bool sparseexists(sparsematrix *s, ae_int_t i, ae_int_t j, ae_state *_state);
double sparsegetdiagonal(sparsematrix *s, ae_int_t i, ae_state *_state);
void sparsemv(sparsematrix *s, RVector *x, RVector *y, ae_state *_state);
void sparsemtv(sparsematrix *s, RVector *x, RVector *y, ae_state *_state);
void sparsegemv(sparsematrix *s, double alpha, ae_int_t ops, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
void sparsemv2(sparsematrix *s, RVector *x, RVector *y0, RVector *y1, ae_state *_state);
void sparsesmv(sparsematrix *s, bool isupper, RVector *x, RVector *y, ae_state *_state);
double sparsevsmv(sparsematrix *s, bool isupper, RVector *x, ae_state *_state);
void sparsemm(sparsematrix *s, RMatrix *a, ae_int_t k, RMatrix *b, ae_state *_state);
void sparsemtm(sparsematrix *s, RMatrix *a, ae_int_t k, RMatrix *b, ae_state *_state);
void sparsemm2(sparsematrix *s, RMatrix *a, ae_int_t k, RMatrix *b0, RMatrix *b1, ae_state *_state);
void sparsesmm(sparsematrix *s, bool isupper, RMatrix *a, ae_int_t k, RMatrix *b, ae_state *_state);
void sparsetrmv(sparsematrix *s, bool isupper, bool isunit, ae_int_t optype, RVector *x, RVector *y, ae_state *_state);
void sparsetrsv(sparsematrix *s, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_state *_state);
void sparsesymmpermtblbuf(sparsematrix *a, bool isupper, ZVector *p, sparsematrix *b, ae_state *_state);
void sparsesymmpermtbl(sparsematrix *a, bool isupper, ZVector *p, sparsematrix *b, ae_state *_state);
double sparsegetaveragelengthofchain(sparsematrix *s, ae_state *_state);
bool sparseenumerate(sparsematrix *s, ae_int_t *t0, ae_int_t *t1, ae_int_t *i, ae_int_t *j, double *v, ae_state *_state);
void sparsegetrow(sparsematrix *s, ae_int_t i, RVector *irow, ae_state *_state);
void sparsegetcompressedrow(sparsematrix *s, ae_int_t i, ZVector *colidx, RVector *vals, ae_int_t *nzcnt, ae_state *_state);
void sparsetransposesks(sparsematrix *s, ae_state *_state);
void sparsetransposecrs(sparsematrix *s, ae_state *_state);
void sparsecopytransposecrsbuf(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytransposecrs(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparseconverttohash(sparsematrix *s, ae_state *_state);
void sparseconverttocrs(sparsematrix *s, ae_state *_state);
void sparseconverttosks(sparsematrix *s, ae_state *_state);
void sparseconvertto(sparsematrix *s0, ae_int_t fmt, ae_state *_state);
void sparsecopytohashbuf(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytohash(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytocrsbuf(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytocrs(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytosksbuf(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytosks(sparsematrix *s0, sparsematrix *s1, ae_state *_state);
void sparsecopytobuf(sparsematrix *s0, ae_int_t fmt, sparsematrix *s1, ae_state *_state);
void sparsecreatecrsinplace(sparsematrix *s, ae_state *_state);
ae_int_t sparsegetmatrixtype(sparsematrix *s, ae_state *_state);
bool sparseishash(sparsematrix *s, ae_state *_state);
bool sparseiscrs(sparsematrix *s, ae_state *_state);
bool sparseissks(sparsematrix *s, ae_state *_state);
void sparsefree(sparsematrix *s, ae_state *_state);
ae_int_t sparsegetncols(sparsematrix *s, ae_state *_state);
ae_int_t sparsegetnrows(sparsematrix *s, ae_state *_state);
ae_int_t sparsegetuppercount(sparsematrix *s, ae_state *_state);
ae_int_t sparsegetlowercount(sparsematrix *s, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(sparsematrix, );
DecClass(sparsebuffers, );
void sparseserialize(sparsematrix &obj, std::string &s_out);
void sparseserialize(sparsematrix &obj, std::ostream &s_out);
void sparseunserialize(const std::string &s_in, sparsematrix &obj);
void sparseunserialize(const std::istream &s_in, sparsematrix &obj);

void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const ae_int_t k, const sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreate(const ae_int_t m, const ae_int_t n, const ae_int_t k, sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreate(const ae_int_t m, const ae_int_t n, sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatecrsbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, const sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatecrs(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatesksbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, const sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatesks(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatesksbandbuf(const ae_int_t m, const ae_int_t n, const ae_int_t bw, const sparsematrix &s, const xparams _xparams = NonTH);
void sparsecreatesksband(const ae_int_t m, const ae_int_t n, const ae_int_t bw, sparsematrix &s, const xparams _xparams = NonTH);
void sparsecopybuf(const sparsematrix &s0, const sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopy(const sparsematrix &s0, sparsematrix &s1, const xparams _xparams = NonTH);
void sparseswap(const sparsematrix &s0, const sparsematrix &s1, const xparams _xparams = NonTH);
bool sparserewriteexisting(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v, const xparams _xparams = NonTH);
void sparseresizematrix(const sparsematrix &s, const xparams _xparams = NonTH);
void sparseset(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v, const xparams _xparams = NonTH);
void sparseadd(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v, const xparams _xparams = NonTH);
double sparseget(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const xparams _xparams = NonTH);
bool sparseexists(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const xparams _xparams = NonTH);
double sparsegetdiagonal(const sparsematrix &s, const ae_int_t i, const xparams _xparams = NonTH);
void sparsemv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void sparsemtv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void sparsegemv(const sparsematrix &s, const double alpha, const ae_int_t ops, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy, const xparams _xparams = NonTH);
void sparsemv2(const sparsematrix &s, const real_1d_array &x, real_1d_array &y0, real_1d_array &y1, const xparams _xparams = NonTH);
void sparsesmv(const sparsematrix &s, const bool isupper, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
double sparsevsmv(const sparsematrix &s, const bool isupper, const real_1d_array &x, const xparams _xparams = NonTH);
void sparsemm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b, const xparams _xparams = NonTH);
void sparsemtm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b, const xparams _xparams = NonTH);
void sparsemm2(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b0, real_2d_array &b1, const xparams _xparams = NonTH);
void sparsesmm(const sparsematrix &s, const bool isupper, const real_2d_array &a, const ae_int_t k, real_2d_array &b, const xparams _xparams = NonTH);
void sparsetrmv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void sparsetrsv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, const xparams _xparams = NonTH);
void sparsesymmpermtblbuf(const sparsematrix &a, const bool isupper, const integer_1d_array &p, const sparsematrix &b, const xparams _xparams = NonTH);
void sparsesymmpermtbl(const sparsematrix &a, const bool isupper, const integer_1d_array &p, sparsematrix &b, const xparams _xparams = NonTH);
bool sparseenumerate(const sparsematrix &s, ae_int_t &t0, ae_int_t &t1, ae_int_t &i, ae_int_t &j, double &v, const xparams _xparams = NonTH);
void sparsegetrow(const sparsematrix &s, const ae_int_t i, real_1d_array &irow, const xparams _xparams = NonTH);
void sparsegetcompressedrow(const sparsematrix &s, const ae_int_t i, integer_1d_array &colidx, real_1d_array &vals, ae_int_t &nzcnt, const xparams _xparams = NonTH);
void sparsetransposesks(const sparsematrix &s, const xparams _xparams = NonTH);
void sparsetransposecrs(const sparsematrix &s, const xparams _xparams = NonTH);
void sparsecopytransposecrsbuf(const sparsematrix &s0, const sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytransposecrs(const sparsematrix &s0, sparsematrix &s1, const xparams _xparams = NonTH);
void sparseconverttohash(const sparsematrix &s, const xparams _xparams = NonTH);
void sparseconverttocrs(const sparsematrix &s, const xparams _xparams = NonTH);
void sparseconverttosks(const sparsematrix &s, const xparams _xparams = NonTH);
void sparseconvertto(const sparsematrix &s0, const ae_int_t fmt, const xparams _xparams = NonTH);
void sparsecopytohashbuf(const sparsematrix &s0, const sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytohash(const sparsematrix &s0, sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytocrsbuf(const sparsematrix &s0, const sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytocrs(const sparsematrix &s0, sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytosksbuf(const sparsematrix &s0, const sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytosks(const sparsematrix &s0, sparsematrix &s1, const xparams _xparams = NonTH);
void sparsecopytobuf(const sparsematrix &s0, const ae_int_t fmt, const sparsematrix &s1, const xparams _xparams = NonTH);
ae_int_t sparsegetmatrixtype(const sparsematrix &s, const xparams _xparams = NonTH);
bool sparseishash(const sparsematrix &s, const xparams _xparams = NonTH);
bool sparseiscrs(const sparsematrix &s, const xparams _xparams = NonTH);
bool sparseissks(const sparsematrix &s, const xparams _xparams = NonTH);
void sparsefree(sparsematrix &s, const xparams _xparams = NonTH);
ae_int_t sparsegetncols(const sparsematrix &s, const xparams _xparams = NonTH);
ae_int_t sparsegetnrows(const sparsematrix &s, const xparams _xparams = NonTH);
ae_int_t sparsegetuppercount(const sparsematrix &s, const xparams _xparams = NonTH);
ae_int_t sparsegetlowercount(const sparsematrix &s, const xparams _xparams = NonTH);
} // end of namespace alglib

// === HSSCHUR Package ===
// Depends on: (AlgLibInternal) BLAS, ROTATIONS
// Depends on: ABLAS
namespace alglib_impl {
void internalschurdecomposition(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info, ae_state *_state);
void rmatrixinternalschurdecomposition(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info, ae_state *_state);
bool upperhessenbergschurdecomposition(RMatrix *h, ae_int_t n, RMatrix *s, ae_state *_state);
} // end of namespace alglib_impl

// === EVD Package ===
// Depends on: (AlgLibInternal) BASICSTATOPS
// Depends on: ORTFAC, MATGEN, SPARSE, HSSCHUR
namespace alglib_impl {
struct eigsubspacestate {
   ae_int_t n;
   ae_int_t k;
   ae_int_t nwork;
   ae_int_t maxits;
   double eps;
   ae_int_t eigenvectorsneeded;
   ae_int_t matrixtype;
   bool usewarmstart;
   bool firstcall;
   hqrndstate rs;
   bool running;
   ae_vector tau;
   ae_matrix q0;
   ae_matrix qcur;
   ae_matrix qnew;
   ae_matrix znew;
   ae_matrix r;
   ae_matrix rz;
   ae_matrix tz;
   ae_matrix rq;
   ae_matrix dummy;
   ae_vector rw;
   ae_vector tw;
   ae_vector wcur;
   ae_vector wprev;
   ae_vector wrank;
   apbuffers buf;
   ae_matrix x;
   ae_matrix ax;
   ae_int_t requesttype;
   ae_int_t requestsize;
   ae_int_t repiterationscount;
   rcommstate rstate;
};
void eigsubspacestate_init(void *_p, ae_state *_state, bool make_automatic);
void eigsubspacestate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void eigsubspacestate_free(void *_p, bool make_automatic);

struct eigsubspacereport {
   ae_int_t iterationscount;
};
void eigsubspacereport_init(void *_p, ae_state *_state, bool make_automatic);
void eigsubspacereport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void eigsubspacereport_free(void *_p, bool make_automatic);

bool rmatrixevd(RMatrix *a, ae_int_t n, ae_int_t vneeded, RVector *wr, RVector *wi, RMatrix *vl, RMatrix *vr, ae_state *_state);
bool smatrixtdevd(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z, ae_state *_state);
bool smatrixtdevdr(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, double a, double b, ae_int_t *m, RMatrix *z, ae_state *_state);
bool smatrixtdevdi(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, ae_int_t i1, ae_int_t i2, RMatrix *z, ae_state *_state);
bool smatrixevd(RMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, RVector *d, RMatrix *z, ae_state *_state);
bool hmatrixevd(CMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, RVector *d, CMatrix *z, ae_state *_state);
bool smatrixevdr(RMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, double b1, double b2, ae_int_t *m, RVector *w, RMatrix *z, ae_state *_state);
bool hmatrixevdr(CMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, double b1, double b2, ae_int_t *m, RVector *w, CMatrix *z, ae_state *_state);
bool smatrixevdi(RMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, ae_int_t i1, ae_int_t i2, RVector *w, RMatrix *z, ae_state *_state);
bool hmatrixevdi(CMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, ae_int_t i1, ae_int_t i2, RVector *w, CMatrix *z, ae_state *_state);
void eigsubspacesetcond(eigsubspacestate *state, double eps, ae_int_t maxits, ae_state *_state);
void eigsubspacecreatebuf(ae_int_t n, ae_int_t k, eigsubspacestate *state, ae_state *_state);
void eigsubspacecreate(ae_int_t n, ae_int_t k, eigsubspacestate *state, ae_state *_state);
void eigsubspacesetwarmstart(eigsubspacestate *state, bool usewarmstart, ae_state *_state);
void eigsubspaceoocstart(eigsubspacestate *state, ae_int_t mtype, ae_state *_state);
bool eigsubspaceiteration(eigsubspacestate *state, ae_state *_state);
bool eigsubspaceooccontinue(eigsubspacestate *state, ae_state *_state);
void eigsubspaceoocgetrequestinfo(eigsubspacestate *state, ae_int_t *requesttype, ae_int_t *requestsize, ae_state *_state);
void eigsubspaceoocgetrequestdata(eigsubspacestate *state, RMatrix *x, ae_state *_state);
void eigsubspaceoocsendresult(eigsubspacestate *state, RMatrix *ax, ae_state *_state);
void eigsubspaceoocstop(eigsubspacestate *state, RVector *w, RMatrix *z, eigsubspacereport *rep, ae_state *_state);
void eigsubspacesolvedenses(eigsubspacestate *state, RMatrix *a, bool isupper, RVector *w, RMatrix *z, eigsubspacereport *rep, ae_state *_state);
void eigsubspacesolvesparses(eigsubspacestate *state, sparsematrix *a, bool isupper, RVector *w, RMatrix *z, eigsubspacereport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(eigsubspacestate, );
DecClass(eigsubspacereport, ae_int_t &iterationscount;);

bool rmatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t vneeded, real_1d_array &wr, real_1d_array &wi, real_2d_array &vl, real_2d_array &vr, const xparams _xparams = NonTH);
bool smatrixtdevd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, real_2d_array &z, const xparams _xparams = NonTH);
bool smatrixtdevdr(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const double a, const double b, ae_int_t &m, real_2d_array &z, const xparams _xparams = NonTH);
bool smatrixtdevdi(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const ae_int_t i1, const ae_int_t i2, real_2d_array &z, const xparams _xparams = NonTH);
bool smatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, real_2d_array &z, const xparams _xparams = NonTH);
bool hmatrixevd(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, complex_2d_array &z, const xparams _xparams = NonTH);
bool smatrixevdr(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, real_2d_array &z, const xparams _xparams = NonTH);
bool hmatrixevdr(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, complex_2d_array &z, const xparams _xparams = NonTH);
bool smatrixevdi(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, real_2d_array &z, const xparams _xparams = NonTH);
bool hmatrixevdi(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, complex_2d_array &z, const xparams _xparams = NonTH);
void eigsubspacesetcond(const eigsubspacestate &state, const double eps, const ae_int_t maxits, const xparams _xparams = NonTH);
void eigsubspacecreatebuf(const ae_int_t n, const ae_int_t k, const eigsubspacestate &state, const xparams _xparams = NonTH);
void eigsubspacecreate(const ae_int_t n, const ae_int_t k, eigsubspacestate &state, const xparams _xparams = NonTH);
void eigsubspacesetwarmstart(const eigsubspacestate &state, const bool usewarmstart, const xparams _xparams = NonTH);
void eigsubspaceoocstart(const eigsubspacestate &state, const ae_int_t mtype, const xparams _xparams = NonTH);
bool eigsubspaceooccontinue(const eigsubspacestate &state, const xparams _xparams = NonTH);
void eigsubspaceoocgetrequestinfo(const eigsubspacestate &state, ae_int_t &requesttype, ae_int_t &requestsize, const xparams _xparams = NonTH);
void eigsubspaceoocgetrequestdata(const eigsubspacestate &state, real_2d_array &x, const xparams _xparams = NonTH);
void eigsubspaceoocsendresult(const eigsubspacestate &state, const real_2d_array &ax, const xparams _xparams = NonTH);
void eigsubspaceoocstop(const eigsubspacestate &state, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep, const xparams _xparams = NonTH);
void eigsubspacesolvedenses(const eigsubspacestate &state, const real_2d_array &a, const bool isupper, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep, const xparams _xparams = NonTH);
void eigsubspacesolvesparses(const eigsubspacestate &state, const sparsematrix &a, const bool isupper, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === DLU Package ===
// Depends on: ABLAS
namespace alglib_impl {
void rmatrixluprec(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, RVector *tmp, ae_state *_state);
void cmatrixluprec(CMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, CVector *tmp, ae_state *_state);
void rmatrixplurec(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, RVector *tmp, ae_state *_state);
void cmatrixplurec(CMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, CVector *tmp, ae_state *_state);
} // end of namespace alglib_impl

// === SPTRF Package ===
// Depends on: SPARSE, DLU
namespace alglib_impl {
struct sluv2list1matrix {
   ae_int_t nfixed;
   ae_int_t ndynamic;
   ae_vector idxfirst;
   ae_vector strgidx;
   ae_vector strgval;
   ae_int_t nallocated;
   ae_int_t nused;
};
void sluv2list1matrix_init(void *_p, ae_state *_state, bool make_automatic);
void sluv2list1matrix_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sluv2list1matrix_free(void *_p, bool make_automatic);

struct sluv2sparsetrail {
   ae_int_t n;
   ae_int_t k;
   ae_vector nzc;
   ae_int_t maxwrkcnt;
   ae_int_t maxwrknz;
   ae_int_t wrkcnt;
   ae_vector wrkset;
   ae_vector colid;
   ae_vector isdensified;
   ae_vector slscolptr;
   ae_vector slsrowptr;
   ae_vector slsidx;
   ae_vector slsval;
   ae_int_t slsused;
   ae_vector tmp0;
};
void sluv2sparsetrail_init(void *_p, ae_state *_state, bool make_automatic);
void sluv2sparsetrail_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sluv2sparsetrail_free(void *_p, bool make_automatic);

struct sluv2densetrail {
   ae_int_t n;
   ae_int_t ndense;
   ae_matrix d;
   ae_vector did;
};
void sluv2densetrail_init(void *_p, ae_state *_state, bool make_automatic);
void sluv2densetrail_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sluv2densetrail_free(void *_p, bool make_automatic);

struct sluv2buffer {
   ae_int_t n;
   sparsematrix sparsel;
   sparsematrix sparseut;
   sluv2list1matrix bleft;
   sluv2list1matrix bupper;
   sluv2sparsetrail strail;
   sluv2densetrail dtrail;
   ae_vector rowpermrawidx;
   ae_matrix dbuf;
   ae_vector v0i;
   ae_vector v1i;
   ae_vector v0r;
   ae_vector v1r;
   ae_vector tmp0;
   ae_vector tmpi;
   ae_vector tmpp;
};
void sluv2buffer_init(void *_p, ae_state *_state, bool make_automatic);
void sluv2buffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sluv2buffer_free(void *_p, bool make_automatic);

bool sptrflu(sparsematrix *a, ae_int_t pivottype, ZVector *pr, ZVector *pc, sluv2buffer *buf, ae_state *_state);
} // end of namespace alglib_impl

// === AMDORDERING Package ===
// Depends on: ABLAS, SPARSE
namespace alglib_impl {
struct amdnset {
   ae_int_t n;
   ae_int_t nstored;
   ae_vector items;
   ae_vector locationof;
   ae_int_t iteridx;
};
void amdnset_init(void *_p, ae_state *_state, bool make_automatic);
void amdnset_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void amdnset_free(void *_p, bool make_automatic);

struct amdknset {
   ae_int_t k;
   ae_int_t n;
   ae_vector flagarray;
   ae_vector vbegin;
   ae_vector vallocated;
   ae_vector vcnt;
   ae_vector data;
   ae_int_t dataused;
   ae_int_t iterrow;
   ae_int_t iteridx;
};
void amdknset_init(void *_p, ae_state *_state, bool make_automatic);
void amdknset_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void amdknset_free(void *_p, bool make_automatic);

struct amdvertexset {
   ae_int_t n;
   bool checkexactdegrees;
   ae_int_t smallestdegree;
   ae_vector approxd;
   ae_vector optionalexactd;
   ae_vector isvertex;
   ae_vector vbegin;
   ae_vector vprev;
   ae_vector vnext;
};
void amdvertexset_init(void *_p, ae_state *_state, bool make_automatic);
void amdvertexset_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void amdvertexset_free(void *_p, bool make_automatic);

struct amdllmatrix {
   ae_int_t n;
   ae_vector vbegin;
   ae_vector vcolcnt;
   ae_vector entries;
   ae_int_t entriesinitialized;
};
void amdllmatrix_init(void *_p, ae_state *_state, bool make_automatic);
void amdllmatrix_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void amdllmatrix_free(void *_p, bool make_automatic);

struct amdbuffer {
   ae_int_t n;
   bool checkexactdegrees;
   ae_vector iseliminated;
   ae_vector issupernode;
   amdknset setsuper;
   amdknset seta;
   amdknset sete;
   amdllmatrix mtxl;
   amdvertexset vertexdegrees;
   amdnset setq;
   ae_vector perm;
   ae_vector invperm;
   ae_vector columnswaps;
   amdnset setp;
   amdnset lp;
   amdnset setrp;
   amdnset ep;
   amdnset adji;
   amdnset adjj;
   ae_vector ls;
   ae_int_t lscnt;
   amdnset setqsupercand;
   amdnset exactdegreetmp0;
   amdknset hashbuckets;
   amdnset nonemptybuckets;
   ae_vector sncandidates;
   ae_vector tmp0;
   ae_vector arrwe;
   ae_matrix dbga;
};
void amdbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void amdbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void amdbuffer_free(void *_p, bool make_automatic);

ae_int_t generateamdpermutationx(sparsematrix *a, ae_int_t n, ZVector *perm, ZVector *invperm, ae_int_t amdtype, amdbuffer *buf, ae_state *_state);
void generateamdpermutation(sparsematrix *a, ae_int_t n, ZVector *perm, ZVector *invperm, amdbuffer *buf, ae_state *_state);
} // end of namespace alglib_impl

// === SPCHOL Package ===
// Depends on: AMDORDERING
namespace alglib_impl {
struct spcholanalysis {
   ae_int_t tasktype;
   ae_int_t n;
   ae_int_t permtype;
   bool unitd;
   ae_int_t modtype;
   double modparam0;
   double modparam1;
   double modparam2;
   double modparam3;
   ae_int_t nsuper;
   ae_vector parentsupernode;
   ae_vector supercolrange;
   ae_vector superrowridx;
   ae_vector superrowidx;
   ae_vector fillinperm;
   ae_vector invfillinperm;
   ae_vector superperm;
   ae_vector invsuperperm;
   ae_vector effectiveperm;
   ae_vector inveffectiveperm;
   bool istopologicalordering;
   bool applypermutationtooutput;
   ae_vector ladjplusr;
   ae_vector ladjplus;
   ae_vector outrowcounts;
   ae_vector inputstorage;
   ae_vector outputstorage;
   ae_vector rowstrides;
   ae_vector rowoffsets;
   ae_vector diagd;
   ae_vector wrkrows;
   ae_vector flagarray;
   ae_vector tmpparent;
   ae_vector node2supernode;
   ae_vector u2smap;
   ae_vector raw2smap;
   amdbuffer amdtmp;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmp2;
   ae_vector tmp3;
   ae_vector tmp4;
   sparsematrix tmpa;
   sparsematrix tmpat;
   sparsematrix tmpa2;
   sparsematrix tmpbottomt;
   sparsematrix tmpupdate;
   sparsematrix tmpupdatet;
   sparsematrix tmpnewtailt;
   ae_vector tmpperm;
   ae_vector invtmpperm;
   ae_vector tmpx;
   ae_vector simdbuf;
};
void spcholanalysis_init(void *_p, ae_state *_state, bool make_automatic);
void spcholanalysis_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spcholanalysis_free(void *_p, bool make_automatic);

ae_int_t spsymmgetmaxfastkernel(ae_state *_state);
bool spsymmanalyze(sparsematrix *a, ae_int_t facttype, ae_int_t permtype, spcholanalysis *analysis, ae_state *_state);
void spsymmsetmodificationstrategy(spcholanalysis *analysis, ae_int_t modstrategy, double p0, double p1, double p2, double p3, ae_state *_state);
void spsymmreload(spcholanalysis *analysis, sparsematrix *a, ae_state *_state);
void spsymmreloaddiagonal(spcholanalysis *analysis, RVector *d, ae_state *_state);
bool spsymmfactorize(spcholanalysis *analysis, ae_state *_state);
void spsymmextract(spcholanalysis *analysis, sparsematrix *a, RVector *d, ZVector *p, ae_state *_state);
void spsymmsolve(spcholanalysis *analysis, RVector *b, ae_state *_state);
void spsymmdiagerr(spcholanalysis *analysis, double *sumsq, double *errsq, ae_state *_state);
} // end of namespace alglib_impl

// === TRFAC Package ===
// Depends on: (AlgLibInternal) ROTATIONS
// Depends on: MATGEN, SPTRF, SPCHOL
namespace alglib_impl {
struct sparsedecompositionanalysis {
   ae_int_t n;
   ae_int_t facttype;
   ae_int_t permtype;
   spcholanalysis analysis;
   sparsematrix wrka;
   sparsematrix wrkat;
   sparsematrix crsa;
   sparsematrix crsat;
};
void sparsedecompositionanalysis_init(void *_p, ae_state *_state, bool make_automatic);
void sparsedecompositionanalysis_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sparsedecompositionanalysis_free(void *_p, bool make_automatic);

void rmatrixlup(RMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
void cmatrixlup(CMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
void rmatrixplu(RMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
void cmatrixplu(CMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
void rmatrixlu(RMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
void cmatrixlu(CMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
bool spdmatrixcholeskyrec(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, RVector *tmp, ae_state *_state);
bool spdmatrixcholesky(RMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
bool hpdmatrixcholesky(CMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
void spdmatrixcholeskyupdateadd1buf(RMatrix *a, ae_int_t n, bool isupper, RVector *u, RVector *bufr, ae_state *_state);
void spdmatrixcholeskyupdateadd1(RMatrix *a, ae_int_t n, bool isupper, RVector *u, ae_state *_state);
void spdmatrixcholeskyupdatefixbuf(RMatrix *a, ae_int_t n, bool isupper, BVector *fix, RVector *bufr, ae_state *_state);
void spdmatrixcholeskyupdatefix(RMatrix *a, ae_int_t n, bool isupper, BVector *fix, ae_state *_state);
bool sparselu(sparsematrix *a, ae_int_t pivottype, ZVector *p, ZVector *q, ae_state *_state);
bool sparsecholeskyskyline(sparsematrix *a, ae_int_t n, bool isupper, ae_state *_state);
bool sparsecholesky(sparsematrix *a, bool isupper, ae_state *_state);
bool sparsecholeskyp(sparsematrix *a, bool isupper, ZVector *p, ae_state *_state);
bool sparsecholeskyanalyze(sparsematrix *a, bool isupper, ae_int_t facttype, ae_int_t permtype, sparsedecompositionanalysis *analysis, ae_state *_state);
void sparsecholeskysetmodtype(sparsedecompositionanalysis *analysis, ae_int_t modstrategy, double p0, double p1, double p2, double p3, ae_state *_state);
bool sparsecholeskyfactorize(sparsedecompositionanalysis *analysis, bool needupper, sparsematrix *a, RVector *d, ZVector *p, ae_state *_state);
void sparsecholeskyreload(sparsedecompositionanalysis *analysis, sparsematrix *a, bool isupper, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(sparsedecompositionanalysis, );

void rmatrixlu(real_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots, const xparams _xparams = NonTH);
void cmatrixlu(complex_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots, const xparams _xparams = NonTH);
bool spdmatrixcholesky(real_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
bool hpdmatrixcholesky(complex_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
void spdmatrixcholeskyupdateadd1buf(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u, real_1d_array &bufr, const xparams _xparams = NonTH);
void spdmatrixcholeskyupdateadd1(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u, const xparams _xparams = NonTH);
void spdmatrixcholeskyupdatefixbuf(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix, real_1d_array &bufr, const xparams _xparams = NonTH);
void spdmatrixcholeskyupdatefix(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix, const xparams _xparams = NonTH);
bool sparselu(const sparsematrix &a, const ae_int_t pivottype, integer_1d_array &p, integer_1d_array &q, const xparams _xparams = NonTH);
bool sparsecholeskyskyline(const sparsematrix &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
bool sparsecholesky(const sparsematrix &a, const bool isupper, const xparams _xparams = NonTH);
bool sparsecholeskyp(const sparsematrix &a, const bool isupper, integer_1d_array &p, const xparams _xparams = NonTH);
bool sparsecholeskyanalyze(const sparsematrix &a, const bool isupper, const ae_int_t facttype, const ae_int_t permtype, sparsedecompositionanalysis &analysis, const xparams _xparams = NonTH);
bool sparsecholeskyfactorize(const sparsedecompositionanalysis &analysis, const bool needupper, sparsematrix &a, real_1d_array &d, integer_1d_array &p, const xparams _xparams = NonTH);
void sparsecholeskyreload(const sparsedecompositionanalysis &analysis, const sparsematrix &a, const bool isupper, const xparams _xparams = NonTH);
} // end of namespace alglib

// === BDSVD Package ===
// Depends on: (AlgLibInternal) ROTATIONS
// Depends on: (AlgLibMisc) HQRND
// Depends on: ABLAS
namespace alglib_impl {
bool rmatrixbdsvd(RVector *d, RVector *e, ae_int_t n, bool isupper, bool isfractionalaccuracyrequired, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt, ae_state *_state);
bool bidiagonalsvddecomposition(RVector *d, RVector *e, ae_int_t n, bool isupper, bool isfractionalaccuracyrequired, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
bool rmatrixbdsvd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const bool isupper, const bool isfractionalaccuracyrequired, real_2d_array &u, const ae_int_t nru, real_2d_array &c, const ae_int_t ncc, real_2d_array &vt, const ae_int_t ncvt, const xparams _xparams = NonTH);
} // end of namespace alglib

// === SVD Package ===
// Depends on: (AlgLibInternal) BLAS
// Depends on: ORTFAC, BDSVD
namespace alglib_impl {
bool rmatrixsvd(RMatrix *a, ae_int_t m, ae_int_t n, ae_int_t uneeded, ae_int_t vtneeded, ae_int_t additionalmemory, RVector *w, RMatrix *u, RMatrix *vt, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
bool rmatrixsvd(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const ae_int_t uneeded, const ae_int_t vtneeded, const ae_int_t additionalmemory, real_1d_array &w, real_2d_array &u, real_2d_array &vt, const xparams _xparams = NonTH);
} // end of namespace alglib

// === RCOND Package ===
// Depends on: (AlgLibInternal) TRLINSOLVE, SAFESOLVE
// Depends on: TRFAC
namespace alglib_impl {
double rcondthreshold(ae_state *_state);
double rmatrixrcond1(RMatrix *a, ae_int_t n, ae_state *_state);
double cmatrixrcond1(CMatrix *a, ae_int_t n, ae_state *_state);
double rmatrixrcondinf(RMatrix *a, ae_int_t n, ae_state *_state);
double cmatrixrcondinf(CMatrix *a, ae_int_t n, ae_state *_state);
double spdmatrixrcond(RMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
double hpdmatrixrcond(CMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
double rmatrixtrrcond1(RMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_state *_state);
double cmatrixtrrcond1(CMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_state *_state);
double rmatrixtrrcondinf(RMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_state *_state);
double cmatrixtrrcondinf(CMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_state *_state);
double rmatrixlurcond1(RMatrix *lua, ae_int_t n, ae_state *_state);
double cmatrixlurcond1(CMatrix *lua, ae_int_t n, ae_state *_state);
double rmatrixlurcondinf(RMatrix *lua, ae_int_t n, ae_state *_state);
double cmatrixlurcondinf(CMatrix *lua, ae_int_t n, ae_state *_state);
double spdmatrixcholeskyrcond(RMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
double hpdmatrixcholeskyrcond(CMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double rmatrixrcond1(const real_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
double cmatrixrcond1(const complex_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
double rmatrixrcondinf(const real_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
double cmatrixrcondinf(const complex_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
double spdmatrixrcond(const real_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
double hpdmatrixrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
double rmatrixtrrcond1(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, const xparams _xparams = NonTH);
double cmatrixtrrcond1(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, const xparams _xparams = NonTH);
double rmatrixtrrcondinf(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, const xparams _xparams = NonTH);
double cmatrixtrrcondinf(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, const xparams _xparams = NonTH);
double rmatrixlurcond1(const real_2d_array &lua, const ae_int_t n, const xparams _xparams = NonTH);
double cmatrixlurcond1(const complex_2d_array &lua, const ae_int_t n, const xparams _xparams = NonTH);
double rmatrixlurcondinf(const real_2d_array &lua, const ae_int_t n, const xparams _xparams = NonTH);
double cmatrixlurcondinf(const complex_2d_array &lua, const ae_int_t n, const xparams _xparams = NonTH);
double spdmatrixcholeskyrcond(const real_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
double hpdmatrixcholeskyrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
} // end of namespace alglib

// === FBLS Package ===
// Depends on: (AlgLibInternal) ROTATIONS
// Depends on: ORTFAC
namespace alglib_impl {
struct fblslincgstate {
   double e1;
   double e2;
   ae_vector x;
   ae_vector ax;
   double xax;
   ae_int_t n;
   ae_vector rk;
   ae_vector rk1;
   ae_vector xk;
   ae_vector xk1;
   ae_vector pk;
   ae_vector pk1;
   ae_vector b;
   rcommstate rstate;
   ae_vector tmp2;
};
void fblslincgstate_init(void *_p, ae_state *_state, bool make_automatic);
void fblslincgstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void fblslincgstate_free(void *_p, bool make_automatic);

struct fblsgmresstate {
   ae_vector b;
   ae_vector x;
   ae_vector ax;
   ae_vector xs;
   ae_matrix qi;
   ae_matrix aqi;
   ae_matrix h;
   ae_matrix hq;
   ae_matrix hr;
   ae_vector hqb;
   ae_vector ys;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_int_t n;
   ae_int_t itscnt;
   double epsort;
   double epsres;
   double epsred;
   double epsdiag;
   ae_int_t itsperformed;
   ae_int_t retcode;
   rcommstate rstate;
};
void fblsgmresstate_init(void *_p, ae_state *_state, bool make_automatic);
void fblsgmresstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void fblsgmresstate_free(void *_p, bool make_automatic);

void fblscholeskysolve(RMatrix *cha, double sqrtscalea, ae_int_t n, bool isupper, RVector *xb, RVector *tmp, ae_state *_state);
void fblssolvecgx(RMatrix *a, ae_int_t m, ae_int_t n, double alpha, RVector *b, RVector *x, RVector *buf, ae_state *_state);
void fblscgcreate(RVector *x, RVector *b, ae_int_t n, fblslincgstate *state, ae_state *_state);
bool fblscgiteration(fblslincgstate *state, ae_state *_state);
void fblsgmrescreate(RVector *b, ae_int_t n, ae_int_t k, fblsgmresstate *state, ae_state *_state);
bool fblsgmresiteration(fblsgmresstate *state, ae_state *_state);
void fblssolvels(RMatrix *a, RVector *b, ae_int_t m, ae_int_t n, RVector *tmp0, RVector *tmp1, RVector *tmp2, ae_state *_state);
} // end of namespace alglib_impl

// === NORMESTIMATOR Package ===
// Depends on: MATGEN, SPARSE
namespace alglib_impl {
struct normestimatorstate {
   ae_int_t n;
   ae_int_t m;
   ae_int_t nstart;
   ae_int_t nits;
   ae_int_t seedval;
   ae_vector x0;
   ae_vector x1;
   ae_vector t;
   ae_vector xbest;
   hqrndstate r;
   ae_vector x;
   ae_vector mv;
   ae_vector mtv;
   bool needmv;
   bool needmtv;
   double repnorm;
   rcommstate rstate;
};
void normestimatorstate_init(void *_p, ae_state *_state, bool make_automatic);
void normestimatorstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void normestimatorstate_free(void *_p, bool make_automatic);

void normestimatorcreate(ae_int_t m, ae_int_t n, ae_int_t nstart, ae_int_t nits, normestimatorstate *state, ae_state *_state);
void normestimatorsetseed(normestimatorstate *state, ae_int_t seedval, ae_state *_state);
bool normestimatoriteration(normestimatorstate *state, ae_state *_state);
void normestimatorrestart(normestimatorstate *state, ae_state *_state);
void normestimatorestimatesparse(normestimatorstate *state, sparsematrix *a, ae_state *_state);
void normestimatorresults(normestimatorstate *state, double *nrm, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(normestimatorstate, );

void normestimatorcreate(const ae_int_t m, const ae_int_t n, const ae_int_t nstart, const ae_int_t nits, normestimatorstate &state, const xparams _xparams = NonTH);
void normestimatorsetseed(const normestimatorstate &state, const ae_int_t seedval, const xparams _xparams = NonTH);
void normestimatorestimatesparse(const normestimatorstate &state, const sparsematrix &a, const xparams _xparams = NonTH);
void normestimatorresults(const normestimatorstate &state, double &nrm, const xparams _xparams = NonTH);
} // end of namespace alglib

// === MATINV Package ===
// Depends on: RCOND
namespace alglib_impl {
struct matinvreport {
   double r1;
   double rinf;
};
void matinvreport_init(void *_p, ae_state *_state, bool make_automatic);
void matinvreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void matinvreport_free(void *_p, bool make_automatic);

void spdmatrixcholeskyinverserec(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, RVector *tmp, ae_state *_state);
void rmatrixluinverse(RMatrix *a, ZVector *pivots, ae_int_t n, ae_int_t *info, matinvreport *rep, ae_state *_state);
void cmatrixluinverse(CMatrix *a, ZVector *pivots, ae_int_t n, ae_int_t *info, matinvreport *rep, ae_state *_state);
void rmatrixinverse(RMatrix *a, ae_int_t n, ae_int_t *info, matinvreport *rep, ae_state *_state);
void cmatrixinverse(CMatrix *a, ae_int_t n, ae_int_t *info, matinvreport *rep, ae_state *_state);
void spdmatrixcholeskyinverse(RMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep, ae_state *_state);
void hpdmatrixcholeskyinverse(CMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep, ae_state *_state);
void spdmatrixinverse(RMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep, ae_state *_state);
void hpdmatrixinverse(CMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep, ae_state *_state);
void rmatrixtrinverse(RMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_int_t *info, matinvreport *rep, ae_state *_state);
void cmatrixtrinverse(CMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_int_t *info, matinvreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(matinvreport, double &r1; double &rinf;);

void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void rmatrixinverse(real_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void rmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void cmatrixinverse(complex_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void cmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void spdmatrixcholeskyinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void spdmatrixcholeskyinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void hpdmatrixcholeskyinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void hpdmatrixcholeskyinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void spdmatrixinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void spdmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void hpdmatrixinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void hpdmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void rmatrixtrinverse(real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void rmatrixtrinverse(real_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void cmatrixtrinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
void cmatrixtrinverse(complex_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === INVERSEUPDATE Package ===
namespace alglib_impl {
void rmatrixinvupdatesimple(RMatrix *inva, ae_int_t n, ae_int_t updrow, ae_int_t updcolumn, double updval, ae_state *_state);
void rmatrixinvupdaterow(RMatrix *inva, ae_int_t n, ae_int_t updrow, RVector *v, ae_state *_state);
void rmatrixinvupdatecolumn(RMatrix *inva, ae_int_t n, ae_int_t updcolumn, RVector *u, ae_state *_state);
void rmatrixinvupdateuv(RMatrix *inva, ae_int_t n, RVector *u, RVector *v, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void rmatrixinvupdatesimple(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const ae_int_t updcolumn, const double updval, const xparams _xparams = NonTH);
void rmatrixinvupdaterow(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const real_1d_array &v, const xparams _xparams = NonTH);
void rmatrixinvupdatecolumn(real_2d_array &inva, const ae_int_t n, const ae_int_t updcolumn, const real_1d_array &u, const xparams _xparams = NonTH);
void rmatrixinvupdateuv(real_2d_array &inva, const ae_int_t n, const real_1d_array &u, const real_1d_array &v, const xparams _xparams = NonTH);
} // end of namespace alglib

// === SCHUR Package ===
// Depends on: ORTFAC, HSSCHUR
namespace alglib_impl {
bool rmatrixschur(RMatrix *a, ae_int_t n, RMatrix *s, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
bool rmatrixschur(real_2d_array &a, const ae_int_t n, real_2d_array &s, const xparams _xparams = NonTH);
} // end of namespace alglib

// === SPDGEVD Package ===
// Depends on: EVD, MATINV
namespace alglib_impl {
bool smatrixgevdreduce(RMatrix *a, ae_int_t n, bool isuppera, RMatrix *b, bool isupperb, ae_int_t problemtype, RMatrix *r, bool *isupperr, ae_state *_state);
bool smatrixgevd(RMatrix *a, ae_int_t n, bool isuppera, RMatrix *b, bool isupperb, ae_int_t zneeded, ae_int_t problemtype, RVector *d, RMatrix *z, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
bool smatrixgevdreduce(real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t problemtype, real_2d_array &r, bool &isupperr, const xparams _xparams = NonTH);
bool smatrixgevd(const real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t zneeded, const ae_int_t problemtype, real_1d_array &d, real_2d_array &z, const xparams _xparams = NonTH);
} // end of namespace alglib

// === MATDET Package ===
// Depends on: TRFAC
namespace alglib_impl {
double rmatrixludet(RMatrix *a, ZVector *pivots, ae_int_t n, ae_state *_state);
complex cmatrixludet(CMatrix *a, ZVector *pivots, ae_int_t n, ae_state *_state);
double rmatrixdet(RMatrix *a, ae_int_t n, ae_state *_state);
complex cmatrixdet(CMatrix *a, ae_int_t n, ae_state *_state);
double spdmatrixcholeskydet(RMatrix *a, ae_int_t n, ae_state *_state);
double spdmatrixdet(RMatrix *a, ae_int_t n, bool isupper, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, const xparams _xparams = NonTH);
double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots, const xparams _xparams = NonTH);
complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, const xparams _xparams = NonTH);
complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots, const xparams _xparams = NonTH);
double rmatrixdet(const real_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
double rmatrixdet(const real_2d_array &a, const xparams _xparams = NonTH);
complex cmatrixdet(const complex_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
complex cmatrixdet(const complex_2d_array &a, const xparams _xparams = NonTH);
double spdmatrixcholeskydet(const real_2d_array &a, const ae_int_t n, const xparams _xparams = NonTH);
double spdmatrixcholeskydet(const real_2d_array &a, const xparams _xparams = NonTH);
double spdmatrixdet(const real_2d_array &a, const ae_int_t n, const bool isupper, const xparams _xparams = NonTH);
double spdmatrixdet(const real_2d_array &a, const xparams _xparams = NonTH);
} // end of namespace alglib

#endif // OnceOnly
