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
#ifndef OnceOnlyAlgLibInternal_h
#define OnceOnlyAlgLibInternal_h

#include "Ap.h"

// === APSERV Package ===
namespace alglib_impl {
struct apbuffers {
   ae_vector ba0;
   ae_vector ia0;
   ae_vector ia1;
   ae_vector ia2;
   ae_vector ia3;
   ae_vector ra0;
   ae_vector ra1;
   ae_vector ra2;
   ae_vector ra3;
   ae_matrix rm0;
   ae_matrix rm1;
};
void apbuffers_init(void *_p, bool make_automatic);
void apbuffers_copy(void *_dst, void *_src, bool make_automatic);
void apbuffers_free(void *_p, bool make_automatic);

// typedef ae_vector RVector; //(@) In Ap.h.
void RVector_init(void *_p, bool make_automatic);
void RVector_copy(void *_dst, void *_src, bool make_automatic);
void RVector_free(void *_p, bool make_automatic);

double coalesce(double a, double b);
ae_int_t coalescei(ae_int_t a, ae_int_t b);
double logbase2(double x);

// Test number A for equality to B (or to 0.0) within the range Tiny; inclusive for *At*().
// Tiny > 0; and for the *At*() cases, Tiny == 0.0 is also allowed.
//(@) Originally: approxequal(double A, double B, double Tiny);
//(@) AlgLib: Copyright 02.12.2009 by Sergey Bochkanov
static inline bool NearAtR(double A, double B, double Tiny) { return fabs(A - B) <= Tiny; }
static inline bool NearAtCR(complex A, double B, double Tiny) { return abscomplex(ae_c_sub_d(A, B)) <= Tiny; }
static inline bool NearAtC(complex A, complex B, double Tiny) { return abscomplex(ae_c_sub(A, B)) <= Tiny; }
static inline bool NearR(double A, double B, double Tiny) { return fabs(A - B) < Tiny; }
static inline bool NearC(complex A, complex B, double Tiny) { return abscomplex(ae_c_sub(A, B)) < Tiny; }
static inline bool SmallAtR(double A, double Tiny) { return fabs(A) <= Tiny; }
static inline bool SmallAtC(complex A, double Tiny) { return abscomplex(A) <= Tiny; }
static inline bool SmallR(double A, double Tiny) { return fabs(A) < Tiny; }
static inline bool SmallC(complex A, double Tiny) { return abscomplex(A) < Tiny; }

void taskgenint1d(double a, double b, ae_int_t n, RVector *x, RVector *y);
void taskgenint1dequidist(double a, double b, ae_int_t n, RVector *x, RVector *y);
void taskgenint1dcheb1(double a, double b, ae_int_t n, RVector *x, RVector *y);
void taskgenint1dcheb2(double a, double b, ae_int_t n, RVector *x, RVector *y);
bool aredistinct(RVector *x, ae_int_t n);
void setlengthzero(RVector *x, ae_int_t n);
void vectorsetlengthatleast(ae_vector *x, ae_int_t n);
void matrixsetlengthatleast(ae_matrix *x, ae_int_t m, ae_int_t n);
void bvectorgrowto(BVector *x, ae_int_t n);
void ivectorgrowto(ZVector *x, ae_int_t n);
void rvectorgrowto(RVector *x, ae_int_t n);
void rmatrixgrowrowsto(RMatrix *a, ae_int_t n, ae_int_t mincols);
void rmatrixgrowcolsto(RMatrix *a, ae_int_t n, ae_int_t minrows);
void ivectorresize(ZVector *x, ae_int_t n);
void rvectorresize(RVector *x, ae_int_t n);
void imatrixresize(ZMatrix *x, ae_int_t m, ae_int_t n);
void rmatrixresize(RMatrix *x, ae_int_t m, ae_int_t n);
bool isfinitevector(RVector *x, ae_int_t n);
bool isfinitecvector(CVector *z, ae_int_t n);
bool apservisfinitematrix(RMatrix *x, ae_int_t m, ae_int_t n);
bool apservisfinitecmatrix(CMatrix *x, ae_int_t m, ae_int_t n);
bool isfinitertrmatrix(RMatrix *x, ae_int_t n, bool isupper);
bool apservisfinitectrmatrix(CMatrix *x, ae_int_t n, bool isupper);
bool apservisfiniteornanmatrix(RMatrix *x, ae_int_t m, ae_int_t n);
double safepythag2(double x, double y);
double safepythag3(double x, double y, double z);
ae_int_t saferdiv(double x, double y, double *r);
double safeminposrv(double x, double y, double v);
void apperiodicmap(double *x, double a, double b, double *k);
double randomnormal();
void randomunit(ae_int_t n, RVector *x);
void swapb(bool *v0, bool *v1);
void swapi(ae_int_t *v0, ae_int_t *v1);
void swapr(double *v0, double *v1);
void swapc(complex *v0, complex *v1);
void swapcols(RMatrix *a, ae_int_t j0, ae_int_t j1, ae_int_t nrows);
void swaprows(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t ncols);
void swapentries(RVector *a, ae_int_t i0, ae_int_t i1, ae_int_t entrywidth);
void swapelementsi(ZVector *a, ae_int_t i0, ae_int_t i1);
void swapelements(RVector *a, ae_int_t i0, ae_int_t i1);
double possign(double x);
void alloccomplex(ae_serializer *s, complex v);
void serializecomplex(ae_serializer *s, complex v);
complex unserializecomplex(ae_serializer *s);
void allocrealarray(ae_serializer *s, RVector *v, ae_int_t n);
void serializerealarray(ae_serializer *s, RVector *v, ae_int_t n);
void unserializerealarray(ae_serializer *s, RVector *v);
void allocintegerarray(ae_serializer *s, ZVector *v, ae_int_t n);
void serializeintegerarray(ae_serializer *s, ZVector *v, ae_int_t n);
void unserializeintegerarray(ae_serializer *s, ZVector *v);
void allocrealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1);
void serializerealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1);
void unserializerealmatrix(ae_serializer *s, RMatrix *v);
void copybooleanarray(BVector *src, BVector *dst);
void copyintegerarray(ZVector *src, ZVector *dst);
void copyrealarray(RVector *src, RVector *dst);
void copyrealmatrix(RMatrix *src, RMatrix *dst);
ae_int_t chunkscount(ae_int_t tasksize, ae_int_t chunksize);
ae_int_t tiledsplit(ae_int_t tasksize, ae_int_t tilesize);
ae_int_t splitlength(ae_int_t tasksize, ae_int_t chunksize);
ae_int_t splitlengtheven(ae_int_t tasksize);
ae_int_t recsearch(ZVector *a, ae_int_t nrec, ae_int_t nheader, ae_int_t i0, ae_int_t i1, ZVector *b);
double sparselevel2density();
ae_int_t matrixtilesizea();
ae_int_t matrixtilesizeb();
double smpactivationlevel();
double spawnlevel();
} // end of namespace alglib_impl

// === ABLASF Package ===
namespace alglib_impl {
double rdotv(ae_int_t n, RVector *x, RVector *y);
double rdotvr(ae_int_t n, RVector *x, RMatrix *a, ae_int_t i);
double rdotrr(ae_int_t n, RMatrix *a, ae_int_t ia, RMatrix *b, ae_int_t ib);
double rdotv2(ae_int_t n, RVector *x);
void raddv(ae_int_t n, double alpha, RVector *y, RVector *x);
void raddvx(ae_int_t n, double alpha, RVector *y, ae_int_t offsy, RVector *x, ae_int_t offsx);
void raddvc(ae_int_t n, double alpha, RVector *y, RMatrix *x, ae_int_t colidx);
void raddvr(ae_int_t n, double alpha, RVector *y, RMatrix *x, ae_int_t rowidx);
void rmergemulvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx);
void rmergemaxrv(ae_int_t n, RMatrix *x, ae_int_t rowidx, RVector *y);
void raddrv(ae_int_t n, double alpha, RMatrix *y, ae_int_t ridx, RVector *x);
void rmulv(ae_int_t n, double v, RVector *x);
void rmulr(ae_int_t n, double v, RMatrix *x, ae_int_t rowidx);
void rmulvx(ae_int_t n, double v, RVector *x, ae_int_t offsx);
double rmaxabsv(ae_int_t n, RVector *x);
void bsetv(ae_int_t n, bool v, BVector *x);
void isetv(ae_int_t n, ae_int_t v, ZVector *x);
void rsetv(ae_int_t n, double v, RVector *x);
void rsetvx(ae_int_t n, double v, RVector *x, ae_int_t offsx);
void rsetm(ae_int_t m, ae_int_t n, double v, RMatrix *a);
void bsetallocv(ae_int_t n, bool v, BVector *x);
void isetallocv(ae_int_t n, ae_int_t v, ZVector *x);
void rsetallocv(ae_int_t n, double v, RVector *x);
void rsetallocm(ae_int_t m, ae_int_t n, double v, RMatrix *a);
void ballocv(ae_int_t n, BVector *x);
void rallocv(ae_int_t n, RVector *x);
void rallocm(ae_int_t m, ae_int_t n, RMatrix *a);
void rsetr(ae_int_t n, double v, RMatrix *a, ae_int_t i);
void rsetc(ae_int_t n, double v, RMatrix *a, ae_int_t j);
void bcopyv(ae_int_t n, BVector *x, BVector *y);
void icopyv(ae_int_t n, ZVector *x, ZVector *y);
void rcopyv(ae_int_t n, RVector *x, RVector *y);
void icopyvx(ae_int_t n, ZVector *x, ae_int_t offsx, ZVector *y, ae_int_t offsy);
void rcopyvx(ae_int_t n, RVector *x, ae_int_t offsx, RVector *y, ae_int_t offsy);
void bcopyallocv(ae_int_t n, BVector *x, BVector *y);
void icopyallocv(ae_int_t n, ZVector *x, ZVector *y);
void rcopyallocv(ae_int_t n, RVector *x, RVector *y);
void rcopyallocm(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y);
void igrowv(ae_int_t newn, ZVector *x);
void rcopymulv(ae_int_t n, double v, RVector *x, RVector *y);
void rcopymulvr(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t ridx);
void rcopyvr(ae_int_t n, RVector *x, RMatrix *a, ae_int_t i);
void rcopyrv(ae_int_t n, RMatrix *a, ae_int_t i, RVector *x);
void rcopyrr(ae_int_t n, RMatrix *a, ae_int_t i, RMatrix *b, ae_int_t k);
void rcopyvc(ae_int_t n, RVector *x, RMatrix *a, ae_int_t j);
void rcopycv(ae_int_t n, RMatrix *a, ae_int_t j, RVector *x);
bool rmatrixgerf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double ralpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv);
bool rmatrixrank1f(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv);
bool cmatrixrank1f(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv);
bool rmatrixlefttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2);
bool cmatrixlefttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2);
bool rmatrixrighttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2);
bool cmatrixrighttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2);
bool rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
bool cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc);
void rmatrixgemmk44v00(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
void rmatrixgemmk44v01(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
void rmatrixgemmk44v10(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
void rmatrixgemmk44v11(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
void rmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
void cmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc);
} // end of namespace alglib_impl

// === HBLAS Package ===
namespace alglib_impl {
void hermitianmatrixvectormultiply(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, complex alpha, CVector *y);
void hermitianrank2update(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, CVector *y, CVector *t, complex alpha);
} // end of namespace alglib_impl

// === CREFLECTIONS Package ===
namespace alglib_impl {
void complexgeneratereflection(CVector *x, ae_int_t n, complex *tau);
void complexapplyreflectionfromtheleft(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work);
void complexapplyreflectionfromtheright(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work);
} // end of namespace alglib_impl

// === SBLAS Package ===
// Depends on: APSERV
namespace alglib_impl {
void symmetricmatrixvectormultiply(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, double alpha, RVector *y);
void symmetricrank2update(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, RVector *y, RVector *t, double alpha);
} // end of namespace alglib_impl

// === ABLASMKL Package ===
namespace alglib_impl {
bool rmatrixgermkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double alpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv);
bool rmatrixrank1mkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv);
bool cmatrixrank1mkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv);
bool rmatrixmvmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, RVector *y, ae_int_t iy);
bool cmatrixmvmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CVector *x, ae_int_t ix, CVector *y, ae_int_t iy);
bool rmatrixgemvmkl(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy);
bool rmatrixtrsvmkl(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix);
bool rmatrixsymvmkl(ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy);
bool rmatrixsyrkmkl(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool cmatrixherkmkl(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool rmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
bool cmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc);
bool rmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2);
bool cmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2);
bool rmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2);
bool cmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2);
bool spdmatrixcholeskymkl(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, bool *cholresult);
bool rmatrixplumkl(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots);
bool rmatrixbdmkl(RMatrix *a, ae_int_t m, ae_int_t n, RVector *d, RVector *e, RVector *tauq, RVector *taup);
bool rmatrixbdmultiplybymkl(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, RVector *taup, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool byq, bool fromtheright, bool dotranspose);
bool rmatrixhessenbergmkl(RMatrix *a, ae_int_t n, RVector *tau);
bool rmatrixhessenbergunpackqmkl(RMatrix *a, ae_int_t n, RVector *tau, RMatrix *q);
bool smatrixtdmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RVector *d, RVector *e);
bool hmatrixtdmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, RVector *d, RVector *e);
bool smatrixtdunpackqmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RMatrix *q);
bool hmatrixtdunpackqmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, CMatrix *q);
bool rmatrixbdsvdmkl(RVector *d, RVector *e, ae_int_t n, bool isupper, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt, bool *svdresult);
bool rmatrixinternalschurdecompositionmkl(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info);
bool rmatrixinternaltrevcmkl(RMatrix *t, ae_int_t n, ae_int_t side, ae_int_t howmny, RMatrix *vl, RMatrix *vr, ae_int_t *m, ae_int_t *info);
bool smatrixtdevdmkl(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z, bool *evdresult);
bool sparsegemvcrsmkl(ae_int_t opa, ae_int_t arows, ae_int_t acols, double alpha, RVector *vals, ZVector *cidx, ZVector *ridx, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy);
} // end of namespace alglib_impl

// === SCODES Package ===
namespace alglib_impl {
ae_int_t getrdfserializationcode();
ae_int_t getkdtreeserializationcode();
ae_int_t getmlpserializationcode();
ae_int_t getmlpeserializationcode();
ae_int_t getrbfserializationcode();
ae_int_t getspline2dserializationcode();
ae_int_t getidwserializationcode();
ae_int_t getknnserializationcode();
} // end of namespace alglib_impl

// === TSORT Package ===
// Depends on: APSERV
namespace alglib_impl {
void tagsortfasti(RVector *a, ZVector *b, RVector *bufa, ZVector *bufb, ae_int_t n);
void tagsortfastr(RVector *a, RVector *b, RVector *bufa, RVector *bufb, ae_int_t n);
void tagsortfast(RVector *a, RVector *bufa, ae_int_t n);
void tagsortmiddleir(ZVector *a, RVector *b, ae_int_t n, ae_int_t offset = 0);
void tagsortmiddlei(ZVector *a, ae_int_t n, ae_int_t offset = 0);
void tagsortbuf(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2, apbuffers *buf);
void tagsort(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2);
void tagheappushi(RVector *a, ZVector *b, ae_int_t *n, double va, ae_int_t vb);
void tagheapreplacetopi(RVector *a, ZVector *b, ae_int_t n, double va, ae_int_t vb);
void tagheappopi(RVector *a, ZVector *b, ae_int_t *n);
ae_int_t lowerbound(RVector *a, ae_int_t n, double t);
ae_int_t upperbound(RVector *a, ae_int_t n, double t);
} // end of namespace alglib_impl

// === BLAS Package ===
namespace alglib_impl {
double vectornorm2(RVector *x, ae_int_t i1, ae_int_t i2);
ae_int_t vectoridxabsmax(RVector *x, ae_int_t i1, ae_int_t i2);
ae_int_t columnidxabsmax(RMatrix *x, ae_int_t i1, ae_int_t i2, ae_int_t j);
ae_int_t rowidxabsmax(RMatrix *x, ae_int_t j1, ae_int_t j2, ae_int_t i);
double upperhessenberg1norm(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work);
void copymatrix(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2);
void inplacetranspose(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work);
void copyandtranspose(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2);
void matrixvectormultiply(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, bool trans, RVector *x, ae_int_t ix1, ae_int_t ix2, double alpha, RVector *y, ae_int_t iy1, ae_int_t iy2, double beta);
void matrixmatrixmultiply(RMatrix *a, ae_int_t ai1, ae_int_t ai2, ae_int_t aj1, ae_int_t aj2, bool transa, RMatrix *b, ae_int_t bi1, ae_int_t bi2, ae_int_t bj1, ae_int_t bj2, bool transb, double alpha, RMatrix *c, ae_int_t ci1, ae_int_t ci2, ae_int_t cj1, ae_int_t cj2, double beta, RVector *work);
} // end of namespace alglib_impl

// === ROTATIONS Package ===
namespace alglib_impl {
void applyrotationsfromtheleft(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work);
void applyrotationsfromtheright(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work);
void generaterotation(double f, double g, double *cs, double *sn, double *r);
} // end of namespace alglib_impl

// === BASICSTATOPS Package ===
// Depends on: TSORT
namespace alglib_impl {
void rankx(RVector *x, ae_int_t n, bool iscentered, apbuffers *buf);
void rankxuntied(RVector *x, ae_int_t n, apbuffers *buf);
} // end of namespace alglib_impl

// === TRLINSOLVE Package ===
namespace alglib_impl {
void safesolvetriangular(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit, bool normin, RVector *cnorm);
void rmatrixtrsafesolve(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit);
} // end of namespace alglib_impl

// === SAFESOLVE Package ===
namespace alglib_impl {
bool rmatrixscaledtrsafesolve(RMatrix *a, double sa, ae_int_t n, RVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth);
bool cmatrixscaledtrsafesolve(CMatrix *a, double sa, ae_int_t n, CVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth);
} // end of namespace alglib_impl

// === XBLAS Package ===
namespace alglib_impl {
void xdot(RVector *a, RVector *b, ae_int_t n, RVector *temp, double *r, double *rerr);
void xcdot(CVector *a, CVector *b, ae_int_t n, RVector *temp, complex *r, double *rerr);
} // end of namespace alglib_impl

// === LINMIN Package ===
namespace alglib_impl {
struct linminstate {
   bool brackt;
   bool stage1;
   ae_int_t infoc;
   double dg;
   double dgm;
   double dginit;
   double dgtest;
   double dgx;
   double dgxm;
   double dgy;
   double dgym;
   double finit;
   double ftest1;
   double fm;
   double fx;
   double fxm;
   double fy;
   double fym;
   double stx;
   double sty;
   double stmin;
   double stmax;
   double width;
   double width1;
   double xtrapf;
};
void linminstate_init(void *_p, bool make_automatic);
void linminstate_copy(void *_dst, void *_src, bool make_automatic);
void linminstate_free(void *_p, bool make_automatic);

struct armijostate {
   bool needf;
   ae_vector x;
   double f;
   ae_int_t n;
   ae_vector xbase;
   ae_vector s;
   double stplen;
   double fcur;
   double stpmax;
   ae_int_t fmax;
   ae_int_t nfev;
   ae_int_t info;
   ae_int_t PQ;
};
void armijostate_init(void *_p, bool make_automatic);
void armijostate_copy(void *_dst, void *_src, bool make_automatic);
void armijostate_free(void *_p, bool make_automatic);

void linminnormalized(RVector *d, double *stp, ae_int_t n);
void mcstep(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, double fp, double dp, bool *brackt, double stmin, double stmax, ae_int_t *info);
bool mcsrch(ae_int_t n, RVector *x, double f, RVector *g, RVector *s, double *stp, double stpmax, double gtol, ae_int_t *info, ae_int_t *nfev, RVector *wa, linminstate *state, ae_int_t *stage);
void armijocreate(ae_int_t n, RVector *x, double f, RVector *s, double stp, double stpmax, ae_int_t fmax, armijostate *state);
bool armijoiteration(armijostate *state);
void armijoresults(armijostate *state, ae_int_t *info, double *stp, double *f);
} // end of namespace alglib_impl

// === NEARUNITYUNIT Package ===
namespace alglib_impl {
double nulog1p(double x);
double nuexpm1(double x);
double nucosm1(double x); //(@) Not used or tested anywhere.
} // end of namespace alglib_impl

// === NTHEORY Package ===
namespace alglib_impl {
void findprimitiverootandinverse(ae_int_t n, ae_int_t *proot, ae_int_t *invproot);
} // end of namespace alglib_impl

// === FTBASE Package ===
// Depends on: APSERV, NTHEORY
namespace alglib_impl {
struct fasttransformplan {
   ae_matrix entries;
   ae_vector buffer;
   ae_vector precr;
   ae_vector preci;
   ae_shared_pool bluesteinpool;
};
void fasttransformplan_init(void *_p, bool make_automatic);
void fasttransformplan_copy(void *_dst, void *_src, bool make_automatic);
void fasttransformplan_free(void *_p, bool make_automatic);

double ftbasegetflopestimate(ae_int_t n);
void ftbasefactorize(ae_int_t n, ae_int_t tasktype, ae_int_t *n1, ae_int_t *n2);
bool ftbaseissmooth(ae_int_t n);
ae_int_t ftbasefindsmooth(ae_int_t n);
ae_int_t ftbasefindsmootheven(ae_int_t n);
void ftapplyplan(fasttransformplan *plan, RVector *a, ae_int_t offsa = 0, ae_int_t repcnt = 1);
void ftcomplexfftplan(ae_int_t n, ae_int_t k, fasttransformplan *plan);
} // end of namespace alglib_impl

// === HPCCORES Package ===
namespace alglib_impl {
struct mlpbuffers {
   ae_int_t chunksize;
   ae_int_t ntotal;
   ae_int_t nin;
   ae_int_t nout;
   ae_int_t wcount;
   ae_vector batch4buf;
   ae_vector hpcbuf;
   ae_matrix xy;
   ae_matrix xy2;
   ae_vector xyrow;
   ae_vector x;
   ae_vector y;
   ae_vector desiredy;
   double e;
   ae_vector g;
   ae_vector tmp0;
};
void mlpbuffers_init(void *_p, bool make_automatic);
void mlpbuffers_copy(void *_dst, void *_src, bool make_automatic);
void mlpbuffers_free(void *_p, bool make_automatic);

void hpcpreparechunkedgradient(RVector *weights, ae_int_t wcount, ae_int_t ntotal, ae_int_t nin, ae_int_t nout, mlpbuffers *buf);
void hpcfinalizechunkedgradient(mlpbuffers *buf, RVector *grad);
bool hpcchunkedgradient(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf, double *e, bool naturalerrorfunc);
bool hpcchunkedprocess(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf);
} // end of namespace alglib_impl

#endif // OnceOnly
