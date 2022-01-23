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
void apbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void apbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void apbuffers_free(void *_p, bool make_automatic);

struct sboolean {
   bool val;
};
void sboolean_init(void *_p, ae_state *_state, bool make_automatic);
void sboolean_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sboolean_free(void *_p, bool make_automatic);

struct sbooleanarray {
   ae_vector val;
};
void sbooleanarray_init(void *_p, ae_state *_state, bool make_automatic);
void sbooleanarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sbooleanarray_free(void *_p, bool make_automatic);

struct sinteger {
   ae_int_t val;
};
void sinteger_init(void *_p, ae_state *_state, bool make_automatic);
void sinteger_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sinteger_free(void *_p, bool make_automatic);

struct sintegerarray {
   ae_vector val;
};
void sintegerarray_init(void *_p, ae_state *_state, bool make_automatic);
void sintegerarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sintegerarray_free(void *_p, bool make_automatic);

struct sreal {
   double val;
};
void sreal_init(void *_p, ae_state *_state, bool make_automatic);
void sreal_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sreal_free(void *_p, bool make_automatic);

struct srealarray {
   ae_vector val;
};
void srealarray_init(void *_p, ae_state *_state, bool make_automatic);
void srealarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void srealarray_free(void *_p, bool make_automatic);

struct scomplex {
   complex val;
};
void scomplex_init(void *_p, ae_state *_state, bool make_automatic);
void scomplex_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void scomplex_free(void *_p, bool make_automatic);

struct scomplexarray {
   ae_vector val;
};
void scomplexarray_init(void *_p, ae_state *_state, bool make_automatic);
void scomplexarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void scomplexarray_free(void *_p, bool make_automatic);

void seterrorflagdiff(bool *flag, double val, double refval, double tol, double s, ae_state *_state);
bool alwaysfalse(ae_state *_state);
void touchint(ae_int_t *a, ae_state *_state);
void touchreal(double *a, ae_state *_state);
double coalesce(double a, double b, ae_state *_state);
ae_int_t coalescei(ae_int_t a, ae_int_t b, ae_state *_state);
double inttoreal(ae_int_t a, ae_state *_state);
double logbase2(double x, ae_state *_state);
bool approxequal(double a, double b, double tol, ae_state *_state);
bool approxequalrel(double a, double b, double tol, ae_state *_state);
void taskgenint1d(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state);
void taskgenint1dequidist(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state);
void taskgenint1dcheb1(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state);
void taskgenint1dcheb2(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state);
bool aredistinct(RVector *x, ae_int_t n, ae_state *_state);
bool aresameboolean(bool v1, bool v2, ae_state *_state);
void setlengthzero(RVector *x, ae_int_t n, ae_state *_state);
void bvectorsetlengthatleast(BVector *x, ae_int_t n, ae_state *_state);
void ivectorsetlengthatleast(ZVector *x, ae_int_t n, ae_state *_state);
void rvectorsetlengthatleast(RVector *x, ae_int_t n, ae_state *_state);
void bmatrixsetlengthatleast(BMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
void rmatrixsetlengthatleast(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
void bvectorgrowto(BVector *x, ae_int_t n, ae_state *_state);
void ivectorgrowto(ZVector *x, ae_int_t n, ae_state *_state);
void rvectorgrowto(RVector *x, ae_int_t n, ae_state *_state);
void rmatrixgrowrowsto(RMatrix *a, ae_int_t n, ae_int_t mincols, ae_state *_state);
void rmatrixgrowcolsto(RMatrix *a, ae_int_t n, ae_int_t minrows, ae_state *_state);
void ivectorresize(ZVector *x, ae_int_t n, ae_state *_state);
void rvectorresize(RVector *x, ae_int_t n, ae_state *_state);
void imatrixresize(ZMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
void rmatrixresize(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
void ivectorappend(ZVector *x, ae_int_t v, ae_state *_state);
bool isfinitevector(RVector *x, ae_int_t n, ae_state *_state);
bool isfinitecvector(CVector *z, ae_int_t n, ae_state *_state);
bool apservisfinitematrix(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
bool apservisfinitecmatrix(CMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
bool isfinitertrmatrix(RMatrix *x, ae_int_t n, bool isupper, ae_state *_state);
bool apservisfinitectrmatrix(CMatrix *x, ae_int_t n, bool isupper, ae_state *_state);
bool apservisfiniteornanmatrix(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state);
double safepythag2(double x, double y, ae_state *_state);
double safepythag3(double x, double y, double z, ae_state *_state);
ae_int_t saferdiv(double x, double y, double *r, ae_state *_state);
double safeminposrv(double x, double y, double v, ae_state *_state);
void apperiodicmap(double *x, double a, double b, double *k, ae_state *_state);
double randomnormal(ae_state *_state);
void randomunit(ae_int_t n, RVector *x, ae_state *_state);
void swapi(ae_int_t *v0, ae_int_t *v1, ae_state *_state);
void swapr(double *v0, double *v1, ae_state *_state);
void swapcols(RMatrix *a, ae_int_t j0, ae_int_t j1, ae_int_t nrows, ae_state *_state);
void swaprows(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t ncols, ae_state *_state);
void swapentries(RVector *a, ae_int_t i0, ae_int_t i1, ae_int_t entrywidth, ae_state *_state);
void swapelementsi(ZVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state);
void swapelements(RVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state);
void inc(ae_int_t *v, ae_state *_state);
void dec(ae_int_t *v, ae_state *_state);
void threadunsafeinc(ae_int_t *v, ae_state *_state);
void threadunsafeincby(ae_int_t *v, ae_int_t k, ae_state *_state);
void countdown(ae_int_t *v, ae_state *_state);
double possign(double x, ae_state *_state);
double rmul2(double v0, double v1, ae_state *_state);
double rmul3(double v0, double v1, double v2, ae_state *_state);
ae_int_t idivup(ae_int_t a, ae_int_t b, ae_state *_state);
ae_int_t imin2(ae_int_t i0, ae_int_t i1, ae_state *_state);
ae_int_t imin3(ae_int_t i0, ae_int_t i1, ae_int_t i2, ae_state *_state);
ae_int_t imax2(ae_int_t i0, ae_int_t i1, ae_state *_state);
ae_int_t imax3(ae_int_t i0, ae_int_t i1, ae_int_t i2, ae_state *_state);
double maxreal3(double v0, double v1, double v2, ae_state *_state);
double rmax3(double r0, double r1, double r2, ae_state *_state);
double rmaxabs3(double r0, double r1, double r2, ae_state *_state);
double boundval(double x, double b1, double b2, ae_state *_state);
ae_int_t iboundval(ae_int_t x, ae_int_t b1, ae_int_t b2, ae_state *_state);
double rboundval(double x, double b1, double b2, ae_state *_state);
ae_int_t countnz1(RVector *v, ae_int_t n, ae_state *_state);
ae_int_t countnz2(RMatrix *v, ae_int_t m, ae_int_t n, ae_state *_state);
void alloccomplex(ae_serializer *s, complex v, ae_state *_state);
void serializecomplex(ae_serializer *s, complex v, ae_state *_state);
complex unserializecomplex(ae_serializer *s, ae_state *_state);
void allocrealarray(ae_serializer *s, RVector *v, ae_int_t n, ae_state *_state);
void serializerealarray(ae_serializer *s, RVector *v, ae_int_t n, ae_state *_state);
void unserializerealarray(ae_serializer *s, RVector *v, ae_state *_state);
void allocintegerarray(ae_serializer *s, ZVector *v, ae_int_t n, ae_state *_state);
void serializeintegerarray(ae_serializer *s, ZVector *v, ae_int_t n, ae_state *_state);
void unserializeintegerarray(ae_serializer *s, ZVector *v, ae_state *_state);
void allocrealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1, ae_state *_state);
void serializerealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1, ae_state *_state);
void unserializerealmatrix(ae_serializer *s, RMatrix *v, ae_state *_state);
void copybooleanarray(BVector *src, BVector *dst, ae_state *_state);
void copyintegerarray(ZVector *src, ZVector *dst, ae_state *_state);
void copyrealarray(RVector *src, RVector *dst, ae_state *_state);
void copyrealmatrix(RMatrix *src, RMatrix *dst, ae_state *_state);
void unsetintegerarray(ZVector *a, ae_state *_state);
void unsetrealarray(RVector *a, ae_state *_state);
void unsetrealmatrix(RMatrix *a, ae_state *_state);
ae_int_t chunkscount(ae_int_t tasksize, ae_int_t chunksize, ae_state *_state);
void tiledsplit(ae_int_t tasksize, ae_int_t tilesize, ae_int_t *task0, ae_int_t *task1, ae_state *_state);
void splitlength(ae_int_t tasksize, ae_int_t chunksize, ae_int_t *task0, ae_int_t *task1, ae_state *_state);
void splitlengtheven(ae_int_t tasksize, ae_int_t *task0, ae_int_t *task1, ae_state *_state);
ae_int_t recsearch(ZVector *a, ae_int_t nrec, ae_int_t nheader, ae_int_t i0, ae_int_t i1, ZVector *b, ae_state *_state);
double sparselevel2density(ae_state *_state);
ae_int_t matrixtilesizea(ae_state *_state);
ae_int_t matrixtilesizeb(ae_state *_state);
double smpactivationlevel(ae_state *_state);
double spawnlevel(ae_state *_state);
void tracevectorautoprec(RVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state);
void tracerowautoprec(RMatrix *a, ae_int_t i, ae_int_t j0, ae_int_t j1, ae_state *_state);
void tracevectorunscaledunshiftedautoprec(RVector *x, ae_int_t n, RVector *scl, bool applyscl, RVector *sft, bool applysft, ae_state *_state);
void tracerownrm1autoprec(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t j0, ae_int_t j1, ae_state *_state);
void tracevectore6(RVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state);
void tracevectore615(RVector *a, ae_int_t i0, ae_int_t i1, bool usee15, ae_state *_state);
void tracerownrm1e6(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t j0, ae_int_t j1, ae_state *_state);
} // end of namespace alglib_impl

// === ABLASF Package ===
namespace alglib_impl {
#ifdef ALGLIB_NO_FAST_KERNELS
double rdotv(ae_int_t n, RVector *x, RVector *y, ae_state *_state);
double rdotvr(ae_int_t n, RVector *x, RMatrix *a, ae_int_t i, ae_state *_state);
double rdotrr(ae_int_t n, RMatrix *a, ae_int_t ia, RMatrix *b, ae_int_t ib, ae_state *_state);
double rdotv2(ae_int_t n, RVector *x, ae_state *_state);
void raddv(ae_int_t n, double alpha, RVector *y, RVector *x, ae_state *_state);
void raddvx(ae_int_t n, double alpha, RVector *y, ae_int_t offsy, RVector *x, ae_int_t offsx, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
void raddvc(ae_int_t n, double alpha, RVector *y, RMatrix *x, ae_int_t colidx, ae_state *_state);
#ifdef ALGLIB_NO_FAST_KERNELS
void raddvr(ae_int_t n, double alpha, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state);
void rmergemulv(ae_int_t n, RVector *y, RVector *x, ae_state *_state);
void rmergemulvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state);
void rmergemulrv(ae_int_t n, RMatrix *y, ae_int_t rowidx, RVector *x, ae_state *_state);
void rmergemaxv(ae_int_t n, RVector *y, RVector *x, ae_state *_state);
void rmergemaxvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state);
void rmergemaxrv(ae_int_t n, RMatrix *x, ae_int_t rowidx, RVector *y, ae_state *_state);
void rmergeminv(ae_int_t n, RVector *y, RVector *x, ae_state *_state);
void rmergeminvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state);
void rmergeminrv(ae_int_t n, RMatrix *x, ae_int_t rowidx, RVector *y, ae_state *_state);
void raddrv(ae_int_t n, double alpha, RMatrix *y, ae_int_t ridx, RVector *x, ae_state *_state);
void raddrr(ae_int_t n, double alpha, RMatrix *y, ae_int_t ridxsrc, RMatrix *x, ae_int_t ridxdst, ae_state *_state);
void rmulv(ae_int_t n, double v, RVector *x, ae_state *_state);
void rmulr(ae_int_t n, double v, RMatrix *x, ae_int_t rowidx, ae_state *_state);
void rmulvx(ae_int_t n, double v, RVector *x, ae_int_t offsx, ae_state *_state);
double rmaxv(ae_int_t n, RVector *x, ae_state *_state);
double rmaxabsv(ae_int_t n, RVector *x, ae_state *_state);
double rmaxr(ae_int_t n, RMatrix *x, ae_int_t rowidx, ae_state *_state);
double rmaxabsr(ae_int_t n, RMatrix *x, ae_int_t rowidx, ae_state *_state);
void bsetv(ae_int_t n, bool v, BVector *x, ae_state *_state);
void isetv(ae_int_t n, ae_int_t v, ZVector *x, ae_state *_state);
void rsetv(ae_int_t n, double v, RVector *x, ae_state *_state);
void rsetvx(ae_int_t n, double v, RVector *x, ae_int_t offsx, ae_state *_state);
void rsetm(ae_int_t m, ae_int_t n, double v, RMatrix *a, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
void bsetallocv(ae_int_t n, bool v, BVector *x, ae_state *_state);
void isetallocv(ae_int_t n, ae_int_t v, ZVector *x, ae_state *_state);
void rsetallocv(ae_int_t n, double v, RVector *x, ae_state *_state);
void rsetallocm(ae_int_t m, ae_int_t n, double v, RMatrix *a, ae_state *_state);
void ballocv(ae_int_t n, BVector *x, ae_state *_state);
void iallocv(ae_int_t n, ZVector *x, ae_state *_state);
void rallocv(ae_int_t n, RVector *x, ae_state *_state);
void rallocm(ae_int_t m, ae_int_t n, RMatrix *a, ae_state *_state);
#ifdef ALGLIB_NO_FAST_KERNELS
void rsetr(ae_int_t n, double v, RMatrix *a, ae_int_t i, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
void rsetc(ae_int_t n, double v, RMatrix *a, ae_int_t j, ae_state *_state);
#ifdef ALGLIB_NO_FAST_KERNELS
void bcopyv(ae_int_t n, BVector *x, BVector *y, ae_state *_state);
void icopyv(ae_int_t n, ZVector *x, ZVector *y, ae_state *_state);
void rcopyv(ae_int_t n, RVector *x, RVector *y, ae_state *_state);
void icopyvx(ae_int_t n, ZVector *x, ae_int_t offsx, ZVector *y, ae_int_t offsy, ae_state *_state);
void rcopyvx(ae_int_t n, RVector *x, ae_int_t offsx, RVector *y, ae_int_t offsy, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
void bcopyallocv(ae_int_t n, BVector *x, BVector *y, ae_state *_state);
void icopyallocv(ae_int_t n, ZVector *x, ZVector *y, ae_state *_state);
void rcopyallocv(ae_int_t n, RVector *x, RVector *y, ae_state *_state);
void rcopym(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y, ae_state *_state);
void rcopyallocm(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y, ae_state *_state);
void igrowv(ae_int_t newn, ZVector *x, ae_state *_state);
void rgrowv(ae_int_t newn, RVector *x, ae_state *_state);
#ifdef ALGLIB_NO_FAST_KERNELS
void rcopymulv(ae_int_t n, double v, RVector *x, RVector *y, ae_state *_state);
void rcopymulvr(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t ridx, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
void rcopymulvc(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t cidx, ae_state *_state);
#ifdef ALGLIB_NO_FAST_KERNELS
void rcopyvr(ae_int_t n, RVector *x, RMatrix *a, ae_int_t i, ae_state *_state);
void rcopyrv(ae_int_t n, RMatrix *a, ae_int_t i, RVector *x, ae_state *_state);
void rcopyrr(ae_int_t n, RMatrix *a, ae_int_t i, RMatrix *b, ae_int_t k, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
void rcopyvc(ae_int_t n, RVector *x, RMatrix *a, ae_int_t j, ae_state *_state);
void rcopycv(ae_int_t n, RMatrix *a, ae_int_t j, RVector *x, ae_state *_state);
#ifdef ALGLIB_NO_FAST_KERNELS
void rgemv(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t opa, RVector *x, double beta, RVector *y, ae_state *_state);
void rgemvx(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
void rger(ae_int_t m, ae_int_t n, double alpha, RVector *u, RVector *v, RMatrix *a, ae_state *_state);
void rtrsvx(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix, ae_state *_state);
#endif // defined ALGLIB_NO_FAST_KERNELS
bool rmatrixgerf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double ralpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state);
bool rmatrixrank1f(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state);
bool cmatrixrank1f(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv, ae_state *_state);
bool rmatrixlefttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool cmatrixlefttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool rmatrixrighttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool cmatrixrighttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
bool cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
bool cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void rmatrixgemmk44v00(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void rmatrixgemmk44v01(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void rmatrixgemmk44v10(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void rmatrixgemmk44v11(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void rmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
void cmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
} // end of namespace alglib_impl

// === HBLAS Package ===
namespace alglib_impl {
void hermitianmatrixvectormultiply(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, complex alpha, CVector *y, ae_state *_state);
void hermitianrank2update(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, CVector *y, CVector *t, complex alpha, ae_state *_state);
} // end of namespace alglib_impl

// === CREFLECTIONS Package ===
namespace alglib_impl {
void complexgeneratereflection(CVector *x, ae_int_t n, complex *tau, ae_state *_state);
void complexapplyreflectionfromtheleft(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work, ae_state *_state);
void complexapplyreflectionfromtheright(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work, ae_state *_state);
} // end of namespace alglib_impl

// === SBLAS Package ===
// Depends on: APSERV
namespace alglib_impl {
void symmetricmatrixvectormultiply(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, double alpha, RVector *y, ae_state *_state);
void symmetricrank2update(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, RVector *y, RVector *t, double alpha, ae_state *_state);
} // end of namespace alglib_impl

// === ABLASMKL Package ===
namespace alglib_impl {
bool rmatrixgermkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double alpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state);
bool rmatrixrank1mkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state);
bool cmatrixrank1mkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv, ae_state *_state);
bool rmatrixmvmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, RVector *y, ae_int_t iy, ae_state *_state);
bool cmatrixmvmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CVector *x, ae_int_t ix, CVector *y, ae_int_t iy, ae_state *_state);
bool rmatrixgemvmkl(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
bool rmatrixtrsvmkl(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix, ae_state *_state);
bool rmatrixsymvmkl(ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
bool rmatrixsyrkmkl(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
bool cmatrixherkmkl(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state);
bool rmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
bool cmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state);
bool rmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool cmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool rmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool cmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state);
bool spdmatrixcholeskymkl(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, bool *cholresult, ae_state *_state);
bool rmatrixplumkl(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state);
bool rmatrixbdmkl(RMatrix *a, ae_int_t m, ae_int_t n, RVector *d, RVector *e, RVector *tauq, RVector *taup, ae_state *_state);
bool rmatrixbdmultiplybymkl(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, RVector *taup, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool byq, bool fromtheright, bool dotranspose, ae_state *_state);
bool rmatrixhessenbergmkl(RMatrix *a, ae_int_t n, RVector *tau, ae_state *_state);
bool rmatrixhessenbergunpackqmkl(RMatrix *a, ae_int_t n, RVector *tau, RMatrix *q, ae_state *_state);
bool smatrixtdmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RVector *d, RVector *e, ae_state *_state);
bool hmatrixtdmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, RVector *d, RVector *e, ae_state *_state);
bool smatrixtdunpackqmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RMatrix *q, ae_state *_state);
bool hmatrixtdunpackqmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, CMatrix *q, ae_state *_state);
bool rmatrixbdsvdmkl(RVector *d, RVector *e, ae_int_t n, bool isupper, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt, bool *svdresult, ae_state *_state);
bool rmatrixinternalschurdecompositionmkl(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info, ae_state *_state);
bool rmatrixinternaltrevcmkl(RMatrix *t, ae_int_t n, ae_int_t side, ae_int_t howmny, RMatrix *vl, RMatrix *vr, ae_int_t *m, ae_int_t *info, ae_state *_state);
bool smatrixtdevdmkl(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z, bool *evdresult, ae_state *_state);
bool sparsegemvcrsmkl(ae_int_t opa, ae_int_t arows, ae_int_t acols, double alpha, RVector *vals, ZVector *cidx, ZVector *ridx, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state);
} // end of namespace alglib_impl

// === SCODES Package ===
namespace alglib_impl {
ae_int_t getrdfserializationcode(ae_state *_state);
ae_int_t getkdtreeserializationcode(ae_state *_state);
ae_int_t getmlpserializationcode(ae_state *_state);
ae_int_t getmlpeserializationcode(ae_state *_state);
ae_int_t getrbfserializationcode(ae_state *_state);
ae_int_t getspline2dserializationcode(ae_state *_state);
ae_int_t getidwserializationcode(ae_state *_state);
ae_int_t getsparsematrixserializationcode(ae_state *_state);
ae_int_t getknnserializationcode(ae_state *_state);
ae_int_t getlptestserializationcode(ae_state *_state);
} // end of namespace alglib_impl

// === TSORT Package ===
// Depends on: APSERV
namespace alglib_impl {
void tagsortfasti(RVector *a, ZVector *b, RVector *bufa, ZVector *bufb, ae_int_t n, ae_state *_state);
void tagsortfastr(RVector *a, RVector *b, RVector *bufa, RVector *bufb, ae_int_t n, ae_state *_state);
void tagsortfast(RVector *a, RVector *bufa, ae_int_t n, ae_state *_state);
void tagsortmiddleir(ZVector *a, RVector *b, ae_int_t offset, ae_int_t n, ae_state *_state);
void tagsortmiddlei(ZVector *a, ae_int_t offset, ae_int_t n, ae_state *_state);
void sortmiddlei(ZVector *a, ae_int_t offset, ae_int_t n, ae_state *_state);
void tagsortbuf(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2, apbuffers *buf, ae_state *_state);
void tagsort(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2, ae_state *_state);
void tagheappushi(RVector *a, ZVector *b, ae_int_t *n, double va, ae_int_t vb, ae_state *_state);
void tagheapreplacetopi(RVector *a, ZVector *b, ae_int_t n, double va, ae_int_t vb, ae_state *_state);
void tagheappopi(RVector *a, ZVector *b, ae_int_t *n, ae_state *_state);
ae_int_t lowerbound(RVector *a, ae_int_t n, double t, ae_state *_state);
ae_int_t upperbound(RVector *a, ae_int_t n, double t, ae_state *_state);
} // end of namespace alglib_impl

// === BLAS Package ===
namespace alglib_impl {
double vectornorm2(RVector *x, ae_int_t i1, ae_int_t i2, ae_state *_state);
ae_int_t vectoridxabsmax(RVector *x, ae_int_t i1, ae_int_t i2, ae_state *_state);
ae_int_t columnidxabsmax(RMatrix *x, ae_int_t i1, ae_int_t i2, ae_int_t j, ae_state *_state);
ae_int_t rowidxabsmax(RMatrix *x, ae_int_t j1, ae_int_t j2, ae_int_t i, ae_state *_state);
double upperhessenberg1norm(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work, ae_state *_state);
void copymatrix(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2, ae_state *_state);
void inplacetranspose(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work, ae_state *_state);
void copyandtranspose(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2, ae_state *_state);
void matrixvectormultiply(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, bool trans, RVector *x, ae_int_t ix1, ae_int_t ix2, double alpha, RVector *y, ae_int_t iy1, ae_int_t iy2, double beta, ae_state *_state);
double pythag2(double x, double y, ae_state *_state);
void matrixmatrixmultiply(RMatrix *a, ae_int_t ai1, ae_int_t ai2, ae_int_t aj1, ae_int_t aj2, bool transa, RMatrix *b, ae_int_t bi1, ae_int_t bi2, ae_int_t bj1, ae_int_t bj2, bool transb, double alpha, RMatrix *c, ae_int_t ci1, ae_int_t ci2, ae_int_t cj1, ae_int_t cj2, double beta, RVector *work, ae_state *_state);
} // end of namespace alglib_impl

// === ROTATIONS Package ===
namespace alglib_impl {
void applyrotationsfromtheleft(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work, ae_state *_state);
void applyrotationsfromtheright(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work, ae_state *_state);
void generaterotation(double f, double g, double *cs, double *sn, double *r, ae_state *_state);
} // end of namespace alglib_impl

// === BASICSTATOPS Package ===
// Depends on: TSORT
namespace alglib_impl {
void rankx(RVector *x, ae_int_t n, bool iscentered, apbuffers *buf, ae_state *_state);
void rankxuntied(RVector *x, ae_int_t n, apbuffers *buf, ae_state *_state);
} // end of namespace alglib_impl

// === TRLINSOLVE Package ===
namespace alglib_impl {
void safesolvetriangular(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit, bool normin, RVector *cnorm, ae_state *_state);
void rmatrixtrsafesolve(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit, ae_state *_state);
} // end of namespace alglib_impl

// === SAFESOLVE Package ===
namespace alglib_impl {
bool rmatrixscaledtrsafesolve(RMatrix *a, double sa, ae_int_t n, RVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth, ae_state *_state);
bool cmatrixscaledtrsafesolve(CMatrix *a, double sa, ae_int_t n, CVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth, ae_state *_state);
} // end of namespace alglib_impl

// === XBLAS Package ===
namespace alglib_impl {
void xdot(RVector *a, RVector *b, ae_int_t n, RVector *temp, double *r, double *rerr, ae_state *_state);
void xcdot(CVector *a, CVector *b, ae_int_t n, RVector *temp, complex *r, double *rerr, ae_state *_state);
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
void linminstate_init(void *_p, ae_state *_state, bool make_automatic);
void linminstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
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
   rcommstate rstate;
};
void armijostate_init(void *_p, ae_state *_state, bool make_automatic);
void armijostate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void armijostate_free(void *_p, bool make_automatic);

void linminnormalized(RVector *d, double *stp, ae_int_t n, ae_state *_state);
void mcsrch(ae_int_t n, RVector *x, double *f, RVector *g, RVector *s, double *stp, double stpmax, double gtol, ae_int_t *info, ae_int_t *nfev, RVector *wa, linminstate *state, ae_int_t *stage, ae_state *_state);
void armijocreate(ae_int_t n, RVector *x, double f, RVector *s, double stp, double stpmax, ae_int_t fmax, armijostate *state, ae_state *_state);
bool armijoiteration(armijostate *state, ae_state *_state);
void armijoresults(armijostate *state, ae_int_t *info, double *stp, double *f, ae_state *_state);
} // end of namespace alglib_impl

// === NEARUNITYUNIT Package ===
namespace alglib_impl {
double nulog1p(double x, ae_state *_state);
double nuexpm1(double x, ae_state *_state);
double nucosm1(double x, ae_state *_state);
} // end of namespace alglib_impl

// === NTHEORY Package ===
namespace alglib_impl {
void findprimitiverootandinverse(ae_int_t n, ae_int_t *proot, ae_int_t *invproot, ae_state *_state);
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
void fasttransformplan_init(void *_p, ae_state *_state, bool make_automatic);
void fasttransformplan_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void fasttransformplan_free(void *_p, bool make_automatic);

double ftbasegetflopestimate(ae_int_t n, ae_state *_state);
void ftbasefactorize(ae_int_t n, ae_int_t tasktype, ae_int_t *n1, ae_int_t *n2, ae_state *_state);
bool ftbaseissmooth(ae_int_t n, ae_state *_state);
ae_int_t ftbasefindsmooth(ae_int_t n, ae_state *_state);
ae_int_t ftbasefindsmootheven(ae_int_t n, ae_state *_state);
void ftapplyplan(fasttransformplan *plan, RVector *a, ae_int_t offsa, ae_int_t repcnt, ae_state *_state);
void ftcomplexfftplan(ae_int_t n, ae_int_t k, fasttransformplan *plan, ae_state *_state);
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
void mlpbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void mlpbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlpbuffers_free(void *_p, bool make_automatic);

void hpcpreparechunkedgradient(RVector *weights, ae_int_t wcount, ae_int_t ntotal, ae_int_t nin, ae_int_t nout, mlpbuffers *buf, ae_state *_state);
void hpcfinalizechunkedgradient(mlpbuffers *buf, RVector *grad, ae_state *_state);
bool hpcchunkedgradient(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf, double *e, bool naturalerrorfunc, ae_state *_state);
bool hpcchunkedprocess(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf, ae_state *_state);
} // end of namespace alglib_impl

#endif // OnceOnly
