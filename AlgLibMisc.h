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
#ifndef OnceOnlyAlgLibMisc_h
#define OnceOnlyAlgLibMisc_h

#include "AlgLibInternal.h"

// === HQRND Package ===
// Depends on: (AlgLibInternal) APSERV, ABLASF
namespace alglib_impl {
struct hqrndstate {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t magicv;
};
void hqrndstate_init(void *_p, bool make_automatic);
void hqrndstate_copy(void *_dst, void *_src, bool make_automatic);
void hqrndstate_free(void *_p, bool make_automatic);

void hqrndrandomize(hqrndstate *state);
void hqrndseed(ae_int_t s1, ae_int_t s2, hqrndstate *state);
double hqrnduniformr(hqrndstate *state);
double hqrndmiduniformr(hqrndstate *state);
ae_int_t hqrnduniformi(hqrndstate *state, ae_int_t n);
void hqrndnormal2(hqrndstate *state, double *x1, double *x2);
double hqrndnormal(hqrndstate *state);
void hqrndnormalv(hqrndstate *state, ae_int_t n, RVector *x);
void hqrndnormalm(hqrndstate *state, ae_int_t m, ae_int_t n, RMatrix *x);
void hqrndunit2(hqrndstate *state, double *x, double *y);
double hqrndexponential(hqrndstate *state, double lambdav);
double hqrnddiscrete(hqrndstate *state, RVector *x, ae_int_t n);
double hqrndcontinuous(hqrndstate *state, RVector *x, ae_int_t n);
} // end of namespace alglib_impl

namespace alglib {
DecClass(hqrndstate, );

void hqrndrandomize(hqrndstate &state);
void hqrndseed(const ae_int_t s1, const ae_int_t s2, hqrndstate &state);
double hqrnduniformr(const hqrndstate &state);
double hqrndmiduniformr(const hqrndstate &state);
ae_int_t hqrnduniformi(const hqrndstate &state, const ae_int_t n);
void hqrndnormal2(const hqrndstate &state, double &x1, double &x2);
double hqrndnormal(const hqrndstate &state);
void hqrndnormalv(const hqrndstate &state, const ae_int_t n, real_1d_array &x);
void hqrndnormalm(const hqrndstate &state, const ae_int_t m, const ae_int_t n, real_2d_array &x);
void hqrndunit2(const hqrndstate &state, double &x, double &y);
double hqrndexponential(const hqrndstate &state, const double lambdav);
double hqrnddiscrete(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
double hqrndcontinuous(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
} // end of namespace alglib

// === XDEBUG Package ===
namespace alglib_impl {
struct xdebugrecord1 {
   ae_int_t i;
   complex c;
   ae_vector a;
};
void xdebugrecord1_init(void *_p, bool make_automatic);
void xdebugrecord1_copy(void *_dst, void *_src, bool make_automatic);
void xdebugrecord1_free(void *_p, bool make_automatic);

void xdebuginitrecord1(xdebugrecord1 *rec1);
ae_int_t xdebugb1count(BVector *a);
void xdebugb1not(BVector *a);
void xdebugb1appendcopy(BVector *a);
void xdebugb1outeven(ae_int_t n, BVector *a);
ae_int_t xdebugi1sum(ZVector *a);
void xdebugi1neg(ZVector *a);
void xdebugi1appendcopy(ZVector *a);
void xdebugi1outeven(ae_int_t n, ZVector *a);
double xdebugr1sum(RVector *a);
void xdebugr1neg(RVector *a);
void xdebugr1appendcopy(RVector *a);
void xdebugr1outeven(ae_int_t n, RVector *a);
complex xdebugc1sum(CVector *a);
void xdebugc1neg(CVector *a);
void xdebugc1appendcopy(CVector *a);
void xdebugc1outeven(ae_int_t n, CVector *a);
ae_int_t xdebugb2count(BMatrix *a);
void xdebugb2not(BMatrix *a);
void xdebugb2transpose(BMatrix *a);
void xdebugb2outsin(ae_int_t m, ae_int_t n, BMatrix *a);
ae_int_t xdebugi2sum(ZMatrix *a);
void xdebugi2neg(ZMatrix *a);
void xdebugi2transpose(ZMatrix *a);
void xdebugi2outsin(ae_int_t m, ae_int_t n, ZMatrix *a);
double xdebugr2sum(RMatrix *a);
void xdebugr2neg(RMatrix *a);
void xdebugr2transpose(RMatrix *a);
void xdebugr2outsin(ae_int_t m, ae_int_t n, RMatrix *a);
complex xdebugc2sum(CMatrix *a);
void xdebugc2neg(CMatrix *a);
void xdebugc2transpose(CMatrix *a);
void xdebugc2outsincos(ae_int_t m, ae_int_t n, CMatrix *a);
double xdebugmaskedbiasedproductsum(ae_int_t m, ae_int_t n, RMatrix *a, RMatrix *b, BMatrix *c);
} // end of namespace alglib_impl

namespace alglib {
DecClass(xdebugrecord1, ae_int_t &i; complex &c; real_1d_array a;);

void xdebuginitrecord1(xdebugrecord1 &rec1);
ae_int_t xdebugb1count(const boolean_1d_array &a);
void xdebugb1not(const boolean_1d_array &a);
void xdebugb1appendcopy(boolean_1d_array &a);
void xdebugb1outeven(const ae_int_t n, boolean_1d_array &a);
ae_int_t xdebugi1sum(const integer_1d_array &a);
void xdebugi1neg(const integer_1d_array &a);
void xdebugi1appendcopy(integer_1d_array &a);
void xdebugi1outeven(const ae_int_t n, integer_1d_array &a);
double xdebugr1sum(const real_1d_array &a);
void xdebugr1neg(const real_1d_array &a);
void xdebugr1appendcopy(real_1d_array &a);
void xdebugr1outeven(const ae_int_t n, real_1d_array &a);
complex xdebugc1sum(const complex_1d_array &a);
void xdebugc1neg(const complex_1d_array &a);
void xdebugc1appendcopy(complex_1d_array &a);
void xdebugc1outeven(const ae_int_t n, complex_1d_array &a);
ae_int_t xdebugb2count(const boolean_2d_array &a);
void xdebugb2not(const boolean_2d_array &a);
void xdebugb2transpose(boolean_2d_array &a);
void xdebugb2outsin(const ae_int_t m, const ae_int_t n, boolean_2d_array &a);
ae_int_t xdebugi2sum(const integer_2d_array &a);
void xdebugi2neg(const integer_2d_array &a);
void xdebugi2transpose(integer_2d_array &a);
void xdebugi2outsin(const ae_int_t m, const ae_int_t n, integer_2d_array &a);
double xdebugr2sum(const real_2d_array &a);
void xdebugr2neg(const real_2d_array &a);
void xdebugr2transpose(real_2d_array &a);
void xdebugr2outsin(const ae_int_t m, const ae_int_t n, real_2d_array &a);
complex xdebugc2sum(const complex_2d_array &a);
void xdebugc2neg(const complex_2d_array &a);
void xdebugc2transpose(complex_2d_array &a);
void xdebugc2outsincos(const ae_int_t m, const ae_int_t n, complex_2d_array &a);
double xdebugmaskedbiasedproductsum(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const real_2d_array &b, const boolean_2d_array &c);
} // end of namespace alglib

// === NEARESTNEIGHBOR Package ===
// Depends on: (AlgLibInternal) SCODES, TSORT
namespace alglib_impl {
struct kdtreerequestbuffer {
   ae_vector x;
   ae_vector boxmin;
   ae_vector boxmax;
   ae_int_t kneeded;
   double rneeded;
   bool selfmatch;
   double approxf;
   ae_int_t kcur;
   ae_vector idx;
   ae_vector r;
   ae_vector buf;
   ae_vector curboxmin;
   ae_vector curboxmax;
   double curdist;
};
void kdtreerequestbuffer_init(void *_p, bool make_automatic);
void kdtreerequestbuffer_copy(void *_dst, void *_src, bool make_automatic);
void kdtreerequestbuffer_free(void *_p, bool make_automatic);

struct kdtree {
   ae_int_t n;
   ae_int_t nx;
   ae_int_t ny;
   ae_int_t normtype;
   ae_matrix xy;
   ae_vector tags;
   ae_vector boxmin;
   ae_vector boxmax;
   ae_vector nodes;
   ae_vector splits;
   kdtreerequestbuffer innerbuf;
   ae_int_t debugcounter;
};
void kdtree_init(void *_p, bool make_automatic);
void kdtree_copy(void *_dst, void *_src, bool make_automatic);
void kdtree_free(void *_p, bool make_automatic);
void kdtreealloc(ae_serializer *s, kdtree *tree);
void kdtreeserialize(ae_serializer *s, kdtree *tree);
void kdtreeunserialize(ae_serializer *s, kdtree *tree);

void kdtreecreaterequestbuffer(kdtree *kdt, kdtreerequestbuffer *buf);
void kdtreebuildtagged(RMatrix *xy, ZVector *tags, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree *kdt);
void kdtreebuild(RMatrix *xy, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree *kdt);
ae_int_t kdtreetsqueryknn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, ae_int_t k, bool selfmatch);
ae_int_t kdtreequeryknn(kdtree *kdt, RVector *x, ae_int_t k, bool selfmatch);
ae_int_t kdtreetsqueryrnn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch);
ae_int_t kdtreequeryrnn(kdtree *kdt, RVector *x, double r, bool selfmatch);
ae_int_t kdtreetsqueryrnnu(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch);
ae_int_t kdtreequeryrnnu(kdtree *kdt, RVector *x, double r, bool selfmatch);
ae_int_t kdtreetsqueryaknn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, ae_int_t k, bool selfmatch, double eps);
ae_int_t kdtreequeryaknn(kdtree *kdt, RVector *x, ae_int_t k, bool selfmatch, double eps);
ae_int_t kdtreetsquerybox(kdtree *kdt, kdtreerequestbuffer *buf, RVector *boxmin, RVector *boxmax);
ae_int_t kdtreequerybox(kdtree *kdt, RVector *boxmin, RVector *boxmax);
void kdtreetsqueryresultsx(kdtree *kdt, kdtreerequestbuffer *buf, RMatrix *x);
void kdtreequeryresultsx(kdtree *kdt, RMatrix *x);
void kdtreetsqueryresultsxy(kdtree *kdt, kdtreerequestbuffer *buf, RMatrix *xy);
void kdtreequeryresultsxy(kdtree *kdt, RMatrix *xy);
void kdtreetsqueryresultstags(kdtree *kdt, kdtreerequestbuffer *buf, ZVector *tags);
void kdtreequeryresultstags(kdtree *kdt, ZVector *tags);
void kdtreetsqueryresultsdistances(kdtree *kdt, kdtreerequestbuffer *buf, RVector *r);
void kdtreequeryresultsdistances(kdtree *kdt, RVector *r);
void kdtreequeryresultsxi(kdtree *kdt, RMatrix *x);
void kdtreequeryresultsxyi(kdtree *kdt, RMatrix *xy);
void kdtreequeryresultstagsi(kdtree *kdt, ZVector *tags);
void kdtreequeryresultsdistancesi(kdtree *kdt, RVector *r);
void kdtreeexplorebox(kdtree *kdt, RVector *boxmin, RVector *boxmax);
void kdtreeexplorenodetype(kdtree *kdt, ae_int_t node, ae_int_t *nodetype);
void kdtreeexploreleaf(kdtree *kdt, ae_int_t node, RMatrix *xy, ae_int_t *k);
void kdtreeexploresplit(kdtree *kdt, ae_int_t node, ae_int_t *d, double *s, ae_int_t *nodele, ae_int_t *nodege);
} // end of namespace alglib_impl

namespace alglib {
DecClass(kdtreerequestbuffer, );
DecClass(kdtree, );
void kdtreeserialize(kdtree &obj, std::string &s_out);
void kdtreeserialize(kdtree &obj, std::ostream &s_out);
void kdtreeunserialize(const std::string &s_in, kdtree &obj);
void kdtreeunserialize(const std::istream &s_in, kdtree &obj);

void kdtreecreaterequestbuffer(const kdtree &kdt, kdtreerequestbuffer &buf);
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuild(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuild(const real_2d_array &xy, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k);
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k);
ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r);
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r);
ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r);
ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r);
ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const double eps);
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const double eps);
ae_int_t kdtreetsquerybox(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &boxmin, const real_1d_array &boxmax);
ae_int_t kdtreequerybox(const kdtree &kdt, const real_1d_array &boxmin, const real_1d_array &boxmax);
void kdtreetsqueryresultsx(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &x);
void kdtreequeryresultsx(const kdtree &kdt, real_2d_array &x);
void kdtreetsqueryresultsxy(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &xy);
void kdtreequeryresultsxy(const kdtree &kdt, real_2d_array &xy);
void kdtreetsqueryresultstags(const kdtree &kdt, const kdtreerequestbuffer &buf, integer_1d_array &tags);
void kdtreequeryresultstags(const kdtree &kdt, integer_1d_array &tags);
void kdtreetsqueryresultsdistances(const kdtree &kdt, const kdtreerequestbuffer &buf, real_1d_array &r);
void kdtreequeryresultsdistances(const kdtree &kdt, real_1d_array &r);
void kdtreequeryresultsxi(const kdtree &kdt, real_2d_array &x);
void kdtreequeryresultsxyi(const kdtree &kdt, real_2d_array &xy);
void kdtreequeryresultstagsi(const kdtree &kdt, integer_1d_array &tags);
void kdtreequeryresultsdistancesi(const kdtree &kdt, real_1d_array &r);
} // end of namespace alglib

#endif // OnceOnly
