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
#ifndef OnceOnlyInterpolation_h
#define OnceOnlyInterpolation_h

#include "Integration.h"
#include "Optimization.h"

// === RATINT Package ===
// Depends on: (AlgLibInternal) TSORT
namespace alglib_impl {
struct barycentricinterpolant {
   ae_int_t n;
   double sy;
   ae_vector x;
   ae_vector y;
   ae_vector w;
};
void barycentricinterpolant_init(void *_p, ae_state *_state, bool make_automatic);
void barycentricinterpolant_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void barycentricinterpolant_free(void *_p, bool make_automatic);

double barycentriccalc(barycentricinterpolant *b, double t, ae_state *_state);
void barycentricdiff1(barycentricinterpolant *b, double t, double *f, double *df, ae_state *_state);
void barycentricdiff2(barycentricinterpolant *b, double t, double *f, double *df, double *d2f, ae_state *_state);
void barycentriclintransx(barycentricinterpolant *b, double ca, double cb, ae_state *_state);
void barycentriclintransy(barycentricinterpolant *b, double ca, double cb, ae_state *_state);
void barycentricunpack(barycentricinterpolant *b, ae_int_t *n, RVector *x, RVector *y, RVector *w, ae_state *_state);
void barycentricbuildxyw(RVector *x, RVector *y, RVector *w, ae_int_t n, barycentricinterpolant *b, ae_state *_state);
void barycentricbuildfloaterhormann(RVector *x, RVector *y, ae_int_t n, ae_int_t d, barycentricinterpolant *b, ae_state *_state);
void barycentriccopy(barycentricinterpolant *b, barycentricinterpolant *b2, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(barycentricinterpolant, );

double barycentriccalc(const barycentricinterpolant &b, const double t, const xparams _xparams = NonTH);
void barycentricdiff1(const barycentricinterpolant &b, const double t, double &f, double &df, const xparams _xparams = NonTH);
void barycentricdiff2(const barycentricinterpolant &b, const double t, double &f, double &df, double &d2f, const xparams _xparams = NonTH);
void barycentriclintransx(const barycentricinterpolant &b, const double ca, const double cb, const xparams _xparams = NonTH);
void barycentriclintransy(const barycentricinterpolant &b, const double ca, const double cb, const xparams _xparams = NonTH);
void barycentricunpack(const barycentricinterpolant &b, ae_int_t &n, real_1d_array &x, real_1d_array &y, real_1d_array &w, const xparams _xparams = NonTH);
void barycentricbuildxyw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, barycentricinterpolant &b, const xparams _xparams = NonTH);
void barycentricbuildfloaterhormann(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t d, barycentricinterpolant &b, const xparams _xparams = NonTH);
} // end of namespace alglib

// === IDW Package ===
// Depends on: (AlgLibMisc) HQRND, NEARESTNEIGHBOR
// Depends on: (LinAlg) ABLAS
namespace alglib_impl {
struct idwcalcbuffer {
   ae_vector x;
   ae_vector y;
   ae_vector tsyw;
   ae_vector tsw;
   ae_matrix tsxy;
   ae_vector tsdist;
   kdtreerequestbuffer requestbuffer;
};
void idwcalcbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void idwcalcbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void idwcalcbuffer_free(void *_p, bool make_automatic);

struct idwmodel {
   ae_int_t nx;
   ae_int_t ny;
   ae_vector globalprior;
   ae_int_t algotype;
   ae_int_t nlayers;
   double r0;
   double rdecay;
   double lambda0;
   double lambdalast;
   double lambdadecay;
   double shepardp;
   kdtree tree;
   ae_int_t npoints;
   ae_vector shepardxy;
   idwcalcbuffer buffer;
};
void idwmodel_init(void *_p, ae_state *_state, bool make_automatic);
void idwmodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void idwmodel_free(void *_p, bool make_automatic);
void idwalloc(ae_serializer *s, idwmodel *model, ae_state *_state);
void idwserialize(ae_serializer *s, idwmodel *model, ae_state *_state);
void idwunserialize(ae_serializer *s, idwmodel *model, ae_state *_state);

struct idwbuilder {
   ae_int_t priortermtype;
   ae_vector priortermval;
   ae_int_t algotype;
   ae_int_t nlayers;
   double r0;
   double rdecay;
   double lambda0;
   double lambdalast;
   double lambdadecay;
   double shepardp;
   ae_vector xy;
   ae_int_t npoints;
   ae_int_t nx;
   ae_int_t ny;
   ae_matrix tmpxy;
   ae_matrix tmplayers;
   ae_vector tmptags;
   ae_vector tmpdist;
   ae_vector tmpx;
   ae_vector tmpwy;
   ae_vector tmpw;
   kdtree tmptree;
   ae_vector tmpmean;
};
void idwbuilder_init(void *_p, ae_state *_state, bool make_automatic);
void idwbuilder_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void idwbuilder_free(void *_p, bool make_automatic);

struct idwreport {
   double rmserror;
   double avgerror;
   double maxerror;
   double r2;
};
void idwreport_init(void *_p, ae_state *_state, bool make_automatic);
void idwreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void idwreport_free(void *_p, bool make_automatic);

void idwcreatecalcbuffer(idwmodel *s, idwcalcbuffer *buf, ae_state *_state);
void idwbuildercreate(ae_int_t nx, ae_int_t ny, idwbuilder *state, ae_state *_state);
void idwbuildersetnlayers(idwbuilder *state, ae_int_t nlayers, ae_state *_state);
void idwbuildersetpoints(idwbuilder *state, RMatrix *xy, ae_int_t n, ae_state *_state);
void idwbuildersetalgomstab(idwbuilder *state, double srad, ae_state *_state);
void idwbuildersetalgotextbookshepard(idwbuilder *state, double p, ae_state *_state);
void idwbuildersetalgotextbookmodshepard(idwbuilder *state, double r, ae_state *_state);
void idwbuildersetuserterm(idwbuilder *state, double v, ae_state *_state);
void idwbuildersetconstterm(idwbuilder *state, ae_state *_state);
void idwbuildersetzeroterm(idwbuilder *state, ae_state *_state);
void idwtscalcbuf(idwmodel *s, idwcalcbuffer *buf, RVector *x, RVector *y, ae_state *_state);
void idwcalc(idwmodel *s, RVector *x, RVector *y, ae_state *_state);
void idwcalcbuf(idwmodel *s, RVector *x, RVector *y, ae_state *_state);
double idwcalc1(idwmodel *s, double x0, ae_state *_state);
double idwcalc2(idwmodel *s, double x0, double x1, ae_state *_state);
double idwcalc3(idwmodel *s, double x0, double x1, double x2, ae_state *_state);
void idwfit(idwbuilder *state, idwmodel *model, idwreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(idwcalcbuffer, );
DecClass(idwmodel, );
DecClass(idwbuilder, );
DecClass(idwreport, double &rmserror; double &avgerror; double &maxerror; double &r2;);
void idwserialize(idwmodel &obj, std::string &s_out);
void idwserialize(idwmodel &obj, std::ostream &s_out);
void idwunserialize(const std::string &s_in, idwmodel &obj);
void idwunserialize(const std::istream &s_in, idwmodel &obj);

void idwcreatecalcbuffer(const idwmodel &s, idwcalcbuffer &buf, const xparams _xparams = NonTH);
void idwbuildercreate(const ae_int_t nx, const ae_int_t ny, idwbuilder &state, const xparams _xparams = NonTH);
void idwbuildersetnlayers(const idwbuilder &state, const ae_int_t nlayers, const xparams _xparams = NonTH);
void idwbuildersetpoints(const idwbuilder &state, const real_2d_array &xy, const ae_int_t n, const xparams _xparams = NonTH);
void idwbuildersetpoints(const idwbuilder &state, const real_2d_array &xy, const xparams _xparams = NonTH);
void idwbuildersetalgomstab(const idwbuilder &state, const double srad, const xparams _xparams = NonTH);
void idwbuildersetalgotextbookshepard(const idwbuilder &state, const double p, const xparams _xparams = NonTH);
void idwbuildersetalgotextbookmodshepard(const idwbuilder &state, const double r, const xparams _xparams = NonTH);
void idwbuildersetuserterm(const idwbuilder &state, const double v, const xparams _xparams = NonTH);
void idwbuildersetconstterm(const idwbuilder &state, const xparams _xparams = NonTH);
void idwbuildersetzeroterm(const idwbuilder &state, const xparams _xparams = NonTH);
void idwtscalcbuf(const idwmodel &s, const idwcalcbuffer &buf, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void idwcalc(const idwmodel &s, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void idwcalcbuf(const idwmodel &s, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
double idwcalc1(const idwmodel &s, const double x0, const xparams _xparams = NonTH);
double idwcalc2(const idwmodel &s, const double x0, const double x1, const xparams _xparams = NonTH);
double idwcalc3(const idwmodel &s, const double x0, const double x1, const double x2, const xparams _xparams = NonTH);
void idwfit(const idwbuilder &state, idwmodel &model, idwreport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === INTFITSERV Package ===
// Depends on: (LinAlg) TRFAC
namespace alglib_impl {
void lsfitscalexy(RVector *x, RVector *y, RVector *w, ae_int_t n, RVector *xc, RVector *yc, ZVector *dc, ae_int_t k, double *xa, double *xb, double *sa, double *sb, RVector *xoriginal, RVector *yoriginal, ae_state *_state);
void buildpriorterm(RMatrix *xy, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t modeltype, double priorval, RMatrix *v, ae_state *_state);
void buildpriorterm1(RVector *xy1, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t modeltype, double priorval, RMatrix *v, ae_state *_state);
} // end of namespace alglib_impl

// === POLINT Package ===
// Depends on: RATINT
namespace alglib_impl {
void polynomialbuild(RVector *x, RVector *y, ae_int_t n, barycentricinterpolant *p, ae_state *_state);
void polynomialbuildeqdist(double a, double b, RVector *y, ae_int_t n, barycentricinterpolant *p, ae_state *_state);
void polynomialbuildcheb1(double a, double b, RVector *y, ae_int_t n, barycentricinterpolant *p, ae_state *_state);
void polynomialbuildcheb2(double a, double b, RVector *y, ae_int_t n, barycentricinterpolant *p, ae_state *_state);
void polynomialbar2cheb(barycentricinterpolant *p, double a, double b, RVector *t, ae_state *_state);
void polynomialcheb2bar(RVector *t, ae_int_t n, double a, double b, barycentricinterpolant *p, ae_state *_state);
void polynomialbar2pow(barycentricinterpolant *p, double c, double s, RVector *a, ae_state *_state);
void polynomialpow2bar(RVector *a, ae_int_t n, double c, double s, barycentricinterpolant *p, ae_state *_state);
double polynomialcalceqdist(double a, double b, RVector *f, ae_int_t n, double t, ae_state *_state);
double polynomialcalccheb1(double a, double b, RVector *f, ae_int_t n, double t, ae_state *_state);
double polynomialcalccheb2(double a, double b, RVector *f, ae_int_t n, double t, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void polynomialbuild(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuild(const real_1d_array &x, const real_1d_array &y, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuildeqdist(const double a, const double b, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuildeqdist(const double a, const double b, const real_1d_array &y, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuildcheb1(const double a, const double b, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuildcheb1(const double a, const double b, const real_1d_array &y, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuildcheb2(const double a, const double b, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbuildcheb2(const double a, const double b, const real_1d_array &y, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbar2cheb(const barycentricinterpolant &p, const double a, const double b, real_1d_array &t, const xparams _xparams = NonTH);
void polynomialcheb2bar(const real_1d_array &t, const ae_int_t n, const double a, const double b, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialcheb2bar(const real_1d_array &t, const double a, const double b, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialbar2pow(const barycentricinterpolant &p, const double c, const double s, real_1d_array &a, const xparams _xparams = NonTH);
void polynomialbar2pow(const barycentricinterpolant &p, real_1d_array &a, const xparams _xparams = NonTH);
void polynomialpow2bar(const real_1d_array &a, const ae_int_t n, const double c, const double s, barycentricinterpolant &p, const xparams _xparams = NonTH);
void polynomialpow2bar(const real_1d_array &a, barycentricinterpolant &p, const xparams _xparams = NonTH);
double polynomialcalceqdist(const double a, const double b, const real_1d_array &f, const ae_int_t n, const double t, const xparams _xparams = NonTH);
double polynomialcalceqdist(const double a, const double b, const real_1d_array &f, const double t, const xparams _xparams = NonTH);
double polynomialcalccheb1(const double a, const double b, const real_1d_array &f, const ae_int_t n, const double t, const xparams _xparams = NonTH);
double polynomialcalccheb1(const double a, const double b, const real_1d_array &f, const double t, const xparams _xparams = NonTH);
double polynomialcalccheb2(const double a, const double b, const real_1d_array &f, const ae_int_t n, const double t, const xparams _xparams = NonTH);
double polynomialcalccheb2(const double a, const double b, const real_1d_array &f, const double t, const xparams _xparams = NonTH);
} // end of namespace alglib

// === SPLINE1D Package ===
// Depends on: (LinAlg) FBLS
// Depends on: (Solvers) LINLSQR
// Depends on: INTFITSERV
namespace alglib_impl {
struct spline1dinterpolant {
   bool periodic;
   ae_int_t n;
   ae_int_t k;
   ae_int_t continuity;
   ae_vector x;
   ae_vector c;
};
void spline1dinterpolant_init(void *_p, ae_state *_state, bool make_automatic);
void spline1dinterpolant_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline1dinterpolant_free(void *_p, bool make_automatic);

struct spline1dfitreport {
   double taskrcond;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double maxerror;
};
void spline1dfitreport_init(void *_p, ae_state *_state, bool make_automatic);
void spline1dfitreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline1dfitreport_free(void *_p, bool make_automatic);

void heapsortdpoints(RVector *x, RVector *y, RVector *d, ae_int_t n, ae_state *_state);
void spline1dbuildlinear(RVector *x, RVector *y, ae_int_t n, spline1dinterpolant *c, ae_state *_state);
void spline1dbuildhermite(RVector *x, RVector *y, RVector *d, ae_int_t n, spline1dinterpolant *c, ae_state *_state);
void spline1dbuildcubic(RVector *x, RVector *y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, spline1dinterpolant *c, ae_state *_state);
void spline1dbuildcatmullrom(RVector *x, RVector *y, ae_int_t n, ae_int_t boundtype, double tension, spline1dinterpolant *c, ae_state *_state);
void spline1dbuildakima(RVector *x, RVector *y, ae_int_t n, spline1dinterpolant *c, ae_state *_state);
void spline1dbuildmonotone(RVector *x, RVector *y, ae_int_t n, spline1dinterpolant *c, ae_state *_state);
void spline1dgriddiffcubic(RVector *x, RVector *y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, RVector *d, ae_state *_state);
void spline1dgriddiff2cubic(RVector *x, RVector *y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, RVector *d1, RVector *d2, ae_state *_state);
void spline1dconvdiffinternal(RVector *xold, RVector *yold, RVector *dold, ae_int_t n, RVector *x2, ae_int_t n2, RVector *y, bool needy, RVector *d1, bool needd1, RVector *d2, bool needd2, ae_state *_state);
void spline1dconvcubic(RVector *x, RVector *y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, RVector *x2, ae_int_t n2, RVector *y2, ae_state *_state);
void spline1dconvdiffcubic(RVector *x, RVector *y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, RVector *x2, ae_int_t n2, RVector *y2, RVector *d2, ae_state *_state);
void spline1dconvdiff2cubic(RVector *x, RVector *y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, RVector *x2, ae_int_t n2, RVector *y2, RVector *d2, RVector *dd2, ae_state *_state);
double spline1dcalc(spline1dinterpolant *c, double x, ae_state *_state);
void spline1ddiff(spline1dinterpolant *c, double x, double *s, double *ds, double *d2s, ae_state *_state);
void spline1dcopy(spline1dinterpolant *c, spline1dinterpolant *cc, ae_state *_state);
void spline1dunpack(spline1dinterpolant *c, ae_int_t *n, RMatrix *tbl, ae_state *_state);
void spline1dlintransx(spline1dinterpolant *c, double a, double b, ae_state *_state);
void spline1dlintransy(spline1dinterpolant *c, double a, double b, ae_state *_state);
double spline1dintegrate(spline1dinterpolant *c, double x, ae_state *_state);
void spline1dfit(RVector *x, RVector *y, ae_int_t n, ae_int_t m, double lambdans, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
void solvepolinom2(double p0, double m0, double p1, double m1, double *x0, double *x1, ae_int_t *nr, ae_state *_state);
ae_int_t bisectmethod(double pa, double ma, double pb, double mb, double a, double b, double *x, ae_state *_state);
void solvecubicpolinom(double pa, double ma, double pb, double mb, double a, double b, double *x0, double *x1, double *x2, double *ex0, double *ex1, ae_int_t *nr, ae_int_t *ne, RVector *tempdata, ae_state *_state);
void spline1drootsandextrema(spline1dinterpolant *c, RVector *r, ae_int_t *nr, bool *dr, RVector *e, ZVector *et, ae_int_t *ne, bool *de, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(spline1dinterpolant, );
DecClass(spline1dfitreport, double &taskrcond; double &rmserror; double &avgerror; double &avgrelerror; double &maxerror;);

void spline1dbuildlinear(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildlinear(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildhermite(const real_1d_array &x, const real_1d_array &y, const real_1d_array &d, const ae_int_t n, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildhermite(const real_1d_array &x, const real_1d_array &y, const real_1d_array &d, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildcubic(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildcatmullrom(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundtype, const double tension, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildcatmullrom(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildakima(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildakima(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildmonotone(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dbuildmonotone(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c, const xparams _xparams = NonTH);
void spline1dgriddiffcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, real_1d_array &d, const xparams _xparams = NonTH);
void spline1dgriddiffcubic(const real_1d_array &x, const real_1d_array &y, real_1d_array &d, const xparams _xparams = NonTH);
void spline1dgriddiff2cubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, real_1d_array &d1, real_1d_array &d2, const xparams _xparams = NonTH);
void spline1dgriddiff2cubic(const real_1d_array &x, const real_1d_array &y, real_1d_array &d1, real_1d_array &d2, const xparams _xparams = NonTH);
void spline1dconvcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y2, const xparams _xparams = NonTH);
void spline1dconvcubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2, const xparams _xparams = NonTH);
void spline1dconvdiffcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y2, real_1d_array &d2, const xparams _xparams = NonTH);
void spline1dconvdiffcubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2, real_1d_array &d2, const xparams _xparams = NonTH);
void spline1dconvdiff2cubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y2, real_1d_array &d2, real_1d_array &dd2, const xparams _xparams = NonTH);
void spline1dconvdiff2cubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2, real_1d_array &d2, real_1d_array &dd2, const xparams _xparams = NonTH);
double spline1dcalc(const spline1dinterpolant &c, const double x, const xparams _xparams = NonTH);
void spline1ddiff(const spline1dinterpolant &c, const double x, double &s, double &ds, double &d2s, const xparams _xparams = NonTH);
void spline1dunpack(const spline1dinterpolant &c, ae_int_t &n, real_2d_array &tbl, const xparams _xparams = NonTH);
void spline1dlintransx(const spline1dinterpolant &c, const double a, const double b, const xparams _xparams = NonTH);
void spline1dlintransy(const spline1dinterpolant &c, const double a, const double b, const xparams _xparams = NonTH);
double spline1dintegrate(const spline1dinterpolant &c, const double x, const xparams _xparams = NonTH);
void spline1dfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, const double lambdans, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, const double lambdans, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === LSFIT Package ===
// Depends on: (Optimization) MINLM
// Depends on: POLINT, SPLINE1D
namespace alglib_impl {
struct polynomialfitreport {
   double taskrcond;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double maxerror;
};
void polynomialfitreport_init(void *_p, ae_state *_state, bool make_automatic);
void polynomialfitreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void polynomialfitreport_free(void *_p, bool make_automatic);

struct barycentricfitreport {
   double taskrcond;
   ae_int_t dbest;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double maxerror;
};
void barycentricfitreport_init(void *_p, ae_state *_state, bool make_automatic);
void barycentricfitreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void barycentricfitreport_free(void *_p, bool make_automatic);

struct lsfitreport {
   double taskrcond;
   ae_int_t iterationscount;
   ae_int_t varidx;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double maxerror;
   double wrmserror;
   ae_matrix covpar;
   ae_vector errpar;
   ae_vector errcurve;
   ae_vector noise;
   double r2;
};
void lsfitreport_init(void *_p, ae_state *_state, bool make_automatic);
void lsfitreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void lsfitreport_free(void *_p, bool make_automatic);

struct lsfitstate {
   ae_int_t optalgo;
   ae_int_t m;
   ae_int_t k;
   double epsx;
   ae_int_t maxits;
   double stpmax;
   bool xrep;
   ae_vector c0;
   ae_vector c1;
   ae_vector s;
   ae_vector bndl;
   ae_vector bndu;
   ae_matrix taskx;
   ae_vector tasky;
   ae_int_t npoints;
   ae_vector taskw;
   ae_int_t nweights;
   ae_int_t wkind;
   ae_int_t wits;
   double diffstep;
   double teststep;
   ae_matrix cleic;
   ae_int_t nec;
   ae_int_t nic;
   bool xupdated;
   bool needf;
   bool needfg;
   bool needfgh;
   ae_int_t pointindex;
   ae_vector x;
   ae_vector c;
   double f;
   ae_vector g;
   ae_matrix h;
   ae_vector wcur;
   ae_vector tmpct;
   ae_vector tmp;
   ae_vector tmpf;
   ae_matrix tmpjac;
   ae_matrix tmpjacw;
   double tmpnoise;
   matinvreport invrep;
   ae_int_t repiterationscount;
   ae_int_t repterminationtype;
   ae_int_t repvaridx;
   double reprmserror;
   double repavgerror;
   double repavgrelerror;
   double repmaxerror;
   double repwrmserror;
   lsfitreport rep;
   minlmstate optstate;
   minlmreport optrep;
   ae_int_t prevnpt;
   ae_int_t prevalgo;
   rcommstate rstate;
};
void lsfitstate_init(void *_p, ae_state *_state, bool make_automatic);
void lsfitstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void lsfitstate_free(void *_p, bool make_automatic);

void lstfitpiecewiselinearrdpfixed(RVector *x, RVector *y, ae_int_t n, ae_int_t m, RVector *x2, RVector *y2, ae_int_t *nsections, ae_state *_state);
void lstfitpiecewiselinearrdp(RVector *x, RVector *y, ae_int_t n, double eps, RVector *x2, RVector *y2, ae_int_t *nsections, ae_state *_state);
void lsfitlinearw(RVector *y, RVector *w, RMatrix *fmatrix, ae_int_t n, ae_int_t m, ae_int_t *info, RVector *c, lsfitreport *rep, ae_state *_state);
void lsfitlinearwc(RVector *y, RVector *w, RMatrix *fmatrix, RMatrix *cmatrix, ae_int_t n, ae_int_t m, ae_int_t k, ae_int_t *info, RVector *c, lsfitreport *rep, ae_state *_state);
void lsfitlinear(RVector *y, RMatrix *fmatrix, ae_int_t n, ae_int_t m, ae_int_t *info, RVector *c, lsfitreport *rep, ae_state *_state);
void lsfitlinearc(RVector *y, RMatrix *fmatrix, RMatrix *cmatrix, ae_int_t n, ae_int_t m, ae_int_t k, ae_int_t *info, RVector *c, lsfitreport *rep, ae_state *_state);
void polynomialfitwc(RVector *x, RVector *y, RVector *w, ae_int_t n, RVector *xc, RVector *yc, ZVector *dc, ae_int_t k, ae_int_t m, ae_int_t *info, barycentricinterpolant *p, polynomialfitreport *rep, ae_state *_state);
void polynomialfit(RVector *x, RVector *y, ae_int_t n, ae_int_t m, ae_int_t *info, barycentricinterpolant *p, polynomialfitreport *rep, ae_state *_state);
double logisticcalc4(double x, double a, double b, double c, double d, ae_state *_state);
double logisticcalc5(double x, double a, double b, double c, double d, double g, ae_state *_state);
void logisticfit45x(RVector *x, RVector *y, ae_int_t n, double cnstrleft, double cnstrright, bool is4pl, double lambdav, double epsx, ae_int_t rscnt, double *a, double *b, double *c, double *d, double *g, lsfitreport *rep, ae_state *_state);
void logisticfit4(RVector *x, RVector *y, ae_int_t n, double *a, double *b, double *c, double *d, lsfitreport *rep, ae_state *_state);
void logisticfit4ec(RVector *x, RVector *y, ae_int_t n, double cnstrleft, double cnstrright, double *a, double *b, double *c, double *d, lsfitreport *rep, ae_state *_state);
void logisticfit5(RVector *x, RVector *y, ae_int_t n, double *a, double *b, double *c, double *d, double *g, lsfitreport *rep, ae_state *_state);
void logisticfit5ec(RVector *x, RVector *y, ae_int_t n, double cnstrleft, double cnstrright, double *a, double *b, double *c, double *d, double *g, lsfitreport *rep, ae_state *_state);
void barycentricfitfloaterhormannwc(RVector *x, RVector *y, RVector *w, ae_int_t n, RVector *xc, RVector *yc, ZVector *dc, ae_int_t k, ae_int_t m, ae_int_t *info, barycentricinterpolant *b, barycentricfitreport *rep, ae_state *_state);
void barycentricfitfloaterhormann(RVector *x, RVector *y, ae_int_t n, ae_int_t m, ae_int_t *info, barycentricinterpolant *b, barycentricfitreport *rep, ae_state *_state);
void spline1dfitcubicwc(RVector *x, RVector *y, RVector *w, ae_int_t n, RVector *xc, RVector *yc, ZVector *dc, ae_int_t k, ae_int_t m, ae_int_t *info, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
void spline1dfitcubic(RVector *x, RVector *y, ae_int_t n, ae_int_t m, ae_int_t *info, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
void spline1dfithermitewc(RVector *x, RVector *y, RVector *w, ae_int_t n, RVector *xc, RVector *yc, ZVector *dc, ae_int_t k, ae_int_t m, ae_int_t *info, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
void spline1dfithermite(RVector *x, RVector *y, ae_int_t n, ae_int_t m, ae_int_t *info, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
void lsfitsetcond(lsfitstate *state, double epsx, ae_int_t maxits, ae_state *_state);
void lsfitsetstpmax(lsfitstate *state, double stpmax, ae_state *_state);
void lsfitsetxrep(lsfitstate *state, bool needxrep, ae_state *_state);
void lsfitsetscale(lsfitstate *state, RVector *s, ae_state *_state);
void lsfitsetbc(lsfitstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void lsfitsetlc(lsfitstate *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void lsfitsetgradientcheck(lsfitstate *state, double teststep, ae_state *_state);
void lsfitcreatewf(RMatrix *x, RVector *y, RVector *w, RVector *c, ae_int_t n, ae_int_t m, ae_int_t k, double diffstep, lsfitstate *state, ae_state *_state);
void lsfitcreatef(RMatrix *x, RVector *y, RVector *c, ae_int_t n, ae_int_t m, ae_int_t k, double diffstep, lsfitstate *state, ae_state *_state);
void lsfitcreatewfg(RMatrix *x, RVector *y, RVector *w, RVector *c, ae_int_t n, ae_int_t m, ae_int_t k, bool cheapfg, lsfitstate *state, ae_state *_state);
void lsfitcreatefg(RMatrix *x, RVector *y, RVector *c, ae_int_t n, ae_int_t m, ae_int_t k, bool cheapfg, lsfitstate *state, ae_state *_state);
void lsfitcreatewfgh(RMatrix *x, RVector *y, RVector *w, RVector *c, ae_int_t n, ae_int_t m, ae_int_t k, lsfitstate *state, ae_state *_state);
void lsfitcreatefgh(RMatrix *x, RVector *y, RVector *c, ae_int_t n, ae_int_t m, ae_int_t k, lsfitstate *state, ae_state *_state);
bool lsfititeration(lsfitstate *state, ae_state *_state);
void lsfitresults(lsfitstate *state, ae_int_t *info, RVector *c, lsfitreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(polynomialfitreport, double &taskrcond; double &rmserror; double &avgerror; double &avgrelerror; double &maxerror;);
DecClass(barycentricfitreport, double &taskrcond; ae_int_t &dbest; double &rmserror; double &avgerror; double &avgrelerror; double &maxerror;);
DecClass(lsfitreport, double &taskrcond; ae_int_t &iterationscount; ae_int_t &varidx; double &rmserror; double &avgerror; double &avgrelerror; double &maxerror; double &wrmserror; real_2d_array covpar; real_1d_array errpar; real_1d_array errcurve; real_1d_array noise; double &r2;);
DecClass(lsfitstate, bool &needf; bool &needfg; bool &needfgh; bool &xupdated; real_1d_array c; double &f; real_1d_array g; real_2d_array h; real_1d_array x;);

void lstfitpiecewiselinearrdpfixed(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, real_1d_array &x2, real_1d_array &y2, ae_int_t &nsections, const xparams _xparams = NonTH);
void lstfitpiecewiselinearrdp(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double eps, real_1d_array &x2, real_1d_array &y2, ae_int_t &nsections, const xparams _xparams = NonTH);
void lsfitlinearw(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const ae_int_t n, const ae_int_t m, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinearw(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinearwc(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const real_2d_array &cmatrix, const ae_int_t n, const ae_int_t m, const ae_int_t k, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinearwc(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const real_2d_array &cmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinear(const real_1d_array &y, const real_2d_array &fmatrix, const ae_int_t n, const ae_int_t m, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinear(const real_1d_array &y, const real_2d_array &fmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinearc(const real_1d_array &y, const real_2d_array &fmatrix, const real_2d_array &cmatrix, const ae_int_t n, const ae_int_t m, const ae_int_t k, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void lsfitlinearc(const real_1d_array &y, const real_2d_array &fmatrix, const real_2d_array &cmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
void polynomialfitwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep, const xparams _xparams = NonTH);
void polynomialfitwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep, const xparams _xparams = NonTH);
void polynomialfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep, const xparams _xparams = NonTH);
void polynomialfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep, const xparams _xparams = NonTH);
double logisticcalc4(const double x, const double a, const double b, const double c, const double d, const xparams _xparams = NonTH);
double logisticcalc5(const double x, const double a, const double b, const double c, const double d, const double g, const xparams _xparams = NonTH);
void logisticfit45x(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double cnstrleft, const double cnstrright, const bool is4pl, const double lambdav, const double epsx, const ae_int_t rscnt, double &a, double &b, double &c, double &d, double &g, lsfitreport &rep, const xparams _xparams = NonTH);
void logisticfit4(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, double &a, double &b, double &c, double &d, lsfitreport &rep, const xparams _xparams = NonTH);
void logisticfit4ec(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double cnstrleft, const double cnstrright, double &a, double &b, double &c, double &d, lsfitreport &rep, const xparams _xparams = NonTH);
void logisticfit5(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, double &a, double &b, double &c, double &d, double &g, lsfitreport &rep, const xparams _xparams = NonTH);
void logisticfit5ec(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double cnstrleft, const double cnstrright, double &a, double &b, double &c, double &d, double &g, lsfitreport &rep, const xparams _xparams = NonTH);
void barycentricfitfloaterhormannwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, barycentricinterpolant &b, barycentricfitreport &rep, const xparams _xparams = NonTH);
void barycentricfitfloaterhormann(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, barycentricinterpolant &b, barycentricfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitcubicwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitcubicwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfithermitewc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfithermitewc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfithermite(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfithermite(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void lsfitsetcond(const lsfitstate &state, const double epsx, const ae_int_t maxits, const xparams _xparams = NonTH);
void lsfitsetstpmax(const lsfitstate &state, const double stpmax, const xparams _xparams = NonTH);
void lsfitsetxrep(const lsfitstate &state, const bool needxrep, const xparams _xparams = NonTH);
void lsfitsetscale(const lsfitstate &state, const real_1d_array &s, const xparams _xparams = NonTH);
void lsfitsetbc(const lsfitstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = NonTH);
void lsfitsetlc(const lsfitstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = NonTH);
void lsfitsetlc(const lsfitstate &state, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = NonTH);
void lsfitsetgradientcheck(const lsfitstate &state, const double teststep, const xparams _xparams = NonTH);
void lsfitcreatewf(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const double diffstep, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatewf(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const double diffstep, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatef(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const double diffstep, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatef(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const double diffstep, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatewfg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const bool cheapfg, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatewfg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const bool cheapfg, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatefg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const bool cheapfg, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatefg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const bool cheapfg, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatewfgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatewfgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatefgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, lsfitstate &state, const xparams _xparams = NonTH);
void lsfitcreatefgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, lsfitstate &state, const xparams _xparams = NonTH);
bool lsfititeration(const lsfitstate &state, const xparams _xparams = NonTH);
void lsfitfit(lsfitstate &state, void (*func)(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr), void (*rep)(const real_1d_array &c, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = NonTH);
void lsfitfit(lsfitstate &state, void (*func)(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr), void (*grad)(const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*rep)(const real_1d_array &c, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = NonTH);
void lsfitfit(lsfitstate &state, void (*func)(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr), void (*grad)(const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*hess)(const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr), void (*rep)(const real_1d_array &c, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = NonTH);
void lsfitresults(const lsfitstate &state, ae_int_t &info, real_1d_array &c, lsfitreport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === FITSPHERE Package ===
// Depends on: (Optimization) MINLM, MINNLC
namespace alglib_impl {
struct fitsphereinternalreport {
   ae_int_t nfev;
   ae_int_t iterationscount;
};
void fitsphereinternalreport_init(void *_p, ae_state *_state, bool make_automatic);
void fitsphereinternalreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void fitsphereinternalreport_free(void *_p, bool make_automatic);

void fitsphereinternal(RMatrix *xy, ae_int_t npoints, ae_int_t nx, ae_int_t problemtype, ae_int_t solvertype, double epsx, ae_int_t aulits, double penalty, RVector *cx, double *rlo, double *rhi, fitsphereinternalreport *rep, ae_state *_state);
void fitspherex(RMatrix *xy, ae_int_t npoints, ae_int_t nx, ae_int_t problemtype, double epsx, ae_int_t aulits, double penalty, RVector *cx, double *rlo, double *rhi, ae_state *_state);
void fitspherels(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *r, ae_state *_state);
void fitspheremc(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *rhi, ae_state *_state);
void fitspheremi(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *rlo, ae_state *_state);
void fitspheremz(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *rlo, double *rhi, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void fitspherex(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, const ae_int_t problemtype, const double epsx, const ae_int_t aulits, const double penalty, real_1d_array &cx, double &rlo, double &rhi, const xparams _xparams = NonTH);
void fitspherels(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &r, const xparams _xparams = NonTH);
void fitspheremc(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &rhi, const xparams _xparams = NonTH);
void fitspheremi(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &rlo, const xparams _xparams = NonTH);
void fitspheremz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &rlo, double &rhi, const xparams _xparams = NonTH);
} // end of namespace alglib

// === PARAMETRIC Package ===
// Depends on: (Integration) AUTOGK
// Depends on: SPLINE1D
namespace alglib_impl {
struct pspline2interpolant {
   ae_int_t n;
   bool periodic;
   ae_vector p;
   spline1dinterpolant x;
   spline1dinterpolant y;
};
void pspline2interpolant_init(void *_p, ae_state *_state, bool make_automatic);
void pspline2interpolant_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void pspline2interpolant_free(void *_p, bool make_automatic);

struct pspline3interpolant {
   ae_int_t n;
   bool periodic;
   ae_vector p;
   spline1dinterpolant x;
   spline1dinterpolant y;
   spline1dinterpolant z;
};
void pspline3interpolant_init(void *_p, ae_state *_state, bool make_automatic);
void pspline3interpolant_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void pspline3interpolant_free(void *_p, bool make_automatic);

void pspline2build(RMatrix *xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline2interpolant *p, ae_state *_state);
void pspline3build(RMatrix *xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline3interpolant *p, ae_state *_state);
void pspline2buildperiodic(RMatrix *xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline2interpolant *p, ae_state *_state);
void pspline3buildperiodic(RMatrix *xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline3interpolant *p, ae_state *_state);
void pspline2parametervalues(pspline2interpolant *p, ae_int_t *n, RVector *t, ae_state *_state);
void pspline3parametervalues(pspline3interpolant *p, ae_int_t *n, RVector *t, ae_state *_state);
void pspline2calc(pspline2interpolant *p, double t, double *x, double *y, ae_state *_state);
void pspline3calc(pspline3interpolant *p, double t, double *x, double *y, double *z, ae_state *_state);
void pspline2diff(pspline2interpolant *p, double t, double *x, double *dx, double *y, double *dy, ae_state *_state);
void pspline3diff(pspline3interpolant *p, double t, double *x, double *dx, double *y, double *dy, double *z, double *dz, ae_state *_state);
void pspline2tangent(pspline2interpolant *p, double t, double *x, double *y, ae_state *_state);
void pspline3tangent(pspline3interpolant *p, double t, double *x, double *y, double *z, ae_state *_state);
void pspline2diff2(pspline2interpolant *p, double t, double *x, double *dx, double *d2x, double *y, double *dy, double *d2y, ae_state *_state);
void pspline3diff2(pspline3interpolant *p, double t, double *x, double *dx, double *d2x, double *y, double *dy, double *d2y, double *z, double *dz, double *d2z, ae_state *_state);
double pspline2arclength(pspline2interpolant *p, double a, double b, ae_state *_state);
double pspline3arclength(pspline3interpolant *p, double a, double b, ae_state *_state);
void parametricrdpfixed(RMatrix *x, ae_int_t n, ae_int_t d, ae_int_t stopm, double stopeps, RMatrix *x2, ZVector *idx2, ae_int_t *nsections, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(pspline2interpolant, );
DecClass(pspline3interpolant, );

void pspline2build(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline2interpolant &p, const xparams _xparams = NonTH);
void pspline3build(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline3interpolant &p, const xparams _xparams = NonTH);
void pspline2buildperiodic(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline2interpolant &p, const xparams _xparams = NonTH);
void pspline3buildperiodic(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline3interpolant &p, const xparams _xparams = NonTH);
void pspline2parametervalues(const pspline2interpolant &p, ae_int_t &n, real_1d_array &t, const xparams _xparams = NonTH);
void pspline3parametervalues(const pspline3interpolant &p, ae_int_t &n, real_1d_array &t, const xparams _xparams = NonTH);
void pspline2calc(const pspline2interpolant &p, const double t, double &x, double &y, const xparams _xparams = NonTH);
void pspline3calc(const pspline3interpolant &p, const double t, double &x, double &y, double &z, const xparams _xparams = NonTH);
void pspline2diff(const pspline2interpolant &p, const double t, double &x, double &dx, double &y, double &dy, const xparams _xparams = NonTH);
void pspline3diff(const pspline3interpolant &p, const double t, double &x, double &dx, double &y, double &dy, double &z, double &dz, const xparams _xparams = NonTH);
void pspline2tangent(const pspline2interpolant &p, const double t, double &x, double &y, const xparams _xparams = NonTH);
void pspline3tangent(const pspline3interpolant &p, const double t, double &x, double &y, double &z, const xparams _xparams = NonTH);
void pspline2diff2(const pspline2interpolant &p, const double t, double &x, double &dx, double &d2x, double &y, double &dy, double &d2y, const xparams _xparams = NonTH);
void pspline3diff2(const pspline3interpolant &p, const double t, double &x, double &dx, double &d2x, double &y, double &dy, double &d2y, double &z, double &dz, double &d2z, const xparams _xparams = NonTH);
double pspline2arclength(const pspline2interpolant &p, const double a, const double b, const xparams _xparams = NonTH);
double pspline3arclength(const pspline3interpolant &p, const double a, const double b, const xparams _xparams = NonTH);
void parametricrdpfixed(const real_2d_array &x, const ae_int_t n, const ae_int_t d, const ae_int_t stopm, const double stopeps, real_2d_array &x2, integer_1d_array &idx2, ae_int_t &nsections, const xparams _xparams = NonTH);
} // end of namespace alglib

// === RBFV1 Package ===
// Depends on: (AlgLibMisc) NEARESTNEIGHBOR
// Depends on: LSFIT
namespace alglib_impl {
struct rbfv1calcbuffer {
   ae_vector calcbufxcx;
   ae_matrix calcbufx;
   ae_vector calcbuftags;
   kdtreerequestbuffer requestbuffer;
};
void rbfv1calcbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv1calcbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv1calcbuffer_free(void *_p, bool make_automatic);

struct rbfv1model {
   ae_int_t ny;
   ae_int_t nx;
   ae_int_t nc;
   ae_int_t nl;
   kdtree tree;
   ae_matrix xc;
   ae_matrix wr;
   double rmax;
   ae_matrix v;
   ae_vector calcbufxcx;
   ae_matrix calcbufx;
   ae_vector calcbuftags;
};
void rbfv1model_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv1model_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv1model_free(void *_p, bool make_automatic);
void rbfv1alloc(ae_serializer *s, rbfv1model *model, ae_state *_state);
void rbfv1serialize(ae_serializer *s, rbfv1model *model, ae_state *_state);
void rbfv1unserialize(ae_serializer *s, rbfv1model *model, ae_state *_state);

struct gridcalc3v1buf {
   ae_vector tx;
   ae_vector cx;
   ae_vector ty;
   ae_vector flag0;
   ae_vector flag1;
   ae_vector flag2;
   ae_vector flag12;
   ae_vector expbuf0;
   ae_vector expbuf1;
   ae_vector expbuf2;
   kdtreerequestbuffer requestbuf;
   ae_matrix calcbufx;
   ae_vector calcbuftags;
};
void gridcalc3v1buf_init(void *_p, ae_state *_state, bool make_automatic);
void gridcalc3v1buf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void gridcalc3v1buf_free(void *_p, bool make_automatic);

struct rbfv1report {
   ae_int_t arows;
   ae_int_t acols;
   ae_int_t annz;
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
};
void rbfv1report_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv1report_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv1report_free(void *_p, bool make_automatic);

void rbfv1create(ae_int_t nx, ae_int_t ny, rbfv1model *s, ae_state *_state);
void rbfv1createcalcbuffer(rbfv1model *s, rbfv1calcbuffer *buf, ae_state *_state);
void rbfv1buildmodel(RMatrix *x, RMatrix *y, ae_int_t n, ae_int_t aterm, ae_int_t algorithmtype, ae_int_t nlayers, double radvalue, double radzvalue, double lambdav, double epsort, double epserr, ae_int_t maxits, rbfv1model *s, rbfv1report *rep, ae_state *_state);
double rbfv1calc2(rbfv1model *s, double x0, double x1, ae_state *_state);
double rbfv1calc3(rbfv1model *s, double x0, double x1, double x2, ae_state *_state);
void rbfv1calcbuf(rbfv1model *s, RVector *x, RVector *y, ae_state *_state);
void rbfv1tscalcbuf(rbfv1model *s, rbfv1calcbuffer *buf, RVector *x, RVector *y, ae_state *_state);
void rbfv1gridcalc2(rbfv1model *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RMatrix *y, ae_state *_state);
void rbfv1gridcalc3vrec(rbfv1model *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *x2, ae_int_t n2, ZVector *blocks0, ae_int_t block0a, ae_int_t block0b, ZVector *blocks1, ae_int_t block1a, ae_int_t block1b, ZVector *blocks2, ae_int_t block2a, ae_int_t block2b, BVector *flagy, bool sparsey, double searchradius, double avgfuncpernode, ae_shared_pool *bufpool, RVector *y, ae_state *_state);
void rbfv1unpack(rbfv1model *s, ae_int_t *nx, ae_int_t *ny, RMatrix *xwr, ae_int_t *nc, RMatrix *v, ae_state *_state);
} // end of namespace alglib_impl

// === SPLINE2D Package ===
// Depends on: SPLINE1D
namespace alglib_impl {
struct spline2dinterpolant {
   ae_int_t stype;
   ae_int_t n;
   ae_int_t m;
   ae_int_t d;
   ae_vector x;
   ae_vector y;
   ae_vector f;
};
void spline2dinterpolant_init(void *_p, ae_state *_state, bool make_automatic);
void spline2dinterpolant_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline2dinterpolant_free(void *_p, bool make_automatic);
void spline2dalloc(ae_serializer *s, spline2dinterpolant *spline, ae_state *_state);
void spline2dserialize(ae_serializer *s, spline2dinterpolant *spline, ae_state *_state);
void spline2dunserialize(ae_serializer *s, spline2dinterpolant *spline, ae_state *_state);

struct spline2dbuilder {
   ae_int_t priorterm;
   double priortermval;
   ae_int_t areatype;
   double xa;
   double xb;
   double ya;
   double yb;
   ae_int_t gridtype;
   ae_int_t kx;
   ae_int_t ky;
   double smoothing;
   ae_int_t nlayers;
   ae_int_t solvertype;
   double lambdabase;
   ae_vector xy;
   ae_int_t npoints;
   ae_int_t d;
   double sx;
   double sy;
   bool adddegreeoffreedom;
   ae_int_t interfacesize;
   ae_int_t lsqrcnt;
   ae_int_t maxcoresize;
};
void spline2dbuilder_init(void *_p, ae_state *_state, bool make_automatic);
void spline2dbuilder_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline2dbuilder_free(void *_p, bool make_automatic);

struct spline2dfitreport {
   double rmserror;
   double avgerror;
   double maxerror;
   double r2;
};
void spline2dfitreport_init(void *_p, ae_state *_state, bool make_automatic);
void spline2dfitreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline2dfitreport_free(void *_p, bool make_automatic);

struct spline2dxdesignmatrix {
   ae_int_t blockwidth;
   ae_int_t kx;
   ae_int_t ky;
   ae_int_t npoints;
   ae_int_t nrows;
   ae_int_t ndenserows;
   ae_int_t ndensebatches;
   ae_int_t d;
   ae_int_t maxbatch;
   ae_matrix vals;
   ae_vector batches;
   ae_vector batchbases;
   double lambdareg;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_matrix tmp2;
};
void spline2dxdesignmatrix_init(void *_p, ae_state *_state, bool make_automatic);
void spline2dxdesignmatrix_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline2dxdesignmatrix_free(void *_p, bool make_automatic);

struct spline2dblockllsbuf {
   linlsqrstate solver;
   linlsqrreport solverrep;
   ae_matrix blockata;
   ae_matrix trsmbuf2;
   ae_matrix cholbuf2;
   ae_vector cholbuf1;
   ae_vector tmp0;
   ae_vector tmp1;
};
void spline2dblockllsbuf_init(void *_p, ae_state *_state, bool make_automatic);
void spline2dblockllsbuf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline2dblockllsbuf_free(void *_p, bool make_automatic);

struct spline2dfastddmbuf {
   spline2dxdesignmatrix xdesignmatrix;
   ae_vector tmp0;
   ae_vector tmpz;
   spline2dfitreport dummyrep;
   spline2dinterpolant localmodel;
   spline2dblockllsbuf blockllsbuf;
};
void spline2dfastddmbuf_init(void *_p, ae_state *_state, bool make_automatic);
void spline2dfastddmbuf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline2dfastddmbuf_free(void *_p, bool make_automatic);

double spline2dcalc(spline2dinterpolant *c, double x, double y, ae_state *_state);
void spline2dcalcvbuf(spline2dinterpolant *c, double x, double y, RVector *f, ae_state *_state);
void spline2dcalcv(spline2dinterpolant *c, double x, double y, RVector *f, ae_state *_state);
double spline2dcalcvi(spline2dinterpolant *c, double x, double y, ae_int_t i, ae_state *_state);
void spline2ddiff(spline2dinterpolant *c, double x, double y, double *f, double *fx, double *fy, double *fxy, ae_state *_state);
void spline2ddiffvi(spline2dinterpolant *c, double x, double y, ae_int_t i, double *f, double *fx, double *fy, double *fxy, ae_state *_state);
void spline2dcopy(spline2dinterpolant *c, spline2dinterpolant *cc, ae_state *_state);
void spline2dresamplebicubic(RMatrix *a, ae_int_t oldheight, ae_int_t oldwidth, RMatrix *b, ae_int_t newheight, ae_int_t newwidth, ae_state *_state);
void spline2dresamplebilinear(RMatrix *a, ae_int_t oldheight, ae_int_t oldwidth, RMatrix *b, ae_int_t newheight, ae_int_t newwidth, ae_state *_state);
void spline2dbuildbilinearv(RVector *x, ae_int_t n, RVector *y, ae_int_t m, RVector *f, ae_int_t d, spline2dinterpolant *c, ae_state *_state);
void spline2dbuildbicubicv(RVector *x, ae_int_t n, RVector *y, ae_int_t m, RVector *f, ae_int_t d, spline2dinterpolant *c, ae_state *_state);
void spline2dbuildbilinear(RVector *x, RVector *y, RMatrix *f, ae_int_t m, ae_int_t n, spline2dinterpolant *c, ae_state *_state);
void spline2dbuildbicubic(RVector *x, RVector *y, RMatrix *f, ae_int_t m, ae_int_t n, spline2dinterpolant *c, ae_state *_state);
void spline2dlintransxy(spline2dinterpolant *c, double ax, double bx, double ay, double by, ae_state *_state);
void spline2dlintransf(spline2dinterpolant *c, double a, double b, ae_state *_state);
void spline2dunpackv(spline2dinterpolant *c, ae_int_t *m, ae_int_t *n, ae_int_t *d, RMatrix *tbl, ae_state *_state);
void spline2dunpack(spline2dinterpolant *c, ae_int_t *m, ae_int_t *n, RMatrix *tbl, ae_state *_state);
void spline2dbuildercreate(ae_int_t d, spline2dbuilder *state, ae_state *_state);
void spline2dbuildersetuserterm(spline2dbuilder *state, double v, ae_state *_state);
void spline2dbuildersetlinterm(spline2dbuilder *state, ae_state *_state);
void spline2dbuildersetconstterm(spline2dbuilder *state, ae_state *_state);
void spline2dbuildersetzeroterm(spline2dbuilder *state, ae_state *_state);
void spline2dbuildersetpoints(spline2dbuilder *state, RMatrix *xy, ae_int_t n, ae_state *_state);
void spline2dbuildersetareaauto(spline2dbuilder *state, ae_state *_state);
void spline2dbuildersetarea(spline2dbuilder *state, double xa, double xb, double ya, double yb, ae_state *_state);
void spline2dbuildersetgrid(spline2dbuilder *state, ae_int_t kx, ae_int_t ky, ae_state *_state);
void spline2dbuildersetalgofastddm(spline2dbuilder *state, ae_int_t nlayers, double lambdav, ae_state *_state);
void spline2dbuildersetalgoblocklls(spline2dbuilder *state, double lambdans, ae_state *_state);
void spline2dbuildersetalgonaivells(spline2dbuilder *state, double lambdans, ae_state *_state);
void spline2dfit(spline2dbuilder *state, spline2dinterpolant *s, spline2dfitreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(spline2dinterpolant, );
DecClass(spline2dbuilder, );
DecClass(spline2dfitreport, double &rmserror; double &avgerror; double &maxerror; double &r2;);
void spline2dserialize(spline2dinterpolant &obj, std::string &s_out);
void spline2dserialize(spline2dinterpolant &obj, std::ostream &s_out);
void spline2dunserialize(const std::string &s_in, spline2dinterpolant &obj);
void spline2dunserialize(const std::istream &s_in, spline2dinterpolant &obj);

double spline2dcalc(const spline2dinterpolant &c, const double x, const double y, const xparams _xparams = NonTH);
void spline2dcalcvbuf(const spline2dinterpolant &c, const double x, const double y, real_1d_array &f, const xparams _xparams = NonTH);
void spline2dcalcv(const spline2dinterpolant &c, const double x, const double y, real_1d_array &f, const xparams _xparams = NonTH);
double spline2dcalcvi(const spline2dinterpolant &c, const double x, const double y, const ae_int_t i, const xparams _xparams = NonTH);
void spline2ddiff(const spline2dinterpolant &c, const double x, const double y, double &f, double &fx, double &fy, double &fxy, const xparams _xparams = NonTH);
void spline2ddiffvi(const spline2dinterpolant &c, const double x, const double y, const ae_int_t i, double &f, double &fx, double &fy, double &fxy, const xparams _xparams = NonTH);
void spline2dcopy(const spline2dinterpolant &c, spline2dinterpolant &cc, const xparams _xparams = NonTH);
void spline2dresamplebicubic(const real_2d_array &a, const ae_int_t oldheight, const ae_int_t oldwidth, real_2d_array &b, const ae_int_t newheight, const ae_int_t newwidth, const xparams _xparams = NonTH);
void spline2dresamplebilinear(const real_2d_array &a, const ae_int_t oldheight, const ae_int_t oldwidth, real_2d_array &b, const ae_int_t newheight, const ae_int_t newwidth, const xparams _xparams = NonTH);
void spline2dbuildbilinearv(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, const real_1d_array &f, const ae_int_t d, spline2dinterpolant &c, const xparams _xparams = NonTH);
void spline2dbuildbicubicv(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, const real_1d_array &f, const ae_int_t d, spline2dinterpolant &c, const xparams _xparams = NonTH);
void spline2dbuildbilinear(const real_1d_array &x, const real_1d_array &y, const real_2d_array &f, const ae_int_t m, const ae_int_t n, spline2dinterpolant &c, const xparams _xparams = NonTH);
void spline2dbuildbicubic(const real_1d_array &x, const real_1d_array &y, const real_2d_array &f, const ae_int_t m, const ae_int_t n, spline2dinterpolant &c, const xparams _xparams = NonTH);
void spline2dlintransxy(const spline2dinterpolant &c, const double ax, const double bx, const double ay, const double by, const xparams _xparams = NonTH);
void spline2dlintransf(const spline2dinterpolant &c, const double a, const double b, const xparams _xparams = NonTH);
void spline2dunpackv(const spline2dinterpolant &c, ae_int_t &m, ae_int_t &n, ae_int_t &d, real_2d_array &tbl, const xparams _xparams = NonTH);
void spline2dunpack(const spline2dinterpolant &c, ae_int_t &m, ae_int_t &n, real_2d_array &tbl, const xparams _xparams = NonTH);
void spline2dbuildercreate(const ae_int_t d, spline2dbuilder &state, const xparams _xparams = NonTH);
void spline2dbuildersetuserterm(const spline2dbuilder &state, const double v, const xparams _xparams = NonTH);
void spline2dbuildersetlinterm(const spline2dbuilder &state, const xparams _xparams = NonTH);
void spline2dbuildersetconstterm(const spline2dbuilder &state, const xparams _xparams = NonTH);
void spline2dbuildersetzeroterm(const spline2dbuilder &state, const xparams _xparams = NonTH);
void spline2dbuildersetpoints(const spline2dbuilder &state, const real_2d_array &xy, const ae_int_t n, const xparams _xparams = NonTH);
void spline2dbuildersetareaauto(const spline2dbuilder &state, const xparams _xparams = NonTH);
void spline2dbuildersetarea(const spline2dbuilder &state, const double xa, const double xb, const double ya, const double yb, const xparams _xparams = NonTH);
void spline2dbuildersetgrid(const spline2dbuilder &state, const ae_int_t kx, const ae_int_t ky, const xparams _xparams = NonTH);
void spline2dbuildersetalgofastddm(const spline2dbuilder &state, const ae_int_t nlayers, const double lambdav, const xparams _xparams = NonTH);
void spline2dbuildersetalgoblocklls(const spline2dbuilder &state, const double lambdans, const xparams _xparams = NonTH);
void spline2dbuildersetalgonaivells(const spline2dbuilder &state, const double lambdans, const xparams _xparams = NonTH);
void spline2dfit(const spline2dbuilder &state, spline2dinterpolant &s, spline2dfitreport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === RBFV2 Package ===
// Depends on: (AlgLibMisc) NEARESTNEIGHBOR
// Depends on: LSFIT
namespace alglib_impl {
struct rbfv2calcbuffer {
   ae_vector x;
   ae_vector curboxmin;
   ae_vector curboxmax;
   double curdist2;
   ae_vector x123;
   ae_vector y123;
};
void rbfv2calcbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv2calcbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv2calcbuffer_free(void *_p, bool make_automatic);

struct rbfv2model {
   ae_int_t ny;
   ae_int_t nx;
   ae_int_t bf;
   ae_int_t nh;
   ae_vector ri;
   ae_vector s;
   ae_vector kdroots;
   ae_vector kdnodes;
   ae_vector kdsplits;
   ae_vector kdboxmin;
   ae_vector kdboxmax;
   ae_vector cw;
   ae_matrix v;
   double lambdareg;
   ae_int_t maxits;
   double supportr;
   ae_int_t basisfunction;
   rbfv2calcbuffer calcbuf;
};
void rbfv2model_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv2model_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv2model_free(void *_p, bool make_automatic);
void rbfv2alloc(ae_serializer *s, rbfv2model *model, ae_state *_state);
void rbfv2serialize(ae_serializer *s, rbfv2model *model, ae_state *_state);
void rbfv2unserialize(ae_serializer *s, rbfv2model *model, ae_state *_state);

struct rbfv2gridcalcbuffer {
   rbfv2calcbuffer calcbuf;
   ae_vector cx;
   ae_vector rx;
   ae_vector ry;
   ae_vector tx;
   ae_vector ty;
   ae_vector rf;
};
void rbfv2gridcalcbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv2gridcalcbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv2gridcalcbuffer_free(void *_p, bool make_automatic);

struct rbfv2report {
   ae_int_t terminationtype;
   double maxerror;
   double rmserror;
};
void rbfv2report_init(void *_p, ae_state *_state, bool make_automatic);
void rbfv2report_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfv2report_free(void *_p, bool make_automatic);

void rbfv2create(ae_int_t nx, ae_int_t ny, rbfv2model *s, ae_state *_state);
void rbfv2createcalcbuffer(rbfv2model *s, rbfv2calcbuffer *buf, ae_state *_state);
double rbfv2farradius(ae_int_t bf, ae_state *_state);
double rbfv2nearradius(ae_int_t bf, ae_state *_state);
double rbfv2basisfunc(ae_int_t bf, double d2, ae_state *_state);
void rbfv2basisfuncdiff2(ae_int_t bf, double d2, double *f, double *df, double *d2f, ae_state *_state);
void rbfv2buildhierarchical(RMatrix *x, RMatrix *y, ae_int_t n, RVector *scalevec, ae_int_t aterm, ae_int_t nh, double rbase, double lambdans, rbfv2model *s, ae_int_t *progress10000, bool *terminationrequest, rbfv2report *rep, ae_state *_state);
void rbfv2tscalcbuf(rbfv2model *s, rbfv2calcbuffer *buf, RVector *x, RVector *y, ae_state *_state);
void rbfv2calcbuf(rbfv2model *s, RVector *x, RVector *y, ae_state *_state);
double rbfv2calc1(rbfv2model *s, double x0, ae_state *_state);
double rbfv2calc2(rbfv2model *s, double x0, double x1, ae_state *_state);
double rbfv2calc3(rbfv2model *s, double x0, double x1, double x2, ae_state *_state);
void rbfv2partialgridcalcrec(rbfv2model *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *x2, ae_int_t n2, RVector *x3, ae_int_t n3, ZVector *blocks0, ae_int_t block0a, ae_int_t block0b, ZVector *blocks1, ae_int_t block1a, ae_int_t block1b, ZVector *blocks2, ae_int_t block2a, ae_int_t block2b, ZVector *blocks3, ae_int_t block3a, ae_int_t block3b, BVector *flagy, bool sparsey, ae_int_t levelidx, double avgfuncpernode, ae_shared_pool *bufpool, RVector *y, ae_state *_state);
void rbfv2gridcalcvx(rbfv2model *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *x2, ae_int_t n2, RVector *x3, ae_int_t n3, BVector *flagy, bool sparsey, RVector *y, ae_state *_state);
void rbfv2gridcalc2(rbfv2model *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RMatrix *y, ae_state *_state);
void rbfv2unpack(rbfv2model *s, ae_int_t *nx, ae_int_t *ny, RMatrix *xwr, ae_int_t *nc, RMatrix *v, ae_state *_state);
} // end of namespace alglib_impl

// === SPLINE3D Package ===
// Depends on: SPLINE1D
namespace alglib_impl {
struct spline3dinterpolant {
   ae_int_t k;
   ae_int_t stype;
   ae_int_t n;
   ae_int_t m;
   ae_int_t l;
   ae_int_t d;
   ae_vector x;
   ae_vector y;
   ae_vector z;
   ae_vector f;
};
void spline3dinterpolant_init(void *_p, ae_state *_state, bool make_automatic);
void spline3dinterpolant_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void spline3dinterpolant_free(void *_p, bool make_automatic);

void spline3dcalcvbuf(spline3dinterpolant *c, double x, double y, double z, RVector *f, ae_state *_state);
void spline3dcalcv(spline3dinterpolant *c, double x, double y, double z, RVector *f, ae_state *_state);
double spline3dcalc(spline3dinterpolant *c, double x, double y, double z, ae_state *_state);
void spline3dresampletrilinear(RVector *a, ae_int_t oldzcount, ae_int_t oldycount, ae_int_t oldxcount, ae_int_t newzcount, ae_int_t newycount, ae_int_t newxcount, RVector *b, ae_state *_state);
void spline3dbuildtrilinearv(RVector *x, ae_int_t n, RVector *y, ae_int_t m, RVector *z, ae_int_t l, RVector *f, ae_int_t d, spline3dinterpolant *c, ae_state *_state);
void spline3dlintransxyz(spline3dinterpolant *c, double ax, double bx, double ay, double by, double az, double bz, ae_state *_state);
void spline3dlintransf(spline3dinterpolant *c, double a, double b, ae_state *_state);
void spline3dcopy(spline3dinterpolant *c, spline3dinterpolant *cc, ae_state *_state);
void spline3dunpackv(spline3dinterpolant *c, ae_int_t *n, ae_int_t *m, ae_int_t *l, ae_int_t *d, ae_int_t *stype, RMatrix *tbl, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(spline3dinterpolant, );

void spline3dcalcvbuf(const spline3dinterpolant &c, const double x, const double y, const double z, real_1d_array &f, const xparams _xparams = NonTH);
void spline3dcalcv(const spline3dinterpolant &c, const double x, const double y, const double z, real_1d_array &f, const xparams _xparams = NonTH);
double spline3dcalc(const spline3dinterpolant &c, const double x, const double y, const double z, const xparams _xparams = NonTH);
void spline3dresampletrilinear(const real_1d_array &a, const ae_int_t oldzcount, const ae_int_t oldycount, const ae_int_t oldxcount, const ae_int_t newzcount, const ae_int_t newycount, const ae_int_t newxcount, real_1d_array &b, const xparams _xparams = NonTH);
void spline3dbuildtrilinearv(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, const real_1d_array &z, const ae_int_t l, const real_1d_array &f, const ae_int_t d, spline3dinterpolant &c, const xparams _xparams = NonTH);
void spline3dlintransxyz(const spline3dinterpolant &c, const double ax, const double bx, const double ay, const double by, const double az, const double bz, const xparams _xparams = NonTH);
void spline3dlintransf(const spline3dinterpolant &c, const double a, const double b, const xparams _xparams = NonTH);
void spline3dunpackv(const spline3dinterpolant &c, ae_int_t &n, ae_int_t &m, ae_int_t &l, ae_int_t &d, ae_int_t &stype, real_2d_array &tbl, const xparams _xparams = NonTH);
} // end of namespace alglib

// === INTCOMP Package ===
// Depends on: SPLINE1D, FITSPHERE
namespace alglib_impl {
void nsfitspherex(RMatrix *xy, ae_int_t npoints, ae_int_t nx, ae_int_t problemtype, double epsx, ae_int_t aulits, double penalty, RVector *cx, double *rlo, double *rhi, ae_state *_state);
void nsfitspheremcc(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *rhi, ae_state *_state);
void nsfitspheremic(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *rlo, ae_state *_state);
void nsfitspheremzc(RMatrix *xy, ae_int_t npoints, ae_int_t nx, RVector *cx, double *rlo, double *rhi, ae_state *_state);
void spline1dfitpenalizedw(RVector *x, RVector *y, RVector *w, ae_int_t n, ae_int_t m, double rho, ae_int_t *info, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
void spline1dfitpenalized(RVector *x, RVector *y, ae_int_t n, ae_int_t m, double rho, ae_int_t *info, spline1dinterpolant *s, spline1dfitreport *rep, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void nsfitspherex(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, const ae_int_t problemtype, const double epsx, const ae_int_t aulits, const double penalty, real_1d_array &cx, double &rlo, double &rhi, const xparams _xparams = NonTH);
void nsfitspheremcc(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &rhi, const xparams _xparams = NonTH);
void nsfitspheremic(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &rlo, const xparams _xparams = NonTH);
void nsfitspheremzc(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nx, real_1d_array &cx, double &rlo, double &rhi, const xparams _xparams = NonTH);
void spline1dfitpenalizedw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitpenalizedw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitpenalized(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
void spline1dfitpenalized(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep, const xparams _xparams = NonTH);
} // end of namespace alglib

// === RBF Package ===
// Depends on: RBFV1, RBFV2
namespace alglib_impl {
struct rbfcalcbuffer {
   ae_int_t modelversion;
   rbfv1calcbuffer bufv1;
   rbfv2calcbuffer bufv2;
};
void rbfcalcbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void rbfcalcbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfcalcbuffer_free(void *_p, bool make_automatic);

struct rbfmodel {
   ae_int_t nx;
   ae_int_t ny;
   ae_int_t modelversion;
   rbfv1model model1;
   rbfv2model model2;
   double lambdav;
   double radvalue;
   double radzvalue;
   ae_int_t nlayers;
   ae_int_t aterm;
   ae_int_t algorithmtype;
   double epsort;
   double epserr;
   ae_int_t maxits;
   ae_int_t nnmaxits;
   ae_int_t n;
   ae_matrix x;
   ae_matrix y;
   bool hasscale;
   ae_vector s;
   ae_int_t progress10000;
   bool terminationrequest;
};
void rbfmodel_init(void *_p, ae_state *_state, bool make_automatic);
void rbfmodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfmodel_free(void *_p, bool make_automatic);
void rbfalloc(ae_serializer *s, rbfmodel *model, ae_state *_state);
void rbfserialize(ae_serializer *s, rbfmodel *model, ae_state *_state);
void rbfunserialize(ae_serializer *s, rbfmodel *model, ae_state *_state);

struct rbfreport {
   double rmserror;
   double maxerror;
   ae_int_t arows;
   ae_int_t acols;
   ae_int_t annz;
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
};
void rbfreport_init(void *_p, ae_state *_state, bool make_automatic);
void rbfreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void rbfreport_free(void *_p, bool make_automatic);

void rbfcreate(ae_int_t nx, ae_int_t ny, rbfmodel *s, ae_state *_state);
void rbfcreatecalcbuffer(rbfmodel *s, rbfcalcbuffer *buf, ae_state *_state);
void rbfsetpoints(rbfmodel *s, RMatrix *xy, ae_int_t n, ae_state *_state);
void rbfsetpointsandscales(rbfmodel *r, RMatrix *xy, ae_int_t n, RVector *s, ae_state *_state);
void rbfsetalgoqnn(rbfmodel *s, double q, double z, ae_state *_state);
void rbfsetalgomultilayer(rbfmodel *s, double rbase, ae_int_t nlayers, double lambdav, ae_state *_state);
void rbfsetalgohierarchical(rbfmodel *s, double rbase, ae_int_t nlayers, double lambdans, ae_state *_state);
void rbfsetlinterm(rbfmodel *s, ae_state *_state);
void rbfsetconstterm(rbfmodel *s, ae_state *_state);
void rbfsetzeroterm(rbfmodel *s, ae_state *_state);
void rbfsetv2bf(rbfmodel *s, ae_int_t bf, ae_state *_state);
void rbfsetv2its(rbfmodel *s, ae_int_t maxits, ae_state *_state);
void rbfsetv2supportr(rbfmodel *s, double r, ae_state *_state);
void rbfsetcond(rbfmodel *s, double epsort, double epserr, ae_int_t maxits, ae_state *_state);
void rbfbuildmodel(rbfmodel *s, rbfreport *rep, ae_state *_state);
double rbfcalc1(rbfmodel *s, double x0, ae_state *_state);
double rbfcalc2(rbfmodel *s, double x0, double x1, ae_state *_state);
double rbfcalc3(rbfmodel *s, double x0, double x1, double x2, ae_state *_state);
void rbftscalcbuf(rbfmodel *s, rbfcalcbuffer *buf, RVector *x, RVector *y, ae_state *_state);
void rbfcalcbuf(rbfmodel *s, RVector *x, RVector *y, ae_state *_state);
void rbfcalc(rbfmodel *s, RVector *x, RVector *y, ae_state *_state);
void rbfgridcalc2(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RMatrix *y, ae_state *_state);
void rbfgridcalc2vx(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, BVector *flagy, bool sparsey, RVector *y, ae_state *_state);
void rbfgridcalc2v(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *y, ae_state *_state);
void rbfgridcalc2vsubset(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, BVector *flagy, RVector *y, ae_state *_state);
void rbfgridcalc3vx(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *x2, ae_int_t n2, BVector *flagy, bool sparsey, RVector *y, ae_state *_state);
void rbfgridcalc3v(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *x2, ae_int_t n2, RVector *y, ae_state *_state);
void rbfgridcalc3vsubset(rbfmodel *s, RVector *x0, ae_int_t n0, RVector *x1, ae_int_t n1, RVector *x2, ae_int_t n2, BVector *flagy, RVector *y, ae_state *_state);
void rbfunpack(rbfmodel *s, ae_int_t *nx, ae_int_t *ny, RMatrix *xwr, ae_int_t *nc, RMatrix *v, ae_int_t *modelversion, ae_state *_state);
ae_int_t rbfgetmodelversion(rbfmodel *s, ae_state *_state);
double rbfpeekprogress(rbfmodel *s, ae_state *_state);
void rbfrequesttermination(rbfmodel *s, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
DecClass(rbfcalcbuffer, );
DecClass(rbfmodel, );
DecClass(rbfreport, double &rmserror; double &maxerror; ae_int_t &arows; ae_int_t &acols; ae_int_t &annz; ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype;);
void rbfserialize(rbfmodel &obj, std::string &s_out);
void rbfserialize(rbfmodel &obj, std::ostream &s_out);
void rbfunserialize(const std::string &s_in, rbfmodel &obj);
void rbfunserialize(const std::istream &s_in, rbfmodel &obj);

void rbfcreate(const ae_int_t nx, const ae_int_t ny, rbfmodel &s, const xparams _xparams = NonTH);
void rbfcreatecalcbuffer(const rbfmodel &s, rbfcalcbuffer &buf, const xparams _xparams = NonTH);
void rbfsetpoints(const rbfmodel &s, const real_2d_array &xy, const ae_int_t n, const xparams _xparams = NonTH);
void rbfsetpoints(const rbfmodel &s, const real_2d_array &xy, const xparams _xparams = NonTH);
void rbfsetpointsandscales(const rbfmodel &r, const real_2d_array &xy, const ae_int_t n, const real_1d_array &s, const xparams _xparams = NonTH);
void rbfsetpointsandscales(const rbfmodel &r, const real_2d_array &xy, const real_1d_array &s, const xparams _xparams = NonTH);
void rbfsetalgoqnn(const rbfmodel &s, const double q, const double z, const xparams _xparams = NonTH);
void rbfsetalgoqnn(const rbfmodel &s, const xparams _xparams = NonTH);
void rbfsetalgomultilayer(const rbfmodel &s, const double rbase, const ae_int_t nlayers, const double lambdav, const xparams _xparams = NonTH);
void rbfsetalgomultilayer(const rbfmodel &s, const double rbase, const ae_int_t nlayers, const xparams _xparams = NonTH);
void rbfsetalgohierarchical(const rbfmodel &s, const double rbase, const ae_int_t nlayers, const double lambdans, const xparams _xparams = NonTH);
void rbfsetlinterm(const rbfmodel &s, const xparams _xparams = NonTH);
void rbfsetconstterm(const rbfmodel &s, const xparams _xparams = NonTH);
void rbfsetzeroterm(const rbfmodel &s, const xparams _xparams = NonTH);
void rbfsetv2bf(const rbfmodel &s, const ae_int_t bf, const xparams _xparams = NonTH);
void rbfsetv2its(const rbfmodel &s, const ae_int_t maxits, const xparams _xparams = NonTH);
void rbfsetv2supportr(const rbfmodel &s, const double r, const xparams _xparams = NonTH);
void rbfbuildmodel(const rbfmodel &s, rbfreport &rep, const xparams _xparams = NonTH);
double rbfcalc1(const rbfmodel &s, const double x0, const xparams _xparams = NonTH);
double rbfcalc2(const rbfmodel &s, const double x0, const double x1, const xparams _xparams = NonTH);
double rbfcalc3(const rbfmodel &s, const double x0, const double x1, const double x2, const xparams _xparams = NonTH);
void rbftscalcbuf(const rbfmodel &s, const rbfcalcbuffer &buf, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void rbfcalcbuf(const rbfmodel &s, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void rbfcalc(const rbfmodel &s, const real_1d_array &x, real_1d_array &y, const xparams _xparams = NonTH);
void rbfgridcalc2(const rbfmodel &s, const real_1d_array &x0, const ae_int_t n0, const real_1d_array &x1, const ae_int_t n1, real_2d_array &y, const xparams _xparams = NonTH);
void rbfgridcalc2v(const rbfmodel &s, const real_1d_array &x0, const ae_int_t n0, const real_1d_array &x1, const ae_int_t n1, real_1d_array &y, const xparams _xparams = NonTH);
void rbfgridcalc2vsubset(const rbfmodel &s, const real_1d_array &x0, const ae_int_t n0, const real_1d_array &x1, const ae_int_t n1, const boolean_1d_array &flagy, real_1d_array &y, const xparams _xparams = NonTH);
void rbfgridcalc3v(const rbfmodel &s, const real_1d_array &x0, const ae_int_t n0, const real_1d_array &x1, const ae_int_t n1, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y, const xparams _xparams = NonTH);
void rbfgridcalc3vsubset(const rbfmodel &s, const real_1d_array &x0, const ae_int_t n0, const real_1d_array &x1, const ae_int_t n1, const real_1d_array &x2, const ae_int_t n2, const boolean_1d_array &flagy, real_1d_array &y, const xparams _xparams = NonTH);
void rbfunpack(const rbfmodel &s, ae_int_t &nx, ae_int_t &ny, real_2d_array &xwr, ae_int_t &nc, real_2d_array &v, ae_int_t &modelversion, const xparams _xparams = NonTH);
ae_int_t rbfgetmodelversion(const rbfmodel &s, const xparams _xparams = NonTH);
double rbfpeekprogress(const rbfmodel &s, const xparams _xparams = NonTH);
void rbfrequesttermination(const rbfmodel &s, const xparams _xparams = NonTH);
} // end of namespace alglib

#endif // OnceOnly
