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
#ifndef OnceOnlyIntegration_h
#define OnceOnlyIntegration_h

#include "LinAlg.h"
#include "SpecialFunctions.h"

// === GQ Package ===
namespace alglib_impl {
void gqgeneraterec(RVector *alpha, RVector *beta, double mu0, ae_int_t n, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
void gqgenerategausslobattorec(RVector *alpha, RVector *beta, double mu0, double a, double b, ae_int_t n, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
void gqgenerategaussradaurec(RVector *alpha, RVector *beta, double mu0, double a, ae_int_t n, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
void gqgenerategausslegendre(ae_int_t n, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
void gqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
void gqgenerategausslaguerre(ae_int_t n, double alpha, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
void gqgenerategausshermite(ae_int_t n, ae_int_t *info, RVector *x, RVector *w, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void gqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
void gqgenerategausslobattorec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const double b, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
void gqgenerategaussradaurec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
void gqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
void gqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
void gqgenerategausslaguerre(const ae_int_t n, const double alpha, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
void gqgenerategausshermite(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w, const xparams _xparams = xdefault);
} // end of namespace alglib

// === GKQ Package ===
namespace alglib_impl {
void gkqgeneraterec(RVector *alpha, RVector *beta, double mu0, ae_int_t n, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss, ae_state *_state);
void gkqgenerategausslegendre(ae_int_t n, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss, ae_state *_state);
void gkqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss, ae_state *_state);
void gkqlegendrecalc(ae_int_t n, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss, ae_state *_state);
void gkqlegendretbl(ae_int_t n, RVector *x, RVector *wkronrod, RVector *wgauss, double *eps, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void gkqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, const xparams _xparams = xdefault);
void gkqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, const xparams _xparams = xdefault);
void gkqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, const xparams _xparams = xdefault);
void gkqlegendrecalc(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, const xparams _xparams = xdefault);
void gkqlegendretbl(const ae_int_t n, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, double &eps, const xparams _xparams = xdefault);
} // end of namespace alglib

// === AUTOGK Package ===
namespace alglib_impl {
typedef struct {
   ae_int_t terminationtype;
   ae_int_t nfev;
   ae_int_t nintervals;
} autogkreport;
typedef struct {
   double a;
   double b;
   double eps;
   double xwidth;
   double x;
   double f;
   ae_int_t info;
   double r;
   ae_matrix heap;
   ae_int_t heapsize;
   ae_int_t heapwidth;
   ae_int_t heapused;
   double sumerr;
   double sumabs;
   ae_vector qn;
   ae_vector wg;
   ae_vector wk;
   ae_vector wr;
   ae_int_t n;
   rcommstate rstate;
} autogkinternalstate;
typedef struct {
   double a;
   double b;
   double alpha;
   double beta;
   double xwidth;
   double x;
   double xminusa;
   double bminusx;
   bool needf;
   double f;
   ae_int_t wrappermode;
   autogkinternalstate internalstate;
   rcommstate rstate;
   double v;
   ae_int_t terminationtype;
   ae_int_t nfev;
   ae_int_t nintervals;
} autogkstate;

void autogksmooth(double a, double b, autogkstate *state, ae_state *_state);
void autogksmoothw(double a, double b, double xwidth, autogkstate *state, ae_state *_state);
void autogksingular(double a, double b, double alpha, double beta, autogkstate *state, ae_state *_state);
bool autogkiteration(autogkstate *state, ae_state *_state);
void autogkresults(autogkstate *state, double *v, autogkreport *rep, ae_state *_state);
void _autogkreport_init(void *_p, ae_state *_state, bool make_automatic);
void _autogkreport_init_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void _autogkreport_clear(void *_p);
void _autogkreport_destroy(void *_p);
void _autogkinternalstate_init(void *_p, ae_state *_state, bool make_automatic);
void _autogkinternalstate_init_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void _autogkinternalstate_clear(void *_p);
void _autogkinternalstate_destroy(void *_p);
void _autogkstate_init(void *_p, ae_state *_state, bool make_automatic);
void _autogkstate_init_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void _autogkstate_clear(void *_p);
void _autogkstate_destroy(void *_p);
} // end of namespace alglib_impl

namespace alglib {
DecClass(autogkreport, ae_int_t &terminationtype; ae_int_t &nfev; ae_int_t &nintervals;);
DecClass(autogkstate, bool &needf; double &x; double &xminusa; double &bminusx; double &f;);

void autogksmooth(const double a, const double b, autogkstate &state, const xparams _xparams = xdefault);
void autogksmoothw(const double a, const double b, const double xwidth, autogkstate &state, const xparams _xparams = xdefault);
void autogksingular(const double a, const double b, const double alpha, const double beta, autogkstate &state, const xparams _xparams = xdefault);
bool autogkiteration(const autogkstate &state, const xparams _xparams = xdefault);
void autogkintegrate(autogkstate &state, void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr), void *ptr = NULL, const xparams _xparams = xdefault);
void autogkresults(const autogkstate &state, double &v, autogkreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

#endif // OnceOnly
