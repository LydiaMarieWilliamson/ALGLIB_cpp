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
// Depends on: (SpecialFunctions) GAMMAFUNC
// Depends on: (LinAlg) EVD
namespace alglib_impl {
void gqgeneraterec(RVector *alpha, RVector *beta, double mu0, ae_int_t n, ae_int_t *info, RVector *x, RVector *w);
void gqgenerategausslobattorec(RVector *alpha, RVector *beta, double mu0, double a, double b, ae_int_t n, ae_int_t *info, RVector *x, RVector *w);
void gqgenerategaussradaurec(RVector *alpha, RVector *beta, double mu0, double a, ae_int_t n, ae_int_t *info, RVector *x, RVector *w);
void gqgenerategausslegendre(ae_int_t n, ae_int_t *info, RVector *x, RVector *w);
void gqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t *info, RVector *x, RVector *w);
void gqgenerategausslaguerre(ae_int_t n, double alpha, ae_int_t *info, RVector *x, RVector *w);
void gqgenerategausshermite(ae_int_t n, ae_int_t *info, RVector *x, RVector *w);
} // end of namespace alglib_impl

namespace alglib {
void gqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
void gqgenerategausslobattorec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const double b, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
void gqgenerategaussradaurec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
void gqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
void gqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &w);
void gqgenerategausslaguerre(const ae_int_t n, const double alpha, ae_int_t &info, real_1d_array &x, real_1d_array &w);
void gqgenerategausshermite(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
} // end of namespace alglib

// === GKQ Package ===
// Depends on: GQ
namespace alglib_impl {
void gkqgeneraterec(RVector *alpha, RVector *beta, double mu0, ae_int_t n, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss);
void gkqlegendrecalc(ae_int_t n, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss);
void gkqlegendretbl(ae_int_t n, RVector *x, RVector *wkronrod, RVector *wgauss, double *eps);
void gkqgenerategausslegendre(ae_int_t n, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss);
void gkqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t *info, RVector *x, RVector *wkronrod, RVector *wgauss);
} // end of namespace alglib_impl

namespace alglib {
void gkqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
void gkqlegendrecalc(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
void gkqlegendretbl(const ae_int_t n, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, double &eps);
void gkqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
void gkqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
} // end of namespace alglib

// === AUTOGK Package ===
// Depends on: GKQ
namespace alglib_impl {
struct autogkreport {
   ae_int_t terminationtype;
   ae_int_t nfev;
   ae_int_t nintervals;
};
void autogkreport_init(void *_p, bool make_automatic);
void autogkreport_copy(void *_dst, const void *_src, bool make_automatic);
void autogkreport_free(void *_p, bool make_automatic);

struct autogkinternalstate {
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
   ae_int_t PQ;
};
void autogkinternalstate_init(void *_p, bool make_automatic);
void autogkinternalstate_copy(void *_dst, const void *_src, bool make_automatic);
void autogkinternalstate_free(void *_p, bool make_automatic);

struct autogkstate {
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
   ae_int_t PQ;
   double v;
   ae_int_t terminationtype;
   ae_int_t nfev;
   ae_int_t nintervals;
};
void autogkstate_init(void *_p, bool make_automatic);
void autogkstate_copy(void *_dst, const void *_src, bool make_automatic);
void autogkstate_free(void *_p, bool make_automatic);

void autogksmoothw(double a, double b, double xwidth, autogkstate *state);
void autogksmooth(double a, double b, autogkstate *state);
void autogksingular(double a, double b, double alpha, double beta, autogkstate *state);
bool autogkiteration(autogkstate *state);
void autogkresults(autogkstate *state, double *v, autogkreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(autogkreport, ae_int_t &terminationtype; ae_int_t &nfev; ae_int_t &nintervals;);
DecClass(autogkstate, bool &needf; double &x; double &xminusa; double &bminusx; double &f;);

void autogksmoothw(const double a, const double b, const double xwidth, autogkstate &state);
void autogksmooth(const double a, const double b, autogkstate &state);
void autogksingular(const double a, const double b, const double alpha, const double beta, autogkstate &state);
bool autogkiteration(const autogkstate &state);
void autogkintegrate(autogkstate &state, void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr), void *ptr = NULL);
void autogkresults(const autogkstate &state, double &v, autogkreport &rep);
} // end of namespace alglib

#endif // OnceOnly
