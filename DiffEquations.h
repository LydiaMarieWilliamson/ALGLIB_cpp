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
#ifndef OnceOnlyDiffEquations_h
#define OnceOnlyDiffEquations_h

#include "AlgLibInternal.h"

// === ODESOLVER Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
struct odesolverstate {
   ae_int_t n;
   ae_int_t m;
   double xscale;
   double h;
   double eps;
   bool fraceps;
   ae_vector yc;
   ae_vector escale;
   ae_vector xg;
   ae_int_t solvertype;
   bool needdy;
   double x;
   ae_vector y;
   ae_vector dy;
   ae_matrix ytbl;
   ae_int_t repterminationtype;
   ae_int_t repnfev;
   ae_vector yn;
   ae_vector yns;
   ae_vector rka;
   ae_vector rkc;
   ae_vector rkcs;
   ae_matrix rkb;
   ae_matrix rkk;
   ae_int_t PQ;
};
void odesolverstate_init(void *_p, bool make_automatic);
void odesolverstate_copy(void *_dst, void *_src, bool make_automatic);
void odesolverstate_free(void *_p, bool make_automatic);

struct odesolverreport {
   ae_int_t nfev;
   ae_int_t terminationtype;
};
void odesolverreport_init(void *_p, bool make_automatic);
void odesolverreport_copy(void *_dst, void *_src, bool make_automatic);
void odesolverreport_free(void *_p, bool make_automatic);

void odesolverrkck(RVector *y, ae_int_t n, RVector *x, ae_int_t m, double eps, double h, odesolverstate *state);
bool odesolveriteration(odesolverstate *state);
void odesolverresults(odesolverstate *state, ae_int_t *m, RVector *xtbl, RMatrix *ytbl, odesolverreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(odesolverstate, bool &needdy; real_1d_array y; real_1d_array dy; double &x;);
DecClass(odesolverreport, ae_int_t &nfev; ae_int_t &terminationtype;);

void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state);
void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state);
bool odesolveriteration(const odesolverstate &state);
void odesolversolve(odesolverstate &state, void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr), void *ptr = NULL);
void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep);
} // end of namespace alglib

#endif // OnceOnly
