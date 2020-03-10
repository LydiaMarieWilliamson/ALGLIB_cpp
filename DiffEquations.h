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
// bool needdy; //(@) Redundant.
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
DecClass(odesolverstate, real_1d_array y; real_1d_array dy; double &x;);
DecClass(odesolverreport, ae_int_t &nfev; ae_int_t &terminationtype;);

// Cash-Karp adaptive ODE solver.
//
// This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
// (here Y may be single variable or vector of N variables).
//
// Inputs:
//     Y       -   initial conditions, array[0..N-1].
//                 contains values of Y[] at X[0]
//     N       -   system size
//     X       -   points at which Y should be tabulated, array[0..M-1]
//                 integrations starts at X[0], ends at X[M-1],  intermediate
//                 values at X[i] are returned too.
//                 SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!
//     M       -   number of intermediate points + first point + last point:
//                 * M > 2 means that you need both Y(X[M-1]) and M-2 values at
//                   intermediate points
//                 * M=2 means that you want just to integrate from  X[0]  to
//                   X[1] and don't interested in intermediate values.
//                 * M=1 means that you don't want to integrate :)
//                   it is degenerate case, but it will be handled correctly.
//                 * M < 1 means error
//     Eps     -   tolerance (absolute/relative error on each  step  will  be
//                 less than Eps). When passing:
//                 * Eps > 0, it means desired ABSOLUTE error
//                 * Eps < 0, it means desired RELATIVE error.  Relative errors
//                   are calculated with respect to maximum values of  Y seen
//                   so far. Be careful to use this criterion  when  starting
//                   from Y[] that are close to zero.
//     H       -   initial  step  lenth,  it  will  be adjusted automatically
//                 after the first  step.  If  H=0,  step  will  be  selected
//                 automatically  (usualy  it  will  be  equal  to  0.001  of
//                 min(x[i]-x[j])).
//
// Outputs:
//     State   -   structure which stores algorithm state between  subsequent
//                 calls of OdeSolverIteration. Used for reverse communication.
//                 This structure should be passed  to the OdeSolverIteration
//                 subroutine.
//
// SEE ALSO
//     AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state);
void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state);

// This function provides reverse communication interface
// Reverse communication interface is not documented or recommended for use.
// See below for functions which provide better documented API
bool odesolveriteration(const odesolverstate &state);

// This function is used to launch iterations of ODE solver
//
// It accepts following parameters:
//     diff    -   callback which calculates dy/dx for given y and x
//     ptr     -   optional pointer which is passed to diff; can be NULL
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolversolve(odesolverstate &state, void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr), void *ptr = NULL);

// ODE solver results
//
// Called after OdeSolverIteration returned False.
//
// Inputs:
//     State   -   algorithm state (used by OdeSolverIteration).
//
// Outputs:
//     M       -   number of tabulated values, M >= 1
//     XTbl    -   array[0..M-1], values of X
//     YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]
//     Rep     -   solver report:
//                 * Rep.TerminationType completetion code:
//                     * -2    X is not ordered  by  ascending/descending  or
//                             there are non-distinct X[],  i.e.  X[i]=X[i+1]
//                     * -1    incorrect parameters were specified
//                     *  1    task has been solved
//                 * Rep.NFEV contains number of function calculations
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep);
} // end of namespace alglib

#endif // OnceOnly
