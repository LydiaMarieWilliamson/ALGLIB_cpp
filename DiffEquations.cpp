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
#define InAlgLib
#include "DiffEquations.h"

// === ODESOLVER Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
static const double odesolver_odesolvermaxgrow = 3.0;
static const double odesolver_odesolvermaxshrink = 10.0;
static void odesolver_odesolverinit(ae_int_t solvertype, RVector y, ae_int_t n, RVector x, ae_int_t m, double eps, double h, odesolverstate *state);

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
//                 * M>2 means that you need both Y(X[M-1]) and M-2 values at
//                   intermediate points
//                 * M=2 means that you want just to integrate from  X[0]  to
//                   X[1] and don't interested in intermediate values.
//                 * M=1 means that you don't want to integrate :)
//                   it is degenerate case, but it will be handled correctly.
//                 * M<1 means error
//     Eps     -   tolerance (absolute/relative error on each  step  will  be
//                 less than Eps). When passing:
//                 * Eps>0, it means desired ABSOLUTE error
//                 * Eps<0, it means desired RELATIVE error.  Relative errors
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
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolverrkck(RVector y, ae_int_t n, RVector x, ae_int_t m, double eps, double h, odesolverstate *state) {

   SetObj(odesolverstate, state);

   ae_assert(n >= 1, "ODESolverRKCK: N<1!");
   ae_assert(m >= 1, "ODESolverRKCK: M<1!");
   ae_assert(y->cnt >= n, "ODESolverRKCK: Length(Y)<N!");
   ae_assert(x->cnt >= m, "ODESolverRKCK: Length(X)<M!");
   ae_assert(isfinitevector(y, n), "ODESolverRKCK: Y contains infinite or NaN values!");
   ae_assert(isfinitevector(x, m), "ODESolverRKCK: Y contains infinite or NaN values!");
   ae_assert(isfinite(eps), "ODESolverRKCK: Eps is not finite!");
   ae_assert(eps != 0.0, "ODESolverRKCK: Eps is zero!");
   ae_assert(isfinite(h), "ODESolverRKCK: H is not finite!");
   odesolver_odesolverinit(0, y, n, x, m, eps, h, state);
}

//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
bool odesolveriteration(odesolverstate *state) {
   AutoS ae_int_t n;
   AutoS ae_int_t m;
   AutoS ae_int_t i;
   AutoS ae_int_t j;
   AutoS ae_int_t k;
   AutoS double xc;
   AutoS double v;
   AutoS double h;
   AutoS double h2;
   AutoS bool gridpoint;
   AutoS double err;
   AutoS double maxgrowpow;
   AutoS ae_int_t klimit;
// Reverse communication preparations
// I know it looks ugly, but it works the same way anywhere from C++ to Python.
// This code initializes locals by:
// * random values determined during code generation - on first subroutine call
// * values from previous call - on subsequent calls
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0;
      default: goto Exit;
   }
Spawn:
   n = 359;
   m = -58;
   i = -919;
   j = -909;
   k = 81;
   klimit = 255;
   gridpoint = false;
   xc = -788;
   v = 809;
   h = 205;
   h2 = -838;
   err = 939;
   maxgrowpow = -526;

// prepare
   if (state->repterminationtype != 0) {
      goto Exit;
   }
   n = state->n;
   m = state->m;
   h = state->h;
   maxgrowpow = pow(odesolver_odesolvermaxgrow, 5.0);
   state->repnfev = 0;

// some preliminary checks for internal errors
// after this we assume that H>0 and M>1
   ae_assert(state->h > 0.0, "ODESolver: internal error");
   ae_assert(m > 1, "ODESolverIteration: internal error");

// choose solver
   if (state->solvertype == 0) {
   // Cask-Karp solver
   // Prepare coefficients table.
   // Check it for errors
      ae_vector_set_length(&state->rka, 6);
      state->rka.ptr.p_double[0] = 0.0;
      state->rka.ptr.p_double[1] = 1.0 / 5.0;
      state->rka.ptr.p_double[2] = 3.0 / 10.0;
      state->rka.ptr.p_double[3] = 3.0 / 5.0;
      state->rka.ptr.p_double[4] = 1.0;
      state->rka.ptr.p_double[5] = 7.0 / 8.0;
      ae_matrix_set_length(&state->rkb, 6, 5);
      state->rkb.ptr.pp_double[1][0] = 1.0 / 5.0;
      state->rkb.ptr.pp_double[2][0] = 3.0 / 40.0;
      state->rkb.ptr.pp_double[2][1] = 9.0 / 40.0;
      state->rkb.ptr.pp_double[3][0] = 3.0 / 10.0;
      state->rkb.ptr.pp_double[3][1] = -9.0 / 10.0;
      state->rkb.ptr.pp_double[3][2] = 6.0 / 5.0;
      state->rkb.ptr.pp_double[4][0] = -11.0 / 54.0;
      state->rkb.ptr.pp_double[4][1] = 5.0 / 2.0;
      state->rkb.ptr.pp_double[4][2] = -70.0 / 27.0;
      state->rkb.ptr.pp_double[4][3] = 35.0 / 27.0;
      state->rkb.ptr.pp_double[5][0] = 1631.0 / 55296.0;
      state->rkb.ptr.pp_double[5][1] = 175.0 / 512.0;
      state->rkb.ptr.pp_double[5][2] = 575.0 / 13824.0;
      state->rkb.ptr.pp_double[5][3] = 44275.0 / 110592.0;
      state->rkb.ptr.pp_double[5][4] = 253.0 / 4096.0;
      ae_vector_set_length(&state->rkc, 6);
      state->rkc.ptr.p_double[0] = 37.0 / 378.0;
      state->rkc.ptr.p_double[1] = 0.0;
      state->rkc.ptr.p_double[2] = 250.0 / 621.0;
      state->rkc.ptr.p_double[3] = 125.0 / 594.0;
      state->rkc.ptr.p_double[4] = 0.0;
      state->rkc.ptr.p_double[5] = 512.0 / 1771.0;
      ae_vector_set_length(&state->rkcs, 6);
      state->rkcs.ptr.p_double[0] = 2825.0 / 27648.0;
      state->rkcs.ptr.p_double[1] = 0.0;
      state->rkcs.ptr.p_double[2] = 18575.0 / 48384.0;
      state->rkcs.ptr.p_double[3] = 13525.0 / 55296.0;
      state->rkcs.ptr.p_double[4] = 277.0 / 14336.0;
      state->rkcs.ptr.p_double[5] = 1.0 / 4.0;
      ae_matrix_set_length(&state->rkk, 6, n);

   // Main cycle consists of two iterations:
   // * outer where we travel from X[i-1] to X[i]
   // * inner where we travel inside [X[i-1],X[i]]
      ae_matrix_set_length(&state->ytbl, m, n);
      ae_vector_set_length(&state->escale, n);
      ae_vector_set_length(&state->yn, n);
      ae_vector_set_length(&state->yns, n);
      xc = state->xg.ptr.p_double[0];
      ae_v_move(state->ytbl.ptr.pp_double[0], 1, state->yc.ptr.p_double, 1, n);
      for (j = 0; j < n; j++) {
         state->escale.ptr.p_double[j] = 0.0;
      }
      for (i = 1; i < m; i++) {
      // begin inner iteration
         while (true) {
         // truncate step if needed (beyond right boundary).
         // determine should we store X or not
            if (xc + h >= state->xg.ptr.p_double[i]) {
               h = state->xg.ptr.p_double[i] - xc;
               gridpoint = true;
            } else {
               gridpoint = false;
            }

         // Update error scale maximums
         //
         // These maximums are initialized by zeros,
         // then updated every iterations.
            for (j = 0; j < n; j++) {
               state->escale.ptr.p_double[j] = ae_maxreal(state->escale.ptr.p_double[j], fabs(state->yc.ptr.p_double[j]));
            }

         // make one step:
         // 1. calculate all info needed to do step
         // 2. update errors scale maximums using values/derivatives
         //    obtained during (1)
         //
         // Take into account that we use scaling of X to reduce task
         // to the form where x[0] < x[1] < ... < x[n-1]. So X is
         // replaced by x=xscale*t, and dy/dx=f(y,x) is replaced
         // by dy/dt=xscale*f(y,xscale*t).
            ae_v_move(state->yn.ptr.p_double, 1, state->yc.ptr.p_double, 1, n);
            ae_v_move(state->yns.ptr.p_double, 1, state->yc.ptr.p_double, 1, n);
            for (k = 0; k < 6; k++) {
            // prepare data for the next update of YN/YNS
               state->x = state->xscale * (xc + state->rka.ptr.p_double[k] * h);
               ae_v_move(state->y.ptr.p_double, 1, state->yc.ptr.p_double, 1, n);
               for (j = 0; j < k; j++) {
                  v = state->rkb.ptr.pp_double[k][j];
                  ae_v_addd(state->y.ptr.p_double, 1, state->rkk.ptr.pp_double[j], 1, n, v);
               }
               state->needdy = true;
               state->PQ = 0; goto Pause; Resume0:
               state->needdy = false;
               state->repnfev = state->repnfev + 1;
               v = h * state->xscale;
               ae_v_moved(state->rkk.ptr.pp_double[k], 1, state->dy.ptr.p_double, 1, n, v);

            // update YN/YNS
               v = state->rkc.ptr.p_double[k];
               ae_v_addd(state->yn.ptr.p_double, 1, state->rkk.ptr.pp_double[k], 1, n, v);
               v = state->rkcs.ptr.p_double[k];
               ae_v_addd(state->yns.ptr.p_double, 1, state->rkk.ptr.pp_double[k], 1, n, v);
            }

         // estimate error
            err = 0.0;
            for (j = 0; j < n; j++) {
               if (!state->fraceps) {

               // absolute error is estimated
                  err = ae_maxreal(err, fabs(state->yn.ptr.p_double[j] - state->yns.ptr.p_double[j]));
               } else {

               // Relative error is estimated
                  v = state->escale.ptr.p_double[j];
                  if (v == 0.0) {
                     v = 1.0;
                  }
                  err = ae_maxreal(err, fabs(state->yn.ptr.p_double[j] - state->yns.ptr.p_double[j]) / v);
               }
            }

         // calculate new step, restart if necessary
            if (maxgrowpow * err <= state->eps) {
               h2 = odesolver_odesolvermaxgrow * h;
            } else {
               h2 = h * pow(state->eps / err, 0.2);
            }
            if (h2 < h / odesolver_odesolvermaxshrink) {
               h2 = h / odesolver_odesolvermaxshrink;
            }
            if (err > state->eps) {
               h = h2;
               continue;
            }
         // advance position
            xc = xc + h;
            ae_v_move(state->yc.ptr.p_double, 1, state->yn.ptr.p_double, 1, n);

         // update H
            h = h2;

         // break on grid point
            if (gridpoint) {
               break;
            }
         }
      // save result
         ae_v_move(state->ytbl.ptr.pp_double[i], 1, state->yc.ptr.p_double, 1, n);
      }
      state->repterminationtype = 1;
      goto Exit;
   }
   goto Exit;

Pause:
   return true;
Exit:
   state->PQ = -1;
   return false;
}

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
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolverresults(odesolverstate *state, ae_int_t *m, RVector xtbl, RMatrix ytbl, odesolverreport *rep) {
   double v;
   ae_int_t i;

   *m = 0;
   SetVector(xtbl);
   SetMatrix(ytbl);
   SetObj(odesolverreport, rep);

   rep->terminationtype = state->repterminationtype;
   if (rep->terminationtype > 0) {
      *m = state->m;
      rep->nfev = state->repnfev;
      ae_vector_set_length(xtbl, state->m);
      v = state->xscale;
      ae_v_moved(xtbl->ptr.p_double, 1, state->xg.ptr.p_double, 1, state->m, v);
      ae_matrix_set_length(ytbl, state->m, state->n);
      for (i = 0; i < state->m; i++) {
         ae_v_move(ytbl->ptr.pp_double[i], 1, state->ytbl.ptr.pp_double[i], 1, state->n);
      }
   } else {
      rep->nfev = 0;
   }
}

// Internal initialization subroutine
static void odesolver_odesolverinit(ae_int_t solvertype, RVector y, ae_int_t n, RVector x, ae_int_t m, double eps, double h, odesolverstate *state) {
   ae_int_t i;
   double v;

   SetObj(odesolverstate, state);

// Prepare RComm
   state->PQ = -1;
   state->needdy = false;

// check parameters.
   if ((n <= 0 || m < 1) || eps == 0.0) {
      state->repterminationtype = -1;
      return;
   }
   if (h < 0.0) {
      h = -h;
   }
// quick exit if necessary.
// after this block we assume that M>1
   if (m == 1) {
      state->repnfev = 0;
      state->repterminationtype = 1;
      ae_matrix_set_length(&state->ytbl, 1, n);
      ae_v_move(state->ytbl.ptr.pp_double[0], 1, y->ptr.p_double, 1, n);
      ae_vector_set_length(&state->xg, m);
      ae_v_move(state->xg.ptr.p_double, 1, x->ptr.p_double, 1, m);
      return;
   }
// check again: correct order of X[]
   if (x->ptr.p_double[1] == x->ptr.p_double[0]) {
      state->repterminationtype = -2;
      return;
   }
   for (i = 1; i < m; i++) {
      if ((x->ptr.p_double[1] > x->ptr.p_double[0] && x->ptr.p_double[i] <= x->ptr.p_double[i - 1]) || (x->ptr.p_double[1] < x->ptr.p_double[0] && x->ptr.p_double[i] >= x->ptr.p_double[i - 1])) {
         state->repterminationtype = -2;
         return;
      }
   }

// auto-select H if necessary
   if (h == 0.0) {
      v = fabs(x->ptr.p_double[1] - x->ptr.p_double[0]);
      for (i = 2; i < m; i++) {
         v = ae_minreal(v, fabs(x->ptr.p_double[i] - x->ptr.p_double[i - 1]));
      }
      h = 0.001 * v;
   }
// store parameters
   state->n = n;
   state->m = m;
   state->h = h;
   state->eps = fabs(eps);
   state->fraceps = eps < 0.0;
   ae_vector_set_length(&state->xg, m);
   ae_v_move(state->xg.ptr.p_double, 1, x->ptr.p_double, 1, m);
   if (x->ptr.p_double[1] > x->ptr.p_double[0]) {
      state->xscale = 1.0;
   } else {
      state->xscale = -1.0;
      ae_v_muld(state->xg.ptr.p_double, 1, m, -1);
   }
   ae_vector_set_length(&state->yc, n);
   ae_v_move(state->yc.ptr.p_double, 1, y->ptr.p_double, 1, n);
   state->solvertype = solvertype;
   state->repterminationtype = 0;

// Allocate arrays
   ae_vector_set_length(&state->y, n);
   ae_vector_set_length(&state->dy, n);
}

void odesolverstate_init(void *_p, bool make_automatic) {
   odesolverstate *p = (odesolverstate *) _p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->yc, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->escale, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xg, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->y, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->dy, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->ytbl, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->yn, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->yns, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rka, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rkc, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rkcs, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->rkb, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->rkk, 0, 0, DT_REAL, make_automatic);
}

void odesolverstate_copy(void *_dst, void *_src, bool make_automatic) {
   odesolverstate *dst = (odesolverstate *) _dst;
   odesolverstate *src = (odesolverstate *) _src;
   dst->n = src->n;
   dst->m = src->m;
   dst->xscale = src->xscale;
   dst->h = src->h;
   dst->eps = src->eps;
   dst->fraceps = src->fraceps;
   ae_vector_copy(&dst->yc, &src->yc, make_automatic);
   ae_vector_copy(&dst->escale, &src->escale, make_automatic);
   ae_vector_copy(&dst->xg, &src->xg, make_automatic);
   dst->solvertype = src->solvertype;
   dst->needdy = src->needdy;
   dst->x = src->x;
   ae_vector_copy(&dst->y, &src->y, make_automatic);
   ae_vector_copy(&dst->dy, &src->dy, make_automatic);
   ae_matrix_copy(&dst->ytbl, &src->ytbl, make_automatic);
   dst->repterminationtype = src->repterminationtype;
   dst->repnfev = src->repnfev;
   ae_vector_copy(&dst->yn, &src->yn, make_automatic);
   ae_vector_copy(&dst->yns, &src->yns, make_automatic);
   ae_vector_copy(&dst->rka, &src->rka, make_automatic);
   ae_vector_copy(&dst->rkc, &src->rkc, make_automatic);
   ae_vector_copy(&dst->rkcs, &src->rkcs, make_automatic);
   ae_matrix_copy(&dst->rkb, &src->rkb, make_automatic);
   ae_matrix_copy(&dst->rkk, &src->rkk, make_automatic);
   dst->PQ = src->PQ;
}

void odesolverstate_free(void *_p, bool make_automatic) {
   odesolverstate *p = (odesolverstate *) _p;
   ae_touch_ptr((void *)p);
   ae_vector_free(&p->yc, make_automatic);
   ae_vector_free(&p->escale, make_automatic);
   ae_vector_free(&p->xg, make_automatic);
   ae_vector_free(&p->y, make_automatic);
   ae_vector_free(&p->dy, make_automatic);
   ae_matrix_free(&p->ytbl, make_automatic);
   ae_vector_free(&p->yn, make_automatic);
   ae_vector_free(&p->yns, make_automatic);
   ae_vector_free(&p->rka, make_automatic);
   ae_vector_free(&p->rkc, make_automatic);
   ae_vector_free(&p->rkcs, make_automatic);
   ae_matrix_free(&p->rkb, make_automatic);
   ae_matrix_free(&p->rkk, make_automatic);
}

void odesolverreport_init(void *_p, bool make_automatic) {
   odesolverreport *p = (odesolverreport *) _p;
   ae_touch_ptr((void *)p);
}

void odesolverreport_copy(void *_dst, void *_src, bool make_automatic) {
   odesolverreport *dst = (odesolverreport *) _dst;
   odesolverreport *src = (odesolverreport *) _src;
   dst->nfev = src->nfev;
   dst->terminationtype = src->terminationtype;
}

void odesolverreport_free(void *_p, bool make_automatic) {
   odesolverreport *p = (odesolverreport *) _p;
   ae_touch_ptr((void *)p);
}
} // end of namespace alglib_impl

namespace alglib {
DefClass(odesolverstate, AndD DecVal(needdy) AndD DecVar(y) AndD DecVar(dy) AndD DecVal(x))
DefClass(odesolverreport, AndD DecVal(nfev) AndD DecVal(terminationtype))

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
//                 * M>2 means that you need both Y(X[M-1]) and M-2 values at
//                   intermediate points
//                 * M=2 means that you want just to integrate from  X[0]  to
//                   X[1] and don't interested in intermediate values.
//                 * M=1 means that you don't want to integrate :)
//                   it is degenerate case, but it will be handled correctly.
//                 * M<1 means error
//     Eps     -   tolerance (absolute/relative error on each  step  will  be
//                 less than Eps). When passing:
//                 * Eps>0, it means desired ABSOLUTE error
//                 * Eps<0, it means desired RELATIVE error.  Relative errors
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
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::odesolverrkck(ConstT(ae_vector, y), n, ConstT(ae_vector, x), m, eps, h, ConstT(odesolverstate, state));
   alglib_impl::ae_state_clear();
}

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
//                 * M>2 means that you need both Y(X[M-1]) and M-2 values at
//                   intermediate points
//                 * M=2 means that you want just to integrate from  X[0]  to
//                   X[1] and don't interested in intermediate values.
//                 * M=1 means that you don't want to integrate :)
//                   it is degenerate case, but it will be handled correctly.
//                 * M<1 means error
//     Eps     -   tolerance (absolute/relative error on each  step  will  be
//                 less than Eps). When passing:
//                 * Eps>0, it means desired ABSOLUTE error
//                 * Eps<0, it means desired RELATIVE error.  Relative errors
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
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
#if !defined AE_NO_EXCEPTIONS
void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state) {
   ae_int_t n = y.length();
   ae_int_t m = x.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::odesolverrkck(ConstT(ae_vector, y), n, ConstT(ae_vector, x), m, eps, h, ConstT(odesolverstate, state));
   alglib_impl::ae_state_clear();
}
#endif

// This function provides reverse communication interface
// Reverse communication interface is not documented or recommended to use.
// See below for functions which provide better documented API
bool odesolveriteration(const odesolverstate &state) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::odesolveriteration(ConstT(odesolverstate, state));
   alglib_impl::ae_state_clear();
   return Ok;
}

void odesolversolve(odesolverstate &state, void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr), void *ptr) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ae_assert(diff != NULL, "ALGLIB: error in 'odesolversolve()' (diff is NULL)");
   while (alglib_impl::odesolveriteration(state.c_ptr()))
   BegPoll
      if (state.needdy) diff(state.y, state.x, state.dy, ptr);
      else alglib_impl::ae_assert(false, "ALGLIB: unexpected error in 'odesolversolve'");
   EndPoll
   alglib_impl::ae_state_clear();
}

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
//
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::odesolverresults(ConstT(odesolverstate, state), &m, ConstT(ae_vector, xtbl), ConstT(ae_matrix, ytbl), ConstT(odesolverreport, rep));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib