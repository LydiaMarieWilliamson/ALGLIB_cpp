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
static double odesolver_odesolvermaxgrow = 3.0;
static double odesolver_odesolvermaxshrink = 10.0;
static double odesolver_odesolverguaranteeddecay = 0.9;
static void odesolver_odesolverinit(ae_int_t solvertype, RVector *y, ae_int_t n, RVector *x, ae_int_t m, double eps, double h, odesolverstate *state, ae_state *_state);

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
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
// API: void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state, const xparams _xparams = xdefault);
// API: void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state, const xparams _xparams = xdefault);
void odesolverrkck(RVector *y, ae_int_t n, RVector *x, ae_int_t m, double eps, double h, odesolverstate *state, ae_state *_state) {

   odesolverstate_free(state, true);

   ae_assert(n >= 1, "ODESolverRKCK: N<1!", _state);
   ae_assert(m >= 1, "ODESolverRKCK: M<1!", _state);
   ae_assert(y->cnt >= n, "ODESolverRKCK: Length(Y)<N!", _state);
   ae_assert(x->cnt >= m, "ODESolverRKCK: Length(X)<M!", _state);
   ae_assert(isfinitevector(y, n, _state), "ODESolverRKCK: Y contains infinite or NaN values!", _state);
   ae_assert(isfinitevector(x, m, _state), "ODESolverRKCK: Y contains infinite or NaN values!", _state);
   ae_assert(ae_isfinite(eps, _state), "ODESolverRKCK: Eps is not finite!", _state);
   ae_assert(ae_fp_neq(eps, (double)(0)), "ODESolverRKCK: Eps is zero!", _state);
   ae_assert(ae_isfinite(h, _state), "ODESolverRKCK: H is not finite!", _state);
   odesolver_odesolverinit(0, y, n, x, m, eps, h, state, _state);
}

// This function provides a reverse communication interface, which is not documented or recommended for use.
// Instead, it is recommended that you use the better-documented API function odesolversolve() listed below.
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
// API: bool odesolveriteration(const odesolverstate &state, const xparams _xparams = xdefault);
// API: void odesolversolve(odesolverstate &state, void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr), void *ptr = NULL, const xparams _xparams = xdefault);
bool odesolveriteration(odesolverstate *state, ae_state *_state) {
   ae_int_t n;
   ae_int_t m;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double xc;
   double v;
   double h;
   double h2;
   bool gridpoint;
   double err;
   double maxgrowpow;
   ae_int_t klimit;
   bool result;

// Reverse communication preparations
// I know it looks ugly, but it works the same way
// anywhere from C++ to Python.
//
// This code initializes locals by:
// * random values determined during code
//   generation - on first subroutine call
// * values from previous call - on subsequent calls
   if (state->rstate.stage >= 0) {
      n = state->rstate.ia.ptr.p_int[0];
      m = state->rstate.ia.ptr.p_int[1];
      i = state->rstate.ia.ptr.p_int[2];
      j = state->rstate.ia.ptr.p_int[3];
      k = state->rstate.ia.ptr.p_int[4];
      klimit = state->rstate.ia.ptr.p_int[5];
      gridpoint = state->rstate.ba.ptr.p_bool[0];
      xc = state->rstate.ra.ptr.p_double[0];
      v = state->rstate.ra.ptr.p_double[1];
      h = state->rstate.ra.ptr.p_double[2];
      h2 = state->rstate.ra.ptr.p_double[3];
      err = state->rstate.ra.ptr.p_double[4];
      maxgrowpow = state->rstate.ra.ptr.p_double[5];
   } else {
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
   }
   if (state->rstate.stage == 0) {
      goto lbl_0;
   }
// Routine body

// prepare
   if (state->repterminationtype != 0) {
      result = false;
      return result;
   }
   n = state->n;
   m = state->m;
   h = state->h;
   maxgrowpow = ae_pow(odesolver_odesolvermaxgrow, (double)(5), _state);
   state->repnfev = 0;

// some preliminary checks for internal errors
// after this we assume that H>0 and M>1
   ae_assert(ae_fp_greater(state->h, (double)(0)), "ODESolver: internal error", _state);
   ae_assert(m > 1, "ODESolverIteration: internal error", _state);

// choose solver
   if (state->solvertype != 0) {
      goto lbl_1;
   }
// Cask-Karp solver
// Prepare coefficients table.
// Check it for errors
   ae_vector_set_length(&state->rka, 6, _state);
   state->rka.ptr.p_double[0] = (double)(0);
   state->rka.ptr.p_double[1] = (double)1 / (double)5;
   state->rka.ptr.p_double[2] = (double)3 / (double)10;
   state->rka.ptr.p_double[3] = (double)3 / (double)5;
   state->rka.ptr.p_double[4] = (double)(1);
   state->rka.ptr.p_double[5] = (double)7 / (double)8;
   ae_matrix_set_length(&state->rkb, 6, 5, _state);
   state->rkb.ptr.pp_double[1][0] = (double)1 / (double)5;
   state->rkb.ptr.pp_double[2][0] = (double)3 / (double)40;
   state->rkb.ptr.pp_double[2][1] = (double)9 / (double)40;
   state->rkb.ptr.pp_double[3][0] = (double)3 / (double)10;
   state->rkb.ptr.pp_double[3][1] = -(double)9 / (double)10;
   state->rkb.ptr.pp_double[3][2] = (double)6 / (double)5;
   state->rkb.ptr.pp_double[4][0] = -(double)11 / (double)54;
   state->rkb.ptr.pp_double[4][1] = (double)5 / (double)2;
   state->rkb.ptr.pp_double[4][2] = -(double)70 / (double)27;
   state->rkb.ptr.pp_double[4][3] = (double)35 / (double)27;
   state->rkb.ptr.pp_double[5][0] = (double)1631 / (double)55296;
   state->rkb.ptr.pp_double[5][1] = (double)175 / (double)512;
   state->rkb.ptr.pp_double[5][2] = (double)575 / (double)13824;
   state->rkb.ptr.pp_double[5][3] = (double)44275 / (double)110592;
   state->rkb.ptr.pp_double[5][4] = (double)253 / (double)4096;
   ae_vector_set_length(&state->rkc, 6, _state);
   state->rkc.ptr.p_double[0] = (double)37 / (double)378;
   state->rkc.ptr.p_double[1] = (double)(0);
   state->rkc.ptr.p_double[2] = (double)250 / (double)621;
   state->rkc.ptr.p_double[3] = (double)125 / (double)594;
   state->rkc.ptr.p_double[4] = (double)(0);
   state->rkc.ptr.p_double[5] = (double)512 / (double)1771;
   ae_vector_set_length(&state->rkcs, 6, _state);
   state->rkcs.ptr.p_double[0] = (double)2825 / (double)27648;
   state->rkcs.ptr.p_double[1] = (double)(0);
   state->rkcs.ptr.p_double[2] = (double)18575 / (double)48384;
   state->rkcs.ptr.p_double[3] = (double)13525 / (double)55296;
   state->rkcs.ptr.p_double[4] = (double)277 / (double)14336;
   state->rkcs.ptr.p_double[5] = (double)1 / (double)4;
   ae_matrix_set_length(&state->rkk, 6, n, _state);

// Main cycle consists of two iterations:
// * outer where we travel from X[i-1] to X[i]
// * inner where we travel inside [X[i-1],X[i]]
   ae_matrix_set_length(&state->ytbl, m, n, _state);
   ae_vector_set_length(&state->escale, n, _state);
   ae_vector_set_length(&state->yn, n, _state);
   ae_vector_set_length(&state->yns, n, _state);
   xc = state->xg.ptr.p_double[0];
   ae_v_move(&state->ytbl.ptr.pp_double[0][0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0, n - 1));
   for (j = 0; j <= n - 1; j++) {
      state->escale.ptr.p_double[j] = (double)(0);
   }
   i = 1;
lbl_3:
   if (i > m - 1) {
      goto lbl_5;
   }
// begin inner iteration
lbl_6:
   if (false) {
      goto lbl_7;
   }
// truncate step if needed (beyond right boundary).
// determine should we store X or not
   if (ae_fp_greater_eq(xc + h, state->xg.ptr.p_double[i])) {
      h = state->xg.ptr.p_double[i] - xc;
      gridpoint = true;
   } else {
      gridpoint = false;
   }

// Update error scale maximums
//
// These maximums are initialized by zeros,
// then updated every iterations.
   for (j = 0; j <= n - 1; j++) {
      state->escale.ptr.p_double[j] = ae_maxreal(state->escale.ptr.p_double[j], ae_fabs(state->yc.ptr.p_double[j], _state), _state);
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
   ae_v_move(&state->yn.ptr.p_double[0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0, n - 1));
   ae_v_move(&state->yns.ptr.p_double[0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0, n - 1));
   k = 0;
lbl_8:
   if (k > 5) {
      goto lbl_10;
   }
// prepare data for the next update of YN/YNS
   state->x = state->xscale * (xc + state->rka.ptr.p_double[k] * h);
   ae_v_move(&state->y.ptr.p_double[0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0, n - 1));
   for (j = 0; j <= k - 1; j++) {
      v = state->rkb.ptr.pp_double[k][j];
      ae_v_addd(&state->y.ptr.p_double[0], 1, &state->rkk.ptr.pp_double[j][0], 1, ae_v_len(0, n - 1), v);
   }
   state->needdy = true;
   state->rstate.stage = 0;
   goto lbl_rcomm;
lbl_0:
   state->needdy = false;
   state->repnfev = state->repnfev + 1;
   v = h * state->xscale;
   ae_v_moved(&state->rkk.ptr.pp_double[k][0], 1, &state->dy.ptr.p_double[0], 1, ae_v_len(0, n - 1), v);

// update YN/YNS
   v = state->rkc.ptr.p_double[k];
   ae_v_addd(&state->yn.ptr.p_double[0], 1, &state->rkk.ptr.pp_double[k][0], 1, ae_v_len(0, n - 1), v);
   v = state->rkcs.ptr.p_double[k];
   ae_v_addd(&state->yns.ptr.p_double[0], 1, &state->rkk.ptr.pp_double[k][0], 1, ae_v_len(0, n - 1), v);
   k = k + 1;
   goto lbl_8;
lbl_10:

// estimate error
   err = (double)(0);
   for (j = 0; j <= n - 1; j++) {
      if (!state->fraceps) {

      // absolute error is estimated
         err = ae_maxreal(err, ae_fabs(state->yn.ptr.p_double[j] - state->yns.ptr.p_double[j], _state), _state);
      } else {

      // Relative error is estimated
         v = state->escale.ptr.p_double[j];
         if (ae_fp_eq(v, (double)(0))) {
            v = (double)(1);
         }
         err = ae_maxreal(err, ae_fabs(state->yn.ptr.p_double[j] - state->yns.ptr.p_double[j], _state) / v, _state);
      }
   }

// calculate new step, restart if necessary
   if (ae_fp_less_eq(maxgrowpow * err, state->eps)) {
      h2 = odesolver_odesolvermaxgrow * h;
   } else {
      h2 = h * ae_pow(state->eps / err, 0.2, _state);
   }
   if (ae_fp_less(h2, h / odesolver_odesolvermaxshrink)) {
      h2 = h / odesolver_odesolvermaxshrink;
   }
   if (ae_fp_greater(err, state->eps)) {
      h = ae_minreal(h2, odesolver_odesolverguaranteeddecay * h, _state);
      goto lbl_6;
   }
// advance position
   xc = xc + h;
   ae_v_move(&state->yc.ptr.p_double[0], 1, &state->yn.ptr.p_double[0], 1, ae_v_len(0, n - 1));

// update H
   h = h2;

// break on grid point
   if (gridpoint) {
      goto lbl_7;
   }
   goto lbl_6;
lbl_7:

// save result
   ae_v_move(&state->ytbl.ptr.pp_double[i][0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0, n - 1));
   i = i + 1;
   goto lbl_3;
lbl_5:
   state->repterminationtype = 1;
   result = false;
   return result;
lbl_1:
   result = false;
   return result;

// Saving state
lbl_rcomm:
   result = true;
   state->rstate.ia.ptr.p_int[0] = n;
   state->rstate.ia.ptr.p_int[1] = m;
   state->rstate.ia.ptr.p_int[2] = i;
   state->rstate.ia.ptr.p_int[3] = j;
   state->rstate.ia.ptr.p_int[4] = k;
   state->rstate.ia.ptr.p_int[5] = klimit;
   state->rstate.ba.ptr.p_bool[0] = gridpoint;
   state->rstate.ra.ptr.p_double[0] = xc;
   state->rstate.ra.ptr.p_double[1] = v;
   state->rstate.ra.ptr.p_double[2] = h;
   state->rstate.ra.ptr.p_double[3] = h2;
   state->rstate.ra.ptr.p_double[4] = err;
   state->rstate.ra.ptr.p_double[5] = maxgrowpow;
   return result;
}

// ODE solver results
//
// Called after OdeSolverIteration returned False.
//
// Inputs:
//     State   -   algorithm state (used by OdeSolverIteration).
//
// Outputs:
//     M       -   number of tabulated values, M>=1
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
// API: void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep, const xparams _xparams = xdefault);
void odesolverresults(odesolverstate *state, ae_int_t *m, RVector *xtbl, RMatrix *ytbl, odesolverreport *rep, ae_state *_state) {
   double v;
   ae_int_t i;

   *m = 0;
   ae_vector_free(xtbl, true);
   ae_matrix_free(ytbl, true);
   odesolverreport_free(rep, true);

   rep->terminationtype = state->repterminationtype;
   if (rep->terminationtype > 0) {
      *m = state->m;
      rep->nfev = state->repnfev;
      ae_vector_set_length(xtbl, state->m, _state);
      v = state->xscale;
      ae_v_moved(&xtbl->ptr.p_double[0], 1, &state->xg.ptr.p_double[0], 1, ae_v_len(0, state->m - 1), v);
      ae_matrix_set_length(ytbl, state->m, state->n, _state);
      for (i = 0; i <= state->m - 1; i++) {
         ae_v_move(&ytbl->ptr.pp_double[i][0], 1, &state->ytbl.ptr.pp_double[i][0], 1, ae_v_len(0, state->n - 1));
      }
   } else {
      rep->nfev = 0;
   }
}

// Internal initialization subroutine
static void odesolver_odesolverinit(ae_int_t solvertype, RVector *y, ae_int_t n, RVector *x, ae_int_t m, double eps, double h, odesolverstate *state, ae_state *_state) {
   ae_int_t i;
   double v;

   odesolverstate_free(state, true);

// Prepare RComm
   ae_vector_set_length(&state->rstate.ia, 5 + 1, _state);
   ae_vector_set_length(&state->rstate.ba, 0 + 1, _state);
   ae_vector_set_length(&state->rstate.ra, 5 + 1, _state);
   state->rstate.stage = -1;
   state->needdy = false;

// check parameters.
   if ((n <= 0 || m < 1) || ae_fp_eq(eps, (double)(0))) {
      state->repterminationtype = -1;
      return;
   }
   if (ae_fp_less(h, (double)(0))) {
      h = -h;
   }
// quick exit if necessary.
// after this block we assume that M>1
   if (m == 1) {
      state->repnfev = 0;
      state->repterminationtype = 1;
      ae_matrix_set_length(&state->ytbl, 1, n, _state);
      ae_v_move(&state->ytbl.ptr.pp_double[0][0], 1, &y->ptr.p_double[0], 1, ae_v_len(0, n - 1));
      ae_vector_set_length(&state->xg, m, _state);
      ae_v_move(&state->xg.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0, m - 1));
      return;
   }
// check again: correct order of X[]
   if (ae_fp_eq(x->ptr.p_double[1], x->ptr.p_double[0])) {
      state->repterminationtype = -2;
      return;
   }
   for (i = 1; i <= m - 1; i++) {
      if ((ae_fp_greater(x->ptr.p_double[1], x->ptr.p_double[0]) && ae_fp_less_eq(x->ptr.p_double[i], x->ptr.p_double[i - 1])) || (ae_fp_less(x->ptr.p_double[1], x->ptr.p_double[0]) && ae_fp_greater_eq(x->ptr.p_double[i], x->ptr.p_double[i - 1]))) {
         state->repterminationtype = -2;
         return;
      }
   }

// auto-select H if necessary
   if (ae_fp_eq(h, (double)(0))) {
      v = ae_fabs(x->ptr.p_double[1] - x->ptr.p_double[0], _state);
      for (i = 2; i <= m - 1; i++) {
         v = ae_minreal(v, ae_fabs(x->ptr.p_double[i] - x->ptr.p_double[i - 1], _state), _state);
      }
      h = 0.001 * v;
   }
// store parameters
   state->n = n;
   state->m = m;
   state->h = h;
   state->eps = ae_fabs(eps, _state);
   state->fraceps = ae_fp_less(eps, (double)(0));
   ae_vector_set_length(&state->xg, m, _state);
   ae_v_move(&state->xg.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0, m - 1));
   if (ae_fp_greater(x->ptr.p_double[1], x->ptr.p_double[0])) {
      state->xscale = (double)(1);
   } else {
      state->xscale = (double)(-1);
      ae_v_muld(&state->xg.ptr.p_double[0], 1, ae_v_len(0, m - 1), -1);
   }
   ae_vector_set_length(&state->yc, n, _state);
   ae_v_move(&state->yc.ptr.p_double[0], 1, &y->ptr.p_double[0], 1, ae_v_len(0, n - 1));
   state->solvertype = solvertype;
   state->repterminationtype = 0;

// Allocate arrays
   ae_vector_set_length(&state->y, n, _state);
   ae_vector_set_length(&state->dy, n, _state);
}

void odesolverstate_init(void *_p, ae_state *_state, bool make_automatic) {
   odesolverstate *p = (odesolverstate *) _p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->yc, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->escale, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->xg, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->dy, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->ytbl, 0, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->yn, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->yns, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->rka, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->rkc, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->rkcs, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->rkb, 0, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->rkk, 0, 0, DT_REAL, _state, make_automatic);
   rcommstate_init(&p->rstate, _state, make_automatic);
}

void odesolverstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   odesolverstate *dst = (odesolverstate *) _dst;
   odesolverstate *src = (odesolverstate *) _src;
   dst->n = src->n;
   dst->m = src->m;
   dst->xscale = src->xscale;
   dst->h = src->h;
   dst->eps = src->eps;
   dst->fraceps = src->fraceps;
   ae_vector_copy(&dst->yc, &src->yc, _state, make_automatic);
   ae_vector_copy(&dst->escale, &src->escale, _state, make_automatic);
   ae_vector_copy(&dst->xg, &src->xg, _state, make_automatic);
   dst->solvertype = src->solvertype;
   dst->needdy = src->needdy;
   dst->x = src->x;
   ae_vector_copy(&dst->y, &src->y, _state, make_automatic);
   ae_vector_copy(&dst->dy, &src->dy, _state, make_automatic);
   ae_matrix_copy(&dst->ytbl, &src->ytbl, _state, make_automatic);
   dst->repterminationtype = src->repterminationtype;
   dst->repnfev = src->repnfev;
   ae_vector_copy(&dst->yn, &src->yn, _state, make_automatic);
   ae_vector_copy(&dst->yns, &src->yns, _state, make_automatic);
   ae_vector_copy(&dst->rka, &src->rka, _state, make_automatic);
   ae_vector_copy(&dst->rkc, &src->rkc, _state, make_automatic);
   ae_vector_copy(&dst->rkcs, &src->rkcs, _state, make_automatic);
   ae_matrix_copy(&dst->rkb, &src->rkb, _state, make_automatic);
   ae_matrix_copy(&dst->rkk, &src->rkk, _state, make_automatic);
   rcommstate_copy(&dst->rstate, &src->rstate, _state, make_automatic);
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
   rcommstate_free(&p->rstate, make_automatic);
}

void odesolverreport_init(void *_p, ae_state *_state, bool make_automatic) {
   odesolverreport *p = (odesolverreport *) _p;
   ae_touch_ptr((void *)p);
}

void odesolverreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
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
//
DefClass(odesolverstate, DecVal(needdy) DecVar(y) DecVar(dy) DecVal(x))

//
DefClass(odesolverreport, DecVal(nfev) DecVal(terminationtype))

void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state, const xparams _xparams) {
   jmp_buf _break_jump;
   alglib_impl::ae_state _alglib_env_state;
   alglib_impl::ae_state_init(&_alglib_env_state);
   if (setjmp(_break_jump)) {
#if !defined AE_NO_EXCEPTIONS
      _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
      _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
      return;
#endif
   }
   ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
   if (_xparams.flags != 0x0)
      ae_state_set_flags(&_alglib_env_state, _xparams.flags);
   alglib_impl::odesolverrkck(const_cast < alglib_impl::ae_vector * >(y.c_ptr()), n, const_cast < alglib_impl::ae_vector * >(x.c_ptr()), m, eps, h, const_cast < alglib_impl::odesolverstate * >(state.c_ptr()), &_alglib_env_state);
   alglib_impl::ae_state_clear(&_alglib_env_state);
   return;
}
#if !defined AE_NO_EXCEPTIONS
void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state, const xparams _xparams) {
   jmp_buf _break_jump;
   alglib_impl::ae_state _alglib_env_state;
   ae_int_t n;
   ae_int_t m;

   n = y.length();
   m = x.length();
   alglib_impl::ae_state_init(&_alglib_env_state);
   if (setjmp(_break_jump))
      _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
   ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
   if (_xparams.flags != 0x0)
      ae_state_set_flags(&_alglib_env_state, _xparams.flags);
   alglib_impl::odesolverrkck(const_cast < alglib_impl::ae_vector * >(y.c_ptr()), n, const_cast < alglib_impl::ae_vector * >(x.c_ptr()), m, eps, h, const_cast < alglib_impl::odesolverstate * >(state.c_ptr()), &_alglib_env_state);

   alglib_impl::ae_state_clear(&_alglib_env_state);
   return;
}
#endif

bool odesolveriteration(const odesolverstate &state, const xparams _xparams) {
   jmp_buf _break_jump;
   alglib_impl::ae_state _alglib_env_state;
   alglib_impl::ae_state_init(&_alglib_env_state);
   if (setjmp(_break_jump)) {
#if !defined AE_NO_EXCEPTIONS
      _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
      _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
      return 0;
#endif
   }
   ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
   if (_xparams.flags != 0x0)
      ae_state_set_flags(&_alglib_env_state, _xparams.flags);
   bool result = alglib_impl::odesolveriteration(const_cast < alglib_impl::odesolverstate * >(state.c_ptr()), &_alglib_env_state);
   alglib_impl::ae_state_clear(&_alglib_env_state);
   return *(reinterpret_cast < bool *>(&result));
}

// This function is used to launcn iterations of ODE solver
//
// It accepts following parameters:
//     diff    -   callback which calculates dy/dx for given y and x
//     ptr     -   optional pointer which is passed to diff; can be NULL
// ALGLIB: Copyright 01.09.2009 by Sergey Bochkanov
void odesolversolve(odesolverstate &state, void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr), void *ptr, const xparams _xparams) {
   jmp_buf _break_jump;
   alglib_impl::ae_state _alglib_env_state;
   alglib_impl::ae_state_init(&_alglib_env_state);
   if (setjmp(_break_jump)) {
#if !defined AE_NO_EXCEPTIONS
      _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
      _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
      return;
#endif
   }
   ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
   if (_xparams.flags != 0x0)
      ae_state_set_flags(&_alglib_env_state, _xparams.flags);
   alglib_impl::ae_assert(diff != NULL, "ALGLIB: error in 'odesolversolve()' (diff is NULL)", &_alglib_env_state);
   while (alglib_impl::odesolveriteration(state.c_ptr(), &_alglib_env_state)) {
      _ALGLIB_CALLBACK_EXCEPTION_GUARD_BEGIN if (state.needdy) {
         diff(state.y, state.x, state.dy, ptr);
         continue;
      }
      goto lbl_no_callback;
      _ALGLIB_CALLBACK_EXCEPTION_GUARD_END lbl_no_callback:alglib_impl::ae_assert(false, "ALGLIB: unexpected error in 'odesolversolve'", &_alglib_env_state);
   }
   alglib_impl::ae_state_clear(&_alglib_env_state);
}

void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep, const xparams _xparams) {
   jmp_buf _break_jump;
   alglib_impl::ae_state _alglib_env_state;
   alglib_impl::ae_state_init(&_alglib_env_state);
   if (setjmp(_break_jump)) {
#if !defined AE_NO_EXCEPTIONS
      _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
      _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
      return;
#endif
   }
   ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
   if (_xparams.flags != 0x0)
      ae_state_set_flags(&_alglib_env_state, _xparams.flags);
   alglib_impl::odesolverresults(const_cast < alglib_impl::odesolverstate * >(state.c_ptr()), &m, const_cast < alglib_impl::ae_vector * >(xtbl.c_ptr()), const_cast < alglib_impl::ae_matrix * >(ytbl.c_ptr()), const_cast < alglib_impl::odesolverreport * >(rep.c_ptr()), &_alglib_env_state);
   alglib_impl::ae_state_clear(&_alglib_env_state);
   return;
}
} // end of namespace alglib
