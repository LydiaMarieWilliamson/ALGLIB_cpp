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
#include "Solvers.h"

// === POLYNOMIALSOLVER Package ===
// Depends on: (LinAlg) TRFAC, EVD
namespace alglib_impl {
// Polynomial root finding.
//
// This function returns all roots of the polynomial
//     P(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
// Both real and complex roots are returned (see below).
//
// Inputs:
//     A       -   array[N+1], polynomial coefficients:
//                 * A[0] is constant term
//                 * A[N] is a coefficient of X^N
//     N       -   polynomial degree
//
// Outputs:
//     X       -   array of complex roots:
//                 * for isolated real root, X[I] is strictly real: IMAGE(X[I])=0
//                 * complex roots are always returned in pairs - roots occupy
//                   positions I and I+1, with:
//                   * X[I+1]=Conj(X[I])
//                   * IMAGE(X[I]) > 0
//                   * IMAGE(X[I+1]) = -IMAGE(X[I]) < 0
//                 * multiple real roots may have non-zero imaginary part due
//                   to roundoff errors. There is no reliable way to distinguish
//                   real root of multiplicity 2 from two  complex  roots  in
//                   the presence of roundoff errors.
//     Rep     -   report, additional information, following fields are set:
//                 * Rep.MaxErr - max( |P(xi)| )  for  i=0..N-1.  This  field
//                   allows to quickly estimate "quality" of the roots  being
//                   returned.
//
// NOTE:   this function uses companion matrix method to find roots. In  case
//         internal EVD  solver  fails  do  find  eigenvalues,  exception  is
//         generated.
//
// NOTE:   roots are not "polished" and  no  matrix  balancing  is  performed
//         for them.
// ALGLIB: Copyright 24.02.2014 by Sergey Bochkanov
// API: void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep);
void polynomialsolve(RVector *a, ae_int_t n, CVector *x, polynomialsolverreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   bool status;
   ae_int_t nz;
   ae_int_t ne;
   complex v;
   complex vv;
   ae_frame_make(&_frame_block);
   DupVector(a);
   SetVector(x);
   SetObj(polynomialsolverreport, rep);
   NewMatrix(c, 0, 0, DT_REAL);
   NewMatrix(vl, 0, 0, DT_REAL);
   NewMatrix(vr, 0, 0, DT_REAL);
   NewVector(wr, 0, DT_REAL);
   NewVector(wi, 0, DT_REAL);
   ae_assert(n > 0, "PolynomialSolve: N <= 0");
   ae_assert(a->cnt >= n + 1, "PolynomialSolve: Length(A)<N+1");
   ae_assert(isfinitevector(a, n + 1), "PolynomialSolve: A contains infitite numbers");
   ae_assert(a->xR[n] != 0.0, "PolynomialSolve: A[N]=0");
// Prepare
   ae_vector_set_length(x, n);
// Normalize A:
// * analytically determine NZ zero roots
// * quick exit for NZ=N
// * make residual NE-th degree polynomial monic
//   (here NE=N-NZ)
   nz = 0;
   while (nz < n && a->xR[nz] == 0.0) {
      nz++;
   }
   ne = n - nz;
   for (i = nz; i <= n; i++) {
      a->xR[i - nz] = a->xR[i] / a->xR[n];
   }
// For NZ<N, build companion matrix and find NE non-zero roots
   if (ne > 0) {
      ae_matrix_set_length(&c, ne, ne);
      for (i = 0; i < ne; i++) {
         for (j = 0; j < ne; j++) {
            c.xyR[i][j] = 0.0;
         }
      }
      c.xyR[0][ne - 1] = -a->xR[0];
      for (i = 1; i < ne; i++) {
         c.xyR[i][i - 1] = 1.0;
         c.xyR[i][ne - 1] = -a->xR[i];
      }
      status = rmatrixevd(&c, ne, 0, &wr, &wi, &vl, &vr);
      ae_assert(status, "PolynomialSolve: inernal error - EVD solver failed");
      for (i = 0; i < ne; i++) {
         x->xC[i] = ae_complex_from_d(wr.xR[i], wi.xR[i]);
      }
   }
// Remaining NZ zero roots
   for (i = ne; i < n; i++) {
      x->xC[i] = ae_complex_from_i(0);
   }
// Rep
   rep->maxerr = 0.0;
   for (i = 0; i < ne; i++) {
      v = ae_complex_from_i(0);
      vv = ae_complex_from_i(1);
      for (j = 0; j <= ne; j++) {
         v = ae_c_add(v, ae_c_mul_d(vv, a->xR[j]));
         vv = ae_c_mul(vv, x->xC[i]);
      }
      rep->maxerr = rmax2(rep->maxerr, abscomplex(v));
   }
   ae_frame_leave();
}

void polynomialsolverreport_init(void *_p, bool make_automatic) {
}

void polynomialsolverreport_copy(void *_dst, void *_src, bool make_automatic) {
   polynomialsolverreport *dst = (polynomialsolverreport *)_dst;
   polynomialsolverreport *src = (polynomialsolverreport *)_src;
   dst->maxerr = src->maxerr;
}

void polynomialsolverreport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
DefClass(polynomialsolverreport, AndD DecVal(maxerr))

void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::polynomialsolve(ConstT(ae_vector, a), n, ConstT(ae_vector, x), ConstT(polynomialsolverreport, rep));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === DIRECTDENSESOLVERS Package ===
// Depends on: (AlgLibInternal) XBLAS
// Depends on: (LinAlg) RCOND, SVD
namespace alglib_impl {
// Dense solver for A*x=b with N*N real matrix A and N*1 real vectorx  x  and
// b. This is "slow-but-feature rich" version of the  linear  solver.  Faster
// version is RMatrixSolveFast() function.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(N^3) complexity
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear  system
//            ! and  performs  iterative   refinement,   which   results   in
//            ! significant performance penalty  when  compared  with  "fast"
//            ! version  which  just  performs  LU  decomposition  and  calls
//            ! triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations. It is also very significant on small matrices.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, RMatrixSolveFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void rmatrixsolve(RMatrix *a, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_REAL);
   NewMatrix(xm, 0, 0, DT_REAL);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_move(bm.xyR[0], bm.stride, b->xR, 1, n);
   rmatrixsolvem(a, n, &bm, 1, true, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_move(x->xR, 1, xm.xyR[0], xm.stride, n);
   ae_frame_leave();
}

// Basic LU solver for PLU*x = y.
//
// This subroutine assumes that:
// * A=PLU is well-conditioned, so no zero divisions or overflow may occur
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_rbasiclusolve(RMatrix *lua, ZVector *p, ae_int_t n, RVector *xb) {
   ae_int_t i;
   double v;
   for (i = 0; i < n; i++) {
      if (p->xZ[i] != i) {
         swapr(&xb->xR[i], &xb->xR[p->xZ[i]]);
      }
   }
   for (i = 1; i < n; i++) {
      v = ae_v_dotproduct(lua->xyR[i], 1, xb->xR, 1, i);
      xb->xR[i] -= v;
   }
   xb->xR[n - 1] /= lua->xyR[n - 1][n - 1];
   for (i = n - 2; i >= 0; i--) {
      v = ae_v_dotproduct(&lua->xyR[i][i + 1], 1, &xb->xR[i + 1], 1, n - i - 1);
      xb->xR[i] = (xb->xR[i] - v) / lua->xyR[i][i];
   }
}

// Dense solver.
//
// This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
// real matrix, x  and  b  are  vectors.  This is a "fast" version of  linear
// solver which does NOT provide  any  additional  functions  like  condition
// number estimation or iterative refinement.
//
// Algorithm features:
// * efficient algorithm O(N^3) complexity
// * no performance overhead from additional functionality
//
// If you need condition number estimation or iterative refinement, use  more
// feature-rich version - RMatrixSolve().
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 16.03.2015 by Sergey Bochkanov
// API: void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
void rmatrixsolvefast(RMatrix *a, ae_int_t n, RVector *b, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   NewVector(p, 0, DT_INT);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   rmatrixlu(a, n, n, &p);
   for (i = 0; i < n; i++) {
      if (a->xyR[i][i] == 0.0) {
         for (j = 0; j < n; j++) {
            b->xR[j] = 0.0;
         }
         *info = -3;
         ae_frame_leave();
         return;
      }
   }
   directdensesolvers_rbasiclusolve(a, &p, n, b);
   *info = 1;
   ae_frame_leave();
}

// Internal subroutine.
// Returns maximum count of RFS iterations as function of:
// 1. machine epsilon
// 2. task size.
// 3. condition number
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static ae_int_t directdensesolvers_densesolverrfsmax(ae_int_t n, double r1, double rinf) {
   ae_int_t result;
   result = 5;
   return result;
}

// Internal LU solver
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_rmatrixlusolveinternal(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *a, bool havea, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t rfs;
   ae_int_t nrfs;
   double v;
   double verr;
   double mxb;
   bool smallerr;
   bool terminatenexttime;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewVector(xc, 0, DT_REAL);
   NewVector(y, 0, DT_REAL);
   NewVector(bc, 0, DT_REAL);
   NewVector(xa, 0, DT_REAL);
   NewVector(xb, 0, DT_REAL);
   NewVector(tx, 0, DT_REAL);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   for (i = 0; i < n; i++) {
      if (p->xZ[i] > n - 1 || p->xZ[i] < i) {
         *info = -1;
         ae_frame_leave();
         return;
      }
   }
   ae_matrix_set_length(x, n, m);
   ae_vector_set_length(&y, n);
   ae_vector_set_length(&xc, n);
   ae_vector_set_length(&bc, n);
   ae_vector_set_length(&tx, n + 1);
   ae_vector_set_length(&xa, n + 1);
   ae_vector_set_length(&xb, n + 1);
// estimate condition number, test for near singularity
   rep->r1 = rmatrixlurcond1(lua, n);
   rep->rinf = rmatrixlurcondinf(lua, n);
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            x->xyR[i][j] = 0.0;
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
   *info = 1;
// First stage of solution: rough solution with TRSM()
   mxb = 0.0;
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
         v = b->xyR[i][j];
         mxb = rmax2(mxb, fabs(v));
         x->xyR[i][j] = v;
      }
   }
   for (i = 0; i < n; i++) {
      if (p->xZ[i] != i) {
         for (j = 0; j < m; j++) {
            swapr(&x->xyR[i][j], &x->xyR[p->xZ[i]][j]);
         }
      }
   }
   rmatrixlefttrsm(n, m, lua, 0, 0, false, true, 0, x, 0, 0);
   rmatrixlefttrsm(n, m, lua, 0, 0, true, false, 0, x, 0, 0);
// Second stage: iterative refinement
   if (havea) {
      for (k = 0; k < m; k++) {
         nrfs = directdensesolvers_densesolverrfsmax(n, rep->r1, rep->rinf);
         terminatenexttime = false;
         for (rfs = 0; rfs < nrfs; rfs++) {
            if (terminatenexttime) {
               break;
            }
         // generate right part
            smallerr = true;
            ae_v_move(xb.xR, 1, &x->xyR[0][k], x->stride, n);
            for (i = 0; i < n; i++) {
               ae_v_move(xa.xR, 1, a->xyR[i], 1, n);
               xa.xR[n] = -1.0;
               xb.xR[n] = b->xyR[i][k];
               xdot(&xa, &xb, n + 1, &tx, &v, &verr);
               y.xR[i] = -v;
               smallerr = smallerr && SmallR(v, 4 * verr);
            }
            if (smallerr) {
               terminatenexttime = true;
            }
         // solve and update
            directdensesolvers_rbasiclusolve(lua, p, n, &y);
            ae_v_add(&x->xyR[0][k], x->stride, y.xR, 1, n);
         }
      }
   }
   ae_frame_leave();
}

// Dense solver.
//
// Similar to RMatrixSolve() but solves task with multiple right parts (where
// b and x are NxM matrices). This is  "slow-but-robust"  version  of  linear
// solver with additional functionality  like  condition  number  estimation.
// There also exists faster version - RMatrixSolveMFast().
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * optional iterative refinement
// * O(N^3+M*N^2) complexity
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear  system
//            ! and  performs  iterative   refinement,   which   results   in
//            ! significant performance penalty  when  compared  with  "fast"
//            ! version  which  just  performs  LU  decomposition  and  calls
//            ! triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations. It also very significant on small matrices.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, RMatrixSolveMFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//     RFS     -   iterative refinement switch:
//                 * True - refinement is used.
//                   Less performance, more precision.
//                 * False - refinement is not used.
//                   More performance, less precision.
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is ill conditioned or singular.
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void rmatrixsolvem(RMatrix *a, ae_int_t n, RMatrix *b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(da, 0, 0, DT_REAL);
   NewMatrix(emptya, 0, 0, DT_REAL);
   NewVector(p, 0, DT_INT);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&da, n, n);
// 1. factorize matrix
// 3. solve
   for (i = 0; i < n; i++) {
      ae_v_move(da.xyR[i], 1, a->xyR[i], 1, n);
   }
   rmatrixlu(&da, n, n, &p);
   if (rfs) {
      directdensesolvers_rmatrixlusolveinternal(&da, &p, n, a, true, b, m, info, rep, x);
   } else {
      directdensesolvers_rmatrixlusolveinternal(&da, &p, n, &emptya, false, b, m, info, rep, x);
   }
   ae_frame_leave();
}

// Dense solver.
//
// Similar to RMatrixSolve() but solves task with multiple right parts (where
// b and x are NxM matrices). This is "fast" version of linear  solver  which
// does NOT offer additional functions like condition  number  estimation  or
// iterative refinement.
//
// Algorithm features:
// * O(N^3+M*N^2) complexity
// * no additional functionality, highest performance
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//     RFS     -   iterative refinement switch:
//                 * True - refinement is used.
//                   Less performance, more precision.
//                 * False - refinement is not used.
//                   More performance, less precision.
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     B       -   array[N]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void rmatrixsolvemfast(RMatrix *a, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   NewVector(p, 0, DT_INT);
// Check for exact degeneracy
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   rmatrixlu(a, n, n, &p);
   for (i = 0; i < n; i++) {
      if (a->xyR[i][i] == 0.0) {
         for (j = 0; j < n; j++) {
            for (k = 0; k < m; k++) {
               b->xyR[j][k] = 0.0;
            }
         }
         *info = -3;
         ae_frame_leave();
         return;
      }
   }
// Solve with TRSM()
   for (i = 0; i < n; i++) {
      if (p.xZ[i] != i) {
         for (j = 0; j < m; j++) {
            swapr(&b->xyR[i][j], &b->xyR[p.xZ[i]][j]);
         }
      }
   }
   rmatrixlefttrsm(n, m, a, 0, 0, false, true, 0, b, 0, 0);
   rmatrixlefttrsm(n, m, a, 0, 0, true, false, 0, b, 0, 0);
   *info = 1;
   ae_frame_leave();
}

// Dense solver.
//
// This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
// real matrix given by its LU decomposition, x and b are real vectors.  This
// is "slow-but-robust" version of the linear LU-based solver. Faster version
// is RMatrixLUSolveFast() function.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(N^2) complexity
// * condition number estimation
//
// No iterative refinement  is provided because exact form of original matrix
// is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which results in 10-15x  performance  penalty  when  compared
//            ! with "fast" version which just calls triangular solver.
//            !
//            ! This performance penalty is insignificant  when compared with
//            ! cost of large LU decomposition.  However,  if you  call  this
//            ! function many times for the same  left  side,  this  overhead
//            ! BECOMES significant. It  also  becomes significant for small-
//            ! scale problems.
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! RMatrixLUSolveFast() function.
//
// Inputs:
//     LUA     -   array[N,N], LU decomposition, RMatrixLU result
//     P       -   array[N], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[N], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void rmatrixlusolve(RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_REAL);
   NewMatrix(xm, 0, 0, DT_REAL);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_move(bm.xyR[0], bm.stride, b->xR, 1, n);
   rmatrixlusolvem(lua, p, n, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_move(x->xR, 1, xm.xyR[0], xm.stride, n);
   ae_frame_leave();
}

// Dense solver.
//
// This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
// real matrix given by its LU decomposition, x and b are real vectors.  This
// is "fast-without-any-checks" version of the linear LU-based solver. Slower
// but more robust version is RMatrixLUSolve() function.
//
// Algorithm features:
// * O(N^2) complexity
// * fast algorithm without ANY additional checks, just triangular solver
//
// Inputs:
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
//     P       -   array[0..N-1], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
// API: void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
void rmatrixlusolvefast(RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info) {
   ae_int_t i;
   ae_int_t j;
   *info = 0;
   if (n <= 0) {
      *info = -1;
      return;
   }
   for (i = 0; i < n; i++) {
      if (lua->xyR[i][i] == 0.0) {
         for (j = 0; j < n; j++) {
            b->xR[j] = 0.0;
         }
         *info = -3;
         return;
      }
   }
   directdensesolvers_rbasiclusolve(lua, p, n, b);
   *info = 1;
}

// Dense solver.
//
// Similar to RMatrixLUSolve() but solves  task  with  multiple  right  parts
// (where b and x are NxM matrices). This  is  "robust-but-slow"  version  of
// LU-based solver which performs additional  checks  for  non-degeneracy  of
// inputs (condition number estimation). If you need  best  performance,  use
// "fast-without-any-checks" version, RMatrixLUSolveMFast().
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(M*N^2) complexity
// * condition number estimation
//
// No iterative refinement  is provided because exact form of original matrix
// is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant  performance   penalty   when
//            ! compared with "fast"  version  which  just  calls  triangular
//            ! solver.
//            !
//            ! This performance penalty is especially apparent when you  use
//            ! ALGLIB parallel capabilities (condition number estimation  is
//            ! inherently  sequential).  It   also   becomes significant for
//            ! small-scale problems.
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! RMatrixLUSolveMFast() function.
//
// Inputs:
//     LUA     -   array[N,N], LU decomposition, RMatrixLU result
//     P       -   array[N], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void rmatrixlusolvem(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(emptya, 0, 0, DT_REAL);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
// solve
   directdensesolvers_rmatrixlusolveinternal(lua, p, n, &emptya, false, b, m, info, rep, x);
   ae_frame_leave();
}

// Dense solver.
//
// Similar to RMatrixLUSolve() but solves  task  with  multiple  right parts,
// where b and x are NxM matrices.  This is "fast-without-any-checks" version
// of LU-based solver. It does not estimate  condition number  of  a  system,
// so it is extremely fast. If you need better detection  of  near-degenerate
// cases, use RMatrixLUSolveM() function.
//
// Algorithm features:
// * O(M*N^2) complexity
// * fast algorithm without ANY additional checks, just triangular solver
//
// Inputs:
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
//     P       -   array[0..N-1], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N,M]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
// API: void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void rmatrixlusolvemfast(RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   *info = 0;
// Check for exact degeneracy
   if (n <= 0 || m <= 0) {
      *info = -1;
      return;
   }
   for (i = 0; i < n; i++) {
      if (lua->xyR[i][i] == 0.0) {
         for (j = 0; j < n; j++) {
            for (k = 0; k < m; k++) {
               b->xyR[j][k] = 0.0;
            }
         }
         *info = -3;
         return;
      }
   }
// Solve with TRSM()
   for (i = 0; i < n; i++) {
      if (p->xZ[i] != i) {
         for (j = 0; j < m; j++) {
            swapr(&b->xyR[i][j], &b->xyR[p->xZ[i]][j]);
         }
      }
   }
   rmatrixlefttrsm(n, m, lua, 0, 0, false, true, 0, b, 0, 0);
   rmatrixlefttrsm(n, m, lua, 0, 0, true, false, 0, b, 0, 0);
   *info = 1;
}

// Dense solver.
//
// This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
// LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
// both A and its LU decomposition.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(N^2) complexity
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
//     P       -   array[0..N-1], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void rmatrixmixedsolve(RMatrix *a, RMatrix *lua, ZVector *p, ae_int_t n, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_REAL);
   NewMatrix(xm, 0, 0, DT_REAL);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_move(bm.xyR[0], bm.stride, b->xR, 1, n);
   rmatrixmixedsolvem(a, lua, p, n, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_move(x->xR, 1, xm.xyR[0], xm.stride, n);
   ae_frame_leave();
}

// Dense solver.
//
// Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
// (where b and x are NxM matrices).
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(M*N^2) complexity
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
//     P       -   array[0..N-1], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void rmatrixmixedsolvem(RMatrix *a, RMatrix *lua, ZVector *p, ae_int_t n, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      return;
   }
// solve
   directdensesolvers_rmatrixlusolveinternal(lua, p, n, a, true, b, m, info, rep, x);
}

// Basic LU solver for ScaleA*PLU*x = y.
//
// This subroutine assumes that:
// * L is well-scaled, and it is U which needs scaling by ScaleA.
// * A=PLU is well-conditioned, so no zero divisions or overflow may occur
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_cbasiclusolve(CMatrix *lua, ZVector *p, ae_int_t n, CVector *xb) {
   ae_int_t i;
   complex v;
   for (i = 0; i < n; i++) {
      if (p->xZ[i] != i) {
         swapc(&xb->xC[i], &xb->xC[p->xZ[i]]);
      }
   }
   for (i = 1; i < n; i++) {
      v = ae_v_cdotproduct(lua->xyC[i], 1, "N", xb->xC, 1, "N", i);
      xb->xC[i] = ae_c_sub(xb->xC[i], v);
   }
   xb->xC[n - 1] = ae_c_div(xb->xC[n - 1], lua->xyC[n - 1][n - 1]);
   for (i = n - 2; i >= 0; i--) {
      v = ae_v_cdotproduct(&lua->xyC[i][i + 1], 1, "N", &xb->xC[i + 1], 1, "N", n - i - 1);
      xb->xC[i] = ae_c_div(ae_c_sub(xb->xC[i], v), lua->xyC[i][i]);
   }
}

// Internal LU solver
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_cmatrixlusolveinternal(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *a, bool havea, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t rfs;
   ae_int_t nrfs;
   complex v;
   double verr;
   bool smallerr;
   bool terminatenexttime;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewVector(xc, 0, DT_COMPLEX);
   NewVector(y, 0, DT_COMPLEX);
   NewVector(bc, 0, DT_COMPLEX);
   NewVector(xa, 0, DT_COMPLEX);
   NewVector(xb, 0, DT_COMPLEX);
   NewVector(tx, 0, DT_COMPLEX);
   NewVector(tmpbuf, 0, DT_REAL);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   for (i = 0; i < n; i++) {
      if (p->xZ[i] > n - 1 || p->xZ[i] < i) {
         *info = -1;
         ae_frame_leave();
         return;
      }
   }
   ae_matrix_set_length(x, n, m);
   ae_vector_set_length(&y, n);
   ae_vector_set_length(&xc, n);
   ae_vector_set_length(&bc, n);
   ae_vector_set_length(&tx, n);
   ae_vector_set_length(&xa, n + 1);
   ae_vector_set_length(&xb, n + 1);
   ae_vector_set_length(&tmpbuf, 2 * n + 2);
// estimate condition number, test for near singularity
   rep->r1 = cmatrixlurcond1(lua, n);
   rep->rinf = cmatrixlurcondinf(lua, n);
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            x->xyC[i][j] = ae_complex_from_i(0);
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
   *info = 1;
// First phase: solve with TRSM()
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
         x->xyC[i][j] = b->xyC[i][j];
      }
   }
   for (i = 0; i < n; i++) {
      if (p->xZ[i] != i) {
         for (j = 0; j < m; j++) {
            swapc(&x->xyC[i][j], &x->xyC[p->xZ[i]][j]);
         }
      }
   }
   cmatrixlefttrsm(n, m, lua, 0, 0, false, true, 0, x, 0, 0);
   cmatrixlefttrsm(n, m, lua, 0, 0, true, false, 0, x, 0, 0);
// solve
   for (k = 0; k < m; k++) {
      ae_v_cmove(bc.xC, 1, &b->xyC[0][k], b->stride, "N", n);
      ae_v_cmove(xc.xC, 1, &x->xyC[0][k], x->stride, "N", n);
   // Iterative refinement of xc:
   // * calculate r = bc-A*xc using extra-precise dot product
   // * solve A*y = r
   // * update x:=x+r
   //
   // This cycle is executed until one of two things happens:
   // 1. maximum number of iterations reached
   // 2. last iteration decreased error to the lower limit
      if (havea) {
         nrfs = directdensesolvers_densesolverrfsmax(n, rep->r1, rep->rinf);
         terminatenexttime = false;
         for (rfs = 0; rfs < nrfs; rfs++) {
            if (terminatenexttime) {
               break;
            }
         // generate right part
            smallerr = true;
            ae_v_cmove(xb.xC, 1, xc.xC, 1, "N", n);
            for (i = 0; i < n; i++) {
               ae_v_cmove(xa.xC, 1, a->xyC[i], 1, "N", n);
               xa.xC[n] = ae_complex_from_i(-1);
               xb.xC[n] = bc.xC[i];
               xcdot(&xa, &xb, n + 1, &tmpbuf, &v, &verr);
               y.xC[i] = ae_c_neg(v);
               smallerr = smallerr && SmallC(v, 4.0 * verr);
            }
            if (smallerr) {
               terminatenexttime = true;
            }
         // solve and update
            directdensesolvers_cbasiclusolve(lua, p, n, &y);
            ae_v_cadd(xc.xC, 1, y.xC, 1, "N", n);
         }
      }
   // Store xc.
   // Post-scale result.
      ae_v_cmove(&x->xyC[0][k], x->stride, xc.xC, 1, "N", n);
   }
   ae_frame_leave();
}

// Complex dense solver for A*X=B with N*N  complex  matrix  A,  N*M  complex
// matrices  X  and  B.  "Slow-but-feature-rich"   version   which   provides
// additional functions, at the cost of slower  performance.  Faster  version
// may be invoked with CMatrixSolveMFast() function.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(N^3+M*N^2) complexity
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear  system
//            ! and  performs  iterative   refinement,   which   results   in
//            ! significant performance penalty  when  compared  with  "fast"
//            ! version  which  just  performs  LU  decomposition  and  calls
//            ! triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, CMatrixSolveMFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//     RFS     -   iterative refinement switch:
//                 * True - refinement is used.
//                   Less performance, more precision.
//                 * False - refinement is not used.
//                   More performance, less precision.
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void cmatrixsolvem(CMatrix *a, ae_int_t n, CMatrix *b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(da, 0, 0, DT_COMPLEX);
   NewMatrix(emptya, 0, 0, DT_COMPLEX);
   NewVector(p, 0, DT_INT);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&da, n, n);
// factorize, solve
   for (i = 0; i < n; i++) {
      ae_v_cmove(da.xyC[i], 1, a->xyC[i], 1, "N", n);
   }
   cmatrixlu(&da, n, n, &p);
   if (rfs) {
      directdensesolvers_cmatrixlusolveinternal(&da, &p, n, a, true, b, m, info, rep, x);
   } else {
      directdensesolvers_cmatrixlusolveinternal(&da, &p, n, &emptya, false, b, m, info, rep, x);
   }
   ae_frame_leave();
}

// Complex dense solver for A*X=B with N*N  complex  matrix  A,  N*M  complex
// matrices  X  and  B.  "Fast-but-lightweight" version which  provides  just
// triangular solver - and no additional functions like iterative  refinement
// or condition number estimation.
//
// Algorithm features:
// * O(N^3+M*N^2) complexity
// * no additional time consuming functions
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N,M]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 16.03.2015 by Sergey Bochkanov
// API: void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void cmatrixsolvemfast(CMatrix *a, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   NewVector(p, 0, DT_INT);
// Check for exact degeneracy
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   cmatrixlu(a, n, n, &p);
   for (i = 0; i < n; i++) {
      if (ae_c_eq_d(a->xyC[i][i], 0.0)) {
         for (j = 0; j < n; j++) {
            for (k = 0; k < m; k++) {
               b->xyC[j][k] = ae_complex_from_d(0.0);
            }
         }
         *info = -3;
         ae_frame_leave();
         return;
      }
   }
// Solve with TRSM()
   for (i = 0; i < n; i++) {
      if (p.xZ[i] != i) {
         for (j = 0; j < m; j++) {
            swapc(&b->xyC[i][j], &b->xyC[p.xZ[i]][j]);
         }
      }
   }
   cmatrixlefttrsm(n, m, a, 0, 0, false, true, 0, b, 0, 0);
   cmatrixlefttrsm(n, m, a, 0, 0, true, false, 0, b, 0, 0);
   *info = 1;
   ae_frame_leave();
}

// Complex dense solver for A*x=B with N*N complex matrix A and  N*1  complex
// vectors x and b. "Slow-but-feature-rich" version of the solver.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(N^3) complexity
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear  system
//            ! and  performs  iterative   refinement,   which   results   in
//            ! significant performance penalty  when  compared  with  "fast"
//            ! version  which  just  performs  LU  decomposition  and  calls
//            ! triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, CMatrixSolveFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void cmatrixsolve(CMatrix *a, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_COMPLEX);
   NewMatrix(xm, 0, 0, DT_COMPLEX);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_cmove(bm.xyC[0], bm.stride, b->xC, 1, "N", n);
   cmatrixsolvem(a, n, &bm, 1, true, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_cmove(x->xC, 1, xm.xyC[0], xm.stride, "N", n);
   ae_frame_leave();
}

// Complex dense solver for A*x=B with N*N complex matrix A and  N*1  complex
// vectors x and b. "Fast-but-lightweight" version of the solver.
//
// Algorithm features:
// * O(N^3) complexity
// * no additional time consuming features, just triangular solver
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
void cmatrixsolvefast(CMatrix *a, ae_int_t n, CVector *b, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   NewVector(p, 0, DT_INT);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   cmatrixlu(a, n, n, &p);
   for (i = 0; i < n; i++) {
      if (ae_c_eq_d(a->xyC[i][i], 0.0)) {
         for (j = 0; j < n; j++) {
            b->xC[j] = ae_complex_from_d(0.0);
         }
         *info = -3;
         ae_frame_leave();
         return;
      }
   }
   directdensesolvers_cbasiclusolve(a, &p, n, b);
   *info = 1;
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N complex A given by its  LU  decomposition,
// and N*M matrices X and B (multiple right sides).   "Slow-but-feature-rich"
// version of the solver.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(M*N^2) complexity
// * condition number estimation
//
// No iterative refinement  is provided because exact form of original matrix
// is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant  performance   penalty   when
//            ! compared with "fast"  version  which  just  calls  triangular
//            ! solver.
//            !
//            ! This performance penalty is especially apparent when you  use
//            ! ALGLIB parallel capabilities (condition number estimation  is
//            ! inherently  sequential).  It   also   becomes significant for
//            ! small-scale problems.
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! CMatrixLUSolveMFast() function.
//
// Inputs:
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
//     P       -   array[0..N-1], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void cmatrixlusolvem(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(emptya, 0, 0, DT_COMPLEX);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
// solve
   directdensesolvers_cmatrixlusolveinternal(lua, p, n, &emptya, false, b, m, info, rep, x);
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N complex A given by its  LU  decomposition,
// and N*M matrices X and B (multiple  right  sides).  "Fast-but-lightweight"
// version of the solver.
//
// Algorithm features:
// * O(M*N^2) complexity
// * no additional time-consuming features
//
// Inputs:
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
//     P       -   array[0..N-1], pivots array, RMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N,M]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void cmatrixlusolvemfast(CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   *info = 0;
// Check for exact degeneracy
   if (n <= 0 || m <= 0) {
      *info = -1;
      return;
   }
   for (i = 0; i < n; i++) {
      if (ae_c_eq_d(lua->xyC[i][i], 0.0)) {
         for (j = 0; j < n; j++) {
            for (k = 0; k < m; k++) {
               b->xyC[j][k] = ae_complex_from_d(0.0);
            }
         }
         *info = -3;
         return;
      }
   }
// Solve with TRSM()
   for (i = 0; i < n; i++) {
      if (p->xZ[i] != i) {
         for (j = 0; j < m; j++) {
            swapc(&b->xyC[i][j], &b->xyC[p->xZ[i]][j]);
         }
      }
   }
   cmatrixlefttrsm(n, m, lua, 0, 0, false, true, 0, b, 0, 0);
   cmatrixlefttrsm(n, m, lua, 0, 0, true, false, 0, b, 0, 0);
   *info = 1;
}

// Complex dense linear solver for A*x=b with complex N*N A  given  by its LU
// decomposition and N*1 vectors x and b. This is  "slow-but-robust"  version
// of  the  complex  linear  solver  with  additional  features   which   add
// significant performance overhead. Faster version  is  CMatrixLUSolveFast()
// function.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(N^2) complexity
// * condition number estimation
//
// No iterative refinement is provided because exact form of original matrix
// is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which results in 10-15x  performance  penalty  when  compared
//            ! with "fast" version which just calls triangular solver.
//            !
//            ! This performance penalty is insignificant  when compared with
//            ! cost of large LU decomposition.  However,  if you  call  this
//            ! function many times for the same  left  side,  this  overhead
//            ! BECOMES significant. It  also  becomes significant for small-
//            ! scale problems.
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! CMatrixLUSolveFast() function.
//
// Inputs:
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
//     P       -   array[0..N-1], pivots array, CMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void cmatrixlusolve(CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_COMPLEX);
   NewMatrix(xm, 0, 0, DT_COMPLEX);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_cmove(bm.xyC[0], bm.stride, b->xC, 1, "N", n);
   cmatrixlusolvem(lua, p, n, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_cmove(x->xC, 1, xm.xyC[0], xm.stride, "N", n);
   ae_frame_leave();
}

// Complex dense linear solver for A*x=b with N*N complex A given by  its  LU
// decomposition and N*1 vectors x and b. This is  fast  lightweight  version
// of solver, which is significantly faster than CMatrixLUSolve(),  but  does
// not provide additional information (like condition numbers).
//
// Algorithm features:
// * O(N^2) complexity
// * no additional time-consuming features, just triangular solver
//
// Inputs:
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
//     P       -   array[0..N-1], pivots array, CMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is exactly singular (ill conditioned matrices
//                         are not recognized).
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N]:
//                 * info > 0  => overwritten by solution
//                 * info = -3 => filled by zeros
//
// NOTE: unlike  CMatrixLUSolve(),  this   function   does   NOT   check  for
//       near-degeneracy of input matrix. It  checks  for  EXACT  degeneracy,
//       because this check is easy to do. However,  very  badly  conditioned
//       matrices may went unnoticed.
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
void cmatrixlusolvefast(CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info) {
   ae_int_t i;
   ae_int_t j;
   *info = 0;
   if (n <= 0) {
      *info = -1;
      return;
   }
   for (i = 0; i < n; i++) {
      if (ae_c_eq_d(lua->xyC[i][i], 0.0)) {
         for (j = 0; j < n; j++) {
            b->xC[j] = ae_complex_from_d(0.0);
         }
         *info = -3;
         return;
      }
   }
   directdensesolvers_cbasiclusolve(lua, p, n, b);
   *info = 1;
}

// Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(M*N^2) complexity
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
//     P       -   array[0..N-1], pivots array, CMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void cmatrixmixedsolvem(CMatrix *a, CMatrix *lua, ZVector *p, ae_int_t n, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      return;
   }
// solve
   directdensesolvers_cmatrixlusolveinternal(lua, p, n, a, true, b, m, info, rep, x);
}

// Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * iterative refinement
// * O(N^2) complexity
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
//     P       -   array[0..N-1], pivots array, CMatrixLU result
//     N       -   size of A
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or exactly singular.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void cmatrixmixedsolve(CMatrix *a, CMatrix *lua, ZVector *p, ae_int_t n, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_COMPLEX);
   NewMatrix(xm, 0, 0, DT_COMPLEX);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_cmove(bm.xyC[0], bm.stride, b->xC, 1, "N", n);
   cmatrixmixedsolvem(a, lua, p, n, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_cmove(x->xC, 1, xm.xyC[0], xm.stride, "N", n);
   ae_frame_leave();
}

// Internal Cholesky solver
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_spdmatrixcholeskysolveinternal(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *a, bool havea, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   ae_int_t i;
   ae_int_t j;
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      return;
   }
   ae_matrix_set_length(x, n, m);
// estimate condition number, test for near singularity
   rep->r1 = spdmatrixcholeskyrcond(cha, n, isupper);
   rep->rinf = rep->r1;
   if (rep->r1 < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            x->xyR[i][j] = 0.0;
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      return;
   }
   *info = 1;
// Solve with TRSM()
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
         x->xyR[i][j] = b->xyR[i][j];
      }
   }
   if (isupper) {
      rmatrixlefttrsm(n, m, cha, 0, 0, true, false, 1, x, 0, 0);
      rmatrixlefttrsm(n, m, cha, 0, 0, true, false, 0, x, 0, 0);
   } else {
      rmatrixlefttrsm(n, m, cha, 0, 0, false, false, 0, x, 0, 0);
      rmatrixlefttrsm(n, m, cha, 0, 0, false, false, 1, x, 0, 0);
   }
}

// Dense solver for A*X=B with N*N symmetric positive definite matrix A,  and
// N*M vectors X and B. It is "slow-but-feature-rich" version of the solver.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * O(N^3+M*N^2) complexity
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant   performance   penalty  when
//            ! compared with "fast" version  which  just  performs  Cholesky
//            ! decomposition and calls triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, SPDMatrixSolveMFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or non-SPD.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void spdmatrixsolvem(RMatrix *a, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(da, 0, 0, DT_REAL);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&da, n, n);
// factorize
// solve
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i;
      }
      ae_v_move(&da.xyR[i][j1], 1, &a->xyR[i][j1], 1, j2 - j1 + 1);
   }
   if (!spdmatrixcholesky(&da, n, isupper)) {
      ae_matrix_set_length(x, n, m);
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            x->xyR[i][j] = 0.0;
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
   *info = 1;
   directdensesolvers_spdmatrixcholeskysolveinternal(&da, n, isupper, a, true, b, m, info, rep, x);
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N symmetric positive definite matrix A,  and
// N*M vectors X and B. It is "fast-but-lightweight" version of the solver.
//
// Algorithm features:
// * O(N^3+M*N^2) complexity
// * matrix is represented by its upper or lower triangle
// * no additional time consuming features
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     B       -   array[N,M], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 17.03.2015 by Sergey Bochkanov
// API: void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void spdmatrixsolvemfast(RMatrix *a, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   if (!spdmatrixcholesky(a, n, isupper)) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            b->xyR[i][j] = 0.0;
         }
      }
      *info = -3;
      ae_frame_leave();
      return;
   }
   if (isupper) {
      rmatrixlefttrsm(n, m, a, 0, 0, true, false, 1, b, 0, 0);
      rmatrixlefttrsm(n, m, a, 0, 0, true, false, 0, b, 0, 0);
   } else {
      rmatrixlefttrsm(n, m, a, 0, 0, false, false, 0, b, 0, 0);
      rmatrixlefttrsm(n, m, a, 0, 0, false, false, 1, b, 0, 0);
   }
   ae_frame_leave();
}

// Dense linear solver for A*x=b with N*N real  symmetric  positive  definite
// matrix A,  N*1 vectors x and b.  "Slow-but-feature-rich"  version  of  the
// solver.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * O(N^3) complexity
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant   performance   penalty  when
//            ! compared with "fast" version  which  just  performs  Cholesky
//            ! decomposition and calls triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, SPDMatrixSolveFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    matrix is very badly conditioned or non-SPD.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void spdmatrixsolve(RMatrix *a, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_REAL);
   NewMatrix(xm, 0, 0, DT_REAL);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_move(bm.xyR[0], bm.stride, b->xR, 1, n);
   spdmatrixsolvem(a, n, isupper, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_move(x->xR, 1, xm.xyR[0], xm.stride, n);
   ae_frame_leave();
}

// Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.
//
// This subroutine assumes that:
// * A*ScaleA is well scaled
// * A is well-conditioned, so no zero divisions or overflow may occur
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_spdbasiccholeskysolve(RMatrix *cha, ae_int_t n, bool isupper, RVector *xb) {
   ae_int_t i;
   double v;
// A = L*L' or A=U'*U
   if (isupper) {
   // Solve U'*y=b first.
      for (i = 0; i < n; i++) {
         xb->xR[i] /= cha->xyR[i][i];
         if (i < n - 1) {
            v = xb->xR[i];
            ae_v_subd(&xb->xR[i + 1], 1, &cha->xyR[i][i + 1], 1, n - i - 1, v);
         }
      }
   // Solve U*x=y then.
      for (i = n - 1; i >= 0; i--) {
         if (i < n - 1) {
            v = ae_v_dotproduct(&cha->xyR[i][i + 1], 1, &xb->xR[i + 1], 1, n - i - 1);
            xb->xR[i] -= v;
         }
         xb->xR[i] /= cha->xyR[i][i];
      }
   } else {
   // Solve L*y=b first
      for (i = 0; i < n; i++) {
         if (i > 0) {
            v = ae_v_dotproduct(cha->xyR[i], 1, xb->xR, 1, i);
            xb->xR[i] -= v;
         }
         xb->xR[i] /= cha->xyR[i][i];
      }
   // Solve L'*x=y then.
      for (i = n - 1; i >= 0; i--) {
         xb->xR[i] /= cha->xyR[i][i];
         if (i > 0) {
            v = xb->xR[i];
            ae_v_subd(xb->xR, 1, cha->xyR[i], 1, i, v);
         }
      }
   }
}

// Dense linear solver for A*x=b with N*N real  symmetric  positive  definite
// matrix A,  N*1 vectors x and  b.  "Fast-but-lightweight"  version  of  the
// solver.
//
// Algorithm features:
// * O(N^3) complexity
// * matrix is represented by its upper or lower triangle
// * no additional time consuming features like condition number estimation
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or non-SPD
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     B       -   array[N], it contains:
//                 * info > 0  => solution
//                 * info = -3 => filled by zeros
// ALGLIB: Copyright 17.03.2015 by Sergey Bochkanov
// API: void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
void spdmatrixsolvefast(RMatrix *a, ae_int_t n, bool isupper, RVector *b, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   if (!spdmatrixcholesky(a, n, isupper)) {
      for (i = 0; i < n; i++) {
         b->xR[i] = 0.0;
      }
      *info = -3;
      ae_frame_leave();
      return;
   }
   directdensesolvers_spdbasiccholeskysolve(a, n, isupper, b);
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N symmetric positive definite matrix A given
// by its Cholesky decomposition, and N*M vectors X and B. It  is  "slow-but-
// feature-rich" version of the solver which estimates  condition  number  of
// the system.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(M*N^2) complexity
// * condition number estimation
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant  performance   penalty   when
//            ! compared with "fast"  version  which  just  calls  triangular
//            ! solver. Amount of  overhead  introduced  depends  on  M  (the
//            ! larger - the more efficient).
//            !
//            ! This performance penalty is insignificant  when compared with
//            ! cost of large LU decomposition.  However,  if you  call  this
//            ! function many times for the same  left  side,  this  overhead
//            ! BECOMES significant. It  also  becomes significant for small-
//            ! scale problems (N < 50).
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! SPDMatrixCholeskySolveMFast() function.
//
// Inputs:
//     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
//                 SPDMatrixCholesky result
//     N       -   size of CHA
//     IsUpper -   what half of CHA is provided
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or badly conditioned
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N]:
//                 * for info > 0 contains solution
//                 * for info = -3 filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void spdmatrixcholeskysolvem(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(emptya, 0, 0, DT_REAL);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
// solve
   directdensesolvers_spdmatrixcholeskysolveinternal(cha, n, isupper, &emptya, false, b, m, info, rep, x);
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N symmetric positive definite matrix A given
// by its Cholesky decomposition, and N*M vectors X and B. It  is  "fast-but-
// lightweight" version of  the  solver  which  just  solves  linear  system,
// without any additional functions.
//
// Algorithm features:
// * O(M*N^2) complexity
// * matrix is represented by its upper or lower triangle
// * no additional functionality
//
// Inputs:
//     CHA     -   array[N,N], Cholesky decomposition,
//                 SPDMatrixCholesky result
//     N       -   size of CHA
//     IsUpper -   what half of CHA is provided
//     B       -   array[N,M], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or badly conditioned
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     B       -   array[N]:
//                 * for info > 0 overwritten by solution
//                 * for info = -3 filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
// API: void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void spdmatrixcholeskysolvemfast(RMatrix *cha, ae_int_t n, bool isupper, RMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      return;
   }
   for (k = 0; k < n; k++) {
      if (cha->xyR[k][k] == 0.0) {
         for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
               b->xyR[i][j] = 0.0;
            }
         }
         *info = -3;
         return;
      }
   }
   if (isupper) {
      rmatrixlefttrsm(n, m, cha, 0, 0, true, false, 1, b, 0, 0);
      rmatrixlefttrsm(n, m, cha, 0, 0, true, false, 0, b, 0, 0);
   } else {
      rmatrixlefttrsm(n, m, cha, 0, 0, false, false, 0, b, 0, 0);
      rmatrixlefttrsm(n, m, cha, 0, 0, false, false, 1, b, 0, 0);
   }
}

// Dense solver for A*x=b with N*N symmetric positive definite matrix A given
// by its Cholesky decomposition, and N*1 real vectors x and b. This is "slow-
// but-feature-rich"  version  of  the  solver  which,  in  addition  to  the
// solution, performs condition number estimation.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(N^2) complexity
// * condition number estimation
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which results in 10-15x  performance  penalty  when  compared
//            ! with "fast" version which just calls triangular solver.
//            !
//            ! This performance penalty is insignificant  when compared with
//            ! cost of large LU decomposition.  However,  if you  call  this
//            ! function many times for the same  left  side,  this  overhead
//            ! BECOMES significant. It  also  becomes significant for small-
//            ! scale problems (N < 50).
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! SPDMatrixCholeskySolveFast() function.
//
// Inputs:
//     CHA     -   array[N,N], Cholesky decomposition,
//                 SPDMatrixCholesky result
//     N       -   size of A
//     IsUpper -   what half of CHA is provided
//     B       -   array[N], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or ill conditioned
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N]:
//                 * for info > 0  - solution
//                 * for info = -3 - filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void spdmatrixcholeskysolve(RMatrix *cha, ae_int_t n, bool isupper, RVector *b, ae_int_t *info, densesolverreport *rep, RVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_REAL);
   NewMatrix(xm, 0, 0, DT_REAL);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_move(bm.xyR[0], bm.stride, b->xR, 1, n);
   spdmatrixcholeskysolvem(cha, n, isupper, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_move(x->xR, 1, xm.xyR[0], xm.stride, n);
   ae_frame_leave();
}

// Dense solver for A*x=b with N*N symmetric positive definite matrix A given
// by its Cholesky decomposition, and N*1 real vectors x and b. This is "fast-
// but-lightweight" version of the solver.
//
// Algorithm features:
// * O(N^2) complexity
// * matrix is represented by its upper or lower triangle
// * no additional features
//
// Inputs:
//     CHA     -   array[N,N], Cholesky decomposition,
//                 SPDMatrixCholesky result
//     N       -   size of A
//     IsUpper -   what half of CHA is provided
//     B       -   array[N], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or ill conditioned
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N]:
//                 * for info > 0  - overwritten by solution
//                 * for info = -3 - filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
void spdmatrixcholeskysolvefast(RMatrix *cha, ae_int_t n, bool isupper, RVector *b, ae_int_t *info) {
   ae_int_t i;
   ae_int_t k;
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      return;
   }
   for (k = 0; k < n; k++) {
      if (cha->xyR[k][k] == 0.0) {
         for (i = 0; i < n; i++) {
            b->xR[i] = 0.0;
         }
         *info = -3;
         return;
      }
   }
   directdensesolvers_spdbasiccholeskysolve(cha, n, isupper, b);
}

// Internal Cholesky solver
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_hpdmatrixcholeskysolveinternal(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *a, bool havea, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewVector(xc, 0, DT_COMPLEX);
   NewVector(y, 0, DT_COMPLEX);
   NewVector(bc, 0, DT_COMPLEX);
   NewVector(xa, 0, DT_COMPLEX);
   NewVector(xb, 0, DT_COMPLEX);
   NewVector(tx, 0, DT_COMPLEX);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(x, n, m);
   ae_vector_set_length(&y, n);
   ae_vector_set_length(&xc, n);
   ae_vector_set_length(&bc, n);
   ae_vector_set_length(&tx, n + 1);
   ae_vector_set_length(&xa, n + 1);
   ae_vector_set_length(&xb, n + 1);
// estimate condition number, test for near singularity
   rep->r1 = hpdmatrixcholeskyrcond(cha, n, isupper);
   rep->rinf = rep->r1;
   if (rep->r1 < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            x->xyC[i][j] = ae_complex_from_i(0);
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
   *info = 1;
// solve
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
         x->xyC[i][j] = b->xyC[i][j];
      }
   }
   if (isupper) {
      cmatrixlefttrsm(n, m, cha, 0, 0, true, false, 2, x, 0, 0);
      cmatrixlefttrsm(n, m, cha, 0, 0, true, false, 0, x, 0, 0);
   } else {
      cmatrixlefttrsm(n, m, cha, 0, 0, false, false, 0, x, 0, 0);
      cmatrixlefttrsm(n, m, cha, 0, 0, false, false, 2, x, 0, 0);
   }
   ae_frame_leave();
}

// Dense solver for A*X=B, with N*N Hermitian positive definite matrix A  and
// N*M  complex  matrices  X  and  B.  "Slow-but-feature-rich" version of the
// solver.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * O(N^3+M*N^2) complexity
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant  performance   penalty   when
//            ! compared with "fast"  version  which  just  calls  triangular
//            ! solver.
//            !
//            ! This performance penalty is especially apparent when you  use
//            ! ALGLIB parallel capabilities (condition number estimation  is
//            ! inherently  sequential).  It   also   becomes significant for
//            ! small-scale problems (N < 100).
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! HPDMatrixSolveMFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   same as in RMatrixSolve.
//                 Returns -3 for non-HPD matrices.
//     Rep     -   same as in RMatrixSolve
//     X       -   same as in RMatrixSolve
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void hpdmatrixsolvem(CMatrix *a, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(da, 0, 0, DT_COMPLEX);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&da, n, n);
// factorize matrix, solve
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i;
      }
      ae_v_cmove(&da.xyC[i][j1], 1, &a->xyC[i][j1], 1, "N", j2 - j1 + 1);
   }
   if (!hpdmatrixcholesky(&da, n, isupper)) {
      ae_matrix_set_length(x, n, m);
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            x->xyC[i][j] = ae_complex_from_i(0);
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
   *info = 1;
   directdensesolvers_hpdmatrixcholeskysolveinternal(&da, n, isupper, a, true, b, m, info, rep, x);
   ae_frame_leave();
}

// Dense solver for A*X=B, with N*N Hermitian positive definite matrix A  and
// N*M complex matrices X and B. "Fast-but-lightweight" version of the solver.
//
// Algorithm features:
// * O(N^3+M*N^2) complexity
// * matrix is represented by its upper or lower triangle
// * no additional time consuming features like condition number estimation
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1,0..M-1], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly  singular or is not positive definite.
//                         B is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[0..N-1]:
//                 * overwritten by solution
//                 * zeros, if problem was not solved
// ALGLIB: Copyright 17.03.2015 by Sergey Bochkanov
// API: void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void hpdmatrixsolvemfast(CMatrix *a, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   if (!hpdmatrixcholesky(a, n, isupper)) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            b->xyC[i][j] = ae_complex_from_d(0.0);
         }
      }
      *info = -3;
      ae_frame_leave();
      return;
   }
   if (isupper) {
      cmatrixlefttrsm(n, m, a, 0, 0, true, false, 2, b, 0, 0);
      cmatrixlefttrsm(n, m, a, 0, 0, true, false, 0, b, 0, 0);
   } else {
      cmatrixlefttrsm(n, m, a, 0, 0, false, false, 0, b, 0, 0);
      cmatrixlefttrsm(n, m, a, 0, 0, false, false, 2, b, 0, 0);
   }
   ae_frame_leave();
}

// Dense solver for A*x=b, with N*N Hermitian positive definite matrix A, and
// N*1 complex vectors  x  and  b.  "Slow-but-feature-rich"  version  of  the
// solver.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * condition number estimation
// * O(N^3) complexity
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant   performance   penalty  when
//            ! compared with "fast" version  which  just  performs  Cholesky
//            ! decomposition and calls triangular solver.
//            !
//            ! This  performance  penalty  is  especially  visible  in   the
//            ! multithreaded mode, because both condition number  estimation
//            ! and   iterative    refinement   are   inherently   sequential
//            ! calculations.
//            !
//            ! Thus, if you need high performance and if you are pretty sure
//            ! that your system is well conditioned, we  strongly  recommend
//            ! you to use faster solver, HPDMatrixSolveFast() function.
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   same as in RMatrixSolve
//                 Returns -3 for non-HPD matrices.
//     Rep     -   same as in RMatrixSolve
//     X       -   same as in RMatrixSolve
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void hpdmatrixsolve(CMatrix *a, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_COMPLEX);
   NewMatrix(xm, 0, 0, DT_COMPLEX);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_cmove(bm.xyC[0], bm.stride, b->xC, 1, "N", n);
   hpdmatrixsolvem(a, n, isupper, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_cmove(x->xC, 1, xm.xyC[0], xm.stride, "N", n);
   ae_frame_leave();
}

// Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.
//
// This subroutine assumes that:
// * A*ScaleA is well scaled
// * A is well-conditioned, so no zero divisions or overflow may occur
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static void directdensesolvers_hpdbasiccholeskysolve(CMatrix *cha, ae_int_t n, bool isupper, CVector *xb) {
   ae_int_t i;
   complex v;
// A = L*L' or A=U'*U
   if (isupper) {
   // Solve U'*y=b first.
      for (i = 0; i < n; i++) {
         xb->xC[i] = ae_c_div(xb->xC[i], conj(cha->xyC[i][i]));
         if (i < n - 1) {
            v = xb->xC[i];
            ae_v_csubc(&xb->xC[i + 1], 1, &cha->xyC[i][i + 1], 1, "Conj", n - i - 1, v);
         }
      }
   // Solve U*x=y then.
      for (i = n - 1; i >= 0; i--) {
         if (i < n - 1) {
            v = ae_v_cdotproduct(&cha->xyC[i][i + 1], 1, "N", &xb->xC[i + 1], 1, "N", n - i - 1);
            xb->xC[i] = ae_c_sub(xb->xC[i], v);
         }
         xb->xC[i] = ae_c_div(xb->xC[i], cha->xyC[i][i]);
      }
   } else {
   // Solve L*y=b first
      for (i = 0; i < n; i++) {
         if (i > 0) {
            v = ae_v_cdotproduct(cha->xyC[i], 1, "N", xb->xC, 1, "N", i);
            xb->xC[i] = ae_c_sub(xb->xC[i], v);
         }
         xb->xC[i] = ae_c_div(xb->xC[i], cha->xyC[i][i]);
      }
   // Solve L'*x=y then.
      for (i = n - 1; i >= 0; i--) {
         xb->xC[i] = ae_c_div(xb->xC[i], conj(cha->xyC[i][i]));
         if (i > 0) {
            v = xb->xC[i];
            ae_v_csubc(xb->xC, 1, cha->xyC[i], 1, "Conj", i, v);
         }
      }
   }
}

// Dense solver for A*x=b, with N*N Hermitian positive definite matrix A, and
// N*1 complex vectors  x  and  b.  "Fast-but-lightweight"  version  of   the
// solver without additional functions.
//
// Algorithm features:
// * O(N^3) complexity
// * matrix is represented by its upper or lower triangle
// * no additional time consuming functions
//
// Inputs:
//     A       -   array[0..N-1,0..N-1], system matrix
//     N       -   size of A
//     IsUpper -   what half of A is provided
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or not positive definite
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     B       -   array[0..N-1]:
//                 * overwritten by solution
//                 * zeros, if A is exactly singular (diagonal of its LU
//                   decomposition has exact zeros).
// ALGLIB: Copyright 17.03.2015 by Sergey Bochkanov
// API: void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
void hpdmatrixsolvefast(CMatrix *a, ae_int_t n, bool isupper, CVector *b, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   if (!hpdmatrixcholesky(a, n, isupper)) {
      for (i = 0; i < n; i++) {
         b->xC[i] = ae_complex_from_d(0.0);
      }
      *info = -3;
      ae_frame_leave();
      return;
   }
   directdensesolvers_hpdbasiccholeskysolve(a, n, isupper, b);
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N Hermitian positive definite matrix A given
// by its Cholesky decomposition and N*M complex matrices X  and  B.  This is
// "slow-but-feature-rich" version of the solver which, in  addition  to  the
// solution, estimates condition number of the system.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(M*N^2) complexity
// * condition number estimation
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which  results  in  significant  performance   penalty   when
//            ! compared with "fast"  version  which  just  calls  triangular
//            ! solver. Amount of  overhead  introduced  depends  on  M  (the
//            ! larger - the more efficient).
//            !
//            ! This performance penalty is insignificant  when compared with
//            ! cost of large Cholesky decomposition.  However,  if  you call
//            ! this  function  many  times  for  the same  left  side,  this
//            ! overhead BECOMES significant. It  also   becomes  significant
//            ! for small-scale problems (N < 50).
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! HPDMatrixCholeskySolveMFast() function.
//
// Inputs:
//     CHA     -   array[N,N], Cholesky decomposition,
//                 HPDMatrixCholesky result
//     N       -   size of CHA
//     IsUpper -   what half of CHA is provided
//     B       -   array[N,M], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is singular, or VERY close to singular.
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N]:
//                 * for info > 0 contains solution
//                 * for info = -3 filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void hpdmatrixcholeskysolvem(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetMatrix(x);
   NewMatrix(emptya, 0, 0, DT_COMPLEX);
// prepare: check inputs, allocate space...
   if (n <= 0 || m <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
// 1. scale matrix, max(|U[i,j]|)
// 2. factorize scaled matrix
// 3. solve
   directdensesolvers_hpdmatrixcholeskysolveinternal(cha, n, isupper, &emptya, false, b, m, info, rep, x);
   ae_frame_leave();
}

// Dense solver for A*X=B with N*N Hermitian positive definite matrix A given
// by its Cholesky decomposition and N*M complex matrices X  and  B.  This is
// "fast-but-lightweight" version of the solver.
//
// Algorithm features:
// * O(M*N^2) complexity
// * matrix is represented by its upper or lower triangle
// * no additional time-consuming features
//
// Inputs:
//     CHA     -   array[N,N], Cholesky decomposition,
//                 HPDMatrixCholesky result
//     N       -   size of CHA
//     IsUpper -   what half of CHA is provided
//     B       -   array[N,M], right part
//     M       -   right part size
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is singular, or VERY close to singular.
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task was solved
//     B       -   array[N]:
//                 * for info > 0 overwritten by solution
//                 * for info = -3 filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
// API: void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void hpdmatrixcholeskysolvemfast(CMatrix *cha, ae_int_t n, bool isupper, CMatrix *b, ae_int_t m, ae_int_t *info) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      return;
   }
   for (k = 0; k < n; k++) {
      if (cha->xyC[k][k].x == 0.0 && cha->xyC[k][k].y == 0.0) {
         for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
               b->xyC[i][j] = ae_complex_from_d(0.0);
            }
         }
         *info = -3;
         return;
      }
   }
   if (isupper) {
      cmatrixlefttrsm(n, m, cha, 0, 0, true, false, 2, b, 0, 0);
      cmatrixlefttrsm(n, m, cha, 0, 0, true, false, 0, b, 0, 0);
   } else {
      cmatrixlefttrsm(n, m, cha, 0, 0, false, false, 0, b, 0, 0);
      cmatrixlefttrsm(n, m, cha, 0, 0, false, false, 2, b, 0, 0);
   }
}

// Dense solver for A*x=b with N*N Hermitian positive definite matrix A given
// by its Cholesky decomposition, and N*1 complex vectors x and  b.  This  is
// "slow-but-feature-rich" version of the solver  which  estimates  condition
// number of the system.
//
// Algorithm features:
// * automatic detection of degenerate cases
// * O(N^2) complexity
// * condition number estimation
// * matrix is represented by its upper or lower triangle
//
// No iterative refinement is provided because such partial representation of
// matrix does not allow efficient calculation of extra-precise  matrix-vector
// products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
// need iterative refinement.
//
// IMPORTANT: ! this function is NOT the most efficient linear solver provided
//            ! by ALGLIB. It estimates condition  number  of  linear system,
//            ! which results in 10-15x  performance  penalty  when  compared
//            ! with "fast" version which just calls triangular solver.
//            !
//            ! This performance penalty is insignificant  when compared with
//            ! cost of large LU decomposition.  However,  if you  call  this
//            ! function many times for the same  left  side,  this  overhead
//            ! BECOMES significant. It  also  becomes significant for small-
//            ! scale problems (N < 50).
//            !
//            ! In such cases we strongly recommend you to use faster solver,
//            ! HPDMatrixCholeskySolveFast() function.
//
// Inputs:
//     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
//                 SPDMatrixCholesky result
//     N       -   size of A
//     IsUpper -   what half of CHA is provided
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or ill conditioned
//                         X is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     Rep     -   additional report, following fields are set:
//                 * rep.r1    condition number in 1-norm
//                 * rep.rinf  condition number in inf-norm
//     X       -   array[N]:
//                 * for info > 0  - solution
//                 * for info = -3 - filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
// API: void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void hpdmatrixcholeskysolve(CMatrix *cha, ae_int_t n, bool isupper, CVector *b, ae_int_t *info, densesolverreport *rep, CVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverreport, rep);
   SetVector(x);
   NewMatrix(bm, 0, 0, DT_COMPLEX);
   NewMatrix(xm, 0, 0, DT_COMPLEX);
   if (n <= 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(&bm, n, 1);
   ae_v_cmove(bm.xyC[0], bm.stride, b->xC, 1, "N", n);
   hpdmatrixcholeskysolvem(cha, n, isupper, &bm, 1, info, rep, &xm);
   ae_vector_set_length(x, n);
   ae_v_cmove(x->xC, 1, xm.xyC[0], xm.stride, "N", n);
   ae_frame_leave();
}

// Dense solver for A*x=b with N*N Hermitian positive definite matrix A given
// by its Cholesky decomposition, and N*1 complex vectors x and  b.  This  is
// "fast-but-lightweight" version of the solver.
//
// Algorithm features:
// * O(N^2) complexity
// * matrix is represented by its upper or lower triangle
// * no additional time-consuming features
//
// Inputs:
//     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
//                 SPDMatrixCholesky result
//     N       -   size of A
//     IsUpper -   what half of CHA is provided
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is is exactly singular or ill conditioned
//                         B is filled by zeros in such cases.
//                 * -1    N <= 0 was passed
//                 *  1    task is solved
//     B       -   array[N]:
//                 * for info > 0  - overwritten by solution
//                 * for info = -3 - filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
// API: void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
void hpdmatrixcholeskysolvefast(CMatrix *cha, ae_int_t n, bool isupper, CVector *b, ae_int_t *info) {
   ae_int_t i;
   ae_int_t k;
   *info = 0;
   *info = 1;
   if (n <= 0) {
      *info = -1;
      return;
   }
   for (k = 0; k < n; k++) {
      if (cha->xyC[k][k].x == 0.0 && cha->xyC[k][k].y == 0.0) {
         for (i = 0; i < n; i++) {
            b->xC[i] = ae_complex_from_d(0.0);
         }
         *info = -3;
         return;
      }
   }
   directdensesolvers_hpdbasiccholeskysolve(cha, n, isupper, b);
}

// Internal subroutine.
// Returns maximum count of RFS iterations as function of:
// 1. machine epsilon
// 2. task size.
// 3. norm-2 condition number
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
static ae_int_t directdensesolvers_densesolverrfsmaxv2(ae_int_t n, double r2) {
   ae_int_t result;
   result = directdensesolvers_densesolverrfsmax(n, 0.0, 0.0);
   return result;
}

// Dense solver.
//
// This subroutine finds solution of the linear system A*X=B with non-square,
// possibly degenerate A.  System  is  solved in the least squares sense, and
// general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
// returned. If A is non-degenerate, solution in the usual sense is returned.
//
// Algorithm features:
// * automatic detection (and correct handling!) of degenerate cases
// * iterative refinement
// * O(N^3) complexity
//
// Inputs:
//     A       -   array[0..NRows-1,0..NCols-1], system matrix
//     NRows   -   vertical size of A
//     NCols   -   horizontal size of A
//     B       -   array[0..NCols-1], right part
//     Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
//                 considered  zero.  Set  it to 0.0, if you don't understand
//                 what it means, so the solver will choose good value on its
//                 own.
//
// Outputs:
//     Info    -   return code:
//                 * -4    SVD subroutine failed
//                 * -1    if NRows <= 0 or NCols <= 0 or Threshold < 0 was passed
//                 *  1    if task is solved
//     Rep     -   solver report, see below for more info
//     X       -   array[0..N-1,0..M-1], it contains:
//                 * solution of A*X=B (even for singular A)
//                 * zeros, if SVD subroutine failed
//
// SOLVER REPORT
//
// Subroutine sets following fields of the Rep structure:
// * R2        reciprocal of condition number: 1/cond(A), 2-norm.
// * N         = NCols
// * K         dim(Null(A))
// * CX        array[0..N-1,0..K-1], kernel of A.
//             Columns of CX store such vectors that A*CX[i]=0.
// ALGLIB: Copyright 24.08.2009 by Sergey Bochkanov
// API: void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);
void rmatrixsolvels(RMatrix *a, ae_int_t nrows, ae_int_t ncols, RVector *b, double threshold, ae_int_t *info, densesolverlsreport *rep, RVector *x) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t nsv;
   ae_int_t kernelidx;
   double v;
   double verr;
   bool svdfailed;
   bool zeroa;
   ae_int_t rfs;
   ae_int_t nrfs;
   bool terminatenexttime;
   bool smallerr;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(densesolverlsreport, rep);
   SetVector(x);
   NewVector(sv, 0, DT_REAL);
   NewMatrix(u, 0, 0, DT_REAL);
   NewMatrix(vt, 0, 0, DT_REAL);
   NewVector(rp, 0, DT_REAL);
   NewVector(utb, 0, DT_REAL);
   NewVector(sutb, 0, DT_REAL);
   NewVector(tmp, 0, DT_REAL);
   NewVector(ta, 0, DT_REAL);
   NewVector(tx, 0, DT_REAL);
   NewVector(buf, 0, DT_REAL);
   NewVector(w, 0, DT_REAL);
   if (nrows <= 0 || ncols <= 0 || threshold < 0.0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   if (threshold == 0.0) {
      threshold = 1000 * machineepsilon;
   }
// Factorize A first
   svdfailed = !rmatrixsvd(a, nrows, ncols, 1, 2, 2, &sv, &u, &vt);
   zeroa = sv.xR[0] == 0.0;
   if (svdfailed || zeroa) {
      if (svdfailed) {
         *info = -4;
      } else {
         *info = 1;
      }
      ae_vector_set_length(x, ncols);
      for (i = 0; i < ncols; i++) {
         x->xR[i] = 0.0;
      }
      rep->n = ncols;
      rep->k = ncols;
      ae_matrix_set_length(&rep->cx, ncols, ncols);
      for (i = 0; i < ncols; i++) {
         for (j = 0; j < ncols; j++) {
            if (i == j) {
               rep->cx.xyR[i][j] = 1.0;
            } else {
               rep->cx.xyR[i][j] = 0.0;
            }
         }
      }
      rep->r2 = 0.0;
      ae_frame_leave();
      return;
   }
   nsv = imin2(ncols, nrows);
   if (nsv == ncols) {
      rep->r2 = sv.xR[nsv - 1] / sv.xR[0];
   } else {
      rep->r2 = 0.0;
   }
   rep->n = ncols;
   *info = 1;
// Iterative refinement of xc combined with solution:
// 1. xc = 0
// 2. calculate r = bc-A*xc using extra-precise dot product
// 3. solve A*y = r
// 4. update x:=x+r
// 5. goto 2
//
// This cycle is executed until one of two things happens:
// 1. maximum number of iterations reached
// 2. last iteration decreased error to the lower limit
   ae_vector_set_length(&utb, nsv);
   ae_vector_set_length(&sutb, nsv);
   ae_vector_set_length(x, ncols);
   ae_vector_set_length(&tmp, ncols);
   ae_vector_set_length(&ta, ncols + 1);
   ae_vector_set_length(&tx, ncols + 1);
   ae_vector_set_length(&buf, ncols + 1);
   for (i = 0; i < ncols; i++) {
      x->xR[i] = 0.0;
   }
   kernelidx = nsv;
   for (i = 0; i < nsv; i++) {
      if (sv.xR[i] <= threshold * sv.xR[0]) {
         kernelidx = i;
         break;
      }
   }
   rep->k = ncols - kernelidx;
   nrfs = directdensesolvers_densesolverrfsmaxv2(ncols, rep->r2);
   terminatenexttime = false;
   ae_vector_set_length(&rp, nrows);
   for (rfs = 0; rfs <= nrfs; rfs++) {
      if (terminatenexttime) {
         break;
      }
   // calculate right part
      if (rfs == 0) {
         ae_v_move(rp.xR, 1, b->xR, 1, nrows);
      } else {
         smallerr = true;
         for (i = 0; i < nrows; i++) {
            ae_v_move(ta.xR, 1, a->xyR[i], 1, ncols);
            ta.xR[ncols] = -1.0;
            ae_v_move(tx.xR, 1, x->xR, 1, ncols);
            tx.xR[ncols] = b->xR[i];
            xdot(&ta, &tx, ncols + 1, &buf, &v, &verr);
            rp.xR[i] = -v;
            smallerr = smallerr && SmallR(v, 4 * verr);
         }
         if (smallerr) {
            terminatenexttime = true;
         }
      }
   // solve A*dx = rp
      for (i = 0; i < ncols; i++) {
         tmp.xR[i] = 0.0;
      }
      for (i = 0; i < nsv; i++) {
         utb.xR[i] = 0.0;
      }
      for (i = 0; i < nrows; i++) {
         v = rp.xR[i];
         ae_v_addd(utb.xR, 1, u.xyR[i], 1, nsv, v);
      }
      for (i = 0; i < nsv; i++) {
         if (i < kernelidx) {
            sutb.xR[i] = utb.xR[i] / sv.xR[i];
         } else {
            sutb.xR[i] = 0.0;
         }
      }
      for (i = 0; i < nsv; i++) {
         v = sutb.xR[i];
         ae_v_addd(tmp.xR, 1, vt.xyR[i], 1, ncols, v);
      }
   // update x:  x:=x+dx
      ae_v_add(x->xR, 1, tmp.xR, 1, ncols);
   }
// fill CX
   if (rep->k > 0) {
      ae_matrix_set_length(&rep->cx, ncols, rep->k);
      for (i = 0; i < rep->k; i++) {
         ae_v_move(&rep->cx.xyR[0][i], rep->cx.stride, vt.xyR[kernelidx + i], 1, ncols);
      }
   }
   ae_frame_leave();
}

void densesolverreport_init(void *_p, bool make_automatic) {
}

void densesolverreport_copy(void *_dst, void *_src, bool make_automatic) {
   densesolverreport *dst = (densesolverreport *)_dst;
   densesolverreport *src = (densesolverreport *)_src;
   dst->r1 = src->r1;
   dst->rinf = src->rinf;
}

void densesolverreport_free(void *_p, bool make_automatic) {
}

void densesolverlsreport_init(void *_p, bool make_automatic) {
   densesolverlsreport *p = (densesolverlsreport *)_p;
   ae_matrix_init(&p->cx, 0, 0, DT_REAL, make_automatic);
}

void densesolverlsreport_copy(void *_dst, void *_src, bool make_automatic) {
   densesolverlsreport *dst = (densesolverlsreport *)_dst;
   densesolverlsreport *src = (densesolverlsreport *)_src;
   dst->r2 = src->r2;
   ae_matrix_copy(&dst->cx, &src->cx, make_automatic);
   dst->n = src->n;
   dst->k = src->k;
}

void densesolverlsreport_free(void *_p, bool make_automatic) {
   densesolverlsreport *p = (densesolverlsreport *)_p;
   ae_matrix_free(&p->cx, make_automatic);
}
} // end of namespace alglib_impl

namespace alglib {
DefClass(densesolverreport, AndD DecVal(r1) AndD DecVal(rinf))
DefClass(densesolverlsreport, AndD DecVal(r2) AndD DecVar(cx) AndD DecVal(n) AndD DecVal(k))

void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsolve(ConstT(ae_matrix, a), n, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsolvefast(ConstT(ae_matrix, a), n, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsolvem(ConstT(ae_matrix, a), n, ConstT(ae_matrix, b), m, rfs, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsolvemfast(ConstT(ae_matrix, a), n, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlusolve(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlusolvefast(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlusolvem(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlusolvemfast(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixmixedsolve(ConstT(ae_matrix, a), ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixmixedsolvem(ConstT(ae_matrix, a), ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixsolvem(ConstT(ae_matrix, a), n, ConstT(ae_matrix, b), m, rfs, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixsolvemfast(ConstT(ae_matrix, a), n, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixsolve(ConstT(ae_matrix, a), n, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixsolvefast(ConstT(ae_matrix, a), n, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlusolvem(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlusolvemfast(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlusolve(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlusolvefast(ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixmixedsolvem(ConstT(ae_matrix, a), ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixmixedsolve(ConstT(ae_matrix, a), ConstT(ae_matrix, lua), ConstT(ae_vector, p), n, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixsolvem(ConstT(ae_matrix, a), n, isupper, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixsolvemfast(ConstT(ae_matrix, a), n, isupper, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixsolve(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixsolvefast(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskysolvem(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskysolvemfast(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskysolve(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskysolvefast(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixsolvem(ConstT(ae_matrix, a), n, isupper, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixsolvemfast(ConstT(ae_matrix, a), n, isupper, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixsolve(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixsolvefast(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixcholeskysolvem(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_matrix, b), m, &info, ConstT(densesolverreport, rep), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixcholeskysolvemfast(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_matrix, b), m, &info);
   alglib_impl::ae_state_clear();
}

void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixcholeskysolve(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_vector, b), &info, ConstT(densesolverreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixcholeskysolvefast(ConstT(ae_matrix, cha), n, isupper, ConstT(ae_vector, b), &info);
   alglib_impl::ae_state_clear();
}

void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsolvels(ConstT(ae_matrix, a), nrows, ncols, ConstT(ae_vector, b), threshold, &info, ConstT(densesolverlsreport, rep), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === DIRECTSPARSESOLVERS Package ===
// Depends on: (LinAlg) TRFAC
namespace alglib_impl {
// Reset report fields
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
static void directsparsesolvers_initreport(sparsesolverreport *rep) {
   SetObj(sparsesolverreport, rep);
   rep->terminationtype = 0;
}

// Sparse linear solver for A*x=b with N*N  sparse  real  symmetric  positive
// definite matrix A, N*1 vectors x and b.
//
// This solver  converts  input  matrix  to  SKS  format,  performs  Cholesky
// factorization using  SKS  Cholesky  subroutine  (works  well  for  limited
// bandwidth matrices) and uses sparse triangular solvers to get solution  of
// the original system.
//
// Inputs:
//     A       -   sparse matrix, must be NxN exactly
//     IsUpper -   which half of A is provided (another half is ignored)
//     B       -   array[0..N-1], right part
//
// Outputs:
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0    =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate or non-SPD system).
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
// API: void sparsespdsolvesks(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsespdsolvesks(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n;
   ae_frame_make(&_frame_block);
   SetVector(x);
   SetObj(sparsesolverreport, rep);
   NewObj(sparsematrix, a2);
   n = sparsegetnrows(a);
   ae_assert(n > 0, "SparseSPDSolveSKS: N <= 0");
   ae_assert(sparsegetnrows(a) == n, "SparseSPDSolveSKS: rows(A) != N");
   ae_assert(sparsegetncols(a) == n, "SparseSPDSolveSKS: cols(A) != N");
   ae_assert(b->cnt >= n, "SparseSPDSolveSKS: length(B) < N");
   ae_assert(isfinitevector(b, n), "SparseSPDSolveSKS: B contains infinities or NANs");
   directsparsesolvers_initreport(rep);
   ae_vector_set_length(x, n);
   sparsecopytosks(a, &a2);
   if (!sparsecholeskyskyline(&a2, n, isupper)) {
      rep->terminationtype = -3;
      for (i = 0; i < n; i++) {
         x->xR[i] = 0.0;
      }
      ae_frame_leave();
      return;
   }
   for (i = 0; i < n; i++) {
      x->xR[i] = b->xR[i];
   }
   if (isupper) {
      sparsetrsv(&a2, isupper, false, 1, x);
      sparsetrsv(&a2, isupper, false, 0, x);
   } else {
      sparsetrsv(&a2, isupper, false, 0, x);
      sparsetrsv(&a2, isupper, false, 1, x);
   }
   rep->terminationtype = 1;
   ae_frame_leave();
}

// Sparse linear solver for A*x=b with N*N  sparse  real  symmetric  positive
// definite matrix A, N*1 vectors x and b.
//
// This solver  converts  input  matrix  to  CRS  format,  performs  Cholesky
// factorization using supernodal Cholesky  decomposition  with  permutation-
// reducing ordering and uses sparse triangular solver to get solution of the
// original system.
//
// Inputs:
//     A       -   sparse matrix, must be NxN exactly
//     IsUpper -   which half of A is provided (another half is ignored)
//     B       -   array[N], right part
//
// Outputs:
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0  =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate or non-SPD system).
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
// API: void sparsespdsolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsespdsolve(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n;
   ae_frame_make(&_frame_block);
   SetVector(x);
   SetObj(sparsesolverreport, rep);
   NewObj(sparsematrix, a2);
   NewVector(p, 0, DT_INT);
   n = sparsegetnrows(a);
   ae_assert(n > 0, "SparseSPDSolve: N <= 0");
   ae_assert(sparsegetnrows(a) == n, "SparseSPDSolve: rows(A) != N");
   ae_assert(sparsegetncols(a) == n, "SparseSPDSolve: cols(A) != N");
   ae_assert(b->cnt >= n, "SparseSPDSolve: length(B)<N");
   ae_assert(isfinitevector(b, n), "SparseSPDSolve: B contains infinities or NANs");
   directsparsesolvers_initreport(rep);
   sparsecopytocrs(a, &a2);
   if (!sparsecholeskyp(&a2, isupper, &p)) {
      rep->terminationtype = -3;
      rsetallocv(n, 0.0, x);
      ae_frame_leave();
      return;
   }
   rcopyallocv(n, b, x);
   for (i = 0; i < n; i++) swapr(&x->xR[i], &x->xR[p.xZ[i]]);
   if (isupper) {
      sparsetrsv(&a2, isupper, false, 1, x);
      sparsetrsv(&a2, isupper, false, 0, x);
   } else {
      sparsetrsv(&a2, isupper, false, 0, x);
      sparsetrsv(&a2, isupper, false, 1, x);
   }
   for (i = n - 1; i >= 0; i--) swapr(&x->xR[i], &x->xR[p.xZ[i]]);
   rep->terminationtype = 1;
   ae_frame_leave();
}

// Sparse linear solver for A*x=b with N*N real  symmetric  positive definite
// matrix A given by its Cholesky decomposition, and N*1 vectors x and b.
//
// IMPORTANT: this solver requires input matrix to be in  the  SKS  (Skyline)
//            or CRS (compressed row storage) format. An  exception  will  be
//            generated if you pass matrix in some other format.
//
// Inputs:
//     A       -   sparse NxN matrix stored in CRS or SKS format, must be NxN
//                 exactly
//     IsUpper -   which half of A is provided (another half is ignored)
//     B       -   array[N], right part
//
// Outputs:
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0    =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate or non-SPD system).
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
// API: void sparsespdcholeskysolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsespdcholeskysolve(sparsematrix *a, bool isupper, RVector *b, RVector *x, sparsesolverreport *rep) {
   ae_int_t i;
   ae_int_t n;
   SetVector(x);
   SetObj(sparsesolverreport, rep);
   n = sparsegetnrows(a);
   ae_assert(n > 0, "SparseSPDCholeskySolve: N <= 0");
   ae_assert(sparsegetnrows(a) == n, "SparseSPDCholeskySolve: rows(A) != N");
   ae_assert(sparsegetncols(a) == n, "SparseSPDCholeskySolve: cols(A) != N");
   ae_assert(sparseissks(a) || sparseiscrs(a), "SparseSPDCholeskySolve: A is not an SKS/CRS matrix");
   ae_assert(b->cnt >= n, "SparseSPDCholeskySolve: length(B)<N");
   ae_assert(isfinitevector(b, n), "SparseSPDCholeskySolve: B contains infinities or NANs");
   directsparsesolvers_initreport(rep);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (sparseget(a, i, i) == 0.0) {
         rep->terminationtype = -3;
         for (i = 0; i < n; i++) {
            x->xR[i] = 0.0;
         }
         return;
      }
   }
   for (i = 0; i < n; i++) {
      x->xR[i] = b->xR[i];
   }
   if (isupper) {
      sparsetrsv(a, isupper, false, 1, x);
      sparsetrsv(a, isupper, false, 0, x);
   } else {
      sparsetrsv(a, isupper, false, 0, x);
      sparsetrsv(a, isupper, false, 1, x);
   }
   rep->terminationtype = 1;
}

// Sparse linear solver for A*x=b with general (nonsymmetric) N*N sparse real
// matrix A, N*1 vectors x and b.
//
// This solver converts input matrix to CRS format, performs LU factorization
// and uses sparse triangular solvers to get solution of the original system.
//
// Inputs:
//     A       -   sparse matrix, must be NxN exactly, any storage format
//     B       -   array[0..N-1], right part
//
// Outputs:
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0    =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate system).
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
// API: void sparsesolve(const sparsematrix &a, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparsesolve(sparsematrix *a, RVector *b, RVector *x, sparsesolverreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n;
   ae_frame_make(&_frame_block);
   SetVector(x);
   SetObj(sparsesolverreport, rep);
   NewObj(sparsematrix, a2);
   NewVector(pivp, 0, DT_INT);
   NewVector(pivq, 0, DT_INT);
   n = sparsegetnrows(a);
   ae_assert(n > 0, "SparseSolve: N <= 0");
   ae_assert(sparsegetnrows(a) == n, "SparseSolve: rows(A) != N");
   ae_assert(sparsegetncols(a) == n, "SparseSolve: cols(A) != N");
   ae_assert(b->cnt >= n, "SparseSolve: length(B)<N");
   ae_assert(isfinitevector(b, n), "SparseSolve: B contains infinities or NANs");
   directsparsesolvers_initreport(rep);
   ae_vector_set_length(x, n);
   sparsecopytocrs(a, &a2);
   if (!sparselu(&a2, 0, &pivp, &pivq)) {
      rep->terminationtype = -3;
      for (i = 0; i < n; i++) {
         x->xR[i] = 0.0;
      }
      ae_frame_leave();
      return;
   }
   for (i = 0; i < n; i++) {
      x->xR[i] = b->xR[i];
   }
   for (i = 0; i < n; i++) {
      swapr(&x->xR[i], &x->xR[pivp.xZ[i]]);
   }
   sparsetrsv(&a2, false, true, 0, x);
   sparsetrsv(&a2, true, false, 0, x);
   for (i = n - 1; i >= 0; i--) {
      swapr(&x->xR[i], &x->xR[pivq.xZ[i]]);
   }
   rep->terminationtype = 1;
   ae_frame_leave();
}

// Sparse linear solver for A*x=b with general (nonsymmetric) N*N sparse real
// matrix A given by its LU factorization, N*1 vectors x and b.
//
// IMPORTANT: this solver requires input matrix  to  be  in  the  CRS  sparse
//            storage format. An exception will  be  generated  if  you  pass
//            matrix in some other format (HASH or SKS).
//
// Inputs:
//     A       -   LU factorization of the sparse matrix, must be NxN exactly
//                 in CRS storage format
//     P, Q    -   pivot indexes from LU factorization
//     B       -   array[0..N-1], right part
//
// Outputs:
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0    =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate system).
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
// API: void sparselusolve(const sparsematrix &a, const integer_1d_array &p, const integer_1d_array &q, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
void sparselusolve(sparsematrix *a, ZVector *p, ZVector *q, RVector *b, RVector *x, sparsesolverreport *rep) {
   ae_int_t i;
   ae_int_t n;
   SetVector(x);
   SetObj(sparsesolverreport, rep);
   n = sparsegetnrows(a);
   ae_assert(n > 0, "SparseLUSolve: N <= 0");
   ae_assert(sparsegetnrows(a) == n, "SparseLUSolve: rows(A) != N");
   ae_assert(sparsegetncols(a) == n, "SparseLUSolve: cols(A) != N");
   ae_assert(sparseiscrs(a), "SparseLUSolve: A is not an SKS matrix");
   ae_assert(b->cnt >= n, "SparseLUSolve: length(B)<N");
   ae_assert(isfinitevector(b, n), "SparseLUSolve: B contains infinities or NANs");
   ae_assert(p->cnt >= n, "SparseLUSolve: length(P)<N");
   ae_assert(q->cnt >= n, "SparseLUSolve: length(Q)<N");
   for (i = 0; i < n; i++) {
      ae_assert(p->xZ[i] >= i && p->xZ[i] < n, "SparseLUSolve: P is corrupted");
      ae_assert(q->xZ[i] >= i && q->xZ[i] < n, "SparseLUSolve: Q is corrupted");
   }
   directsparsesolvers_initreport(rep);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (a->didx.xZ[i] == a->uidx.xZ[i] || a->vals.xR[a->didx.xZ[i]] == 0.0) {
         rep->terminationtype = -3;
         for (i = 0; i < n; i++) {
            x->xR[i] = 0.0;
         }
         return;
      }
   }
   for (i = 0; i < n; i++) {
      x->xR[i] = b->xR[i];
   }
   for (i = 0; i < n; i++) {
      swapr(&x->xR[i], &x->xR[p->xZ[i]]);
   }
   sparsetrsv(a, false, true, 0, x);
   sparsetrsv(a, true, false, 0, x);
   for (i = n - 1; i >= 0; i--) {
      swapr(&x->xR[i], &x->xR[q->xZ[i]]);
   }
   rep->terminationtype = 1;
}

void sparsesolverreport_init(void *_p, bool make_automatic) {
}

void sparsesolverreport_copy(void *_dst, void *_src, bool make_automatic) {
   sparsesolverreport *dst = (sparsesolverreport *)_dst;
   sparsesolverreport *src = (sparsesolverreport *)_src;
   dst->terminationtype = src->terminationtype;
}

void sparsesolverreport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
// This structure is a sparse solver report.
//
// Following fields can be accessed by users:
DefClass(sparsesolverreport, AndD DecVal(terminationtype))

void sparsespdsolvesks(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsespdsolvesks(ConstT(sparsematrix, a), isupper, ConstT(ae_vector, b), ConstT(ae_vector, x), ConstT(sparsesolverreport, rep));
   alglib_impl::ae_state_clear();
}

void sparsespdsolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsespdsolve(ConstT(sparsematrix, a), isupper, ConstT(ae_vector, b), ConstT(ae_vector, x), ConstT(sparsesolverreport, rep));
   alglib_impl::ae_state_clear();
}

void sparsespdcholeskysolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsespdcholeskysolve(ConstT(sparsematrix, a), isupper, ConstT(ae_vector, b), ConstT(ae_vector, x), ConstT(sparsesolverreport, rep));
   alglib_impl::ae_state_clear();
}

void sparsesolve(const sparsematrix &a, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsesolve(ConstT(sparsematrix, a), ConstT(ae_vector, b), ConstT(ae_vector, x), ConstT(sparsesolverreport, rep));
   alglib_impl::ae_state_clear();
}

void sparselusolve(const sparsematrix &a, const integer_1d_array &p, const integer_1d_array &q, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparselusolve(ConstT(sparsematrix, a), ConstT(ae_vector, p), ConstT(ae_vector, q), ConstT(ae_vector, b), ConstT(ae_vector, x), ConstT(sparsesolverreport, rep));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === LINCG Package ===
// Depends on: (LinAlg) SPARSE, MATGEN
namespace alglib_impl {
static const double lincg_defaultprecision = 1.0E-6;

// Clears request fields (to be sure that we don't forgot to clear something)
static void lincg_updateitersdata(lincgstate *state) {
   state->repiterationscount = 0;
   state->repnmv = 0;
   state->repterminationtype = 0;
}

// This function initializes linear CG Solver. This solver is used  to  solve
// symmetric positive definite problems. If you want  to  solve  nonsymmetric
// (or non-positive definite) problem you may use LinLSQR solver provided  by
// ALGLIB.
//
// USAGE:
// 1. User initializes algorithm state with LinCGCreate() call
// 2. User tunes solver parameters with  LinCGSetCond() and other functions
// 3. Optionally, user sets starting point with LinCGSetStartingPoint()
// 4. User  calls LinCGSolveSparse() function which takes algorithm state and
//    SparseMatrix object.
// 5. User calls LinCGResults() to get solution
// 6. Optionally, user may call LinCGSolveSparse()  again  to  solve  another
//    problem  with different matrix and/or right part without reinitializing
//    LinCGState structure.
//
// Inputs:
//     N       -   problem dimension, N > 0
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgcreate(const ae_int_t n, lincgstate &state);
void lincgcreate(ae_int_t n, lincgstate *state) {
   ae_int_t i;
   SetObj(lincgstate, state);
   ae_assert(n > 0, "LinCGCreate: N <= 0");
   state->n = n;
   state->prectype = 0;
   state->itsbeforerestart = n;
   state->itsbeforerupdate = 10;
   state->epsf = lincg_defaultprecision;
   state->maxits = 0;
   state->xrep = false;
   state->running = false;
// * allocate arrays
// * set RX to NAN (just for the case user calls Results() without
//   calling SolveSparse()
// * set starting point to zero
// * we do NOT initialize B here because we assume that user should
//   initializate it using LinCGSetB() function. In case he forgets
//   to do so, exception will be thrown in the LinCGIteration().
   ae_vector_set_length(&state->rx, state->n);
   ae_vector_set_length(&state->startx, state->n);
   ae_vector_set_length(&state->b, state->n);
   for (i = 0; i < state->n; i++) {
      state->rx.xR[i] = NAN;
      state->startx.xR[i] = 0.0;
      state->b.xR[i] = 0.0;
   }
   ae_vector_set_length(&state->cx, state->n);
   ae_vector_set_length(&state->p, state->n);
   ae_vector_set_length(&state->r, state->n);
   ae_vector_set_length(&state->cr, state->n);
   ae_vector_set_length(&state->z, state->n);
   ae_vector_set_length(&state->cz, state->n);
   ae_vector_set_length(&state->x, state->n);
   ae_vector_set_length(&state->mv, state->n);
   ae_vector_set_length(&state->pv, state->n);
   lincg_updateitersdata(state);
   state->PQ = -1;
}

// This function sets starting point.
// By default, zero starting point is used.
//
// Inputs:
//     X       -   starting point, array[N]
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x);
void lincgsetstartingpoint(lincgstate *state, RVector *x) {
   ae_assert(!state->running, "LinCGSetStartingPoint: you can not change starting point because LinCGIteration() function is running");
   ae_assert(state->n <= x->cnt, "LinCGSetStartingPoint: Length(X)<N");
   ae_assert(isfinitevector(x, state->n), "LinCGSetStartingPoint: X contains infinite or NaN values!");
   ae_v_move(state->startx.xR, 1, x->xR, 1, state->n);
}

// This function sets right part. By default, right part is zero.
//
// Inputs:
//     B       -   right part, array[N].
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsetb(const lincgstate &state, const real_1d_array &b);
void lincgsetb(lincgstate *state, RVector *b) {
   ae_assert(!state->running, "LinCGSetB: you can not set B, because function LinCGIteration is running!");
   ae_assert(b->cnt >= state->n, "LinCGSetB: Length(B)<N");
   ae_assert(isfinitevector(b, state->n), "LinCGSetB: B contains infinite or NaN values!");
   ae_v_move(state->b.xR, 1, b->xR, 1, state->n);
}

// This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
// function. By default, SolveSparse() uses diagonal preconditioner,  but  if
// you want to use solver without preconditioning, you can call this function
// which forces solver to use unit matrix for preconditioning.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
// API: void lincgsetprecunit(const lincgstate &state);
void lincgsetprecunit(lincgstate *state) {
   ae_assert(!state->running, "LinCGSetPrecUnit: you can not change preconditioner, because function LinCGIteration is running!");
   state->prectype = -1;
}

// This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
// function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
// preconditioner. This preconditioning mode is active by default.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
// API: void lincgsetprecdiag(const lincgstate &state);
void lincgsetprecdiag(lincgstate *state) {
   ae_assert(!state->running, "LinCGSetPrecDiag: you can not change preconditioner, because function LinCGIteration is running!");
   state->prectype = 0;
}

// This function sets stopping criteria.
//
// Inputs:
//     EpsF    -   algorithm will be stopped if norm of residual is less than
//                 EpsF*||b||.
//     MaxIts  -   algorithm will be stopped if number of iterations is  more
//                 than MaxIts.
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// NOTES:
// If  both  EpsF  and  MaxIts  are  zero then small EpsF will be set to small
// value.
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits);
void lincgsetcond(lincgstate *state, double epsf, ae_int_t maxits) {
   ae_assert(!state->running, "LinCGSetCond: you can not change stopping criteria when LinCGIteration() is running");
   ae_assert(isfinite(epsf) && epsf >= 0.0, "LinCGSetCond: EpsF is negative or contains infinite or NaN values");
   ae_assert(maxits >= 0, "LinCGSetCond: MaxIts is negative");
   if (epsf == 0.0 && maxits == 0) {
      state->epsf = lincg_defaultprecision;
      state->maxits = maxits;
   } else {
      state->epsf = epsf;
      state->maxits = maxits;
   }
}

// Reverse communication version of linear CG.
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
bool lincgiteration(lincgstate *state) {
   AutoS ae_int_t i;
   AutoS double uvar;
   AutoS double bnorm;
   AutoS double v;
// Manually threaded two-way signalling.
// Locals are set arbitrarily the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0; case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      case 4: goto Resume4; case 5: goto Resume5; case 6: goto Resume6; case 7: goto Resume7;
      default: goto Exit;
   }
Spawn:
   i = 359;
   uvar = -58;
   bnorm = -919;
   v = -909;
   ae_assert(state->b.cnt > 0, "LinCGIteration: B is not initialized (you must initialize B by LinCGSetB() call");
   state->needprec = state->needvmv = state->needmv = state->xupdated = false;
   state->running = true;
   state->repnmv = 0;
   lincg_updateitersdata(state);
// Start 0-th iteration
   ae_v_move(state->rx.xR, 1, state->startx.xR, 1, state->n);
   ae_v_move(state->x.xR, 1, state->rx.xR, 1, state->n);
   state->repnmv++, state->needvmv = true, state->PQ = 0; goto Pause; Resume0: state->needvmv = false;
   bnorm = 0.0;
   state->r2 = 0.0;
   state->meritfunction = 0.0;
   for (i = 0; i < state->n; i++) {
      state->r.xR[i] = state->b.xR[i] - state->mv.xR[i];
      state->r2 += state->r.xR[i] * state->r.xR[i];
      state->meritfunction += state->mv.xR[i] * state->rx.xR[i] - 2 * state->b.xR[i] * state->rx.xR[i];
      bnorm += state->b.xR[i] * state->b.xR[i];
   }
   bnorm = sqrt(bnorm);
// Output first report
   if (state->xrep) {
      ae_v_move(state->x.xR, 1, state->rx.xR, 1, state->n);
      state->xupdated = true, state->PQ = 1; goto Pause; Resume1: state->xupdated = false;
   }
// Is x0 a solution?
   if (!isfinite(state->r2) || sqrt(state->r2) <= state->epsf * bnorm) {
      state->running = false;
      if (isfinite(state->r2)) {
         state->repterminationtype = 1;
      } else {
         state->repterminationtype = -4;
      }
      goto Exit;
   }
// Calculate Z and P
   ae_v_move(state->x.xR, 1, state->r.xR, 1, state->n);
   state->repnmv++, state->needprec = true, state->PQ = 2; goto Pause; Resume2: state->needprec = false;
   for (i = 0; i < state->n; i++) {
      state->z.xR[i] = state->pv.xR[i];
      state->p.xR[i] = state->z.xR[i];
   }
// Other iterations(1..N)
   state->repiterationscount = 0;
   while (true) {
      state->repiterationscount++;
   // Calculate Alpha
      ae_v_move(state->x.xR, 1, state->p.xR, 1, state->n);
      state->repnmv++, state->needvmv = true, state->PQ = 3; goto Pause; Resume3: state->needvmv = false;
      if (!isfinite(state->vmv) || state->vmv <= 0.0) {
      // a) Overflow when calculating VMV
      // b) non-positive VMV (non-SPD matrix)
         state->running = false;
         if (isfinite(state->vmv)) {
            state->repterminationtype = -5;
         } else {
            state->repterminationtype = -4;
         }
         goto Exit;
      }
      state->alpha = 0.0;
      for (i = 0; i < state->n; i++) {
         state->alpha += state->r.xR[i] * state->z.xR[i];
      }
      state->alpha /= state->vmv;
      if (!isfinite(state->alpha)) {
      // Overflow when calculating Alpha
         state->running = false;
         state->repterminationtype = -4;
         goto Exit;
      }
   // Next step toward solution
      for (i = 0; i < state->n; i++) {
         state->cx.xR[i] = state->rx.xR[i] + state->alpha * state->p.xR[i];
      }
   // Calculate R:
   // * use recurrent relation to update R
   // * at every ItsBeforeRUpdate-th iteration recalculate it from scratch, using matrix-vector product
   //   in case R grows instead of decreasing, algorithm is terminated with positive completion code
      if (state->itsbeforerupdate == 0 || state->repiterationscount % state->itsbeforerupdate != 0) {
      // Calculate R using recurrent formula
         for (i = 0; i < state->n; i++) {
            state->cr.xR[i] = state->r.xR[i] - state->alpha * state->mv.xR[i];
            state->x.xR[i] = state->cr.xR[i];
         }
      } else {
      // Calculate R using matrix-vector multiplication
         ae_v_move(state->x.xR, 1, state->cx.xR, 1, state->n);
         state->repnmv++, state->needmv = true, state->PQ = 4; goto Pause; Resume4: state->needmv = false;
         for (i = 0; i < state->n; i++) {
            state->cr.xR[i] = state->b.xR[i] - state->mv.xR[i];
            state->x.xR[i] = state->cr.xR[i];
         }
      // Calculating merit function
      // Check emergency stopping criterion
         v = 0.0;
         for (i = 0; i < state->n; i++) {
            v += state->mv.xR[i] * state->cx.xR[i] - 2 * state->b.xR[i] * state->cx.xR[i];
         }
         if (v >= state->meritfunction) {
            for (i = 0; i < state->n; i++) {
               if (!isfinite(state->rx.xR[i])) {
                  state->running = false;
                  state->repterminationtype = -4;
                  goto Exit;
               }
            }
         // output last report
            if (state->xrep) {
               ae_v_move(state->x.xR, 1, state->rx.xR, 1, state->n);
               state->xupdated = true, state->PQ = 5; goto Pause; Resume5: state->xupdated = false;
            }
            state->running = false;
            state->repterminationtype = 7;
            goto Exit;
         }
         state->meritfunction = v;
      }
      ae_v_move(state->rx.xR, 1, state->cx.xR, 1, state->n);
   // calculating RNorm
   //
   // NOTE: monotonic decrease of R2 is not guaranteed by algorithm.
      state->r2 = 0.0;
      for (i = 0; i < state->n; i++) {
         state->r2 += state->cr.xR[i] * state->cr.xR[i];
      }
   // output report
      if (state->xrep) {
         ae_v_move(state->x.xR, 1, state->rx.xR, 1, state->n);
         state->xupdated = true, state->PQ = 6; goto Pause; Resume6: state->xupdated = false;
      }
   // stopping criterion
   // achieved the required precision
      if (!isfinite(state->r2) || sqrt(state->r2) <= state->epsf * bnorm) {
         state->running = false;
         if (isfinite(state->r2)) {
            state->repterminationtype = 1;
         } else {
            state->repterminationtype = -4;
         }
         goto Exit;
      }
      if (state->repiterationscount >= state->maxits && state->maxits > 0) {
         for (i = 0; i < state->n; i++) {
            if (!isfinite(state->rx.xR[i])) {
               state->running = false;
               state->repterminationtype = -4;
               goto Exit;
            }
         }
      // if X is finite number
         state->running = false;
         state->repterminationtype = 5;
         goto Exit;
      }
      ae_v_move(state->x.xR, 1, state->cr.xR, 1, state->n);
   // prepere of parameters for next iteration
      state->repnmv++, state->needprec = true, state->PQ = 7; goto Pause; Resume7: state->needprec = false;
      ae_v_move(state->cz.xR, 1, state->pv.xR, 1, state->n);
      if (state->repiterationscount % state->itsbeforerestart != 0) {
         state->beta = 0.0;
         uvar = 0.0;
         for (i = 0; i < state->n; i++) {
            state->beta += state->cz.xR[i] * state->cr.xR[i];
            uvar += state->z.xR[i] * state->r.xR[i];
         }
      // check that UVar is't INF or is't zero
         if (!isfinite(uvar) || uvar == 0.0) {
            state->running = false;
            state->repterminationtype = -4;
            goto Exit;
         }
      // calculate .BETA
         state->beta /= uvar;
      // check that .BETA neither INF nor NaN
         if (!isfinite(state->beta)) {
            state->running = false;
            state->repterminationtype = -1;
            goto Exit;
         }
         for (i = 0; i < state->n; i++) {
            state->p.xR[i] = state->cz.xR[i] + state->beta * state->p.xR[i];
         }
      } else {
         ae_v_move(state->p.xR, 1, state->cz.xR, 1, state->n);
      }
   // prepere data for next iteration
      for (i = 0; i < state->n; i++) {
      // write (k+1)th iteration to (k )th iteration
         state->r.xR[i] = state->cr.xR[i];
         state->z.xR[i] = state->cz.xR[i];
      }
   }
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// Procedure for solution of A*x=b with sparse A.
//
// Inputs:
//     State   -   algorithm state
//     A       -   sparse matrix in the CRS format (you MUST contvert  it  to
//                 CRS format by calling SparseConvertToCRS() function).
//     IsUpper -   whether upper or lower triangle of A is used:
//                 * IsUpper=True  => only upper triangle is used and lower
//                                    triangle is not referenced at all
//                 * IsUpper=False => only lower triangle is used and upper
//                                    triangle is not referenced at all
//     B       -   right part, array[N]
//
// Result:
//     This function returns no result.
//     You can get solution by calling LinCGResults()
//
// NOTE: this function uses lightweight preconditioning -  multiplication  by
//       inverse of diag(A). If you want, you can turn preconditioning off by
//       calling LinCGSetPrecUnit(). However, preconditioning cost is low and
//       preconditioner  is  very  important  for  solution  of  badly scaled
//       problems.
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b);
void lincgsolvesparse(lincgstate *state, sparsematrix *a, bool isupper, RVector *b) {
   ae_int_t n;
   ae_int_t i;
   double v;
   double vmv;
   n = state->n;
   ae_assert(b->cnt >= state->n, "LinCGSetB: Length(B)<N");
   ae_assert(isfinitevector(b, state->n), "LinCGSetB: B contains infinite or NaN values!");
// Allocate temporaries
   vectorsetlengthatleast(&state->tmpd, n);
// Compute diagonal scaling matrix D
   if (state->prectype == 0) {
   // Default preconditioner - inverse of matrix diagonal
      for (i = 0; i < n; i++) {
         v = sparsegetdiagonal(a, i);
         if (v > 0.0) {
            state->tmpd.xR[i] = 1 / sqrt(v);
         } else {
            state->tmpd.xR[i] = 1.0;
         }
      }
   } else {
   // No diagonal scaling
      for (i = 0; i < n; i++) {
         state->tmpd.xR[i] = 1.0;
      }
   }
// Solve
   lincgrestart(state);
   lincgsetb(state, b);
   while (lincgiteration(state)) {
   // Process different requests from optimizer
      if (state->needmv) {
         sparsesmv(a, isupper, &state->x, &state->mv);
      } else if (state->needvmv) {
         sparsesmv(a, isupper, &state->x, &state->mv);
         vmv = ae_v_dotproduct(state->x.xR, 1, state->mv.xR, 1, state->n);
         state->vmv = vmv;
      } else if (state->needprec) {
         for (i = 0; i < n; i++) {
            state->pv.xR[i] = state->x.xR[i] * sqr(state->tmpd.xR[i]);
         }
      }
   }
}

// CG-solver: results.
//
// This function must be called after LinCGSolve
//
// Inputs:
//     State   -   algorithm state
//
// Outputs:
//     X       -   array[N], solution
//     Rep     -   optimization report:
//                 * Rep.TerminationType completetion code:
//                     * -5    input matrix is either not positive definite,
//                             too large or too small
//                     * -4    overflow/underflow during solution
//                             (ill conditioned problem)
//                     *  1    ||residual|| <= EpsF*||b||
//                     *  5    MaxIts steps was taken
//                     *  7    rounding errors prevent further progress,
//                             best point found is returned
//                 * Rep.IterationsCount contains iterations count
//                 * NMV countains number of matrix-vector calculations
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep);
void lincgresults(lincgstate *state, RVector *x, lincgreport *rep) {
   SetVector(x);
   SetObj(lincgreport, rep);
   ae_assert(!state->running, "LinCGResult: you can not get result, because function LinCGIteration has been launched!");
   if (x->cnt < state->n) {
      ae_vector_set_length(x, state->n);
   }
   ae_v_move(x->xR, 1, state->rx.xR, 1, state->n);
   rep->iterationscount = state->repiterationscount;
   rep->nmv = state->repnmv;
   rep->terminationtype = state->repterminationtype;
   rep->r2 = state->r2;
}

// This function sets restart frequency. By default, algorithm  is  restarted
// after N subsequent iterations.
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf);
void lincgsetrestartfreq(lincgstate *state, ae_int_t srf) {
   ae_assert(!state->running, "LinCGSetRestartFreq: you can not change restart frequency when LinCGIteration() is running");
   ae_assert(srf > 0, "LinCGSetRestartFreq: non-positive SRF");
   state->itsbeforerestart = srf;
}

// This function sets frequency of residual recalculations.
//
// Algorithm updates residual r_k using iterative formula,  but  recalculates
// it from scratch after each 10 iterations. It is done to avoid accumulation
// of numerical errors and to stop algorithm when r_k starts to grow.
//
// Such low update frequence (1/10) gives very  little  overhead,  but  makes
// algorithm a bit more robust against numerical errors. However, you may
// change it
//
// Inputs:
//     Freq    -   desired update frequency, Freq >= 0.
//                 Zero value means that no updates will be done.
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq);
void lincgsetrupdatefreq(lincgstate *state, ae_int_t freq) {
   ae_assert(!state->running, "LinCGSetRUpdateFreq: you can not change update frequency when LinCGIteration() is running");
   ae_assert(freq >= 0, "LinCGSetRUpdateFreq: non-positive Freq");
   state->itsbeforerupdate = freq;
}

// This function turns on/off reporting.
//
// Inputs:
//     State   -   structure which stores algorithm state
//     NeedXRep-   whether iteration reports are needed or not
//
// If NeedXRep is True, algorithm will call rep() callback function if  it is
// provided to MinCGOptimize().
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
// API: void lincgsetxrep(const lincgstate &state, const bool needxrep);
void lincgsetxrep(lincgstate *state, bool needxrep) {
   state->xrep = needxrep;
}

// Procedure for restart function LinCGIteration
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
void lincgrestart(lincgstate *state) {
   state->PQ = -1;
}

void lincgstate_init(void *_p, bool make_automatic) {
   lincgstate *p = (lincgstate *)_p;
   ae_vector_init(&p->rx, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->b, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->cx, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->cr, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->cz, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->p, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->r, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->z, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->mv, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->pv, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->startx, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmpd, 0, DT_REAL, make_automatic);
}

void lincgstate_copy(void *_dst, void *_src, bool make_automatic) {
   lincgstate *dst = (lincgstate *)_dst;
   lincgstate *src = (lincgstate *)_src;
   ae_vector_copy(&dst->rx, &src->rx, make_automatic);
   ae_vector_copy(&dst->b, &src->b, make_automatic);
   dst->n = src->n;
   dst->prectype = src->prectype;
   ae_vector_copy(&dst->cx, &src->cx, make_automatic);
   ae_vector_copy(&dst->cr, &src->cr, make_automatic);
   ae_vector_copy(&dst->cz, &src->cz, make_automatic);
   ae_vector_copy(&dst->p, &src->p, make_automatic);
   ae_vector_copy(&dst->r, &src->r, make_automatic);
   ae_vector_copy(&dst->z, &src->z, make_automatic);
   dst->alpha = src->alpha;
   dst->beta = src->beta;
   dst->r2 = src->r2;
   dst->meritfunction = src->meritfunction;
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   ae_vector_copy(&dst->mv, &src->mv, make_automatic);
   ae_vector_copy(&dst->pv, &src->pv, make_automatic);
   dst->vmv = src->vmv;
   ae_vector_copy(&dst->startx, &src->startx, make_automatic);
   dst->epsf = src->epsf;
   dst->maxits = src->maxits;
   dst->itsbeforerestart = src->itsbeforerestart;
   dst->itsbeforerupdate = src->itsbeforerupdate;
   dst->xrep = src->xrep;
   dst->xupdated = src->xupdated;
   dst->needmv = src->needmv;
   dst->needvmv = src->needvmv;
   dst->needprec = src->needprec;
   dst->repiterationscount = src->repiterationscount;
   dst->repnmv = src->repnmv;
   dst->repterminationtype = src->repterminationtype;
   dst->running = src->running;
   ae_vector_copy(&dst->tmpd, &src->tmpd, make_automatic);
   dst->PQ = src->PQ;
}

void lincgstate_free(void *_p, bool make_automatic) {
   lincgstate *p = (lincgstate *)_p;
   ae_vector_free(&p->rx, make_automatic);
   ae_vector_free(&p->b, make_automatic);
   ae_vector_free(&p->cx, make_automatic);
   ae_vector_free(&p->cr, make_automatic);
   ae_vector_free(&p->cz, make_automatic);
   ae_vector_free(&p->p, make_automatic);
   ae_vector_free(&p->r, make_automatic);
   ae_vector_free(&p->z, make_automatic);
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->mv, make_automatic);
   ae_vector_free(&p->pv, make_automatic);
   ae_vector_free(&p->startx, make_automatic);
   ae_vector_free(&p->tmpd, make_automatic);
}

void lincgreport_init(void *_p, bool make_automatic) {
}

void lincgreport_copy(void *_dst, void *_src, bool make_automatic) {
   lincgreport *dst = (lincgreport *)_dst;
   lincgreport *src = (lincgreport *)_src;
   dst->iterationscount = src->iterationscount;
   dst->nmv = src->nmv;
   dst->terminationtype = src->terminationtype;
   dst->r2 = src->r2;
}

void lincgreport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the linear CG method.
// You should use ALGLIB functions to work with this object.
// Never try to access its fields directly!
DefClass(lincgstate, EndD)
DefClass(lincgreport, AndD DecVal(iterationscount) AndD DecVal(nmv) AndD DecVal(terminationtype) AndD DecVal(r2))

void lincgcreate(const ae_int_t n, lincgstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgcreate(n, ConstT(lincgstate, state));
   alglib_impl::ae_state_clear();
}

void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetstartingpoint(ConstT(lincgstate, state), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void lincgsetb(const lincgstate &state, const real_1d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetb(ConstT(lincgstate, state), ConstT(ae_vector, b));
   alglib_impl::ae_state_clear();
}

void lincgsetprecunit(const lincgstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetprecunit(ConstT(lincgstate, state));
   alglib_impl::ae_state_clear();
}

void lincgsetprecdiag(const lincgstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetprecdiag(ConstT(lincgstate, state));
   alglib_impl::ae_state_clear();
}

void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetcond(ConstT(lincgstate, state), epsf, maxits);
   alglib_impl::ae_state_clear();
}

void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsolvesparse(ConstT(lincgstate, state), ConstT(sparsematrix, a), isupper, ConstT(ae_vector, b));
   alglib_impl::ae_state_clear();
}

void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgresults(ConstT(lincgstate, state), ConstT(ae_vector, x), ConstT(lincgreport, rep));
   alglib_impl::ae_state_clear();
}

void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetrestartfreq(ConstT(lincgstate, state), srf);
   alglib_impl::ae_state_clear();
}

void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetrupdatefreq(ConstT(lincgstate, state), freq);
   alglib_impl::ae_state_clear();
}

void lincgsetxrep(const lincgstate &state, const bool needxrep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::lincgsetxrep(ConstT(lincgstate, state), needxrep);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === LINLSQR Package ===
// Depends on: (LinAlg) SVD, NORMESTIMATOR
namespace alglib_impl {
static const double linlsqr_atol = 1.0E-6;
static const double linlsqr_btol = 1.0E-6;

// This function initializes linear LSQR Solver. This solver is used to solve
// non-symmetric (and, possibly, non-square) problems. Least squares solution
// is returned for non-compatible systems.
//
// USAGE:
// 1. User initializes algorithm state with LinLSQRCreate() call
// 2. User tunes solver parameters with  LinLSQRSetCond() and other functions
// 3. User  calls  LinLSQRSolveSparse()  function which takes algorithm state
//    and SparseMatrix object.
// 4. User calls LinLSQRResults() to get solution
// 5. Optionally, user may call LinLSQRSolveSparse() again to  solve  another
//    problem  with different matrix and/or right part without reinitializing
//    LinLSQRState structure.
//
// Inputs:
//     M       -   number of rows in A
//     N       -   number of variables, N > 0
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// NOTE: see also linlsqrcreatebuf()  for  version  which  reuses  previously
//       allocated place as much as possible.
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
// API: void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state);
void linlsqrcreate(ae_int_t m, ae_int_t n, linlsqrstate *state) {
   SetObj(linlsqrstate, state);
   ae_assert(m > 0, "LinLSQRCreate: M <= 0");
   ae_assert(n > 0, "LinLSQRCreate: N <= 0");
   linlsqrcreatebuf(m, n, state);
}

// This function initializes linear LSQR Solver.  It  provides  exactly  same
// functionality as linlsqrcreate(), but reuses  previously  allocated  space
// as much as possible.
//
// Inputs:
//     M       -   number of rows in A
//     N       -   number of variables, N > 0
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 14.11.2018 by Sergey Bochkanov
// API: void linlsqrcreatebuf(const ae_int_t m, const ae_int_t n, const linlsqrstate &state);
void linlsqrcreatebuf(ae_int_t m, ae_int_t n, linlsqrstate *state) {
   ae_int_t i;
   ae_assert(m > 0, "LinLSQRCreateBuf: M <= 0");
   ae_assert(n > 0, "LinLSQRCreateBuf: N <= 0");
   state->m = m;
   state->n = n;
   state->prectype = 0;
   state->epsa = linlsqr_atol;
   state->epsb = linlsqr_btol;
   state->epsc = 1 / sqrt(machineepsilon);
   state->maxits = 0;
   state->lambdai = 0.0;
   state->xrep = false;
   state->running = false;
   state->repiterationscount = 0;
// * allocate arrays
// * set RX to NAN (just for the case user calls Results() without
//   calling SolveSparse()
// * set B to zero
   normestimatorcreate(m, n, 2, 2, &state->nes);
   ae_vector_set_length(&state->rx, state->n);
   ae_vector_set_length(&state->ui, state->m + state->n);
   ae_vector_set_length(&state->uip1, state->m + state->n);
   ae_vector_set_length(&state->vip1, state->n);
   ae_vector_set_length(&state->vi, state->n);
   ae_vector_set_length(&state->omegai, state->n);
   ae_vector_set_length(&state->omegaip1, state->n);
   ae_vector_set_length(&state->d, state->n);
   ae_vector_set_length(&state->x, state->m + state->n);
   ae_vector_set_length(&state->mv, state->m + state->n);
   ae_vector_set_length(&state->mtv, state->n);
   ae_vector_set_length(&state->b, state->m);
   for (i = 0; i < n; i++) {
      state->rx.xR[i] = NAN;
   }
   for (i = 0; i < m; i++) {
      state->b.xR[i] = 0.0;
   }
   state->PQ = -1;
}

// This function sets right part. By default, right part is zero.
//
// Inputs:
//     B       -   right part, array[N].
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
void linlsqrsetb(linlsqrstate *state, RVector *b) {
   ae_int_t i;
   ae_assert(!state->running, "LinLSQRSetB: you can not change B when LinLSQRIteration is running");
   ae_assert(state->m <= b->cnt, "LinLSQRSetB: Length(B)<M");
   ae_assert(isfinitevector(b, state->m), "LinLSQRSetB: B contains infinite or NaN values");
   state->bnorm2 = 0.0;
   for (i = 0; i < state->m; i++) {
      state->b.xR[i] = b->xR[i];
      state->bnorm2 += b->xR[i] * b->xR[i];
   }
}

// This  function  changes  preconditioning  settings of LinLSQQSolveSparse()
// function. By default, SolveSparse() uses diagonal preconditioner,  but  if
// you want to use solver without preconditioning, you can call this function
// which forces solver to use unit matrix for preconditioning.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
// API: void linlsqrsetprecunit(const linlsqrstate &state);
void linlsqrsetprecunit(linlsqrstate *state) {
   ae_assert(!state->running, "LinLSQRSetPrecUnit: you can not change preconditioner, because function LinLSQRIteration is running!");
   state->prectype = -1;
}

// This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
// function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
// preconditioner. This preconditioning mode is active by default.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
// API: void linlsqrsetprecdiag(const linlsqrstate &state);
void linlsqrsetprecdiag(linlsqrstate *state) {
   ae_assert(!state->running, "LinLSQRSetPrecDiag: you can not change preconditioner, because function LinCGIteration is running!");
   state->prectype = 0;
}

// This function sets optional Tikhonov regularization coefficient.
// It is zero by default.
//
// Inputs:
//     LambdaI -   regularization factor, LambdaI >= 0
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
// API: void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai);
void linlsqrsetlambdai(linlsqrstate *state, double lambdai) {
   ae_assert(!state->running, "LinLSQRSetLambdaI: you can not set LambdaI, because function LinLSQRIteration is running");
   ae_assert(isfinite(lambdai) && lambdai >= 0.0, "LinLSQRSetLambdaI: LambdaI is infinite or NaN");
   state->lambdai = lambdai;
}

// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
bool linlsqriteration(linlsqrstate *state) {
   AutoS ae_int_t summn;
   AutoS double bnorm;
   AutoS ae_int_t i;
// Manually threaded two-way signalling.
// Locals are set arbitrarily the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0; case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      case 4: goto Resume4; case 5: goto Resume5; case 6: goto Resume6;
      default: goto Exit;
   }
Spawn:
   i = -58;
   ae_assert(state->b.cnt > 0, "LinLSQRIteration: using non-allocated array B");
   summn = state->m + state->n;
   bnorm = sqrt(state->bnorm2);
   state->needmtv = state->needmv = state->xupdated = false;
   state->userterminationneeded = false;
   state->running = true;
   state->repnmv = 0;
   state->repiterationscount = 0;
   state->r2 = state->bnorm2;
// estimate for ANorm
   for (normestimatorrestart(&state->nes); normestimatoriteration(&state->nes); ) {
      if (state->nes.needmv) {
         ae_v_move(state->x.xR, 1, state->nes.x.xR, 1, state->n);
         state->repnmv++, state->needmv = true, state->PQ = 0; goto Pause; Resume0: state->needmv = false;
         ae_v_move(state->nes.mv.xR, 1, state->mv.xR, 1, state->m);
      } else if (state->nes.needmtv) {
         ae_v_move(state->x.xR, 1, state->nes.x.xR, 1, state->m);
      // matrix-vector multiplication
         state->repnmv++, state->needmtv = true, state->PQ = 1; goto Pause; Resume1: state->needmtv = false;
         ae_v_move(state->nes.mtv.xR, 1, state->mtv.xR, 1, state->n);
      }
   }
   normestimatorresults(&state->nes, &state->anorm);
// initialize .RX by zeros
   for (i = 0; i < state->n; i++) {
      state->rx.xR[i] = 0.0;
   }
// output first report
   if (state->xrep) {
      ae_v_move(state->x.xR, 1, state->rx.xR, 1, state->n);
      state->xupdated = true, state->PQ = 2; goto Pause; Resume2: state->xupdated = false;
   }
// LSQR, Step 0.
//
// Algorithm outline corresponds to one which was described at p.50 of
// "LSQR - an algorithm for sparse linear equations and sparse least
// squares" by C.Paige and M.Saunders with one small addition - we
// explicitly extend system matrix by additional N lines in order
// to handle non-zero lambda, i.e. original A is replaced by
//         [ A        ]
// A_mod = [          ]
//         [ lambda*I ].
//
// Step 0:
//     x[0]          = 0
//     beta[1]*u[1]  = b
//     alpha[1]*v[1] = A_mod'*u[1]
//     w[1]          = v[1]
//     phiBar[1]     = beta[1]
//     rhoBar[1]     = alpha[1]
//     d[0]          = 0
//
// NOTE:
// There are three criteria for stopping:
// (S0) maximum number of iterations
// (S1) ||Rk|| <= EpsB*||B||;
// (S2) ||A^T*Rk||/(||A||*||Rk||) <= EpsA.
// It is very important that S2 always checked AFTER S1. It is necessary
// to avoid division by zero when Rk=0.
   state->betai = bnorm;
   if (state->betai == 0.0) {
   // Zero right part
      state->running = false;
      state->repterminationtype = 1;
      goto Exit;
   }
   for (i = 0; i < summn; i++) {
      if (i < state->m) {
         state->ui.xR[i] = state->b.xR[i] / state->betai;
      } else {
         state->ui.xR[i] = 0.0;
      }
      state->x.xR[i] = state->ui.xR[i];
   }
   state->repnmv++, state->needmtv = true, state->PQ = 3; goto Pause; Resume3: state->needmtv = false;
   for (i = 0; i < state->n; i++) {
      state->mtv.xR[i] += state->lambdai * state->ui.xR[state->m + i];
   }
   state->alphai = 0.0;
   for (i = 0; i < state->n; i++) {
      state->alphai += state->mtv.xR[i] * state->mtv.xR[i];
   }
   state->alphai = sqrt(state->alphai);
   if (state->alphai == 0.0) {
   // Orthogonality stopping criterion is met
      state->running = false;
      state->repterminationtype = 4;
      goto Exit;
   }
   for (i = 0; i < state->n; i++) {
      state->vi.xR[i] = state->mtv.xR[i] / state->alphai;
      state->omegai.xR[i] = state->vi.xR[i];
   }
   state->phibari = state->betai;
   state->rhobari = state->alphai;
   for (i = 0; i < state->n; i++) {
      state->d.xR[i] = 0.0;
   }
   state->dnorm = 0.0;
// Steps I=1, 2, ...
   while (true) {
   // At I-th step State.RepIterationsCount=I.
      state->repiterationscount++;
   // Bidiagonalization part:
   //     beta[i+1]*u[i+1]  = A_mod*v[i]-alpha[i]*u[i]
   //     alpha[i+1]*v[i+1] = A_mod'*u[i+1] - beta[i+1]*v[i]
   //
   // NOTE:  beta[i+1]=0 or alpha[i+1]=0 will lead to successful termination
   //        in the end of the current iteration. In this case u/v are zero.
   // NOTE2: algorithm won't fail on zero alpha or beta (there will be no
   //        division by zero because it will be stopped BEFORE division
   //        occurs). However, near-zero alpha and beta won't stop algorithm
   //        and, although no division by zero will happen, orthogonality
   //        in U and V will be lost.
      ae_v_move(state->x.xR, 1, state->vi.xR, 1, state->n);
      state->repnmv++, state->needmv = true, state->PQ = 4; goto Pause; Resume4: state->needmv = false;
      for (i = 0; i < state->n; i++) {
         state->mv.xR[state->m + i] = state->lambdai * state->vi.xR[i];
      }
      state->betaip1 = 0.0;
      for (i = 0; i < summn; i++) {
         state->uip1.xR[i] = state->mv.xR[i] - state->alphai * state->ui.xR[i];
         state->betaip1 += state->uip1.xR[i] * state->uip1.xR[i];
      }
      if (state->betaip1 != 0.0) {
         state->betaip1 = sqrt(state->betaip1);
         for (i = 0; i < summn; i++) {
            state->uip1.xR[i] /= state->betaip1;
         }
      }
      ae_v_move(state->x.xR, 1, state->uip1.xR, 1, state->m);
      state->repnmv++, state->needmtv = true, state->PQ = 5; goto Pause; Resume5: state->needmtv = false;
      for (i = 0; i < state->n; i++) {
         state->mtv.xR[i] += state->lambdai * state->uip1.xR[state->m + i];
      }
      state->alphaip1 = 0.0;
      for (i = 0; i < state->n; i++) {
         state->vip1.xR[i] = state->mtv.xR[i] - state->betaip1 * state->vi.xR[i];
         state->alphaip1 += state->vip1.xR[i] * state->vip1.xR[i];
      }
      if (state->alphaip1 != 0.0) {
         state->alphaip1 = sqrt(state->alphaip1);
         for (i = 0; i < state->n; i++) {
            state->vip1.xR[i] /= state->alphaip1;
         }
      }
   // Build next orthogonal transformation
      state->rhoi = safepythag2(state->rhobari, state->betaip1);
      state->ci = state->rhobari / state->rhoi;
      state->si = state->betaip1 / state->rhoi;
      state->theta = state->si * state->alphaip1;
      state->rhobarip1 = -state->ci * state->alphaip1;
      state->phii = state->ci * state->phibari;
      state->phibarip1 = state->si * state->phibari;
   // Update .RNorm
   //
   // This tricky  formula  is  necessary  because  simply  writing
   // State.R2:=State.PhiBarIP1*State.PhiBarIP1 does NOT guarantees
   // monotonic decrease of R2. Roundoff error combined with 80-bit
   // precision used internally by Intel chips allows R2 to increase
   // slightly in some rare, but possible cases. This property is
   // undesirable, so we prefer to guard against R increase.
      state->r2 = rmin2(state->r2, state->phibarip1 * state->phibarip1);
   // Update d and DNorm, check condition-related stopping criteria
      for (i = 0; i < state->n; i++) {
         state->d.xR[i] = 1 / state->rhoi * (state->vi.xR[i] - state->theta * state->d.xR[i]);
         state->dnorm += state->d.xR[i] * state->d.xR[i];
      }
      if (sqrt(state->dnorm) * state->anorm >= state->epsc) {
         state->running = false;
         state->repterminationtype = 7;
         goto Exit;
      }
   // Update x, output report
      for (i = 0; i < state->n; i++) {
         state->rx.xR[i] += state->phii / state->rhoi * state->omegai.xR[i];
      }
      if (state->xrep) {
         ae_v_move(state->x.xR, 1, state->rx.xR, 1, state->n);
         state->xupdated = true, state->PQ = 6; goto Pause; Resume6: state->xupdated = false;
      }
   // Check stopping criteria
   // 1. achieved required number of iterations;
   // 2. ||Rk|| <= EpsB*||B||;
   // 3. ||A^T*Rk||/(||A||*||Rk||) <= EpsA;
      if (state->maxits > 0 && state->repiterationscount >= state->maxits) {
      // Achieved required number of iterations
         state->running = false;
         state->repterminationtype = 5;
         goto Exit;
      }
      if (state->phibarip1 <= state->epsb * bnorm) {
      // ||Rk|| <= EpsB*||B||, here ||Rk||=PhiBar
         state->running = false;
         state->repterminationtype = 1;
         goto Exit;
      }
      if (state->alphaip1 * fabs(state->ci) / state->anorm <= state->epsa) {
      // ||A^T*Rk||/(||A||*||Rk||) <= EpsA, here ||A^T*Rk||=PhiBar*Alpha[i+1]*|.C|
         state->running = false;
         state->repterminationtype = 4;
         goto Exit;
      }
      if (state->userterminationneeded) {
      // User requested termination
         state->running = false;
         state->repterminationtype = 8;
         goto Exit;
      }
   // Update omega
      for (i = 0; i < state->n; i++) {
         state->omegaip1.xR[i] = state->vip1.xR[i] - state->theta / state->rhoi * state->omegai.xR[i];
      }
   // Prepare for the next iteration - rename variables:
   // u[i]   := u[i+1]
   // v[i]   := v[i+1]
   // rho[i] := rho[i+1]
   // ...
      ae_v_move(state->ui.xR, 1, state->uip1.xR, 1, summn);
      ae_v_move(state->vi.xR, 1, state->vip1.xR, 1, state->n);
      ae_v_move(state->omegai.xR, 1, state->omegaip1.xR, 1, state->n);
      state->alphai = state->alphaip1;
      state->betai = state->betaip1;
      state->phibari = state->phibarip1;
      state->rhobari = state->rhobarip1;
   }
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// Procedure for solution of A*x=b with sparse A.
//
// Inputs:
//     State   -   algorithm state
//     A       -   sparse M*N matrix in the CRS format (you MUST contvert  it
//                 to CRS format  by  calling  SparseConvertToCRS()  function
//                 BEFORE you pass it to this function).
//     B       -   right part, array[M]
//
// Result:
//     This function returns no result.
//     You can get solution by calling LinCGResults()
//
// NOTE: this function uses lightweight preconditioning -  multiplication  by
//       inverse of diag(A). If you want, you can turn preconditioning off by
//       calling LinLSQRSetPrecUnit(). However, preconditioning cost is   low
//       and preconditioner is very important for solution  of  badly  scaled
//       problems.
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
// API: void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b);
void linlsqrsolvesparse(linlsqrstate *state, sparsematrix *a, RVector *b) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t j;
   ae_int_t t0;
   ae_int_t t1;
   double v;
   n = state->n;
   ae_assert(!state->running, "LinLSQRSolveSparse: you can not call this function when LinLSQRIteration is running");
   ae_assert(b->cnt >= state->m, "LinLSQRSolveSparse: Length(B)<M");
   ae_assert(isfinitevector(b, state->m), "LinLSQRSolveSparse: B contains infinite or NaN values");
// Allocate temporaries
   vectorsetlengthatleast(&state->tmpd, n);
   vectorsetlengthatleast(&state->tmpx, n);
// Compute diagonal scaling matrix D
   if (state->prectype == 0) {
   // Default preconditioner - inverse of column norms
      for (i = 0; i < n; i++) {
         state->tmpd.xR[i] = 0.0;
      }
      t0 = 0;
      t1 = 0;
      while (sparseenumerate(a, &t0, &t1, &i, &j, &v)) {
         state->tmpd.xR[j] += sqr(v);
      }
      for (i = 0; i < n; i++) {
         if (state->tmpd.xR[i] > 0.0) {
            state->tmpd.xR[i] = 1 / sqrt(state->tmpd.xR[i]);
         } else {
            state->tmpd.xR[i] = 1.0;
         }
      }
   } else {
   // No diagonal scaling
      for (i = 0; i < n; i++) {
         state->tmpd.xR[i] = 1.0;
      }
   }
// Solve.
//
// Instead of solving A*x=b we solve preconditioned system (A*D)*(inv(D)*x)=b.
// Transformed A is not calculated explicitly, we just modify multiplication
// by A or A'. After solution we modify State.RX so it will store untransformed
// variables
   linlsqrsetb(state, b);
   linlsqrrestart(state);
   while (linlsqriteration(state)) {
      if (state->needmv) {
         for (i = 0; i < n; i++) {
            state->tmpx.xR[i] = state->tmpd.xR[i] * state->x.xR[i];
         }
         sparsemv(a, &state->tmpx, &state->mv);
      } else if (state->needmtv) {
         sparsemtv(a, &state->x, &state->mtv);
         for (i = 0; i < n; i++) {
            state->mtv.xR[i] *= state->tmpd.xR[i];
         }
      }
   }
   for (i = 0; i < n; i++) {
      state->rx.xR[i] *= state->tmpd.xR[i];
   }
}

// This function sets stopping criteria.
//
// Inputs:
//     EpsA    -   algorithm will be stopped if ||A^T*Rk||/(||A||*||Rk||) <= EpsA.
//     EpsB    -   algorithm will be stopped if ||Rk|| <= EpsB*||B||
//     MaxIts  -   algorithm will be stopped if number of iterations
//                 more than MaxIts.
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// NOTE: if EpsA,EpsB,EpsC and MaxIts are zero then these variables will
// be setted as default values.
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
// API: void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits);
void linlsqrsetcond(linlsqrstate *state, double epsa, double epsb, ae_int_t maxits) {
   ae_assert(!state->running, "LinLSQRSetCond: you can not call this function when LinLSQRIteration is running");
   ae_assert(isfinite(epsa) && epsa >= 0.0, "LinLSQRSetCond: EpsA is negative, INF or NAN");
   ae_assert(isfinite(epsb) && epsb >= 0.0, "LinLSQRSetCond: EpsB is negative, INF or NAN");
   ae_assert(maxits >= 0, "LinLSQRSetCond: MaxIts is negative");
   if (epsa == 0.0 && epsb == 0.0 && maxits == 0) {
      state->epsa = linlsqr_atol;
      state->epsb = linlsqr_btol;
      state->maxits = state->n;
   } else {
      state->epsa = epsa;
      state->epsb = epsb;
      state->maxits = maxits;
   }
}

// LSQR solver: results.
//
// This function must be called after LinLSQRSolve
//
// Inputs:
//     State   -   algorithm state
//
// Outputs:
//     X       -   array[N], solution
//     Rep     -   optimization report:
//                 * Rep.TerminationType completetion code:
//                     *  1    ||Rk|| <= EpsB*||B||
//                     *  4    ||A^T*Rk||/(||A||*||Rk||) <= EpsA
//                     *  5    MaxIts steps was taken
//                     *  7    rounding errors prevent further progress,
//                             X contains best point found so far.
//                             (sometimes returned on singular systems)
//                     *  8    user requested termination via calling
//                             linlsqrrequesttermination()
//                 * Rep.IterationsCount contains iterations count
//                 * NMV countains number of matrix-vector calculations
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
// API: void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep);
void linlsqrresults(linlsqrstate *state, RVector *x, linlsqrreport *rep) {
   SetVector(x);
   SetObj(linlsqrreport, rep);
   ae_assert(!state->running, "LinLSQRResult: you can not call this function when LinLSQRIteration is running");
   if (x->cnt < state->n) {
      ae_vector_set_length(x, state->n);
   }
   ae_v_move(x->xR, 1, state->rx.xR, 1, state->n);
   rep->iterationscount = state->repiterationscount;
   rep->nmv = state->repnmv;
   rep->terminationtype = state->repterminationtype;
}

// This function turns on/off reporting.
//
// Inputs:
//     State   -   structure which stores algorithm state
//     NeedXRep-   whether iteration reports are needed or not
//
// If NeedXRep is True, algorithm will call rep() callback function if  it is
// provided to MinCGOptimize().
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
// API: void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep);
void linlsqrsetxrep(linlsqrstate *state, bool needxrep) {
   state->xrep = needxrep;
}

// This function restarts LinLSQRIteration
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
void linlsqrrestart(linlsqrstate *state) {
   state->PQ = -1;
   state->repiterationscount = 0;
}

// This function is used to peek into LSQR solver and get  current  iteration
// counter. You can safely "peek" into the solver from another thread.
//
// Inputs:
//     S           -   solver object
//
// Result:
//     iteration counter, in [0,INF)
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
// API: ae_int_t linlsqrpeekiterationscount(const linlsqrstate &s);
ae_int_t linlsqrpeekiterationscount(linlsqrstate *s) {
   ae_int_t result;
   result = s->repiterationscount;
   return result;
}

// This subroutine submits request for termination of the running solver.  It
// can be called from some other thread which wants LSQR solver to  terminate
// (obviously, the  thread  running  LSQR  solver can not request termination
// because it is already busy working on LSQR).
//
// As result, solver  stops  at  point  which  was  "current  accepted"  when
// termination  request  was  submitted  and returns error code 8 (successful
// termination).  Such   termination   is  a smooth  process  which  properly
// deallocates all temporaries.
//
// Inputs:
//     State   -   solver structure
//
// NOTE: calling this function on solver which is NOT running  will  have  no
//       effect.
//
// NOTE: multiple calls to this function are possible. First call is counted,
//       subsequent calls are silently ignored.
//
// NOTE: solver clears termination flag on its start, it means that  if  some
//       other thread will request termination too soon, its request will went
//       unnoticed.
// ALGLIB: Copyright 08.10.2014 by Sergey Bochkanov
// API: void linlsqrrequesttermination(const linlsqrstate &state);
void linlsqrrequesttermination(linlsqrstate *state) {
   state->userterminationneeded = true;
}

void linlsqrstate_init(void *_p, bool make_automatic) {
   linlsqrstate *p = (linlsqrstate *)_p;
   normestimatorstate_init(&p->nes, make_automatic);
   ae_vector_init(&p->rx, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->b, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->ui, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->uip1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->vi, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->vip1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->omegai, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->omegaip1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->d, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->mv, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->mtv, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmpd, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmpx, 0, DT_REAL, make_automatic);
}

void linlsqrstate_copy(void *_dst, void *_src, bool make_automatic) {
   linlsqrstate *dst = (linlsqrstate *)_dst;
   linlsqrstate *src = (linlsqrstate *)_src;
   normestimatorstate_copy(&dst->nes, &src->nes, make_automatic);
   ae_vector_copy(&dst->rx, &src->rx, make_automatic);
   ae_vector_copy(&dst->b, &src->b, make_automatic);
   dst->n = src->n;
   dst->m = src->m;
   dst->prectype = src->prectype;
   ae_vector_copy(&dst->ui, &src->ui, make_automatic);
   ae_vector_copy(&dst->uip1, &src->uip1, make_automatic);
   ae_vector_copy(&dst->vi, &src->vi, make_automatic);
   ae_vector_copy(&dst->vip1, &src->vip1, make_automatic);
   ae_vector_copy(&dst->omegai, &src->omegai, make_automatic);
   ae_vector_copy(&dst->omegaip1, &src->omegaip1, make_automatic);
   dst->alphai = src->alphai;
   dst->alphaip1 = src->alphaip1;
   dst->betai = src->betai;
   dst->betaip1 = src->betaip1;
   dst->phibari = src->phibari;
   dst->phibarip1 = src->phibarip1;
   dst->phii = src->phii;
   dst->rhobari = src->rhobari;
   dst->rhobarip1 = src->rhobarip1;
   dst->rhoi = src->rhoi;
   dst->ci = src->ci;
   dst->si = src->si;
   dst->theta = src->theta;
   dst->lambdai = src->lambdai;
   ae_vector_copy(&dst->d, &src->d, make_automatic);
   dst->anorm = src->anorm;
   dst->bnorm2 = src->bnorm2;
   dst->dnorm = src->dnorm;
   dst->r2 = src->r2;
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   ae_vector_copy(&dst->mv, &src->mv, make_automatic);
   ae_vector_copy(&dst->mtv, &src->mtv, make_automatic);
   dst->epsa = src->epsa;
   dst->epsb = src->epsb;
   dst->epsc = src->epsc;
   dst->maxits = src->maxits;
   dst->xrep = src->xrep;
   dst->xupdated = src->xupdated;
   dst->needmv = src->needmv;
   dst->needmtv = src->needmtv;
   dst->repiterationscount = src->repiterationscount;
   dst->repnmv = src->repnmv;
   dst->repterminationtype = src->repterminationtype;
   dst->running = src->running;
   dst->userterminationneeded = src->userterminationneeded;
   ae_vector_copy(&dst->tmpd, &src->tmpd, make_automatic);
   ae_vector_copy(&dst->tmpx, &src->tmpx, make_automatic);
   dst->PQ = src->PQ;
}

void linlsqrstate_free(void *_p, bool make_automatic) {
   linlsqrstate *p = (linlsqrstate *)_p;
   normestimatorstate_free(&p->nes, make_automatic);
   ae_vector_free(&p->rx, make_automatic);
   ae_vector_free(&p->b, make_automatic);
   ae_vector_free(&p->ui, make_automatic);
   ae_vector_free(&p->uip1, make_automatic);
   ae_vector_free(&p->vi, make_automatic);
   ae_vector_free(&p->vip1, make_automatic);
   ae_vector_free(&p->omegai, make_automatic);
   ae_vector_free(&p->omegaip1, make_automatic);
   ae_vector_free(&p->d, make_automatic);
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->mv, make_automatic);
   ae_vector_free(&p->mtv, make_automatic);
   ae_vector_free(&p->tmpd, make_automatic);
   ae_vector_free(&p->tmpx, make_automatic);
}

void linlsqrreport_init(void *_p, bool make_automatic) {
}

void linlsqrreport_copy(void *_dst, void *_src, bool make_automatic) {
   linlsqrreport *dst = (linlsqrreport *)_dst;
   linlsqrreport *src = (linlsqrreport *)_src;
   dst->iterationscount = src->iterationscount;
   dst->nmv = src->nmv;
   dst->terminationtype = src->terminationtype;
}

void linlsqrreport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the LinLSQR method.
// You should use ALGLIB functions to work with this object.
DefClass(linlsqrstate, EndD)
DefClass(linlsqrreport, AndD DecVal(iterationscount) AndD DecVal(nmv) AndD DecVal(terminationtype))

void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrcreate(m, n, ConstT(linlsqrstate, state));
   alglib_impl::ae_state_clear();
}

void linlsqrcreatebuf(const ae_int_t m, const ae_int_t n, const linlsqrstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrcreatebuf(m, n, ConstT(linlsqrstate, state));
   alglib_impl::ae_state_clear();
}

void linlsqrsetprecunit(const linlsqrstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrsetprecunit(ConstT(linlsqrstate, state));
   alglib_impl::ae_state_clear();
}

void linlsqrsetprecdiag(const linlsqrstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrsetprecdiag(ConstT(linlsqrstate, state));
   alglib_impl::ae_state_clear();
}

void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrsetlambdai(ConstT(linlsqrstate, state), lambdai);
   alglib_impl::ae_state_clear();
}

void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrsolvesparse(ConstT(linlsqrstate, state), ConstT(sparsematrix, a), ConstT(ae_vector, b));
   alglib_impl::ae_state_clear();
}

void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrsetcond(ConstT(linlsqrstate, state), epsa, epsb, maxits);
   alglib_impl::ae_state_clear();
}

void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrresults(ConstT(linlsqrstate, state), ConstT(ae_vector, x), ConstT(linlsqrreport, rep));
   alglib_impl::ae_state_clear();
}

void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrsetxrep(ConstT(linlsqrstate, state), needxrep);
   alglib_impl::ae_state_clear();
}

ae_int_t linlsqrpeekiterationscount(const linlsqrstate &s) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::linlsqrpeekiterationscount(ConstT(linlsqrstate, s));
   alglib_impl::ae_state_clear();
   return Z;
}

void linlsqrrequesttermination(const linlsqrstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::linlsqrrequesttermination(ConstT(linlsqrstate, state));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === NLEQ Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: (LinAlg) FBLS
namespace alglib_impl {
// LEVENBERG-MARQUARDT-LIKE NONLINEAR SOLVER
// This algorithm solves a system of nonlinear equations
//     F[0](x[0], ..., x[N-1])   = 0
//     F[1](x[0], ..., x[N-1])   = 0
//     ...
//     F[M-1](x[0], ..., x[N-1]) = 0
// where M/N do not necessarily coincide. The algorithm converges quadratically
// under following conditions:
//     * the solution set XS is nonempty
//     * for some xs in XS there exist such neighbourhood N(xs) that:
//       * vector function F(x) and its Jacobian J(x) are continuously
//         differentiable on N
//       * ||F(x)|| provides local error bound on N, i.e. there  exists  such
//         c1, that ||F(x)|| > c1*distance(x,XS)
// Note that these conditions are much more weaker than usual non-singularity
// conditions. For example, algorithm will converge for any  affine  function
// F (whether its Jacobian singular or not).
//
// REQUIREMENTS:
// Algorithm will request following information during its operation:
// * function vector F[] and Jacobian matrix at given point X
// * value of merit function f(x)=F[0]^2(x)+...+F[M-1]^2(x) at given point X
//
// USAGE:
// 1. User initializes algorithm state with NLEQCreateLM() call
// 2. User tunes solver parameters with  NLEQSetCond(),  NLEQSetStpMax()  and
//    other functions
// 3. User  calls  NLEQSolve()  function  which  takes  algorithm  state  and
//    pointers (delegates, etc.) to callback functions which calculate  merit
//    function value and Jacobian.
// 4. User calls NLEQResults() to get solution
// 5. Optionally, user may call NLEQRestartFrom() to  solve  another  problem
//    with same parameters (N/M) but another starting  point  and/or  another
//    function vector. NLEQRestartFrom() allows to reuse already  initialized
//    structure.
//
// Inputs:
//     N       -   space dimension, N > 1:
//                 * if provided, only leading N elements of X are used
//                 * if not provided, determined automatically from size of X
//     M       -   system size
//     X       -   starting point
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// NOTES:
// 1. you may tune stopping conditions with NLEQSetCond() function
// 2. if target function contains exp() or other fast growing functions,  and
//    optimization algorithm makes too large steps which leads  to  overflow,
//    use NLEQSetStpMax() function to bound algorithm's steps.
// 3. this  algorithm  is  a  slightly  modified implementation of the method
//    described  in  'Levenberg-Marquardt  method  for constrained  nonlinear
//    equations with strong local convergence properties' by Christian Kanzow
//    Nobuo Yamashita and Masao Fukushima and further  developed  in  'On the
//    convergence of a New Levenberg-Marquardt Method'  by  Jin-yan  Fan  and
//    Ya-Xiang Yuan.
// ALGLIB: Copyright 20.08.2009 by Sergey Bochkanov
// API: void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state);
// API: void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state);
void nleqcreatelm(ae_int_t n, ae_int_t m, RVector *x, nleqstate *state) {
   SetObj(nleqstate, state);
   ae_assert(n >= 1, "NLEQCreateLM: N<1!");
   ae_assert(m >= 1, "NLEQCreateLM: M < 1!");
   ae_assert(x->cnt >= n, "NLEQCreateLM: Length(X)<N!");
   ae_assert(isfinitevector(x, n), "NLEQCreateLM: X contains infinite or NaN values!");
// Initialize
   state->n = n;
   state->m = m;
   nleqsetcond(state, 0.0, 0);
   nleqsetxrep(state, false);
   nleqsetstpmax(state, 0.0);
   ae_vector_set_length(&state->x, n);
   ae_vector_set_length(&state->xbase, n);
   ae_matrix_set_length(&state->j, m, n);
   ae_vector_set_length(&state->fi, m);
   ae_vector_set_length(&state->rightpart, n);
   ae_vector_set_length(&state->candstep, n);
   nleqrestartfrom(state, x);
}

// This function sets stopping conditions for the nonlinear solver
//
// Inputs:
//     State   -   structure which stores algorithm state
//     EpsF    - >= 0
//                 The subroutine finishes  its work if on k+1-th iteration
//                 the condition ||F|| <= EpsF is satisfied
//     MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
//                 iterations is unlimited.
//
// Passing EpsF=0 and MaxIts=0 simultaneously will lead to  automatic
// stopping criterion selection (small EpsF).
//
// NOTES:
// ALGLIB: Copyright 20.08.2010 by Sergey Bochkanov
// API: void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits);
void nleqsetcond(nleqstate *state, double epsf, ae_int_t maxits) {
   ae_assert(isfinite(epsf), "NLEQSetCond: EpsF is not finite number!");
   ae_assert(epsf >= 0.0, "NLEQSetCond: negative EpsF!");
   ae_assert(maxits >= 0, "NLEQSetCond: negative MaxIts!");
   if (epsf == 0.0 && maxits == 0) {
      epsf = 1.0E-6;
   }
   state->epsf = epsf;
   state->maxits = maxits;
}

// This function turns on/off reporting.
//
// Inputs:
//     State   -   structure which stores algorithm state
//     NeedXRep-   whether iteration reports are needed or not
//
// If NeedXRep is True, algorithm will call rep() callback function if  it is
// provided to NLEQSolve().
// ALGLIB: Copyright 20.08.2010 by Sergey Bochkanov
// API: void nleqsetxrep(const nleqstate &state, const bool needxrep);
void nleqsetxrep(nleqstate *state, bool needxrep) {
   state->xrep = needxrep;
}

// This function sets maximum step length
//
// Inputs:
//     State   -   structure which stores algorithm state
//     StpMax  -   maximum step length, >= 0. Set StpMax to 0.0,  if you don't
//                 want to limit step length.
//
// Use this subroutine when target function  contains  exp()  or  other  fast
// growing functions, and algorithm makes  too  large  steps  which  lead  to
// overflow. This function allows us to reject steps that are too large  (and
// therefore expose us to the possible overflow) without actually calculating
// function value at the x+stp*d.
// ALGLIB: Copyright 20.08.2010 by Sergey Bochkanov
// API: void nleqsetstpmax(const nleqstate &state, const double stpmax);
void nleqsetstpmax(nleqstate *state, double stpmax) {
   ae_assert(isfinite(stpmax), "NLEQSetStpMax: StpMax is not finite!");
   ae_assert(stpmax >= 0.0, "NLEQSetStpMax: StpMax<0!");
   state->stpmax = stpmax;
}

// Increases lambda, returns False when there is a danger of overflow
static bool nleq_increaselambda(double *lambdav, double *nu, double lambdaup) {
   double lnlambda;
   double lnnu;
   double lnlambdaup;
   double lnmax;
   bool result;
   result = false;
   lnlambda = log(*lambdav);
   lnlambdaup = log(lambdaup);
   lnnu = log(*nu);
   lnmax = 0.5 * log(maxrealnumber);
   if (lnlambda + lnlambdaup + lnnu > lnmax) {
      return result;
   }
   if (lnnu + log(2.0) > lnmax) {
      return result;
   }
   *lambdav *= lambdaup * (*nu);
   *nu *= 2;
   result = true;
   return result;
}

// Decreases lambda, but leaves it unchanged when there is danger of underflow.
static void nleq_decreaselambda(double *lambdav, double *nu, double lambdadown) {
   *nu = 1.0;
   if (log(*lambdav) + log(lambdadown) < log(minrealnumber)) {
      *lambdav = minrealnumber;
   } else {
      *lambdav *= lambdadown;
   }
}

// This function provides a reverse communication interface, which is not documented or recommended for use.
// Instead, it is recommended that you use the better-documented API function nleqsolve() listed below.
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
// API: bool nleqiteration(const nleqstate &state);
// API: void nleqsolve(nleqstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL);
bool nleqiteration(nleqstate *state) {
   AutoS ae_int_t n;
   AutoS ae_int_t m;
   AutoS ae_int_t i;
   AutoS double lambdaup;
   AutoS double lambdadown;
   AutoS double lambdav;
   AutoS double rho;
   AutoS double mu;
   AutoS double stepnorm;
   AutoS bool b;
// Manually threaded two-way signalling.
// Locals are set arbitrarily the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0; case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3; case 4: goto Resume4;
      default: goto Exit;
   }
Spawn:
   i = -919;
   b = true;
   lambdaup = 81;
   lambdadown = 255;
   lambdav = 74;
   rho = -788;
   mu = 809;
   stepnorm = 205;
// Prepare
   n = state->n;
   m = state->m;
   state->xupdated = state->needfij = state->needf = false;
   state->repterminationtype = 0;
   state->repiterationscount = 0;
   state->repnfunc = 0;
   state->repnjac = 0;
// Calculate F/G, initialize algorithm
   state->needf = true, state->PQ = 0; goto Pause; Resume0: state->needf = false, state->repnfunc++;
   ae_v_move(state->xbase.xR, 1, state->x.xR, 1, n);
   state->fbase = state->f;
   state->fprev = maxrealnumber;
   if (state->xrep) {
   // progress report
      state->xupdated = true, state->PQ = 1; goto Pause; Resume1: state->xupdated = false;
   }
   if (state->f <= sqr(state->epsf)) {
      state->repterminationtype = 1;
      goto Exit;
   }
// Main cycle
   lambdaup = 10.0;
   lambdadown = 0.3;
   lambdav = 0.001;
   rho = 1.0;
   do {
   // Get Jacobian;
   // before we get to this point we already have State.XBase filled
   // with current point and State.FBase filled with function value
   // at XBase
      ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
      state->needfij = true, state->PQ = 2; goto Pause; Resume2: state->needfij = false, state->repnfunc++, state->repnjac++;
      rmatrixmv(n, m, &state->j, 0, 0, 1, &state->fi, 0, &state->rightpart, 0);
      ae_v_muld(state->rightpart.xR, 1, n, -1);
   // Inner cycle: find good lambda
      while (true) {
      // Solve (J^T*J + (Lambda+Mu)*I)*y = J^T*F
      // to get step d=-y where:
      // * Mu=||F|| - is damping parameter for nonlinear system
      // * Lambda   - is additional Levenberg-Marquardt parameter
      //              for better convergence when far away from minimum
         for (i = 0; i < n; i++) {
            state->candstep.xR[i] = 0.0;
         }
         fblssolvecgx(&state->j, m, n, lambdav, &state->rightpart, &state->candstep, &state->cgbuf);
      // Normalize step (it must be no more than StpMax)
         stepnorm = 0.0;
         for (i = 0; i < n; i++) {
            if (state->candstep.xR[i] != 0.0) {
               stepnorm = 1.0;
               break;
            }
         }
         linminnormalized(&state->candstep, &stepnorm, n);
         if (state->stpmax != 0.0) {
            stepnorm = rmin2(stepnorm, state->stpmax);
         }
      // Test new step - is it good enough?
      // * if not, Lambda is increased and we try again.
      // * if step is good, we decrease Lambda and move on.
      //
      // We can break this cycle on two occasions:
      // * step is so small that x + step == x (in floating point arithmetics)
      // * lambda is so large
         ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
         ae_v_addd(state->x.xR, 1, state->candstep.xR, 1, n, stepnorm);
         b = true;
         for (i = 0; i < n; i++) {
            if (state->x.xR[i] != state->xbase.xR[i]) {
               b = false;
               break;
            }
         }
         if (b) {
         // Step is too small, force zero step and break
            stepnorm = 0.0;
            ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
            state->f = state->fbase;
            break;
         }
         state->needf = true, state->PQ = 3; goto Pause; Resume3: state->needf = false, state->repnfunc++;
         if (state->f < state->fbase) {
         // function value decreased, move on
            nleq_decreaselambda(&lambdav, &rho, lambdadown);
            break;
         }
         if (!nleq_increaselambda(&lambdav, &rho, lambdaup)) {
         // Lambda is too large (near overflow), force zero step and break
            stepnorm = 0.0;
            ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
            state->f = state->fbase;
            break;
         }
      }
   // Accept step:
   // * new position
   // * new function value
      state->fbase = state->f;
      ae_v_addd(state->xbase.xR, 1, state->candstep.xR, 1, n, stepnorm);
      state->repiterationscount++;
   // Report new iteration
      if (state->xrep) {
         state->f = state->fbase;
         ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
         state->xupdated = true, state->PQ = 4; goto Pause; Resume4: state->xupdated = false;
      }
   // Test stopping conditions on F, step (zero/non-zero) and MaxIts;
   // If one of the conditions is met, RepTerminationType is changed.
      if (sqrt(state->f) <= state->epsf) {
         state->repterminationtype = 1;
      }
      if (stepnorm == 0.0 && state->repterminationtype == 0) {
         state->repterminationtype = -4;
      }
      if (state->repiterationscount >= state->maxits && state->maxits > 0) {
         state->repterminationtype = 5;
      }
   // Now, iteration is finally over
   } while (state->repterminationtype == 0);
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// NLEQ solver results
//
// Inputs:
//     State   -   algorithm state.
//
// Outputs:
//     X       -   array[0..N-1], solution
//     Rep     -   optimization report:
//                 * Rep.TerminationType completetion code:
//                     * -4    ERROR:  algorithm   has   converged   to   the
//                             stationary point Xf which is local minimum  of
//                             f=F[0]^2+...+F[m-1]^2, but is not solution  of
//                             nonlinear system.
//                     *  1    sqrt(f) <= EpsF.
//                     *  5    MaxIts steps was taken
//                     *  7    stopping conditions are too stringent,
//                             further improvement is impossible
//                 * Rep.IterationsCount contains iterations count
//                 * NFEV countains number of function calculations
//                 * ActiveConstraints contains number of active constraints
// ALGLIB: Copyright 20.08.2009 by Sergey Bochkanov
// API: void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep);
void nleqresults(nleqstate *state, RVector *x, nleqreport *rep) {
   SetVector(x);
   SetObj(nleqreport, rep);
   nleqresultsbuf(state, x, rep);
}

// NLEQ solver results
//
// Buffered implementation of NLEQResults(), which uses pre-allocated  buffer
// to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
// intended to be used in the inner cycles of performance critical algorithms
// where array reallocation penalty is too large to be ignored.
// ALGLIB: Copyright 20.08.2009 by Sergey Bochkanov
// API: void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep);
void nleqresultsbuf(nleqstate *state, RVector *x, nleqreport *rep) {
   if (x->cnt < state->n) {
      ae_vector_set_length(x, state->n);
   }
   ae_v_move(x->xR, 1, state->xbase.xR, 1, state->n);
   rep->iterationscount = state->repiterationscount;
   rep->nfunc = state->repnfunc;
   rep->njac = state->repnjac;
   rep->terminationtype = state->repterminationtype;
}

// This  subroutine  restarts  CG  algorithm from new point. All optimization
// parameters are left unchanged.
//
// This  function  allows  to  solve multiple  optimization  problems  (which
// must have same number of dimensions) without object reallocation penalty.
//
// Inputs:
//     State   -   structure used for reverse communication previously
//                 allocated with MinCGCreate call.
//     X       -   new starting point.
//     BndL    -   new lower bounds
//     BndU    -   new upper bounds
// ALGLIB: Copyright 30.07.2010 by Sergey Bochkanov
// API: void nleqrestartfrom(const nleqstate &state, const real_1d_array &x);
void nleqrestartfrom(nleqstate *state, RVector *x) {
   ae_assert(x->cnt >= state->n, "NLEQRestartFrom: Length(X)<N!");
   ae_assert(isfinitevector(x, state->n), "NLEQRestartFrom: X contains infinite or NaN values!");
   ae_v_move(state->x.xR, 1, x->xR, 1, state->n);
   state->PQ = -1;
}

void nleqstate_init(void *_p, bool make_automatic) {
   nleqstate *p = (nleqstate *)_p;
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->fi, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->j, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xbase, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->candstep, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rightpart, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->cgbuf, 0, DT_REAL, make_automatic);
}

void nleqstate_copy(void *_dst, void *_src, bool make_automatic) {
   nleqstate *dst = (nleqstate *)_dst;
   nleqstate *src = (nleqstate *)_src;
   dst->n = src->n;
   dst->m = src->m;
   dst->epsf = src->epsf;
   dst->maxits = src->maxits;
   dst->xrep = src->xrep;
   dst->stpmax = src->stpmax;
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   dst->f = src->f;
   ae_vector_copy(&dst->fi, &src->fi, make_automatic);
   ae_matrix_copy(&dst->j, &src->j, make_automatic);
   dst->needf = src->needf;
   dst->needfij = src->needfij;
   dst->xupdated = src->xupdated;
   dst->PQ = src->PQ;
   dst->repiterationscount = src->repiterationscount;
   dst->repnfunc = src->repnfunc;
   dst->repnjac = src->repnjac;
   dst->repterminationtype = src->repterminationtype;
   ae_vector_copy(&dst->xbase, &src->xbase, make_automatic);
   dst->fbase = src->fbase;
   dst->fprev = src->fprev;
   ae_vector_copy(&dst->candstep, &src->candstep, make_automatic);
   ae_vector_copy(&dst->rightpart, &src->rightpart, make_automatic);
   ae_vector_copy(&dst->cgbuf, &src->cgbuf, make_automatic);
}

void nleqstate_free(void *_p, bool make_automatic) {
   nleqstate *p = (nleqstate *)_p;
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->fi, make_automatic);
   ae_matrix_free(&p->j, make_automatic);
   ae_vector_free(&p->xbase, make_automatic);
   ae_vector_free(&p->candstep, make_automatic);
   ae_vector_free(&p->rightpart, make_automatic);
   ae_vector_free(&p->cgbuf, make_automatic);
}

void nleqreport_init(void *_p, bool make_automatic) {
}

void nleqreport_copy(void *_dst, void *_src, bool make_automatic) {
   nleqreport *dst = (nleqreport *)_dst;
   nleqreport *src = (nleqreport *)_src;
   dst->iterationscount = src->iterationscount;
   dst->nfunc = src->nfunc;
   dst->njac = src->njac;
   dst->terminationtype = src->terminationtype;
}

void nleqreport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
DefClass(nleqstate, AndD DecVal(needf) AndD DecVal(needfij) AndD DecVal(xupdated) AndD DecVal(f) AndD DecVar(fi) AndD DecVar(j) AndD DecVar(x))
DefClass(nleqreport, AndD DecVal(iterationscount) AndD DecVal(nfunc) AndD DecVal(njac) AndD DecVal(terminationtype))

void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqcreatelm(n, m, ConstT(ae_vector, x), ConstT(nleqstate, state));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state) {
   ae_int_t n = x.length();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqcreatelm(n, m, ConstT(ae_vector, x), ConstT(nleqstate, state));
   alglib_impl::ae_state_clear();
}
#endif

void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqsetcond(ConstT(nleqstate, state), epsf, maxits);
   alglib_impl::ae_state_clear();
}

void nleqsetxrep(const nleqstate &state, const bool needxrep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqsetxrep(ConstT(nleqstate, state), needxrep);
   alglib_impl::ae_state_clear();
}

void nleqsetstpmax(const nleqstate &state, const double stpmax) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqsetstpmax(ConstT(nleqstate, state), stpmax);
   alglib_impl::ae_state_clear();
}

bool nleqiteration(const nleqstate &state) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::nleqiteration(ConstT(nleqstate, state));
   alglib_impl::ae_state_clear();
   return Ok;
}

// This family of functions is used to launch iterations of nonlinear solver
//
// These functions accept following parameters:
//     state   -   algorithm state
//     func    -   callback which calculates function (or merit function)
//                 value func at given point x
//     jac     -   callback which calculates function vector fi[]
//                 and Jacobian jac at given point x
//     rep     -   optional callback which is called after each iteration
//                 can be NULL
//     ptr     -   optional pointer which is passed to func/grad/hess/jac/rep
//                 can be NULL
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void nleqsolve(nleqstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr)/* = NULL*/, void *ptr/* = NULL*/) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ae_assert(func != NULL, "nleqsolve: func is NULL");
   alglib_impl::ae_assert(jac != NULL, "nleqsolve: jac is NULL");
   while (alglib_impl::nleqiteration(state.c_ptr()))
   BegPoll
      if (state.needf) func(state.x, state.f, ptr);
      else if (state.needfij) jac(state.x, state.fi, state.j, ptr);
      else if (state.xupdated) { if (rep != NULL) rep(state.x, state.f, ptr); }
      else alglib_impl::ae_assert(false, "nleqsolve: some derivatives were not provided?");
   EndPoll
   alglib_impl::ae_state_clear();
}

void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqresults(ConstT(nleqstate, state), ConstT(ae_vector, x), ConstT(nleqreport, rep));
   alglib_impl::ae_state_clear();
}

void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqresultsbuf(ConstT(nleqstate, state), ConstT(ae_vector, x), ConstT(nleqreport, rep));
   alglib_impl::ae_state_clear();
}

void nleqrestartfrom(const nleqstate &state, const real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::nleqrestartfrom(ConstT(nleqstate, state), ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib
