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
#ifndef OnceOnlySolvers_h
#define OnceOnlySolvers_h

#include "LinAlg.h"

// === DIRECTDENSESOLVERS Package ===
// Depends on: (AlgLibInternal) XBLAS
// Depends on: (LinAlg) RCOND, SVD
namespace alglib_impl {
typedef struct {
   double r1;
   double rinf;
} densesolverreport;
void densesolverreport_init(void *_p, bool make_automatic);
void densesolverreport_copy(void *_dst, void *_src, bool make_automatic);
void densesolverreport_free(void *_p, bool make_automatic);

typedef struct {
   double r2;
   ae_matrix cx;
   ae_int_t n;
   ae_int_t k;
} densesolverlsreport;
void densesolverlsreport_init(void *_p, bool make_automatic);
void densesolverlsreport_copy(void *_dst, void *_src, bool make_automatic);
void densesolverlsreport_free(void *_p, bool make_automatic);

void rmatrixsolve(RMatrix a, ae_int_t n, RVector b, ae_int_t *info, densesolverreport *rep, RVector x);
void rmatrixsolvefast(RMatrix a, ae_int_t n, RVector b, ae_int_t *info);
void rmatrixsolvem(RMatrix a, ae_int_t n, RMatrix b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, RMatrix x);
void rmatrixsolvemfast(RMatrix a, ae_int_t n, RMatrix b, ae_int_t m, ae_int_t *info);
void rmatrixlusolve(RMatrix lua, ZVector p, ae_int_t n, RVector b, ae_int_t *info, densesolverreport *rep, RVector x);
void rmatrixlusolvefast(RMatrix lua, ZVector p, ae_int_t n, RVector b, ae_int_t *info);
void rmatrixlusolvem(RMatrix lua, ZVector p, ae_int_t n, RMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix x);
void rmatrixlusolvemfast(RMatrix lua, ZVector p, ae_int_t n, RMatrix b, ae_int_t m, ae_int_t *info);
void rmatrixmixedsolve(RMatrix a, RMatrix lua, ZVector p, ae_int_t n, RVector b, ae_int_t *info, densesolverreport *rep, RVector x);
void rmatrixmixedsolvem(RMatrix a, RMatrix lua, ZVector p, ae_int_t n, RMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix x);
void cmatrixsolvem(CMatrix a, ae_int_t n, CMatrix b, ae_int_t m, bool rfs, ae_int_t *info, densesolverreport *rep, CMatrix x);
void cmatrixsolvemfast(CMatrix a, ae_int_t n, CMatrix b, ae_int_t m, ae_int_t *info);
void cmatrixsolve(CMatrix a, ae_int_t n, CVector b, ae_int_t *info, densesolverreport *rep, CVector x);
void cmatrixsolvefast(CMatrix a, ae_int_t n, CVector b, ae_int_t *info);
void cmatrixlusolvem(CMatrix lua, ZVector p, ae_int_t n, CMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix x);
void cmatrixlusolvemfast(CMatrix lua, ZVector p, ae_int_t n, CMatrix b, ae_int_t m, ae_int_t *info);
void cmatrixlusolve(CMatrix lua, ZVector p, ae_int_t n, CVector b, ae_int_t *info, densesolverreport *rep, CVector x);
void cmatrixlusolvefast(CMatrix lua, ZVector p, ae_int_t n, CVector b, ae_int_t *info);
void cmatrixmixedsolvem(CMatrix a, CMatrix lua, ZVector p, ae_int_t n, CMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix x);
void cmatrixmixedsolve(CMatrix a, CMatrix lua, ZVector p, ae_int_t n, CVector b, ae_int_t *info, densesolverreport *rep, CVector x);
void spdmatrixsolvem(RMatrix a, ae_int_t n, bool isupper, RMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix x);
void spdmatrixsolvemfast(RMatrix a, ae_int_t n, bool isupper, RMatrix b, ae_int_t m, ae_int_t *info);
void spdmatrixsolve(RMatrix a, ae_int_t n, bool isupper, RVector b, ae_int_t *info, densesolverreport *rep, RVector x);
void spdmatrixsolvefast(RMatrix a, ae_int_t n, bool isupper, RVector b, ae_int_t *info);
void spdmatrixcholeskysolvem(RMatrix cha, ae_int_t n, bool isupper, RMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, RMatrix x);
void spdmatrixcholeskysolvemfast(RMatrix cha, ae_int_t n, bool isupper, RMatrix b, ae_int_t m, ae_int_t *info);
void spdmatrixcholeskysolve(RMatrix cha, ae_int_t n, bool isupper, RVector b, ae_int_t *info, densesolverreport *rep, RVector x);
void spdmatrixcholeskysolvefast(RMatrix cha, ae_int_t n, bool isupper, RVector b, ae_int_t *info);
void hpdmatrixsolvem(CMatrix a, ae_int_t n, bool isupper, CMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix x);
void hpdmatrixsolvemfast(CMatrix a, ae_int_t n, bool isupper, CMatrix b, ae_int_t m, ae_int_t *info);
void hpdmatrixsolve(CMatrix a, ae_int_t n, bool isupper, CVector b, ae_int_t *info, densesolverreport *rep, CVector x);
void hpdmatrixsolvefast(CMatrix a, ae_int_t n, bool isupper, CVector b, ae_int_t *info);
void hpdmatrixcholeskysolvem(CMatrix cha, ae_int_t n, bool isupper, CMatrix b, ae_int_t m, ae_int_t *info, densesolverreport *rep, CMatrix x);
void hpdmatrixcholeskysolvemfast(CMatrix cha, ae_int_t n, bool isupper, CMatrix b, ae_int_t m, ae_int_t *info);
void hpdmatrixcholeskysolve(CMatrix cha, ae_int_t n, bool isupper, CVector b, ae_int_t *info, densesolverreport *rep, CVector x);
void hpdmatrixcholeskysolvefast(CMatrix cha, ae_int_t n, bool isupper, CVector b, ae_int_t *info);
void rmatrixsolvels(RMatrix a, ae_int_t nrows, ae_int_t ncols, RVector b, double threshold, ae_int_t *info, densesolverlsreport *rep, RVector x);
} // end of namespace alglib_impl

namespace alglib {
DecClass(densesolverreport, double &r1; double &rinf;);
DecClass(densesolverlsreport, double &r2; real_2d_array cx; ae_int_t &n; ae_int_t &k;);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 16.03.2015 by Sergey Bochkanov
void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
//
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
//
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
//
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
//
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 16.03.2015 by Sergey Bochkanov
void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
//
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);

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
//                 * info > 0   =>  overwritten by solution
//                 * info=-3   =>  filled by zeros
//
// NOTE: unlike  CMatrixLUSolve(),  this   function   does   NOT   check  for
//       near-degeneracy of input matrix. It  checks  for  EXACT  degeneracy,
//       because this check is easy to do. However,  very  badly  conditioned
//       matrices may went unnoticed.
//
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 17.03.2015 by Sergey Bochkanov
void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);

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
//                 * info > 0   =>  solution
//                 * info=-3   =>  filled by zeros
// ALGLIB: Copyright 17.03.2015 by Sergey Bochkanov
void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);

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
//                 * for info > 0contains solution
//                 * for info=-3 filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);

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
//                 * for info > 0overwritten by solution
//                 * for info=-3 filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * for info > 0 - solution
//                 * for info=-3 - filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);

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
//                 * for info > 0 - overwritten by solution
//                 * for info=-3 - filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);

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
void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);

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
void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);

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
void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);

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
void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);

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
//                 * for info > 0contains solution
//                 * for info=-3 filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);

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
//                 * for info > 0overwritten by solution
//                 * for info=-3 filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);

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
//                 * for info > 0 - solution
//                 * for info=-3 - filled by zeros
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);

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
//                 * for info > 0 - overwritten by solution
//                 * for info=-3 - filled by zeros
// ALGLIB: Copyright 18.03.2015 by Sergey Bochkanov
void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);

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
void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);
} // end of namespace alglib

// === LINLSQR Package ===
// Depends on: (LinAlg) SVD, NORMESTIMATOR
namespace alglib_impl {
typedef struct {
   normestimatorstate nes;
   ae_vector rx;
   ae_vector b;
   ae_int_t n;
   ae_int_t m;
   ae_int_t prectype;
   ae_vector ui;
   ae_vector uip1;
   ae_vector vi;
   ae_vector vip1;
   ae_vector omegai;
   ae_vector omegaip1;
   double alphai;
   double alphaip1;
   double betai;
   double betaip1;
   double phibari;
   double phibarip1;
   double phii;
   double rhobari;
   double rhobarip1;
   double rhoi;
   double ci;
   double si;
   double theta;
   double lambdai;
   ae_vector d;
   double anorm;
   double bnorm2;
   double dnorm;
   double r2;
   ae_vector x;
   ae_vector mv;
   ae_vector mtv;
   double epsa;
   double epsb;
   double epsc;
   ae_int_t maxits;
   bool xrep;
   bool xupdated;
   bool needmv;
   bool needmtv;
// bool needmv2; //(@) Not used.
// bool needvmv; //(@) Not used.
// bool needprec; //(@) Not used.
   ae_int_t repiterationscount;
   ae_int_t repnmv;
   ae_int_t repterminationtype;
   bool running;
   bool userterminationneeded;
   ae_vector tmpd;
   ae_vector tmpx;
   ae_int_t PQ;
} linlsqrstate;
void linlsqrstate_init(void *_p, bool make_automatic);
void linlsqrstate_copy(void *_dst, void *_src, bool make_automatic);
void linlsqrstate_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
} linlsqrreport;
void linlsqrreport_init(void *_p, bool make_automatic);
void linlsqrreport_copy(void *_dst, void *_src, bool make_automatic);
void linlsqrreport_free(void *_p, bool make_automatic);

void linlsqrcreate(ae_int_t m, ae_int_t n, linlsqrstate *state);
void linlsqrcreatebuf(ae_int_t m, ae_int_t n, linlsqrstate *state);
void linlsqrsetb(linlsqrstate *state, RVector b);
void linlsqrsetprecunit(linlsqrstate *state);
void linlsqrsetprecdiag(linlsqrstate *state);
void linlsqrsetlambdai(linlsqrstate *state, double lambdai);
bool linlsqriteration(linlsqrstate *state);
void linlsqrsolvesparse(linlsqrstate *state, sparsematrix *a, RVector b);
void linlsqrsetcond(linlsqrstate *state, double epsa, double epsb, ae_int_t maxits);
void linlsqrresults(linlsqrstate *state, RVector x, linlsqrreport *rep);
void linlsqrsetxrep(linlsqrstate *state, bool needxrep);
void linlsqrrestart(linlsqrstate *state);
ae_int_t linlsqrpeekiterationscount(linlsqrstate *s);
void linlsqrrequesttermination(linlsqrstate *state);
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the LinLSQR method.
//
// You should use ALGLIB functions to work with this object.
DecClass(linlsqrstate, EndD);
DecClass(linlsqrreport, ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype;);

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
void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state);

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
void linlsqrcreatebuf(const ae_int_t m, const ae_int_t n, const linlsqrstate &state);

// This  function  changes  preconditioning  settings of LinLSQQSolveSparse()
// function. By default, SolveSparse() uses diagonal preconditioner,  but  if
// you want to use solver without preconditioning, you can call this function
// which forces solver to use unit matrix for preconditioning.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
void linlsqrsetprecunit(const linlsqrstate &state);

// This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
// function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
// preconditioner. This preconditioning mode is active by default.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
void linlsqrsetprecdiag(const linlsqrstate &state);

// This function sets optional Tikhonov regularization coefficient.
// It is zero by default.
//
// Inputs:
//     LambdaI -   regularization factor, LambdaI >= 0
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai);

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
void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b);

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
void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits);

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
void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep);

// This function turns on/off reporting.
//
// Inputs:
//     State   -   structure which stores algorithm state
//     NeedXRep-   whether iteration reports are needed or not
//
// If NeedXRep is True, algorithm will call rep() callback function if  it is
// provided to MinCGOptimize().
// ALGLIB: Copyright 30.11.2011 by Sergey Bochkanov
void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep);

// This function is used to peek into LSQR solver and get  current  iteration
// counter. You can safely "peek" into the solver from another thread.
//
// Inputs:
//     S           -   solver object
//
// Result:
//     iteration counter, in [0,INF)
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
ae_int_t linlsqrpeekiterationscount(const linlsqrstate &s);

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
void linlsqrrequesttermination(const linlsqrstate &state);
} // end of namespace alglib

// === POLYNOMIALSOLVER Package ===
// Depends on: (LinAlg) TRFAC, EVD
namespace alglib_impl {
typedef struct {
   double maxerr;
} polynomialsolverreport;
void polynomialsolverreport_init(void *_p, bool make_automatic);
void polynomialsolverreport_copy(void *_dst, void *_src, bool make_automatic);
void polynomialsolverreport_free(void *_p, bool make_automatic);

void polynomialsolve(RVector a, ae_int_t n, CVector x, polynomialsolverreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(polynomialsolverreport, double &maxerr;);

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
void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep);
} // end of namespace alglib

// === NLEQ Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: (LinAlg) FBLS
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t m;
   double epsf;
   ae_int_t maxits;
   bool xrep;
   double stpmax;
   ae_vector x;
   double f;
   ae_vector fi;
   ae_matrix j;
   bool needf;
   bool needfij;
   bool xupdated;
   ae_int_t PQ;
   ae_int_t repiterationscount;
   ae_int_t repnfunc;
   ae_int_t repnjac;
   ae_int_t repterminationtype;
   ae_vector xbase;
   double fbase;
   double fprev;
   ae_vector candstep;
   ae_vector rightpart;
   ae_vector cgbuf;
} nleqstate;
void nleqstate_init(void *_p, bool make_automatic);
void nleqstate_copy(void *_dst, void *_src, bool make_automatic);
void nleqstate_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfunc;
   ae_int_t njac;
   ae_int_t terminationtype;
} nleqreport;
void nleqreport_init(void *_p, bool make_automatic);
void nleqreport_copy(void *_dst, void *_src, bool make_automatic);
void nleqreport_free(void *_p, bool make_automatic);

void nleqcreatelm(ae_int_t n, ae_int_t m, RVector x, nleqstate *state);
void nleqsetcond(nleqstate *state, double epsf, ae_int_t maxits);
void nleqsetxrep(nleqstate *state, bool needxrep);
void nleqsetstpmax(nleqstate *state, double stpmax);
bool nleqiteration(nleqstate *state);
void nleqresults(nleqstate *state, RVector x, nleqreport *rep);
void nleqresultsbuf(nleqstate *state, RVector x, nleqreport *rep);
void nleqrestartfrom(nleqstate *state, RVector x);
} // end of namespace alglib_impl

namespace alglib {
DecClass(nleqstate, bool &needf; bool &needfij; bool &xupdated; double &f; real_1d_array fi; real_2d_array j; real_1d_array x;);
DecClass(nleqreport, ae_int_t &iterationscount; ae_int_t &nfunc; ae_int_t &njac; ae_int_t &terminationtype;);

//                 LEVENBERG-MARQUARDT-LIKE NONLINEAR SOLVER
//
// DESCRIPTION:
// This algorithm solves system of nonlinear equations
//     F[0](x[0], ..., x[n-1])   = 0
//     F[1](x[0], ..., x[n-1])   = 0
//     ...
//     F[M-1](x[0], ..., x[n-1]) = 0
// with M/N do not necessarily coincide.  Algorithm  converges  quadratically
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
//
// REQUIREMENTS:
// Algorithm will request following information during its operation:
// * function vector F[] and Jacobian matrix at given point X
// * value of merit function f(x)=F[0]^2(x)+...+F[M-1]^2(x) at given point X
//
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
//
// Inputs:
//     N       -   space dimension, N > 1:
//                 * if provided, only leading N elements of X are used
//                 * if not provided, determined automatically from size of X
//     M       -   system size
//     X       -   starting point
//
//
// Outputs:
//     State   -   structure which stores algorithm state
//
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
//
// ALGLIB: Copyright 20.08.2009 by Sergey Bochkanov
void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state);
void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state);

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
void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits);

// This function turns on/off reporting.
//
// Inputs:
//     State   -   structure which stores algorithm state
//     NeedXRep-   whether iteration reports are needed or not
//
// If NeedXRep is True, algorithm will call rep() callback function if  it is
// provided to NLEQSolve().
// ALGLIB: Copyright 20.08.2010 by Sergey Bochkanov
void nleqsetxrep(const nleqstate &state, const bool needxrep);

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
void nleqsetstpmax(const nleqstate &state, const double stpmax);

// This function provides reverse communication interface
// Reverse communication interface is not documented or recommended for use.
// See below for functions which provide better documented API
bool nleqiteration(const nleqstate &state);

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
//
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void nleqsolve(nleqstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL);

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
void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep);

// NLEQ solver results
//
// Buffered implementation of NLEQResults(), which uses pre-allocated  buffer
// to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
// intended to be used in the inner cycles of performance critical algorithms
// where array reallocation penalty is too large to be ignored.
// ALGLIB: Copyright 20.08.2009 by Sergey Bochkanov
void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep);

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
void nleqrestartfrom(const nleqstate &state, const real_1d_array &x);
} // end of namespace alglib

// === DIRECTSPARSESOLVERS Package ===
// Depends on: (LinAlg) TRFAC
namespace alglib_impl {
typedef struct {
   ae_int_t terminationtype;
} sparsesolverreport;
void sparsesolverreport_init(void *_p, bool make_automatic);
void sparsesolverreport_copy(void *_dst, void *_src, bool make_automatic);
void sparsesolverreport_free(void *_p, bool make_automatic);

void sparsesolvesks(sparsematrix *a, ae_int_t n, bool isupper, RVector b, sparsesolverreport *rep, RVector x);
void sparsecholeskysolvesks(sparsematrix *a, ae_int_t n, bool isupper, RVector b, sparsesolverreport *rep, RVector x);
void sparsesolve(sparsematrix *a, ae_int_t n, RVector b, RVector x, sparsesolverreport *rep);
void sparselusolve(sparsematrix *a, ZVector p, ZVector q, ae_int_t n, RVector b, RVector x, sparsesolverreport *rep);
} // end of namespace alglib_impl

namespace alglib {
// This structure is a sparse solver report.
//
// Following fields can be accessed by users:
DecClass(sparsesolverreport, ae_int_t &terminationtype;);

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
//     N       -   size of A, N > 0
//     IsUpper -   which half of A is provided (another half is ignored)
//     B       -   array[0..N-1], right part
//
// Outputs:
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate or non-SPD system).
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0    =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
void sparsesolvesks(const sparsematrix &a, const ae_int_t n, const bool isupper, const real_1d_array &b, sparsesolverreport &rep, real_1d_array &x);

// Sparse linear solver for A*x=b with N*N real  symmetric  positive definite
// matrix A given by its Cholesky decomposition, and N*1 vectors x and b.
//
// IMPORTANT: this solver requires input matrix to be in  the  SKS  (Skyline)
//            sparse storage format. An exception will be  generated  if  you
//            pass matrix in some other format (HASH or CRS).
//
// Inputs:
//     A       -   sparse NxN matrix stored in SKS format, must be NxN exactly
//     N       -   size of A, N > 0
//     IsUpper -   which half of A is provided (another half is ignored)
//     B       -   array[N], right part
//
// Outputs:
//     Rep     -   solver report, following fields are set:
//                 * rep.terminationtype - solver status; > 0 for success,
//                   set to -3 on failure (degenerate or non-SPD system).
//     X       -   array[N], it contains:
//                 * rep.terminationtype > 0    =>  solution
//                 * rep.terminationtype=-3   =>  filled by zeros
// ALGLIB: Copyright 26.12.2017 by Sergey Bochkanov
void sparsecholeskysolvesks(const sparsematrix &a, const ae_int_t n, const bool isupper, const real_1d_array &b, sparsesolverreport &rep, real_1d_array &x);

// Sparse linear solver for A*x=b with general (nonsymmetric) N*N sparse real
// matrix A, N*1 vectors x and b.
//
// This solver converts input matrix to CRS format, performs LU factorization
// and uses sparse triangular solvers to get solution of the original system.
//
// Inputs:
//     A       -   sparse matrix, must be NxN exactly, any storage format
//     N       -   size of A, N > 0
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
void sparsesolve(const sparsematrix &a, const ae_int_t n, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);

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
//     N       -   size of A, N > 0
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
void sparselusolve(const sparsematrix &a, const integer_1d_array &p, const integer_1d_array &q, const ae_int_t n, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep);
} // end of namespace alglib

// === LINCG Package ===
// Depends on: (LinAlg) SPARSE, MATGEN
namespace alglib_impl {
typedef struct {
   ae_vector rx;
   ae_vector b;
   ae_int_t n;
   ae_int_t prectype;
   ae_vector cx;
   ae_vector cr;
   ae_vector cz;
   ae_vector p;
   ae_vector r;
   ae_vector z;
   double alpha;
   double beta;
   double r2;
   double meritfunction;
   ae_vector x;
   ae_vector mv;
   ae_vector pv;
   double vmv;
   ae_vector startx;
   double epsf;
   ae_int_t maxits;
   ae_int_t itsbeforerestart;
   ae_int_t itsbeforerupdate;
   bool xrep;
   bool xupdated;
   bool needmv;
// bool needmtv; //(@) Not used.
// bool needmv2; //(@) Not used.
   bool needvmv;
   bool needprec;
   ae_int_t repiterationscount;
   ae_int_t repnmv;
   ae_int_t repterminationtype;
   bool running;
   ae_vector tmpd;
   ae_int_t PQ;
} lincgstate;
void lincgstate_init(void *_p, bool make_automatic);
void lincgstate_copy(void *_dst, void *_src, bool make_automatic);
void lincgstate_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t iterationscount;
   ae_int_t nmv;
   ae_int_t terminationtype;
   double r2;
} lincgreport;
void lincgreport_init(void *_p, bool make_automatic);
void lincgreport_copy(void *_dst, void *_src, bool make_automatic);
void lincgreport_free(void *_p, bool make_automatic);

void lincgcreate(ae_int_t n, lincgstate *state);
void lincgsetstartingpoint(lincgstate *state, RVector x);
void lincgsetb(lincgstate *state, RVector b);
void lincgsetprecunit(lincgstate *state);
void lincgsetprecdiag(lincgstate *state);
void lincgsetcond(lincgstate *state, double epsf, ae_int_t maxits);
bool lincgiteration(lincgstate *state);
void lincgsolvesparse(lincgstate *state, sparsematrix *a, bool isupper, RVector b);
void lincgresults(lincgstate *state, RVector x, lincgreport *rep);
void lincgsetrestartfreq(lincgstate *state, ae_int_t srf);
void lincgsetrupdatefreq(lincgstate *state, ae_int_t freq);
void lincgsetxrep(lincgstate *state, bool needxrep);
void lincgrestart(lincgstate *state);
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the linear CG method.
//
// You should use ALGLIB functions to work with this object.
// Never try to access its fields directly!
DecClass(lincgstate, EndD);
DecClass(lincgreport, ae_int_t &iterationscount; ae_int_t &nmv; ae_int_t &terminationtype; double &r2;);

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
void lincgcreate(const ae_int_t n, lincgstate &state);

// This function sets starting point.
// By default, zero starting point is used.
//
// Inputs:
//     X       -   starting point, array[N]
//
// Outputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x);

// This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
// function. By default, SolveSparse() uses diagonal preconditioner,  but  if
// you want to use solver without preconditioning, you can call this function
// which forces solver to use unit matrix for preconditioning.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
void lincgsetprecunit(const lincgstate &state);

// This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
// function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
// preconditioner. This preconditioning mode is active by default.
//
// Inputs:
//     State   -   structure which stores algorithm state
// ALGLIB: Copyright 19.11.2012 by Sergey Bochkanov
void lincgsetprecdiag(const lincgstate &state);

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
void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits);

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
void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b);

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
void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep);

// This function sets restart frequency. By default, algorithm  is  restarted
// after N subsequent iterations.
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf);

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
void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq);

// This function turns on/off reporting.
//
// Inputs:
//     State   -   structure which stores algorithm state
//     NeedXRep-   whether iteration reports are needed or not
//
// If NeedXRep is True, algorithm will call rep() callback function if  it is
// provided to MinCGOptimize().
// ALGLIB: Copyright 14.11.2011 by Sergey Bochkanov
void lincgsetxrep(const lincgstate &state, const bool needxrep);
} // end of namespace alglib

#endif // OnceOnly
