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
void gqgeneraterec(RVector alpha, RVector beta, double mu0, ae_int_t n, ae_int_t *info, RVector x, RVector w);
void gqgenerategausslobattorec(RVector alpha, RVector beta, double mu0, double a, double b, ae_int_t n, ae_int_t *info, RVector x, RVector w);
void gqgenerategaussradaurec(RVector alpha, RVector beta, double mu0, double a, ae_int_t n, ae_int_t *info, RVector x, RVector w);
void gqgenerategausslegendre(ae_int_t n, ae_int_t *info, RVector x, RVector w);
void gqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t *info, RVector x, RVector w);
void gqgenerategausslaguerre(ae_int_t n, double alpha, ae_int_t *info, RVector x, RVector w);
void gqgenerategausshermite(ae_int_t n, ae_int_t *info, RVector x, RVector w);
} // end of namespace alglib_impl

namespace alglib {
// Computation of nodes and weights for a Gauss quadrature formula
//
// The algorithm generates the N-point Gauss quadrature formula  with  weight
// function given by coefficients alpha and beta  of  a  recurrence  relation
// which generates a system of orthogonal polynomials:
//
// P-1(x)   =  0
// P0(x)    =  1
// Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)
//
// and zeroth moment Mu0
//
// Mu0 = integral(W(x)dx,a,b)
//
// Inputs:
//     Alpha   -   array[0..N-1], alpha coefficients
//     Beta    -   array[0..N-1], beta coefficients
//                 Zero-indexed element is not used and may be arbitrary.
//                 Beta[I]>0.
//     Mu0     -   zeroth moment of the weight function.
//     N       -   number of nodes of the quadrature formula, N >= 1
//
// Outputs:
//     Info    -   error code:
//                 * -3    internal eigenproblem solver hasn't converged
//                 * -2    Beta[i] <= 0
//                 * -1    incorrect N was passed
//                 *  1    OK
//     X       -   array[0..N-1] - array of quadrature nodes,
//                 in ascending order.
//     W       -   array[0..N-1] - array of quadrature weights.
//
// ALGLIB: Copyright 2005-2009 by Sergey Bochkanov
void gqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);

// Computation of nodes and weights for a Gauss-Lobatto quadrature formula
//
// The algorithm generates the N-point Gauss-Lobatto quadrature formula  with
// weight function given by coefficients alpha and beta of a recurrence which
// generates a system of orthogonal polynomials.
//
// P-1(x)   =  0
// P0(x)    =  1
// Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)
//
// and zeroth moment Mu0
//
// Mu0 = integral(W(x)dx,a,b)
//
// Inputs:
//     Alpha   -   array[0..N-2], alpha coefficients
//     Beta    -   array[0..N-2], beta coefficients.
//                 Zero-indexed element is not used, may be arbitrary.
//                 Beta[I]>0
//     Mu0     -   zeroth moment of the weighting function.
//     A       -   left boundary of the integration interval.
//     B       -   right boundary of the integration interval.
//     N       -   number of nodes of the quadrature formula, N >= 3
//                 (including the left and right boundary nodes).
//
// Outputs:
//     Info    -   error code:
//                 * -3    internal eigenproblem solver hasn't converged
//                 * -2    Beta[i] <= 0
//                 * -1    incorrect N was passed
//                 *  1    OK
//     X       -   array[0..N-1] - array of quadrature nodes,
//                 in ascending order.
//     W       -   array[0..N-1] - array of quadrature weights.
//
// ALGLIB: Copyright 2005-2009 by Sergey Bochkanov
void gqgenerategausslobattorec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const double b, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);

// Computation of nodes and weights for a Gauss-Radau quadrature formula
//
// The algorithm generates the N-point Gauss-Radau  quadrature  formula  with
// weight function given by the coefficients alpha and  beta  of a recurrence
// which generates a system of orthogonal polynomials.
//
// P-1(x)   =  0
// P0(x)    =  1
// Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)
//
// and zeroth moment Mu0
//
// Mu0 = integral(W(x)dx,a,b)
//
// Inputs:
//     Alpha   -   array[0..N-2], alpha coefficients.
//     Beta    -   array[0..N-1], beta coefficients
//                 Zero-indexed element is not used.
//                 Beta[I]>0
//     Mu0     -   zeroth moment of the weighting function.
//     A       -   left boundary of the integration interval.
//     N       -   number of nodes of the quadrature formula, N >= 2
//                 (including the left boundary node).
//
// Outputs:
//     Info    -   error code:
//                 * -3    internal eigenproblem solver hasn't converged
//                 * -2    Beta[i] <= 0
//                 * -1    incorrect N was passed
//                 *  1    OK
//     X       -   array[0..N-1] - array of quadrature nodes,
//                 in ascending order.
//     W       -   array[0..N-1] - array of quadrature weights.
//
//
// ALGLIB: Copyright 2005-2009 by Sergey Bochkanov
void gqgenerategaussradaurec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);

// Returns nodes/weights for Gauss-Legendre quadrature on [-1,1] with N
// nodes.
//
// Inputs:
//     N           -   number of nodes, >= 1
//
// Outputs:
//     Info        -   error code:
//                     * -4    an  error   was   detected   when  calculating
//                             weights/nodes.  N  is  too  large   to  obtain
//                             weights/nodes  with  high   enough   accuracy.
//                             Try  to   use   multiple   precision  version.
//                     * -3    internal eigenproblem solver hasn't  converged
//                     * -1    incorrect N was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes,
//                     in ascending order.
//     W           -   array[0..N-1] - array of quadrature weights.
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);

// Returns  nodes/weights  for  Gauss-Jacobi quadrature on [-1,1] with weight
// function W(x)=power(1-x,Alpha)*power(1+x,Beta).
//
// Inputs:
//     N           -   number of nodes, >= 1
//     Alpha       -   power-law coefficient, Alpha>-1
//     Beta        -   power-law coefficient, Beta>-1
//
// Outputs:
//     Info        -   error code:
//                     * -4    an  error  was   detected   when   calculating
//                             weights/nodes. Alpha or  Beta  are  too  close
//                             to -1 to obtain weights/nodes with high enough
//                             accuracy, or, may be, N is too large.  Try  to
//                             use multiple precision version.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -1    incorrect N/Alpha/Beta was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes,
//                     in ascending order.
//     W           -   array[0..N-1] - array of quadrature weights.
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &w);

// Returns  nodes/weights  for  Gauss-Laguerre  quadrature  on  [0,+inf) with
// weight function W(x)=power(x,Alpha)*exp(-x)
//
// Inputs:
//     N           -   number of nodes, >= 1
//     Alpha       -   power-law coefficient, Alpha>-1
//
// Outputs:
//     Info        -   error code:
//                     * -4    an  error  was   detected   when   calculating
//                             weights/nodes. Alpha is too  close  to  -1  to
//                             obtain weights/nodes with high enough accuracy
//                             or, may  be,  N  is  too  large.  Try  to  use
//                             multiple precision version.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -1    incorrect N/Alpha was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes,
//                     in ascending order.
//     W           -   array[0..N-1] - array of quadrature weights.
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gqgenerategausslaguerre(const ae_int_t n, const double alpha, ae_int_t &info, real_1d_array &x, real_1d_array &w);

// Returns  nodes/weights  for  Gauss-Hermite  quadrature on (-inf,+inf) with
// weight function W(x)=exp(-x*x)
//
// Inputs:
//     N           -   number of nodes, >= 1
//
// Outputs:
//     Info        -   error code:
//                     * -4    an  error  was   detected   when   calculating
//                             weights/nodes.  May be, N is too large. Try to
//                             use multiple precision version.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -1    incorrect N/Alpha was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes,
//                     in ascending order.
//     W           -   array[0..N-1] - array of quadrature weights.
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gqgenerategausshermite(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
} // end of namespace alglib

// === GKQ Package ===
// Depends on: GQ
namespace alglib_impl {
void gkqgeneraterec(RVector alpha, RVector beta, double mu0, ae_int_t n, ae_int_t *info, RVector x, RVector wkronrod, RVector wgauss);
void gkqgenerategausslegendre(ae_int_t n, ae_int_t *info, RVector x, RVector wkronrod, RVector wgauss);
void gkqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t *info, RVector x, RVector wkronrod, RVector wgauss);
void gkqlegendrecalc(ae_int_t n, ae_int_t *info, RVector x, RVector wkronrod, RVector wgauss);
void gkqlegendretbl(ae_int_t n, RVector x, RVector wkronrod, RVector wgauss, double *eps);
} // end of namespace alglib_impl

namespace alglib {
// Computation of nodes and weights of a Gauss-Kronrod quadrature formula
//
// The algorithm generates the N-point Gauss-Kronrod quadrature formula  with
// weight  function  given  by  coefficients  alpha  and beta of a recurrence
// relation which generates a system of orthogonal polynomials:
//
//     P-1(x)   =  0
//     P0(x)    =  1
//     Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)
//
// and zero moment Mu0
//
//     Mu0 = integral(W(x)dx,a,b)
//
//
// Inputs:
//     Alpha       -   alpha coefficients, array[0..floor(3*K/2)].
//     Beta        -   beta coefficients,  array[0..ceil(3*K/2)].
//                     Beta[0] is not used and may be arbitrary.
//                     Beta[I]>0.
//     Mu0         -   zeroth moment of the weight function.
//     N           -   number of nodes of the Gauss-Kronrod quadrature formula,
//                     N >= 3,
//                     N =  2*K+1.
//
// Outputs:
//     Info        -   error code:
//                     * -5    no real and positive Gauss-Kronrod formula can
//                             be created for such a weight function  with  a
//                             given number of nodes.
//                     * -4    N is too large, task may be ill  conditioned -
//                             x[i]=x[i+1] found.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -2    Beta[i] <= 0
//                     * -1    incorrect N was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes,
//                     in ascending order.
//     WKronrod    -   array[0..N-1] - Kronrod weights
//     WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
//                     corresponding to extended Kronrod nodes).
//
// ALGLIB: Copyright 08.05.2009 by Sergey Bochkanov
void gkqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);

// Returns   Gauss   and   Gauss-Kronrod   nodes/weights  for  Gauss-Legendre
// quadrature with N points.
//
// GKQLegendreCalc (calculation) or  GKQLegendreTbl  (precomputed  table)  is
// used depending on machine precision and number of nodes.
//
// Inputs:
//     N           -   number of Kronrod nodes, must be odd number, >= 3.
//
// Outputs:
//     Info        -   error code:
//                     * -4    an  error   was   detected   when  calculating
//                             weights/nodes.  N  is  too  large   to  obtain
//                             weights/nodes  with  high   enough   accuracy.
//                             Try  to   use   multiple   precision  version.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -1    incorrect N was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes, ordered in
//                     ascending order.
//     WKronrod    -   array[0..N-1] - Kronrod weights
//     WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
//                     corresponding to extended Kronrod nodes).
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gkqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);

// Returns   Gauss   and   Gauss-Kronrod   nodes/weights   for   Gauss-Jacobi
// quadrature on [-1,1] with weight function
//
//     W(x)=power(1-x,Alpha)*power(1+x,Beta).
//
// Inputs:
//     N           -   number of Kronrod nodes, must be odd number, >= 3.
//     Alpha       -   power-law coefficient, Alpha>-1
//     Beta        -   power-law coefficient, Beta>-1
//
// Outputs:
//     Info        -   error code:
//                     * -5    no real and positive Gauss-Kronrod formula can
//                             be created for such a weight function  with  a
//                             given number of nodes.
//                     * -4    an  error  was   detected   when   calculating
//                             weights/nodes. Alpha or  Beta  are  too  close
//                             to -1 to obtain weights/nodes with high enough
//                             accuracy, or, may be, N is too large.  Try  to
//                             use multiple precision version.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -1    incorrect N was passed
//                     * +1    OK
//                     * +2    OK, but quadrature rule have exterior  nodes,
//                             x[0]<-1 or x[n-1]>+1
//     X           -   array[0..N-1] - array of quadrature nodes, ordered in
//                     ascending order.
//     WKronrod    -   array[0..N-1] - Kronrod weights
//     WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
//                     corresponding to extended Kronrod nodes).
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gkqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);

// Returns Gauss and Gauss-Kronrod nodes for quadrature with N points.
//
// Reduction to tridiagonal eigenproblem is used.
//
// Inputs:
//     N           -   number of Kronrod nodes, must be odd number, >= 3.
//
// Outputs:
//     Info        -   error code:
//                     * -4    an  error   was   detected   when  calculating
//                             weights/nodes.  N  is  too  large   to  obtain
//                             weights/nodes  with  high   enough   accuracy.
//                             Try  to   use   multiple   precision  version.
//                     * -3    internal eigenproblem solver hasn't converged
//                     * -1    incorrect N was passed
//                     * +1    OK
//     X           -   array[0..N-1] - array of quadrature nodes, ordered in
//                     ascending order.
//     WKronrod    -   array[0..N-1] - Kronrod weights
//     WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
//                     corresponding to extended Kronrod nodes).
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gkqlegendrecalc(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);

// Returns Gauss and Gauss-Kronrod nodes for quadrature with N  points  using
// pre-calculated table. Nodes/weights were  computed  with  accuracy  up  to
// 1.0E-32 (if MPFR version of ALGLIB is used). In standard double  precision
// accuracy reduces to something about 2.0E-16 (depending  on your compiler's
// handling of long floating point constants).
//
// Inputs:
//     N           -   number of Kronrod nodes.
//                     N can be 15, 21, 31, 41, 51, 61.
//
// Outputs:
//     X           -   array[0..N-1] - array of quadrature nodes, ordered in
//                     ascending order.
//     WKronrod    -   array[0..N-1] - Kronrod weights
//     WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
//                     corresponding to extended Kronrod nodes).
//
//
// ALGLIB: Copyright 12.05.2009 by Sergey Bochkanov
void gkqlegendretbl(const ae_int_t n, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss, double &eps);
} // end of namespace alglib

// === AUTOGK Package ===
// Depends on: GKQ
namespace alglib_impl {
typedef struct {
   ae_int_t terminationtype;
   ae_int_t nfev;
   ae_int_t nintervals;
} autogkreport;
void autogkreport_init(void *_p, bool make_automatic);
void autogkreport_copy(void *_dst, void *_src, bool make_automatic);
void autogkreport_free(void *_p, bool make_automatic);

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
   ae_int_t PQ;
} autogkinternalstate;
void autogkinternalstate_init(void *_p, bool make_automatic);
void autogkinternalstate_copy(void *_dst, void *_src, bool make_automatic);
void autogkinternalstate_free(void *_p, bool make_automatic);

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
   ae_int_t PQ;
   double v;
   ae_int_t terminationtype;
   ae_int_t nfev;
   ae_int_t nintervals;
} autogkstate;
void autogkstate_init(void *_p, bool make_automatic);
void autogkstate_copy(void *_dst, void *_src, bool make_automatic);
void autogkstate_free(void *_p, bool make_automatic);

void autogksmooth(double a, double b, autogkstate *state);
void autogksmoothw(double a, double b, double xwidth, autogkstate *state);
void autogksingular(double a, double b, double alpha, double beta, autogkstate *state);
bool autogkiteration(autogkstate *state);
void autogkresults(autogkstate *state, double *v, autogkreport *rep);
} // end of namespace alglib_impl

namespace alglib {
// Integration report:
// * TerminationType = completetion code:
//     * -5    non-convergence of Gauss-Kronrod nodes
//             calculation subroutine.
//     * -1    incorrect parameters were specified
//     *  1    OK
// * Rep.NFEV countains number of function calculations
// * Rep.NIntervals contains number of intervals [a,b]
//   was partitioned into.
DecClass(autogkreport, ae_int_t &terminationtype; ae_int_t &nfev; ae_int_t &nintervals;);

// This structure stores state of the integration algorithm.
//
// Although this class has public fields,  they are not intended for external
// use. You should use ALGLIB functions to work with this class:
// * autogksmooth()/AutoGKSmoothW()/... to create objects
// * autogkintegrate() to begin integration
// * autogkresults() to get results
DecClass(autogkstate, bool &needf; double &x; double &xminusa; double &bminusx; double &f;);

// Integration of a smooth function F(x) on a finite interval [a,b].
//
// Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
// is calculated with accuracy close to the machine precision.
//
// Algorithm works well only with smooth integrands.  It  may  be  used  with
// continuous non-smooth integrands, but with  less  performance.
//
// It should never be used with integrands which have integrable singularities
// at lower or upper limits - algorithm may crash. Use AutoGKSingular in such
// cases.
//
// Inputs:
//     A, B    -   interval boundaries (A<B, A=B or A>B)
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// SEE ALSO
//     AutoGKSmoothW, AutoGKSingular, AutoGKResults.
//
//
// ALGLIB: Copyright 06.05.2009 by Sergey Bochkanov
void autogksmooth(const double a, const double b, autogkstate &state);

// Integration of a smooth function F(x) on a finite interval [a,b].
//
// This subroutine is same as AutoGKSmooth(), but it guarantees that interval
// [a,b] is partitioned into subintervals which have width at most XWidth.
//
// Subroutine  can  be  used  when  integrating nearly-constant function with
// narrow "bumps" (about XWidth wide). If "bumps" are too narrow, AutoGKSmooth
// subroutine can overlook them.
//
// Inputs:
//     A, B    -   interval boundaries (A<B, A=B or A>B)
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// SEE ALSO
//     AutoGKSmooth, AutoGKSingular, AutoGKResults.
//
//
// ALGLIB: Copyright 06.05.2009 by Sergey Bochkanov
void autogksmoothw(const double a, const double b, const double xwidth, autogkstate &state);

// Integration on a finite interval [A,B].
// Integrand have integrable singularities at A/B.
//
// F(X) must diverge as "(x-A)^alpha" at A, as "(B-x)^beta" at B,  with known
// alpha/beta (alpha>-1, beta>-1).  If alpha/beta  are  not known,  estimates
// from below can be used (but these estimates should be greater than -1 too).
//
// One  of  alpha/beta variables (or even both alpha/beta) may be equal to 0,
// which means than function F(x) is non-singular at A/B. Anyway (singular at
// bounds or not), function F(x) is supposed to be continuous on (A,B).
//
// Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
// is calculated with accuracy close to the machine precision.
//
// Inputs:
//     A, B    -   interval boundaries (A<B, A=B or A>B)
//     Alpha   -   power-law coefficient of the F(x) at A,
//                 Alpha>-1
//     Beta    -   power-law coefficient of the F(x) at B,
//                 Beta>-1
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// SEE ALSO
//     AutoGKSmooth, AutoGKSmoothW, AutoGKResults.
//
//
// ALGLIB: Copyright 06.05.2009 by Sergey Bochkanov
void autogksingular(const double a, const double b, const double alpha, const double beta, autogkstate &state);

// This function provides reverse communication interface
// Reverse communication interface is not documented or recommended to use.
// See below for functions which provide better documented API
bool autogkiteration(const autogkstate &state);

// This function is used to launcn iterations of the 1-dimensional integrator
//
// It accepts following parameters:
//     func    -   callback which calculates f(x) for given x
//     ptr     -   optional pointer which is passed to func; can be NULL
//
//
// ALGLIB: Copyright 07.05.2009 by Sergey Bochkanov
//
void autogkintegrate(autogkstate &state, void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr), void *ptr = NULL);

// Adaptive integration results
//
// Called after AutoGKIteration returned False.
//
// Inputs:
//     State   -   algorithm state (used by AutoGKIteration).
//
// Outputs:
//     V       -   integral(f(x)dx,a,b)
//     Rep     -   optimization report (see AutoGKReport description)
//
// ALGLIB: Copyright 14.11.2007 by Sergey Bochkanov
void autogkresults(const autogkstate &state, double &v, autogkreport &rep);
} // end of namespace alglib

#endif // OnceOnly
