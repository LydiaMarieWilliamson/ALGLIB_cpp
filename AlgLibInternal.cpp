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
// Must be defined before we include anything depending on any kernel headers.
#define _ALGLIB_IMPL_DEFINES
#define _ALGLIB_INTEGRITY_CHECKS_ONCE
#include "AlgLibInternal.h"

// === APSERV Package ===
namespace alglib_impl {
static const ae_int_t apserv_maxtemporariesinnpool = 1000;

// The function performs zero-coalescing on real value.
//
// NOTE: no check is performed for B != 0
// ALGLIB: Copyright 18.05.2015 by Sergey Bochkanov
double coalesce(double a, double b) {
   double result;
   result = a;
   if (a == 0.0) {
      result = b;
   }
   return result;
}

// The function performs zero-coalescing on integer value.
//
// NOTE: no check is performed for B != 0
// ALGLIB: Copyright 18.05.2015 by Sergey Bochkanov
ae_int_t coalescei(ae_int_t a, ae_int_t b) {
   ae_int_t result;
   result = a;
   if (a == 0) {
      result = b;
   }
   return result;
}

// The function calculates binary logarithm.
//
// NOTE: it costs twice as much as Ln(x)
// ALGLIB: Copyright 17.09.2012 by Sergey Bochkanov
double logbase2(double x) {
   const double log2e = 1.0 / log(2.0);
   double result;
   result = log(x) * log2e;
   return result;
}

// This  function  generates  1-dimensional  general  interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N == 1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1d(double a, double b, ae_int_t n, RVector *x, RVector *y) {
   ae_int_t i;
   double h;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "taskgenint1d: n < 1!");
   ae_vector_set_length(x, n);
   ae_vector_set_length(y, n);
   if (n > 1) {
      x->xR[0] = a;
      y->xR[0] = randommid();
      h = (b - a) / (n - 1);
      for (i = 1; i < n; i++) {
         if (i != n - 1) {
            x->xR[i] = a + (i + 0.2 * randommid()) * h;
         } else {
            x->xR[i] = b;
         }
         y->xR[i] = y->xR[i - 1] + randommid() * (x->xR[i] - x->xR[i - 1]);
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = randommid();
   }
}

// This function generates  1-dimensional equidistant interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N == 1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1dequidist(double a, double b, ae_int_t n, RVector *x, RVector *y) {
   ae_int_t i;
   double h;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "taskgenint1dequidist: n < 1!");
   ae_vector_set_length(x, n);
   ae_vector_set_length(y, n);
   if (n > 1) {
      x->xR[0] = a;
      y->xR[0] = randommid();
      h = (b - a) / (n - 1);
      for (i = 1; i < n; i++) {
         x->xR[i] = a + i * h;
         y->xR[i] = y->xR[i - 1] + randommid() * h;
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = randommid();
   }
}

// This function generates  1-dimensional Chebyshev-1 interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N == 1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1dcheb1(double a, double b, ae_int_t n, RVector *x, RVector *y) {
   ae_int_t i;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "taskgenint1dcheb1: n < 1!");
   ae_vector_set_length(x, n);
   ae_vector_set_length(y, n);
   if (n > 1) {
      for (i = 0; i < n; i++) {
         x->xR[i] = 0.5 * (b + a) + 0.5 * (b - a) * cos(pi * (2 * i + 1) / (2 * n));
         if (i == 0) {
            y->xR[i] = randommid();
         } else {
            y->xR[i] = y->xR[i - 1] + randommid() * (x->xR[i] - x->xR[i - 1]);
         }
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = randommid();
   }
}

// This function generates  1-dimensional Chebyshev-2 interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N == 1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1dcheb2(double a, double b, ae_int_t n, RVector *x, RVector *y) {
   ae_int_t i;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "taskgenint1dcheb2: n < 1!");
   ae_vector_set_length(x, n);
   ae_vector_set_length(y, n);
   if (n > 1) {
      for (i = 0; i < n; i++) {
         x->xR[i] = 0.5 * (b + a) + 0.5 * (b - a) * cos(pi * i / (n - 1));
         if (i == 0) {
            y->xR[i] = randommid();
         } else {
            y->xR[i] = y->xR[i - 1] + randommid() * (x->xR[i] - x->xR[i - 1]);
         }
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = randommid();
   }
}

// This function checks that all values from X[] are distinct. It does more
// than just usual floating point comparison:
// * first, it calculates max(X) and min(X)
// * second, it maps X[] from [min,max] to [0,1]
// * only at this stage actual comparison is done
//
// The meaning of such check is to ensure that all values are "distinct enough"
// and will not cause interpolation subroutine to fail.
//
// NOTE:
//     X[] must be sorted by ascending (subroutine ASSERT's it)
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
bool aredistinct(RVector *x, ae_int_t n) {
   double a;
   double b;
   ae_int_t i;
   bool nonsorted;
   bool result;
   ae_assert(n >= 1, "aredistinct: internal error (n < 1)");
   if (n == 1) {
   // everything is alright, it is up to caller to decide whether it
   // can interpolate something with just one point
      result = true;
      return result;
   }
   a = x->xR[0];
   b = x->xR[0];
   nonsorted = false;
   for (i = 1; i < n; i++) {
      a = rmin2(a, x->xR[i]);
      b = rmax2(b, x->xR[i]);
      nonsorted = nonsorted || x->xR[i - 1] >= x->xR[i];
   }
   ae_assert(!nonsorted, "aredistinct: internal error (not sorted)");
   for (i = 1; i < n; i++) {
      if ((x->xR[i] - a) / (b - a) == (x->xR[i - 1] - a) / (b - a)) {
         result = false;
         return result;
      }
   }
   result = true;
   return result;
}

// Resizes X and fills by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void setlengthzero(RVector *x, ae_int_t n) {
   ae_int_t i;
   ae_assert(n >= 0, "setlengthzero: n < 0");
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      x->xR[i] = 0.0;
   }
}

// If Length(X) < N, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void vectorsetlengthatleast(ae_vector *x, ae_int_t n) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n);
   }
}

// If Cols(X) < N or Rows(X) < M, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void matrixsetlengthatleast(ae_matrix *x, ae_int_t m, ae_int_t n) {
   if (x->rows < m || x->cols < n) {
      if (m < x->rows) m = x->rows; else if (n < x->cols) n = x->cols;
      ae_matrix_set_length(x, m, n);
   }
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) new size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void bvectorgrowto(BVector *x, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewVector(oldx, 0, DT_BOOL);
// Enough place
   if (x->cnt >= n) {
      ae_frame_leave();
      return;
   }
// Choose new size
   n = imax2(n, iround(1.8 * x->cnt + 1.0));
// Grow
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xB[i] = oldx.xB[i];
      } else {
         x->xB[i] = false;
      }
   }
   ae_frame_leave();
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) new size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void ivectorgrowto(ZVector *x, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewVector(oldx, 0, DT_INT);
// Enough place
   if (x->cnt >= n) {
      ae_frame_leave();
      return;
   }
// Choose new size
   n = imax2(n, iround(1.8 * x->cnt + 1.0));
// Grow
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xZ[i] = oldx.xZ[i];
      } else {
         x->xZ[i] = 0;
      }
   }
   ae_frame_leave();
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) new size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rvectorgrowto(RVector *x, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewVector(oldx, 0, DT_REAL);
// Enough place
   if (x->cnt >= n) {
      ae_frame_leave();
      return;
   }
// Choose new size
   n = imax2(n, iround(1.8 * x->cnt + 1.0));
// Grow
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xR[i] = oldx.xR[i];
      } else {
         x->xR[i] = 0.0;
      }
   }
   ae_frame_leave();
}

// Grows X, i.e. appends cols in such a way that:
// a) contents is preserved
// b) new col count is at least N
// c) new col count can be larger than N, so subsequent grow() calls can return
//    without reallocation
// d) new matrix has at least MinRows row (if less than specified amount
//    of rows is present, new rows are added with undefined contents);
//    MinRows can be 0 or negative value = ignored
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rmatrixgrowcolsto(RMatrix *a, ae_int_t n, ae_int_t minrows) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t n2;
   ae_int_t m;
   ae_frame_make(&_frame_block);
   NewMatrix(olda, 0, 0, DT_REAL);
// Enough place?
   if (a->cols >= n && a->rows >= minrows) {
      ae_frame_leave();
      return;
   }
// Sizes and metrics
   if (a->cols < n) {
      n = imax2(n, iround(1.8 * a->cols + 1.0));
   }
   n2 = imin2(a->cols, n);
   m = a->rows;
// Grow
   ae_swap_matrices(a, &olda);
   ae_matrix_set_length(a, imax2(m, minrows), n);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n2; j++) {
         a->xyR[i][j] = olda.xyR[i][j];
      }
   }
   ae_frame_leave();
}

// Grows X, i.e. appends rows in such a way that:
// a) contents is preserved
// b) new row count is at least N
// c) new row count can be larger than N, so subsequent grow() calls can return
//    without reallocation
// d) new matrix has at least MinCols columns (if less than specified amount
//    of columns is present, new columns are added with undefined contents);
//    MinCols can be 0 or negative value = ignored
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rmatrixgrowrowsto(RMatrix *a, ae_int_t n, ae_int_t mincols) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t n2;
   ae_int_t m;
   ae_frame_make(&_frame_block);
   NewMatrix(olda, 0, 0, DT_REAL);
// Enough place?
   if (a->rows >= n && a->cols >= mincols) {
      ae_frame_leave();
      return;
   }
// Sizes and metrics
   if (a->rows < n) {
      n = imax2(n, iround(1.8 * a->rows + 1.0));
   }
   n2 = imin2(a->rows, n);
   m = a->cols;
// Grow
   ae_swap_matrices(a, &olda);
   ae_matrix_set_length(a, n, imax2(m, mincols));
   for (i = 0; i < n2; i++) {
      for (j = 0; j < m; j++) {
         a->xyR[i][j] = olda.xyR[i][j];
      }
   }
   ae_frame_leave();
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void ivectorresize(ZVector *x, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewVector(oldx, 0, DT_INT);
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xZ[i] = oldx.xZ[i];
      } else {
         x->xZ[i] = 0;
      }
   }
   ae_frame_leave();
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rvectorresize(RVector *x, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewVector(oldx, 0, DT_REAL);
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xR[i] = oldx.xR[i];
      } else {
         x->xR[i] = 0.0;
      }
   }
   ae_frame_leave();
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void imatrixresize(ZMatrix *x, ae_int_t m, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t m2;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewMatrix(oldx, 0, 0, DT_INT);
   m2 = x->rows;
   n2 = x->cols;
   ae_swap_matrices(x, &oldx);
   ae_matrix_set_length(x, m, n);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (i < m2 && j < n2) {
            x->xyZ[i][j] = oldx.xyZ[i][j];
         } else {
            x->xyZ[i][j] = 0;
         }
      }
   }
   ae_frame_leave();
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rmatrixresize(RMatrix *x, ae_int_t m, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t m2;
   ae_int_t n2;
   ae_frame_make(&_frame_block);
   NewMatrix(oldx, 0, 0, DT_REAL);
   m2 = x->rows;
   n2 = x->cols;
   ae_swap_matrices(x, &oldx);
   ae_matrix_set_length(x, m, n);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (i < m2 && j < n2) {
            x->xyR[i][j] = oldx.xyR[i][j];
         } else {
            x->xyR[i][j] = 0.0;
         }
      }
   }
   ae_frame_leave();
}

// This function checks that length(X) is at least N and first N values  from
// X[] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool isfinitevector(RVector *x, ae_int_t n) {
   ae_int_t i;
   double v;
   bool result;
   ae_assert(n >= 0, "isfinitevector: internal error (n < 0)");
   if (n == 0) {
      result = true;
      return result;
   }
   if (x->cnt < n) {
      result = false;
      return result;
   }
   v = 0.0;
   for (i = 0; i < n; i++) {
      v = 0.01 * v + x->xR[i];
   }
   result = isfinite(v);
   return result;
}

// This function checks that first N values from X[] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool isfinitecvector(CVector *z, ae_int_t n) {
   ae_int_t i;
   bool result;
   ae_assert(n >= 0, "isfinitecvector: internal error (n < 0)");
   for (i = 0; i < n; i++) {
      if (!isfinite(z->xC[i].x) || !isfinite(z->xC[i].y)) {
         result = false;
         return result;
      }
   }
   result = true;
   return result;
}

// This function checks that length(X) is at least N and first N values  from
// X[] are finite or NANs
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool isfiniteornanvector(RVector *x, ae_int_t n) {
   ae_int_t i;
   double v;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteVector: internal error (N < 0)");
   if (n == 0) {
      result = true;
      return result;
   }
   if (x->cnt < n) {
      result = false;
      return result;
   }
// Is it entirely finite?
   v = 0.0;
   for (i = 0; i < n; i++) {
      v = 0.01 * v + x->xR[i];
   }
   if (isfinite(v)) {
      result = true;
      return result;
   }
// OK, check that either finite or nan
   for (i = 0; i < n; i++) {
      if (!isfinite(x->xR[i]) && !isnan(x->xR[i])) {
         result = false;
         return result;
      }
   }
   result = true;
   return result;
}

// This function checks that size of X is at least MxN and values from
// X[0..M-1,0..N-1] are finite.
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool apservisfinitematrix(RMatrix *x, ae_int_t m, ae_int_t n) {
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "apservisfinitematrix: internal error (n < 0)");
   ae_assert(m >= 0, "apservisfinitematrix: internal error (m < 0)");
   if (m == 0 || n == 0) {
      result = true;
      return result;
   }
   if (x->rows < m || x->cols < n) {
      result = false;
      return result;
   }
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (!isfinite(x->xyR[i][j])) {
            result = false;
            return result;
         }
      }
   }
   result = true;
   return result;
}

// This function checks that all values from X[0..M-1,0..N-1] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool apservisfinitecmatrix(CMatrix *x, ae_int_t m, ae_int_t n) {
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "apservisfinitecmatrix: internal error (n < 0)");
   ae_assert(m >= 0, "apservisfinitecmatrix: internal error (m < 0)");
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (!isfinite(x->xyC[i][j].x) || !isfinite(x->xyC[i][j].y)) {
            result = false;
            return result;
         }
      }
   }
   result = true;
   return result;
}

// This function checks that size of X is at least NxN and all values from
// upper/lower triangle of X[0..N-1,0..N-1] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool isfinitertrmatrix(RMatrix *x, ae_int_t n, bool isupper) {
   ae_int_t i;
   ae_int_t j1;
   ae_int_t j2;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "isfinitertrmatrix: internal error (n < 0)");
   if (n == 0) {
      result = true;
      return result;
   }
   if (x->rows < n || x->cols < n) {
      result = false;
      return result;
   }
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i;
      }
      for (j = j1; j <= j2; j++) {
         if (!isfinite(x->xyR[i][j])) {
            result = false;
            return result;
         }
      }
   }
   result = true;
   return result;
}

// This function checks that all values from upper/lower triangle of
// X[0..N-1,0..N-1] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool apservisfinitectrmatrix(CMatrix *x, ae_int_t n, bool isupper) {
   ae_int_t i;
   ae_int_t j1;
   ae_int_t j2;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "apservisfinitectrmatrix: internal error (n < 0)");
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i;
      }
      for (j = j1; j <= j2; j++) {
         if (!isfinite(x->xyC[i][j].x) || !isfinite(x->xyC[i][j].y)) {
            result = false;
            return result;
         }
      }
   }
   result = true;
   return result;
}

// This function checks that all values from X[0..M-1,0..N-1] are  finite  or
// NaN's.
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool apservisfiniteornanmatrix(RMatrix *x, ae_int_t m, ae_int_t n) {
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "apservisfiniteornanmatrix: internal error (n < 0)");
   ae_assert(m >= 0, "apservisfiniteornanmatrix: internal error (m < 0)");
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (!(isfinite(x->xyR[i][j]) || isnan(x->xyR[i][j]))) {
            result = false;
            return result;
         }
      }
   }
   result = true;
   return result;
}

// Safe sqrt(x^2+y^2)
// ALGLIB: Copyright by Sergey Bochkanov
double safepythag2(double x, double y) {
   double w;
   double xabs;
   double yabs;
   double z;
   double result;
   xabs = fabs(x);
   yabs = fabs(y);
   w = rmax2(xabs, yabs);
   z = rmin2(xabs, yabs);
   if (z == 0.0) {
      result = w;
   } else {
      result = w * sqrt(1.0 + sqr(z / w));
   }
   return result;
}

// Safe sqrt(x^2+y^2+z^2)
// ALGLIB: Copyright by Sergey Bochkanov
double safepythag3(double x, double y, double z) {
   double w;
   double result;
   w = rmax2(fabs(x), rmax2(fabs(y), fabs(z)));
   if (w == 0.0) {
      result = 0.0;
      return result;
   }
   x /= w;
   y /= w;
   z /= w;
   result = w * sqrt(sqr(x) + sqr(y) + sqr(z));
   return result;
}

// Safe division.
//
// This function attempts to calculate R == X/Y without overflow.
//
// It returns:
// * +1, if abs(X/Y) >= MaxRealNumber or undefined - overflow-like situation
//       (no overlfow is generated, R is either NAN, PosINF, NegINF)
// *  0, if MinRealNumber < abs(X/Y) < MaxRealNumber or X == 0, Y != 0
//       (R contains result, may be zero)
// * -1, if 0 < abs(X/Y) < MinRealNumber - underflow-like situation
//       (R contains zero; it corresponds to underflow)
//
// No overflow is generated in any case.
// ALGLIB: Copyright by Sergey Bochkanov
ae_int_t saferdiv(double x, double y, double *r) {
   ae_int_t result;
   *r = 0.0;
// Two special cases:
// * Y == 0
// * X == 0 and Y != 0
   if (y == 0.0) {
      result = 1;
      if (x == 0.0) {
         *r = NAN;
      }
      if (x > 0.0) {
         *r = +INFINITY;
      }
      if (x < 0.0) {
         *r = -INFINITY;
      }
      return result;
   }
   if (x == 0.0) {
      *r = 0.0;
      result = 0;
      return result;
   }
// make Y == 0
   if (y < 0.0) {
      x = -x;
      y = -y;
   }
//
   if (y >= 1.0) {
      *r = x / y;
      if (SmallAtR(*r, minrealnumber)) {
         result = -1;
         *r = 0.0;
      } else {
         result = 0;
      }
   } else {
      if (!SmallR(x, maxrealnumber * y)) {
         if (x > 0.0) {
            *r = +INFINITY;
         } else {
            *r = -INFINITY;
         }
         result = 1;
      } else {
         *r = x / y;
         result = 0;
      }
   }
   return result;
}

// This function calculates "safe" min(X/Y,V) for positive finite X, Y, V.
// No overflow is generated in any case.
// ALGLIB: Copyright by Sergey Bochkanov
double safeminposrv(double x, double y, double v) {
   double r;
   double result;
   if (y >= 1.0) {
   // Y >= 1, we can safely divide by Y
      r = x / y;
      result = v;
      if (v > r) {
         result = r;
      } else {
         result = v;
      }
   } else {
   // Y < 1, we can safely multiply by Y
      if (x < v * y) {
         result = x / y;
      } else {
         result = v;
      }
   }
   return result;
}

// This function makes periodic mapping of X to [A,B].
//
// It accepts X, A, B (A > B). It returns T which lies in  [A,B] and integer K,
// such that X = T + K*(B-A).
//
// NOTES:
// * K is represented as real value, although actually it is integer
// * T is guaranteed to be in [A,B]
// * T replaces X
// ALGLIB: Copyright by Sergey Bochkanov
void apperiodicmap(double *x, double a, double b, double *k) {
   *k = 0.0;
   ae_assert(a < b, "apperiodicmap: internal error!");
   *k = floor((*x - a) / (b - a));
   *x -= *k * (b - a);
   while (*x < a) {
      *x += b - a;
      --*k;
   }
   while (*x > b) {
      *x -= b - a;
      ++*k;
   }
   *x = rmax2(*x, a);
   *x = rmin2(*x, b);
}

// Returns random normal number using low-quality system-provided generator
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
double randomnormal() {
   double u;
   double v;
   double s;
   double result;
   while (true) {
      u = randommid();
      v = randommid();
      s = sqr(u) + sqr(v);
      if (s > 0.0 && s < 1.0) {
      // two sqrt's instead of one to
      // avoid overflow when S is too small
         s = sqrt(-2.0 * log(s)) / sqrt(s);
         result = u * s;
         break;
      }
   }
   return result;
}

// Generates random unit vector using low-quality system-provided generator.
// Reallocates array if its size is too short.
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void randomunit(ae_int_t n, RVector *x) {
   ae_int_t i;
   double v;
   double vv;
   ae_assert(n > 0, "randomunit: N <= 0");
   if (x->cnt < n) {
      ae_vector_set_length(x, n);
   }
   do {
      v = 0.0;
      for (i = 0; i < n; i++) {
         vv = randomnormal();
         x->xR[i] = vv;
         v += vv * vv;
      }
   } while (v <= 0.0);
   v = 1.0 / sqrt(v);
   for (i = 0; i < n; i++) {
      x->xR[i] *= v;
   }
}

// This function is used to swap two boolean values.
void swapb(bool *v0, bool *v1) {
   bool v;
   v = *v0;
   *v0 = *v1;
   *v1 = v;
}

// This function is used to swap two integer values.
void swapi(ae_int_t *v0, ae_int_t *v1) {
   ae_int_t v;
   v = *v0;
   *v0 = *v1;
   *v1 = v;
}

// This function is used to swap two real values.
void swapr(double *v0, double *v1) {
   double v;
   v = *v0;
   *v0 = *v1;
   *v1 = v;
}

// This function is used to swap two complex values.
void swapc(complex *v0, complex *v1) {
   swapr(&v0->x, &v1->x);
   swapr(&v0->y, &v1->y);
}

// This function is used to swap two cols of the matrix; if NRows < 0, automatically
// determined from the matrix size.
void swapcols(RMatrix *a, ae_int_t j0, ae_int_t j1, ae_int_t nrows) {
   ae_int_t i;
   if (j0 == j1) {
      return;
   }
   if (nrows < 0) {
      nrows = a->rows;
   }
   for (i = 0; i < nrows; i++) {
      swapr(&a->xyR[i][j0], &a->xyR[i][j1]);
   }
}

// This function is used to swap two rows of the matrix; if NCols < 0, automatically
// determined from the matrix size.
void swaprows(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t ncols) {
   ae_int_t j;
   if (i0 == i1) {
      return;
   }
   if (ncols < 0) {
      ncols = a->cols;
   }
   for (j = 0; j < ncols; j++) {
      swapr(&a->xyR[i0][j], &a->xyR[i1][j]);
   }
}

// This function is used to swap two "entries" in 1-dimensional array composed
// from D-element entries
void swapentries(RVector *a, ae_int_t i0, ae_int_t i1, ae_int_t entrywidth) {
   ae_int_t offs0;
   ae_int_t offs1;
   ae_int_t j;
   if (i0 == i1) {
      return;
   }
   offs0 = i0 * entrywidth;
   offs1 = i1 * entrywidth;
   for (j = 0; j < entrywidth; j++) {
      swapr(&a->xR[offs0 + j], &a->xR[offs1 + j]);
   }
}

// This function is used to swap two "entries" in 1-dimensional array composed
// from D-element entries
void swapentriesb(BVector *a, ae_int_t i0, ae_int_t i1, ae_int_t entrywidth) {
   ae_int_t offs0;
   ae_int_t offs1;
   ae_int_t j;
   if (i0 == i1) {
      return;
   }
   offs0 = i0 * entrywidth;
   offs1 = i1 * entrywidth;
   for (j = 0; j < entrywidth; j++) {
      swapb(&a->xB[offs0 + j], &a->xB[offs1 + j]);
   }
}

// This function is used to swap two elements of the vector
void swapelementsi(ZVector *a, ae_int_t i0, ae_int_t i1) {
   if (i0 == i1) {
      return;
   }
   swapi(&a->xZ[i0], &a->xZ[i1]);
}

// This function is used to swap two elements of the vector
void swapelements(RVector *a, ae_int_t i0, ae_int_t i1) {
   if (i0 == i1) {
      return;
   }
   swapr(&a->xR[i0], &a->xR[i1]);
}

// This function returns +1 or -1 depending on sign of X.
// x == 0 results in +1 being returned.
double possign(double x) {
   double result;
   if (x >= 0.0) {
      result = 1.0;
   } else {
      result = -1.0;
   }
   return result;
}

// Return ceil(a/b).
// It expects that a > 0, b > 0, but does not check it.
ae_int_t idivup(ae_int_t a, ae_int_t b) {
   return (a + b - 1) / b;
}

// Allocation of serializer: complex value
void alloccomplex(ae_serializer *s, complex v) {
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
}

// Serialization: complex value
void serializecomplex(ae_serializer *s, complex v) {
   ae_serializer_serialize_double(s, v.x);
   ae_serializer_serialize_double(s, v.y);
}

// Unserialization: complex value
complex unserializecomplex(ae_serializer *s) {
   complex result;
   result.x = ae_serializer_unserialize_double(s);
   result.y = ae_serializer_unserialize_double(s);
   return result;
}

// Allocation of serializer: boolean array
void allocbooleanarray(ae_serializer *s, BVector *v, ae_int_t n) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_alloc_entry(s);
   for (i = 0; i < n; i++) {
      ae_serializer_alloc_entry(s);
   }
}

// Serialization: boolean array
void serializebooleanarray(ae_serializer *s, BVector *v, ae_int_t n) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_serialize_int(s, n);
   for (i = 0; i < n; i++) {
      ae_serializer_serialize_bool(s, v->xB[i]);
   }
}

// Unserialization: boolean value
void unserializebooleanarray(ae_serializer *s, BVector *v) {
   ae_int_t n;
   ae_int_t i;
   SetVector(v);
   n = ae_serializer_unserialize_int(s);
   if (n == 0) {
      return;
   }
   ae_vector_set_length(v, n);
   for (i = 0; i < n; i++) {
      v->xB[i] = ae_serializer_unserialize_bool(s);
   }
}

// Allocation of serializer: real array
void allocrealarray(ae_serializer *s, RVector *v, ae_int_t n) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_alloc_entry(s);
   for (i = 0; i < n; i++) {
      ae_serializer_alloc_entry(s);
   }
}

// Serialization: complex value
void serializerealarray(ae_serializer *s, RVector *v, ae_int_t n) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_serialize_int(s, n);
   for (i = 0; i < n; i++) {
      ae_serializer_serialize_double(s, v->xR[i]);
   }
}

// Unserialization: complex value
void unserializerealarray(ae_serializer *s, RVector *v) {
   ae_int_t n;
   ae_int_t i;
   SetVector(v);
   n = ae_serializer_unserialize_int(s);
   if (n == 0) {
      return;
   }
   ae_vector_set_length(v, n);
   for (i = 0; i < n; i++) {
      v->xR[i] = ae_serializer_unserialize_double(s);
   }
}

// Allocation of serializer: Integer array
void allocintegerarray(ae_serializer *s, ZVector *v, ae_int_t n) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_alloc_entry(s);
   for (i = 0; i < n; i++) {
      ae_serializer_alloc_entry(s);
   }
}

// Serialization: Integer array
void serializeintegerarray(ae_serializer *s, ZVector *v, ae_int_t n) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_serialize_int(s, n);
   for (i = 0; i < n; i++) {
      ae_serializer_serialize_int(s, v->xZ[i]);
   }
}

// Unserialization: complex value
void unserializeintegerarray(ae_serializer *s, ZVector *v) {
   ae_int_t n;
   ae_int_t i;
   SetVector(v);
   n = ae_serializer_unserialize_int(s);
   if (n == 0) {
      return;
   }
   ae_vector_set_length(v, n);
   for (i = 0; i < n; i++) {
      v->xZ[i] = ae_serializer_unserialize_int(s);
   }
}

// Allocation of serializer: real matrix
void allocrealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1) {
   ae_int_t i;
   ae_int_t j;
   if (n0 < 0) {
      n0 = v->rows;
   }
   if (n1 < 0) {
      n1 = v->cols;
   }
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
   for (i = 0; i < n0; i++) {
      for (j = 0; j < n1; j++) {
         ae_serializer_alloc_entry(s);
      }
   }
}

// Serialization: complex value
void serializerealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1) {
   ae_int_t i;
   ae_int_t j;
   if (n0 < 0) {
      n0 = v->rows;
   }
   if (n1 < 0) {
      n1 = v->cols;
   }
   ae_serializer_serialize_int(s, n0);
   ae_serializer_serialize_int(s, n1);
   for (i = 0; i < n0; i++) {
      for (j = 0; j < n1; j++) {
         ae_serializer_serialize_double(s, v->xyR[i][j]);
      }
   }
}

// Unserialization: complex value
void unserializerealmatrix(ae_serializer *s, RMatrix *v) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t n0;
   ae_int_t n1;
   SetMatrix(v);
   n0 = ae_serializer_unserialize_int(s);
   n1 = ae_serializer_unserialize_int(s);
   if (n0 == 0 || n1 == 0) {
      return;
   }
   ae_matrix_set_length(v, n0, n1);
   for (i = 0; i < n0; i++) {
      for (j = 0; j < n1; j++) {
         v->xyR[i][j] = ae_serializer_unserialize_double(s);
      }
   }
}

// Copy boolean array
void copybooleanarray(BVector *src, BVector *dst) {
   ae_int_t i;
   SetVector(dst);
   if (src->cnt > 0) {
      ae_vector_set_length(dst, src->cnt);
      for (i = 0; i < src->cnt; i++) {
         dst->xB[i] = src->xB[i];
      }
   }
}

// Copy integer array
void copyintegerarray(ZVector *src, ZVector *dst) {
   ae_int_t i;
   SetVector(dst);
   if (src->cnt > 0) {
      ae_vector_set_length(dst, src->cnt);
      for (i = 0; i < src->cnt; i++) {
         dst->xZ[i] = src->xZ[i];
      }
   }
}

// Copy real array
void copyrealarray(RVector *src, RVector *dst) {
   ae_int_t i;
   SetVector(dst);
   if (src->cnt > 0) {
      ae_vector_set_length(dst, src->cnt);
      for (i = 0; i < src->cnt; i++) {
         dst->xR[i] = src->xR[i];
      }
   }
}

// Copy real matrix
void copyrealmatrix(RMatrix *src, RMatrix *dst) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(dst);
   if (src->rows > 0 && src->cols > 0) {
      ae_matrix_set_length(dst, src->rows, src->cols);
      for (i = 0; i < src->rows; i++) {
         for (j = 0; j < src->cols; j++) {
            dst->xyR[i][j] = src->xyR[i][j];
         }
      }
   }
}

// Clears integer array
void unsetintegerarray(ZVector *a) {
   SetVector(a), ae_vector_set_length(a, 0);
}

// Clears real array
void unsetrealarray(RVector *a) {
   SetVector(a), ae_vector_set_length(a, 0);
}

// Clears real matrix
void unsetrealmatrix(RMatrix *a) {
   SetMatrix(a), ae_matrix_set_length(a, 0, 0);
}

// Initialize nbPool - prepare it to store N-length arrays, N >= 0.
// Tries to reuse previously allocated memory as much as possible.
void nbpoolinit(nbpool *pool, ae_int_t n) {
   ae_assert(n >= 0, "niPoolInit: N < 0");
   pool->n = n;
   pool->temporariescount = 0;
   if (n == 0) {
      return;
   }
   if (pool->seed0.cnt != 0) {
      ae_vector_set_length(&pool->seed0, 0);
   }
   if (pool->seedn.cnt != n) {
      ae_vector_set_length(&pool->seedn, n);
   }
   ae_shared_pool_set_seed(&pool->sourcepool, &pool->seedn, sizeof(pool->seedn), BVector_init, BVector_copy, BVector_free);
   ae_shared_pool_set_seed(&pool->temporarypool, &pool->seed0, sizeof(pool->seed0), BVector_init, BVector_copy, BVector_free);
}

// Thread-safe retrieval of array from the nbPool. If there are enough arrays
// in the pool, it is performed without additional dynamic allocations.
//
// Inputs:
//     Pool        -   nbPool properly initialized with nbPoolInit
//     A           -   array[0], must have exactly zero length (exception will
//                     be generated if length is different from zero)
//
// Outputs:
//     A           -   array[N], contents undefined
void nbpoolretrieve(nbpool *pool, BVector *a) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(BVector, tmp);
   ae_assert(a->cnt == 0, "nbPoolRetrieve: A has non-zero length on entry");
   if (pool->n == 0) {
      ae_frame_leave();
      return;
   }
   ae_shared_pool_retrieve(&pool->sourcepool, &_tmp);
   ae_swap_vectors(tmp, a);
   ae_shared_pool_recycle(&pool->temporarypool, &_tmp);
   pool->temporariescount++;
   if (pool->temporariescount > apserv_maxtemporariesinnpool) {
      pool->temporariescount = 0;
      ae_shared_pool_clear_recycled(&pool->temporarypool, true);
   }
   ae_frame_leave();
}

// Thread-safe recycling of N-length array to the nbPool.
//
// Inputs:
//     Pool        -   nbPool properly initialized with nbPoolInit
//     A           -   array[N], length must be N exactly (exception will
//                     be generated if length is different from N)
//
// Outputs:
//     A           -   array[0], length is exactly zero on exit
void nbpoolrecycle(nbpool *pool, BVector *a) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(BVector, tmp);
   ae_assert(a->cnt == pool->n, "nbPoolRecycle: A has length != N on entry");
   if (pool->n == 0) {
      ae_frame_leave();
      return;
   }
   ae_shared_pool_retrieve(&pool->temporarypool, &_tmp);
   ae_swap_vectors(tmp, a);
   ae_shared_pool_recycle(&pool->sourcepool, &_tmp);
   pool->temporariescount--;
   if (pool->temporariescount < 0) {
      pool->temporariescount = 0;
   }
   ae_frame_leave();
}

// Initialize niPool - prepare it to store N-length arrays, N >= 0.
// Tries to reuse previously allocated memory as much as possible.
void nipoolinit(nipool *pool, ae_int_t n) {
   ae_assert(n >= 0, "niPoolInit: N < 0");
   pool->n = n;
   pool->temporariescount = 0;
   if (n == 0) {
      return;
   }
   if (pool->seed0.cnt != 0) {
      ae_vector_set_length(&pool->seed0, 0);
   }
   if (pool->seedn.cnt != n) {
      ae_vector_set_length(&pool->seedn, n);
   }
   ae_shared_pool_set_seed(&pool->sourcepool, &pool->seedn, sizeof(pool->seedn), ZVector_init, ZVector_copy, ZVector_free);
   ae_shared_pool_set_seed(&pool->temporarypool, &pool->seed0, sizeof(pool->seed0), ZVector_init, ZVector_copy, ZVector_free);
}

// Thread-safe retrieval of array from the nrPool. If there are enough arrays
// in the pool, it is performed without additional dynamic allocations.
//
// Inputs:
//     Pool        -   niPool properly initialized with niPoolInit
//     A           -   array[0], must have exactly zero length (exception will
//                     be generated if length is different from zero)
//
// Outputs:
//     A           -   array[N], contents undefined
void nipoolretrieve(nipool *pool, ZVector *a) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(ZVector, tmp);
   ae_assert(a->cnt == 0, "niPoolRetrieve: A has non-zero length on entry");
   if (pool->n == 0) {
      ae_frame_leave();
      return;
   }
   ae_shared_pool_retrieve(&pool->sourcepool, &_tmp);
   ae_swap_vectors(tmp, a);
   ae_shared_pool_recycle(&pool->temporarypool, &_tmp);
   pool->temporariescount++;
   if (pool->temporariescount > apserv_maxtemporariesinnpool) {
      pool->temporariescount = 0;
      ae_shared_pool_clear_recycled(&pool->temporarypool, true);
   }
   ae_frame_leave();
}

// Thread-safe recycling of N-length array to the niPool.
//
// Inputs:
//     Pool        -   niPool properly initialized with niPoolInit
//     A           -   array[N], length must be N exactly (exception will
//                     be generated if length is different from N)
//
// Outputs:
//     A           -   array[0], length is exactly zero on exit
void nipoolrecycle(nipool *pool, ZVector *a) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(ZVector, tmp);
   ae_assert(a->cnt == pool->n, "niPoolRecycle: A has length != N on entry");
   if (pool->n == 0) {
      ae_frame_leave();
      return;
   }
   ae_shared_pool_retrieve(&pool->temporarypool, &_tmp);
   ae_swap_vectors(tmp, a);
   ae_shared_pool_recycle(&pool->sourcepool, &_tmp);
   pool->temporariescount--;
   if (pool->temporariescount < 0) {
      pool->temporariescount = 0;
   }
   ae_frame_leave();
}

// Initialize nrPool - prepare it to store N-length arrays, N >= 0.
// Tries to reuse previously allocated memory as much as possible.
void nrpoolinit(nrpool *pool, ae_int_t n) {
   ae_assert(n >= 0, "nrPoolInit: N < 0");
   pool->n = n;
   pool->temporariescount = 0;
   if (n == 0) {
      return;
   }
   if (pool->seed0.cnt != 0) {
      ae_vector_set_length(&pool->seed0, 0);
   }
   if (pool->seedn.cnt != n) {
      ae_vector_set_length(&pool->seedn, n);
   }
   ae_shared_pool_set_seed(&pool->sourcepool, &pool->seedn, sizeof(pool->seedn), RVector_init, RVector_copy, RVector_free);
   ae_shared_pool_set_seed(&pool->temporarypool, &pool->seed0, sizeof(pool->seed0), RVector_init, RVector_copy, RVector_free);
}

// Thread-safe retrieval of array from the nrPool. If there are enough arrays
// in the pool, it is performed without additional dynamic allocations.
//
// Inputs:
//     Pool        -   nrPool properly initialized with nrPoolInit
//     A           -   array[0], must have exactly zero length (exception will
//                     be generated if length is different from zero)
//
// Outputs:
//     A           -   array[N], contents undefined
void nrpoolretrieve(nrpool *pool, RVector *a) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(RVector, tmp);
   ae_assert(a->cnt == 0, "nrPoolRetrieve: A has non-zero length on entry");
   if (pool->n == 0) {
      ae_frame_leave();
      return;
   }
   ae_shared_pool_retrieve(&pool->sourcepool, &_tmp);
   ae_swap_vectors(tmp, a);
   ae_shared_pool_recycle(&pool->temporarypool, &_tmp);
   pool->temporariescount++;
   if (pool->temporariescount > apserv_maxtemporariesinnpool) {
      pool->temporariescount = 0;
      ae_shared_pool_clear_recycled(&pool->temporarypool, true);
   }
   ae_frame_leave();
}

// Thread-safe recycling of N-length array to the nrPool.
//
// Inputs:
//     Pool        -   nrPool properly initialized with nrPoolInit
//     A           -   array[N], length must be N exactly (exception will
//                     be generated if length is different from N)
//
// Outputs:
//     A           -   array[0], length is exactly zero on exit
void nrpoolrecycle(nrpool *pool, RVector *a) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(RVector, tmp);
   ae_assert(a->cnt == pool->n, "nrPoolRecycle: A has length != N on entry");
   if (pool->n == 0) {
      ae_frame_leave();
      return;
   }
   ae_shared_pool_retrieve(&pool->temporarypool, &_tmp);
   ae_swap_vectors(tmp, a);
   ae_shared_pool_recycle(&pool->sourcepool, &_tmp);
   pool->temporariescount--;
   if (pool->temporariescount < 0) {
      pool->temporariescount = 0;
   }
   ae_frame_leave();
}

// This function is used to calculate number of chunks (including partial,
// non-complete chunks) in some set. It expects that ChunkSize >= 1, TaskSize >= 0.
// Assertion is thrown otherwise.
//
// Function result is equivalent to ceil(TaskSize/ChunkSize), but with guarantees
// that rounding errors won't ruin results.
// ALGLIB: Copyright 21.01.2015 by Sergey Bochkanov
ae_int_t chunkscount(ae_int_t tasksize, ae_int_t chunksize) {
   ae_int_t result;
   ae_assert(tasksize >= 0, "chunkscount: tasksize < 0");
   ae_assert(chunksize >= 1, "chunkscount: chunksize < 1");
   result = tasksize / chunksize;
   if (tasksize % chunksize != 0) {
      result++;
   }
   return result;
}

// This function is used in parallel functions for recurrent division of large
// task into two smaller tasks.
//
// Return: Task0, with the following properties:
// * it works only for TaskSize >= 2 and TaskSize > TileSize (assertion is thrown otherwise)
// * 0 < Task0 < TaskSize
// * Task0 ~ TaskSize/2
// * Task0 >= TaskSize/2
// * Task0 is always divisible by TileSize
// ALGLIB: Copyright 07.04.2013 by Sergey Bochkanov
ae_int_t tiledsplit(ae_int_t tasksize, ae_int_t tilesize) {
   ae_assert(tasksize >= 2, "tiledsplit: tasksize < 2");
   ae_assert(tasksize > tilesize, "tiledsplit: tasksize <= tilesize");
   ae_int_t cc = chunkscount(tasksize, tilesize);
   ae_assert(cc >= 2, "tiledsplit: integrity check failed");
   ae_int_t task0 = (cc + 1) / 2 * tilesize;
   ae_assert(task0 > 0, "tiledsplit: internal error");
   ae_assert(task0 < tasksize, "tiledsplit: internal error");
   ae_assert(task0 % tilesize == 0, "tiledsplit: internal error");
   ae_assert(task0 >= tasksize - task0, "tiledsplit: internal error");
   return task0;
}

// --- OBSOLETE FUNCTION, USE TILED SPLIT INSTEAD ---
//
// This function is used in parallel functions for recurrent division of large
// task into two smaller tasks.
//
// Return: Task0, with the following properties:
// * it works only for TaskSize >= 2 and ChunkSize >= 2 (assertion is thrown otherwise)
// * 0 < Task0 < TaskSize
// * Task0 ~ TaskSize/2
// * in case TaskSize > ChunkSize, Task0 is always divisible by ChunkSize
// ALGLIB: Copyright 07.04.2013 by Sergey Bochkanov
ae_int_t splitlength(ae_int_t tasksize, ae_int_t chunksize) {
   ae_assert(chunksize >= 2, "splitlength: chunksize < 2");
   ae_assert(tasksize >= 2, "splitlength: tasksize < 2");
   ae_int_t task0 = tasksize / 2;
   if (task0 > chunksize && task0 % chunksize != 0) task0 -= task0 % chunksize;
   ae_assert(task0 > 0, "splitlength: internal error");
   ae_assert(task0 < tasksize, "splitlength: internal error");
   return task0;
}

// This function is used in parallel functions for recurrent division of large
// task into two smaller tasks.
//
// Return: Task0, with the following properties:
// * it works only for TaskSize >= 2 (assertion is thrown otherwise)
// * for TaskSize == 2, it returns Task0 = 1
// * in case TaskSize is odd, Task0 = TaskSize - 1
// * in case TaskSize is even, Task0 ~ TaskSize/2 and Task0 is even, Task0 >= TaskSize/2
// ALGLIB: Copyright 07.04.2013 by Sergey Bochkanov
ae_int_t splitlengtheven(ae_int_t tasksize) {
   ae_assert(tasksize >= 2, "splitlengtheven: tasksize < 2");
   if (tasksize == 2) return 1;
   ae_int_t task0 = tasksize;
// For odd task sizes, split the trailing odd part from it.
   if (tasksize % 2 != 0) task0--;
// For even task sizes, divide evenly.
   else {
      task0 /= 2;
      if (task0 % 2 != 0) task0++;
   }
   ae_assert(task0 > 0, "splitlengtheven: internal error");
   ae_assert(task0 < tasksize, "splitlengtheven: internal error");
   return task0;
}

// This function searches integer array. Elements in this array are actually
// records, each NRec elements wide. Each record has unique header - NHeader
// integer values, which identify it. Records are lexicographically sorted by
// header.
//
// Records are identified by their index, not offset (offset = NRec*index).
//
// This function searches A (records with indices [I0,I1)) for a record with
// header B. It returns index of this record (not offset!), or -1 on failure.
// ALGLIB: Copyright 28.03.2011 by Sergey Bochkanov
ae_int_t recsearch(ZVector *a, ae_int_t nrec, ae_int_t nheader, ae_int_t i0, ae_int_t i1, ZVector *b) {
   ae_int_t mididx;
   ae_int_t cflag;
   ae_int_t k;
   ae_int_t offs;
   ae_int_t result;
   result = -1;
   while (true) {
      if (i0 >= i1) {
         break;
      }
      mididx = (i0 + i1) / 2;
      offs = nrec * mididx;
      cflag = 0;
      for (k = 0; k < nheader; k++) {
         if (a->xZ[offs + k] < b->xZ[k]) {
            cflag = -1;
            break;
         }
         if (a->xZ[offs + k] > b->xZ[k]) {
            cflag = 1;
            break;
         }
      }
      if (cflag == 0) {
         result = mididx;
         return result;
      }
      if (cflag < 0) {
         i0 = mididx + 1;
      } else {
         i1 = mididx;
      }
   }
   return result;
}

// Returns maximum density for level 2 sparse/dense functions. Density values
// below one returned by this function are better to handle via sparse Level 2
// functionality.
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
double sparselevel2density() {
   double result;
   result = 0.1;
   return result;
}

// Returns A-tile size for a matrix.
//
// A-tiles are smallest tiles (32x32), suitable for processing by ALGLIB  own
// implementation of Level 3 linear algebra.
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
ae_int_t matrixtilesizea() {
   ae_int_t result;
   result = 32;
   return result;
}

// Returns B-tile size for a matrix.
//
// B-tiles are larger  tiles (64x64), suitable for parallel execution or for
// processing by vendor's implementation of Level 3 linear algebra.
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
ae_int_t matrixtilesizeb() {
#ifndef ALGLIB_INTERCEPTS_MKL
   ae_int_t result;
   result = 64;
   return result;
#else
   return _ialglib_i_matrixtilesizeb();
#endif
}

// This function returns minimum cost of task which is feasible for
// multithreaded processing. It returns real number in order to avoid overflow
// problems.
// ALGLIB: Copyright 10.01.2018 by Sergey Bochkanov
double smpactivationlevel() {
   double nn;
   double result;
   nn = 2.0 * matrixtilesizeb();
   result = rmax2(0.95 * 2.0 * nn * nn * nn, 10000000.0);
   return result;
}

// This function returns minimum cost of task which is feasible for
// spawn (given that multithreading is active).
//
// It returns real number in order to avoid overflow problems.
// ALGLIB: Copyright 10.01.2018 by Sergey Bochkanov
double spawnlevel() {
   double nn;
   double result;
   nn = 2.0 * matrixtilesizea();
   result = 0.95 * 2.0 * nn * nn * nn;
   return result;
}

// Minimum speedup feasible for multithreading
double minspeedup() {
   double result;
   result = 1.5;
   return result;
}

// Initialize SAvgCounter
//
// Prior value is a value that is returned when no values are in the buffer
void savgcounterinit(savgcounter *c, double priorvalue) {
   c->rsum = 0.0;
   c->rcnt = 0.0;
   c->prior = priorvalue;
}

// Enqueue value into SAvgCounter
void savgcounterenqueue(savgcounter *c, double v) {
   c->rsum += v;
   c->rcnt++;
}

// Enqueue value into SAvgCounter
double savgcounterget(savgcounter *c) {
   double result;
   if (c->rcnt == 0.0) {
      result = c->prior;
   } else {
      result = c->rsum / c->rcnt;
   }
   return result;
}

void apbuffers_init(void *_p, bool make_automatic) {
   apbuffers *p = (apbuffers *)_p;
   ae_vector_init(&p->ba0, 0, DT_BOOL, make_automatic);
   ae_vector_init(&p->ia0, 0, DT_INT, make_automatic);
   ae_vector_init(&p->ia1, 0, DT_INT, make_automatic);
   ae_vector_init(&p->ia2, 0, DT_INT, make_automatic);
   ae_vector_init(&p->ia3, 0, DT_INT, make_automatic);
   ae_vector_init(&p->ra0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->ra1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->ra2, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->ra3, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->rm0, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->rm1, 0, 0, DT_REAL, make_automatic);
}

void apbuffers_copy(void *_dst, const void *_src, bool make_automatic) {
   apbuffers *dst = (apbuffers *)_dst;
   const apbuffers *src = (const apbuffers *)_src;
   ae_vector_copy(&dst->ba0, &src->ba0, make_automatic);
   ae_vector_copy(&dst->ia0, &src->ia0, make_automatic);
   ae_vector_copy(&dst->ia1, &src->ia1, make_automatic);
   ae_vector_copy(&dst->ia2, &src->ia2, make_automatic);
   ae_vector_copy(&dst->ia3, &src->ia3, make_automatic);
   ae_vector_copy(&dst->ra0, &src->ra0, make_automatic);
   ae_vector_copy(&dst->ra1, &src->ra1, make_automatic);
   ae_vector_copy(&dst->ra2, &src->ra2, make_automatic);
   ae_vector_copy(&dst->ra3, &src->ra3, make_automatic);
   ae_matrix_copy(&dst->rm0, &src->rm0, make_automatic);
   ae_matrix_copy(&dst->rm1, &src->rm1, make_automatic);
}

void apbuffers_free(void *_p, bool make_automatic) {
   apbuffers *p = (apbuffers *)_p;
   ae_vector_free(&p->ba0, make_automatic);
   ae_vector_free(&p->ia0, make_automatic);
   ae_vector_free(&p->ia1, make_automatic);
   ae_vector_free(&p->ia2, make_automatic);
   ae_vector_free(&p->ia3, make_automatic);
   ae_vector_free(&p->ra0, make_automatic);
   ae_vector_free(&p->ra1, make_automatic);
   ae_vector_free(&p->ra2, make_automatic);
   ae_vector_free(&p->ra3, make_automatic);
   ae_matrix_free(&p->rm0, make_automatic);
   ae_matrix_free(&p->rm1, make_automatic);
}

void BVector_init(void *_p, bool make_automatic) {
   BVector *p = (BVector *)_p;
   ae_vector_init(p, 0, DT_BOOL, make_automatic);
}

void BVector_copy(void *_dst, const void *_src, bool make_automatic) {
   BVector *dst = (BVector *)_dst;
   const BVector *src = (const BVector *)_src;
   ae_vector_copy(dst, src, make_automatic);
}

void BVector_free(void *_p, bool make_automatic) {
   BVector *p = (BVector *)_p;
   ae_vector_free(p, make_automatic);
}

void ZVector_init(void *_p, bool make_automatic) {
   ZVector *p = (ZVector *)_p;
   ae_vector_init(p, 0, DT_INT, make_automatic);
}

void ZVector_copy(void *_dst, const void *_src, bool make_automatic) {
   ZVector *dst = (ZVector *)_dst;
   const ZVector *src = (const ZVector *)_src;
   ae_vector_copy(dst, src, make_automatic);
}

void ZVector_free(void *_p, bool make_automatic) {
   ZVector *p = (ZVector *)_p;
   ae_vector_free(p, make_automatic);
}

void RVector_init(void *_p, bool make_automatic) {
   RVector *p = (RVector *)_p;
   ae_vector_init(p, 0, DT_REAL, make_automatic);
}

void RVector_copy(void *_dst, const void *_src, bool make_automatic) {
   RVector *dst = (RVector *)_dst;
   const RVector *src = (const RVector *)_src;
   ae_vector_copy(dst, src, make_automatic);
}

void RVector_free(void *_p, bool make_automatic) {
   RVector *p = (RVector *)_p;
   ae_vector_free(p, make_automatic);
}

void nbpool_init(void *_p, bool make_automatic) {
   nbpool *p = (nbpool *)_p;
   ae_shared_pool_init(&p->sourcepool, make_automatic);
   ae_shared_pool_init(&p->temporarypool, make_automatic);
   BVector_init(&p->seed0, make_automatic);
   BVector_init(&p->seedn, make_automatic);
}

void nbpool_copy(void *_dst, const void *_src, bool make_automatic) {
   nbpool *dst = (nbpool *)_dst;
   const nbpool *src = (const nbpool *)_src;
   dst->n = src->n;
   dst->temporariescount = src->temporariescount;
   ae_shared_pool_copy(&dst->sourcepool, &src->sourcepool, make_automatic);
   ae_shared_pool_copy(&dst->temporarypool, &src->temporarypool, make_automatic);
   BVector_copy(&dst->seed0, &src->seed0, make_automatic);
   BVector_copy(&dst->seedn, &src->seedn, make_automatic);
}

void nbpool_free(void *_p, bool make_automatic) {
   nbpool *p = (nbpool *)_p;
   ae_shared_pool_free(&p->sourcepool, make_automatic);
   ae_shared_pool_free(&p->temporarypool, make_automatic);
   BVector_free(&p->seed0, make_automatic);
   BVector_free(&p->seedn, make_automatic);
}

void nipool_init(void *_p, bool make_automatic) {
   nipool *p = (nipool *)_p;
   ae_shared_pool_init(&p->sourcepool, make_automatic);
   ae_shared_pool_init(&p->temporarypool, make_automatic);
   ZVector_init(&p->seed0, make_automatic);
   ZVector_init(&p->seedn, make_automatic);
}

void nipool_copy(void *_dst, const void *_src, bool make_automatic) {
   nipool *dst = (nipool *)_dst;
   const nipool *src = (const nipool *)_src;
   dst->n = src->n;
   dst->temporariescount = src->temporariescount;
   ae_shared_pool_copy(&dst->sourcepool, &src->sourcepool, make_automatic);
   ae_shared_pool_copy(&dst->temporarypool, &src->temporarypool, make_automatic);
   ZVector_copy(&dst->seed0, &src->seed0, make_automatic);
   ZVector_copy(&dst->seedn, &src->seedn, make_automatic);
}

void nipool_free(void *_p, bool make_automatic) {
   nipool *p = (nipool *)_p;
   ae_shared_pool_free(&p->sourcepool, make_automatic);
   ae_shared_pool_free(&p->temporarypool, make_automatic);
   ZVector_free(&p->seed0, make_automatic);
   ZVector_free(&p->seedn, make_automatic);
}

void nrpool_init(void *_p, bool make_automatic) {
   nrpool *p = (nrpool *)_p;
   ae_shared_pool_init(&p->sourcepool, make_automatic);
   ae_shared_pool_init(&p->temporarypool, make_automatic);
   RVector_init(&p->seed0, make_automatic);
   RVector_init(&p->seedn, make_automatic);
}

void nrpool_copy(void *_dst, const void *_src, bool make_automatic) {
   nrpool *dst = (nrpool *)_dst;
   const nrpool *src = (const nrpool *)_src;
   dst->n = src->n;
   dst->temporariescount = src->temporariescount;
   ae_shared_pool_copy(&dst->sourcepool, &src->sourcepool, make_automatic);
   ae_shared_pool_copy(&dst->temporarypool, &src->temporarypool, make_automatic);
   RVector_copy(&dst->seed0, &src->seed0, make_automatic);
   RVector_copy(&dst->seedn, &src->seedn, make_automatic);
}

void nrpool_free(void *_p, bool make_automatic) {
   nrpool *p = (nrpool *)_p;
   ae_shared_pool_free(&p->sourcepool, make_automatic);
   ae_shared_pool_free(&p->temporarypool, make_automatic);
   RVector_free(&p->seed0, make_automatic);
   RVector_free(&p->seedn, make_automatic);
}

void savgcounter_init(void *_p, bool make_automatic) {
}

void savgcounter_copy(void *_dst, const void *_src, bool make_automatic) {
   savgcounter *dst = (savgcounter *)_dst;
   const savgcounter *src = (const savgcounter *)_src;
   dst->rsum = src->rsum;
   dst->rcnt = src->rcnt;
   dst->prior = src->prior;
}

void savgcounter_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

// === ABLASF Package ===
namespace alglib_impl {
// The dot product of elements [0,n) of vectors x and y.
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to process.
//	y:	An n-vector to process.
// Result:
//	The inner product (x,y).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotv(ae_int_t n, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
// Use KerSub* for a kernel that does not return a result.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2Fma(rdotv(n, x->xR, y->xR))
#endif
// Original generic C implementation.
   double dot = 0.0;
   for (ae_int_t i = 0; i < n; i++) dot += x->xR[i] * y->xR[i];
   return dot;
}

// The dot product of elements [0,n) of vectors x and y[iy,_].
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to process.
//	y:	A ? x n matrix to process.
//	iy:	A row index in y.
// Result:
//	The inner product (x,y[iy,_]).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotvr(ae_int_t n, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2Fma(rdotv(n, x->xR, y->xyR[iy]))
#endif
   double dot = 0.0;
   for (ae_int_t j = 0; j < n; j++) dot += x->xR[j] * y->xyR[iy][j];
   return dot;
}

// The dot product of elements [0,n) of vectors x[ix,_] and y[iy,_].
// Inputs:
//	n:	The vector length.
//	x:	An n x n matrix to process.
//	ix:	A row index in x.
//	y:	A ? x n matrix to process.
//	iy:	A row index in y.
// Result:
//	The inner product (x[ix,_],y[iy,_]).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotrr(ae_int_t n, RMatrix *x, ae_int_t ix, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2Fma(rdotv(n, x->xyR[ix], y->xyR[iy]))
#endif
   double dot = 0.0;
   for (ae_int_t j = 0; j < n; j++) dot += x->xyR[ix][j] * y->xyR[iy][j];
   return dot;
}

// The dot product square of elements [0,n) of the vector x.
// Inputs:
//	n:	The vector length.
//	x:	The n-vector to process.
// Result:
//	The inner product square (x,x).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotv2(ae_int_t n, RVector *x) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2Fma(rdotv2(n, x->xR))
#endif
   double dot = 0.0;
   for (ae_int_t i = 0; i < n; i++) dot += sqr(x->xR[i]);
   return dot;
}

// Add elements [0,n) of the vector x to y.
// Inputs:
//	n:	The vector length.
//	alpha:	The multiplier.
//	x:	An n-vector to process.
//	y:	An n-vector to process.
// Result:
//	y += alpha x.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddv(ae_int_t n, double alpha, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2Fma(raddv(n, alpha, x->xR, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] += alpha * x->xR[i];
}

// Add elements [0,n) of the vector x to y[iy,_].
// Inputs:
//	n:	The vector length.
//	alpha:	The multiplier.
//	x:	An n-vector to add.
//	y:	The destination matrix.
//	iy:	The destination row index in y.
// Result:
//	y[iy,_] += alpha x.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddvr(ae_int_t n, double alpha, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2Fma(raddv(n, alpha, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] += alpha * x->xR[i];
}

// Add elements [0,n) of the vector x[ix,_] to y.
// Inputs:
//	n:	The vector length.
//	alpha:	The multiplier.
//	x:	The source ? x n matrix.
//	ix:	A row index in x.
//	y:	The n-vector to process.
// Result:
//	y += alpha x[ix,_].
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddrv(ae_int_t n, double alpha, RMatrix *x, ae_int_t ix, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2Fma(raddv(n, alpha, x->xyR[ix], y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] += alpha * x->xyR[ix][i];
}

// Add elements [0,n) of the vector x[ix,_] to y[iy,_].
// Inputs:
//	n:	The vector length.
//	alpha:	A multiplier.
//	x:	The source ? x n matrix.
//	ix:	The source row index in x.
//	y:	The destination ? x n matrix.
//	iy:	The destination row index in y.
// Result:
//	y[iy,_] += alpha x[ix,_].
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddrr(ae_int_t n, double alpha, RMatrix *x, ae_int_t ix, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2Fma(raddv(n, alpha, x->xyR[ix], y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] += alpha * x->xyR[ix][i];
}

// Add elements [x0,x0+n) of the vector x to elements [y0,y0+n) of y.
// Inputs:
//	n:	The vector length.
//	alpha:	The multiplier.
//	x:	The source vector.
//	x0:	The source offset in x.
//	y:	The destination vector.
//	y0:	The destination offset in y.
// Result:
//	y[y0+i] += alpha x[x0+i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddvx(ae_int_t n, double alpha, RVector *x, ae_int_t x0, RVector *y, ae_int_t y0) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2Fma(raddvx(n, alpha, x->xR + x0, y->xR + y0))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[y0 + i] += alpha * x->xR[x0 + i];
}

// Add elements [0,n) of the vector x to y[_,jy].
// Inputs:
//	n:	The vector length.
//	alpha:	The multiplier.
//	x:	The vector to add.
//	y:	The target matrix.
//	jy:	The target column index in y.
// Result:
//	y[i,jy] += alpha x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddvc(ae_int_t n, double alpha, RVector *x, RMatrix *y, ae_int_t jy) {
   for (ae_int_t i = 0; i < n; i++) y->xyR[i][jy] += alpha * x->xR[i];
}

// Component-wise multiply elements [0,n) of the vector y by x.
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to multiply by.
//	y:	The target vector.
// Result:
//	y[i] *= x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemulv(ae_int_t n, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergemulv(n, x->xR, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] *= x->xR[i];
}

// Component-wise multiply elements [0,n) of the vector y[iy,_] by x.
// Inputs:
//	n:	The vector length.
//	x:	A vector to multiply by.
//	y:	The target matrix.
//	iy:	The target row index in y.
// Result:
//	y[iy,i] *= x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemulvr(ae_int_t n, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergemulv(n, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] *= x->xR[i];
}

// Component-wise multiply elements [0,n) of the vector y by x[ix,_].
// Inputs:
//	n:	The vector length.
//	x:	A matrix to multiply from.
//	ix:	The row index in x.
//	y:	The target vector.
// Result:
//	y[i] *= x[ix,i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemulrv(ae_int_t n, RMatrix *x, ae_int_t ix, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergemulv(n, x->xyR[ix], y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] *= x->xyR[ix][i];
}

// Component-wise divide elements [0,n) of vector y by vector x.
// Inputs:
//	n:	The vector length.
//	x:	The vector to divide by.
//	y:	The target vector.
// Result:
//	y[i] /= x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergedivv(ae_int_t n, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubAvx2(rmergedivv(n, x->xR, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] /= x->xR[i];
}

// Component-wise divide elements [0,n) of the vector y[iy,_] by x.
// Inputs:
//	n:	The vector length.
//	x:	The vector divide by.
//	y:	The target matrix.
//	iy:	The row index in y.
// Result:
//	y[iy,i] /= x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergedivvr(ae_int_t n, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubAvx2(rmergedivv(n, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] /= x->xR[i];
}

// Component-wise divide elements [0,n) of vector x[ix,_] by y.
// Inputs:
//	n:	The vector length.
//	x:	The source matrix.
//	ix:	The row index in x.
//	y:	The target vector.
// Result:
//	y[i] /= x[ix,i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergedivrv(ae_int_t n, RMatrix *x, ae_int_t ix, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubAvx2(rmergedivv(n, x->xyR[ix], y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] /= x->xyR[ix][i];
}

// Component-wise max elements [0,n) of the vector x into y.
// Inputs:
//	n:	The vector length.
//	x:	The vector to max from.
//	y:	The target vector.
// Result:
//	y[i] = max(y[i],x[i]), for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemaxv(ae_int_t n, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergemaxv(n, x->xR, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = rmax2(y->xR[i], x->xR[i]);
}

// Component-wise max elements [0,n) of the vector x into y[iy,_].
// Inputs:
//	n:	The vector length.
//	x:	The vector to max from.
//	y:	The matrix to max into.
//	iy:	A row index in y.
// Result:
//	y[iy,i] = max(y[iy,i],x[i]), for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemaxvr(ae_int_t n, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergemaxv(n, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] = rmax2(y->xyR[iy][i], x->xR[i]);
}

// Component-wise max elements [0,n) of the vector x[ix,_] into y.
// Inputs:
//	n:	The vector length.
//	x:	The source matrix.
//	ix:	A row index in x.
//	y:	The destination vector.
// Result:
//	y[i] = max(y[i],x[ix,i]), for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemaxrv(ae_int_t n, RMatrix *x, ae_int_t ix, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergemaxv(n, x->xyR[ix], y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = rmax2(y->xR[i], x->xyR[ix][i]);
}

// Component-wise min elements [0,n) of the vector x into y.
// Inputs:
//	n:	The vector length.
//	x:	The source vector.
//	y:	The target vector.
// Result:
//	y[i] = min(y[i],x[i]), for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergeminv(ae_int_t n, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergeminv(n, x->xR, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = rmin2(y->xR[i], x->xR[i]);
}

// Component-wise min elements [0,n) of the vector x into y[iy,_].
// Inputs:
//	n:	The vector length.
//	x:	The vector to min from.
//	y:	The target matrix.
//	iy:	A row index in y.
// Result:
//	y[iy,i] = min(y[iy,i],x[i]), for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergeminvr(ae_int_t n, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergeminv(n, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] = rmin2(y->xyR[iy][i], x->xR[i]);
}

// Component-wise min elements [0,n) of the vector x[ix,_] into y.
// Inputs:
//	n:	The vector length.
//	x:	The source matrix.
//	ix:	A row index in x.
//	y:	The target vector.
// Result:
//	y[i] = min(y[i],x[ix,i]), for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergeminrv(ae_int_t n, RMatrix *x, ae_int_t ix, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmergeminv(n, x->xyR[ix], y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = rmin2(y->xR[i], x->xyR[ix][i]);
}

// Square root elements [0,n) of vector y.
// Inputs:
//	n:	The vector length.
//	y:	An n-vector to process.
// Output:
//	y:	Elements [0,n) of y square root'ed.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsqrtv(ae_int_t n, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubAvx2(rsqrtv(n, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = sqrt(y->xR[i]);
}

// Square root elements [0,n) of vector y[iy,_].
// Inputs:
//	n:	The vector length.
//	y:	A ? y n matrix to process.
//	iy:	A row index in y.
// Output:
//	y:	Elements [0,n) of y[iy,_] square root'ed.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsqrtr(ae_int_t n, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubAvx2(rsqrtv(n, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] = sqrt(y->xyR[iy][i]);
}

// Multiply elements [0,n) of the vector y by the scalar v.
// Inputs:
//	n:	The vector length.
//	v:	The multiplier.
//	y:	The n-vector to process.
// Output:
//	y:	Elements [0,n) of y *= v.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmulv(ae_int_t n, double v, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmulv(n, v, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] *= v;
}

// Multiply elements [0,n) of the vector y[iy,_] by the scalar v.
// Inputs:
//	n:	The row length.
//	v:	The multiplier.
//	y:	The ? x n matrix to process.
//	iy:	The row of y to process.
// Output:
//	y:	Elements [0,n) of y[iy,_] *= v.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmulr(ae_int_t n, double v, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmulv(n, v, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] *= v;
}

// Multiply elements [y0,y0+n) of the vector y by the scalar v.
// Inputs:
//	n:	The sub-vector length.
//	v:	The multiplier.
//	y:	The vector to process.
//	y0:	An offset location in y.
// Output:
//	y:	y[y0+i] *= v, for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmulvx(ae_int_t n, double v, RVector *y, ae_int_t y0) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rmulvx(n, v, y->xR + y0))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[y0 + i] *= v;
}

// Component-wise multiply elements [0,n) of the vectors x and y and add into z.
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to process.
//	y:	An n-vector to process.
//	z:	An n-vector to process.
// Result:
//	z[i] += x[i] y[i], for i in [0,n).
// ALGLIB: Copyright 29.10.2021 by Sergey Bochkanov
void rmuladdv(ae_int_t n, RVector *x, RVector *y, RVector *z) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubFma(rmuladdv(n, x->xR, y->xR, z->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) z->xR[i] += x->xR[i] * y->xR[i];
}

// Component-wise multiply elements [0,n) of the vectors x and y and subtract into z.
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to process.
//	y:	An n-vector to process.
//	z:	An n-vector to process.
// Result:
//	z[i] -= x[i] y[i], for i in [0,n).
// ALGLIB: Copyright 29.10.2021 by Sergey Bochkanov
void rnegmuladdv(ae_int_t n, RVector *x, RVector *y, RVector *z) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubFma(rnegmuladdv(n, x->xR, y->xR, z->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) z->xR[i] -= x->xR[i] * y->xR[i];
}

// The maximum over elements [0,n) of the vector x.
// Inputs:
//	n:	The vector length.
//	x:	The n-vector to process.
// Output:
//	max(x[i], i in [0,n)), or 0 for n == 0.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxv(ae_int_t n, RVector *x) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2(rmaxv(n, x->xR))
#endif
   if (n <= 0) return 0.0;
   double max = x->xR[0];
   for (ae_int_t i = 1; i < n; i++) max = rmax2(max, x->xR[i]);
   return max;
}

// The maximum over elements [0,n) of the vector x[ix,_].
// Inputs:
//	n:	The vector length.
//	x:	The matrix to process.
//	ix:	A row index in x.
// Output:
//	max(x[ix,i], i in [0,n)), or 0 for n == 0.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxr(ae_int_t n, RMatrix *x, ae_int_t ix) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2(rmaxv(n, x->xyR[ix]))
#endif
   if (n <= 0) return 0.0;
   double max = x->xyR[ix][0];
   for (ae_int_t i = 1; i < n; i++) max = rmax2(max, x->xyR[ix][i]);
   return max;
}

// The maximum over elements [0,n) of the vector |x[_]|.
// Inputs:
//	n:	The vector length.
//	x:	The n-vector to process.
// Output:
//	max(|x[i]|, i in [0,n)), or 0 for n == 0.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxabsv(ae_int_t n, RVector *x) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2(rmaxabsv(n, x->xR))
#endif
   double max = 0.0;
   for (ae_int_t i = 0; i < n; i++) max = rmax2(max, fabs(x->xR[i]));
   return max;
}

// The maximum over elements [0,n) of the vector |x[ix,_]|.
// Inputs:
//	n:	The vector length.
//	x:	The matrix to process.
//	ix:	A row index in x.
// Output:
//	max(|x[ix,i]|, i in [0,n)), or 0 for n == 0.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxabsr(ae_int_t n, RMatrix *x, ae_int_t ix) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerFunSse2Avx2(rmaxabsv(n, x->xyR[ix]))
#endif
   double max = 0.0;
   for (ae_int_t i = 0; i < n; i++) max = rmax2(max, fabs(x->xyR[ix][i]));
   return max;
}

// Set elements [0,n) of the vector y to v.
// Inputs:
//	n:	The vector length.
//	v:	The value to set.
//	y:	An n-vector.
// Output:
//	y:	Elements [0,n) of y are set to v.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
// Boolean:
void bsetv(ae_int_t n, bool v, BVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= 8 * _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(bsetv(n, v, y->xB))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xB[j] = v;
}

// Integer:
void isetv(ae_int_t n, ae_int_t v, ZVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(isetv(n, v, y->xZ))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xZ[j] = v;
}

// Real:
void rsetv(ae_int_t n, double v, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rsetv(n, v, y->xR))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xR[j] = v;
}

// Complex:
void csetv(ae_int_t n, complex v, CVector *y) {
   for (ae_int_t j = 0; j < n; j++) y->xC[j] = v;
}

// Set elements [y0,y0+n) of the vector y to v.
// Inputs:
//	n:	The sub-vector length.
//	v:	The value to set.
//	y:	A vector to process.
//	y0:	The offset of y.
// Output:
//	y:	y[y0+i] = v, for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetvx(ae_int_t n, double v, RVector *y, ae_int_t y0) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rsetvx(n, v, y->xR + y0))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xR[y0 + j] = v;
}

#if !defined ALGLIB_NO_FAST_KERNELS
// Set elements [0,n) of the vector yp[_] to v.
// Inputs:
//	n:	The vector size.
//	v:	The value to set.
//	yp:	An n-vector.
// Output:
//	yp:	Elements [0,n) of yp are set to v.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
static void rsetm_simd(const ae_int_t n, const double v, double *yp) {
   KerSubSse2Avx2(rsetv(n, v, yp))
   for (ae_int_t j = 0; j < n; j++) yp[j] = v;
}
#endif

// Set elements [0,m) x [0,n) of the matrix y to v.
// Inputs:
//	m, n:	The row and column counts.
//	v:	The value to set.
//	y:	The m x n matrix.
// Output:
//	y:	Elements [0,m) x [0,n) of y are set to v.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetm(ae_int_t m, ae_int_t n, double v, RMatrix *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) {
      for (ae_int_t i = 0; i < m; i++) rsetm_simd(n, v, y->xyR[i]);
      return;
   }
#endif
   for (ae_int_t i = 0; i < m; i++) for (ae_int_t j = 0; j < n; j++) y->xyR[i][j] = v;
}

// Set elements [0,n) of the vector y[iy,_] to v.
// Inputs:
//	n:	The vector length.
//	v:	The value to set.
//	y:	A matrix of size n x n or larger.
//	iy:	A row index in y.
// Output:
//	y:	y[iy,j] = v, for j in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetr(ae_int_t n, double v, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rsetv(n, v, y->xyR[iy]))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xyR[iy][j] = v;
}

// Set elements [0,n) of the vector y[_,jy] to v.
// Inputs:
//	n:	The vector length.
//	v:	The value to set.
//	y:	A matrix of size n x n or larger.
//	jy:	A column index in y.
// Output:
//	y:	y[i,jy] = v, for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetc(ae_int_t n, double v, RMatrix *y, ae_int_t jy) {
   for (ae_int_t i = 0; i < n; i++) y->xyR[i][jy] = v;
}

// Set elements [0,n) of the vector y to v, reallocating y if it's too small.
// Inputs:
//	n:	The vector length.
//	v:	The value to set.
//	y:	A possibly-preallocated vector.
// Output:
//	y:	Elements [0,n) of y are set to v; y is reallocated, if it is shorter than n.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
// Boolean:
void bsetallocv(ae_int_t n, bool v, BVector *y) {
   vectorsetlengthatleast(y, n);
   bsetv(n, v, y);
}
// Integer:
void isetallocv(ae_int_t n, ae_int_t v, ZVector *y) {
   vectorsetlengthatleast(y, n);
   isetv(n, v, y);
}
// Real:
void rsetallocv(ae_int_t n, double v, RVector *y) {
   vectorsetlengthatleast(y, n);
   rsetv(n, v, y);
}
// Complex:
void csetallocv(ae_int_t n, complex v, CVector *y) {
   vectorsetlengthatleast(y, n);
   csetv(n, v, y);
}

// Set elements [0,m) x [0,n) of the matrix y to v, reallocating y if it's too small.
// Inputs:
//	m:	The row count.
//	n:	The column count.
//	v:	The value to set.
//	y:	A possibly-preallocated matrix.
// Output:
//	y:	Elements [0,m) x [0,n) of y are set to v; y is reallocated, if its row or column counts are too small.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetallocm(ae_int_t m, ae_int_t n, double v, RMatrix *y) {
   matrixsetlengthatleast(y, m, n);
   rsetm(m, n, v, y);
}

// Reallocate the vector x, if its length is less than n.
// Leave it unchanged, if it is large enough.
// Inputs:
//	n:	The desired vector length.
//	x:	A possibly-preallocated vector.
// Output:
//	x:	A vector, at least of length n.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void allocv(ae_int_t n, ae_vector *x) {
   vectorsetlengthatleast(x, n);
}

// Reallocate the matrix x, if its size is less than m x n.
// Leave it unchanged, if it is large enough.
// Inputs:
//	m:	The row count.
//	n:	The column count.
//	x:	A possibly-preallocated matrix.
// Output:
//	x:	A matrix, at least of size m x n.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void allocm(ae_int_t m, ae_int_t n, ae_matrix *x) {
   matrixsetlengthatleast(x, m, n);
}

// Copy elements [0,n) of the vector x into y.
// Inputs:
//	n:	The vector length.
//	x:	The source n-vector.
//	y:	A preallocated n-vector.
// Output:
//	y:	A copy of elements [0,n) of x into y.
// Note:
// *	The destination and source should NOT overlap.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
// Boolean:
void bcopyv(ae_int_t n, BVector *x, BVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= 8 * _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(bcopyv(n, x->xB, y->xB))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xB[j] = x->xB[j];
}
// Integer:
void icopyv(ae_int_t n, ZVector *x, ZVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(icopyv(n, x->xZ, y->xZ))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xZ[j] = x->xZ[j];
}
// Real:
void rcopyv(ae_int_t n, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopyv(n, x->xR, y->xR))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xR[j] = x->xR[j];
}

// Copy elements [0,n) of the vector x into y[iy,_].
// Inputs:
//	n:	The vector length.
//	x:	The source n-vector.
//	y:	A preallocated matrix large enough to store the result.
//	iy:	The destination row index in y.
// Output:
//	y:	A copy of elements [0,n) of x into y[iy,_].
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyvr(ae_int_t n, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopyv(n, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xyR[iy][j] = x->xR[j];
}

// Copy elements [0,n) of the vector x[ix,_] into y.
// Inputs:
//	n:	The vector length.
//	x:	The source matrix.
//	ix:	The source row index in x.
//	y:	A preallocated destination vector.
// Output:
//	y:	A copy of elements [0,n) of x[ix,_] into y.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyrv(ae_int_t n, RMatrix *x, ae_int_t ix, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopyv(n, x->xyR[ix], y->xR))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xR[j] = x->xyR[ix][j];
}

// Copy elements [0,n) of the vector x[ix,_] into y[iy,_].
// Inputs:
//	n:	The vector length.
//	x:	The source n x n matrix.
//	ix:	The source row index in x.
//	y:	A preallocated destination matrix.
//	iy:	The destination row index in y.
// Output:
//	y:	A copy of elements [0,n) of x[ix,_] into y[iy,_].
// Note:
// *	x[ix,_] and y[iy,_] may not overlap. //(@) The original said they may overlap.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyrr(ae_int_t n, RMatrix *x, ae_int_t ix, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopyv(n, x->xyR[ix], y->xyR[iy]))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xyR[iy][j] = x->xyR[ix][j];
}

// Copy elements [x0,x0+n) of the vector x into elements [y0,y0+n) of the vector y, extended version.
// Inputs:
//	n:	The vector length.
//	x:	The source vector.
//	x0:	The source offset in x.
//	y:	A preallocated n-vector.
//	y0:	The destination offset in y.
// Output:
//	y:	y[y0+j] = x[x0+j], for j in [0,n).
// Note:
// *	The destination and source should NOT overlap.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
// Integer:
void icopyvx(ae_int_t n, ZVector *x, ae_int_t x0, ZVector *y, ae_int_t y0) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(icopyvx(n, x->xZ + x0, y->xZ + y0))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xZ[y0 + j] = x->xZ[x0 + j];
}
// Real:
void rcopyvx(ae_int_t n, RVector *x, ae_int_t x0, RVector *y, ae_int_t y0) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopyvx(n, x->xR + x0, y->xR + y0))
#endif
   for (ae_int_t j = 0; j < n; j++) y->xR[y0 + j] = x->xR[x0 + j];
}

// Copy elements [0,n) of the vector x into y[_,jy].
// Inputs:
//	n:	The vector length.
//	x:	The source n-vector.
//	y:	A preallocated matrix large enough to store the result.
//	jy:	The destination column index in y.
// Output:
//	y:	A copy of elements [0,n) of x into y[_,jy].
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyvc(ae_int_t n, RVector *x, RMatrix *y, ae_int_t jy) {
   for (ae_int_t i = 0; i < n; i++) y->xyR[i][jy] = x->xR[i];
}

// Copy elements [0,n) of the vector x[_,jx] into y.
// Inputs:
//	n:	The vector length.
//	x:	The source matrix.
//	jx:	The source column index in x.
//	y:	A preallocated destination n-vector.
// Output:
//	y:	A copy of elements [0,n) of x[_,jx] into y.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopycv(ae_int_t n, RMatrix *x, ae_int_t jx, RVector *y) {
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = x->xyR[i][jx];
}

// Copy elements [0,m) x [0,n) of the matrix x into y, resizing y if needed.
// On resize, the dimensions of y are increased - but not decreased.
// Inputs:
//	m:	The row count.
//	n:	The column count.
//	x:	The source m x n matrix.
//	y:	A possibly-preallocated m x n matrix, resized if needed.
// Output:
//	y:	Elements [0,m) x [0,n) of x copied into y.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopym(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y) {
   if (m == 0 || n == 0) return;
   for (ae_int_t i = 0; i < m; i++) for (ae_int_t j = 0; j < n; j++) y->xyR[i][j] = x->xyR[i][j];
}

// Multiply elements [0,n) of the vector x by v into y.
// Inputs:
//	n:	The vector length.
//	v:	The multiplier.
//	x:	The source n-vector.
//	y:	A preallocated n-vector.
// Output:
//	y:	y[i] = v x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopymulv(ae_int_t n, double v, RVector *x, RVector *y) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopymulv(n, v, x->xR, y->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xR[i] = v * x->xR[i];
}

// Multiply elements [0,n) of the vector x by v into y[iy,_].
// Inputs:
//	n:	The vector length.
//	v:	The multiplier.
//	x:	The source n-vector.
//	y:	A preallocated ? x n matrix.
//	iy:	The destination row index in y.
// Output:
//	y:	y[iy,i] = v x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopymulvr(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t iy) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubSse2Avx2(rcopymulv(n, v, x->xR, y->xyR[iy]))
#endif
   for (ae_int_t i = 0; i < n; i++) y->xyR[iy][i] = v * x->xR[i];
}

// Multiply elements [0,n) of the vector x by v into y[_,jy].
// Inputs:
//	n:	The vector length.
//	v:	The multiplier.
//	x:	The source n-vector.
//	y:	A preallocated n x ? matrix.
//	jy:	The destination column index in y.
// Output:
//	y:	y[i,jy] = v x[i], for i in [0,n).
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopymulvc(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t jy) {
   for (ae_int_t i = 0; i < n; i++) y->xyR[i][jy] = v * x->xR[i];
}

// Component-wise multiply elements [0,n) of the vectors x and y adding to z into w.
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to process.
//	y:	An n-vector to process.
//	z:	An n-vector to process.
//	w:	An n-vector to process.
// Result:
//	w[i] = z[i] + x[i] y[i], for i in [0,n).
// ALGLIB: Copyright 29.10.2021 by Sergey Bochkanov
void rcopymuladdv(ae_int_t n, RVector *x, RVector *y, RVector *z, RVector *w) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubFma(rcopymuladdv(n, x->xR, y->xR, z->xR, w->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) w->xR[i] = z->xR[i] + x->xR[i] * y->xR[i];
}

// Component-wise multiply elements [0,n) of the vectors x and y subtracting from z into w.
// Inputs:
//	n:	The vector length.
//	x:	An n-vector to process.
//	y:	An n-vector to process.
//	z:	An n-vector to process.
//	w:	An n-vector to process.
// Result:
//	w[i] = z[i] - x[i] y[i], for i in [0,n).
// ALGLIB: Copyright 29.10.2021 by Sergey Bochkanov
void rcopynegmuladdv(ae_int_t n, RVector *x, RVector *y, RVector *z, RVector *w) {
#if !defined ALGLIB_NO_FAST_KERNELS
// Try fast kernels.
// On success this macro will return, on failure to find kernel it will pass execution to the generic C implementation.
   if (n >= _ABLASF_KERNEL_SIZE1) KerSubFma(rcopynegmuladdv(n, x->xR, y->xR, z->xR, w->xR))
#endif
   for (ae_int_t i = 0; i < n; i++) w->xR[i] = z->xR[i] - x->xR[i] * y->xR[i];
}

// Copy elements [0,n) of the vector x into y, resizing y if needed.
// Inputs:
//	n:	The vector length.
//	x:	The source n-vector.
//	y:	A possibly-preallocated n-vector, resized if needed.
// Output:
//	y:	Elements [0,n) of x are copied into y.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
// Boolean:
void bcopyallocv(ae_int_t n, BVector *x, BVector *y) {
   vectorsetlengthatleast(y, n);
   bcopyv(n, x, y);
}
// Integer:
void icopyallocv(ae_int_t n, ZVector *x, ZVector *y) {
   vectorsetlengthatleast(y, n);
   icopyv(n, x, y);
}
// Real:
void rcopyallocv(ae_int_t n, RVector *x, RVector *y) {
   vectorsetlengthatleast(y, n);
   rcopyv(n, x, y);
}

// Copy elements [0,m) x [0,n) of the matrix x into y, resizing y if needed.
// On resize, the dimensions of y are increased - but not decreased.
// Inputs:
//	m:	The row count.
//	n:	The column count.
//	x:	The m x n source matrix.
//	y:	A possibly-preallocated m x n matrix, resized if needed.
// Output:
//	y:	Elements [0,m) x [0,n) of x are copied into y.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyallocm(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y) {
   if (m == 0 || n == 0) return;
   matrixsetlengthatleast(y, m, n);
   rcopym(m, n, x, y);
}

// Grow the vector x, i.e. change its size in such a way that:
// *	its contents is preserved,
// *	its new size is at least newn,
// *	its actual size may be larger than newn, so that subsequent grow() calls can return without reallocation.
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
// Integer:
void igrowv(ae_int_t newn, ZVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   if (x->cnt >= newn) {
      ae_frame_leave();
      return;
   }
   ae_int_t oldn = x->cnt;
   newn = imax2(newn, iround(1.8 * oldn + 1.0));
   NewVector(oldx, 0, DT_INT);
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, newn);
   icopyv(oldn, &oldx, x);
   ae_frame_leave();
}
// Real:
void rgrowv(ae_int_t newn, RVector *x) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   if (x->cnt >= newn) {
      ae_frame_leave();
      return;
   }
   ae_int_t oldn = x->cnt;
   newn = imax2(newn, iround(1.8 * oldn + 1.0));
   NewVector(oldx, 0, DT_REAL);
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, newn);
   rcopyv(oldn, &oldx, x);
   ae_frame_leave();
}

// Matrix-vector product: y = alpha a' x + beta y.
// Note:
// *	This function expects y to be large enough to store the result.
//	No automatic preallocation happens for smaller arrays.
//	No integrity checks are performed on the sizes of a, x, y.
// Inputs:
//	m:	The number of rows of a'.
//	n:	The number of columns of a'.
//	alpha:	The coefficient.
//	a:	The source matrix.
//	opa:	The operation type:
//		0:	a' = a,
//		1:	a' = a^T.
//	x:	An input vector, at least of size n.
//	beta:	A coefficient.
//	y:	A preallocated output vector, at least of size m.
// Output:
//	y:	The vector which stores the result.
//
// Handling of Special Cases:
// --------------------------
// *	If m == 0, then the subroutine does nothing.
//	It does not even touch the arrays.
// *	If n == 0 or alpha == 0, then:
//	-	If beta == 0, then y is filled by 0's.
//		a and x are not referenced at all.
//		The initial values of y are ignored (we do not multiply y by 0, we just rewrite it by 0's).
//	-	If beta != 0, then y is multiplied by beta.
// *	If m > 0, n > 0, alpha != 0, but beta == 0, then y is replaced by a x;
//	the initial state of y is ignored (rewritten by a x, without initial multiplication by 0's).
// ALGLIB Routine: Copyright 01.09.2021 by Sergey Bochkanov
void rgemv(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t opa, RVector *x, double beta, RVector *y) {
// Properly premultiply y by beta.
// Quick exit for m == 0, n == 0 or alpha == 0.
// After this block we have m > 0, n > 0, alpha != 0.
   if (m <= 0) return;
   if (beta != 0.0) rmulv(m, beta, y); else rsetv(m, 0.0, y);
   if (n <= 0 || alpha == 0.0) return;
// Straight or transposed?
   switch (opa) {
   // y += a x.
      case 0:
#if !defined ALGLIB_NO_FAST_KERNELS
      // Try SIMD code.
         if (n >= _ABLASF_KERNEL_SIZE2) KerSubAvx2Fma(rgemv_straight(m, n, alpha, a, x->xR, y->xR))
#endif
      // Generic C version.
         for (ae_int_t i = 0; i < m; i++) {
            double v = 0.0;
            for (ae_int_t j = 0; j < n; j++) v += a->xyR[i][j] * x->xR[j];
            y->xR[i] += alpha * v;
         }
      break;
   // y += a^T x.
      case 1:
#if !defined ALGLIB_NO_FAST_KERNELS
      // Try SIMD code.
         if (m >= _ABLASF_KERNEL_SIZE2) KerSubAvx2Fma(rgemv_transposed(m, n, alpha, a, x->xR, y->xR))
#endif
      // Generic C version.
         for (ae_int_t i = 0; i < n; i++) {
            double v = alpha * x->xR[i];
            for (ae_int_t j = 0; j < m; j++) y->xR[j] += v * a->xyR[i][j];
         }
      break;
   }
}

// Matrix-vector product: y = alpha a' x + beta y.
// Here x, y, a are sub-vectors/sub-matrices of larger vectors/matrices.
// Note:
// *	This function expects y to be large enough to store the result.
//	No automatic preallocation happens for smaller arrays.
//	No integrity checks are performed on the sizes of a, x, y.
// Inputs:
//	m:	The number of rows of a'.
//	n:	The number of columns of a'.
//	alpha:	A coefficient.
//	a:	The source matrix.
//	ia:	The sub-matrix offset (row index).
//	ja:	The sub-matrix offset (column index).
//	opa:	The operation type:
//		0:	a' = a,
//		1:	a' = a^T.
//	x:	The input vector, at least of size n+ix.
//	ix:	The sub-vector offset in x.
//	beta:	A coefficient.
//	y:	A preallocated output array, at least of size m+iy.
//	iy:	The sub-vector offset in y.
// Output:
//	y:	The vector which stores the result.
//
// Handling of Special Cases:
// --------------------------
// *	If m == 0, then the subroutine does nothing.
//	It does not even touch the arrays.
// *	If n == 0 or alpha == 0, then:
//	-	If beta == 0, then y is filled by 0's.
//		a and x are not referenced at all.
//		The initial values of y are ignored (we do not multiply y by 0, we just rewrite it by 0's).
//	-	If beta != 0, then y is multiplied by beta.
// *	If m > 0, n > 0, alpha != 0, but beta == 0, then y is replaced by a x;
//	the initial state of y is ignored (rewritten by a x, without initial multiplication by 0's).
// ALGLIB Routine: Copyright 01.09.2021 by Sergey Bochkanov
void rgemvx(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
// Properly premultiply y by beta.
// Quick exit for m == 0, n == 0 or alpha == 0.
// After this block we have m > 0, n > 0, alpha != 0.
   if (m <= 0) return;
   if (beta != 0.0) rmulvx(m, beta, y, iy); else rsetvx(m, 0.0, y, iy);
   if (n <= 0 || alpha == 0.0) return;
// Straight or transposed?
   switch (opa) {
   // y += a x.
      case 0:
#if !defined ALGLIB_NO_FAST_KERNELS
      // Try SIMD code.
         if (n >= _ABLASF_KERNEL_SIZE2) KerSubAvx2Fma(rgemvx_straight(m, n, alpha, a, ia, ja, x->xR + ix, y->xR + iy))
#endif
      // Generic C code.
         for (ae_int_t i = 0; i < m; i++) {
            double v = 0.0;
            for (ae_int_t j = 0; j < n; j++) v += a->xyR[ia + i][ja + j] * x->xR[ix + j];
            y->xR[iy + i] += alpha * v;
         }
      break;
   // y += a^T x.
      case 1:
#if !defined ALGLIB_NO_FAST_KERNELS
      // Try SIMD code.
         if (m >= _ABLASF_KERNEL_SIZE2) KerSubAvx2Fma(rgemvx_transposed(m, n, alpha, a, ia, ja, x->xR + ix, y->xR + iy))
#endif
      // Generic C code.
         for (ae_int_t i = 0; i < n; i++) {
            double v = alpha * x->xR[ix + i];
            for (ae_int_t j = 0; j < m; j++) y->xR[iy + j] += v * a->xyR[ia + i][ja + j];
         }
      break;
   }
}

// Rank-1 correction: a += alpha u v^T.
// Note:
// *	This function expects a to be large enough to store the result.
//	No automatic preallocation happens for smaller arrays.
//	No integrity checks are performed on the sizes of a, u, v.
// Inputs:
//	m:	The number of rows.
//	n:	The number of columns.
//	alpha:	A coefficient.
//	u:	Vector #1.
//	v:	Vector #2.
//	a:	The target m x n matrix.
// ALGLIB Routine: Copyright 07.09.2021 by Sergey Bochkanov
void rger(ae_int_t m, ae_int_t n, double alpha, RVector *u, RVector *v, RMatrix *a) {
   if (m <= 0 || n <= 0 || alpha == 0.0) return;
   for (ae_int_t i = 0; i < m; i++) {
      double s = alpha * u->xR[i];
      for (ae_int_t j = 0; j < n; j++) a->xyR[i][j] += s * v->xR[j];
   }
}

// Solve the linear system a' x == b for x, where:
// *	a is an n x n upper/lower triangular/unitriangular matrix,
// *	x and b are n x 1 vectors,
// *	the operation a' may be the identity operation a |-> a, or transposition a |-> a^T.
// The solution replaces x.
// IMPORTANT:
// *	No overflow/underflow/denegeracy tests are performed.
// *	No integrity checks for operand sizes, out-of-bounds accesses and so on are performed.
// Inputs:
//	n:	The matrix size; n >= 0.
//	a:	A matrix containing the n x n matrix in elements [ia,ia+n) x [ja,ja+n).
//	ia:	The row sub-matrix offset.
//	ja:	The column sub-matrix offset.
//	isupper: Whether the matrix is upper triangular.
//	isunit:	Whether the matrix is unitriangular.
//	opa:	The transformation type:
//		0:	no transformation: a' = a,
//		1:	transposition: a' = a^T.
//	x:	A vector containing the right-hand side n-vector in elements [ix,ix+n).
//	ix:	The offset in x.
// Output:
//	x:	The solution n-vector in elements [ix,ix+n) of x.
// ALGLIB Routine: Copyright 07.09.2021 by Sergey Bochkanov
void rtrsvx(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t opa, RVector *x, ae_int_t ix) {
   if (n <= 0) return;
   else switch (opa) {
      case 0:
         if (isupper) for (ae_int_t i = n - 1; i >= 0; i--) {
            double v = x->xR[ix + i];
            for (ae_int_t j = i + 1; j < n; j++) v -= a->xyR[ia + i][ja + j] * x->xR[ix + j];
            if (!isunit) v /= a->xyR[ia + i][ja + i];
            x->xR[ix + i] = v;
         } else for (ae_int_t i = 0; i < n; i++) {
            double v = x->xR[ix + i];
            for (ae_int_t j = 0; j < i; j++) v -= a->xyR[ia + i][ja + j] * x->xR[ix + j];
            if (!isunit) v /= a->xyR[ia + i][ja + i];
            x->xR[ix + i] = v;
         }
      break;
      case 1:
         if (isupper) for (ae_int_t i = 0; i < n; i++) {
            double v = x->xR[ix + i];
            if (!isunit) v /= a->xyR[ia + i][ja + i];
            x->xR[ix + i] = v;
            if (v == 0.0) continue;
            for (ae_int_t j = i + 1; j < n; j++) x->xR[ix + j] -= v * a->xyR[ia + i][ja + j];
         } else for (ae_int_t i = n - 1; i >= 0; i--) {
            double v = x->xR[ix + i];
            if (!isunit) v /= a->xyR[ia + i][ja + i];
            x->xR[ix + i] = v;
            if (v == 0.0) continue;
            for (ae_int_t j = 0; j < i; j++) x->xR[ix + j] -= v * a->xyR[ia + i][ja + j];
         }
      break;
      default: ae_assert(false, "rtrsvx: unexpected operation type"); break;
   }
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixgerf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double ralpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_rmatrixgerf(m, n, a, ia, ja, ralpha, u, iu, v, iv);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixrank1f(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixrank1f(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixlefttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, RMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixlefttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, CMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixrighttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, RMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixrighttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, CMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_rmatrixsyrkf(n, k, alpha, a, ia, ja, opa, beta, c, ic, jc, isupper);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_cmatrixherkf(n, k, alpha, a, ia, ja, opa, beta, c, ic, jc, isupper);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_rmatrixgemmf(m, n, k, alpha, a, ia, ja, opa, b, ib, jb, opb, beta, c, ic, jc);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   return false;
#else
   return _ialglib_i_cmatrixgemmf(m, n, k, alpha, a, ia, ja, opa, b, ib, jb, opb, beta, c, ic, jc);
#endif
}

// rmatrixgemm() kernels: the base case code for rmatrixgemm(), and its specializations to (opa, opb) in {0,1} x {0,1}.
//	c = alpha a' b' + beta c
// where:
// *	c is an m x n general matrix, a' is an m x k matrix, b' is a k x n matrix, with:
// *	a' = a if opa == 0, a' = a^T if opa == 1,
// *	b' = b if opb == 0, b' = b^T if opb == 1.
// Additional information:
// *	The multiplication result replaces c.
// *	If beta == 0, c is not used in the calculations (not multiplied by 0 - just not referenced).
// *	If alpha == 0, a is not used (not multiplied by 0 - just not referenced).
// *	If both beta == 0 and alpha == 0, c is filled by 0's.
// *	The specialized functions require that alpha != 0; an exception will be thrown otherwise.
// Important:
// *	This function requires the output matrix c to be preallocated.
//	An exception will be thrown, if it is not large enough to store the result.
// Inputs:
//	m, n, k:	The matrix sizes; m, n, k > 0.
//	alpha, beta:	The coefficients.
//	a, ia, ja:	A matrix and its sub-matrix offsets.
//	opa:		The transformation type for a:
//			0:	no transformation: a' = a,
//			1:	transposition: a' = a^T.
//	b, ib, jb:	A matrix and its sub-matrix offsets.
//	opb:		The transformation type for b:
//			0:	no transformation: b' = b,
//			1:	transposition: b' = b^T.
//	c, ic, jc:	The m x n PREALLOCATED output matrix and its sub-matrix offsets.
// ALGLIB Routines: Copyright 27.03.2013 by Sergey Bochkanov

// rmatrixgemm() kernel: the base case code for rmatrixgemm(), specialized to opa == 0 and opb == 0.
void rmatrixgemmk44v00(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_assert(alpha != 0.0, "rmatrixgemmk44v00: internal error (alpha == 0)");
// Size-0 matrices require nothing to be done.
   if (m == 0 || n == 0) return;
// a b:
   for (ae_int_t i = 0; i < m; i += 4) for (ae_int_t j = 0; j < n; j += 4) {
   // Choose between the specialized 4 x 4 code and the general code.
      if (i + 4 <= m && j + 4 <= n) {
      // Specialized 4 x 4 code for the sub-matrix of c at elements [i,i+4) x [j,j+4):
      // a sum of k rank-1 products, with operands cached in local variables in order to speed up array operations.
         ae_int_t ax0 = ia + i, ax1 = ax0 + 1, ax2 = ax0 + 2, ax3 = ax0 + 3, ay = ja;
         ae_int_t bx0 = jb + j, bx1 = bx0 + 1, bx2 = bx0 + 2, bx3 = bx0 + 3, by = ib;
         double v00 = 0.0, v01 = 0.0, v02 = 0.0, v03 = 0.0;
         double v10 = 0.0, v11 = 0.0, v12 = 0.0, v13 = 0.0;
         double v20 = 0.0, v21 = 0.0, v22 = 0.0, v23 = 0.0;
         double v30 = 0.0, v31 = 0.0, v32 = 0.0, v33 = 0.0;
      // Different variants of the internal loop.
         for (ae_int_t t = 0; t < k; ay++, by++, t++) {
            double a0 = a->xyR[ax0][ay], a1 = a->xyR[ax1][ay], a2 = a->xyR[ax2][ay], a3 = a->xyR[ax3][ay];
            double b0 = b->xyR[by][bx0], b1 = b->xyR[by][bx1], b2 = b->xyR[by][bx2], b3 = b->xyR[by][bx3];
            v00 += a0 * b0, v01 += a0 * b1, v02 += a0 * b2, v03 += a0 * b3;
            v10 += a1 * b0, v11 += a1 * b1, v12 += a1 * b2, v13 += a1 * b3;
            v20 += a2 * b0, v21 += a2 * b1, v22 += a2 * b2, v23 += a2 * b3;
            v30 += a3 * b0, v31 += a3 * b1, v32 += a3 * b2, v33 += a3 * b3;
         }
         v00 *= alpha, v01 *= alpha, v02 *= alpha, v03 *= alpha;
         v10 *= alpha, v11 *= alpha, v12 *= alpha, v13 *= alpha;
         v20 *= alpha, v21 *= alpha, v22 *= alpha, v23 *= alpha;
         v30 *= alpha, v31 *= alpha, v32 *= alpha, v33 *= alpha;
         ae_int_t ic0 = ic + i, ic1 = ic0 + 1, ic2 = ic0 + 2, ic3 = ic0 + 3;
         ae_int_t jc0 = jc + j, jc1 = jc0 + 1, jc2 = jc0 + 2, jc3 = jc0 + 3;
         if (beta == 0.0) {
            c->xyR[ic0][jc0] = v00, c->xyR[ic0][jc1] = v01, c->xyR[ic0][jc2] = v02, c->xyR[ic0][jc3] = v03;
            c->xyR[ic1][jc0] = v10, c->xyR[ic1][jc1] = v11, c->xyR[ic1][jc2] = v12, c->xyR[ic1][jc3] = v13;
            c->xyR[ic2][jc0] = v20, c->xyR[ic2][jc1] = v21, c->xyR[ic2][jc2] = v22, c->xyR[ic2][jc3] = v23;
            c->xyR[ic3][jc0] = v30, c->xyR[ic3][jc1] = v31, c->xyR[ic3][jc2] = v32, c->xyR[ic3][jc3] = v33;
         } else {
            c->xyR[ic0][jc0] = beta * c->xyR[ic0][jc0] + v00, c->xyR[ic0][jc1] = beta * c->xyR[ic0][jc1] + v01;
            c->xyR[ic0][jc2] = beta * c->xyR[ic0][jc2] + v02, c->xyR[ic0][jc3] = beta * c->xyR[ic0][jc3] + v03;
            c->xyR[ic1][jc0] = beta * c->xyR[ic1][jc0] + v10, c->xyR[ic1][jc1] = beta * c->xyR[ic1][jc1] + v11;
            c->xyR[ic1][jc2] = beta * c->xyR[ic1][jc2] + v12, c->xyR[ic1][jc3] = beta * c->xyR[ic1][jc3] + v13;
            c->xyR[ic2][jc0] = beta * c->xyR[ic2][jc0] + v20, c->xyR[ic2][jc1] = beta * c->xyR[ic2][jc1] + v21;
            c->xyR[ic2][jc2] = beta * c->xyR[ic2][jc2] + v22, c->xyR[ic2][jc3] = beta * c->xyR[ic2][jc3] + v23;
            c->xyR[ic3][jc0] = beta * c->xyR[ic3][jc0] + v30, c->xyR[ic3][jc1] = beta * c->xyR[ic3][jc1] + v31;
            c->xyR[ic3][jc2] = beta * c->xyR[ic3][jc2] + v32, c->xyR[ic3][jc3] = beta * c->xyR[ic3][jc3] + v33;
         }
      } else {
      // Determine the sub-matrix of c of which elements [i0,i1) x [j0,j1) to process.
         ae_int_t i0 = i, i1 = imin2(i + 4, m);
         ae_int_t j0 = j, j1 = imin2(j + 4, n);
      // Process the sub-matrix.
         for (ae_int_t ik = i0; ik < i1; ik++) for (ae_int_t jk = j0; jk < j1; jk++) {
            double v = k == 0 || alpha == 0.0 ? 0.0 : alpha * ae_v_dotproduct(&a->xyR[ia + ik][ja], 1, &b->xyR[ib][jb + jk], b->stride, k);
            c->xyR[ic + ik][jc + jk] = beta == 0.0 ? v : beta * c->xyR[ic + ik][jc + jk] + v;
         }
      }
   }
}

// rmatrixgemm() kernel: the base case code for rmatrixgemm(), specialized to opa == 0 and opb == 1.
void rmatrixgemmk44v01(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_assert(alpha != 0.0, "rmatrixgemmk44v01: internal error (alpha == 0)");
// Size-0 matrices require nothing to be done.
   if (m == 0 || n == 0) return;
// a b':
   for (ae_int_t i = 0; i < m; i += 4) for (ae_int_t j = 0; j < n; j += 4) {
   // Choose between the specialized 4 x 4 code and the general code.
      if (i + 4 <= m && j + 4 <= n) {
      // Specialized 4 x 4 code for the sub-matrix of c at elements [i,i+4) x [j,j+4):
      // a sum of k rank-1 products, with operands cached in local variables in order to speed up array operations.
         ae_int_t ax0 = ia + i, ax1 = ax0 + 1, ax2 = ax0 + 2, ax3 = ax0 + 3, ay = ja;
         ae_int_t bx0 = ib + j, bx1 = bx0 + 1, bx2 = bx0 + 2, bx3 = bx0 + 3, by = jb;
         double v00 = 0.0, v01 = 0.0, v02 = 0.0, v03 = 0.0;
         double v10 = 0.0, v11 = 0.0, v12 = 0.0, v13 = 0.0;
         double v20 = 0.0, v21 = 0.0, v22 = 0.0, v23 = 0.0;
         double v30 = 0.0, v31 = 0.0, v32 = 0.0, v33 = 0.0;
      // Different variants of the internal loop.
         for (ae_int_t t = 0; t < k; ay++, by++, t++) {
            double a0 = a->xyR[ax0][ay], a1 = a->xyR[ax1][ay], a2 = a->xyR[ax2][ay], a3 = a->xyR[ax3][ay];
            double b0 = b->xyR[bx0][by], b1 = b->xyR[bx1][by], b2 = b->xyR[bx2][by], b3 = b->xyR[bx3][by];
            v00 += a0 * b0, v01 += a0 * b1, v02 += a0 * b2, v03 += a0 * b3;
            v10 += a1 * b0, v11 += a1 * b1, v12 += a1 * b2, v13 += a1 * b3;
            v20 += a2 * b0, v21 += a2 * b1, v22 += a2 * b2, v23 += a2 * b3;
            v30 += a3 * b0, v31 += a3 * b1, v32 += a3 * b2, v33 += a3 * b3;
         }
         v00 *= alpha, v01 *= alpha, v02 *= alpha, v03 *= alpha;
         v10 *= alpha, v11 *= alpha, v12 *= alpha, v13 *= alpha;
         v20 *= alpha, v21 *= alpha, v22 *= alpha, v23 *= alpha;
         v30 *= alpha, v31 *= alpha, v32 *= alpha, v33 *= alpha;
         ae_int_t ic0 = ic + i, ic1 = ic0 + 1, ic2 = ic0 + 2, ic3 = ic0 + 3;
         ae_int_t jc0 = jc + j, jc1 = jc0 + 1, jc2 = jc0 + 2, jc3 = jc0 + 3;
         if (beta == 0.0) {
            c->xyR[ic0][jc0] = v00, c->xyR[ic0][jc1] = v01, c->xyR[ic0][jc2] = v02, c->xyR[ic0][jc3] = v03;
            c->xyR[ic1][jc0] = v10, c->xyR[ic1][jc1] = v11, c->xyR[ic1][jc2] = v12, c->xyR[ic1][jc3] = v13;
            c->xyR[ic2][jc0] = v20, c->xyR[ic2][jc1] = v21, c->xyR[ic2][jc2] = v22, c->xyR[ic2][jc3] = v23;
            c->xyR[ic3][jc0] = v30, c->xyR[ic3][jc1] = v31, c->xyR[ic3][jc2] = v32, c->xyR[ic3][jc3] = v33;
         } else {
            c->xyR[ic0][jc0] = beta * c->xyR[ic0][jc0] + v00, c->xyR[ic0][jc1] = beta * c->xyR[ic0][jc1] + v01;
            c->xyR[ic0][jc2] = beta * c->xyR[ic0][jc2] + v02, c->xyR[ic0][jc3] = beta * c->xyR[ic0][jc3] + v03;
            c->xyR[ic1][jc0] = beta * c->xyR[ic1][jc0] + v10, c->xyR[ic1][jc1] = beta * c->xyR[ic1][jc1] + v11;
            c->xyR[ic1][jc2] = beta * c->xyR[ic1][jc2] + v12, c->xyR[ic1][jc3] = beta * c->xyR[ic1][jc3] + v13;
            c->xyR[ic2][jc0] = beta * c->xyR[ic2][jc0] + v20, c->xyR[ic2][jc1] = beta * c->xyR[ic2][jc1] + v21;
            c->xyR[ic2][jc2] = beta * c->xyR[ic2][jc2] + v22, c->xyR[ic2][jc3] = beta * c->xyR[ic2][jc3] + v23;
            c->xyR[ic3][jc0] = beta * c->xyR[ic3][jc0] + v30, c->xyR[ic3][jc1] = beta * c->xyR[ic3][jc1] + v31;
            c->xyR[ic3][jc2] = beta * c->xyR[ic3][jc2] + v32, c->xyR[ic3][jc3] = beta * c->xyR[ic3][jc3] + v33;
         }
      } else {
      // Determine the sub-matrix of c of which elements [i0,i1) x [j0,j1) to process.
         ae_int_t i0 = i, i1 = imin2(i + 4, m);
         ae_int_t j0 = j, j1 = imin2(j + 4, n);
      // Process the sub-matrix.
         for (ae_int_t ik = i0; ik < i1; ik++) for (ae_int_t jk = j0; jk < j1; jk++) {
            double v = k == 0 || alpha == 0.0 ? 0.0 : alpha * ae_v_dotproduct(&a->xyR[ia + ik][ja], 1, &b->xyR[ib + jk][jb], 1, k);
            c->xyR[ic + ik][jc + jk] = beta == 0.0 ? v : beta * c->xyR[ic + ik][jc + jk] + v;
         }
      }
   }
}

// rmatrixgemm() kernel: the base case code for rmatrixgemm(), specialized to opa == 1 and opb == 0.
void rmatrixgemmk44v10(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_assert(alpha != 0.0, "rmatrixgemmk44v10: internal error (alpha == 0)");
// Size-0 matrices require nothing to be done.
   if (m == 0 || n == 0) return;
// a' b:
   for (ae_int_t i = 0; i < m; i += 4) for (ae_int_t j = 0; j < n; j += 4) {
   // Choose between the specialized 4 x 4 code and the general code.
      if (i + 4 <= m && j + 4 <= n) {
      // Specialized 4 x 4 code for the sub-matrix of c at elements [i,i+4) x [j,j+4):
      // a sum of k rank-1 products, with operands cached in local variables in order to speed up array operations.
         ae_int_t ax0 = ja + i, ax1 = ax0 + 1, ax2 = ax0 + 2, ax3 = ax0 + 3, ay = ia;
         ae_int_t bx0 = jb + j, bx1 = bx0 + 1, bx2 = bx0 + 2, bx3 = bx0 + 3, by = ib;
         double v00 = 0.0, v01 = 0.0, v02 = 0.0, v03 = 0.0;
         double v10 = 0.0, v11 = 0.0, v12 = 0.0, v13 = 0.0;
         double v20 = 0.0, v21 = 0.0, v22 = 0.0, v23 = 0.0;
         double v30 = 0.0, v31 = 0.0, v32 = 0.0, v33 = 0.0;
      // Different variants of the internal loop.
         for (ae_int_t t = 0; t < k; ay++, by++, t++) {
            double a0 = a->xyR[ay][ax0], a1 = a->xyR[ay][ax1], a2 = a->xyR[ay][ax2], a3 = a->xyR[ay][ax3];
            double b0 = b->xyR[by][bx0], b1 = b->xyR[by][bx1], b2 = b->xyR[by][bx2], b3 = b->xyR[by][bx3];
            v00 += a0 * b0, v01 += a0 * b1, v02 += a0 * b2, v03 += a0 * b3;
            v10 += a1 * b0, v11 += a1 * b1, v12 += a1 * b2, v13 += a1 * b3;
            v20 += a2 * b0, v21 += a2 * b1, v22 += a2 * b2, v23 += a2 * b3;
            v30 += a3 * b0, v31 += a3 * b1, v32 += a3 * b2, v33 += a3 * b3;
         }
         v00 *= alpha, v01 *= alpha, v02 *= alpha, v03 *= alpha;
         v10 *= alpha, v11 *= alpha, v12 *= alpha, v13 *= alpha;
         v20 *= alpha, v21 *= alpha, v22 *= alpha, v23 *= alpha;
         v30 *= alpha, v31 *= alpha, v32 *= alpha, v33 *= alpha;
         ae_int_t ic0 = ic + i, ic1 = ic0 + 1, ic2 = ic0 + 2, ic3 = ic0 + 3;
         ae_int_t jc0 = jc + j, jc1 = jc0 + 1, jc2 = jc0 + 2, jc3 = jc0 + 3;
         if (beta == 0.0) {
            c->xyR[ic0][jc0] = v00, c->xyR[ic0][jc1] = v01, c->xyR[ic0][jc2] = v02, c->xyR[ic0][jc3] = v03;
            c->xyR[ic1][jc0] = v10, c->xyR[ic1][jc1] = v11, c->xyR[ic1][jc2] = v12, c->xyR[ic1][jc3] = v13;
            c->xyR[ic2][jc0] = v20, c->xyR[ic2][jc1] = v21, c->xyR[ic2][jc2] = v22, c->xyR[ic2][jc3] = v23;
            c->xyR[ic3][jc0] = v30, c->xyR[ic3][jc1] = v31, c->xyR[ic3][jc2] = v32, c->xyR[ic3][jc3] = v33;
         } else {
            c->xyR[ic0][jc0] = beta * c->xyR[ic0][jc0] + v00, c->xyR[ic0][jc1] = beta * c->xyR[ic0][jc1] + v01;
            c->xyR[ic0][jc2] = beta * c->xyR[ic0][jc2] + v02, c->xyR[ic0][jc3] = beta * c->xyR[ic0][jc3] + v03;
            c->xyR[ic1][jc0] = beta * c->xyR[ic1][jc0] + v10, c->xyR[ic1][jc1] = beta * c->xyR[ic1][jc1] + v11;
            c->xyR[ic1][jc2] = beta * c->xyR[ic1][jc2] + v12, c->xyR[ic1][jc3] = beta * c->xyR[ic1][jc3] + v13;
            c->xyR[ic2][jc0] = beta * c->xyR[ic2][jc0] + v20, c->xyR[ic2][jc1] = beta * c->xyR[ic2][jc1] + v21;
            c->xyR[ic2][jc2] = beta * c->xyR[ic2][jc2] + v22, c->xyR[ic2][jc3] = beta * c->xyR[ic2][jc3] + v23;
            c->xyR[ic3][jc0] = beta * c->xyR[ic3][jc0] + v30, c->xyR[ic3][jc1] = beta * c->xyR[ic3][jc1] + v31;
            c->xyR[ic3][jc2] = beta * c->xyR[ic3][jc2] + v32, c->xyR[ic3][jc3] = beta * c->xyR[ic3][jc3] + v33;
         }
      } else {
      // Determine the sub-matrix of c of which elements [i0,i1) x [j0,j1) to process.
         ae_int_t i0 = i, i1 = imin2(i + 4, m);
         ae_int_t j0 = j, j1 = imin2(j + 4, n);
      // Process the sub-matrix.
         for (ae_int_t ik = i0; ik < i1; ik++) for (ae_int_t jk = j0; jk < j1; jk++) {
            double v = k == 0 || alpha == 0.0 ? 0.0 : alpha * ae_v_dotproduct(&a->xyR[ia][ja + ik], a->stride, &b->xyR[ib][jb + jk], b->stride, k);
            c->xyR[ic + ik][jc + jk] = beta == 0.0 ? v : beta * c->xyR[ic + ik][jc + jk] + v;
         }
      }
   }
}

// rmatrixgemm() kernel: the base case code for rmatrixgemm(), specialized to opa == 1 and opb == 1.
void rmatrixgemmk44v11(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_assert(alpha != 0.0, "rmatrixgemmk44v11: internal error (alpha == 0)");
// Size-0 matrices require nothing to be done.
   if (m == 0 || n == 0) return;
// a' b':
   for (ae_int_t i = 0; i < m; i += 4) for (ae_int_t j = 0; j < n; j += 4) {
   // Choose between the specialized 4 x 4 code and the general code.
      if (i + 4 <= m && j + 4 <= n) {
      // Specialized 4 x 4 code for the sub-matrix of c at elements [i,i+4) x [j,j+4):
      // a sum of k rank-1 products, with operands cached in local variables in order to speed up array operations.
         ae_int_t ax0 = ja + i, ax1 = ax0 + 1, ax2 = ax0 + 2, ax3 = ax0 + 3, ay = ia;
         ae_int_t bx0 = ib + j, bx1 = bx0 + 1, bx2 = bx0 + 2, bx3 = bx0 + 3, by = jb;
         double v00 = 0.0, v01 = 0.0, v02 = 0.0, v03 = 0.0;
         double v10 = 0.0, v11 = 0.0, v12 = 0.0, v13 = 0.0;
         double v20 = 0.0, v21 = 0.0, v22 = 0.0, v23 = 0.0;
         double v30 = 0.0, v31 = 0.0, v32 = 0.0, v33 = 0.0;
      // Different variants of the internal loop.
         for (ae_int_t t = 0; t < k; ay++, by++, t++) {
            double a0 = a->xyR[ay][ax0], a1 = a->xyR[ay][ax1], a2 = a->xyR[ay][ax2], a3 = a->xyR[ay][ax3];
            double b0 = b->xyR[bx0][by], b1 = b->xyR[bx1][by], b2 = b->xyR[bx2][by], b3 = b->xyR[bx3][by];
            v00 += a0 * b0, v01 += a0 * b1, v02 += a0 * b2, v03 += a0 * b3;
            v10 += a1 * b0, v11 += a1 * b1, v12 += a1 * b2, v13 += a1 * b3;
            v20 += a2 * b0, v21 += a2 * b1, v22 += a2 * b2, v23 += a2 * b3;
            v30 += a3 * b0, v31 += a3 * b1, v32 += a3 * b2, v33 += a3 * b3;
         }
         v00 *= alpha, v01 *= alpha, v02 *= alpha, v03 *= alpha;
         v10 *= alpha, v11 *= alpha, v12 *= alpha, v13 *= alpha;
         v20 *= alpha, v21 *= alpha, v22 *= alpha, v23 *= alpha;
         v30 *= alpha, v31 *= alpha, v32 *= alpha, v33 *= alpha;
         ae_int_t ic0 = ic + i, ic1 = ic0 + 1, ic2 = ic0 + 2, ic3 = ic0 + 3;
         ae_int_t jc0 = jc + j, jc1 = jc0 + 1, jc2 = jc0 + 2, jc3 = jc0 + 3;
         if (beta == 0.0) {
            c->xyR[ic0][jc0] = v00, c->xyR[ic0][jc1] = v01, c->xyR[ic0][jc2] = v02, c->xyR[ic0][jc3] = v03;
            c->xyR[ic1][jc0] = v10, c->xyR[ic1][jc1] = v11, c->xyR[ic1][jc2] = v12, c->xyR[ic1][jc3] = v13;
            c->xyR[ic2][jc0] = v20, c->xyR[ic2][jc1] = v21, c->xyR[ic2][jc2] = v22, c->xyR[ic2][jc3] = v23;
            c->xyR[ic3][jc0] = v30, c->xyR[ic3][jc1] = v31, c->xyR[ic3][jc2] = v32, c->xyR[ic3][jc3] = v33;
         } else {
            c->xyR[ic0][jc0] = beta * c->xyR[ic0][jc0] + v00, c->xyR[ic0][jc1] = beta * c->xyR[ic0][jc1] + v01;
            c->xyR[ic0][jc2] = beta * c->xyR[ic0][jc2] + v02, c->xyR[ic0][jc3] = beta * c->xyR[ic0][jc3] + v03;
            c->xyR[ic1][jc0] = beta * c->xyR[ic1][jc0] + v10, c->xyR[ic1][jc1] = beta * c->xyR[ic1][jc1] + v11;
            c->xyR[ic1][jc2] = beta * c->xyR[ic1][jc2] + v12, c->xyR[ic1][jc3] = beta * c->xyR[ic1][jc3] + v13;
            c->xyR[ic2][jc0] = beta * c->xyR[ic2][jc0] + v20, c->xyR[ic2][jc1] = beta * c->xyR[ic2][jc1] + v21;
            c->xyR[ic2][jc2] = beta * c->xyR[ic2][jc2] + v22, c->xyR[ic2][jc3] = beta * c->xyR[ic2][jc3] + v23;
            c->xyR[ic3][jc0] = beta * c->xyR[ic3][jc0] + v30, c->xyR[ic3][jc1] = beta * c->xyR[ic3][jc1] + v31;
            c->xyR[ic3][jc2] = beta * c->xyR[ic3][jc2] + v32, c->xyR[ic3][jc3] = beta * c->xyR[ic3][jc3] + v33;
         }
      } else {
      // Determine the sub-matrix of c of which elements [i0,i1) x [j0,j1) to process.
         ae_int_t i0 = i, i1 = imin2(i + 4, m);
         ae_int_t j0 = j, j1 = imin2(j + 4, n);
      // Process the sub-matrix.
         for (ae_int_t ik = i0; ik < i1; ik++) for (ae_int_t jk = j0; jk < j1; jk++) {
            double v = k == 0 || alpha == 0.0 ? 0.0 : alpha * ae_v_dotproduct(&a->xyR[ia][ja + ik], a->stride, &b->xyR[ib + jk][jb], 1, k);
            c->xyR[ic + ik][jc + jk] = beta == 0.0 ? v : beta * c->xyR[ic + ik][jc + jk] + v;
         }
      }
   }
}

#if 0
//(@) Returned to Ap.cpp, because merely putting it here, instead of in Ap.cpp, seems to degrade the speed of the GEMM routine!
//(@) The slow-down is apparently the same issue, inherited from ALGLIB, that brought about the splitting of rmatrixgemmk() in ALGLIB.
// Fast rmatrixgemm() kernel: with AVX2/FMA support.
static bool ablasf_rgemm32basecase(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
#endif

// rmatrixgemm() kernel: the base case code for rmatrixgemm().
void rmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
// Size-0 matrices require nothing to be done.
   if (m == 0 || n == 0) return;
// Try the optimized code.
   if (ablasf_rgemm32basecase(m, n, k, alpha, a, ia, ja, opa, b, ib, jb, opb, beta, c, ic, jc)) return;
// If k == 0 or alpha == 0, then c *= beta.
   if (k == 0 || alpha == 0.0) {
      if (beta == 0.0)
         for (ae_int_t i = 0; i < m; i++) for (ae_int_t j = 0; j < n; j++) c->xyR[ic + i][jc + j] = 0.0;
      else if (beta != 1.0)
         for (ae_int_t i = 0; i < m; i++) for (ae_int_t j = 0; j < n; j++) c->xyR[ic + i][jc + j] *= beta;
      return;
   }
// Call specialized code.
// Note:
// *	The specialized code was moved to separate functions because of strange issues with the instruction cache on some systems.
//	Having too-long functions significantly slows down the internal loop of the algorithm.
   switch (opa) {
      case 0: switch (opb) {
         case 0: rmatrixgemmk44v00(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc); break;
         default:
         case 1: rmatrixgemmk44v01(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc); break;
      }
      break;
      default:
      case 1: switch (opb) {
         case 0: rmatrixgemmk44v10(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc); break;
         default:
         case 1: rmatrixgemmk44v11(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc); break;
      }
      break;
   }
}

// cmatrixgemm() kernels: the base case code for cmatrixgemm().
//	c = alpha a' b' + beta c
// where:
// *	c is an m x n general matrix, a' is an m x k matrix, b' is a k x n matrix, with:
// *	a' = a if opa == 0, a' = a^T if opa == 1, a' = a^+ if opa == 2,
// *	b' = b if opb == 0, b' = b^T if opb == 1, b' = b^+ if opb == 2.
// Additional information:
// *	The multiplication result replaces c.
// *	If beta == 0, c is not used in the calculations (not multiplied by 0 - just not referenced).
// *	If alpha == 0, a is not used (not multiplied by 0 - just not referenced).
// *	If both beta == 0 and alpha == 0, c is filled by 0's.
// Important:
// *	This function requires the output matrix c to be preallocated.
//	An exception will be thrown, if it is not large enough to store the result.
// Inputs:
//	m, n, k:	The matrix sizes; m, n, k > 0.
//	alpha, beta:	The coefficients.
//	a, ia, ja:	A matrix and its sub-matrix offsets.
//	opa:		The transformation type for a:
//			0:	no transformation: a' = a,
//			1:	transposition: a' = a^T,
//			2:	conjugate transposition: a' = a^+.
//	b, ib, jb:	A matrix and its sub-matrix offsets.
//	opb:		The transformation type for b:
//			0:	no transformation: b' = b,
//			1:	transposition: b' = b^T,
//			2:	conjugate transposition: b' = b^+.
//	c, ic, jc:	The m x n PREALLOCATED output matrix and its sub-matrix offsets.
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void cmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc) {
// Size-0 matrices require nothing to be done.
   if (m == 0 || n == 0) return;
// Try the optimized code.
   if (cmatrixgemmf(m, n, k, alpha, a, ia, ja, opa, b, ib, jb, opb, beta, c, ic, jc)) return;
// If k == 0 or alpha == 0, then c *= beta.
   if (k == 0 || ae_c_eq_d(alpha, 0.0)) {
      if (ae_c_eq_d(beta, 0.0))
         for (ae_int_t i = 0; i < m; i++) for (ae_int_t j = 0; j < n; j++) c->xyC[ic + i][jc + j] = complex_from_i(0);
      else if (ae_c_neq_d(beta, 1.0))
         for (ae_int_t i = 0; i < m; i++) for (ae_int_t j = 0; j < n; j++) c->xyC[ic + i][jc + j] = ae_c_mul(beta, c->xyC[ic + i][jc + j]);
      return;
   }
// This phase is not really necessary, but the compiler complains about "possibly uninitialized variables".
   double a0x = 0.0, a0y = 0.0, a1x = 0.0, a1y = 0.0;
   double b0x = 0.0, b0y = 0.0, b1x = 0.0, b1y = 0.0;
// The general case.
   for (ae_int_t i = 0; i < m; i += 2) for (ae_int_t j = 0; j < n; j += 2) {
   // Choose between the specialized 4 x 4 code and the general code.
      if (i + 2 <= m && j + 2 <= n) {
      // Specialized 4 x 4 code for the sub-matrix of c at elements [i,i+4) x [j,j+4):
      // a sum of k rank-1 products, with operands cached in local variables in order to speed up array operations.
         double v00x = 0.0, v00y = 0.0, v01x = 0.0, v01y = 0.0;
         double v10x = 0.0, v10y = 0.0, v11x = 0.0, v11y = 0.0;
         ae_int_t xa = opa == 0 ? ia : ja, ay = opa == 0 ? ja : ia, ax0 = xa + i, ax1 = ax0 + 1;
         ae_int_t xb = opb == 0 ? jb : ib, by = opb == 0 ? ib : jb, bx0 = xb + j, bx1 = bx0 + 1;
      // Different variants of the internal loop.
         for (ae_int_t t = 0; t < k; ay++, by++, t++) {
            switch (opa) {
               case 0: a0x = a->xyC[ax0][ay].x, a0y = a->xyC[ax0][ay].y, a1x = a->xyC[ax1][ay].x, a1y = a->xyC[ax1][ay].y; break;
               case 1: a0x = a->xyC[ay][ax0].x, a0y = a->xyC[ay][ax0].y, a1x = a->xyC[ay][ax1].x, a1y = a->xyC[ay][ax1].y; break;
               case 2: a0x = a->xyC[ay][ax0].x, a0y = -a->xyC[ay][ax0].y, a1x = a->xyC[ay][ax1].x, a1y = -a->xyC[ay][ax1].y; break;
            }
            switch (opb) {
               case 0: b0x = b->xyC[by][bx0].x, b0y = b->xyC[by][bx0].y, b1x = b->xyC[by][bx1].x, b1y = b->xyC[by][bx1].y; break;
               case 1: b0x = b->xyC[bx0][by].x, b0y = b->xyC[bx0][by].y, b1x = b->xyC[bx1][by].x, b1y = b->xyC[bx1][by].y; break;
               case 2: b0x = b->xyC[bx0][by].x, b0y = -b->xyC[bx0][by].y, b1x = b->xyC[bx1][by].x, b1y = -b->xyC[bx1][by].y; break;
            }
            v00x += a0x * b0x - a0y * b0y, v00y += a0x * b0y + a0y * b0x;
            v01x += a0x * b1x - a0y * b1y, v01y += a0x * b1y + a0y * b1x;
            v10x += a1x * b0x - a1y * b0y, v10y += a1x * b0y + a1y * b0x;
            v11x += a1x * b1x - a1y * b1y, v11y += a1x * b1y + a1y * b1x;
         }
         complex v00 = ae_c_mul(alpha, complex_from_d(v00x, v00y)), v01 = ae_c_mul(alpha, complex_from_d(v01x, v01y));
         complex v10 = ae_c_mul(alpha, complex_from_d(v10x, v10y)), v11 = ae_c_mul(alpha, complex_from_d(v11x, v11y));
         ae_int_t ic0 = ic + i, ic1 = ic0 + 1;
         ae_int_t jc0 = jc + j, jc1 = jc0 + 1;
         if (ae_c_eq_d(beta, 0.0)) {
            c->xyC[ic0][jc0] = v00, c->xyC[ic0][jc1] = v01;
            c->xyC[ic1][jc0] = v10, c->xyC[ic1][jc1] = v11;
         } else {
            c->xyC[ic0][jc0] = ae_c_add(ae_c_mul(beta, c->xyC[ic0][jc0]), v00);
            c->xyC[ic0][jc1] = ae_c_add(ae_c_mul(beta, c->xyC[ic0][jc1]), v01);
            c->xyC[ic1][jc0] = ae_c_add(ae_c_mul(beta, c->xyC[ic1][jc0]), v10);
            c->xyC[ic1][jc1] = ae_c_add(ae_c_mul(beta, c->xyC[ic1][jc1]), v11);
         }
      } else {
      // Determine the sub-matrix of c of which elements [i0,i1) x [j0,j1) to process.
         ae_int_t i0 = i, i1 = imin2(i + 2, m);
         ae_int_t j0 = j, j1 = imin2(j + 2, n);
      // Process the sub-matrix.
         for (ae_int_t ik = i0; ik < i1; ik++) for (ae_int_t jk = j0; jk < j1; jk++) {
            complex v = complex_from_i(0);
            if (k > 0 && ae_c_neq_d(alpha, 0.0)) switch (opa) {
               case 0: switch (opb) {
                  case 0: v = ae_v_cdotproduct(&a->xyC[ia + ik][ja], 1, "N", &b->xyC[ib][jb + jk], b->stride, "N", k); break;
                  case 1: v = ae_v_cdotproduct(&a->xyC[ia + ik][ja], 1, "N", &b->xyC[ib + jk][jb], 1, "N", k); break;
                  case 2: v = ae_v_cdotproduct(&a->xyC[ia + ik][ja], 1, "N", &b->xyC[ib + jk][jb], 1, "Conj", k); break;
               }
               break;
               case 1: switch (opb) {
                  case 0: v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "N", &b->xyC[ib][jb + jk], b->stride, "N", k); break;
                  case 1: v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "N", &b->xyC[ib + jk][jb], 1, "N", k); break;
                  case 2: v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "N", &b->xyC[ib + jk][jb], 1, "Conj", k); break;
               }
               break;
               case 2: switch (opb) {
                  case 0: v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "Conj", &b->xyC[ib][jb + jk], b->stride, "N", k); break;
                  case 1: v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "Conj", &b->xyC[ib + jk][jb], 1, "N", k); break;
                  case 2: v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "Conj", &b->xyC[ib + jk][jb], 1, "Conj", k); break;
               }
               break;
            }
            v = ae_c_mul(alpha, v);
            c->xyC[ic + ik][jc + jk] = ae_c_eq_d(beta, 0.0) ? v : ae_c_add(ae_c_mul(beta, c->xyC[ic + ik][jc + jk]), v);
         }
      }
   }
}
} // end of namespace alglib_impl

// === HBLAS Package ===
namespace alglib_impl {
void hermitianmatrixvectormultiply(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, complex alpha, CVector *y) {
   ae_int_t i;
   ae_int_t ba1;
   ae_int_t by1;
   ae_int_t by2;
   ae_int_t bx1;
   ae_int_t bx2;
   ae_int_t n;
   complex v;
   n = i2 - i1 + 1;
   if (n <= 0) {
      return;
   }
// Let A = L + D + U, where
//	L is strictly lower triangular (main diagonal is zero)
//	D is diagonal
//	U is strictly upper triangular (main diagonal is zero)
//
//	A*x = L*x + D*x + U*x
//
// Calculate D*x first
   for (i = i1; i <= i2; i++) {
      y->xC[i - i1 + 1] = ae_c_mul(a->xyC[i][i], x->xC[i - i1 + 1]);
   }
// Add L*x + U*x
   if (isupper) {
      for (i = i1; i < i2; i++) {
      // Add L*x to the result
         v = x->xC[i - i1 + 1];
         by1 = i - i1 + 2;
         by2 = n;
         ba1 = i + 1;
         ae_v_caddc(&y->xC[by1], 1, &a->xyC[i][ba1], 1, "Conj", by2 - by1 + 1, v);
      // Add U*x to the result
         bx1 = i - i1 + 2;
         bx2 = n;
         ba1 = i + 1;
         v = ae_v_cdotproduct(&x->xC[bx1], 1, "N", &a->xyC[i][ba1], 1, "N", bx2 - bx1 + 1);
         y->xC[i - i1 + 1] = ae_c_add(y->xC[i - i1 + 1], v);
      }
   } else {
      for (i = i1 + 1; i <= i2; i++) {
      // Add L*x to the result
         bx1 = 1;
         bx2 = i - i1;
         ba1 = i1;
         v = ae_v_cdotproduct(&x->xC[bx1], 1, "N", &a->xyC[i][ba1], 1, "N", bx2 - bx1 + 1);
         y->xC[i - i1 + 1] = ae_c_add(y->xC[i - i1 + 1], v);
      // Add U*x to the result
         v = x->xC[i - i1 + 1];
         by1 = 1;
         by2 = i - i1;
         ba1 = i1;
         ae_v_caddc(&y->xC[by1], 1, &a->xyC[i][ba1], 1, "Conj", by2 - by1 + 1, v);
      }
   }
   ae_v_cmulc(&y->xC[1], 1, n, alpha);
}

void hermitianrank2update(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, CVector *y, CVector *t, complex alpha) {
   ae_int_t i;
   ae_int_t tp1;
   ae_int_t tp2;
   complex v;
   if (isupper) {
      for (i = i1; i <= i2; i++) {
         tp1 = i + 1 - i1;
         tp2 = i2 - i1 + 1;
         v = ae_c_mul(alpha, x->xC[i + 1 - i1]);
         ae_v_cmovec(&t->xC[tp1], 1, &y->xC[tp1], 1, "Conj", tp2 - tp1 + 1, v);
         v = ae_c_mul(conj(alpha), y->xC[i + 1 - i1]);
         ae_v_caddc(&t->xC[tp1], 1, &x->xC[tp1], 1, "Conj", tp2 - tp1 + 1, v);
         ae_v_cadd(&a->xyC[i][i], 1, &t->xC[tp1], 1, "N", i2 - i + 1);
      }
   } else {
      for (i = i1; i <= i2; i++) {
         tp1 = 1;
         tp2 = i + 1 - i1;
         v = ae_c_mul(alpha, x->xC[i + 1 - i1]);
         ae_v_cmovec(&t->xC[tp1], 1, &y->xC[tp1], 1, "Conj", tp2 - tp1 + 1, v);
         v = ae_c_mul(conj(alpha), y->xC[i + 1 - i1]);
         ae_v_caddc(&t->xC[tp1], 1, &x->xC[tp1], 1, "Conj", tp2 - tp1 + 1, v);
         ae_v_cadd(&a->xyC[i][i1], 1, &t->xC[tp1], 1, "N", i - i1 + 1);
      }
   }
}
} // end of namespace alglib_impl

// === CREFLECTIONS Package ===
namespace alglib_impl {
// Generation of an elementary complex reflection transformation
//
// The subroutine generates elementary complex reflection H of  order  N,  so
// that, for a given X, the following equality holds true:
//
//      ( X(1) )   ( Beta )
// H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number
//      ( X(n) )   (  0   )
//
// where
//
//               ( V(1) )
// H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )
//               ( V(n) )
//
// where the first component of vector V equals 1.
//
// Inputs:
//     X   -   vector. Array with elements [1..N].
//     N   -   reflection order.
//
// Outputs:
//     X   -   components from 2 to N are replaced by vector V.
//             The first component is replaced with parameter Beta.
//     Tau -   scalar value Tau.
//
// This subroutine is the modification of CLARFG subroutines  from the LAPACK
// library. It has similar functionality except for the fact that it  doesn't
// handle errors when intermediate results cause an overflow.
//
// LAPACK auxiliary routine (version 3.0):
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void complexgeneratereflection(CVector *x, ae_int_t n, complex *tau) {
   ae_int_t j;
   complex alpha;
   double alphi;
   double alphr;
   double beta;
   double xnorm;
   double mx;
   complex t;
   double s;
   complex v;
   tau->x = 0.0;
   tau->y = 0.0;
   if (n <= 0) {
      *tau = complex_from_i(0);
      return;
   }
// Scale if needed (to avoid overflow/underflow during intermediate
// calculations).
   mx = 0.0;
   for (j = 1; j <= n; j++) {
      mx = rmax2(abscomplex(x->xC[j]), mx);
   }
   s = 1.0;
   if (mx != 0.0) {
      if (mx < 1.0) {
         s = sqrt(minrealnumber);
         v = complex_from_d(1.0 / s);
         ae_v_cmulc(&x->xC[1], 1, n, v);
      } else {
         s = sqrt(maxrealnumber);
         v = complex_from_d(1.0 / s);
         ae_v_cmulc(&x->xC[1], 1, n, v);
      }
   }
// calculate
   alpha = x->xC[1];
   mx = 0.0;
   for (j = 2; j <= n; j++) {
      mx = rmax2(abscomplex(x->xC[j]), mx);
   }
   xnorm = 0.0;
   if (mx != 0.0) {
      for (j = 2; j <= n; j++) {
         t = ae_c_div_d(x->xC[j], mx);
         xnorm += ae_c_mul(t, conj(t)).x;
      }
      xnorm = sqrt(xnorm) * mx;
   }
   alphr = alpha.x;
   alphi = alpha.y;
   if (xnorm == 0.0 && alphi == 0.0) {
      *tau = complex_from_i(0);
      x->xC[1] = ae_c_mul_d(x->xC[1], s);
      return;
   }
   mx = rmax2(fabs(alphr), fabs(alphi));
   mx = rmax2(mx, fabs(xnorm));
   beta = -mx * sqrt(sqr(alphr / mx) + sqr(alphi / mx) + sqr(xnorm / mx));
   if (alphr < 0.0) {
      beta = -beta;
   }
   tau->x = (beta - alphr) / beta;
   tau->y = -alphi / beta;
   alpha = ae_c_d_div(1.0, ae_c_sub_d(alpha, beta));
   if (n > 1) {
      ae_v_cmulc(&x->xC[2], 1, n - 1, alpha);
   }
   alpha = complex_from_d(beta);
   x->xC[1] = alpha;
// Scale back
   x->xC[1] = ae_c_mul_d(x->xC[1], s);
}

// Application of an elementary reflection to a rectangular matrix of size MxN
//
// The  algorithm  pre-multiplies  the  matrix  by  an  elementary reflection
// transformation  which  is  given  by  column  V  and  scalar  Tau (see the
// description of the GenerateReflection). Not the whole matrix  but  only  a
// part of it is transformed (rows from M1 to M2, columns from N1 to N2). Only
// the elements of this submatrix are changed.
//
// Note: the matrix is multiplied by H, not by H'.   If  it  is  required  to
// multiply the matrix by H', it is necessary to pass Conj(Tau) instead of Tau.
//
// Inputs:
//     C       -   matrix to be transformed.
//     Tau     -   scalar defining transformation.
//     V       -   column defining transformation.
//                 Array whose index ranges within [1..M2-M1+1]
//     M1, M2  -   range of rows to be transformed.
//     N1, N2  -   range of columns to be transformed.
//     WORK    -   working array whose index goes from N1 to N2.
//
// Outputs:
//     C       -   the result of multiplying the input matrix C by the
//                 transformation matrix which is given by Tau and V.
//                 If N1 > N2 or M1 > M2, C is not modified.
//
// LAPACK auxiliary routine (version 3.0):
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void complexapplyreflectionfromtheleft(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work) {
   complex t;
   ae_int_t i;
   if (ae_c_eq_d(tau, 0.0) || n1 > n2 || m1 > m2) {
      return;
   }
// w = C^T * conj(v)
   for (i = n1; i <= n2; i++) {
      work->xC[i] = complex_from_i(0);
   }
   for (i = m1; i <= m2; i++) {
      t = conj(v->xC[i + 1 - m1]);
      ae_v_caddc(&work->xC[n1], 1, &c->xyC[i][n1], 1, "N", n2 - n1 + 1, t);
   }
// C = C - tau * v * w^T
   for (i = m1; i <= m2; i++) {
      t = ae_c_mul(v->xC[i - m1 + 1], tau);
      ae_v_csubc(&c->xyC[i][n1], 1, &work->xC[n1], 1, "N", n2 - n1 + 1, t);
   }
}

// Application of an elementary reflection to a rectangular matrix of size MxN
//
// The  algorithm  post-multiplies  the  matrix  by  an elementary reflection
// transformation  which  is  given  by  column  V  and  scalar  Tau (see the
// description  of  the  GenerateReflection). Not the whole matrix but only a
// part  of  it  is  transformed (rows from M1 to M2, columns from N1 to N2).
// Only the elements of this submatrix are changed.
//
// Inputs:
//     C       -   matrix to be transformed.
//     Tau     -   scalar defining transformation.
//     V       -   column defining transformation.
//                 Array whose index ranges within [1..N2-N1+1]
//     M1, M2  -   range of rows to be transformed.
//     N1, N2  -   range of columns to be transformed.
//     WORK    -   working array whose index goes from M1 to M2.
//
// Outputs:
//     C       -   the result of multiplying the input matrix C by the
//                 transformation matrix which is given by Tau and V.
//                 If N1 > N2 or M1 > M2, C is not modified.
//
// LAPACK auxiliary routine (version 3.0):
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void complexapplyreflectionfromtheright(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work) {
   complex t;
   ae_int_t i;
   ae_int_t vm;
   if (ae_c_eq_d(tau, 0.0) || n1 > n2 || m1 > m2) {
      return;
   }
// w = C * v
   vm = n2 - n1 + 1;
   for (i = m1; i <= m2; i++) {
      t = ae_v_cdotproduct(&c->xyC[i][n1], 1, "N", &v->xC[1], 1, "N", n2 - n1 + 1);
      work->xC[i] = t;
   }
// C = C - w * conj(v^T)
   ae_v_cmove(&v->xC[1], 1, &v->xC[1], 1, "Conj", vm);
   for (i = m1; i <= m2; i++) {
      t = ae_c_mul(work->xC[i], tau);
      ae_v_csubc(&c->xyC[i][n1], 1, &v->xC[1], 1, "N", n2 - n1 + 1, t);
   }
   ae_v_cmove(&v->xC[1], 1, &v->xC[1], 1, "Conj", vm);
}
} // end of namespace alglib_impl

// === SBLAS Package ===
// Depends on: APSERV
namespace alglib_impl {
void symmetricmatrixvectormultiply(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, double alpha, RVector *y) {
   ae_int_t i;
   ae_int_t ba1;
   ae_int_t by1;
   ae_int_t by2;
   ae_int_t bx1;
   ae_int_t bx2;
   ae_int_t n;
   double v;
   n = i2 - i1 + 1;
   if (n <= 0) {
      return;
   }
// Let A = L + D + U, where
//	L is strictly lower triangular (main diagonal is zero)
//	D is diagonal
//	U is strictly upper triangular (main diagonal is zero)
//
//	A*x = L*x + D*x + U*x
//
// Calculate D*x first
   for (i = i1; i <= i2; i++) {
      y->xR[i - i1 + 1] = a->xyR[i][i] * x->xR[i - i1 + 1];
   }
// Add L*x + U*x
   if (isupper) {
      for (i = i1; i < i2; i++) {
      // Add L*x to the result
         v = x->xR[i - i1 + 1];
         by1 = i - i1 + 2;
         by2 = n;
         ba1 = i + 1;
         ae_v_addd(&y->xR[by1], 1, &a->xyR[i][ba1], 1, by2 - by1 + 1, v);
      // Add U*x to the result
         bx1 = i - i1 + 2;
         bx2 = n;
         ba1 = i + 1;
         v = ae_v_dotproduct(&x->xR[bx1], 1, &a->xyR[i][ba1], 1, bx2 - bx1 + 1);
         y->xR[i - i1 + 1] += v;
      }
   } else {
      for (i = i1 + 1; i <= i2; i++) {
      // Add L*x to the result
         bx1 = 1;
         bx2 = i - i1;
         ba1 = i1;
         v = ae_v_dotproduct(&x->xR[bx1], 1, &a->xyR[i][ba1], 1, bx2 - bx1 + 1);
         y->xR[i - i1 + 1] += v;
      // Add U*x to the result
         v = x->xR[i - i1 + 1];
         by1 = 1;
         by2 = i - i1;
         ba1 = i1;
         ae_v_addd(&y->xR[by1], 1, &a->xyR[i][ba1], 1, by2 - by1 + 1, v);
      }
   }
   ae_v_muld(&y->xR[1], 1, n, alpha);
}

void symmetricrank2update(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, RVector *y, RVector *t, double alpha) {
   ae_int_t i;
   ae_int_t tp1;
   ae_int_t tp2;
   double v;
   if (isupper) {
      for (i = i1; i <= i2; i++) {
         tp1 = i + 1 - i1;
         tp2 = i2 - i1 + 1;
         v = x->xR[i + 1 - i1];
         ae_v_moved(&t->xR[tp1], 1, &y->xR[tp1], 1, tp2 - tp1 + 1, v);
         v = y->xR[i + 1 - i1];
         ae_v_addd(&t->xR[tp1], 1, &x->xR[tp1], 1, tp2 - tp1 + 1, v);
         ae_v_muld(&t->xR[tp1], 1, tp2 - tp1 + 1, alpha);
         ae_v_add(&a->xyR[i][i], 1, &t->xR[tp1], 1, i2 - i + 1);
      }
   } else {
      for (i = i1; i <= i2; i++) {
         tp1 = 1;
         tp2 = i + 1 - i1;
         v = x->xR[i + 1 - i1];
         ae_v_moved(&t->xR[tp1], 1, &y->xR[tp1], 1, tp2 - tp1 + 1, v);
         v = y->xR[i + 1 - i1];
         ae_v_addd(&t->xR[tp1], 1, &x->xR[tp1], 1, tp2 - tp1 + 1, v);
         ae_v_muld(&t->xR[tp1], 1, tp2 - tp1 + 1, alpha);
         ae_v_add(&a->xyR[i][i1], 1, &t->xR[tp1], 1, i - i1 + 1);
      }
   }
}
} // end of namespace alglib_impl

// === ABLASMKL Package ===
namespace alglib_impl {
// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool rmatrixgermkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double alpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixgermkl(m, n, a, ia, ja, alpha, u, iu, v, iv);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool rmatrixrank1mkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixrank1mkl(m, n, a, ia, ja, u, iu, v, iv);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool cmatrixrank1mkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixrank1mkl(m, n, a, ia, ja, u, iu, v, iv);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool rmatrixmvmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, RVector *y, ae_int_t iy) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixmvmkl(m, n, a, ia, ja, opa, x, ix, y, iy);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool cmatrixmvmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CVector *x, ae_int_t ix, CVector *y, ae_int_t iy) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixmvmkl(m, n, a, ia, ja, opa, x, ix, y, iy);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool rmatrixgemvmkl(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixgemvmkl(m, n, alpha, a, ia, ja, opa, x, ix, beta, y, iy);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 12.10.2017 by Sergey Bochkanov
bool rmatrixtrsvmkl(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t opa, RVector *x, ae_int_t ix) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixtrsvmkl(n, a, ia, ja, isupper, isunit, opa, x, ix);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2017 by Sergey Bochkanov
bool rmatrixsymvmkl(ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixsymvmkl(n, alpha, a, ia, ja, isupper, x, ix, beta, y, iy);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2013 by Sergey Bochkanov
bool rmatrixsyrkmkl(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixsyrkmkl(n, k, alpha, a, ia, ja, opa, beta, c, ic, jc, isupper);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2013 by Sergey Bochkanov
bool cmatrixherkmkl(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixherkmkl(n, k, alpha, a, ia, ja, opa, beta, c, ic, jc, isupper);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2013 by Sergey Bochkanov
bool rmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixgemmmkl(m, n, k, alpha, a, ia, ja, opa, b, ib, jb, opb, beta, c, ic, jc);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool cmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixgemmmkl(m, n, k, alpha, a, ia, ja, opa, b, ib, jb, opb, beta, c, ic, jc);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool rmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, RMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool cmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, CMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool rmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, RMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool cmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, CMatrix *x, ae_int_t i2, ae_int_t j2) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, opa, x, i2, j2);
#endif
}

// MKL-based kernel.
//
// NOTE:
//
// if function returned False, CholResult is NOT modified. Not ever referenced!
// if function returned True, CholResult is set to status of Cholesky decomposition
// (True on succeess).
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool spdmatrixcholeskymkl(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, bool *cholresult) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_spdmatrixcholeskymkl(a, offs, n, isupper, cholresult);
#endif
}

// MKL-based kernel.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixplumkl(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixplumkl(a, offs, m, n, pivots);
#endif
}

// MKL-based kernel.
//
// NOTE: this function needs preallocated output/temporary arrays.
//       D and E must be at least max(M,N)-wide.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixbdmkl(RMatrix *a, ae_int_t m, ae_int_t n, RVector *d, RVector *e, RVector *tauq, RVector *taup) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixbdmkl(a, m, n, d, e, tauq, taup);
#endif
}

// MKL-based kernel.
//
// If ByQ is True,  TauP is not used (can be empty array).
// If ByQ is False, TauQ is not used (can be empty array).
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixbdmultiplybymkl(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, RVector *taup, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool byq, bool fromtheright, bool dotranspose) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixbdmultiplybymkl(qp, m, n, tauq, taup, z, zrows, zcolumns, byq, fromtheright, dotranspose);
#endif
}

// MKL-based kernel.
//
// NOTE: Tau must be preallocated array with at least N-1 elements.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixhessenbergmkl(RMatrix *a, ae_int_t n, RVector *tau) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixhessenbergmkl(a, n, tau);
#endif
}

// MKL-based kernel.
//
// NOTE: Q must be preallocated N*N array
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixhessenbergunpackqmkl(RMatrix *a, ae_int_t n, RVector *tau, RMatrix *q) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixhessenbergunpackqmkl(a, n, tau, q);
#endif
}

// MKL-based kernel.
//
// NOTE: Tau, D, E must be preallocated arrays;
//       length(E) == length(Tau) == N-1 (or larger)
//       length(D) == N (or larger)
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool smatrixtdmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RVector *d, RVector *e) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_smatrixtdmkl(a, n, isupper, tau, d, e);
#endif
}

// MKL-based kernel.
//
// NOTE: Tau, D, E must be preallocated arrays;
//       length(E) == length(Tau) == N-1 (or larger)
//       length(D) == N (or larger)
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool hmatrixtdmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, RVector *d, RVector *e) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hmatrixtdmkl(a, n, isupper, tau, d, e);
#endif
}

// MKL-based kernel.
//
// NOTE: Q must be preallocated N*N array
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool smatrixtdunpackqmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RMatrix *q) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_smatrixtdunpackqmkl(a, n, isupper, tau, q);
#endif
}

// MKL-based kernel.
//
// NOTE: Q must be preallocated N*N array
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool hmatrixtdunpackqmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, CMatrix *q) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hmatrixtdunpackqmkl(a, n, isupper, tau, q);
#endif
}

// MKL-based kernel.
//
// Returns True if MKL was present and handled request (MKL  completion  code
// is returned as separate output parameter).
//
// D and E are preallocated arrays with length N (both of them!). On output,
// D constraints singular values, and E is destroyed.
//
// SVDResult is modified if and only if MKL is present.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixbdsvdmkl(RVector *d, RVector *e, ae_int_t n, bool isupper, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt, bool *svdresult) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixbdsvdmkl(d, e, n, isupper, u, nru, c, ncc, vt, ncvt, svdresult);
#endif
}

// MKL-based DHSEQR kernel.
//
// Returns True if MKL was present and handled request.
//
// WR and WI are preallocated arrays with length N.
// Z is preallocated array[N,N].
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixinternalschurdecompositionmkl(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixinternalschurdecompositionmkl(h, n, tneeded, zneeded, wr, wi, z, info);
#endif
}

// MKL-based DTREVC kernel.
//
// Returns True if MKL was present and handled request.
//
// NOTE: this function does NOT support HOWMNY == 3!!!!
//
// VL and VR are preallocated arrays with length N*N, if required. If particalar
// variables is not required, it can be dummy (empty) array.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixinternaltrevcmkl(RMatrix *t, ae_int_t n, ae_int_t side, ae_int_t howmny, RMatrix *vl, RMatrix *vr, ae_int_t *m, ae_int_t *info) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixinternaltrevcmkl(t, n, side, howmny, vl, vr, m, info);
#endif
}

// MKL-based kernel.
//
// Returns True if MKL was present and handled request (MKL  completion  code
// is returned as separate output parameter).
//
// D and E are preallocated arrays with length N (both of them!). On output,
// D constraints eigenvalues, and E is destroyed.
//
// Z is preallocated array[N,N] for ZNeeded != 0; ignored for ZNeeded == 0.
//
// EVDResult is modified if and only if MKL is present.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool smatrixtdevdmkl(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z, bool *evdresult) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_smatrixtdevdmkl(d, e, n, zneeded, z, evdresult);
#endif
}

// MKL-based kernel.
//
// Returns True if MKL was present and handled request (MKL  completion  code
// is returned as separate output parameter).
//
// D and E are preallocated arrays with length N (both of them!). On output,
// D constraints eigenvalues, and E is destroyed.
//
// Z is preallocated array[N,N] for ZNeeded != 0; ignored for ZNeeded == 0.
//
// EVDResult is modified if and only if MKL is present.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool sparsegemvcrsmkl(ae_int_t opa, ae_int_t arows, ae_int_t acols, double alpha, RVector *vals, ZVector *cidx, ZVector *ridx, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_sparsegemvcrsmkl(opa, arows, acols, alpha, vals, cidx, ridx, x, ix, beta, y, iy);
#endif
}
} // end of namespace alglib_impl

// === SCODES Package ===
namespace alglib_impl {
ae_int_t getrdfserializationcode() {
   ae_int_t result;
   result = 1;
   return result;
}

ae_int_t getkdtreeserializationcode() {
   ae_int_t result;
   result = 2;
   return result;
}

ae_int_t getmlpserializationcode() {
   ae_int_t result;
   result = 3;
   return result;
}

ae_int_t getmlpeserializationcode() {
   ae_int_t result;
   result = 4;
   return result;
}

ae_int_t getrbfserializationcode() {
   ae_int_t result;
   result = 5;
   return result;
}

ae_int_t getspline2dserializationcode() {
   ae_int_t result;
   result = 6;
   return result;
}

ae_int_t getidwserializationcode() {
   ae_int_t result;
   result = 7;
   return result;
}

ae_int_t getsparsematrixserializationcode() {
   ae_int_t result;
   result = 8;
   return result;
}

ae_int_t getspline2dwithmissingnodesserializationcode() {
   ae_int_t result;
   result = 9;
   return result;
}

ae_int_t getknnserializationcode() {
   ae_int_t result;
   result = 108;
   return result;
}

ae_int_t getlptestserializationcode() {
   ae_int_t result;
   result = 200;
   return result;
}
} // end of namespace alglib_impl

// === TSORT Package ===
// Depends on: APSERV
namespace alglib_impl {
// Internal versions respectively of tagsortfasti(), tagsortfastr() and tagsortfast():
// sort the n-vector ap and apply the same permutations to abp (and to bp and bbp, if present).
// ALGLIB: Copyright 06.09.2010 by Sergey Bochkanov
static void tsort_tagsortfastirec(double *ap, ae_int_t *bp, double *abp, ae_int_t *bbp, ae_int_t n) {
   while (n > 16) {
   // Quicksort: choose the pivot.
      double al = ap[0], am = ap[(n - 1) / 2], ah = ap[n - 1];
      if (al > am) swapr(&al, &am);
      if (am > ah) swapr(&am, &ah);
      if (al > am) swapr(&al, &am);
   // Pass through ap/bp:
   // * everything < am goes to the left of ap/bp,
   // * everything == am goes to the right of abp/bbp (in the reverse order),
   // * everything > am goes to the left of abp/bbp (in the normal order,
   // * everything from the right of abp/bbp to the middle goes to the right of ap/bp back in normal order,
   // * everything from the left of abp/bbp to the middle goes to the right of ap/bp.
      ae_int_t lts = 0, eqs = 0, gts = 0;
      for (ae_int_t i = 0; i < n; i++) {
         double ai = ap[i];
         ae_int_t bi = bp[i];
         ae_int_t k = ai < am ? lts++ : ai > am ? gts++ : n - ++eqs;
         if (ai >= am)
            abp[k] = ai, bbp[k] = bi;
         else if (i != k)
            ap[k] = ai, bp[k] = bi;
      }
      ae_int_t les = lts + eqs;
      for (ae_int_t k = n - eqs, j = les - 1; k < n; k++, j--)
         ap[j] = abp[k], bp[j] = bbp[k];
      for (ae_int_t k = 0, j = les; k < gts; k++, j++)
         ap[j] = abp[k], bp[j] = bbp[k];
   // Sort the left part of the vector, and skip over the middle part.
      tsort_tagsortfastirec(ap, bp, abp, bbp, lts);
      ap += les, bp += les, abp += les, bbp += les, n -= les;
   }
// The non-recursive sort for small arrays.
   for (ae_int_t j = 1; j < n; j++) {
   // Find a place, if any, in [0,j) to insert ap[j] at.
   // Skip of ap[j] can be left intact; i.e. if all elements have the same value of ap[j] or larger than any of them.
      double aj = ap[j];
      ae_int_t k = j;
      while (--k >= 0 && ap[k] > aj);
      if (++k >= j) continue;
   // Move element j down to position k.
      ae_int_t bj = bp[j];
      for (ae_int_t i = j; i > k; i--) ap[i] = ap[i - 1], bp[i] = bp[i - 1];
      ap[k] = aj, bp[k] = bj;
   }
}

static void tsort_tagsortfastrrec(double *ap, double *bp, double *abp, double *bbp, ae_int_t n) {
   while (n > 16) {
   // Quicksort: choose the pivot.
      double al = ap[0], am = ap[(n - 1) / 2], ah = ap[n - 1];
      if (al > am) swapr(&al, &am);
      if (am > ah) swapr(&am, &ah);
      if (al > am) swapr(&al, &am);
   // Pass through ap/bp:
   // * everything < am goes to the left of ap/bp,
   // * everything == am goes to the right of abp/bbp (in the reverse order),
   // * everything > am goes to the left of abp/bbp (in the normal order,
   // * everything from the right of abp/bbp to the middle goes to the right of ap/bp back in normal order,
   // * everything from the left of abp/bbp to the middle goes to the right of ap/bp.
      ae_int_t lts = 0, eqs = 0, gts = 0;
      for (ae_int_t i = 0; i < n; i++) {
         double ai = ap[i];
         ae_int_t k = ai < am ? lts++ : ai > am ? gts++ : n - ++eqs;
         if (ai >= am)
            abp[k] = ai, bbp[k] = bp[i];
         else if (i != k)
            ap[k] = ai, bp[k] = bp[i];
      }
      ae_int_t les = lts + eqs;
      for (ae_int_t k = n - eqs, j = les - 1; k < n; k++, j--)
         ap[j] = abp[k], bp[j] = bbp[k];
      for (ae_int_t k = 0, j = les; k < gts; k++, j++)
         ap[j] = abp[k], bp[j] = bbp[k];
   // Sort the left part of the vector, and skip over the middle part.
      tsort_tagsortfastrrec(ap, bp, abp, bbp, lts);
      ap += les, bp += les, abp += les, bbp += les, n -= les;
   }
// The non-recursive sort for small arrays.
   for (ae_int_t j = 1; j < n; j++) {
   // Find a place, if any, in [0,j) to insert ap[j] at.
   // Skip of ap[j] can be left intact; i.e. if all elements have the same value of ap[j] or larger than any of them.
      double aj = ap[j];
      ae_int_t k = j;
      while (--k >= 0 && ap[k] > aj);
      if (++k >= j) continue;
   // Move element j down to position k.
      double bj = bp[j];
      for (ae_int_t i = j; i > k; i--) ap[i] = ap[i - 1], bp[i] = bp[i - 1];
      ap[k] = aj, bp[k] = bj;
   }
}

static void tsort_tagsortfastrec(double *ap, double *abp, ae_int_t n) {
   while (n > 16) {
   // Quicksort: choose the pivot.
      double al = ap[0], am = ap[(n - 1) / 2], ah = ap[n - 1];
      if (al > am) swapr(&al, &am);
      if (am > ah) swapr(&am, &ah);
      if (al > am) swapr(&al, &am);
   // Pass through ap:
   // * everything < am goes to the left of ap,
   // * everything == am goes to the right of abp (in the reverse order),
   // * everything > am goes to the left of abp (in the normal order,
   // * everything from the right of abp to the middle goes to the right of ap back in normal order,
   // * everything from the left of abp to the middle goes to the right of ap.
      ae_int_t lts = 0, eqs = 0, gts = 0;
      for (ae_int_t i = 0; i < n; i++) {
         double ai = ap[i];
         ae_int_t k = ai < am ? lts++ : ai > am ? gts++ : n - ++eqs;
         if (ai >= am)
            abp[k] = ai;
         else if (i != k)
            ap[k] = ai;
      }
      ae_int_t les = lts + eqs;
      for (ae_int_t k = n - eqs, j = les - 1; k < n; k++, j--)
         ap[j] = abp[k];
      for (ae_int_t k = 0, j = les; k < gts; k++, j++)
         ap[j] = abp[k];
   // Sort the left part of the vector, and skip over the middle part.
      tsort_tagsortfastrec(ap, abp, lts);
      ap += les, abp += les, n -= les;
   }
// The non-recursive sort for small arrays.
   for (ae_int_t j = 1; j < n; j++) {
   // Find a place, if any, in [0,j) to insert ap[j] at.
   // Skip of ap[j] can be left intact; i.e. if all elements have the same value of ap[j] or larger than any of them.
      double aj = ap[j];
      ae_int_t k = j;
      while (--k >= 0 && ap[k] > aj);
      if (++k >= j) continue;
   // Move element j down to position k.
      for (ae_int_t i = j; i > k; i--) ap[i] = ap[i - 1];
      ap[k] = aj;
   }
}

// Versions of tagsort() optimized for real keys and {integer, real, no} labels respectively.
// The vector a is sorted, and same permutations are applied to the vector b, if present.
// NOTES:
// 1.	these functions assume that the vector a is finite
//	and don't check for it or for other conditions (size of input arrays, etc.).
// 2.	these functions use buffer(s), bufa and bufb, that will be reallocated to size n, if necessary.
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortfasti(RVector *a, ZVector *b, RVector *bufa, ZVector *bufb, ae_int_t n) {
// Special case.
   if (n <= 1) return;
// Test for already sorted set
   bool isascending = true, isdescending = true;
   for (ae_int_t i = 1; i < n; i++) {
      isascending = isascending && a->xR[i] >= a->xR[i - 1];
      isdescending = isdescending && a->xR[i] <= a->xR[i - 1];
   }
   if (isascending) return;
   else if (isdescending) {
      for (ae_int_t i = 0; i < n; i++) {
         ae_int_t j = n - 1 - i;
         if (j <= i) break;
         swapr(&a->xR[i], &a->xR[j]);
         swapi(&b->xZ[i], &b->xZ[j]);
      }
   } else {
   // General case.
      if (bufa->cnt < n) ae_vector_set_length(bufa, n);
      if (bufb->cnt < n) ae_vector_set_length(bufb, n);
      tsort_tagsortfastirec(a->xR, b->xZ, bufa->xR, bufb->xZ, n);
   }
}

void tagsortfastr(RVector *a, RVector *b, RVector *bufa, RVector *bufb, ae_int_t n) {
// Special case.
   if (n <= 1) return;
// Test for already sorted set
   bool isascending = true, isdescending = true;
   for (ae_int_t i = 1; i < n; i++) {
      isascending = isascending && a->xR[i] >= a->xR[i - 1];
      isdescending = isdescending && a->xR[i] <= a->xR[i - 1];
   }
   if (isascending) return;
   else if (isdescending) {
      for (ae_int_t i = 0; i < n; i++) {
         ae_int_t j = n - 1 - i;
         if (j <= i) break;
         swapr(&a->xR[i], &a->xR[j]);
         swapr(&b->xR[i], &b->xR[j]);
      }
   } else {
   // General case.
      if (bufa->cnt < n) ae_vector_set_length(bufa, n);
      if (bufb->cnt < n) ae_vector_set_length(bufb, n);
      tsort_tagsortfastrrec(a->xR, b->xR, bufa->xR, bufb->xR, n);
   }
}

void tagsortfast(RVector *a, RVector *bufa, ae_int_t n) {
// Special case.
   if (n <= 1) return;
// Test for already sorted set
   bool isascending = true, isdescending = true;
   for (ae_int_t i = 1; i < n; i++) {
      isascending = isascending && a->xR[i] >= a->xR[i - 1];
      isdescending = isdescending && a->xR[i] <= a->xR[i - 1];
   }
   if (isascending) return;
   else if (isdescending) {
      for (ae_int_t i = 0; i < n; i++) {
         ae_int_t j = n - 1 - i;
         if (j <= i) break;
         swapr(&a->xR[i], &a->xR[j]);
      }
   } else {
   // General case.
      if (bufa->cnt < n) ae_vector_set_length(bufa, n);
      tsort_tagsortfastrec(a->xR, bufa->xR, n);
   }
}

// Versions of tagsort*(), optimized for integer keys optionally with integer or real labels,
// suitable for sorting the middle of the vector.
// Sort the vector a, starting at offset, and apply the same permutations to the vector b, if it is specified.
// NOTE:
//	These functions assume that the vector a is finite
//	and do not check for it or for other conditions (size of input arrays, etc.).
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortmiddleii(ZVector *a, ZVector *b, ae_int_t n, ae_int_t offset/* = 0*/) {
// Special case.
   if (n <= 1) return;
// General case, n > 1: sort, update bp.
   ae_int_t *ap = a->xZ + offset;
   ae_int_t *bp = b->xZ + offset;
   for (ae_int_t nh = 1; nh < n; nh++)
      for (ae_int_t nm = nh, nl = (nh - 1) / 2; nm > 0 && ap[nl] < ap[nm]; nm = nl, nl = (nl - 1) / 2) {
         swapi(ap + nl, ap + nm);
         swapi(bp + nl, bp + nm);
      }
   for (ae_int_t nh = n - 1; nh > 0; nh--) {
      swapi(ap, ap + nh);
      swapi(bp, bp + nh);
      for (ae_int_t nl = 0, nm = 1; nm < nh; nl = nm, nm = 2 * nm + 1) {
         if (nm + 1 < nh && ap[nm + 1] > ap[nm]) nm++;
         if (ap[nl] >= ap[nm]) break;
         swapi(ap + nl, ap + nm);
         swapi(bp + nl, bp + nm);
      }
   }
}

void tagsortmiddleir(ZVector *a, RVector *b, ae_int_t n, ae_int_t offset/* = 0*/) {
// Special case.
   if (n <= 1) return;
// General case, n > 1: sort, update bp.
   ae_int_t *ap = a->xZ + offset;
   double *bp = b->xR + offset;
   for (ae_int_t nh = 1; nh < n; nh++)
      for (ae_int_t nm = nh, nl = (nh - 1) / 2; nm > 0 && ap[nl] < ap[nm]; nm = nl, nl = (nl - 1) / 2) {
         swapi(ap + nl, ap + nm);
         swapr(bp + nl, bp + nm);
      }
   for (ae_int_t nh = n - 1; nh > 0; nh--) {
      swapi(ap, ap + nh);
      swapr(bp, bp + nh);
      for (ae_int_t nl = 0, nm = 1; nm < nh; nl = nm, nm = 2 * nm + 1) {
         if (nm + 1 < nh && ap[nm + 1] > ap[nm]) nm++;
         if (ap[nl] >= ap[nm]) break;
         swapi(ap + nl, ap + nm);
         swapr(bp + nl, bp + nm);
      }
   }
}

void tagsortmiddleri(RVector *a, ZVector *b, ae_int_t n, ae_int_t offset/* = 0*/) {
// Special case.
   if (n <= 1) return;
// General case, n > 1: sort, update bp.
   double *ap = a->xR + offset;
   ae_int_t *bp = b->xZ + offset;
   for (ae_int_t nh = 1; nh < n; nh++)
      for (ae_int_t nm = nh, nl = (nh - 1) / 2; nm > 0 && ap[nl] < ap[nm]; nm = nl, nl = (nl - 1) / 2) {
         swapr(ap + nl, ap + nm);
         swapi(bp + nl, bp + nm);
      }
   for (ae_int_t nh = n - 1; nh > 0; nh--) {
      swapr(ap, ap + nh);
      swapi(bp, bp + nh);
      for (ae_int_t nl = 0, nm = 1; nm < nh; nl = nm, nm = 2 * nm + 1) {
         if (nm + 1 < nh && ap[nm + 1] > ap[nm]) nm++;
         if (ap[nl] >= ap[nm]) break;
         swapr(ap + nl, ap + nm);
         swapi(bp + nl, bp + nm);
      }
   }
}

void tagsortmiddlei(ZVector *a, ae_int_t n, ae_int_t offset/* = 0*/) {
// Special case.
   if (n <= 1) return;
// General case, n > 1: sort.
   ae_int_t *ap = a->xZ + offset;
   for (ae_int_t nh = 1; nh < n; nh++)
      for (ae_int_t nm = nh, nl = (nh - 1) / 2; nm > 0 && ap[nl] < ap[nm]; nm = nl, nl = (nl - 1) / 2)
         swapi(ap + nl, ap + nm);
   for (ae_int_t nh = n - 1; nh > 0; nh--) {
      swapi(ap, ap + nh);
      for (ae_int_t nl = 0, nm = 1; nm < nh; nl = nm, nm = 2 * nm + 1) {
         if (nm + 1 < nh && ap[nm + 1] > ap[nm]) nm++;
         if (ap[nl] >= ap[nm]) break;
         swapi(ap + nl, ap + nm);
      }
   }
}

// Buffered variant of TagSort, which accepts preallocated output arrays as
// well as special structure for buffered allocations. If arrays are too
// short, they are reallocated. If they are large enough, no memory
// allocation is done.
//
// It is intended to be used in the performance-critical parts of code, where
// additional allocations can lead to severe performance degradation
// ALGLIB: Copyright 14.05.2008 by Sergey Bochkanov
void tagsortbuf(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2, apbuffers *buf) {
   ae_int_t i;
   ae_int_t lv;
   ae_int_t lp;
   ae_int_t rv;
   ae_int_t rp;
// Special cases.
   if (n <= 0) {
      return;
   }
   if (n == 1) {
      vectorsetlengthatleast(p1, 1);
      vectorsetlengthatleast(p2, 1);
      p1->xZ[0] = 0;
      p2->xZ[0] = 0;
      return;
   }
// General case, N > 1: prepare permutations table P1
   vectorsetlengthatleast(p1, n);
   for (i = 0; i < n; i++) {
      p1->xZ[i] = i;
   }
// General case, N > 1: sort, update P1
   vectorsetlengthatleast(&buf->ra0, n);
   vectorsetlengthatleast(&buf->ia0, n);
   tagsortfasti(a, p1, &buf->ra0, &buf->ia0, n);
// General case, N > 1: fill permutations table P2
//
// To fill P2 we maintain two arrays:
// * PV (Buf.IA0), Position(Value). PV[i] contains position of I-th key at the moment
// * VP (Buf.IA1), Value(Position). VP[i] contains key which has position I at the moment
//
// At each step we making permutation of two items:
//   Left, which is given by position/value pair LP/LV
//   and Right, which is given by RP/RV
// and updating PV[] and VP[] correspondingly.
   vectorsetlengthatleast(&buf->ia0, n);
   vectorsetlengthatleast(&buf->ia1, n);
   vectorsetlengthatleast(p2, n);
   for (i = 0; i < n; i++) {
      buf->ia0.xZ[i] = i;
      buf->ia1.xZ[i] = i;
   }
   for (i = 0; i < n; i++) {
   // calculate LP, LV, RP, RV
      lp = i;
      lv = buf->ia1.xZ[lp];
      rv = p1->xZ[i];
      rp = buf->ia0.xZ[rv];
   // Fill P2
      p2->xZ[i] = rp;
   // update PV and VP
      buf->ia1.xZ[lp] = rv;
      buf->ia1.xZ[rp] = lv;
      buf->ia0.xZ[lv] = rp;
      buf->ia0.xZ[rv] = lp;
   }
}

// This function sorts array of real keys by ascending.
//
// Its results are:
// * sorted array A
// * permutation tables P1, P2
//
// Algorithm outputs permutation tables using two formats:
// * as usual permutation of [0..N-1]. If P1[i] == j, then sorted A[i]  contains
//   value which was moved there from J-th position.
// * as a sequence of pairwise permutations. Sorted A[] may  be  obtained  by
//   swaping A[i] and A[P2[i]] for all i from 0 to N-1.
//
// Inputs:
//     A       -   unsorted array
//     N       -   array size
//
// Outputs:
//     A       -   sorted array
//     P1, P2  -   permutation tables, array[N]
//
// NOTES:
//     this function assumes that A[] is finite; it doesn't checks that
//     condition. All other conditions (size of input arrays, etc.) are not
//     checked too.
// ALGLIB: Copyright 14.05.2008 by Sergey Bochkanov
void tagsort(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   SetVector(p1);
   SetVector(p2);
   NewObj(apbuffers, buf);
   tagsortbuf(a, n, p1, p2, &buf);
   ae_frame_leave();
}

// Heap operations: adds element to the heap
//
// Parameters:
//     A       -   heap itself, must be at least array[0..N]
//     B       -   array of integer tags, which are updated according to
//                 permutations in the heap
//     N       -   size of the heap (without new element).
//                 updated on output
//     VA      -   value of the element being added
//     VB      -   value of the tag
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void tagheappushi(RVector *a, ZVector *b, ae_int_t *n, double va, ae_int_t vb) {
   ae_int_t j;
   ae_int_t k;
   double v;
   if (*n < 0) {
      return;
   }
// N == 0 is a special case
   if (*n == 0) {
      a->xR[0] = va;
      b->xZ[0] = vb;
      ++*n;
      return;
   }
// add current point to the heap
// (add to the bottom, then move up)
//
// we don't write point to the heap
// until its final position is determined
// (it allow us to reduce number of array access operations)
   j = *n;
   ++*n;
   while (j > 0) {
      k = (j - 1) / 2;
      v = a->xR[k];
      if (v < va) {
      // swap with higher element
         a->xR[j] = v;
         b->xZ[j] = b->xZ[k];
         j = k;
      } else {
      // element in its place. terminate.
         break;
      }
   }
   a->xR[j] = va;
   b->xZ[j] = vb;
}

// Heap operations: replaces top element with new element
// (which is moved down)
//
// Parameters:
//     A       -   heap itself, must be at least array[0..N-1]
//     B       -   array of integer tags, which are updated according to
//                 permutations in the heap
//     N       -   size of the heap
//     VA      -   value of the element which replaces top element
//     VB      -   value of the tag
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void tagheapreplacetopi(RVector *a, ZVector *b, ae_int_t n, double va, ae_int_t vb) {
   ae_int_t j;
   ae_int_t k1;
   ae_int_t k2;
   double v;
   double v1;
   double v2;
   if (n < 1) {
      return;
   }
// N == 1 is a special case
   if (n == 1) {
      a->xR[0] = va;
      b->xZ[0] = vb;
      return;
   }
// move down through heap:
// * J  -   current element
// * K1 -   first child (always exists)
// * K2 -   second child (may not exists)
//
// we don't write point to the heap
// until its final position is determined
// (it allow us to reduce number of array access operations)
   j = 0;
   k1 = 1;
   k2 = 2;
   while (k1 < n) {
      if (k2 >= n) {
      // only one child.
      //
      // swap and terminate (because this child
      // have no siblings due to heap structure)
         v = a->xR[k1];
         if (v > va) {
            a->xR[j] = v;
            b->xZ[j] = b->xZ[k1];
            j = k1;
         }
         break;
      } else {
      // two childs
         v1 = a->xR[k1];
         v2 = a->xR[k2];
         if (v1 > v2) {
            if (va < v1) {
               a->xR[j] = v1;
               b->xZ[j] = b->xZ[k1];
               j = k1;
            } else {
               break;
            }
         } else {
            if (va < v2) {
               a->xR[j] = v2;
               b->xZ[j] = b->xZ[k2];
               j = k2;
            } else {
               break;
            }
         }
         k1 = 2 * j + 1;
         k2 = 2 * j + 2;
      }
   }
   a->xR[j] = va;
   b->xZ[j] = vb;
}

// Heap operations: pops top element from the heap
//
// Parameters:
//     A       -   heap itself, must be at least array[0..N-1]
//     B       -   array of integer tags, which are updated according to
//                 permutations in the heap
//     N       -   size of the heap, N >= 1
//
// On output top element is moved to A[N-1], B[N-1], heap is reordered, N is
// decreased by 1.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void tagheappopi(RVector *a, ZVector *b, ae_int_t *n) {
   double va;
   ae_int_t vb;
   if (*n < 1) {
      return;
   }
// N == 1 is a special case
   if (*n == 1) {
      *n = 0;
      return;
   }
// swap top element and last element,
// then reorder heap
   va = a->xR[*n - 1];
   vb = b->xZ[*n - 1];
   a->xR[*n - 1] = a->xR[0];
   b->xZ[*n - 1] = b->xZ[0];
   --*n;
   tagheapreplacetopi(a, b, *n, va, vb);
}

// Search first element less than T in sorted array.
//
// Parameters:
//     A - sorted array by ascending from 0 to N-1
//     N - number of elements in array
//     T - the desired element
//
// Result:
//     The very first element's index, which isn't less than T.
// In the case when there aren't such elements, returns N.
ae_int_t lowerbound(RVector *a, ae_int_t n, double t) {
   ae_int_t l;
   ae_int_t half;
   ae_int_t first;
   ae_int_t middle;
   ae_int_t result;
   l = n;
   first = 0;
   while (l > 0) {
      half = l / 2;
      middle = first + half;
      if (a->xR[middle] < t) {
         first = middle + 1;
         l -= half + 1;
      } else {
         l = half;
      }
   }
   result = first;
   return result;
}

// Search first element more than T in sorted array.
//
// Parameters:
//     A - sorted array by ascending from 0 to N-1
//     N - number of elements in array
//     T - the desired element
//
// Result:
//     The very first element's index, which more than T.
// In the case when there aren't such elements, returns N.
ae_int_t upperbound(RVector *a, ae_int_t n, double t) {
   ae_int_t l;
   ae_int_t half;
   ae_int_t first;
   ae_int_t middle;
   ae_int_t result;
   l = n;
   first = 0;
   while (l > 0) {
      half = l / 2;
      middle = first + half;
      if (t < a->xR[middle]) {
         l = half;
      } else {
         first = middle + 1;
         l -= half + 1;
      }
   }
   result = first;
   return result;
}
} // end of namespace alglib_impl

// === BLAS Package ===
namespace alglib_impl {
double vectornorm2(RVector *x, ae_int_t i1, ae_int_t i2) {
   ae_int_t n;
   ae_int_t ix;
   double absxi;
   double scl;
   double ssq;
   double result;
   n = i2 - i1 + 1;
   if (n < 1) {
      result = 0.0;
      return result;
   }
   if (n == 1) {
      result = fabs(x->xR[i1]);
      return result;
   }
   scl = 0.0;
   ssq = 1.0;
   for (ix = i1; ix <= i2; ix++) {
      if (x->xR[ix] != 0.0) {
         absxi = fabs(x->xR[ix]);
         if (scl < absxi) {
            ssq = 1.0 + ssq * sqr(scl / absxi);
            scl = absxi;
         } else {
            ssq += sqr(absxi / scl);
         }
      }
   }
   result = scl * sqrt(ssq);
   return result;
}

ae_int_t vectoridxabsmax(RVector *x, ae_int_t i1, ae_int_t i2) {
   ae_int_t i;
   ae_int_t result;
   result = i1;
   for (i = i1 + 1; i <= i2; i++) {
      if (fabs(x->xR[i]) > fabs(x->xR[result])) {
         result = i;
      }
   }
   return result;
}

ae_int_t columnidxabsmax(RMatrix *x, ae_int_t i1, ae_int_t i2, ae_int_t j) {
   ae_int_t i;
   ae_int_t result;
   result = i1;
   for (i = i1 + 1; i <= i2; i++) {
      if (fabs(x->xyR[i][j]) > fabs(x->xyR[result][j])) {
         result = i;
      }
   }
   return result;
}

ae_int_t rowidxabsmax(RMatrix *x, ae_int_t j1, ae_int_t j2, ae_int_t i) {
   ae_int_t j;
   ae_int_t result;
   result = j1;
   for (j = j1 + 1; j <= j2; j++) {
      if (fabs(x->xyR[i][j]) > fabs(x->xyR[i][result])) {
         result = j;
      }
   }
   return result;
}

double upperhessenberg1norm(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work) {
   ae_int_t i;
   ae_int_t j;
   double result;
   ae_assert(i2 - i1 == j2 - j1, "upperhessenberg1norm: i2 - i1 != j2 - j1!");
   for (j = j1; j <= j2; j++) {
      work->xR[j] = 0.0;
   }
   for (i = i1; i <= i2; i++) {
      for (j = imax2(j1, j1 + i - i1 - 1); j <= j2; j++) {
         work->xR[j] += fabs(a->xyR[i][j]);
      }
   }
   result = 0.0;
   for (j = j1; j <= j2; j++) {
      result = rmax2(result, work->xR[j]);
   }
   return result;
}

void copymatrix(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2) {
   ae_int_t isrc;
   ae_int_t idst;
   if (is1 > is2 || js1 > js2) {
      return;
   }
   ae_assert(is2 - is1 == id2 - id1, "copymatrix: different sizes!");
   ae_assert(js2 - js1 == jd2 - jd1, "copymatrix: different sizes!");
   for (isrc = is1; isrc <= is2; isrc++) {
      idst = isrc - is1 + id1;
      ae_v_move(&b->xyR[idst][jd1], 1, &a->xyR[isrc][js1], 1, jd2 - jd1 + 1);
   }
}

void inplacetranspose(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t ips;
   ae_int_t jps;
   ae_int_t l;
   if (i1 > i2 || j1 > j2) {
      return;
   }
   ae_assert(i1 - i2 == j1 - j2, "inplacetranspose: error: incorrect array size!");
   for (i = i1; i < i2; i++) {
      j = j1 + i - i1;
      ips = i + 1;
      jps = j1 + ips - i1;
      l = i2 - i;
      ae_v_move(&work->xR[1], 1, &a->xyR[ips][j], a->stride, l);
      ae_v_move(&a->xyR[ips][j], a->stride, &a->xyR[i][jps], 1, i2 - ips + 1);
      ae_v_move(&a->xyR[i][jps], 1, &work->xR[1], 1, j2 - jps + 1);
   }
}

void copyandtranspose(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2) {
   ae_int_t isrc;
   ae_int_t jdst;
   if (is1 > is2 || js1 > js2) {
      return;
   }
   ae_assert(is2 - is1 == jd2 - jd1, "copyandtranspose: different sizes!");
   ae_assert(js2 - js1 == id2 - id1, "copyandtranspose: different sizes!");
   for (isrc = is1; isrc <= is2; isrc++) {
      jdst = isrc - is1 + jd1;
      ae_v_move(&b->xyR[id1][jdst], b->stride, &a->xyR[isrc][js1], 1, id2 - id1 + 1);
   }
}

void matrixvectormultiply(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, bool trans, RVector *x, ae_int_t ix1, ae_int_t ix2, double alpha, RVector *y, ae_int_t iy1, ae_int_t iy2, double beta) {
   ae_int_t i;
   double v;
   if (!trans) {
   // y = alpha*A*x + beta*y;
      if (i1 > i2 || j1 > j2) {
         return;
      }
      ae_assert(j2 - j1 == ix2 - ix1, "matrixvectormultiply: a and x dont match!");
      ae_assert(i2 - i1 == iy2 - iy1, "matrixvectormultiply: a and y dont match!");
   // beta*y
      if (beta == 0.0) {
         for (i = iy1; i <= iy2; i++) {
            y->xR[i] = 0.0;
         }
      } else {
         ae_v_muld(&y->xR[iy1], 1, iy2 - iy1 + 1, beta);
      }
   // alpha*A*x
      for (i = i1; i <= i2; i++) {
         v = ae_v_dotproduct(&a->xyR[i][j1], 1, &x->xR[ix1], 1, j2 - j1 + 1);
         y->xR[iy1 + i - i1] += alpha * v;
      }
   } else {
   // y = alpha*A'*x + beta*y;
      if (i1 > i2 || j1 > j2) {
         return;
      }
      ae_assert(i2 - i1 == ix2 - ix1, "matrixvectormultiply: a and x dont match!");
      ae_assert(j2 - j1 == iy2 - iy1, "matrixvectormultiply: a and y dont match!");
   // beta*y
      if (beta == 0.0) {
         for (i = iy1; i <= iy2; i++) {
            y->xR[i] = 0.0;
         }
      } else {
         ae_v_muld(&y->xR[iy1], 1, iy2 - iy1 + 1, beta);
      }
   // alpha*A'*x
      for (i = i1; i <= i2; i++) {
         v = alpha * x->xR[ix1 + i - i1];
         ae_v_addd(&y->xR[iy1], 1, &a->xyR[i][j1], 1, iy2 - iy1 + 1, v);
      }
   }
}

void matrixmatrixmultiply(RMatrix *a, ae_int_t ai1, ae_int_t ai2, ae_int_t aj1, ae_int_t aj2, bool transa, RMatrix *b, ae_int_t bi1, ae_int_t bi2, ae_int_t bj1, ae_int_t bj2, bool transb, double alpha, RMatrix *c, ae_int_t ci1, ae_int_t ci2, ae_int_t cj1, ae_int_t cj2, double beta, RVector *work) {
   ae_int_t arows;
   ae_int_t acols;
   ae_int_t brows;
   ae_int_t bcols;
   ae_int_t crows;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t l;
   ae_int_t r;
   double v;
// Setup
   if (!transa) {
      arows = ai2 - ai1 + 1;
      acols = aj2 - aj1 + 1;
   } else {
      arows = aj2 - aj1 + 1;
      acols = ai2 - ai1 + 1;
   }
   if (!transb) {
      brows = bi2 - bi1 + 1;
      bcols = bj2 - bj1 + 1;
   } else {
      brows = bj2 - bj1 + 1;
      bcols = bi2 - bi1 + 1;
   }
   ae_assert(acols == brows, "matrixmatrixmultiply: incorrect matrix sizes!");
   if (arows <= 0 || acols <= 0 || brows <= 0 || bcols <= 0) {
      return;
   }
   crows = arows;
// Test WORK
   i = imax2(arows, acols);
   i = imax2(brows, i);
   i = imax2(i, bcols);
   work->xR[1] = 0.0;
   work->xR[i] = 0.0;
// Prepare C
   if (beta == 0.0) {
      for (i = ci1; i <= ci2; i++) {
         for (j = cj1; j <= cj2; j++) {
            c->xyR[i][j] = 0.0;
         }
      }
   } else {
      for (i = ci1; i <= ci2; i++) {
         ae_v_muld(&c->xyR[i][cj1], 1, cj2 - cj1 + 1, beta);
      }
   }
// A*B
   if (!transa && !transb) {
      for (l = ai1; l <= ai2; l++) {
         for (r = bi1; r <= bi2; r++) {
            v = alpha * a->xyR[l][aj1 + r - bi1];
            k = ci1 + l - ai1;
            ae_v_addd(&c->xyR[k][cj1], 1, &b->xyR[r][bj1], 1, cj2 - cj1 + 1, v);
         }
      }
      return;
   }
// A*B'
   if (!transa && transb) {
      if (arows * acols < brows * bcols) {
         for (r = bi1; r <= bi2; r++) {
            for (l = ai1; l <= ai2; l++) {
               v = ae_v_dotproduct(&a->xyR[l][aj1], 1, &b->xyR[r][bj1], 1, aj2 - aj1 + 1);
               c->xyR[ci1 + l - ai1][cj1 + r - bi1] += alpha * v;
            }
         }
         return;
      } else {
         for (l = ai1; l <= ai2; l++) {
            for (r = bi1; r <= bi2; r++) {
               v = ae_v_dotproduct(&a->xyR[l][aj1], 1, &b->xyR[r][bj1], 1, aj2 - aj1 + 1);
               c->xyR[ci1 + l - ai1][cj1 + r - bi1] += alpha * v;
            }
         }
         return;
      }
   }
// A'*B
   if (transa && !transb) {
      for (l = aj1; l <= aj2; l++) {
         for (r = bi1; r <= bi2; r++) {
            v = alpha * a->xyR[ai1 + r - bi1][l];
            k = ci1 + l - aj1;
            ae_v_addd(&c->xyR[k][cj1], 1, &b->xyR[r][bj1], 1, cj2 - cj1 + 1, v);
         }
      }
      return;
   }
// A'*B'
   if (transa && transb) {
      if (arows * acols < brows * bcols) {
         for (r = bi1; r <= bi2; r++) {
            k = cj1 + r - bi1;
            for (i = 1; i <= crows; i++) {
               work->xR[i] = 0.0;
            }
            for (l = ai1; l <= ai2; l++) {
               v = alpha * b->xyR[r][bj1 + l - ai1];
               ae_v_addd(&work->xR[1], 1, &a->xyR[l][aj1], 1, crows, v);
            }
            ae_v_add(&c->xyR[ci1][k], c->stride, &work->xR[1], 1, ci2 - ci1 + 1);
         }
         return;
      } else {
         for (l = aj1; l <= aj2; l++) {
            k = ai2 - ai1 + 1;
            ae_v_move(&work->xR[1], 1, &a->xyR[ai1][l], a->stride, k);
            for (r = bi1; r <= bi2; r++) {
               v = ae_v_dotproduct(&work->xR[1], 1, &b->xyR[r][bj1], 1, k);
               c->xyR[ci1 + l - aj1][cj1 + r - bi1] += alpha * v;
            }
         }
         return;
      }
   }
}
} // end of namespace alglib_impl

// === ROTATIONS Package ===
namespace alglib_impl {
// Application of a sequence of  elementary rotations to a matrix
//
// The algorithm pre-multiplies the matrix by a sequence of rotation
// transformations which is given by arrays C and S. Depending on the value
// of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward == true)
// rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.
//
// Not the whole matrix but only a part of it is transformed (rows from M1 to
// M2, columns from N1 to N2). Only the elements of this submatrix are changed.
//
// Inputs:
//     IsForward   -   the sequence of the rotation application.
//     M1,M2       -   the range of rows to be transformed.
//     N1, N2      -   the range of columns to be transformed.
//     C,S         -   transformation coefficients.
//                     Array whose index ranges within [1..M2-M1].
//     A           -   processed matrix.
//     WORK        -   working array whose index ranges within [N1..N2].
//
// Outputs:
//     A           -   transformed matrix.
//
// Utility subroutine.
void applyrotationsfromtheleft(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work) {
   ae_int_t j;
   ae_int_t jp1;
   double ctemp;
   double stemp;
   double temp;
   if (m1 > m2 || n1 > n2) {
      return;
   }
// Form  P * A
   if (isforward) {
      if (n1 != n2) {
      // Common case: N1 != N2
         for (j = m1; j < m2; j++) {
            ctemp = c->xR[j - m1 + 1];
            stemp = s->xR[j - m1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               jp1 = j + 1;
               ae_v_moved(&work->xR[n1], 1, &a->xyR[jp1][n1], 1, n2 - n1 + 1, ctemp);
               ae_v_subd(&work->xR[n1], 1, &a->xyR[j][n1], 1, n2 - n1 + 1, stemp);
               ae_v_muld(&a->xyR[j][n1], 1, n2 - n1 + 1, ctemp);
               ae_v_addd(&a->xyR[j][n1], 1, &a->xyR[jp1][n1], 1, n2 - n1 + 1, stemp);
               ae_v_move(&a->xyR[jp1][n1], 1, &work->xR[n1], 1, n2 - n1 + 1);
            }
         }
      } else {
      // Special case: N1 == N2
         for (j = m1; j < m2; j++) {
            ctemp = c->xR[j - m1 + 1];
            stemp = s->xR[j - m1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               temp = a->xyR[j + 1][n1];
               a->xyR[j + 1][n1] = ctemp * temp - stemp * a->xyR[j][n1];
               a->xyR[j][n1] = stemp * temp + ctemp * a->xyR[j][n1];
            }
         }
      }
   } else {
      if (n1 != n2) {
      // Common case: N1 != N2
         for (j = m2 - 1; j >= m1; j--) {
            ctemp = c->xR[j - m1 + 1];
            stemp = s->xR[j - m1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               jp1 = j + 1;
               ae_v_moved(&work->xR[n1], 1, &a->xyR[jp1][n1], 1, n2 - n1 + 1, ctemp);
               ae_v_subd(&work->xR[n1], 1, &a->xyR[j][n1], 1, n2 - n1 + 1, stemp);
               ae_v_muld(&a->xyR[j][n1], 1, n2 - n1 + 1, ctemp);
               ae_v_addd(&a->xyR[j][n1], 1, &a->xyR[jp1][n1], 1, n2 - n1 + 1, stemp);
               ae_v_move(&a->xyR[jp1][n1], 1, &work->xR[n1], 1, n2 - n1 + 1);
            }
         }
      } else {
      // Special case: N1 == N2
         for (j = m2 - 1; j >= m1; j--) {
            ctemp = c->xR[j - m1 + 1];
            stemp = s->xR[j - m1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               temp = a->xyR[j + 1][n1];
               a->xyR[j + 1][n1] = ctemp * temp - stemp * a->xyR[j][n1];
               a->xyR[j][n1] = stemp * temp + ctemp * a->xyR[j][n1];
            }
         }
      }
   }
}

// Application of a sequence of  elementary rotations to a matrix
//
// The algorithm post-multiplies the matrix by a sequence of rotation
// transformations which is given by arrays C and S. Depending on the value
// of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward == true)
// rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.
//
// Not the whole matrix but only a part of it is transformed (rows from M1
// to M2, columns from N1 to N2). Only the elements of this submatrix are changed.
//
// Inputs:
//     IsForward   -   the sequence of the rotation application.
//     M1,M2       -   the range of rows to be transformed.
//     N1, N2      -   the range of columns to be transformed.
//     C,S         -   transformation coefficients.
//                     Array whose index ranges within [1..N2-N1].
//     A           -   processed matrix.
//     WORK        -   working array whose index ranges within [M1..M2].
//
// Outputs:
//     A           -   transformed matrix.
//
// Utility subroutine.
void applyrotationsfromtheright(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work) {
   ae_int_t j;
   ae_int_t jp1;
   double ctemp;
   double stemp;
   double temp;
// Form A * P'
   if (isforward) {
      if (m1 != m2) {
      // Common case: M1 != M2
         for (j = n1; j < n2; j++) {
            ctemp = c->xR[j - n1 + 1];
            stemp = s->xR[j - n1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               jp1 = j + 1;
               ae_v_moved(&work->xR[m1], 1, &a->xyR[m1][jp1], a->stride, m2 - m1 + 1, ctemp);
               ae_v_subd(&work->xR[m1], 1, &a->xyR[m1][j], a->stride, m2 - m1 + 1, stemp);
               ae_v_muld(&a->xyR[m1][j], a->stride, m2 - m1 + 1, ctemp);
               ae_v_addd(&a->xyR[m1][j], a->stride, &a->xyR[m1][jp1], a->stride, m2 - m1 + 1, stemp);
               ae_v_move(&a->xyR[m1][jp1], a->stride, &work->xR[m1], 1, m2 - m1 + 1);
            }
         }
      } else {
      // Special case: M1 == M2
         for (j = n1; j < n2; j++) {
            ctemp = c->xR[j - n1 + 1];
            stemp = s->xR[j - n1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               temp = a->xyR[m1][j + 1];
               a->xyR[m1][j + 1] = ctemp * temp - stemp * a->xyR[m1][j];
               a->xyR[m1][j] = stemp * temp + ctemp * a->xyR[m1][j];
            }
         }
      }
   } else {
      if (m1 != m2) {
      // Common case: M1 != M2
         for (j = n2 - 1; j >= n1; j--) {
            ctemp = c->xR[j - n1 + 1];
            stemp = s->xR[j - n1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               jp1 = j + 1;
               ae_v_moved(&work->xR[m1], 1, &a->xyR[m1][jp1], a->stride, m2 - m1 + 1, ctemp);
               ae_v_subd(&work->xR[m1], 1, &a->xyR[m1][j], a->stride, m2 - m1 + 1, stemp);
               ae_v_muld(&a->xyR[m1][j], a->stride, m2 - m1 + 1, ctemp);
               ae_v_addd(&a->xyR[m1][j], a->stride, &a->xyR[m1][jp1], a->stride, m2 - m1 + 1, stemp);
               ae_v_move(&a->xyR[m1][jp1], a->stride, &work->xR[m1], 1, m2 - m1 + 1);
            }
         }
      } else {
      // Special case: M1 == M2
         for (j = n2 - 1; j >= n1; j--) {
            ctemp = c->xR[j - n1 + 1];
            stemp = s->xR[j - n1 + 1];
            if (ctemp != 1.0 || stemp != 0.0) {
               temp = a->xyR[m1][j + 1];
               a->xyR[m1][j + 1] = ctemp * temp - stemp * a->xyR[m1][j];
               a->xyR[m1][j] = stemp * temp + ctemp * a->xyR[m1][j];
            }
         }
      }
   }
}

// The subroutine generates the elementary rotation, so that:
//
// [  CS  SN  ]  .  [ F ]  =  [ R ]
// [ -SN  CS  ]     [ G ]     [ 0 ]
//
// CS**2 + SN**2 = 1
void generaterotation(double f, double g, double *cs, double *sn, double *r) {
   double f1;
   double g1;
   *cs = 0.0;
   *sn = 0.0;
   *r = 0.0;
   if (g == 0.0) {
      *cs = 1.0;
      *sn = 0.0;
      *r = f;
   } else {
      if (f == 0.0) {
         *cs = 0.0;
         *sn = 1.0;
         *r = g;
      } else {
         f1 = f;
         g1 = g;
         if (fabs(f1) > fabs(g1)) {
            *r = fabs(f1) * sqrt(1.0 + sqr(g1 / f1));
         } else {
            *r = fabs(g1) * sqrt(1.0 + sqr(f1 / g1));
         }
         *cs = f1 / *r;
         *sn = g1 / *r;
         if (fabs(f) > fabs(g) && *cs < 0.0) {
            *cs = -*cs;
            *sn = -*sn;
            *r = -*r;
         }
      }
   }
}
} // end of namespace alglib_impl

// === BASICSTATOPS Package ===
// Depends on: TSORT
namespace alglib_impl {
// Internal tied ranking subroutine.
//
// Inputs:
//     X       -   array to rank
//     N       -   array size
//     IsCentered- whether ranks are centered or not:
//                 * True      -   ranks are centered in such way that  their
//                                 sum is zero
//                 * False     -   ranks are not centered
//     Buf     -   temporary buffers
//
// NOTE: when IsCentered is True and all X[] are equal, this  function  fills
//       X by zeros (exact zeros are used, not sum which is only approximately
//       equal to zero).
void rankx(RVector *x, ae_int_t n, bool iscentered, apbuffers *buf) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double tmp;
   double voffs;
// Prepare
   if (n < 1) {
      return;
   }
   if (n == 1) {
      x->xR[0] = 0.0;
      return;
   }
   if (buf->ra1.cnt < n) {
      ae_vector_set_length(&buf->ra1, n);
   }
   if (buf->ia1.cnt < n) {
      ae_vector_set_length(&buf->ia1, n);
   }
   for (i = 0; i < n; i++) {
      buf->ra1.xR[i] = x->xR[i];
      buf->ia1.xZ[i] = i;
   }
   tagsortfasti(&buf->ra1, &buf->ia1, &buf->ra2, &buf->ia2, n);
// Special test for all values being equal
   if (buf->ra1.xR[0] == buf->ra1.xR[n - 1]) {
      if (iscentered) {
         tmp = 0.0;
      } else {
         tmp = (n - 1) / 2.0;
      }
      for (i = 0; i < n; i++) {
         x->xR[i] = tmp;
      }
      return;
   }
// compute tied ranks
   i = 0;
   while (i < n) {
      j = i + 1;
      while (j < n) {
         if (buf->ra1.xR[j] != buf->ra1.xR[i]) {
            break;
         }
         j++;
      }
      for (k = i; k < j; k++) {
         buf->ra1.xR[k] = (i + j - 1) / 2.0;
      }
      i = j;
   }
// back to x
   if (iscentered) {
      voffs = (n - 1) / 2.0;
   } else {
      voffs = 0.0;
   }
   for (i = 0; i < n; i++) {
      x->xR[buf->ia1.xZ[i]] = buf->ra1.xR[i] - voffs;
   }
}

// Internal untied ranking subroutine.
//
// Inputs:
//     X       -   array to rank
//     N       -   array size
//     Buf     -   temporary buffers
//
// Returns untied ranks (in case of a tie ranks are resolved arbitrarily).
void rankxuntied(RVector *x, ae_int_t n, apbuffers *buf) {
   ae_int_t i;
// Prepare
   if (n < 1) {
      return;
   }
   if (n == 1) {
      x->xR[0] = 0.0;
      return;
   }
   if (buf->ra1.cnt < n) {
      ae_vector_set_length(&buf->ra1, n);
   }
   if (buf->ia1.cnt < n) {
      ae_vector_set_length(&buf->ia1, n);
   }
   for (i = 0; i < n; i++) {
      buf->ra1.xR[i] = x->xR[i];
      buf->ia1.xZ[i] = i;
   }
   tagsortfasti(&buf->ra1, &buf->ia1, &buf->ra2, &buf->ia2, n);
   for (i = 0; i < n; i++) {
      x->xR[buf->ia1.xZ[i]] = i;
   }
}
} // end of namespace alglib_impl

// === APSTRUCT Package ===
// Depends on: ABLASF
namespace alglib_impl {
// Initializes n-set by empty structure.
//
// IMPORTANT: this function need O(N) time for initialization. It is recommended
//            to reduce its usage as much as possible, and use nisClear()
//            where possible.
//
// Inputs:
//     N           -   possible set size
//
// Outputs:
//     SA          -   empty N-set
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
void nisinitemptyslow(ae_int_t n, niset *sa) {
   sa->n = n;
   sa->nstored = 0;
   isetallocv(n, -999999999, &sa->locationof);
   isetallocv(n, -999999999, &sa->items);
}

// Clears set
//
// Inputs:
//     SA          -   set to be cleared
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
void nisclear(niset *sa) {
   ae_int_t i;
   ae_int_t ns;
   ns = sa->nstored;
   for (i = 0; i < ns; i++) {
      sa->locationof.xZ[sa->items.xZ[i]] = -1;
   }
   sa->nstored = 0;
}

// Copies n-set to properly initialized target set. The target set has to  be
// properly initialized, and it can be non-empty. If  it  is  non-empty,  its
// contents is quickly erased before copying.
//
// The cost of this function is O(max(SrcSize,DstSize))
//
// Inputs:
//     SSrc        -   source N-set
//     SDst        -   destination N-set (has same size as SSrc)
//
// Outputs:
//     SDst        -   copy of SSrc
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
void niscopy(niset *ssrc, niset *sdst) {
   ae_int_t ns;
   ae_int_t i;
   ae_int_t k;
   nisclear(sdst);
   ns = ssrc->nstored;
   for (i = 0; i < ns; i++) {
      k = ssrc->items.xZ[i];
      sdst->items.xZ[i] = k;
      sdst->locationof.xZ[k] = i;
   }
   sdst->nstored = ns;
}

// Add K-th element to the set. The element may already exist in the set.
//
// Inputs:
//     SA          -   set
//     K           -   element to add, 0 <= K < N.
//
// Outputs:
//     SA          -   modified SA
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
void nisaddelement(niset *sa, ae_int_t k) {
   ae_int_t ns;
   if (sa->locationof.xZ[k] >= 0) {
      return;
   }
   ns = sa->nstored;
   sa->locationof.xZ[k] = ns;
   sa->items.xZ[ns] = k;
   sa->nstored = ns + 1;
}

// Subtracts K-th set from the source structure
//
// Inputs:
//     SA          -   set
//     Src, K      -   source kn-set and set index K
//
// Outputs:
//     SA          -   modified SA
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
void nissubtract1(niset *sa, niset *src) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t loc;
   ae_int_t item;
   ae_int_t ns;
   ae_int_t ss;
   ns = sa->nstored;
   ss = src->nstored;
   if (ss < ns) {
      for (i = 0; i < ss; i++) {
         j = src->items.xZ[i];
         loc = sa->locationof.xZ[j];
         if (loc >= 0) {
            item = sa->items.xZ[ns - 1];
            sa->items.xZ[loc] = item;
            sa->locationof.xZ[item] = loc;
            sa->locationof.xZ[j] = -1;
            ns--;
         }
      }
   } else {
      i = 0;
      while (i < ns) {
         j = sa->items.xZ[i];
         loc = src->locationof.xZ[j];
         if (loc >= 0) {
            item = sa->items.xZ[ns - 1];
            sa->items.xZ[i] = item;
            sa->locationof.xZ[item] = i;
            sa->locationof.xZ[j] = -1;
            ns--;
         } else {
            i++;
         }
      }
   }
   sa->nstored = ns;
}

// Counts set elements
//
// Inputs:
//     SA          -   set
//
// Result:
//     number of elements in SA
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
ae_int_t niscount(niset *sa) {
   ae_int_t result;
   result = sa->nstored;
   return result;
}

// Compare two sets, returns True for equal sets
//
// Inputs:
//     S0          -   set 0
//     S1          -   set 1, must have same parameter N as set 0
//
// Result:
//     True, if sets are equal
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
bool nisequal(niset *s0, niset *s1) {
   ae_int_t i;
   ae_int_t ns0;
   ae_int_t ns1;
   bool result;
   result = false;
   if (s0->n != s1->n) {
      return result;
   }
   if (s0->nstored != s1->nstored) {
      return result;
   }
   ns0 = s0->nstored;
   ns1 = s1->nstored;
   for (i = 0; i < ns0; i++) {
      if (s1->locationof.xZ[s0->items.xZ[i]] < 0) {
         return result;
      }
   }
   for (i = 0; i < ns1; i++) {
      if (s0->locationof.xZ[s1->items.xZ[i]] < 0) {
         return result;
      }
   }
   result = true;
   return result;
}

// Prepares iteration over set
//
// Inputs:
//     SA          -   set
//
// Outputs:
//     SA          -   SA ready for repeated calls of nisEnumerate()
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
void nisstartenumeration(niset *sa) {
   sa->iteridx = 0;
}

// Iterates over the set. Subsequent calls return True and set J to  new  set
// item until iteration stops and False is returned.
//
// Inputs:
//     SA          -   n-set
//
// Outputs:
//     J           -   if:
//                     * Result=True - index of element in the set
//                     * Result=False - not set
// ALGLIB Project: Copyright 05.10.2020 by Sergey Bochkanov.
bool nisenumerate(niset *sa, ae_int_t *i) {
   ae_int_t k;
   bool result;
   *i = 0;
   k = sa->iteridx;
   if (k >= sa->nstored) {
      result = false;
      return result;
   }
   *i = sa->items.xZ[k];
   sa->iteridx = k + 1;
   result = true;
   return result;
}

void niset_init(void *_p, bool make_automatic) {
   niset *p = (niset *)_p;
   ae_vector_init(&p->items, 0, DT_INT, make_automatic);
   ae_vector_init(&p->locationof, 0, DT_INT, make_automatic);
}

void niset_copy(void *_dst, const void *_src, bool make_automatic) {
   niset *dst = (niset *)_dst;
   const niset *src = (const niset *)_src;
   dst->n = src->n;
   dst->nstored = src->nstored;
   ae_vector_copy(&dst->items, &src->items, make_automatic);
   ae_vector_copy(&dst->locationof, &src->locationof, make_automatic);
   dst->iteridx = src->iteridx;
}

void niset_free(void *_p, bool make_automatic) {
   niset *p = (niset *)_p;
   ae_vector_free(&p->items, make_automatic);
   ae_vector_free(&p->locationof, make_automatic);
}
} // end of namespace alglib_impl

// === TRLINSOLVE Package ===
namespace alglib_impl {
// Obsolete 1-based subroutine.
// See RMatrixTRSafeSolve for 0-based replacement.
void safesolvetriangular(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit, bool normin, RVector *cnorm) {
   ae_int_t i;
   ae_int_t imax;
   ae_int_t j;
   ae_int_t jfirst;
   ae_int_t jinc;
   ae_int_t jlast;
   ae_int_t jm1;
   ae_int_t jp1;
   ae_int_t ip1;
   ae_int_t im1;
   ae_int_t k;
   ae_int_t flg;
   double v;
   double vd;
   double bignum;
   double grow;
   double rec;
   double smlnum;
   double sumj;
   double tjj;
   double tjjs;
   double tmax;
   double tscal;
   double uscal;
   double xbnd;
   double xj;
   double xmax;
   bool notran;
   bool upper;
   bool nounit;
   *s = 0.0;
   upper = isupper;
   notran = !istrans;
   nounit = !isunit;
// these initializers are not really necessary,
// but without them compiler complains about uninitialized locals
   tjjs = 0.0;
// Quick return if possible
   if (n == 0) {
      return;
   }
// Determine machine dependent parameters to control overflow.
   smlnum = minrealnumber / (machineepsilon * 2.0);
   bignum = 1.0 / smlnum;
   *s = 1.0;
   if (!normin) {
      ae_vector_set_length(cnorm, n + 1);
   // Compute the 1-norm of each column, not including the diagonal.
      if (upper) {
      // A is upper triangular.
         for (j = 1; j <= n; j++) {
            v = 0.0;
            for (k = 1; k < j; k++) {
               v += fabs(a->xyR[k][j]);
            }
            cnorm->xR[j] = v;
         }
      } else {
      // A is lower triangular.
         for (j = 1; j < n; j++) {
            v = 0.0;
            for (k = j + 1; k <= n; k++) {
               v += fabs(a->xyR[k][j]);
            }
            cnorm->xR[j] = v;
         }
         cnorm->xR[n] = 0.0;
      }
   }
// Scale the column norms by TSCAL if the maximum element in CNORM is
// greater than BIGNUM.
   imax = 1;
   for (k = 2; k <= n; k++) {
      if (cnorm->xR[k] > cnorm->xR[imax]) {
         imax = k;
      }
   }
   tmax = cnorm->xR[imax];
   if (tmax <= bignum) {
      tscal = 1.0;
   } else {
      tscal = 1.0 / (smlnum * tmax);
      ae_v_muld(&cnorm->xR[1], 1, n, tscal);
   }
// Compute a bound on the computed solution vector to see if the
// Level 2 BLAS routine DTRSV can be used.
   j = 1;
   for (k = 2; k <= n; k++) {
      if (fabs(x->xR[k]) > fabs(x->xR[j])) {
         j = k;
      }
   }
   xmax = fabs(x->xR[j]);
   xbnd = xmax;
   if (notran) {
   // Compute the growth in A * x = b.
      if (upper) {
         jfirst = n;
         jlast = 1;
         jinc = -1;
      } else {
         jfirst = 1;
         jlast = n;
         jinc = 1;
      }
      if (tscal != 1.0) {
         grow = 0.0;
      } else {
         if (nounit) {
         // A is non-unit triangular.
         //
         // Compute GROW = 1/G(j) and XBND = 1/M(j).
         // Initially, G(0) = max{x(i), i = 1,...,n}.
            grow = 1.0 / rmax2(xbnd, smlnum);
            xbnd = grow;
            j = jfirst;
            while (jinc > 0 && j <= jlast || jinc < 0 && j >= jlast) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // M(j) = G(j-1) / abs(A(j,j))
               tjj = fabs(a->xyR[j][j]);
               xbnd = rmin2(xbnd, rmin2(1.0, tjj) * grow);
               if (tjj + cnorm->xR[j] >= smlnum) {
               // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                  grow *= tjj / (tjj + cnorm->xR[j]);
               } else {
               // G(j) could overflow, set GROW to 0.
                  grow = 0.0;
               }
               if (j == jlast) {
                  grow = xbnd;
               }
               j += jinc;
            }
         } else {
         // A is unit triangular.
         //
         // Compute GROW = 1/G(j), where G(0) = max{x(i), i = 1,...,n}.
            grow = rmin2(1.0, 1.0 / rmax2(xbnd, smlnum));
            j = jfirst;
            while (jinc > 0 && j <= jlast || jinc < 0 && j >= jlast) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // G(j) = G(j-1)*( 1 + CNORM(j) )
               grow /= 1.0 + cnorm->xR[j];
               j += jinc;
            }
         }
      }
   } else {
   // Compute the growth in A' * x = b.
      if (upper) {
         jfirst = 1;
         jlast = n;
         jinc = 1;
      } else {
         jfirst = n;
         jlast = 1;
         jinc = -1;
      }
      if (tscal != 1.0) {
         grow = 0.0;
      } else {
         if (nounit) {
         // A is non-unit triangular.
         //
         // Compute GROW = 1/G(j) and XBND = 1/M(j).
         // Initially, M(0) = max{x(i), i = 1,...,n}.
            grow = 1.0 / rmax2(xbnd, smlnum);
            xbnd = grow;
            j = jfirst;
            while (jinc > 0 && j <= jlast || jinc < 0 && j >= jlast) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
               xj = 1.0 + cnorm->xR[j];
               grow = rmin2(grow, xbnd / xj);
            // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
               tjj = fabs(a->xyR[j][j]);
               if (xj > tjj) {
                  xbnd *= tjj / xj;
               }
               if (j == jlast) {
                  grow = rmin2(grow, xbnd);
               }
               j += jinc;
            }
         } else {
         // A is unit triangular.
         //
         // Compute GROW = 1/G(j), where G(0) = max{x(i), i = 1,...,n}.
            grow = rmin2(1.0, 1.0 / rmax2(xbnd, smlnum));
            j = jfirst;
            while (jinc > 0 && j <= jlast || jinc < 0 && j >= jlast) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // G(j) = ( 1 + CNORM(j) )*G(j-1)
               xj = 1.0 + cnorm->xR[j];
               grow /= xj;
               j += jinc;
            }
         }
      }
   }
   if (grow * tscal > smlnum) {
   // Use the Level 2 BLAS solve if the reciprocal of the bound on
   // elements of X is not too small.
      if (upper == notran) {
         if (nounit) {
            vd = a->xyR[n][n];
         } else {
            vd = 1.0;
         }
         x->xR[n] /= vd;
         for (i = n - 1; i >= 1; i--) {
            ip1 = i + 1;
            if (upper) {
               v = ae_v_dotproduct(&a->xyR[i][ip1], 1, &x->xR[ip1], 1, n - ip1 + 1);
            } else {
               v = ae_v_dotproduct(&a->xyR[ip1][i], a->stride, &x->xR[ip1], 1, n - ip1 + 1);
            }
            if (nounit) {
               vd = a->xyR[i][i];
            } else {
               vd = 1.0;
            }
            x->xR[i] = (x->xR[i] - v) / vd;
         }
      } else {
         if (nounit) {
            vd = a->xyR[1][1];
         } else {
            vd = 1.0;
         }
         x->xR[1] /= vd;
         for (i = 2; i <= n; i++) {
            im1 = i - 1;
            if (upper) {
               v = ae_v_dotproduct(&a->xyR[1][i], a->stride, &x->xR[1], 1, im1);
            } else {
               v = ae_v_dotproduct(&a->xyR[i][1], 1, &x->xR[1], 1, im1);
            }
            if (nounit) {
               vd = a->xyR[i][i];
            } else {
               vd = 1.0;
            }
            x->xR[i] = (x->xR[i] - v) / vd;
         }
      }
   } else {
   // Use a Level 1 BLAS solve, scaling intermediate results.
      if (xmax > bignum) {
      // Scale X so that its components are less than or equal to
      // BIGNUM in absolute value.
         *s = bignum / xmax;
         ae_v_muld(&x->xR[1], 1, n, *s);
         xmax = bignum;
      }
      if (notran) {
      // Solve A * x = b
         j = jfirst;
         while (jinc > 0 && j <= jlast || jinc < 0 && j >= jlast) {
         // Compute x(j) = b(j) / A(j,j), scaling x if necessary.
            xj = fabs(x->xR[j]);
            flg = 0;
            if (nounit) {
               tjjs = a->xyR[j][j] * tscal;
            } else {
               tjjs = tscal;
               if (tscal == 1.0) {
                  flg = 100;
               }
            }
            if (flg != 100) {
               tjj = fabs(tjjs);
               if (tjj > smlnum) {
               // abs(A(j,j)) > SMLNUM:
                  if (tjj < 1.0) {
                     if (xj > tjj * bignum) {
                     // Scale x by 1/b(j).
                        rec = 1.0 / xj;
                        ae_v_muld(&x->xR[1], 1, n, rec);
                        *s *= rec;
                        xmax *= rec;
                     }
                  }
                  x->xR[j] /= tjjs;
                  xj = fabs(x->xR[j]);
               } else {
                  if (tjj > 0.0) {
                  // 0 < abs(A(j,j)) <= SMLNUM:
                     if (xj > tjj * bignum) {
                     // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                     // to avoid overflow when dividing by A(j,j).
                        rec = tjj * bignum / xj;
                        if (cnorm->xR[j] > 1.0) {
                        // Scale by 1/CNORM(j) to avoid overflow when
                        // multiplying x(j) times column j.
                           rec /= cnorm->xR[j];
                        }
                        ae_v_muld(&x->xR[1], 1, n, rec);
                        *s *= rec;
                        xmax *= rec;
                     }
                     x->xR[j] /= tjjs;
                     xj = fabs(x->xR[j]);
                  } else {
                  // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                  // scale = 0, and compute a solution to A*x = 0.
                     for (i = 1; i <= n; i++) {
                        x->xR[i] = 0.0;
                     }
                     x->xR[j] = 1.0;
                     xj = 1.0;
                     *s = 0.0;
                     xmax = 0.0;
                  }
               }
            }
         // Scale x if necessary to avoid overflow when adding a
         // multiple of column j of A.
            if (xj > 1.0) {
               rec = 1.0 / xj;
               if (cnorm->xR[j] > (bignum - xmax) * rec) {
               // Scale x by 1/(2*abs(x(j))).
                  rec *= 0.5;
                  ae_v_muld(&x->xR[1], 1, n, rec);
                  *s *= rec;
               }
            } else {
               if (xj * cnorm->xR[j] > bignum - xmax) {
               // Scale x by 1/2.
                  ae_v_muld(&x->xR[1], 1, n, 0.5);
                  *s *= 0.5;
               }
            }
            if (upper) {
               if (j > 1) {
               // Compute the update
               // x(1:j-1) = x(1:j-1) - x(j) * A(1:j-1,j)
                  v = x->xR[j] * tscal;
                  jm1 = j - 1;
                  ae_v_subd(&x->xR[1], 1, &a->xyR[1][j], a->stride, jm1, v);
                  i = 1;
                  for (k = 2; k < j; k++) {
                     if (fabs(x->xR[k]) > fabs(x->xR[i])) {
                        i = k;
                     }
                  }
                  xmax = fabs(x->xR[i]);
               }
            } else {
               if (j < n) {
               // Compute the update
               // x(j+1:n) = x(j+1:n) - x(j) * A(j+1:n,j)
                  jp1 = j + 1;
                  v = x->xR[j] * tscal;
                  ae_v_subd(&x->xR[jp1], 1, &a->xyR[jp1][j], a->stride, n - jp1 + 1, v);
                  i = j + 1;
                  for (k = j + 2; k <= n; k++) {
                     if (fabs(x->xR[k]) > fabs(x->xR[i])) {
                        i = k;
                     }
                  }
                  xmax = fabs(x->xR[i]);
               }
            }
            j += jinc;
         }
      } else {
      // Solve A' * x = b
         j = jfirst;
         while (jinc > 0 && j <= jlast || jinc < 0 && j >= jlast) {
         // Compute x(j) = b(j) - sum A(k,j)*x(k).
         //   k != j
            xj = fabs(x->xR[j]);
            uscal = tscal;
            rec = 1.0 / rmax2(xmax, 1.0);
            if (cnorm->xR[j] > (bignum - xj) * rec) {
            // If x(j) could overflow, scale x by 1/(2*XMAX).
               rec *= 0.5;
               if (nounit) {
                  tjjs = a->xyR[j][j] * tscal;
               } else {
                  tjjs = tscal;
               }
               tjj = fabs(tjjs);
               if (tjj > 1.0) {
               // Divide by A(j,j) when scaling x if A(j,j) > 1.
                  rec = rmin2(1.0, rec * tjj);
                  uscal /= tjjs;
               }
               if (rec < 1.0) {
                  ae_v_muld(&x->xR[1], 1, n, rec);
                  *s *= rec;
                  xmax *= rec;
               }
            }
            sumj = 0.0;
            if (uscal == 1.0) {
            // If the scaling needed for A in the dot product is 1,
            // call DDOT to perform the dot product.
               if (upper) {
                  if (j > 1) {
                     jm1 = j - 1;
                     sumj = ae_v_dotproduct(&a->xyR[1][j], a->stride, &x->xR[1], 1, jm1);
                  } else {
                     sumj = 0.0;
                  }
               } else {
                  if (j < n) {
                     jp1 = j + 1;
                     sumj = ae_v_dotproduct(&a->xyR[jp1][j], a->stride, &x->xR[jp1], 1, n - jp1 + 1);
                  }
               }
            } else {
            // Otherwise, use in-line code for the dot product.
               if (upper) {
                  for (i = 1; i < j; i++) {
                     v = a->xyR[i][j] * uscal;
                     sumj += v * x->xR[i];
                  }
               } else {
                  if (j < n) {
                     for (i = j + 1; i <= n; i++) {
                        v = a->xyR[i][j] * uscal;
                        sumj += v * x->xR[i];
                     }
                  }
               }
            }
            if (uscal == tscal) {
            // Compute x(j) = ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
            // was not used to scale the dotproduct.
               x->xR[j] -= sumj;
               xj = fabs(x->xR[j]);
               flg = 0;
               if (nounit) {
                  tjjs = a->xyR[j][j] * tscal;
               } else {
                  tjjs = tscal;
                  if (tscal == 1.0) {
                     flg = 150;
                  }
               }
            // Compute x(j) = x(j) / A(j,j), scaling if necessary.
               if (flg != 150) {
                  tjj = fabs(tjjs);
                  if (tjj > smlnum) {
                  // abs(A(j,j)) > SMLNUM:
                     if (tjj < 1.0) {
                        if (xj > tjj * bignum) {
                        // Scale X by 1/abs(x(j)).
                           rec = 1.0 / xj;
                           ae_v_muld(&x->xR[1], 1, n, rec);
                           *s *= rec;
                           xmax *= rec;
                        }
                     }
                     x->xR[j] /= tjjs;
                  } else {
                     if (tjj > 0.0) {
                     // 0 < abs(A(j,j)) <= SMLNUM:
                        if (xj > tjj * bignum) {
                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                           rec = tjj * bignum / xj;
                           ae_v_muld(&x->xR[1], 1, n, rec);
                           *s *= rec;
                           xmax *= rec;
                        }
                        x->xR[j] /= tjjs;
                     } else {
                     // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     // scale = 0, and compute a solution to A'*x = 0.
                        for (i = 1; i <= n; i++) {
                           x->xR[i] = 0.0;
                        }
                        x->xR[j] = 1.0;
                        *s = 0.0;
                        xmax = 0.0;
                     }
                  }
               }
            } else {
            // Compute x(j) = x(j) / A(j,j)  - sumj if the dot
            // product has already been divided by 1/A(j,j).
               x->xR[j] = x->xR[j] / tjjs - sumj;
            }
            xmax = rmax2(xmax, fabs(x->xR[j]));
            j += jinc;
         }
      }
      *s /= tscal;
   }
// Scale the column norms by 1/TSCAL for return.
   if (tscal != 1.0) {
      v = 1.0 / tscal;
      ae_v_muld(&cnorm->xR[1], 1, n, v);
   }
}

// Utility subroutine performing the "safe" solution of system of linear
// equations with triangular coefficient matrices.
//
// The subroutine uses scaling and solves the scaled system A*x == s*b (where  s
// is  a  scalar  value)  instead  of  A*x == b,  choosing  s  so  that x can be
// represented by a floating-point number. The closer the system  gets  to  a
// singular, the less s is. If the system is singular, s == 0 and x contains the
// non-trivial solution of equation A*x == 0.
//
// The feature of an algorithm is that it could not cause an  overflow  or  a
// division by zero regardless of the matrix used as the input.
//
// The algorithm can solve systems of equations with  upper/lower  triangular
// matrices,  with/without unit diagonal, and systems of type A*x == b or A'*x == b
// (where A' is a transposed matrix A).
//
// Inputs:
//     A       -   system matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     X       -   right-hand member of a system.
//                 Array whose index ranges within [0..N-1].
//     IsUpper -   matrix type. If it is True, the system matrix is the upper
//                 triangular and is located in  the  corresponding  part  of
//                 matrix A.
//     Trans   -   problem type. If it is True, the problem to be  solved  is
//                 A'*x == b, otherwise it is A*x == b.
//     Isunit  -   matrix type. If it is True, the system matrix has  a  unit
//                 diagonal (the elements on the main diagonal are  not  used
//                 in the calculation process), otherwise the matrix is considered
//                 to be a general triangular matrix.
//
// Outputs:
//     X       -   solution. Array whose index ranges within [0..N-1].
//     S       -   scaling factor.
//
// LAPACK auxiliary routine (version 3.0):
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      June 30, 1992
void rmatrixtrsafesolve(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit) {
   ae_frame _frame_block;
   bool normin;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   *s = 0.0;
   NewVector(cnorm, 0, DT_REAL);
   NewMatrix(a1, 0, 0, DT_REAL);
   NewVector(x1, 0, DT_REAL);
// From 0-based to 1-based
   normin = false;
   ae_matrix_set_length(&a1, n + 1, n + 1);
   ae_vector_set_length(&x1, n + 1);
   for (i = 1; i <= n; i++) {
      ae_v_move(&a1.xyR[i][1], 1, a->xyR[i - 1], 1, n);
   }
   ae_v_move(&x1.xR[1], 1, x->xR, 1, n);
// Solve 1-based
   safesolvetriangular(&a1, n, &x1, s, isupper, istrans, isunit, normin, &cnorm);
// From 1-based to 0-based
   ae_v_move(x->xR, 1, &x1.xR[1], 1, n);
   ae_frame_leave();
}
} // end of namespace alglib_impl

// === SAFESOLVE Package ===
namespace alglib_impl {
// complex basic solver-updater for reduced linear system
//
//     alpha*x[i] = beta
//
// solves this equation and updates it in overlfow-safe manner (keeping track
// of relative growth of solution).
//
// Parameters:
//     Alpha   -   alpha
//     Beta    -   beta
//     LnMax   -   precomputed Ln(MaxRealNumber)
//     BNorm   -   inf-norm of b (right part of original system)
//     MaxGrowth-  maximum growth of norm(x) relative to norm(b)
//     XNorm   -   inf-norm of other components of X (which are already processed)
//                 it is updated by CBasicSolveAndUpdate.
//     X       -   solution
// ALGLIB Routine: Copyright 26.01.2009 by Sergey Bochkanov
static bool safesolve_cbasicsolveandupdate(complex alpha, complex beta, double lnmax, double bnorm, double maxgrowth, double *xnorm, complex *x) {
   double v;
   bool result;
   x->x = 0.0;
   x->y = 0.0;
   result = false;
   if (ae_c_eq_d(alpha, 0.0)) {
      return result;
   }
   if (ae_c_neq_d(beta, 0.0)) {
   // alpha*x[i] == beta
      v = log(abscomplex(beta)) - log(abscomplex(alpha));
      if (v > lnmax) {
         return result;
      }
      *x = ae_c_div(beta, alpha);
   } else {
   // alpha*x[i] == 0
      *x = complex_from_i(0);
   }
// update NrmX, test growth limit
   *xnorm = rmax2(*xnorm, abscomplex(*x));
   if (*xnorm > maxgrowth * bnorm) {
      return result;
   }
   result = true;
   return result;
}

// Real implementation of CMatrixScaledTRSafeSolve
// ALGLIB Routine: Copyright 21.01.2010 by Sergey Bochkanov
bool rmatrixscaledtrsafesolve(RMatrix *a, double sa, ae_int_t n, RVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth) {
   ae_frame _frame_block;
   double lnmax;
   double nrmb;
   double nrmx;
   ae_int_t i;
   complex alpha;
   complex beta;
   double vr;
   complex cx;
   bool result;
   ae_frame_make(&_frame_block);
   NewVector(tmp, 0, DT_REAL);
   ae_assert(n > 0, "rmatrixscaledtrsafesolve: incorrect n!");
   ae_assert(trans == 0 || trans == 1, "rmatrixscaledtrsafesolve: incorrect trans!");
   result = true;
   lnmax = log(maxrealnumber);
// Quick return if possible
   if (n <= 0) {
      ae_frame_leave();
      return result;
   }
// Load norms: right part and X
   nrmb = 0.0;
   for (i = 0; i < n; i++) {
      nrmb = rmax2(nrmb, fabs(x->xR[i]));
   }
   nrmx = 0.0;
// Solve
   ae_vector_set_length(&tmp, n);
   result = true;
   if (isupper && trans == 0) {
   // U*x = b
      for (i = n - 1; i >= 0; i--) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = complex_from_d(a->xyR[i][i] * sa);
         }
         if (i < n - 1) {
            ae_v_moved(&tmp.xR[i + 1], 1, &a->xyR[i][i + 1], 1, n - i - 1, sa);
            vr = ae_v_dotproduct(&tmp.xR[i + 1], 1, &x->xR[i + 1], 1, n - i - 1);
            beta = complex_from_d(x->xR[i] - vr);
         } else {
            beta = complex_from_d(x->xR[i]);
         }
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xR[i] = cx.x;
      }
      ae_frame_leave();
      return result;
   }
   if (!isupper && trans == 0) {
   // L*x = b
      for (i = 0; i < n; i++) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = complex_from_d(a->xyR[i][i] * sa);
         }
         if (i > 0) {
            ae_v_moved(tmp.xR, 1, a->xyR[i], 1, i, sa);
            vr = ae_v_dotproduct(tmp.xR, 1, x->xR, 1, i);
            beta = complex_from_d(x->xR[i] - vr);
         } else {
            beta = complex_from_d(x->xR[i]);
         }
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xR[i] = cx.x;
      }
      ae_frame_leave();
      return result;
   }
   if (isupper && trans == 1) {
   // U^T*x = b
      for (i = 0; i < n; i++) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = complex_from_d(a->xyR[i][i] * sa);
         }
         beta = complex_from_d(x->xR[i]);
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xR[i] = cx.x;
      // update the rest of right part
         if (i < n - 1) {
            vr = cx.x;
            ae_v_moved(&tmp.xR[i + 1], 1, &a->xyR[i][i + 1], 1, n - i - 1, sa);
            ae_v_subd(&x->xR[i + 1], 1, &tmp.xR[i + 1], 1, n - i - 1, vr);
         }
      }
      ae_frame_leave();
      return result;
   }
   if (!isupper && trans == 1) {
   // L^T*x = b
      for (i = n - 1; i >= 0; i--) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = complex_from_d(a->xyR[i][i] * sa);
         }
         beta = complex_from_d(x->xR[i]);
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xR[i] = cx.x;
      // update the rest of right part
         if (i > 0) {
            vr = cx.x;
            ae_v_moved(tmp.xR, 1, a->xyR[i], 1, i, sa);
            ae_v_subd(x->xR, 1, tmp.xR, 1, i, vr);
         }
      }
      ae_frame_leave();
      return result;
   }
   result = false;
   ae_frame_leave();
   return result;
}

// Internal subroutine for safe solution of
//
//     SA*op(A) == b
//
// where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
// either identity transform, transposition or Hermitian transposition, SA is
// a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.
//
// This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
// MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
// near-degenerate matrices are handled correctly (False is returned) as long
// as MaxGrowth is significantly less than MaxRealNumber/norm(b).
// ALGLIB Routine: Copyright 21.01.2010 by Sergey Bochkanov
bool cmatrixscaledtrsafesolve(CMatrix *a, double sa, ae_int_t n, CVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth) {
   ae_frame _frame_block;
   double lnmax;
   double nrmb;
   double nrmx;
   ae_int_t i;
   complex alpha;
   complex beta;
   complex vc;
   bool result;
   ae_frame_make(&_frame_block);
   NewVector(tmp, 0, DT_COMPLEX);
   ae_assert(n > 0, "cmatrixscaledtrsafesolve: incorrect n!");
   ae_assert(trans == 0 || trans == 1 || trans == 2, "cmatrixscaledtrsafesolve: incorrect trans!");
   result = true;
   lnmax = log(maxrealnumber);
// Quick return if possible
   if (n <= 0) {
      ae_frame_leave();
      return result;
   }
// Load norms: right part and X
   nrmb = 0.0;
   for (i = 0; i < n; i++) {
      nrmb = rmax2(nrmb, abscomplex(x->xC[i]));
   }
   nrmx = 0.0;
// Solve
   ae_vector_set_length(&tmp, n);
   result = true;
   if (isupper && trans == 0) {
   // U*x = b
      for (i = n - 1; i >= 0; i--) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(a->xyC[i][i], sa);
         }
         if (i < n - 1) {
            ae_v_cmoved(&tmp.xC[i + 1], 1, &a->xyC[i][i + 1], 1, "N", n - i - 1, sa);
            vc = ae_v_cdotproduct(&tmp.xC[i + 1], 1, "N", &x->xC[i + 1], 1, "N", n - i - 1);
            beta = ae_c_sub(x->xC[i], vc);
         } else {
            beta = x->xC[i];
         }
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xC[i] = vc;
      }
      ae_frame_leave();
      return result;
   }
   if (!isupper && trans == 0) {
   // L*x = b
      for (i = 0; i < n; i++) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(a->xyC[i][i], sa);
         }
         if (i > 0) {
            ae_v_cmoved(tmp.xC, 1, a->xyC[i], 1, "N", i, sa);
            vc = ae_v_cdotproduct(tmp.xC, 1, "N", x->xC, 1, "N", i);
            beta = ae_c_sub(x->xC[i], vc);
         } else {
            beta = x->xC[i];
         }
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xC[i] = vc;
      }
      ae_frame_leave();
      return result;
   }
   if (isupper && trans == 1) {
   // U^T*x = b
      for (i = 0; i < n; i++) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(a->xyC[i][i], sa);
         }
         beta = x->xC[i];
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i < n - 1) {
            ae_v_cmoved(&tmp.xC[i + 1], 1, &a->xyC[i][i + 1], 1, "N", n - i - 1, sa);
            ae_v_csubc(&x->xC[i + 1], 1, &tmp.xC[i + 1], 1, "N", n - i - 1, vc);
         }
      }
      ae_frame_leave();
      return result;
   }
   if (!isupper && trans == 1) {
   // L^T*x = b
      for (i = n - 1; i >= 0; i--) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(a->xyC[i][i], sa);
         }
         beta = x->xC[i];
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i > 0) {
            ae_v_cmoved(tmp.xC, 1, a->xyC[i], 1, "N", i, sa);
            ae_v_csubc(x->xC, 1, tmp.xC, 1, "N", i, vc);
         }
      }
      ae_frame_leave();
      return result;
   }
   if (isupper && trans == 2) {
   // U^H*x = b
      for (i = 0; i < n; i++) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(conj(a->xyC[i][i]), sa);
         }
         beta = x->xC[i];
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i < n - 1) {
            ae_v_cmoved(&tmp.xC[i + 1], 1, &a->xyC[i][i + 1], 1, "Conj", n - i - 1, sa);
            ae_v_csubc(&x->xC[i + 1], 1, &tmp.xC[i + 1], 1, "N", n - i - 1, vc);
         }
      }
      ae_frame_leave();
      return result;
   }
   if (!isupper && trans == 2) {
   // L^T*x = b
      for (i = n - 1; i >= 0; i--) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(conj(a->xyC[i][i]), sa);
         }
         beta = x->xC[i];
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc);
         if (!result) {
            ae_frame_leave();
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i > 0) {
            ae_v_cmoved(tmp.xC, 1, a->xyC[i], 1, "Conj", i, sa);
            ae_v_csubc(x->xC, 1, tmp.xC, 1, "N", i, vc);
         }
      }
      ae_frame_leave();
      return result;
   }
   result = false;
   ae_frame_leave();
   return result;
}
} // end of namespace alglib_impl

// === XBLAS Package ===
namespace alglib_impl {
// Fast pow
// ALGLIB: Copyright 24.08.2009 by Sergey Bochkanov
static double xblas_xfastpow(double r, ae_int_t n) {
   double result;
   result = 0.0;
   if (n > 0) {
      if (n % 2 == 0) {
         result = sqr(xblas_xfastpow(r, n / 2));
      } else {
         result = r * xblas_xfastpow(r, n - 1);
      }
      return result;
   }
   if (n == 0) {
      result = 1.0;
   }
   if (n < 0) {
      result = xblas_xfastpow(1.0 / r, -n);
   }
   return result;
}

// Internal subroutine for extra-precise calculation of SUM(w[i]).
//
// Inputs:
//     W   -   array[0..N-1], values to be added
//             W is modified during calculations.
//     MX  -   max(W[i])
//     N   -   array size
//
// Outputs:
//     R   -   SUM(w[i])
//     RErr-   error estimate for R
// ALGLIB: Copyright 24.08.2009 by Sergey Bochkanov
static void xblas_xsum(RVector *w, double mx, ae_int_t n, double *r, double *rerr) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t ks;
   double v;
   double s;
   double ln2;
   double chunk;
   double invchunk;
   bool allzeros;
   *r = 0.0;
   *rerr = 0.0;
// special cases:
// * N == 0
// * N is too large to use integer arithmetics
   if (n == 0) {
      *r = 0.0;
      *rerr = 0.0;
      return;
   }
   if (mx == 0.0) {
      *r = 0.0;
      *rerr = 0.0;
      return;
   }
   ae_assert(n < 0x20000000, "xblas_xsum: n is too large!");
// Prepare
   ln2 = log(2.0);
   *rerr = mx * machineepsilon;
// 1. find S such that 0.5 <= S*MX < 1
// 2. multiply W by S, so task is normalized in some sense
// 3. S = 1/S so we can obtain original vector multiplying by S
   k = iround(log(mx) / ln2);
   s = xblas_xfastpow(2.0, -k);
   if (!isfinite(s)) {
   // Overflow or underflow during evaluation of S; fallback low-precision code
      *r = 0.0;
      *rerr = mx * machineepsilon;
      for (i = 0; i < n; i++) {
         *r += w->xR[i];
      }
      return;
   }
   while (s * mx >= 1.0) {
      s *= 0.5;
   }
   while (s * mx < 0.5) {
      s *= 2.0;
   }
   ae_v_muld(w->xR, 1, n, s);
   s = 1.0 / s;
// find Chunk == 2^M such that N*Chunk < 2^29
//
// we have chosen upper limit (2^29) with enough space left
// to tolerate possible problems with rounding and N's close
// to the limit, so we don't want to be very strict here.
   k = itrunc(log(536870912.0 / n) / ln2);
   chunk = xblas_xfastpow(2.0, k);
   if (chunk < 2.0) {
      chunk = 2.0;
   }
   invchunk = 1.0 / chunk;
// calculate result
   *r = 0.0;
   ae_v_muld(w->xR, 1, n, chunk);
   while (true) {
      s *= invchunk;
      allzeros = true;
      ks = 0;
      for (i = 0; i < n; i++) {
         v = w->xR[i];
         k = itrunc(v);
         if (v != k) {
            allzeros = false;
         }
         w->xR[i] = chunk * (v - k);
         ks += k;
      }
      *r += s * ks;
      if (allzeros || s * n + mx == mx) {
         break;
      }
   }
// correct error
   *rerr = rmax2(*rerr, fabs(*r) * machineepsilon);
}

// More precise dot-product. Absolute error of  subroutine  result  is  about
// 1 ulp of max(MX,V), where:
//     MX = max( |a[i]*b[i]| )
//     V  = |(a,b)|
//
// Inputs:
//     A       -   array[0..N-1], vector 1
//     B       -   array[0..N-1], vector 2
//     N       -   vectors length, N < 2^29.
//     Temp    -   array[0..N-1], preallocated temporary storage
//
// Outputs:
//     R       -   (A,B)
//     RErr    -   estimate of error. This estimate accounts for both  errors
//                 during  calculation  of  (A,B)  and  errors  introduced by
//                 rounding of A and B to fit in double (about 1 ulp).
// ALGLIB: Copyright 24.08.2009 by Sergey Bochkanov
void xdot(RVector *a, RVector *b, ae_int_t n, RVector *temp, double *r, double *rerr) {
   ae_int_t i;
   double mx;
   double v;
   *r = 0.0;
   *rerr = 0.0;
// special cases:
// * N == 0
   if (n == 0) {
      *r = 0.0;
      *rerr = 0.0;
      return;
   }
   mx = 0.0;
   for (i = 0; i < n; i++) {
      v = a->xR[i] * b->xR[i];
      temp->xR[i] = v;
      mx = rmax2(mx, fabs(v));
   }
   if (mx == 0.0) {
      *r = 0.0;
      *rerr = 0.0;
      return;
   }
   xblas_xsum(temp, mx, n, r, rerr);
}

// More precise complex dot-product. Absolute error of  subroutine  result is
// about 1 ulp of max(MX,V), where:
//     MX = max( |a[i]*b[i]| )
//     V  = |(a,b)|
//
// Inputs:
//     A       -   array[0..N-1], vector 1
//     B       -   array[0..N-1], vector 2
//     N       -   vectors length, N < 2^29.
//     Temp    -   array[0..2*N-1], preallocated temporary storage
//
// Outputs:
//     R       -   (A,B)
//     RErr    -   estimate of error. This estimate accounts for both  errors
//                 during  calculation  of  (A,B)  and  errors  introduced by
//                 rounding of A and B to fit in double (about 1 ulp).
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void xcdot(CVector *a, CVector *b, ae_int_t n, RVector *temp, complex *r, double *rerr) {
   ae_int_t i;
   double mx;
   double v;
   double rerrx;
   double rerry;
   r->x = 0.0;
   r->y = 0.0;
   *rerr = 0.0;
// special cases:
// * N == 0
   if (n == 0) {
      *r = complex_from_i(0);
      *rerr = 0.0;
      return;
   }
// calculate real part
   mx = 0.0;
   for (i = 0; i < n; i++) {
      v = a->xC[i].x * b->xC[i].x;
      temp->xR[2 * i] = v;
      mx = rmax2(mx, fabs(v));
      v = -a->xC[i].y * b->xC[i].y;
      temp->xR[2 * i + 1] = v;
      mx = rmax2(mx, fabs(v));
   }
   if (mx == 0.0) {
      r->x = 0.0;
      rerrx = 0.0;
   } else {
      xblas_xsum(temp, mx, 2 * n, &r->x, &rerrx);
   }
// calculate imaginary part
   mx = 0.0;
   for (i = 0; i < n; i++) {
      v = a->xC[i].x * b->xC[i].y;
      temp->xR[2 * i] = v;
      mx = rmax2(mx, fabs(v));
      v = a->xC[i].y * b->xC[i].x;
      temp->xR[2 * i + 1] = v;
      mx = rmax2(mx, fabs(v));
   }
   if (mx == 0.0) {
      r->y = 0.0;
      rerry = 0.0;
   } else {
      xblas_xsum(temp, mx, 2 * n, &r->y, &rerry);
   }
// total error
   if (rerrx == 0.0 && rerry == 0.0) {
      *rerr = 0.0;
   } else {
      *rerr = rmax2(rerrx, rerry) * sqrt(1.0 + sqr(rmin2(rerrx, rerry) / rmax2(rerrx, rerry)));
   }
}
} // end of namespace alglib_impl

// === LINMIN Package ===
namespace alglib_impl {
static const double linmin_stpmin = 1.0E-50;

// Normalizes direction/step pair: makes |D| == 1, scales Stp.
// If |D| == 0, it returns, leavind D/Stp unchanged.
// ALGLIB: Copyright 01.04.2010 by Sergey Bochkanov
void linminnormalized(RVector *d, double *stp, ae_int_t n) {
   double mx;
   double s;
   ae_int_t i;
// first, scale D to avoid underflow/overflow durng squaring
   mx = 0.0;
   for (i = 0; i < n; i++) {
      mx = rmax2(mx, fabs(d->xR[i]));
   }
   if (mx == 0.0) {
      return;
   }
   s = 1.0 / mx;
   ae_v_muld(d->xR, 1, n, s);
   *stp /= s;
// normalize D
   s = ae_v_dotproduct(d->xR, 1, d->xR, 1, n);
   s = 1.0 / sqrt(s);
   ae_v_muld(d->xR, 1, n, s);
   *stp /= s;
}

void mcstep(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, double fp, double dp, bool *brackt, double stmin, double stmax, ae_int_t *info) {
   bool bracktbound;
   double stpf;
   *info = 0;
// Check the inputs for errors.
   if (*brackt && (*stp <= rmin2(*stx, *sty) || *stp >= rmax2(*stx, *sty)) || *dx * (*stp - *stx) >= 0.0 || stmax < stmin) return;
// Determine if the derivatives have opposite sign.
   double sgnd = dp * (*dx / fabs(*dx));
// The first case: a higher function value.
// The minimum is bracketed.
// If the cubic step is closer to *stx than the quadratic step, the cubic step is taken,
// else the average of the cubic and quadratic steps is taken.
   if (fp > *fx) {
      *info = 1;
      double theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
      double s = rmax2(fabs(theta), rmax2(fabs(*dx), fabs(dp)));
      double gamma = s * sqrt(sqr(theta / s) - *dx / s * (dp / s));
      if (*stp < *stx) gamma = -gamma;
      double r = (gamma - *dx + theta) / (gamma - *dx + gamma + dp);
      double stpc = *stx + r * (*stp - *stx);
      double stpq = *stx + *dx / ((*fx - fp) / (*stp - *stx) + *dx) / 2 * (*stp - *stx);
      stpf = fabs(stpc - *stx) < fabs(stpq - *stx) ? stpc : stpc + (stpq - stpc) / 2;
      bracktbound = *brackt = true;
   } else if (sgnd < 0.0) {
   // The second case: a lower function value and derivatives of opposite sign.
   // The minimum is bracketed.
   // If the cubic step is closer to *stx than the quadratic (secant) step, the cubic step is taken, else the quadratic step is taken.
      *info = 2;
      bracktbound = false;
      double theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
      double s = rmax2(fabs(theta), rmax2(fabs(*dx), fabs(dp)));
      double gamma = s * sqrt(sqr(theta / s) - *dx / s * (dp / s));
      if (*stp > *stx) gamma = -gamma;
      double r = (gamma - dp + theta) / (gamma - dp + gamma + *dx);
      double stpc = *stp + r * (*stx - *stp);
      double stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
      stpf = fabs(stpc - *stp) > fabs(stpq - *stp) ? stpc : stpq;
      *brackt = true;
   } else if (fabs(dp) < fabs(*dx)) {
   // The third case: a lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
   // The cubic step is only used if the cubic tends to infinity in the direction of the step or if the minimum of the cubic is beyond *stp.
   // Otherwise the cubic step is defined to be either stpmin or stpmax.
   // The quadratic (secant) step is also computed and if the minimum is bracketed
   // then the step closest to *stx is taken, else the step farthest away is taken.
      *info = 3;
      bracktbound = *brackt;
      double theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
      double s = rmax2(fabs(theta), rmax2(fabs(*dx), fabs(dp)));
   // The case gamma == 0.0 only arises if the cubic does not tend to infinity in the direction of the step.
      double gamma = s * sqrt(rmax2(0.0, sqr(theta / s) - *dx / s * (dp / s)));
      if (*stp > *stx) gamma = -gamma;
      double r = (gamma - dp + theta) / (gamma - dp + gamma + *dx);
      double stpc = r < 0.0 && gamma != 0.0 ? *stp + r * (*stx - *stp) : *stp > *stx ? stmax : stmin;
      double stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
      stpf = (*brackt ? fabs(*stp - stpc) < fabs(*stp - stpq) : fabs(*stp - stpc) > fabs(*stp - stpq)) ? stpc : stpq;
   } else {
   // The fourth case: a lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
   // If the minimum is not bracketed, the step is either stpmin or stpmax, else the cubic step is taken.
      *info = 4;
      bracktbound = false;
      if (*brackt) {
         double theta = 3.0 * (fp - *fy) / (*sty - *stp) + *dy + dp;
         double s = rmax2(fabs(theta), rmax2(fabs(*dy), fabs(dp)));
         double gamma = s * sqrt(sqr(theta / s) - *dy / s * (dp / s));
         if (*stp > *sty) gamma = -gamma;
         double r = (gamma - dp + theta) / (gamma - dp + gamma + *dy);
         stpf = *stp + r * (*sty - *stp);
      } else stpf = *stp > *stx ? stmax : stmin;
   }
// Update the uncertainty interval.
// This update does not depend on the new step or the case analysis above.
   if (fp > *fx) {
      *sty = *stp;
      *fy = fp;
      *dy = dp;
   } else {
      if (sgnd < 0.0) {
         *sty = *stx;
         *fy = *fx;
         *dy = *dx;
      }
      *stx = *stp;
      *fx = fp;
      *dx = dp;
   }
// Compute the new step and safeguard it.
   *stp = rmax2(stmin, rmin2(stmax, stpf));
   if (bracktbound) *stp = *sty > *stx ? rmin2(*stx + 0.66 * (*sty - *stx), *stp) : rmax2(*stx + 0.66 * (*sty - *stx), *stp);
}

// The purpose of mcsrch() is to find a step which satisfies a sufficient decrease condition and a curvature condition.
// At each stage the subroutine updates an uncertainty interval with the state members stx and sty as the endpoints.
// The uncertainty interval is initially chosen so that it contains a minimizer of the modified function
//	F(x + *stp s) - F(x) - ftol *stp (F'(x)^T s).
// If a step is obtained for which the modified function has a non-positive function value and non-negative derivative,
// then the uncertainty interval is chosen so that it contains a minimizer of F(x + *stp s).
//
// The algorithm is designed to find a step which satisfies the sufficient decrease condition
//	F(x + *stp s) <= F(x) + ftol *stp (F'(x)^T s),
// and the curvature condition
//	|F'(x + *stp s)^T s| <= gtol |F'(x)^T s|.
// If ftol < gtol and if, for example, the function is bounded below, then there is always a step which satisfies both conditions.
// If no step can be found which satisfies both conditions,
// then the algorithm usually stops when rounding errors prevent further progress.
// In this case *stp only satisfies the sufficient decrease condition.
//
// IMPORTANT NOTES:
// *	This routine guarantees that it will stop at the last point where function value was calculated.
//	It won't make several additional function evaluations after finding a good point.
//	So if you store function evaluations requested by this routine,
//	you can be sure that last one is the point where we've stopped.
// *	When 0 < stpmax < stpmin, the algorithm will terminate with *info == 5 and *stp == stpmax.
// *	This algorithm guarantees that, if *info == 1 or *info == 5, then:
//	*	F(final_point) < F(initial_point) - strict inequality,
//	*	final_point != initial_point - after rounding to machine precision.
// *	When a non-descent direction is specified, algorithm stops with *info == 0, *stp == 0 and initial point at X[].
//
// Parameters and Inputs:
// *	n:	The number of variables; n > 0.
// *	x:	An n-vector for the base point for the line search, updated to x + *stp s.
// *	f:	The value, set to F(x) and updated to F(x + *stp s).
// *	g:	An n-vector, set to F'(x) and updated to F'(x + *stp s).
// *	s:	An n-vector indicating the search direction.
// *	*stp:	The step estimate; *stp >= 0; updated on output; accessed via the pointer stp.
// *	stpmin:	The minimum step size; stpmin >= 0 (represented here as the constant linmin_stpmin).
// *	stpmax:	The maximum step size; stpmax >= 0.
// *	xtol:	The tolerance for the relative width of the uncertainty interval; xtol >= 0.
// *	ftol:	The tolerance for sufficient decrease; ftol >= 0.
// *	gtol:	The tolerance for the directional derivative curvature condition; gtol >= 0.
// *	*info:	The return code; accessed via the pointer info:
//		0:	Improper inputs or parameters.
//		1:	The sufficient decrease condition and the directional derivative condition hold.
//		2:	The relative width of the uncertainty interval is at most xtol.
//		3:	The number of function calls has reached maxfev.
//		4:	The step is at the lower bound stpmin.
//		5:	The step is at the upper bound stpmax.
//		6:	Rounding errors prevent further progress.
//			There may not be a step which satisfies the sufficient decrease and curvature conditions.
//			The tolerances may be too small.
// *	*nfev:	The number of function calls; accessed via the pointer nfev.
// *	maxfev:	The number of function calls allowed for the algorithm; maxfev > 0.
// *	wa:	An n-vector for work space.
// *	state:	The algorithm state.
// *	*stage:	The algorithm stage; accessed via the pointer stage.
// Return Value:
// *	The condition (*stage != 0), indicating that iteration is in progress.
// Argonne National Laboratory. MINPACK Project. 1983 June.
// Jorge J. More', David J. Thuente.
bool mcsrch(ae_int_t n, RVector *x, double f, RVector *g, RVector *s, double *stp, double stpmax, double gtol, ae_int_t *info, ae_int_t *nfev, RVector *wa, linminstate *state, ae_int_t *stage) {
   const double ftol = 0.001, xtol = 100.0 * machineepsilon;
   const ae_int_t maxfev = 20;
   const double defstpmax = 1.0E+50;
   ae_int_t i;
   double v;
// Initialize.
   const double p5 = 0.5, p66 = 0.66;
   state->xtrapf = 4.0;
   const double zero = 0.0;
   if (stpmax == 0.0) {
      stpmax = defstpmax;
   }
   if (*stp < linmin_stpmin) {
      *stp = linmin_stpmin;
   }
   if (*stp > stpmax) {
      *stp = stpmax;
   }
// Manually threaded two-way signalling.
// Locals are zeroed out the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (*stage > 0) switch (*stage) {
   // case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      case 4: goto Resume4;
      default: goto Exit;
   }
Spawn:
// The main cycle.
#if 0
// Next.
   *stage = 2;
   Resume2:
#endif
   state->infoc = 1;
   *info = 0;
// Check the inputs and parameters for errors.
   if (stpmax < linmin_stpmin && stpmax > 0.0) {
      *info = 5;
      *stp = stpmax;
      goto Exit;
   }
   if (n <= 0 || *stp <= 0.0 || ftol < 0.0 || gtol < zero || xtol < zero || linmin_stpmin < zero || stpmax < linmin_stpmin || maxfev <= 0) {
      goto Exit;
   }
// Compute the initial gradient in the search direction and check that s is a descent direction.
   v = ae_v_dotproduct(g->xR, 1, s->xR, 1, n);
   state->dginit = v;
   if (state->dginit >= 0.0) {
      *stp = 0.0;
      goto Exit;
   }
// Initialize the local variables.
   state->brackt = false;
   state->stage1 = true;
   *nfev = 0;
   state->finit = f;
   state->dgtest = ftol * state->dginit;
   state->width = stpmax - linmin_stpmin;
   state->width1 = state->width / p5;
   ae_v_move(wa->xR, 1, x->xR, 1, n);
// The members stx, fx, dgx contain the values of the step, function, and directional derivative at the best step.
// The members sty, fy, dgy contain the values of the step, function, and derivative at the other endpoint of the uncertainty interval.
// The variables *stp, f and member dg contain the values of the step, function, and derivative at the current step.
   state->stx = 0.0;
   state->fx = state->finit;
   state->dgx = state->dginit;
   state->sty = 0.0;
   state->fy = state->finit;
   state->dgy = state->dginit;
   while (true) {
#if 0
   // Next.
      *stage = 3;
      Resume3:
#endif
   // Start the iteration.
   // Set the minimum and maximum steps to correspond to the present uncertainty interval.
      if (state->brackt) {
         if (state->stx < state->sty) {
            state->stmin = state->stx;
            state->stmax = state->sty;
         } else {
            state->stmin = state->sty;
            state->stmax = state->stx;
         }
      } else {
         state->stmin = state->stx;
         state->stmax = *stp + state->xtrapf * (*stp - state->stx);
      }
   // Force the step to be within the bounds stpmax and stpmin.
      if (*stp > stpmax) {
         *stp = stpmax;
      }
      if (*stp < linmin_stpmin) {
         *stp = linmin_stpmin;
      }
   // If an unusual termination is to occur then let *stp be the lowest point obtained so far.
      if (state->brackt && (*stp <= state->stmin || *stp >= state->stmax) || *nfev >= maxfev - 1 || state->infoc == 0 || state->brackt && state->stmax - state->stmin <= xtol * state->stmax) {
         *stp = state->stx;
      }
   // Evaluate the function and gradient at *stp and compute the directional derivative.
      ae_v_move(x->xR, 1, wa->xR, 1, n);
      ae_v_addd(x->xR, 1, s->xR, 1, n, *stp);
   // Next.
      *stage = 4; goto Pause; Resume4: ++*nfev;
      *info = 0;
      v = ae_v_dotproduct(g->xR, 1, s->xR, 1, n);
      state->dg = v;
      state->ftest1 = state->finit + *stp * state->dgtest;
   // Test for convergence.
      if (state->brackt && (*stp <= state->stmin || *stp >= state->stmax) || state->infoc == 0) {
         *info = 6;
      }
      if (*stp == stpmax && f < state->finit && f <= state->ftest1 && state->dg <= state->dgtest) {
         *info = 5;
      }
      if (*stp == linmin_stpmin && (f >= state->finit || f > state->ftest1 || state->dg >= state->dgtest)) {
         *info = 4;
      }
      if (*nfev >= maxfev) {
         *info = 3;
      }
      if (state->brackt && state->stmax - state->stmin <= xtol * state->stmax) {
         *info = 2;
      }
      if (f < state->finit && f <= state->ftest1 && SmallAtR(state->dg, -gtol * state->dginit)) {
         *info = 1;
      }
   // Check for termination.
      if (*info != 0) {
      // Check the guarantees provided by the function for *info == 1 or *info == 5.
         if (*info == 1 || *info == 5) {
            v = 0.0;
            for (i = 0; i < n; i++) {
               v += (wa->xR[i] - x->xR[i]) * (wa->xR[i] - x->xR[i]);
            }
            if (f >= state->finit || v == 0.0) {
               *info = 6;
            }
         }
         goto Exit;
      }
   // In the first stage we seek a step for which the modified function has a non-positive value and non-negative derivative.
      if (state->stage1 && f <= state->ftest1 && state->dg >= rmin2(ftol, gtol) * state->dginit) {
         state->stage1 = false;
      }
   // A modified function is used to predict the step only if we have not obtained a step
   // for which the modified function has a non-positive function value and non-negative derivative,
   // and if a lower function value has been obtained but the decrease is not sufficient.
      if (state->stage1 && f <= state->fx && f > state->ftest1) {
      // Define the modified function and derivative values.
         state->fm = f - *stp * state->dgtest;
         state->fxm = state->fx - state->stx * state->dgtest;
         state->fym = state->fy - state->sty * state->dgtest;
         state->dgm = state->dg - state->dgtest;
         state->dgxm = state->dgx - state->dgtest;
         state->dgym = state->dgy - state->dgtest;
      // Update the uncertainty interval and compute the new step.
         mcstep(&state->stx, &state->fxm, &state->dgxm, &state->sty, &state->fym, &state->dgym, stp, state->fm, state->dgm, &state->brackt, state->stmin, state->stmax, &state->infoc);
      // Reset the function and gradient values for f.
         state->fx = state->fxm + state->stx * state->dgtest;
         state->fy = state->fym + state->sty * state->dgtest;
         state->dgx = state->dgxm + state->dgtest;
         state->dgy = state->dgym + state->dgtest;
      } else {
      // Update the uncertainty interval and compute the new step.
         mcstep(&state->stx, &state->fx, &state->dgx, &state->sty, &state->fy, &state->dgy, stp, f, state->dg, &state->brackt, state->stmin, state->stmax, &state->infoc);
      }
   // Force a sufficient decrease in the size of the uncertainty interval.
      if (state->brackt) {
         if (!NearR(state->sty, state->stx, p66 * state->width1)) {
            *stp = state->stx + p5 * (state->sty - state->stx);
         }
         state->width1 = state->width;
         state->width = fabs(state->sty - state->stx);
      }
   }
Exit:
   *stage = 0;
   return false;
Pause:
   return true;
}

// These functions perform Armijo line search using  at  most  FMAX  function
// evaluations.  It  doesn't  enforce  some  kind  of  " sufficient decrease"
// criterion - it just tries different Armijo steps and returns optimum found
// so far.
//
// Optimization is done using F-RComm interface:
// * ArmijoCreate initializes State structure
//   (reusing previously allocated buffers)
// * ArmijoIteration is subsequently called
// * ArmijoResults returns results
//
// Inputs:
//     N       -   problem size
//     X       -   array[N], starting point
//     F       -   F(X+S*STP)
//     S       -   step direction, S > 0
//     STP     -   step length
//     STPMAX  -   maximum value for STP or zero (if no limit is imposed)
//     FMAX    -   maximum number of function evaluations
//     State   -   optimization state
// ALGLIB: Copyright 05.10.2010 by Sergey Bochkanov
void armijocreate(ae_int_t n, RVector *x, double f, RVector *s, double stp, double stpmax, ae_int_t fmax, armijostate *state) {
   if (state->x.cnt < n) {
      ae_vector_set_length(&state->x, n);
   }
   if (state->xbase.cnt < n) {
      ae_vector_set_length(&state->xbase, n);
   }
   if (state->s.cnt < n) {
      ae_vector_set_length(&state->s, n);
   }
   state->stpmax = stpmax;
   state->fmax = fmax;
   state->stplen = stp;
   state->fcur = f;
   state->n = n;
   ae_v_move(state->xbase.xR, 1, x->xR, 1, n);
   ae_v_move(state->s.xR, 1, s->xR, 1, n);
   state->PQ = -1;
}

// This is RComm-based search function
// ALGLIB: Copyright 05.10.2010 by Sergey Bochkanov
bool armijoiteration(armijostate *state) {
   const double armijofactor = 1.3;
   AutoS double v;
   AutoS ae_int_t n;
// Manually threaded two-way signalling.
// Locals are zeroed out the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0; case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      default: goto Exit;
   }
Spawn:
   if (state->stplen <= 0.0 || state->stpmax < 0.0 || state->fmax < 2) {
      state->info = 0;
      goto Exit;
   }
   if (state->stplen <= linmin_stpmin) {
      state->info = 4;
      goto Exit;
   }
   n = state->n;
   state->nfev = 0;
// We always need F
   state->needf = true;
// Bound StpLen
   if (state->stplen > state->stpmax && state->stpmax != 0.0) {
      state->stplen = state->stpmax;
   }
// Increase length
   v = state->stplen * armijofactor;
   if (v > state->stpmax && state->stpmax != 0.0) {
      v = state->stpmax;
   }
   ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
   ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
   state->PQ = 0; goto Pause; Resume0: state->nfev++;
   if (state->f < state->fcur) {
      state->stplen = v;
      state->fcur = state->f;
      while (true) {
      // test stopping conditions
         if (state->nfev >= state->fmax) {
            state->info = 3;
            goto Exit;
         }
         if (state->stplen >= state->stpmax) {
            state->info = 5;
            goto Exit;
         }
      // evaluate F
         v = state->stplen * armijofactor;
         if (v > state->stpmax && state->stpmax != 0.0) {
            v = state->stpmax;
         }
         ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
         ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
         state->PQ = 1; goto Pause; Resume1: state->nfev++;
      // make decision
         if (state->f < state->fcur) {
            state->stplen = v;
            state->fcur = state->f;
         } else {
            state->info = 1;
            goto Exit;
         }
      }
   }
// Decrease length
   v = state->stplen / armijofactor;
   ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
   ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
   state->PQ = 2; goto Pause; Resume2: state->nfev++;
   if (state->f < state->fcur) {
      state->stplen /= armijofactor;
      state->fcur = state->f;
      while (true) {
      // test stopping conditions
         if (state->nfev >= state->fmax) {
            state->info = 3;
            goto Exit;
         }
         if (state->stplen <= linmin_stpmin) {
            state->info = 4;
            goto Exit;
         }
      // evaluate F
         v = state->stplen / armijofactor;
         ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
         ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
         state->PQ = 3; goto Pause; Resume3: state->nfev++;
      // make decision
         if (state->f < state->fcur) {
            state->stplen /= armijofactor;
            state->fcur = state->f;
         } else {
            state->info = 1;
            goto Exit;
         }
      }
   }
// Nothing to be done
   state->info = 1;
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// Results of Armijo search
//
// Outputs:
//     INFO    -   on output it is set to one of the return codes:
//                 * 0     improper input params
//                 * 1     optimum step is found with at most FMAX evaluations
//                 * 3     FMAX evaluations were used,
//                         X contains optimum found so far
//                 * 4     step is at lower bound STPMIN
//                 * 5     step is at upper bound
//     STP     -   step length (in case of failure it is still returned)
//     F       -   function value (in case of failure it is still returned)
// ALGLIB: Copyright 05.10.2010 by Sergey Bochkanov
void armijoresults(armijostate *state, ae_int_t *info, double *stp, double *f) {
   *info = state->info;
   *stp = state->stplen;
   *f = state->fcur;
}

void linminstate_init(void *_p, bool make_automatic) {
}

void linminstate_copy(void *_dst, const void *_src, bool make_automatic) {
   linminstate *dst = (linminstate *)_dst;
   const linminstate *src = (const linminstate *)_src;
   dst->brackt = src->brackt;
   dst->stage1 = src->stage1;
   dst->infoc = src->infoc;
   dst->dg = src->dg;
   dst->dgm = src->dgm;
   dst->dginit = src->dginit;
   dst->dgtest = src->dgtest;
   dst->dgx = src->dgx;
   dst->dgxm = src->dgxm;
   dst->dgy = src->dgy;
   dst->dgym = src->dgym;
   dst->finit = src->finit;
   dst->ftest1 = src->ftest1;
   dst->fm = src->fm;
   dst->fx = src->fx;
   dst->fxm = src->fxm;
   dst->fy = src->fy;
   dst->fym = src->fym;
   dst->stx = src->stx;
   dst->sty = src->sty;
   dst->stmin = src->stmin;
   dst->stmax = src->stmax;
   dst->width = src->width;
   dst->width1 = src->width1;
   dst->xtrapf = src->xtrapf;
}

void linminstate_free(void *_p, bool make_automatic) {
}

void armijostate_init(void *_p, bool make_automatic) {
   armijostate *p = (armijostate *)_p;
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xbase, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->s, 0, DT_REAL, make_automatic);
}

void armijostate_copy(void *_dst, const void *_src, bool make_automatic) {
   armijostate *dst = (armijostate *)_dst;
   const armijostate *src = (const armijostate *)_src;
   dst->needf = src->needf;
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   dst->f = src->f;
   dst->n = src->n;
   ae_vector_copy(&dst->xbase, &src->xbase, make_automatic);
   ae_vector_copy(&dst->s, &src->s, make_automatic);
   dst->stplen = src->stplen;
   dst->fcur = src->fcur;
   dst->stpmax = src->stpmax;
   dst->fmax = src->fmax;
   dst->nfev = src->nfev;
   dst->info = src->info;
   dst->PQ = src->PQ;
}

void armijostate_free(void *_p, bool make_automatic) {
   armijostate *p = (armijostate *)_p;
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->xbase, make_automatic);
   ae_vector_free(&p->s, make_automatic);
}
} // end of namespace alglib_impl

// === NEARUNITYUNIT Package ===
namespace alglib_impl {
double nulog1p(double x) {
   double z;
   double lp;
   double lq;
   double result;
   z = 1.0 + x;
   if (z < 0.70710678118654752440 || z > 1.41421356237309504880) {
      result = log(z);
      return result;
   }
   z = x * x;
   lp = 4.5270000862445199635215E-5;
   lp = lp * x + 4.9854102823193375972212E-1;
   lp = lp * x + 6.5787325942061044846969E0;
   lp = lp * x + 2.9911919328553073277375E1;
   lp = lp * x + 6.0949667980987787057556E1;
   lp = lp * x + 5.7112963590585538103336E1;
   lp = lp * x + 2.0039553499201281259648E1;
   lq = 1.0000000000000000000000E0;
   lq = lq * x + 1.5062909083469192043167E1;
   lq = lq * x + 8.3047565967967209469434E1;
   lq = lq * x + 2.2176239823732856465394E2;
   lq = lq * x + 3.0909872225312059774938E2;
   lq = lq * x + 2.1642788614495947685003E2;
   lq = lq * x + 6.0118660497603843919306E1;
   z = -0.5 * z + x * (z * lp / lq);
   result = x + z;
   return result;
}

double nuexpm1(double x) {
   double r;
   double xx;
   double ep;
   double eq;
   double result;
   if (x < -0.5 || x > 0.5) {
      result = exp(x) - 1.0;
      return result;
   }
   xx = x * x;
   ep = 1.2617719307481059087798E-4;
   ep = ep * xx + 3.0299440770744196129956E-2;
   ep = ep * xx + 9.9999999999999999991025E-1;
   eq = 3.0019850513866445504159E-6;
   eq = eq * xx + 2.5244834034968410419224E-3;
   eq = eq * xx + 2.2726554820815502876593E-1;
   eq = eq * xx + 2.0000000000000000000897E0;
   r = x * ep;
   r /= eq - r;
   result = r + r;
   return result;
}

double nucosm1(double x) {
   double xx;
   double c;
   double result;
   if (x < -0.25 * pi || x > 0.25 * pi) {
      result = cos(x) - 1.0;
      return result;
   }
   xx = x * x;
   c = 4.7377507964246204691685E-14;
   c = c * xx - 1.1470284843425359765671E-11;
   c = c * xx + 2.0876754287081521758361E-9;
   c = c * xx - 2.7557319214999787979814E-7;
   c = c * xx + 2.4801587301570552304991E-5;
   c = c * xx - 1.3888888888888872993737E-3;
   c = c * xx + 4.1666666666666666609054E-2;
   result = -0.5 * xx + xx * xx * c;
   return result;
}
} // end of namespace alglib_impl

// === NTHEORY Package ===
namespace alglib_impl {
static bool ntheory_isprime(ae_int_t n) {
   ae_int_t p;
   bool result;
   result = false;
   p = 2;
   while (p * p <= n) {
      if (n % p == 0) {
         return result;
      }
      p++;
   }
   result = true;
   return result;
}

static ae_int_t ntheory_modmul(ae_int_t a, ae_int_t b, ae_int_t n) {
   ae_int_t t;
   double ra;
   double rb;
   ae_int_t result;
   ae_assert(a >= 0 && a < n, "ntheory_modmul: a < 0 or a >= n");
   ae_assert(b >= 0 && b < n, "ntheory_modmul: b < 0 or b >= n");
// Base cases
   ra = a;
   rb = b;
   if (b == 0 || a == 0) {
      result = 0;
      return result;
   }
   if (b == 1 || a == 1) {
      result = a * b;
      return result;
   }
   if (ra * rb == a * b) {
      result = a * b % n;
      return result;
   }
// Non-base cases
   if (b % 2 == 0) {
   // A*B = (A*(B/2)) * 2
   //
   // Product T == A*(B/2) is calculated recursively, product T*2 is
   // calculated as follows:
   // * result = T-N
   // * result = result+T
   // * if result < 0 then result = result+N
   //
   // In case integer result overflows, we generate exception
      t = ntheory_modmul(a, b / 2, n);
      result = t - n;
      result += t;
      if (result < 0) {
         result += n;
      }
   } else {
   // A*B = (A*(B div 2)) * 2 + A
   //
   // Product T == A*(B/2) is calculated recursively, product T*2 is
   // calculated as follows:
   // * result = T-N
   // * result = result+T
   // * if result < 0 then result = result+N
   //
   // In case integer result overflows, we generate exception
      t = ntheory_modmul(a, b / 2, n);
      result = t - n;
      result += t;
      if (result < 0) {
         result += n;
      }
      result -= n;
      result += a;
      if (result < 0) {
         result += n;
      }
   }
   return result;
}

static ae_int_t ntheory_modexp(ae_int_t a, ae_int_t b, ae_int_t n) {
   ae_int_t t;
   ae_int_t result;
   ae_assert(a >= 0 && a < n, "ntheory_modexp: a < 0 or a >= n");
   ae_assert(b >= 0, "ntheory_modexp: b < 0");
// Base cases
   if (b == 0) {
      result = 1;
      return result;
   }
   if (b == 1) {
      result = a;
      return result;
   }
// Non-base cases
   if (b % 2 == 0) {
      t = ntheory_modmul(a, a, n);
      result = ntheory_modexp(t, b / 2, n);
   } else {
      t = ntheory_modmul(a, a, n);
      result = ntheory_modexp(t, b / 2, n);
      result = ntheory_modmul(result, a, n);
   }
   return result;
}

void findprimitiverootandinverse(ae_int_t n, ae_int_t *proot, ae_int_t *invproot) {
   ae_int_t candroot;
   ae_int_t phin;
   ae_int_t q;
   ae_int_t f;
   bool allnonone;
   ae_int_t x;
   ae_int_t lastx;
   ae_int_t y;
   ae_int_t lasty;
   ae_int_t a;
   ae_int_t b;
   ae_int_t t;
   ae_int_t n2;
   *proot = 0;
   *invproot = 0;
   ae_assert(n >= 3, "findprimitiverootandinverse: n < 3");
   *proot = 0;
   *invproot = 0;
// check that N is prime
   ae_assert(ntheory_isprime(n), "findprimitiverootandinverse: n is not prime");
// Because N is prime, Euler totient function is equal to N-1
   phin = n - 1;
// Test different values of PRoot - from 2 to N-1.
// One of these values MUST be primitive root.
//
// For testing we use algorithm from Wiki (Primitive root modulo n):
// * compute phi(N)
// * determine the different prime factors of phi(N), say p1, ..., pk
// * for every element m of Zn*, compute m^(phi(N)/pi) mod N for i = 1..k
//   using a fast algorithm for modular exponentiation.
// * a number m for which these k results are all different from 1 is a
//   primitive root.
   for (candroot = 2; candroot < n; candroot++) {
   // We have current candidate root in CandRoot.
   //
   // Scan different prime factors of PhiN. Here:
   // * F is a current candidate factor
   // * Q is a current quotient - amount which was left after dividing PhiN
   //   by all previous factors
   //
   // For each factor, perform test mentioned above.
      q = phin;
      f = 2;
      allnonone = true;
      while (q > 1) {
         if (q % f == 0) {
            t = ntheory_modexp(candroot, phin / f, n);
            if (t == 1) {
               allnonone = false;
               break;
            }
            while (q % f == 0) {
               q /= f;
            }
         }
         f++;
      }
      if (allnonone) {
         *proot = candroot;
         break;
      }
   }
   ae_assert(*proot >= 2, "findprimitiverootandinverse: internal error (root not found)");
// Use extended Euclidean algorithm to find multiplicative inverse of primitive root
   x = 0;
   lastx = 1;
   y = 1;
   lasty = 0;
   a = *proot;
   b = n;
   while (b != 0) {
      q = a / b;
      t = a % b;
      a = b;
      b = t;
      t = lastx - q * x;
      lastx = x;
      x = t;
      t = lasty - q * y;
      lasty = y;
      y = t;
   }
   while (lastx < 0) {
      lastx += n;
   }
   *invproot = lastx;
// Check that it is safe to perform multiplication modulo N.
// Check results for consistency.
   n2 = (n - 1) * (n - 1);
   ae_assert(n2 / (n - 1) == n - 1, "findprimitiverootandinverse: internal error");
   ae_assert(*proot * *invproot / *proot == *invproot, "findprimitiverootandinverse: internal error");
   ae_assert(*proot * *invproot / *invproot == *proot, "findprimitiverootandinverse: internal error");
   ae_assert(*proot * *invproot % n == 1, "findprimitiverootandinverse: internal error");
}
} // end of namespace alglib_impl

// === FTBASE Package ===
// Depends on: APSERV, NTHEORY
namespace alglib_impl {
static const ae_int_t ftbase_coltype = 0, ftbase_coloperandscnt = 1, ftbase_coloperandsize = 2;
static const ae_int_t ftbase_colmicrovectorsize = 3, ftbase_colparam0 = 4, ftbase_colparam1 = 5;
static const ae_int_t ftbase_colparam2 = 6, ftbase_colparam3 = 7, ftbase_colscnt = 8;
static const ae_int_t ftbase_opend = 0, ftbase_opcomplexreffft = 1, ftbase_opbluesteinsfft = 2;
static const ae_int_t ftbase_opcomplexcodeletfft = 3, ftbase_opcomplexcodelettwfft = 4, ftbase_opradersfft = 5;
static const ae_int_t ftbase_opcomplextranspose = -1, ftbase_opcomplexfftfactors = -2;
static const ae_int_t ftbase_opstart = -3, ftbase_opjmp = -4, ftbase_opparallelcall = -5;
static const ae_int_t ftbase_maxradix = 6, ftbase_updatetw = 16;
static const ae_int_t ftbase_recursivethreshold = 0x400, ftbase_raderthreshold = 19;
static const ae_int_t ftbase_ftbasemaxsmoothfactor = 5;

// An optimistic cost estimate of the size n > 0 FFT, in units of 100 KFlops rounded down to nearest integer.
// Inputs:
//     N       -   the task size; N > 0.
// Result:
//     The cost in 100-KFlop units, rounded down to the nearest integer.
// NOTE:
//	A cost of less than 1 unit is rounded down to 0.
// ALGLIB: Copyright 08.04.2013 by Sergey Bochkanov
static ae_int_t ftbase_ftoptimisticestimate(ae_int_t n) {
   ae_assert(n > 0, "ftbase_ftoptimisticestimate: n <= 0");
   return ifloor(0.00005 * n * logbase2(n));
}

// Apply a complex reference FFT to an input/output vector.
// This is the VERY SLOW baseline operation, and is not meant for use in real life plans. :)
// Inputs:
//     ap   - the vector to be tranformed; assumed to be large enough for the result.
//     args - the operands count (see the description of fasttransformplan).
//     n    - the operand size (see the description of fasttransformplan).
//     buf  - a temporary vector of size at least 2 n args.
// Outputs:
//     ap   - the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplycomplexreffft(double *ap, ae_int_t args, ae_int_t n, RVector *buf) {
   const double twopi = 2.0 * pi;
   ae_assert(args >= 1, "ftbase_ftapplycomplexreffft: args < 1");
   ae_assert(n >= 1, "ftbase_ftapplycomplexreffft: n < 1");
   ae_int_t astride = 2 * n;
   double *bp = buf->xR;
   for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
      for (ae_int_t i = 0; i < n; i++) {
         double bx = 0.0, by = 0.0;
         for (ae_int_t k = 0; k < n; k++) {
            double ax = ap[2 * k], ay = ap[2 * k + 1];
            double omega = twopi * k * i / n, x = cos(omega), y = sin(omega);
            bx += x * ax + y * ay, by += x * ay - y * ax;
         }
         bp[2 * i] = bx, bp[2 * i + 1] = by;
      }
      for (ae_int_t i = 0; i < astride; i++) ap[i] = bp[i];
   }
}

// Apply a complex codelet FFT to an input/output vector.
// Inputs:
//     ap   - the vector to be transformed; assumed to be large enough for the result.
//     args - the operands count (see the description of fasttransformplan).
//     n    - the operand size (see the description of fasttransformplan).
// Outputs:
//     ap   - the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplycomplexcodeletfft(double *ap, ae_int_t args, ae_int_t n) {
   const double twopi = 2.0 * pi;
   const double sin30 = 0.5, sin60 = sqrt(0.75);
   const double root5 = sqrt(5.0), sin36 = sqrt((5.0 - root5)/8.0), sin72 = sqrt((5.0 + root5)/8.0);
   const double ca = 0.25, cb = 0.25 * root5, cc = sin72 + sin36, cd = sin72 - sin36;
   ae_assert(args >= 1, "ftbase_ftapplycomplexcodeletfft: args < 1");
   ae_assert(n >= 1, "ftbase_ftapplycomplexcodeletfft: n < 1");
// Hard-coded transforms for different n's.
   ae_assert(n <= ftbase_maxradix, "ftbase_ftapplycomplexcodeletfft: n > ftbase_maxradix");
   ae_int_t astride = 2 * n;
   switch (n) {
      case 2:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double ax = ap[0], ay = ap[1];
            double bx = ap[2], by = ap[3];
            ap[0] = ax + bx, ap[1] = ay + by;
            ap[2] = ax - bx, ap[3] = ay - by;
         }
      break;
      case 3:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double ax = ap[0], ay = ap[1];
            double bx = ap[2], by = ap[3];
            double cx = ap[4], cy = ap[5];
            double dx = bx + cx, dy = by + cy;
            double ex = sin60 * (by - cy), ey = sin60 * (cx - bx);
            double fx = ax - sin30 * dx, fy = ay - sin30 * dy;
            ap[0] = ax + dx, ap[1] = ay + dy;
            ap[2] = fx + ex, ap[3] = fy + ey;
            ap[4] = fx - ex, ap[5] = fy - ey;
         }
      break;
      case 4:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double ax = ap[0], ay = ap[1];
            double bx = ap[2], by = ap[3];
            double cx = ap[4], cy = ap[5];
            double dx = ap[6], dy = ap[7];
            double ex = ax + cx, ey = ay + cy;
            double fx = ax - cx, fy = ay - cy;
            double gx = bx + dx, gy = by + dy;
            double hx = by - dy, hy = dx - bx;
            ap[0] = ex + gx, ap[1] = ey + gy;
            ap[2] = fx + hx, ap[3] = fy + hy;
            ap[4] = ex - gx, ap[5] = ey - gy;
            ap[6] = fx - hx, ap[7] = fy - hy;
         }
      break;
      case 5:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double ax = ap[0], ay = ap[1];
            double bx = ap[2], by = ap[3];
            double cx = ap[4], cy = ap[5];
            double dx = ap[6], dy = ap[7];
            double ex = ap[8], ey = ap[9];
            double fx = bx + ex, fy = by + ey;
            double gx = bx - ex, gy = by - ey;
            double hx = dx + cx, hy = dy + cy;
            double ix = dx - cx, iy = dy - cy;
            double jx = fx + hx, jy = fy + hy;
            double kx = cb * (fx - hx), ky = cb * (fy - hy);
            double lx = sin72 * (gy + iy), ly = -sin72 * (gx + ix);
            double mx = lx - cc * iy, my = ly + cc * ix;
            double nx = lx - cd * gy, ny = ly + cd * gx;
            double ox = ax - ca * jx, oy = ay - ca * jy;
            double px = ox + kx, py = oy + ky;
            double qx = ox - kx, qy = oy - ky;
            ap[0] = ax + jx, ap[1] = ay + jy;
            ap[2] = px + mx, ap[3] = py + my;
            ap[4] = qx + nx, ap[5] = qy + ny;
            ap[6] = qx - nx, ap[7] = qy - ny;
            ap[8] = px - mx, ap[9] = py - my;
         }
      break;
      case 6:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double ax = ap[0], ay = ap[1];
            double bx = ap[2], by = ap[3];
            double cx = ap[4], cy = ap[5];
            double dx = ap[6], dy = ap[7];
            double ex = ap[8], ey = ap[9];
            double fx = ap[10], fy = ap[11];
            double gx = ax + dx, gy = ay + dy;
            double hx = ax - dx, hy = ay - dy;
            double ix = bx + ex, iy = by + ey;
            double jx = bx - ex, jy = by - ey;
            double kx = cx + fx, ky = cy + fy;
            double lx = cx - fx, ly = cy - fy;
            double mx = jx * sin30 + jy * sin60, my = -jx * sin60 + jy * sin30;
            double nx = -lx * sin30 + ly * sin60, ny = -lx * sin60 - ly * sin30;
            double ox = ix + kx, oy = iy + ky;
            double px = sin60 * (iy - ky), py = sin60 * (kx - ix);
            double qx = mx + nx, qy = my + ny;
            double rx = sin60 * (my - ny), ry = sin60 * (nx - mx);
            double sx = gx - sin30 * ox, sy = gy - sin30 * oy;
            double tx = hx - sin30 * qx, ty = hy - sin30 * qy;
            ap[0] = gx + ox, ap[1] = gy + oy;
            ap[2] = hx + qx, ap[3] = hy + qy;
            ap[4] = sx + px, ap[5] = sy + py;
            ap[6] = tx + rx, ap[7] = ty + ry;
            ap[8] = sx - px, ap[9] = sy - py;
            ap[10] = tx - rx, ap[11] = ty - ry;
         }
      break;
   }
}

// Apply the complex "integrated" codelet FFT to an input/output vector.
// The "integrated" codelet differs from the "normal" one in the following ways:
// * it can work with mun > 1,
// * hence, it can be used in the Cooley-Tukey FFT without transpositions,
// * it performs inlined multiplication by twiddle factors of Cooley-Tukey FFT with n2 == mun/2.
// Inputs:
//     ap   - the vector to be transformed; assumed to be large enough for the result.
//     args - the operand count (see the description of fasttransformplan).
//     n    - the operand size (see the description of fasttransformplan).
//     mun  - the microvector size; assumed to be at least 1.
// Outputs:
//     ap   - the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplycomplexcodelettwfft(double *ap, ae_int_t args, ae_int_t n, ae_int_t mun) {
   const double twopi = 2.0 * pi;
   const double sin30 = 0.5, sin60 = sqrt(0.75);
   const double root5 = sqrt(5.0), sin36 = sqrt((5.0 - root5) / 8.0), sin72 = sqrt((5.0 + root5) / 8.0);
   const double ca = 0.25, cb = 0.25 * root5, cc = sin72 + sin36, cd = sin72 - sin36;
   ae_assert(args >= 1, "ftbase_ftapplycomplexcodelettwfft: args < 1");
   ae_assert(n >= 1, "ftbase_ftapplycomplexcodelettwfft: n < 1");
   ae_assert(mun >= 1, "ftbase_ftapplycomplexcodelettwfft: mun < 1");
   ae_assert(mun % 2 == 0, "ftbase_ftapplycomplexcodelettwfft: mun is not even");
   ae_int_t m = mun / 2, astride = n * mun;
// Hard-coded transforms for different n's.
   ae_assert(n <= ftbase_maxradix, "ftbase_ftapplycomplexcodelettwfft: n > ftbase_maxradix");
   double omega = -twopi / (n * m), twx = cos(omega), twy = sin(omega);
   switch (n) {
      case 2:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double *ap0 = ap, *ap1 = ap0 + mun;
            double tw1x = 1.0, tw1y = 0.0;
            for (ae_int_t mvidx = 0; mvidx < m; ap0 += 2, ap1 += 2, mvidx++) {
               double ax = ap0[0], ay = ap0[1];
               double bx = ap1[0], by = ap1[1];
               double cx = ax - bx, cy = ay - by;
               ap0[0] = ax + bx, ap0[1] = ay + by;
               ap1[0] = cx * tw1x - cy * tw1y, ap1[1] = cy * tw1x + cx * tw1y;
               if ((mvidx + 1) % ftbase_updatetw == 0) {
                  double omegamu = -twopi * (mvidx + 1) / (n * m);
                  tw1x = cos(omegamu), tw1y = sin(omegamu);
               } else {
                  double x = tw1x * twx - tw1y * twy, y = tw1x * twy + tw1y * twx;
                  tw1x = x, tw1y = y;
               }
            }
         }
      break;
      case 3:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double *ap0 = ap, *ap1 = ap0 + mun, *ap2 = ap1 + mun;
            double tw1x = 1.0, tw1y = 0.0;
            for (ae_int_t mvidx = 0; mvidx < m; ap0 += 2, ap1 += 2, ap2 += 2, mvidx++) {
               double ax = ap0[0], ay = ap0[1];
               double bx = ap1[0], by = ap1[1];
               double cx = ap2[0], cy = ap2[1];
               double dx = bx + cx, dy = by + cy;
               double ex = sin60 * (by - cy), ey = sin60 * (cx - bx);
               double fx = ax - sin30 * dx, fy = ay - sin30 * dy;
               double tw2x = tw1x * tw1x - tw1y * tw1y, tw2y = 2.0 * tw1x * tw1y;
               bx = fx + ex, by = fy + ey;
               cx = fx - ex, cy = fy - ey;
               ap0[0] = ax + dx, ap0[1] = ay + dy;
               ap1[0] = bx * tw1x - by * tw1y, ap1[1] = by * tw1x + bx * tw1y;
               ap2[0] = cx * tw2x - cy * tw2y, ap2[1] = cy * tw2x + cx * tw2y;
               if ((mvidx + 1) % ftbase_updatetw == 0) {
                  double omegamu = -twopi * (mvidx + 1) / (n * m);
                  tw1x = cos(omegamu), tw1y = sin(omegamu);
               } else {
                  double x = tw1x * twx - tw1y * twy, y = tw1x * twy + tw1y * twx;
                  tw1x = x, tw1y = y;
               }
            }
         }
      break;
      case 4:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double *ap0 = ap, *ap1 = ap0 + mun, *ap2 = ap1 + mun, *ap3 = ap2 + mun;
            double tw1x = 1.0, tw1y = 0.0;
            for (ae_int_t mvidx = 0; mvidx < m; ap0 += 2, ap1 += 2, ap2 += 2, ap3 += 2, mvidx++) {
               double ax = ap0[0], ay = ap0[1];
               double bx = ap1[0], by = ap1[1];
               double cx = ap2[0], cy = ap2[1];
               double dx = ap3[0], dy = ap3[1];
               double ex = ax + cx, ey = ay + cy;
               double fx = ax - cx, fy = ay - cy;
               double gx = bx + dx, gy = by + dy;
               double hx = by - dy, hy = dx - bx;
               double tw2x = tw1x * tw1x - tw1y * tw1y, tw2y = 2.0 * tw1x * tw1y;
               double tw3x = tw1x * tw2x - tw1y * tw2y, tw3y = tw1x * tw2y + tw1y * tw2x;
               bx = fx + hx, by = fy + hy;
               dx = fx - hx, dy = fy - hy;
               cx = ex - gx, cy = ey - gy;
               ap0[0] = ex + gx, ap0[1] = ey + gy;
               ap1[0] = bx * tw1x - by * tw1y, ap1[1] = by * tw1x + bx * tw1y;
               ap2[0] = cx * tw2x - cy * tw2y, ap2[1] = cy * tw2x + cx * tw2y;
               ap3[0] = dx * tw3x - dy * tw3y, ap3[1] = dy * tw3x + dx * tw3y;
               if ((mvidx + 1) % ftbase_updatetw == 0) {
                  double omegamu = -twopi * (mvidx + 1) / (n * m);
                  tw1x = cos(omegamu), tw1y = sin(omegamu);
               } else {
                  double x = tw1x * twx - tw1y * twy, y = tw1x * twy + tw1y * twx;
                  tw1x = x, tw1y = y;
               }
            }
         }
      break;
      case 5:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double *ap0 = ap, *ap1 = ap0 + mun, *ap2 = ap1 + mun, *ap3 = ap2 + mun, *ap4 = ap3 + mun;
            double tw1x = 1.0, tw1y = 0.0;
            for (ae_int_t mvidx = 0; mvidx < m; ap0 += 2, ap1 += 2, ap2 += 2, ap3 += 2, ap4 += 2, mvidx++) {
               double ax = ap0[0], ay = ap0[1];
               double bx = ap1[0], by = ap1[1];
               double cx = ap2[0], cy = ap2[1];
               double dx = ap3[0], dy = ap3[1];
               double ex = ap4[0], ey = ap4[1];
               double fx = bx + ex, fy = by + ey;
               double gx = bx - ex, gy = by - ey;
               double hx = dx + cx, hy = dy + cy;
               double ix = dx - cx, iy = dy - cy;
               double jx = fx + hx, jy = fy + hy;
               double kx = cb * (fx - hx), ky = cb * (fy - hy);
               double lx = sin72 * (gy + iy), ly = -sin72 * (gx + ix);
               double mx = lx - cc * iy, my = ly + cc * ix;
               double nx = lx - cd * gy, ny = ly + cd * gx;
               double ox = ax - ca * jx, oy = ay - ca * jy;
               double px = ox + kx, py = oy + ky;
               double qx = ox - kx, qy = oy - ky;
               double tw2x = tw1x * tw1x - tw1y * tw1y, tw2y = 2.0 * tw1x * tw1y;
               double tw3x = tw1x * tw2x - tw1y * tw2y, tw3y = tw1x * tw2y + tw1y * tw2x;
               double tw4x = tw2x * tw2x - tw2y * tw2y, tw4y = tw2x * tw2y + tw2y * tw2x;
               bx = px + mx, by = py + my;
               ex = px - mx, ey = py - my;
               cx = qx + nx, cy = qy + ny;
               dx = qx - nx, dy = qy - ny;
               ap0[0] = ax + jx, ap0[1] = ay + jy;
               ap1[0] = bx * tw1x - by * tw1y, ap1[1] = bx * tw1y + by * tw1x;
               ap2[0] = cx * tw2x - cy * tw2y, ap2[1] = cx * tw2y + cy * tw2x;
               ap3[0] = dx * tw3x - dy * tw3y, ap3[1] = dx * tw3y + dy * tw3x;
               ap4[0] = ex * tw4x - ey * tw4y, ap4[1] = ex * tw4y + ey * tw4x;
               if ((mvidx + 1) % ftbase_updatetw == 0) {
                  double omegamu = -twopi * (mvidx + 1) / (n * m);
                  tw1x = cos(omegamu), tw1y = sin(omegamu);
               } else {
                  double x = tw1x * twx - tw1y * twy, y = tw1x * twy + tw1y * twx;
                  tw1x = x, tw1y = y;
               }
            }
         }
      break;
      case 6:
         for (ae_int_t arg = 0; arg < args; ap += astride, arg++) {
            double *ap0 = ap, *ap1 = ap0 + mun, *ap2 = ap1 + mun, *ap3 = ap2 + mun, *ap4 = ap3 + mun, *ap5 = ap4 + mun;
            double tw1x = 1.0, tw1y = 0.0;
            for (ae_int_t mvidx = 0; mvidx < m; ap0 += 2, ap1 += 2, ap2 += 2, ap3 += 2, ap4 += 2, ap5 += 2, mvidx++) {
               double ax = ap0[0], ay = ap0[1];
               double bx = ap1[0], by = ap1[1];
               double cx = ap2[0], cy = ap2[1];
               double dx = ap3[0], dy = ap3[1];
               double ex = ap4[0], ey = ap4[1];
               double fx = ap5[0], fy = ap5[1];
               double gx = ax + dx, gy = ay + dy;
               double hx = ax - dx, hy = ay - dy;
               double ix = bx + ex, iy = by + ey;
               double jx = bx - ex, jy = by - ey;
               double kx = cx + fx, ky = cy + fy;
               double lx = cx - fx, ly = cy - fy;
               double mx = jx * sin30 + jy * sin60, my = -jx * sin60 + jy * sin30;
               double nx = -lx * sin30 + ly * sin60, ny = -lx * sin60 - ly * sin30;
               double ox = ix + kx, oy = iy + ky;
               double px = sin60 * (iy - ky), py = sin60 * (kx - ix);
               double qx = mx + nx, qy = my + ny;
               double rx = sin60 * (my - ny), ry = sin60 * (nx - mx);
               double sx = gx - sin30 * ox, sy = gy - sin30 * oy;
               double tx = hx - sin30 * qx, ty = hy - sin30 * qy;
               double tw2x = tw1x * tw1x - tw1y * tw1y, tw2y = 2.0 * tw1x * tw1y;
               double tw3x = tw1x * tw2x - tw1y * tw2y, tw3y = tw1x * tw2y + tw1y * tw2x;
               double tw4x = tw2x * tw2x - tw2y * tw2y, tw4y = 2.0 * tw2x * tw2y;
               double tw5x = tw3x * tw2x - tw3y * tw2y, tw5y = tw3x * tw2y + tw3y * tw2x;
               cx = sx + px, cy = sy + py;
               ex = sx - px, ey = sy - py;
               bx = hx + qx, by = hy + qy;
               dx = tx + rx, dy = ty + ry;
               fx = tx - rx, fy = ty - ry;
               ap0[0] = gx + ox, ap0[1] = gy + oy;
               ap1[0] = bx * tw1x - by * tw1y, ap1[1] = by * tw1x + bx * tw1y;
               ap2[0] = cx * tw2x - cy * tw2y, ap2[1] = cy * tw2x + cx * tw2y;
               ap3[0] = dx * tw3x - dy * tw3y, ap3[1] = dy * tw3x + dx * tw3y;
               ap4[0] = ex * tw4x - ey * tw4y, ap4[1] = ey * tw4x + ex * tw4y;
               ap5[0] = fx * tw5x - fy * tw5y, ap5[1] = fy * tw5x + fx * tw5y;
               if ((mvidx + 1) % ftbase_updatetw == 0) {
                  double omegamu = -twopi * (mvidx + 1) / (n * m);
                  tw1x = cos(omegamu), tw1y = sin(omegamu);
               } else {
                  double x = tw1x * twx - tw1y * twy, y = tw1x * twy + tw1y * twx;
                  tw1x = x, tw1y = y;
               }
            }
         }
      break;
   }
}

// Twiddle factors calculation.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_ffttwcalc(double *ap, ae_int_t n1, ae_int_t n2) {
// Multiplication by twiddle factors for the complex Cooley-Tukey FFT with the factoring n = n1 n2.
// In the following description, 1^x stands for exp(2 pi i x) = cos(2 pi x) + i sin(2 pi x).
// The naive solution to this problem is:
//	for (k in [1,n2), j in [1,n1)) {
//		n = k n1 + j;
//		x + iy = ap[2n] + i ap[2n + 1];
//		twx + i twy = 1^{-k j/(n1 n2)};
//		ap[2n] + i ap[2n + 1] = (x + iy) (twx + i twy) /* (x twx - y twy) + i (x twy + y twx).
//	}
// However, there are more efficient solutions.
//
// Each pass of the inner multiplies an entry of ap by 1^{-kj/n} = (1^{-k/n})^j.
// Exponentiation thus becomes multiplication 1^{-(j + 1) k/n} = 1^{-jk/n} 1^{-k/n}, with 1^{-k/n} tabulated only once.
// Moreover, since 1^{-k/n} = (1^{-1/n})^k, we also have 1^{-(k + 1)/n} = 1^{-k/n} 1^{-1/n}.
//
// In our loop, we use the following variables:
// * twbase = twbasex + i twbasey = 1^{1/n} - the base twiddle factor,
// * twrow = twrowx + i twrowy = 1^{k/n} - the base row twiddle,
// * tw = twx + i twy = 1^{jk/N} - the current twiddle factor.
// The inner loop multiplies tw by twrow, the outer loop multiplies twrow by twbase.
   ae_assert(ftbase_updatetw >= 2, "ftbase_ffttwcalc: internal error - ftbase_updatetw < 2");
   ae_int_t updatetw2 = ftbase_updatetw / 2, n1q = n1 / 2, n1r = n1 % 2, n = n1 * n2;
   const double twopi = 2.0 * pi;
   double omega = -twopi / n, twbasex = cos(omega), twbasey = sin(omega);
   double twrowx = 1.0, twrowy = 0.0;
   for (ae_int_t k = 0; k < n2; k++) {
   // Initialize tw = 1.
      double twx = 1.0, twy = 0.0;
   // The n1-point block is separated into 2-point chunks and, for odd n1, a residual 1-point chunk.
   // The unrolled loop is several times faster.
      for (ae_int_t j2 = 0; j2 < n1q; j2++) {
      // ap[0] + i ap[1] *= tw;
         double x = ap[0], y = ap[1];
         ap[0] = x * twx - y * twy, ap[1] = x * twy + y * twx;
      // tw *= twrow;
         x = twx * twrowx - twy * twrowy, y = twx * twrowy + twy * twrowx;
         twx = x, twy = y;
      // ap[2] + i ap[3] *= tw;
         x = ap[2], y = ap[3];
         ap[2] = x * twx - y * twy, ap[3] = x * twy + y * twx;
         ap += 4;
      // Conditionally update tw.
         if ((j2 + 1) % updatetw2 == 0 && j2 < n1q - 1) {
         // Re-twiddle: tw = 1^{-(2 k (j2 + 1)/n)};
            double omega = -twopi * 2.0 * k * (j2 + 1) / n;
            twx = cos(omega), twy = sin(omega);
         } else {
         // Update: tw *= twrow;
            double x = twx * twrowx - twy * twrowy, y = twx * twrowy + twy * twrowx;
            twx = x, twy = y;
         }
      }
      if (n1r == 1) {
      // Handle the residual chunk.
         double x = ap[0], y = ap[1];
         ap[0] = x * twx - y * twy, ap[1] = x * twy + y * twx;
         ap += 2;
      }
      if (k < n2 - 1) {
      // Update twrow.
         if ((k + 1) % ftbase_updatetw == 0) {
         // Re-twiddle: twrow = 1^{-(k + 1)/n};
            double omega = -twopi * (k + 1) / n;
            twrowx = cos(omega), twrowy = sin(omega);
         } else {
         // Update: twrow *= tw;
            double x = twrowx * twbasex - twrowy * twbasey, y = twrowx * twbasey + twrowy * twbasex;
            twrowx = x, twrowy = y;
         }
      }
   }
}

// An estimate of the FLOP count for the FFT, only really accurate when n is a power of 2;
// based on the operations count for the *perfect* FFT and the relative inefficiency of the algorithm actually used.
// The estimates are badly wrong for non-power-of-2 n's.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
double ftbasegetflopestimate(ae_int_t n) {
   const double ftbaseinefficiencyfactor = 1.3;
   return ftbaseinefficiencyfactor * (4.0 * n * logbase2(n) - 6.0 * n + 8.0);
}

#if 0
// A recurrent subroutine for the (non-existent) routine ftbase_internalreallintranspose().
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_fftirltrec(double *ap, ae_int_t astride, double *bp, ae_int_t bstride, ae_int_t m, ae_int_t n) {
   if (m == 0 || n == 0) return;
   else if (imax2(m, n) <= 8) {
      for (ae_int_t i = 0; i < m; i++) {
         ae_int_t idx1 = i, idx2 = i * astride;
         for (ae_int_t j = 0; j < n; j++) {
            bp[idx1] = ap[idx2];
            idx1 += bstride;
            idx2++;
         }
      }
   } else if (n > m) {
   // New partition:
   // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
   //                                  ( B2 )
      ae_int_t n1 = n / 2;
      if (n - n1 >= 8 && n1 % 8 != 0) n1 += 8 - n1 % 8;
      ae_assert(n - n1 > 0, "ftbase_fftirltrec: Assertion n > n1 failed");
      ftbase_fftirltrec(ap, astride, bp, bstride, m, n1);
      ftbase_fftirltrec(ap + n1, astride, bp + n1 * bstride, bstride, m, n - n1);
   } else {
   // New partition:
   // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
   //                     ( A2 )
      ae_int_t m1 = m / 2;
      if (m - m1 >= 8 && m1 % 8 != 0) m1 += 8 - m1 % 8;
      ae_assert(m - m1 > 0, "ftbase_fftirltrec: Assertion m > m1 failed");
      ftbase_fftirltrec(ap, astride, bp, bstride, m1, n);
      ftbase_fftirltrec(ap + m1 * astride, astride, bp + m1, bstride, m - m1, n);
   }
}
#endif

// A recurrent subroutine for ftbase_internalcomplexlintranspose().
// Write A^T to B, where:
// * A is formatted as an m x n complex matrix ap in real/imaginary value pairs, with stride astride,
// * B is formatted as an n x m complex matrix bp in real/imaginary value pairs, with stride bstride.
// Strides are in complex number units, i.e. in real/imaginary pairs.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_ffticltrec(double *ap, ae_int_t astride, double *bp, ae_int_t bstride, ae_int_t m, ae_int_t n) {
   if (m == 0 || n == 0) return;
   else if (imax2(m, n) <= 8) {
      for (ae_int_t i = 0, m2 = 2 * bstride; i < m; i++) {
         double *ai = ap + 2 * i * astride, *bi = bp + 2 * i;
         for (ae_int_t j = 0; j < n; ai += 2, bi += m2, j++) {
            bi[0] = ai[0], bi[1] = ai[1];
         }
      }
   } else if (n > m) {
   // New partition:
   // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
   //                                  ( B2 )
      ae_int_t n1 = n / 2;
      if (n - n1 >= 8 && n1 % 8 != 0) n1 += 8 - n1 % 8;
      ae_assert(n - n1 > 0, "ftbase_ffticltrec: Assertion n > n1 failed");
      ftbase_ffticltrec(ap, astride, bp, bstride, m, n1);
      ftbase_ffticltrec(ap + 2 * n1, astride, bp + 2 * n1 * bstride, bstride, m, n - n1);
   } else {
   // New partition:
   // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
   //                     ( A2 )
      ae_int_t m1 = m / 2;
      if (m - m1 >= 8 && m1 % 8 != 0) m1 += 8 - m1 % 8;
      ae_assert(m - m1 > 0, "ftbase_ffticltrec: Assertion m > m1 failed");
      ftbase_ffticltrec(ap, astride, bp, bstride, m1, n);
      ftbase_ffticltrec(ap + 2 * m1 * astride, astride, bp + 2 * m1, bstride, m - m1, n);
   }
}

// Linear transpose: transpose a complex matrix stored as a vector.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_internalcomplexlintranspose(double *ap, ae_int_t m, ae_int_t n, RVector *buf) {
   double *bp = buf->xR;
   ftbase_ffticltrec(ap, n, bp, m, m, n);
   ae_v_move(ap, 1, bp, 1, 2 * m * n);
}

// Factor the task size n into tasks of sizes n1 and n2, such that:
// * n1 == 0 == n2, if n is prime or composite of size <= ftbase_maxradix,
// * 1 <= n1 <= n2, with n1 n2 = n, if n is composite or of size > ftbase_maxradix.
// Inputs:
//     n       -   the task size; n > 0.
//     isroot  -   whether or not the task is the root task: the first one in a sequence.
// Outputs:
//     n1, n2  -   the sizes of the smaller tasks, returned respectively as *n1p and *n2p.
// ALGLIB: Copyright 08.04.2013 by Sergey Bochkanov
static void ftbase_ftfactorize(ae_int_t n, bool isroot, ae_int_t *n1p, ae_int_t *n2p) {
   ae_assert(n > 0, "ftbase_ftfactorize: n <= 0");
// Small n.
   if (n <= ftbase_maxradix) {
      *n2p = *n1p = 0;
      return;
   }
// Large n: recursive split.
   if (n > ftbase_recursivethreshold) {
      ae_int_t k = iceil(sqrt(n)) + 1;
      ae_assert(k * k >= n, "ftbase_ftfactorize: internal error during recursive factorization");
      for (ae_int_t j = k; j >= 2; j--) if (n % j == 0) {
         *n1p = imin2(n / j, j), *n2p = imax2(n / j, j);
         return;
      }
   }
   ae_int_t n1 = 0, n2 = 0;
// n > ftbase_maxradix: try to find a good codelet.
   for (ae_int_t j = ftbase_maxradix; j >= 2; j--) if (n % j == 0) {
      n1 = j, n2 = n / j;
      break;
   }
// If there are no good codelets, try to factor n into product of ANY primes.
   if (n1 * n2 != n)
      for (ae_int_t j = 2; j < n; j++)
         if (n % j == 0) {
            n1 = j, n2 = n / j;
            break;
         } else if (j * j > n) break;
// Normalize.
   if (n1 > n2)
      *n1p = n2, *n2p = n1;
   else
      *n1p = n1, *n2p = n2;
}

// Find a good factoring n = n1 n2; returned as *n1p = n1 and *n2p = n2.
// Usually n1 <= n2, but not always: small n's may be exception.
// if n1 != 1 then n2 != 1.
// Factoring is chosen depending on the task type and codelets we have.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
void ftbasefactorize(ae_int_t n, ae_int_t tasktype, ae_int_t *n1p, ae_int_t *n2p) {
   const ae_int_t ftbasecodeletrecommended = 5;
   ae_int_t n1 = 0, n2 = 0;
// Look for a good codelet size.
   if (n1 * n2 != n)
      for (ae_int_t j = ftbasecodeletrecommended; j >= 2; j--) if (n % j == 0) {
         n1 = j, n2 = n / j;
         break;
      }
// Try to factor n.
   if (n1 * n2 != n)
      for (ae_int_t j = ftbasecodeletrecommended + 1; j < n; j++) if (n % j == 0) {
         n1 = j, n2 = n / j;
         break;
      }
// n is prime. :(
   if (n1 * n2 != n) n1 = 1, n2 = n;
// Normalize.
   if (n2 == 1 && n1 != 1) n2 = n1, n1 = 1;
   *n1p = n1, *n2p = n2;
}

// Is the number smooth?
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
bool ftbaseissmooth(ae_int_t n) {
   for (ae_int_t i = 2; i <= ftbase_ftbasemaxsmoothfactor; i++)
      for (; n % i == 0; n /= i);
   return n == 1;
}

// A recurrent subroutine for ftbasefindsmooth().
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_ftbasefindsmoothrec(ae_int_t n, ae_int_t seed, ae_int_t leastfactor, ae_int_t *bestp) {
   ae_assert(ftbase_ftbasemaxsmoothfactor <= 5, "ftbase_ftbasefindsmoothrec: internal error!");
   if (seed >= n) {
      *bestp = imin2(*bestp, seed);
      return;
   }
   if (leastfactor <= 2) ftbase_ftbasefindsmoothrec(n, seed * 2, 2, bestp);
   if (leastfactor <= 3) ftbase_ftbasefindsmoothrec(n, seed * 3, 3, bestp);
   if (leastfactor <= 5) ftbase_ftbasefindsmoothrec(n, seed * 5, 5, bestp);
}

// The smallest smooth (divisible only by 2, 3, 5) number at least as large as n and 2.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
ae_int_t ftbasefindsmooth(ae_int_t n) {
   ae_int_t best = 2;
   for (; best < n; best *= 2);
   ftbase_ftbasefindsmoothrec(n, 1, 2, &best);
   return best;
}

// The smallest smooth (divisible only by 2, 3, 5) even number at least as large as n and 2.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
ae_int_t ftbasefindsmootheven(ae_int_t n) {
   ae_int_t best = 2;
   for (; best < n; best *= 2);
   ftbase_ftbasefindsmoothrec(n, 2, 2, &best);
   return best;
}

// The space requirements for the n-point FFT: highly dependent on details of the different FFTs used by this unit,
// So, every time the algorithm is changed this function has to be rewritten.
// Inputs:
//     n          - the transform length.
//     *precrsize - must be set to zero.
//     *precisize - must be set to zero.
// Outputs:
//     *precrsize - the number of real temporaries required for transformation.
//     *precisize - the number of integer temporaries required for transformation.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftdeterminespacerequirements(ae_int_t n, ae_int_t *precrsize, ae_int_t *precisize) {
// Determine the required sizes of the precomputed real and integer buffers.
// This stage of the code is highly dependent on the internals ftbase_ftcomplexfftplanrec()
// and must be kept in sync with possible changes in the internals of the plan generation function.
//
// The buffer size is determined as follows:
// *	n is factored,
// *	we factor out anything which is less or equal to ftbase_maxradix,
// *	a prime factor f > ftbase_raderthreshold requires 4 ftbasefindsmooth(2f - 1) real entries
//	for the precomputed Quantities for Bluestein's transformation,
// *	a prime factor f <= ftbase_raderthreshold requires 2(f - 1) + ESTIMATE(F - 1) precomputed storage.
   ae_int_t ncur = n;
   for (ae_int_t i = 2; i <= ftbase_maxradix; i++)
      for (; ncur % i == 0; ncur /= i);
   for (ae_int_t f = 2; f <= ncur; f++)
      for (; ncur % f == 0; ncur /= f)
         if (f > ftbase_raderthreshold)
            *precrsize += 4 * ftbasefindsmooth(2 * f - 1);
         else {
            *precrsize += 2 * (f - 1);
            ftbase_ftdeterminespacerequirements(f - 1, precrsize, precisize);
         }
}

// Forward declaration for indirect recursion.
static void ftbase_ftapplysubplan(fasttransformplan *plan, ae_int_t subplan, double *ap, ae_int_t aoffset, RVector *buf, ae_int_t repcnt);

// Apply the complex Bluestein's FFT to an input/output vector.
// Inputs:
//     plan        -   the transformation plan.
//     ap          -   the vector to be transformed; assumed to be large enough for the plan to work.
//     args        -   number of repeated operands (length n each).
//     n           -   the original complex data length.
//     m           -   the padded complex data length.
//     precoffs    -   the offset of the precomputed data for the plan.
//     subplan     -   the position of the length-m FFT subplan which is used by the transformation.
//     bufb        -   a temporary buffer, at least 2m elements.
//     bufc        -   a temporary buffer, at least 2m elements.
// Outputs:
//     ap          -   the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftbluesteinsfft(fasttransformplan *plan, double *ap, ae_int_t args, ae_int_t n, ae_int_t m, ae_int_t precoffs, ae_int_t subplan, RVector *bufb, RVector *bufc) {
   double *bp = bufb->xR;
   for (ae_int_t arg = 0; arg < args; ap += 2 * n, arg++) {
   // bp = ap zp*, padding bp with zeros.
   // NOTE:
   //	z[k] == 1^{ik^2/2n}
      double *ai = ap, *bi = bp, *zi = plan->precr.xR + precoffs;
      for (ae_int_t i = 0; i < n; ai += 2, bi += 2, zi += 2, i++) {
         double ax = ai[0], ay = ai[1];
         double zx = zi[0], zy = zi[1];
         bi[0] = zx * ax + zy * ay, bi[1] = zx * ay - zy * ax;
      }
      for (ae_int_t i = 2 * n; i < 2 * m; i++) bp[i] = 0.0;
   // bp = convolution of bp and zp, using the precomputed DFT of zp in plan.
      ftbase_ftapplysubplan(plan, subplan, bp, 0, bufc, 1);
      bi = bp, zi = plan->precr.xR + precoffs + 2 * m;
      for (ae_int_t i = 0; i < m; bi += 2, zi += 2, i++) {
         double ax = bi[0], ay = bi[1];
         double zx = zi[0], zy = zi[1];
         bi[0] = ax * zx - ay * zy, bi[1] = -(ax * zy + ay * zx);
      }
      ftbase_ftapplysubplan(plan, subplan, bp, 0, bufc, 1);
   // Post-processing:
   //     ap = zp* bp*/m = zp* (bp*/m)
   // where zp* comes from the Bluestein's FFT algorithm and bp*/m corresponds to the last stage of the inverse DFT,
      ai = ap, bi = bp, zi = plan->precr.xR + precoffs;
      for (ae_int_t i = 0; i < n; ai += 2, bi += 2, zi += 2, i++) {
         double bx = bi[0], by = bi[1];
         double zx = zi[0], zy = zi[1];
         ai[0] = (zx * bx - zy * by) / m, ai[1] = -(zx * by + zy * bx) / m;
      }
   }
}

// Apply the complex Rader's FFT to an input/output vector.
// Inputs:
//     plan        -   the transformation plan.
//     ap          -   the vector to be transformed; assumed to be large enough for the plan to work.
//     aoffset     -   an offset in [0, PlanLength), within the large PlanLength-subarray of the chunk to process.
//     args        -   the number of length-n repeated operands.
//     n           -   the original complex data length.
//     subplan     -   the position of the (n-1)-point FFT subplan to be used by transformation.
//     rq          -   a primitive root modulo n.
//     riq         -   the inverse of the primitive root modulo n.
//     precoffs    -   the offset of the precomputed data for the plan.
//     buf         -   a temporary vector.
// Outputs:
//     ap          -   the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftradersfft(fasttransformplan *plan, double *ap, ae_int_t aoffset, ae_int_t args, ae_int_t n, ae_int_t subplan, ae_int_t rq, ae_int_t riq, ae_int_t precoffs, RVector *buf) {
   ae_assert(args >= 1, "ftbase_ftradersfft: args < 1");
   double *bp = buf->xR + aoffset;
// Process the operands.
   for (ae_int_t arg = 0; arg < args; ap += 2 * n, bp += 2 * n, arg++) {
   // Fill QA.
      double rx = ap[0], ry = ap[1];
      double x0 = rx, y0 = ry;
      double *ai = ap, *bi = bp;
      for (ae_int_t q = 0, kq = 1; q < n - 1; bi += 2, q++, kq = kq * rq % n) {
         rx += bi[0] = ai[2 * kq], ry += bi[1] = ai[2 * kq + 1];
      }
      ai = ap, bi = bp;
      for (ae_int_t q = 0; q < n - 1; ai += 2, bi += 2, q++) {
         ai[0] = bi[0], ai[1] = bi[1];
      }
   // Convolution.
      ftbase_ftapplysubplan(plan, subplan, ap, aoffset, buf, 1);
      double *zi = plan->precr.xR + precoffs;
      ai = ap;
      for (ae_int_t i = 0; i < n - 1; ai += 2, zi += 2, i++) {
         double ax = ai[0], ay = ai[1];
         double zx = zi[0], zy = zi[1];
         ai[0] = ax * zx - ay * zy, ai[1] = -(ax * zy + ay * zx);
      }
      ftbase_ftapplysubplan(plan, subplan, ap, aoffset, buf, 1);
      ai = ap;
      for (ae_int_t i = 0; i < n - 1; ai += 2, i++) {
         ai[0] /= n - 1, ai[1] /= 1 - n;
      }
   // Result.
      ai = ap, bi = bp;
      bp[0] = rx, bp[1] = ry;
      for (ae_int_t q = 0, kiq = 1; q < n - 1; ai += 2, q++, kiq = kiq * riq % n) {
         bi[2 * kiq] = x0 + ai[0], bi[2 * kiq + 1] = y0 + ai[1];
      }
      ai = ap, bi = bp;
      for (ae_int_t q = 0; q < n; ai += 2, bi += 2, q++) {
         ai[0] = bi[0], ai[1] = bi[1];
      }
   }
}

// Apply a subplan to an input/output vector.
// Inputs:
//     plan    - the transformation plan.
//     subplan - the subplan index.
//     ap      - the vector to be transformed; assumed to be large enough for the plan to work.
//     aoffset - an offset in [0, PlanLength), within the large PlanLength-subarray of the chunk to process.
//     buf     - a temporary buffer of length at least the plan length (without taking into account repcnt).
//     repcnt  - the repetition count (transformation is repeatedly applied to subsequent subvectors).
// Outputs:
//     plan    - the transformation plan: temporary buffers can be modified, though the plan itself is unchanged and can be reused.
//     ap      - the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplysubplan(fasttransformplan *plan, ae_int_t subplan, double *ap, ae_int_t aoffset, RVector *buf, ae_int_t repcnt) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   RefObj(RVector, bufb);
   RefObj(RVector, bufc);
   ae_int_t **plantab = plan->entries.xyZ;
   ae_assert(plantab[subplan][ftbase_coltype] == ftbase_opstart, "ftbase_ftapplysubplan: incorrect subplan header");
   for (ae_int_t rowidx = subplan + 1; plantab[rowidx][ftbase_coltype] != ftbase_opend; ) {
      ae_int_t *planrow = plantab[rowidx];
      ae_int_t operation = planrow[ftbase_coltype], args = repcnt * planrow[ftbase_coloperandscnt];
      ae_int_t n = planrow[ftbase_coloperandsize], mun = planrow[ftbase_colmicrovectorsize];
      ae_int_t param0 = planrow[ftbase_colparam0], param1 = planrow[ftbase_colparam1];
      ae_int_t param2 = planrow[ftbase_colparam2], param3 = planrow[ftbase_colparam3];
      switch (operation) {
      // Jump:
         case ftbase_opjmp: rowidx += param0; break;
      // Parallel Call:
         case ftbase_opparallelcall: {
         // * initially check for consistency between the sup and sub plans,
         // * call ftbase_ftapplysubplan() to split the parallel plan into several parallel tasks.
            ae_int_t *subrow = plantab[rowidx + param0];
            ae_int_t supn = n * mun;
            ae_int_t subn = subrow[ftbase_coloperandscnt] * subrow[ftbase_coloperandsize] * subrow[ftbase_colmicrovectorsize];
            ae_assert(subrow[ftbase_coltype] == ftbase_opstart, "ftbase_ftapplysubplan: incorrect child subplan header");
            ae_assert(supn == subn, "ftbase_ftapplysubplan: incorrect child subplan header");
            ae_int_t chunksize = imax2(ftbase_recursivethreshold / subn, 1);
            for (ae_int_t i = 0; i < args; i += chunksize) {
               chunksize = imin2(chunksize, args - i);
               ftbase_ftapplysubplan(plan, rowidx + param0, ap + i * subn, aoffset + i * subn, buf, chunksize);
            }
            rowidx++;
         }
         break;
      // Reference Complex FFT:
         case ftbase_opcomplexreffft:
            ae_assert(mun == 2, "ftbase_ftapplysubplan: mun != 2 for complex FFT");
            ftbase_ftapplycomplexreffft(ap, args, n, buf);
            rowidx++;
         break;
      // Codelet FFT:
         case ftbase_opcomplexcodeletfft:
            ae_assert(mun == 2, "ftbase_ftapplysubplan: mun != 2 for codelet FFT");
            ftbase_ftapplycomplexcodeletfft(ap, args, n);
            rowidx++;
         break;
      // Integrated Codelet FFT:
         case ftbase_opcomplexcodelettwfft:
            ftbase_ftapplycomplexcodelettwfft(ap, args, n, mun);
            rowidx++;
         break;
      // Bluestein's FFT:
         case ftbase_opbluesteinsfft:
            ae_assert(mun == 2, "ftbase_ftapplysubplan: mun != 2 for Bluestein's FFT");
            ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufb);
            ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufc);
            ftbase_ftbluesteinsfft(plan, ap, args, n, param0, param2, rowidx + param1, bufb, bufc);
            ae_shared_pool_recycle(&plan->bluesteinpool, &_bufb);
            ae_shared_pool_recycle(&plan->bluesteinpool, &_bufc);
            rowidx++;
         break;
      // Rader's FFT:
         case ftbase_opradersfft:
            ftbase_ftradersfft(plan, ap, aoffset, args, n, rowidx + param0, param1, param2, param3, buf);
            rowidx++;
         break;
      // Complex Twiddle Factors:
         case ftbase_opcomplexfftfactors: {
            ae_assert(mun == 2, "ftbase_ftapplysubplan: mun != 2 for twiddle FFT");
            ae_int_t n1 = param0, n2 = n / n1;
            for (ae_int_t i = 0; i < args; i++) ftbase_ffttwcalc(ap + i * n * 2, n1, n2);
            rowidx++;
         }
         break;
      // Complex Transposition:
         case ftbase_opcomplextranspose: {
            ae_assert(mun == 2, "ftbase_ftapplysubplan: mun != 2 for complex transposition");
            ae_int_t n1 = param0, n2 = n / n1;
            for (ae_int_t i = 0; i < args; i++) ftbase_internalcomplexlintranspose(ap + 2 * n * i, n1, n2, buf);
            rowidx++;
         }
         break;
      // Error:
         default: ae_assert(false, "ftbase_ftapplysubplan: unexpected plan type");
      }
   }
   ae_frame_leave();
}

// Apply a transformation plan to an input/output vector.
// Inputs:
//     plan   - the transformation plan
//     a      - the vector to be transformed; assumed to be large enough for plan to work
//     offsa  - the offset of the subvector to process
//     repcnt - the repetition count; the transformation is repeatedly applied to subsequent subvectors.
// Outputs:
//     plan   - the transformation plan: temporary buffers may be modified, though the plan itself is unchanged and can be reused.
//     a      - the transformed vector.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
void ftapplyplan(fasttransformplan *plan, RVector *a, ae_int_t offsa/* = 0*/, ae_int_t repcnt/* = 1*/) {
   double *ap = a->xR + offsa;
   ae_int_t *planrow = plan->entries.xyZ[0];
   ae_int_t plansize = planrow[ftbase_coloperandscnt] * planrow[ftbase_coloperandsize] * planrow[ftbase_colmicrovectorsize];
   for (ae_int_t i = 0; i < repcnt; i++) ftbase_ftapplysubplan(plan, 0, ap + plansize * i, 0, &plan->buffer, 1);
}

// Push an entry with up to 4 parameters to a plan grid, resizing it, if necessary.
// Inputs:
//     gridp       -   the plan grid to be updated.
//     *rowptr     -   the index which points past the last entry generated so far; accessed via the pointer rowptr.
//     etype       -   the entry type.
//     eopcnt      -   the operand count.
//     eopsize     -   the operand size.
//     emcvsize    -   the microvector size.
//     eparam0     -   parameter 0.
//     eparam1     -   parameter 1.
//     eparam2     -   parameter 2.
//     eparam3     -   parameter 3.
// Outputs:
//     *rowptr     -   the updated index.
//     gridp       -   the updated plan grid.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftpushentry(ZMatrix *gridp, ae_int_t *rowptr, ae_int_t etype, ae_int_t eopcnt, ae_int_t eopsize, ae_int_t emcvsize, ae_int_t eparam0, ae_int_t eparam1 = -1, ae_int_t eparam2 = 0, ae_int_t eparam3 = 0) {
   ae_int_t row = (*rowptr)++;
   if (row >= gridp->rows) imatrixresize(gridp, imax2(2 * gridp->rows, 1), ftbase_colscnt);
   ae_int_t *planrow = gridp->xyZ[row];
   planrow[ftbase_coltype] = etype;
   planrow[ftbase_coloperandscnt] = eopcnt;
   planrow[ftbase_coloperandsize] = eopsize;
   planrow[ftbase_colmicrovectorsize] = emcvsize;
   planrow[ftbase_colparam0] = eparam0;
   planrow[ftbase_colparam1] = eparam1;
   planrow[ftbase_colparam2] = eparam2;
   planrow[ftbase_colparam3] = eparam3;
}

#if 0
// Forward reference to an indirect recursive call. //(@) Already declared externally.
void ftcomplexfftplan(ae_int_t n, ae_int_t k, fasttransformplan *plan);
#endif

// Precompute the complex Bluestein's FFT and write the data to a vector.
// The caller must ensure that the vector is large enough.
// Inputs:
//     n  - the original size of the transform.
//     m  - the size of the "padded" Bluestein's transform.
//     zp - the preallocated array; assumed to be of length at least 4m.
// Outputs:
//     zp - the 4M-vector, modified so that:
//          * zp[0:2m-1] = (1^{i k^2/2n}: 0 <= i < 2m)
//          * zp[2m:4m-1] = the DFT of Z.
// NOTE:
//	This function performs an internal M-point FFT.
//	It allocates a temporary plan which is freed after leaving this function.
// ALGLIB: Copyright 08.05.2013 by Sergey Bochkanov
static void ftbase_ftprecomputebluesteinsfft(ae_int_t n, ae_int_t m, double *zp) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   NewObj(fasttransformplan, plan);
// Fill first half of PrecR with b[k] = exp(i*pi*k^2/N)
   for (ae_int_t i = 0; i < 2 * m; i++) zp[i] = 0.0;
   for (ae_int_t i = 0; i < n; i++) {
      double omega = pi / n * i * i;
      ae_int_t i1 = (m - i) % m;
      zp[2 * i1] = zp[2 * i] = cos(omega), zp[2 * i1 + 1] = zp[2 * i + 1] = sin(omega);
   }
// Precomputed the FFT.
   ftcomplexfftplan(m, 1, &plan);
   for (ae_int_t i = 0; i < 2 * m; i++) zp[2 * m + i] = zp[i];
   ftbase_ftapplysubplan(&plan, 0, zp + 2 * m, 0, &plan.buffer, 1);
   ae_frame_leave();
}

// Precompute the complex Rader's FFT and write the data to a vector.
// The caller must ensure that the vector is large enough for the result.
// Inputs:
//     n   - the original size of the transform; before reduction to n - 1.
//     rq  - a primitive root modulo n.
//     riq - the inverse of the primitive root modulo n.
//     zp  - a preallocated vector of size at least 2 (n - 1).
// Outputs:
//     zp  - a 2(n-1)-vector for the FFT of Rader's factors.
// NOTE:
//	This function performs an internal (N-1)-point FFT.
//	It allocates a temporary plan which is freed after leaving this function.
// ALGLIB: Copyright 08.05.2013 by Sergey Bochkanov
static void ftbase_ftprecomputeradersfft(ae_int_t n, ae_int_t rq, ae_int_t riq, double *zp) {
   const double twopi = 2.0 * pi;
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   NewObj(fasttransformplan, plan);
// Fill zp with Rader factors, perform an FFT.
   for (ae_int_t q = 0, kiq = 1; q < n - 1; q++, kiq = kiq * riq % n) {
      double omega = -twopi * kiq / n;
      zp[2 * q] = cos(omega), zp[2 * q + 1] = sin(omega);
   }
   ftcomplexfftplan(n - 1, 1, &plan);
   ftbase_ftapplysubplan(&plan, 0, zp, 0, &plan.buffer, 1);
   ae_frame_leave();
}

// A recurrent function called by ftcomplexfftplan() and other functions.
// It recursively builds a transformation plan.
// Inputs:
//     N          - the FFT length (in complex numbers); N >= 1.
//     K          - the number of repetitions; K >= 1.
//     issub      - if true, the plan generator inserts OpStart/OpEnd in the plan header/footer.
//     istop      - if true, the plan generator assumes that it is the topmost plan:
//                  * it may use a global buffer for transpositions and there is no other plan which executes in parallel.
//     rowptr     - an index which points past the last entry generated so far.
//     bluesteinn - the amount of storage (in real numbers) required for Bluestein buffer; stored as *bluesteinnp.
//     precrptr   - a pointer to the unused part of the precomputed real buffer (plan->precr):
//                  * when this function stores some data to the precomputed buffer, it advances the pointer,
//                  * the function will assert that plan->precr has enough space for the data before actually writing to the buffer.
//                  * the caller must allocate enough space before calling this function.
//     preciptr   - a pointer to the unused part of the precomputed integer buffer (plan->preci):
//                  * when this function stores some data to the precomputed buffer, it advances the pointer,
//                  * the function will assert that plan->preci has enough space for the data before actually writing to the buffer,
//                  * the caller must allocate enough space before calling this function.
//     plan       - the plan; generated so far.
// Outputs:
//     rowptr     - the updated pointer: advanced by the number of entries generated by the function.
//     bluesteinn - the updated amount: may be increased, but may never be decreased; returned as *bluesteinnp.
// NOTE: in case istop is true, issub is also must be true.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftcomplexfftplanrec(ae_int_t n, ae_int_t k, bool issub, bool istop, ae_int_t *rowptr, ae_int_t *bluesteinnp, ae_int_t *precrptr, ae_int_t *preciptr, fasttransformplan *plan) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   ae_assert(n > 0, "ftbase_ftcomplexfftplanrec: n <= 0");
   ae_assert(k > 0, "ftbase_ftcomplexfftplanrec: k <= 0");
   ae_assert(!istop || issub, "ftbase_ftcomplexfftplanrec: issub is inconsistent with istop");
   ae_int_t n1, n2;
   ftbase_ftfactorize(n, false, &n1, &n2);
   ZMatrix *gridp = &plan->entries;
   if (istop && n > ftbase_recursivethreshold) { // Try to generate the "topmost" plan.
      if (n1 * n2 == 0) { // Prime-factor DFT with Bluestein's FFT.
      // Determine the size of the Bluestein buffer.
         ae_int_t m = ftbasefindsmooth(2 * n - 1);
         *bluesteinnp = imax2(2 * m, *bluesteinnp);
      // Generate the plan.
         ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opbluesteinsfft, k, n, 2, m, 2, *precrptr, 0);
         ae_int_t row0 = *rowptr;
         ftbase_ftpushentry(gridp, rowptr, ftbase_opjmp, 0, 0, 0, 0);
         ftbase_ftcomplexfftplanrec(m, 1, true, true, rowptr, bluesteinnp, precrptr, preciptr, plan);
         gridp->xyZ[row0][ftbase_colparam0] = *rowptr - row0;
         ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
      // Fill the precomputed buffer.
         ftbase_ftprecomputebluesteinsfft(n, m, plan->precr.xR + *precrptr);
      // Update the pointer to the precomputed area.
         *precrptr += 4 * m;
      } else { // Composite FFT: done recursively with Cooley-Tukey, which uses the global buffer instead of a local one.
         ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         ae_int_t row0 = *rowptr;
         ftbase_ftpushentry(gridp, rowptr, ftbase_opparallelcall, k * n2, n1, 2, 0, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n2);
         ae_int_t row2 = *rowptr;
         ftbase_ftpushentry(gridp, rowptr, ftbase_opparallelcall, k * n1, n2, 2, 0, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
         ae_int_t row1 = *rowptr;
         ftbase_ftcomplexfftplanrec(n1, 1, true, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         gridp->xyZ[row0][ftbase_colparam0] = row1 - row0;
         ae_int_t row3 = *rowptr;
         ftbase_ftcomplexfftplanrec(n2, 1, true, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         gridp->xyZ[row2][ftbase_colparam0] = row3 - row2;
      }
   } else {
   // Generate subordinate plans.
   // A local (shared) buffer is needed and used, and the buffer size is updated:
   // ANY plan will need at least 2 N temporaries, additional requirements can be applied later.
      if (n1 * n2 == 0) { // Small-N or prime-factor DFT.
         if (n <= ftbase_maxradix) { // Small-N DFT.
            if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
            ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplexcodeletfft, k, n, 2, 0);
            if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
         } else if (n <= ftbase_raderthreshold) { // Prime-factor DFT with Rader's FFT.
            ae_int_t m = n - 1;
            if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
            ae_int_t gq, giq;
            findprimitiverootandinverse(n, &gq, &giq);
            ftbase_ftpushentry(gridp, rowptr, ftbase_opradersfft, k, n, 2, 2, gq, giq, *precrptr);
            ftbase_ftprecomputeradersfft(n, gq, giq, plan->precr.xR + *precrptr);
            *precrptr += 2 * (n - 1);
            ae_int_t row0 = *rowptr;
            ftbase_ftpushentry(gridp, rowptr, ftbase_opjmp, 0, 0, 0, 0);
            ftbase_ftcomplexfftplanrec(m, 1, true, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
            gridp->xyZ[row0][ftbase_colparam0] = *rowptr - row0;
            if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
         } else { // Prime-factor DFT with Bluestein's FFT.
            ae_int_t m = ftbasefindsmooth(2 * n - 1);
            *bluesteinnp = imax2(2 * m, *bluesteinnp);
            if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
            ftbase_ftpushentry(gridp, rowptr, ftbase_opbluesteinsfft, k, n, 2, m, 2, *precrptr, 0);
            ftbase_ftprecomputebluesteinsfft(n, m, plan->precr.xR + *precrptr);
            *precrptr += 4 * m;
            ae_int_t row0 = *rowptr;
            ftbase_ftpushentry(gridp, rowptr, ftbase_opjmp, 0, 0, 0, 0);
            ftbase_ftcomplexfftplanrec(m, 1, true, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
            gridp->xyZ[row0][ftbase_colparam0] = *rowptr - row0;
            if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
         }
      } else if (n1 <= ftbase_maxradix) { // Small n1 DFT with Handle Cooley-Tukey FFT.
      // Specialized transformation for small n1:
      // * n2 short in-place FFT's, each n1-point, with integrated twiddle factors,
      // * n1 long FFT's,
      // * a final transposition.
         if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplexcodelettwfft, k, n1, 2 * n2, 0);
         ftbase_ftcomplexfftplanrec(n2, k * n1, false, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
      } else if (n <= ftbase_recursivethreshold) { // Small n DFT: general "flat" Cooley-Tukey FFT.
      // General code for large n1/n2, "flat" version without explicit recurrence
      // (nested subplans are inserted directly into the body of the plan).
         if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         ftbase_ftcomplexfftplanrec(n1, k * n2, false, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n2);
         ftbase_ftcomplexfftplanrec(n2, k * n1, false, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
      } else { // Large n DFT: general "recursive" Cooley-Tukey DFT - nested subplans are separated from the plan body.
      // Generate the sup-plan.
         if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         ae_int_t row0 = *rowptr;
         ftbase_ftpushentry(gridp, rowptr, ftbase_opparallelcall, k * n2, n1, 2, 0, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1);
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n2);
         ae_int_t row2 = *rowptr;
         ftbase_ftpushentry(gridp, rowptr, ftbase_opparallelcall, k * n1, n2, 2, 0, ftbase_ftoptimisticestimate(n));
         ftbase_ftpushentry(gridp, rowptr, ftbase_opcomplextranspose, k, n, 2, n1);
         if (issub) ftbase_ftpushentry(gridp, rowptr, ftbase_opend, k, n, 2, 0);
      // Generate the sub-plans, linking them to their sup-plans.
         ae_int_t row1 = *rowptr;
         ftbase_ftcomplexfftplanrec(n1, 1, true, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         gridp->xyZ[row0][ftbase_colparam0] = row1 - row0;
         ae_int_t row3 = *rowptr;
         ftbase_ftcomplexfftplanrec(n2, 1, true, false, rowptr, bluesteinnp, precrptr, preciptr, plan);
         gridp->xyZ[row2][ftbase_colparam0] = row3 - row2;
      }
   }
   ae_frame_leave();
}

// Generate an FFT plan for k length-n complex DFT's.
// Inputs:
//     n    - the size of the complex vector to be processed by the FFT; n >= 1.
//     k    - the number of repetitions; k >= 1.
// Outputs:
//     plan - the plan.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
void ftcomplexfftplan(ae_int_t n, ae_int_t k, fasttransformplan *plan) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   SetObj(fasttransformplan, plan);
   NewObj(RVector, bluesteinbuf);
// Verify the parameters.
   ae_assert(n > 0, "ftcomplexfftplan: n <= 0");
   ae_assert(k > 0, "ftcomplexfftplan: k <= 0");
// Determine the required sizes of the precomputed real and integer buffers.
// This stage of the code is highly dependent on the internals of ftbase_ftcomplexfftplanrec()
// and must be kept in sync with possible changes in the internals of the plan generation function.
//
// The buffer size is determined as follows:
// *	n is factored,
// *	we factor out anything which is less or equal to ftbase_maxradix,
// *	a prime factor f > ftbase_raderthreshold requires 4 ftbasefindsmooth(2f - 1) real entries
//	for the precomputed Quantities for Bluestein's transformation,
// *	a prime factor f <= ftbase_raderthreshold does NOT require precomputed storage.
   ae_int_t precrsize = 0, precisize = 0;
   ftbase_ftdeterminespacerequirements(n, &precrsize, &precisize);
   if (precrsize > 0) ae_vector_set_length(&plan->precr, precrsize);
   if (precisize > 0) ae_vector_set_length(&plan->preci, precisize);
// Generate the plan.
   ae_int_t row = 0, precr = 0, preci = 0, bluesteinsize = 1;
   ae_vector_set_length(&plan->buffer, 2 * n * k);
   ftbase_ftcomplexfftplanrec(n, k, true, true, &row, &bluesteinsize, &precr, &preci, plan);
   ae_vector_set_length(&bluesteinbuf, bluesteinsize);
   ae_shared_pool_set_seed(&plan->bluesteinpool, &bluesteinbuf, sizeof bluesteinbuf, RVector_init, RVector_copy, RVector_free);
// Check that the actual amount of precomputed space used by the transformation plan
// EXACTLY matches the amount of space allocated by us.
   ae_assert(precr == precrsize, "ftcomplexfftplan: internal error (precr != precrsize)");
   ae_assert(preci == precisize, "ftcomplexfftplan: internal error (preci != precisize)");
   ae_frame_leave();
}

void fasttransformplan_init(void *_p, bool make_automatic) {
   fasttransformplan *p = (fasttransformplan *)_p;
   ae_matrix_init(&p->entries, 0, 0, DT_INT, make_automatic);
   ae_vector_init(&p->buffer, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->precr, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->preci, 0, DT_REAL, make_automatic);
   ae_shared_pool_init(&p->bluesteinpool, make_automatic);
}

void fasttransformplan_copy(void *_dst, const void *_src, bool make_automatic) {
   fasttransformplan *dst = (fasttransformplan *)_dst;
   const fasttransformplan *src = (const fasttransformplan *)_src;
   ae_matrix_copy(&dst->entries, &src->entries, make_automatic);
   ae_vector_copy(&dst->buffer, &src->buffer, make_automatic);
   ae_vector_copy(&dst->precr, &src->precr, make_automatic);
   ae_vector_copy(&dst->preci, &src->preci, make_automatic);
   ae_shared_pool_copy(&dst->bluesteinpool, &src->bluesteinpool, make_automatic);
}

void fasttransformplan_free(void *_p, bool make_automatic) {
   fasttransformplan *p = (fasttransformplan *)_p;
   ae_matrix_free(&p->entries, make_automatic);
   ae_vector_free(&p->buffer, make_automatic);
   ae_vector_free(&p->precr, make_automatic);
   ae_vector_free(&p->preci, make_automatic);
   ae_shared_pool_free(&p->bluesteinpool, make_automatic);
}
} // end of namespace alglib_impl

// === HPCCORES Package ===
namespace alglib_impl {
// Stub function.
// ALGLIB Routine: Copyright 14.06.2013 by Sergey Bochkanov
static bool hpccores_hpcpreparechunkedgradientx(RVector *weights, ae_int_t wcount, RVector *hpcbuf) {
#ifndef ALGLIB_INTERCEPTS_SSE2
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hpcpreparechunkedgradientx(weights, wcount, hpcbuf);
#endif
}

// Stub function.
// ALGLIB Routine: Copyright 14.06.2013 by Sergey Bochkanov
static bool hpccores_hpcfinalizechunkedgradientx(RVector *buf, ae_int_t wcount, RVector *grad) {
#ifndef ALGLIB_INTERCEPTS_SSE2
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hpcfinalizechunkedgradientx(buf, wcount, grad);
#endif
}

// Prepares HPC compuations  of  chunked  gradient with HPCChunkedGradient().
// You  have to call this function  before  calling  HPCChunkedGradient() for
// a new set of weights. You have to call it only once, see example below:
//
// HOW TO PROCESS DATASET WITH THIS FUNCTION:
//     Grad = 0
//     HPCPrepareChunkedGradient(Weights, WCount, NTotal, NOut, Buf)
//     foreach chunk-of-dataset do
//         HPCChunkedGradient(...)
//     HPCFinalizeChunkedGradient(Buf, Grad)
//
void hpcpreparechunkedgradient(RVector *weights, ae_int_t wcount, ae_int_t ntotal, ae_int_t nin, ae_int_t nout, mlpbuffers *buf) {
   ae_int_t i;
   ae_int_t batch4size;
   ae_int_t chunksize;
   chunksize = 4;
   batch4size = 3 * chunksize * ntotal + chunksize * (2 * nout + 1);
   matrixsetlengthatleast(&buf->xy, chunksize, nin + nout);
   matrixsetlengthatleast(&buf->xy2, chunksize, nin + nout);
   vectorsetlengthatleast(&buf->xyrow, nin + nout);
   vectorsetlengthatleast(&buf->x, nin);
   vectorsetlengthatleast(&buf->y, nout);
   vectorsetlengthatleast(&buf->desiredy, nout);
   vectorsetlengthatleast(&buf->batch4buf, batch4size);
   vectorsetlengthatleast(&buf->hpcbuf, wcount);
   vectorsetlengthatleast(&buf->g, wcount);
   if (!hpccores_hpcpreparechunkedgradientx(weights, wcount, &buf->hpcbuf)) {
      for (i = 0; i < wcount; i++) {
         buf->hpcbuf.xR[i] = 0.0;
      }
   }
   buf->wcount = wcount;
   buf->ntotal = ntotal;
   buf->nin = nin;
   buf->nout = nout;
   buf->chunksize = chunksize;
}

// Finalizes HPC compuations  of  chunked gradient with HPCChunkedGradient().
// You  have to call this function  after  calling  HPCChunkedGradient()  for
// a new set of weights. You have to call it only once, see example below:
//
// HOW TO PROCESS DATASET WITH THIS FUNCTION:
//     Grad = 0
//     HPCPrepareChunkedGradient(Weights, WCount, NTotal, NOut, Buf)
//     foreach chunk-of-dataset do
//         HPCChunkedGradient(...)
//     HPCFinalizeChunkedGradient(Buf, Grad)
//
void hpcfinalizechunkedgradient(mlpbuffers *buf, RVector *grad) {
   ae_int_t i;
   if (!hpccores_hpcfinalizechunkedgradientx(&buf->hpcbuf, buf->wcount, grad)) {
      for (i = 0; i < buf->wcount; i++) {
         grad->xR[i] += buf->hpcbuf.xR[i];
      }
   }
}

// Fast kernel for chunked gradient.
//
bool hpcchunkedgradient(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf, double *e, bool naturalerrorfunc) {
#ifndef ALGLIB_INTERCEPTS_SSE2
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hpcchunkedgradient(weights, structinfo, columnmeans, columnsigmas, xy, cstart, csize, batch4buf, hpcbuf, e, naturalerrorfunc);
#endif
}

// Fast kernel for chunked processing.
//
bool hpcchunkedprocess(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf) {
#ifndef ALGLIB_INTERCEPTS_SSE2
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hpcchunkedprocess(weights, structinfo, columnmeans, columnsigmas, xy, cstart, csize, batch4buf, hpcbuf);
#endif
}

void mlpbuffers_init(void *_p, bool make_automatic) {
   mlpbuffers *p = (mlpbuffers *)_p;
   ae_vector_init(&p->batch4buf, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->hpcbuf, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->xy, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->xy2, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xyrow, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->y, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->desiredy, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->g, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmp0, 0, DT_REAL, make_automatic);
}

void mlpbuffers_copy(void *_dst, const void *_src, bool make_automatic) {
   mlpbuffers *dst = (mlpbuffers *)_dst;
   const mlpbuffers *src = (const mlpbuffers *)_src;
   dst->chunksize = src->chunksize;
   dst->ntotal = src->ntotal;
   dst->nin = src->nin;
   dst->nout = src->nout;
   dst->wcount = src->wcount;
   ae_vector_copy(&dst->batch4buf, &src->batch4buf, make_automatic);
   ae_vector_copy(&dst->hpcbuf, &src->hpcbuf, make_automatic);
   ae_matrix_copy(&dst->xy, &src->xy, make_automatic);
   ae_matrix_copy(&dst->xy2, &src->xy2, make_automatic);
   ae_vector_copy(&dst->xyrow, &src->xyrow, make_automatic);
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   ae_vector_copy(&dst->y, &src->y, make_automatic);
   ae_vector_copy(&dst->desiredy, &src->desiredy, make_automatic);
   dst->e = src->e;
   ae_vector_copy(&dst->g, &src->g, make_automatic);
   ae_vector_copy(&dst->tmp0, &src->tmp0, make_automatic);
}

void mlpbuffers_free(void *_p, bool make_automatic) {
   mlpbuffers *p = (mlpbuffers *)_p;
   ae_vector_free(&p->batch4buf, make_automatic);
   ae_vector_free(&p->hpcbuf, make_automatic);
   ae_matrix_free(&p->xy, make_automatic);
   ae_matrix_free(&p->xy2, make_automatic);
   ae_vector_free(&p->xyrow, make_automatic);
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->y, make_automatic);
   ae_vector_free(&p->desiredy, make_automatic);
   ae_vector_free(&p->g, make_automatic);
   ae_vector_free(&p->tmp0, make_automatic);
}
} // end of namespace alglib_impl
