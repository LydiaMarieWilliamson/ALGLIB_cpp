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
#include "AlgLibInternal.h"

// === APSERV Package ===
namespace alglib_impl {
// The function always returns False.
// It may be used sometimes to prevent spurious warnings.
// ALGLIB: Copyright 17.09.2012 by Sergey Bochkanov
bool alwaysfalse(ae_state *_state) {
   bool result;
   result = false;
   return result;
}

// The function "touches" integer - it is used  to  avoid  compiler  messages
// about unused variables (in rare cases when we do NOT want to remove  these
// variables).
// ALGLIB: Copyright 17.09.2012 by Sergey Bochkanov
void touchint(ae_int_t *a, ae_state *_state) {
}

// The function "touches" real   -  it is used  to  avoid  compiler  messages
// about unused variables (in rare cases when we do NOT want to remove  these
// variables).
// ALGLIB: Copyright 17.09.2012 by Sergey Bochkanov
void touchreal(double *a, ae_state *_state) {
}

// The function performs zero-coalescing on real value.
//
// NOTE: no check is performed for B != 0
// ALGLIB: Copyright 18.05.2015 by Sergey Bochkanov
double coalesce(double a, double b, ae_state *_state) {
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
ae_int_t coalescei(ae_int_t a, ae_int_t b, ae_state *_state) {
   ae_int_t result;
   result = a;
   if (a == 0) {
      result = b;
   }
   return result;
}

// The function convert integer value to real value.
// ALGLIB: Copyright 17.09.2012 by Sergey Bochkanov
double inttoreal(ae_int_t a, ae_state *_state) {
   double result;
   result = (double)a;
   return result;
}

// The function calculates binary logarithm.
//
// NOTE: it costs twice as much as Ln(x)
// ALGLIB: Copyright 17.09.2012 by Sergey Bochkanov
double logbase2(double x, ae_state *_state) {
   double result;
   result = ae_log(x, _state) / ae_log(2.0, _state);
   return result;
}

// This function compares two numbers for approximate equality, with tolerance
// to errors as large as tol.
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
bool approxequal(double a, double b, double tol, ae_state *_state) {
   bool result;
   result = ae_fabs(a - b, _state) <= tol;
   return result;
}

// This function compares two numbers for approximate equality, with tolerance
// to errors as large as max(|a|,|b|)*tol.
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
bool approxequalrel(double a, double b, double tol, ae_state *_state) {
   bool result;
   result = ae_fabs(a - b, _state) <= maxreal(ae_fabs(a, _state), ae_fabs(b, _state), _state) * tol;
   return result;
}

// This  function  generates  1-dimensional  general  interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N=1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1d(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t i;
   double h;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "TaskGenInterpolationEqdist1D: N<1!", _state);
   ae_vector_set_length(x, n, _state);
   ae_vector_set_length(y, n, _state);
   if (n > 1) {
      x->xR[0] = a;
      y->xR[0] = 2 * randomreal(_state) - 1;
      h = (b - a) / (n - 1);
      for (i = 1; i < n; i++) {
         if (i != n - 1) {
            x->xR[i] = a + (i + 0.2 * (2 * randomreal(_state) - 1)) * h;
         } else {
            x->xR[i] = b;
         }
         y->xR[i] = y->xR[i - 1] + (2 * randomreal(_state) - 1) * (x->xR[i] - x->xR[i - 1]);
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = 2 * randomreal(_state) - 1;
   }
}

// This function generates  1-dimensional equidistant interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N=1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1dequidist(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t i;
   double h;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "TaskGenInterpolationEqdist1D: N<1!", _state);
   ae_vector_set_length(x, n, _state);
   ae_vector_set_length(y, n, _state);
   if (n > 1) {
      x->xR[0] = a;
      y->xR[0] = 2 * randomreal(_state) - 1;
      h = (b - a) / (n - 1);
      for (i = 1; i < n; i++) {
         x->xR[i] = a + i * h;
         y->xR[i] = y->xR[i - 1] + (2 * randomreal(_state) - 1) * h;
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = 2 * randomreal(_state) - 1;
   }
}

// This function generates  1-dimensional Chebyshev-1 interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N=1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1dcheb1(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t i;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "TaskGenInterpolation1DCheb1: N<1!", _state);
   ae_vector_set_length(x, n, _state);
   ae_vector_set_length(y, n, _state);
   if (n > 1) {
      for (i = 0; i < n; i++) {
         x->xR[i] = 0.5 * (b + a) + 0.5 * (b - a) * ae_cos(pi * (2 * i + 1) / (2 * n), _state);
         if (i == 0) {
            y->xR[i] = 2 * randomreal(_state) - 1;
         } else {
            y->xR[i] = y->xR[i - 1] + (2 * randomreal(_state) - 1) * (x->xR[i] - x->xR[i - 1]);
         }
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = 2 * randomreal(_state) - 1;
   }
}

// This function generates  1-dimensional Chebyshev-2 interpolation task with
// moderate Lipshitz constant (close to 1.0)
//
// If N=1 then suborutine generates only one point at the middle of [A,B]
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void taskgenint1dcheb2(double a, double b, ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t i;
   SetVector(x);
   SetVector(y);
   ae_assert(n >= 1, "TaskGenInterpolation1DCheb2: N<1!", _state);
   ae_vector_set_length(x, n, _state);
   ae_vector_set_length(y, n, _state);
   if (n > 1) {
      for (i = 0; i < n; i++) {
         x->xR[i] = 0.5 * (b + a) + 0.5 * (b - a) * ae_cos(pi * i / (n - 1), _state);
         if (i == 0) {
            y->xR[i] = 2 * randomreal(_state) - 1;
         } else {
            y->xR[i] = y->xR[i - 1] + (2 * randomreal(_state) - 1) * (x->xR[i] - x->xR[i - 1]);
         }
      }
   } else {
      x->xR[0] = 0.5 * (a + b);
      y->xR[0] = 2 * randomreal(_state) - 1;
   }
}

// This function checks that all values from X[] are distinct. It does more
// than just usual floating point comparison:
// * first, it calculates max(X) and min(X)
// * second, it maps X[] from [min,max] to [1,2]
// * only at this stage actual comparison is done
//
// The meaning of such check is to ensure that all values are "distinct enough"
// and will not cause interpolation subroutine to fail.
//
// NOTE:
//     X[] must be sorted by ascending (subroutine ASSERT's it)
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
bool aredistinct(RVector *x, ae_int_t n, ae_state *_state) {
   double a;
   double b;
   ae_int_t i;
   bool nonsorted;
   bool result;
   ae_assert(n >= 1, "APSERVAreDistinct: internal error (N<1)", _state);
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
      a = minreal(a, x->xR[i], _state);
      b = maxreal(b, x->xR[i], _state);
      nonsorted = nonsorted || x->xR[i - 1] >= x->xR[i];
   }
   ae_assert(!nonsorted, "APSERVAreDistinct: internal error (not sorted)", _state);
   for (i = 1; i < n; i++) {
      if ((x->xR[i] - a) / (b - a) + 1 == (x->xR[i - 1] - a) / (b - a) + 1) {
         result = false;
         return result;
      }
   }
   result = true;
   return result;
}

// This function checks that two boolean values are the same (both  are  True
// or both are False).
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
bool aresameboolean(bool v1, bool v2, ae_state *_state) {
   bool result;
   result = (v1 && v2) || (!v1 && !v2);
   return result;
}

// Resizes X and fills by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void setlengthzero(RVector *x, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_assert(n >= 0, "SetLengthZero: N<0", _state);
   ae_vector_set_length(x, n, _state);
   for (i = 0; i < n; i++) {
      x->xR[i] = 0.0;
   }
}

// If Length(X)<N, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void bvectorsetlengthatleast(BVector *x, ae_int_t n, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
}

// If Length(X)<N, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void ivectorsetlengthatleast(ZVector *x, ae_int_t n, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
}

// If Length(X)<N, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rvectorsetlengthatleast(RVector *x, ae_int_t n, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
}

// If Cols(X)<N or Rows(X)<M, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void bmatrixsetlengthatleast(BMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   if (m > 0 && n > 0) {
      if (x->rows < m || x->cols < n) {
         ae_matrix_set_length(x, m, n, _state);
      }
   }
}

// If Cols(X)<N or Rows(X)<M, resizes X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rmatrixsetlengthatleast(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   if (m > 0 && n > 0) {
      if (x->rows < m || x->cols < n) {
         ae_matrix_set_length(x, m, n, _state);
      }
   }
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) new size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void bvectorgrowto(BVector *x, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_BOOL, _state);
// Enough place
   if (x->cnt >= n) {
      ae_frame_leave(_state);
      return;
   }
// Choose new size
   n = maxint(n, iround(1.8 * x->cnt + 1, _state), _state);
// Grow
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n, _state);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xB[i] = oldx.xB[i];
      } else {
         x->xB[i] = false;
      }
   }
   ae_frame_leave(_state);
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) new size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void ivectorgrowto(ZVector *x, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_INT, _state);
// Enough place
   if (x->cnt >= n) {
      ae_frame_leave(_state);
      return;
   }
// Choose new size
   n = maxint(n, iround(1.8 * x->cnt + 1, _state), _state);
// Grow
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n, _state);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xZ[i] = oldx.xZ[i];
      } else {
         x->xZ[i] = 0;
      }
   }
   ae_frame_leave(_state);
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) new size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rvectorgrowto(RVector *x, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_REAL, _state);
// Enough place
   if (x->cnt >= n) {
      ae_frame_leave(_state);
      return;
   }
// Choose new size
   n = maxint(n, iround(1.8 * x->cnt + 1, _state), _state);
// Grow
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n, _state);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xR[i] = oldx.xR[i];
      } else {
         x->xR[i] = 0.0;
      }
   }
   ae_frame_leave(_state);
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
void rmatrixgrowrowsto(RMatrix *a, ae_int_t n, ae_int_t mincols, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t n2;
   ae_int_t m;
   ae_frame_make(_state, &_frame_block);
   NewMatrix(olda, 0, 0, DT_REAL, _state);
// Enough place?
   if (a->rows >= n && a->cols >= mincols) {
      ae_frame_leave(_state);
      return;
   }
// Sizes and metrics
   if (a->rows < n) {
      n = maxint(n, iround(1.8 * a->rows + 1, _state), _state);
   }
   n2 = minint(a->rows, n, _state);
   m = a->cols;
// Grow
   ae_swap_matrices(a, &olda);
   ae_matrix_set_length(a, n, maxint(m, mincols, _state), _state);
   for (i = 0; i < n2; i++) {
      for (j = 0; j < m; j++) {
         a->xyR[i][j] = olda.xyR[i][j];
      }
   }
   ae_frame_leave(_state);
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
void rmatrixgrowcolsto(RMatrix *a, ae_int_t n, ae_int_t minrows, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t n2;
   ae_int_t m;
   ae_frame_make(_state, &_frame_block);
   NewMatrix(olda, 0, 0, DT_REAL, _state);
// Enough place?
   if (a->cols >= n && a->rows >= minrows) {
      ae_frame_leave(_state);
      return;
   }
// Sizes and metrics
   if (a->cols < n) {
      n = maxint(n, iround(1.8 * a->cols + 1, _state), _state);
   }
   n2 = minint(a->cols, n, _state);
   m = a->rows;
// Grow
   ae_swap_matrices(a, &olda);
   ae_matrix_set_length(a, maxint(m, minrows, _state), n, _state);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n2; j++) {
         a->xyR[i][j] = olda.xyR[i][j];
      }
   }
   ae_frame_leave(_state);
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void ivectorresize(ZVector *x, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_INT, _state);
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n, _state);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xZ[i] = oldx.xZ[i];
      } else {
         x->xZ[i] = 0;
      }
   }
   ae_frame_leave(_state);
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rvectorresize(RVector *x, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_REAL, _state);
   n2 = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n, _state);
   for (i = 0; i < n; i++) {
      if (i < n2) {
         x->xR[i] = oldx.xR[i];
      } else {
         x->xR[i] = 0.0;
      }
   }
   ae_frame_leave(_state);
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void imatrixresize(ZMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t m2;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewMatrix(oldx, 0, 0, DT_INT, _state);
   m2 = x->rows;
   n2 = x->cols;
   ae_swap_matrices(x, &oldx);
   ae_matrix_set_length(x, m, n, _state);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (i < m2 && j < n2) {
            x->xyZ[i][j] = oldx.xyZ[i][j];
         } else {
            x->xyZ[i][j] = 0;
         }
      }
   }
   ae_frame_leave(_state);
}

// Resizes X and:
// * preserves old contents of X
// * fills new elements by zeros
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rmatrixresize(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t m2;
   ae_int_t n2;
   ae_frame_make(_state, &_frame_block);
   NewMatrix(oldx, 0, 0, DT_REAL, _state);
   m2 = x->rows;
   n2 = x->cols;
   ae_swap_matrices(x, &oldx);
   ae_matrix_set_length(x, m, n, _state);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (i < m2 && j < n2) {
            x->xyR[i][j] = oldx.xyR[i][j];
         } else {
            x->xyR[i][j] = 0.0;
         }
      }
   }
   ae_frame_leave(_state);
}

// appends element to X
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void ivectorappend(ZVector *x, ae_int_t v, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t n;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_INT, _state);
   n = x->cnt;
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, n + 1, _state);
   for (i = 0; i < n; i++) {
      x->xZ[i] = oldx.xZ[i];
   }
   x->xZ[n] = v;
   ae_frame_leave(_state);
}

// This function checks that length(X) is at least N and first N values  from
// X[] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool isfinitevector(RVector *x, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   double v;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteVector: internal error (N<0)", _state);
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
   result = ae_isfinite(v, _state);
   return result;
}

// This function checks that first N values from X[] are finite
// ALGLIB: Copyright 18.06.2010 by Sergey Bochkanov
bool isfinitecvector(CVector *z, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteCVector: internal error (N<0)", _state);
   for (i = 0; i < n; i++) {
      if (!ae_isfinite(z->xC[i].x, _state) || !ae_isfinite(z->xC[i].y, _state)) {
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
bool apservisfinitematrix(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteMatrix: internal error (N<0)", _state);
   ae_assert(m >= 0, "APSERVIsFiniteMatrix: internal error (M<0)", _state);
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
         if (!ae_isfinite(x->xyR[i][j], _state)) {
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
bool apservisfinitecmatrix(CMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteCMatrix: internal error (N<0)", _state);
   ae_assert(m >= 0, "APSERVIsFiniteCMatrix: internal error (M<0)", _state);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (!ae_isfinite(x->xyC[i][j].x, _state) || !ae_isfinite(x->xyC[i][j].y, _state)) {
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
bool isfinitertrmatrix(RMatrix *x, ae_int_t n, bool isupper, ae_state *_state) {
   ae_int_t i;
   ae_int_t j1;
   ae_int_t j2;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteRTRMatrix: internal error (N<0)", _state);
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
         if (!ae_isfinite(x->xyR[i][j], _state)) {
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
bool apservisfinitectrmatrix(CMatrix *x, ae_int_t n, bool isupper, ae_state *_state) {
   ae_int_t i;
   ae_int_t j1;
   ae_int_t j2;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteCTRMatrix: internal error (N<0)", _state);
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i;
      }
      for (j = j1; j <= j2; j++) {
         if (!ae_isfinite(x->xyC[i][j].x, _state) || !ae_isfinite(x->xyC[i][j].y, _state)) {
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
bool apservisfiniteornanmatrix(RMatrix *x, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_assert(n >= 0, "APSERVIsFiniteOrNaNMatrix: internal error (N<0)", _state);
   ae_assert(m >= 0, "APSERVIsFiniteOrNaNMatrix: internal error (M<0)", _state);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (!(ae_isfinite(x->xyR[i][j], _state) || ae_isnan(x->xyR[i][j], _state))) {
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
double safepythag2(double x, double y, ae_state *_state) {
   double w;
   double xabs;
   double yabs;
   double z;
   double result;
   xabs = ae_fabs(x, _state);
   yabs = ae_fabs(y, _state);
   w = maxreal(xabs, yabs, _state);
   z = minreal(xabs, yabs, _state);
   if (z == 0.0) {
      result = w;
   } else {
      result = w * ae_sqrt(1 + sqr(z / w, _state), _state);
   }
   return result;
}

// Safe sqrt(x^2+y^2)
// ALGLIB: Copyright by Sergey Bochkanov
double safepythag3(double x, double y, double z, ae_state *_state) {
   double w;
   double result;
   w = maxreal(ae_fabs(x, _state), maxreal(ae_fabs(y, _state), ae_fabs(z, _state), _state), _state);
   if (w == 0.0) {
      result = 0.0;
      return result;
   }
   x = x / w;
   y = y / w;
   z = z / w;
   result = w * ae_sqrt(sqr(x, _state) + sqr(y, _state) + sqr(z, _state), _state);
   return result;
}

// Safe division.
//
// This function attempts to calculate R=X/Y without overflow.
//
// It returns:
// * +1, if abs(X/Y) >= MaxRealNumber or undefined - overflow-like situation
//       (no overlfow is generated, R is either NAN, PosINF, NegINF)
// *  0, if MinRealNumber<abs(X/Y)<MaxRealNumber or X=0, Y != 0
//       (R contains result, may be zero)
// * -1, if 0<abs(X/Y)<MinRealNumber - underflow-like situation
//       (R contains zero; it corresponds to underflow)
//
// No overflow is generated in any case.
// ALGLIB: Copyright by Sergey Bochkanov
ae_int_t saferdiv(double x, double y, double *r, ae_state *_state) {
   ae_int_t result;
   *r = 0;
// Two special cases:
// * Y=0
// * X=0 and Y != 0
   if (y == 0.0) {
      result = 1;
      if (x == 0.0) {
         *r = _state->v_nan;
      }
      if (x > 0.0) {
         *r = _state->v_posinf;
      }
      if (x < 0.0) {
         *r = _state->v_neginf;
      }
      return result;
   }
   if (x == 0.0) {
      *r = 0.0;
      result = 0;
      return result;
   }
// make Y>0
   if (y < 0.0) {
      x = -x;
      y = -y;
   }
//
   if (y >= 1.0) {
      *r = x / y;
      if (ae_fabs(*r, _state) <= minrealnumber) {
         result = -1;
         *r = 0.0;
      } else {
         result = 0;
      }
   } else {
      if (ae_fabs(x, _state) >= maxrealnumber * y) {
         if (x > 0.0) {
            *r = _state->v_posinf;
         } else {
            *r = _state->v_neginf;
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
double safeminposrv(double x, double y, double v, ae_state *_state) {
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
   // Y<1, we can safely multiply by Y
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
// It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
// such that X = T + K*(B-A).
//
// NOTES:
// * K is represented as real value, although actually it is integer
// * T is guaranteed to be in [A,B]
// * T replaces X
// ALGLIB: Copyright by Sergey Bochkanov
void apperiodicmap(double *x, double a, double b, double *k, ae_state *_state) {
   *k = 0;
   ae_assert(a < b, "APPeriodicMap: internal error!", _state);
   *k = floor((*x - a) / (b - a));
   *x = *x - *k * (b - a);
   while (*x < a) {
      *x = *x + (b - a);
      *k = *k - 1;
   }
   while (*x > b) {
      *x = *x - (b - a);
      *k = *k + 1;
   }
   *x = maxreal(*x, a, _state);
   *x = minreal(*x, b, _state);
}

// Returns random normal number using low-quality system-provided generator
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
double randomnormal(ae_state *_state) {
   double u;
   double v;
   double s;
   double result;
   while (true) {
      u = 2 * randomreal(_state) - 1;
      v = 2 * randomreal(_state) - 1;
      s = sqr(u, _state) + sqr(v, _state);
      if (s > 0.0 && s < 1.0) {
      // two Sqrt's instead of one to
      // avoid overflow when S is too small
         s = ae_sqrt(-2 * ae_log(s, _state), _state) / ae_sqrt(s, _state);
         result = u * s;
         break;
      }
   }
   return result;
}

// Generates random unit vector using low-quality system-provided generator.
// Reallocates array if its size is too short.
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void randomunit(ae_int_t n, RVector *x, ae_state *_state) {
   ae_int_t i;
   double v;
   double vv;
   ae_assert(n > 0, "RandomUnit: N <= 0", _state);
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
   do {
      v = 0.0;
      for (i = 0; i < n; i++) {
         vv = randomnormal(_state);
         x->xR[i] = vv;
         v = v + vv * vv;
      }
   } while (v <= 0.0);
   v = 1 / ae_sqrt(v, _state);
   for (i = 0; i < n; i++) {
      x->xR[i] = x->xR[i] * v;
   }
}

// This function is used to swap two integer values
void swapi(ae_int_t *v0, ae_int_t *v1, ae_state *_state) {
   ae_int_t v;
   v = *v0;
   *v0 = *v1;
   *v1 = v;
}

// This function is used to swap two real values
void swapr(double *v0, double *v1, ae_state *_state) {
   double v;
   v = *v0;
   *v0 = *v1;
   *v1 = v;
}

// This function is used to swap two cols of the matrix; if NRows<0, automatically
// determined from the matrix size.
void swapcols(RMatrix *a, ae_int_t j0, ae_int_t j1, ae_int_t nrows, ae_state *_state) {
   ae_int_t i;
   double v;
   if (j0 == j1) {
      return;
   }
   if (nrows < 0) {
      nrows = a->rows;
   }
   for (i = 0; i < nrows; i++) {
      v = a->xyR[i][j0];
      a->xyR[i][j0] = a->xyR[i][j1];
      a->xyR[i][j1] = v;
   }
}

// This function is used to swap two rows of the matrix; if NCols<0, automatically
// determined from the matrix size.
void swaprows(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t ncols, ae_state *_state) {
   ae_int_t j;
   double v;
   if (i0 == i1) {
      return;
   }
   if (ncols < 0) {
      ncols = a->cols;
   }
   for (j = 0; j < ncols; j++) {
      v = a->xyR[i0][j];
      a->xyR[i0][j] = a->xyR[i1][j];
      a->xyR[i1][j] = v;
   }
}

// This function is used to swap two "entries" in 1-dimensional array composed
// from D-element entries
void swapentries(RVector *a, ae_int_t i0, ae_int_t i1, ae_int_t entrywidth, ae_state *_state) {
   ae_int_t offs0;
   ae_int_t offs1;
   ae_int_t j;
   double v;
   if (i0 == i1) {
      return;
   }
   offs0 = i0 * entrywidth;
   offs1 = i1 * entrywidth;
   for (j = 0; j < entrywidth; j++) {
      v = a->xR[offs0 + j];
      a->xR[offs0 + j] = a->xR[offs1 + j];
      a->xR[offs1 + j] = v;
   }
}

// This function is used to swap two elements of the vector
void swapelementsi(ZVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state) {
   ae_int_t v;
   if (i0 == i1) {
      return;
   }
   v = a->xZ[i0];
   a->xZ[i0] = a->xZ[i1];
   a->xZ[i1] = v;
}

// This function is used to swap two elements of the vector
void swapelements(RVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state) {
   double v;
   if (i0 == i1) {
      return;
   }
   v = a->xR[i0];
   a->xR[i0] = a->xR[i1];
   a->xR[i1] = v;
}

// This function is used to increment value of integer variable
void inc(ae_int_t *v, ae_state *_state) {
   *v = *v + 1;
}

// This function is used to decrement value of integer variable
void dec(ae_int_t *v, ae_state *_state) {
   *v = *v - 1;
}

// This function is used to increment value of integer variable; name of  the
// function suggests that increment is done in multithreaded setting  in  the
// thread-unsafe manner (optional progress reports which do not need guaranteed
// correctness)
void threadunsafeinc(ae_int_t *v, ae_state *_state) {
   *v = *v + 1;
}

// This function is used to increment value of integer variable; name of  the
// function suggests that increment is done in multithreaded setting  in  the
// thread-unsafe manner (optional progress reports which do not need guaranteed
// correctness)
void threadunsafeincby(ae_int_t *v, ae_int_t k, ae_state *_state) {
   *v = *v + k;
}

// This function performs two operations:
// 1. decrements value of integer variable, if it is positive
// 2. explicitly sets variable to zero if it is non-positive
// It is used by some algorithms to decrease value of internal counters.
void countdown(ae_int_t *v, ae_state *_state) {
   if (*v > 0) {
      *v = *v - 1;
   } else {
      *v = 0;
   }
}

// This function returns +1 or -1 depending on sign of X.
// x=0 results in +1 being returned.
double possign(double x, ae_state *_state) {
   double result;
   if (x >= 0.0) {
      result = 1.0;
   } else {
      result = -1.0;
   }
   return result;
}

// This function returns product of two real numbers. It is convenient when
// you have to perform typecast-and-product of two INTEGERS.
double rmul2(double v0, double v1, ae_state *_state) {
   double result;
   result = v0 * v1;
   return result;
}

// This function returns product of three real numbers. It is convenient when
// you have to perform typecast-and-product of two INTEGERS.
double rmul3(double v0, double v1, double v2, ae_state *_state) {
   double result;
   result = v0 * v1 * v2;
   return result;
}

// This function returns (A div B) rounded up; it expects that A>0, B>0, but
// does not check it.
ae_int_t idivup(ae_int_t a, ae_int_t b, ae_state *_state) {
   ae_int_t result;
   result = a / b;
   if (a % b > 0) {
      result = result + 1;
   }
   return result;
}

// This function returns min(i0,i1)
ae_int_t imin2(ae_int_t i0, ae_int_t i1, ae_state *_state) {
   ae_int_t result;
   result = i0;
   if (i1 < result) {
      result = i1;
   }
   return result;
}

// This function returns min(i0,i1,i2)
ae_int_t imin3(ae_int_t i0, ae_int_t i1, ae_int_t i2, ae_state *_state) {
   ae_int_t result;
   result = i0;
   if (i1 < result) {
      result = i1;
   }
   if (i2 < result) {
      result = i2;
   }
   return result;
}

// This function returns max(i0,i1)
ae_int_t imax2(ae_int_t i0, ae_int_t i1, ae_state *_state) {
   ae_int_t result;
   result = i0;
   if (i1 > result) {
      result = i1;
   }
   return result;
}

// This function returns max(i0,i1,i2)
ae_int_t imax3(ae_int_t i0, ae_int_t i1, ae_int_t i2, ae_state *_state) {
   ae_int_t result;
   result = i0;
   if (i1 > result) {
      result = i1;
   }
   if (i2 > result) {
      result = i2;
   }
   return result;
}

// This function is used to return maximum of three real values
double maxreal3(double v0, double v1, double v2, ae_state *_state) {
   double result;
   result = v0;
   if (result < v1) {
      result = v1;
   }
   if (result < v2) {
      result = v2;
   }
   return result;
}

// This function returns max(r0,r1,r2)
double rmax3(double r0, double r1, double r2, ae_state *_state) {
   double result;
   result = r0;
   if (r1 > result) {
      result = r1;
   }
   if (r2 > result) {
      result = r2;
   }
   return result;
}

// This function returns max(|r0|,|r1|,|r2|)
double rmaxabs3(double r0, double r1, double r2, ae_state *_state) {
   double result;
   r0 = ae_fabs(r0, _state);
   r1 = ae_fabs(r1, _state);
   r2 = ae_fabs(r2, _state);
   result = r0;
   if (r1 > result) {
      result = r1;
   }
   if (r2 > result) {
      result = r2;
   }
   return result;
}

// 'bounds' value: maps X to [B1,B2]
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
double boundval(double x, double b1, double b2, ae_state *_state) {
   double result;
   if (x <= b1) {
      result = b1;
      return result;
   }
   if (x >= b2) {
      result = b2;
      return result;
   }
   result = x;
   return result;
}

// 'bounds' value: maps X to [B1,B2]
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
ae_int_t iboundval(ae_int_t x, ae_int_t b1, ae_int_t b2, ae_state *_state) {
   ae_int_t result;
   if (x <= b1) {
      result = b1;
      return result;
   }
   if (x >= b2) {
      result = b2;
      return result;
   }
   result = x;
   return result;
}

// 'bounds' value: maps X to [B1,B2]
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
double rboundval(double x, double b1, double b2, ae_state *_state) {
   double result;
   if (x <= b1) {
      result = b1;
      return result;
   }
   if (x >= b2) {
      result = b2;
      return result;
   }
   result = x;
   return result;
}

// Returns number of non-zeros
ae_int_t countnz1(RVector *v, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t result;
   result = 0;
   for (i = 0; i < n; i++) {
      if (!(v->xR[i] == 0)) {
         result = result + 1;
      }
   }
   return result;
}

// Returns number of non-zeros
ae_int_t countnz2(RMatrix *v, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t result;
   result = 0;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         if (!(v->xyR[i][j] == 0)) {
            result = result + 1;
         }
      }
   }
   return result;
}

// Allocation of serializer: complex value
void alloccomplex(ae_serializer *s, complex v, ae_state *_state) {
   ae_serializer_alloc_entry(s);
   ae_serializer_alloc_entry(s);
}

// Serialization: complex value
void serializecomplex(ae_serializer *s, complex v, ae_state *_state) {
   ae_serializer_serialize_double(s, v.x, _state);
   ae_serializer_serialize_double(s, v.y, _state);
}

// Unserialization: complex value
complex unserializecomplex(ae_serializer *s, ae_state *_state) {
   complex result;
   ae_serializer_unserialize_double(s, &result.x, _state);
   ae_serializer_unserialize_double(s, &result.y, _state);
   return result;
}

// Allocation of serializer: real array
void allocrealarray(ae_serializer *s, RVector *v, ae_int_t n, ae_state *_state) {
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
void serializerealarray(ae_serializer *s, RVector *v, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_serialize_int(s, n, _state);
   for (i = 0; i < n; i++) {
      ae_serializer_serialize_double(s, v->xR[i], _state);
   }
}

// Unserialization: complex value
void unserializerealarray(ae_serializer *s, RVector *v, ae_state *_state) {
   ae_int_t n;
   ae_int_t i;
   double t;
   SetVector(v);
   ae_serializer_unserialize_int(s, &n, _state);
   if (n == 0) {
      return;
   }
   ae_vector_set_length(v, n, _state);
   for (i = 0; i < n; i++) {
      ae_serializer_unserialize_double(s, &t, _state);
      v->xR[i] = t;
   }
}

// Allocation of serializer: Integer array
void allocintegerarray(ae_serializer *s, ZVector *v, ae_int_t n, ae_state *_state) {
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
void serializeintegerarray(ae_serializer *s, ZVector *v, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   if (n < 0) {
      n = v->cnt;
   }
   ae_serializer_serialize_int(s, n, _state);
   for (i = 0; i < n; i++) {
      ae_serializer_serialize_int(s, v->xZ[i], _state);
   }
}

// Unserialization: complex value
void unserializeintegerarray(ae_serializer *s, ZVector *v, ae_state *_state) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t t;
   SetVector(v);
   ae_serializer_unserialize_int(s, &n, _state);
   if (n == 0) {
      return;
   }
   ae_vector_set_length(v, n, _state);
   for (i = 0; i < n; i++) {
      ae_serializer_unserialize_int(s, &t, _state);
      v->xZ[i] = t;
   }
}

// Allocation of serializer: real matrix
void allocrealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1, ae_state *_state) {
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
void serializerealmatrix(ae_serializer *s, RMatrix *v, ae_int_t n0, ae_int_t n1, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   if (n0 < 0) {
      n0 = v->rows;
   }
   if (n1 < 0) {
      n1 = v->cols;
   }
   ae_serializer_serialize_int(s, n0, _state);
   ae_serializer_serialize_int(s, n1, _state);
   for (i = 0; i < n0; i++) {
      for (j = 0; j < n1; j++) {
         ae_serializer_serialize_double(s, v->xyR[i][j], _state);
      }
   }
}

// Unserialization: complex value
void unserializerealmatrix(ae_serializer *s, RMatrix *v, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t n0;
   ae_int_t n1;
   double t;
   SetMatrix(v);
   ae_serializer_unserialize_int(s, &n0, _state);
   ae_serializer_unserialize_int(s, &n1, _state);
   if (n0 == 0 || n1 == 0) {
      return;
   }
   ae_matrix_set_length(v, n0, n1, _state);
   for (i = 0; i < n0; i++) {
      for (j = 0; j < n1; j++) {
         ae_serializer_unserialize_double(s, &t, _state);
         v->xyR[i][j] = t;
      }
   }
}

// Copy boolean array
void copybooleanarray(BVector *src, BVector *dst, ae_state *_state) {
   ae_int_t i;
   SetVector(dst);
   if (src->cnt > 0) {
      ae_vector_set_length(dst, src->cnt, _state);
      for (i = 0; i < src->cnt; i++) {
         dst->xB[i] = src->xB[i];
      }
   }
}

// Copy integer array
void copyintegerarray(ZVector *src, ZVector *dst, ae_state *_state) {
   ae_int_t i;
   SetVector(dst);
   if (src->cnt > 0) {
      ae_vector_set_length(dst, src->cnt, _state);
      for (i = 0; i < src->cnt; i++) {
         dst->xZ[i] = src->xZ[i];
      }
   }
}

// Copy real array
void copyrealarray(RVector *src, RVector *dst, ae_state *_state) {
   ae_int_t i;
   SetVector(dst);
   if (src->cnt > 0) {
      ae_vector_set_length(dst, src->cnt, _state);
      for (i = 0; i < src->cnt; i++) {
         dst->xR[i] = src->xR[i];
      }
   }
}

// Copy real matrix
void copyrealmatrix(RMatrix *src, RMatrix *dst, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(dst);
   if (src->rows > 0 && src->cols > 0) {
      ae_matrix_set_length(dst, src->rows, src->cols, _state);
      for (i = 0; i < src->rows; i++) {
         for (j = 0; j < src->cols; j++) {
            dst->xyR[i][j] = src->xyR[i][j];
         }
      }
   }
}

// Clears integer array
void unsetintegerarray(ZVector *a, ae_state *_state) {
   SetVector(a);
}

// Clears real array
void unsetrealarray(RVector *a, ae_state *_state) {
   SetVector(a);
}

// Clears real matrix
void unsetrealmatrix(RMatrix *a, ae_state *_state) {
   SetMatrix(a);
}

// This function is used to calculate number of chunks (including partial,
// non-complete chunks) in some set. It expects that ChunkSize >= 1, TaskSize >= 0.
// Assertion is thrown otherwise.
//
// Function result is equivalent to Ceil(TaskSize/ChunkSize), but with guarantees
// that rounding errors won't ruin results.
// ALGLIB: Copyright 21.01.2015 by Sergey Bochkanov
ae_int_t chunkscount(ae_int_t tasksize, ae_int_t chunksize, ae_state *_state) {
   ae_int_t result;
   ae_assert(tasksize >= 0, "ChunksCount: TaskSize<0", _state);
   ae_assert(chunksize >= 1, "ChunksCount: ChunkSize<1", _state);
   result = tasksize / chunksize;
   if (tasksize % chunksize != 0) {
      result = result + 1;
   }
   return result;
}

// This function is used in parallel functions for recurrent division of large
// task into two smaller tasks.
//
// It has following properties:
// * it works only for TaskSize >= 2 and TaskSize>TileSize (assertion is thrown otherwise)
// * Task0+Task1=TaskSize, Task0>0, Task1>0
// * Task0 and Task1 are close to each other
// * Task0 >= Task1
// * Task0 is always divisible by TileSize
// ALGLIB: Copyright 07.04.2013 by Sergey Bochkanov
void tiledsplit(ae_int_t tasksize, ae_int_t tilesize, ae_int_t *task0, ae_int_t *task1, ae_state *_state) {
   ae_int_t cc;
   *task0 = 0;
   *task1 = 0;
   ae_assert(tasksize >= 2, "TiledSplit: TaskSize<2", _state);
   ae_assert(tasksize > tilesize, "TiledSplit: TaskSize <= TileSize", _state);
   cc = chunkscount(tasksize, tilesize, _state);
   ae_assert(cc >= 2, "TiledSplit: integrity check failed", _state);
   *task0 = idivup(cc, 2, _state) * tilesize;
   *task1 = tasksize - (*task0);
   ae_assert(*task0 >= 1, "TiledSplit: internal error", _state);
   ae_assert(*task1 >= 1, "TiledSplit: internal error", _state);
   ae_assert(*task0 % tilesize == 0, "TiledSplit: internal error", _state);
   ae_assert(*task0 >= (*task1), "TiledSplit: internal error", _state);
}

// --- OBSOLETE FUNCTION, USE TILED SPLIT INSTEAD ---
//
// This function is used in parallel functions for recurrent division of large
// task into two smaller tasks.
//
// It has following properties:
// * it works only for TaskSize >= 2 and ChunkSize >= 2
//   (assertion is thrown otherwise)
// * Task0+Task1=TaskSize, Task0>0, Task1>0
// * Task0 and Task1 are close to each other
// * in case TaskSize>ChunkSize, Task0 is always divisible by ChunkSize
// ALGLIB: Copyright 07.04.2013 by Sergey Bochkanov
void splitlength(ae_int_t tasksize, ae_int_t chunksize, ae_int_t *task0, ae_int_t *task1, ae_state *_state) {
   *task0 = 0;
   *task1 = 0;
   ae_assert(chunksize >= 2, "SplitLength: ChunkSize<2", _state);
   ae_assert(tasksize >= 2, "SplitLength: TaskSize<2", _state);
   *task0 = tasksize / 2;
   if (*task0 > chunksize && *task0 % chunksize != 0) {
      *task0 = *task0 - *task0 % chunksize;
   }
   *task1 = tasksize - (*task0);
   ae_assert(*task0 >= 1, "SplitLength: internal error", _state);
   ae_assert(*task1 >= 1, "SplitLength: internal error", _state);
}

// This function is used in parallel functions for recurrent division of large
// task into two smaller tasks.
//
// It has following properties:
// * it works only for TaskSize >= 2 (assertion is thrown otherwise)
// * for TaskSize=2, it returns Task0=1, Task1=1
// * in case TaskSize is odd,  Task0=TaskSize-1, Task1=1
// * in case TaskSize is even, Task0 and Task1 are approximately TaskSize/2
//   and both Task0 and Task1 are even, Task0 >= Task1
// ALGLIB: Copyright 07.04.2013 by Sergey Bochkanov
void splitlengtheven(ae_int_t tasksize, ae_int_t *task0, ae_int_t *task1, ae_state *_state) {
   *task0 = 0;
   *task1 = 0;
   ae_assert(tasksize >= 2, "SplitLengthEven: TaskSize<2", _state);
   if (tasksize == 2) {
      *task0 = 1;
      *task1 = 1;
      return;
   }
   if (tasksize % 2 == 0) {
   // Even division
      *task0 = tasksize / 2;
      *task1 = tasksize / 2;
      if (*task0 % 2 != 0) {
         *task0 = *task0 + 1;
         *task1 = *task1 - 1;
      }
   } else {
   // Odd task size, split trailing odd part from it.
      *task0 = tasksize - 1;
      *task1 = 1;
   }
   ae_assert(*task0 >= 1, "SplitLengthEven: internal error", _state);
   ae_assert(*task1 >= 1, "SplitLengthEven: internal error", _state);
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
ae_int_t recsearch(ZVector *a, ae_int_t nrec, ae_int_t nheader, ae_int_t i0, ae_int_t i1, ZVector *b, ae_state *_state) {
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
double sparselevel2density(ae_state *_state) {
   double result;
   result = 0.1;
   return result;
}

// Returns A-tile size for a matrix.
//
// A-tiles are smallest tiles (32x32), suitable for processing by ALGLIB  own
// implementation of Level 3 linear algebra.
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
ae_int_t matrixtilesizea(ae_state *_state) {
   ae_int_t result;
   result = 32;
   return result;
}

// Returns B-tile size for a matrix.
//
// B-tiles are larger  tiles (64x64), suitable for parallel execution or for
// processing by vendor's implementation of Level 3 linear algebra.
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
ae_int_t matrixtilesizeb(ae_state *_state) {
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
double smpactivationlevel(ae_state *_state) {
   double nn;
   double result;
   nn = (double)(2 * matrixtilesizeb(_state));
   result = maxreal(0.95 * 2 * nn * nn * nn, 1.0E7, _state);
   return result;
}

// This function returns minimum cost of task which is feasible for
// spawn (given that multithreading is active).
//
// It returns real number in order to avoid overflow problems.
// ALGLIB: Copyright 10.01.2018 by Sergey Bochkanov
double spawnlevel(ae_state *_state) {
   double nn;
   double result;
   nn = (double)(2 * matrixtilesizea(_state));
   result = 0.95 * 2 * nn * nn * nn;
   return result;
}

// Outputs vector A[I0,I1-1] to trace log using either:
// a)  6-digit exponential format (no trace flags is set)
// b) 15-ditit exponential format ('PREC.E15' trace flag is set)
// b)  6-ditit fixed-point format ('PREC.F6' trace flag is set)
//
// This function checks trace flags every time it is called.
void tracevectorautoprec(RVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state) {
   ae_int_t i;
   ae_int_t prectouse;
// Determine precision to use
   prectouse = 0;
   if (ae_is_trace_enabled("PREC.E15")) {
      prectouse = 1;
   }
   if (ae_is_trace_enabled("PREC.F6")) {
      prectouse = 2;
   }
// Output
   ae_trace("[ ");
   for (i = i0; i < i1; i++) {
      if (prectouse == 0) {
         ae_trace("%14.6e", a->xR[i]);
      }
      if (prectouse == 1) {
         ae_trace("%23.15e", a->xR[i]);
      }
      if (prectouse == 2) {
         ae_trace("%13.6f", a->xR[i]);
      }
      if (i < i1 - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

// Outputs row A[I,J0..J1-1] to trace log using either:
// a)  6-digit exponential format (no trace flags is set)
// b) 15-ditit exponential format ('PREC.E15' trace flag is set)
// b)  6-ditit fixed-point format ('PREC.F6' trace flag is set)
//
// This function checks trace flags every time it is called.
void tracerowautoprec(RMatrix *a, ae_int_t i, ae_int_t j0, ae_int_t j1, ae_state *_state) {
   ae_int_t j;
   ae_int_t prectouse;
// Determine precision to use
   prectouse = 0;
   if (ae_is_trace_enabled("PREC.E15")) {
      prectouse = 1;
   }
   if (ae_is_trace_enabled("PREC.F6")) {
      prectouse = 2;
   }
// Output
   ae_trace("[ ");
   for (j = j0; j < j1; j++) {
      if (prectouse == 0) {
         ae_trace("%14.6e", a->xyR[i][j]);
      }
      if (prectouse == 1) {
         ae_trace("%23.15e", a->xyR[i][j]);
      }
      if (prectouse == 2) {
         ae_trace("%13.6f", a->xyR[i][j]);
      }
      if (j < j1 - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

// Unscales/unshifts vector A[N] by computing A*Scl+Sft and outputs result to
// trace log using either:
// a)  6-digit exponential format (no trace flags is set)
// b) 15-ditit exponential format ('PREC.E15' trace flag is set)
// b)  6-ditit fixed-point format ('PREC.F6' trace flag is set)
//
// This function checks trace flags every time it is called.
// Both Scl and Sft can be omitted.
void tracevectorunscaledunshiftedautoprec(RVector *x, ae_int_t n, RVector *scl, bool applyscl, RVector *sft, bool applysft, ae_state *_state) {
   ae_int_t i;
   ae_int_t prectouse;
   double v;
// Determine precision to use
   prectouse = 0;
   if (ae_is_trace_enabled("PREC.E15")) {
      prectouse = 1;
   }
   if (ae_is_trace_enabled("PREC.F6")) {
      prectouse = 2;
   }
// Output
   ae_trace("[ ");
   for (i = 0; i < n; i++) {
      v = x->xR[i];
      if (applyscl) {
         v = v * scl->xR[i];
      }
      if (applysft) {
         v = v + sft->xR[i];
      }
      if (prectouse == 0) {
         ae_trace("%14.6e", v);
      }
      if (prectouse == 1) {
         ae_trace("%23.15e", v);
      }
      if (prectouse == 2) {
         ae_trace("%13.6f", v);
      }
      if (i < n - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

// Outputs vector of 1-norms of rows [I0,I1-1] of A[I0...I1-1,J0...J1-1]   to
// trace log using either:
// a)  6-digit exponential format (no trace flags is set)
// b) 15-ditit exponential format ('PREC.E15' trace flag is set)
// b)  6-ditit fixed-point format ('PREC.F6' trace flag is set)
//
// This function checks trace flags every time it is called.
void tracerownrm1autoprec(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t j0, ae_int_t j1, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_int_t prectouse;
// Determine precision to use
   prectouse = 0;
   if (ae_is_trace_enabled("PREC.E15")) {
      prectouse = 1;
   }
   if (ae_is_trace_enabled("PREC.F6")) {
      prectouse = 2;
   }
// Output
   ae_trace("[ ");
   for (i = i0; i < i1; i++) {
      v = 0.0;
      for (j = j0; j < j1; j++) {
         v = maxreal(v, ae_fabs(a->xyR[i][j], _state), _state);
      }
      if (prectouse == 0) {
         ae_trace("%14.6e", v);
      }
      if (prectouse == 1) {
         ae_trace("%23.15e", v);
      }
      if (prectouse == 2) {
         ae_trace("%13.6f", v);
      }
      if (i < i1 - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

// Outputs vector A[I0,I1-1] to trace log using E8 precision
void tracevectore6(RVector *a, ae_int_t i0, ae_int_t i1, ae_state *_state) {
   ae_int_t i;
   ae_trace("[ ");
   for (i = i0; i < i1; i++) {
      ae_trace("%14.6e", a->xR[i]);
      if (i < i1 - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

// Outputs vector A[I0,I1-1] to trace log using E8 or E15 precision
void tracevectore615(RVector *a, ae_int_t i0, ae_int_t i1, bool usee15, ae_state *_state) {
   ae_int_t i;
   ae_trace("[ ");
   for (i = i0; i < i1; i++) {
      if (usee15) {
         ae_trace("%23.15e", a->xR[i]);
      } else {
         ae_trace("%14.6e", a->xR[i]);
      }
      if (i < i1 - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

// Outputs vector of 1-norms of rows [I0,I1-1] of A[I0...I1-1,J0...J1-1]   to
// trace log using E8 precision
void tracerownrm1e6(RMatrix *a, ae_int_t i0, ae_int_t i1, ae_int_t j0, ae_int_t j1, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_trace("[ ");
   for (i = i0; i < i1; i++) {
      v = 0.0;
      for (j = j0; j < j1; j++) {
         v = maxreal(v, ae_fabs(a->xyR[i][j], _state), _state);
      }
      ae_trace("%14.6e", v);
      if (i < i1 - 1) {
         ae_trace(" ");
      }
   }
   ae_trace(" ]");
}

void apbuffers_init(void *_p, ae_state *_state, bool make_automatic) {
   apbuffers *p = (apbuffers *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->ba0, 0, DT_BOOL, _state, make_automatic);
   ae_vector_init(&p->ia0, 0, DT_INT, _state, make_automatic);
   ae_vector_init(&p->ia1, 0, DT_INT, _state, make_automatic);
   ae_vector_init(&p->ia2, 0, DT_INT, _state, make_automatic);
   ae_vector_init(&p->ia3, 0, DT_INT, _state, make_automatic);
   ae_vector_init(&p->ra0, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->ra1, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->ra2, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->ra3, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->rm0, 0, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->rm1, 0, 0, DT_REAL, _state, make_automatic);
}

void apbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   apbuffers *dst = (apbuffers *)_dst;
   apbuffers *src = (apbuffers *)_src;
   ae_vector_copy(&dst->ba0, &src->ba0, _state, make_automatic);
   ae_vector_copy(&dst->ia0, &src->ia0, _state, make_automatic);
   ae_vector_copy(&dst->ia1, &src->ia1, _state, make_automatic);
   ae_vector_copy(&dst->ia2, &src->ia2, _state, make_automatic);
   ae_vector_copy(&dst->ia3, &src->ia3, _state, make_automatic);
   ae_vector_copy(&dst->ra0, &src->ra0, _state, make_automatic);
   ae_vector_copy(&dst->ra1, &src->ra1, _state, make_automatic);
   ae_vector_copy(&dst->ra2, &src->ra2, _state, make_automatic);
   ae_vector_copy(&dst->ra3, &src->ra3, _state, make_automatic);
   ae_matrix_copy(&dst->rm0, &src->rm0, _state, make_automatic);
   ae_matrix_copy(&dst->rm1, &src->rm1, _state, make_automatic);
}

void apbuffers_free(void *_p, bool make_automatic) {
   apbuffers *p = (apbuffers *)_p;
   ae_touch_ptr((void *)p);
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

void sboolean_init(void *_p, ae_state *_state, bool make_automatic) {
   sboolean *p = (sboolean *)_p;
   ae_touch_ptr((void *)p);
}

void sboolean_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   sboolean *dst = (sboolean *)_dst;
   sboolean *src = (sboolean *)_src;
   dst->val = src->val;
}

void sboolean_free(void *_p, bool make_automatic) {
   sboolean *p = (sboolean *)_p;
   ae_touch_ptr((void *)p);
}

void sbooleanarray_init(void *_p, ae_state *_state, bool make_automatic) {
   sbooleanarray *p = (sbooleanarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->val, 0, DT_BOOL, _state, make_automatic);
}

void sbooleanarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   sbooleanarray *dst = (sbooleanarray *)_dst;
   sbooleanarray *src = (sbooleanarray *)_src;
   ae_vector_copy(&dst->val, &src->val, _state, make_automatic);
}

void sbooleanarray_free(void *_p, bool make_automatic) {
   sbooleanarray *p = (sbooleanarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_free(&p->val, make_automatic);
}

void sinteger_init(void *_p, ae_state *_state, bool make_automatic) {
   sinteger *p = (sinteger *)_p;
   ae_touch_ptr((void *)p);
}

void sinteger_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   sinteger *dst = (sinteger *)_dst;
   sinteger *src = (sinteger *)_src;
   dst->val = src->val;
}

void sinteger_free(void *_p, bool make_automatic) {
   sinteger *p = (sinteger *)_p;
   ae_touch_ptr((void *)p);
}

void sintegerarray_init(void *_p, ae_state *_state, bool make_automatic) {
   sintegerarray *p = (sintegerarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->val, 0, DT_INT, _state, make_automatic);
}

void sintegerarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   sintegerarray *dst = (sintegerarray *)_dst;
   sintegerarray *src = (sintegerarray *)_src;
   ae_vector_copy(&dst->val, &src->val, _state, make_automatic);
}

void sintegerarray_free(void *_p, bool make_automatic) {
   sintegerarray *p = (sintegerarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_free(&p->val, make_automatic);
}

void sreal_init(void *_p, ae_state *_state, bool make_automatic) {
   sreal *p = (sreal *)_p;
   ae_touch_ptr((void *)p);
}

void sreal_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   sreal *dst = (sreal *)_dst;
   sreal *src = (sreal *)_src;
   dst->val = src->val;
}

void sreal_free(void *_p, bool make_automatic) {
   sreal *p = (sreal *)_p;
   ae_touch_ptr((void *)p);
}

void srealarray_init(void *_p, ae_state *_state, bool make_automatic) {
   srealarray *p = (srealarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->val, 0, DT_REAL, _state, make_automatic);
}

void srealarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   srealarray *dst = (srealarray *)_dst;
   srealarray *src = (srealarray *)_src;
   ae_vector_copy(&dst->val, &src->val, _state, make_automatic);
}

void srealarray_free(void *_p, bool make_automatic) {
   srealarray *p = (srealarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_free(&p->val, make_automatic);
}

void scomplex_init(void *_p, ae_state *_state, bool make_automatic) {
   scomplex *p = (scomplex *)_p;
   ae_touch_ptr((void *)p);
}

void scomplex_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   scomplex *dst = (scomplex *)_dst;
   scomplex *src = (scomplex *)_src;
   dst->val = src->val;
}

void scomplex_free(void *_p, bool make_automatic) {
   scomplex *p = (scomplex *)_p;
   ae_touch_ptr((void *)p);
}

void scomplexarray_init(void *_p, ae_state *_state, bool make_automatic) {
   scomplexarray *p = (scomplexarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->val, 0, DT_COMPLEX, _state, make_automatic);
}

void scomplexarray_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   scomplexarray *dst = (scomplexarray *)_dst;
   scomplexarray *src = (scomplexarray *)_src;
   ae_vector_copy(&dst->val, &src->val, _state, make_automatic);
}

void scomplexarray_free(void *_p, bool make_automatic) {
   scomplexarray *p = (scomplexarray *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_free(&p->val, make_automatic);
}
} // end of namespace alglib_impl

// === ABLASF Package ===
namespace alglib_impl {
#ifdef ALGLIB_NO_FAST_KERNELS
// Computes dot product (X,Y) for elements [0,N) of X[] and Y[]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//     Y       -   array[N], vector to process
//
// Result:
//     (X,Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotv(ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t i;
   double result;
   result = 0.0;
   for (i = 0; i < n; i++) {
      result = result + x->xR[i] * y->xR[i];
   }
   return result;
}

// Computes dot product (X,A[i]) for elements [0,N) of vector X[] and row A[i,*]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//     A       -   array[?,N], matrix to process
//     I       -   row index
//
// Result:
//     (X,Ai)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotvr(ae_int_t n, RVector *x, RMatrix *a, ae_int_t i, ae_state *_state) {
   ae_int_t j;
   double result;
   result = 0.0;
   for (j = 0; j < n; j++) {
      result = result + x->xR[j] * a->xyR[i][j];
   }
   return result;
}

// Computes dot product (X,A[i]) for rows A[ia,*] and B[ib,*]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//     A       -   array[?,N], matrix to process
//     I       -   row index
//
// Result:
//     (X,Ai)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotrr(ae_int_t n, RMatrix *a, ae_int_t ia, RMatrix *b, ae_int_t ib, ae_state *_state) {
   ae_int_t j;
   double result;
   result = 0.0;
   for (j = 0; j < n; j++) {
      result = result + a->xyR[ia][j] * b->xyR[ib][j];
   }
   return result;
}

// Computes dot product (X,X) for elements [0,N) of X[]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//
// Result:
//     (X,X)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rdotv2(ae_int_t n, RVector *x, ae_state *_state) {
   ae_int_t i;
   double v;
   double result;
   result = 0.0;
   for (i = 0; i < n; i++) {
      v = x->xR[i];
      result = result + v * v;
   }
   return result;
}

// Performs inplace addition of Y[] to X[]
//
// Inputs:
//     N       -   vector length
//     Alpha   -   multiplier
//     Y       -   array[N], vector to process
//     X       -   array[N], vector to process
//
// Result:
//     X := X + alpha*Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddv(ae_int_t n, double alpha, RVector *y, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = x->xR[i] + alpha * y->xR[i];
   }
}

// Performs inplace addition of Y[] to X[]
//
// Inputs:
//     N       -   vector length
//     Alpha   -   multiplier
//     Y       -   source vector
//     OffsY   -   source offset
//     X       -   destination vector
//     OffsX   -   destination offset
//
// Result:
//     X := X + alpha*Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddvx(ae_int_t n, double alpha, RVector *y, ae_int_t offsy, RVector *x, ae_int_t offsx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[offsx + i] = x->xR[offsx + i] + alpha * y->xR[offsy + i];
   }
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Performs inplace addition of vector Y[] to column X[]
//
// Inputs:
//     N       -   vector length
//     Alpha   -   multiplier
//     Y       -   vector to add
//     X       -   target column ColIdx
//
// Result:
//     X := X + alpha*Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddvc(ae_int_t n, double alpha, RVector *y, RMatrix *x, ae_int_t colidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[i][colidx] = x->xyR[i][colidx] + alpha * y->xR[i];
   }
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Performs inplace addition of vector Y[] to row X[]
//
// Inputs:
//     N       -   vector length
//     Alpha   -   multiplier
//     Y       -   vector to add
//     X       -   target row RowIdx
//
// Result:
//     X := X + alpha*Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddvr(ae_int_t n, double alpha, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[rowidx][i] = x->xyR[rowidx][i] + alpha * y->xR[i];
   }
}

// Performs componentwise multiplication of vector X[] by vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target vector
//
// Result:
//     X := componentwise(X*Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemulv(ae_int_t n, RVector *y, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = x->xR[i] * y->xR[i];
   }
}

// Performs componentwise multiplication of row X[] by vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target row RowIdx
//
// Result:
//     X := componentwise(X*Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemulvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[rowidx][i] = x->xyR[rowidx][i] * y->xR[i];
   }
}

// Performs componentwise multiplication of row X[] by vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target row RowIdx
//
// Result:
//     X := componentwise(X*Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemulrv(ae_int_t n, RMatrix *y, ae_int_t rowidx, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = x->xR[i] * y->xyR[rowidx][i];
   }
}

// Performs componentwise max of vector X[] and vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target vector
//
// Result:
//     X := componentwise_max(X,Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemaxv(ae_int_t n, RVector *y, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = maxreal(x->xR[i], y->xR[i], _state);
   }
}

// Performs componentwise max of row X[] and vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target row RowIdx
//
// Result:
//     X := componentwise_max(X,Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemaxvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[rowidx][i] = maxreal(x->xyR[rowidx][i], y->xR[i], _state);
   }
}

// Performs componentwise max of row X[I] and vector Y[]
//
// Inputs:
//     N       -   vector length
//     X       -   matrix, I-th row is source
//     X       -   target row RowIdx
//
// Result:
//     Y := componentwise_max(Y,X)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergemaxrv(ae_int_t n, RMatrix *x, ae_int_t rowidx, RVector *y, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      y->xR[i] = maxreal(y->xR[i], x->xyR[rowidx][i], _state);
   }
}

// Performs componentwise max of vector X[] and vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target vector
//
// Result:
//     X := componentwise_max(X,Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergeminv(ae_int_t n, RVector *y, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = minreal(x->xR[i], y->xR[i], _state);
   }
}

// Performs componentwise max of row X[] and vector Y[]
//
// Inputs:
//     N       -   vector length
//     Y       -   vector to multiply by
//     X       -   target row RowIdx
//
// Result:
//     X := componentwise_max(X,Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergeminvr(ae_int_t n, RVector *y, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[rowidx][i] = minreal(x->xyR[rowidx][i], y->xR[i], _state);
   }
}

// Performs componentwise max of row X[I] and vector Y[]
//
// Inputs:
//     N       -   vector length
//     X       -   matrix, I-th row is source
//     X       -   target row RowIdx
//
// Result:
//     X := componentwise_max(X,Y)
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmergeminrv(ae_int_t n, RMatrix *x, ae_int_t rowidx, RVector *y, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      y->xR[i] = minreal(y->xR[i], x->xyR[rowidx][i], _state);
   }
}

// Performs inplace addition of Y[RIdx,...] to X[]
//
// Inputs:
//     N       -   vector length
//     Alpha   -   multiplier
//     Y       -   array[?,N], matrix whose RIdx-th row is added
//     RIdx    -   row index
//     X       -   array[N], vector to process
//
// Result:
//     X := X + alpha*Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddrv(ae_int_t n, double alpha, RMatrix *y, ae_int_t ridx, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = x->xR[i] + alpha * y->xyR[ridx][i];
   }
}

// Performs inplace addition of Y[RIdx,...] to X[RIdxDst]
//
// Inputs:
//     N       -   vector length
//     Alpha   -   multiplier
//     Y       -   array[?,N], matrix whose RIdxSrc-th row is added
//     RIdxSrc -   source row index
//     X       -   array[?,N], matrix whose RIdxDst-th row is target
//     RIdxDst -   destination row index
//
// Result:
//     X := X + alpha*Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void raddrr(ae_int_t n, double alpha, RMatrix *y, ae_int_t ridxsrc, RMatrix *x, ae_int_t ridxdst, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[ridxdst][i] = x->xyR[ridxdst][i] + alpha * y->xyR[ridxsrc][i];
   }
}

// Performs inplace multiplication of X[] by V
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//     V       -   multiplier
//
// Outputs:
//     X       -   elements 0...N-1 multiplied by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmulv(ae_int_t n, double v, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = x->xR[i] * v;
   }
}

// Performs inplace multiplication of X[] by V
//
// Inputs:
//     N       -   row length
//     X       -   array[?,N], row to process
//     V       -   multiplier
//
// Outputs:
//     X       -   elements 0...N-1 of row RowIdx are multiplied by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmulr(ae_int_t n, double v, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xyR[rowidx][i] = x->xyR[rowidx][i] * v;
   }
}

// Performs inplace multiplication of X[OffsX:OffsX+N-1] by V
//
// Inputs:
//     N       -   subvector length
//     X       -   vector to process
//     V       -   multiplier
//
// Outputs:
//     X       -   elements OffsX:OffsX+N-1 multiplied by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rmulvx(ae_int_t n, double v, RVector *x, ae_int_t offsx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[offsx + i] = x->xR[offsx + i] * v;
   }
}

// Returns maximum X
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//
// Outputs:
//     max(X[i])
//     zero for N=0
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxv(ae_int_t n, RVector *x, ae_state *_state) {
   ae_int_t i;
   double v;
   double result;
   if (n <= 0) {
      result = 0.0;
      return result;
   }
   result = x->xR[0];
   for (i = 1; i < n; i++) {
      v = x->xR[i];
      if (v > result) {
         result = v;
      }
   }
   return result;
}

// Returns maximum |X|
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], vector to process
//
// Outputs:
//     max(|X[i]|)
//     zero for N=0
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxabsv(ae_int_t n, RVector *x, ae_state *_state) {
   ae_int_t i;
   double v;
   double result;
   result = 0.0;
   for (i = 0; i < n; i++) {
      v = ae_fabs(x->xR[i], _state);
      if (v > result) {
         result = v;
      }
   }
   return result;
}

// Returns maximum X
//
// Inputs:
//     N       -   vector length
//     X       -   matrix to process, RowIdx-th row is processed
//
// Outputs:
//     max(X[RowIdx,i])
//     zero for N=0
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxr(ae_int_t n, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   double v;
   double result;
   if (n <= 0) {
      result = 0.0;
      return result;
   }
   result = x->xyR[rowidx][0];
   for (i = 1; i < n; i++) {
      v = x->xyR[rowidx][i];
      if (v > result) {
         result = v;
      }
   }
   return result;
}

// Returns maximum |X|
//
// Inputs:
//     N       -   vector length
//     X       -   matrix to process, RowIdx-th row is processed
//
// Outputs:
//     max(|X[RowIdx,i]|)
//     zero for N=0
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
double rmaxabsr(ae_int_t n, RMatrix *x, ae_int_t rowidx, ae_state *_state) {
   ae_int_t i;
   double v;
   double result;
   result = 0.0;
   for (i = 0; i < n; i++) {
      v = ae_fabs(x->xyR[rowidx][i], _state);
      if (v > result) {
         result = v;
      }
   }
   return result;
}

// Sets vector X[] to V
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     X       -   array[N]
//
// Outputs:
//     X       -   leading N elements are replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void bsetv(ae_int_t n, bool v, BVector *x, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      x->xB[j] = v;
   }
}

// Sets vector X[] to V
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     X       -   array[N]
//
// Outputs:
//     X       -   leading N elements are replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void isetv(ae_int_t n, ae_int_t v, ZVector *x, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      x->xZ[j] = v;
   }
}

// Sets vector X[] to V
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     X       -   array[N]
//
// Outputs:
//     X       -   leading N elements are replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetv(ae_int_t n, double v, RVector *x, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      x->xR[j] = v;
   }
}

// Sets X[OffsX:OffsX+N-1] to V
//
// Inputs:
//     N       -   subvector length
//     V       -   value to set
//     X       -   array[N]
//
// Outputs:
//     X       -   X[OffsX:OffsX+N-1] is replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetvx(ae_int_t n, double v, RVector *x, ae_int_t offsx, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      x->xR[offsx + j] = v;
   }
}

// Sets matrix A[] to V
//
// Inputs:
//     M, N    -   rows/cols count
//     V       -   value to set
//     A       -   array[M,N]
//
// Outputs:
//     A       -   leading M rows, N cols are replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetm(ae_int_t m, ae_int_t n, double v, RMatrix *a, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         a->xyR[i][j] = v;
      }
   }
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Sets vector X[] to V, reallocating X[] if too small
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     X       -   possibly preallocated array
//
// Outputs:
//     X       -   leading N elements are replaced by V; array is reallocated
//                 if its length is less than N.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void bsetallocv(ae_int_t n, bool v, BVector *x, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
   bsetv(n, v, x, _state);
}

// Sets vector X[] to V, reallocating X[] if too small
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     X       -   possibly preallocated array
//
// Outputs:
//     X       -   leading N elements are replaced by V; array is reallocated
//                 if its length is less than N.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void isetallocv(ae_int_t n, ae_int_t v, ZVector *x, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
   isetv(n, v, x, _state);
}

// Sets vector X[] to V, reallocating X[] if too small
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     X       -   possibly preallocated array
//
// Outputs:
//     X       -   leading N elements are replaced by V; array is reallocated
//                 if its length is less than N.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetallocv(ae_int_t n, double v, RVector *x, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
   rsetv(n, v, x, _state);
}

// Sets vector A[] to V, reallocating A[] if too small.
//
// Inputs:
//     M       -   rows count
//     N       -   cols count
//     V       -   value to set
//     A       -   possibly preallocated matrix
//
// Outputs:
//     A       -   leading M rows, N cols are replaced by V; the matrix is
//                 reallocated if its rows/cols count is less than M/N.
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetallocm(ae_int_t m, ae_int_t n, double v, RMatrix *a, ae_state *_state) {
   if (a->rows < m || a->cols < n) {
      ae_matrix_set_length(a, m, n, _state);
   }
   rsetm(m, n, v, a, _state);
}

// Reallocates X[] if its length is less than required value. Does not change
// its length and contents if it is large enough.
//
// Inputs:
//     N       -   desired vector length
//     X       -   possibly preallocated array
//
// Outputs:
//     X       -   length(X) >= N
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void ballocv(ae_int_t n, BVector *x, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
}

// Reallocates X[] if its length is less than required value. Does not change
// its length and contents if it is large enough.
//
// Inputs:
//     N       -   desired vector length
//     X       -   possibly preallocated array
//
// Outputs:
//     X       -   length(X) >= N
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void iallocv(ae_int_t n, ZVector *x, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
}

// Reallocates X[] if its length is less than required value. Does not change
// its length and contents if it is large enough.
//
// Inputs:
//     N       -   desired vector length
//     X       -   possibly preallocated array
//
// Outputs:
//     X       -   length(X) >= N
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rallocv(ae_int_t n, RVector *x, ae_state *_state) {
   if (x->cnt < n) {
      ae_vector_set_length(x, n, _state);
   }
}

// Reallocates matrix if its rows or cols count is less than  required.  Does
// not change its size if it is exactly that size or larger.
//
// Inputs:
//     M       -   rows count
//     N       -   cols count
//     A       -   possibly preallocated matrix
//
// Outputs:
//     A       -   size is at least M*N
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rallocm(ae_int_t m, ae_int_t n, RMatrix *a, ae_state *_state) {
   if (a->rows < m || a->cols < n) {
      ae_matrix_set_length(a, m, n, _state);
   }
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Sets row I of A[,] to V
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     A       -   array[N,N] or larger
//     I       -   row index
//
// Outputs:
//     A       -   leading N elements of I-th row are replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetr(ae_int_t n, double v, RMatrix *a, ae_int_t i, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      a->xyR[i][j] = v;
   }
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Sets col J of A[,] to V
//
// Inputs:
//     N       -   vector length
//     V       -   value to set
//     A       -   array[N,N] or larger
//     J       -   col index
//
// Outputs:
//     A       -   leading N elements of I-th col are replaced by V
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rsetc(ae_int_t n, double v, RMatrix *a, ae_int_t j, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      a->xyR[i][j] = v;
   }
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Copies vector X[] to Y[]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     Y       -   preallocated array[N]
//
// Outputs:
//     Y       -   leading N elements are replaced by X
//
//
// NOTE: destination and source should NOT overlap
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void bcopyv(ae_int_t n, BVector *x, BVector *y, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      y->xB[j] = x->xB[j];
   }
}

// Copies vector X[] to Y[]
//
// Inputs:
//     N       -   vector length
//     X       -   source array
//     Y       -   preallocated array[N]
//
// Outputs:
//     Y       -   X copied to Y
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void icopyv(ae_int_t n, ZVector *x, ZVector *y, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      y->xZ[j] = x->xZ[j];
   }
}

// Copies vector X[] to Y[]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     Y       -   preallocated array[N]
//
// Outputs:
//     Y       -   leading N elements are replaced by X
//
//
// NOTE: destination and source should NOT overlap
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyv(ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      y->xR[j] = x->xR[j];
   }
}

// Copies vector X[] to Y[], extended version
//
// Inputs:
//     N       -   vector length
//     X       -   source array
//     OffsX   -   source offset
//     Y       -   preallocated array[N]
//     OffsY   -   destination offset
//
// Outputs:
//     Y       -   N elements starting from OffsY are replaced by X[OffsX:OffsX+N-1]
//
// NOTE: destination and source should NOT overlap
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void icopyvx(ae_int_t n, ZVector *x, ae_int_t offsx, ZVector *y, ae_int_t offsy, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      y->xZ[offsy + j] = x->xZ[offsx + j];
   }
}

// Copies vector X[] to Y[], extended version
//
// Inputs:
//     N       -   vector length
//     X       -   source array
//     OffsX   -   source offset
//     Y       -   preallocated array[N]
//     OffsY   -   destination offset
//
// Outputs:
//     Y       -   N elements starting from OffsY are replaced by X[OffsX:OffsX+N-1]
//
// NOTE: destination and source should NOT overlap
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyvx(ae_int_t n, RVector *x, ae_int_t offsx, RVector *y, ae_int_t offsy, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      y->xR[offsy + j] = x->xR[offsx + j];
   }
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Copies vector X[] to Y[], resizing Y[] if needed.
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     Y       -   possibly preallocated array[N] (resized if needed)
//
// Outputs:
//     Y       -   leading N elements are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void bcopyallocv(ae_int_t n, BVector *x, BVector *y, ae_state *_state) {
   if (y->cnt < n) {
      ae_vector_set_length(y, n, _state);
   }
   bcopyv(n, x, y, _state);
}

// Copies vector X[] to Y[], resizing Y[] if needed.
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     Y       -   possibly preallocated array[N] (resized if needed)
//
// Outputs:
//     Y       -   leading N elements are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void icopyallocv(ae_int_t n, ZVector *x, ZVector *y, ae_state *_state) {
   if (y->cnt < n) {
      ae_vector_set_length(y, n, _state);
   }
   icopyv(n, x, y, _state);
}

// Copies vector X[] to Y[], resizing Y[] if needed.
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     Y       -   possibly preallocated array[N] (resized if needed)
//
// Outputs:
//     Y       -   leading N elements are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyallocv(ae_int_t n, RVector *x, RVector *y, ae_state *_state) {
   if (y->cnt < n) {
      ae_vector_set_length(y, n, _state);
   }
   rcopyv(n, x, y, _state);
}

// Copies matrix X[] to Y[], resizing Y[] if needed. On resize, dimensions of
// Y[] are increased - but not decreased.
//
// Inputs:
//     M       -   rows count
//     N       -   cols count
//     X       -   array[M,N], source
//     Y       -   possibly preallocated array[M,N] (resized if needed)
//
// Outputs:
//     Y       -   leading [M,N] elements are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopym(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   if (m == 0 || n == 0) {
      return;
   }
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         y->xyR[i][j] = x->xyR[i][j];
      }
   }
}

// Copies matrix X[] to Y[], resizing Y[] if needed. On resize, dimensions of
// Y[] are increased - but not decreased.
//
// Inputs:
//     M       -   rows count
//     N       -   cols count
//     X       -   array[M,N], source
//     Y       -   possibly preallocated array[M,N] (resized if needed)
//
// Outputs:
//     Y       -   leading [M,N] elements are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyallocm(ae_int_t m, ae_int_t n, RMatrix *x, RMatrix *y, ae_state *_state) {
   if (m == 0 || n == 0) {
      return;
   }
   if (y->rows < m || y->cols < n) {
      ae_matrix_set_length(y, maxint(m, y->rows, _state), maxint(n, y->cols, _state), _state);
   }
   rcopym(m, n, x, y, _state);
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) actual size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void igrowv(ae_int_t newn, ZVector *x, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t oldn;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_INT, _state);
   if (x->cnt >= newn) {
      ae_frame_leave(_state);
      return;
   }
   oldn = x->cnt;
   newn = maxint(newn, iround(1.8 * oldn + 1, _state), _state);
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, newn, _state);
   icopyv(oldn, &oldx, x, _state);
   ae_frame_leave(_state);
}

// Grows X, i.e. changes its size in such a way that:
// a) contents is preserved
// b) new size is at least N
// c) actual size can be larger than N, so subsequent grow() calls can return
//    without reallocation
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
void rgrowv(ae_int_t newn, RVector *x, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t oldn;
   ae_frame_make(_state, &_frame_block);
   NewVector(oldx, 0, DT_REAL, _state);
   if (x->cnt >= newn) {
      ae_frame_leave(_state);
      return;
   }
   oldn = x->cnt;
   newn = maxint(newn, iround(1.8 * oldn + 1, _state), _state);
   ae_swap_vectors(x, &oldx);
   ae_vector_set_length(x, newn, _state);
   rcopyv(oldn, &oldx, x, _state);
   ae_frame_leave(_state);
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Performs copying with multiplication of V*X[] to Y[]
//
// Inputs:
//     N       -   vector length
//     V       -   multiplier
//     X       -   array[N], source
//     Y       -   preallocated array[N]
//
// Outputs:
//     Y       -   array[N], Y = V*X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopymulv(ae_int_t n, double v, RVector *x, RVector *y, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      y->xR[i] = v * x->xR[i];
   }
}

// Performs copying with multiplication of V*X[] to Y[I,*]
//
// Inputs:
//     N       -   vector length
//     V       -   multiplier
//     X       -   array[N], source
//     Y       -   preallocated array[?,N]
//     RIdx    -   destination row index
//
// Outputs:
//     Y       -   Y[RIdx,...] = V*X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopymulvr(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t ridx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      y->xyR[ridx][i] = v * x->xR[i];
   }
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Performs copying with multiplication of V*X[] to Y[*,J]
//
// Inputs:
//     N       -   vector length
//     V       -   multiplier
//     X       -   array[N], source
//     Y       -   preallocated array[N,?]
//     CIdx    -   destination rocol index
//
// Outputs:
//     Y       -   Y[RIdx,...] = V*X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopymulvc(ae_int_t n, double v, RVector *x, RMatrix *y, ae_int_t cidx, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      y->xyR[i][cidx] = v * x->xR[i];
   }
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Copies vector X[] to row I of A[,]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     A       -   preallocated 2D array large enough to store result
//     I       -   destination row index
//
// Outputs:
//     A       -   leading N elements of I-th row are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyvr(ae_int_t n, RVector *x, RMatrix *a, ae_int_t i, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      a->xyR[i][j] = x->xR[j];
   }
}

// Copies row I of A[,] to vector X[]
//
// Inputs:
//     N       -   vector length
//     A       -   2D array, source
//     I       -   source row index
//     X       -   preallocated destination
//
// Outputs:
//     X       -   array[N], destination
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyrv(ae_int_t n, RMatrix *a, ae_int_t i, RVector *x, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      x->xR[j] = a->xyR[i][j];
   }
}

// Copies row I of A[,] to row K of B[,].
//
// A[i,...] and B[k,...] may overlap.
//
// Inputs:
//     N       -   vector length
//     A       -   2D array, source
//     I       -   source row index
//     B       -   preallocated destination
//     K       -   destination row index
//
// Outputs:
//     B       -   row K overwritten
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyrr(ae_int_t n, RMatrix *a, ae_int_t i, RMatrix *b, ae_int_t k, ae_state *_state) {
   ae_int_t j;
   for (j = 0; j < n; j++) {
      b->xyR[k][j] = a->xyR[i][j];
   }
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Copies vector X[] to column J of A[,]
//
// Inputs:
//     N       -   vector length
//     X       -   array[N], source
//     A       -   preallocated 2D array large enough to store result
//     J       -   destination col index
//
// Outputs:
//     A       -   leading N elements of J-th column are replaced by X
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopyvc(ae_int_t n, RVector *x, RMatrix *a, ae_int_t j, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      a->xyR[i][j] = x->xR[i];
   }
}

// Copies column J of A[,] to vector X[]
//
// Inputs:
//     N       -   vector length
//     A       -   source 2D array
//     J       -   source col index
//
// Outputs:
//     X       -   preallocated array[N], destination
// ALGLIB: Copyright 20.01.2020 by Sergey Bochkanov
void rcopycv(ae_int_t n, RMatrix *a, ae_int_t j, RVector *x, ae_state *_state) {
   ae_int_t i;
   for (i = 0; i < n; i++) {
      x->xR[i] = a->xyR[i][j];
   }
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Matrix-vector product: y := alpha*op(A)*x + beta*y
//
// NOTE: this  function  expects  Y  to  be  large enough to store result. No
//       automatic preallocation happens for  smaller  arrays.  No  integrity
//       checks is performed for sizes of A, x, y.
//
// Inputs:
//     M   -   number of rows of op(A)
//     N   -   number of columns of op(A)
//     Alpha-  coefficient
//     A   -   source matrix
//     OpA -   operation type:
//             * OpA=0     =>  op(A) = A
//             * OpA=1     =>  op(A) = A^T
//     X   -   input vector, has at least N elements
//     Beta-   coefficient
//     Y   -   preallocated output array, has at least M elements
//
// Outputs:
//     Y   -   vector which stores result
//
// HANDLING OF SPECIAL CASES:
//     * if M=0, then subroutine does nothing. It does not even touch arrays.
//     * if N=0 or Alpha=0.0, then:
//       * if Beta=0, then Y is filled by zeros. A and X are  not  referenced
//         at all. Initial values of Y are ignored (we do not  multiply  Y by
//         zero, we just rewrite it by zeros)
//       * if Beta != 0, then Y is replaced by Beta*Y
//     * if M>0, N>0, Alpha != 0, but  Beta=0,  then  Y  is  replaced  by  A*x;
//        initial state of Y is ignored (rewritten by  A*x,  without  initial
//        multiplication by zeros).
// ALGLIB Routine: Copyright 01.09.2021 by Sergey Bochkanov
void rgemv(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t opa, RVector *x, double beta, RVector *y, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
// Properly premultiply Y by Beta.
//
// Quick exit for M=0, N=0 or Alpha=0.
// After this block we have M>0, N>0, Alpha != 0.
   if (m <= 0) {
      return;
   }
   if (beta != 0.0) {
      rmulv(m, beta, y, _state);
   } else {
      rsetv(m, 0.0, y, _state);
   }
   if (n <= 0 || alpha == 0.0) {
      return;
   }
// Generic code
   if (opa == 0) {
   // y += A*x
      for (i = 0; i < m; i++) {
         v = 0.0;
         for (j = 0; j < n; j++) {
            v = v + a->xyR[i][j] * x->xR[j];
         }
         y->xR[i] = alpha * v + y->xR[i];
      }
      return;
   }
   if (opa == 1) {
   // y += A^T*x
      for (i = 0; i < n; i++) {
         v = alpha * x->xR[i];
         for (j = 0; j < m; j++) {
            y->xR[j] = y->xR[j] + v * a->xyR[i][j];
         }
      }
      return;
   }
}

// Matrix-vector product: y := alpha*op(A)*x + beta*y
//
// Here x, y, A are subvectors/submatrices of larger vectors/matrices.
//
// NOTE: this  function  expects  Y  to  be  large enough to store result. No
//       automatic preallocation happens for  smaller  arrays.  No  integrity
//       checks is performed for sizes of A, x, y.
//
// Inputs:
//     M   -   number of rows of op(A)
//     N   -   number of columns of op(A)
//     Alpha-  coefficient
//     A   -   source matrix
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     OpA -   operation type:
//             * OpA=0     =>  op(A) = A
//             * OpA=1     =>  op(A) = A^T
//     X   -   input vector, has at least N+IX elements
//     IX  -   subvector offset
//     Beta-   coefficient
//     Y   -   preallocated output array, has at least M+IY elements
//     IY  -   subvector offset
//
// Outputs:
//     Y   -   vector which stores result
//
// HANDLING OF SPECIAL CASES:
//     * if M=0, then subroutine does nothing. It does not even touch arrays.
//     * if N=0 or Alpha=0.0, then:
//       * if Beta=0, then Y is filled by zeros. A and X are  not  referenced
//         at all. Initial values of Y are ignored (we do not  multiply  Y by
//         zero, we just rewrite it by zeros)
//       * if Beta != 0, then Y is replaced by Beta*Y
//     * if M>0, N>0, Alpha != 0, but  Beta=0,  then  Y  is  replaced  by  A*x;
//        initial state of Y is ignored (rewritten by  A*x,  without  initial
//        multiplication by zeros).
// ALGLIB Routine: Copyright 01.09.2021 by Sergey Bochkanov
void rgemvx(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
// Properly premultiply Y by Beta.
//
// Quick exit for M=0, N=0 or Alpha=0.
// After this block we have M>0, N>0, Alpha != 0.
   if (m <= 0) {
      return;
   }
   if (beta != 0.0) {
      rmulvx(m, beta, y, iy, _state);
   } else {
      rsetvx(m, 0.0, y, iy, _state);
   }
   if (n <= 0 || alpha == 0.0) {
      return;
   }
// Generic code
   if (opa == 0) {
   // y += A*x
      for (i = 0; i < m; i++) {
         v = 0.0;
         for (j = 0; j < n; j++) {
            v = v + a->xyR[ia + i][ja + j] * x->xR[ix + j];
         }
         y->xR[iy + i] = alpha * v + y->xR[iy + i];
      }
      return;
   }
   if (opa == 1) {
   // y += A^T*x
      for (i = 0; i < n; i++) {
         v = alpha * x->xR[ix + i];
         for (j = 0; j < m; j++) {
            y->xR[iy + j] = y->xR[iy + j] + v * a->xyR[ia + i][ja + j];
         }
      }
      return;
   }
}

// Rank-1 correction: A := A + alpha*u*v'
//
// NOTE: this  function  expects  A  to  be  large enough to store result. No
//       automatic preallocation happens for  smaller  arrays.  No  integrity
//       checks is performed for sizes of A, u, v.
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   target MxN matrix
//     Alpha-  coefficient
//     U   -   vector #1
//     V   -   vector #2
// ALGLIB Routine: Copyright 07.09.2021 by Sergey Bochkanov
void rger(ae_int_t m, ae_int_t n, double alpha, RVector *u, RVector *v, RMatrix *a, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double s;
   if ((m <= 0 || n <= 0) || alpha == 0.0) {
      return;
   }
   for (i = 0; i < m; i++) {
      s = alpha * u->xR[i];
      for (j = 0; j < n; j++) {
         a->xyR[i][j] = a->xyR[i][j] + s * v->xR[j];
      }
   }
}

// This subroutine solves linear system op(A)*x=b where:
// * A is NxN upper/lower triangular/unitriangular matrix
// * X and B are Nx1 vectors
// * "op" may be identity transformation or transposition
//
// Solution replaces X.
//
// IMPORTANT: * no overflow/underflow/denegeracy tests is performed.
//            * no integrity checks for operand sizes, out-of-bounds accesses
//              and so on is performed
//
// Inputs:
//     N   -   matrix size, N >= 0
//     A       -   matrix, actial matrix is stored in A[IA:IA+N-1,JA:JA+N-1]
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     IsUpper -   whether matrix is upper triangular
//     IsUnit  -   whether matrix is unitriangular
//     OpType  -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//     X       -   right part, actual vector is stored in X[IX:IX+N-1]
//     IX      -   offset
//
// Outputs:
//     X       -   solution replaces elements X[IX:IX+N-1]
// ALGLIB Routine: Copyright 07.09.2021 by Sergey Bochkanov
void rtrsvx(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   if (n <= 0) {
      return;
   }
   if (optype == 0 && isupper) {
      for (i = n - 1; i >= 0; i--) {
         v = x->xR[ix + i];
         for (j = i + 1; j < n; j++) {
            v = v - a->xyR[ia + i][ja + j] * x->xR[ix + j];
         }
         if (!isunit) {
            v = v / a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
      }
      return;
   }
   if (optype == 0 && !isupper) {
      for (i = 0; i < n; i++) {
         v = x->xR[ix + i];
         for (j = 0; j < i; j++) {
            v = v - a->xyR[ia + i][ja + j] * x->xR[ix + j];
         }
         if (!isunit) {
            v = v / a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
      }
      return;
   }
   if (optype == 1 && isupper) {
      for (i = 0; i < n; i++) {
         v = x->xR[ix + i];
         if (!isunit) {
            v = v / a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
         if (v == 0) {
            continue;
         }
         for (j = i + 1; j < n; j++) {
            x->xR[ix + j] = x->xR[ix + j] - v * a->xyR[ia + i][ja + j];
         }
      }
      return;
   }
   if (optype == 1 && !isupper) {
      for (i = n - 1; i >= 0; i--) {
         v = x->xR[ix + i];
         if (!isunit) {
            v = v / a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
         if (v == 0) {
            continue;
         }
         for (j = 0; j < i; j++) {
            x->xR[ix + j] = x->xR[ix + j] - v * a->xyR[ia + i][ja + j];
         }
      }
      return;
   }
   ae_assert(false, "rTRSVX: unexpected operation type", _state);
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixgerf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double ralpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixgerf(m, n, a, ia, ja, ralpha, u, iu, v, iv);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixrank1f(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixrank1f(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixlefttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixlefttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixrighttrsmf(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixrighttrsmf(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixherkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}

// Fast kernel
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
bool cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_ABLAS
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}

// RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
// with OpTypeA=0 and OpTypeB=0.
//
// Additional info:
// * this function requires that Alpha != 0 (assertion is thrown otherwise)
//
// Inputs:
//     M       -   matrix size, M>0
//     N       -   matrix size, N>0
//     K       -   matrix size, K>0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     B       -   matrix
//     IB      -   submatrix offset
//     JB      -   submatrix offset
//     Beta    -   coefficient
//     C       -   PREALLOCATED output matrix
//     IC      -   submatrix offset
//     JC      -   submatrix offset
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void rmatrixgemmk44v00(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   double v00;
   double v01;
   double v02;
   double v03;
   double v10;
   double v11;
   double v12;
   double v13;
   double v20;
   double v21;
   double v22;
   double v23;
   double v30;
   double v31;
   double v32;
   double v33;
   double a0;
   double a1;
   double a2;
   double a3;
   double b0;
   double b1;
   double b2;
   double b3;
   ae_int_t idxa0;
   ae_int_t idxa1;
   ae_int_t idxa2;
   ae_int_t idxa3;
   ae_int_t idxb0;
   ae_int_t idxb1;
   ae_int_t idxb2;
   ae_int_t idxb3;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t ik;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t jk;
   ae_int_t t;
   ae_int_t offsa;
   ae_int_t offsb;
   ae_assert(alpha != 0.0, "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
// if matrix size is zero
   if (m == 0 || n == 0) {
      return;
   }
// A*B
   i = 0;
   while (i < m) {
      j = 0;
      while (j < n) {
      // Choose between specialized 4x4 code and general code
         if (i + 4 <= m && j + 4 <= n) {
         // Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
         //
         // This submatrix is calculated as sum of K rank-1 products,
         // with operands cached in local variables in order to speed
         // up operations with arrays.
            idxa0 = ia + i + 0;
            idxa1 = ia + i + 1;
            idxa2 = ia + i + 2;
            idxa3 = ia + i + 3;
            offsa = ja;
            idxb0 = jb + j + 0;
            idxb1 = jb + j + 1;
            idxb2 = jb + j + 2;
            idxb3 = jb + j + 3;
            offsb = ib;
            v00 = 0.0;
            v01 = 0.0;
            v02 = 0.0;
            v03 = 0.0;
            v10 = 0.0;
            v11 = 0.0;
            v12 = 0.0;
            v13 = 0.0;
            v20 = 0.0;
            v21 = 0.0;
            v22 = 0.0;
            v23 = 0.0;
            v30 = 0.0;
            v31 = 0.0;
            v32 = 0.0;
            v33 = 0.0;
         // Different variants of internal loop
            for (t = 0; t < k; t++) {
               a0 = a->xyR[idxa0][offsa];
               a1 = a->xyR[idxa1][offsa];
               b0 = b->xyR[offsb][idxb0];
               b1 = b->xyR[offsb][idxb1];
               v00 = v00 + a0 * b0;
               v01 = v01 + a0 * b1;
               v10 = v10 + a1 * b0;
               v11 = v11 + a1 * b1;
               a2 = a->xyR[idxa2][offsa];
               a3 = a->xyR[idxa3][offsa];
               v20 = v20 + a2 * b0;
               v21 = v21 + a2 * b1;
               v30 = v30 + a3 * b0;
               v31 = v31 + a3 * b1;
               b2 = b->xyR[offsb][idxb2];
               b3 = b->xyR[offsb][idxb3];
               v22 = v22 + a2 * b2;
               v23 = v23 + a2 * b3;
               v32 = v32 + a3 * b2;
               v33 = v33 + a3 * b3;
               v02 = v02 + a0 * b2;
               v03 = v03 + a0 * b3;
               v12 = v12 + a1 * b2;
               v13 = v13 + a1 * b3;
               offsa = offsa + 1;
               offsb = offsb + 1;
            }
            if (beta == 0.0) {
               c->xyR[ic + i + 0][jc + j + 0] = alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = alpha * v33;
            } else {
               c->xyR[ic + i + 0][jc + j + 0] = beta * c->xyR[ic + i + 0][jc + j + 0] + alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = beta * c->xyR[ic + i + 0][jc + j + 1] + alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = beta * c->xyR[ic + i + 0][jc + j + 2] + alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = beta * c->xyR[ic + i + 0][jc + j + 3] + alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = beta * c->xyR[ic + i + 1][jc + j + 0] + alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = beta * c->xyR[ic + i + 1][jc + j + 1] + alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = beta * c->xyR[ic + i + 1][jc + j + 2] + alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = beta * c->xyR[ic + i + 1][jc + j + 3] + alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = beta * c->xyR[ic + i + 2][jc + j + 0] + alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = beta * c->xyR[ic + i + 2][jc + j + 1] + alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = beta * c->xyR[ic + i + 2][jc + j + 2] + alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = beta * c->xyR[ic + i + 2][jc + j + 3] + alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = beta * c->xyR[ic + i + 3][jc + j + 0] + alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = beta * c->xyR[ic + i + 3][jc + j + 1] + alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = beta * c->xyR[ic + i + 3][jc + j + 2] + alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = beta * c->xyR[ic + i + 3][jc + j + 3] + alpha * v33;
            }
         } else {
         // Determine submatrix [I0..I1]x[J0..J1] to process
            i0 = i;
            i1 = minint(i + 3, m - 1, _state);
            j0 = j;
            j1 = minint(j + 3, n - 1, _state);
         // Process submatrix
            for (ik = i0; ik <= i1; ik++) {
               for (jk = j0; jk <= j1; jk++) {
                  if (k == 0 || alpha == 0.0) {
                     v = 0.0;
                  } else {
                     v = ae_v_dotproduct(&a->xyR[ia + ik][ja], 1, &b->xyR[ib][jb + jk], b->stride, k);
                  }
                  if (beta == 0.0) {
                     c->xyR[ic + ik][jc + jk] = alpha * v;
                  } else {
                     c->xyR[ic + ik][jc + jk] = beta * c->xyR[ic + ik][jc + jk] + alpha * v;
                  }
               }
            }
         }
         j = j + 4;
      }
      i = i + 4;
   }
}

// RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
// with OpTypeA=0 and OpTypeB=1.
//
// Additional info:
// * this function requires that Alpha != 0 (assertion is thrown otherwise)
//
// Inputs:
//     M       -   matrix size, M>0
//     N       -   matrix size, N>0
//     K       -   matrix size, K>0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     B       -   matrix
//     IB      -   submatrix offset
//     JB      -   submatrix offset
//     Beta    -   coefficient
//     C       -   PREALLOCATED output matrix
//     IC      -   submatrix offset
//     JC      -   submatrix offset
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void rmatrixgemmk44v01(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   double v00;
   double v01;
   double v02;
   double v03;
   double v10;
   double v11;
   double v12;
   double v13;
   double v20;
   double v21;
   double v22;
   double v23;
   double v30;
   double v31;
   double v32;
   double v33;
   double a0;
   double a1;
   double a2;
   double a3;
   double b0;
   double b1;
   double b2;
   double b3;
   ae_int_t idxa0;
   ae_int_t idxa1;
   ae_int_t idxa2;
   ae_int_t idxa3;
   ae_int_t idxb0;
   ae_int_t idxb1;
   ae_int_t idxb2;
   ae_int_t idxb3;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t ik;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t jk;
   ae_int_t t;
   ae_int_t offsa;
   ae_int_t offsb;
   ae_assert(alpha != 0.0, "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
// if matrix size is zero
   if (m == 0 || n == 0) {
      return;
   }
// A*B'
   i = 0;
   while (i < m) {
      j = 0;
      while (j < n) {
      // Choose between specialized 4x4 code and general code
         if (i + 4 <= m && j + 4 <= n) {
         // Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
         //
         // This submatrix is calculated as sum of K rank-1 products,
         // with operands cached in local variables in order to speed
         // up operations with arrays.
            idxa0 = ia + i + 0;
            idxa1 = ia + i + 1;
            idxa2 = ia + i + 2;
            idxa3 = ia + i + 3;
            offsa = ja;
            idxb0 = ib + j + 0;
            idxb1 = ib + j + 1;
            idxb2 = ib + j + 2;
            idxb3 = ib + j + 3;
            offsb = jb;
            v00 = 0.0;
            v01 = 0.0;
            v02 = 0.0;
            v03 = 0.0;
            v10 = 0.0;
            v11 = 0.0;
            v12 = 0.0;
            v13 = 0.0;
            v20 = 0.0;
            v21 = 0.0;
            v22 = 0.0;
            v23 = 0.0;
            v30 = 0.0;
            v31 = 0.0;
            v32 = 0.0;
            v33 = 0.0;
            for (t = 0; t < k; t++) {
               a0 = a->xyR[idxa0][offsa];
               a1 = a->xyR[idxa1][offsa];
               b0 = b->xyR[idxb0][offsb];
               b1 = b->xyR[idxb1][offsb];
               v00 = v00 + a0 * b0;
               v01 = v01 + a0 * b1;
               v10 = v10 + a1 * b0;
               v11 = v11 + a1 * b1;
               a2 = a->xyR[idxa2][offsa];
               a3 = a->xyR[idxa3][offsa];
               v20 = v20 + a2 * b0;
               v21 = v21 + a2 * b1;
               v30 = v30 + a3 * b0;
               v31 = v31 + a3 * b1;
               b2 = b->xyR[idxb2][offsb];
               b3 = b->xyR[idxb3][offsb];
               v22 = v22 + a2 * b2;
               v23 = v23 + a2 * b3;
               v32 = v32 + a3 * b2;
               v33 = v33 + a3 * b3;
               v02 = v02 + a0 * b2;
               v03 = v03 + a0 * b3;
               v12 = v12 + a1 * b2;
               v13 = v13 + a1 * b3;
               offsa = offsa + 1;
               offsb = offsb + 1;
            }
            if (beta == 0.0) {
               c->xyR[ic + i + 0][jc + j + 0] = alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = alpha * v33;
            } else {
               c->xyR[ic + i + 0][jc + j + 0] = beta * c->xyR[ic + i + 0][jc + j + 0] + alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = beta * c->xyR[ic + i + 0][jc + j + 1] + alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = beta * c->xyR[ic + i + 0][jc + j + 2] + alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = beta * c->xyR[ic + i + 0][jc + j + 3] + alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = beta * c->xyR[ic + i + 1][jc + j + 0] + alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = beta * c->xyR[ic + i + 1][jc + j + 1] + alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = beta * c->xyR[ic + i + 1][jc + j + 2] + alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = beta * c->xyR[ic + i + 1][jc + j + 3] + alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = beta * c->xyR[ic + i + 2][jc + j + 0] + alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = beta * c->xyR[ic + i + 2][jc + j + 1] + alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = beta * c->xyR[ic + i + 2][jc + j + 2] + alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = beta * c->xyR[ic + i + 2][jc + j + 3] + alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = beta * c->xyR[ic + i + 3][jc + j + 0] + alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = beta * c->xyR[ic + i + 3][jc + j + 1] + alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = beta * c->xyR[ic + i + 3][jc + j + 2] + alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = beta * c->xyR[ic + i + 3][jc + j + 3] + alpha * v33;
            }
         } else {
         // Determine submatrix [I0..I1]x[J0..J1] to process
            i0 = i;
            i1 = minint(i + 3, m - 1, _state);
            j0 = j;
            j1 = minint(j + 3, n - 1, _state);
         // Process submatrix
            for (ik = i0; ik <= i1; ik++) {
               for (jk = j0; jk <= j1; jk++) {
                  if (k == 0 || alpha == 0.0) {
                     v = 0.0;
                  } else {
                     v = ae_v_dotproduct(&a->xyR[ia + ik][ja], 1, &b->xyR[ib + jk][jb], 1, k);
                  }
                  if (beta == 0.0) {
                     c->xyR[ic + ik][jc + jk] = alpha * v;
                  } else {
                     c->xyR[ic + ik][jc + jk] = beta * c->xyR[ic + ik][jc + jk] + alpha * v;
                  }
               }
            }
         }
         j = j + 4;
      }
      i = i + 4;
   }
}

// RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
// with OpTypeA=1 and OpTypeB=0.
//
// Additional info:
// * this function requires that Alpha != 0 (assertion is thrown otherwise)
//
// Inputs:
//     M       -   matrix size, M>0
//     N       -   matrix size, N>0
//     K       -   matrix size, K>0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     B       -   matrix
//     IB      -   submatrix offset
//     JB      -   submatrix offset
//     Beta    -   coefficient
//     C       -   PREALLOCATED output matrix
//     IC      -   submatrix offset
//     JC      -   submatrix offset
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void rmatrixgemmk44v10(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   double v00;
   double v01;
   double v02;
   double v03;
   double v10;
   double v11;
   double v12;
   double v13;
   double v20;
   double v21;
   double v22;
   double v23;
   double v30;
   double v31;
   double v32;
   double v33;
   double a0;
   double a1;
   double a2;
   double a3;
   double b0;
   double b1;
   double b2;
   double b3;
   ae_int_t idxa0;
   ae_int_t idxa1;
   ae_int_t idxa2;
   ae_int_t idxa3;
   ae_int_t idxb0;
   ae_int_t idxb1;
   ae_int_t idxb2;
   ae_int_t idxb3;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t ik;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t jk;
   ae_int_t t;
   ae_int_t offsa;
   ae_int_t offsb;
   ae_assert(alpha != 0.0, "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
// if matrix size is zero
   if (m == 0 || n == 0) {
      return;
   }
// A'*B
   i = 0;
   while (i < m) {
      j = 0;
      while (j < n) {
      // Choose between specialized 4x4 code and general code
         if (i + 4 <= m && j + 4 <= n) {
         // Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
         //
         // This submatrix is calculated as sum of K rank-1 products,
         // with operands cached in local variables in order to speed
         // up operations with arrays.
            idxa0 = ja + i + 0;
            idxa1 = ja + i + 1;
            idxa2 = ja + i + 2;
            idxa3 = ja + i + 3;
            offsa = ia;
            idxb0 = jb + j + 0;
            idxb1 = jb + j + 1;
            idxb2 = jb + j + 2;
            idxb3 = jb + j + 3;
            offsb = ib;
            v00 = 0.0;
            v01 = 0.0;
            v02 = 0.0;
            v03 = 0.0;
            v10 = 0.0;
            v11 = 0.0;
            v12 = 0.0;
            v13 = 0.0;
            v20 = 0.0;
            v21 = 0.0;
            v22 = 0.0;
            v23 = 0.0;
            v30 = 0.0;
            v31 = 0.0;
            v32 = 0.0;
            v33 = 0.0;
            for (t = 0; t < k; t++) {
               a0 = a->xyR[offsa][idxa0];
               a1 = a->xyR[offsa][idxa1];
               b0 = b->xyR[offsb][idxb0];
               b1 = b->xyR[offsb][idxb1];
               v00 = v00 + a0 * b0;
               v01 = v01 + a0 * b1;
               v10 = v10 + a1 * b0;
               v11 = v11 + a1 * b1;
               a2 = a->xyR[offsa][idxa2];
               a3 = a->xyR[offsa][idxa3];
               v20 = v20 + a2 * b0;
               v21 = v21 + a2 * b1;
               v30 = v30 + a3 * b0;
               v31 = v31 + a3 * b1;
               b2 = b->xyR[offsb][idxb2];
               b3 = b->xyR[offsb][idxb3];
               v22 = v22 + a2 * b2;
               v23 = v23 + a2 * b3;
               v32 = v32 + a3 * b2;
               v33 = v33 + a3 * b3;
               v02 = v02 + a0 * b2;
               v03 = v03 + a0 * b3;
               v12 = v12 + a1 * b2;
               v13 = v13 + a1 * b3;
               offsa = offsa + 1;
               offsb = offsb + 1;
            }
            if (beta == 0.0) {
               c->xyR[ic + i + 0][jc + j + 0] = alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = alpha * v33;
            } else {
               c->xyR[ic + i + 0][jc + j + 0] = beta * c->xyR[ic + i + 0][jc + j + 0] + alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = beta * c->xyR[ic + i + 0][jc + j + 1] + alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = beta * c->xyR[ic + i + 0][jc + j + 2] + alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = beta * c->xyR[ic + i + 0][jc + j + 3] + alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = beta * c->xyR[ic + i + 1][jc + j + 0] + alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = beta * c->xyR[ic + i + 1][jc + j + 1] + alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = beta * c->xyR[ic + i + 1][jc + j + 2] + alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = beta * c->xyR[ic + i + 1][jc + j + 3] + alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = beta * c->xyR[ic + i + 2][jc + j + 0] + alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = beta * c->xyR[ic + i + 2][jc + j + 1] + alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = beta * c->xyR[ic + i + 2][jc + j + 2] + alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = beta * c->xyR[ic + i + 2][jc + j + 3] + alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = beta * c->xyR[ic + i + 3][jc + j + 0] + alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = beta * c->xyR[ic + i + 3][jc + j + 1] + alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = beta * c->xyR[ic + i + 3][jc + j + 2] + alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = beta * c->xyR[ic + i + 3][jc + j + 3] + alpha * v33;
            }
         } else {
         // Determine submatrix [I0..I1]x[J0..J1] to process
            i0 = i;
            i1 = minint(i + 3, m - 1, _state);
            j0 = j;
            j1 = minint(j + 3, n - 1, _state);
         // Process submatrix
            for (ik = i0; ik <= i1; ik++) {
               for (jk = j0; jk <= j1; jk++) {
                  if (k == 0 || alpha == 0.0) {
                     v = 0.0;
                  } else {
                     v = 0.0;
                     v = ae_v_dotproduct(&a->xyR[ia][ja + ik], a->stride, &b->xyR[ib][jb + jk], b->stride, k);
                  }
                  if (beta == 0.0) {
                     c->xyR[ic + ik][jc + jk] = alpha * v;
                  } else {
                     c->xyR[ic + ik][jc + jk] = beta * c->xyR[ic + ik][jc + jk] + alpha * v;
                  }
               }
            }
         }
         j = j + 4;
      }
      i = i + 4;
   }
}

// RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
// with OpTypeA=1 and OpTypeB=1.
//
// Additional info:
// * this function requires that Alpha != 0 (assertion is thrown otherwise)
//
// Inputs:
//     M       -   matrix size, M>0
//     N       -   matrix size, N>0
//     K       -   matrix size, K>0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     B       -   matrix
//     IB      -   submatrix offset
//     JB      -   submatrix offset
//     Beta    -   coefficient
//     C       -   PREALLOCATED output matrix
//     IC      -   submatrix offset
//     JC      -   submatrix offset
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void rmatrixgemmk44v11(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double v;
   double v00;
   double v01;
   double v02;
   double v03;
   double v10;
   double v11;
   double v12;
   double v13;
   double v20;
   double v21;
   double v22;
   double v23;
   double v30;
   double v31;
   double v32;
   double v33;
   double a0;
   double a1;
   double a2;
   double a3;
   double b0;
   double b1;
   double b2;
   double b3;
   ae_int_t idxa0;
   ae_int_t idxa1;
   ae_int_t idxa2;
   ae_int_t idxa3;
   ae_int_t idxb0;
   ae_int_t idxb1;
   ae_int_t idxb2;
   ae_int_t idxb3;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t ik;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t jk;
   ae_int_t t;
   ae_int_t offsa;
   ae_int_t offsb;
   ae_assert(alpha != 0.0, "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
// if matrix size is zero
   if (m == 0 || n == 0) {
      return;
   }
// A'*B'
   i = 0;
   while (i < m) {
      j = 0;
      while (j < n) {
      // Choose between specialized 4x4 code and general code
         if (i + 4 <= m && j + 4 <= n) {
         // Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
         //
         // This submatrix is calculated as sum of K rank-1 products,
         // with operands cached in local variables in order to speed
         // up operations with arrays.
            idxa0 = ja + i + 0;
            idxa1 = ja + i + 1;
            idxa2 = ja + i + 2;
            idxa3 = ja + i + 3;
            offsa = ia;
            idxb0 = ib + j + 0;
            idxb1 = ib + j + 1;
            idxb2 = ib + j + 2;
            idxb3 = ib + j + 3;
            offsb = jb;
            v00 = 0.0;
            v01 = 0.0;
            v02 = 0.0;
            v03 = 0.0;
            v10 = 0.0;
            v11 = 0.0;
            v12 = 0.0;
            v13 = 0.0;
            v20 = 0.0;
            v21 = 0.0;
            v22 = 0.0;
            v23 = 0.0;
            v30 = 0.0;
            v31 = 0.0;
            v32 = 0.0;
            v33 = 0.0;
            for (t = 0; t < k; t++) {
               a0 = a->xyR[offsa][idxa0];
               a1 = a->xyR[offsa][idxa1];
               b0 = b->xyR[idxb0][offsb];
               b1 = b->xyR[idxb1][offsb];
               v00 = v00 + a0 * b0;
               v01 = v01 + a0 * b1;
               v10 = v10 + a1 * b0;
               v11 = v11 + a1 * b1;
               a2 = a->xyR[offsa][idxa2];
               a3 = a->xyR[offsa][idxa3];
               v20 = v20 + a2 * b0;
               v21 = v21 + a2 * b1;
               v30 = v30 + a3 * b0;
               v31 = v31 + a3 * b1;
               b2 = b->xyR[idxb2][offsb];
               b3 = b->xyR[idxb3][offsb];
               v22 = v22 + a2 * b2;
               v23 = v23 + a2 * b3;
               v32 = v32 + a3 * b2;
               v33 = v33 + a3 * b3;
               v02 = v02 + a0 * b2;
               v03 = v03 + a0 * b3;
               v12 = v12 + a1 * b2;
               v13 = v13 + a1 * b3;
               offsa = offsa + 1;
               offsb = offsb + 1;
            }
            if (beta == 0.0) {
               c->xyR[ic + i + 0][jc + j + 0] = alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = alpha * v33;
            } else {
               c->xyR[ic + i + 0][jc + j + 0] = beta * c->xyR[ic + i + 0][jc + j + 0] + alpha * v00;
               c->xyR[ic + i + 0][jc + j + 1] = beta * c->xyR[ic + i + 0][jc + j + 1] + alpha * v01;
               c->xyR[ic + i + 0][jc + j + 2] = beta * c->xyR[ic + i + 0][jc + j + 2] + alpha * v02;
               c->xyR[ic + i + 0][jc + j + 3] = beta * c->xyR[ic + i + 0][jc + j + 3] + alpha * v03;
               c->xyR[ic + i + 1][jc + j + 0] = beta * c->xyR[ic + i + 1][jc + j + 0] + alpha * v10;
               c->xyR[ic + i + 1][jc + j + 1] = beta * c->xyR[ic + i + 1][jc + j + 1] + alpha * v11;
               c->xyR[ic + i + 1][jc + j + 2] = beta * c->xyR[ic + i + 1][jc + j + 2] + alpha * v12;
               c->xyR[ic + i + 1][jc + j + 3] = beta * c->xyR[ic + i + 1][jc + j + 3] + alpha * v13;
               c->xyR[ic + i + 2][jc + j + 0] = beta * c->xyR[ic + i + 2][jc + j + 0] + alpha * v20;
               c->xyR[ic + i + 2][jc + j + 1] = beta * c->xyR[ic + i + 2][jc + j + 1] + alpha * v21;
               c->xyR[ic + i + 2][jc + j + 2] = beta * c->xyR[ic + i + 2][jc + j + 2] + alpha * v22;
               c->xyR[ic + i + 2][jc + j + 3] = beta * c->xyR[ic + i + 2][jc + j + 3] + alpha * v23;
               c->xyR[ic + i + 3][jc + j + 0] = beta * c->xyR[ic + i + 3][jc + j + 0] + alpha * v30;
               c->xyR[ic + i + 3][jc + j + 1] = beta * c->xyR[ic + i + 3][jc + j + 1] + alpha * v31;
               c->xyR[ic + i + 3][jc + j + 2] = beta * c->xyR[ic + i + 3][jc + j + 2] + alpha * v32;
               c->xyR[ic + i + 3][jc + j + 3] = beta * c->xyR[ic + i + 3][jc + j + 3] + alpha * v33;
            }
         } else {
         // Determine submatrix [I0..I1]x[J0..J1] to process
            i0 = i;
            i1 = minint(i + 3, m - 1, _state);
            j0 = j;
            j1 = minint(j + 3, n - 1, _state);
         // Process submatrix
            for (ik = i0; ik <= i1; ik++) {
               for (jk = j0; jk <= j1; jk++) {
                  if (k == 0 || alpha == 0.0) {
                     v = 0.0;
                  } else {
                     v = 0.0;
                     v = ae_v_dotproduct(&a->xyR[ia][ja + ik], a->stride, &b->xyR[ib + jk][jb], 1, k);
                  }
                  if (beta == 0.0) {
                     c->xyR[ic + ik][jc + jk] = alpha * v;
                  } else {
                     c->xyR[ic + ik][jc + jk] = beta * c->xyR[ic + ik][jc + jk] + alpha * v;
                  }
               }
            }
         }
         j = j + 4;
      }
      i = i + 4;
   }
}

#ifdef ALGLIB_NO_FAST_KERNELS
// Fast kernel (new version with AVX2/SSE2)
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
static bool ablasf_rgemm32basecase(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   bool result;
   result = false;
   return result;
}
#endif // defined ALGLIB_NO_FAST_KERNELS

// RMatrixGEMM kernel, basecase code for RMatrixGEMM.
//
// This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
// * C is MxN general matrix
// * op1(A) is MxK matrix
// * op2(B) is KxN matrix
// * "op" may be identity transformation, transposition
//
// Additional info:
// * multiplication result replaces C. If Beta=0, C elements are not used in
//   calculations (not multiplied by zero - just not referenced)
// * if Alpha=0, A is not used (not multiplied by zero - just not referenced)
// * if both Beta and Alpha are zero, C is filled by zeros.
//
// IMPORTANT:
//
// This function does NOT preallocate output matrix C, it MUST be preallocated
// by caller prior to calling this function. In case C does not have  enough
// space to store result, exception will be generated.
//
// Inputs:
//     M       -   matrix size, M>0
//     N       -   matrix size, N>0
//     K       -   matrix size, K>0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     OpTypeA -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//     B       -   matrix
//     IB      -   submatrix offset
//     JB      -   submatrix offset
//     OpTypeB -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//     Beta    -   coefficient
//     C       -   PREALLOCATED output matrix
//     IC      -   submatrix offset
//     JC      -   submatrix offset
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void rmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
// if matrix size is zero
   if (m == 0 || n == 0) {
      return;
   }
// Try optimized code
   if (ablasf_rgemm32basecase(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc, _state)) {
      return;
   }
// if K=0 or Alpha=0, then C=Beta*C
   if (k == 0 || alpha == 0.0) {
      if (beta != 1.0) {
         if (beta != 0.0) {
            for (i = 0; i < m; i++) {
               for (j = 0; j < n; j++) {
                  c->xyR[ic + i][jc + j] = beta * c->xyR[ic + i][jc + j];
               }
            }
         } else {
            for (i = 0; i < m; i++) {
               for (j = 0; j < n; j++) {
                  c->xyR[ic + i][jc + j] = 0.0;
               }
            }
         }
      }
      return;
   }
// Call specialized code.
//
// NOTE: specialized code was moved to separate function because of strange
//       issues with instructions cache on some systems; Having too long
//       functions significantly slows down internal loop of the algorithm.
   if (optypea == 0 && optypeb == 0) {
      rmatrixgemmk44v00(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
   }
   if (optypea == 0 && optypeb != 0) {
      rmatrixgemmk44v01(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
   }
   if (optypea != 0 && optypeb == 0) {
      rmatrixgemmk44v10(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
   }
   if (optypea != 0 && optypeb != 0) {
      rmatrixgemmk44v11(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
   }
}

// CMatrixGEMM kernel, basecase code for CMatrixGEMM.
//
// This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
// * C is MxN general matrix
// * op1(A) is MxK matrix
// * op2(B) is KxN matrix
// * "op" may be identity transformation, transposition, conjugate transposition
//
// Additional info:
// * multiplication result replaces C. If Beta=0, C elements are not used in
//   calculations (not multiplied by zero - just not referenced)
// * if Alpha=0, A is not used (not multiplied by zero - just not referenced)
// * if both Beta and Alpha are zero, C is filled by zeros.
//
// IMPORTANT:
//
// This function does NOT preallocate output matrix C, it MUST be preallocated
// by caller prior to calling this function. In case C does not have  enough
// space to store result, exception will be generated.
//
// Inputs:
//     M       -   matrix size, M>0
//     N       -   matrix size, N>0
//     K       -   matrix size, K>0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset
//     JA      -   submatrix offset
//     OpTypeA -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//                 * 2 - conjugate transposition
//     B       -   matrix
//     IB      -   submatrix offset
//     JB      -   submatrix offset
//     OpTypeB -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//                 * 2 - conjugate transposition
//     Beta    -   coefficient
//     C       -   PREALLOCATED output matrix
//     IC      -   submatrix offset
//     JC      -   submatrix offset
// ALGLIB Routine: Copyright 27.03.2013 by Sergey Bochkanov
void cmatrixgemmk(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   complex v;
   complex v00;
   complex v01;
   complex v10;
   complex v11;
   double v00x;
   double v00y;
   double v01x;
   double v01y;
   double v10x;
   double v10y;
   double v11x;
   double v11y;
   double a0x;
   double a0y;
   double a1x;
   double a1y;
   double b0x;
   double b0y;
   double b1x;
   double b1y;
   ae_int_t idxa0;
   ae_int_t idxa1;
   ae_int_t idxb0;
   ae_int_t idxb1;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t ik;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t jk;
   ae_int_t t;
   ae_int_t offsa;
   ae_int_t offsb;
// if matrix size is zero
   if (m == 0 || n == 0) {
      return;
   }
// Try optimized code
   if (cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc, _state)) {
      return;
   }
// if K=0 or Alpha=0, then C=Beta*C
   if (k == 0 || ae_c_eq_d(alpha, 0.0)) {
      if (ae_c_neq_d(beta, 1.0)) {
         if (ae_c_neq_d(beta, 0.0)) {
            for (i = 0; i < m; i++) {
               for (j = 0; j < n; j++) {
                  c->xyC[ic + i][jc + j] = ae_c_mul(beta, c->xyC[ic + i][jc + j]);
               }
            }
         } else {
            for (i = 0; i < m; i++) {
               for (j = 0; j < n; j++) {
                  c->xyC[ic + i][jc + j] = complex_from_i(0);
               }
            }
         }
      }
      return;
   }
// This phase is not really necessary, but compiler complains
// about "possibly uninitialized variables"
   a0x = 0.0;
   a0y = 0.0;
   a1x = 0.0;
   a1y = 0.0;
   b0x = 0.0;
   b0y = 0.0;
   b1x = 0.0;
   b1y = 0.0;
// General case
   i = 0;
   while (i < m) {
      j = 0;
      while (j < n) {
      // Choose between specialized 4x4 code and general code
         if (i + 2 <= m && j + 2 <= n) {
         // Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
         //
         // This submatrix is calculated as sum of K rank-1 products,
         // with operands cached in local variables in order to speed
         // up operations with arrays.
            v00x = 0.0;
            v00y = 0.0;
            v01x = 0.0;
            v01y = 0.0;
            v10x = 0.0;
            v10y = 0.0;
            v11x = 0.0;
            v11y = 0.0;
            if (optypea == 0) {
               idxa0 = ia + i + 0;
               idxa1 = ia + i + 1;
               offsa = ja;
            } else {
               idxa0 = ja + i + 0;
               idxa1 = ja + i + 1;
               offsa = ia;
            }
            if (optypeb == 0) {
               idxb0 = jb + j + 0;
               idxb1 = jb + j + 1;
               offsb = ib;
            } else {
               idxb0 = ib + j + 0;
               idxb1 = ib + j + 1;
               offsb = jb;
            }
            for (t = 0; t < k; t++) {
               if (optypea == 0) {
                  a0x = a->xyC[idxa0][offsa].x;
                  a0y = a->xyC[idxa0][offsa].y;
                  a1x = a->xyC[idxa1][offsa].x;
                  a1y = a->xyC[idxa1][offsa].y;
               }
               if (optypea == 1) {
                  a0x = a->xyC[offsa][idxa0].x;
                  a0y = a->xyC[offsa][idxa0].y;
                  a1x = a->xyC[offsa][idxa1].x;
                  a1y = a->xyC[offsa][idxa1].y;
               }
               if (optypea == 2) {
                  a0x = a->xyC[offsa][idxa0].x;
                  a0y = -a->xyC[offsa][idxa0].y;
                  a1x = a->xyC[offsa][idxa1].x;
                  a1y = -a->xyC[offsa][idxa1].y;
               }
               if (optypeb == 0) {
                  b0x = b->xyC[offsb][idxb0].x;
                  b0y = b->xyC[offsb][idxb0].y;
                  b1x = b->xyC[offsb][idxb1].x;
                  b1y = b->xyC[offsb][idxb1].y;
               }
               if (optypeb == 1) {
                  b0x = b->xyC[idxb0][offsb].x;
                  b0y = b->xyC[idxb0][offsb].y;
                  b1x = b->xyC[idxb1][offsb].x;
                  b1y = b->xyC[idxb1][offsb].y;
               }
               if (optypeb == 2) {
                  b0x = b->xyC[idxb0][offsb].x;
                  b0y = -b->xyC[idxb0][offsb].y;
                  b1x = b->xyC[idxb1][offsb].x;
                  b1y = -b->xyC[idxb1][offsb].y;
               }
               v00x = v00x + a0x * b0x - a0y * b0y;
               v00y = v00y + a0x * b0y + a0y * b0x;
               v01x = v01x + a0x * b1x - a0y * b1y;
               v01y = v01y + a0x * b1y + a0y * b1x;
               v10x = v10x + a1x * b0x - a1y * b0y;
               v10y = v10y + a1x * b0y + a1y * b0x;
               v11x = v11x + a1x * b1x - a1y * b1y;
               v11y = v11y + a1x * b1y + a1y * b1x;
               offsa = offsa + 1;
               offsb = offsb + 1;
            }
            v00 = complex_from_d(v00x, v00y);
            v10 = complex_from_d(v10x, v10y);
            v01 = complex_from_d(v01x, v01y);
            v11 = complex_from_d(v11x, v11y);
            if (ae_c_eq_d(beta, 0.0)) {
               c->xyC[ic + i + 0][jc + j + 0] = ae_c_mul(alpha, v00);
               c->xyC[ic + i + 0][jc + j + 1] = ae_c_mul(alpha, v01);
               c->xyC[ic + i + 1][jc + j + 0] = ae_c_mul(alpha, v10);
               c->xyC[ic + i + 1][jc + j + 1] = ae_c_mul(alpha, v11);
            } else {
               c->xyC[ic + i + 0][jc + j + 0] = ae_c_add(ae_c_mul(beta, c->xyC[ic + i + 0][jc + j + 0]), ae_c_mul(alpha, v00));
               c->xyC[ic + i + 0][jc + j + 1] = ae_c_add(ae_c_mul(beta, c->xyC[ic + i + 0][jc + j + 1]), ae_c_mul(alpha, v01));
               c->xyC[ic + i + 1][jc + j + 0] = ae_c_add(ae_c_mul(beta, c->xyC[ic + i + 1][jc + j + 0]), ae_c_mul(alpha, v10));
               c->xyC[ic + i + 1][jc + j + 1] = ae_c_add(ae_c_mul(beta, c->xyC[ic + i + 1][jc + j + 1]), ae_c_mul(alpha, v11));
            }
         } else {
         // Determine submatrix [I0..I1]x[J0..J1] to process
            i0 = i;
            i1 = minint(i + 1, m - 1, _state);
            j0 = j;
            j1 = minint(j + 1, n - 1, _state);
         // Process submatrix
            for (ik = i0; ik <= i1; ik++) {
               for (jk = j0; jk <= j1; jk++) {
                  if (k == 0 || ae_c_eq_d(alpha, 0.0)) {
                     v = complex_from_i(0);
                  } else {
                     v = complex_from_d(0.0);
                     if (optypea == 0 && optypeb == 0) {
                        v = ae_v_cdotproduct(&a->xyC[ia + ik][ja], 1, "N", &b->xyC[ib][jb + jk], b->stride, "N", k);
                     }
                     if (optypea == 0 && optypeb == 1) {
                        v = ae_v_cdotproduct(&a->xyC[ia + ik][ja], 1, "N", &b->xyC[ib + jk][jb], 1, "N", k);
                     }
                     if (optypea == 0 && optypeb == 2) {
                        v = ae_v_cdotproduct(&a->xyC[ia + ik][ja], 1, "N", &b->xyC[ib + jk][jb], 1, "Conj", k);
                     }
                     if (optypea == 1 && optypeb == 0) {
                        v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "N", &b->xyC[ib][jb + jk], b->stride, "N", k);
                     }
                     if (optypea == 1 && optypeb == 1) {
                        v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "N", &b->xyC[ib + jk][jb], 1, "N", k);
                     }
                     if (optypea == 1 && optypeb == 2) {
                        v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "N", &b->xyC[ib + jk][jb], 1, "Conj", k);
                     }
                     if (optypea == 2 && optypeb == 0) {
                        v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "Conj", &b->xyC[ib][jb + jk], b->stride, "N", k);
                     }
                     if (optypea == 2 && optypeb == 1) {
                        v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "Conj", &b->xyC[ib + jk][jb], 1, "N", k);
                     }
                     if (optypea == 2 && optypeb == 2) {
                        v = ae_v_cdotproduct(&a->xyC[ia][ja + ik], a->stride, "Conj", &b->xyC[ib + jk][jb], 1, "Conj", k);
                     }
                  }
                  if (ae_c_eq_d(beta, 0.0)) {
                     c->xyC[ic + ik][jc + jk] = ae_c_mul(alpha, v);
                  } else {
                     c->xyC[ic + ik][jc + jk] = ae_c_add(ae_c_mul(beta, c->xyC[ic + ik][jc + jk]), ae_c_mul(alpha, v));
                  }
               }
            }
         }
         j = j + 2;
      }
      i = i + 2;
   }
}
} // end of namespace alglib_impl

// === HBLAS Package ===
namespace alglib_impl {
void hermitianmatrixvectormultiply(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, complex alpha, CVector *y, ae_state *_state) {
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
//  L is strictly lower triangular (main diagonal is zero)
//  D is diagonal
//  U is strictly upper triangular (main diagonal is zero)
//
// A*x = L*x + D*x + U*x
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

void hermitianrank2update(CMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, CVector *x, CVector *y, CVector *t, complex alpha, ae_state *_state) {
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
         v = ae_c_mul(conj(alpha, _state), y->xC[i + 1 - i1]);
         ae_v_caddc(&t->xC[tp1], 1, &x->xC[tp1], 1, "Conj", tp2 - tp1 + 1, v);
         ae_v_cadd(&a->xyC[i][i], 1, &t->xC[tp1], 1, "N", i2 - i + 1);
      }
   } else {
      for (i = i1; i <= i2; i++) {
         tp1 = 1;
         tp2 = i + 1 - i1;
         v = ae_c_mul(alpha, x->xC[i + 1 - i1]);
         ae_v_cmovec(&t->xC[tp1], 1, &y->xC[tp1], 1, "Conj", tp2 - tp1 + 1, v);
         v = ae_c_mul(conj(alpha, _state), y->xC[i + 1 - i1]);
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
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void complexgeneratereflection(CVector *x, ae_int_t n, complex *tau, ae_state *_state) {
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
   tau->x = 0;
   tau->y = 0;
   if (n <= 0) {
      *tau = complex_from_i(0);
      return;
   }
// Scale if needed (to avoid overflow/underflow during intermediate
// calculations).
   mx = 0.0;
   for (j = 1; j <= n; j++) {
      mx = maxreal(abscomplex(x->xC[j], _state), mx, _state);
   }
   s = 1.0;
   if (mx != 0.0) {
      if (mx < 1.0) {
         s = ae_sqrt(minrealnumber, _state);
         v = complex_from_d(1 / s);
         ae_v_cmulc(&x->xC[1], 1, n, v);
      } else {
         s = ae_sqrt(maxrealnumber, _state);
         v = complex_from_d(1 / s);
         ae_v_cmulc(&x->xC[1], 1, n, v);
      }
   }
// calculate
   alpha = x->xC[1];
   mx = 0.0;
   for (j = 2; j <= n; j++) {
      mx = maxreal(abscomplex(x->xC[j], _state), mx, _state);
   }
   xnorm = 0.0;
   if (mx != 0.0) {
      for (j = 2; j <= n; j++) {
         t = ae_c_div_d(x->xC[j], mx);
         xnorm = xnorm + ae_c_mul(t, conj(t, _state)).x;
      }
      xnorm = ae_sqrt(xnorm, _state) * mx;
   }
   alphr = alpha.x;
   alphi = alpha.y;
   if (xnorm == 0.0 && alphi == 0.0) {
      *tau = complex_from_i(0);
      x->xC[1] = ae_c_mul_d(x->xC[1], s);
      return;
   }
   mx = maxreal(ae_fabs(alphr, _state), ae_fabs(alphi, _state), _state);
   mx = maxreal(mx, ae_fabs(xnorm, _state), _state);
   beta = -mx * ae_sqrt(sqr(alphr / mx, _state) + sqr(alphi / mx, _state) + sqr(xnorm / mx, _state), _state);
   if (alphr < 0.0) {
      beta = -beta;
   }
   tau->x = (beta - alphr) / beta;
   tau->y = -alphi / beta;
   alpha = ae_c_d_div(1, ae_c_sub_d(alpha, beta));
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
//                 If N1>N2 or M1>M2, C is not modified.
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void complexapplyreflectionfromtheleft(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work, ae_state *_state) {
   complex t;
   ae_int_t i;
   if ((ae_c_eq_d(tau, 0.0) || n1 > n2) || m1 > m2) {
      return;
   }
// w := C^T * conj(v)
   for (i = n1; i <= n2; i++) {
      work->xC[i] = complex_from_i(0);
   }
   for (i = m1; i <= m2; i++) {
      t = conj(v->xC[i + 1 - m1], _state);
      ae_v_caddc(&work->xC[n1], 1, &c->xyC[i][n1], 1, "N", n2 - n1 + 1, t);
   }
// C := C - tau * v * w^T
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
//                 If N1>N2 or M1>M2, C is not modified.
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void complexapplyreflectionfromtheright(CMatrix *c, complex tau, CVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, CVector *work, ae_state *_state) {
   complex t;
   ae_int_t i;
   ae_int_t vm;
   if ((ae_c_eq_d(tau, 0.0) || n1 > n2) || m1 > m2) {
      return;
   }
// w := C * v
   vm = n2 - n1 + 1;
   for (i = m1; i <= m2; i++) {
      t = ae_v_cdotproduct(&c->xyC[i][n1], 1, "N", &v->xC[1], 1, "N", n2 - n1 + 1);
      work->xC[i] = t;
   }
// C := C - w * conj(v^T)
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
void symmetricmatrixvectormultiply(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, double alpha, RVector *y, ae_state *_state) {
   ae_int_t i;
   ae_int_t ba1;
   ae_int_t ba2;
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
//  L is strictly lower triangular (main diagonal is zero)
//  D is diagonal
//  U is strictly upper triangular (main diagonal is zero)
//
// A*x = L*x + D*x + U*x
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
         ba2 = i2;
         ae_v_addd(&y->xR[by1], 1, &a->xyR[i][ba1], 1, by2 - by1 + 1, v);
      // Add U*x to the result
         bx1 = i - i1 + 2;
         bx2 = n;
         ba1 = i + 1;
         ba2 = i2;
         v = ae_v_dotproduct(&x->xR[bx1], 1, &a->xyR[i][ba1], 1, bx2 - bx1 + 1);
         y->xR[i - i1 + 1] = y->xR[i - i1 + 1] + v;
      }
   } else {
      for (i = i1 + 1; i <= i2; i++) {
      // Add L*x to the result
         bx1 = 1;
         bx2 = i - i1;
         ba1 = i1;
         ba2 = i - 1;
         v = ae_v_dotproduct(&x->xR[bx1], 1, &a->xyR[i][ba1], 1, bx2 - bx1 + 1);
         y->xR[i - i1 + 1] = y->xR[i - i1 + 1] + v;
      // Add U*x to the result
         v = x->xR[i - i1 + 1];
         by1 = 1;
         by2 = i - i1;
         ba1 = i1;
         ba2 = i - 1;
         ae_v_addd(&y->xR[by1], 1, &a->xyR[i][ba1], 1, by2 - by1 + 1, v);
      }
   }
   ae_v_muld(&y->xR[1], 1, n, alpha);
   touchint(&ba2, _state);
}

void symmetricrank2update(RMatrix *a, bool isupper, ae_int_t i1, ae_int_t i2, RVector *x, RVector *y, RVector *t, double alpha, ae_state *_state) {
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
bool rmatrixgermkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double alpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state) {
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
bool rmatrixrank1mkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv, ae_state *_state) {
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
bool cmatrixrank1mkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv, ae_state *_state) {
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
bool rmatrixmvmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, RVector *y, ae_int_t iy, ae_state *_state) {
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
bool cmatrixmvmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CVector *x, ae_int_t ix, CVector *y, ae_int_t iy, ae_state *_state) {
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
bool rmatrixgemvmkl(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state) {
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
bool rmatrixtrsvmkl(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixtrsvmkl(n, a, ia, ja, isupper, isunit, optype, x, ix);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2017 by Sergey Bochkanov
bool rmatrixsymvmkl(ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state) {
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
bool rmatrixsyrkmkl(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixsyrkmkl(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2013 by Sergey Bochkanov
bool cmatrixherkmkl(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixherkmkl(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 01.10.2013 by Sergey Bochkanov
bool rmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixgemmmkl(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool cmatrixgemmmkl(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixgemmmkl(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool rmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool cmatrixlefttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool rmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_rmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}

// MKL-based kernel
// ALGLIB Routine: Copyright 16.10.2014 by Sergey Bochkanov
bool cmatrixrighttrsmmkl(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_MKL
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_cmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
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
bool spdmatrixcholeskymkl(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, bool *cholresult, ae_state *_state) {
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
bool rmatrixplumkl(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, ae_state *_state) {
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
bool rmatrixbdmkl(RMatrix *a, ae_int_t m, ae_int_t n, RVector *d, RVector *e, RVector *tauq, RVector *taup, ae_state *_state) {
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
bool rmatrixbdmultiplybymkl(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, RVector *taup, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool byq, bool fromtheright, bool dotranspose, ae_state *_state) {
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
bool rmatrixhessenbergmkl(RMatrix *a, ae_int_t n, RVector *tau, ae_state *_state) {
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
bool rmatrixhessenbergunpackqmkl(RMatrix *a, ae_int_t n, RVector *tau, RMatrix *q, ae_state *_state) {
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
//       length(E)=length(Tau)=N-1 (or larger)
//       length(D)=N (or larger)
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool smatrixtdmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RVector *d, RVector *e, ae_state *_state) {
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
//       length(E)=length(Tau)=N-1 (or larger)
//       length(D)=N (or larger)
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool hmatrixtdmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, RVector *d, RVector *e, ae_state *_state) {
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
bool smatrixtdunpackqmkl(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RMatrix *q, ae_state *_state) {
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
bool hmatrixtdunpackqmkl(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, CMatrix *q, ae_state *_state) {
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
// D and E are pre-allocated arrays with length N (both of them!). On output,
// D constraints singular values, and E is destroyed.
//
// SVDResult is modified if and only if MKL is present.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixbdsvdmkl(RVector *d, RVector *e, ae_int_t n, bool isupper, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt, bool *svdresult, ae_state *_state) {
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
// WR and WI are pre-allocated arrays with length N.
// Z is pre-allocated array[N,N].
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixinternalschurdecompositionmkl(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info, ae_state *_state) {
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
// NOTE: this function does NOT support HOWMNY=3!!!!
//
// VL and VR are pre-allocated arrays with length N*N, if required. If particalar
// variables is not required, it can be dummy (empty) array.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool rmatrixinternaltrevcmkl(RMatrix *t, ae_int_t n, ae_int_t side, ae_int_t howmny, RMatrix *vl, RMatrix *vr, ae_int_t *m, ae_int_t *info, ae_state *_state) {
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
// D and E are pre-allocated arrays with length N (both of them!). On output,
// D constraints eigenvalues, and E is destroyed.
//
// Z is preallocated array[N,N] for ZNeeded != 0; ignored for ZNeeded=0.
//
// EVDResult is modified if and only if MKL is present.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool smatrixtdevdmkl(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z, bool *evdresult, ae_state *_state) {
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
// D and E are pre-allocated arrays with length N (both of them!). On output,
// D constraints eigenvalues, and E is destroyed.
//
// Z is preallocated array[N,N] for ZNeeded != 0; ignored for ZNeeded=0.
//
// EVDResult is modified if and only if MKL is present.
// ALGLIB Routine: Copyright 20.10.2014 by Sergey Bochkanov
bool sparsegemvcrsmkl(ae_int_t opa, ae_int_t arows, ae_int_t acols, double alpha, RVector *vals, ZVector *cidx, ZVector *ridx, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy, ae_state *_state) {
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
ae_int_t getrdfserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 1;
   return result;
}

ae_int_t getkdtreeserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 2;
   return result;
}

ae_int_t getmlpserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 3;
   return result;
}

ae_int_t getmlpeserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 4;
   return result;
}

ae_int_t getrbfserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 5;
   return result;
}

ae_int_t getspline2dserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 6;
   return result;
}

ae_int_t getidwserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 7;
   return result;
}

ae_int_t getsparsematrixserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 8;
   return result;
}

ae_int_t getknnserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 108;
   return result;
}

ae_int_t getlptestserializationcode(ae_state *_state) {
   ae_int_t result;
   result = 200;
   return result;
}
} // end of namespace alglib_impl

// === TSORT Package ===
// Depends on: APSERV
namespace alglib_impl {
// Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
// applies same permutations to B.
// ALGLIB: Copyright 06.09.2010 by Sergey Bochkanov
static void tsort_tagsortfastirec(RVector *a, ZVector *b, RVector *bufa, ZVector *bufb, ae_int_t i1, ae_int_t i2, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t cntless;
   ae_int_t cnteq;
   ae_int_t cntgreater;
   double tmpr;
   ae_int_t tmpi;
   double v0;
   double v1;
   double v2;
   double vp;
// Fast exit
   if (i2 <= i1) {
      return;
   }
// Non-recursive sort for small arrays
   if (i2 - i1 <= 16) {
      for (j = i1 + 1; j <= i2; j++) {
      // Search elements [I1..J-1] for place to insert Jth element.
      //
      // This code stops immediately if we can leave A[J] at J-th position
      // (all elements have same value of A[J] larger than any of them)
         tmpr = a->xR[j];
         tmpi = j;
         for (k = j - 1; k >= i1; k--) {
            if (a->xR[k] <= tmpr) {
               break;
            }
            tmpi = k;
         }
         k = tmpi;
      // Insert Jth element into Kth position
         if (k != j) {
            tmpr = a->xR[j];
            tmpi = b->xZ[j];
            for (i = j - 1; i >= k; i--) {
               a->xR[i + 1] = a->xR[i];
               b->xZ[i + 1] = b->xZ[i];
            }
            a->xR[k] = tmpr;
            b->xZ[k] = tmpi;
         }
      }
      return;
   }
// Quicksort: choose pivot
// Here we assume that I2-I1 >= 2
   v0 = a->xR[i1];
   v1 = a->xR[i1 + (i2 - i1) / 2];
   v2 = a->xR[i2];
   if (v0 > v1) {
      tmpr = v1;
      v1 = v0;
      v0 = tmpr;
   }
   if (v1 > v2) {
      tmpr = v2;
      v2 = v1;
      v1 = tmpr;
   }
   if (v0 > v1) {
      tmpr = v1;
      v1 = v0;
      v0 = tmpr;
   }
   vp = v1;
// now pass through A/B and:
// * move elements that are LESS than VP to the left of A/B
// * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
// * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
// * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
// * move elements from the left of BufA/BufB to the end of A/B
   cntless = 0;
   cnteq = 0;
   cntgreater = 0;
   for (i = i1; i <= i2; i++) {
      v0 = a->xR[i];
      if (v0 < vp) {
      // LESS
         k = i1 + cntless;
         if (i != k) {
            a->xR[k] = v0;
            b->xZ[k] = b->xZ[i];
         }
         cntless = cntless + 1;
         continue;
      }
      if (v0 == vp) {
      // EQUAL
         k = i2 - cnteq;
         bufa->xR[k] = v0;
         bufb->xZ[k] = b->xZ[i];
         cnteq = cnteq + 1;
         continue;
      }
   // GREATER
      k = i1 + cntgreater;
      bufa->xR[k] = v0;
      bufb->xZ[k] = b->xZ[i];
      cntgreater = cntgreater + 1;
   }
   for (i = 0; i < cnteq; i++) {
      j = i1 + cntless + cnteq - 1 - i;
      k = i2 + i - (cnteq - 1);
      a->xR[j] = bufa->xR[k];
      b->xZ[j] = bufb->xZ[k];
   }
   for (i = 0; i < cntgreater; i++) {
      j = i1 + cntless + cnteq + i;
      k = i1 + i;
      a->xR[j] = bufa->xR[k];
      b->xZ[j] = bufb->xZ[k];
   }
// Sort left and right parts of the array (ignoring middle part)
   tsort_tagsortfastirec(a, b, bufa, bufb, i1, i1 + cntless - 1, _state);
   tsort_tagsortfastirec(a, b, bufa, bufb, i1 + cntless + cnteq, i2, _state);
}

// Internal TagSortFastR: sorts A[I1...I2] (both bounds are included),
// applies same permutations to B.
// ALGLIB: Copyright 06.09.2010 by Sergey Bochkanov
static void tsort_tagsortfastrrec(RVector *a, RVector *b, RVector *bufa, RVector *bufb, ae_int_t i1, ae_int_t i2, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double tmpr;
   double tmpr2;
   ae_int_t tmpi;
   ae_int_t cntless;
   ae_int_t cnteq;
   ae_int_t cntgreater;
   double v0;
   double v1;
   double v2;
   double vp;
// Fast exit
   if (i2 <= i1) {
      return;
   }
// Non-recursive sort for small arrays
   if (i2 - i1 <= 16) {
      for (j = i1 + 1; j <= i2; j++) {
      // Search elements [I1..J-1] for place to insert Jth element.
      //
      // This code stops immediatly if we can leave A[J] at J-th position
      // (all elements have same value of A[J] larger than any of them)
         tmpr = a->xR[j];
         tmpi = j;
         for (k = j - 1; k >= i1; k--) {
            if (a->xR[k] <= tmpr) {
               break;
            }
            tmpi = k;
         }
         k = tmpi;
      // Insert Jth element into Kth position
         if (k != j) {
            tmpr = a->xR[j];
            tmpr2 = b->xR[j];
            for (i = j - 1; i >= k; i--) {
               a->xR[i + 1] = a->xR[i];
               b->xR[i + 1] = b->xR[i];
            }
            a->xR[k] = tmpr;
            b->xR[k] = tmpr2;
         }
      }
      return;
   }
// Quicksort: choose pivot
// Here we assume that I2-I1 >= 16
   v0 = a->xR[i1];
   v1 = a->xR[i1 + (i2 - i1) / 2];
   v2 = a->xR[i2];
   if (v0 > v1) {
      tmpr = v1;
      v1 = v0;
      v0 = tmpr;
   }
   if (v1 > v2) {
      tmpr = v2;
      v2 = v1;
      v1 = tmpr;
   }
   if (v0 > v1) {
      tmpr = v1;
      v1 = v0;
      v0 = tmpr;
   }
   vp = v1;
// now pass through A/B and:
// * move elements that are LESS than VP to the left of A/B
// * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
// * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
// * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
// * move elements from the left of BufA/BufB to the end of A/B
   cntless = 0;
   cnteq = 0;
   cntgreater = 0;
   for (i = i1; i <= i2; i++) {
      v0 = a->xR[i];
      if (v0 < vp) {
      // LESS
         k = i1 + cntless;
         if (i != k) {
            a->xR[k] = v0;
            b->xR[k] = b->xR[i];
         }
         cntless = cntless + 1;
         continue;
      }
      if (v0 == vp) {
      // EQUAL
         k = i2 - cnteq;
         bufa->xR[k] = v0;
         bufb->xR[k] = b->xR[i];
         cnteq = cnteq + 1;
         continue;
      }
   // GREATER
      k = i1 + cntgreater;
      bufa->xR[k] = v0;
      bufb->xR[k] = b->xR[i];
      cntgreater = cntgreater + 1;
   }
   for (i = 0; i < cnteq; i++) {
      j = i1 + cntless + cnteq - 1 - i;
      k = i2 + i - (cnteq - 1);
      a->xR[j] = bufa->xR[k];
      b->xR[j] = bufb->xR[k];
   }
   for (i = 0; i < cntgreater; i++) {
      j = i1 + cntless + cnteq + i;
      k = i1 + i;
      a->xR[j] = bufa->xR[k];
      b->xR[j] = bufb->xR[k];
   }
// Sort left and right parts of the array (ignoring middle part)
   tsort_tagsortfastrrec(a, b, bufa, bufb, i1, i1 + cntless - 1, _state);
   tsort_tagsortfastrrec(a, b, bufa, bufb, i1 + cntless + cnteq, i2, _state);
}

// Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
// applies same permutations to B.
// ALGLIB: Copyright 06.09.2010 by Sergey Bochkanov
static void tsort_tagsortfastrec(RVector *a, RVector *bufa, ae_int_t i1, ae_int_t i2, ae_state *_state) {
   ae_int_t cntless;
   ae_int_t cnteq;
   ae_int_t cntgreater;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double tmpr;
   ae_int_t tmpi;
   double v0;
   double v1;
   double v2;
   double vp;
// Fast exit
   if (i2 <= i1) {
      return;
   }
// Non-recursive sort for small arrays
   if (i2 - i1 <= 16) {
      for (j = i1 + 1; j <= i2; j++) {
      // Search elements [I1..J-1] for place to insert Jth element.
      //
      // This code stops immediatly if we can leave A[J] at J-th position
      // (all elements have same value of A[J] larger than any of them)
         tmpr = a->xR[j];
         tmpi = j;
         for (k = j - 1; k >= i1; k--) {
            if (a->xR[k] <= tmpr) {
               break;
            }
            tmpi = k;
         }
         k = tmpi;
      // Insert Jth element into Kth position
         if (k != j) {
            tmpr = a->xR[j];
            for (i = j - 1; i >= k; i--) {
               a->xR[i + 1] = a->xR[i];
            }
            a->xR[k] = tmpr;
         }
      }
      return;
   }
// Quicksort: choose pivot
// Here we assume that I2-I1 >= 16
   v0 = a->xR[i1];
   v1 = a->xR[i1 + (i2 - i1) / 2];
   v2 = a->xR[i2];
   if (v0 > v1) {
      tmpr = v1;
      v1 = v0;
      v0 = tmpr;
   }
   if (v1 > v2) {
      tmpr = v2;
      v2 = v1;
      v1 = tmpr;
   }
   if (v0 > v1) {
      tmpr = v1;
      v1 = v0;
      v0 = tmpr;
   }
   vp = v1;
// now pass through A/B and:
// * move elements that are LESS than VP to the left of A/B
// * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
// * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
// * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
// * move elements from the left of BufA/BufB to the end of A/B
   cntless = 0;
   cnteq = 0;
   cntgreater = 0;
   for (i = i1; i <= i2; i++) {
      v0 = a->xR[i];
      if (v0 < vp) {
      // LESS
         k = i1 + cntless;
         if (i != k) {
            a->xR[k] = v0;
         }
         cntless = cntless + 1;
         continue;
      }
      if (v0 == vp) {
      // EQUAL
         k = i2 - cnteq;
         bufa->xR[k] = v0;
         cnteq = cnteq + 1;
         continue;
      }
   // GREATER
      k = i1 + cntgreater;
      bufa->xR[k] = v0;
      cntgreater = cntgreater + 1;
   }
   for (i = 0; i < cnteq; i++) {
      j = i1 + cntless + cnteq - 1 - i;
      k = i2 + i - (cnteq - 1);
      a->xR[j] = bufa->xR[k];
   }
   for (i = 0; i < cntgreater; i++) {
      j = i1 + cntless + cnteq + i;
      k = i1 + i;
      a->xR[j] = bufa->xR[k];
   }
// Sort left and right parts of the array (ignoring middle part)
   tsort_tagsortfastrec(a, bufa, i1, i1 + cntless - 1, _state);
   tsort_tagsortfastrec(a, bufa, i1 + cntless + cnteq, i2, _state);
}

// Same as TagSort, but optimized for real keys and integer labels.
//
// A is sorted, and same permutations are applied to B.
//
// NOTES:
// 1.  this function assumes that A[] is finite; it doesn't checks that
//     condition. All other conditions (size of input arrays, etc.) are not
//     checked too.
// 2.  this function uses two buffers, BufA and BufB, each is N elements large.
//     They may be preallocated (which will save some time) or not, in which
//     case function will automatically allocate memory.
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortfasti(RVector *a, ZVector *b, RVector *bufa, ZVector *bufb, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   bool isascending;
   bool isdescending;
   double tmpr;
   ae_int_t tmpi;
// Special case
   if (n <= 1) {
      return;
   }
// Test for already sorted set
   isascending = true;
   isdescending = true;
   for (i = 1; i < n; i++) {
      isascending = isascending && a->xR[i] >= a->xR[i - 1];
      isdescending = isdescending && a->xR[i] <= a->xR[i - 1];
   }
   if (isascending) {
      return;
   }
   if (isdescending) {
      for (i = 0; i < n; i++) {
         j = n - 1 - i;
         if (j <= i) {
            break;
         }
         tmpr = a->xR[i];
         a->xR[i] = a->xR[j];
         a->xR[j] = tmpr;
         tmpi = b->xZ[i];
         b->xZ[i] = b->xZ[j];
         b->xZ[j] = tmpi;
      }
      return;
   }
// General case
   if (bufa->cnt < n) {
      ae_vector_set_length(bufa, n, _state);
   }
   if (bufb->cnt < n) {
      ae_vector_set_length(bufb, n, _state);
   }
   tsort_tagsortfastirec(a, b, bufa, bufb, 0, n - 1, _state);
}

// Same as TagSort, but optimized for real keys and real labels.
//
// A is sorted, and same permutations are applied to B.
//
// NOTES:
// 1.  this function assumes that A[] is finite; it doesn't checks that
//     condition. All other conditions (size of input arrays, etc.) are not
//     checked too.
// 2.  this function uses two buffers, BufA and BufB, each is N elements large.
//     They may be preallocated (which will save some time) or not, in which
//     case function will automatically allocate memory.
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortfastr(RVector *a, RVector *b, RVector *bufa, RVector *bufb, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   bool isascending;
   bool isdescending;
   double tmpr;
// Special case
   if (n <= 1) {
      return;
   }
// Test for already sorted set
   isascending = true;
   isdescending = true;
   for (i = 1; i < n; i++) {
      isascending = isascending && a->xR[i] >= a->xR[i - 1];
      isdescending = isdescending && a->xR[i] <= a->xR[i - 1];
   }
   if (isascending) {
      return;
   }
   if (isdescending) {
      for (i = 0; i < n; i++) {
         j = n - 1 - i;
         if (j <= i) {
            break;
         }
         tmpr = a->xR[i];
         a->xR[i] = a->xR[j];
         a->xR[j] = tmpr;
         tmpr = b->xR[i];
         b->xR[i] = b->xR[j];
         b->xR[j] = tmpr;
      }
      return;
   }
// General case
   if (bufa->cnt < n) {
      ae_vector_set_length(bufa, n, _state);
   }
   if (bufb->cnt < n) {
      ae_vector_set_length(bufb, n, _state);
   }
   tsort_tagsortfastrrec(a, b, bufa, bufb, 0, n - 1, _state);
}

// Same as TagSort, but optimized for real keys without labels.
//
// A is sorted, and that's all.
//
// NOTES:
// 1.  this function assumes that A[] is finite; it doesn't checks that
//     condition. All other conditions (size of input arrays, etc.) are not
//     checked too.
// 2.  this function uses buffer, BufA, which is N elements large. It may be
//     preallocated (which will save some time) or not, in which case
//     function will automatically allocate memory.
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortfast(RVector *a, RVector *bufa, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   bool isascending;
   bool isdescending;
   double tmpr;
// Special case
   if (n <= 1) {
      return;
   }
// Test for already sorted set
   isascending = true;
   isdescending = true;
   for (i = 1; i < n; i++) {
      isascending = isascending && a->xR[i] >= a->xR[i - 1];
      isdescending = isdescending && a->xR[i] <= a->xR[i - 1];
   }
   if (isascending) {
      return;
   }
   if (isdescending) {
      for (i = 0; i < n; i++) {
         j = n - 1 - i;
         if (j <= i) {
            break;
         }
         tmpr = a->xR[i];
         a->xR[i] = a->xR[j];
         a->xR[j] = tmpr;
      }
      return;
   }
// General case
   if (bufa->cnt < n) {
      ae_vector_set_length(bufa, n, _state);
   }
   tsort_tagsortfastrec(a, bufa, 0, n - 1, _state);
}

// Sorting function optimized for integer keys and real labels, can be used
// to sort middle of the array
//
// A is sorted, and same permutations are applied to B.
//
// NOTES:
//     this function assumes that A[] is finite; it doesn't checks that
//     condition. All other conditions (size of input arrays, etc.) are not
//     checked too.
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortmiddleir(ZVector *a, RVector *b, ae_int_t offset, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t t;
   ae_int_t tmp;
   double tmpr;
   ae_int_t p0;
   ae_int_t p1;
   ae_int_t at;
   ae_int_t ak;
   ae_int_t ak1;
   double bt;
// Special cases
   if (n <= 1) {
      return;
   }
// General case, N>1: sort, update B
   for (i = 2; i <= n; i++) {
      t = i;
      while (t != 1) {
         k = t / 2;
         p0 = offset + k - 1;
         p1 = offset + t - 1;
         ak = a->xZ[p0];
         at = a->xZ[p1];
         if (ak >= at) {
            break;
         }
         a->xZ[p0] = at;
         a->xZ[p1] = ak;
         tmpr = b->xR[p0];
         b->xR[p0] = b->xR[p1];
         b->xR[p1] = tmpr;
         t = k;
      }
   }
   for (i = n - 1; i >= 1; i--) {
      p0 = offset + 0;
      p1 = offset + i;
      tmp = a->xZ[p1];
      a->xZ[p1] = a->xZ[p0];
      a->xZ[p0] = tmp;
      at = tmp;
      tmpr = b->xR[p1];
      b->xR[p1] = b->xR[p0];
      b->xR[p0] = tmpr;
      bt = tmpr;
      t = 0;
      while (true) {
         k = 2 * t + 1;
         if (k + 1 > i) {
            break;
         }
         p0 = offset + t;
         p1 = offset + k;
         ak = a->xZ[p1];
         if (k + 1 < i) {
            ak1 = a->xZ[p1 + 1];
            if (ak1 > ak) {
               ak = ak1;
               p1 = p1 + 1;
               k = k + 1;
            }
         }
         if (at >= ak) {
            break;
         }
         a->xZ[p1] = at;
         a->xZ[p0] = ak;
         b->xR[p0] = b->xR[p1];
         b->xR[p1] = bt;
         t = k;
      }
   }
}

// Sorting function optimized for integer keys and real labels, can be used
// to sort middle of the array
//
// A is sorted, and same permutations are applied to B.
//
// NOTES:
//     this function assumes that A[] is finite; it doesn't checks that
//     condition. All other conditions (size of input arrays, etc.) are not
//     checked too.
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void tagsortmiddlei(ZVector *a, ae_int_t offset, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t t;
   ae_int_t tmp;
   ae_int_t p0;
   ae_int_t p1;
   ae_int_t at;
   ae_int_t ak;
   ae_int_t ak1;
// Special cases
   if (n <= 1) {
      return;
   }
// General case, N>1: sort, update B
   for (i = 2; i <= n; i++) {
      t = i;
      while (t != 1) {
         k = t / 2;
         p0 = offset + k - 1;
         p1 = offset + t - 1;
         ak = a->xZ[p0];
         at = a->xZ[p1];
         if (ak >= at) {
            break;
         }
         a->xZ[p0] = at;
         a->xZ[p1] = ak;
         t = k;
      }
   }
   for (i = n - 1; i >= 1; i--) {
      p0 = offset + 0;
      p1 = offset + i;
      tmp = a->xZ[p1];
      a->xZ[p1] = a->xZ[p0];
      a->xZ[p0] = tmp;
      at = tmp;
      t = 0;
      while (true) {
         k = 2 * t + 1;
         if (k + 1 > i) {
            break;
         }
         p0 = offset + t;
         p1 = offset + k;
         ak = a->xZ[p1];
         if (k + 1 < i) {
            ak1 = a->xZ[p1 + 1];
            if (ak1 > ak) {
               ak = ak1;
               p1 = p1 + 1;
               k = k + 1;
            }
         }
         if (at >= ak) {
            break;
         }
         a->xZ[p1] = at;
         a->xZ[p0] = ak;
         t = k;
      }
   }
}

// Sorting function optimized for integer values (only keys, no labels),  can
// be used to sort middle of the array
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void sortmiddlei(ZVector *a, ae_int_t offset, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t t;
   ae_int_t tmp;
   ae_int_t p0;
   ae_int_t p1;
   ae_int_t at;
   ae_int_t ak;
   ae_int_t ak1;
// Special cases
   if (n <= 1) {
      return;
   }
// General case, N>1: sort, update B
   for (i = 2; i <= n; i++) {
      t = i;
      while (t != 1) {
         k = t / 2;
         p0 = offset + k - 1;
         p1 = offset + t - 1;
         ak = a->xZ[p0];
         at = a->xZ[p1];
         if (ak >= at) {
            break;
         }
         a->xZ[p0] = at;
         a->xZ[p1] = ak;
         t = k;
      }
   }
   for (i = n - 1; i >= 1; i--) {
      p0 = offset + 0;
      p1 = offset + i;
      tmp = a->xZ[p1];
      a->xZ[p1] = a->xZ[p0];
      a->xZ[p0] = tmp;
      at = tmp;
      t = 0;
      while (true) {
         k = 2 * t + 1;
         if (k + 1 > i) {
            break;
         }
         p0 = offset + t;
         p1 = offset + k;
         ak = a->xZ[p1];
         if (k + 1 < i) {
            ak1 = a->xZ[p1 + 1];
            if (ak1 > ak) {
               ak = ak1;
               p1 = p1 + 1;
               k = k + 1;
            }
         }
         if (at >= ak) {
            break;
         }
         a->xZ[p1] = at;
         a->xZ[p0] = ak;
         t = k;
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
void tagsortbuf(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2, apbuffers *buf, ae_state *_state) {
   ae_int_t i;
   ae_int_t lv;
   ae_int_t lp;
   ae_int_t rv;
   ae_int_t rp;
// Special cases
   if (n <= 0) {
      return;
   }
   if (n == 1) {
      ivectorsetlengthatleast(p1, 1, _state);
      ivectorsetlengthatleast(p2, 1, _state);
      p1->xZ[0] = 0;
      p2->xZ[0] = 0;
      return;
   }
// General case, N>1: prepare permutations table P1
   ivectorsetlengthatleast(p1, n, _state);
   for (i = 0; i < n; i++) {
      p1->xZ[i] = i;
   }
// General case, N>1: sort, update P1
   rvectorsetlengthatleast(&buf->ra0, n, _state);
   ivectorsetlengthatleast(&buf->ia0, n, _state);
   tagsortfasti(a, p1, &buf->ra0, &buf->ia0, n, _state);
// General case, N>1: fill permutations table P2
//
// To fill P2 we maintain two arrays:
// * PV (Buf.IA0), Position(Value). PV[i] contains position of I-th key at the moment
// * VP (Buf.IA1), Value(Position). VP[i] contains key which has position I at the moment
//
// At each step we making permutation of two items:
//   Left, which is given by position/value pair LP/LV
//   and Right, which is given by RP/RV
// and updating PV[] and VP[] correspondingly.
   ivectorsetlengthatleast(&buf->ia0, n, _state);
   ivectorsetlengthatleast(&buf->ia1, n, _state);
   ivectorsetlengthatleast(p2, n, _state);
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
// * as usual permutation of [0..N-1]. If P1[i]=j, then sorted A[i]  contains
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
void tagsort(RVector *a, ae_int_t n, ZVector *p1, ZVector *p2, ae_state *_state) {
   ae_frame _frame_block;
   ae_frame_make(_state, &_frame_block);
   SetVector(p1);
   SetVector(p2);
   NewObj(apbuffers, buf, _state);
   tagsortbuf(a, n, p1, p2, &buf, _state);
   ae_frame_leave(_state);
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
void tagheappushi(RVector *a, ZVector *b, ae_int_t *n, double va, ae_int_t vb, ae_state *_state) {
   ae_int_t j;
   ae_int_t k;
   double v;
   if (*n < 0) {
      return;
   }
// N=0 is a special case
   if (*n == 0) {
      a->xR[0] = va;
      b->xZ[0] = vb;
      *n = *n + 1;
      return;
   }
// add current point to the heap
// (add to the bottom, then move up)
//
// we don't write point to the heap
// until its final position is determined
// (it allow us to reduce number of array access operations)
   j = *n;
   *n = *n + 1;
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
void tagheapreplacetopi(RVector *a, ZVector *b, ae_int_t n, double va, ae_int_t vb, ae_state *_state) {
   ae_int_t j;
   ae_int_t k1;
   ae_int_t k2;
   double v;
   double v1;
   double v2;
   if (n < 1) {
      return;
   }
// N=1 is a special case
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
void tagheappopi(RVector *a, ZVector *b, ae_int_t *n, ae_state *_state) {
   double va;
   ae_int_t vb;
   if (*n < 1) {
      return;
   }
// N=1 is a special case
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
   *n = *n - 1;
   tagheapreplacetopi(a, b, *n, va, vb, _state);
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
ae_int_t lowerbound(RVector *a, ae_int_t n, double t, ae_state *_state) {
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
         l = l - half - 1;
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
ae_int_t upperbound(RVector *a, ae_int_t n, double t, ae_state *_state) {
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
         l = l - half - 1;
      }
   }
   result = first;
   return result;
}
} // end of namespace alglib_impl

// === BLAS Package ===
namespace alglib_impl {
double vectornorm2(RVector *x, ae_int_t i1, ae_int_t i2, ae_state *_state) {
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
      result = ae_fabs(x->xR[i1], _state);
      return result;
   }
   scl = 0.0;
   ssq = 1.0;
   for (ix = i1; ix <= i2; ix++) {
      if (x->xR[ix] != 0.0) {
         absxi = ae_fabs(x->xR[ix], _state);
         if (scl < absxi) {
            ssq = 1 + ssq * sqr(scl / absxi, _state);
            scl = absxi;
         } else {
            ssq = ssq + sqr(absxi / scl, _state);
         }
      }
   }
   result = scl * ae_sqrt(ssq, _state);
   return result;
}

ae_int_t vectoridxabsmax(RVector *x, ae_int_t i1, ae_int_t i2, ae_state *_state) {
   ae_int_t i;
   ae_int_t result;
   result = i1;
   for (i = i1 + 1; i <= i2; i++) {
      if (ae_fabs(x->xR[i], _state) > ae_fabs(x->xR[result], _state)) {
         result = i;
      }
   }
   return result;
}

ae_int_t columnidxabsmax(RMatrix *x, ae_int_t i1, ae_int_t i2, ae_int_t j, ae_state *_state) {
   ae_int_t i;
   ae_int_t result;
   result = i1;
   for (i = i1 + 1; i <= i2; i++) {
      if (ae_fabs(x->xyR[i][j], _state) > ae_fabs(x->xyR[result][j], _state)) {
         result = i;
      }
   }
   return result;
}

ae_int_t rowidxabsmax(RMatrix *x, ae_int_t j1, ae_int_t j2, ae_int_t i, ae_state *_state) {
   ae_int_t j;
   ae_int_t result;
   result = j1;
   for (j = j1 + 1; j <= j2; j++) {
      if (ae_fabs(x->xyR[i][j], _state) > ae_fabs(x->xyR[i][result], _state)) {
         result = j;
      }
   }
   return result;
}

double upperhessenberg1norm(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   double result;
   ae_assert(i2 - i1 == j2 - j1, "UpperHessenberg1Norm: I2-I1 != J2-J1!", _state);
   for (j = j1; j <= j2; j++) {
      work->xR[j] = 0.0;
   }
   for (i = i1; i <= i2; i++) {
      for (j = maxint(j1, j1 + i - i1 - 1, _state); j <= j2; j++) {
         work->xR[j] = work->xR[j] + ae_fabs(a->xyR[i][j], _state);
      }
   }
   result = 0.0;
   for (j = j1; j <= j2; j++) {
      result = maxreal(result, work->xR[j], _state);
   }
   return result;
}

void copymatrix(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2, ae_state *_state) {
   ae_int_t isrc;
   ae_int_t idst;
   if (is1 > is2 || js1 > js2) {
      return;
   }
   ae_assert(is2 - is1 == id2 - id1, "CopyMatrix: different sizes!", _state);
   ae_assert(js2 - js1 == jd2 - jd1, "CopyMatrix: different sizes!", _state);
   for (isrc = is1; isrc <= is2; isrc++) {
      idst = isrc - is1 + id1;
      ae_v_move(&b->xyR[idst][jd1], 1, &a->xyR[isrc][js1], 1, jd2 - jd1 + 1);
   }
}

void inplacetranspose(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, RVector *work, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t ips;
   ae_int_t jps;
   ae_int_t l;
   if (i1 > i2 || j1 > j2) {
      return;
   }
   ae_assert(i1 - i2 == j1 - j2, "InplaceTranspose error: incorrect array size!", _state);
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

void copyandtranspose(RMatrix *a, ae_int_t is1, ae_int_t is2, ae_int_t js1, ae_int_t js2, RMatrix *b, ae_int_t id1, ae_int_t id2, ae_int_t jd1, ae_int_t jd2, ae_state *_state) {
   ae_int_t isrc;
   ae_int_t jdst;
   if (is1 > is2 || js1 > js2) {
      return;
   }
   ae_assert(is2 - is1 == jd2 - jd1, "CopyAndTranspose: different sizes!", _state);
   ae_assert(js2 - js1 == id2 - id1, "CopyAndTranspose: different sizes!", _state);
   for (isrc = is1; isrc <= is2; isrc++) {
      jdst = isrc - is1 + jd1;
      ae_v_move(&b->xyR[id1][jdst], b->stride, &a->xyR[isrc][js1], 1, id2 - id1 + 1);
   }
}

void matrixvectormultiply(RMatrix *a, ae_int_t i1, ae_int_t i2, ae_int_t j1, ae_int_t j2, bool trans, RVector *x, ae_int_t ix1, ae_int_t ix2, double alpha, RVector *y, ae_int_t iy1, ae_int_t iy2, double beta, ae_state *_state) {
   ae_int_t i;
   double v;
   if (!trans) {
   // y := alpha*A*x + beta*y;
      if (i1 > i2 || j1 > j2) {
         return;
      }
      ae_assert(j2 - j1 == ix2 - ix1, "MatrixVectorMultiply: A and X dont match!", _state);
      ae_assert(i2 - i1 == iy2 - iy1, "MatrixVectorMultiply: A and Y dont match!", _state);
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
         y->xR[iy1 + i - i1] = y->xR[iy1 + i - i1] + alpha * v;
      }
   } else {
   // y := alpha*A'*x + beta*y;
      if (i1 > i2 || j1 > j2) {
         return;
      }
      ae_assert(i2 - i1 == ix2 - ix1, "MatrixVectorMultiply: A and X dont match!", _state);
      ae_assert(j2 - j1 == iy2 - iy1, "MatrixVectorMultiply: A and Y dont match!", _state);
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

double pythag2(double x, double y, ae_state *_state) {
   double w;
   double xabs;
   double yabs;
   double z;
   double result;
   xabs = ae_fabs(x, _state);
   yabs = ae_fabs(y, _state);
   w = maxreal(xabs, yabs, _state);
   z = minreal(xabs, yabs, _state);
   if (z == 0.0) {
      result = w;
   } else {
      result = w * ae_sqrt(1 + sqr(z / w, _state), _state);
   }
   return result;
}

void matrixmatrixmultiply(RMatrix *a, ae_int_t ai1, ae_int_t ai2, ae_int_t aj1, ae_int_t aj2, bool transa, RMatrix *b, ae_int_t bi1, ae_int_t bi2, ae_int_t bj1, ae_int_t bj2, bool transb, double alpha, RMatrix *c, ae_int_t ci1, ae_int_t ci2, ae_int_t cj1, ae_int_t cj2, double beta, RVector *work, ae_state *_state) {
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
   ae_assert(acols == brows, "MatrixMatrixMultiply: incorrect matrix sizes!", _state);
   if (((arows <= 0 || acols <= 0) || brows <= 0) || bcols <= 0) {
      return;
   }
   crows = arows;
// Test WORK
   i = maxint(arows, acols, _state);
   i = maxint(brows, i, _state);
   i = maxint(i, bcols, _state);
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
               c->xyR[ci1 + l - ai1][cj1 + r - bi1] = c->xyR[ci1 + l - ai1][cj1 + r - bi1] + alpha * v;
            }
         }
         return;
      } else {
         for (l = ai1; l <= ai2; l++) {
            for (r = bi1; r <= bi2; r++) {
               v = ae_v_dotproduct(&a->xyR[l][aj1], 1, &b->xyR[r][bj1], 1, aj2 - aj1 + 1);
               c->xyR[ci1 + l - ai1][cj1 + r - bi1] = c->xyR[ci1 + l - ai1][cj1 + r - bi1] + alpha * v;
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
               c->xyR[ci1 + l - aj1][cj1 + r - bi1] = c->xyR[ci1 + l - aj1][cj1 + r - bi1] + alpha * v;
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
// of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
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
void applyrotationsfromtheleft(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work, ae_state *_state) {
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
      // Special case: N1=N2
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
      // Special case: N1=N2
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
// of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
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
void applyrotationsfromtheright(bool isforward, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *c, RVector *s, RMatrix *a, RVector *work, ae_state *_state) {
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
      // Special case: M1=M2
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
      // Special case: M1=M2
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
void generaterotation(double f, double g, double *cs, double *sn, double *r, ae_state *_state) {
   double f1;
   double g1;
   *cs = 0;
   *sn = 0;
   *r = 0;
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
         if (ae_fabs(f1, _state) > ae_fabs(g1, _state)) {
            *r = ae_fabs(f1, _state) * ae_sqrt(1 + sqr(g1 / f1, _state), _state);
         } else {
            *r = ae_fabs(g1, _state) * ae_sqrt(1 + sqr(f1 / g1, _state), _state);
         }
         *cs = f1 / (*r);
         *sn = g1 / (*r);
         if (ae_fabs(f, _state) > ae_fabs(g, _state) && *cs < 0.0) {
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
void rankx(RVector *x, ae_int_t n, bool iscentered, apbuffers *buf, ae_state *_state) {
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
      ae_vector_set_length(&buf->ra1, n, _state);
   }
   if (buf->ia1.cnt < n) {
      ae_vector_set_length(&buf->ia1, n, _state);
   }
   for (i = 0; i < n; i++) {
      buf->ra1.xR[i] = x->xR[i];
      buf->ia1.xZ[i] = i;
   }
   tagsortfasti(&buf->ra1, &buf->ia1, &buf->ra2, &buf->ia2, n, _state);
// Special test for all values being equal
   if (buf->ra1.xR[0] == buf->ra1.xR[n - 1]) {
      if (iscentered) {
         tmp = 0.0;
      } else {
         tmp = (double)(n - 1) / 2.0;
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
         j = j + 1;
      }
      for (k = i; k < j; k++) {
         buf->ra1.xR[k] = (double)(i + j - 1) / 2.0;
      }
      i = j;
   }
// back to x
   if (iscentered) {
      voffs = (double)(n - 1) / 2.0;
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
void rankxuntied(RVector *x, ae_int_t n, apbuffers *buf, ae_state *_state) {
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
      ae_vector_set_length(&buf->ra1, n, _state);
   }
   if (buf->ia1.cnt < n) {
      ae_vector_set_length(&buf->ia1, n, _state);
   }
   for (i = 0; i < n; i++) {
      buf->ra1.xR[i] = x->xR[i];
      buf->ia1.xZ[i] = i;
   }
   tagsortfasti(&buf->ra1, &buf->ia1, &buf->ra2, &buf->ia2, n, _state);
   for (i = 0; i < n; i++) {
      x->xR[buf->ia1.xZ[i]] = (double)i;
   }
}
} // end of namespace alglib_impl

// === TRLINSOLVE Package ===
namespace alglib_impl {
// Obsolete 1-based subroutine.
// See RMatrixTRSafeSolve for 0-based replacement.
void safesolvetriangular(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit, bool normin, RVector *cnorm, ae_state *_state) {
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
   *s = 0;
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
   smlnum = minrealnumber / (machineepsilon * 2);
   bignum = 1 / smlnum;
   *s = 1.0;
   if (!normin) {
      ae_vector_set_length(cnorm, n + 1, _state);
   // Compute the 1-norm of each column, not including the diagonal.
      if (upper) {
      // A is upper triangular.
         for (j = 1; j <= n; j++) {
            v = 0.0;
            for (k = 1; k < j; k++) {
               v = v + ae_fabs(a->xyR[k][j], _state);
            }
            cnorm->xR[j] = v;
         }
      } else {
      // A is lower triangular.
         for (j = 1; j < n; j++) {
            v = 0.0;
            for (k = j + 1; k <= n; k++) {
               v = v + ae_fabs(a->xyR[k][j], _state);
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
      tscal = 1 / (smlnum * tmax);
      ae_v_muld(&cnorm->xR[1], 1, n, tscal);
   }
// Compute a bound on the computed solution vector to see if the
// Level 2 BLAS routine DTRSV can be used.
   j = 1;
   for (k = 2; k <= n; k++) {
      if (ae_fabs(x->xR[k], _state) > ae_fabs(x->xR[j], _state)) {
         j = k;
      }
   }
   xmax = ae_fabs(x->xR[j], _state);
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
         // Initially, G(0) = max{x(i), i=1,...,n}.
            grow = 1 / maxreal(xbnd, smlnum, _state);
            xbnd = grow;
            j = jfirst;
            while ((jinc > 0 && j <= jlast) || (jinc < 0 && j >= jlast)) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // M(j) = G(j-1) / abs(A(j,j))
               tjj = ae_fabs(a->xyR[j][j], _state);
               xbnd = minreal(xbnd, minreal(1.0, tjj, _state) * grow, _state);
               if (tjj + cnorm->xR[j] >= smlnum) {
               // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                  grow = grow * (tjj / (tjj + cnorm->xR[j]));
               } else {
               // G(j) could overflow, set GROW to 0.
                  grow = 0.0;
               }
               if (j == jlast) {
                  grow = xbnd;
               }
               j = j + jinc;
            }
         } else {
         // A is unit triangular.
         //
         // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
            grow = minreal(1.0, 1 / maxreal(xbnd, smlnum, _state), _state);
            j = jfirst;
            while ((jinc > 0 && j <= jlast) || (jinc < 0 && j >= jlast)) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // G(j) = G(j-1)*( 1 + CNORM(j) )
               grow = grow * (1 / (1 + cnorm->xR[j]));
               j = j + jinc;
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
         // Initially, M(0) = max{x(i), i=1,...,n}.
            grow = 1 / maxreal(xbnd, smlnum, _state);
            xbnd = grow;
            j = jfirst;
            while ((jinc > 0 && j <= jlast) || (jinc < 0 && j >= jlast)) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
               xj = 1 + cnorm->xR[j];
               grow = minreal(grow, xbnd / xj, _state);
            // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
               tjj = ae_fabs(a->xyR[j][j], _state);
               if (xj > tjj) {
                  xbnd = xbnd * (tjj / xj);
               }
               if (j == jlast) {
                  grow = minreal(grow, xbnd, _state);
               }
               j = j + jinc;
            }
         } else {
         // A is unit triangular.
         //
         // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
            grow = minreal(1.0, 1 / maxreal(xbnd, smlnum, _state), _state);
            j = jfirst;
            while ((jinc > 0 && j <= jlast) || (jinc < 0 && j >= jlast)) {
            // Exit the loop if the growth factor is too small.
               if (grow <= smlnum) {
                  break;
               }
            // G(j) = ( 1 + CNORM(j) )*G(j-1)
               xj = 1 + cnorm->xR[j];
               grow = grow / xj;
               j = j + jinc;
            }
         }
      }
   }
   if (grow * tscal > smlnum) {
   // Use the Level 2 BLAS solve if the reciprocal of the bound on
   // elements of X is not too small.
      if ((upper && notran) || (!upper && !notran)) {
         if (nounit) {
            vd = a->xyR[n][n];
         } else {
            vd = 1.0;
         }
         x->xR[n] = x->xR[n] / vd;
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
         x->xR[1] = x->xR[1] / vd;
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
         while ((jinc > 0 && j <= jlast) || (jinc < 0 && j >= jlast)) {
         // Compute x(j) = b(j) / A(j,j), scaling x if necessary.
            xj = ae_fabs(x->xR[j], _state);
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
               tjj = ae_fabs(tjjs, _state);
               if (tjj > smlnum) {
               // abs(A(j,j)) > SMLNUM:
                  if (tjj < 1.0) {
                     if (xj > tjj * bignum) {
                     // Scale x by 1/b(j).
                        rec = 1 / xj;
                        ae_v_muld(&x->xR[1], 1, n, rec);
                        *s = *s * rec;
                        xmax = xmax * rec;
                     }
                  }
                  x->xR[j] = x->xR[j] / tjjs;
                  xj = ae_fabs(x->xR[j], _state);
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
                           rec = rec / cnorm->xR[j];
                        }
                        ae_v_muld(&x->xR[1], 1, n, rec);
                        *s = *s * rec;
                        xmax = xmax * rec;
                     }
                     x->xR[j] = x->xR[j] / tjjs;
                     xj = ae_fabs(x->xR[j], _state);
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
               rec = 1 / xj;
               if (cnorm->xR[j] > (bignum - xmax) * rec) {
               // Scale x by 1/(2*abs(x(j))).
                  rec = rec * 0.5;
                  ae_v_muld(&x->xR[1], 1, n, rec);
                  *s = *s * rec;
               }
            } else {
               if (xj * cnorm->xR[j] > bignum - xmax) {
               // Scale x by 1/2.
                  ae_v_muld(&x->xR[1], 1, n, 0.5);
                  *s = *s * 0.5;
               }
            }
            if (upper) {
               if (j > 1) {
               // Compute the update
               // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                  v = x->xR[j] * tscal;
                  jm1 = j - 1;
                  ae_v_subd(&x->xR[1], 1, &a->xyR[1][j], a->stride, jm1, v);
                  i = 1;
                  for (k = 2; k < j; k++) {
                     if (ae_fabs(x->xR[k], _state) > ae_fabs(x->xR[i], _state)) {
                        i = k;
                     }
                  }
                  xmax = ae_fabs(x->xR[i], _state);
               }
            } else {
               if (j < n) {
               // Compute the update
               // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                  jp1 = j + 1;
                  v = x->xR[j] * tscal;
                  ae_v_subd(&x->xR[jp1], 1, &a->xyR[jp1][j], a->stride, n - jp1 + 1, v);
                  i = j + 1;
                  for (k = j + 2; k <= n; k++) {
                     if (ae_fabs(x->xR[k], _state) > ae_fabs(x->xR[i], _state)) {
                        i = k;
                     }
                  }
                  xmax = ae_fabs(x->xR[i], _state);
               }
            }
            j = j + jinc;
         }
      } else {
      // Solve A' * x = b
         j = jfirst;
         while ((jinc > 0 && j <= jlast) || (jinc < 0 && j >= jlast)) {
         // Compute x(j) = b(j) - sum A(k,j)*x(k).
         //   k != j
            xj = ae_fabs(x->xR[j], _state);
            uscal = tscal;
            rec = 1 / maxreal(xmax, 1.0, _state);
            if (cnorm->xR[j] > (bignum - xj) * rec) {
            // If x(j) could overflow, scale x by 1/(2*XMAX).
               rec = rec * 0.5;
               if (nounit) {
                  tjjs = a->xyR[j][j] * tscal;
               } else {
                  tjjs = tscal;
               }
               tjj = ae_fabs(tjjs, _state);
               if (tjj > 1.0) {
               // Divide by A(j,j) when scaling x if A(j,j) > 1.
                  rec = minreal(1.0, rec * tjj, _state);
                  uscal = uscal / tjjs;
               }
               if (rec < 1.0) {
                  ae_v_muld(&x->xR[1], 1, n, rec);
                  *s = *s * rec;
                  xmax = xmax * rec;
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
                     sumj = sumj + v * x->xR[i];
                  }
               } else {
                  if (j < n) {
                     for (i = j + 1; i <= n; i++) {
                        v = a->xyR[i][j] * uscal;
                        sumj = sumj + v * x->xR[i];
                     }
                  }
               }
            }
            if (uscal == tscal) {
            // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
            // was not used to scale the dotproduct.
               x->xR[j] = x->xR[j] - sumj;
               xj = ae_fabs(x->xR[j], _state);
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
                  tjj = ae_fabs(tjjs, _state);
                  if (tjj > smlnum) {
                  // abs(A(j,j)) > SMLNUM:
                     if (tjj < 1.0) {
                        if (xj > tjj * bignum) {
                        // Scale X by 1/abs(x(j)).
                           rec = 1 / xj;
                           ae_v_muld(&x->xR[1], 1, n, rec);
                           *s = *s * rec;
                           xmax = xmax * rec;
                        }
                     }
                     x->xR[j] = x->xR[j] / tjjs;
                  } else {
                     if (tjj > 0.0) {
                     // 0 < abs(A(j,j)) <= SMLNUM:
                        if (xj > tjj * bignum) {
                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                           rec = tjj * bignum / xj;
                           ae_v_muld(&x->xR[1], 1, n, rec);
                           *s = *s * rec;
                           xmax = xmax * rec;
                        }
                        x->xR[j] = x->xR[j] / tjjs;
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
            // Compute x(j) := x(j) / A(j,j)  - sumj if the dot
            // product has already been divided by 1/A(j,j).
               x->xR[j] = x->xR[j] / tjjs - sumj;
            }
            xmax = maxreal(xmax, ae_fabs(x->xR[j], _state), _state);
            j = j + jinc;
         }
      }
      *s = *s / tscal;
   }
// Scale the column norms by 1/TSCAL for return.
   if (tscal != 1.0) {
      v = 1 / tscal;
      ae_v_muld(&cnorm->xR[1], 1, n, v);
   }
}

// Utility subroutine performing the "safe" solution of system of linear
// equations with triangular coefficient matrices.
//
// The subroutine uses scaling and solves the scaled system A*x=s*b (where  s
// is  a  scalar  value)  instead  of  A*x=b,  choosing  s  so  that x can be
// represented by a floating-point number. The closer the system  gets  to  a
// singular, the less s is. If the system is singular, s=0 and x contains the
// non-trivial solution of equation A*x=0.
//
// The feature of an algorithm is that it could not cause an  overflow  or  a
// division by zero regardless of the matrix used as the input.
//
// The algorithm can solve systems of equations with  upper/lower  triangular
// matrices,  with/without unit diagonal, and systems of type A*x=b or A'*x=b
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
//                 A'*x=b, otherwise it is A*x=b.
//     Isunit  -   matrix type. If it is True, the system matrix has  a  unit
//                 diagonal (the elements on the main diagonal are  not  used
//                 in the calculation process), otherwise the matrix is considered
//                 to be a general triangular matrix.
//
// Outputs:
//     X       -   solution. Array whose index ranges within [0..N-1].
//     S       -   scaling factor.
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      June 30, 1992
void rmatrixtrsafesolve(RMatrix *a, ae_int_t n, RVector *x, double *s, bool isupper, bool istrans, bool isunit, ae_state *_state) {
   ae_frame _frame_block;
   bool normin;
   ae_int_t i;
   ae_frame_make(_state, &_frame_block);
   *s = 0;
   NewVector(cnorm, 0, DT_REAL, _state);
   NewMatrix(a1, 0, 0, DT_REAL, _state);
   NewVector(x1, 0, DT_REAL, _state);
// From 0-based to 1-based
   normin = false;
   ae_matrix_set_length(&a1, n + 1, n + 1, _state);
   ae_vector_set_length(&x1, n + 1, _state);
   for (i = 1; i <= n; i++) {
      ae_v_move(&a1.xyR[i][1], 1, a->xyR[i - 1], 1, n);
   }
   ae_v_move(&x1.xR[1], 1, x->xR, 1, n);
// Solve 1-based
   safesolvetriangular(&a1, n, &x1, s, isupper, istrans, isunit, normin, &cnorm, _state);
// From 1-based to 0-based
   ae_v_move(x->xR, 1, &x1.xR[1], 1, n);
   ae_frame_leave(_state);
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
static bool safesolve_cbasicsolveandupdate(complex alpha, complex beta, double lnmax, double bnorm, double maxgrowth, double *xnorm, complex *x, ae_state *_state) {
   double v;
   bool result;
   x->x = 0;
   x->y = 0;
   result = false;
   if (ae_c_eq_d(alpha, 0.0)) {
      return result;
   }
   if (ae_c_neq_d(beta, 0.0)) {
   // alpha*x[i]=beta
      v = ae_log(abscomplex(beta, _state), _state) - ae_log(abscomplex(alpha, _state), _state);
      if (v > lnmax) {
         return result;
      }
      *x = ae_c_div(beta, alpha);
   } else {
   // alpha*x[i]=0
      *x = complex_from_i(0);
   }
// update NrmX, test growth limit
   *xnorm = maxreal(*xnorm, abscomplex(*x, _state), _state);
   if (*xnorm > maxgrowth * bnorm) {
      return result;
   }
   result = true;
   return result;
}

// Real implementation of CMatrixScaledTRSafeSolve
// ALGLIB Routine: Copyright 21.01.2010 by Sergey Bochkanov
bool rmatrixscaledtrsafesolve(RMatrix *a, double sa, ae_int_t n, RVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth, ae_state *_state) {
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
   ae_frame_make(_state, &_frame_block);
   NewVector(tmp, 0, DT_REAL, _state);
   ae_assert(n > 0, "RMatrixTRSafeSolve: incorrect N!", _state);
   ae_assert(trans == 0 || trans == 1, "RMatrixTRSafeSolve: incorrect Trans!", _state);
   result = true;
   lnmax = ae_log(maxrealnumber, _state);
// Quick return if possible
   if (n <= 0) {
      ae_frame_leave(_state);
      return result;
   }
// Load norms: right part and X
   nrmb = 0.0;
   for (i = 0; i < n; i++) {
      nrmb = maxreal(nrmb, ae_fabs(x->xR[i], _state), _state);
   }
   nrmx = 0.0;
// Solve
   ae_vector_set_length(&tmp, n, _state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xR[i] = cx.x;
      }
      ae_frame_leave(_state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xR[i] = cx.x;
      }
      ae_frame_leave(_state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
         if (!result) {
            ae_frame_leave(_state);
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
      ae_frame_leave(_state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
         if (!result) {
            ae_frame_leave(_state);
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
      ae_frame_leave(_state);
      return result;
   }
   result = false;
   ae_frame_leave(_state);
   return result;
}

// Internal subroutine for safe solution of
//
//     SA*op(A)=b
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
bool cmatrixscaledtrsafesolve(CMatrix *a, double sa, ae_int_t n, CVector *x, bool isupper, ae_int_t trans, bool isunit, double maxgrowth, ae_state *_state) {
   ae_frame _frame_block;
   double lnmax;
   double nrmb;
   double nrmx;
   ae_int_t i;
   complex alpha;
   complex beta;
   complex vc;
   bool result;
   ae_frame_make(_state, &_frame_block);
   NewVector(tmp, 0, DT_COMPLEX, _state);
   ae_assert(n > 0, "CMatrixTRSafeSolve: incorrect N!", _state);
   ae_assert((trans == 0 || trans == 1) || trans == 2, "CMatrixTRSafeSolve: incorrect Trans!", _state);
   result = true;
   lnmax = ae_log(maxrealnumber, _state);
// Quick return if possible
   if (n <= 0) {
      ae_frame_leave(_state);
      return result;
   }
// Load norms: right part and X
   nrmb = 0.0;
   for (i = 0; i < n; i++) {
      nrmb = maxreal(nrmb, abscomplex(x->xC[i], _state), _state);
   }
   nrmx = 0.0;
// Solve
   ae_vector_set_length(&tmp, n, _state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xC[i] = vc;
      }
      ae_frame_leave(_state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xC[i] = vc;
      }
      ae_frame_leave(_state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i < n - 1) {
            ae_v_cmoved(&tmp.xC[i + 1], 1, &a->xyC[i][i + 1], 1, "N", n - i - 1, sa);
            ae_v_csubc(&x->xC[i + 1], 1, &tmp.xC[i + 1], 1, "N", n - i - 1, vc);
         }
      }
      ae_frame_leave(_state);
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
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i > 0) {
            ae_v_cmoved(tmp.xC, 1, a->xyC[i], 1, "N", i, sa);
            ae_v_csubc(x->xC, 1, tmp.xC, 1, "N", i, vc);
         }
      }
      ae_frame_leave(_state);
      return result;
   }
   if (isupper && trans == 2) {
   // U^H*x = b
      for (i = 0; i < n; i++) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(conj(a->xyC[i][i], _state), sa);
         }
         beta = x->xC[i];
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i < n - 1) {
            ae_v_cmoved(&tmp.xC[i + 1], 1, &a->xyC[i][i + 1], 1, "Conj", n - i - 1, sa);
            ae_v_csubc(&x->xC[i + 1], 1, &tmp.xC[i + 1], 1, "N", n - i - 1, vc);
         }
      }
      ae_frame_leave(_state);
      return result;
   }
   if (!isupper && trans == 2) {
   // L^T*x = b
      for (i = n - 1; i >= 0; i--) {
      // Task is reduced to alpha*x[i] = beta
         if (isunit) {
            alpha = complex_from_d(sa);
         } else {
            alpha = ae_c_mul_d(conj(a->xyC[i][i], _state), sa);
         }
         beta = x->xC[i];
      // solve alpha*x[i] = beta
         result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
         if (!result) {
            ae_frame_leave(_state);
            return result;
         }
         x->xC[i] = vc;
      // update the rest of right part
         if (i > 0) {
            ae_v_cmoved(tmp.xC, 1, a->xyC[i], 1, "Conj", i, sa);
            ae_v_csubc(x->xC, 1, tmp.xC, 1, "N", i, vc);
         }
      }
      ae_frame_leave(_state);
      return result;
   }
   result = false;
   ae_frame_leave(_state);
   return result;
}
} // end of namespace alglib_impl

// === XBLAS Package ===
namespace alglib_impl {
// Fast Pow
// ALGLIB: Copyright 24.08.2009 by Sergey Bochkanov
static double xblas_xfastpow(double r, ae_int_t n, ae_state *_state) {
   double result;
   result = 0.0;
   if (n > 0) {
      if (n % 2 == 0) {
         result = sqr(xblas_xfastpow(r, n / 2, _state), _state);
      } else {
         result = r * xblas_xfastpow(r, n - 1, _state);
      }
      return result;
   }
   if (n == 0) {
      result = 1.0;
   }
   if (n < 0) {
      result = xblas_xfastpow(1 / r, -n, _state);
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
static void xblas_xsum(RVector *w, double mx, ae_int_t n, double *r, double *rerr, ae_state *_state) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t ks;
   double v;
   double s;
   double ln2;
   double chunk;
   double invchunk;
   bool allzeros;
   *r = 0;
   *rerr = 0;
// special cases:
// * N=0
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
   ae_assert(n < 536870912, "XDot: N is too large!", _state);
// Prepare
   ln2 = ae_log(2.0, _state);
   *rerr = mx * machineepsilon;
// 1. find S such that 0.5 <= S*MX<1
// 2. multiply W by S, so task is normalized in some sense
// 3. S:=1/S so we can obtain original vector multiplying by S
   k = iround(ae_log(mx, _state) / ln2, _state);
   s = xblas_xfastpow(2.0, -k, _state);
   if (!ae_isfinite(s, _state)) {
   // Overflow or underflow during evaluation of S; fallback low-precision code
      *r = 0.0;
      *rerr = mx * machineepsilon;
      for (i = 0; i < n; i++) {
         *r = *r + w->xR[i];
      }
      return;
   }
   while (s * mx >= 1.0) {
      s = 0.5 * s;
   }
   while (s * mx < 0.5) {
      s = 2 * s;
   }
   ae_v_muld(w->xR, 1, n, s);
   s = 1 / s;
// find Chunk=2^M such that N*Chunk<2^29
//
// we have chosen upper limit (2^29) with enough space left
// to tolerate possible problems with rounding and N's close
// to the limit, so we don't want to be very strict here.
   k = itrunc(ae_log(536870912.0 / n, _state) / ln2, _state);
   chunk = xblas_xfastpow(2.0, k, _state);
   if (chunk < 2.0) {
      chunk = 2.0;
   }
   invchunk = 1 / chunk;
// calculate result
   *r = 0.0;
   ae_v_muld(w->xR, 1, n, chunk);
   while (true) {
      s = s * invchunk;
      allzeros = true;
      ks = 0;
      for (i = 0; i < n; i++) {
         v = w->xR[i];
         k = itrunc(v, _state);
         if (v != (double)k) {
            allzeros = false;
         }
         w->xR[i] = chunk * (v - k);
         ks = ks + k;
      }
      *r = *r + s * ks;
      v = ae_fabs(*r, _state);
      if (allzeros || s * n + mx == mx) {
         break;
      }
   }
// correct error
   *rerr = maxreal(*rerr, ae_fabs(*r, _state) * machineepsilon, _state);
}

// More precise dot-product. Absolute error of  subroutine  result  is  about
// 1 ulp of max(MX,V), where:
//     MX = max( |a[i]*b[i]| )
//     V  = |(a,b)|
//
// Inputs:
//     A       -   array[0..N-1], vector 1
//     B       -   array[0..N-1], vector 2
//     N       -   vectors length, N<2^29.
//     Temp    -   array[0..N-1], pre-allocated temporary storage
//
// Outputs:
//     R       -   (A,B)
//     RErr    -   estimate of error. This estimate accounts for both  errors
//                 during  calculation  of  (A,B)  and  errors  introduced by
//                 rounding of A and B to fit in double (about 1 ulp).
// ALGLIB: Copyright 24.08.2009 by Sergey Bochkanov
void xdot(RVector *a, RVector *b, ae_int_t n, RVector *temp, double *r, double *rerr, ae_state *_state) {
   ae_int_t i;
   double mx;
   double v;
   *r = 0;
   *rerr = 0;
// special cases:
// * N=0
   if (n == 0) {
      *r = 0.0;
      *rerr = 0.0;
      return;
   }
   mx = 0.0;
   for (i = 0; i < n; i++) {
      v = a->xR[i] * b->xR[i];
      temp->xR[i] = v;
      mx = maxreal(mx, ae_fabs(v, _state), _state);
   }
   if (mx == 0.0) {
      *r = 0.0;
      *rerr = 0.0;
      return;
   }
   xblas_xsum(temp, mx, n, r, rerr, _state);
}

// More precise complex dot-product. Absolute error of  subroutine  result is
// about 1 ulp of max(MX,V), where:
//     MX = max( |a[i]*b[i]| )
//     V  = |(a,b)|
//
// Inputs:
//     A       -   array[0..N-1], vector 1
//     B       -   array[0..N-1], vector 2
//     N       -   vectors length, N<2^29.
//     Temp    -   array[0..2*N-1], pre-allocated temporary storage
//
// Outputs:
//     R       -   (A,B)
//     RErr    -   estimate of error. This estimate accounts for both  errors
//                 during  calculation  of  (A,B)  and  errors  introduced by
//                 rounding of A and B to fit in double (about 1 ulp).
// ALGLIB: Copyright 27.01.2010 by Sergey Bochkanov
void xcdot(CVector *a, CVector *b, ae_int_t n, RVector *temp, complex *r, double *rerr, ae_state *_state) {
   ae_int_t i;
   double mx;
   double v;
   double rerrx;
   double rerry;
   r->x = 0;
   r->y = 0;
   *rerr = 0;
// special cases:
// * N=0
   if (n == 0) {
      *r = complex_from_i(0);
      *rerr = 0.0;
      return;
   }
// calculate real part
   mx = 0.0;
   for (i = 0; i < n; i++) {
      v = a->xC[i].x * b->xC[i].x;
      temp->xR[2 * i + 0] = v;
      mx = maxreal(mx, ae_fabs(v, _state), _state);
      v = -a->xC[i].y * b->xC[i].y;
      temp->xR[2 * i + 1] = v;
      mx = maxreal(mx, ae_fabs(v, _state), _state);
   }
   if (mx == 0.0) {
      r->x = 0.0;
      rerrx = 0.0;
   } else {
      xblas_xsum(temp, mx, 2 * n, &r->x, &rerrx, _state);
   }
// calculate imaginary part
   mx = 0.0;
   for (i = 0; i < n; i++) {
      v = a->xC[i].x * b->xC[i].y;
      temp->xR[2 * i + 0] = v;
      mx = maxreal(mx, ae_fabs(v, _state), _state);
      v = a->xC[i].y * b->xC[i].x;
      temp->xR[2 * i + 1] = v;
      mx = maxreal(mx, ae_fabs(v, _state), _state);
   }
   if (mx == 0.0) {
      r->y = 0.0;
      rerry = 0.0;
   } else {
      xblas_xsum(temp, mx, 2 * n, &r->y, &rerry, _state);
   }
// total error
   if (rerrx == 0.0 && rerry == 0.0) {
      *rerr = 0.0;
   } else {
      *rerr = maxreal(rerrx, rerry, _state) * ae_sqrt(1 + sqr(minreal(rerrx, rerry, _state) / maxreal(rerrx, rerry, _state), _state), _state);
   }
}
} // end of namespace alglib_impl

// === LINMIN Package ===
namespace alglib_impl {
static const double linmin_ftol = 0.001;
static const double linmin_xtol = 100 * machineepsilon;
static const ae_int_t linmin_maxfev = 20;
static const double linmin_stpmin = 1.0E-50;
static const double linmin_defstpmax = 1.0E+50;
static const double linmin_armijofactor = 1.3;

// Normalizes direction/step pair: makes |D|=1, scales Stp.
// If |D|=0, it returns, leavind D/Stp unchanged.
// ALGLIB: Copyright 01.04.2010 by Sergey Bochkanov
void linminnormalized(RVector *d, double *stp, ae_int_t n, ae_state *_state) {
   double mx;
   double s;
   ae_int_t i;
// first, scale D to avoid underflow/overflow durng squaring
   mx = 0.0;
   for (i = 0; i < n; i++) {
      mx = maxreal(mx, ae_fabs(d->xR[i], _state), _state);
   }
   if (mx == 0.0) {
      return;
   }
   s = 1 / mx;
   ae_v_muld(d->xR, 1, n, s);
   *stp = *stp / s;
// normalize D
   s = ae_v_dotproduct(d->xR, 1, d->xR, 1, n);
   s = 1 / ae_sqrt(s, _state);
   ae_v_muld(d->xR, 1, n, s);
   *stp = *stp / s;
}

static void linmin_mcstep(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, double fp, double dp, bool *brackt, double stmin, double stmax, ae_int_t *info, ae_state *_state) {
   bool bound;
   double gamma;
   double p;
   double q;
   double r;
   double s;
   double sgnd;
   double stpc;
   double stpf;
   double stpq;
   double theta;
   *info = 0;
//     CHECK THE INPUT PARAMETERS FOR ERRORS.
   if (((*brackt && (*stp <= minreal(*stx, *sty, _state) || *stp >= maxreal(*stx, *sty, _state))) || *dx * (*stp - (*stx)) >= 0.0) || stmax < stmin) {
      return;
   }
//     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
   sgnd = dp * (*dx / ae_fabs(*dx, _state));
//     FIRST CASE. A HIGHER FUNCTION VALUE.
//     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
//     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
//     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
   if (fp > *fx) {
      *info = 1;
      bound = true;
      theta = 3 * (*fx - fp) / (*stp - (*stx)) + (*dx) + dp;
      s = maxreal(ae_fabs(theta, _state), maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
      gamma = s * ae_sqrt(sqr(theta / s, _state) - *dx / s * (dp / s), _state);
      if (*stp < *stx) {
         gamma = -gamma;
      }
      p = gamma - (*dx) + theta;
      q = gamma - (*dx) + gamma + dp;
      r = p / q;
      stpc = *stx + r * (*stp - (*stx));
      stpq = *stx + *dx / ((*fx - fp) / (*stp - (*stx)) + (*dx)) / 2 * (*stp - (*stx));
      if (ae_fabs(stpc - (*stx), _state) < ae_fabs(stpq - (*stx), _state)) {
         stpf = stpc;
      } else {
         stpf = stpc + (stpq - stpc) / 2;
      }
      *brackt = true;
   } else {
      if (sgnd < 0.0) {
      //     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
      //     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
      //     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
      //     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
         *info = 2;
         bound = false;
         theta = 3 * (*fx - fp) / (*stp - (*stx)) + (*dx) + dp;
         s = maxreal(ae_fabs(theta, _state), maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
         gamma = s * ae_sqrt(sqr(theta / s, _state) - *dx / s * (dp / s), _state);
         if (*stp > *stx) {
            gamma = -gamma;
         }
         p = gamma - dp + theta;
         q = gamma - dp + gamma + (*dx);
         r = p / q;
         stpc = *stp + r * (*stx - (*stp));
         stpq = *stp + dp / (dp - (*dx)) * (*stx - (*stp));
         if (ae_fabs(stpc - (*stp), _state) > ae_fabs(stpq - (*stp), _state)) {
            stpf = stpc;
         } else {
            stpf = stpq;
         }
         *brackt = true;
      } else {
         if (ae_fabs(dp, _state) < ae_fabs(*dx, _state)) {
         //     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
         //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
         //     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
         //     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
         //     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
         //     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
         //     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
         //     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
            *info = 3;
            bound = true;
            theta = 3 * (*fx - fp) / (*stp - (*stx)) + (*dx) + dp;
            s = maxreal(ae_fabs(theta, _state), maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
         //        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
         //        TO INFINITY IN THE DIRECTION OF THE STEP.
            gamma = s * ae_sqrt(maxreal(0.0, sqr(theta / s, _state) - *dx / s * (dp / s), _state), _state);
            if (*stp > *stx) {
               gamma = -gamma;
            }
            p = gamma - dp + theta;
            q = gamma + (*dx - dp) + gamma;
            r = p / q;
            if (r < 0.0 && gamma != 0.0) {
               stpc = *stp + r * (*stx - (*stp));
            } else {
               if (*stp > *stx) {
                  stpc = stmax;
               } else {
                  stpc = stmin;
               }
            }
            stpq = *stp + dp / (dp - (*dx)) * (*stx - (*stp));
            if (*brackt) {
               if (ae_fabs(*stp - stpc, _state) < ae_fabs(*stp - stpq, _state)) {
                  stpf = stpc;
               } else {
                  stpf = stpq;
               }
            } else {
               if (ae_fabs(*stp - stpc, _state) > ae_fabs(*stp - stpq, _state)) {
                  stpf = stpc;
               } else {
                  stpf = stpq;
               }
            }
         } else {
         //     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
         //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
         //     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
         //     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
            *info = 4;
            bound = false;
            if (*brackt) {
               theta = 3 * (fp - (*fy)) / (*sty - (*stp)) + (*dy) + dp;
               s = maxreal(ae_fabs(theta, _state), maxreal(ae_fabs(*dy, _state), ae_fabs(dp, _state), _state), _state);
               gamma = s * ae_sqrt(sqr(theta / s, _state) - *dy / s * (dp / s), _state);
               if (*stp > *sty) {
                  gamma = -gamma;
               }
               p = gamma - dp + theta;
               q = gamma - dp + gamma + (*dy);
               r = p / q;
               stpc = *stp + r * (*sty - (*stp));
               stpf = stpc;
            } else {
               if (*stp > *stx) {
                  stpf = stmax;
               } else {
                  stpf = stmin;
               }
            }
         }
      }
   }
//     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
//     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
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
//     COMPUTE THE NEW STEP AND SAFEGUARD IT.
   stpf = minreal(stmax, stpf, _state);
   stpf = maxreal(stmin, stpf, _state);
   *stp = stpf;
   if (*brackt && bound) {
      if (*sty > *stx) {
         *stp = minreal(*stx + 0.66 * (*sty - (*stx)), *stp, _state);
      } else {
         *stp = maxreal(*stx + 0.66 * (*sty - (*stx)), *stp, _state);
      }
   }
}

// THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
// DECREASE CONDITION AND A CURVATURE CONDITION.
//
// AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
// ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
// SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION
//
//     F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
//
// IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
// FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
// UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).
//
// THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
// DECREASE CONDITION
//
//     F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
//
// AND THE CURVATURE CONDITION
//
//     ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
//
// IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
// BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
// IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
// ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
// IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.
//
//
// :::::::::::::IMPORTANT NOTES:::::::::::::
//
// NOTE 1:
//
// This routine  guarantees that it will stop at the last point where function
// value was calculated. It won't make several additional function evaluations
// after finding good point. So if you store function evaluations requested by
// this routine, you can be sure that last one is the point where we've stopped.
//
// NOTE 2:
//
// when 0<StpMax<StpMin, algorithm will terminate with INFO=5 and Stp=StpMax
//
// NOTE 3:
//
// this algorithm guarantees that, if MCINFO=1 or MCINFO=5, then:
// * F(final_point)<F(initial_point) - strict inequality
// * final_point != initial_point - after rounding to machine precision
//
// NOTE 4:
//
// when non-descent direction is specified, algorithm stops with MCINFO=0,
// Stp=0 and initial point at X[].
// :::::::::::::::::::::::::::::::::::::::::
//
//
// PARAMETERS DESCRIPRION
//
// STAGE IS ZERO ON FIRST CALL, ZERO ON FINAL EXIT
//
// N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.
//
// X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
// THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.
//
// F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
// IT CONTAINS THE VALUE OF F AT X + STP*S.
//
// G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
// ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.
//
// S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.
//
// STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
// OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.
//
// FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
// SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
// SATISFIED.
//
// XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
// WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.
//
// STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
// UPPER BOUNDS FOR THE STEP.
//
// MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
// NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.
//
// INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
//     INFO = 0  IMPROPER INPUT PARAMETERS.
//
//     INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
//               DIRECTIONAL DERIVATIVE CONDITION HOLD.
//
//     INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
//               IS AT MOST XTOL.
//
//     INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
//
//     INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
//
//     INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
//
//     INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
//               THERE MAY NOT BE A STEP WHICH SATISFIES THE
//               SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
//               TOLERANCES MAY BE TOO SMALL.
//
// NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.
//
// WA IS A WORK ARRAY OF LENGTH N.
//
// ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
// JORGE J. MORE', DAVID J. THUENTE
void mcsrch(ae_int_t n, RVector *x, double *f, RVector *g, RVector *s, double *stp, double stpmax, double gtol, ae_int_t *info, ae_int_t *nfev, RVector *wa, linminstate *state, ae_int_t *stage, ae_state *_state) {
   ae_int_t i;
   double v;
   double p5;
   double p66;
   double zero;
// init
   p5 = 0.5;
   p66 = 0.66;
   state->xtrapf = 4.0;
   zero = 0.0;
   if (stpmax == 0.0) {
      stpmax = linmin_defstpmax;
   }
   if (*stp < linmin_stpmin) {
      *stp = linmin_stpmin;
   }
   if (*stp > stpmax) {
      *stp = stpmax;
   }
// Main cycle
   while (true) {
      if (*stage == 0) {
      // NEXT
         *stage = 2;
         continue;
      }
      if (*stage == 2) {
         state->infoc = 1;
         *info = 0;
      //     CHECK THE INPUT PARAMETERS FOR ERRORS.
         if (stpmax < linmin_stpmin && stpmax > 0.0) {
            *info = 5;
            *stp = stpmax;
            *stage = 0;
            return;
         }
         if (((((((n <= 0 || *stp <= 0.0) || linmin_ftol < 0.0) || gtol < zero) || linmin_xtol < zero) || linmin_stpmin < zero) || stpmax < linmin_stpmin) || linmin_maxfev <= 0) {
            *stage = 0;
            return;
         }
      //     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
      //     AND CHECK THAT S IS A DESCENT DIRECTION.
         v = ae_v_dotproduct(g->xR, 1, s->xR, 1, n);
         state->dginit = v;
         if (state->dginit >= 0.0) {
            *stage = 0;
            *stp = 0.0;
            return;
         }
      //     INITIALIZE LOCAL VARIABLES.
         state->brackt = false;
         state->stage1 = true;
         *nfev = 0;
         state->finit = *f;
         state->dgtest = linmin_ftol * state->dginit;
         state->width = stpmax - linmin_stpmin;
         state->width1 = state->width / p5;
         ae_v_move(wa->xR, 1, x->xR, 1, n);
      //     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
      //     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
      //     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
      //     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
      //     THE INTERVAL OF UNCERTAINTY.
      //     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
      //     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
         state->stx = 0.0;
         state->fx = state->finit;
         state->dgx = state->dginit;
         state->sty = 0.0;
         state->fy = state->finit;
         state->dgy = state->dginit;
      // NEXT
         *stage = 3;
         continue;
      }
      if (*stage == 3) {
      //     START OF ITERATION.
      //
      //     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
      //     TO THE PRESENT INTERVAL OF UNCERTAINTY.
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
      //        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
         if (*stp > stpmax) {
            *stp = stpmax;
         }
         if (*stp < linmin_stpmin) {
            *stp = linmin_stpmin;
         }
      //        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
      //        STP BE THE LOWEST POINT OBTAINED SO FAR.
         if ((((state->brackt && (*stp <= state->stmin || *stp >= state->stmax)) || *nfev >= linmin_maxfev - 1) || state->infoc == 0) || (state->brackt && state->stmax - state->stmin <= linmin_xtol * state->stmax)) {
            *stp = state->stx;
         }
      //        EVALUATE THE FUNCTION AND GRADIENT AT STP
      //        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
         ae_v_move(x->xR, 1, wa->xR, 1, n);
         ae_v_addd(x->xR, 1, s->xR, 1, n, *stp);
      // NEXT
         *stage = 4;
         return;
      }
      if (*stage == 4) {
         *info = 0;
         *nfev = *nfev + 1;
         v = ae_v_dotproduct(g->xR, 1, s->xR, 1, n);
         state->dg = v;
         state->ftest1 = state->finit + *stp * state->dgtest;
      //        TEST FOR CONVERGENCE.
         if ((state->brackt && (*stp <= state->stmin || *stp >= state->stmax)) || state->infoc == 0) {
            *info = 6;
         }
         if (((*stp == stpmax && *f < state->finit) && *f <= state->ftest1) && state->dg <= state->dgtest) {
            *info = 5;
         }
         if (*stp == linmin_stpmin && ((*f >= state->finit || *f > state->ftest1) || state->dg >= state->dgtest)) {
            *info = 4;
         }
         if (*nfev >= linmin_maxfev) {
            *info = 3;
         }
         if (state->brackt && state->stmax - state->stmin <= linmin_xtol * state->stmax) {
            *info = 2;
         }
         if ((*f < state->finit && *f <= state->ftest1) && ae_fabs(state->dg, _state) <= -gtol * state->dginit) {
            *info = 1;
         }
      //        CHECK FOR TERMINATION.
         if (*info != 0) {
         // Check guarantees provided by the function for INFO=1 or INFO=5
            if (*info == 1 || *info == 5) {
               v = 0.0;
               for (i = 0; i < n; i++) {
                  v = v + (wa->xR[i] - x->xR[i]) * (wa->xR[i] - x->xR[i]);
               }
               if (*f >= state->finit || v == 0.0) {
                  *info = 6;
               }
            }
            *stage = 0;
            return;
         }
      //        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
      //        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
         if ((state->stage1 && *f <= state->ftest1) && state->dg >= minreal(linmin_ftol, gtol, _state) * state->dginit) {
            state->stage1 = false;
         }
      //        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
      //        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
      //        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
      //        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
      //        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
         if ((state->stage1 && *f <= state->fx) && *f > state->ftest1) {
         //           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
            state->fm = *f - *stp * state->dgtest;
            state->fxm = state->fx - state->stx * state->dgtest;
            state->fym = state->fy - state->sty * state->dgtest;
            state->dgm = state->dg - state->dgtest;
            state->dgxm = state->dgx - state->dgtest;
            state->dgym = state->dgy - state->dgtest;
         //           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
         //           AND TO COMPUTE THE NEW STEP.
            linmin_mcstep(&state->stx, &state->fxm, &state->dgxm, &state->sty, &state->fym, &state->dgym, stp, state->fm, state->dgm, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
         //           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
            state->fx = state->fxm + state->stx * state->dgtest;
            state->fy = state->fym + state->sty * state->dgtest;
            state->dgx = state->dgxm + state->dgtest;
            state->dgy = state->dgym + state->dgtest;
         } else {
         //           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
         //           AND TO COMPUTE THE NEW STEP.
            linmin_mcstep(&state->stx, &state->fx, &state->dgx, &state->sty, &state->fy, &state->dgy, stp, *f, state->dg, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
         }
      //        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
      //        INTERVAL OF UNCERTAINTY.
         if (state->brackt) {
            if (ae_fabs(state->sty - state->stx, _state) >= p66 * state->width1) {
               *stp = state->stx + p5 * (state->sty - state->stx);
            }
            state->width1 = state->width;
            state->width = ae_fabs(state->sty - state->stx, _state);
         }
      //  NEXT.
         *stage = 3;
         continue;
      }
   }
}

// These functions perform Armijo line search using  at  most  FMAX  function
// evaluations.  It  doesn't  enforce  some  kind  of  " sufficient decrease"
// criterion - it just tries different Armijo steps and returns optimum found
// so far.
//
// Optimization is done using F-rcomm interface:
// * ArmijoCreate initializes State structure
//   (reusing previously allocated buffers)
// * ArmijoIteration is subsequently called
// * ArmijoResults returns results
//
// Inputs:
//     N       -   problem size
//     X       -   array[N], starting point
//     F       -   F(X+S*STP)
//     S       -   step direction, S>0
//     STP     -   step length
//     STPMAX  -   maximum value for STP or zero (if no limit is imposed)
//     FMAX    -   maximum number of function evaluations
//     State   -   optimization state
// ALGLIB: Copyright 05.10.2010 by Sergey Bochkanov
void armijocreate(ae_int_t n, RVector *x, double f, RVector *s, double stp, double stpmax, ae_int_t fmax, armijostate *state, ae_state *_state) {
   if (state->x.cnt < n) {
      ae_vector_set_length(&state->x, n, _state);
   }
   if (state->xbase.cnt < n) {
      ae_vector_set_length(&state->xbase, n, _state);
   }
   if (state->s.cnt < n) {
      ae_vector_set_length(&state->s, n, _state);
   }
   state->stpmax = stpmax;
   state->fmax = fmax;
   state->stplen = stp;
   state->fcur = f;
   state->n = n;
   ae_v_move(state->xbase.xR, 1, x->xR, 1, n);
   ae_v_move(state->s.xR, 1, s->xR, 1, n);
   ae_vector_set_length(&state->rstate.ia, 0 + 1, _state);
   ae_vector_set_length(&state->rstate.ra, 0 + 1, _state);
   state->rstate.stage = -1;
}

// This is rcomm-based search function
// ALGLIB: Copyright 05.10.2010 by Sergey Bochkanov
bool armijoiteration(armijostate *state, ae_state *_state) {
   double v;
   ae_int_t n;
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
      n = state->rstate.ia.xZ[0];
      v = state->rstate.ra.xR[0];
   } else {
      n = 359;
      v = -58;
   }
   if (state->rstate.stage == 0) {
      goto lbl_0;
   }
   if (state->rstate.stage == 1) {
      goto lbl_1;
   }
   if (state->rstate.stage == 2) {
      goto lbl_2;
   }
   if (state->rstate.stage == 3) {
      goto lbl_3;
   }
// Routine body
   if ((state->stplen <= 0.0 || state->stpmax < 0.0) || state->fmax < 2) {
      state->info = 0;
      result = false;
      return result;
   }
   if (state->stplen <= linmin_stpmin) {
      state->info = 4;
      result = false;
      return result;
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
   v = state->stplen * linmin_armijofactor;
   if (v > state->stpmax && state->stpmax != 0.0) {
      v = state->stpmax;
   }
   ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
   ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
   state->rstate.stage = 0;
   goto lbl_rcomm;
lbl_0:
   state->nfev = state->nfev + 1;
   if (state->f >= state->fcur) {
      goto lbl_4;
   }
   state->stplen = v;
   state->fcur = state->f;
lbl_6:
   if (false) {
      goto lbl_7;
   }
// test stopping conditions
   if (state->nfev >= state->fmax) {
      state->info = 3;
      result = false;
      return result;
   }
   if (state->stplen >= state->stpmax) {
      state->info = 5;
      result = false;
      return result;
   }
// evaluate F
   v = state->stplen * linmin_armijofactor;
   if (v > state->stpmax && state->stpmax != 0.0) {
      v = state->stpmax;
   }
   ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
   ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
   state->rstate.stage = 1;
   goto lbl_rcomm;
lbl_1:
   state->nfev = state->nfev + 1;
// make decision
   if (state->f < state->fcur) {
      state->stplen = v;
      state->fcur = state->f;
   } else {
      state->info = 1;
      result = false;
      return result;
   }
   goto lbl_6;
lbl_7:
lbl_4:
// Decrease length
   v = state->stplen / linmin_armijofactor;
   ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
   ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
   state->rstate.stage = 2;
   goto lbl_rcomm;
lbl_2:
   state->nfev = state->nfev + 1;
   if (state->f >= state->fcur) {
      goto lbl_8;
   }
   state->stplen = state->stplen / linmin_armijofactor;
   state->fcur = state->f;
lbl_10:
   if (false) {
      goto lbl_11;
   }
// test stopping conditions
   if (state->nfev >= state->fmax) {
      state->info = 3;
      result = false;
      return result;
   }
   if (state->stplen <= linmin_stpmin) {
      state->info = 4;
      result = false;
      return result;
   }
// evaluate F
   v = state->stplen / linmin_armijofactor;
   ae_v_move(state->x.xR, 1, state->xbase.xR, 1, n);
   ae_v_addd(state->x.xR, 1, state->s.xR, 1, n, v);
   state->rstate.stage = 3;
   goto lbl_rcomm;
lbl_3:
   state->nfev = state->nfev + 1;
// make decision
   if (state->f < state->fcur) {
      state->stplen = state->stplen / linmin_armijofactor;
      state->fcur = state->f;
   } else {
      state->info = 1;
      result = false;
      return result;
   }
   goto lbl_10;
lbl_11:
lbl_8:
// Nothing to be done
   state->info = 1;
   result = false;
   return result;
// Saving state
lbl_rcomm:
   result = true;
   state->rstate.ia.xZ[0] = n;
   state->rstate.ra.xR[0] = v;
   return result;
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
void armijoresults(armijostate *state, ae_int_t *info, double *stp, double *f, ae_state *_state) {
   *info = state->info;
   *stp = state->stplen;
   *f = state->fcur;
}

void linminstate_init(void *_p, ae_state *_state, bool make_automatic) {
   linminstate *p = (linminstate *)_p;
   ae_touch_ptr((void *)p);
}

void linminstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   linminstate *dst = (linminstate *)_dst;
   linminstate *src = (linminstate *)_src;
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
   linminstate *p = (linminstate *)_p;
   ae_touch_ptr((void *)p);
}

void armijostate_init(void *_p, ae_state *_state, bool make_automatic) {
   armijostate *p = (armijostate *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->xbase, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->s, 0, DT_REAL, _state, make_automatic);
   rcommstate_init(&p->rstate, _state, make_automatic);
}

void armijostate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   armijostate *dst = (armijostate *)_dst;
   armijostate *src = (armijostate *)_src;
   dst->needf = src->needf;
   ae_vector_copy(&dst->x, &src->x, _state, make_automatic);
   dst->f = src->f;
   dst->n = src->n;
   ae_vector_copy(&dst->xbase, &src->xbase, _state, make_automatic);
   ae_vector_copy(&dst->s, &src->s, _state, make_automatic);
   dst->stplen = src->stplen;
   dst->fcur = src->fcur;
   dst->stpmax = src->stpmax;
   dst->fmax = src->fmax;
   dst->nfev = src->nfev;
   dst->info = src->info;
   rcommstate_copy(&dst->rstate, &src->rstate, _state, make_automatic);
}

void armijostate_free(void *_p, bool make_automatic) {
   armijostate *p = (armijostate *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->xbase, make_automatic);
   ae_vector_free(&p->s, make_automatic);
   rcommstate_free(&p->rstate, make_automatic);
}
} // end of namespace alglib_impl

// === NEARUNITYUNIT Package ===
namespace alglib_impl {
double nulog1p(double x, ae_state *_state) {
   double z;
   double lp;
   double lq;
   double result;
   z = 1.0 + x;
   if (z < 0.70710678118654752440 || z > 1.41421356237309504880) {
      result = ae_log(z, _state);
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

double nuexpm1(double x, ae_state *_state) {
   double r;
   double xx;
   double ep;
   double eq;
   double result;
   if (x < -0.5 || x > 0.5) {
      result = ae_exp(x, _state) - 1.0;
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
   r = r / (eq - r);
   result = r + r;
   return result;
}

double nucosm1(double x, ae_state *_state) {
   double xx;
   double c;
   double result;
   if (x < -0.25 * pi || x > 0.25 * pi) {
      result = ae_cos(x, _state) - 1;
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
static bool ntheory_isprime(ae_int_t n, ae_state *_state) {
   ae_int_t p;
   bool result;
   result = false;
   p = 2;
   while (p * p <= n) {
      if (n % p == 0) {
         return result;
      }
      p = p + 1;
   }
   result = true;
   return result;
}

static ae_int_t ntheory_modmul(ae_int_t a, ae_int_t b, ae_int_t n, ae_state *_state) {
   ae_int_t t;
   double ra;
   double rb;
   ae_int_t result;
   ae_assert(a >= 0 && a < n, "ModMul: A<0 or A >= N", _state);
   ae_assert(b >= 0 && b < n, "ModMul: B<0 or B >= N", _state);
// Base cases
   ra = (double)a;
   rb = (double)b;
   if (b == 0 || a == 0) {
      result = 0;
      return result;
   }
   if (b == 1 || a == 1) {
      result = a * b;
      return result;
   }
   if (ra * rb == (double)(a * b)) {
      result = a * b % n;
      return result;
   }
// Non-base cases
   if (b % 2 == 0) {
   // A*B = (A*(B/2)) * 2
   //
   // Product T=A*(B/2) is calculated recursively, product T*2 is
   // calculated as follows:
   // * result:=T-N
   // * result:=result+T
   // * if result<0 then result:=result+N
   //
   // In case integer result overflows, we generate exception
      t = ntheory_modmul(a, b / 2, n, _state);
      result = t - n;
      result = result + t;
      if (result < 0) {
         result = result + n;
      }
   } else {
   // A*B = (A*(B div 2)) * 2 + A
   //
   // Product T=A*(B/2) is calculated recursively, product T*2 is
   // calculated as follows:
   // * result:=T-N
   // * result:=result+T
   // * if result<0 then result:=result+N
   //
   // In case integer result overflows, we generate exception
      t = ntheory_modmul(a, b / 2, n, _state);
      result = t - n;
      result = result + t;
      if (result < 0) {
         result = result + n;
      }
      result = result - n;
      result = result + a;
      if (result < 0) {
         result = result + n;
      }
   }
   return result;
}

static ae_int_t ntheory_modexp(ae_int_t a, ae_int_t b, ae_int_t n, ae_state *_state) {
   ae_int_t t;
   ae_int_t result;
   ae_assert(a >= 0 && a < n, "ModExp: A<0 or A >= N", _state);
   ae_assert(b >= 0, "ModExp: B<0", _state);
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
      t = ntheory_modmul(a, a, n, _state);
      result = ntheory_modexp(t, b / 2, n, _state);
   } else {
      t = ntheory_modmul(a, a, n, _state);
      result = ntheory_modexp(t, b / 2, n, _state);
      result = ntheory_modmul(result, a, n, _state);
   }
   return result;
}

void findprimitiverootandinverse(ae_int_t n, ae_int_t *proot, ae_int_t *invproot, ae_state *_state) {
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
   ae_assert(n >= 3, "FindPrimitiveRootAndInverse: N<3", _state);
   *proot = 0;
   *invproot = 0;
// check that N is prime
   ae_assert(ntheory_isprime(n, _state), "FindPrimitiveRoot: N is not prime", _state);
// Because N is prime, Euler totient function is equal to N-1
   phin = n - 1;
// Test different values of PRoot - from 2 to N-1.
// One of these values MUST be primitive root.
//
// For testing we use algorithm from Wiki (Primitive root modulo n):
// * compute phi(N)
// * determine the different prime factors of phi(N), say p1, ..., pk
// * for every element m of Zn*, compute m^(phi(N)/pi) mod N for i=1..k
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
            t = ntheory_modexp(candroot, phin / f, n, _state);
            if (t == 1) {
               allnonone = false;
               break;
            }
            while (q % f == 0) {
               q = q / f;
            }
         }
         f = f + 1;
      }
      if (allnonone) {
         *proot = candroot;
         break;
      }
   }
   ae_assert(*proot >= 2, "FindPrimitiveRoot: internal error (root not found)", _state);
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
      lastx = lastx + n;
   }
   *invproot = lastx;
// Check that it is safe to perform multiplication modulo N.
// Check results for consistency.
   n2 = (n - 1) * (n - 1);
   ae_assert(n2 / (n - 1) == n - 1, "FindPrimitiveRoot: internal error", _state);
   ae_assert(*proot * (*invproot) / (*proot) == (*invproot), "FindPrimitiveRoot: internal error", _state);
   ae_assert(*proot * (*invproot) / (*invproot) == (*proot), "FindPrimitiveRoot: internal error", _state);
   ae_assert(*proot * (*invproot) % n == 1, "FindPrimitiveRoot: internal error", _state);
}
} // end of namespace alglib_impl

// === FTBASE Package ===
// Depends on: APSERV, NTHEORY
namespace alglib_impl {
static const ae_int_t ftbase_coltype = 0;
static const ae_int_t ftbase_coloperandscnt = 1;
static const ae_int_t ftbase_coloperandsize = 2;
static const ae_int_t ftbase_colmicrovectorsize = 3;
static const ae_int_t ftbase_colparam0 = 4;
static const ae_int_t ftbase_colparam1 = 5;
static const ae_int_t ftbase_colparam2 = 6;
static const ae_int_t ftbase_colparam3 = 7;
static const ae_int_t ftbase_colscnt = 8;
static const ae_int_t ftbase_opend = 0;
static const ae_int_t ftbase_opcomplexreffft = 1;
static const ae_int_t ftbase_opbluesteinsfft = 2;
static const ae_int_t ftbase_opcomplexcodeletfft = 3;
static const ae_int_t ftbase_opcomplexcodelettwfft = 4;
static const ae_int_t ftbase_opradersfft = 5;
static const ae_int_t ftbase_opcomplextranspose = -1;
static const ae_int_t ftbase_opcomplexfftfactors = -2;
static const ae_int_t ftbase_opstart = -3;
static const ae_int_t ftbase_opjmp = -4;
static const ae_int_t ftbase_opparallelcall = -5;
static const ae_int_t ftbase_maxradix = 6;
static const ae_int_t ftbase_updatetw = 16;
static const ae_int_t ftbase_recursivethreshold = 1024;
static const ae_int_t ftbase_raderthreshold = 19;
static const ae_int_t ftbase_ftbasecodeletrecommended = 5;
static const double ftbase_ftbaseinefficiencyfactor = 1.3;
static const ae_int_t ftbase_ftbasemaxsmoothfactor = 5;

// Returns optimistic estimate of the FFT cost, in UNITs (1 UNIT = 100 KFLOPs)
//
// Inputs:
//     N       -   task size, N>0
//
// Result:
//     cost in UNITs, rounded down to nearest integer
//
// NOTE: If FFT cost is less than 1 UNIT, it will return 0 as result.
// ALGLIB: Copyright 08.04.2013 by Sergey Bochkanov
static ae_int_t ftbase_ftoptimisticestimate(ae_int_t n, ae_state *_state) {
   ae_int_t result;
   ae_assert(n > 0, "FTOptimisticEstimate: N <= 0", _state);
   result = ifloor(1.0E-5 * 5 * n * ae_log((double)n, _state) / ae_log(2.0, _state), _state);
   return result;
}

// This subroutine applies complex reference FFT to input/output array A.
//
// VERY SLOW OPERATION, do not use it in real life plans :)
//
// Inputs:
//     A           -   array, must be large enough for plan to work
//     Offs        -   offset of the subarray to process
//     OperandsCnt -   operands count (see description of FastTransformPlan)
//     OperandSize -   operand size (see description of FastTransformPlan)
//     MicrovectorSize-microvector size (see description of FastTransformPlan)
//     Buf         -   temporary array, must be at least OperandsCnt*OperandSize*MicrovectorSize
//
// Outputs:
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplycomplexreffft(RVector *a, ae_int_t offs, ae_int_t operandscnt, ae_int_t operandsize, ae_int_t microvectorsize, RVector *buf, ae_state *_state) {
   ae_int_t opidx;
   ae_int_t i;
   ae_int_t k;
   double hre;
   double him;
   double c;
   double s;
   double re;
   double im;
   ae_int_t n;
   ae_assert(operandscnt >= 1, "FTApplyComplexRefFFT: OperandsCnt<1", _state);
   ae_assert(operandsize >= 1, "FTApplyComplexRefFFT: OperandSize<1", _state);
   ae_assert(microvectorsize == 2, "FTApplyComplexRefFFT: MicrovectorSize != 2", _state);
   n = operandsize;
   for (opidx = 0; opidx < operandscnt; opidx++) {
      for (i = 0; i < n; i++) {
         hre = 0.0;
         him = 0.0;
         for (k = 0; k < n; k++) {
            re = a->xR[offs + opidx * operandsize * 2 + 2 * k + 0];
            im = a->xR[offs + opidx * operandsize * 2 + 2 * k + 1];
            c = ae_cos(-2 * pi * k * i / n, _state);
            s = ae_sin(-2 * pi * k * i / n, _state);
            hre = hre + c * re - s * im;
            him = him + c * im + s * re;
         }
         buf->xR[2 * i + 0] = hre;
         buf->xR[2 * i + 1] = him;
      }
      for (i = 0; i < operandsize * 2; i++) {
         a->xR[offs + opidx * operandsize * 2 + i] = buf->xR[i];
      }
   }
}

// This subroutine applies complex codelet FFT to input/output array A.
//
// Inputs:
//     A           -   array, must be large enough for plan to work
//     Offs        -   offset of the subarray to process
//     OperandsCnt -   operands count (see description of FastTransformPlan)
//     OperandSize -   operand size (see description of FastTransformPlan)
//     MicrovectorSize-microvector size, must be 2
//
// Outputs:
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplycomplexcodeletfft(RVector *a, ae_int_t offs, ae_int_t operandscnt, ae_int_t operandsize, ae_int_t microvectorsize, ae_state *_state) {
   ae_int_t opidx;
   ae_int_t n;
   ae_int_t aoffset;
   double a0x;
   double a0y;
   double a1x;
   double a1y;
   double a2x;
   double a2y;
   double a3x;
   double a3y;
   double a4x;
   double a4y;
   double a5x;
   double a5y;
   double v0;
   double v1;
   double v2;
   double v3;
   double t1x;
   double t1y;
   double t2x;
   double t2y;
   double t3x;
   double t3y;
   double t4x;
   double t4y;
   double t5x;
   double t5y;
   double m1x;
   double m1y;
   double m2x;
   double m2y;
   double m3x;
   double m3y;
   double m4x;
   double m4y;
   double m5x;
   double m5y;
   double s1x;
   double s1y;
   double s2x;
   double s2y;
   double s3x;
   double s3y;
   double s4x;
   double s4y;
   double s5x;
   double s5y;
   double c1;
   double c2;
   double c3;
   double c4;
   double c5;
   double v;
   ae_assert(operandscnt >= 1, "FTApplyComplexCodeletFFT: OperandsCnt<1", _state);
   ae_assert(operandsize >= 1, "FTApplyComplexCodeletFFT: OperandSize<1", _state);
   ae_assert(microvectorsize == 2, "FTApplyComplexCodeletFFT: MicrovectorSize != 2", _state);
   n = operandsize;
// Hard-coded transforms for different N's
   ae_assert(n <= ftbase_maxradix, "FTApplyComplexCodeletFFT: N>MaxRadix", _state);
   if (n == 2) {
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset = offs + opidx * operandsize * 2;
         a0x = a->xR[aoffset + 0];
         a0y = a->xR[aoffset + 1];
         a1x = a->xR[aoffset + 2];
         a1y = a->xR[aoffset + 3];
         v0 = a0x + a1x;
         v1 = a0y + a1y;
         v2 = a0x - a1x;
         v3 = a0y - a1y;
         a->xR[aoffset + 0] = v0;
         a->xR[aoffset + 1] = v1;
         a->xR[aoffset + 2] = v2;
         a->xR[aoffset + 3] = v3;
      }
      return;
   }
   if (n == 3) {
      c1 = ae_cos(2 * pi / 3, _state) - 1;
      c2 = ae_sin(2 * pi / 3, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset = offs + opidx * operandsize * 2;
         a0x = a->xR[aoffset + 0];
         a0y = a->xR[aoffset + 1];
         a1x = a->xR[aoffset + 2];
         a1y = a->xR[aoffset + 3];
         a2x = a->xR[aoffset + 4];
         a2y = a->xR[aoffset + 5];
         t1x = a1x + a2x;
         t1y = a1y + a2y;
         a0x = a0x + t1x;
         a0y = a0y + t1y;
         m1x = c1 * t1x;
         m1y = c1 * t1y;
         m2x = c2 * (a1y - a2y);
         m2y = c2 * (a2x - a1x);
         s1x = a0x + m1x;
         s1y = a0y + m1y;
         a1x = s1x + m2x;
         a1y = s1y + m2y;
         a2x = s1x - m2x;
         a2y = s1y - m2y;
         a->xR[aoffset + 0] = a0x;
         a->xR[aoffset + 1] = a0y;
         a->xR[aoffset + 2] = a1x;
         a->xR[aoffset + 3] = a1y;
         a->xR[aoffset + 4] = a2x;
         a->xR[aoffset + 5] = a2y;
      }
      return;
   }
   if (n == 4) {
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset = offs + opidx * operandsize * 2;
         a0x = a->xR[aoffset + 0];
         a0y = a->xR[aoffset + 1];
         a1x = a->xR[aoffset + 2];
         a1y = a->xR[aoffset + 3];
         a2x = a->xR[aoffset + 4];
         a2y = a->xR[aoffset + 5];
         a3x = a->xR[aoffset + 6];
         a3y = a->xR[aoffset + 7];
         t1x = a0x + a2x;
         t1y = a0y + a2y;
         t2x = a1x + a3x;
         t2y = a1y + a3y;
         m2x = a0x - a2x;
         m2y = a0y - a2y;
         m3x = a1y - a3y;
         m3y = a3x - a1x;
         a->xR[aoffset + 0] = t1x + t2x;
         a->xR[aoffset + 1] = t1y + t2y;
         a->xR[aoffset + 4] = t1x - t2x;
         a->xR[aoffset + 5] = t1y - t2y;
         a->xR[aoffset + 2] = m2x + m3x;
         a->xR[aoffset + 3] = m2y + m3y;
         a->xR[aoffset + 6] = m2x - m3x;
         a->xR[aoffset + 7] = m2y - m3y;
      }
      return;
   }
   if (n == 5) {
      v = 2 * pi / 5;
      c1 = (ae_cos(v, _state) + ae_cos(2 * v, _state)) / 2 - 1;
      c2 = (ae_cos(v, _state) - ae_cos(2 * v, _state)) / 2;
      c3 = -ae_sin(v, _state);
      c4 = -(ae_sin(v, _state) + ae_sin(2 * v, _state));
      c5 = ae_sin(v, _state) - ae_sin(2 * v, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset = offs + opidx * operandsize * 2;
         t1x = a->xR[aoffset + 2] + a->xR[aoffset + 8];
         t1y = a->xR[aoffset + 3] + a->xR[aoffset + 9];
         t2x = a->xR[aoffset + 4] + a->xR[aoffset + 6];
         t2y = a->xR[aoffset + 5] + a->xR[aoffset + 7];
         t3x = a->xR[aoffset + 2] - a->xR[aoffset + 8];
         t3y = a->xR[aoffset + 3] - a->xR[aoffset + 9];
         t4x = a->xR[aoffset + 6] - a->xR[aoffset + 4];
         t4y = a->xR[aoffset + 7] - a->xR[aoffset + 5];
         t5x = t1x + t2x;
         t5y = t1y + t2y;
         a->xR[aoffset + 0] = a->xR[aoffset + 0] + t5x;
         a->xR[aoffset + 1] = a->xR[aoffset + 1] + t5y;
         m1x = c1 * t5x;
         m1y = c1 * t5y;
         m2x = c2 * (t1x - t2x);
         m2y = c2 * (t1y - t2y);
         m3x = -c3 * (t3y + t4y);
         m3y = c3 * (t3x + t4x);
         m4x = -c4 * t4y;
         m4y = c4 * t4x;
         m5x = -c5 * t3y;
         m5y = c5 * t3x;
         s3x = m3x - m4x;
         s3y = m3y - m4y;
         s5x = m3x + m5x;
         s5y = m3y + m5y;
         s1x = a->xR[aoffset + 0] + m1x;
         s1y = a->xR[aoffset + 1] + m1y;
         s2x = s1x + m2x;
         s2y = s1y + m2y;
         s4x = s1x - m2x;
         s4y = s1y - m2y;
         a->xR[aoffset + 2] = s2x + s3x;
         a->xR[aoffset + 3] = s2y + s3y;
         a->xR[aoffset + 4] = s4x + s5x;
         a->xR[aoffset + 5] = s4y + s5y;
         a->xR[aoffset + 6] = s4x - s5x;
         a->xR[aoffset + 7] = s4y - s5y;
         a->xR[aoffset + 8] = s2x - s3x;
         a->xR[aoffset + 9] = s2y - s3y;
      }
      return;
   }
   if (n == 6) {
      c1 = ae_cos(2 * pi / 3, _state) - 1;
      c2 = ae_sin(2 * pi / 3, _state);
      c3 = ae_cos(-pi / 3, _state);
      c4 = ae_sin(-pi / 3, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset = offs + opidx * operandsize * 2;
         a0x = a->xR[aoffset + 0];
         a0y = a->xR[aoffset + 1];
         a1x = a->xR[aoffset + 2];
         a1y = a->xR[aoffset + 3];
         a2x = a->xR[aoffset + 4];
         a2y = a->xR[aoffset + 5];
         a3x = a->xR[aoffset + 6];
         a3y = a->xR[aoffset + 7];
         a4x = a->xR[aoffset + 8];
         a4y = a->xR[aoffset + 9];
         a5x = a->xR[aoffset + 10];
         a5y = a->xR[aoffset + 11];
         v0 = a0x;
         v1 = a0y;
         a0x = a0x + a3x;
         a0y = a0y + a3y;
         a3x = v0 - a3x;
         a3y = v1 - a3y;
         v0 = a1x;
         v1 = a1y;
         a1x = a1x + a4x;
         a1y = a1y + a4y;
         a4x = v0 - a4x;
         a4y = v1 - a4y;
         v0 = a2x;
         v1 = a2y;
         a2x = a2x + a5x;
         a2y = a2y + a5y;
         a5x = v0 - a5x;
         a5y = v1 - a5y;
         t4x = a4x * c3 - a4y * c4;
         t4y = a4x * c4 + a4y * c3;
         a4x = t4x;
         a4y = t4y;
         t5x = -a5x * c3 - a5y * c4;
         t5y = a5x * c4 - a5y * c3;
         a5x = t5x;
         a5y = t5y;
         t1x = a1x + a2x;
         t1y = a1y + a2y;
         a0x = a0x + t1x;
         a0y = a0y + t1y;
         m1x = c1 * t1x;
         m1y = c1 * t1y;
         m2x = c2 * (a1y - a2y);
         m2y = c2 * (a2x - a1x);
         s1x = a0x + m1x;
         s1y = a0y + m1y;
         a1x = s1x + m2x;
         a1y = s1y + m2y;
         a2x = s1x - m2x;
         a2y = s1y - m2y;
         t1x = a4x + a5x;
         t1y = a4y + a5y;
         a3x = a3x + t1x;
         a3y = a3y + t1y;
         m1x = c1 * t1x;
         m1y = c1 * t1y;
         m2x = c2 * (a4y - a5y);
         m2y = c2 * (a5x - a4x);
         s1x = a3x + m1x;
         s1y = a3y + m1y;
         a4x = s1x + m2x;
         a4y = s1y + m2y;
         a5x = s1x - m2x;
         a5y = s1y - m2y;
         a->xR[aoffset + 0] = a0x;
         a->xR[aoffset + 1] = a0y;
         a->xR[aoffset + 2] = a3x;
         a->xR[aoffset + 3] = a3y;
         a->xR[aoffset + 4] = a1x;
         a->xR[aoffset + 5] = a1y;
         a->xR[aoffset + 6] = a4x;
         a->xR[aoffset + 7] = a4y;
         a->xR[aoffset + 8] = a2x;
         a->xR[aoffset + 9] = a2y;
         a->xR[aoffset + 10] = a5x;
         a->xR[aoffset + 11] = a5y;
      }
      return;
   }
}

// This subroutine applies complex "integrated" codelet FFT  to  input/output
// array A. "Integrated" codelet differs from "normal" one in following ways:
// * it can work with MicrovectorSize>1
// * hence, it can be used in Cooley-Tukey FFT without transpositions
// * it performs inlined multiplication by twiddle factors of Cooley-Tukey
//   FFT with N2=MicrovectorSize/2.
//
// Inputs:
//     A           -   array, must be large enough for plan to work
//     Offs        -   offset of the subarray to process
//     OperandsCnt -   operands count (see description of FastTransformPlan)
//     OperandSize -   operand size (see description of FastTransformPlan)
//     MicrovectorSize-microvector size, must be 1
//
// Outputs:
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplycomplexcodelettwfft(RVector *a, ae_int_t offs, ae_int_t operandscnt, ae_int_t operandsize, ae_int_t microvectorsize, ae_state *_state) {
   ae_int_t opidx;
   ae_int_t mvidx;
   ae_int_t n;
   ae_int_t m;
   ae_int_t aoffset0;
   ae_int_t aoffset2;
   ae_int_t aoffset4;
   ae_int_t aoffset6;
   ae_int_t aoffset8;
   ae_int_t aoffset10;
   double a0x;
   double a0y;
   double a1x;
   double a1y;
   double a2x;
   double a2y;
   double a3x;
   double a3y;
   double a4x;
   double a4y;
   double a5x;
   double a5y;
   double v0;
   double v1;
   double v2;
   double v3;
   double q0x;
   double q0y;
   double t1x;
   double t1y;
   double t2x;
   double t2y;
   double t3x;
   double t3y;
   double t4x;
   double t4y;
   double t5x;
   double t5y;
   double m1x;
   double m1y;
   double m2x;
   double m2y;
   double m3x;
   double m3y;
   double m4x;
   double m4y;
   double m5x;
   double m5y;
   double s1x;
   double s1y;
   double s2x;
   double s2y;
   double s3x;
   double s3y;
   double s4x;
   double s4y;
   double s5x;
   double s5y;
   double c1;
   double c2;
   double c3;
   double c4;
   double c5;
   double v;
   double tw0;
   double tw1;
   double twx;
   double twxm1;
   double twy;
   double tw2x;
   double tw2y;
   double tw3x;
   double tw3y;
   double tw4x;
   double tw4y;
   double tw5x;
   double tw5y;
   ae_assert(operandscnt >= 1, "FTApplyComplexCodeletFFT: OperandsCnt<1", _state);
   ae_assert(operandsize >= 1, "FTApplyComplexCodeletFFT: OperandSize<1", _state);
   ae_assert(microvectorsize >= 1, "FTApplyComplexCodeletFFT: MicrovectorSize != 1", _state);
   ae_assert(microvectorsize % 2 == 0, "FTApplyComplexCodeletFFT: MicrovectorSize is not even", _state);
   n = operandsize;
   m = microvectorsize / 2;
// Hard-coded transforms for different N's
   ae_assert(n <= ftbase_maxradix, "FTApplyComplexCodeletTwFFT: N>MaxRadix", _state);
   if (n == 2) {
      v = -2 * pi / (n * m);
      tw0 = -2 * sqr(ae_sin(0.5 * v, _state), _state);
      tw1 = ae_sin(v, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset0 = offs + opidx * operandsize * microvectorsize;
         aoffset2 = aoffset0 + microvectorsize;
         twxm1 = 0.0;
         twy = 0.0;
         for (mvidx = 0; mvidx < m; mvidx++) {
            a0x = a->xR[aoffset0];
            a0y = a->xR[aoffset0 + 1];
            a1x = a->xR[aoffset2];
            a1y = a->xR[aoffset2 + 1];
            v0 = a0x + a1x;
            v1 = a0y + a1y;
            v2 = a0x - a1x;
            v3 = a0y - a1y;
            a->xR[aoffset0] = v0;
            a->xR[aoffset0 + 1] = v1;
            a->xR[aoffset2] = v2 * (1 + twxm1) - v3 * twy;
            a->xR[aoffset2 + 1] = v3 * (1 + twxm1) + v2 * twy;
            aoffset0 = aoffset0 + 2;
            aoffset2 = aoffset2 + 2;
            if ((mvidx + 1) % ftbase_updatetw == 0) {
               v = -2 * pi * (mvidx + 1) / (n * m);
               twxm1 = ae_sin(0.5 * v, _state);
               twxm1 = -2 * twxm1 * twxm1;
               twy = ae_sin(v, _state);
            } else {
               v = twxm1 + tw0 + twxm1 * tw0 - twy * tw1;
               twy = twy + tw1 + twxm1 * tw1 + twy * tw0;
               twxm1 = v;
            }
         }
      }
      return;
   }
   if (n == 3) {
      v = -2 * pi / (n * m);
      tw0 = -2 * sqr(ae_sin(0.5 * v, _state), _state);
      tw1 = ae_sin(v, _state);
      c1 = ae_cos(2 * pi / 3, _state) - 1;
      c2 = ae_sin(2 * pi / 3, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset0 = offs + opidx * operandsize * microvectorsize;
         aoffset2 = aoffset0 + microvectorsize;
         aoffset4 = aoffset2 + microvectorsize;
         twx = 1.0;
         twxm1 = 0.0;
         twy = 0.0;
         for (mvidx = 0; mvidx < m; mvidx++) {
            a0x = a->xR[aoffset0];
            a0y = a->xR[aoffset0 + 1];
            a1x = a->xR[aoffset2];
            a1y = a->xR[aoffset2 + 1];
            a2x = a->xR[aoffset4];
            a2y = a->xR[aoffset4 + 1];
            t1x = a1x + a2x;
            t1y = a1y + a2y;
            a0x = a0x + t1x;
            a0y = a0y + t1y;
            m1x = c1 * t1x;
            m1y = c1 * t1y;
            m2x = c2 * (a1y - a2y);
            m2y = c2 * (a2x - a1x);
            s1x = a0x + m1x;
            s1y = a0y + m1y;
            a1x = s1x + m2x;
            a1y = s1y + m2y;
            a2x = s1x - m2x;
            a2y = s1y - m2y;
            tw2x = twx * twx - twy * twy;
            tw2y = 2 * twx * twy;
            a->xR[aoffset0] = a0x;
            a->xR[aoffset0 + 1] = a0y;
            a->xR[aoffset2] = a1x * twx - a1y * twy;
            a->xR[aoffset2 + 1] = a1y * twx + a1x * twy;
            a->xR[aoffset4] = a2x * tw2x - a2y * tw2y;
            a->xR[aoffset4 + 1] = a2y * tw2x + a2x * tw2y;
            aoffset0 = aoffset0 + 2;
            aoffset2 = aoffset2 + 2;
            aoffset4 = aoffset4 + 2;
            if ((mvidx + 1) % ftbase_updatetw == 0) {
               v = -2 * pi * (mvidx + 1) / (n * m);
               twxm1 = ae_sin(0.5 * v, _state);
               twxm1 = -2 * twxm1 * twxm1;
               twy = ae_sin(v, _state);
               twx = twxm1 + 1;
            } else {
               v = twxm1 + tw0 + twxm1 * tw0 - twy * tw1;
               twy = twy + tw1 + twxm1 * tw1 + twy * tw0;
               twxm1 = v;
               twx = v + 1;
            }
         }
      }
      return;
   }
   if (n == 4) {
      v = -2 * pi / (n * m);
      tw0 = -2 * sqr(ae_sin(0.5 * v, _state), _state);
      tw1 = ae_sin(v, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset0 = offs + opidx * operandsize * microvectorsize;
         aoffset2 = aoffset0 + microvectorsize;
         aoffset4 = aoffset2 + microvectorsize;
         aoffset6 = aoffset4 + microvectorsize;
         twx = 1.0;
         twxm1 = 0.0;
         twy = 0.0;
         for (mvidx = 0; mvidx < m; mvidx++) {
            a0x = a->xR[aoffset0];
            a0y = a->xR[aoffset0 + 1];
            a1x = a->xR[aoffset2];
            a1y = a->xR[aoffset2 + 1];
            a2x = a->xR[aoffset4];
            a2y = a->xR[aoffset4 + 1];
            a3x = a->xR[aoffset6];
            a3y = a->xR[aoffset6 + 1];
            t1x = a0x + a2x;
            t1y = a0y + a2y;
            t2x = a1x + a3x;
            t2y = a1y + a3y;
            m2x = a0x - a2x;
            m2y = a0y - a2y;
            m3x = a1y - a3y;
            m3y = a3x - a1x;
            tw2x = twx * twx - twy * twy;
            tw2y = 2 * twx * twy;
            tw3x = twx * tw2x - twy * tw2y;
            tw3y = twx * tw2y + twy * tw2x;
            a1x = m2x + m3x;
            a1y = m2y + m3y;
            a2x = t1x - t2x;
            a2y = t1y - t2y;
            a3x = m2x - m3x;
            a3y = m2y - m3y;
            a->xR[aoffset0] = t1x + t2x;
            a->xR[aoffset0 + 1] = t1y + t2y;
            a->xR[aoffset2] = a1x * twx - a1y * twy;
            a->xR[aoffset2 + 1] = a1y * twx + a1x * twy;
            a->xR[aoffset4] = a2x * tw2x - a2y * tw2y;
            a->xR[aoffset4 + 1] = a2y * tw2x + a2x * tw2y;
            a->xR[aoffset6] = a3x * tw3x - a3y * tw3y;
            a->xR[aoffset6 + 1] = a3y * tw3x + a3x * tw3y;
            aoffset0 = aoffset0 + 2;
            aoffset2 = aoffset2 + 2;
            aoffset4 = aoffset4 + 2;
            aoffset6 = aoffset6 + 2;
            if ((mvidx + 1) % ftbase_updatetw == 0) {
               v = -2 * pi * (mvidx + 1) / (n * m);
               twxm1 = ae_sin(0.5 * v, _state);
               twxm1 = -2 * twxm1 * twxm1;
               twy = ae_sin(v, _state);
               twx = twxm1 + 1;
            } else {
               v = twxm1 + tw0 + twxm1 * tw0 - twy * tw1;
               twy = twy + tw1 + twxm1 * tw1 + twy * tw0;
               twxm1 = v;
               twx = v + 1;
            }
         }
      }
      return;
   }
   if (n == 5) {
      v = -2 * pi / (n * m);
      tw0 = -2 * sqr(ae_sin(0.5 * v, _state), _state);
      tw1 = ae_sin(v, _state);
      v = 2 * pi / 5;
      c1 = (ae_cos(v, _state) + ae_cos(2 * v, _state)) / 2 - 1;
      c2 = (ae_cos(v, _state) - ae_cos(2 * v, _state)) / 2;
      c3 = -ae_sin(v, _state);
      c4 = -(ae_sin(v, _state) + ae_sin(2 * v, _state));
      c5 = ae_sin(v, _state) - ae_sin(2 * v, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset0 = offs + opidx * operandsize * microvectorsize;
         aoffset2 = aoffset0 + microvectorsize;
         aoffset4 = aoffset2 + microvectorsize;
         aoffset6 = aoffset4 + microvectorsize;
         aoffset8 = aoffset6 + microvectorsize;
         twx = 1.0;
         twxm1 = 0.0;
         twy = 0.0;
         for (mvidx = 0; mvidx < m; mvidx++) {
            a0x = a->xR[aoffset0];
            a0y = a->xR[aoffset0 + 1];
            a1x = a->xR[aoffset2];
            a1y = a->xR[aoffset2 + 1];
            a2x = a->xR[aoffset4];
            a2y = a->xR[aoffset4 + 1];
            a3x = a->xR[aoffset6];
            a3y = a->xR[aoffset6 + 1];
            a4x = a->xR[aoffset8];
            a4y = a->xR[aoffset8 + 1];
            t1x = a1x + a4x;
            t1y = a1y + a4y;
            t2x = a2x + a3x;
            t2y = a2y + a3y;
            t3x = a1x - a4x;
            t3y = a1y - a4y;
            t4x = a3x - a2x;
            t4y = a3y - a2y;
            t5x = t1x + t2x;
            t5y = t1y + t2y;
            q0x = a0x + t5x;
            q0y = a0y + t5y;
            m1x = c1 * t5x;
            m1y = c1 * t5y;
            m2x = c2 * (t1x - t2x);
            m2y = c2 * (t1y - t2y);
            m3x = -c3 * (t3y + t4y);
            m3y = c3 * (t3x + t4x);
            m4x = -c4 * t4y;
            m4y = c4 * t4x;
            m5x = -c5 * t3y;
            m5y = c5 * t3x;
            s3x = m3x - m4x;
            s3y = m3y - m4y;
            s5x = m3x + m5x;
            s5y = m3y + m5y;
            s1x = q0x + m1x;
            s1y = q0y + m1y;
            s2x = s1x + m2x;
            s2y = s1y + m2y;
            s4x = s1x - m2x;
            s4y = s1y - m2y;
            tw2x = twx * twx - twy * twy;
            tw2y = 2 * twx * twy;
            tw3x = twx * tw2x - twy * tw2y;
            tw3y = twx * tw2y + twy * tw2x;
            tw4x = tw2x * tw2x - tw2y * tw2y;
            tw4y = tw2x * tw2y + tw2y * tw2x;
            a1x = s2x + s3x;
            a1y = s2y + s3y;
            a2x = s4x + s5x;
            a2y = s4y + s5y;
            a3x = s4x - s5x;
            a3y = s4y - s5y;
            a4x = s2x - s3x;
            a4y = s2y - s3y;
            a->xR[aoffset0] = q0x;
            a->xR[aoffset0 + 1] = q0y;
            a->xR[aoffset2] = a1x * twx - a1y * twy;
            a->xR[aoffset2 + 1] = a1x * twy + a1y * twx;
            a->xR[aoffset4] = a2x * tw2x - a2y * tw2y;
            a->xR[aoffset4 + 1] = a2x * tw2y + a2y * tw2x;
            a->xR[aoffset6] = a3x * tw3x - a3y * tw3y;
            a->xR[aoffset6 + 1] = a3x * tw3y + a3y * tw3x;
            a->xR[aoffset8] = a4x * tw4x - a4y * tw4y;
            a->xR[aoffset8 + 1] = a4x * tw4y + a4y * tw4x;
            aoffset0 = aoffset0 + 2;
            aoffset2 = aoffset2 + 2;
            aoffset4 = aoffset4 + 2;
            aoffset6 = aoffset6 + 2;
            aoffset8 = aoffset8 + 2;
            if ((mvidx + 1) % ftbase_updatetw == 0) {
               v = -2 * pi * (mvidx + 1) / (n * m);
               twxm1 = ae_sin(0.5 * v, _state);
               twxm1 = -2 * twxm1 * twxm1;
               twy = ae_sin(v, _state);
               twx = twxm1 + 1;
            } else {
               v = twxm1 + tw0 + twxm1 * tw0 - twy * tw1;
               twy = twy + tw1 + twxm1 * tw1 + twy * tw0;
               twxm1 = v;
               twx = v + 1;
            }
         }
      }
      return;
   }
   if (n == 6) {
      c1 = ae_cos(2 * pi / 3, _state) - 1;
      c2 = ae_sin(2 * pi / 3, _state);
      c3 = ae_cos(-pi / 3, _state);
      c4 = ae_sin(-pi / 3, _state);
      v = -2 * pi / (n * m);
      tw0 = -2 * sqr(ae_sin(0.5 * v, _state), _state);
      tw1 = ae_sin(v, _state);
      for (opidx = 0; opidx < operandscnt; opidx++) {
         aoffset0 = offs + opidx * operandsize * microvectorsize;
         aoffset2 = aoffset0 + microvectorsize;
         aoffset4 = aoffset2 + microvectorsize;
         aoffset6 = aoffset4 + microvectorsize;
         aoffset8 = aoffset6 + microvectorsize;
         aoffset10 = aoffset8 + microvectorsize;
         twx = 1.0;
         twxm1 = 0.0;
         twy = 0.0;
         for (mvidx = 0; mvidx < m; mvidx++) {
            a0x = a->xR[aoffset0 + 0];
            a0y = a->xR[aoffset0 + 1];
            a1x = a->xR[aoffset2 + 0];
            a1y = a->xR[aoffset2 + 1];
            a2x = a->xR[aoffset4 + 0];
            a2y = a->xR[aoffset4 + 1];
            a3x = a->xR[aoffset6 + 0];
            a3y = a->xR[aoffset6 + 1];
            a4x = a->xR[aoffset8 + 0];
            a4y = a->xR[aoffset8 + 1];
            a5x = a->xR[aoffset10 + 0];
            a5y = a->xR[aoffset10 + 1];
            v0 = a0x;
            v1 = a0y;
            a0x = a0x + a3x;
            a0y = a0y + a3y;
            a3x = v0 - a3x;
            a3y = v1 - a3y;
            v0 = a1x;
            v1 = a1y;
            a1x = a1x + a4x;
            a1y = a1y + a4y;
            a4x = v0 - a4x;
            a4y = v1 - a4y;
            v0 = a2x;
            v1 = a2y;
            a2x = a2x + a5x;
            a2y = a2y + a5y;
            a5x = v0 - a5x;
            a5y = v1 - a5y;
            t4x = a4x * c3 - a4y * c4;
            t4y = a4x * c4 + a4y * c3;
            a4x = t4x;
            a4y = t4y;
            t5x = -a5x * c3 - a5y * c4;
            t5y = a5x * c4 - a5y * c3;
            a5x = t5x;
            a5y = t5y;
            t1x = a1x + a2x;
            t1y = a1y + a2y;
            a0x = a0x + t1x;
            a0y = a0y + t1y;
            m1x = c1 * t1x;
            m1y = c1 * t1y;
            m2x = c2 * (a1y - a2y);
            m2y = c2 * (a2x - a1x);
            s1x = a0x + m1x;
            s1y = a0y + m1y;
            a1x = s1x + m2x;
            a1y = s1y + m2y;
            a2x = s1x - m2x;
            a2y = s1y - m2y;
            t1x = a4x + a5x;
            t1y = a4y + a5y;
            a3x = a3x + t1x;
            a3y = a3y + t1y;
            m1x = c1 * t1x;
            m1y = c1 * t1y;
            m2x = c2 * (a4y - a5y);
            m2y = c2 * (a5x - a4x);
            s1x = a3x + m1x;
            s1y = a3y + m1y;
            a4x = s1x + m2x;
            a4y = s1y + m2y;
            a5x = s1x - m2x;
            a5y = s1y - m2y;
            tw2x = twx * twx - twy * twy;
            tw2y = 2 * twx * twy;
            tw3x = twx * tw2x - twy * tw2y;
            tw3y = twx * tw2y + twy * tw2x;
            tw4x = tw2x * tw2x - tw2y * tw2y;
            tw4y = 2 * tw2x * tw2y;
            tw5x = tw3x * tw2x - tw3y * tw2y;
            tw5y = tw3x * tw2y + tw3y * tw2x;
            a->xR[aoffset0 + 0] = a0x;
            a->xR[aoffset0 + 1] = a0y;
            a->xR[aoffset2 + 0] = a3x * twx - a3y * twy;
            a->xR[aoffset2 + 1] = a3y * twx + a3x * twy;
            a->xR[aoffset4 + 0] = a1x * tw2x - a1y * tw2y;
            a->xR[aoffset4 + 1] = a1y * tw2x + a1x * tw2y;
            a->xR[aoffset6 + 0] = a4x * tw3x - a4y * tw3y;
            a->xR[aoffset6 + 1] = a4y * tw3x + a4x * tw3y;
            a->xR[aoffset8 + 0] = a2x * tw4x - a2y * tw4y;
            a->xR[aoffset8 + 1] = a2y * tw4x + a2x * tw4y;
            a->xR[aoffset10 + 0] = a5x * tw5x - a5y * tw5y;
            a->xR[aoffset10 + 1] = a5y * tw5x + a5x * tw5y;
            aoffset0 = aoffset0 + 2;
            aoffset2 = aoffset2 + 2;
            aoffset4 = aoffset4 + 2;
            aoffset6 = aoffset6 + 2;
            aoffset8 = aoffset8 + 2;
            aoffset10 = aoffset10 + 2;
            if ((mvidx + 1) % ftbase_updatetw == 0) {
               v = -2 * pi * (mvidx + 1) / (n * m);
               twxm1 = ae_sin(0.5 * v, _state);
               twxm1 = -2 * twxm1 * twxm1;
               twy = ae_sin(v, _state);
               twx = twxm1 + 1;
            } else {
               v = twxm1 + tw0 + twxm1 * tw0 - twy * tw1;
               twy = twy + tw1 + twxm1 * tw1 + twy * tw0;
               twxm1 = v;
               twx = v + 1;
            }
         }
      }
      return;
   }
}

// Twiddle factors calculation
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_ffttwcalc(RVector *a, ae_int_t aoffset, ae_int_t n1, ae_int_t n2, ae_state *_state) {
   ae_int_t i;
   ae_int_t j2;
   ae_int_t n;
   ae_int_t halfn1;
   ae_int_t offs;
   double x;
   double y;
   double twxm1;
   double twy;
   double twbasexm1;
   double twbasey;
   double twrowxm1;
   double twrowy;
   double tmpx;
   double tmpy;
   double v;
   ae_int_t updatetw2;
// Multiplication by twiddle factors for complex Cooley-Tukey FFT
// with N factorized as N1*N2.
//
// Naive solution to this problem is given below:
//
//     > for K:=1 to N2-1 do
//     >     for J:=1 to N1-1 do
//     >     begin
//     >         Idx:=K*N1+J;
//     >         X:=A[AOffset+2*Idx+0];
//     >         Y:=A[AOffset+2*Idx+1];
//     >         TwX:=Cos(-2*Pi()*K*J/(N1*N2));
//     >         TwY:=Sin(-2*Pi()*K*J/(N1*N2));
//     >         A[AOffset+2*Idx+0]:=X*TwX-Y*TwY;
//     >         A[AOffset+2*Idx+1]:=X*TwY+Y*TwX;
//     >     end;
//
// However, there are exist more efficient solutions.
//
// Each pass of the inner cycle corresponds to multiplication of one
// entry of A by W[k,j]=exp(-I*2*pi*k*j/N). This factor can be rewritten
// as exp(-I*2*pi*k/N)^j. So we can replace costly exponentiation by
// repeated multiplication: W[k,j+1]=W[k,j]*exp(-I*2*pi*k/N), with
// second factor being computed once in the beginning of the iteration.
//
// Also, exp(-I*2*pi*k/N) can be represented as exp(-I*2*pi/N)^k, i.e.
// we have W[K+1,1]=W[K,1]*W[1,1].
//
// In our loop we use following variables:
// * [TwBaseXM1,TwBaseY] =   [cos(2*pi/N)-1,     sin(2*pi/N)]
// * [TwRowXM1, TwRowY]  =   [cos(2*pi*I/N)-1,   sin(2*pi*I/N)]
// * [TwXM1,    TwY]     =   [cos(2*pi*I*J/N)-1, sin(2*pi*I*J/N)]
//
// Meaning of the variables:
// * [TwXM1,TwY] is current twiddle factor W[I,J]
// * [TwRowXM1, TwRowY] is W[I,1]
// * [TwBaseXM1,TwBaseY] is W[1,1]
//
// During inner loop we multiply current twiddle factor by W[I,1],
// during outer loop we update W[I,1].
//
   ae_assert(ftbase_updatetw >= 2, "FFTTwCalc: internal error - UpdateTw<2", _state);
   updatetw2 = ftbase_updatetw / 2;
   halfn1 = n1 / 2;
   n = n1 * n2;
   v = -2 * pi / n;
   twbasexm1 = -2 * sqr(ae_sin(0.5 * v, _state), _state);
   twbasey = ae_sin(v, _state);
   twrowxm1 = 0.0;
   twrowy = 0.0;
   offs = aoffset;
   for (i = 0; i < n2; i++) {
   // Initialize twiddle factor for current row
      twxm1 = 0.0;
      twy = 0.0;
   // N1-point block is separated into 2-point chunks and residual 1-point chunk
   // (in case N1 is odd). Unrolled loop is several times faster.
      for (j2 = 0; j2 < halfn1; j2++) {
      // Processing:
      // * process first element in a chunk.
      // * update twiddle factor (unconditional update)
      // * process second element
      // * conditional update of the twiddle factor
         x = a->xR[offs + 0];
         y = a->xR[offs + 1];
         tmpx = x * (1 + twxm1) - y * twy;
         tmpy = x * twy + y * (1 + twxm1);
         a->xR[offs + 0] = tmpx;
         a->xR[offs + 1] = tmpy;
         tmpx = (1 + twxm1) * twrowxm1 - twy * twrowy;
         twy = twy + (1 + twxm1) * twrowy + twy * twrowxm1;
         twxm1 = twxm1 + tmpx;
         x = a->xR[offs + 2];
         y = a->xR[offs + 3];
         tmpx = x * (1 + twxm1) - y * twy;
         tmpy = x * twy + y * (1 + twxm1);
         a->xR[offs + 2] = tmpx;
         a->xR[offs + 3] = tmpy;
         offs = offs + 4;
         if ((j2 + 1) % updatetw2 == 0 && j2 < halfn1 - 1) {
         // Recalculate twiddle factor
            v = -2 * pi * i * 2 * (j2 + 1) / n;
            twxm1 = ae_sin(0.5 * v, _state);
            twxm1 = -2 * twxm1 * twxm1;
            twy = ae_sin(v, _state);
         } else {
         // Update twiddle factor
            tmpx = (1 + twxm1) * twrowxm1 - twy * twrowy;
            twy = twy + (1 + twxm1) * twrowy + twy * twrowxm1;
            twxm1 = twxm1 + tmpx;
         }
      }
      if (n1 % 2 == 1) {
      // Handle residual chunk
         x = a->xR[offs + 0];
         y = a->xR[offs + 1];
         tmpx = x * (1 + twxm1) - y * twy;
         tmpy = x * twy + y * (1 + twxm1);
         a->xR[offs + 0] = tmpx;
         a->xR[offs + 1] = tmpy;
         offs = offs + 2;
      }
   // update TwRow: TwRow(new) = TwRow(old)*TwBase
      if (i < n2 - 1) {
         if ((i + 1) % ftbase_updatetw == 0) {
            v = -2 * pi * (i + 1) / n;
            twrowxm1 = ae_sin(0.5 * v, _state);
            twrowxm1 = -2 * twrowxm1 * twrowxm1;
            twrowy = ae_sin(v, _state);
         } else {
            tmpx = twbasexm1 + twrowxm1 * twbasexm1 - twrowy * twbasey;
            tmpy = twbasey + twrowxm1 * twbasey + twrowy * twbasexm1;
            twrowxm1 = twrowxm1 + tmpx;
            twrowy = twrowy + tmpy;
         }
      }
   }
}

// Returns estimate of FLOP count for the FFT.
//
// It is only an estimate based on operations count for the PERFECT FFT
// and relative inefficiency of the algorithm actually used.
//
// N should be power of 2, estimates are badly wrong for non-power-of-2 N's.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
double ftbasegetflopestimate(ae_int_t n, ae_state *_state) {
   double result;
   result = ftbase_ftbaseinefficiencyfactor * (4 * n * ae_log((double)n, _state) / ae_log(2.0, _state) - 6 * n + 8);
   return result;
}

#if 0
// Recurrent subroutine for a (non-existent) InternalRealLinTranspose
//
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_fftirltrec(RVector *a, ae_int_t astart, ae_int_t astride, RVector *b, ae_int_t bstart, ae_int_t bstride, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t idx1;
   ae_int_t idx2;
   ae_int_t m1;
   ae_int_t n1;
   if (m == 0 || n == 0) {
      return;
   }
   if (maxint(m, n, _state) <= 8) {
      for (i = 0; i < m; i++) {
         idx1 = bstart + i;
         idx2 = astart + i * astride;
         for (j = 0; j < n; j++) {
            b->xR[idx1] = a->xR[idx2];
            idx1 = idx1 + bstride;
            idx2 = idx2 + 1;
         }
      }
      return;
   }
   if (n > m) {
   // New partition:
   //
   // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
   //                                  ( B2 )
      n1 = n / 2;
      if (n - n1 >= 8 && n1 % 8 != 0) {
         n1 = n1 + (8 - n1 % 8);
      }
      ae_assert(n - n1 > 0, "Assertion failed", _state);
      ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
      ftbase_fftirltrec(a, astart + n1, astride, b, bstart + n1 * bstride, bstride, m, n - n1, _state);
   } else {
   // New partition:
   //
   // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
   //                     ( A2 )
      m1 = m / 2;
      if (m - m1 >= 8 && m1 % 8 != 0) {
         m1 = m1 + (8 - m1 % 8);
      }
      ae_assert(m - m1 > 0, "Assertion failed", _state);
      ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
      ftbase_fftirltrec(a, astart + m1 * astride, astride, b, bstart + m1, bstride, m - m1, n, _state);
   }
}
#endif

// Recurrent subroutine for a InternalComplexLinTranspose
//
// Write A^T to B, where:
// * A is m*n complex matrix stored in array A as pairs of real/image values,
//   beginning from AStart position, with AStride stride
// * B is n*m complex matrix stored in array B as pairs of real/image values,
//   beginning from BStart position, with BStride stride
// stride is measured in complex numbers, i.e. in real/image pairs.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_ffticltrec(RVector *a, ae_int_t astart, ae_int_t astride, RVector *b, ae_int_t bstart, ae_int_t bstride, ae_int_t m, ae_int_t n, ae_state *_state) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t idx1;
   ae_int_t idx2;
   ae_int_t m2;
   ae_int_t m1;
   ae_int_t n1;
   if (m == 0 || n == 0) {
      return;
   }
   if (maxint(m, n, _state) <= 8) {
      m2 = 2 * bstride;
      for (i = 0; i < m; i++) {
         idx1 = bstart + 2 * i;
         idx2 = astart + 2 * i * astride;
         for (j = 0; j < n; j++) {
            b->xR[idx1 + 0] = a->xR[idx2 + 0];
            b->xR[idx1 + 1] = a->xR[idx2 + 1];
            idx1 = idx1 + m2;
            idx2 = idx2 + 2;
         }
      }
      return;
   }
   if (n > m) {
   // New partition:
   //
   // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
   //                                  ( B2 )
      n1 = n / 2;
      if (n - n1 >= 8 && n1 % 8 != 0) {
         n1 = n1 + (8 - n1 % 8);
      }
      ae_assert(n - n1 > 0, "Assertion failed", _state);
      ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
      ftbase_ffticltrec(a, astart + 2 * n1, astride, b, bstart + 2 * n1 * bstride, bstride, m, n - n1, _state);
   } else {
   // New partition:
   //
   // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
   //                     ( A2 )
      m1 = m / 2;
      if (m - m1 >= 8 && m1 % 8 != 0) {
         m1 = m1 + (8 - m1 % 8);
      }
      ae_assert(m - m1 > 0, "Assertion failed", _state);
      ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
      ftbase_ffticltrec(a, astart + 2 * m1 * astride, astride, b, bstart + 2 * m1, bstride, m - m1, n, _state);
   }
}

// Linear transpose: transpose complex matrix stored in 1-dimensional array
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_internalcomplexlintranspose(RVector *a, ae_int_t m, ae_int_t n, ae_int_t astart, RVector *buf, ae_state *_state) {
   ftbase_ffticltrec(a, astart, n, buf, 0, m, m, n, _state);
   ae_v_move(&a->xR[astart], 1, buf->xR, 1, 2 * m * n);
}

// Factorizes task size N into product of two smaller sizes N1 and N2
//
// Inputs:
//     N       -   task size, N>0
//     IsRoot  -   whether taks is root task (first one in a sequence)
//
// Outputs:
//     N1, N2  -   such numbers that:
//                 * for prime N:                  N1=N2=0
//                 * for composite N <= MaxRadix:    N1=N2=0
//                 * for composite N>MaxRadix:     1 <= N1 <= N2, N1*N2=N
// ALGLIB: Copyright 08.04.2013 by Sergey Bochkanov
static void ftbase_ftfactorize(ae_int_t n, bool isroot, ae_int_t *n1, ae_int_t *n2, ae_state *_state) {
   ae_int_t j;
   ae_int_t k;
   *n1 = 0;
   *n2 = 0;
   ae_assert(n > 0, "FTFactorize: N <= 0", _state);
   *n1 = 0;
   *n2 = 0;
// Small N
   if (n <= ftbase_maxradix) {
      return;
   }
// Large N, recursive split
   if (n > ftbase_recursivethreshold) {
      k = iceil(ae_sqrt((double)n, _state), _state) + 1;
      ae_assert(k * k >= n, "FTFactorize: internal error during recursive factorization", _state);
      for (j = k; j >= 2; j--) {
         if (n % j == 0) {
            *n1 = minint(n / j, j, _state);
            *n2 = maxint(n / j, j, _state);
            return;
         }
      }
   }
// N>MaxRadix, try to find good codelet
   for (j = ftbase_maxradix; j >= 2; j--) {
      if (n % j == 0) {
         *n1 = j;
         *n2 = n / j;
         break;
      }
   }
// In case no good codelet was found,
// try to factorize N into product of ANY primes.
   if (*n1 * (*n2) != n) {
      for (j = 2; j < n; j++) {
         if (n % j == 0) {
            *n1 = j;
            *n2 = n / j;
            break;
         }
         if (j * j > n) {
            break;
         }
      }
   }
// normalize
   if (*n1 > (*n2)) {
      j = *n1;
      *n1 = *n2;
      *n2 = j;
   }
}

// Returns good factorization N=N1*N2.
//
// Usually N1 <= N2 (but not always - small N's may be exception).
// if N1 != 1 then N2 != 1.
//
// Factorization is chosen depending on task type and codelets we have.
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
void ftbasefactorize(ae_int_t n, ae_int_t tasktype, ae_int_t *n1, ae_int_t *n2, ae_state *_state) {
   ae_int_t j;
   *n1 = 0;
   *n2 = 0;
   *n1 = 0;
   *n2 = 0;
// try to find good codelet
   if (*n1 * (*n2) != n) {
      for (j = ftbase_ftbasecodeletrecommended; j >= 2; j--) {
         if (n % j == 0) {
            *n1 = j;
            *n2 = n / j;
            break;
         }
      }
   }
// try to factorize N
   if (*n1 * (*n2) != n) {
      for (j = ftbase_ftbasecodeletrecommended + 1; j < n; j++) {
         if (n % j == 0) {
            *n1 = j;
            *n2 = n / j;
            break;
         }
      }
   }
// looks like N is prime :(
   if (*n1 * (*n2) != n) {
      *n1 = 1;
      *n2 = n;
   }
// normalize
   if (*n2 == 1 && *n1 != 1) {
      *n2 = *n1;
      *n1 = 1;
   }
}

// Is number smooth?
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
bool ftbaseissmooth(ae_int_t n, ae_state *_state) {
   ae_int_t i;
   bool result;
   for (i = 2; i <= ftbase_ftbasemaxsmoothfactor; i++) {
      while (n % i == 0) {
         n = n / i;
      }
   }
   result = n == 1;
   return result;
}

// recurrent subroutine for FFTFindSmoothRec
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
static void ftbase_ftbasefindsmoothrec(ae_int_t n, ae_int_t seed, ae_int_t leastfactor, ae_int_t *best, ae_state *_state) {
   ae_assert(ftbase_ftbasemaxsmoothfactor <= 5, "FTBaseFindSmoothRec: internal error!", _state);
   if (seed >= n) {
      *best = minint(*best, seed, _state);
      return;
   }
   if (leastfactor <= 2) {
      ftbase_ftbasefindsmoothrec(n, seed * 2, 2, best, _state);
   }
   if (leastfactor <= 3) {
      ftbase_ftbasefindsmoothrec(n, seed * 3, 3, best, _state);
   }
   if (leastfactor <= 5) {
      ftbase_ftbasefindsmoothrec(n, seed * 5, 5, best, _state);
   }
}

// Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
// than or equal to max(N,2)
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
ae_int_t ftbasefindsmooth(ae_int_t n, ae_state *_state) {
   ae_int_t best;
   ae_int_t result;
   best = 2;
   while (best < n) {
      best = 2 * best;
   }
   ftbase_ftbasefindsmoothrec(n, 1, 2, &best, _state);
   result = best;
   return result;
}

// Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
// greater than or equal to max(N,2)
// ALGLIB: Copyright 01.05.2009 by Sergey Bochkanov
ae_int_t ftbasefindsmootheven(ae_int_t n, ae_state *_state) {
   ae_int_t best;
   ae_int_t result;
   best = 2;
   while (best < n) {
      best = 2 * best;
   }
   ftbase_ftbasefindsmoothrec(n, 2, 2, &best, _state);
   result = best;
   return result;
}

// This function returns EXACT estimate of the space requirements for N-point
// FFT. Internals of this function are highly dependent on details of different
// FFTs employed by this unit, so every time algorithm is changed this function
// has to be rewritten.
//
// Inputs:
//     N           -   transform length
//     PrecRSize   -   must be set to zero
//     PrecISize   -   must be set to zero
//
// Outputs:
//     PrecRSize   -   number of real temporaries required for transformation
//     PrecISize   -   number of integer temporaries required for transformation
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftdeterminespacerequirements(ae_int_t n, ae_int_t *precrsize, ae_int_t *precisize, ae_state *_state) {
   ae_int_t ncur;
   ae_int_t f;
   ae_int_t i;
// Determine required sizes of precomputed real and integer
// buffers. This stage of code is highly dependent on internals
// of FTComplexFFTPlanRec() and must be kept synchronized with
// possible changes in internals of plan generation function.
//
// Buffer size is determined as follows:
// * N is factorized
// * we factor out anything which is less or equal to MaxRadix
// * prime factor F>RaderThreshold requires 4*FTBaseFindSmooth(2*F-1)
//   real entries to store precomputed Quantities for Bluestein's
//   transformation
// * prime factor F <= RaderThreshold requires 2*(F-1)+ESTIMATE(F-1)
//   precomputed storage
   ncur = n;
   for (i = 2; i <= ftbase_maxradix; i++) {
      while (ncur % i == 0) {
         ncur = ncur / i;
      }
   }
   f = 2;
   while (f <= ncur) {
      while (ncur % f == 0) {
         if (f > ftbase_raderthreshold) {
            *precrsize = *precrsize + 4 * ftbasefindsmooth(2 * f - 1, _state);
         } else {
            *precrsize = *precrsize + 2 * (f - 1);
            ftbase_ftdeterminespacerequirements(f - 1, precrsize, precisize, _state);
         }
         ncur = ncur / f;
      }
      f = f + 1;
   }
}

// Forward reference to an indirect recursive call.
static void ftbase_ftapplysubplan(fasttransformplan *plan, ae_int_t subplan, RVector *a, ae_int_t abase, ae_int_t aoffset, RVector *buf, ae_int_t repcnt, ae_state *_state);

// This subroutine applies complex Bluestein's FFT to input/output array A.
//
// Inputs:
//     Plan        -   transformation plan
//     A           -   array, must be large enough for plan to work
//     ABase       -   base offset in array A, this value points to start of
//                     subarray whose length is equal to length of the plan
//     AOffset     -   offset with respect to ABase, 0 <= AOffset<PlanLength.
//                     This is an offset within large PlanLength-subarray of
//                     the chunk to process.
//     OperandsCnt -   number of repeated operands (length N each)
//     N           -   original data length (measured in complex numbers)
//     M           -   padded data length (measured in complex numbers)
//     PrecOffs    -   offset of the precomputed data for the plan
//     SubPlan     -   position of the length-M FFT subplan which is used by
//                     transformation
//     BufA        -   temporary buffer, at least 2*M elements
//     BufB        -   temporary buffer, at least 2*M elements
//     BufC        -   temporary buffer, at least 2*M elements
//     BufD        -   temporary buffer, at least 2*M elements
//
// Outputs:
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftbluesteinsfft(fasttransformplan *plan, RVector *a, ae_int_t abase, ae_int_t aoffset, ae_int_t operandscnt, ae_int_t n, ae_int_t m, ae_int_t precoffs, ae_int_t subplan, RVector *bufa, RVector *bufb, RVector *bufc, RVector *bufd, ae_state *_state) {
   ae_int_t op;
   ae_int_t i;
   double x;
   double y;
   double bx;
   double by;
   double ax;
   double ay;
   double rx;
   double ry;
   ae_int_t p0;
   ae_int_t p1;
   ae_int_t p2;
   for (op = 0; op < operandscnt; op++) {
   // Multiply A by conj(Z), store to buffer.
   // Pad A by zeros.
   //
   // NOTE: Z[k]=exp(i*pi*k^2/N)
      p0 = abase + aoffset + op * 2 * n;
      p1 = precoffs;
      for (i = 0; i < n; i++) {
         x = a->xR[p0 + 0];
         y = a->xR[p0 + 1];
         bx = plan->precr.xR[p1 + 0];
         by = -plan->precr.xR[p1 + 1];
         bufa->xR[2 * i + 0] = x * bx - y * by;
         bufa->xR[2 * i + 1] = x * by + y * bx;
         p0 = p0 + 2;
         p1 = p1 + 2;
      }
      for (i = 2 * n; i < 2 * m; i++) {
         bufa->xR[i] = 0.0;
      }
   // Perform convolution of A and Z (using precomputed
   // FFT of Z stored in Plan structure).
      ftbase_ftapplysubplan(plan, subplan, bufa, 0, 0, bufc, 1, _state);
      p0 = 0;
      p1 = precoffs + 2 * m;
      for (i = 0; i < m; i++) {
         ax = bufa->xR[p0 + 0];
         ay = bufa->xR[p0 + 1];
         bx = plan->precr.xR[p1 + 0];
         by = plan->precr.xR[p1 + 1];
         bufa->xR[p0 + 0] = ax * bx - ay * by;
         bufa->xR[p0 + 1] = -(ax * by + ay * bx);
         p0 = p0 + 2;
         p1 = p1 + 2;
      }
      ftbase_ftapplysubplan(plan, subplan, bufa, 0, 0, bufc, 1, _state);
   // Post processing:
   //     A:=conj(Z)*conj(A)/M
   // Here conj(A)/M corresponds to last stage of inverse DFT,
   // and conj(Z) comes from Bluestein's FFT algorithm.
      p0 = precoffs;
      p1 = 0;
      p2 = abase + aoffset + op * 2 * n;
      for (i = 0; i < n; i++) {
         bx = plan->precr.xR[p0 + 0];
         by = plan->precr.xR[p0 + 1];
         rx = bufa->xR[p1 + 0] / m;
         ry = -bufa->xR[p1 + 1] / m;
         a->xR[p2 + 0] = rx * bx - ry * (-by);
         a->xR[p2 + 1] = rx * (-by) + ry * bx;
         p0 = p0 + 2;
         p1 = p1 + 2;
         p2 = p2 + 2;
      }
   }
}

// This subroutine applies complex Rader's FFT to input/output array A.
//
// Inputs:
//     A           -   array, must be large enough for plan to work
//     ABase       -   base offset in array A, this value points to start of
//                     subarray whose length is equal to length of the plan
//     AOffset     -   offset with respect to ABase, 0 <= AOffset<PlanLength.
//                     This is an offset within large PlanLength-subarray of
//                     the chunk to process.
//     OperandsCnt -   number of repeated operands (length N each)
//     N           -   original data length (measured in complex numbers)
//     SubPlan     -   position of the (N-1)-point FFT subplan which is used
//                     by transformation
//     RQ          -   primitive root modulo N
//     RIQ         -   inverse of primitive root modulo N
//     PrecOffs    -   offset of the precomputed data for the plan
//     Buf         -   temporary array
//
// Outputs:
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftradersfft(fasttransformplan *plan, RVector *a, ae_int_t abase, ae_int_t aoffset, ae_int_t operandscnt, ae_int_t n, ae_int_t subplan, ae_int_t rq, ae_int_t riq, ae_int_t precoffs, RVector *buf, ae_state *_state) {
   ae_int_t opidx;
   ae_int_t i;
   ae_int_t q;
   ae_int_t kq;
   ae_int_t kiq;
   double x0;
   double y0;
   ae_int_t p0;
   ae_int_t p1;
   double ax;
   double ay;
   double bx;
   double by;
   double rx;
   double ry;
   ae_assert(operandscnt >= 1, "FTApplyComplexRefFFT: OperandsCnt<1", _state);
// Process operands
   for (opidx = 0; opidx < operandscnt; opidx++) {
   // fill QA
      kq = 1;
      p0 = abase + aoffset + opidx * n * 2;
      p1 = aoffset + opidx * n * 2;
      rx = a->xR[p0 + 0];
      ry = a->xR[p0 + 1];
      x0 = rx;
      y0 = ry;
      for (q = 0; q < n - 1; q++) {
         ax = a->xR[p0 + 2 * kq + 0];
         ay = a->xR[p0 + 2 * kq + 1];
         buf->xR[p1 + 0] = ax;
         buf->xR[p1 + 1] = ay;
         rx = rx + ax;
         ry = ry + ay;
         kq = kq * rq % n;
         p1 = p1 + 2;
      }
      p0 = abase + aoffset + opidx * n * 2;
      p1 = aoffset + opidx * n * 2;
      for (q = 0; q < n - 1; q++) {
         a->xR[p0] = buf->xR[p1];
         a->xR[p0 + 1] = buf->xR[p1 + 1];
         p0 = p0 + 2;
         p1 = p1 + 2;
      }
   // Convolution
      ftbase_ftapplysubplan(plan, subplan, a, abase, aoffset + opidx * n * 2, buf, 1, _state);
      p0 = abase + aoffset + opidx * n * 2;
      p1 = precoffs;
      for (i = 0; i < n - 1; i++) {
         ax = a->xR[p0 + 0];
         ay = a->xR[p0 + 1];
         bx = plan->precr.xR[p1 + 0];
         by = plan->precr.xR[p1 + 1];
         a->xR[p0 + 0] = ax * bx - ay * by;
         a->xR[p0 + 1] = -(ax * by + ay * bx);
         p0 = p0 + 2;
         p1 = p1 + 2;
      }
      ftbase_ftapplysubplan(plan, subplan, a, abase, aoffset + opidx * n * 2, buf, 1, _state);
      p0 = abase + aoffset + opidx * n * 2;
      for (i = 0; i < n - 1; i++) {
         a->xR[p0 + 0] = a->xR[p0 + 0] / (n - 1);
         a->xR[p0 + 1] = -a->xR[p0 + 1] / (n - 1);
         p0 = p0 + 2;
      }
   // Result
      buf->xR[aoffset + opidx * n * 2 + 0] = rx;
      buf->xR[aoffset + opidx * n * 2 + 1] = ry;
      kiq = 1;
      p0 = aoffset + opidx * n * 2;
      p1 = abase + aoffset + opidx * n * 2;
      for (q = 0; q < n - 1; q++) {
         buf->xR[p0 + 2 * kiq + 0] = x0 + a->xR[p1 + 0];
         buf->xR[p0 + 2 * kiq + 1] = y0 + a->xR[p1 + 1];
         kiq = kiq * riq % n;
         p1 = p1 + 2;
      }
      p0 = abase + aoffset + opidx * n * 2;
      p1 = aoffset + opidx * n * 2;
      for (q = 0; q < n; q++) {
         a->xR[p0] = buf->xR[p1];
         a->xR[p0 + 1] = buf->xR[p1 + 1];
         p0 = p0 + 2;
         p1 = p1 + 2;
      }
   }
}

// This subroutine applies subplan to input/output array A.
//
// Inputs:
//     Plan        -   transformation plan
//     SubPlan     -   subplan index
//     A           -   array, must be large enough for plan to work
//     ABase       -   base offset in array A, this value points to start of
//                     subarray whose length is equal to length of the plan
//     AOffset     -   offset with respect to ABase, 0 <= AOffset<PlanLength.
//                     This is an offset within large PlanLength-subarray of
//                     the chunk to process.
//     Buf         -   temporary buffer whose length is equal to plan length
//                     (without taking into account RepCnt) or larger.
//     OffsBuf     -   offset in the buffer array
//     RepCnt      -   repetition count (transformation is repeatedly applied
//                     to subsequent subarrays)
//
// Outputs:
//     Plan        -   plan (temporary buffers can be modified, plan itself
//                     is unchanged and can be reused)
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftapplysubplan(fasttransformplan *plan, ae_int_t subplan, RVector *a, ae_int_t abase, ae_int_t aoffset, RVector *buf, ae_int_t repcnt, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t rowidx;
   ae_int_t i;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t operation;
   ae_int_t operandscnt;
   ae_int_t operandsize;
   ae_int_t microvectorsize;
   ae_int_t param0;
   ae_int_t param1;
   ae_int_t parentsize;
   ae_int_t childsize;
   ae_int_t chunksize;
   ae_int_t lastchunksize;
   ae_frame_make(_state, &_frame_block);
   RefObj(srealarray, bufa, _state);
   RefObj(srealarray, bufb, _state);
   RefObj(srealarray, bufc, _state);
   RefObj(srealarray, bufd, _state);
   ae_assert(plan->entries.xyZ[subplan][ftbase_coltype] == ftbase_opstart, "FTApplySubPlan: incorrect subplan header", _state);
   rowidx = subplan + 1;
   while (plan->entries.xyZ[rowidx][ftbase_coltype] != ftbase_opend) {
      operation = plan->entries.xyZ[rowidx][ftbase_coltype];
      operandscnt = repcnt * plan->entries.xyZ[rowidx][ftbase_coloperandscnt];
      operandsize = plan->entries.xyZ[rowidx][ftbase_coloperandsize];
      microvectorsize = plan->entries.xyZ[rowidx][ftbase_colmicrovectorsize];
      param0 = plan->entries.xyZ[rowidx][ftbase_colparam0];
      param1 = plan->entries.xyZ[rowidx][ftbase_colparam1];
      touchint(&param1, _state);
   // Process "jump" operation
      if (operation == ftbase_opjmp) {
         rowidx = rowidx + plan->entries.xyZ[rowidx][ftbase_colparam0];
         continue;
      }
   // Process "parallel call" operation:
   // * we perform initial check for consistency between parent and child plans
   // * we call FTSplitAndApplyParallelPlan(), which splits parallel plan into
   //   several parallel tasks
      if (operation == ftbase_opparallelcall) {
         parentsize = operandsize * microvectorsize;
         childsize = plan->entries.xyZ[rowidx + param0][ftbase_coloperandscnt] * plan->entries.xyZ[rowidx + param0][ftbase_coloperandsize] * plan->entries.xyZ[rowidx + param0][ftbase_colmicrovectorsize];
         ae_assert(plan->entries.xyZ[rowidx + param0][ftbase_coltype] == ftbase_opstart, "FTApplySubPlan: incorrect child subplan header", _state);
         ae_assert(parentsize == childsize, "FTApplySubPlan: incorrect child subplan header", _state);
         chunksize = maxint(ftbase_recursivethreshold / childsize, 1, _state);
         lastchunksize = operandscnt % chunksize;
         if (lastchunksize == 0) {
            lastchunksize = chunksize;
         }
         i = 0;
         while (i < operandscnt) {
            chunksize = minint(chunksize, operandscnt - i, _state);
            ftbase_ftapplysubplan(plan, rowidx + param0, a, abase, aoffset + i * childsize, buf, chunksize, _state);
            i = i + chunksize;
         }
         rowidx = rowidx + 1;
         continue;
      }
   // Process "reference complex FFT" operation
      if (operation == ftbase_opcomplexreffft) {
         ftbase_ftapplycomplexreffft(a, abase + aoffset, operandscnt, operandsize, microvectorsize, buf, _state);
         rowidx = rowidx + 1;
         continue;
      }
   // Process "codelet FFT" operation
      if (operation == ftbase_opcomplexcodeletfft) {
         ftbase_ftapplycomplexcodeletfft(a, abase + aoffset, operandscnt, operandsize, microvectorsize, _state);
         rowidx = rowidx + 1;
         continue;
      }
   // Process "integrated codelet FFT" operation
      if (operation == ftbase_opcomplexcodelettwfft) {
         ftbase_ftapplycomplexcodelettwfft(a, abase + aoffset, operandscnt, operandsize, microvectorsize, _state);
         rowidx = rowidx + 1;
         continue;
      }
   // Process Bluestein's FFT operation
      if (operation == ftbase_opbluesteinsfft) {
         ae_assert(microvectorsize == 2, "FTApplySubPlan: microvectorsize != 2 for Bluesteins FFT", _state);
         ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufa, _state);
         ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufb, _state);
         ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufc, _state);
         ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufd, _state);
         ftbase_ftbluesteinsfft(plan, a, abase, aoffset, operandscnt, operandsize, plan->entries.xyZ[rowidx][ftbase_colparam0], plan->entries.xyZ[rowidx][ftbase_colparam2], rowidx + plan->entries.xyZ[rowidx][ftbase_colparam1], &bufa->val, &bufb->val, &bufc->val, &bufd->val, _state);
         ae_shared_pool_recycle(&plan->bluesteinpool, &_bufa, _state);
         ae_shared_pool_recycle(&plan->bluesteinpool, &_bufb, _state);
         ae_shared_pool_recycle(&plan->bluesteinpool, &_bufc, _state);
         ae_shared_pool_recycle(&plan->bluesteinpool, &_bufd, _state);
         rowidx = rowidx + 1;
         continue;
      }
   // Process Rader's FFT
      if (operation == ftbase_opradersfft) {
         ftbase_ftradersfft(plan, a, abase, aoffset, operandscnt, operandsize, rowidx + plan->entries.xyZ[rowidx][ftbase_colparam0], plan->entries.xyZ[rowidx][ftbase_colparam1], plan->entries.xyZ[rowidx][ftbase_colparam2], plan->entries.xyZ[rowidx][ftbase_colparam3], buf, _state);
         rowidx = rowidx + 1;
         continue;
      }
   // Process "complex twiddle factors" operation
      if (operation == ftbase_opcomplexfftfactors) {
         ae_assert(microvectorsize == 2, "FTApplySubPlan: MicrovectorSize != 1", _state);
         n1 = plan->entries.xyZ[rowidx][ftbase_colparam0];
         n2 = operandsize / n1;
         for (i = 0; i < operandscnt; i++) {
            ftbase_ffttwcalc(a, abase + aoffset + i * operandsize * 2, n1, n2, _state);
         }
         rowidx = rowidx + 1;
         continue;
      }
   // Process "complex transposition" operation
      if (operation == ftbase_opcomplextranspose) {
         ae_assert(microvectorsize == 2, "FTApplySubPlan: MicrovectorSize != 1", _state);
         n1 = plan->entries.xyZ[rowidx][ftbase_colparam0];
         n2 = operandsize / n1;
         for (i = 0; i < operandscnt; i++) {
            ftbase_internalcomplexlintranspose(a, n1, n2, abase + aoffset + i * operandsize * 2, buf, _state);
         }
         rowidx = rowidx + 1;
         continue;
      }
   // Error
      ae_assert(false, "FTApplySubPlan: unexpected plan type", _state);
   }
   ae_frame_leave(_state);
}

// This subroutine applies transformation plan to input/output array A.
//
// Inputs:
//     Plan        -   transformation plan
//     A           -   array, must be large enough for plan to work
//     OffsA       -   offset of the subarray to process
//     RepCnt      -   repetition count (transformation is repeatedly applied
//                     to subsequent subarrays)
//
// Outputs:
//     Plan        -   plan (temporary buffers can be modified, plan itself
//                     is unchanged and can be reused)
//     A           -   transformed array
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
void ftapplyplan(fasttransformplan *plan, RVector *a, ae_int_t offsa, ae_int_t repcnt, ae_state *_state) {
   ae_int_t plansize;
   ae_int_t i;
   plansize = plan->entries.xyZ[0][ftbase_coloperandscnt] * plan->entries.xyZ[0][ftbase_coloperandsize] * plan->entries.xyZ[0][ftbase_colmicrovectorsize];
   for (i = 0; i < repcnt; i++) {
      ftbase_ftapplysubplan(plan, 0, a, offsa + plansize * i, 0, &plan->buffer, 1, _state);
   }
}

// Same as FTPushEntry(), but sets Param0 AND Param1.
// This function pushes one more entry to the plan. It resized Entries matrix
// if needed.
//
// Inputs:
//     Plan        -   plan (generated so far)
//     RowPtr      -   index which points to past-the-last entry generated so far
//     EType       -   entry type
//     EOpCnt      -   operands count
//     EOpSize     -   operand size
//     EMcvSize    -   microvector size
//     EParam0     -   parameter 0
//     EParam1     -   parameter 1
//
// Outputs:
//     Plan        -   updated plan
//     RowPtr      -   updated pointer
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftpushentry2(fasttransformplan *plan, ae_int_t *rowptr, ae_int_t etype, ae_int_t eopcnt, ae_int_t eopsize, ae_int_t emcvsize, ae_int_t eparam0, ae_int_t eparam1, ae_state *_state) {
   if (*rowptr >= plan->entries.rows) {
      imatrixresize(&plan->entries, maxint(2 * plan->entries.rows, 1, _state), ftbase_colscnt, _state);
   }
   plan->entries.xyZ[*rowptr][ftbase_coltype] = etype;
   plan->entries.xyZ[*rowptr][ftbase_coloperandscnt] = eopcnt;
   plan->entries.xyZ[*rowptr][ftbase_coloperandsize] = eopsize;
   plan->entries.xyZ[*rowptr][ftbase_colmicrovectorsize] = emcvsize;
   plan->entries.xyZ[*rowptr][ftbase_colparam0] = eparam0;
   plan->entries.xyZ[*rowptr][ftbase_colparam1] = eparam1;
   plan->entries.xyZ[*rowptr][ftbase_colparam2] = 0;
   plan->entries.xyZ[*rowptr][ftbase_colparam3] = 0;
   *rowptr = *rowptr + 1;
}

// Same as FTPushEntry(), but sets Param0, Param1, Param2 and Param3.
// This function pushes one more entry to the plan. It resized Entries matrix
// if needed.
//
// Inputs:
//     Plan        -   plan (generated so far)
//     RowPtr      -   index which points to past-the-last entry generated so far
//     EType       -   entry type
//     EOpCnt      -   operands count
//     EOpSize     -   operand size
//     EMcvSize    -   microvector size
//     EParam0     -   parameter 0
//     EParam1     -   parameter 1
//     EParam2     -   parameter 2
//     EParam3     -   parameter 3
//
// Outputs:
//     Plan        -   updated plan
//     RowPtr      -   updated pointer
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftpushentry4(fasttransformplan *plan, ae_int_t *rowptr, ae_int_t etype, ae_int_t eopcnt, ae_int_t eopsize, ae_int_t emcvsize, ae_int_t eparam0, ae_int_t eparam1, ae_int_t eparam2, ae_int_t eparam3, ae_state *_state) {
   if (*rowptr >= plan->entries.rows) {
      imatrixresize(&plan->entries, maxint(2 * plan->entries.rows, 1, _state), ftbase_colscnt, _state);
   }
   plan->entries.xyZ[*rowptr][ftbase_coltype] = etype;
   plan->entries.xyZ[*rowptr][ftbase_coloperandscnt] = eopcnt;
   plan->entries.xyZ[*rowptr][ftbase_coloperandsize] = eopsize;
   plan->entries.xyZ[*rowptr][ftbase_colmicrovectorsize] = emcvsize;
   plan->entries.xyZ[*rowptr][ftbase_colparam0] = eparam0;
   plan->entries.xyZ[*rowptr][ftbase_colparam1] = eparam1;
   plan->entries.xyZ[*rowptr][ftbase_colparam2] = eparam2;
   plan->entries.xyZ[*rowptr][ftbase_colparam3] = eparam3;
   *rowptr = *rowptr + 1;
}

// This function pushes one more entry to the plan. It resizes Entries matrix
// if needed.
//
// Inputs:
//     Plan        -   plan (generated so far)
//     RowPtr      -   index which points to past-the-last entry generated so far
//     EType       -   entry type
//     EOpCnt      -   operands count
//     EOpSize     -   operand size
//     EMcvSize    -   microvector size
//     EParam0     -   parameter 0
//
// Outputs:
//     Plan        -   updated plan
//     RowPtr      -   updated pointer
//
// NOTE: Param1 is set to -1.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftpushentry(fasttransformplan *plan, ae_int_t *rowptr, ae_int_t etype, ae_int_t eopcnt, ae_int_t eopsize, ae_int_t emcvsize, ae_int_t eparam0, ae_state *_state) {
   ftbase_ftpushentry2(plan, rowptr, etype, eopcnt, eopsize, emcvsize, eparam0, -1, _state);
}

#if 0
// Forward reference to an indirect recursive call. //(@) Already declared externally.
void ftcomplexfftplan(ae_int_t n, ae_int_t k, fasttransformplan *plan, ae_state *_state);
#endif

// This subroutine precomputes data for complex Bluestein's  FFT  and  writes
// them to array PrecR[] at specified offset. It  is  responsibility  of  the
// caller to make sure that PrecR[] is large enough.
//
// Inputs:
//     N           -   original size of the transform
//     M           -   size of the "padded" Bluestein's transform
//     PrecR       -   preallocated array
//     Offs        -   offset
//
// Outputs:
//     PrecR       -   data at Offs:Offs+4*M-1 are modified:
//                     * PrecR[Offs:Offs+2*M-1] stores Z[k]=exp(i*pi*k^2/N)
//                     * PrecR[Offs+2*M:Offs+4*M-1] stores FFT of the Z
//                     Other parts of PrecR are unchanged.
//
// NOTE: this function performs internal M-point FFT. It allocates temporary
//       plan which is destroyed after leaving this function.
// ALGLIB: Copyright 08.05.2013 by Sergey Bochkanov
static void ftbase_ftprecomputebluesteinsfft(ae_int_t n, ae_int_t m, RVector *precr, ae_int_t offs, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t i;
   double bx;
   double by;
   ae_frame_make(_state, &_frame_block);
   NewObj(fasttransformplan, plan, _state);
// Fill first half of PrecR with b[k] = exp(i*pi*k^2/N)
   for (i = 0; i < 2 * m; i++) {
      precr->xR[offs + i] = 0.0;
   }
   for (i = 0; i < n; i++) {
      bx = ae_cos(pi / n * i * i, _state);
      by = ae_sin(pi / n * i * i, _state);
      precr->xR[offs + 2 * i + 0] = bx;
      precr->xR[offs + 2 * i + 1] = by;
      precr->xR[offs + 2 * ((m - i) % m) + 0] = bx;
      precr->xR[offs + 2 * ((m - i) % m) + 1] = by;
   }
// Precomputed FFT
   ftcomplexfftplan(m, 1, &plan, _state);
   for (i = 0; i < 2 * m; i++) {
      precr->xR[offs + 2 * m + i] = precr->xR[offs + i];
   }
   ftbase_ftapplysubplan(&plan, 0, precr, offs + 2 * m, 0, &plan.buffer, 1, _state);
   ae_frame_leave(_state);
}

// This subroutine precomputes data for complex Rader's FFT and  writes  them
// to array PrecR[] at specified offset. It  is  responsibility of the caller
// to make sure that PrecR[] is large enough.
//
// Inputs:
//     N           -   original size of the transform (before reduction to N-1)
//     RQ          -   primitive root modulo N
//     RIQ         -   inverse of primitive root modulo N
//     PrecR       -   preallocated array
//     Offs        -   offset
//
// Outputs:
//     PrecR       -   data at Offs:Offs+2*(N-1)-1 store FFT of Rader's factors,
//                     other parts of PrecR are unchanged.
//
// NOTE: this function performs internal (N-1)-point FFT. It allocates temporary
//       plan which is destroyed after leaving this function.
// ALGLIB: Copyright 08.05.2013 by Sergey Bochkanov
static void ftbase_ftprecomputeradersfft(ae_int_t n, ae_int_t rq, ae_int_t riq, RVector *precr, ae_int_t offs, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t q;
   ae_int_t kiq;
   double v;
   ae_frame_make(_state, &_frame_block);
   NewObj(fasttransformplan, plan, _state);
// Fill PrecR with Rader factors, perform FFT
   kiq = 1;
   for (q = 0; q < n - 1; q++) {
      v = -2 * pi * kiq / n;
      precr->xR[offs + 2 * q + 0] = ae_cos(v, _state);
      precr->xR[offs + 2 * q + 1] = ae_sin(v, _state);
      kiq = kiq * riq % n;
   }
   ftcomplexfftplan(n - 1, 1, &plan, _state);
   ftbase_ftapplysubplan(&plan, 0, precr, offs, 0, &plan.buffer, 1, _state);
   ae_frame_leave(_state);
}

// Recurrent function called by FTComplexFFTPlan() and other functions. It
// recursively builds transformation plan
//
// Inputs:
//     N           -   FFT length (in complex numbers), N >= 1
//     K           -   number of repetitions, K >= 1
//     ChildPlan   -   if True, plan generator inserts OpStart/opEnd in the
//                     plan header/footer.
//     TopmostPlan -   if True, plan generator assumes that it is topmost plan:
//                     * it may use global buffer for transpositions
//                     and there is no other plan which executes in parallel
//     RowPtr      -   index which points to past-the-last entry generated so far
//     BluesteinSize-  amount of storage (in real numbers) required for Bluestein buffer
//     PrecRPtr    -   pointer to unused part of precomputed real buffer (Plan.PrecR):
//                     * when this function stores some data to precomputed buffer,
//                       it advances pointer.
//                     * it is responsibility of the function to assert that
//                       Plan.PrecR has enough space to store data before actually
//                       writing to buffer.
//                     * it is responsibility of the caller to allocate enough
//                       space before calling this function
//     PrecIPtr    -   pointer to unused part of precomputed integer buffer (Plan.PrecI):
//                     * when this function stores some data to precomputed buffer,
//                       it advances pointer.
//                     * it is responsibility of the function to assert that
//                       Plan.PrecR has enough space to store data before actually
//                       writing to buffer.
//                     * it is responsibility of the caller to allocate enough
//                       space before calling this function
//     Plan        -   plan (generated so far)
//
// Outputs:
//     RowPtr      -   updated pointer (advanced by number of entries generated
//                     by function)
//     BluesteinSize-  updated amount
//                     (may be increased, but may never be decreased)
//
// NOTE: in case TopmostPlan is True, ChildPlan is also must be True.
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
static void ftbase_ftcomplexfftplanrec(ae_int_t n, ae_int_t k, bool childplan, bool topmostplan, ae_int_t *rowptr, ae_int_t *bluesteinsize, ae_int_t *precrptr, ae_int_t *preciptr, fasttransformplan *plan, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t m;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t gq;
   ae_int_t giq;
   ae_int_t row0;
   ae_int_t row1;
   ae_int_t row2;
   ae_int_t row3;
   ae_frame_make(_state, &_frame_block);
   NewObj(srealarray, localbuf, _state);
   ae_assert(n > 0, "FTComplexFFTPlan: N <= 0", _state);
   ae_assert(k > 0, "FTComplexFFTPlan: K <= 0", _state);
   ae_assert(!topmostplan || childplan, "FTComplexFFTPlan: ChildPlan is inconsistent with TopmostPlan", _state);
// Try to generate "topmost" plan
   if (topmostplan && n > ftbase_recursivethreshold) {
      ftbase_ftfactorize(n, false, &n1, &n2, _state);
      if (n1 * n2 == 0) {
      // Handle prime-factor FFT with Bluestein's FFT.
      // Determine size of Bluestein's buffer.
         m = ftbasefindsmooth(2 * n - 1, _state);
         *bluesteinsize = maxint(2 * m, *bluesteinsize, _state);
      // Generate plan
         ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
         ftbase_ftpushentry4(plan, rowptr, ftbase_opbluesteinsfft, k, n, 2, m, 2, *precrptr, 0, _state);
         row0 = *rowptr;
         ftbase_ftpushentry(plan, rowptr, ftbase_opjmp, 0, 0, 0, 0, _state);
         ftbase_ftcomplexfftplanrec(m, 1, true, true, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
         row1 = *rowptr;
         plan->entries.xyZ[row0][ftbase_colparam0] = row1 - row0;
         ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
      // Fill precomputed buffer
         ftbase_ftprecomputebluesteinsfft(n, m, &plan->precr, *precrptr, _state);
      // Update pointer to the precomputed area
         *precrptr = *precrptr + 4 * m;
      } else {
      // Handle composite FFT with recursive Cooley-Tukey which
      // uses global buffer instead of local one.
         ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
         ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
         row0 = *rowptr;
         ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k * n2, n1, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
         ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1, _state);
         ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n2, _state);
         row2 = *rowptr;
         ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k * n1, n2, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
         ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
         ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
         row1 = *rowptr;
         ftbase_ftcomplexfftplanrec(n1, 1, true, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
         plan->entries.xyZ[row0][ftbase_colparam0] = row1 - row0;
         row3 = *rowptr;
         ftbase_ftcomplexfftplanrec(n2, 1, true, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
         plan->entries.xyZ[row2][ftbase_colparam0] = row3 - row2;
      }
      ae_frame_leave(_state);
      return;
   }
// Prepare "non-topmost" plan:
// * calculate factorization
// * use local (shared) buffer
// * update buffer size - ANY plan will need at least
//   2*N temporaries, additional requirements can be
//   applied later
   ftbase_ftfactorize(n, false, &n1, &n2, _state);
// Handle FFT's with N1*N2=0: either small-N or prime-factor
   if (n1 * n2 == 0) {
      if (n <= ftbase_maxradix) {
      // Small-N FFT
         if (childplan) {
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
         }
         ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexcodeletfft, k, n, 2, 0, _state);
         if (childplan) {
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
         }
         ae_frame_leave(_state);
         return;
      }
      if (n <= ftbase_raderthreshold) {
      // Handle prime-factor FFT's with Rader's FFT
         m = n - 1;
         if (childplan) {
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
         }
         findprimitiverootandinverse(n, &gq, &giq, _state);
         ftbase_ftpushentry4(plan, rowptr, ftbase_opradersfft, k, n, 2, 2, gq, giq, *precrptr, _state);
         ftbase_ftprecomputeradersfft(n, gq, giq, &plan->precr, *precrptr, _state);
         *precrptr = *precrptr + 2 * (n - 1);
         row0 = *rowptr;
         ftbase_ftpushentry(plan, rowptr, ftbase_opjmp, 0, 0, 0, 0, _state);
         ftbase_ftcomplexfftplanrec(m, 1, true, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
         row1 = *rowptr;
         plan->entries.xyZ[row0][ftbase_colparam0] = row1 - row0;
         if (childplan) {
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
         }
      } else {
      // Handle prime-factor FFT's with Bluestein's FFT
         m = ftbasefindsmooth(2 * n - 1, _state);
         *bluesteinsize = maxint(2 * m, *bluesteinsize, _state);
         if (childplan) {
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
         }
         ftbase_ftpushentry4(plan, rowptr, ftbase_opbluesteinsfft, k, n, 2, m, 2, *precrptr, 0, _state);
         ftbase_ftprecomputebluesteinsfft(n, m, &plan->precr, *precrptr, _state);
         *precrptr = *precrptr + 4 * m;
         row0 = *rowptr;
         ftbase_ftpushentry(plan, rowptr, ftbase_opjmp, 0, 0, 0, 0, _state);
         ftbase_ftcomplexfftplanrec(m, 1, true, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
         row1 = *rowptr;
         plan->entries.xyZ[row0][ftbase_colparam0] = row1 - row0;
         if (childplan) {
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
         }
      }
      ae_frame_leave(_state);
      return;
   }
// Handle Cooley-Tukey FFT with small N1
   if (n1 <= ftbase_maxradix) {
   // Specialized transformation for small N1:
   // * N2 short inplace FFT's, each N1-point, with integrated twiddle factors
   // * N1 long FFT's
   // * final transposition
      if (childplan) {
         ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
      }
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexcodelettwfft, k, n1, 2 * n2, 0, _state);
      ftbase_ftcomplexfftplanrec(n2, k * n1, false, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
      if (childplan) {
         ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
      }
      ae_frame_leave(_state);
      return;
   }
// Handle general Cooley-Tukey FFT, either "flat" or "recursive"
   if (n <= ftbase_recursivethreshold) {
   // General code for large N1/N2, "flat" version without explicit recurrence
   // (nested subplans are inserted directly into the body of the plan)
      if (childplan) {
         ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
      }
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
      ftbase_ftcomplexfftplanrec(n1, k * n2, false, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1, _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n2, _state);
      ftbase_ftcomplexfftplanrec(n2, k * n1, false, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
      if (childplan) {
         ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
      }
   } else {
   // General code for large N1/N2, "recursive" version - nested subplans
   // are separated from the plan body.
   //
   // Generate parent plan.
      if (childplan) {
         ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
      }
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
      row0 = *rowptr;
      ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k * n2, n1, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1, _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n2, _state);
      row2 = *rowptr;
      ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k * n1, n2, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
      ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
      if (childplan) {
         ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
      }
   // Generate child subplans, insert refence to parent plans
      row1 = *rowptr;
      ftbase_ftcomplexfftplanrec(n1, 1, true, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
      plan->entries.xyZ[row0][ftbase_colparam0] = row1 - row0;
      row3 = *rowptr;
      ftbase_ftcomplexfftplanrec(n2, 1, true, false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
      plan->entries.xyZ[row2][ftbase_colparam0] = row3 - row2;
   }
   ae_frame_leave(_state);
}

// This subroutine generates FFT plan for K complex FFT's with length N each.
//
// Inputs:
//     N           -   FFT length (in complex numbers), N >= 1
//     K           -   number of repetitions, K >= 1
//
// Outputs:
//     Plan        -   plan
// ALGLIB: Copyright 05.04.2013 by Sergey Bochkanov
void ftcomplexfftplan(ae_int_t n, ae_int_t k, fasttransformplan *plan, ae_state *_state) {
   ae_frame _frame_block;
   ae_int_t rowptr;
   ae_int_t bluesteinsize;
   ae_int_t precrptr;
   ae_int_t preciptr;
   ae_int_t precrsize;
   ae_int_t precisize;
   ae_frame_make(_state, &_frame_block);
   SetObj(fasttransformplan, plan);
   NewObj(srealarray, bluesteinbuf, _state);
// Initial check for parameters
   ae_assert(n > 0, "FTComplexFFTPlan: N <= 0", _state);
   ae_assert(k > 0, "FTComplexFFTPlan: K <= 0", _state);
// Determine required sizes of precomputed real and integer
// buffers. This stage of code is highly dependent on internals
// of FTComplexFFTPlanRec() and must be kept synchronized with
// possible changes in internals of plan generation function.
//
// Buffer size is determined as follows:
// * N is factorized
// * we factor out anything which is less or equal to MaxRadix
// * prime factor F>RaderThreshold requires 4*FTBaseFindSmooth(2*F-1)
//   real entries to store precomputed Quantities for Bluestein's
//   transformation
// * prime factor F <= RaderThreshold does NOT require
//   precomputed storage
   precrsize = 0;
   precisize = 0;
   ftbase_ftdeterminespacerequirements(n, &precrsize, &precisize, _state);
   if (precrsize > 0) {
      ae_vector_set_length(&plan->precr, precrsize, _state);
   }
   if (precisize > 0) {
      ae_vector_set_length(&plan->preci, precisize, _state);
   }
// Generate plan
   rowptr = 0;
   precrptr = 0;
   preciptr = 0;
   bluesteinsize = 1;
   ae_vector_set_length(&plan->buffer, 2 * n * k, _state);
   ftbase_ftcomplexfftplanrec(n, k, true, true, &rowptr, &bluesteinsize, &precrptr, &preciptr, plan, _state);
   ae_vector_set_length(&bluesteinbuf.val, bluesteinsize, _state);
   ae_shared_pool_set_seed(&plan->bluesteinpool, &bluesteinbuf, sizeof(bluesteinbuf), srealarray_init, srealarray_copy, srealarray_free, _state);
// Check that actual amount of precomputed space used by transformation
// plan is EXACTLY equal to amount of space allocated by us.
   ae_assert(precrptr == precrsize, "FTComplexFFTPlan: internal error (PrecRPtr != PrecRSize)", _state);
   ae_assert(preciptr == precisize, "FTComplexFFTPlan: internal error (PrecRPtr != PrecRSize)", _state);
   ae_frame_leave(_state);
}

void fasttransformplan_init(void *_p, ae_state *_state, bool make_automatic) {
   fasttransformplan *p = (fasttransformplan *)_p;
   ae_touch_ptr((void *)p);
   ae_matrix_init(&p->entries, 0, 0, DT_INT, _state, make_automatic);
   ae_vector_init(&p->buffer, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->precr, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->preci, 0, DT_REAL, _state, make_automatic);
   ae_shared_pool_init(&p->bluesteinpool, _state, make_automatic);
}

void fasttransformplan_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   fasttransformplan *dst = (fasttransformplan *)_dst;
   fasttransformplan *src = (fasttransformplan *)_src;
   ae_matrix_copy(&dst->entries, &src->entries, _state, make_automatic);
   ae_vector_copy(&dst->buffer, &src->buffer, _state, make_automatic);
   ae_vector_copy(&dst->precr, &src->precr, _state, make_automatic);
   ae_vector_copy(&dst->preci, &src->preci, _state, make_automatic);
   ae_shared_pool_copy(&dst->bluesteinpool, &src->bluesteinpool, _state, make_automatic);
}

void fasttransformplan_free(void *_p, bool make_automatic) {
   fasttransformplan *p = (fasttransformplan *)_p;
   ae_touch_ptr((void *)p);
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
static bool hpccores_hpcpreparechunkedgradientx(RVector *weights, ae_int_t wcount, RVector *hpcbuf, ae_state *_state) {
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
static bool hpccores_hpcfinalizechunkedgradientx(RVector *buf, ae_int_t wcount, RVector *grad, ae_state *_state) {
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
//     Grad:=0
//     HPCPrepareChunkedGradient(Weights, WCount, NTotal, NOut, Buf)
//     foreach chunk-of-dataset do
//         HPCChunkedGradient(...)
//     HPCFinalizeChunkedGradient(Buf, Grad)
//
void hpcpreparechunkedgradient(RVector *weights, ae_int_t wcount, ae_int_t ntotal, ae_int_t nin, ae_int_t nout, mlpbuffers *buf, ae_state *_state) {
   ae_int_t i;
   ae_int_t batch4size;
   ae_int_t chunksize;
   chunksize = 4;
   batch4size = 3 * chunksize * ntotal + chunksize * (2 * nout + 1);
   if (buf->xy.rows < chunksize || buf->xy.cols < nin + nout) {
      ae_matrix_set_length(&buf->xy, chunksize, nin + nout, _state);
   }
   if (buf->xy2.rows < chunksize || buf->xy2.cols < nin + nout) {
      ae_matrix_set_length(&buf->xy2, chunksize, nin + nout, _state);
   }
   if (buf->xyrow.cnt < nin + nout) {
      ae_vector_set_length(&buf->xyrow, nin + nout, _state);
   }
   if (buf->x.cnt < nin) {
      ae_vector_set_length(&buf->x, nin, _state);
   }
   if (buf->y.cnt < nout) {
      ae_vector_set_length(&buf->y, nout, _state);
   }
   if (buf->desiredy.cnt < nout) {
      ae_vector_set_length(&buf->desiredy, nout, _state);
   }
   if (buf->batch4buf.cnt < batch4size) {
      ae_vector_set_length(&buf->batch4buf, batch4size, _state);
   }
   if (buf->hpcbuf.cnt < wcount) {
      ae_vector_set_length(&buf->hpcbuf, wcount, _state);
   }
   if (buf->g.cnt < wcount) {
      ae_vector_set_length(&buf->g, wcount, _state);
   }
   if (!hpccores_hpcpreparechunkedgradientx(weights, wcount, &buf->hpcbuf, _state)) {
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
//     Grad:=0
//     HPCPrepareChunkedGradient(Weights, WCount, NTotal, NOut, Buf)
//     foreach chunk-of-dataset do
//         HPCChunkedGradient(...)
//     HPCFinalizeChunkedGradient(Buf, Grad)
//
void hpcfinalizechunkedgradient(mlpbuffers *buf, RVector *grad, ae_state *_state) {
   ae_int_t i;
   if (!hpccores_hpcfinalizechunkedgradientx(&buf->hpcbuf, buf->wcount, grad, _state)) {
      for (i = 0; i < buf->wcount; i++) {
         grad->xR[i] = grad->xR[i] + buf->hpcbuf.xR[i];
      }
   }
}

// Fast kernel for chunked gradient.
//
bool hpcchunkedgradient(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf, double *e, bool naturalerrorfunc, ae_state *_state) {
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
bool hpcchunkedprocess(RVector *weights, ZVector *structinfo, RVector *columnmeans, RVector *columnsigmas, RMatrix *xy, ae_int_t cstart, ae_int_t csize, RVector *batch4buf, RVector *hpcbuf, ae_state *_state) {
#ifndef ALGLIB_INTERCEPTS_SSE2
   bool result;
   result = false;
   return result;
#else
   return _ialglib_i_hpcchunkedprocess(weights, structinfo, columnmeans, columnsigmas, xy, cstart, csize, batch4buf, hpcbuf);
#endif
}

void mlpbuffers_init(void *_p, ae_state *_state, bool make_automatic) {
   mlpbuffers *p = (mlpbuffers *)_p;
   ae_touch_ptr((void *)p);
   ae_vector_init(&p->batch4buf, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->hpcbuf, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->xy, 0, 0, DT_REAL, _state, make_automatic);
   ae_matrix_init(&p->xy2, 0, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->xyrow, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->desiredy, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic);
   ae_vector_init(&p->tmp0, 0, DT_REAL, _state, make_automatic);
}

void mlpbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic) {
   mlpbuffers *dst = (mlpbuffers *)_dst;
   mlpbuffers *src = (mlpbuffers *)_src;
   dst->chunksize = src->chunksize;
   dst->ntotal = src->ntotal;
   dst->nin = src->nin;
   dst->nout = src->nout;
   dst->wcount = src->wcount;
   ae_vector_copy(&dst->batch4buf, &src->batch4buf, _state, make_automatic);
   ae_vector_copy(&dst->hpcbuf, &src->hpcbuf, _state, make_automatic);
   ae_matrix_copy(&dst->xy, &src->xy, _state, make_automatic);
   ae_matrix_copy(&dst->xy2, &src->xy2, _state, make_automatic);
   ae_vector_copy(&dst->xyrow, &src->xyrow, _state, make_automatic);
   ae_vector_copy(&dst->x, &src->x, _state, make_automatic);
   ae_vector_copy(&dst->y, &src->y, _state, make_automatic);
   ae_vector_copy(&dst->desiredy, &src->desiredy, _state, make_automatic);
   dst->e = src->e;
   ae_vector_copy(&dst->g, &src->g, _state, make_automatic);
   ae_vector_copy(&dst->tmp0, &src->tmp0, _state, make_automatic);
}

void mlpbuffers_free(void *_p, bool make_automatic) {
   mlpbuffers *p = (mlpbuffers *)_p;
   ae_touch_ptr((void *)p);
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
