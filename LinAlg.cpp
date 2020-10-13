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
#include "LinAlg.h"

// === SPARSE Package ===
// Depends on: (AlgLibInternal) ABLASMKL, TSORT
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
static const double sparse_desiredloadfactor = 0.66;
static const double sparse_maxloadfactor = 0.75;
static const double sparse_growfactor = 2.00;
static const ae_int_t sparse_additional = 10;
static const ae_int_t sparse_linalgswitch = 16;

// This function creates sparse matrix in a Hash-Table format.
//
// This function creates Hast-Table matrix, which can be  converted  to  CRS
// format after its initialization is over. Typical  usage  scenario  for  a
// sparse matrix is:
// 1. creation in a Hash-Table format
// 2. insertion of the matrix elements
// 3. conversion to the CRS representation
// 4. matrix is passed to some linear algebra algorithm
//
// Some  information  about  different matrix formats can be found below, in
// the "NOTES" section.
//
// Inputs:
//     M           -   number of rows in a matrix, M >= 1
//     N           -   number of columns in a matrix, N >= 1
//     K           -   K >= 0, expected number of non-zero elements in a matrix.
//                     K can be inexact approximation, can be less than actual
//                     number  of  elements  (table will grow when needed) or
//                     even zero).
//                     It is important to understand that although hash-table
//                     may grow automatically, it is better to  provide  good
//                     estimate of data size.
//
// Outputs:
//     S           -   sparse M*N matrix in Hash-Table representation.
//                     All elements of the matrix are zero.
//
// NOTE 1
//
// Hash-tables use memory inefficiently, and they have to keep  some  amount
// of the "spare memory" in order to have good performance. Hash  table  for
// matrix with K non-zero elements will  need  C*K*(8+2*sizeof(int))  bytes,
// where C is a small constant, about 1.5-2 in magnitude.
//
// CRS storage, from the other side, is  more  memory-efficient,  and  needs
// just K*(8+sizeof(int))+M*sizeof(int) bytes, where M is a number  of  rows
// in a matrix.
//
// When you convert from the Hash-Table to CRS  representation, all unneeded
// memory will be freed.
//
// NOTE 2
//
// Comments of SparseMatrix structure outline  information  about  different
// sparse storage formats. We recommend you to read them before starting  to
// use ALGLIB sparse matrices.
//
// NOTE 3
//
// This function completely  overwrites S with new sparse matrix. Previously
// allocated storage is NOT reused. If you  want  to reuse already allocated
// memory, call SparseCreateBuf function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsecreate(const ae_int_t m, const ae_int_t n, const ae_int_t k, sparsematrix &s);
// API: void sparsecreate(const ae_int_t m, const ae_int_t n, sparsematrix &s);
void sparsecreate(ae_int_t m, ae_int_t n, ae_int_t k, sparsematrix *s) {
   SetObj(sparsematrix, s);
   sparsecreatebuf(m, n, k, s);
}

// This version of SparseCreate function creates sparse matrix in Hash-Table
// format, reusing previously allocated storage as much  as  possible.  Read
// comments for SparseCreate() for more information.
//
// Inputs:
//     M           -   number of rows in a matrix, M >= 1
//     N           -   number of columns in a matrix, N >= 1
//     K           -   K >= 0, expected number of non-zero elements in a matrix.
//                     K can be inexact approximation, can be less than actual
//                     number  of  elements  (table will grow when needed) or
//                     even zero).
//                     It is important to understand that although hash-table
//                     may grow automatically, it is better to  provide  good
//                     estimate of data size.
//     S           -   SparseMatrix structure which MAY contain some  already
//                     allocated storage.
//
// Outputs:
//     S           -   sparse M*N matrix in Hash-Table representation.
//                     All elements of the matrix are zero.
//                     Previously allocated storage is reused, if  its  size
//                     is compatible with expected number of non-zeros K.
//
// ALGLIB Project: Copyright 14.01.2014 by Sergey Bochkanov
// API: void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const ae_int_t k, const sparsematrix &s);
// API: void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const sparsematrix &s);
void sparsecreatebuf(ae_int_t m, ae_int_t n, ae_int_t k, sparsematrix *s) {
   ae_int_t i;
   ae_assert(m > 0, "SparseCreateBuf: M <= 0");
   ae_assert(n > 0, "SparseCreateBuf: N <= 0");
   ae_assert(k >= 0, "SparseCreateBuf: K<0");
// Hash-table size is max(existing_size,requested_size)
//
// NOTE: it is important to use ALL available memory for hash table
//       because it is impossible to efficiently reallocate table
//       without temporary storage. So, if we want table with up to
//       1.000.000 elements, we have to create such table from the
//       very beginning. Otherwise, the very idea of memory reuse
//       will be compromised.
   s->tablesize = RoundZ(k / sparse_desiredloadfactor + sparse_additional);
   vectorsetlengthatleast(&s->vals, s->tablesize);
   s->tablesize = s->vals.cnt;
// Initialize other fields
   s->matrixtype = 0;
   s->m = m;
   s->n = n;
   s->nfree = s->tablesize;
   vectorsetlengthatleast(&s->idx, 2 * s->tablesize);
   for (i = 0; i < s->tablesize; i++) {
      s->idx.xZ[2 * i] = -1;
   }
}

// This function creates sparse matrix in a CRS format (expert function for
// situations when you are running out of memory).
//
// This function creates CRS matrix. Typical usage scenario for a CRS matrix
// is:
// 1. creation (you have to tell number of non-zero elements at each row  at
//    this moment)
// 2. insertion of the matrix elements (row by row, from left to right)
// 3. matrix is passed to some linear algebra algorithm
//
// This function is a memory-efficient alternative to SparseCreate(), but it
// is more complex because it requires you to know in advance how large your
// matrix is. Some  information about  different matrix formats can be found
// in comments on SparseMatrix structure.  We recommend  you  to  read  them
// before starting to use ALGLIB sparse matrices..
//
// Inputs:
//     M           -   number of rows in a matrix, M >= 1
//     N           -   number of columns in a matrix, N >= 1
//     NER         -   number of elements at each row, array[M], NER[I] >= 0
//
// Outputs:
//     S           -   sparse M*N matrix in CRS representation.
//                     You have to fill ALL non-zero elements by calling
//                     SparseSet() BEFORE you try to use this matrix.
//
// NOTE: this function completely  overwrites  S  with  new  sparse  matrix.
//       Previously allocated storage is NOT reused. If you  want  to  reuse
//       already allocated memory, call SparseCreateCRSBuf function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsecreatecrs(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, sparsematrix &s);
void sparsecreatecrs(ae_int_t m, ae_int_t n, ZVector *ner, sparsematrix *s) {
   ae_int_t i;
   SetObj(sparsematrix, s);
   ae_assert(m > 0, "SparseCreateCRS: M <= 0");
   ae_assert(n > 0, "SparseCreateCRS: N <= 0");
   ae_assert(ner->cnt >= m, "SparseCreateCRS: Length(NER)<M");
   for (i = 0; i < m; i++) {
      ae_assert(ner->xZ[i] >= 0, "SparseCreateCRS: NER[] contains negative elements");
   }
   sparsecreatecrsbuf(m, n, ner, s);
}

// This function creates sparse matrix in a CRS format (expert function  for
// situations when you are running out  of  memory).  This  version  of  CRS
// matrix creation function may reuse memory already allocated in S.
//
// This function creates CRS matrix. Typical usage scenario for a CRS matrix
// is:
// 1. creation (you have to tell number of non-zero elements at each row  at
//    this moment)
// 2. insertion of the matrix elements (row by row, from left to right)
// 3. matrix is passed to some linear algebra algorithm
//
// This function is a memory-efficient alternative to SparseCreate(), but it
// is more complex because it requires you to know in advance how large your
// matrix is. Some  information about  different matrix formats can be found
// in comments on SparseMatrix structure.  We recommend  you  to  read  them
// before starting to use ALGLIB sparse matrices..
//
// Inputs:
//     M           -   number of rows in a matrix, M >= 1
//     N           -   number of columns in a matrix, N >= 1
//     NER         -   number of elements at each row, array[M], NER[I] >= 0
//     S           -   sparse matrix structure with possibly preallocated
//                     memory.
//
// Outputs:
//     S           -   sparse M*N matrix in CRS representation.
//                     You have to fill ALL non-zero elements by calling
//                     SparseSet() BEFORE you try to use this matrix.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsecreatecrsbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, const sparsematrix &s);
void sparsecreatecrsbuf(ae_int_t m, ae_int_t n, ZVector *ner, sparsematrix *s) {
   ae_int_t i;
   ae_int_t noe;
   ae_assert(m > 0, "SparseCreateCRSBuf: M <= 0");
   ae_assert(n > 0, "SparseCreateCRSBuf: N <= 0");
   ae_assert(ner->cnt >= m, "SparseCreateCRSBuf: Length(NER)<M");
   noe = 0;
   s->matrixtype = 1;
   s->ninitialized = 0;
   s->m = m;
   s->n = n;
   vectorsetlengthatleast(&s->ridx, s->m + 1);
   s->ridx.xZ[0] = 0;
   for (i = 0; i < s->m; i++) {
      ae_assert(ner->xZ[i] >= 0, "SparseCreateCRSBuf: NER[] contains negative elements");
      noe += ner->xZ[i];
      s->ridx.xZ[i + 1] = s->ridx.xZ[i] + ner->xZ[i];
   }
   vectorsetlengthatleast(&s->vals, noe);
   vectorsetlengthatleast(&s->idx, noe);
   if (noe == 0) {
      sparseinitduidx(s);
   }
}

// This function creates sparse matrix in  a  SKS  format  (skyline  storage
// format). In most cases you do not need this function - CRS format  better
// suits most use cases.
//
// Inputs:
//     M, N        -   number of rows(M) and columns (N) in a matrix:
//                     * M=N (as for now, ALGLIB supports only square SKS)
//                     * N >= 1
//                     * M >= 1
//     D           -   "bottom" bandwidths, array[M], D[I] >= 0.
//                     I-th element stores number of non-zeros at I-th  row,
//                     below the diagonal (diagonal itself is not  included)
//     U           -   "top" bandwidths, array[N], U[I] >= 0.
//                     I-th element stores number of non-zeros  at I-th row,
//                     above the diagonal (diagonal itself  is not included)
//
// Outputs:
//     S           -   sparse M*N matrix in SKS representation.
//                     All elements are filled by zeros.
//                     You may use sparseset() to change their values.
//
// NOTE: this function completely  overwrites  S  with  new  sparse  matrix.
//       Previously allocated storage is NOT reused. If you  want  to  reuse
//       already allocated memory, call SparseCreateSKSBuf function.
//
// ALGLIB Project: Copyright 13.01.2014 by Sergey Bochkanov
// API: void sparsecreatesks(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, sparsematrix &s);
void sparsecreatesks(ae_int_t m, ae_int_t n, ZVector *d, ZVector *u, sparsematrix *s) {
   ae_int_t i;
   SetObj(sparsematrix, s);
   ae_assert(m > 0, "SparseCreateSKS: M <= 0");
   ae_assert(n > 0, "SparseCreateSKS: N <= 0");
   ae_assert(m == n, "SparseCreateSKS: M != N");
   ae_assert(d->cnt >= m, "SparseCreateSKS: Length(D)<M");
   ae_assert(u->cnt >= n, "SparseCreateSKS: Length(U)<N");
   for (i = 0; i < m; i++) {
      ae_assert(d->xZ[i] >= 0, "SparseCreateSKS: D[] contains negative elements");
      ae_assert(d->xZ[i] <= i, "SparseCreateSKS: D[I]>I for some I");
   }
   for (i = 0; i < n; i++) {
      ae_assert(u->xZ[i] >= 0, "SparseCreateSKS: U[] contains negative elements");
      ae_assert(u->xZ[i] <= i, "SparseCreateSKS: U[I]>I for some I");
   }
   sparsecreatesksbuf(m, n, d, u, s);
}

// This is "buffered"  version  of  SparseCreateSKS()  which  reuses  memory
// previously allocated in S (of course, memory is reallocated if needed).
//
// This function creates sparse matrix in  a  SKS  format  (skyline  storage
// format). In most cases you do not need this function - CRS format  better
// suits most use cases.
//
// Inputs:
//     M, N        -   number of rows(M) and columns (N) in a matrix:
//                     * M=N (as for now, ALGLIB supports only square SKS)
//                     * N >= 1
//                     * M >= 1
//     D           -   "bottom" bandwidths, array[M], 0 <= D[I] <= I.
//                     I-th element stores number of non-zeros at I-th row,
//                     below the diagonal (diagonal itself is not included)
//     U           -   "top" bandwidths, array[N], 0 <= U[I] <= I.
//                     I-th element stores number of non-zeros at I-th row,
//                     above the diagonal (diagonal itself is not included)
//
// Outputs:
//     S           -   sparse M*N matrix in SKS representation.
//                     All elements are filled by zeros.
//                     You may use sparseset() to change their values.
//
// ALGLIB Project: Copyright 13.01.2014 by Sergey Bochkanov
// API: void sparsecreatesksbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, const sparsematrix &s);
void sparsecreatesksbuf(ae_int_t m, ae_int_t n, ZVector *d, ZVector *u, sparsematrix *s) {
   ae_int_t i;
   ae_int_t minmn;
   ae_int_t nz;
   ae_int_t mxd;
   ae_int_t mxu;
   ae_assert(m > 0, "SparseCreateSKSBuf: M <= 0");
   ae_assert(n > 0, "SparseCreateSKSBuf: N <= 0");
   ae_assert(m == n, "SparseCreateSKSBuf: M != N");
   ae_assert(d->cnt >= m, "SparseCreateSKSBuf: Length(D)<M");
   ae_assert(u->cnt >= n, "SparseCreateSKSBuf: Length(U)<N");
   for (i = 0; i < m; i++) {
      ae_assert(d->xZ[i] >= 0, "SparseCreateSKSBuf: D[] contains negative elements");
      ae_assert(d->xZ[i] <= i, "SparseCreateSKSBuf: D[I]>I for some I");
   }
   for (i = 0; i < n; i++) {
      ae_assert(u->xZ[i] >= 0, "SparseCreateSKSBuf: U[] contains negative elements");
      ae_assert(u->xZ[i] <= i, "SparseCreateSKSBuf: U[I]>I for some I");
   }
   minmn = imin2(m, n);
   s->matrixtype = 2;
   s->ninitialized = 0;
   s->m = m;
   s->n = n;
   vectorsetlengthatleast(&s->ridx, minmn + 1);
   s->ridx.xZ[0] = 0;
   nz = 0;
   for (i = 0; i < minmn; i++) {
      nz += 1 + d->xZ[i] + u->xZ[i];
      s->ridx.xZ[i + 1] = s->ridx.xZ[i] + 1 + d->xZ[i] + u->xZ[i];
   }
   vectorsetlengthatleast(&s->vals, nz);
   for (i = 0; i < nz; i++) {
      s->vals.xR[i] = 0.0;
   }
   vectorsetlengthatleast(&s->didx, m + 1);
   mxd = 0;
   for (i = 0; i < m; i++) {
      s->didx.xZ[i] = d->xZ[i];
      mxd = imax2(mxd, d->xZ[i]);
   }
   s->didx.xZ[m] = mxd;
   vectorsetlengthatleast(&s->uidx, n + 1);
   mxu = 0;
   for (i = 0; i < n; i++) {
      s->uidx.xZ[i] = u->xZ[i];
      mxu = imax2(mxu, u->xZ[i]);
   }
   s->uidx.xZ[n] = mxu;
}

// This function creates sparse matrix in  a  SKS  format  (skyline  storage
// format). Unlike more general  sparsecreatesks(),  this  function  creates
// sparse matrix with constant bandwidth.
//
// You may want to use this function instead of sparsecreatesks() when  your
// matrix has  constant  or  nearly-constant  bandwidth,  and  you  want  to
// simplify source code.
//
// Inputs:
//     M, N        -   number of rows(M) and columns (N) in a matrix:
//                     * M=N (as for now, ALGLIB supports only square SKS)
//                     * N >= 1
//                     * M >= 1
//     BW          -   matrix bandwidth, BW >= 0
//
// Outputs:
//     S           -   sparse M*N matrix in SKS representation.
//                     All elements are filled by zeros.
//                     You may use sparseset() to  change  their values.
//
// NOTE: this function completely  overwrites  S  with  new  sparse  matrix.
//       Previously allocated storage is NOT reused. If you  want  to  reuse
//       already allocated memory, call sparsecreatesksbandbuf function.
//
// ALGLIB Project: Copyright 25.12.2017 by Sergey Bochkanov
// API: void sparsecreatesksband(const ae_int_t m, const ae_int_t n, const ae_int_t bw, sparsematrix &s);
void sparsecreatesksband(ae_int_t m, ae_int_t n, ae_int_t bw, sparsematrix *s) {
   SetObj(sparsematrix, s);
   ae_assert(m > 0, "SparseCreateSKSBand: M <= 0");
   ae_assert(n > 0, "SparseCreateSKSBand: N <= 0");
   ae_assert(bw >= 0, "SparseCreateSKSBand: BW<0");
   ae_assert(m == n, "SparseCreateSKSBand: M != N");
   sparsecreatesksbandbuf(m, n, bw, s);
}

// This is "buffered" version  of  sparsecreatesksband() which reuses memory
// previously allocated in S (of course, memory is reallocated if needed).
//
// You may want to use this function instead  of  sparsecreatesksbuf()  when
// your matrix has  constant or nearly-constant  bandwidth,  and you want to
// simplify source code.
//
// Inputs:
//     M, N        -   number of rows(M) and columns (N) in a matrix:
//                     * M=N (as for now, ALGLIB supports only square SKS)
//                     * N >= 1
//                     * M >= 1
//     BW          -   bandwidth, BW >= 0
//
// Outputs:
//     S           -   sparse M*N matrix in SKS representation.
//                     All elements are filled by zeros.
//                     You may use sparseset() to change their values.
//
// ALGLIB Project: Copyright 13.01.2014 by Sergey Bochkanov
// API: void sparsecreatesksbandbuf(const ae_int_t m, const ae_int_t n, const ae_int_t bw, const sparsematrix &s);
void sparsecreatesksbandbuf(ae_int_t m, ae_int_t n, ae_int_t bw, sparsematrix *s) {
   ae_int_t i;
   ae_int_t minmn;
   ae_int_t nz;
   ae_int_t mxd;
   ae_int_t mxu;
   ae_int_t dui;
   ae_assert(m > 0, "SparseCreateSKSBandBuf: M <= 0");
   ae_assert(n > 0, "SparseCreateSKSBandBuf: N <= 0");
   ae_assert(m == n, "SparseCreateSKSBandBuf: M != N");
   ae_assert(bw >= 0, "SparseCreateSKSBandBuf: BW<0");
   minmn = imin2(m, n);
   s->matrixtype = 2;
   s->ninitialized = 0;
   s->m = m;
   s->n = n;
   vectorsetlengthatleast(&s->ridx, minmn + 1);
   s->ridx.xZ[0] = 0;
   nz = 0;
   for (i = 0; i < minmn; i++) {
      dui = imin2(i, bw);
      nz += 1 + 2 * dui;
      s->ridx.xZ[i + 1] = s->ridx.xZ[i] + 1 + 2 * dui;
   }
   vectorsetlengthatleast(&s->vals, nz);
   for (i = 0; i < nz; i++) {
      s->vals.xR[i] = 0.0;
   }
   vectorsetlengthatleast(&s->didx, m + 1);
   mxd = 0;
   for (i = 0; i < m; i++) {
      dui = imin2(i, bw);
      s->didx.xZ[i] = dui;
      mxd = imax2(mxd, dui);
   }
   s->didx.xZ[m] = mxd;
   vectorsetlengthatleast(&s->uidx, n + 1);
   mxu = 0;
   for (i = 0; i < n; i++) {
      dui = imin2(i, bw);
      s->uidx.xZ[i] = dui;
      mxu = imax2(mxu, dui);
   }
   s->uidx.xZ[n] = mxu;
}

// This function copies S0 to S1.
// This function completely deallocates memory owned by S1 before creating a
// copy of S0. If you want to reuse memory, use SparseCopyBuf.
//
// NOTE:  this  function  does  not verify its arguments, it just copies all
// fields of the structure.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsecopy(const sparsematrix &s0, sparsematrix &s1);
void sparsecopy(sparsematrix *s0, sparsematrix *s1) {
   SetObj(sparsematrix, s1);
   sparsecopybuf(s0, s1);
}

// This function copies S0 to S1.
// Memory already allocated in S1 is reused as much as possible.
//
// NOTE:  this  function  does  not verify its arguments, it just copies all
// fields of the structure.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsecopybuf(const sparsematrix &s0, const sparsematrix &s1);
void sparsecopybuf(sparsematrix *s0, sparsematrix *s1) {
   ae_int_t l;
   ae_int_t i;
   s1->matrixtype = s0->matrixtype;
   s1->m = s0->m;
   s1->n = s0->n;
   s1->nfree = s0->nfree;
   s1->ninitialized = s0->ninitialized;
   s1->tablesize = s0->tablesize;
// Initialization for arrays
   l = s0->vals.cnt;
   vectorsetlengthatleast(&s1->vals, l);
   for (i = 0; i < l; i++) {
      s1->vals.xR[i] = s0->vals.xR[i];
   }
   l = s0->ridx.cnt;
   vectorsetlengthatleast(&s1->ridx, l);
   for (i = 0; i < l; i++) {
      s1->ridx.xZ[i] = s0->ridx.xZ[i];
   }
   l = s0->idx.cnt;
   vectorsetlengthatleast(&s1->idx, l);
   for (i = 0; i < l; i++) {
      s1->idx.xZ[i] = s0->idx.xZ[i];
   }
// Initalization for CRS-parameters
   l = s0->uidx.cnt;
   vectorsetlengthatleast(&s1->uidx, l);
   for (i = 0; i < l; i++) {
      s1->uidx.xZ[i] = s0->uidx.xZ[i];
   }
   l = s0->didx.cnt;
   vectorsetlengthatleast(&s1->didx, l);
   for (i = 0; i < l; i++) {
      s1->didx.xZ[i] = s0->didx.xZ[i];
   }
}

// This function efficiently swaps contents of S0 and S1.
//
// ALGLIB Project: Copyright 16.01.2014 by Sergey Bochkanov
// API: void sparseswap(const sparsematrix &s0, const sparsematrix &s1);
void sparseswap(sparsematrix *s0, sparsematrix *s1) {
   swapi(&s1->matrixtype, &s0->matrixtype);
   swapi(&s1->m, &s0->m);
   swapi(&s1->n, &s0->n);
   swapi(&s1->nfree, &s0->nfree);
   swapi(&s1->ninitialized, &s0->ninitialized);
   swapi(&s1->tablesize, &s0->tablesize);
   ae_swap_vectors(&s1->vals, &s0->vals);
   ae_swap_vectors(&s1->ridx, &s0->ridx);
   ae_swap_vectors(&s1->idx, &s0->idx);
   ae_swap_vectors(&s1->uidx, &s0->uidx);
   ae_swap_vectors(&s1->didx, &s0->didx);
}

// This is hash function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
static ae_int_t sparse_hash(ae_int_t i, ae_int_t j, ae_int_t tabsize) {
   ae_frame _frame_block;
   ae_int_t result;
   ae_frame_make(&_frame_block);
   NewObj(hqrndstate, r);
   hqrndseed(i, j, &r);
   result = hqrnduniformi(&r, tabsize);
   ae_frame_leave();
   return result;
}

// This function adds value to S[i,j] - element of the sparse matrix. Matrix
// must be in a Hash-Table mode.
//
// In case S[i,j] already exists in the table, V i added to  its  value.  In
// case  S[i,j]  is  non-existent,  it  is  inserted  in  the  table.  Table
// automatically grows when necessary.
//
// Inputs:
//     S           -   sparse M*N matrix in Hash-Table representation.
//                     Exception will be thrown for CRS matrix.
//     I           -   row index of the element to modify, 0 <= I < M
//     J           -   column index of the element to modify, 0 <= J < N
//     V           -   value to add, must be finite number
//
// Outputs:
//     S           -   modified matrix
//
// NOTE 1:  when  S[i,j]  is exactly zero after modification, it is  deleted
// from the table.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparseadd(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v);
void sparseadd(sparsematrix *s, ae_int_t i, ae_int_t j, double v) {
   ae_int_t hashcode;
   ae_int_t tcode;
   ae_int_t k;
   ae_assert(s->matrixtype == 0, "SparseAdd: matrix must be in the Hash-Table mode to do this operation");
   ae_assert(i >= 0, "SparseAdd: I<0");
   ae_assert(i < s->m, "SparseAdd: I >= M");
   ae_assert(j >= 0, "SparseAdd: J<0");
   ae_assert(j < s->n, "SparseAdd: J >= N");
   ae_assert(isfinite(v), "SparseAdd: V is not finite number");
   if (v == 0.0) {
      return;
   }
   tcode = -1;
   k = s->tablesize;
   if ((1 - sparse_maxloadfactor) * k >= (double)(s->nfree)) {
      sparseresizematrix(s);
      k = s->tablesize;
   }
   hashcode = sparse_hash(i, j, k);
   while (true) {
      if (s->idx.xZ[2 * hashcode] == -1) {
         if (tcode != -1) {
            hashcode = tcode;
         }
         s->vals.xR[hashcode] = v;
         s->idx.xZ[2 * hashcode] = i;
         s->idx.xZ[2 * hashcode + 1] = j;
         if (tcode == -1) {
            s->nfree--;
         }
         return;
      } else {
         if (s->idx.xZ[2 * hashcode] == i && s->idx.xZ[2 * hashcode + 1] == j) {
            s->vals.xR[hashcode] += v;
            if (s->vals.xR[hashcode] == 0.0) {
               s->idx.xZ[2 * hashcode] = -2;
            }
            return;
         }
      // Is it deleted element?
         if (tcode == -1 && s->idx.xZ[2 * hashcode] == -2) {
            tcode = hashcode;
         }
      // Next step
         hashcode = (hashcode + 1) % k;
      }
   }
}

// This function modifies S[i,j] - element of the sparse matrix.
//
// For Hash-based storage format:
// * this function can be called at any moment - during matrix initialization
//   or later
// * new value can be zero or non-zero.  In case new value of S[i,j] is zero,
//   this element is deleted from the table.
// * this  function  has  no  effect when called with zero V for non-existent
//   element.
//
// For CRS-bases storage format:
// * this function can be called ONLY DURING MATRIX INITIALIZATION
// * zero values are stored in the matrix similarly to non-zero ones
// * elements must be initialized in correct order -  from top row to bottom,
//   within row - from left to right.
//
// For SKS storage:
// * this function can be called at any moment - during matrix initialization
//   or later
// * zero values are stored in the matrix similarly to non-zero ones
// * this function CAN NOT be called for non-existent (outside  of  the  band
//   specified during SKS matrix creation) elements. Say, if you created  SKS
//   matrix  with  bandwidth=2  and  tried to call sparseset(s,0,10,VAL),  an
//   exception will be generated.
//
// Inputs:
//     S           -   sparse M*N matrix in Hash-Table, SKS or CRS format.
//     I           -   row index of the element to modify, 0 <= I < M
//     J           -   column index of the element to modify, 0 <= J < N
//     V           -   value to set, must be finite number, can be zero
//
// Outputs:
//     S           -   modified matrix
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparseset(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v);
void sparseset(sparsematrix *s, ae_int_t i, ae_int_t j, double v) {
   ae_int_t hashcode;
   ae_int_t tcode;
   ae_int_t k;
   bool b;
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseSet: unsupported matrix storage format");
   ae_assert(i >= 0, "SparseSet: I<0");
   ae_assert(i < s->m, "SparseSet: I >= M");
   ae_assert(j >= 0, "SparseSet: J<0");
   ae_assert(j < s->n, "SparseSet: J >= N");
   ae_assert(isfinite(v), "SparseSet: V is not finite number");
// Hash-table matrix
   if (s->matrixtype == 0) {
      tcode = -1;
      k = s->tablesize;
      if ((1 - sparse_maxloadfactor) * k >= (double)(s->nfree)) {
         sparseresizematrix(s);
         k = s->tablesize;
      }
      hashcode = sparse_hash(i, j, k);
      while (true) {
         if (s->idx.xZ[2 * hashcode] == -1) {
            if (v != 0.0) {
               if (tcode != -1) {
                  hashcode = tcode;
               }
               s->vals.xR[hashcode] = v;
               s->idx.xZ[2 * hashcode] = i;
               s->idx.xZ[2 * hashcode + 1] = j;
               if (tcode == -1) {
                  s->nfree--;
               }
            }
            return;
         } else {
            if (s->idx.xZ[2 * hashcode] == i && s->idx.xZ[2 * hashcode + 1] == j) {
               if (v == 0.0) {
                  s->idx.xZ[2 * hashcode] = -2;
               } else {
                  s->vals.xR[hashcode] = v;
               }
               return;
            }
            if (tcode == -1 && s->idx.xZ[2 * hashcode] == -2) {
               tcode = hashcode;
            }
         // Next step
            hashcode = (hashcode + 1) % k;
         }
      }
   }
// CRS matrix
   if (s->matrixtype == 1) {
      ae_assert(s->ridx.xZ[i] <= s->ninitialized, "SparseSet: too few initialized elements at some row (you have promised more when called SparceCreateCRS)");
      ae_assert(s->ridx.xZ[i + 1] > s->ninitialized, "SparseSet: too many initialized elements at some row (you have promised less when called SparceCreateCRS)");
      ae_assert(s->ninitialized == s->ridx.xZ[i] || s->idx.xZ[s->ninitialized - 1] < j, "SparseSet: incorrect column order (you must fill every row from left to right)");
      s->vals.xR[s->ninitialized] = v;
      s->idx.xZ[s->ninitialized] = j;
      s->ninitialized++;
   // If matrix has been created then
   // initiale 'S.UIdx' and 'S.DIdx'
      if (s->ninitialized == s->ridx.xZ[s->m]) {
         sparseinitduidx(s);
      }
      return;
   }
// SKS matrix
   if (s->matrixtype == 2) {
      b = sparserewriteexisting(s, i, j, v);
      ae_assert(b, "SparseSet: an attempt to initialize out-of-band element of the SKS matrix");
      return;
   }
}

// This function returns S[i,j] - element of the sparse matrix.  Matrix  can
// be in any mode (Hash-Table, CRS, SKS), but this function is less efficient
// for CRS matrices. Hash-Table and SKS matrices can find  element  in  O(1)
// time, while  CRS  matrices need O(log(RS)) time, where RS is an number of
// non-zero elements in a row.
//
// Inputs:
//     S           -   sparse M*N matrix in Hash-Table representation.
//                     Exception will be thrown for CRS matrix.
//     I           -   row index of the element to modify, 0 <= I < M
//     J           -   column index of the element to modify, 0 <= J < N
//
// Result:
//     value of S[I,J] or zero (in case no element with such index is found)
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: double sparseget(const sparsematrix &s, const ae_int_t i, const ae_int_t j);
double sparseget(sparsematrix *s, ae_int_t i, ae_int_t j) {
   ae_int_t hashcode;
   ae_int_t k;
   ae_int_t k0;
   ae_int_t k1;
   double result;
   ae_assert(i >= 0, "SparseGet: I<0");
   ae_assert(i < s->m, "SparseGet: I >= M");
   ae_assert(j >= 0, "SparseGet: J<0");
   ae_assert(j < s->n, "SparseGet: J >= N");
   result = 0.0;
   if (s->matrixtype == 0) {
   // Hash-based storage
      result = 0.0;
      k = s->tablesize;
      hashcode = sparse_hash(i, j, k);
      while (true) {
         if (s->idx.xZ[2 * hashcode] == -1) {
            return result;
         }
         if (s->idx.xZ[2 * hashcode] == i && s->idx.xZ[2 * hashcode + 1] == j) {
            result = s->vals.xR[hashcode];
            return result;
         }
         hashcode = (hashcode + 1) % k;
      }
   }
   if (s->matrixtype == 1) {
   // CRS
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseGet: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      k0 = s->ridx.xZ[i];
      k1 = s->ridx.xZ[i + 1] - 1;
      result = 0.0;
      while (k0 <= k1) {
         k = (k0 + k1) / 2;
         if (s->idx.xZ[k] == j) {
            result = s->vals.xR[k];
            return result;
         }
         if (s->idx.xZ[k] < j) {
            k0 = k + 1;
         } else {
            k1 = k - 1;
         }
      }
      return result;
   }
   if (s->matrixtype == 2) {
   // SKS
      ae_assert(s->m == s->n, "SparseGet: non-square SKS matrix not supported");
      result = 0.0;
      if (i == j) {
      // Return diagonal element
         result = s->vals.xR[s->ridx.xZ[i] + s->didx.xZ[i]];
         return result;
      }
      if (j < i) {
      // Return subdiagonal element at I-th "skyline block"
         k = s->didx.xZ[i];
         if (i - j <= k) {
            result = s->vals.xR[s->ridx.xZ[i] + k + j - i];
         }
      } else {
      // Return superdiagonal element at J-th "skyline block"
         k = s->uidx.xZ[j];
         if (j - i <= k) {
            result = s->vals.xR[s->ridx.xZ[j + 1] - (j - i)];
         }
         return result;
      }
      return result;
   }
   ae_assert(false, "SparseGet: unexpected matrix type");
   return result;
}

// This function returns I-th diagonal element of the sparse matrix.
//
// Matrix can be in any mode (Hash-Table or CRS storage), but this  function
// is most efficient for CRS matrices - it requires less than 50 CPU  cycles
// to extract diagonal element. For Hash-Table matrices we still  have  O(1)
// query time, but function is many times slower.
//
// Inputs:
//     S           -   sparse M*N matrix in Hash-Table representation.
//                     Exception will be thrown for CRS matrix.
//     I           -   index of the element to modify, 0 <= I < min(M,N)
//
// Result:
//     value of S[I,I] or zero (in case no element with such index is found)
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: double sparsegetdiagonal(const sparsematrix &s, const ae_int_t i);
double sparsegetdiagonal(sparsematrix *s, ae_int_t i) {
   double result;
   ae_assert(i >= 0, "SparseGetDiagonal: I<0");
   ae_assert(i < s->m, "SparseGetDiagonal: I >= M");
   ae_assert(i < s->n, "SparseGetDiagonal: I >= N");
   result = 0.0;
   if (s->matrixtype == 0) {
      result = sparseget(s, i, i);
      return result;
   }
   if (s->matrixtype == 1) {
      if (s->didx.xZ[i] != s->uidx.xZ[i]) {
         result = s->vals.xR[s->didx.xZ[i]];
      }
      return result;
   }
   if (s->matrixtype == 2) {
      ae_assert(s->m == s->n, "SparseGetDiagonal: non-square SKS matrix not supported");
      result = s->vals.xR[s->ridx.xZ[i] + s->didx.xZ[i]];
      return result;
   }
   ae_assert(false, "SparseGetDiagonal: unexpected matrix type");
   return result;
}

// This function calculates matrix-vector product  S*x.  Matrix  S  must  be
// stored in CRS or SKS format (exception will be thrown otherwise).
//
// Inputs:
//     S           -   sparse M*N matrix in CRS or SKS format.
//     X           -   array[N], input vector. For  performance  reasons  we
//                     make only quick checks - we check that array size  is
//                     at least N, but we do not check for NAN's or INF's.
//     Y           -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     Y           -   array[M], S*x
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsemv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y);
void sparsemv(sparsematrix *s, RVector *x, RVector *y) {
   double tval;
   double v;
   double vv;
   ae_int_t i;
   ae_int_t j;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_int_t n;
   ae_int_t m;
   ae_int_t d;
   ae_int_t u;
   ae_int_t ri;
   ae_int_t ri1;
   ae_assert(x->cnt >= s->n, "SparseMV: length(X)<N");
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseMV: incorrect matrix type (convert your matrix to CRS/SKS)");
   vectorsetlengthatleast(y, s->m);
   n = s->n;
   m = s->m;
   if (s->matrixtype == 1) {
   // CRS format.
   // Perform integrity check.
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseMV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
   // Try vendor kernels
      if (sparsegemvcrsmkl(0, s->m, s->n, 1.0, &s->vals, &s->idx, &s->ridx, x, 0, 0.0, y, 0)) {
         return;
      }
   // Our own implementation
      for (i = 0; i < m; i++) {
         tval = 0.0;
         lt = s->ridx.xZ[i];
         rt = s->ridx.xZ[i + 1] - 1;
         for (j = lt; j <= rt; j++) {
            tval += x->xR[s->idx.xZ[j]] * s->vals.xR[j];
         }
         y->xR[i] = tval;
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(s->m == s->n, "SparseMV: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         v = s->vals.xR[ri + d] * x->xR[i];
         if (d > 0) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
            v += vv;
         }
         y->xR[i] = v;
         if (u > 0) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            v = x->xR[i];
            ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
         }
      }
      return;
   }
}

// This function calculates matrix-vector product  S^T*x. Matrix S  must  be
// stored in CRS or SKS format (exception will be thrown otherwise).
//
// Inputs:
//     S           -   sparse M*N matrix in CRS or SKS format.
//     X           -   array[M], input vector. For  performance  reasons  we
//                     make only quick checks - we check that array size  is
//                     at least M, but we do not check for NAN's or INF's.
//     Y           -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     Y           -   array[N], S^T*x
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsemtv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y);
void sparsemtv(sparsematrix *s, RVector *x, RVector *y) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t ct;
   ae_int_t lt1;
   ae_int_t rt1;
   double v;
   double vv;
   ae_int_t n;
   ae_int_t m;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t d;
   ae_int_t u;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseMTV: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(x->cnt >= s->m, "SparseMTV: Length(X)<M");
   n = s->n;
   m = s->m;
   vectorsetlengthatleast(y, n);
   for (i = 0; i < n; i++) {
      y->xR[i] = 0.0;
   }
   if (s->matrixtype == 1) {
   // CRS format
   // Perform integrity check.
      ae_assert(s->ninitialized == s->ridx.xZ[m], "SparseMTV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
   // Try vendor kernels
      if (sparsegemvcrsmkl(1, s->m, s->n, 1.0, &s->vals, &s->idx, &s->ridx, x, 0, 0.0, y, 0)) {
         return;
      }
   // Our own implementation
      for (i = 0; i < m; i++) {
         lt = s->ridx.xZ[i];
         rt = s->ridx.xZ[i + 1];
         v = x->xR[i];
         for (j = lt; j < rt; j++) {
            ct = s->idx.xZ[j];
            y->xR[ct] += v * s->vals.xR[j];
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(s->m == s->n, "SparseMV: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         if (d > 0) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            v = x->xR[i];
            ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
         }
         v = s->vals.xR[ri + d] * x->xR[i];
         if (u > 0) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
            v += vv;
         }
         y->xR[i] = v;
      }
      return;
   }
}

// This function calculates generalized sparse matrix-vector product
//
//     y := alpha*op(S)*x + beta*y
//
// Matrix S must be stored in CRS or SKS format (exception  will  be  thrown
// otherwise). op(S) can be either S or S^T.
//
// NOTE: this  function  expects  Y  to  be  large enough to store result. No
//       automatic preallocation happens for smaller arrays.
//
// Inputs:
//     S           -   sparse matrix in CRS or SKS format.
//     Alpha       -   source coefficient
//     OpS         -   operation type:
//                     * OpS=0     =>  op(S) = S
//                     * OpS=1     =>  op(S) = S^T
//     X           -   input vector, must have at least Cols(op(S))+IX elements
//     IX          -   subvector offset
//     Beta        -   destination coefficient
//     Y           -   preallocated output array, must have at least Rows(op(S))+IY elements
//     IY          -   subvector offset
//
// Outputs:
//     Y           -   elements [IY...IY+Rows(op(S))-1] are replaced by result,
//                     other elements are not modified
//
// HANDLING OF SPECIAL CASES:
// * below M=Rows(op(S)) and N=Cols(op(S)). Although current  ALGLIB  version
//   does not allow you to  create  zero-sized  sparse  matrices,  internally
//   ALGLIB  can  deal  with  such matrices. So, comments for M or N equal to
//   zero are for internal use only.
// * if M=0, then subroutine does nothing. It does not even touch arrays.
// * if N=0 or Alpha=0.0, then:
//   * if Beta=0, then Y is filled by zeros. S and X are  not  referenced  at
//     all. Initial values of Y are ignored (we do not  multiply  Y by  zero,
//     we just rewrite it by zeros)
//   * if Beta != 0, then Y is replaced by Beta*Y
// * if M > 0, N > 0, Alpha != 0, but  Beta=0, then  Y is replaced by alpha*op(S)*x
//   initial state of Y  is ignored (rewritten without initial multiplication
//   by zeros).
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 10.12.2019 by Sergey Bochkanov
// API: void sparsegemv(const sparsematrix &s, const double alpha, const ae_int_t ops, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy);
void sparsegemv(sparsematrix *s, double alpha, ae_int_t ops, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
   ae_int_t opm;
   ae_int_t opn;
   ae_int_t rawm;
   ae_int_t rawn;
   ae_int_t i;
   ae_int_t j;
   double tval;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t ct;
   ae_int_t d;
   ae_int_t u;
   ae_int_t ri;
   ae_int_t ri1;
   double v;
   double vv;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_assert(ops == 0 || ops == 1, "SparseGEMV: incorrect OpS");
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseGEMV: incorrect matrix type (convert your matrix to CRS/SKS)");
   if (ops == 0) {
      opm = s->m;
      opn = s->n;
   } else {
      opm = s->n;
      opn = s->m;
   }
   ae_assert(opm >= 0 && opn >= 0, "SparseGEMV: op(S) has negative size");
   ae_assert(opn == 0 || x->cnt + ix >= opn, "SparseGEMV: X is too short");
   ae_assert(opm == 0 || y->cnt + iy >= opm, "SparseGEMV: X is too short");
   rawm = s->m;
   rawn = s->n;
// Quick exit strategies
   if (opm == 0) {
      return;
   }
   if (beta != 0.0) {
      for (i = 0; i < opm; i++) {
         y->xR[iy + i] *= beta;
      }
   } else {
      for (i = 0; i < opm; i++) {
         y->xR[iy + i] = 0.0;
      }
   }
   if (opn == 0 || alpha == 0.0) {
      return;
   }
// Now we have OpM >= 1, OpN >= 1, Alpha != 0
   if (ops == 0) {
   // Compute generalized product y := alpha*S*x + beta*y
   // (with "beta*y" part already computed).
      if (s->matrixtype == 1) {
      // CRS format.
      // Perform integrity check.
         ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseGEMV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      // Try vendor kernels
         if (sparsegemvcrsmkl(0, s->m, s->n, alpha, &s->vals, &s->idx, &s->ridx, x, ix, 1.0, y, iy)) {
            return;
         }
      // Our own implementation
         for (i = 0; i < rawm; i++) {
            tval = 0.0;
            lt = s->ridx.xZ[i];
            rt = s->ridx.xZ[i + 1] - 1;
            for (j = lt; j <= rt; j++) {
               tval += x->xR[s->idx.xZ[j] + ix] * s->vals.xR[j];
            }
            y->xR[i + iy] += alpha * tval;
         }
         return;
      }
      if (s->matrixtype == 2) {
      // SKS format
         ae_assert(s->m == s->n, "SparseMV: non-square SKS matrices are not supported");
         for (i = 0; i < rawn; i++) {
            ri = s->ridx.xZ[i];
            ri1 = s->ridx.xZ[i + 1];
            d = s->didx.xZ[i];
            u = s->uidx.xZ[i];
            v = s->vals.xR[ri + d] * x->xR[i + ix];
            if (d > 0) {
               lt = ri;
               rt = ri + d - 1;
               lt1 = i - d + ix;
               rt1 = i - 1 + ix;
               vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
               v += vv;
            }
            y->xR[i + iy] += alpha * v;
            if (u > 0) {
               lt = ri1 - u;
               rt = ri1 - 1;
               lt1 = i - u + iy;
               rt1 = i - 1 + iy;
               v = alpha * x->xR[i + ix];
               ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            }
         }
         return;
      }
   } else {
   // Compute generalized product y := alpha*S^T*x + beta*y
   // (with "beta*y" part already computed).
      if (s->matrixtype == 1) {
      // CRS format
      // Perform integrity check.
         ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseGEMV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      // Try vendor kernels
         if (sparsegemvcrsmkl(1, s->m, s->n, alpha, &s->vals, &s->idx, &s->ridx, x, ix, 1.0, y, iy)) {
            return;
         }
      // Our own implementation
         for (i = 0; i < rawm; i++) {
            lt = s->ridx.xZ[i];
            rt = s->ridx.xZ[i + 1];
            v = alpha * x->xR[i + ix];
            for (j = lt; j < rt; j++) {
               ct = s->idx.xZ[j] + iy;
               y->xR[ct] += v * s->vals.xR[j];
            }
         }
         return;
      }
      if (s->matrixtype == 2) {
      // SKS format
         ae_assert(s->m == s->n, "SparseGEMV: non-square SKS matrices are not supported");
         for (i = 0; i < rawn; i++) {
            ri = s->ridx.xZ[i];
            ri1 = s->ridx.xZ[i + 1];
            d = s->didx.xZ[i];
            u = s->uidx.xZ[i];
            if (d > 0) {
               lt = ri;
               rt = ri + d - 1;
               lt1 = i - d + iy;
               rt1 = i - 1 + iy;
               v = alpha * x->xR[i + ix];
               ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            }
            v = alpha * s->vals.xR[ri + d] * x->xR[i + ix];
            if (u > 0) {
               lt = ri1 - u;
               rt = ri1 - 1;
               lt1 = i - u + ix;
               rt1 = i - 1 + ix;
               vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
               v += alpha * vv;
            }
            y->xR[i + iy] += v;
         }
         return;
      }
   }
}

// This function simultaneously calculates two matrix-vector products:
//     S*x and S^T*x.
// S must be square (non-rectangular) matrix stored in  CRS  or  SKS  format
// (exception will be thrown otherwise).
//
// Inputs:
//     S           -   sparse N*N matrix in CRS or SKS format.
//     X           -   array[N], input vector. For  performance  reasons  we
//                     make only quick checks - we check that array size  is
//                     at least N, but we do not check for NAN's or INF's.
//     Y0          -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//     Y1          -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     Y0          -   array[N], S*x
//     Y1          -   array[N], S^T*x
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsemv2(const sparsematrix &s, const real_1d_array &x, real_1d_array &y0, real_1d_array &y1);
void sparsemv2(sparsematrix *s, RVector *x, RVector *y0, RVector *y1) {
   ae_int_t l;
   double tval;
   ae_int_t i;
   ae_int_t j;
   double vx;
   double vs;
   double v;
   double vv;
   double vd0;
   double vd1;
   ae_int_t vi;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t n;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t d;
   ae_int_t u;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseMV2: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(s->m == s->n, "SparseMV2: matrix is non-square");
   l = x->cnt;
   ae_assert(l >= s->n, "SparseMV2: Length(X)<N");
   n = s->n;
   vectorsetlengthatleast(y0, l);
   vectorsetlengthatleast(y1, l);
   for (i = 0; i < n; i++) {
      y0->xR[i] = 0.0;
      y1->xR[i] = 0.0;
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseMV2: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      for (i = 0; i < s->m; i++) {
         tval = 0.0;
         vx = x->xR[i];
         j0 = s->ridx.xZ[i];
         j1 = s->ridx.xZ[i + 1] - 1;
         for (j = j0; j <= j1; j++) {
            vi = s->idx.xZ[j];
            vs = s->vals.xR[j];
            tval += x->xR[vi] * vs;
            y1->xR[vi] += vx * vs;
         }
         y0->xR[i] = tval;
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         vd0 = s->vals.xR[ri + d] * x->xR[i];
         vd1 = vd0;
         if (d > 0) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            v = x->xR[i];
            ae_v_addd(&y1->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
            vd0 += vv;
         }
         if (u > 0) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            v = x->xR[i];
            ae_v_addd(&y0->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
            vd1 += vv;
         }
         y0->xR[i] = vd0;
         y1->xR[i] = vd1;
      }
      return;
   }
}

// This function calculates matrix-vector product  S*x, when S is  symmetric
// matrix. Matrix S  must be stored in CRS or SKS format  (exception will be
// thrown otherwise).
//
// Inputs:
//     S           -   sparse M*M matrix in CRS or SKS format.
//     IsUpper     -   whether upper or lower triangle of S is given:
//                     * if upper triangle is given,  only   S[i,j] for j >= i
//                       are used, and lower triangle is ignored (it can  be
//                       empty - these elements are not referenced at all).
//                     * if lower triangle is given,  only   S[i,j] for j <= i
//                       are used, and upper triangle is ignored.
//     X           -   array[N], input vector. For  performance  reasons  we
//                     make only quick checks - we check that array size  is
//                     at least N, but we do not check for NAN's or INF's.
//     Y           -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     Y           -   array[M], S*x
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsesmv(const sparsematrix &s, const bool isupper, const real_1d_array &x, real_1d_array &y);
void sparsesmv(sparsematrix *s, bool isupper, RVector *x, RVector *y) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t j;
   ae_int_t id;
   ae_int_t lt;
   ae_int_t rt;
   double v;
   double vv;
   double vy;
   double vx;
   double vd;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t d;
   ae_int_t u;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseSMV: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(x->cnt >= s->n, "SparseSMV: length(X)<N");
   ae_assert(s->m == s->n, "SparseSMV: non-square matrix");
   n = s->n;
   vectorsetlengthatleast(y, n);
   for (i = 0; i < n; i++) {
      y->xR[i] = 0.0;
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseSMV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      for (i = 0; i < n; i++) {
         if (s->didx.xZ[i] != s->uidx.xZ[i]) {
            y->xR[i] += s->vals.xR[s->didx.xZ[i]] * x->xR[s->idx.xZ[s->didx.xZ[i]]];
         }
         if (isupper) {
            lt = s->uidx.xZ[i];
            rt = s->ridx.xZ[i + 1];
            vy = 0.0;
            vx = x->xR[i];
            for (j = lt; j < rt; j++) {
               id = s->idx.xZ[j];
               v = s->vals.xR[j];
               vy += x->xR[id] * v;
               y->xR[id] += vx * v;
            }
            y->xR[i] += vy;
         } else {
            lt = s->ridx.xZ[i];
            rt = s->didx.xZ[i];
            vy = 0.0;
            vx = x->xR[i];
            for (j = lt; j < rt; j++) {
               id = s->idx.xZ[j];
               v = s->vals.xR[j];
               vy += x->xR[id] * v;
               y->xR[id] += vx * v;
            }
            y->xR[i] += vy;
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         vd = s->vals.xR[ri + d] * x->xR[i];
         if (d > 0 && !isupper) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            v = x->xR[i];
            ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
            vd += vv;
         }
         if (u > 0 && isupper) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            v = x->xR[i];
            ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            vv = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
            vd += vv;
         }
         y->xR[i] = vd;
      }
      return;
   }
}

// This function calculates vector-matrix-vector product x'*S*x, where  S is
// symmetric matrix. Matrix S must be stored in CRS or SKS format (exception
// will be thrown otherwise).
//
// Inputs:
//     S           -   sparse M*M matrix in CRS or SKS format.
//     IsUpper     -   whether upper or lower triangle of S is given:
//                     * if upper triangle is given,  only   S[i,j] for j >= i
//                       are used, and lower triangle is ignored (it can  be
//                       empty - these elements are not referenced at all).
//                     * if lower triangle is given,  only   S[i,j] for j <= i
//                       are used, and upper triangle is ignored.
//     X           -   array[N], input vector. For  performance  reasons  we
//                     make only quick checks - we check that array size  is
//                     at least N, but we do not check for NAN's or INF's.
//
// Result:
//     x'*S*x
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 27.01.2014 by Sergey Bochkanov
// API: double sparsevsmv(const sparsematrix &s, const bool isupper, const real_1d_array &x);
double sparsevsmv(sparsematrix *s, bool isupper, RVector *x) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t id;
   ae_int_t lt;
   ae_int_t rt;
   double v;
   double v0;
   double v1;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t d;
   ae_int_t u;
   ae_int_t lt1;
   double result;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseVSMV: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(x->cnt >= s->n, "SparseVSMV: length(X)<N");
   ae_assert(s->m == s->n, "SparseVSMV: non-square matrix");
   n = s->n;
   result = 0.0;
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseVSMV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      for (i = 0; i < n; i++) {
         if (s->didx.xZ[i] != s->uidx.xZ[i]) {
            v = x->xR[s->idx.xZ[s->didx.xZ[i]]];
            result += v * s->vals.xR[s->didx.xZ[i]] * v;
         }
         if (isupper) {
            lt = s->uidx.xZ[i];
            rt = s->ridx.xZ[i + 1];
         } else {
            lt = s->ridx.xZ[i];
            rt = s->didx.xZ[i];
         }
         v0 = x->xR[i];
         for (j = lt; j < rt; j++) {
            id = s->idx.xZ[j];
            v1 = x->xR[id];
            v = s->vals.xR[j];
            result += 2 * v0 * v1 * v;
         }
      }
      return result;
   }
   if (s->matrixtype == 2) {
   // SKS format
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         v = x->xR[i];
         result += v * s->vals.xR[ri + d] * v;
         if (d > 0 && !isupper) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            k = d - 1;
            v0 = x->xR[i];
            v = 0.0;
            for (j = 0; j <= k; j++) {
               v += x->xR[lt1 + j] * s->vals.xR[lt + j];
            }
            result += 2 * v0 * v;
         }
         if (u > 0 && isupper) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            k = u - 1;
            v0 = x->xR[i];
            v = 0.0;
            for (j = 0; j <= k; j++) {
               v += x->xR[lt1 + j] * s->vals.xR[lt + j];
            }
            result += 2 * v0 * v;
         }
      }
      return result;
   }
   return result;
}

// This function calculates matrix-matrix product  S*A.  Matrix  S  must  be
// stored in CRS or SKS format (exception will be thrown otherwise).
//
// Inputs:
//     S           -   sparse M*N matrix in CRS or SKS format.
//     A           -   array[N][K], input dense matrix. For  performance reasons
//                     we make only quick checks - we check that array size
//                     is at least N, but we do not check for NAN's or INF's.
//     K           -   number of columns of matrix (A).
//     B           -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     B           -   array[M][K], S*A
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsemm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b);
void sparsemm(sparsematrix *s, RMatrix *a, ae_int_t k, RMatrix *b) {
   double tval;
   double v;
   ae_int_t id;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k0;
   ae_int_t k1;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t m;
   ae_int_t n;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_int_t d;
   ae_int_t u;
   double vd;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseMM: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(a->rows >= s->n, "SparseMM: Rows(A)<N");
   ae_assert(k > 0, "SparseMM: K <= 0");
   m = s->m;
   n = s->n;
   k1 = k - 1;
   matrixsetlengthatleast(b, m, k);
   for (i = 0; i < m; i++) {
      for (j = 0; j < k; j++) {
         b->xyR[i][j] = 0.0;
      }
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[m], "SparseMM: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      if (k < sparse_linalgswitch) {
         for (i = 0; i < m; i++) {
            for (j = 0; j < k; j++) {
               tval = 0.0;
               lt = s->ridx.xZ[i];
               rt = s->ridx.xZ[i + 1];
               for (k0 = lt; k0 < rt; k0++) {
                  tval += s->vals.xR[k0] * a->xyR[s->idx.xZ[k0]][j];
               }
               b->xyR[i][j] = tval;
            }
         }
      } else {
         for (i = 0; i < m; i++) {
            lt = s->ridx.xZ[i];
            rt = s->ridx.xZ[i + 1];
            for (j = lt; j < rt; j++) {
               id = s->idx.xZ[j];
               v = s->vals.xR[j];
               ae_v_addd(b->xyR[i], 1, a->xyR[id], 1, k, v);
            }
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(m == n, "SparseMM: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         if (d > 0) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b->xyR[i][k0] += v * a->xyR[j][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b->xyR[i], 1, a->xyR[j], 1, k, v);
               }
            }
         }
         if (u > 0) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b->xyR[j][k0] += v * a->xyR[i][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b->xyR[j], 1, a->xyR[i], 1, k, v);
               }
            }
         }
         vd = s->vals.xR[ri + d];
         ae_v_addd(b->xyR[i], 1, a->xyR[i], 1, k, vd);
      }
      return;
   }
}

// This function calculates matrix-matrix product  S^T*A. Matrix S  must  be
// stored in CRS or SKS format (exception will be thrown otherwise).
//
// Inputs:
//     S           -   sparse M*N matrix in CRS or SKS format.
//     A           -   array[M][K], input dense matrix. For performance reasons
//                     we make only quick checks - we check that array size  is
//                     at least M, but we do not check for NAN's or INF's.
//     K           -   number of columns of matrix (A).
//     B           -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     B           -   array[N][K], S^T*A
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsemtm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b);
void sparsemtm(sparsematrix *s, RMatrix *a, ae_int_t k, RMatrix *b) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k0;
   ae_int_t k1;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t ct;
   double v;
   ae_int_t m;
   ae_int_t n;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_int_t d;
   ae_int_t u;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseMTM: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(a->rows >= s->m, "SparseMTM: Rows(A)<M");
   ae_assert(k > 0, "SparseMTM: K <= 0");
   m = s->m;
   n = s->n;
   k1 = k - 1;
   matrixsetlengthatleast(b, n, k);
   for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
         b->xyR[i][j] = 0.0;
      }
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[m], "SparseMTM: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      if (k < sparse_linalgswitch) {
         for (i = 0; i < m; i++) {
            lt = s->ridx.xZ[i];
            rt = s->ridx.xZ[i + 1];
            for (k0 = lt; k0 < rt; k0++) {
               v = s->vals.xR[k0];
               ct = s->idx.xZ[k0];
               for (j = 0; j < k; j++) {
                  b->xyR[ct][j] += v * a->xyR[i][j];
               }
            }
         }
      } else {
         for (i = 0; i < m; i++) {
            lt = s->ridx.xZ[i];
            rt = s->ridx.xZ[i + 1];
            for (j = lt; j < rt; j++) {
               v = s->vals.xR[j];
               ct = s->idx.xZ[j];
               ae_v_addd(b->xyR[ct], 1, a->xyR[i], 1, k, v);
            }
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(m == n, "SparseMTM: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         if (d > 0) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b->xyR[j][k0] += v * a->xyR[i][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b->xyR[j], 1, a->xyR[i], 1, k, v);
               }
            }
         }
         if (u > 0) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b->xyR[i][k0] += v * a->xyR[j][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b->xyR[i], 1, a->xyR[j], 1, k, v);
               }
            }
         }
         v = s->vals.xR[ri + d];
         ae_v_addd(b->xyR[i], 1, a->xyR[i], 1, k, v);
      }
      return;
   }
}

// This function simultaneously calculates two matrix-matrix products:
//     S*A and S^T*A.
// S  must  be  square (non-rectangular) matrix stored in CRS or  SKS  format
// (exception will be thrown otherwise).
//
// Inputs:
//     S           -   sparse N*N matrix in CRS or SKS format.
//     A           -   array[N][K], input dense matrix. For performance reasons
//                     we make only quick checks - we check that array size  is
//                     at least N, but we do not check for NAN's or INF's.
//     K           -   number of columns of matrix (A).
//     B0          -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//     B1          -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     B0          -   array[N][K], S*A
//     B1          -   array[N][K], S^T*A
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsemm2(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b0, real_2d_array &b1);
void sparsemm2(sparsematrix *s, RMatrix *a, ae_int_t k, RMatrix *b0, RMatrix *b1) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k0;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t ct;
   double v;
   double tval;
   ae_int_t n;
   ae_int_t k1;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_int_t d;
   ae_int_t u;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseMM2: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(s->m == s->n, "SparseMM2: matrix is non-square");
   ae_assert(a->rows >= s->n, "SparseMM2: Rows(A)<N");
   ae_assert(k > 0, "SparseMM2: K <= 0");
   n = s->n;
   k1 = k - 1;
   matrixsetlengthatleast(b0, n, k);
   matrixsetlengthatleast(b1, n, k);
   for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
         b1->xyR[i][j] = 0.0;
         b0->xyR[i][j] = 0.0;
      }
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseMM2: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      if (k < sparse_linalgswitch) {
         for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
               tval = 0.0;
               lt = s->ridx.xZ[i];
               rt = s->ridx.xZ[i + 1];
               v = a->xyR[i][j];
               for (k0 = lt; k0 < rt; k0++) {
                  ct = s->idx.xZ[k0];
                  b1->xyR[ct][j] += s->vals.xR[k0] * v;
                  tval += s->vals.xR[k0] * a->xyR[ct][j];
               }
               b0->xyR[i][j] = tval;
            }
         }
      } else {
         for (i = 0; i < n; i++) {
            lt = s->ridx.xZ[i];
            rt = s->ridx.xZ[i + 1];
            for (j = lt; j < rt; j++) {
               v = s->vals.xR[j];
               ct = s->idx.xZ[j];
               ae_v_addd(b0->xyR[i], 1, a->xyR[ct], 1, k, v);
               ae_v_addd(b1->xyR[ct], 1, a->xyR[i], 1, k, v);
            }
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(s->m == s->n, "SparseMM2: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         if (d > 0) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b0->xyR[i][k0] += v * a->xyR[j][k0];
                     b1->xyR[j][k0] += v * a->xyR[i][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b0->xyR[i], 1, a->xyR[j], 1, k, v);
                  ae_v_addd(b1->xyR[j], 1, a->xyR[i], 1, k, v);
               }
            }
         }
         if (u > 0) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b0->xyR[j][k0] += v * a->xyR[i][k0];
                     b1->xyR[i][k0] += v * a->xyR[j][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b0->xyR[j], 1, a->xyR[i], 1, k, v);
                  ae_v_addd(b1->xyR[i], 1, a->xyR[j], 1, k, v);
               }
            }
         }
         v = s->vals.xR[ri + d];
         ae_v_addd(b0->xyR[i], 1, a->xyR[i], 1, k, v);
         ae_v_addd(b1->xyR[i], 1, a->xyR[i], 1, k, v);
      }
      return;
   }
}

// This function calculates matrix-matrix product  S*A, when S  is  symmetric
// matrix. Matrix S must be stored in CRS or SKS format  (exception  will  be
// thrown otherwise).
//
// Inputs:
//     S           -   sparse M*M matrix in CRS or SKS format.
//     IsUpper     -   whether upper or lower triangle of S is given:
//                     * if upper triangle is given,  only   S[i,j] for j >= i
//                       are used, and lower triangle is ignored (it can  be
//                       empty - these elements are not referenced at all).
//                     * if lower triangle is given,  only   S[i,j] for j <= i
//                       are used, and upper triangle is ignored.
//     A           -   array[M][K], input dense matrix. For performance reasons
//                     we make only quick checks - we check that array size is
//                     at least M, but we do not check for NAN's or INF's.
//     K           -   number of columns of matrix (A).
//     B           -   output buffer, possibly preallocated. In case  buffer
//                     size is too small to store  result,  this  buffer  is
//                     automatically resized.
//
// Outputs:
//     B           -   array[M][K], S*A
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparsesmm(const sparsematrix &s, const bool isupper, const real_2d_array &a, const ae_int_t k, real_2d_array &b);
void sparsesmm(sparsematrix *s, bool isupper, RMatrix *a, ae_int_t k, RMatrix *b) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k0;
   ae_int_t id;
   ae_int_t k1;
   ae_int_t lt;
   ae_int_t rt;
   double v;
   double vb;
   double va;
   ae_int_t n;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_int_t d;
   ae_int_t u;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseSMM: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(a->rows >= s->n, "SparseSMM: Rows(X)<N");
   ae_assert(s->m == s->n, "SparseSMM: matrix is non-square");
   n = s->n;
   k1 = k - 1;
   matrixsetlengthatleast(b, n, k);
   for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
         b->xyR[i][j] = 0.0;
      }
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseSMM: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      if (k > sparse_linalgswitch) {
         for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
               if (s->didx.xZ[i] != s->uidx.xZ[i]) {
                  id = s->didx.xZ[i];
                  b->xyR[i][j] += s->vals.xR[id] * a->xyR[s->idx.xZ[id]][j];
               }
               if (isupper) {
                  lt = s->uidx.xZ[i];
                  rt = s->ridx.xZ[i + 1];
                  vb = 0.0;
                  va = a->xyR[i][j];
                  for (k0 = lt; k0 < rt; k0++) {
                     id = s->idx.xZ[k0];
                     v = s->vals.xR[k0];
                     vb += a->xyR[id][j] * v;
                     b->xyR[id][j] += va * v;
                  }
                  b->xyR[i][j] += vb;
               } else {
                  lt = s->ridx.xZ[i];
                  rt = s->didx.xZ[i];
                  vb = 0.0;
                  va = a->xyR[i][j];
                  for (k0 = lt; k0 < rt; k0++) {
                     id = s->idx.xZ[k0];
                     v = s->vals.xR[k0];
                     vb += a->xyR[id][j] * v;
                     b->xyR[id][j] += va * v;
                  }
                  b->xyR[i][j] += vb;
               }
            }
         }
      } else {
         for (i = 0; i < n; i++) {
            if (s->didx.xZ[i] != s->uidx.xZ[i]) {
               id = s->didx.xZ[i];
               v = s->vals.xR[id];
               ae_v_addd(b->xyR[i], 1, a->xyR[s->idx.xZ[id]], 1, k, v);
            }
            if (isupper) {
               lt = s->uidx.xZ[i];
               rt = s->ridx.xZ[i + 1];
               for (j = lt; j < rt; j++) {
                  id = s->idx.xZ[j];
                  v = s->vals.xR[j];
                  ae_v_addd(b->xyR[i], 1, a->xyR[id], 1, k, v);
                  ae_v_addd(b->xyR[id], 1, a->xyR[i], 1, k, v);
               }
            } else {
               lt = s->ridx.xZ[i];
               rt = s->didx.xZ[i];
               for (j = lt; j < rt; j++) {
                  id = s->idx.xZ[j];
                  v = s->vals.xR[j];
                  ae_v_addd(b->xyR[i], 1, a->xyR[id], 1, k, v);
                  ae_v_addd(b->xyR[id], 1, a->xyR[i], 1, k, v);
               }
            }
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(s->m == s->n, "SparseMM2: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         if (d > 0 && !isupper) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b->xyR[i][k0] += v * a->xyR[j][k0];
                     b->xyR[j][k0] += v * a->xyR[i][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b->xyR[i], 1, a->xyR[j], 1, k, v);
                  ae_v_addd(b->xyR[j], 1, a->xyR[i], 1, k, v);
               }
            }
         }
         if (u > 0 && isupper) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            for (j = lt1; j <= rt1; j++) {
               v = s->vals.xR[lt + (j - lt1)];
               if (k < sparse_linalgswitch) {
               // Use loop
                  for (k0 = 0; k0 <= k1; k0++) {
                     b->xyR[j][k0] += v * a->xyR[i][k0];
                     b->xyR[i][k0] += v * a->xyR[j][k0];
                  }
               } else {
               // Use vector operation
                  ae_v_addd(b->xyR[j], 1, a->xyR[i], 1, k, v);
                  ae_v_addd(b->xyR[i], 1, a->xyR[j], 1, k, v);
               }
            }
         }
         v = s->vals.xR[ri + d];
         ae_v_addd(b->xyR[i], 1, a->xyR[i], 1, k, v);
      }
      return;
   }
}

// This function calculates matrix-vector product op(S)*x, when x is  vector,
// S is symmetric triangular matrix, op(S) is transposition or no  operation.
// Matrix S must be stored in CRS or SKS format  (exception  will  be  thrown
// otherwise).
//
// Inputs:
//     S           -   sparse square matrix in CRS or SKS format.
//     IsUpper     -   whether upper or lower triangle of S is used:
//                     * if upper triangle is given,  only   S[i,j] for  j >= i
//                       are used, and lower triangle is  ignored (it can  be
//                       empty - these elements are not referenced at all).
//                     * if lower triangle is given,  only   S[i,j] for  j <= i
//                       are used, and upper triangle is ignored.
//     IsUnit      -   unit or non-unit diagonal:
//                     * if True, diagonal elements of triangular matrix  are
//                       considered equal to 1.0. Actual elements  stored  in
//                       S are not referenced at all.
//                     * if False, diagonal stored in S is used
//     OpType      -   operation type:
//                     * if 0, S*x is calculated
//                     * if 1, (S^T)*x is calculated (transposition)
//     X           -   array[N] which stores input  vector.  For  performance
//                     reasons we make only quick  checks  -  we  check  that
//                     array  size  is  at  least  N, but we do not check for
//                     NAN's or INF's.
//     Y           -   possibly  preallocated  input   buffer.  Automatically
//                     resized if its size is too small.
//
// Outputs:
//     Y           -   array[N], op(S)*x
//
// NOTE: this function throws exception when called for non-CRS/SKS  matrix.
// You must convert your matrix with SparseConvertToCRS/SKS()  before  using
// this function.
//
// ALGLIB Project: Copyright 20.01.2014 by Sergey Bochkanov
// API: void sparsetrmv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, real_1d_array &y);
void sparsetrmv(sparsematrix *s, bool isupper, bool isunit, ae_int_t optype, RVector *x, RVector *y) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t j0;
   ae_int_t j1;
   double v;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t d;
   ae_int_t u;
   ae_int_t lt;
   ae_int_t rt;
   ae_int_t lt1;
   ae_int_t rt1;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseTRMV: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(optype == 0 || optype == 1, "SparseTRMV: incorrect operation type (must be 0 or 1)");
   ae_assert(x->cnt >= s->n, "SparseTRMV: Length(X)<N");
   ae_assert(s->m == s->n, "SparseTRMV: matrix is non-square");
   n = s->n;
   vectorsetlengthatleast(y, n);
   if (isunit) {
   // Set initial value of y to x
      for (i = 0; i < n; i++) {
         y->xR[i] = x->xR[i];
      }
   } else {
   // Set initial value of y to 0
      for (i = 0; i < n; i++) {
         y->xR[i] = 0.0;
      }
   }
   if (s->matrixtype == 1) {
   // CRS format
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseTRMV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      for (i = 0; i < n; i++) {
      // Depending on IsUpper/IsUnit, select range of indexes to process
         if (isupper) {
            if (isunit || s->didx.xZ[i] == s->uidx.xZ[i]) {
               j0 = s->uidx.xZ[i];
            } else {
               j0 = s->didx.xZ[i];
            }
            j1 = s->ridx.xZ[i + 1] - 1;
         } else {
            j0 = s->ridx.xZ[i];
            if (isunit || s->didx.xZ[i] == s->uidx.xZ[i]) {
               j1 = s->didx.xZ[i] - 1;
            } else {
               j1 = s->didx.xZ[i];
            }
         }
      // Depending on OpType, process subset of I-th row of input matrix
         if (optype == 0) {
            v = 0.0;
            for (j = j0; j <= j1; j++) {
               v += s->vals.xR[j] * x->xR[s->idx.xZ[j]];
            }
            y->xR[i] += v;
         } else {
            v = x->xR[i];
            for (j = j0; j <= j1; j++) {
               k = s->idx.xZ[j];
               y->xR[k] += v * s->vals.xR[j];
            }
         }
      }
      return;
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(s->m == s->n, "SparseTRMV: non-square SKS matrices are not supported");
      for (i = 0; i < n; i++) {
         ri = s->ridx.xZ[i];
         ri1 = s->ridx.xZ[i + 1];
         d = s->didx.xZ[i];
         u = s->uidx.xZ[i];
         if (!isunit) {
            y->xR[i] += s->vals.xR[ri + d] * x->xR[i];
         }
         if (d > 0 && !isupper) {
            lt = ri;
            rt = ri + d - 1;
            lt1 = i - d;
            rt1 = i - 1;
            if (optype == 0) {
               v = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
               y->xR[i] += v;
            } else {
               v = x->xR[i];
               ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            }
         }
         if (u > 0 && isupper) {
            lt = ri1 - u;
            rt = ri1 - 1;
            lt1 = i - u;
            rt1 = i - 1;
            if (optype == 0) {
               v = x->xR[i];
               ae_v_addd(&y->xR[lt1], 1, &s->vals.xR[lt], 1, rt1 - lt1 + 1, v);
            } else {
               v = ae_v_dotproduct(&s->vals.xR[lt], 1, &x->xR[lt1], 1, rt - lt + 1);
               y->xR[i] += v;
            }
         }
      }
      return;
   }
}

// This function solves linear system op(S)*y=x  where  x  is  vector,  S  is
// symmetric  triangular  matrix,  op(S)  is  transposition  or no operation.
// Matrix S must be stored in CRS or SKS format  (exception  will  be  thrown
// otherwise).
//
// Inputs:
//     S           -   sparse square matrix in CRS or SKS format.
//     IsUpper     -   whether upper or lower triangle of S is used:
//                     * if upper triangle is given,  only   S[i,j] for  j >= i
//                       are used, and lower triangle is  ignored (it can  be
//                       empty - these elements are not referenced at all).
//                     * if lower triangle is given,  only   S[i,j] for  j <= i
//                       are used, and upper triangle is ignored.
//     IsUnit      -   unit or non-unit diagonal:
//                     * if True, diagonal elements of triangular matrix  are
//                       considered equal to 1.0. Actual elements  stored  in
//                       S are not referenced at all.
//                     * if False, diagonal stored in S is used. It  is  your
//                       responsibility  to  make  sure  that   diagonal   is
//                       non-zero.
//     OpType      -   operation type:
//                     * if 0, S*x is calculated
//                     * if 1, (S^T)*x is calculated (transposition)
//     X           -   array[N] which stores input  vector.  For  performance
//                     reasons we make only quick  checks  -  we  check  that
//                     array  size  is  at  least  N, but we do not check for
//                     NAN's or INF's.
//
// Outputs:
//     X           -   array[N], inv(op(S))*x
//
// NOTE: this function throws exception when called for  non-CRS/SKS  matrix.
//       You must convert your matrix  with  SparseConvertToCRS/SKS()  before
//       using this function.
//
// NOTE: no assertion or tests are done during algorithm  operation.   It  is
//       your responsibility to provide invertible matrix to algorithm.
//
// ALGLIB Project: Copyright 20.01.2014 by Sergey Bochkanov
// API: void sparsetrsv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x);
void sparsetrsv(sparsematrix *s, bool isupper, bool isunit, ae_int_t optype, RVector *x) {
   ae_int_t n;
   ae_int_t fst;
   ae_int_t lst;
   ae_int_t stp;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double v;
   double vd;
   double v0;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t ri;
   ae_int_t ri1;
   ae_int_t d;
   ae_int_t u;
   ae_int_t lt;
   ae_int_t lt1;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseTRSV: incorrect matrix type (convert your matrix to CRS/SKS)");
   ae_assert(optype == 0 || optype == 1, "SparseTRSV: incorrect operation type (must be 0 or 1)");
   ae_assert(x->cnt >= s->n, "SparseTRSV: Length(X)<N");
   ae_assert(s->m == s->n, "SparseTRSV: matrix is non-square");
   n = s->n;
   if (s->matrixtype == 1) {
   // CRS format.
   //
   // Several branches for different combinations of IsUpper and OpType
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseTRSV: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      if (optype == 0) {
      // No transposition.
      //
      // S*x=y with upper or lower triangular S.
         v0 = 0.0;
         if (isupper) {
            fst = n - 1;
            lst = 0;
            stp = -1;
         } else {
            fst = 0;
            lst = n - 1;
            stp = 1;
         }
         i = fst;
         while (stp > 0 && i <= lst || stp < 0 && i >= lst) {
         // Select range of indexes to process
            if (isupper) {
               j0 = s->uidx.xZ[i];
               j1 = s->ridx.xZ[i + 1] - 1;
            } else {
               j0 = s->ridx.xZ[i];
               j1 = s->didx.xZ[i] - 1;
            }
         // Calculate X[I]
            v = 0.0;
            for (j = j0; j <= j1; j++) {
               v += s->vals.xR[j] * x->xR[s->idx.xZ[j]];
            }
            if (!isunit) {
               if (s->didx.xZ[i] == s->uidx.xZ[i]) {
                  vd = 0.0;
               } else {
                  vd = s->vals.xR[s->didx.xZ[i]];
               }
            } else {
               vd = 1.0;
            }
            v = (x->xR[i] - v) / vd;
            x->xR[i] = v;
            v0 = 0.25 * v0 + v;
         // Next I
            i += stp;
         }
         ae_assert(isfinite(v0), "SparseTRSV: overflow or division by exact zero");
         return;
      }
      if (optype == 1) {
      // Transposition.
      //
      // (S^T)*x=y with upper or lower triangular S.
         if (isupper) {
            fst = 0;
            lst = n - 1;
            stp = 1;
         } else {
            fst = n - 1;
            lst = 0;
            stp = -1;
         }
         i = fst;
         v0 = 0.0;
         while (stp > 0 && i <= lst || stp < 0 && i >= lst) {
            v = x->xR[i];
            if (v != 0.0) {
            // X[i] already stores A[i,i]*Y[i], the only thing left
            // is to divide by diagonal element.
               if (!isunit) {
                  if (s->didx.xZ[i] == s->uidx.xZ[i]) {
                     vd = 0.0;
                  } else {
                     vd = s->vals.xR[s->didx.xZ[i]];
                  }
               } else {
                  vd = 1.0;
               }
               v /= vd;
               x->xR[i] = v;
               v0 = 0.25 * v0 + v;
            // For upper triangular case:
            //     subtract X[i]*Ai from X[i+1:N-1]
            //
            // For lower triangular case:
            //     subtract X[i]*Ai from X[0:i-1]
            //
            // (here Ai is I-th row of original, untransposed A).
               if (isupper) {
                  j0 = s->uidx.xZ[i];
                  j1 = s->ridx.xZ[i + 1] - 1;
               } else {
                  j0 = s->ridx.xZ[i];
                  j1 = s->didx.xZ[i] - 1;
               }
               for (j = j0; j <= j1; j++) {
                  k = s->idx.xZ[j];
                  x->xR[k] -= s->vals.xR[j] * v;
               }
            }
         // Next I
            i += stp;
         }
         ae_assert(isfinite(v0), "SparseTRSV: overflow or division by exact zero");
         return;
      }
      ae_assert(false, "SparseTRSV: internal error");
   }
   if (s->matrixtype == 2) {
   // SKS format
      ae_assert(s->m == s->n, "SparseTRSV: non-square SKS matrices are not supported");
      if (optype == (isupper? 1: 0)) {
      // Lower triangular op(S) (matrix itself can be upper triangular).
         v0 = 0.0;
         for (i = 0; i < n; i++) {
         // Select range of indexes to process
            ri = s->ridx.xZ[i];
            ri1 = s->ridx.xZ[i + 1];
            d = s->didx.xZ[i];
            u = s->uidx.xZ[i];
            if (isupper) {
               lt = i - u;
               lt1 = ri1 - u;
               k = u - 1;
            } else {
               lt = i - d;
               lt1 = ri;
               k = d - 1;
            }
         // Calculate X[I]
            v = 0.0;
            for (j = 0; j <= k; j++) {
               v += s->vals.xR[lt1 + j] * x->xR[lt + j];
            }
            if (isunit) {
               vd = 1.0;
            } else {
               vd = s->vals.xR[ri + d];
            }
            v = (x->xR[i] - v) / vd;
            x->xR[i] = v;
            v0 = 0.25 * v0 + v;
         }
         ae_assert(isfinite(v0), "SparseTRSV: overflow or division by exact zero");
         return;
      }
      if (optype == (isupper? 0: 1)) {
      // Upper triangular op(S) (matrix itself can be lower triangular).
         v0 = 0.0;
         for (i = n - 1; i >= 0; i--) {
            ri = s->ridx.xZ[i];
            ri1 = s->ridx.xZ[i + 1];
            d = s->didx.xZ[i];
            u = s->uidx.xZ[i];
         // X[i] already stores A[i,i]*Y[i], the only thing left
         // is to divide by diagonal element.
            if (isunit) {
               vd = 1.0;
            } else {
               vd = s->vals.xR[ri + d];
            }
            v = x->xR[i] / vd;
            x->xR[i] = v;
            v0 = 0.25 * v0 + v;
         // Subtract product of X[i] and I-th column of "effective" A from
         // unprocessed variables.
            v = x->xR[i];
            if (isupper) {
               lt = i - u;
               lt1 = ri1 - u;
               k = u - 1;
            } else {
               lt = i - d;
               lt1 = ri;
               k = d - 1;
            }
            for (j = 0; j <= k; j++) {
               x->xR[lt + j] -= v * s->vals.xR[lt1 + j];
            }
         }
         ae_assert(isfinite(v0), "SparseTRSV: overflow or division by exact zero");
         return;
      }
      ae_assert(false, "SparseTRSV: internal error");
   }
   ae_assert(false, "SparseTRSV: internal error");
}

// This procedure resizes Hash-Table matrix. It can be called when you  have
// deleted too many elements from the matrix, and you want to  free unneeded
// memory.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparseresizematrix(const sparsematrix &s);
void sparseresizematrix(sparsematrix *s) {
   ae_frame _frame_block;
   ae_int_t k;
   ae_int_t k1;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(tvals, 0, DT_REAL);
   NewVector(tidx, 0, DT_INT);
   ae_assert(s->matrixtype == 0, "SparseResizeMatrix: incorrect matrix type");
// Initialization for length and number of non-null elementd
   k = s->tablesize;
   k1 = 0;
// Calculating number of non-null elements
   for (i = 0; i < k; i++) {
      if (s->idx.xZ[2 * i] >= 0) {
         k1++;
      }
   }
// Initialization value for free space
   s->tablesize = RoundZ(k1 / sparse_desiredloadfactor * sparse_growfactor + sparse_additional);
   s->nfree = s->tablesize - k1;
   ae_vector_set_length(&tvals, s->tablesize);
   ae_vector_set_length(&tidx, 2 * s->tablesize);
   ae_swap_vectors(&s->vals, &tvals);
   ae_swap_vectors(&s->idx, &tidx);
   for (i = 0; i < s->tablesize; i++) {
      s->idx.xZ[2 * i] = -1;
   }
   for (i = 0; i < k; i++) {
      if (tidx.xZ[2 * i] >= 0) {
         sparseset(s, tidx.xZ[2 * i], tidx.xZ[2 * i + 1], tvals.xR[i]);
      }
   }
   ae_frame_leave();
}

// Procedure for initialization 'S.DIdx' and 'S.UIdx'
//
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
void sparseinitduidx(sparsematrix *s) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t lt;
   ae_int_t rt;
   ae_assert(s->matrixtype == 1, "SparseInitDUIdx: internal error, incorrect matrix type");
   vectorsetlengthatleast(&s->didx, s->m);
   vectorsetlengthatleast(&s->uidx, s->m);
   for (i = 0; i < s->m; i++) {
      s->uidx.xZ[i] = -1;
      s->didx.xZ[i] = -1;
      lt = s->ridx.xZ[i];
      rt = s->ridx.xZ[i + 1];
      for (j = lt; j < rt; j++) {
         k = s->idx.xZ[j];
         if (k == i) {
            s->didx.xZ[i] = j;
         } else {
            if (k > i && s->uidx.xZ[i] == -1) {
               s->uidx.xZ[i] = j;
               break;
            }
         }
      }
      if (s->uidx.xZ[i] == -1) {
         s->uidx.xZ[i] = s->ridx.xZ[i + 1];
      }
      if (s->didx.xZ[i] == -1) {
         s->didx.xZ[i] = s->uidx.xZ[i];
      }
   }
}

// This function return average length of chain at hash-table.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
double sparsegetaveragelengthofchain(sparsematrix *s) {
   ae_int_t nchains;
   ae_int_t talc;
   ae_int_t l;
   ae_int_t i;
   ae_int_t ind0;
   ae_int_t ind1;
   ae_int_t hashcode;
   double result;
// If matrix represent in CRS then return zero and exit
   if (s->matrixtype != 0) {
      result = 0.0;
      return result;
   }
   nchains = 0;
   talc = 0;
   l = s->tablesize;
   for (i = 0; i < l; i++) {
      ind0 = 2 * i;
      if (s->idx.xZ[ind0] != -1) {
         nchains++;
         hashcode = sparse_hash(s->idx.xZ[ind0], s->idx.xZ[ind0 + 1], l);
         while (true) {
            talc++;
            ind1 = 2 * hashcode;
            if (s->idx.xZ[ind0] == s->idx.xZ[ind1] && s->idx.xZ[ind0 + 1] == s->idx.xZ[ind1 + 1]) {
               break;
            }
            hashcode = (hashcode + 1) % l;
         }
      }
   }
   if (nchains == 0) {
      result = 0.0;
   } else {
      result = (double)talc / (double)nchains;
   }
   return result;
}

// This  function  is  used  to enumerate all elements of the sparse matrix.
// Before  first  call  user  initializes  T0 and T1 counters by zero. These
// counters are used to remember current position in a  matrix;  after  each
// call they are updated by the function.
//
// Subsequent calls to this function return non-zero elements of the  sparse
// matrix, one by one. If you enumerate CRS matrix, matrix is traversed from
// left to right, from top to bottom. In case you enumerate matrix stored as
// Hash table, elements are returned in random order.
//
// EXAMPLE
//     > T0=0
//     > T1=0
//     > while SparseEnumerate(S,T0,T1,I,J,V) do
//     >     ....do something with I,J,V
//
// Inputs:
//     S           -   sparse M*N matrix in Hash-Table or CRS representation.
//     T0          -   internal counter
//     T1          -   internal counter
//
// Outputs:
//     T0          -   new value of the internal counter
//     T1          -   new value of the internal counter
//     I           -   row index of non-zero element, 0 <= I < M.
//     J           -   column index of non-zero element, 0 <= J < N
//     V           -   value of the T-th element
//
// Result:
//     True in case of success (next non-zero element was retrieved)
//     False in case all non-zero elements were enumerated
//
// NOTE: you may call SparseRewriteExisting() during enumeration, but it  is
//       THE  ONLY  matrix  modification  function  you  can  call!!!  Other
//       matrix modification functions should not be called during enumeration!
//
// ALGLIB Project: Copyright 14.03.2012 by Sergey Bochkanov
// API: bool sparseenumerate(const sparsematrix &s, ae_int_t &t0, ae_int_t &t1, ae_int_t &i, ae_int_t &j, double &v);
bool sparseenumerate(sparsematrix *s, ae_int_t *t0, ae_int_t *t1, ae_int_t *i, ae_int_t *j, double *v) {
   ae_int_t sz;
   ae_int_t i0;
   bool result;
   *i = 0;
   *j = 0;
   *v = 0;
   result = false;
   if (*t0 < 0 || s->matrixtype != 0 && *t1 < 0) {
   // Incorrect T0/T1, terminate enumeration
      result = false;
      return result;
   }
   if (s->matrixtype == 0) {
   // Hash-table matrix
      sz = s->tablesize;
      for (i0 = *t0; i0 < sz; i0++) {
         if (s->idx.xZ[2 * i0] == -1 || s->idx.xZ[2 * i0] == -2) {
            continue;
         } else {
            *i = s->idx.xZ[2 * i0];
            *j = s->idx.xZ[2 * i0 + 1];
            *v = s->vals.xR[i0];
            *t0 = i0 + 1;
            result = true;
            return result;
         }
      }
      *t0 = 0;
      *t1 = 0;
      result = false;
      return result;
   }
   if (s->matrixtype == 1) {
   // CRS matrix
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseEnumerate: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      if (*t0 >= s->ninitialized) {
         *t0 = 0;
         *t1 = 0;
         result = false;
         return result;
      }
      while (*t0 > s->ridx.xZ[*t1 + 1] - 1 && *t1 < s->m) {
         ++*t1;
      }
      *i = *t1;
      *j = s->idx.xZ[*t0];
      *v = s->vals.xR[*t0];
      ++*t0;
      result = true;
      return result;
   }
   if (s->matrixtype == 2) {
   // SKS matrix:
   // * T0 stores current offset in Vals[] array
   // * T1 stores index of the diagonal block
      ae_assert(s->m == s->n, "SparseEnumerate: non-square SKS matrices are not supported");
      if (*t0 >= s->ridx.xZ[s->m]) {
         *t0 = 0;
         *t1 = 0;
         result = false;
         return result;
      }
      while (*t0 > s->ridx.xZ[*t1 + 1] - 1 && *t1 < s->m) {
         ++*t1;
      }
      i0 = *t0 - s->ridx.xZ[*t1];
      if (i0 < s->didx.xZ[*t1] + 1) {
      // subdiagonal or diagonal element, row index is T1.
         *i = *t1;
         *j = *t1 - s->didx.xZ[*t1] + i0;
      } else {
      // superdiagonal element, column index is T1.
         *i = *t1 - (s->ridx.xZ[*t1 + 1] - (*t0));
         *j = *t1;
      }
      *v = s->vals.xR[*t0];
      ++*t0;
      result = true;
      return result;
   }
   ae_assert(false, "SparseEnumerate: unexpected matrix type");
   return result;
}

// This function rewrites existing (non-zero) element. It  returns  True   if
// element  exists  or  False,  when  it  is  called for non-existing  (zero)
// element.
//
// This function works with any kind of the matrix.
//
// The purpose of this function is to provide convenient thread-safe  way  to
// modify  sparse  matrix.  Such  modification  (already  existing element is
// rewritten) is guaranteed to be thread-safe without any synchronization, as
// long as different threads modify different elements.
//
// Inputs:
//     S           -   sparse M*N matrix in any kind of representation
//                     (Hash, SKS, CRS).
//     I           -   row index of non-zero element to modify, 0 <= I < M
//     J           -   column index of non-zero element to modify, 0 <= J < N
//     V           -   value to rewrite, must be finite number
//
// Outputs:
//     S           -   modified matrix
// Result:
//     True in case when element exists
//     False in case when element doesn't exist or it is zero
//
// ALGLIB Project: Copyright 14.03.2012 by Sergey Bochkanov
// API: bool sparserewriteexisting(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v);
bool sparserewriteexisting(sparsematrix *s, ae_int_t i, ae_int_t j, double v) {
   ae_int_t hashcode;
   ae_int_t k;
   ae_int_t k0;
   ae_int_t k1;
   bool result;
   ae_assert(0 <= i && i < s->m, "SparseRewriteExisting: invalid argument I(either I<0 or I >= S.M)");
   ae_assert(0 <= j && j < s->n, "SparseRewriteExisting: invalid argument J(either J<0 or J >= S.N)");
   ae_assert(isfinite(v), "SparseRewriteExisting: invalid argument V(either V is infinite or V is NaN)");
   result = false;
// Hash-table matrix
   if (s->matrixtype == 0) {
      k = s->tablesize;
      hashcode = sparse_hash(i, j, k);
      while (true) {
         if (s->idx.xZ[2 * hashcode] == -1) {
            return result;
         }
         if (s->idx.xZ[2 * hashcode] == i && s->idx.xZ[2 * hashcode + 1] == j) {
            s->vals.xR[hashcode] = v;
            result = true;
            return result;
         }
         hashcode = (hashcode + 1) % k;
      }
   }
// CRS matrix
   if (s->matrixtype == 1) {
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseRewriteExisting: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      k0 = s->ridx.xZ[i];
      k1 = s->ridx.xZ[i + 1] - 1;
      while (k0 <= k1) {
         k = (k0 + k1) / 2;
         if (s->idx.xZ[k] == j) {
            s->vals.xR[k] = v;
            result = true;
            return result;
         }
         if (s->idx.xZ[k] < j) {
            k0 = k + 1;
         } else {
            k1 = k - 1;
         }
      }
   }
// SKS
   if (s->matrixtype == 2) {
      ae_assert(s->m == s->n, "SparseRewriteExisting: non-square SKS matrix not supported");
      if (i == j) {
      // Rewrite diagonal element
         result = true;
         s->vals.xR[s->ridx.xZ[i] + s->didx.xZ[i]] = v;
         return result;
      }
      if (j < i) {
      // Return subdiagonal element at I-th "skyline block"
         k = s->didx.xZ[i];
         if (i - j <= k) {
            s->vals.xR[s->ridx.xZ[i] + k + j - i] = v;
            result = true;
         }
      } else {
      // Return superdiagonal element at J-th "skyline block"
         k = s->uidx.xZ[j];
         if (j - i <= k) {
            s->vals.xR[s->ridx.xZ[j + 1] - (j - i)] = v;
            result = true;
         }
      }
      return result;
   }
   return result;
}

// This function returns I-th row of the sparse matrix. Matrix must be stored
// in CRS or SKS format.
//
// Inputs:
//     S           -   sparse M*N matrix in CRS format
//     I           -   row index, 0 <= I < M
//     IRow        -   output buffer, can be  preallocated.  In  case  buffer
//                     size  is  too  small  to  store  I-th   row,   it   is
//                     automatically reallocated.
//
// Outputs:
//     IRow        -   array[M], I-th row.
//
// NOTE: this function has O(N) running time, where N is a  column  count. It
//       allocates and fills N-element  array,  even  although  most  of  its
//       elemets are zero.
//
// NOTE: If you have O(non-zeros-per-row) time and memory  requirements,  use
//       SparseGetCompressedRow() function. It  returns  data  in  compressed
//       format.
//
// NOTE: when  incorrect  I  (outside  of  [0,M-1]) or  matrix (non  CRS/SKS)
//       is passed, this function throws exception.
//
// ALGLIB Project: Copyright 10.12.2014 by Sergey Bochkanov
// API: void sparsegetrow(const sparsematrix &s, const ae_int_t i, real_1d_array &irow);
void sparsegetrow(sparsematrix *s, ae_int_t i, RVector *irow) {
   ae_int_t i0;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t j;
   ae_int_t upperprofile;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseGetRow: S must be CRS/SKS-based matrix");
   ae_assert(i >= 0 && i < s->m, "SparseGetRow: I<0 or I >= M");
// Prepare output buffer
   vectorsetlengthatleast(irow, s->n);
   for (i0 = 0; i0 < s->n; i0++) {
      irow->xR[i0] = 0.0;
   }
// Output
   if (s->matrixtype == 1) {
      for (i0 = s->ridx.xZ[i]; i0 < s->ridx.xZ[i + 1]; i0++) {
         irow->xR[s->idx.xZ[i0]] = s->vals.xR[i0];
      }
      return;
   }
   if (s->matrixtype == 2) {
   // Copy subdiagonal and diagonal parts
      ae_assert(s->n == s->m, "SparseGetRow: non-square SKS matrices are not supported");
      j0 = i - s->didx.xZ[i];
      i0 = -j0 + s->ridx.xZ[i];
      for (j = j0; j <= i; j++) {
         irow->xR[j] = s->vals.xR[j + i0];
      }
   // Copy superdiagonal part
      upperprofile = s->uidx.xZ[s->n];
      j0 = i + 1;
      j1 = imin2(s->n - 1, i + upperprofile);
      for (j = j0; j <= j1; j++) {
         if (j - i <= s->uidx.xZ[j]) {
            irow->xR[j] = s->vals.xR[s->ridx.xZ[j + 1] - (j - i)];
         }
      }
      return;
   }
}

// This function returns I-th row of the sparse matrix IN COMPRESSED FORMAT -
// only non-zero elements are returned (with their indexes). Matrix  must  be
// stored in CRS or SKS format.
//
// Inputs:
//     S           -   sparse M*N matrix in CRS format
//     I           -   row index, 0 <= I < M
//     ColIdx      -   output buffer for column indexes, can be preallocated.
//                     In case buffer size is too small to store I-th row, it
//                     is automatically reallocated.
//     Vals        -   output buffer for values, can be preallocated. In case
//                     buffer size is too small to  store  I-th  row,  it  is
//                     automatically reallocated.
//
// Outputs:
//     ColIdx      -   column   indexes   of  non-zero  elements,  sorted  by
//                     ascending. Symbolically non-zero elements are  counted
//                     (i.e. if you allocated place for element, but  it  has
//                     zero numerical value - it is counted).
//     Vals        -   values. Vals[K] stores value of  matrix  element  with
//                     indexes (I,ColIdx[K]). Symbolically non-zero  elements
//                     are counted (i.e. if you allocated place for  element,
//                     but it has zero numerical value - it is counted).
//     NZCnt       -   number of symbolically non-zero elements per row.
//
// NOTE: when  incorrect  I  (outside  of  [0,M-1]) or  matrix (non  CRS/SKS)
//       is passed, this function throws exception.
//
// NOTE: this function may allocate additional, unnecessary place for  ColIdx
//       and Vals arrays. It is dictated by  performance  reasons  -  on  SKS
//       matrices it is faster  to  allocate  space  at  the  beginning  with
//       some "extra"-space, than performing two passes over matrix  -  first
//       time to calculate exact space required for data, second  time  -  to
//       store data itself.
//
// ALGLIB Project: Copyright 10.12.2014 by Sergey Bochkanov
// API: void sparsegetcompressedrow(const sparsematrix &s, const ae_int_t i, integer_1d_array &colidx, real_1d_array &vals, ae_int_t &nzcnt);
void sparsegetcompressedrow(sparsematrix *s, ae_int_t i, ZVector *colidx, RVector *vals, ae_int_t *nzcnt) {
   ae_int_t k;
   ae_int_t k0;
   ae_int_t j;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t i0;
   ae_int_t upperprofile;
   *nzcnt = 0;
   ae_assert(s->matrixtype == 1 || s->matrixtype == 2, "SparseGetRow: S must be CRS/SKS-based matrix");
   ae_assert(i >= 0 && i < s->m, "SparseGetRow: I<0 or I >= M");
// Initialize NZCnt
   *nzcnt = 0;
// CRS matrix - just copy data
   if (s->matrixtype == 1) {
      *nzcnt = s->ridx.xZ[i + 1] - s->ridx.xZ[i];
      vectorsetlengthatleast(colidx, *nzcnt);
      vectorsetlengthatleast(vals, *nzcnt);
      k0 = s->ridx.xZ[i];
      for (k = 0; k < *nzcnt; k++) {
         colidx->xZ[k] = s->idx.xZ[k0 + k];
         vals->xR[k] = s->vals.xR[k0 + k];
      }
      return;
   }
// SKS matrix - a bit more complex sequence
   if (s->matrixtype == 2) {
      ae_assert(s->n == s->m, "SparseGetCompressedRow: non-square SKS matrices are not supported");
   // Allocate enough place for storage
      upperprofile = s->uidx.xZ[s->n];
      vectorsetlengthatleast(colidx, s->didx.xZ[i] + 1 + upperprofile);
      vectorsetlengthatleast(vals, s->didx.xZ[i] + 1 + upperprofile);
   // Copy subdiagonal and diagonal parts
      j0 = i - s->didx.xZ[i];
      i0 = -j0 + s->ridx.xZ[i];
      for (j = j0; j <= i; j++) {
         colidx->xZ[*nzcnt] = j;
         vals->xR[*nzcnt] = s->vals.xR[j + i0];
         ++*nzcnt;
      }
   // Copy superdiagonal part
      j0 = i + 1;
      j1 = imin2(s->n - 1, i + upperprofile);
      for (j = j0; j <= j1; j++) {
         if (j - i <= s->uidx.xZ[j]) {
            colidx->xZ[*nzcnt] = j;
            vals->xR[*nzcnt] = s->vals.xR[s->ridx.xZ[j + 1] - (j - i)];
            ++*nzcnt;
         }
      }
      return;
   }
}

// This function performs efficient in-place  transpose  of  SKS  matrix.  No
// additional memory is allocated during transposition.
//
// This function supports only skyline storage format (SKS).
//
// Inputs:
//     S       -   sparse matrix in SKS format.
//
// Outputs:
//     S           -   sparse matrix, transposed.
//
// ALGLIB Project: Copyright 16.01.2014 by Sergey Bochkanov
// API: void sparsetransposesks(const sparsematrix &s);
void sparsetransposesks(sparsematrix *s) {
   ae_int_t n;
   ae_int_t d;
   ae_int_t u;
   ae_int_t i;
   ae_int_t k;
   ae_int_t t0;
   ae_int_t t1;
   ae_assert(s->matrixtype == 2, "SparseTransposeSKS: only SKS matrices are supported");
   ae_assert(s->m == s->n, "SparseTransposeSKS: non-square SKS matrices are not supported");
   n = s->n;
   for (i = 1; i < n; i++) {
      d = s->didx.xZ[i];
      u = s->uidx.xZ[i];
      swapi(&s->uidx.xZ[i], &s->didx.xZ[i]);
      if (d == u) {
      // Upper skyline height equal to lower skyline height,
      // simple exchange is needed for transposition
         t0 = s->ridx.xZ[i];
         for (k = 0; k < d; k++) {
            swapr(&s->vals.xR[t0 + k], &s->vals.xR[t0 + d + 1 + k]);
         }
      }
      if (d > u) {
      // Upper skyline height is less than lower skyline height.
      //
      // Transposition becomes a bit tricky: we have to rearrange
      // "L0 L1 D U" to "U D L0 L1", where |L0|=|U|=u, |L1|=d-u.
      //
      // In order to do this we perform a sequence of swaps and
      // in-place reversals:
      // * swap(L0,U)         =>  "U   L1  D   L0"
      // * reverse("L1 D L0") =>  "U   L0~ D   L1~" (where X~ is a reverse of X)
      // * reverse("L0~ D")   =>  "U   D   L0  L1~"
      // * reverse("L1")      =>  "U   D   L0  L1"
         t0 = s->ridx.xZ[i];
         t1 = s->ridx.xZ[i] + d + 1;
         for (k = 0; k < u; k++) {
            swapr(&s->vals.xR[t0 + k], &s->vals.xR[t1 + k]);
         }
         t0 = s->ridx.xZ[i] + u;
         t1 = s->ridx.xZ[i + 1] - 1;
         while (t1 > t0) {
            swapr(&s->vals.xR[t0], &s->vals.xR[t1]);
            t0++;
            t1--;
         }
         t0 = s->ridx.xZ[i] + u;
         t1 = s->ridx.xZ[i] + u + u;
         while (t1 > t0) {
            swapr(&s->vals.xR[t0], &s->vals.xR[t1]);
            t0++;
            t1--;
         }
         t0 = s->ridx.xZ[i + 1] - (d - u);
         t1 = s->ridx.xZ[i + 1] - 1;
         while (t1 > t0) {
            swapr(&s->vals.xR[t0], &s->vals.xR[t1]);
            t0++;
            t1--;
         }
      }
      if (d < u) {
      // Upper skyline height is greater than lower skyline height.
      //
      // Transposition becomes a bit tricky: we have to rearrange
      // "L D U0 U1" to "U0 U1 D L", where |U1|=|L|=d, |U0|=u-d.
      //
      // In order to do this we perform a sequence of swaps and
      // in-place reversals:
      // * swap(L,U1)         =>  "U1  D   U0  L"
      // * reverse("U1 D U0") =>  "U0~ D   U1~ L" (where X~ is a reverse of X)
      // * reverse("U0~")     =>  "U0  D   U1~ L"
      // * reverse("D U1~")   =>  "U0  U1  D   L"
         t0 = s->ridx.xZ[i];
         t1 = s->ridx.xZ[i + 1] - d;
         for (k = 0; k < d; k++) {
            swapr(&s->vals.xR[t0 + k], &s->vals.xR[t1 + k]);
         }
         t0 = s->ridx.xZ[i];
         t1 = s->ridx.xZ[i] + u;
         while (t1 > t0) {
            swapr(&s->vals.xR[t0], &s->vals.xR[t1]);
            t0++;
            t1--;
         }
         t0 = s->ridx.xZ[i];
         t1 = s->ridx.xZ[i] + u - d - 1;
         while (t1 > t0) {
            swapr(&s->vals.xR[t0], &s->vals.xR[t1]);
            t0++;
            t1--;
         }
         t0 = s->ridx.xZ[i] + u - d;
         t1 = s->ridx.xZ[i + 1] - d - 1;
         while (t1 > t0) {
            swapr(&s->vals.xR[t0], &s->vals.xR[t1]);
            t0++;
            t1--;
         }
      }
   }
   swapi(&s->uidx.xZ[n], &s->didx.xZ[n]);
}

// This function performs transpose of CRS matrix.
//
// Inputs:
//     S       -   sparse matrix in CRS format.
//
// Outputs:
//     S           -   sparse matrix, transposed.
//
// NOTE: internal  temporary  copy  is  allocated   for   the   purposes   of
//       transposition. It is deallocated after transposition.
//
// ALGLIB Project: Copyright 30.01.2018 by Sergey Bochkanov
// API: void sparsetransposecrs(const sparsematrix &s);
void sparsetransposecrs(sparsematrix *s) {
   ae_frame _frame_block;
   ae_int_t oldn;
   ae_int_t oldm;
   ae_int_t newn;
   ae_int_t newm;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t nonne;
   ae_frame_make(&_frame_block);
   NewVector(oldvals, 0, DT_REAL);
   NewVector(oldidx, 0, DT_INT);
   NewVector(oldridx, 0, DT_INT);
   NewVector(counts, 0, DT_INT);
   ae_assert(s->matrixtype == 1, "SparseTransposeCRS: only CRS matrices are supported");
   ae_swap_vectors(&s->vals, &oldvals);
   ae_swap_vectors(&s->idx, &oldidx);
   ae_swap_vectors(&s->ridx, &oldridx);
   oldn = s->n;
   oldm = s->m;
   newn = oldm;
   newm = oldn;
// Update matrix size
   s->n = newn;
   s->m = newm;
// Fill RIdx by number of elements per row:
// RIdx[I+1] stores number of elements in I-th row.
//
// Convert RIdx from row sizes to row offsets.
// Set NInitialized
   nonne = 0;
   vectorsetlengthatleast(&s->ridx, newm + 1);
   for (i = 0; i <= newm; i++) {
      s->ridx.xZ[i] = 0;
   }
   for (i = 0; i < oldm; i++) {
      for (j = oldridx.xZ[i]; j < oldridx.xZ[i + 1]; j++) {
         k = oldidx.xZ[j] + 1;
         s->ridx.xZ[k]++;
         nonne++;
      }
   }
   for (i = 0; i < newm; i++) {
      s->ridx.xZ[i + 1] += s->ridx.xZ[i];
   }
   s->ninitialized = s->ridx.xZ[newm];
// Allocate memory and move elements to Vals/Idx.
   ae_vector_set_length(&counts, newm);
   for (i = 0; i < newm; i++) {
      counts.xZ[i] = 0;
   }
   vectorsetlengthatleast(&s->vals, nonne);
   vectorsetlengthatleast(&s->idx, nonne);
   for (i = 0; i < oldm; i++) {
      for (j = oldridx.xZ[i]; j < oldridx.xZ[i + 1]; j++) {
         k = oldidx.xZ[j];
         k = s->ridx.xZ[k] + counts.xZ[k];
         s->idx.xZ[k] = i;
         s->vals.xR[k] = oldvals.xR[j];
         k = oldidx.xZ[j];
         counts.xZ[k]++;
      }
   }
// Initialization 'S.UIdx' and 'S.DIdx'
   sparseinitduidx(s);
   ae_frame_leave();
}

// This function performs copying with transposition of CRS matrix  (buffered
// version which reuses memory already allocated by  the  target as  much  as
// possible).
//
// Inputs:
//     S0      -   sparse matrix in CRS format.
//
// Outputs:
//     S1      -   sparse matrix, transposed; previously allocated memory  is
//                 reused if possible.
//
// ALGLIB Project: Copyright 23.07.2018 by Sergey Bochkanov
// API: void sparsecopytransposecrsbuf(const sparsematrix &s0, const sparsematrix &s1);
void sparsecopytransposecrsbuf(sparsematrix *s0, sparsematrix *s1) {
   ae_frame _frame_block;
   ae_int_t oldn;
   ae_int_t oldm;
   ae_int_t newn;
   ae_int_t newm;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t nonne;
   ae_frame_make(&_frame_block);
   NewVector(counts, 0, DT_INT);
   ae_assert(s0->matrixtype == 1, "SparseCopyTransposeCRSBuf: only CRS matrices are supported");
   oldn = s0->n;
   oldm = s0->m;
   newn = oldm;
   newm = oldn;
// Update matrix size
   s1->matrixtype = 1;
   s1->n = newn;
   s1->m = newm;
// Fill RIdx by number of elements per row:
// RIdx[I+1] stores number of elements in I-th row.
//
// Convert RIdx from row sizes to row offsets.
// Set NInitialized
   nonne = 0;
   vectorsetlengthatleast(&s1->ridx, newm + 1);
   for (i = 0; i <= newm; i++) {
      s1->ridx.xZ[i] = 0;
   }
   for (i = 0; i < oldm; i++) {
      for (j = s0->ridx.xZ[i]; j < s0->ridx.xZ[i + 1]; j++) {
         k = s0->idx.xZ[j] + 1;
         s1->ridx.xZ[k]++;
         nonne++;
      }
   }
   for (i = 0; i < newm; i++) {
      s1->ridx.xZ[i + 1] += s1->ridx.xZ[i];
   }
   s1->ninitialized = s1->ridx.xZ[newm];
// Allocate memory and move elements to Vals/Idx.
   ae_vector_set_length(&counts, newm);
   for (i = 0; i < newm; i++) {
      counts.xZ[i] = 0;
   }
   vectorsetlengthatleast(&s1->vals, nonne);
   vectorsetlengthatleast(&s1->idx, nonne);
   for (i = 0; i < oldm; i++) {
      for (j = s0->ridx.xZ[i]; j < s0->ridx.xZ[i + 1]; j++) {
         k = s0->idx.xZ[j];
         k = s1->ridx.xZ[k] + counts.xZ[k];
         s1->idx.xZ[k] = i;
         s1->vals.xR[k] = s0->vals.xR[j];
         k = s0->idx.xZ[j];
         counts.xZ[k]++;
      }
   }
// Initialization 'S.UIdx' and 'S.DIdx'
   sparseinitduidx(s1);
   ae_frame_leave();
}

// This function performs copying with transposition of CRS matrix.
//
// Inputs:
//     S0      -   sparse matrix in CRS format.
//
// Outputs:
//     S1      -   sparse matrix, transposed
//
// ALGLIB Project: Copyright 23.07.2018 by Sergey Bochkanov
// API: void sparsecopytransposecrs(const sparsematrix &s0, sparsematrix &s1);
void sparsecopytransposecrs(sparsematrix *s0, sparsematrix *s1) {
   SetObj(sparsematrix, s1);
   sparsecopytransposecrsbuf(s0, s1);
}

// This  function  performs  in-place  conversion  to  desired sparse storage
// format.
//
// Inputs:
//     S0      -   sparse matrix in any format.
//     Fmt     -   desired storage format  of  the  output,  as  returned  by
//                 SparseGetMatrixType() function:
//                 * 0 for hash-based storage
//                 * 1 for CRS
//                 * 2 for SKS
//
// Outputs:
//     S0          -   sparse matrix in requested format.
//
// NOTE: in-place conversion wastes a lot of memory which is  used  to  store
//       temporaries.  If  you  perform  a  lot  of  repeated conversions, we
//       recommend to use out-of-place buffered  conversion  functions,  like
//       SparseCopyToBuf(), which can reuse already allocated memory.
//
// ALGLIB Project: Copyright 16.01.2014 by Sergey Bochkanov
// API: void sparseconvertto(const sparsematrix &s0, const ae_int_t fmt);
void sparseconvertto(sparsematrix *s0, ae_int_t fmt) {
   ae_assert(fmt == 0 || fmt == 1 || fmt == 2, "SparseConvertTo: invalid fmt parameter");
   if (fmt == 0) {
      sparseconverttohash(s0);
      return;
   }
   if (fmt == 1) {
      sparseconverttocrs(s0);
      return;
   }
   if (fmt == 2) {
      sparseconverttosks(s0);
      return;
   }
   ae_assert(false, "SparseConvertTo: invalid matrix type");
}

// This  function  performs out-of-place conversion to desired sparse storage
// format. S0 is copied to S1 and converted on-the-fly. Memory  allocated  in
// S1 is reused to maximum extent possible.
//
// Inputs:
//     S0      -   sparse matrix in any format.
//     Fmt     -   desired storage format  of  the  output,  as  returned  by
//                 SparseGetMatrixType() function:
//                 * 0 for hash-based storage
//                 * 1 for CRS
//                 * 2 for SKS
//
// Outputs:
//     S1          -   sparse matrix in requested format.
//
// ALGLIB Project: Copyright 16.01.2014 by Sergey Bochkanov
// API: void sparsecopytobuf(const sparsematrix &s0, const ae_int_t fmt, const sparsematrix &s1);
void sparsecopytobuf(sparsematrix *s0, ae_int_t fmt, sparsematrix *s1) {
   ae_assert(fmt == 0 || fmt == 1 || fmt == 2, "SparseCopyToBuf: invalid fmt parameter");
   if (fmt == 0) {
      sparsecopytohashbuf(s0, s1);
      return;
   }
   if (fmt == 1) {
      sparsecopytocrsbuf(s0, s1);
      return;
   }
   if (fmt == 2) {
      sparsecopytosksbuf(s0, s1);
      return;
   }
   ae_assert(false, "SparseCopyToBuf: invalid matrix type");
}

// This function performs in-place conversion to Hash table storage.
//
// Inputs:
//     S           -   sparse matrix in CRS format.
//
// Outputs:
//     S           -   sparse matrix in Hash table format.
//
// NOTE: this  function  has   no  effect  when  called with matrix which  is
//       already in Hash table mode.
//
// NOTE: in-place conversion involves allocation of temporary arrays. If  you
//       perform a lot of repeated in- place  conversions,  it  may  lead  to
//       memory fragmentation. Consider using out-of-place SparseCopyToHashBuf()
//       function in this case.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparseconverttohash(const sparsematrix &s);
void sparseconverttohash(sparsematrix *s) {
   ae_frame _frame_block;
   ae_int_t n;
   ae_int_t m;
   ae_int_t offs0;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_frame_make(&_frame_block);
   NewVector(tidx, 0, DT_INT);
   NewVector(tridx, 0, DT_INT);
   NewVector(tdidx, 0, DT_INT);
   NewVector(tuidx, 0, DT_INT);
   NewVector(tvals, 0, DT_REAL);
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseConvertToHash: invalid matrix type");
   if (s->matrixtype == 0) {
   // Already in Hash mode
      ae_frame_leave();
      return;
   }
   if (s->matrixtype == 1) {
   // From CRS to Hash
      s->matrixtype = 0;
      m = s->m;
      n = s->n;
      ae_swap_vectors(&s->idx, &tidx);
      ae_swap_vectors(&s->ridx, &tridx);
      ae_swap_vectors(&s->vals, &tvals);
      sparsecreatebuf(m, n, tridx.xZ[m], s);
      for (i = 0; i < m; i++) {
         for (j = tridx.xZ[i]; j < tridx.xZ[i + 1]; j++) {
            sparseset(s, i, tidx.xZ[j], tvals.xR[j]);
         }
      }
      ae_frame_leave();
      return;
   }
   if (s->matrixtype == 2) {
   // From SKS to Hash
      s->matrixtype = 0;
      m = s->m;
      n = s->n;
      ae_swap_vectors(&s->ridx, &tridx);
      ae_swap_vectors(&s->didx, &tdidx);
      ae_swap_vectors(&s->uidx, &tuidx);
      ae_swap_vectors(&s->vals, &tvals);
      sparsecreatebuf(m, n, tridx.xZ[m], s);
      for (i = 0; i < m; i++) {
      // copy subdiagonal and diagonal parts of I-th block
         offs0 = tridx.xZ[i];
         k = tdidx.xZ[i] + 1;
         for (j = 0; j < k; j++) {
            sparseset(s, i, i - tdidx.xZ[i] + j, tvals.xR[offs0 + j]);
         }
      // Copy superdiagonal part of I-th block
         offs0 = tridx.xZ[i] + tdidx.xZ[i] + 1;
         k = tuidx.xZ[i];
         for (j = 0; j < k; j++) {
            sparseset(s, i - k + j, i, tvals.xR[offs0 + j]);
         }
      }
      ae_frame_leave();
      return;
   }
   ae_assert(false, "SparseConvertToHash: invalid matrix type");
   ae_frame_leave();
}

// This  function  performs  out-of-place  conversion  to  Hash table storage
// format. S0 is copied to S1 and converted on-the-fly.
//
// Inputs:
//     S0          -   sparse matrix in any format.
//
// Outputs:
//     S1          -   sparse matrix in Hash table format.
//
// NOTE: if S0 is stored as Hash-table, it is just copied without conversion.
//
// NOTE: this function de-allocates memory  occupied  by  S1 before  starting
//       conversion. If you perform a  lot  of  repeated  conversions, it may
//       lead to memory fragmentation. In this case we recommend you  to  use
//       SparseCopyToHashBuf() function which re-uses memory in S1 as much as
//       possible.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparsecopytohash(const sparsematrix &s0, sparsematrix &s1);
void sparsecopytohash(sparsematrix *s0, sparsematrix *s1) {
   SetObj(sparsematrix, s1);
   ae_assert(s0->matrixtype == 0 || s0->matrixtype == 1 || s0->matrixtype == 2, "SparseCopyToHash: invalid matrix type");
   sparsecopytohashbuf(s0, s1);
}

// This  function  performs  out-of-place  conversion  to  Hash table storage
// format. S0 is copied to S1 and converted on-the-fly. Memory  allocated  in
// S1 is reused to maximum extent possible.
//
// Inputs:
//     S0          -   sparse matrix in any format.
//
// Outputs:
//     S1          -   sparse matrix in Hash table format.
//
// NOTE: if S0 is stored as Hash-table, it is just copied without conversion.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparsecopytohashbuf(const sparsematrix &s0, const sparsematrix &s1);
void sparsecopytohashbuf(sparsematrix *s0, sparsematrix *s1) {
   double val;
   ae_int_t t0;
   ae_int_t t1;
   ae_int_t i;
   ae_int_t j;
   ae_assert(s0->matrixtype == 0 || s0->matrixtype == 1 || s0->matrixtype == 2, "SparseCopyToHashBuf: invalid matrix type");
   if (s0->matrixtype == 0) {
   // Already hash, just copy
      sparsecopybuf(s0, s1);
      return;
   }
   if (s0->matrixtype == 1) {
   // CRS storage
      t0 = 0;
      t1 = 0;
      sparsecreatebuf(s0->m, s0->n, s0->ridx.xZ[s0->m], s1);
      while (sparseenumerate(s0, &t0, &t1, &i, &j, &val)) {
         sparseset(s1, i, j, val);
      }
      return;
   }
   if (s0->matrixtype == 2) {
   // SKS storage
      t0 = 0;
      t1 = 0;
      sparsecreatebuf(s0->m, s0->n, s0->ridx.xZ[s0->m], s1);
      while (sparseenumerate(s0, &t0, &t1, &i, &j, &val)) {
         sparseset(s1, i, j, val);
      }
      return;
   }
   ae_assert(false, "SparseCopyToHashBuf: invalid matrix type");
}

// This function converts matrix to CRS format.
//
// Some  algorithms  (linear  algebra ones, for example) require matrices in
// CRS format. This function allows to perform in-place conversion.
//
// Inputs:
//     S           -   sparse M*N matrix in any format
//
// Outputs:
//     S           -   matrix in CRS format
//
// NOTE: this   function  has  no  effect  when  called with matrix which is
//       already in CRS mode.
//
// NOTE: this function allocates temporary memory to store a   copy  of  the
//       matrix. If you perform a lot of repeated conversions, we  recommend
//       you  to  use  SparseCopyToCRSBuf()  function,   which   can   reuse
//       previously allocated memory.
//
// ALGLIB Project: Copyright 14.10.2011 by Sergey Bochkanov
// API: void sparseconverttocrs(const sparsematrix &s);
void sparseconverttocrs(sparsematrix *s) {
   ae_frame _frame_block;
   ae_int_t m;
   ae_int_t i;
   ae_int_t j;
   ae_int_t nonne;
   ae_int_t k;
   ae_int_t offs0;
   ae_int_t offs1;
   ae_frame_make(&_frame_block);
   NewVector(tvals, 0, DT_REAL);
   NewVector(tidx, 0, DT_INT);
   NewVector(temp, 0, DT_INT);
   NewVector(tridx, 0, DT_INT);
   m = s->m;
   if (s->matrixtype == 0) {
   // From Hash-table to CRS.
   // First, create local copy of the hash table.
      s->matrixtype = 1;
      k = s->tablesize;
      ae_swap_vectors(&s->vals, &tvals);
      ae_swap_vectors(&s->idx, &tidx);
   // Fill RIdx by number of elements per row:
   // RIdx[I+1] stores number of elements in I-th row.
   //
   // Convert RIdx from row sizes to row offsets.
   // Set NInitialized
      nonne = 0;
      vectorsetlengthatleast(&s->ridx, s->m + 1);
      for (i = 0; i <= s->m; i++) {
         s->ridx.xZ[i] = 0;
      }
      for (i = 0; i < k; i++) {
         if (tidx.xZ[2 * i] >= 0) {
            s->ridx.xZ[tidx.xZ[2 * i] + 1]++;
            nonne++;
         }
      }
      for (i = 0; i < s->m; i++) {
         s->ridx.xZ[i + 1] += s->ridx.xZ[i];
      }
      s->ninitialized = s->ridx.xZ[s->m];
   // Allocate memory and move elements to Vals/Idx.
   // Initially, elements are sorted by rows, but unsorted within row.
   // After initial insertion we sort elements within row.
      ae_vector_set_length(&temp, s->m);
      for (i = 0; i < s->m; i++) {
         temp.xZ[i] = 0;
      }
      vectorsetlengthatleast(&s->vals, nonne);
      vectorsetlengthatleast(&s->idx, nonne);
      for (i = 0; i < k; i++) {
         if (tidx.xZ[2 * i] >= 0) {
            s->vals.xR[s->ridx.xZ[tidx.xZ[2 * i]] + temp.xZ[tidx.xZ[2 * i]]] = tvals.xR[i];
            s->idx.xZ[s->ridx.xZ[tidx.xZ[2 * i]] + temp.xZ[tidx.xZ[2 * i]]] = tidx.xZ[2 * i + 1];
            temp.xZ[tidx.xZ[2 * i]]++;
         }
      }
      for (i = 0; i < s->m; i++) {
         tagsortmiddleir(&s->idx, &s->vals, s->ridx.xZ[i + 1] - s->ridx.xZ[i], s->ridx.xZ[i]);
      }
   // Initialization 'S.UIdx' and 'S.DIdx'
      sparseinitduidx(s);
      ae_frame_leave();
      return;
   }
   if (s->matrixtype == 1) {
   // Already CRS
      ae_frame_leave();
      return;
   }
   if (s->matrixtype == 2) {
      ae_assert(s->m == s->n, "SparseConvertToCRS: non-square SKS matrices are not supported");
   // From SKS to CRS.
   //
   // First, create local copy of the SKS matrix (Vals,
   // Idx, RIdx are stored; DIdx/UIdx for some time are
   // left in the SparseMatrix structure).
      s->matrixtype = 1;
      ae_swap_vectors(&s->vals, &tvals);
      ae_swap_vectors(&s->idx, &tidx);
      ae_swap_vectors(&s->ridx, &tridx);
   // Fill RIdx by number of elements per row:
   // RIdx[I+1] stores number of elements in I-th row.
   //
   // Convert RIdx from row sizes to row offsets.
   // Set NInitialized
      vectorsetlengthatleast(&s->ridx, m + 1);
      s->ridx.xZ[0] = 0;
      for (i = 1; i <= m; i++) {
         s->ridx.xZ[i] = 1;
      }
      nonne = 0;
      for (i = 0; i < m; i++) {
         s->ridx.xZ[i + 1] += s->didx.xZ[i];
         for (j = i - s->uidx.xZ[i]; j < i; j++) {
            s->ridx.xZ[j + 1]++;
         }
         nonne += s->didx.xZ[i] + 1 + s->uidx.xZ[i];
      }
      for (i = 0; i < s->m; i++) {
         s->ridx.xZ[i + 1] += s->ridx.xZ[i];
      }
      s->ninitialized = s->ridx.xZ[s->m];
   // Allocate memory and move elements to Vals/Idx.
   // Initially, elements are sorted by rows, and are sorted within row too.
   // No additional post-sorting is required.
      ae_vector_set_length(&temp, m);
      for (i = 0; i < m; i++) {
         temp.xZ[i] = 0;
      }
      vectorsetlengthatleast(&s->vals, nonne);
      vectorsetlengthatleast(&s->idx, nonne);
      for (i = 0; i < m; i++) {
      // copy subdiagonal and diagonal parts of I-th block
         offs0 = tridx.xZ[i];
         offs1 = s->ridx.xZ[i] + temp.xZ[i];
         k = s->didx.xZ[i] + 1;
         for (j = 0; j < k; j++) {
            s->vals.xR[offs1 + j] = tvals.xR[offs0 + j];
            s->idx.xZ[offs1 + j] = i - s->didx.xZ[i] + j;
         }
         temp.xZ[i] += s->didx.xZ[i] + 1;
      // Copy superdiagonal part of I-th block
         offs0 = tridx.xZ[i] + s->didx.xZ[i] + 1;
         k = s->uidx.xZ[i];
         for (j = 0; j < k; j++) {
            offs1 = s->ridx.xZ[i - k + j] + temp.xZ[i - k + j];
            s->vals.xR[offs1] = tvals.xR[offs0 + j];
            s->idx.xZ[offs1] = i;
            temp.xZ[i - k + j]++;
         }
      }
   // Initialization 'S.UIdx' and 'S.DIdx'
      sparseinitduidx(s);
      ae_frame_leave();
      return;
   }
   ae_assert(false, "SparseConvertToCRS: invalid matrix type");
   ae_frame_leave();
}

// This  function  performs  out-of-place  conversion  to  CRS format.  S0 is
// copied to S1 and converted on-the-fly.
//
// Inputs:
//     S0          -   sparse matrix in any format.
//
// Outputs:
//     S1          -   sparse matrix in CRS format.
//
// NOTE: if S0 is stored as CRS, it is just copied without conversion.
//
// NOTE: this function de-allocates memory occupied by S1 before starting CRS
//       conversion. If you perform a lot of repeated CRS conversions, it may
//       lead to memory fragmentation. In this case we recommend you  to  use
//       SparseCopyToCRSBuf() function which re-uses memory in S1 as much  as
//       possible.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparsecopytocrs(const sparsematrix &s0, sparsematrix &s1);
void sparsecopytocrs(sparsematrix *s0, sparsematrix *s1) {
   SetObj(sparsematrix, s1);
   ae_assert(s0->matrixtype == 0 || s0->matrixtype == 1 || s0->matrixtype == 2, "SparseCopyToCRS: invalid matrix type");
   sparsecopytocrsbuf(s0, s1);
}

// This  function  performs  out-of-place  conversion  to  CRS format.  S0 is
// copied to S1 and converted on-the-fly. Memory allocated in S1 is reused to
// maximum extent possible.
//
// Inputs:
//     S0          -   sparse matrix in any format.
//     S1          -   matrix which may contain some pre-allocated memory, or
//                     can be just uninitialized structure.
//
// Outputs:
//     S1          -   sparse matrix in CRS format.
//
// NOTE: if S0 is stored as CRS, it is just copied without conversion.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparsecopytocrsbuf(const sparsematrix &s0, const sparsematrix &s1);
void sparsecopytocrsbuf(sparsematrix *s0, sparsematrix *s1) {
   ae_frame _frame_block;
   ae_int_t nonne;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t offs0;
   ae_int_t offs1;
   ae_int_t m;
   ae_frame_make(&_frame_block);
   NewVector(temp, 0, DT_INT);
   ae_assert(s0->matrixtype == 0 || s0->matrixtype == 1 || s0->matrixtype == 2, "SparseCopyToCRSBuf: invalid matrix type");
   m = s0->m;
   if (s0->matrixtype == 0) {
   // Convert from hash-table to CRS
   // Done like ConvertToCRS function
      s1->matrixtype = 1;
      s1->m = s0->m;
      s1->n = s0->n;
      s1->nfree = s0->nfree;
      nonne = 0;
      k = s0->tablesize;
      vectorsetlengthatleast(&s1->ridx, s1->m + 1);
      for (i = 0; i <= s1->m; i++) {
         s1->ridx.xZ[i] = 0;
      }
      ae_vector_set_length(&temp, s1->m);
      for (i = 0; i < s1->m; i++) {
         temp.xZ[i] = 0;
      }
   // Number of elements per row
      for (i = 0; i < k; i++) {
         if (s0->idx.xZ[2 * i] >= 0) {
            s1->ridx.xZ[s0->idx.xZ[2 * i] + 1]++;
            nonne++;
         }
      }
   // Fill RIdx (offsets of rows)
      for (i = 0; i < s1->m; i++) {
         s1->ridx.xZ[i + 1] += s1->ridx.xZ[i];
      }
   // Allocate memory
      vectorsetlengthatleast(&s1->vals, nonne);
      vectorsetlengthatleast(&s1->idx, nonne);
      for (i = 0; i < k; i++) {
         if (s0->idx.xZ[2 * i] >= 0) {
            s1->vals.xR[s1->ridx.xZ[s0->idx.xZ[2 * i]] + temp.xZ[s0->idx.xZ[2 * i]]] = s0->vals.xR[i];
            s1->idx.xZ[s1->ridx.xZ[s0->idx.xZ[2 * i]] + temp.xZ[s0->idx.xZ[2 * i]]] = s0->idx.xZ[2 * i + 1];
            temp.xZ[s0->idx.xZ[2 * i]]++;
         }
      }
   // Set NInitialized
      s1->ninitialized = s1->ridx.xZ[s1->m];
   // Sorting of elements
      for (i = 0; i < s1->m; i++) {
         tagsortmiddleir(&s1->idx, &s1->vals, s1->ridx.xZ[i + 1] - s1->ridx.xZ[i], s1->ridx.xZ[i]);
      }
   // Initialization 'S.UIdx' and 'S.DIdx'
      sparseinitduidx(s1);
      ae_frame_leave();
      return;
   }
   if (s0->matrixtype == 1) {
   // Already CRS, just copy
      sparsecopybuf(s0, s1);
      ae_frame_leave();
      return;
   }
   if (s0->matrixtype == 2) {
      ae_assert(s0->m == s0->n, "SparseCopyToCRS: non-square SKS matrices are not supported");
   // From SKS to CRS.
      s1->m = s0->m;
      s1->n = s0->n;
      s1->matrixtype = 1;
   // Fill RIdx by number of elements per row:
   // RIdx[I+1] stores number of elements in I-th row.
   //
   // Convert RIdx from row sizes to row offsets.
   // Set NInitialized
      vectorsetlengthatleast(&s1->ridx, m + 1);
      s1->ridx.xZ[0] = 0;
      for (i = 1; i <= m; i++) {
         s1->ridx.xZ[i] = 1;
      }
      nonne = 0;
      for (i = 0; i < m; i++) {
         s1->ridx.xZ[i + 1] += s0->didx.xZ[i];
         for (j = i - s0->uidx.xZ[i]; j < i; j++) {
            s1->ridx.xZ[j + 1]++;
         }
         nonne += s0->didx.xZ[i] + 1 + s0->uidx.xZ[i];
      }
      for (i = 0; i < m; i++) {
         s1->ridx.xZ[i + 1] += s1->ridx.xZ[i];
      }
      s1->ninitialized = s1->ridx.xZ[m];
   // Allocate memory and move elements to Vals/Idx.
   // Initially, elements are sorted by rows, and are sorted within row too.
   // No additional post-sorting is required.
      ae_vector_set_length(&temp, m);
      for (i = 0; i < m; i++) {
         temp.xZ[i] = 0;
      }
      vectorsetlengthatleast(&s1->vals, nonne);
      vectorsetlengthatleast(&s1->idx, nonne);
      for (i = 0; i < m; i++) {
      // copy subdiagonal and diagonal parts of I-th block
         offs0 = s0->ridx.xZ[i];
         offs1 = s1->ridx.xZ[i] + temp.xZ[i];
         k = s0->didx.xZ[i] + 1;
         for (j = 0; j < k; j++) {
            s1->vals.xR[offs1 + j] = s0->vals.xR[offs0 + j];
            s1->idx.xZ[offs1 + j] = i - s0->didx.xZ[i] + j;
         }
         temp.xZ[i] += s0->didx.xZ[i] + 1;
      // Copy superdiagonal part of I-th block
         offs0 = s0->ridx.xZ[i] + s0->didx.xZ[i] + 1;
         k = s0->uidx.xZ[i];
         for (j = 0; j < k; j++) {
            offs1 = s1->ridx.xZ[i - k + j] + temp.xZ[i - k + j];
            s1->vals.xR[offs1] = s0->vals.xR[offs0 + j];
            s1->idx.xZ[offs1] = i;
            temp.xZ[i - k + j]++;
         }
      }
   // Initialization 'S.UIdx' and 'S.DIdx'
      sparseinitduidx(s1);
      ae_frame_leave();
      return;
   }
   ae_assert(false, "SparseCopyToCRSBuf: unexpected matrix type");
   ae_frame_leave();
}

// This function performs in-place conversion to SKS format.
//
// Inputs:
//     S           -   sparse matrix in any format.
//
// Outputs:
//     S           -   sparse matrix in SKS format.
//
// NOTE: this  function  has   no  effect  when  called with matrix which  is
//       already in SKS mode.
//
// NOTE: in-place conversion involves allocation of temporary arrays. If  you
//       perform a lot of repeated in- place  conversions,  it  may  lead  to
//       memory fragmentation. Consider using out-of-place SparseCopyToSKSBuf()
//       function in this case.
//
// ALGLIB Project: Copyright 15.01.2014 by Sergey Bochkanov
// API: void sparseconverttosks(const sparsematrix &s);
void sparseconverttosks(sparsematrix *s) {
   ae_frame _frame_block;
   ae_int_t n;
   ae_int_t t0;
   ae_int_t t1;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double v;
   ae_frame_make(&_frame_block);
   NewVector(tridx, 0, DT_INT);
   NewVector(tdidx, 0, DT_INT);
   NewVector(tuidx, 0, DT_INT);
   NewVector(tvals, 0, DT_REAL);
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseConvertToSKS: invalid matrix type");
   ae_assert(s->m == s->n, "SparseConvertToSKS: rectangular matrices are not supported");
   n = s->n;
   if (s->matrixtype == 2) {
   // Already in SKS mode
      ae_frame_leave();
      return;
   }
// Generate internal copy of SKS matrix
   vectorsetlengthatleast(&tdidx, n + 1);
   vectorsetlengthatleast(&tuidx, n + 1);
   for (i = 0; i <= n; i++) {
      tdidx.xZ[i] = 0;
      tuidx.xZ[i] = 0;
   }
   t0 = 0;
   t1 = 0;
   while (sparseenumerate(s, &t0, &t1, &i, &j, &v)) {
      if (j < i) {
         tdidx.xZ[i] = imax2(tdidx.xZ[i], i - j);
      } else {
         tuidx.xZ[j] = imax2(tuidx.xZ[j], j - i);
      }
   }
   vectorsetlengthatleast(&tridx, n + 1);
   tridx.xZ[0] = 0;
   for (i = 1; i <= n; i++) {
      tridx.xZ[i] = tridx.xZ[i - 1] + tdidx.xZ[i - 1] + 1 + tuidx.xZ[i - 1];
   }
   vectorsetlengthatleast(&tvals, tridx.xZ[n]);
   k = tridx.xZ[n];
   for (i = 0; i < k; i++) {
      tvals.xR[i] = 0.0;
   }
   t0 = 0;
   t1 = 0;
   while (sparseenumerate(s, &t0, &t1, &i, &j, &v)) {
      if (j <= i) {
         tvals.xR[tridx.xZ[i] + tdidx.xZ[i] - (i - j)] = v;
      } else {
         tvals.xR[tridx.xZ[j + 1] - (j - i)] = v;
      }
   }
   for (i = 0; i < n; i++) {
      tdidx.xZ[n] = imax2(tdidx.xZ[n], tdidx.xZ[i]);
      tuidx.xZ[n] = imax2(tuidx.xZ[n], tuidx.xZ[i]);
   }
   s->matrixtype = 2;
   s->ninitialized = 0;
   s->nfree = 0;
   s->m = n;
   s->n = n;
   ae_swap_vectors(&s->didx, &tdidx);
   ae_swap_vectors(&s->uidx, &tuidx);
   ae_swap_vectors(&s->ridx, &tridx);
   ae_swap_vectors(&s->vals, &tvals);
   ae_frame_leave();
}

// This  function  performs  out-of-place  conversion  to SKS storage format.
// S0 is copied to S1 and converted on-the-fly.
//
// Inputs:
//     S0          -   sparse matrix in any format.
//
// Outputs:
//     S1          -   sparse matrix in SKS format.
//
// NOTE: if S0 is stored as SKS, it is just copied without conversion.
//
// NOTE: this function de-allocates memory  occupied  by  S1 before  starting
//       conversion. If you perform a  lot  of  repeated  conversions, it may
//       lead to memory fragmentation. In this case we recommend you  to  use
//       SparseCopyToSKSBuf() function which re-uses memory in S1 as much  as
//       possible.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparsecopytosks(const sparsematrix &s0, sparsematrix &s1);
void sparsecopytosks(sparsematrix *s0, sparsematrix *s1) {
   SetObj(sparsematrix, s1);
   ae_assert(s0->matrixtype == 0 || s0->matrixtype == 1 || s0->matrixtype == 2, "SparseCopyToSKS: invalid matrix type");
   sparsecopytosksbuf(s0, s1);
}

// This  function  performs  out-of-place  conversion  to SKS format.  S0  is
// copied to S1 and converted on-the-fly. Memory  allocated  in S1 is  reused
// to maximum extent possible.
//
// Inputs:
//     S0          -   sparse matrix in any format.
//
// Outputs:
//     S1          -   sparse matrix in SKS format.
//
// NOTE: if S0 is stored as SKS, it is just copied without conversion.
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: void sparsecopytosksbuf(const sparsematrix &s0, const sparsematrix &s1);
void sparsecopytosksbuf(sparsematrix *s0, sparsematrix *s1) {
   double v;
   ae_int_t n;
   ae_int_t t0;
   ae_int_t t1;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_assert(s0->matrixtype == 0 || s0->matrixtype == 1 || s0->matrixtype == 2, "SparseCopyToSKSBuf: invalid matrix type");
   ae_assert(s0->m == s0->n, "SparseCopyToSKSBuf: rectangular matrices are not supported");
   n = s0->n;
   if (s0->matrixtype == 2) {
   // Already SKS, just copy
      sparsecopybuf(s0, s1);
      return;
   }
// Generate copy of matrix in the SKS format
   vectorsetlengthatleast(&s1->didx, n + 1);
   vectorsetlengthatleast(&s1->uidx, n + 1);
   for (i = 0; i <= n; i++) {
      s1->didx.xZ[i] = 0;
      s1->uidx.xZ[i] = 0;
   }
   t0 = 0;
   t1 = 0;
   while (sparseenumerate(s0, &t0, &t1, &i, &j, &v)) {
      if (j < i) {
         s1->didx.xZ[i] = imax2(s1->didx.xZ[i], i - j);
      } else {
         s1->uidx.xZ[j] = imax2(s1->uidx.xZ[j], j - i);
      }
   }
   vectorsetlengthatleast(&s1->ridx, n + 1);
   s1->ridx.xZ[0] = 0;
   for (i = 1; i <= n; i++) {
      s1->ridx.xZ[i] = s1->ridx.xZ[i - 1] + s1->didx.xZ[i - 1] + 1 + s1->uidx.xZ[i - 1];
   }
   vectorsetlengthatleast(&s1->vals, s1->ridx.xZ[n]);
   k = s1->ridx.xZ[n];
   for (i = 0; i < k; i++) {
      s1->vals.xR[i] = 0.0;
   }
   t0 = 0;
   t1 = 0;
   while (sparseenumerate(s0, &t0, &t1, &i, &j, &v)) {
      if (j <= i) {
         s1->vals.xR[s1->ridx.xZ[i] + s1->didx.xZ[i] - (i - j)] = v;
      } else {
         s1->vals.xR[s1->ridx.xZ[j + 1] - (j - i)] = v;
      }
   }
   for (i = 0; i < n; i++) {
      s1->didx.xZ[n] = imax2(s1->didx.xZ[n], s1->didx.xZ[i]);
      s1->uidx.xZ[n] = imax2(s1->uidx.xZ[n], s1->uidx.xZ[i]);
   }
   s1->matrixtype = 2;
   s1->ninitialized = 0;
   s1->nfree = 0;
   s1->m = n;
   s1->n = n;
}

// This non-accessible to user function performs  in-place  creation  of  CRS
// matrix. It is expected that:
// * S.M and S.N are initialized
// * S.RIdx, S.Idx and S.Vals are loaded with values in CRS  format  used  by
//   ALGLIB, with elements of S.Idx/S.Vals  possibly  being  unsorted  within
//   each row (this constructor function may post-sort matrix,  assuming that
//   it is sorted by rows).
//
// Only 5 fields should be set by caller. Other fields will be  rewritten  by
// this constructor function.
//
// This function performs integrity check on user-specified values, with  the
// only exception being Vals[] array:
// * it does not require values to be non-zero
// * it does not checks for element of Vals[] being finite IEEE-754 values
//
// Inputs:
//     S   -   sparse matrix with corresponding fields set by caller
//
// Outputs:
//     S   -   sparse matrix in CRS format.
//
// ALGLIB Project: Copyright 20.08.2016 by Sergey Bochkanov
void sparsecreatecrsinplace(sparsematrix *s) {
   ae_int_t m;
   ae_int_t n;
   ae_int_t i;
   ae_int_t j;
   ae_int_t j0;
   ae_int_t j1;
   m = s->m;
   n = s->n;
// Quick exit for M=0 or N=0
   ae_assert(s->m >= 0, "SparseCreateCRSInplace: integrity check failed");
   ae_assert(s->n >= 0, "SparseCreateCRSInplace: integrity check failed");
   if (m == 0 || n == 0) {
      s->matrixtype = 1;
      s->ninitialized = 0;
      vectorsetlengthatleast(&s->ridx, s->m + 1);
      vectorsetlengthatleast(&s->didx, s->m);
      vectorsetlengthatleast(&s->uidx, s->m);
      for (i = 0; i < s->m; i++) {
         s->ridx.xZ[i] = 0;
         s->uidx.xZ[i] = 0;
         s->didx.xZ[i] = 0;
      }
      s->ridx.xZ[s->m] = 0;
      return;
   }
// Perform integrity check
   ae_assert(s->m > 0, "SparseCreateCRSInplace: integrity check failed");
   ae_assert(s->n > 0, "SparseCreateCRSInplace: integrity check failed");
   ae_assert(s->ridx.cnt >= m + 1, "SparseCreateCRSInplace: integrity check failed");
   for (i = 0; i < m; i++) {
      ae_assert(s->ridx.xZ[i] >= 0 && s->ridx.xZ[i] <= s->ridx.xZ[i + 1], "SparseCreateCRSInplace: integrity check failed");
   }
   ae_assert(s->ridx.xZ[m] <= s->idx.cnt, "SparseCreateCRSInplace: integrity check failed");
   ae_assert(s->ridx.xZ[m] <= s->vals.cnt, "SparseCreateCRSInplace: integrity check failed");
   for (i = 0; i < m; i++) {
      j0 = s->ridx.xZ[i];
      j1 = s->ridx.xZ[i + 1] - 1;
      for (j = j0; j <= j1; j++) {
         ae_assert(s->idx.xZ[j] >= 0 && s->idx.xZ[j] < n, "SparseCreateCRSInplace: integrity check failed");
      }
   }
// Initialize
   s->matrixtype = 1;
   s->ninitialized = s->ridx.xZ[m];
   for (i = 0; i < m; i++) {
      tagsortmiddleir(&s->idx, &s->vals, s->ridx.xZ[i + 1] - s->ridx.xZ[i], s->ridx.xZ[i]);
   }
   sparseinitduidx(s);
}

// This function returns type of the matrix storage format.
//
// Inputs:
//     S           -   sparse matrix.
//
// Result:
//     sparse storage format used by matrix:
//         0   -   Hash-table
//         1   -   CRS (compressed row storage)
//         2   -   SKS (skyline)
//
// NOTE: future  versions  of  ALGLIB  may  include additional sparse storage
//       formats.
//
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: ae_int_t sparsegetmatrixtype(const sparsematrix &s);
ae_int_t sparsegetmatrixtype(sparsematrix *s) {
   ae_int_t result;
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseGetMatrixType: invalid matrix type");
   result = s->matrixtype;
   return result;
}

// This function checks matrix storage format and returns True when matrix is
// stored using Hash table representation.
//
// Inputs:
//     S   -   sparse matrix.
//
// Result:
//     True if matrix type is Hash table
//     False if matrix type is not Hash table
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: bool sparseishash(const sparsematrix &s);
bool sparseishash(sparsematrix *s) {
   bool result;
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseIsHash: invalid matrix type");
   result = s->matrixtype == 0;
   return result;
}

// This function checks matrix storage format and returns True when matrix is
// stored using CRS representation.
//
// Inputs:
//     S   -   sparse matrix.
//
// Result:
//     True if matrix type is CRS
//     False if matrix type is not CRS
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: bool sparseiscrs(const sparsematrix &s);
bool sparseiscrs(sparsematrix *s) {
   bool result;
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseIsCRS: invalid matrix type");
   result = s->matrixtype == 1;
   return result;
}

// This function checks matrix storage format and returns True when matrix is
// stored using SKS representation.
//
// Inputs:
//     S   -   sparse matrix.
//
// Result:
//     True if matrix type is SKS
//     False if matrix type is not SKS
//
// ALGLIB Project: Copyright 20.07.2012 by Sergey Bochkanov
// API: bool sparseissks(const sparsematrix &s);
bool sparseissks(sparsematrix *s) {
   bool result;
   ae_assert(s->matrixtype == 0 || s->matrixtype == 1 || s->matrixtype == 2, "SparseIsSKS: invalid matrix type");
   result = s->matrixtype == 2;
   return result;
}

// The function frees all memory occupied by  sparse  matrix.  Sparse  matrix
// structure becomes unusable after this call.
//
// Outputs:
//     S   -   sparse matrix to delete
//
// ALGLIB Project: Copyright 24.07.2012 by Sergey Bochkanov
// API: void sparsefree(sparsematrix &s);
void sparsefree(sparsematrix *s) {
   SetObj(sparsematrix, s);
   s->matrixtype = -1;
   s->m = 0;
   s->n = 0;
   s->nfree = 0;
   s->ninitialized = 0;
   s->tablesize = 0;
}

// The function returns number of columns of a sparse matrix.
//
// Result: number of columns of a sparse matrix.
//
// ALGLIB Project: Copyright 23.08.2012 by Sergey Bochkanov
// API: ae_int_t sparsegetncols(const sparsematrix &s);
ae_int_t sparsegetncols(sparsematrix *s) {
   ae_int_t result;
   result = s->n;
   return result;
}

// The function returns number of rows of a sparse matrix.
//
// Result: number of rows of a sparse matrix.
//
// ALGLIB Project: Copyright 23.08.2012 by Sergey Bochkanov
// API: ae_int_t sparsegetnrows(const sparsematrix &s);
ae_int_t sparsegetnrows(sparsematrix *s) {
   ae_int_t result;
   result = s->m;
   return result;
}

// The function returns number of strictly upper triangular non-zero elements
// in  the  matrix.  It  counts  SYMBOLICALLY non-zero elements, i.e. entries
// in the sparse matrix data structure. If some element  has  zero  numerical
// value, it is still counted.
//
// This function has different cost for different types of matrices:
// * for hash-based matrices it involves complete pass over entire hash-table
//   with O(NNZ) cost, where NNZ is number of non-zero elements
// * for CRS and SKS matrix types cost of counting is O(N) (N - matrix size).
//
// Result: number of non-zero elements strictly above main diagonal
//
// ALGLIB Project: Copyright 12.02.2014 by Sergey Bochkanov
// API: ae_int_t sparsegetuppercount(const sparsematrix &s);
ae_int_t sparsegetuppercount(sparsematrix *s) {
   ae_int_t sz;
   ae_int_t i0;
   ae_int_t i;
   ae_int_t result;
   result = -1;
   if (s->matrixtype == 0) {
   // Hash-table matrix
      result = 0;
      sz = s->tablesize;
      for (i0 = 0; i0 < sz; i0++) {
         i = s->idx.xZ[2 * i0];
         if (i >= 0 && s->idx.xZ[2 * i0 + 1] > i) {
            result++;
         }
      }
      return result;
   }
   if (s->matrixtype == 1) {
   // CRS matrix
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseGetUpperCount: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      result = 0;
      sz = s->m;
      for (i = 0; i < sz; i++) {
         result += s->ridx.xZ[i + 1] - s->uidx.xZ[i];
      }
      return result;
   }
   if (s->matrixtype == 2) {
   // SKS matrix
      ae_assert(s->m == s->n, "SparseGetUpperCount: non-square SKS matrices are not supported");
      result = 0;
      sz = s->m;
      for (i = 0; i < sz; i++) {
         result += s->uidx.xZ[i];
      }
      return result;
   }
   ae_assert(false, "SparseGetUpperCount: internal error");
   return result;
}

// The function returns number of strictly lower triangular non-zero elements
// in  the  matrix.  It  counts  SYMBOLICALLY non-zero elements, i.e. entries
// in the sparse matrix data structure. If some element  has  zero  numerical
// value, it is still counted.
//
// This function has different cost for different types of matrices:
// * for hash-based matrices it involves complete pass over entire hash-table
//   with O(NNZ) cost, where NNZ is number of non-zero elements
// * for CRS and SKS matrix types cost of counting is O(N) (N - matrix size).
//
// Result: number of non-zero elements strictly below main diagonal
//
// ALGLIB Project: Copyright 12.02.2014 by Sergey Bochkanov
// API: ae_int_t sparsegetlowercount(const sparsematrix &s);
ae_int_t sparsegetlowercount(sparsematrix *s) {
   ae_int_t sz;
   ae_int_t i0;
   ae_int_t i;
   ae_int_t result;
   result = -1;
   if (s->matrixtype == 0) {
   // Hash-table matrix
      result = 0;
      sz = s->tablesize;
      for (i0 = 0; i0 < sz; i0++) {
         i = s->idx.xZ[2 * i0];
         if (i >= 0 && s->idx.xZ[2 * i0 + 1] < i) {
            result++;
         }
      }
      return result;
   }
   if (s->matrixtype == 1) {
   // CRS matrix
      ae_assert(s->ninitialized == s->ridx.xZ[s->m], "SparseGetUpperCount: some rows/elements of the CRS matrix were not initialized (you must initialize everything you promised to SparseCreateCRS)");
      result = 0;
      sz = s->m;
      for (i = 0; i < sz; i++) {
         result += s->didx.xZ[i] - s->ridx.xZ[i];
      }
      return result;
   }
   if (s->matrixtype == 2) {
   // SKS matrix
      ae_assert(s->m == s->n, "SparseGetUpperCount: non-square SKS matrices are not supported");
      result = 0;
      sz = s->m;
      for (i = 0; i < sz; i++) {
         result += s->didx.xZ[i];
      }
      return result;
   }
   ae_assert(false, "SparseGetUpperCount: internal error");
   return result;
}

void sparsematrix_init(void *_p, bool make_automatic) {
   sparsematrix *p = (sparsematrix *)_p;
   ae_vector_init(&p->vals, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->idx, 0, DT_INT, make_automatic);
   ae_vector_init(&p->ridx, 0, DT_INT, make_automatic);
   ae_vector_init(&p->didx, 0, DT_INT, make_automatic);
   ae_vector_init(&p->uidx, 0, DT_INT, make_automatic);
}

void sparsematrix_copy(void *_dst, void *_src, bool make_automatic) {
   sparsematrix *dst = (sparsematrix *)_dst;
   sparsematrix *src = (sparsematrix *)_src;
   ae_vector_copy(&dst->vals, &src->vals, make_automatic);
   ae_vector_copy(&dst->idx, &src->idx, make_automatic);
   ae_vector_copy(&dst->ridx, &src->ridx, make_automatic);
   ae_vector_copy(&dst->didx, &src->didx, make_automatic);
   ae_vector_copy(&dst->uidx, &src->uidx, make_automatic);
   dst->matrixtype = src->matrixtype;
   dst->m = src->m;
   dst->n = src->n;
   dst->nfree = src->nfree;
   dst->ninitialized = src->ninitialized;
   dst->tablesize = src->tablesize;
}

void sparsematrix_free(void *_p, bool make_automatic) {
   sparsematrix *p = (sparsematrix *)_p;
   ae_vector_free(&p->vals, make_automatic);
   ae_vector_free(&p->idx, make_automatic);
   ae_vector_free(&p->ridx, make_automatic);
   ae_vector_free(&p->didx, make_automatic);
   ae_vector_free(&p->uidx, make_automatic);
}

void sparsebuffers_init(void *_p, bool make_automatic) {
   sparsebuffers *p = (sparsebuffers *)_p;
   ae_vector_init(&p->d, 0, DT_INT, make_automatic);
   ae_vector_init(&p->u, 0, DT_INT, make_automatic);
   sparsematrix_init(&p->s, make_automatic);
}

void sparsebuffers_copy(void *_dst, void *_src, bool make_automatic) {
   sparsebuffers *dst = (sparsebuffers *)_dst;
   sparsebuffers *src = (sparsebuffers *)_src;
   ae_vector_copy(&dst->d, &src->d, make_automatic);
   ae_vector_copy(&dst->u, &src->u, make_automatic);
   sparsematrix_copy(&dst->s, &src->s, make_automatic);
}

void sparsebuffers_free(void *_p, bool make_automatic) {
   sparsebuffers *p = (sparsebuffers *)_p;
   ae_vector_free(&p->d, make_automatic);
   ae_vector_free(&p->u, make_automatic);
   sparsematrix_free(&p->s, make_automatic);
}
} // end of namespace alglib_impl

namespace alglib {
// Sparse matrix structure.
// You should use ALGLIB functions to work with sparse matrix. Never  try  to
// access its fields directly!
//
// NOTES ON THE SPARSE STORAGE FORMATS
//
// Sparse matrices can be stored using several formats:
// * Hash-Table representation
// * Compressed Row Storage (CRS)
// * Skyline matrix storage (SKS)
//
// Each of the formats has benefits and drawbacks:
// * Hash-table is good for dynamic operations (insertion of new elements),
//   but does not support linear algebra operations
// * CRS is good for operations like matrix-vector or matrix-matrix products,
//   but its initialization is less convenient - you have to tell row   sizes
//   at the initialization, and you have to fill  matrix  only  row  by  row,
//   from left to right.
// * SKS is a special format which is used to store triangular  factors  from
//   Cholesky factorization. It does not support  dynamic  modification,  and
//   support for linear algebra operations is very limited.
//
// Tables below outline information about these two formats:
//
//     OPERATIONS WITH MATRIX      HASH        CRS         SKS
//     creation                    +           +           +
//     SparseGet                   +           +           +
//     SparseRewriteExisting       +           +           +
//     SparseSet                   +           +           +
//     SparseAdd                   +
//     SparseGetRow                            +           +
//     SparseGetCompressedRow                  +           +
//     sparse-dense linear algebra             +           +
DefClass(sparsematrix, EndD)

// Temporary buffers for sparse matrix operations.
//
// You should pass an instance of this structure to factorization  functions.
// It allows to reuse memory during repeated sparse  factorizations.  You  do
// not have to call some initialization function - simply passing an instance
// to factorization function is enough.
DefClass(sparsebuffers, EndD)

void sparsecreate(const ae_int_t m, const ae_int_t n, const ae_int_t k, sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreate(m, n, k, ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void sparsecreate(const ae_int_t m, const ae_int_t n, sparsematrix &s) {
   ae_int_t k = 0;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreate(m, n, k, ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}
#endif

void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const ae_int_t k, const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatebuf(m, n, k, ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const sparsematrix &s) {
   ae_int_t k = 0;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatebuf(m, n, k, ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}
#endif

void sparsecreatecrs(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatecrs(m, n, ConstT(ae_vector, ner), ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecreatecrsbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatecrsbuf(m, n, ConstT(ae_vector, ner), ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecreatesks(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatesks(m, n, ConstT(ae_vector, d), ConstT(ae_vector, u), ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecreatesksbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatesksbuf(m, n, ConstT(ae_vector, d), ConstT(ae_vector, u), ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecreatesksband(const ae_int_t m, const ae_int_t n, const ae_int_t bw, sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatesksband(m, n, bw, ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecreatesksbandbuf(const ae_int_t m, const ae_int_t n, const ae_int_t bw, const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecreatesksbandbuf(m, n, bw, ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecopy(const sparsematrix &s0, sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopy(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparsecopybuf(const sparsematrix &s0, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopybuf(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparseswap(const sparsematrix &s0, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseswap(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparseadd(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseadd(ConstT(sparsematrix, s), i, j, v);
   alglib_impl::ae_state_clear();
}

void sparseset(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseset(ConstT(sparsematrix, s), i, j, v);
   alglib_impl::ae_state_clear();
}

double sparseget(const sparsematrix &s, const ae_int_t i, const ae_int_t j) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::sparseget(ConstT(sparsematrix, s), i, j);
   alglib_impl::ae_state_clear();
   return D;
}

double sparsegetdiagonal(const sparsematrix &s, const ae_int_t i) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::sparsegetdiagonal(ConstT(sparsematrix, s), i);
   alglib_impl::ae_state_clear();
   return D;
}

void sparsemv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsemv(ConstT(sparsematrix, s), ConstT(ae_vector, x), ConstT(ae_vector, y));
   alglib_impl::ae_state_clear();
}

void sparsemtv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsemtv(ConstT(sparsematrix, s), ConstT(ae_vector, x), ConstT(ae_vector, y));
   alglib_impl::ae_state_clear();
}

void sparsegemv(const sparsematrix &s, const double alpha, const ae_int_t ops, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsegemv(ConstT(sparsematrix, s), alpha, ops, ConstT(ae_vector, x), ix, beta, ConstT(ae_vector, y), iy);
   alglib_impl::ae_state_clear();
}

void sparsemv2(const sparsematrix &s, const real_1d_array &x, real_1d_array &y0, real_1d_array &y1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsemv2(ConstT(sparsematrix, s), ConstT(ae_vector, x), ConstT(ae_vector, y0), ConstT(ae_vector, y1));
   alglib_impl::ae_state_clear();
}

void sparsesmv(const sparsematrix &s, const bool isupper, const real_1d_array &x, real_1d_array &y) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsesmv(ConstT(sparsematrix, s), isupper, ConstT(ae_vector, x), ConstT(ae_vector, y));
   alglib_impl::ae_state_clear();
}

double sparsevsmv(const sparsematrix &s, const bool isupper, const real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::sparsevsmv(ConstT(sparsematrix, s), isupper, ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
   return D;
}

void sparsemm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsemm(ConstT(sparsematrix, s), ConstT(ae_matrix, a), k, ConstT(ae_matrix, b));
   alglib_impl::ae_state_clear();
}

void sparsemtm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsemtm(ConstT(sparsematrix, s), ConstT(ae_matrix, a), k, ConstT(ae_matrix, b));
   alglib_impl::ae_state_clear();
}

void sparsemm2(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b0, real_2d_array &b1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsemm2(ConstT(sparsematrix, s), ConstT(ae_matrix, a), k, ConstT(ae_matrix, b0), ConstT(ae_matrix, b1));
   alglib_impl::ae_state_clear();
}

void sparsesmm(const sparsematrix &s, const bool isupper, const real_2d_array &a, const ae_int_t k, real_2d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsesmm(ConstT(sparsematrix, s), isupper, ConstT(ae_matrix, a), k, ConstT(ae_matrix, b));
   alglib_impl::ae_state_clear();
}

void sparsetrmv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, real_1d_array &y) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsetrmv(ConstT(sparsematrix, s), isupper, isunit, optype, ConstT(ae_vector, x), ConstT(ae_vector, y));
   alglib_impl::ae_state_clear();
}

void sparsetrsv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsetrsv(ConstT(sparsematrix, s), isupper, isunit, optype, ConstT(ae_vector, x));
   alglib_impl::ae_state_clear();
}

void sparseresizematrix(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseresizematrix(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

bool sparseenumerate(const sparsematrix &s, ae_int_t &t0, ae_int_t &t1, ae_int_t &i, ae_int_t &j, double &v) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparseenumerate(ConstT(sparsematrix, s), &t0, &t1, &i, &j, &v);
   alglib_impl::ae_state_clear();
   return Ok;
}

bool sparserewriteexisting(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparserewriteexisting(ConstT(sparsematrix, s), i, j, v);
   alglib_impl::ae_state_clear();
   return Ok;
}

void sparsegetrow(const sparsematrix &s, const ae_int_t i, real_1d_array &irow) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsegetrow(ConstT(sparsematrix, s), i, ConstT(ae_vector, irow));
   alglib_impl::ae_state_clear();
}

void sparsegetcompressedrow(const sparsematrix &s, const ae_int_t i, integer_1d_array &colidx, real_1d_array &vals, ae_int_t &nzcnt) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsegetcompressedrow(ConstT(sparsematrix, s), i, ConstT(ae_vector, colidx), ConstT(ae_vector, vals), &nzcnt);
   alglib_impl::ae_state_clear();
}

void sparsetransposesks(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsetransposesks(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsetransposecrs(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsetransposecrs(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecopytransposecrsbuf(const sparsematrix &s0, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytransposecrsbuf(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparsecopytransposecrs(const sparsematrix &s0, sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytransposecrs(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparseconvertto(const sparsematrix &s0, const ae_int_t fmt) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseconvertto(ConstT(sparsematrix, s0), fmt);
   alglib_impl::ae_state_clear();
}

void sparsecopytobuf(const sparsematrix &s0, const ae_int_t fmt, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytobuf(ConstT(sparsematrix, s0), fmt, ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparseconverttohash(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseconverttohash(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecopytohash(const sparsematrix &s0, sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytohash(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparsecopytohashbuf(const sparsematrix &s0, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytohashbuf(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparseconverttocrs(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseconverttocrs(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecopytocrs(const sparsematrix &s0, sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytocrs(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparsecopytocrsbuf(const sparsematrix &s0, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytocrsbuf(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparseconverttosks(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparseconverttosks(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

void sparsecopytosks(const sparsematrix &s0, sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytosks(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

void sparsecopytosksbuf(const sparsematrix &s0, const sparsematrix &s1) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsecopytosksbuf(ConstT(sparsematrix, s0), ConstT(sparsematrix, s1));
   alglib_impl::ae_state_clear();
}

ae_int_t sparsegetmatrixtype(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::sparsegetmatrixtype(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Z;
}

bool sparseishash(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparseishash(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool sparseiscrs(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparseiscrs(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool sparseissks(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparseissks(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Ok;
}

void sparsefree(sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sparsefree(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
}

ae_int_t sparsegetncols(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::sparsegetncols(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Z;
}

ae_int_t sparsegetnrows(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::sparsegetnrows(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Z;
}

ae_int_t sparsegetuppercount(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::sparsegetuppercount(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Z;
}

ae_int_t sparsegetlowercount(const sparsematrix &s) {
   alglib_impl::ae_state_init();
   TryCatch(0)
   ae_int_t Z = alglib_impl::sparsegetlowercount(ConstT(sparsematrix, s));
   alglib_impl::ae_state_clear();
   return Z;
}
} // end of namespace alglib

// === ABLAS Package ===
// Depends on: (AlgLibInternal) APSERV, ABLASMKL, ABLASF
namespace alglib_impl {
static const ae_int_t ablas_blas2minvendorkernelsize = 8;

// Complex ABLASSplitLength
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
static ae_int_t ablas_ablasinternalsplitlength(ae_int_t n, ae_int_t nb) {
   if (n <= nb) ; // Block size: no further splitting.
// Greater than the block size.
   else if (n % nb != 0) n -= n % nb; // Split the remainder.
   else { // Split on the block boundaries.
      n -= n / 2;
      if (n % nb != 0) n += nb - n % nb;
   }
   return n;
}

// Splits matrix length in two parts, left part should match ABLAS block size
//
// Inputs:
//     A   -   real matrix, is passed to ensure that we didn't split
//             complex matrix using real splitting subroutine.
//             matrix itself is not changed.
//     N   -   length, N > 0
//
// Output:
//     N1  -   length of the left half
//
// N1 >= N/2, N1 may be N
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
ae_int_t ablassplitlength(RMatrix *a, ae_int_t n) {
   return ablas_ablasinternalsplitlength(n, n > ablasblocksize(a) ? ablasblocksize(a) : ablasmicroblocksize());
}

// Complex ABLASSplitLength
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
ae_int_t ablascomplexsplitlength(CMatrix *a, ae_int_t n) {
   return ablas_ablasinternalsplitlength(n, n > ablascomplexblocksize(a) ? ablascomplexblocksize(a) : ablasmicroblocksize());
}

// Returns switch point for parallelism.
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
ae_int_t gemmparallelsize() {
   ae_int_t result;
   result = 64;
   return result;
}

// Returns block size - subdivision size where  cache-oblivious  soubroutines
// switch to the optimized kernel.
//
// Inputs:
//     A   -   real matrix, is passed to ensure that we didn't split
//             complex matrix using real splitting subroutine.
//             matrix itself is not changed.
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
ae_int_t ablasblocksize(RMatrix *a) {
   ae_int_t result;
   result = 32;
   return result;
}

// Block size for complex subroutines.
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
ae_int_t ablascomplexblocksize(CMatrix *a) {
   ae_int_t result;
   result = 24;
   return result;
}

// Microblock size
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
ae_int_t ablasmicroblocksize() {
   ae_int_t result;
   result = 8;
   return result;
}

// Generation of an elementary reflection transformation
//
// The subroutine generates elementary reflection H of order N, so that, for
// a given X, the following equality holds true:
//
//     ( X(1) )   ( Beta )
// H * (  ..  ) = (  0   )
//     ( X(n) )   (  0   )
//
// where
//               ( V(1) )
// H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
//               ( V(n) )
//
// where the first component of vector V equals 1.
//
// Inputs:
//     X   -   vector. Array whose index ranges within [1..N].
//     N   -   reflection order.
//
// Outputs:
//     X   -   components from 2 to N are replaced with vector V.
//             The first component is replaced with parameter Beta.
//     Tau -   scalar value Tau. If X is a null vector, Tau equals 0,
//             otherwise 1 <= Tau <= 2.
//
// This subroutine is the modification of the DLARFG subroutines from
// the LAPACK library.
//
// MODIFICATIONS:
//     24.12.2005 sign(Alpha) was replaced with an analogous to the Fortran SIGN code.
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
void generatereflection(RVector *x, ae_int_t n, double *tau) {
   ae_int_t j;
   double alpha;
   double xnorm;
   double v;
   double beta;
   double mx;
   double s;
   *tau = 0;
   if (n <= 1) {
      *tau = 0.0;
      return;
   }
// Scale if needed (to avoid overflow/underflow during intermediate
// calculations).
   mx = 0.0;
   for (j = 1; j <= n; j++) {
      mx = rmax2(fabs(x->xR[j]), mx);
   }
   s = 1.0;
   if (mx != 0.0) {
      if (mx <= ae_minrealnumber / ae_machineepsilon) {
         s = ae_minrealnumber / ae_machineepsilon;
         v = 1 / s;
         ae_v_muld(&x->xR[1], 1, n, v);
         mx *= v;
      } else {
         if (mx >= ae_maxrealnumber * ae_machineepsilon) {
            s = ae_maxrealnumber * ae_machineepsilon;
            v = 1 / s;
            ae_v_muld(&x->xR[1], 1, n, v);
            mx *= v;
         }
      }
   }
// XNORM = DNRM2( N-1, X, INCX )
   alpha = x->xR[1];
   xnorm = 0.0;
   if (mx != 0.0) {
      for (j = 2; j <= n; j++) {
         xnorm += ae_sqr(x->xR[j] / mx);
      }
      xnorm = sqrt(xnorm) * mx;
   }
   if (xnorm == 0.0) {
   // H  =  I
      *tau = 0.0;
      x->xR[1] *= s;
      return;
   }
// general case
   mx = rmax2(fabs(alpha), fabs(xnorm));
   beta = -mx * sqrt(ae_sqr(alpha / mx) + ae_sqr(xnorm / mx));
   if (alpha < 0.0) {
      beta = -beta;
   }
   *tau = (beta - alpha) / beta;
   v = 1 / (alpha - beta);
   ae_v_muld(&x->xR[2], 1, n - 1, v);
   x->xR[1] = beta;
// Scale back outputs
   x->xR[1] *= s;
}

// Application of an elementary reflection to a rectangular matrix of size MxN
//
// The algorithm pre-multiplies the matrix by an elementary reflection transformation
// which is given by column V and scalar Tau (see the description of the
// GenerateReflection procedure). Not the whole matrix but only a part of it
// is transformed (rows from M1 to M2, columns from N1 to N2). Only the elements
// of this submatrix are changed.
//
// Inputs:
//     C       -   matrix to be transformed.
//     Tau     -   scalar defining the transformation.
//     V       -   column defining the transformation.
//                 Array whose index ranges within [1..M2-M1+1].
//     M1, M2  -   range of rows to be transformed.
//     N1, N2  -   range of columns to be transformed.
//     WORK    -   working array whose indexes goes from N1 to N2.
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
void applyreflectionfromtheleft(RMatrix *c, double tau, RVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *work) {
   if (tau == 0.0 || n1 > n2 || m1 > m2) {
      return;
   }
   vectorsetlengthatleast(work, n2 - n1 + 1);
   rmatrixgemv(n2 - n1 + 1, m2 - m1 + 1, 1.0, c, m1, n1, 1, v, 1, 0.0, work, 0);
   rmatrixger(m2 - m1 + 1, n2 - n1 + 1, c, m1, n1, -tau, v, 1, work, 0);
}

// Application of an elementary reflection to a rectangular matrix of size MxN
//
// The algorithm post-multiplies the matrix by an elementary reflection transformation
// which is given by column V and scalar Tau (see the description of the
// GenerateReflection procedure). Not the whole matrix but only a part of it
// is transformed (rows from M1 to M2, columns from N1 to N2). Only the
// elements of this submatrix are changed.
//
// Inputs:
//     C       -   matrix to be transformed.
//     Tau     -   scalar defining the transformation.
//     V       -   column defining the transformation.
//                 Array whose index ranges within [1..N2-N1+1].
//     M1, M2  -   range of rows to be transformed.
//     N1, N2  -   range of columns to be transformed.
//     WORK    -   working array whose indexes goes from M1 to M2.
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
void applyreflectionfromtheright(RMatrix *c, double tau, RVector *v, ae_int_t m1, ae_int_t m2, ae_int_t n1, ae_int_t n2, RVector *work) {
   if (tau == 0.0 || n1 > n2 || m1 > m2) {
      return;
   }
   vectorsetlengthatleast(work, m2 - m1 + 1);
   rmatrixgemv(m2 - m1 + 1, n2 - n1 + 1, 1.0, c, m1, n1, 0, v, 1, 0.0, work, 0);
   rmatrixger(m2 - m1 + 1, n2 - n1 + 1, c, m1, n1, -tau, work, 0, v, 1);
}

// Cache-oblivous real "copy-and-transpose"
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   source matrix, MxN submatrix is copied and transposed
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     B   -   destination matrix, must be large enough to store result
//     IB  -   submatrix offset (row index)
//     JB  -   submatrix offset (column index)
// API: void rmatrixtranspose(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_2d_array &b, const ae_int_t ib, const ae_int_t jb);
void rmatrixtranspose(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb) {
   ae_int_t i;
   ae_int_t s1;
   ae_int_t s2;
   if (m <= 2 * ablasblocksize(a) && n <= 2 * ablasblocksize(a)) {
   // base case
      for (i = 0; i < m; i++) {
         ae_v_move(&b->xyR[ib][jb + i], b->stride, &a->xyR[ia + i][ja], 1, n);
      }
   } else {
   // Cache-oblivious recursion
      if (m > n) {
         s1 = ablassplitlength(a, m), s2 = m - s1;
         rmatrixtranspose(s1, n, a, ia, ja, b, ib, jb);
         rmatrixtranspose(s2, n, a, ia + s1, ja, b, ib, jb + s1);
      } else {
         s1 = ablassplitlength(a, n), s2 = n - s1;
         rmatrixtranspose(m, s1, a, ia, ja, b, ib, jb);
         rmatrixtranspose(m, s2, a, ia, ja + s1, b, ib + s1, jb);
      }
   }
}

// Cache-oblivous complex "copy-and-transpose"
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   source matrix, MxN submatrix is copied and transposed
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     B   -   destination matrix, must be large enough to store result
//     IB  -   submatrix offset (row index)
//     JB  -   submatrix offset (column index)
// API: void cmatrixtranspose(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb);
void cmatrixtranspose(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CMatrix *b, ae_int_t ib, ae_int_t jb) {
   ae_int_t i;
   ae_int_t s1;
   ae_int_t s2;
   if (m <= 2 * ablascomplexblocksize(a) && n <= 2 * ablascomplexblocksize(a)) {
   // base case
      for (i = 0; i < m; i++) {
         ae_v_cmove(&b->xyC[ib][jb + i], b->stride, &a->xyC[ia + i][ja], 1, "N", n);
      }
   } else {
   // Cache-oblivious recursion
      if (m > n) {
         s1 = ablascomplexsplitlength(a, m), s2 = m - s1;
         cmatrixtranspose(s1, n, a, ia, ja, b, ib, jb);
         cmatrixtranspose(s2, n, a, ia + s1, ja, b, ib, jb + s1);
      } else {
         s1 = ablascomplexsplitlength(a, n), s2 = n - s1;
         cmatrixtranspose(m, s1, a, ia, ja, b, ib, jb);
         cmatrixtranspose(m, s2, a, ia, ja + s1, b, ib + s1, jb);
      }
   }
}

// This code enforces symmetricy of the matrix by copying Upper part to lower
// one (or vice versa).
//
// Inputs:
//     A   -   matrix
//     N   -   number of rows/columns
//     IsUpper - whether we want to copy upper triangle to lower one (True)
//             or vice versa (False).
// API: void rmatrixenforcesymmetricity(const real_2d_array &a, const ae_int_t n, const bool isupper);
void rmatrixenforcesymmetricity(RMatrix *a, ae_int_t n, bool isupper) {
   ae_int_t i;
   ae_int_t j;
   if (isupper) {
      for (i = 0; i < n; i++) {
         for (j = i + 1; j < n; j++) {
            a->xyR[j][i] = a->xyR[i][j];
         }
      }
   } else {
      for (i = 0; i < n; i++) {
         for (j = i + 1; j < n; j++) {
            a->xyR[i][j] = a->xyR[j][i];
         }
      }
   }
}

// Copy
//
// Inputs:
//     N   -   subvector size
//     A   -   source vector, N elements are copied
//     IA  -   source offset (first element index)
//     B   -   destination vector, must be large enough to store result
//     IB  -   destination offset (first element index)
// API: void rvectorcopy(const ae_int_t n, const real_1d_array &a, const ae_int_t ia, const real_1d_array &b, const ae_int_t ib);
void rvectorcopy(ae_int_t n, RVector *a, ae_int_t ia, RVector *b, ae_int_t ib) {
   ae_int_t i;
   if (n == 0) {
      return;
   }
   for (i = 0; i < n; i++) {
      b->xR[ib + i] = a->xR[ia + i];
   }
}

// Copy
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   source matrix, MxN submatrix is copied and transposed
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     B   -   destination matrix, must be large enough to store result
//     IB  -   submatrix offset (row index)
//     JB  -   submatrix offset (column index)
// API: void rmatrixcopy(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_2d_array &b, const ae_int_t ib, const ae_int_t jb);
void rmatrixcopy(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RMatrix *b, ae_int_t ib, ae_int_t jb) {
   ae_int_t i;
   if (m == 0 || n == 0) {
      return;
   }
   for (i = 0; i < m; i++) {
      ae_v_move(&b->xyR[ib + i][jb], 1, &a->xyR[ia + i][ja], 1, n);
   }
}

// Copy
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   source matrix, MxN submatrix is copied and transposed
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     B   -   destination matrix, must be large enough to store result
//     IB  -   submatrix offset (row index)
//     JB  -   submatrix offset (column index)
// API: void cmatrixcopy(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb);
void cmatrixcopy(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CMatrix *b, ae_int_t ib, ae_int_t jb) {
   ae_int_t i;
   if (m == 0 || n == 0) {
      return;
   }
   for (i = 0; i < m; i++) {
      ae_v_cmove(&b->xyC[ib + i][jb], 1, &a->xyC[ia + i][ja], 1, "N", n);
   }
}

// Performs generalized copy: B := Beta*B + Alpha*A.
//
// If Beta=0, then previous contents of B is simply ignored. If Alpha=0, then
// A is ignored and not referenced. If both Alpha and Beta  are  zero,  B  is
// filled by zeros.
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     Alpha-  coefficient
//     A   -   source matrix, MxN submatrix is copied and transposed
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     Beta-   coefficient
//     B   -   destination matrix, must be large enough to store result
//     IB  -   submatrix offset (row index)
//     JB  -   submatrix offset (column index)
// API: void rmatrixgencopy(const ae_int_t m, const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const double beta, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb);
void rmatrixgencopy(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, double beta, RMatrix *b, ae_int_t ib, ae_int_t jb) {
   ae_int_t i;
   ae_int_t j;
   if (m == 0 || n == 0) {
      return;
   }
// Zero-fill
   if (alpha == 0.0 && beta == 0.0) {
      for (i = 0; i < m; i++) {
         for (j = 0; j < n; j++) {
            b->xyR[ib + i][jb + j] = 0.0;
         }
      }
      return;
   }
// Inplace multiply
   if (alpha == 0.0) {
      for (i = 0; i < m; i++) {
         for (j = 0; j < n; j++) {
            b->xyR[ib + i][jb + j] *= beta;
         }
      }
      return;
   }
// Multiply and copy
   if (beta == 0.0) {
      for (i = 0; i < m; i++) {
         for (j = 0; j < n; j++) {
            b->xyR[ib + i][jb + j] = alpha * a->xyR[ia + i][ja + j];
         }
      }
      return;
   }
// Generic
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         b->xyR[ib + i][jb + j] = alpha * a->xyR[ia + i][ja + j] + beta * b->xyR[ib + i][jb + j];
      }
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
//     A   -   target matrix, MxN submatrix is updated
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     Alpha-  coefficient
//     U   -   vector #1
//     IU  -   subvector offset
//     V   -   vector #2
//     IV  -   subvector offset
//
//
// ALGLIB Routine: Copyright 16.10.2017 by Sergey Bochkanov
// API: void rmatrixger(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const double alpha, const real_1d_array &u, const ae_int_t iu, const real_1d_array &v, const ae_int_t iv);
void rmatrixger(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, double alpha, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv) {
   ae_int_t i;
   double s;
// Quick exit
   if (m <= 0 || n <= 0) {
      return;
   }
// Try fast kernels:
// * vendor kernel
// * internal kernel
   if (m > ablas_blas2minvendorkernelsize && n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel first
      if (rmatrixgermkl(m, n, a, ia, ja, alpha, u, iu, v, iv)) {
         return;
      }
   }
   if (rmatrixgerf(m, n, a, ia, ja, alpha, u, iu, v, iv)) {
      return;
   }
// Generic code
   for (i = 0; i < m; i++) {
      s = alpha * u->xR[iu + i];
      ae_v_addd(&a->xyR[ia + i][ja], 1, &v->xR[iv], 1, n, s);
   }
}

// IMPORTANT: this function is deprecated since ALGLIB 3.13. Use RMatrixGER()
//            which is more generic version of this function.
//
// Rank-1 correction: A := A + u*v'
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   target matrix, MxN submatrix is updated
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     U   -   vector #1
//     IU  -   subvector offset
//     V   -   vector #2
//     IV  -   subvector offset
// API: void rmatrixrank1(const ae_int_t m, const ae_int_t n, real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_1d_array &u, const ae_int_t iu, real_1d_array &v, const ae_int_t iv);
void rmatrixrank1(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, RVector *u, ae_int_t iu, RVector *v, ae_int_t iv) {
   ae_int_t i;
   double s;
// Quick exit
   if (m <= 0 || n <= 0) {
      return;
   }
// Try fast kernels:
// * vendor kernel
// * internal kernel
   if (m > ablas_blas2minvendorkernelsize && n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel first
      if (rmatrixrank1mkl(m, n, a, ia, ja, u, iu, v, iv)) {
         return;
      }
   }
   if (rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv)) {
      return;
   }
// Generic code
   for (i = 0; i < m; i++) {
      s = u->xR[iu + i];
      ae_v_addd(&a->xyR[ia + i][ja], 1, &v->xR[iv], 1, n, s);
   }
}

// Rank-1 correction: A := A + u*v'
//
// Inputs:
//     M   -   number of rows
//     N   -   number of columns
//     A   -   target matrix, MxN submatrix is updated
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     U   -   vector #1
//     IU  -   subvector offset
//     V   -   vector #2
//     IV  -   subvector offset
// API: void cmatrixrank1(const ae_int_t m, const ae_int_t n, complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_1d_array &u, const ae_int_t iu, complex_1d_array &v, const ae_int_t iv);
void cmatrixrank1(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, CVector *u, ae_int_t iu, CVector *v, ae_int_t iv) {
   ae_int_t i;
   ae_complex s;
// Quick exit
   if (m <= 0 || n <= 0) {
      return;
   }
// Try fast kernels:
// * vendor kernel
// * internal kernel
   if (m > ablas_blas2minvendorkernelsize && n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel first
      if (cmatrixrank1mkl(m, n, a, ia, ja, u, iu, v, iv)) {
         return;
      }
   }
   if (cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv)) {
      return;
   }
// Generic code
   for (i = 0; i < m; i++) {
      s = u->xC[iu + i];
      ae_v_caddc(&a->xyC[ia + i][ja], 1, &v->xC[iv], 1, "N", n, s);
   }
}

// Scaled matrix-vector addition: y += alpha a x.
// API: void rmatrixgemv(const ae_int_t m, const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy);
void rmatrixgemv(ae_int_t m, ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
   ae_int_t i;
   double v;
// Quick exit for M=0, N=0 or Alpha=0.
//
// After this block we have M > 0, N > 0, Alpha != 0.
   if (m <= 0) {
      return;
   }
   if (n <= 0 || alpha == 0.0) {
      if (beta != 0.0) {
         for (i = 0; i < m; i++) {
            y->xR[iy + i] *= beta;
         }
      } else {
         for (i = 0; i < m; i++) {
            y->xR[iy + i] = 0.0;
         }
      }
      return;
   }
// Try fast kernels
   if (m > ablas_blas2minvendorkernelsize && n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel
      if (rmatrixgemvmkl(m, n, alpha, a, ia, ja, opa, x, ix, beta, y, iy)) {
         return;
      }
   }
// Generic code
   if (opa == 0) {
   // y = A*x
      for (i = 0; i < m; i++) {
         v = ae_v_dotproduct(&a->xyR[ia + i][ja], 1, &x->xR[ix], 1, n);
         if (beta == 0.0) {
            y->xR[iy + i] = alpha * v;
         } else {
            y->xR[iy + i] = alpha * v + beta * y->xR[iy + i];
         }
      }
      return;
   }
   if (opa == 1) {
   // Prepare output array
      if (beta == 0.0) {
         for (i = 0; i < m; i++) {
            y->xR[iy + i] = 0.0;
         }
      } else {
         for (i = 0; i < m; i++) {
            y->xR[iy + i] *= beta;
         }
      }
   // y += A^T*x
      for (i = 0; i < n; i++) {
         v = alpha * x->xR[ix + i];
         ae_v_addd(&y->xR[iy], 1, &a->xyR[ia + i][ja], 1, m, v);
      }
      return;
   }
}

// IMPORTANT: this function is deprecated since ALGLIB 3.13. Use RMatrixGEMV()
//            which is more generic version of this function.
//
// Matrix-vector product: y := op(A)*x
//
// Inputs:
//     M   -   number of rows of op(A)
//     N   -   number of columns of op(A)
//     A   -   target matrix
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     OpA -   operation type:
//             * OpA=0     =>  op(A) = A
//             * OpA=1     =>  op(A) = A^T
//     X   -   input vector
//     IX  -   subvector offset
//     IY  -   subvector offset
//     Y   -   preallocated matrix, must be large enough to store result
//
// Outputs:
//     Y   -   vector which stores result
//
// if M=0, then subroutine does nothing.
// if N=0, Y is filled by zeros.
//
//
// ALGLIB Routine: Copyright 28.01.2010 by Sergey Bochkanov
// API: void rmatrixmv(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, real_1d_array &y, const ae_int_t iy);
void rmatrixmv(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RVector *x, ae_int_t ix, RVector *y, ae_int_t iy) {
   ae_int_t i;
   double v;
// Quick exit
   if (m == 0) {
      return;
   }
   if (n == 0) {
      for (i = 0; i < m; i++) {
         y->xR[iy + i] = 0.0;
      }
      return;
   }
// Try fast kernels
   if (m > ablas_blas2minvendorkernelsize && n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel
      if (rmatrixmvmkl(m, n, a, ia, ja, opa, x, ix, y, iy)) {
         return;
      }
   }
// Generic code
   if (opa == 0) {
   // y = A*x
      for (i = 0; i < m; i++) {
         v = ae_v_dotproduct(&a->xyR[ia + i][ja], 1, &x->xR[ix], 1, n);
         y->xR[iy + i] = v;
      }
      return;
   }
   if (opa == 1) {
   // y = A^T*x
      for (i = 0; i < m; i++) {
         y->xR[iy + i] = 0.0;
      }
      for (i = 0; i < n; i++) {
         v = x->xR[ix + i];
         ae_v_addd(&y->xR[iy], 1, &a->xyR[ia + i][ja], 1, m, v);
      }
      return;
   }
}

// Matrix-vector product: y := op(A)*x
//
// Inputs:
//     M   -   number of rows of op(A)
//             M >= 0
//     N   -   number of columns of op(A)
//             N >= 0
//     A   -   target matrix
//     IA  -   submatrix offset (row index)
//     JA  -   submatrix offset (column index)
//     OpA -   operation type:
//             * OpA=0     =>  op(A) = A
//             * OpA=1     =>  op(A) = A^T
//             * OpA=2     =>  op(A) = A^H
//     X   -   input vector
//     IX  -   subvector offset
//     IY  -   subvector offset
//     Y   -   preallocated matrix, must be large enough to store result
//
// Outputs:
//     Y   -   vector which stores result
//
// if M=0, then subroutine does nothing.
// if N=0, Y is filled by zeros.
//
//
// ALGLIB Routine: Copyright 28.01.2010 by Sergey Bochkanov
// API: void cmatrixmv(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const complex_1d_array &x, const ae_int_t ix, complex_1d_array &y, const ae_int_t iy);
void cmatrixmv(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, CVector *x, ae_int_t ix, CVector *y, ae_int_t iy) {
   ae_int_t i;
   ae_complex v;
// Quick exit
   if (m == 0) {
      return;
   }
   if (n == 0) {
      for (i = 0; i < m; i++) {
         y->xC[iy + i] = ae_complex_from_i(0);
      }
      return;
   }
// Try fast kernels
   if (m > ablas_blas2minvendorkernelsize && n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel
      if (cmatrixmvmkl(m, n, a, ia, ja, opa, x, ix, y, iy)) {
         return;
      }
   }
// Generic code
   if (opa == 0) {
   // y = A*x
      for (i = 0; i < m; i++) {
         v = ae_v_cdotproduct(&a->xyC[ia + i][ja], 1, "N", &x->xC[ix], 1, "N", n);
         y->xC[iy + i] = v;
      }
      return;
   }
   if (opa == 1) {
   // y = A^T*x
      for (i = 0; i < m; i++) {
         y->xC[iy + i] = ae_complex_from_i(0);
      }
      for (i = 0; i < n; i++) {
         v = x->xC[ix + i];
         ae_v_caddc(&y->xC[iy], 1, &a->xyC[ia + i][ja], 1, "N", m, v);
      }
      return;
   }
   if (opa == 2) {
   // y = A^H*x
      for (i = 0; i < m; i++) {
         y->xC[iy + i] = ae_complex_from_i(0);
      }
      for (i = 0; i < n; i++) {
         v = x->xC[ix + i];
         ae_v_caddc(&y->xC[iy], 1, &a->xyC[ia + i][ja], 1, "Conj", m, v);
      }
      return;
   }
}

// Scaled vector-matrix-vector addition: y = alpha a x + beta y.
// API: void rmatrixsymv(const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy);
void rmatrixsymv(ae_int_t n, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, double beta, RVector *y, ae_int_t iy) {
   ae_int_t i;
   ae_int_t j;
   double v;
   double vr;
   double vx;
// Quick exit for M=0, N=0 or Alpha=0.
//
// After this block we have M > 0, N > 0, Alpha != 0.
   if (n <= 0) {
      return;
   }
   if (alpha == 0.0) {
      if (beta != 0.0) {
         for (i = 0; i < n; i++) {
            y->xR[iy + i] *= beta;
         }
      } else {
         for (i = 0; i < n; i++) {
            y->xR[iy + i] = 0.0;
         }
      }
      return;
   }
// Try fast kernels
   if (n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel
      if (rmatrixsymvmkl(n, alpha, a, ia, ja, isupper, x, ix, beta, y, iy)) {
         return;
      }
   }
// Generic code
   if (beta != 0.0) {
      for (i = 0; i < n; i++) {
         y->xR[iy + i] *= beta;
      }
   } else {
      for (i = 0; i < n; i++) {
         y->xR[iy + i] = 0.0;
      }
   }
   if (isupper) {
   // Upper triangle of A is stored
      for (i = 0; i < n; i++) {
      // Process diagonal element
         v = alpha * a->xyR[ia + i][ja + i];
         y->xR[iy + i] += v * x->xR[ix + i];
      // Process off-diagonal elements
         vr = 0.0;
         vx = x->xR[ix + i];
         for (j = i + 1; j < n; j++) {
            v = alpha * a->xyR[ia + i][ja + j];
            y->xR[iy + j] += v * vx;
            vr += v * x->xR[ix + j];
         }
         y->xR[iy + i] += vr;
      }
   } else {
   // Lower triangle of A is stored
      for (i = 0; i < n; i++) {
      // Process diagonal element
         v = alpha * a->xyR[ia + i][ja + i];
         y->xR[iy + i] += v * x->xR[ix + i];
      // Process off-diagonal elements
         vr = 0.0;
         vx = x->xR[ix + i];
         for (j = 0; j < i; j++) {
            v = alpha * a->xyR[ia + i][ja + j];
            y->xR[iy + j] += v * vx;
            vr += v * x->xR[ix + j];
         }
         y->xR[iy + i] += vr;
      }
   }
}

// Vector-matrix-vector multiplication: x^T a x (with tmp = a x).
// API: double rmatrixsyvmv(const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const real_1d_array &x, const ae_int_t ix, const real_1d_array &tmp);
double rmatrixsyvmv(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, RVector *x, ae_int_t ix, RVector *tmp) {
   ae_int_t i;
   double result;
// Quick exit for N=0
   if (n <= 0) {
      result = 0.0;
      return result;
   }
// Generic code
   rmatrixsymv(n, 1.0, a, ia, ja, isupper, x, ix, 0.0, tmp, 0);
   result = 0.0;
   for (i = 0; i < n; i++) {
      result += x->xR[ix + i] * tmp->xR[i];
   }
   return result;
}

// This subroutine solves linear system op(A)*x=b where:
// * A is NxN upper/lower triangular/unitriangular matrix
// * X and B are Nx1 vectors
// * "op" may be identity transformation, transposition, conjugate transposition
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
//
// ALGLIB Routine: Copyright (c) 2017 by Sergey Bochkanov - converted to ALGLIB, remastered from LAPACK's DTRSV.
// Copyright (c) 2016 Reference BLAS level1 routine (LAPACK version 3.7.0)
//      Reference BLAS is a software package provided by Univ. of Tennessee,
//      Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
// API: void rmatrixtrsv(const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, const ae_int_t ix);
void rmatrixtrsv(ae_int_t n, RMatrix *a, ae_int_t ia, ae_int_t ja, bool isupper, bool isunit, ae_int_t optype, RVector *x, ae_int_t ix) {
   ae_int_t i;
   ae_int_t j;
   double v;
// Quick exit
   if (n <= 0) {
      return;
   }
// Try fast kernels
   if (n > ablas_blas2minvendorkernelsize) {
   // Try MKL kernel
      if (rmatrixtrsvmkl(n, a, ia, ja, isupper, isunit, optype, x, ix)) {
         return;
      }
   }
// Generic code
   if (optype == 0 && isupper) {
      for (i = n - 1; i >= 0; i--) {
         v = x->xR[ix + i];
         for (j = i + 1; j < n; j++) {
            v -= a->xyR[ia + i][ja + j] * x->xR[ix + j];
         }
         if (!isunit) {
            v /= a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
      }
      return;
   }
   if (optype == 0 && !isupper) {
      for (i = 0; i < n; i++) {
         v = x->xR[ix + i];
         for (j = 0; j < i; j++) {
            v -= a->xyR[ia + i][ja + j] * x->xR[ix + j];
         }
         if (!isunit) {
            v /= a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
      }
      return;
   }
   if (optype == 1 && isupper) {
      for (i = 0; i < n; i++) {
         v = x->xR[ix + i];
         if (!isunit) {
            v /= a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
         if (v == 0) {
            continue;
         }
         for (j = i + 1; j < n; j++) {
            x->xR[ix + j] -= v * a->xyR[ia + i][ja + j];
         }
      }
      return;
   }
   if (optype == 1 && !isupper) {
      for (i = n - 1; i >= 0; i--) {
         v = x->xR[ix + i];
         if (!isunit) {
            v /= a->xyR[ia + i][ja + i];
         }
         x->xR[ix + i] = v;
         if (v == 0) {
            continue;
         }
         for (j = 0; j < i; j++) {
            x->xR[ix + j] -= v * a->xyR[ia + i][ja + j];
         }
      }
      return;
   }
   ae_assert(false, "RMatrixTRSV: unexpected operation type");
}

// Level 2 subroutine
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
static void ablas_rmatrixrighttrsm2(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t i;
   ae_int_t j;
   double vr;
   double vd;
// Special case
   if (n * m == 0) {
      return;
   }
// Try to use "fast" code
   if (rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
      return;
   }
// General case
   if (isupper) {
   // Upper triangular matrix
      if (optype == 0) {
      // X*A^(-1)
         for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
               if (isunit) {
                  vd = 1.0;
               } else {
                  vd = a->xyR[i1 + j][j1 + j];
               }
               x->xyR[i2 + i][j2 + j] /= vd;
               if (j < n - 1) {
                  vr = x->xyR[i2 + i][j2 + j];
                  ae_v_subd(&x->xyR[i2 + i][j2 + j + 1], 1, &a->xyR[i1 + j][j1 + j + 1], 1, n - j - 1, vr);
               }
            }
         }
         return;
      }
      if (optype == 1) {
      // X*A^(-T)
         for (i = 0; i < m; i++) {
            for (j = n - 1; j >= 0; j--) {
               vr = 0.0;
               vd = 1.0;
               if (j < n - 1) {
                  vr = ae_v_dotproduct(&x->xyR[i2 + i][j2 + j + 1], 1, &a->xyR[i1 + j][j1 + j + 1], 1, n - j - 1);
               }
               if (!isunit) {
                  vd = a->xyR[i1 + j][j1 + j];
               }
               x->xyR[i2 + i][j2 + j] = (x->xyR[i2 + i][j2 + j] - vr) / vd;
            }
         }
         return;
      }
   } else {
   // Lower triangular matrix
      if (optype == 0) {
      // X*A^(-1)
         for (i = 0; i < m; i++) {
            for (j = n - 1; j >= 0; j--) {
               if (isunit) {
                  vd = 1.0;
               } else {
                  vd = a->xyR[i1 + j][j1 + j];
               }
               x->xyR[i2 + i][j2 + j] /= vd;
               if (j > 0) {
                  vr = x->xyR[i2 + i][j2 + j];
                  ae_v_subd(&x->xyR[i2 + i][j2], 1, &a->xyR[i1 + j][j1], 1, j, vr);
               }
            }
         }
         return;
      }
      if (optype == 1) {
      // X*A^(-T)
         for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
               vr = 0.0;
               vd = 1.0;
               if (j > 0) {
                  vr = ae_v_dotproduct(&x->xyR[i2 + i][j2], 1, &a->xyR[i1 + j][j1], 1, j);
               }
               if (!isunit) {
                  vd = a->xyR[i1 + j][j1 + j];
               }
               x->xyR[i2 + i][j2 + j] = (x->xyR[i2 + i][j2 + j] - vr) / vd;
            }
         }
         return;
      }
   }
}

// Level 2 variant of CMatrixRightTRSM
static void ablas_cmatrixrighttrsm2(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t i;
   ae_int_t j;
   ae_complex vc;
   ae_complex vd;
// Special case
   if (n * m == 0) {
      return;
   }
// Try to call fast TRSM
   if (cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
      return;
   }
// General case
   if (isupper) {
   // Upper triangular matrix
      if (optype == 0) {
      // X*A^(-1)
         for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
               if (isunit) {
                  vd = ae_complex_from_i(1);
               } else {
                  vd = a->xyC[i1 + j][j1 + j];
               }
               x->xyC[i2 + i][j2 + j] = ae_c_div(x->xyC[i2 + i][j2 + j], vd);
               if (j < n - 1) {
                  vc = x->xyC[i2 + i][j2 + j];
                  ae_v_csubc(&x->xyC[i2 + i][j2 + j + 1], 1, &a->xyC[i1 + j][j1 + j + 1], 1, "N", n - j - 1, vc);
               }
            }
         }
         return;
      }
      if (optype == 1) {
      // X*A^(-T)
         for (i = 0; i < m; i++) {
            for (j = n - 1; j >= 0; j--) {
               vc = ae_complex_from_i(0);
               vd = ae_complex_from_i(1);
               if (j < n - 1) {
                  vc = ae_v_cdotproduct(&x->xyC[i2 + i][j2 + j + 1], 1, "N", &a->xyC[i1 + j][j1 + j + 1], 1, "N", n - j - 1);
               }
               if (!isunit) {
                  vd = a->xyC[i1 + j][j1 + j];
               }
               x->xyC[i2 + i][j2 + j] = ae_c_div(ae_c_sub(x->xyC[i2 + i][j2 + j], vc), vd);
            }
         }
         return;
      }
      if (optype == 2) {
      // X*A^(-H)
         for (i = 0; i < m; i++) {
            for (j = n - 1; j >= 0; j--) {
               vc = ae_complex_from_i(0);
               vd = ae_complex_from_i(1);
               if (j < n - 1) {
                  vc = ae_v_cdotproduct(&x->xyC[i2 + i][j2 + j + 1], 1, "N", &a->xyC[i1 + j][j1 + j + 1], 1, "Conj", n - j - 1);
               }
               if (!isunit) {
                  vd = ae_c_conj(a->xyC[i1 + j][j1 + j]);
               }
               x->xyC[i2 + i][j2 + j] = ae_c_div(ae_c_sub(x->xyC[i2 + i][j2 + j], vc), vd);
            }
         }
         return;
      }
   } else {
   // Lower triangular matrix
      if (optype == 0) {
      // X*A^(-1)
         for (i = 0; i < m; i++) {
            for (j = n - 1; j >= 0; j--) {
               if (isunit) {
                  vd = ae_complex_from_i(1);
               } else {
                  vd = a->xyC[i1 + j][j1 + j];
               }
               x->xyC[i2 + i][j2 + j] = ae_c_div(x->xyC[i2 + i][j2 + j], vd);
               if (j > 0) {
                  vc = x->xyC[i2 + i][j2 + j];
                  ae_v_csubc(&x->xyC[i2 + i][j2], 1, &a->xyC[i1 + j][j1], 1, "N", j, vc);
               }
            }
         }
         return;
      }
      if (optype == 1) {
      // X*A^(-T)
         for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
               vc = ae_complex_from_i(0);
               vd = ae_complex_from_i(1);
               if (j > 0) {
                  vc = ae_v_cdotproduct(&x->xyC[i2 + i][j2], 1, "N", &a->xyC[i1 + j][j1], 1, "N", j);
               }
               if (!isunit) {
                  vd = a->xyC[i1 + j][j1 + j];
               }
               x->xyC[i2 + i][j2 + j] = ae_c_div(ae_c_sub(x->xyC[i2 + i][j2 + j], vc), vd);
            }
         }
         return;
      }
      if (optype == 2) {
      // X*A^(-H)
         for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
               vc = ae_complex_from_i(0);
               vd = ae_complex_from_i(1);
               if (j > 0) {
                  vc = ae_v_cdotproduct(&x->xyC[i2 + i][j2], 1, "N", &a->xyC[i1 + j][j1], 1, "Conj", j);
               }
               if (!isunit) {
                  vd = ae_c_conj(a->xyC[i1 + j][j1 + j]);
               }
               x->xyC[i2 + i][j2 + j] = ae_c_div(ae_c_sub(x->xyC[i2 + i][j2 + j], vc), vd);
            }
         }
         return;
      }
   }
}

// This subroutine calculates X*op(A^-1) where:
// * X is MxN general matrix
// * A is NxN upper/lower triangular/unitriangular matrix
// * "op" may be identity transformation, transposition
// Multiplication result replaces X.
//
// Inputs:
//     N   -   matrix size, N >= 0
//     M   -   matrix size, N >= 0
//     A       -   matrix, actial matrix is stored in A[I1:I1+N-1,J1:J1+N-1]
//     I1      -   submatrix offset
//     J1      -   submatrix offset
//     IsUpper -   whether matrix is upper triangular
//     IsUnit  -   whether matrix is unitriangular
//     OpType  -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//     X   -   matrix, actial matrix is stored in X[I2:I2+M-1,J2:J2+N-1]
//     I2  -   submatrix offset
//     J2  -   submatrix offset
//
// ALGLIB Routine: Copyright 15.12.2009-22.01.2018 by Sergey Bochkanov
// API: void rmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2);
void rmatrixrighttrsm(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (imax2(m, n) <= tsb) {
      tscur = tsa;
   }
   ae_assert(tscur >= 1, "RMatrixRightTRSM: integrity check failed");
// Upper level parallelization:
// * decide whether it is feasible to activate multithreading
// * perform optionally parallelized splits on M
// Parallelism was activated if: m >= 2 * tsb && (double)m * n * n >= smpactivationlevel()
   if (m >= 2 * tsb) {
   // Split X: X*A = (X1 X2)^T*A
      s1 = tiledsplit(m, tsb), s2 = m - s1;
      rmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      rmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, x, i2 + s1, j2);
      return;
   }
// Basecase: MKL or ALGLIB code
   if (imax2(m, n) <= tsb) {
      if (rmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
         return;
      }
   }
   if (imax2(m, n) <= tsa) {
      ablas_rmatrixrighttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      return;
   }
// Recursive subdivision
   if (m >= n) {
   // Split X: X*A = (X1 X2)^T*A
      s1 = tiledsplit(m, tscur), s2 = m - s1;
      rmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      rmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, x, i2 + s1, j2);
   } else {
   // Split A:
   //               (A1  A12)
   // X*op(A) = X*op(       )
   //               (     A2)
   //
   // Different variants depending on
   // IsUpper/OpType combinations
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      if (isupper && optype == 0) {
      //                  (A1  A12)-1
      // X*A^-1 = (X1 X2)*(       )
      //                  (     A2)
         rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         rmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1, j1 + s1, 0, 1.0, x, i2, j2 + s1);
         rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
      }
      if (isupper && optype != 0) {
      //                  (A1'     )-1
      // X*A^-1 = (X1 X2)*(        )
      //                  (A12' A2')
         rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
         rmatrixgemm(m, s1, s2, -1.0, x, i2, j2 + s1, 0, a, i1, j1 + s1, optype, 1.0, x, i2, j2);
         rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
      if (!isupper && optype == 0) {
      //                  (A1     )-1
      // X*A^-1 = (X1 X2)*(       )
      //                  (A21  A2)
         rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
         rmatrixgemm(m, s1, s2, -1.0, x, i2, j2 + s1, 0, a, i1 + s1, j1, 0, 1.0, x, i2, j2);
         rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
      if (!isupper && optype != 0) {
      //                  (A1' A21')-1
      // X*A^-1 = (X1 X2)*(        )
      //                  (     A2')
         rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         rmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1 + s1, j1, optype, 1.0, x, i2, j2 + s1);
         rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
      }
   }
}

// This subroutine calculates X*op(A^-1) where:
// * X is MxN general matrix
// * A is NxN upper/lower triangular/unitriangular matrix
// * "op" may be identity transformation, transposition, conjugate transposition
// Multiplication result replaces X.
//
// Inputs:
//     N   -   matrix size, N >= 0
//     M   -   matrix size, N >= 0
//     A       -   matrix, actial matrix is stored in A[I1:I1+N-1,J1:J1+N-1]
//     I1      -   submatrix offset
//     J1      -   submatrix offset
//     IsUpper -   whether matrix is upper triangular
//     IsUnit  -   whether matrix is unitriangular
//     OpType  -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//                 * 2 - conjugate transposition
//     X   -   matrix, actial matrix is stored in X[I2:I2+M-1,J2:J2+N-1]
//     I2  -   submatrix offset
//     J2  -   submatrix offset
//
// ALGLIB Routine: Copyright 20.01.2018 by Sergey Bochkanov
// API: void cmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2);
void cmatrixrighttrsm(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (imax2(m, n) <= tsb) {
      tscur = tsa;
   }
   ae_assert(tscur >= 1, "CMatrixRightTRSM: integrity check failed");
// Upper level parallelization:
// * decide whether it is feasible to activate multithreading
// * perform optionally parallelized splits on M
// Parallelism was activated if: m >= 2 * tsb && 4.0 * m * n * n >= smpactivationlevel()
   if (m >= 2 * tsb) {
   // Split X: X*A = (X1 X2)^T*A
      s1 = tiledsplit(m, tsb), s2 = m - s1;
      cmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      cmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, x, i2 + s1, j2);
      return;
   }
// Basecase: either MKL-supported code or ALGLIB basecase code
   if (imax2(m, n) <= tsb) {
      if (cmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
         return;
      }
   }
   if (imax2(m, n) <= tsa) {
      ablas_cmatrixrighttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      return;
   }
// Recursive subdivision
   if (m >= n) {
   // Split X: X*A = (X1 X2)^T*A
      s1 = tiledsplit(m, tscur), s2 = m - s1;
      cmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      cmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, x, i2 + s1, j2);
   } else {
   // Split A:
   //               (A1  A12)
   // X*op(A) = X*op(       )
   //               (     A2)
   //
   // Different variants depending on
   // IsUpper/OpType combinations
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      if (isupper && optype == 0) {
      //                  (A1  A12)-1
      // X*A^-1 = (X1 X2)*(       )
      //                  (     A2)
         cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         cmatrixgemm(m, s2, s1, ae_complex_from_d(-1.0), x, i2, j2, 0, a, i1, j1 + s1, 0, ae_complex_from_d(1.0), x, i2, j2 + s1);
         cmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
      }
      if (isupper && optype != 0) {
      //                  (A1'     )-1
      // X*A^-1 = (X1 X2)*(        )
      //                  (A12' A2')
         cmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
         cmatrixgemm(m, s1, s2, ae_complex_from_d(-1.0), x, i2, j2 + s1, 0, a, i1, j1 + s1, optype, ae_complex_from_d(1.0), x, i2, j2);
         cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
      if (!isupper && optype == 0) {
      //                  (A1     )-1
      // X*A^-1 = (X1 X2)*(       )
      //                  (A21  A2)
         cmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
         cmatrixgemm(m, s1, s2, ae_complex_from_d(-1.0), x, i2, j2 + s1, 0, a, i1 + s1, j1, 0, ae_complex_from_d(1.0), x, i2, j2);
         cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
      if (!isupper && optype != 0) {
      //                  (A1' A21')-1
      // X*A^-1 = (X1 X2)*(        )
      //                  (     A2')
         cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         cmatrixgemm(m, s2, s1, ae_complex_from_d(-1.0), x, i2, j2, 0, a, i1 + s1, j1, optype, ae_complex_from_d(1.0), x, i2, j2 + s1);
         cmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2, j2 + s1);
      }
   }
}

// Level 2 subroutine
static void ablas_rmatrixlefttrsm2(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t i;
   ae_int_t j;
   double vr;
   double vd;
// Special case
   if (n == 0 || m == 0) {
      return;
   }
// Try fast code
   if (rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
      return;
   }
// General case
   if (isupper) {
   // Upper triangular matrix
      if (optype == 0) {
      // A^(-1)*X
         for (i = m - 1; i >= 0; i--) {
            for (j = i + 1; j < m; j++) {
               vr = a->xyR[i1 + i][j1 + j];
               ae_v_subd(&x->xyR[i2 + i][j2], 1, &x->xyR[i2 + j][j2], 1, n, vr);
            }
            if (!isunit) {
               vd = 1 / a->xyR[i1 + i][j1 + i];
               ae_v_muld(&x->xyR[i2 + i][j2], 1, n, vd);
            }
         }
         return;
      }
      if (optype == 1) {
      // A^(-T)*X
         for (i = 0; i < m; i++) {
            if (isunit) {
               vd = 1.0;
            } else {
               vd = 1 / a->xyR[i1 + i][j1 + i];
            }
            ae_v_muld(&x->xyR[i2 + i][j2], 1, n, vd);
            for (j = i + 1; j < m; j++) {
               vr = a->xyR[i1 + i][j1 + j];
               ae_v_subd(&x->xyR[i2 + j][j2], 1, &x->xyR[i2 + i][j2], 1, n, vr);
            }
         }
         return;
      }
   } else {
   // Lower triangular matrix
      if (optype == 0) {
      // A^(-1)*X
         for (i = 0; i < m; i++) {
            for (j = 0; j < i; j++) {
               vr = a->xyR[i1 + i][j1 + j];
               ae_v_subd(&x->xyR[i2 + i][j2], 1, &x->xyR[i2 + j][j2], 1, n, vr);
            }
            if (isunit) {
               vd = 1.0;
            } else {
               vd = 1 / a->xyR[i1 + j][j1 + j];
            }
            ae_v_muld(&x->xyR[i2 + i][j2], 1, n, vd);
         }
         return;
      }
      if (optype == 1) {
      // A^(-T)*X
         for (i = m - 1; i >= 0; i--) {
            if (isunit) {
               vd = 1.0;
            } else {
               vd = 1 / a->xyR[i1 + i][j1 + i];
            }
            ae_v_muld(&x->xyR[i2 + i][j2], 1, n, vd);
            for (j = i - 1; j >= 0; j--) {
               vr = a->xyR[i1 + i][j1 + j];
               ae_v_subd(&x->xyR[i2 + j][j2], 1, &x->xyR[i2 + i][j2], 1, n, vr);
            }
         }
         return;
      }
   }
}

// Level-2 subroutine
static void ablas_cmatrixlefttrsm2(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t i;
   ae_int_t j;
   ae_complex vc;
   ae_complex vd;
// Special case
   if (n * m == 0) {
      return;
   }
// Try to call fast TRSM
   if (cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
      return;
   }
// General case
   if (isupper) {
   // Upper triangular matrix
      if (optype == 0) {
      // A^(-1)*X
         for (i = m - 1; i >= 0; i--) {
            for (j = i + 1; j < m; j++) {
               vc = a->xyC[i1 + i][j1 + j];
               ae_v_csubc(&x->xyC[i2 + i][j2], 1, &x->xyC[i2 + j][j2], 1, "N", n, vc);
            }
            if (!isunit) {
               vd = ae_c_d_div(1, a->xyC[i1 + i][j1 + i]);
               ae_v_cmulc(&x->xyC[i2 + i][j2], 1, n, vd);
            }
         }
         return;
      }
      if (optype == 1) {
      // A^(-T)*X
         for (i = 0; i < m; i++) {
            if (isunit) {
               vd = ae_complex_from_i(1);
            } else {
               vd = ae_c_d_div(1, a->xyC[i1 + i][j1 + i]);
            }
            ae_v_cmulc(&x->xyC[i2 + i][j2], 1, n, vd);
            for (j = i + 1; j < m; j++) {
               vc = a->xyC[i1 + i][j1 + j];
               ae_v_csubc(&x->xyC[i2 + j][j2], 1, &x->xyC[i2 + i][j2], 1, "N", n, vc);
            }
         }
         return;
      }
      if (optype == 2) {
      // A^(-H)*X
         for (i = 0; i < m; i++) {
            if (isunit) {
               vd = ae_complex_from_i(1);
            } else {
               vd = ae_c_d_div(1, ae_c_conj(a->xyC[i1 + i][j1 + i]));
            }
            ae_v_cmulc(&x->xyC[i2 + i][j2], 1, n, vd);
            for (j = i + 1; j < m; j++) {
               vc = ae_c_conj(a->xyC[i1 + i][j1 + j]);
               ae_v_csubc(&x->xyC[i2 + j][j2], 1, &x->xyC[i2 + i][j2], 1, "N", n, vc);
            }
         }
         return;
      }
   } else {
   // Lower triangular matrix
      if (optype == 0) {
      // A^(-1)*X
         for (i = 0; i < m; i++) {
            for (j = 0; j < i; j++) {
               vc = a->xyC[i1 + i][j1 + j];
               ae_v_csubc(&x->xyC[i2 + i][j2], 1, &x->xyC[i2 + j][j2], 1, "N", n, vc);
            }
            if (isunit) {
               vd = ae_complex_from_i(1);
            } else {
               vd = ae_c_d_div(1, a->xyC[i1 + j][j1 + j]);
            }
            ae_v_cmulc(&x->xyC[i2 + i][j2], 1, n, vd);
         }
         return;
      }
      if (optype == 1) {
      // A^(-T)*X
         for (i = m - 1; i >= 0; i--) {
            if (isunit) {
               vd = ae_complex_from_i(1);
            } else {
               vd = ae_c_d_div(1, a->xyC[i1 + i][j1 + i]);
            }
            ae_v_cmulc(&x->xyC[i2 + i][j2], 1, n, vd);
            for (j = i - 1; j >= 0; j--) {
               vc = a->xyC[i1 + i][j1 + j];
               ae_v_csubc(&x->xyC[i2 + j][j2], 1, &x->xyC[i2 + i][j2], 1, "N", n, vc);
            }
         }
         return;
      }
      if (optype == 2) {
      // A^(-H)*X
         for (i = m - 1; i >= 0; i--) {
            if (isunit) {
               vd = ae_complex_from_i(1);
            } else {
               vd = ae_c_d_div(1, ae_c_conj(a->xyC[i1 + i][j1 + i]));
            }
            ae_v_cmulc(&x->xyC[i2 + i][j2], 1, n, vd);
            for (j = i - 1; j >= 0; j--) {
               vc = ae_c_conj(a->xyC[i1 + i][j1 + j]);
               ae_v_csubc(&x->xyC[i2 + j][j2], 1, &x->xyC[i2 + i][j2], 1, "N", n, vc);
            }
         }
         return;
      }
   }
}

// This subroutine calculates op(A^-1)*X where:
// * X is MxN general matrix
// * A is MxM upper/lower triangular/unitriangular matrix
// * "op" may be identity transformation, transposition
// Multiplication result replaces X.
//
// Inputs:
//     N   -   matrix size, N >= 0
//     M   -   matrix size, N >= 0
//     A       -   matrix, actial matrix is stored in A[I1:I1+M-1,J1:J1+M-1]
//     I1      -   submatrix offset
//     J1      -   submatrix offset
//     IsUpper -   whether matrix is upper triangular
//     IsUnit  -   whether matrix is unitriangular
//     OpType  -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//     X   -   matrix, actial matrix is stored in X[I2:I2+M-1,J2:J2+N-1]
//     I2  -   submatrix offset
//     J2  -   submatrix offset
//
// ALGLIB Routine: Copyright 15.12.2009-22.01.2018 by Sergey Bochkanov
// API: void rmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2);
void rmatrixlefttrsm(ae_int_t m, ae_int_t n, RMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, RMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (imax2(m, n) <= tsb) {
      tscur = tsa;
   }
   ae_assert(tscur >= 1, "RMatrixLeftTRSMRec: integrity check failed");
// Upper level parallelization:
// * decide whether it is feasible to activate multithreading
// * perform optionally parallelized splits on N
// Parallelism was activated if: n >= 2 * tsb && (double)n * m * m >= smpactivationlevel()
   if (n >= 2 * tsb) {
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      rmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2 + s1);
      rmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      return;
   }
// Basecase: MKL or ALGLIB code
   if (imax2(m, n) <= tsb) {
      if (rmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
         return;
      }
   }
   if (imax2(m, n) <= tsa) {
      ablas_rmatrixlefttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      return;
   }
// Recursive subdivision
   if (n >= m) {
   // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      rmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      rmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2 + s1);
   } else {
   // Split A
      s1 = tiledsplit(m, tscur), s2 = m - s1;
      if (isupper && optype == 0) {
      //           (A1  A12)-1  ( X1 )
      // A^-1*X* = (       )   *(    )
      //           (     A2)    ( X2 )
         rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
         rmatrixgemm(s1, n, s2, -1.0, a, i1, j1 + s1, 0, x, i2 + s1, j2, 0, 1.0, x, i2, j2);
         rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
      if (isupper && optype != 0) {
      //          (A1'     )-1 ( X1 )
      // A^-1*X = (        )  *(    )
      //          (A12' A2')   ( X2 )
         rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         rmatrixgemm(s2, n, s1, -1.0, a, i1, j1 + s1, optype, x, i2, j2, 0, 1.0, x, i2 + s1, j2);
         rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
      }
      if (!isupper && optype == 0) {
      //          (A1     )-1 ( X1 )
      // A^-1*X = (       )  *(    )
      //          (A21  A2)   ( X2 )
         rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         rmatrixgemm(s2, n, s1, -1.0, a, i1 + s1, j1, 0, x, i2, j2, 0, 1.0, x, i2 + s1, j2);
         rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
      }
      if (!isupper && optype != 0) {
      //          (A1' A21')-1 ( X1 )
      // A^-1*X = (        )  *(    )
      //          (     A2')   ( X2 )
         rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
         rmatrixgemm(s1, n, s2, -1.0, a, i1 + s1, j1, optype, x, i2 + s1, j2, 0, 1.0, x, i2, j2);
         rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
   }
}

// This subroutine calculates op(A^-1)*X where:
// * X is MxN general matrix
// * A is MxM upper/lower triangular/unitriangular matrix
// * "op" may be identity transformation, transposition, conjugate transposition
// Multiplication result replaces X.
//
// Inputs:
//     N   -   matrix size, N >= 0
//     M   -   matrix size, N >= 0
//     A       -   matrix, actial matrix is stored in A[I1:I1+M-1,J1:J1+M-1]
//     I1      -   submatrix offset
//     J1      -   submatrix offset
//     IsUpper -   whether matrix is upper triangular
//     IsUnit  -   whether matrix is unitriangular
//     OpType  -   transformation type:
//                 * 0 - no transformation
//                 * 1 - transposition
//                 * 2 - conjugate transposition
//     X   -   matrix, actial matrix is stored in X[I2:I2+M-1,J2:J2+N-1]
//     I2  -   submatrix offset
//     J2  -   submatrix offset
//
// ALGLIB Routine: Copyright 15.12.2009-22.01.2018 by Sergey Bochkanov
// API: void cmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2);
void cmatrixlefttrsm(ae_int_t m, ae_int_t n, CMatrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, CMatrix *x, ae_int_t i2, ae_int_t j2) {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (imax2(m, n) <= tsb) {
      tscur = tsa;
   }
   ae_assert(tscur >= 1, "CMatrixLeftTRSM: integrity check failed");
// Upper level parallelization:
// * decide whether it is feasible to activate multithreading
// * perform optionally parallelized splits on N
// Parallelism was activated if: n >= 2 * tsb && 4.0 * n * m * m >= smpactivationlevel()
   if (n >= 2 * tsb) {
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      cmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2 + s1);
      cmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      return;
   }
// Basecase: either MKL-supported code or ALGLIB basecase code
   if (imax2(m, n) <= tsb) {
      if (cmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2)) {
         return;
      }
   }
   if (imax2(m, n) <= tsa) {
      ablas_cmatrixlefttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      return;
   }
// Recursive subdivision
   if (n >= m) {
   // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      cmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      cmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2 + s1);
   } else {
   // Split A
      s1 = tiledsplit(m, tscur), s2 = m - s1;
      if (isupper && optype == 0) {
      //           (A1  A12)-1  ( X1 )
      // A^-1*X* = (       )   *(    )
      //           (     A2)    ( X2 )
         cmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
         cmatrixgemm(s1, n, s2, ae_complex_from_d(-1.0), a, i1, j1 + s1, 0, x, i2 + s1, j2, 0, ae_complex_from_d(1.0), x, i2, j2);
         cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
      if (isupper && optype != 0) {
      //          (A1'     )-1 ( X1 )
      // A^-1*X = (        )  *(    )
      //          (A12' A2')   ( X2 )
         cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         cmatrixgemm(s2, n, s1, ae_complex_from_d(-1.0), a, i1, j1 + s1, optype, x, i2, j2, 0, ae_complex_from_d(1.0), x, i2 + s1, j2);
         cmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
      }
      if (!isupper && optype == 0) {
      //          (A1     )-1 ( X1 )
      // A^-1*X = (       )  *(    )
      //          (A21  A2)   ( X2 )
         cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
         cmatrixgemm(s2, n, s1, ae_complex_from_d(-1.0), a, i1 + s1, j1, 0, x, i2, j2, 0, ae_complex_from_d(1.0), x, i2 + s1, j2);
         cmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
      }
      if (!isupper && optype != 0) {
      //          (A1' A21')-1 ( X1 )
      // A^-1*X = (        )  *(    )
      //          (     A2')   ( X2 )
         cmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, x, i2 + s1, j2);
         cmatrixgemm(s1, n, s2, ae_complex_from_d(-1.0), a, i1 + s1, j1, optype, x, i2 + s1, j2, 0, ae_complex_from_d(1.0), x, i2, j2);
         cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
      }
   }
}

// Level 2 subrotuine
static void ablas_rmatrixsyrk2(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   double v;
// Fast exit (nothing to be done)
   if ((alpha == 0.0 || k == 0) && beta == 1.0) {
      return;
   }
// Try to call fast SYRK
   if (rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper)) {
      return;
   }
// SYRK
   if (optypea == 0) {
   // C=alpha*A*A^H+beta*C
      for (i = 0; i < n; i++) {
         if (isupper) {
            j1 = i;
            j2 = n - 1;
         } else {
            j1 = 0;
            j2 = i;
         }
         for (j = j1; j <= j2; j++) {
            if (alpha != 0.0 && k > 0) {
               v = ae_v_dotproduct(&a->xyR[ia + i][ja], 1, &a->xyR[ia + j][ja], 1, k);
            } else {
               v = 0.0;
            }
            if (beta == 0.0) {
               c->xyR[ic + i][jc + j] = alpha * v;
            } else {
               c->xyR[ic + i][jc + j] = beta * c->xyR[ic + i][jc + j] + alpha * v;
            }
         }
      }
      return;
   } else {
   // C=alpha*A^H*A+beta*C
      for (i = 0; i < n; i++) {
         if (isupper) {
            j1 = i;
            j2 = n - 1;
         } else {
            j1 = 0;
            j2 = i;
         }
         if (beta == 0.0) {
            for (j = j1; j <= j2; j++) {
               c->xyR[ic + i][jc + j] = 0.0;
            }
         } else {
            ae_v_muld(&c->xyR[ic + i][jc + j1], 1, j2 - j1 + 1, beta);
         }
      }
      if (alpha != 0.0 && k > 0) {
         for (i = 0; i < k; i++) {
            for (j = 0; j < n; j++) {
               if (isupper) {
                  j1 = j;
                  j2 = n - 1;
               } else {
                  j1 = 0;
                  j2 = j;
               }
               v = alpha * a->xyR[ia + i][ja + j];
               ae_v_addd(&c->xyR[ic + j][jc + j1], 1, &a->xyR[ia + i][ja + j1], 1, j2 - j1 + 1, v);
            }
         }
      }
      return;
   }
}

// Level 2 subroutine
static void ablas_cmatrixherk2(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   ae_complex v;
// Fast exit (nothing to be done)
   if ((alpha == 0.0 || k == 0) && beta == 1.0) {
      return;
   }
// Try to call fast SYRK
   if (cmatrixherkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper)) {
      return;
   }
// SYRK
   if (optypea == 0) {
   // C=alpha*A*A^H+beta*C
      for (i = 0; i < n; i++) {
         if (isupper) {
            j1 = i;
            j2 = n - 1;
         } else {
            j1 = 0;
            j2 = i;
         }
         for (j = j1; j <= j2; j++) {
            if (alpha != 0.0 && k > 0) {
               v = ae_v_cdotproduct(&a->xyC[ia + i][ja], 1, "N", &a->xyC[ia + j][ja], 1, "Conj", k);
            } else {
               v = ae_complex_from_i(0);
            }
            if (beta == 0.0) {
               c->xyC[ic + i][jc + j] = ae_c_mul_d(v, alpha);
            } else {
               c->xyC[ic + i][jc + j] = ae_c_add(ae_c_mul_d(c->xyC[ic + i][jc + j], beta), ae_c_mul_d(v, alpha));
            }
         }
      }
      return;
   } else {
   // C=alpha*A^H*A+beta*C
      for (i = 0; i < n; i++) {
         if (isupper) {
            j1 = i;
            j2 = n - 1;
         } else {
            j1 = 0;
            j2 = i;
         }
         if (beta == 0.0) {
            for (j = j1; j <= j2; j++) {
               c->xyC[ic + i][jc + j] = ae_complex_from_i(0);
            }
         } else {
            ae_v_cmuld(&c->xyC[ic + i][jc + j1], 1, j2 - j1 + 1, beta);
         }
      }
      if (alpha != 0.0 && k > 0) {
         for (i = 0; i < k; i++) {
            for (j = 0; j < n; j++) {
               if (isupper) {
                  j1 = j;
                  j2 = n - 1;
               } else {
                  j1 = 0;
                  j2 = j;
               }
               v = ae_c_mul_d(ae_c_conj(a->xyC[ia + i][ja + j]), alpha);
               ae_v_caddc(&c->xyC[ic + j][jc + j1], 1, &a->xyC[ia + i][ja + j1], 1, "N", j2 - j1 + 1, v);
            }
         }
      }
      return;
   }
}

// This subroutine calculates  C=alpha*A*A^T+beta*C  or  C=alpha*A^T*A+beta*C
// where:
// * C is NxN symmetric matrix given by its upper/lower triangle
// * A is NxK matrix when A*A^T is calculated, KxN matrix otherwise
//
// Additional info:
// * multiplication result replaces C. If Beta=0, C elements are not used in
//   calculations (not multiplied by zero - just not referenced)
// * if Alpha=0, A is not used (not multiplied by zero - just not referenced)
// * if both Beta and Alpha are zero, C is filled by zeros.
//
// Inputs:
//     N       -   matrix size, N >= 0
//     K       -   matrix size, K >= 0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset (row index)
//     JA      -   submatrix offset (column index)
//     OpTypeA -   multiplication type:
//                 * 0 - A*A^T is calculated
//                 * 2 - A^T*A is calculated
//     Beta    -   coefficient
//     C       -   preallocated input/output matrix
//     IC      -   submatrix offset (row index)
//     JC      -   submatrix offset (column index)
//     IsUpper -   whether C is upper triangular or lower triangular
//
// ALGLIB Routine: Copyright 16.12.2009-22.01.2018 by Sergey Bochkanov
// API: void rmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
void rmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (imax2(n, k) <= tsb) {
      tscur = tsa;
   }
   ae_assert(tscur >= 1, "RMatrixSYRK: integrity check failed");
// Decide whether it is feasible to activate multithreading
// Parallelism was activated if: n >= 2 * tsb && 2.0 * k * n * n / 2 >= smpactivationlevel()
// Use MKL or generic basecase code
   if (imax2(n, k) <= tsb) {
      if (rmatrixsyrkmkl(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper)) {
         return;
      }
   }
   if (imax2(n, k) <= tsa) {
      ablas_rmatrixsyrk2(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
      return;
   }
// Recursive subdivision of the problem
   if (k >= n) {
   // Split K
      s1 = tiledsplit(k, tscur), s2 = k - s1;
      if (optypea == 0) {
         rmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         rmatrixsyrk(n, s2, alpha, a, ia, ja + s1, optypea, 1.0, c, ic, jc, isupper);
      } else {
         rmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         rmatrixsyrk(n, s2, alpha, a, ia + s1, ja, optypea, 1.0, c, ic, jc, isupper);
      }
   } else {
   // Split N
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      if (optypea == 0 && isupper) {
         rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         rmatrixsyrk(s2, k, alpha, a, ia + s1, ja, optypea, beta, c, ic + s1, jc + s1, isupper);
         rmatrixgemm(s1, s2, k, alpha, a, ia, ja, 0, a, ia + s1, ja, 1, beta, c, ic, jc + s1);
      }
      if (optypea == 0 && !isupper) {
         rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         rmatrixsyrk(s2, k, alpha, a, ia + s1, ja, optypea, beta, c, ic + s1, jc + s1, isupper);
         rmatrixgemm(s2, s1, k, alpha, a, ia + s1, ja, 0, a, ia, ja, 1, beta, c, ic + s1, jc);
      }
      if (optypea != 0 && isupper) {
         rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         rmatrixsyrk(s2, k, alpha, a, ia, ja + s1, optypea, beta, c, ic + s1, jc + s1, isupper);
         rmatrixgemm(s1, s2, k, alpha, a, ia, ja, 1, a, ia, ja + s1, 0, beta, c, ic, jc + s1);
      }
      if (optypea != 0 && !isupper) {
         rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         rmatrixsyrk(s2, k, alpha, a, ia, ja + s1, optypea, beta, c, ic + s1, jc + s1, isupper);
         rmatrixgemm(s2, s1, k, alpha, a, ia, ja + s1, 1, a, ia, ja, 0, beta, c, ic + s1, jc);
      }
   }
}

// This subroutine calculates  C=alpha*A*A^H+beta*C  or  C=alpha*A^H*A+beta*C
// where:
// * C is NxN Hermitian matrix given by its upper/lower triangle
// * A is NxK matrix when A*A^H is calculated, KxN matrix otherwise
//
// Additional info:
// * multiplication result replaces C. If Beta=0, C elements are not used in
//   calculations (not multiplied by zero - just not referenced)
// * if Alpha=0, A is not used (not multiplied by zero - just not referenced)
// * if both Beta and Alpha are zero, C is filled by zeros.
//
// Inputs:
//     N       -   matrix size, N >= 0
//     K       -   matrix size, K >= 0
//     Alpha   -   coefficient
//     A       -   matrix
//     IA      -   submatrix offset (row index)
//     JA      -   submatrix offset (column index)
//     OpTypeA -   multiplication type:
//                 * 0 - A*A^H is calculated
//                 * 2 - A^H*A is calculated
//     Beta    -   coefficient
//     C       -   preallocated input/output matrix
//     IC      -   submatrix offset (row index)
//     JC      -   submatrix offset (column index)
//     IsUpper -   whether upper or lower triangle of C is updated;
//                 this function updates only one half of C, leaving
//                 other half unchanged (not referenced at all).
//
// ALGLIB Routine: Copyright 16.12.2009-22.01.2018 by Sergey Bochkanov
// API: void cmatrixherk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
void cmatrixherk(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (imax2(n, k) <= tsb) {
      tscur = tsa;
   }
   ae_assert(tscur >= 1, "CMatrixHERK: integrity check failed");
// Decide whether it is feasible to activate multithreading
// Parallelism was activated if: n >= 2 * tsb && 8.0 * k * n * n / 2 >= smpactivationlevel()
// Use MKL or ALGLIB basecase code
   if (imax2(n, k) <= tsb) {
      if (cmatrixherkmkl(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper)) {
         return;
      }
   }
   if (imax2(n, k) <= tsa) {
      ablas_cmatrixherk2(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
      return;
   }
// Recursive division of the problem
   if (k >= n) {
   // Split K
      s1 = tiledsplit(k, tscur), s2 = k - s1;
      if (optypea == 0) {
         cmatrixherk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         cmatrixherk(n, s2, alpha, a, ia, ja + s1, optypea, 1.0, c, ic, jc, isupper);
      } else {
         cmatrixherk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         cmatrixherk(n, s2, alpha, a, ia + s1, ja, optypea, 1.0, c, ic, jc, isupper);
      }
   } else {
   // Split N
      s1 = tiledsplit(n, tscur), s2 = n - s1;
      if (optypea == 0 && isupper) {
         cmatrixherk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         cmatrixherk(s2, k, alpha, a, ia + s1, ja, optypea, beta, c, ic + s1, jc + s1, isupper);
         cmatrixgemm(s1, s2, k, ae_complex_from_d(alpha), a, ia, ja, 0, a, ia + s1, ja, 2, ae_complex_from_d(beta), c, ic, jc + s1);
      }
      if (optypea == 0 && !isupper) {
         cmatrixherk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         cmatrixherk(s2, k, alpha, a, ia + s1, ja, optypea, beta, c, ic + s1, jc + s1, isupper);
         cmatrixgemm(s2, s1, k, ae_complex_from_d(alpha), a, ia + s1, ja, 0, a, ia, ja, 2, ae_complex_from_d(beta), c, ic + s1, jc);
      }
      if (optypea != 0 && isupper) {
         cmatrixherk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         cmatrixherk(s2, k, alpha, a, ia, ja + s1, optypea, beta, c, ic + s1, jc + s1, isupper);
         cmatrixgemm(s1, s2, k, ae_complex_from_d(alpha), a, ia, ja, 2, a, ia, ja + s1, 0, ae_complex_from_d(beta), c, ic, jc + s1);
      }
      if (optypea != 0 && !isupper) {
         cmatrixherk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
         cmatrixherk(s2, k, alpha, a, ia, ja + s1, optypea, beta, c, ic + s1, jc + s1, isupper);
         cmatrixgemm(s2, s1, k, ae_complex_from_d(alpha), a, ia, ja + s1, 2, a, ia, ja, 0, ae_complex_from_d(beta), c, ic + s1, jc);
      }
   }
}

// This subroutine is an older version of CMatrixHERK(), one with wrong  name
// (it is HErmitian update, not SYmmetric). It  is  left  here  for  backward
// compatibility.
//
// ALGLIB Routine: Copyright 16.12.2009 by Sergey Bochkanov
// API: void cmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
void cmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, CMatrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   cmatrixherk(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
}

// This subroutine is an actual implementation of RMatrixGEMM.  It  does  not
// perform some integrity checks performed in the  driver  function,  and  it
// does not activate multithreading  framework  (driver  decides  whether  to
// activate workers or not).
//
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
static void ablas_rmatrixgemmrec(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_int_t mnk = imax3(m, n, k);
   ae_int_t tsa = matrixtilesizea(), tsb = matrixtilesizeb(), tscur = mnk <= tsb? tsa: tsb;
   ae_assert(tscur >= 1, "RMatrixGEMMRec: integrity check failed");
// Use MKL or ALGLIB basecase code
   if (mnk <= tsb && rmatrixgemmmkl(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc)) return;
   else if (mnk <= tsa) rmatrixgemmk(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
// Recursive algorithm: split on M or N
   else if (mnk <= m) { // A*B = (A1 A2)^T*B
      ae_int_t m0 = tiledsplit(m, tscur), m1 = m - m0;
      ablas_rmatrixgemmrec(m0, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
      ic += m0;
      if (optypea == 0) ia += m0; else ja += m0;
      ablas_rmatrixgemmrec(m1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
   } else if (mnk <= n) { // A*B = A*(B1 B2)
      ae_int_t n0 = tiledsplit(n, tscur), n1 = n - n0;
      ablas_rmatrixgemmrec(m, n0, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
      jc += n0;
      if (optypeb == 0) jb += n0; else ib += n0;
      ablas_rmatrixgemmrec(m, n1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
   } else { // Recursive algorithm: split on K: A*B = (A1 A2)*(B1 B2)^T
      ae_int_t k0 = tiledsplit(k, tscur), k1 = k - k0;
      ablas_rmatrixgemmrec(m, n, k0, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
      if (optypea == 0) ja += k0; else ia += k0;
      if (optypeb == 0) ib += k0; else jb += k0;
      ablas_rmatrixgemmrec(m, n, k1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, 1.0, c, ic, jc);
   }
}

// This subroutine is an actual implementation of CMatrixGEMM.  It  does  not
// perform some integrity checks performed in the  driver  function,  and  it
// does not activate multithreading  framework  (driver  decides  whether  to
// activate workers or not).
//
// ALGLIB Routine: Copyright 10.01.2019 by Sergey Bochkanov
static void ablas_cmatrixgemmrec(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_int_t mnk = imax3(m, n, k);
// Tile hierarchy: B -> A -> A/2
   ae_int_t tsa = matrixtilesizea() / 2, tsb = matrixtilesizeb(), tscur = mnk <= tsb? tsa: tsb;
   ae_assert(tscur >= 1, "CMatrixGEMMRec: integrity check failed");
// Use MKL or ALGLIB basecase code
   if (mnk <= tsb && cmatrixgemmmkl(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc)) return;
   else if (mnk <= tsa) cmatrixgemmk(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
// Recursive algorithm: parallel splitting on M/N
   else if (mnk <= m) { // A*B = (A1 A2)^T*B
      ae_int_t m0 = tiledsplit(m, tscur), m1 = m - m0;
      ablas_cmatrixgemmrec(m0, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
      ic += m0;
      if (optypea == 0) ia += m0; else ja += m0;
      ablas_cmatrixgemmrec(m1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
   } else if (mnk <= n) { // A*B = A*(B1 B2)
      ae_int_t n0 = tiledsplit(n, tscur), n1 = n - n0;
      ablas_cmatrixgemmrec(m, n0, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
      jc += n0;
      if (optypeb == 0) jb += n0; else ib += n0;
      ablas_cmatrixgemmrec(m, n1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
   } else { // Recursive algorithm: serial splitting on K: A*B = (A1 A2)*(B1 B2)^T
      ae_int_t k0 = tiledsplit(k, tscur), k1 = k - k0;
      ablas_cmatrixgemmrec(m, n, k0, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
      if (optypea == 0) ja += k0; else ia += k0;
      if (optypeb == 0) ib += k0; else jb += k0;
      ablas_cmatrixgemmrec(m, n, k1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, ae_complex_from_d(1.0), c, ic, jc);
   }
}

// This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
// * C is MxN general matrix
// * op1(A) is MxK matrix
// * op2(B) is KxN matrix
// * "op" may be identity transformation, transposition
//
// Additional info:
// * cache-oblivious algorithm is used.
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
//     M       -   matrix size, M > 0
//     N       -   matrix size, N > 0
//     K       -   matrix size, K > 0
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
//     C       -   PREALLOCATED output matrix, large enough to store result
//     IC      -   submatrix offset
//     JC      -   submatrix offset
//
// ALGLIB Routine: Copyright 2009-2019 by Sergey Bochkanov
// API: void rmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc);
void rmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_int_t ts;
   ts = matrixtilesizeb();
// Check input sizes for correctness
   ae_assert(optypea == 0 || optypea == 1, "RMatrixGEMM: incorrect OpTypeA (must be 0 or 1)");
   ae_assert(optypeb == 0 || optypeb == 1, "RMatrixGEMM: incorrect OpTypeB (must be 0 or 1)");
   ae_assert(ic + m <= c->rows, "RMatrixGEMM: incorect size of output matrix C");
   ae_assert(jc + n <= c->cols, "RMatrixGEMM: incorect size of output matrix C");
// Decide whether it is feasible to activate multithreading
// Parallelism was activated if: (m >= 2 * ts || n >= 2 * ts) && 2.0 * m * n * k >= smpactivationlevel()
// Start actual work
   ablas_rmatrixgemmrec(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
}

// This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
// * C is MxN general matrix
// * op1(A) is MxK matrix
// * op2(B) is KxN matrix
// * "op" may be identity transformation, transposition, conjugate transposition
//
// Additional info:
// * cache-oblivious algorithm is used.
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
//     M       -   matrix size, M > 0
//     N       -   matrix size, N > 0
//     K       -   matrix size, K > 0
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
//     C       -   matrix (PREALLOCATED, large enough to store result)
//     IC      -   submatrix offset
//     JC      -   submatrix offset
//
// ALGLIB Routine: Copyright 2009-2019 by Sergey Bochkanov
// API: void cmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const complex alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const complex beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc);
void cmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, CMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, CMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, CMatrix *c, ae_int_t ic, ae_int_t jc) {
   ae_int_t ts;
   ts = matrixtilesizeb();
// Check input sizes for correctness
   ae_assert(optypea == 0 || optypea == 1 || optypea == 2, "CMatrixGEMM: incorrect OpTypeA (must be 0 or 1 or 2)");
   ae_assert(optypeb == 0 || optypeb == 1 || optypeb == 2, "CMatrixGEMM: incorrect OpTypeB (must be 0 or 1 or 2)");
   ae_assert(ic + m <= c->rows, "CMatrixGEMM: incorect size of output matrix C");
   ae_assert(jc + n <= c->cols, "CMatrixGEMM: incorect size of output matrix C");
// Decide whether it is feasible to activate multithreading
// Parallelism was activated if: (m >= 2 * ts || n >= 2 * ts) && 8.0 * m * n * k >= smpactivationlevel()
// Start actual work
   ablas_cmatrixgemmrec(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
}
} // end of namespace alglib_impl

namespace alglib {
void rmatrixtranspose(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_2d_array &b, const ae_int_t ib, const ae_int_t jb) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixtranspose(m, n, ConstT(ae_matrix, a), ia, ja, ConstT(ae_matrix, b), ib, jb);
   alglib_impl::ae_state_clear();
}

void cmatrixtranspose(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixtranspose(m, n, ConstT(ae_matrix, a), ia, ja, ConstT(ae_matrix, b), ib, jb);
   alglib_impl::ae_state_clear();
}

void rmatrixenforcesymmetricity(const real_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixenforcesymmetricity(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
}

void rvectorcopy(const ae_int_t n, const real_1d_array &a, const ae_int_t ia, const real_1d_array &b, const ae_int_t ib) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rvectorcopy(n, ConstT(ae_vector, a), ia, ConstT(ae_vector, b), ib);
   alglib_impl::ae_state_clear();
}

void rmatrixcopy(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_2d_array &b, const ae_int_t ib, const ae_int_t jb) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixcopy(m, n, ConstT(ae_matrix, a), ia, ja, ConstT(ae_matrix, b), ib, jb);
   alglib_impl::ae_state_clear();
}

void cmatrixcopy(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixcopy(m, n, ConstT(ae_matrix, a), ia, ja, ConstT(ae_matrix, b), ib, jb);
   alglib_impl::ae_state_clear();
}

void rmatrixgencopy(const ae_int_t m, const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const double beta, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixgencopy(m, n, alpha, ConstT(ae_matrix, a), ia, ja, beta, ConstT(ae_matrix, b), ib, jb);
   alglib_impl::ae_state_clear();
}

void rmatrixger(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const double alpha, const real_1d_array &u, const ae_int_t iu, const real_1d_array &v, const ae_int_t iv) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixger(m, n, ConstT(ae_matrix, a), ia, ja, alpha, ConstT(ae_vector, u), iu, ConstT(ae_vector, v), iv);
   alglib_impl::ae_state_clear();
}

void rmatrixrank1(const ae_int_t m, const ae_int_t n, real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_1d_array &u, const ae_int_t iu, real_1d_array &v, const ae_int_t iv) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixrank1(m, n, ConstT(ae_matrix, a), ia, ja, ConstT(ae_vector, u), iu, ConstT(ae_vector, v), iv);
   alglib_impl::ae_state_clear();
}

void cmatrixrank1(const ae_int_t m, const ae_int_t n, complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_1d_array &u, const ae_int_t iu, complex_1d_array &v, const ae_int_t iv) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixrank1(m, n, ConstT(ae_matrix, a), ia, ja, ConstT(ae_vector, u), iu, ConstT(ae_vector, v), iv);
   alglib_impl::ae_state_clear();
}

void rmatrixgemv(const ae_int_t m, const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixgemv(m, n, alpha, ConstT(ae_matrix, a), ia, ja, opa, ConstT(ae_vector, x), ix, beta, ConstT(ae_vector, y), iy);
   alglib_impl::ae_state_clear();
}

void rmatrixmv(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, real_1d_array &y, const ae_int_t iy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixmv(m, n, ConstT(ae_matrix, a), ia, ja, opa, ConstT(ae_vector, x), ix, ConstT(ae_vector, y), iy);
   alglib_impl::ae_state_clear();
}

void cmatrixmv(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const complex_1d_array &x, const ae_int_t ix, complex_1d_array &y, const ae_int_t iy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixmv(m, n, ConstT(ae_matrix, a), ia, ja, opa, ConstT(ae_vector, x), ix, ConstT(ae_vector, y), iy);
   alglib_impl::ae_state_clear();
}

void rmatrixsymv(const ae_int_t n, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const real_1d_array &x, const ae_int_t ix, const double beta, const real_1d_array &y, const ae_int_t iy) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsymv(n, alpha, ConstT(ae_matrix, a), ia, ja, isupper, ConstT(ae_vector, x), ix, beta, ConstT(ae_vector, y), iy);
   alglib_impl::ae_state_clear();
}

double rmatrixsyvmv(const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const real_1d_array &x, const ae_int_t ix, const real_1d_array &tmp) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixsyvmv(n, ConstT(ae_matrix, a), ia, ja, isupper, ConstT(ae_vector, x), ix, ConstT(ae_vector, tmp));
   alglib_impl::ae_state_clear();
   return D;
}

void rmatrixtrsv(const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, const ae_int_t ix) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixtrsv(n, ConstT(ae_matrix, a), ia, ja, isupper, isunit, optype, ConstT(ae_vector, x), ix);
   alglib_impl::ae_state_clear();
}

void rmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixrighttrsm(m, n, ConstT(ae_matrix, a), i1, j1, isupper, isunit, optype, ConstT(ae_matrix, x), i2, j2);
   alglib_impl::ae_state_clear();
}

void cmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixrighttrsm(m, n, ConstT(ae_matrix, a), i1, j1, isupper, isunit, optype, ConstT(ae_matrix, x), i2, j2);
   alglib_impl::ae_state_clear();
}

void rmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlefttrsm(m, n, ConstT(ae_matrix, a), i1, j1, isupper, isunit, optype, ConstT(ae_matrix, x), i2, j2);
   alglib_impl::ae_state_clear();
}

void cmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlefttrsm(m, n, ConstT(ae_matrix, a), i1, j1, isupper, isunit, optype, ConstT(ae_matrix, x), i2, j2);
   alglib_impl::ae_state_clear();
}

void rmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixsyrk(n, k, alpha, ConstT(ae_matrix, a), ia, ja, optypea, beta, ConstT(ae_matrix, c), ic, jc, isupper);
   alglib_impl::ae_state_clear();
}

void cmatrixherk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixherk(n, k, alpha, ConstT(ae_matrix, a), ia, ja, optypea, beta, ConstT(ae_matrix, c), ic, jc, isupper);
   alglib_impl::ae_state_clear();
}

void cmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixsyrk(n, k, alpha, ConstT(ae_matrix, a), ia, ja, optypea, beta, ConstT(ae_matrix, c), ic, jc, isupper);
   alglib_impl::ae_state_clear();
}

void rmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixgemm(m, n, k, alpha, ConstT(ae_matrix, a), ia, ja, optypea, ConstT(ae_matrix, b), ib, jb, optypeb, beta, ConstT(ae_matrix, c), ic, jc);
   alglib_impl::ae_state_clear();
}

void cmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const complex alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const complex beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixgemm(m, n, k, *alpha.c_ptr(), ConstT(ae_matrix, a), ia, ja, optypea, ConstT(ae_matrix, b), ib, jb, optypeb, *beta.c_ptr(), ConstT(ae_matrix, c), ic, jc);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === DLU Package ===
// Depends on: ABLAS
namespace alglib_impl {
// Real LUP kernel
//
// ALGLIB Routine: Copyright 10.01.2010 by Sergey Bochkanov
static void dlu_rmatrixlup2(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, RVector *tmp) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t jp;
   double s;
   if (m == 0 || n == 0) {
      return;
   }
   for (j = 0; j < imin2(m, n); j++) {
      jp = j;
      for (i = j + 1; i < n; i++) {
         if (fabs(a->xyR[offs + j][offs + i]) > fabs(a->xyR[offs + j][offs + jp])) {
            jp = i;
         }
      }
      pivots->xZ[offs + j] = offs + jp;
      if (jp != j) {
         ae_v_move(tmp->xR, 1, &a->xyR[offs][offs + j], a->stride, m);
         ae_v_move(&a->xyR[offs][offs + j], a->stride, &a->xyR[offs][offs + jp], a->stride, m);
         ae_v_move(&a->xyR[offs][offs + jp], a->stride, tmp->xR, 1, m);
      }
      if (a->xyR[offs + j][offs + j] != 0.0 && j + 1 < n) {
         s = 1 / a->xyR[offs + j][offs + j];
         ae_v_muld(&a->xyR[offs + j][offs + j + 1], 1, n - j - 1, s);
      }
      if (j < imin2(m - 1, n - 1)) {
         ae_v_move(tmp->xR, 1, &a->xyR[offs + j + 1][offs + j], a->stride, m - j - 1);
         ae_v_moveneg(&tmp->xR[m], 1, &a->xyR[offs + j][offs + j + 1], 1, n - j - 1);
         rmatrixrank1(m - j - 1, n - j - 1, a, offs + j + 1, offs + j + 1, tmp, 0, tmp, m);
      }
   }
}

// Complex LUP kernel
//
// ALGLIB Routine: Copyright 10.01.2010 by Sergey Bochkanov
static void dlu_cmatrixlup2(CMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, CVector *tmp) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t jp;
   ae_complex s;
   if (m == 0 || n == 0) {
      return;
   }
   for (j = 0; j < imin2(m, n); j++) {
      jp = j;
      for (i = j + 1; i < n; i++) {
         if (ae_c_abs(a->xyC[offs + j][offs + i]) > ae_c_abs(a->xyC[offs + j][offs + jp])) {
            jp = i;
         }
      }
      pivots->xZ[offs + j] = offs + jp;
      if (jp != j) {
         ae_v_cmove(tmp->xC, 1, &a->xyC[offs][offs + j], a->stride, "N", m);
         ae_v_cmove(&a->xyC[offs][offs + j], a->stride, &a->xyC[offs][offs + jp], a->stride, "N", m);
         ae_v_cmove(&a->xyC[offs][offs + jp], a->stride, tmp->xC, 1, "N", m);
      }
      if (ae_c_neq_d(a->xyC[offs + j][offs + j], 0.0) && j + 1 < n) {
         s = ae_c_d_div(1, a->xyC[offs + j][offs + j]);
         ae_v_cmulc(&a->xyC[offs + j][offs + j + 1], 1, n - j - 1, s);
      }
      if (j < imin2(m - 1, n - 1)) {
         ae_v_cmove(tmp->xC, 1, &a->xyC[offs + j + 1][offs + j], a->stride, "N", m - j - 1);
         ae_v_cmoveneg(&tmp->xC[m], 1, &a->xyC[offs + j][offs + j + 1], 1, "N", n - j - 1);
         cmatrixrank1(m - j - 1, n - j - 1, a, offs + j + 1, offs + j + 1, tmp, 0, tmp, m);
      }
   }
}

// Recurrent real LU subroutine.
// Never call it directly.
//
// ALGLIB Routine: Copyright 04.01.2010 by Sergey Bochkanov
void rmatrixluprec(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, RVector *tmp) {
   ae_int_t i;
   ae_int_t m1;
   ae_int_t m2;
   if (imin2(m, n) <= ablasblocksize(a)) {
      dlu_rmatrixlup2(a, offs, m, n, pivots, tmp);
      return;
   }
   if (m > n) {
      rmatrixluprec(a, offs, n, n, pivots, tmp);
      for (i = 0; i < n; i++) {
         if (offs + i != pivots->xZ[offs + i]) {
            ae_v_move(tmp->xR, 1, &a->xyR[offs + n][offs + i], a->stride, m - n);
            ae_v_move(&a->xyR[offs + n][offs + i], a->stride, &a->xyR[offs + n][pivots->xZ[offs + i]], a->stride, m - n);
            ae_v_move(&a->xyR[offs + n][pivots->xZ[offs + i]], a->stride, tmp->xR, 1, m - n);
         }
      }
      rmatrixrighttrsm(m - n, n, a, offs, offs, true, true, 0, a, offs + n, offs);
      return;
   }
   m1 = ablassplitlength(a, m), m2 = m - m1;
   rmatrixluprec(a, offs, m1, n, pivots, tmp);
   if (m2 > 0) {
      for (i = 0; i < m1; i++) {
         if (offs + i != pivots->xZ[offs + i]) {
            ae_v_move(tmp->xR, 1, &a->xyR[offs + m1][offs + i], a->stride, m2);
            ae_v_move(&a->xyR[offs + m1][offs + i], a->stride, &a->xyR[offs + m1][pivots->xZ[offs + i]], a->stride, m - m1);
            ae_v_move(&a->xyR[offs + m1][pivots->xZ[offs + i]], a->stride, tmp->xR, 1, m - m1);
         }
      }
      rmatrixrighttrsm(m2, m1, a, offs, offs, true, true, 0, a, offs + m1, offs);
      rmatrixgemm(m - m1, n - m1, m1, -1.0, a, offs + m1, offs, 0, a, offs, offs + m1, 0, 1.0, a, offs + m1, offs + m1);
      rmatrixluprec(a, offs + m1, m - m1, n - m1, pivots, tmp);
      for (i = 0; i < m2; i++) {
         if (offs + m1 + i != pivots->xZ[offs + m1 + i]) {
            ae_v_move(tmp->xR, 1, &a->xyR[offs][offs + m1 + i], a->stride, m1);
            ae_v_move(&a->xyR[offs][offs + m1 + i], a->stride, &a->xyR[offs][pivots->xZ[offs + m1 + i]], a->stride, m1);
            ae_v_move(&a->xyR[offs][pivots->xZ[offs + m1 + i]], a->stride, tmp->xR, 1, m1);
         }
      }
   }
}

// Recurrent complex LU subroutine.
// Never call it directly.
//
// ALGLIB Routine: Copyright 04.01.2010 by Sergey Bochkanov
void cmatrixluprec(CMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, CVector *tmp) {
   ae_int_t i;
   ae_int_t m1;
   ae_int_t m2;
   if (imin2(m, n) <= ablascomplexblocksize(a)) {
      dlu_cmatrixlup2(a, offs, m, n, pivots, tmp);
      return;
   }
   if (m > n) {
      cmatrixluprec(a, offs, n, n, pivots, tmp);
      for (i = 0; i < n; i++) {
         ae_v_cmove(tmp->xC, 1, &a->xyC[offs + n][offs + i], a->stride, "N", m - n);
         ae_v_cmove(&a->xyC[offs + n][offs + i], a->stride, &a->xyC[offs + n][pivots->xZ[offs + i]], a->stride, "N", m - n);
         ae_v_cmove(&a->xyC[offs + n][pivots->xZ[offs + i]], a->stride, tmp->xC, 1, "N", m - n);
      }
      cmatrixrighttrsm(m - n, n, a, offs, offs, true, true, 0, a, offs + n, offs);
      return;
   }
   m1 = ablascomplexsplitlength(a, m), m2 = m - m1;
   cmatrixluprec(a, offs, m1, n, pivots, tmp);
   if (m2 > 0) {
      for (i = 0; i < m1; i++) {
         if (offs + i != pivots->xZ[offs + i]) {
            ae_v_cmove(tmp->xC, 1, &a->xyC[offs + m1][offs + i], a->stride, "N", m2);
            ae_v_cmove(&a->xyC[offs + m1][offs + i], a->stride, &a->xyC[offs + m1][pivots->xZ[offs + i]], a->stride, "N", m - m1);
            ae_v_cmove(&a->xyC[offs + m1][pivots->xZ[offs + i]], a->stride, tmp->xC, 1, "N", m - m1);
         }
      }
      cmatrixrighttrsm(m2, m1, a, offs, offs, true, true, 0, a, offs + m1, offs);
      cmatrixgemm(m - m1, n - m1, m1, ae_complex_from_d(-1.0), a, offs + m1, offs, 0, a, offs, offs + m1, 0, ae_complex_from_d(1.0), a, offs + m1, offs + m1);
      cmatrixluprec(a, offs + m1, m - m1, n - m1, pivots, tmp);
      for (i = 0; i < m2; i++) {
         if (offs + m1 + i != pivots->xZ[offs + m1 + i]) {
            ae_v_cmove(tmp->xC, 1, &a->xyC[offs][offs + m1 + i], a->stride, "N", m1);
            ae_v_cmove(&a->xyC[offs][offs + m1 + i], a->stride, &a->xyC[offs][pivots->xZ[offs + m1 + i]], a->stride, "N", m1);
            ae_v_cmove(&a->xyC[offs][pivots->xZ[offs + m1 + i]], a->stride, tmp->xC, 1, "N", m1);
         }
      }
   }
}

// Real PLU kernel
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      June 30, 1992
static void dlu_rmatrixplu2(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, RVector *tmp) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t jp;
   double s;
   if (m == 0 || n == 0) {
      return;
   }
   for (j = 0; j < imin2(m, n); j++) {
      jp = j;
      for (i = j + 1; i < m; i++) {
         if (fabs(a->xyR[offs + i][offs + j]) > fabs(a->xyR[offs + jp][offs + j])) {
            jp = i;
         }
      }
      pivots->xZ[offs + j] = offs + jp;
      if (a->xyR[offs + jp][offs + j] != 0.0) {
         if (jp != j) {
            for (i = 0; i < n; i++) {
               swapr(&a->xyR[offs + j][offs + i], &a->xyR[offs + jp][offs + i]);
            }
         }
         if (j + 1 < m) {
            s = 1 / a->xyR[offs + j][offs + j];
            ae_v_muld(&a->xyR[offs + j + 1][offs + j], a->stride, m - j - 1, s);
         }
      }
      if (j < imin2(m, n) - 1) {
         ae_v_move(tmp->xR, 1, &a->xyR[offs + j + 1][offs + j], a->stride, m - j - 1);
         ae_v_moveneg(&tmp->xR[m], 1, &a->xyR[offs + j][offs + j + 1], 1, n - j - 1);
         rmatrixrank1(m - j - 1, n - j - 1, a, offs + j + 1, offs + j + 1, tmp, 0, tmp, m);
      }
   }
}

// Complex PLU kernel
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      June 30, 1992
static void dlu_cmatrixplu2(CMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, CVector *tmp) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t jp;
   ae_complex s;
   if (m == 0 || n == 0) {
      return;
   }
   for (j = 0; j < imin2(m, n); j++) {
      jp = j;
      for (i = j + 1; i < m; i++) {
         if (ae_c_abs(a->xyC[offs + i][offs + j]) > ae_c_abs(a->xyC[offs + jp][offs + j])) {
            jp = i;
         }
      }
      pivots->xZ[offs + j] = offs + jp;
      if (ae_c_neq_d(a->xyC[offs + jp][offs + j], 0.0)) {
         if (jp != j) {
            for (i = 0; i < n; i++) {
               swapc(&a->xyC[offs + j][offs + i], &a->xyC[offs + jp][offs + i]);
            }
         }
         if (j + 1 < m) {
            s = ae_c_d_div(1, a->xyC[offs + j][offs + j]);
            ae_v_cmulc(&a->xyC[offs + j + 1][offs + j], a->stride, m - j - 1, s);
         }
      }
      if (j < imin2(m, n) - 1) {
         ae_v_cmove(tmp->xC, 1, &a->xyC[offs + j + 1][offs + j], a->stride, "N", m - j - 1);
         ae_v_cmoveneg(&tmp->xC[m], 1, &a->xyC[offs + j][offs + j + 1], 1, "N", n - j - 1);
         cmatrixrank1(m - j - 1, n - j - 1, a, offs + j + 1, offs + j + 1, tmp, 0, tmp, m);
      }
   }
}

// Recurrent real LU subroutine.
// Never call it directly.
//
// ALGLIB Routine: Copyright 04.01.2010 by Sergey Bochkanov
void rmatrixplurec(RMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, RVector *tmp) {
   ae_int_t i;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   if (n <= tsb) {
      if (rmatrixplumkl(a, offs, m, n, pivots)) {
         return;
      }
   }
   if (n <= tsa) {
      dlu_rmatrixplu2(a, offs, m, n, pivots, tmp);
      return;
   }
   if (n > m) {
      rmatrixplurec(a, offs, m, m, pivots, tmp);
      for (i = 0; i < m; i++) {
         ae_v_move(tmp->xR, 1, &a->xyR[offs + i][offs + m], 1, n - m);
         ae_v_move(&a->xyR[offs + i][offs + m], 1, &a->xyR[pivots->xZ[offs + i]][offs + m], 1, n - m);
         ae_v_move(&a->xyR[pivots->xZ[offs + i]][offs + m], 1, tmp->xR, 1, n - m);
      }
      rmatrixlefttrsm(m, n - m, a, offs, offs, false, true, 0, a, offs, offs + m);
      return;
   }
   if (n > tsb) {
      n1 = tsb;
      n2 = n - n1;
   } else {
      n1 = tiledsplit(n, tsa), n2 = n - n1;
   }
   rmatrixplurec(a, offs, m, n1, pivots, tmp);
   if (n2 > 0) {
      for (i = 0; i < n1; i++) {
         if (offs + i != pivots->xZ[offs + i]) {
            ae_v_move(tmp->xR, 1, &a->xyR[offs + i][offs + n1], 1, n2);
            ae_v_move(&a->xyR[offs + i][offs + n1], 1, &a->xyR[pivots->xZ[offs + i]][offs + n1], 1, n - n1);
            ae_v_move(&a->xyR[pivots->xZ[offs + i]][offs + n1], 1, tmp->xR, 1, n - n1);
         }
      }
      rmatrixlefttrsm(n1, n2, a, offs, offs, false, true, 0, a, offs, offs + n1);
      rmatrixgemm(m - n1, n - n1, n1, -1.0, a, offs + n1, offs, 0, a, offs, offs + n1, 0, 1.0, a, offs + n1, offs + n1);
      rmatrixplurec(a, offs + n1, m - n1, n - n1, pivots, tmp);
      for (i = 0; i < n2; i++) {
         if (offs + n1 + i != pivots->xZ[offs + n1 + i]) {
            ae_v_move(tmp->xR, 1, &a->xyR[offs + n1 + i][offs], 1, n1);
            ae_v_move(&a->xyR[offs + n1 + i][offs], 1, &a->xyR[pivots->xZ[offs + n1 + i]][offs], 1, n1);
            ae_v_move(&a->xyR[pivots->xZ[offs + n1 + i]][offs], 1, tmp->xR, 1, n1);
         }
      }
   }
}

// Recurrent complex LU subroutine.
// Never call it directly.
//
// ALGLIB Routine: Copyright 04.01.2010 by Sergey Bochkanov
void cmatrixplurec(CMatrix *a, ae_int_t offs, ae_int_t m, ae_int_t n, ZVector *pivots, CVector *tmp) {
   ae_int_t i;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   if (n <= tsa) {
      dlu_cmatrixplu2(a, offs, m, n, pivots, tmp);
      return;
   }
   if (n > m) {
      cmatrixplurec(a, offs, m, m, pivots, tmp);
      for (i = 0; i < m; i++) {
         ae_v_cmove(tmp->xC, 1, &a->xyC[offs + i][offs + m], 1, "N", n - m);
         ae_v_cmove(&a->xyC[offs + i][offs + m], 1, &a->xyC[pivots->xZ[offs + i]][offs + m], 1, "N", n - m);
         ae_v_cmove(&a->xyC[pivots->xZ[offs + i]][offs + m], 1, tmp->xC, 1, "N", n - m);
      }
      cmatrixlefttrsm(m, n - m, a, offs, offs, false, true, 0, a, offs, offs + m);
      return;
   }
   if (n > tsb) {
      n1 = tsb;
      n2 = n - n1;
   } else {
      n1 = tiledsplit(n, tsa), n2 = n - n1;
   }
   cmatrixplurec(a, offs, m, n1, pivots, tmp);
   if (n2 > 0) {
      for (i = 0; i < n1; i++) {
         if (offs + i != pivots->xZ[offs + i]) {
            ae_v_cmove(tmp->xC, 1, &a->xyC[offs + i][offs + n1], 1, "N", n2);
            ae_v_cmove(&a->xyC[offs + i][offs + n1], 1, &a->xyC[pivots->xZ[offs + i]][offs + n1], 1, "N", n - n1);
            ae_v_cmove(&a->xyC[pivots->xZ[offs + i]][offs + n1], 1, tmp->xC, 1, "N", n - n1);
         }
      }
      cmatrixlefttrsm(n1, n2, a, offs, offs, false, true, 0, a, offs, offs + n1);
      cmatrixgemm(m - n1, n - n1, n1, ae_complex_from_d(-1.0), a, offs + n1, offs, 0, a, offs, offs + n1, 0, ae_complex_from_d(1.0), a, offs + n1, offs + n1);
      cmatrixplurec(a, offs + n1, m - n1, n - n1, pivots, tmp);
      for (i = 0; i < n2; i++) {
         if (offs + n1 + i != pivots->xZ[offs + n1 + i]) {
            ae_v_cmove(tmp->xC, 1, &a->xyC[offs + n1 + i][offs], 1, "N", n1);
            ae_v_cmove(&a->xyC[offs + n1 + i][offs], 1, &a->xyC[pivots->xZ[offs + n1 + i]][offs], 1, "N", n1);
            ae_v_cmove(&a->xyC[pivots->xZ[offs + n1 + i]][offs], 1, tmp->xC, 1, "N", n1);
         }
      }
   }
}
} // end of namespace alglib_impl

// === SPTRF Package ===
// Depends on: SPARSE, DLU
namespace alglib_impl {
static const double sptrf_densebnd = 0.10;
static const ae_int_t sptrf_slswidth = 8;

// This function drops sequence #I from the structure
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sluv2list1dropsequence(sluv2list1matrix *a, ae_int_t i) {
   a->idxfirst.xZ[i] = -1;
}

// This function appends column with id=ID to the dense trail (column IDs are
// integer numbers in [0,N) which can be used to track column permutations).
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_densetrailappendcolumn(sluv2densetrail *d, RVector *x, ae_int_t id) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t targetidx;
   n = d->n;
// Reallocate storage
   rmatrixgrowcolsto(&d->d, d->ndense + 1, n);
// Copy to dense storage:
// * BUpper
// * BTrail
// Remove from sparse storage
   targetidx = d->ndense;
   for (i = 0; i < n; i++) {
      d->d.xyR[i][targetidx] = x->xR[i];
   }
   d->did.xZ[targetidx] = id;
   d->ndense = targetidx + 1;
}

// This function densifies I1-th column of the sparse trail.
//
// PARAMETERS:
//     A           -   sparse trail
//     I1          -   column index
//     BUpper      -   upper rectangular submatrix, updated during densification
//                     of the columns (densified columns are removed)
//     DTrail      -   dense trail, receives densified columns from sparse
//                     trail and BUpper
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sparsetraildensify(sluv2sparsetrail *a, ae_int_t i1, sluv2list1matrix *bupper, sluv2densetrail *dtrail) {
   ae_int_t n;
   ae_int_t k;
   ae_int_t i;
   ae_int_t jp;
   ae_int_t entry;
   ae_int_t pprev;
   ae_int_t pnext;
   n = a->n;
   k = a->k;
   ae_assert(k < n, "SparseTrailDensify: integrity check failed");
   ae_assert(k <= i1, "SparseTrailDensify: integrity check failed");
   ae_assert(!a->isdensified.xB[i1], "SparseTrailDensify: integrity check failed");
// Offload items [0,K) of densified column from BUpper
   for (i = 0; i < n; i++) {
      a->tmp0.xR[i] = 0.0;
   }
   jp = bupper->idxfirst.xZ[i1];
   while (jp >= 0) {
      a->tmp0.xR[bupper->strgidx.xZ[2 * jp + 1]] = bupper->strgval.xR[jp];
      jp = bupper->strgidx.xZ[2 * jp];
   }
   sptrf_sluv2list1dropsequence(bupper, i1);
// Offload items [K,N) of densified column from BLeft
   entry = a->slscolptr.xZ[i1];
   while (entry >= 0) {
   // Offload element
      i = a->slsidx.xZ[entry * sptrf_slswidth + 4];
      a->tmp0.xR[i] = a->slsval.xR[entry];
   // Remove element from the row list
      pprev = a->slsidx.xZ[entry * sptrf_slswidth + 2];
      pnext = a->slsidx.xZ[entry * sptrf_slswidth + 3];
      if (pprev >= 0) {
         a->slsidx.xZ[pprev * sptrf_slswidth + 3] = pnext;
      } else {
         a->slsrowptr.xZ[i] = pnext;
      }
      if (pnext >= 0) {
         a->slsidx.xZ[pnext * sptrf_slswidth + 2] = pprev;
      }
   // Select next entry
      entry = a->slsidx.xZ[entry * sptrf_slswidth + 1];
   }
// Densify
   a->nzc.xZ[i1] = 0;
   a->isdensified.xB[i1] = true;
   a->slscolptr.xZ[i1] = -1;
   sptrf_densetrailappendcolumn(dtrail, &a->tmp0, a->colid.xZ[i1]);
}

// This function appends rank-1 update to the sparse trail.  Dense  trail  is
// not  updated  here,  but  we  may  move some columns to dense trail during
// update (i.e. densify them). Thus, you have to update  dense  trail  BEFORE
// you start updating sparse one (otherwise, recently densified columns  will
// be updated twice).
//
// PARAMETERS:
//     A           -   sparse trail
//     V0I, V0R    -   update column returned by SparseTrailPivotOut (MUST be
//                     array[N] independently of the NZ0).
//     NZ0         -   non-zero count for update column
//     V1I, V1R    -   update row returned by SparseTrailPivotOut
//     NZ1         -   non-zero count for update row
//     BUpper      -   upper rectangular submatrix, updated during densification
//                     of the columns (densified columns are removed)
//     DTrail      -   dense trail, receives densified columns from sparse
//                     trail and BUpper
//     DensificationSupported- if False, no densification is performed
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sparsetrailupdate(sluv2sparsetrail *a, ZVector *v0i, RVector *v0r, ae_int_t nz0, ZVector *v1i, RVector *v1r, ae_int_t nz1, sluv2list1matrix *bupper, sluv2densetrail *dtrail, bool densificationsupported) {
   ae_int_t n;
   ae_int_t k;
   ae_int_t i;
   ae_int_t j;
   ae_int_t i0;
   ae_int_t i1;
   double v1;
   ae_int_t densifyabove;
   ae_int_t nnz;
   ae_int_t entry;
   ae_int_t newentry;
   ae_int_t pprev;
   ae_int_t pnext;
   ae_int_t p;
   ae_int_t nexti;
   ae_int_t newoffs;
   n = a->n;
   k = a->k;
   ae_assert(k < n, "SparseTrailPivotOut: integrity check failed");
   densifyabove = RoundZ(sptrf_densebnd * (n - k)) + 1;
   ae_assert(v0i->cnt >= nz0 + 1, "SparseTrailUpdate: integrity check failed");
   ae_assert(v0r->cnt >= nz0 + 1, "SparseTrailUpdate: integrity check failed");
   v0i->xZ[nz0] = -1;
   v0r->xR[nz0] = 0.0;
// Update sparse representation
   ivectorgrowto(&a->slsidx, (a->slsused + nz0 * nz1) * sptrf_slswidth);
   rvectorgrowto(&a->slsval, a->slsused + nz0 * nz1);
   for (j = 0; j < nz1; j++) {
      if (nz0 == 0) {
         continue;
      }
      i1 = v1i->xZ[j];
      v1 = v1r->xR[j];
   // Update column #I1
      nnz = a->nzc.xZ[i1];
      i = 0;
      i0 = v0i->xZ[i];
      entry = a->slscolptr.xZ[i1];
      pprev = -1;
      while (i < nz0) {
      // Handle possible fill-in happening BEFORE already existing
      // entry of the column list (or simply fill-in, if no entry
      // is present).
         pnext = entry;
         if (entry >= 0) {
            nexti = a->slsidx.xZ[entry * sptrf_slswidth + 4];
         } else {
            nexti = n + 1;
         }
         while (i < nz0) {
            if (i0 >= nexti) {
               break;
            }
         // Allocate new entry, store column/row/value
            newentry = a->slsused;
            a->slsused = newentry + 1;
            nnz++;
            newoffs = newentry * sptrf_slswidth;
            a->slsidx.xZ[newoffs + 4] = i0;
            a->slsidx.xZ[newoffs + 5] = i1;
            a->slsval.xR[newentry] = -v1 * v0r->xR[i];
         // Insert entry into column list
            a->slsidx.xZ[newoffs] = pprev;
            a->slsidx.xZ[newoffs + 1] = pnext;
            if (pprev >= 0) {
               a->slsidx.xZ[pprev * sptrf_slswidth + 1] = newentry;
            } else {
               a->slscolptr.xZ[i1] = newentry;
            }
            if (entry >= 0) {
               a->slsidx.xZ[entry * sptrf_slswidth] = newentry;
            }
         // Insert entry into row list
            p = a->slsrowptr.xZ[i0];
            a->slsidx.xZ[newoffs + 2] = -1;
            a->slsidx.xZ[newoffs + 3] = p;
            if (p >= 0) {
               a->slsidx.xZ[p * sptrf_slswidth + 2] = newentry;
            }
            a->slsrowptr.xZ[i0] = newentry;
         // Advance pointers
            pprev = newentry;
            i++;
            i0 = v0i->xZ[i];
         }
         if (i >= nz0) {
            break;
         }
      // Update already existing entry of the column list, if needed
         if (entry >= 0) {
            if (i0 == nexti) {
               a->slsval.xR[entry] -= v1 * v0r->xR[i];
               i++;
               i0 = v0i->xZ[i];
            }
            pprev = entry;
         }
      // Advance to the next pre-existing entry (if present)
         if (entry >= 0) {
            entry = a->slsidx.xZ[entry * sptrf_slswidth + 1];
         }
      }
      a->nzc.xZ[i1] = nnz;
   // Densify column if needed
      if (densificationsupported && nnz > densifyabove && !a->isdensified.xB[i1]) {
         sptrf_sparsetraildensify(a, i1, bupper, dtrail);
      }
   }
}

// This function initialized rectangular submatrix structure.
//
// After initialization this structure stores  matrix[N,0],  which contains N
// rows (sequences), stored as single-linked lists.
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sluv2list1init(ae_int_t n, sluv2list1matrix *a) {
   ae_int_t i;
   ae_assert(n >= 1, "SLUV2List1Init: N<1");
   a->nfixed = n;
   a->ndynamic = 0;
   a->nallocated = n;
   a->nused = 0;
   ivectorgrowto(&a->idxfirst, n);
   ivectorgrowto(&a->strgidx, 2 * a->nallocated);
   rvectorgrowto(&a->strgval, a->nallocated);
   for (i = 0; i < n; i++) {
      a->idxfirst.xZ[i] = -1;
   }
}

// This function swaps sequences #I and #J stored by the structure
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sluv2list1swap(sluv2list1matrix *a, ae_int_t i, ae_int_t j) {
   swapi(&a->idxfirst.xZ[i], &a->idxfirst.xZ[j]);
}

// This function appends sequence from the structure to the sparse matrix.
//
// It is assumed that S is a lower triangular  matrix,  and A stores strictly
// lower triangular elements (no diagonal ones!). You can explicitly  control
// whether you want to add diagonal elements or not.
//
// Output matrix is assumed to be stored in CRS format and  to  be  partially
// initialized (up to, but not including, Dst-th row). DIdx and UIdx are  NOT
// updated by this function as well as NInitialized.
//
// Inputs:
//     A           -   rectangular matrix structure
//     Src         -   sequence (row or column) index in the structure
//     HasDiagonal -   whether we want to add diagonal element
//     D           -   diagonal element, if HasDiagonal=True
//     NZMAX       -   maximum estimated number of non-zeros in the row,
//                     this function will preallocate storage in the output
//                     matrix.
//     S           -   destination matrix in CRS format, partially initialized
//     Dst         -   destination row index
//
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sluv2list1appendsequencetomatrix(sluv2list1matrix *a, ae_int_t src, bool hasdiagonal, double d, ae_int_t nzmax, sparsematrix *s, ae_int_t dst) {
   ae_int_t i;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t jp;
   ae_int_t nnz;
   i0 = s->ridx.xZ[dst];
   ivectorgrowto(&s->idx, i0 + nzmax);
   rvectorgrowto(&s->vals, i0 + nzmax);
   if (hasdiagonal) {
      i1 = i0 + nzmax - 1;
      s->idx.xZ[i1] = dst;
      s->vals.xR[i1] = d;
      nnz = 1;
   } else {
      i1 = i0 + nzmax;
      nnz = 0;
   }
   jp = a->idxfirst.xZ[src];
   while (jp >= 0) {
      i1--;
      s->idx.xZ[i1] = a->strgidx.xZ[2 * jp + 1];
      s->vals.xR[i1] = a->strgval.xR[jp];
      nnz++;
      jp = a->strgidx.xZ[2 * jp];
   }
   for (i = 0; i < nnz; i++) {
      s->idx.xZ[i0 + i] = s->idx.xZ[i1 + i];
      s->vals.xR[i0 + i] = s->vals.xR[i1 + i];
   }
   s->ridx.xZ[dst + 1] = s->ridx.xZ[dst] + nnz;
}

// This function appends sparse column to the  matrix,  increasing  its  size
// from [N,K] to [N,K+1]
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sluv2list1pushsparsevector(sluv2list1matrix *a, ZVector *si, RVector *sv, ae_int_t nz) {
   ae_int_t idx;
   ae_int_t i;
   ae_int_t k;
   ae_int_t nused;
   double v;
// Fetch matrix size, increase
   k = a->ndynamic;
   ae_assert(k < a->nfixed, "Assertion failed");
   a->ndynamic = k + 1;
// Allocate new storage if needed
   nused = a->nused;
   a->nallocated = imax2(a->nallocated, nused + nz);
   ivectorgrowto(&a->strgidx, 2 * a->nallocated);
   rvectorgrowto(&a->strgval, a->nallocated);
// Append to list
   for (idx = 0; idx < nz; idx++) {
      i = si->xZ[idx];
      v = sv->xR[idx];
      a->strgidx.xZ[2 * nused] = a->idxfirst.xZ[i];
      a->strgidx.xZ[2 * nused + 1] = k;
      a->strgval.xR[nused] = v;
      a->idxfirst.xZ[i] = nused;
      nused++;
   }
   a->nused = nused;
}

// This function initializes dense trail, by default it is matrix[N,0]
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_densetrailinit(sluv2densetrail *d, ae_int_t n) {
   ae_int_t excessivesize;
// Note: excessive rows are allocated to accomodate for situation when
//       this buffer is used to solve successive problems with increasing
//       sizes.
   excessivesize = imax2(RoundZ(1.333 * n), n);
   d->n = n;
   d->ndense = 0;
   vectorsetlengthatleast(&d->did, n);
   if (d->d.rows <= excessivesize) {
      matrixsetlengthatleast(&d->d, n, 1);
   } else {
      ae_matrix_set_length(&d->d, excessivesize, 1);
   }
}

// This function initializes sparse trail from the sparse matrix. By default,
// sparse trail spans columns and rows in [0,N)  range.  Subsequent  pivoting
// out of rows/columns changes its range to [K,N), [K+1,N) and so on.
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sparsetrailinit(sparsematrix *s, sluv2sparsetrail *a) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t n;
   ae_int_t j0;
   ae_int_t j1;
   ae_int_t jj;
   ae_int_t p;
   ae_int_t slsused;
   ae_assert(s->m == s->n, "SparseTrailInit: M != N");
   ae_assert(s->matrixtype == 1, "SparseTrailInit: non-CRS input");
   n = s->n;
   a->n = s->n;
   a->k = 0;
   vectorsetlengthatleast(&a->nzc, n);
   vectorsetlengthatleast(&a->colid, n);
   vectorsetlengthatleast(&a->tmp0, n);
   for (i = 0; i < n; i++) {
      a->colid.xZ[i] = i;
   }
   vectorsetlengthatleast(&a->isdensified, n);
   for (i = 0; i < n; i++) {
      a->isdensified.xB[i] = false;
   }
// Working set of columns
   a->maxwrkcnt = iboundval(RoundZ(1 + (double)n / 3.0), 1, imin2(n, 50));
   a->wrkcnt = 0;
   vectorsetlengthatleast(&a->wrkset, a->maxwrkcnt);
// Sparse linked storage (SLS). Store CRS matrix to SLS format,
// row by row, starting from the last one.
   vectorsetlengthatleast(&a->slscolptr, n);
   vectorsetlengthatleast(&a->slsrowptr, n);
   vectorsetlengthatleast(&a->slsidx, s->ridx.xZ[n] * sptrf_slswidth);
   vectorsetlengthatleast(&a->slsval, s->ridx.xZ[n]);
   for (i = 0; i < n; i++) {
      a->nzc.xZ[i] = 0;
   }
   for (i = 0; i < n; i++) {
      a->slscolptr.xZ[i] = -1;
      a->slsrowptr.xZ[i] = -1;
   }
   slsused = 0;
   for (i = n - 1; i >= 0; i--) {
      j0 = s->ridx.xZ[i];
      j1 = s->ridx.xZ[i + 1] - 1;
      for (jj = j1; jj >= j0; jj--) {
         j = s->idx.xZ[jj];
      // Update non-zero counts for columns
         a->nzc.xZ[j]++;
      // Insert into column list
         p = a->slscolptr.xZ[j];
         if (p >= 0) {
            a->slsidx.xZ[p * sptrf_slswidth] = slsused;
         }
         a->slsidx.xZ[slsused * sptrf_slswidth] = -1;
         a->slsidx.xZ[slsused * sptrf_slswidth + 1] = p;
         a->slscolptr.xZ[j] = slsused;
      // Insert into row list
         p = a->slsrowptr.xZ[i];
         if (p >= 0) {
            a->slsidx.xZ[p * sptrf_slswidth + 2] = slsused;
         }
         a->slsidx.xZ[slsused * sptrf_slswidth + 2] = -1;
         a->slsidx.xZ[slsused * sptrf_slswidth + 3] = p;
         a->slsrowptr.xZ[i] = slsused;
      // Store index and value
         a->slsidx.xZ[slsused * sptrf_slswidth + 4] = i;
         a->slsidx.xZ[slsused * sptrf_slswidth + 5] = j;
         a->slsval.xR[slsused] = s->vals.xR[jj];
         slsused++;
      }
   }
   a->slsused = slsused;
}

// This function searches for a appropriate pivot column/row.
//
// If there exists non-densified column, it returns indexes of  pivot  column
// and row, with most sparse column selected for column pivoting, and largest
// element selected for row pivoting. Function result is True.
//
// PivotType=1 means that no column pivoting is performed
// PivotType=2 means that both column and row pivoting are supported
//
// If all columns were densified, False is returned.
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static bool sptrf_sparsetrailfindpivot(sluv2sparsetrail *a, ae_int_t pivottype, ae_int_t *ipiv, ae_int_t *jpiv) {
   ae_int_t n;
   ae_int_t k;
   ae_int_t j;
   ae_int_t jp;
   ae_int_t entry;
   ae_int_t nz;
   ae_int_t maxwrknz;
   ae_int_t nnzbest;
   double s;
   double bbest;
   ae_int_t wrk0;
   ae_int_t wrk1;
   bool result;
   *ipiv = 0;
   *jpiv = 0;
   n = a->n;
   k = a->k;
   nnzbest = n + 1;
   *jpiv = -1;
   *ipiv = -1;
   result = true;
// Select pivot column
   if (pivottype == 1) {
   // No column pivoting
      ae_assert(!a->isdensified.xB[k], "SparseTrailFindPivot: integrity check failed");
      *jpiv = k;
   } else {
   // Find pivot column
      while (true) {
      // Scan working set (if non-empty) for good columns
         maxwrknz = a->maxwrknz;
         for (j = 0; j < a->wrkcnt; j++) {
            jp = a->wrkset.xZ[j];
            if (jp < k) {
               continue;
            }
            if (a->isdensified.xB[jp]) {
               continue;
            }
            nz = a->nzc.xZ[jp];
            if (nz > maxwrknz) {
               continue;
            }
            if (*jpiv < 0 || nz < nnzbest) {
               nnzbest = nz;
               *jpiv = jp;
            }
         }
         if (*jpiv >= 0) {
            break;
         }
      // Well, nothing found. Recompute working set:
      // * determine most sparse unprocessed yet column
      // * gather all columns with density in [Wrk0,Wrk1) range,
      //   increase range, repeat, until working set is full
         a->wrkcnt = 0;
         a->maxwrknz = 0;
         wrk0 = n + 1;
         for (jp = k; jp < n; jp++) {
            if (!a->isdensified.xB[jp] && a->nzc.xZ[jp] < wrk0) {
               wrk0 = a->nzc.xZ[jp];
            }
         }
         if (wrk0 > n) {
         // Only densified columns are present, exit.
            result = false;
            return result;
         }
         wrk1 = wrk0 + 1;
         while (a->wrkcnt < a->maxwrkcnt && wrk0 <= n) {
         // Find columns with non-zero count in [Wrk0,Wrk1) range
            for (jp = k; jp < n; jp++) {
               if (a->wrkcnt == a->maxwrkcnt) {
                  break;
               }
               if (a->isdensified.xB[jp]) {
                  continue;
               }
               if (a->nzc.xZ[jp] >= wrk0 && a->nzc.xZ[jp] < wrk1) {
                  a->wrkset.xZ[a->wrkcnt] = jp;
                  a->wrkcnt++;
                  a->maxwrknz = imax2(a->maxwrknz, a->nzc.xZ[jp]);
               }
            }
         // Advance scan range
            jp = RoundZ(1.41 * (wrk1 - wrk0)) + 1;
            wrk0 = wrk1;
            wrk1 = wrk0 + jp;
         }
      }
   }
// Select pivot row
   bbest = 0.0;
   entry = a->slscolptr.xZ[*jpiv];
   while (entry >= 0) {
      s = fabs(a->slsval.xR[entry]);
      if (*ipiv < 0 || s > bbest) {
         bbest = s;
         *ipiv = a->slsidx.xZ[entry * sptrf_slswidth + 4];
      }
      entry = a->slsidx.xZ[entry * sptrf_slswidth + 1];
   }
   if (*ipiv < 0) {
      *ipiv = k;
   }
   return result;
}

// This function pivots out specified row and column.
//
// Sparse trail range changes from [K,N) to [K+1,N).
//
// V0I, V0R, V1I, V1R must be preallocated arrays[N].
//
// Following data are returned:
// * UU - diagonal element (pivoted out), can be zero
// * V0I, V0R, NZ0 - sparse column pivoted out to the left (after permutation
//   is applied to its elements) and divided by UU.
//   V0I is array[NZ0] which stores row indexes in [K+1,N) range, V0R  stores
//   values.
// * V1I, V1R, NZ1 - sparse row pivoted out to the top.
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
static void sptrf_sparsetrailpivotout(sluv2sparsetrail *a, ae_int_t ipiv, ae_int_t jpiv, double *uu, ZVector *v0i, RVector *v0r, ae_int_t *nz0, ZVector *v1i, RVector *v1r, ae_int_t *nz1) {
   ae_int_t n;
   ae_int_t k;
   ae_int_t i;
   ae_int_t j;
   ae_int_t entry;
   double s;
   ae_int_t pos0k;
   ae_int_t pos0piv;
   ae_int_t pprev;
   ae_int_t pnext;
   ae_int_t pnextnext;
   *uu = 0;
   *nz0 = 0;
   *nz1 = 0;
   n = a->n;
   k = a->k;
   ae_assert(k < n, "SparseTrailPivotOut: integrity check failed");
// Pivot out column JPiv from the sparse linked storage:
// * remove column JPiv from the matrix
// * update column K:
//   * change element indexes after it is permuted to JPiv
//   * resort rows affected by move K->JPiv
//
// NOTE: this code leaves V0I/V0R/NZ0 in the unfinalized state,
//       i.e. these arrays do not account for pivoting performed
//       on rows. They will be post-processed later.
   *nz0 = 0;
   pos0k = -1;
   pos0piv = -1;
   entry = a->slscolptr.xZ[jpiv];
   while (entry >= 0) {
   // Offload element
      i = a->slsidx.xZ[entry * sptrf_slswidth + 4];
      v0i->xZ[*nz0] = i;
      v0r->xR[*nz0] = a->slsval.xR[entry];
      if (i == k) {
         pos0k = *nz0;
      }
      if (i == ipiv) {
         pos0piv = *nz0;
      }
      ++*nz0;
   // Remove element from the row list
      pprev = a->slsidx.xZ[entry * sptrf_slswidth + 2];
      pnext = a->slsidx.xZ[entry * sptrf_slswidth + 3];
      if (pprev >= 0) {
         a->slsidx.xZ[pprev * sptrf_slswidth + 3] = pnext;
      } else {
         a->slsrowptr.xZ[i] = pnext;
      }
      if (pnext >= 0) {
         a->slsidx.xZ[pnext * sptrf_slswidth + 2] = pprev;
      }
   // Select next entry
      entry = a->slsidx.xZ[entry * sptrf_slswidth + 1];
   }
   entry = a->slscolptr.xZ[k];
   a->slscolptr.xZ[jpiv] = entry;
   while (entry >= 0) {
   // Change column index
      a->slsidx.xZ[entry * sptrf_slswidth + 5] = jpiv;
   // Next entry
      entry = a->slsidx.xZ[entry * sptrf_slswidth + 1];
   }
// Post-process V0, account for pivoting.
// Compute diagonal element UU.
   *uu = 0.0;
   if (pos0k >= 0 || pos0piv >= 0) {
   // Apply permutation to rows of pivoted out column, specific
   // implementation depends on the sparsity at locations #Pos0K
   // and #Pos0Piv of the V0 array.
      if (pos0k >= 0 && pos0piv >= 0) {
      // Obtain diagonal element
         *uu = v0r->xR[pos0piv];
         if (*uu != 0) {
            s = 1 / (*uu);
         } else {
            s = 1.0;
         }
      // Move pivoted out element, shift array by one in order
      // to remove heading diagonal element (not needed here
      // anymore).
         v0r->xR[pos0piv] = v0r->xR[pos0k];
         for (i = 0; i < *nz0 - 1; i++) {
            v0i->xZ[i] = v0i->xZ[i + 1];
            v0r->xR[i] = v0r->xR[i + 1] * s;
         }
         --*nz0;
      }
      if (pos0k >= 0 && pos0piv < 0) {
      // Diagonal element is zero
         *uu = 0.0;
      // Pivot out element, reorder array
         v0i->xZ[pos0k] = ipiv;
         for (i = pos0k; i < *nz0 - 1; i++) {
            if (v0i->xZ[i] < v0i->xZ[i + 1]) {
               break;
            }
            swapi(&v0i->xZ[i], &v0i->xZ[i + 1]);
            swapr(&v0r->xR[i], &v0r->xR[i + 1]);
         }
      }
      if (pos0k < 0 && pos0piv >= 0) {
      // Get diagonal element
         *uu = v0r->xR[pos0piv];
         if (*uu != 0) {
            s = 1 / (*uu);
         } else {
            s = 1.0;
         }
      // Shift array past the pivoted in element by one
      // in order to remove pivot
         for (i = 0; i < pos0piv; i++) {
            v0r->xR[i] *= s;
         }
         for (i = pos0piv; i < *nz0 - 1; i++) {
            v0i->xZ[i] = v0i->xZ[i + 1];
            v0r->xR[i] = v0r->xR[i + 1] * s;
         }
         --*nz0;
      }
   }
// Pivot out row IPiv from the sparse linked storage:
// * remove row IPiv from the matrix
// * reindex elements of row K after it is permuted to IPiv
// * apply permutation to the cols of the pivoted out row,
//   resort columns
   *nz1 = 0;
   entry = a->slsrowptr.xZ[ipiv];
   while (entry >= 0) {
   // Offload element
      j = a->slsidx.xZ[entry * sptrf_slswidth + 5];
      v1i->xZ[*nz1] = j;
      v1r->xR[*nz1] = a->slsval.xR[entry];
      ++*nz1;
   // Remove element from the column list
      pprev = a->slsidx.xZ[entry * sptrf_slswidth];
      pnext = a->slsidx.xZ[entry * sptrf_slswidth + 1];
      if (pprev >= 0) {
         a->slsidx.xZ[pprev * sptrf_slswidth + 1] = pnext;
      } else {
         a->slscolptr.xZ[j] = pnext;
      }
      if (pnext >= 0) {
         a->slsidx.xZ[pnext * sptrf_slswidth] = pprev;
      }
   // Select next entry
      entry = a->slsidx.xZ[entry * sptrf_slswidth + 3];
   }
   a->slsrowptr.xZ[ipiv] = a->slsrowptr.xZ[k];
   entry = a->slsrowptr.xZ[ipiv];
   while (entry >= 0) {
   // Change row index
      a->slsidx.xZ[entry * sptrf_slswidth + 4] = ipiv;
   // Resort column affected by row pivoting
      j = a->slsidx.xZ[entry * sptrf_slswidth + 5];
      pprev = a->slsidx.xZ[entry * sptrf_slswidth];
      pnext = a->slsidx.xZ[entry * sptrf_slswidth + 1];
      while (pnext >= 0 && a->slsidx.xZ[pnext * sptrf_slswidth + 4] < ipiv) {
         pnextnext = a->slsidx.xZ[pnext * sptrf_slswidth + 1];
      // prev->next
         if (pprev >= 0) {
            a->slsidx.xZ[pprev * sptrf_slswidth + 1] = pnext;
         } else {
            a->slscolptr.xZ[j] = pnext;
         }
      // entry->prev, entry->next
         a->slsidx.xZ[entry * sptrf_slswidth] = pnext;
         a->slsidx.xZ[entry * sptrf_slswidth + 1] = pnextnext;
      // next->prev, next->next
         a->slsidx.xZ[pnext * sptrf_slswidth] = pprev;
         a->slsidx.xZ[pnext * sptrf_slswidth + 1] = entry;
      // nextnext->prev
         if (pnextnext >= 0) {
            a->slsidx.xZ[pnextnext * sptrf_slswidth] = entry;
         }
      // PPrev, Item, PNext
         pprev = pnext;
         pnext = pnextnext;
      }
   // Next entry
      entry = a->slsidx.xZ[entry * sptrf_slswidth + 3];
   }
// Reorder other structures
   swapi(&a->nzc.xZ[k], &a->nzc.xZ[jpiv]);
   swapi(&a->colid.xZ[k], &a->colid.xZ[jpiv]);
   swapb(&a->isdensified.xB[k], &a->isdensified.xB[jpiv]);
// Handle removal of col/row #K
   for (i = 0; i < *nz1; i++) {
      j = v1i->xZ[i];
      a->nzc.xZ[j]--;
   }
   a->k++;
}

// Sparse LU for square NxN CRS matrix with both row and column permutations.
//
// Represents A as Pr*L*U*Pc, where:
// * Pr is a product of row permutations Pr=Pr(0)*Pr(1)*...*Pr(n-2)*Pr(n-1)
// * Pc is a product of col permutations Pc=Pc(n-1)*Pc(n-2)*...*Pc(1)*Pc(0)
// * L is lower unitriangular
// * U is upper triangular
//
// Inputs:
//     A           -   sparse square matrix in CRS format
//     PivotType   -   pivot type:
//                     * 0 - for best pivoting available
//                     * 1 - row-only pivoting
//                     * 2 - row and column greedy pivoting  algorithm  (most
//                           sparse pivot column is selected from the trailing
//                           matrix at each step)
//     Buf         -   temporary buffer, previously allocated memory is
//                     reused as much as possible
//
// Outputs:
//     A           -   LU decomposition of A
//     PR          -   array[N], row pivots
//     PC          -   array[N], column pivots
//     Buf         -   following fields of Buf are set:
//                     * Buf.RowPermRawIdx[] - contains row permutation, with
//                       RawIdx[I]=J meaning that J-th row  of  the  original
//                       input matrix was moved to Ith position of the output
//                       factorization
//
// This function always succeeds  i.e. it ALWAYS returns valid factorization,
// but for your convenience it also  returns boolean  value  which  helps  to
// detect symbolically degenerate matrix:
// * function returns TRUE if the matrix was factorized AND symbolically
//   non-degenerate
// * function returns FALSE if the matrix was factorized but U has strictly
//   zero elements at the diagonal (the factorization is returned anyway).
//
// ALGLIB Routine: Copyright 15.01.2019 by Sergey Bochkanov
bool sptrflu(sparsematrix *a, ae_int_t pivottype, ZVector *pr, ZVector *pc, sluv2buffer *buf) {
   ae_int_t n;
   ae_int_t k;
   ae_int_t i;
   ae_int_t j;
   ae_int_t jp;
   ae_int_t i0;
   ae_int_t i1;
   ae_int_t ibest;
   ae_int_t jbest;
   double v;
   double v0;
   ae_int_t nz0;
   ae_int_t nz1;
   double uu;
   ae_int_t offs;
   ae_int_t tmpndense;
   bool densificationsupported;
   ae_int_t densifyabove;
   bool result;
   ae_assert(sparseiscrs(a), "SparseLU: A is not stored in CRS format");
   ae_assert(sparsegetnrows(a) == sparsegetncols(a), "SparseLU: non-square A");
   ae_assert(pivottype == 0 || pivottype == 1 || pivottype == 2, "SparseLU: unexpected pivot type");
   result = true;
   n = sparsegetnrows(a);
   if (pivottype == 0) {
      pivottype = 2;
   }
   densificationsupported = pivottype == 2;
//
   buf->n = n;
   vectorsetlengthatleast(&buf->rowpermrawidx, n);
   for (i = 0; i < n; i++) {
      buf->rowpermrawidx.xZ[i] = i;
   }
// Allocate storage for sparse L and U factors
//
// NOTE: SparseMatrix structure for these factors is only
//       partially initialized; we use it just as a temporary
//       storage and do not intend to use facilities of the
//       'sparse' subpackage to work with these objects.
   buf->sparsel.matrixtype = 1;
   buf->sparsel.m = n;
   buf->sparsel.n = n;
   vectorsetlengthatleast(&buf->sparsel.ridx, n + 1);
   buf->sparsel.ridx.xZ[0] = 0;
   buf->sparseut.matrixtype = 1;
   buf->sparseut.m = n;
   buf->sparseut.n = n;
   vectorsetlengthatleast(&buf->sparseut.ridx, n + 1);
   buf->sparseut.ridx.xZ[0] = 0;
// Allocate unprocessed yet part of the matrix,
// two submatrices:
// * BU, upper J rows of columns [J,N), upper submatrix
// * BL, left J  cols of rows [J,N), left submatrix
// * B1, (N-J)*(N-J) square submatrix
   sptrf_sluv2list1init(n, &buf->bleft);
   sptrf_sluv2list1init(n, &buf->bupper);
   vectorsetlengthatleast(pr, n);
   vectorsetlengthatleast(pc, n);
   vectorsetlengthatleast(&buf->v0i, n);
   vectorsetlengthatleast(&buf->v1i, n);
   vectorsetlengthatleast(&buf->v0r, n);
   vectorsetlengthatleast(&buf->v1r, n);
   sptrf_sparsetrailinit(a, &buf->strail);
// Prepare dense trail, initial densification
   sptrf_densetrailinit(&buf->dtrail, n);
   densifyabove = RoundZ(sptrf_densebnd * n) + 1;
   if (densificationsupported) {
      for (i = 0; i < n; i++) {
         if (buf->strail.nzc.xZ[i] > densifyabove) {
            sptrf_sparsetraildensify(&buf->strail, i, &buf->bupper, &buf->dtrail);
         }
      }
   }
// Process sparse part
   for (k = 0; k < n; k++) {
   // Find pivot column and pivot row
      if (!sptrf_sparsetrailfindpivot(&buf->strail, pivottype, &ibest, &jbest)) {
      // Only densified columns are left, break sparse iteration
         ae_assert(buf->dtrail.ndense + k == n, "SPTRF: integrity check failed (35741)");
         break;
      }
      pc->xZ[k] = jbest;
      pr->xZ[k] = ibest;
      swapi(&buf->rowpermrawidx.xZ[k], &buf->rowpermrawidx.xZ[ibest]);
   // Apply pivoting to BL and BU
      sptrf_sluv2list1swap(&buf->bleft, k, ibest);
      sptrf_sluv2list1swap(&buf->bupper, k, jbest);
   // Apply pivoting to sparse trail, pivot out
      sptrf_sparsetrailpivotout(&buf->strail, ibest, jbest, &uu, &buf->v0i, &buf->v0r, &nz0, &buf->v1i, &buf->v1r, &nz1);
      result = result && uu != 0;
   // Pivot dense trail
      tmpndense = buf->dtrail.ndense;
      for (i = 0; i < tmpndense; i++) {
         swapr(&buf->dtrail.d.xyR[k][i], &buf->dtrail.d.xyR[ibest][i]);
      }
   // Output to LU matrix
      sptrf_sluv2list1appendsequencetomatrix(&buf->bupper, k, true, uu, n, &buf->sparseut, k);
      sptrf_sluv2list1appendsequencetomatrix(&buf->bleft, k, false, 0.0, n, &buf->sparsel, k);
   // Extract K-th col/row of B1, generate K-th col/row of BL/BU, update NZC
      sptrf_sluv2list1pushsparsevector(&buf->bleft, &buf->v0i, &buf->v0r, nz0);
      sptrf_sluv2list1pushsparsevector(&buf->bupper, &buf->v1i, &buf->v1r, nz1);
   // Update the rest of the matrix
      if (nz0 * (nz1 + buf->dtrail.ndense) > 0) {
      // Update dense trail
      //
      // NOTE: this update MUST be performed before we update sparse trail,
      //       because sparse update may move columns to dense storage after
      //       update is performed on them. Thus, we have to avoid applying
      //       same update twice.
         if (buf->dtrail.ndense > 0) {
            tmpndense = buf->dtrail.ndense;
            for (i = 0; i < nz0; i++) {
               i0 = buf->v0i.xZ[i];
               v0 = buf->v0r.xR[i];
               for (j = 0; j < tmpndense; j++) {
                  buf->dtrail.d.xyR[i0][j] -= v0 * buf->dtrail.d.xyR[k][j];
               }
            }
         }
      // Update sparse trail
         sptrf_sparsetrailupdate(&buf->strail, &buf->v0i, &buf->v0r, nz0, &buf->v1i, &buf->v1r, nz1, &buf->bupper, &buf->dtrail, densificationsupported);
      }
   }
// Process densified trail
   if (buf->dtrail.ndense > 0) {
      tmpndense = buf->dtrail.ndense;
   // Generate column pivots to bring actual order of columns in the
   // working part of the matrix to one used for dense storage
      for (i = n - tmpndense; i < n; i++) {
         k = buf->dtrail.did.xZ[i - (n - tmpndense)];
         jp = -1;
         for (j = i; j < n; j++) {
            if (buf->strail.colid.xZ[j] == k) {
               jp = j;
               break;
            }
         }
         ae_assert(jp >= 0, "SPTRF: integrity check failed during reordering");
         swapi(&buf->strail.colid.xZ[i], &buf->strail.colid.xZ[jp]);
         pc->xZ[i] = jp;
      }
   // Perform dense LU decomposition on dense trail
      matrixsetlengthatleast(&buf->dbuf, buf->dtrail.ndense, buf->dtrail.ndense);
      for (i = 0; i < tmpndense; i++) {
         for (j = 0; j < tmpndense; j++) {
            buf->dbuf.xyR[i][j] = buf->dtrail.d.xyR[i + (n - tmpndense)][j];
         }
      }
      vectorsetlengthatleast(&buf->tmp0, 2 * n);
      vectorsetlengthatleast(&buf->tmpp, n);
      rmatrixplurec(&buf->dbuf, 0, tmpndense, tmpndense, &buf->tmpp, &buf->tmp0);
   // Convert indexes of rows pivots, swap elements of BLeft
      for (i = 0; i < tmpndense; i++) {
         pr->xZ[i + (n - tmpndense)] = buf->tmpp.xZ[i] + (n - tmpndense);
         sptrf_sluv2list1swap(&buf->bleft, i + (n - tmpndense), pr->xZ[i + (n - tmpndense)]);
         swapi(&buf->rowpermrawidx.xZ[i + (n - tmpndense)], &buf->rowpermrawidx.xZ[pr->xZ[i + (n - tmpndense)]]);
      }
   // Convert U-factor
      ivectorgrowto(&buf->sparseut.idx, buf->sparseut.ridx.xZ[n - tmpndense] + n * tmpndense);
      rvectorgrowto(&buf->sparseut.vals, buf->sparseut.ridx.xZ[n - tmpndense] + n * tmpndense);
      for (j = 0; j < tmpndense; j++) {
         offs = buf->sparseut.ridx.xZ[j + (n - tmpndense)];
         k = n - tmpndense;
      // Convert leading N-NDense columns
         for (i = 0; i < k; i++) {
            v = buf->dtrail.d.xyR[i][j];
            if (v != 0) {
               buf->sparseut.idx.xZ[offs] = i;
               buf->sparseut.vals.xR[offs] = v;
               offs++;
            }
         }
      // Convert upper diagonal elements
         for (i = 0; i < j; i++) {
            v = buf->dbuf.xyR[i][j];
            if (v != 0) {
               buf->sparseut.idx.xZ[offs] = i + (n - tmpndense);
               buf->sparseut.vals.xR[offs] = v;
               offs++;
            }
         }
      // Convert diagonal element (always stored)
         v = buf->dbuf.xyR[j][j];
         buf->sparseut.idx.xZ[offs] = j + (n - tmpndense);
         buf->sparseut.vals.xR[offs] = v;
         offs++;
         result = result && v != 0;
      // Column is done
         buf->sparseut.ridx.xZ[j + (n - tmpndense) + 1] = offs;
      }
   // Convert L-factor
      ivectorgrowto(&buf->sparsel.idx, buf->sparsel.ridx.xZ[n - tmpndense] + n * tmpndense);
      rvectorgrowto(&buf->sparsel.vals, buf->sparsel.ridx.xZ[n - tmpndense] + n * tmpndense);
      for (i = 0; i < tmpndense; i++) {
         sptrf_sluv2list1appendsequencetomatrix(&buf->bleft, i + (n - tmpndense), false, 0.0, n, &buf->sparsel, i + (n - tmpndense));
         offs = buf->sparsel.ridx.xZ[i + (n - tmpndense) + 1];
         for (j = 0; j < i; j++) {
            v = buf->dbuf.xyR[i][j];
            if (v != 0) {
               buf->sparsel.idx.xZ[offs] = j + (n - tmpndense);
               buf->sparsel.vals.xR[offs] = v;
               offs++;
            }
         }
         buf->sparsel.ridx.xZ[i + (n - tmpndense) + 1] = offs;
      }
   }
// Allocate output
   vectorsetlengthatleast(&buf->tmpi, n);
   for (i = 0; i < n; i++) {
      buf->tmpi.xZ[i] = buf->sparsel.ridx.xZ[i + 1] - buf->sparsel.ridx.xZ[i];
   }
   for (i = 0; i < n; i++) {
      i0 = buf->sparseut.ridx.xZ[i];
      i1 = buf->sparseut.ridx.xZ[i + 1] - 1;
      for (j = i0; j <= i1; j++) {
         k = buf->sparseut.idx.xZ[j];
         buf->tmpi.xZ[k]++;
      }
   }
   a->matrixtype = 1;
   a->ninitialized = buf->sparsel.ridx.xZ[n] + buf->sparseut.ridx.xZ[n];
   a->m = n;
   a->n = n;
   vectorsetlengthatleast(&a->ridx, n + 1);
   vectorsetlengthatleast(&a->idx, a->ninitialized);
   vectorsetlengthatleast(&a->vals, a->ninitialized);
   a->ridx.xZ[0] = 0;
   for (i = 0; i < n; i++) {
      a->ridx.xZ[i + 1] = a->ridx.xZ[i] + buf->tmpi.xZ[i];
   }
   for (i = 0; i < n; i++) {
      i0 = buf->sparsel.ridx.xZ[i];
      i1 = buf->sparsel.ridx.xZ[i + 1] - 1;
      jp = a->ridx.xZ[i];
      for (j = i0; j <= i1; j++) {
         a->idx.xZ[jp + (j - i0)] = buf->sparsel.idx.xZ[j];
         a->vals.xR[jp + (j - i0)] = buf->sparsel.vals.xR[j];
      }
      buf->tmpi.xZ[i] = buf->sparsel.ridx.xZ[i + 1] - buf->sparsel.ridx.xZ[i];
   }
   vectorsetlengthatleast(&a->didx, n);
   vectorsetlengthatleast(&a->uidx, n);
   for (i = 0; i < n; i++) {
      a->didx.xZ[i] = a->ridx.xZ[i] + buf->tmpi.xZ[i];
      a->uidx.xZ[i] = a->didx.xZ[i] + 1;
      buf->tmpi.xZ[i] = a->didx.xZ[i];
   }
   for (i = 0; i < n; i++) {
      i0 = buf->sparseut.ridx.xZ[i];
      i1 = buf->sparseut.ridx.xZ[i + 1] - 1;
      for (j = i0; j <= i1; j++) {
         k = buf->sparseut.idx.xZ[j];
         offs = buf->tmpi.xZ[k];
         a->idx.xZ[offs] = i;
         a->vals.xR[offs] = buf->sparseut.vals.xR[j];
         buf->tmpi.xZ[k] = offs + 1;
      }
   }
   return result;
}

void sluv2list1matrix_init(void *_p, bool make_automatic) {
   sluv2list1matrix *p = (sluv2list1matrix *)_p;
   ae_vector_init(&p->idxfirst, 0, DT_INT, make_automatic);
   ae_vector_init(&p->strgidx, 0, DT_INT, make_automatic);
   ae_vector_init(&p->strgval, 0, DT_REAL, make_automatic);
}

void sluv2list1matrix_copy(void *_dst, void *_src, bool make_automatic) {
   sluv2list1matrix *dst = (sluv2list1matrix *)_dst;
   sluv2list1matrix *src = (sluv2list1matrix *)_src;
   dst->nfixed = src->nfixed;
   dst->ndynamic = src->ndynamic;
   ae_vector_copy(&dst->idxfirst, &src->idxfirst, make_automatic);
   ae_vector_copy(&dst->strgidx, &src->strgidx, make_automatic);
   ae_vector_copy(&dst->strgval, &src->strgval, make_automatic);
   dst->nallocated = src->nallocated;
   dst->nused = src->nused;
}

void sluv2list1matrix_free(void *_p, bool make_automatic) {
   sluv2list1matrix *p = (sluv2list1matrix *)_p;
   ae_vector_free(&p->idxfirst, make_automatic);
   ae_vector_free(&p->strgidx, make_automatic);
   ae_vector_free(&p->strgval, make_automatic);
}

void sluv2sparsetrail_init(void *_p, bool make_automatic) {
   sluv2sparsetrail *p = (sluv2sparsetrail *)_p;
   ae_vector_init(&p->nzc, 0, DT_INT, make_automatic);
   ae_vector_init(&p->wrkset, 0, DT_INT, make_automatic);
   ae_vector_init(&p->colid, 0, DT_INT, make_automatic);
   ae_vector_init(&p->isdensified, 0, DT_BOOL, make_automatic);
   ae_vector_init(&p->slscolptr, 0, DT_INT, make_automatic);
   ae_vector_init(&p->slsrowptr, 0, DT_INT, make_automatic);
   ae_vector_init(&p->slsidx, 0, DT_INT, make_automatic);
   ae_vector_init(&p->slsval, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmp0, 0, DT_REAL, make_automatic);
}

void sluv2sparsetrail_copy(void *_dst, void *_src, bool make_automatic) {
   sluv2sparsetrail *dst = (sluv2sparsetrail *)_dst;
   sluv2sparsetrail *src = (sluv2sparsetrail *)_src;
   dst->n = src->n;
   dst->k = src->k;
   ae_vector_copy(&dst->nzc, &src->nzc, make_automatic);
   dst->maxwrkcnt = src->maxwrkcnt;
   dst->maxwrknz = src->maxwrknz;
   dst->wrkcnt = src->wrkcnt;
   ae_vector_copy(&dst->wrkset, &src->wrkset, make_automatic);
   ae_vector_copy(&dst->colid, &src->colid, make_automatic);
   ae_vector_copy(&dst->isdensified, &src->isdensified, make_automatic);
   ae_vector_copy(&dst->slscolptr, &src->slscolptr, make_automatic);
   ae_vector_copy(&dst->slsrowptr, &src->slsrowptr, make_automatic);
   ae_vector_copy(&dst->slsidx, &src->slsidx, make_automatic);
   ae_vector_copy(&dst->slsval, &src->slsval, make_automatic);
   dst->slsused = src->slsused;
   ae_vector_copy(&dst->tmp0, &src->tmp0, make_automatic);
}

void sluv2sparsetrail_free(void *_p, bool make_automatic) {
   sluv2sparsetrail *p = (sluv2sparsetrail *)_p;
   ae_vector_free(&p->nzc, make_automatic);
   ae_vector_free(&p->wrkset, make_automatic);
   ae_vector_free(&p->colid, make_automatic);
   ae_vector_free(&p->isdensified, make_automatic);
   ae_vector_free(&p->slscolptr, make_automatic);
   ae_vector_free(&p->slsrowptr, make_automatic);
   ae_vector_free(&p->slsidx, make_automatic);
   ae_vector_free(&p->slsval, make_automatic);
   ae_vector_free(&p->tmp0, make_automatic);
}

void sluv2densetrail_init(void *_p, bool make_automatic) {
   sluv2densetrail *p = (sluv2densetrail *)_p;
   ae_matrix_init(&p->d, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->did, 0, DT_INT, make_automatic);
}

void sluv2densetrail_copy(void *_dst, void *_src, bool make_automatic) {
   sluv2densetrail *dst = (sluv2densetrail *)_dst;
   sluv2densetrail *src = (sluv2densetrail *)_src;
   dst->n = src->n;
   dst->ndense = src->ndense;
   ae_matrix_copy(&dst->d, &src->d, make_automatic);
   ae_vector_copy(&dst->did, &src->did, make_automatic);
}

void sluv2densetrail_free(void *_p, bool make_automatic) {
   sluv2densetrail *p = (sluv2densetrail *)_p;
   ae_matrix_free(&p->d, make_automatic);
   ae_vector_free(&p->did, make_automatic);
}

void sluv2buffer_init(void *_p, bool make_automatic) {
   sluv2buffer *p = (sluv2buffer *)_p;
   sparsematrix_init(&p->sparsel, make_automatic);
   sparsematrix_init(&p->sparseut, make_automatic);
   sluv2list1matrix_init(&p->bleft, make_automatic);
   sluv2list1matrix_init(&p->bupper, make_automatic);
   sluv2sparsetrail_init(&p->strail, make_automatic);
   sluv2densetrail_init(&p->dtrail, make_automatic);
   ae_vector_init(&p->rowpermrawidx, 0, DT_INT, make_automatic);
   ae_matrix_init(&p->dbuf, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->v0i, 0, DT_INT, make_automatic);
   ae_vector_init(&p->v1i, 0, DT_INT, make_automatic);
   ae_vector_init(&p->v0r, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->v1r, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmp0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmpi, 0, DT_INT, make_automatic);
   ae_vector_init(&p->tmpp, 0, DT_INT, make_automatic);
}

void sluv2buffer_copy(void *_dst, void *_src, bool make_automatic) {
   sluv2buffer *dst = (sluv2buffer *)_dst;
   sluv2buffer *src = (sluv2buffer *)_src;
   dst->n = src->n;
   sparsematrix_copy(&dst->sparsel, &src->sparsel, make_automatic);
   sparsematrix_copy(&dst->sparseut, &src->sparseut, make_automatic);
   sluv2list1matrix_copy(&dst->bleft, &src->bleft, make_automatic);
   sluv2list1matrix_copy(&dst->bupper, &src->bupper, make_automatic);
   sluv2sparsetrail_copy(&dst->strail, &src->strail, make_automatic);
   sluv2densetrail_copy(&dst->dtrail, &src->dtrail, make_automatic);
   ae_vector_copy(&dst->rowpermrawidx, &src->rowpermrawidx, make_automatic);
   ae_matrix_copy(&dst->dbuf, &src->dbuf, make_automatic);
   ae_vector_copy(&dst->v0i, &src->v0i, make_automatic);
   ae_vector_copy(&dst->v1i, &src->v1i, make_automatic);
   ae_vector_copy(&dst->v0r, &src->v0r, make_automatic);
   ae_vector_copy(&dst->v1r, &src->v1r, make_automatic);
   ae_vector_copy(&dst->tmp0, &src->tmp0, make_automatic);
   ae_vector_copy(&dst->tmpi, &src->tmpi, make_automatic);
   ae_vector_copy(&dst->tmpp, &src->tmpp, make_automatic);
}

void sluv2buffer_free(void *_p, bool make_automatic) {
   sluv2buffer *p = (sluv2buffer *)_p;
   sparsematrix_free(&p->sparsel, make_automatic);
   sparsematrix_free(&p->sparseut, make_automatic);
   sluv2list1matrix_free(&p->bleft, make_automatic);
   sluv2list1matrix_free(&p->bupper, make_automatic);
   sluv2sparsetrail_free(&p->strail, make_automatic);
   sluv2densetrail_free(&p->dtrail, make_automatic);
   ae_vector_free(&p->rowpermrawidx, make_automatic);
   ae_matrix_free(&p->dbuf, make_automatic);
   ae_vector_free(&p->v0i, make_automatic);
   ae_vector_free(&p->v1i, make_automatic);
   ae_vector_free(&p->v0r, make_automatic);
   ae_vector_free(&p->v1r, make_automatic);
   ae_vector_free(&p->tmp0, make_automatic);
   ae_vector_free(&p->tmpi, make_automatic);
   ae_vector_free(&p->tmpp, make_automatic);
}
} // end of namespace alglib_impl

// === MATGEN Package ===
// Depends on: (AlgLibInternal) CREFLECTIONS
// Depends on: (AlgLibMisc) HQRND
// Depends on: ABLAS
namespace alglib_impl {
// Generation of a random uniformly distributed (Haar) orthogonal matrix
//
// Inputs:
//     N   -   matrix size, N >= 1
//
// Outputs:
//     A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]
//
// NOTE: this function uses algorithm  described  in  Stewart, G. W.  (1980),
//       "The Efficient Generation of  Random  Orthogonal  Matrices  with  an
//       Application to Condition Estimators".
//
//       Speaking short, to generate an (N+1)x(N+1) orthogonal matrix, it:
//       * takes an NxN one
//       * takes uniformly distributed unit vector of dimension N+1.
//       * constructs a Householder reflection from the vector, then applies
//         it to the smaller matrix (embedded in the larger size with a 1 at
//         the bottom right corner).
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void rmatrixrndorthogonal(const ae_int_t n, real_2d_array &a);
void rmatrixrndorthogonal(ae_int_t n, RMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(a);
   ae_assert(n >= 1, "RMatrixRndOrthogonal: N<1!");
   ae_matrix_set_length(a, n, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            a->xyR[i][j] = 1.0;
         } else {
            a->xyR[i][j] = 0.0;
         }
      }
   }
   rmatrixrndorthogonalfromtheright(a, n, n);
}

// Generation of random NxN matrix with given condition number and norm2(A)=1
//
// Inputs:
//     N   -   matrix size
//     C   -   condition number (in 2-norm)
//
// Outputs:
//     A   -   random matrix with norm2(A)=1 and cond(A)=C
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void rmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a);
void rmatrixrndcond(ae_int_t n, double c, RMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double l1;
   double l2;
   ae_frame_make(&_frame_block);
   SetMatrix(a);
   NewObj(hqrndstate, rs);
   ae_assert(n >= 1 && c >= 1.0, "RMatrixRndCond: N<1 or C<1!");
   ae_matrix_set_length(a, n, n);
   if (n == 1) {
   // special case
      a->xyR[0][0] = (double)(2 * ae_randominteger(2) - 1);
      ae_frame_leave();
      return;
   }
   hqrndrandomize(&rs);
   l1 = 0.0;
   l2 = log(1 / c);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a->xyR[i][j] = 0.0;
      }
   }
   a->xyR[0][0] = exp(l1);
   for (i = 1; i < n - 1; i++) {
      a->xyR[i][i] = exp(hqrnduniformr(&rs) * (l2 - l1) + l1);
   }
   a->xyR[n - 1][n - 1] = exp(l2);
   rmatrixrndorthogonalfromtheleft(a, n, n);
   rmatrixrndorthogonalfromtheright(a, n, n);
   ae_frame_leave();
}

// Generation of a random Haar distributed orthogonal complex matrix
//
// Inputs:
//     N   -   matrix size, N >= 1
//
// Outputs:
//     A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]
//
// NOTE: this function uses algorithm  described  in  Stewart, G. W.  (1980),
//       "The Efficient Generation of  Random  Orthogonal  Matrices  with  an
//       Application to Condition Estimators".
//
//       Speaking short, to generate an (N+1)x(N+1) orthogonal matrix, it:
//       * takes an NxN one
//       * takes uniformly distributed unit vector of dimension N+1.
//       * constructs a Householder reflection from the vector, then applies
//         it to the smaller matrix (embedded in the larger size with a 1 at
//         the bottom right corner).
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void cmatrixrndorthogonal(const ae_int_t n, complex_2d_array &a);
void cmatrixrndorthogonal(ae_int_t n, CMatrix *a) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(a);
   ae_assert(n >= 1, "CMatrixRndOrthogonal: N<1!");
   ae_matrix_set_length(a, n, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            a->xyC[i][j] = ae_complex_from_i(1);
         } else {
            a->xyC[i][j] = ae_complex_from_i(0);
         }
      }
   }
   cmatrixrndorthogonalfromtheright(a, n, n);
}

// Generation of random NxN complex matrix with given condition number C and
// norm2(A)=1
//
// Inputs:
//     N   -   matrix size
//     C   -   condition number (in 2-norm)
//
// Outputs:
//     A   -   random matrix with norm2(A)=1 and cond(A)=C
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void cmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a);
void cmatrixrndcond(ae_int_t n, double c, CMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double l1;
   double l2;
   ae_complex v;
   ae_frame_make(&_frame_block);
   SetMatrix(a);
   NewObj(hqrndstate, state);
   ae_assert(n >= 1 && c >= 1.0, "CMatrixRndCond: N<1 or C<1!");
   ae_matrix_set_length(a, n, n);
   if (n == 1) {
   // special case
      hqrndrandomize(&state);
      hqrndunit2(&state, &v.x, &v.y);
      a->xyC[0][0] = v;
      ae_frame_leave();
      return;
   }
   hqrndrandomize(&state);
   l1 = 0.0;
   l2 = log(1 / c);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a->xyC[i][j] = ae_complex_from_i(0);
      }
   }
   a->xyC[0][0] = ae_complex_from_d(exp(l1));
   for (i = 1; i < n - 1; i++) {
      a->xyC[i][i] = ae_complex_from_d(exp(hqrnduniformr(&state) * (l2 - l1) + l1));
   }
   a->xyC[n - 1][n - 1] = ae_complex_from_d(exp(l2));
   cmatrixrndorthogonalfromtheleft(a, n, n);
   cmatrixrndorthogonalfromtheright(a, n, n);
   ae_frame_leave();
}

// Generation of random NxN symmetric matrix with given condition number  and
// norm2(A)=1
//
// Inputs:
//     N   -   matrix size
//     C   -   condition number (in 2-norm)
//
// Outputs:
//     A   -   random matrix with norm2(A)=1 and cond(A)=C
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void smatrixrndcond(const ae_int_t n, const double c, real_2d_array &a);
void smatrixrndcond(ae_int_t n, double c, RMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double l1;
   double l2;
   ae_frame_make(&_frame_block);
   SetMatrix(a);
   NewObj(hqrndstate, rs);
   ae_assert(n >= 1 && c >= 1.0, "SMatrixRndCond: N<1 or C<1!");
   ae_matrix_set_length(a, n, n);
   if (n == 1) {
   // special case
      a->xyR[0][0] = (double)(2 * ae_randominteger(2) - 1);
      ae_frame_leave();
      return;
   }
// Prepare matrix
   hqrndrandomize(&rs);
   l1 = 0.0;
   l2 = log(1 / c);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a->xyR[i][j] = 0.0;
      }
   }
   a->xyR[0][0] = exp(l1);
   for (i = 1; i < n - 1; i++) {
      a->xyR[i][i] = (2 * hqrnduniformi(&rs, 2) - 1) * exp(hqrnduniformr(&rs) * (l2 - l1) + l1);
   }
   a->xyR[n - 1][n - 1] = exp(l2);
// Multiply
   smatrixrndmultiply(a, n);
   ae_frame_leave();
}

// Generation of random NxN symmetric positive definite matrix with given
// condition number and norm2(A)=1
//
// Inputs:
//     N   -   matrix size
//     C   -   condition number (in 2-norm)
//
// Outputs:
//     A   -   random SPD matrix with norm2(A)=1 and cond(A)=C
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void spdmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a);
void spdmatrixrndcond(ae_int_t n, double c, RMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double l1;
   double l2;
   ae_frame_make(&_frame_block);
   SetMatrix(a);
   NewObj(hqrndstate, rs);
// Special cases
   if (n <= 0 || c < 1.0) {
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(a, n, n);
   if (n == 1) {
      a->xyR[0][0] = 1.0;
      ae_frame_leave();
      return;
   }
// Prepare matrix
   hqrndrandomize(&rs);
   l1 = 0.0;
   l2 = log(1 / c);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a->xyR[i][j] = 0.0;
      }
   }
   a->xyR[0][0] = exp(l1);
   for (i = 1; i < n - 1; i++) {
      a->xyR[i][i] = exp(hqrnduniformr(&rs) * (l2 - l1) + l1);
   }
   a->xyR[n - 1][n - 1] = exp(l2);
// Multiply
   smatrixrndmultiply(a, n);
   ae_frame_leave();
}

// Generation of random NxN Hermitian matrix with given condition number  and
// norm2(A)=1
//
// Inputs:
//     N   -   matrix size
//     C   -   condition number (in 2-norm)
//
// Outputs:
//     A   -   random matrix with norm2(A)=1 and cond(A)=C
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void hmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a);
void hmatrixrndcond(ae_int_t n, double c, CMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double l1;
   double l2;
   ae_frame_make(&_frame_block);
   SetMatrix(a);
   NewObj(hqrndstate, rs);
   ae_assert(n >= 1 && c >= 1.0, "HMatrixRndCond: N<1 or C<1!");
   ae_matrix_set_length(a, n, n);
   if (n == 1) {
   // special case
      a->xyC[0][0] = ae_complex_from_i(2 * ae_randominteger(2) - 1);
      ae_frame_leave();
      return;
   }
// Prepare matrix
   hqrndrandomize(&rs);
   l1 = 0.0;
   l2 = log(1 / c);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a->xyC[i][j] = ae_complex_from_i(0);
      }
   }
   a->xyC[0][0] = ae_complex_from_d(exp(l1));
   for (i = 1; i < n - 1; i++) {
      a->xyC[i][i] = ae_complex_from_d((2 * hqrnduniformi(&rs, 2) - 1) * exp(hqrnduniformr(&rs) * (l2 - l1) + l1));
   }
   a->xyC[n - 1][n - 1] = ae_complex_from_d(exp(l2));
// Multiply
   hmatrixrndmultiply(a, n);
// post-process to ensure that matrix diagonal is real
   for (i = 0; i < n; i++) {
      a->xyC[i][i].y = 0.0;
   }
   ae_frame_leave();
}

// Generation of random NxN Hermitian positive definite matrix with given
// condition number and norm2(A)=1
//
// Inputs:
//     N   -   matrix size
//     C   -   condition number (in 2-norm)
//
// Outputs:
//     A   -   random HPD matrix with norm2(A)=1 and cond(A)=C
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void hpdmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a);
void hpdmatrixrndcond(ae_int_t n, double c, CMatrix *a) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double l1;
   double l2;
   ae_frame_make(&_frame_block);
   SetMatrix(a);
   NewObj(hqrndstate, rs);
// Special cases
   if (n <= 0 || c < 1.0) {
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(a, n, n);
   if (n == 1) {
      a->xyC[0][0] = ae_complex_from_i(1);
      ae_frame_leave();
      return;
   }
// Prepare matrix
   hqrndrandomize(&rs);
   l1 = 0.0;
   l2 = log(1 / c);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a->xyC[i][j] = ae_complex_from_i(0);
      }
   }
   a->xyC[0][0] = ae_complex_from_d(exp(l1));
   for (i = 1; i < n - 1; i++) {
      a->xyC[i][i] = ae_complex_from_d(exp(hqrnduniformr(&rs) * (l2 - l1) + l1));
   }
   a->xyC[n - 1][n - 1] = ae_complex_from_d(exp(l2));
// Multiply
   hmatrixrndmultiply(a, n);
// post-process to ensure that matrix diagonal is real
   for (i = 0; i < n; i++) {
      a->xyC[i][i].y = 0.0;
   }
   ae_frame_leave();
}

// Multiplication of MxN matrix by NxN random Haar distributed orthogonal matrix
//
// Inputs:
//     A   -   matrix, array[0..M-1, 0..N-1]
//     M, N-   matrix size
//
// Outputs:
//     A   -   A*Q, where Q is random NxN orthogonal matrix
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void rmatrixrndorthogonalfromtheright(real_2d_array &a, const ae_int_t m, const ae_int_t n);
void rmatrixrndorthogonalfromtheright(RMatrix *a, ae_int_t m, ae_int_t n) {
   ae_frame _frame_block;
   double tau;
   double lambdav;
   ae_int_t s;
   ae_int_t i;
   double u1;
   double u2;
   ae_frame_make(&_frame_block);
   NewVector(w, 0, DT_REAL);
   NewVector(v, 0, DT_REAL);
   NewObj(hqrndstate, state);
   ae_assert(n >= 1 && m >= 1, "RMatrixRndOrthogonalFromTheRight: N<1 or M < 1!");
   if (n == 1) {
   // Special case
      tau = (double)(2 * ae_randominteger(2) - 1);
      for (i = 0; i < m; i++) {
         a->xyR[i][0] *= tau;
      }
      ae_frame_leave();
      return;
   }
// General case.
// First pass.
   ae_vector_set_length(&w, m);
   ae_vector_set_length(&v, n + 1);
   hqrndrandomize(&state);
   for (s = 2; s <= n; s++) {
   // Prepare random normal v
      do {
         i = 1;
         while (i <= s) {
            hqrndnormal2(&state, &u1, &u2);
            v.xR[i] = u1;
            if (i + 1 <= s) {
               v.xR[i + 1] = u2;
            }
            i += 2;
         }
         lambdav = ae_v_dotproduct(&v.xR[1], 1, &v.xR[1], 1, s);
      } while (lambdav == 0.0);
   // Prepare and apply reflection
      generatereflection(&v, s, &tau);
      v.xR[1] = 1.0;
      applyreflectionfromtheright(a, tau, &v, 0, m - 1, n - s, n - 1, &w);
   }
// Second pass.
   for (i = 0; i < n; i++) {
      tau = (double)(2 * hqrnduniformi(&state, 2) - 1);
      ae_v_muld(&a->xyR[0][i], a->stride, m, tau);
   }
   ae_frame_leave();
}

// Multiplication of MxN matrix by MxM random Haar distributed orthogonal matrix
//
// Inputs:
//     A   -   matrix, array[0..M-1, 0..N-1]
//     M, N-   matrix size
//
// Outputs:
//     A   -   Q*A, where Q is random MxM orthogonal matrix
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void rmatrixrndorthogonalfromtheleft(real_2d_array &a, const ae_int_t m, const ae_int_t n);
void rmatrixrndorthogonalfromtheleft(RMatrix *a, ae_int_t m, ae_int_t n) {
   ae_frame _frame_block;
   double tau;
   double lambdav;
   ae_int_t s;
   ae_int_t i;
   ae_int_t j;
   double u1;
   double u2;
   ae_frame_make(&_frame_block);
   NewVector(w, 0, DT_REAL);
   NewVector(v, 0, DT_REAL);
   NewObj(hqrndstate, state);
   ae_assert(n >= 1 && m >= 1, "RMatrixRndOrthogonalFromTheRight: N<1 or M < 1!");
   if (m == 1) {
   // special case
      tau = (double)(2 * ae_randominteger(2) - 1);
      for (j = 0; j < n; j++) {
         a->xyR[0][j] *= tau;
      }
      ae_frame_leave();
      return;
   }
// General case.
// First pass.
   ae_vector_set_length(&w, n);
   ae_vector_set_length(&v, m + 1);
   hqrndrandomize(&state);
   for (s = 2; s <= m; s++) {
   // Prepare random normal v
      do {
         i = 1;
         while (i <= s) {
            hqrndnormal2(&state, &u1, &u2);
            v.xR[i] = u1;
            if (i + 1 <= s) {
               v.xR[i + 1] = u2;
            }
            i += 2;
         }
         lambdav = ae_v_dotproduct(&v.xR[1], 1, &v.xR[1], 1, s);
      } while (lambdav == 0.0);
   // Prepare and apply reflection
      generatereflection(&v, s, &tau);
      v.xR[1] = 1.0;
      applyreflectionfromtheleft(a, tau, &v, m - s, m - 1, 0, n - 1, &w);
   }
// Second pass.
   for (i = 0; i < m; i++) {
      tau = (double)(2 * hqrnduniformi(&state, 2) - 1);
      ae_v_muld(a->xyR[i], 1, n, tau);
   }
   ae_frame_leave();
}

// Multiplication of MxN complex matrix by NxN random Haar distributed
// complex orthogonal matrix
//
// Inputs:
//     A   -   matrix, array[0..M-1, 0..N-1]
//     M, N-   matrix size
//
// Outputs:
//     A   -   A*Q, where Q is random NxN orthogonal matrix
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void cmatrixrndorthogonalfromtheright(complex_2d_array &a, const ae_int_t m, const ae_int_t n);
void cmatrixrndorthogonalfromtheright(CMatrix *a, ae_int_t m, ae_int_t n) {
   ae_frame _frame_block;
   ae_complex lambdav;
   ae_complex tau;
   ae_int_t s;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(w, 0, DT_COMPLEX);
   NewVector(v, 0, DT_COMPLEX);
   NewObj(hqrndstate, state);
   ae_assert(n >= 1 && m >= 1, "CMatrixRndOrthogonalFromTheRight: N<1 or M < 1!");
   if (n == 1) {
   // Special case
      hqrndrandomize(&state);
      hqrndunit2(&state, &tau.x, &tau.y);
      for (i = 0; i < m; i++) {
         a->xyC[i][0] = ae_c_mul(a->xyC[i][0], tau);
      }
      ae_frame_leave();
      return;
   }
// General case.
// First pass.
   ae_vector_set_length(&w, m);
   ae_vector_set_length(&v, n + 1);
   hqrndrandomize(&state);
   for (s = 2; s <= n; s++) {
   // Prepare random normal v
      do {
         for (i = 1; i <= s; i++) {
            hqrndnormal2(&state, &tau.x, &tau.y);
            v.xC[i] = tau;
         }
         lambdav = ae_v_cdotproduct(&v.xC[1], 1, "N", &v.xC[1], 1, "Conj", s);
      } while (ae_c_eq_d(lambdav, 0.0));
   // Prepare and apply reflection
      complexgeneratereflection(&v, s, &tau);
      v.xC[1] = ae_complex_from_i(1);
      complexapplyreflectionfromtheright(a, tau, &v, 0, m - 1, n - s, n - 1, &w);
   }
// Second pass.
   for (i = 0; i < n; i++) {
      hqrndunit2(&state, &tau.x, &tau.y);
      ae_v_cmulc(&a->xyC[0][i], a->stride, m, tau);
   }
   ae_frame_leave();
}

// Multiplication of MxN complex matrix by MxM random Haar distributed
// complex orthogonal matrix
//
// Inputs:
//     A   -   matrix, array[0..M-1, 0..N-1]
//     M, N-   matrix size
//
// Outputs:
//     A   -   Q*A, where Q is random MxM orthogonal matrix
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void cmatrixrndorthogonalfromtheleft(complex_2d_array &a, const ae_int_t m, const ae_int_t n);
void cmatrixrndorthogonalfromtheleft(CMatrix *a, ae_int_t m, ae_int_t n) {
   ae_frame _frame_block;
   ae_complex tau;
   ae_complex lambdav;
   ae_int_t s;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   NewVector(w, 0, DT_COMPLEX);
   NewVector(v, 0, DT_COMPLEX);
   NewObj(hqrndstate, state);
   ae_assert(n >= 1 && m >= 1, "CMatrixRndOrthogonalFromTheRight: N<1 or M < 1!");
   if (m == 1) {
   // special case
      hqrndrandomize(&state);
      hqrndunit2(&state, &tau.x, &tau.y);
      for (j = 0; j < n; j++) {
         a->xyC[0][j] = ae_c_mul(a->xyC[0][j], tau);
      }
      ae_frame_leave();
      return;
   }
// General case.
// First pass.
   ae_vector_set_length(&w, n);
   ae_vector_set_length(&v, m + 1);
   hqrndrandomize(&state);
   for (s = 2; s <= m; s++) {
   // Prepare random normal v
      do {
         for (i = 1; i <= s; i++) {
            hqrndnormal2(&state, &tau.x, &tau.y);
            v.xC[i] = tau;
         }
         lambdav = ae_v_cdotproduct(&v.xC[1], 1, "N", &v.xC[1], 1, "Conj", s);
      } while (ae_c_eq_d(lambdav, 0.0));
   // Prepare and apply reflection
      complexgeneratereflection(&v, s, &tau);
      v.xC[1] = ae_complex_from_i(1);
      complexapplyreflectionfromtheleft(a, tau, &v, m - s, m - 1, 0, n - 1, &w);
   }
// Second pass.
   for (i = 0; i < m; i++) {
      hqrndunit2(&state, &tau.x, &tau.y);
      ae_v_cmulc(a->xyC[i], 1, n, tau);
   }
   ae_frame_leave();
}

// Symmetric multiplication of NxN matrix by random Haar distributed
// orthogonal  matrix
//
// Inputs:
//     A   -   matrix, array[0..N-1, 0..N-1]
//     N   -   matrix size
//
// Outputs:
//     A   -   Q'*A*Q, where Q is random NxN orthogonal matrix
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void smatrixrndmultiply(real_2d_array &a, const ae_int_t n);
void smatrixrndmultiply(RMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   double tau;
   double lambdav;
   ae_int_t s;
   ae_int_t i;
   double u1;
   double u2;
   ae_frame_make(&_frame_block);
   NewVector(w, 0, DT_REAL);
   NewVector(v, 0, DT_REAL);
   NewObj(hqrndstate, state);
// General case.
   ae_vector_set_length(&w, n);
   ae_vector_set_length(&v, n + 1);
   hqrndrandomize(&state);
   for (s = 2; s <= n; s++) {
   // Prepare random normal v
      do {
         i = 1;
         while (i <= s) {
            hqrndnormal2(&state, &u1, &u2);
            v.xR[i] = u1;
            if (i + 1 <= s) {
               v.xR[i + 1] = u2;
            }
            i += 2;
         }
         lambdav = ae_v_dotproduct(&v.xR[1], 1, &v.xR[1], 1, s);
      } while (lambdav == 0.0);
   // Prepare and apply reflection
      generatereflection(&v, s, &tau);
      v.xR[1] = 1.0;
      applyreflectionfromtheright(a, tau, &v, 0, n - 1, n - s, n - 1, &w);
      applyreflectionfromtheleft(a, tau, &v, n - s, n - 1, 0, n - 1, &w);
   }
// Second pass.
   for (i = 0; i < n; i++) {
      tau = (double)(2 * hqrnduniformi(&state, 2) - 1);
      ae_v_muld(&a->xyR[0][i], a->stride, n, tau);
      ae_v_muld(a->xyR[i], 1, n, tau);
   }
// Copy upper triangle to lower
   for (i = 0; i < n - 1; i++) {
      ae_v_move(&a->xyR[i + 1][i], a->stride, &a->xyR[i][i + 1], 1, n - i - 1);
   }
   ae_frame_leave();
}

// Hermitian multiplication of NxN matrix by random Haar distributed
// complex orthogonal matrix
//
// Inputs:
//     A   -   matrix, array[0..N-1, 0..N-1]
//     N   -   matrix size
//
// Outputs:
//     A   -   Q^H*A*Q, where Q is random NxN orthogonal matrix
//
// ALGLIB Routine: Copyright 04.12.2009 by Sergey Bochkanov
// API: void hmatrixrndmultiply(complex_2d_array &a, const ae_int_t n);
void hmatrixrndmultiply(CMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_complex tau;
   ae_complex lambdav;
   ae_int_t s;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   NewVector(w, 0, DT_COMPLEX);
   NewVector(v, 0, DT_COMPLEX);
   NewObj(hqrndstate, state);
// General case.
   ae_vector_set_length(&w, n);
   ae_vector_set_length(&v, n + 1);
   hqrndrandomize(&state);
   for (s = 2; s <= n; s++) {
   // Prepare random normal v
      do {
         for (i = 1; i <= s; i++) {
            hqrndnormal2(&state, &tau.x, &tau.y);
            v.xC[i] = tau;
         }
         lambdav = ae_v_cdotproduct(&v.xC[1], 1, "N", &v.xC[1], 1, "Conj", s);
      } while (ae_c_eq_d(lambdav, 0.0));
   // Prepare and apply reflection
      complexgeneratereflection(&v, s, &tau);
      v.xC[1] = ae_complex_from_i(1);
      complexapplyreflectionfromtheright(a, tau, &v, 0, n - 1, n - s, n - 1, &w);
      complexapplyreflectionfromtheleft(a, ae_c_conj(tau), &v, n - s, n - 1, 0, n - 1, &w);
   }
// Second pass.
   for (i = 0; i < n; i++) {
      hqrndunit2(&state, &tau.x, &tau.y);
      ae_v_cmulc(&a->xyC[0][i], a->stride, n, tau);
      tau = ae_c_conj(tau);
      ae_v_cmulc(a->xyC[i], 1, n, tau);
   }
// Change all values from lower triangle by complex-conjugate values
// from upper one
   for (i = 0; i < n - 1; i++) {
      ae_v_cmove(&a->xyC[i + 1][i], a->stride, &a->xyC[i][i + 1], 1, "N", n - i - 1);
   }
   for (s = 0; s < n - 1; s++) {
      for (i = s + 1; i < n; i++) {
         a->xyC[i][s].y = -a->xyC[i][s].y;
      }
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void rmatrixrndorthogonal(const ae_int_t n, real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixrndorthogonal(n, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void rmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixrndcond(n, c, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void cmatrixrndorthogonal(const ae_int_t n, complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixrndorthogonal(n, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void cmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixrndcond(n, c, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void smatrixrndcond(const ae_int_t n, const double c, real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::smatrixrndcond(n, c, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void spdmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixrndcond(n, c, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void hmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hmatrixrndcond(n, c, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void hpdmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixrndcond(n, c, ConstT(ae_matrix, a));
   alglib_impl::ae_state_clear();
}

void rmatrixrndorthogonalfromtheright(real_2d_array &a, const ae_int_t m, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixrndorthogonalfromtheright(ConstT(ae_matrix, a), m, n);
   alglib_impl::ae_state_clear();
}

void rmatrixrndorthogonalfromtheleft(real_2d_array &a, const ae_int_t m, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixrndorthogonalfromtheleft(ConstT(ae_matrix, a), m, n);
   alglib_impl::ae_state_clear();
}

void cmatrixrndorthogonalfromtheright(complex_2d_array &a, const ae_int_t m, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixrndorthogonalfromtheright(ConstT(ae_matrix, a), m, n);
   alglib_impl::ae_state_clear();
}

void cmatrixrndorthogonalfromtheleft(complex_2d_array &a, const ae_int_t m, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixrndorthogonalfromtheleft(ConstT(ae_matrix, a), m, n);
   alglib_impl::ae_state_clear();
}

void smatrixrndmultiply(real_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::smatrixrndmultiply(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
}

void hmatrixrndmultiply(complex_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hmatrixrndmultiply(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === TRFAC Package ===
// Depends on: (AlgLibInternal) ROTATIONS
// Depends on: SPTRF, MATGEN
namespace alglib_impl {
// LU decomposition of a general real matrix with row pivoting
//
// A is represented as A = P*L*U, where:
// * L is lower unitriangular matrix
// * U is upper triangular matrix
// * P = P0*P1*...*PK, K=min(M,N)-1,
//   Pi - permutation matrix for I and Pivots[I]
//
// Inputs:
//     A       -   array[0..M-1, 0..N-1].
//     M       -   number of rows in matrix A.
//     N       -   number of columns in matrix A.
//
//
// Outputs:
//     A       -   matrices L and U in compact form:
//                 * L is stored under main diagonal
//                 * U is stored on and above main diagonal
//     Pivots  -   permutation matrix in compact form.
//                 array[0..Min(M-1,N-1)].
//
// ALGLIB Routine: Copyright 10.01.2010 by Sergey Bochkanov
// API: void rmatrixlu(real_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots);
void rmatrixlu(RMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots) {
   SetVector(pivots);
   ae_assert(m > 0, "RMatrixLU: incorrect M!");
   ae_assert(n > 0, "RMatrixLU: incorrect N!");
   rmatrixplu(a, m, n, pivots);
}

// LU decomposition of a general complex matrix with row pivoting
//
// A is represented as A = P*L*U, where:
// * L is lower unitriangular matrix
// * U is upper triangular matrix
// * P = P0*P1*...*PK, K=min(M,N)-1,
//   Pi - permutation matrix for I and Pivots[I]
//
// Inputs:
//     A       -   array[0..M-1, 0..N-1].
//     M       -   number of rows in matrix A.
//     N       -   number of columns in matrix A.
//
//
// Outputs:
//     A       -   matrices L and U in compact form:
//                 * L is stored under main diagonal
//                 * U is stored on and above main diagonal
//     Pivots  -   permutation matrix in compact form.
//                 array[0..Min(M-1,N-1)].
//
// ALGLIB Routine: Copyright 10.01.2010 by Sergey Bochkanov
// API: void cmatrixlu(complex_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots);
void cmatrixlu(CMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots) {
   SetVector(pivots);
   ae_assert(m > 0, "CMatrixLU: incorrect M!");
   ae_assert(n > 0, "CMatrixLU: incorrect N!");
   cmatrixplu(a, m, n, pivots);
}

// Level-2 Hermitian Cholesky subroutine.
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static bool trfac_hpdmatrixcholesky2(CMatrix *aaa, ae_int_t offs, ae_int_t n, bool isupper, CVector *tmp) {
   ae_int_t i;
   ae_int_t j;
   double ajj;
   ae_complex v;
   double r;
   bool result;
   result = true;
   if (n < 0) {
      result = false;
      return result;
   }
// Quick return if possible
   if (n == 0) {
      return result;
   }
   if (isupper) {
   // Compute the Cholesky factorization A = U'*U.
      for (j = 0; j < n; j++) {
      // Compute U(J,J) and test for non-positive-definiteness.
         v = ae_v_cdotproduct(&aaa->xyC[offs][offs + j], aaa->stride, "Conj", &aaa->xyC[offs][offs + j], aaa->stride, "N", j);
         ajj = ae_c_sub(aaa->xyC[offs + j][offs + j], v).x;
         if (ajj <= 0.0) {
            aaa->xyC[offs + j][offs + j] = ae_complex_from_d(ajj);
            result = false;
            return result;
         }
         ajj = sqrt(ajj);
         aaa->xyC[offs + j][offs + j] = ae_complex_from_d(ajj);
      // Compute elements J+1:N-1 of row J.
         if (j < n - 1) {
            if (j > 0) {
               ae_v_cmoveneg(tmp->xC, 1, &aaa->xyC[offs][offs + j], aaa->stride, "Conj", j);
               cmatrixmv(n - j - 1, j, aaa, offs, offs + j + 1, 1, tmp, 0, tmp, n);
               ae_v_cadd(&aaa->xyC[offs + j][offs + j + 1], 1, &tmp->xC[n], 1, "N", n - j - 1);
            }
            r = 1 / ajj;
            ae_v_cmuld(&aaa->xyC[offs + j][offs + j + 1], 1, n - j - 1, r);
         }
      }
   } else {
   // Compute the Cholesky factorization A = L*L'.
      for (j = 0; j < n; j++) {
      // Compute L(J+1,J+1) and test for non-positive-definiteness.
         v = ae_v_cdotproduct(&aaa->xyC[offs + j][offs], 1, "Conj", &aaa->xyC[offs + j][offs], 1, "N", j);
         ajj = ae_c_sub(aaa->xyC[offs + j][offs + j], v).x;
         if (ajj <= 0.0) {
            aaa->xyC[offs + j][offs + j] = ae_complex_from_d(ajj);
            result = false;
            return result;
         }
         ajj = sqrt(ajj);
         aaa->xyC[offs + j][offs + j] = ae_complex_from_d(ajj);
      // Compute elements J+1:N of column J.
         if (j < n - 1) {
            r = 1 / ajj;
            if (j > 0) {
               ae_v_cmove(tmp->xC, 1, &aaa->xyC[offs + j][offs], 1, "Conj", j);
               cmatrixmv(n - j - 1, j, aaa, offs + j + 1, offs, 0, tmp, 0, tmp, n);
               for (i = 0; i < n - j - 1; i++) {
                  aaa->xyC[offs + j + 1 + i][offs + j] = ae_c_mul_d(ae_c_sub(aaa->xyC[offs + j + 1 + i][offs + j], tmp->xC[n + i]), r);
               }
            } else {
               for (i = 0; i < n - j - 1; i++) {
                  aaa->xyC[offs + j + 1 + i][offs + j] = ae_c_mul_d(aaa->xyC[offs + j + 1 + i][offs + j], r);
               }
            }
         }
      }
   }
   return result;
}

// Recursive computational subroutine for HPDMatrixCholesky
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
static bool trfac_hpdmatrixcholeskyrec(CMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, CVector *tmp) {
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   bool result;
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
// check N
   if (n < 1) {
      result = false;
      return result;
   }
// Prepare buffer
   if (tmp->cnt < 2 * n) {
      ae_vector_set_length(tmp, 2 * n);
   }
// Basecases
//
// NOTE: we do not use MKL for basecases because their price is only
//       minor part of overall running time for N>256.
   if (n == 1) {
      if (a->xyC[offs][offs].x > 0.0) {
         a->xyC[offs][offs] = ae_complex_from_d(sqrt(a->xyC[offs][offs].x));
         result = true;
      } else {
         result = false;
      }
      return result;
   }
   if (n <= tsa) {
      result = trfac_hpdmatrixcholesky2(a, offs, n, isupper, tmp);
      return result;
   }
// Split task into smaller ones
   if (n > tsb) {
   // Split leading B-sized block from the beginning (block-matrix approach)
      n1 = tsb;
      n2 = n - n1;
   } else {
   // Smaller than B-size, perform cache-oblivious split
      n1 = tiledsplit(n, tsa), n2 = n - n1;
   }
   result = trfac_hpdmatrixcholeskyrec(a, offs, n1, isupper, tmp);
   if (!result) {
      return result;
   }
   if (n2 > 0) {
      if (isupper) {
         cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 2, a, offs, offs + n1);
         cmatrixherk(n2, n1, -1.0, a, offs, offs + n1, 2, 1.0, a, offs + n1, offs + n1, isupper);
      } else {
         cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 2, a, offs + n1, offs);
         cmatrixherk(n2, n1, -1.0, a, offs + n1, offs, 0, 1.0, a, offs + n1, offs + n1, isupper);
      }
      result = trfac_hpdmatrixcholeskyrec(a, offs + n1, n2, isupper, tmp);
      if (!result) {
         return result;
      }
   }
   return result;
}

// Cache-oblivious Cholesky decomposition
//
// The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
// definite matrix. The result of an algorithm is a representation  of  A  as
// A=U^T*U  or A=L*L^T
//
// Inputs:
//     A       -   upper or lower triangle of a factorized matrix.
//                 array with elements [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   if IsUpper=True, then A contains an upper triangle of
//                 a symmetric matrix, otherwise A contains a lower one.
//
// Outputs:
//     A       -   the result of factorization. If IsUpper=True, then
//                 the upper triangle contains matrix U, so that A = U^T*U,
//                 and the elements below the main diagonal are not modified.
//                 Similarly, if IsUpper = False.
//
// Result:
//     If  the  matrix  is  positive-definite,  the  function  returns  True.
//     Otherwise, the function returns False. Contents of A is not determined
//     in such case.
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
// API: bool spdmatrixcholesky(real_2d_array &a, const ae_int_t n, const bool isupper);
bool spdmatrixcholesky(RMatrix *a, ae_int_t n, bool isupper) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   NewVector(tmp, 0, DT_REAL);
   if (n < 1) {
      result = false;
      ae_frame_leave();
      return result;
   }
   result = spdmatrixcholeskyrec(a, 0, n, isupper, &tmp);
   ae_frame_leave();
   return result;
}

// Cache-oblivious Cholesky decomposition
//
// The algorithm computes Cholesky decomposition  of  a  Hermitian  positive-
// definite matrix. The result of an algorithm is a representation  of  A  as
// A=U'*U  or A=L*L' (here X' denotes conj(X^T)).
//
// Inputs:
//     A       -   upper or lower triangle of a factorized matrix.
//                 array with elements [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   if IsUpper=True, then A contains an upper triangle of
//                 a symmetric matrix, otherwise A contains a lower one.
//
// Outputs:
//     A       -   the result of factorization. If IsUpper=True, then
//                 the upper triangle contains matrix U, so that A = U'*U,
//                 and the elements below the main diagonal are not modified.
//                 Similarly, if IsUpper = False.
//
// Result:
//     If  the  matrix  is  positive-definite,  the  function  returns  True.
//     Otherwise, the function returns False. Contents of A is not determined
//     in such case.
//
// ALGLIB Routine: Copyright 15.12.2009-22.01.2018 by Sergey Bochkanov
// API: bool hpdmatrixcholesky(complex_2d_array &a, const ae_int_t n, const bool isupper);
bool hpdmatrixcholesky(CMatrix *a, ae_int_t n, bool isupper) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   NewVector(tmp, 0, DT_COMPLEX);
   if (n < 1) {
      result = false;
      ae_frame_leave();
      return result;
   }
   result = trfac_hpdmatrixcholeskyrec(a, 0, n, isupper, &tmp);
   ae_frame_leave();
   return result;
}

// Update of Cholesky decomposition: rank-1 update to original A.  "Buffered"
// version which uses preallocated buffer which is saved  between  subsequent
// function calls.
//
// This function uses internally allocated buffer which is not saved  between
// subsequent  calls.  So,  if  you  perform  a lot  of  subsequent  updates,
// we  recommend   you   to   use   "buffered"   version   of  this function:
// SPDMatrixCholeskyUpdateAdd1Buf().
//
// Inputs:
//     A       -   upper or lower Cholesky factor.
//                 array with elements [0..N-1, 0..N-1].
//                 Exception is thrown if array size is too small.
//     N       -   size of matrix A, N > 0
//     IsUpper -   if IsUpper=True, then A contains  upper  Cholesky  factor;
//                 otherwise A contains a lower one.
//     U       -   array[N], rank-1 update to A: A_mod = A + u*u'
//                 Exception is thrown if array size is too small.
//     BufR    -   possibly preallocated  buffer;  automatically  resized  if
//                 needed. It is recommended to  reuse  this  buffer  if  you
//                 perform a lot of subsequent decompositions.
//
// Outputs:
//     A       -   updated factorization.  If  IsUpper=True,  then  the  upper
//                 triangle contains matrix U, and the elements below the main
//                 diagonal are not modified. Similarly, if IsUpper = False.
//
// NOTE: this function always succeeds, so it does not return completion code
//
// NOTE: this function checks sizes of input arrays, but it does  NOT  checks
//       for presence of infinities or NAN's.
// ALGLIB: Copyright 03.02.2014 by Sergey Bochkanov
// API: void spdmatrixcholeskyupdateadd1(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u);
void spdmatrixcholeskyupdateadd1(RMatrix *a, ae_int_t n, bool isupper, RVector *u) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   NewVector(bufr, 0, DT_REAL);
   ae_assert(n > 0, "SPDMatrixCholeskyUpdateAdd1: N <= 0");
   ae_assert(a->rows >= n, "SPDMatrixCholeskyUpdateAdd1: Rows(A)<N");
   ae_assert(a->cols >= n, "SPDMatrixCholeskyUpdateAdd1: Cols(A)<N");
   ae_assert(u->cnt >= n, "SPDMatrixCholeskyUpdateAdd1: Length(U)<N");
   spdmatrixcholeskyupdateadd1buf(a, n, isupper, u, &bufr);
   ae_frame_leave();
}

// Update of Cholesky decomposition: "fixing" some variables.
//
// This function uses internally allocated buffer which is not saved  between
// subsequent  calls.  So,  if  you  perform  a lot  of  subsequent  updates,
// we  recommend   you   to   use   "buffered"   version   of  this function:
// SPDMatrixCholeskyUpdateFixBuf().
//
// "FIXING" EXPLAINED:
//
//     Suppose we have N*N positive definite matrix A. "Fixing" some variable
//     means filling corresponding row/column of  A  by  zeros,  and  setting
//     diagonal element to 1.
//
//     For example, if we fix 2nd variable in 4*4 matrix A, it becomes Af:
//
//         ( A00  A01  A02  A03 )      ( Af00  0   Af02 Af03 )
//         ( A10  A11  A12  A13 )      (  0    1    0    0   )
//         ( A20  A21  A22  A23 )  =>  ( Af20  0   Af22 Af23 )
//         ( A30  A31  A32  A33 )      ( Af30  0   Af32 Af33 )
//
//     If we have Cholesky decomposition of A, it must be recalculated  after
//     variables were  fixed.  However,  it  is  possible  to  use  efficient
//     algorithm, which needs O(K*N^2)  time  to  "fix"  K  variables,  given
//     Cholesky decomposition of original, "unfixed" A.
//
// Inputs:
//     A       -   upper or lower Cholesky factor.
//                 array with elements [0..N-1, 0..N-1].
//                 Exception is thrown if array size is too small.
//     N       -   size of matrix A, N > 0
//     IsUpper -   if IsUpper=True, then A contains  upper  Cholesky  factor;
//                 otherwise A contains a lower one.
//     Fix     -   array[N], I-th element is True if I-th  variable  must  be
//                 fixed. Exception is thrown if array size is too small.
//     BufR    -   possibly preallocated  buffer;  automatically  resized  if
//                 needed. It is recommended to  reuse  this  buffer  if  you
//                 perform a lot of subsequent decompositions.
//
// Outputs:
//     A       -   updated factorization.  If  IsUpper=True,  then  the  upper
//                 triangle contains matrix U, and the elements below the main
//                 diagonal are not modified. Similarly, if IsUpper = False.
//
// NOTE: this function always succeeds, so it does not return completion code
//
// NOTE: this function checks sizes of input arrays, but it does  NOT  checks
//       for presence of infinities or NAN's.
//
// NOTE: this  function  is  efficient  only  for  moderate amount of updated
//       variables - say, 0.1*N or 0.3*N. For larger amount of  variables  it
//       will  still  work,  but  you  may  get   better   performance   with
//       straightforward Cholesky.
// ALGLIB: Copyright 03.02.2014 by Sergey Bochkanov
// API: void spdmatrixcholeskyupdatefix(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix);
void spdmatrixcholeskyupdatefix(RMatrix *a, ae_int_t n, bool isupper, BVector *fix) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   NewVector(bufr, 0, DT_REAL);
   ae_assert(n > 0, "SPDMatrixCholeskyUpdateFix: N <= 0");
   ae_assert(a->rows >= n, "SPDMatrixCholeskyUpdateFix: Rows(A)<N");
   ae_assert(a->cols >= n, "SPDMatrixCholeskyUpdateFix: Cols(A)<N");
   ae_assert(fix->cnt >= n, "SPDMatrixCholeskyUpdateFix: Length(Fix)<N");
   spdmatrixcholeskyupdatefixbuf(a, n, isupper, fix, &bufr);
   ae_frame_leave();
}

// Update of Cholesky decomposition: rank-1 update to original A.  "Buffered"
// version which uses preallocated buffer which is saved  between  subsequent
// function calls.
//
// See comments for SPDMatrixCholeskyUpdateAdd1() for more information.
//
// Inputs:
//     A       -   upper or lower Cholesky factor.
//                 array with elements [0..N-1, 0..N-1].
//                 Exception is thrown if array size is too small.
//     N       -   size of matrix A, N > 0
//     IsUpper -   if IsUpper=True, then A contains  upper  Cholesky  factor;
//                 otherwise A contains a lower one.
//     U       -   array[N], rank-1 update to A: A_mod = A + u*u'
//                 Exception is thrown if array size is too small.
//     BufR    -   possibly preallocated  buffer;  automatically  resized  if
//                 needed. It is recommended to  reuse  this  buffer  if  you
//                 perform a lot of subsequent decompositions.
//
// Outputs:
//     A       -   updated factorization.  If  IsUpper=True,  then  the  upper
//                 triangle contains matrix U, and the elements below the main
//                 diagonal are not modified. Similarly, if IsUpper = False.
// ALGLIB: Copyright 03.02.2014 by Sergey Bochkanov
// API: void spdmatrixcholeskyupdateadd1buf(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u, real_1d_array &bufr);
void spdmatrixcholeskyupdateadd1buf(RMatrix *a, ae_int_t n, bool isupper, RVector *u, RVector *bufr) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t nz;
   double cs;
   double sn;
   double v;
   double vv;
   ae_assert(n > 0, "SPDMatrixCholeskyUpdateAdd1Buf: N <= 0");
   ae_assert(a->rows >= n, "SPDMatrixCholeskyUpdateAdd1Buf: Rows(A)<N");
   ae_assert(a->cols >= n, "SPDMatrixCholeskyUpdateAdd1Buf: Cols(A)<N");
   ae_assert(u->cnt >= n, "SPDMatrixCholeskyUpdateAdd1Buf: Length(U)<N");
// Find index of first non-zero entry in U
   nz = n;
   for (i = 0; i < n; i++) {
      if (u->xR[i] != 0.0) {
         nz = i;
         break;
      }
   }
   if (nz == n) {
   // Nothing to update
      return;
   }
// If working with upper triangular matrix
   if (isupper) {
   // Perform a sequence of updates which fix variables one by one.
   // This approach is different from one which is used when we work
   // with lower triangular matrix.
      vectorsetlengthatleast(bufr, n);
      for (j = nz; j < n; j++) {
         bufr->xR[j] = u->xR[j];
      }
      for (i = nz; i < n; i++) {
         if (bufr->xR[i] != 0.0) {
            generaterotation(a->xyR[i][i], bufr->xR[i], &cs, &sn, &v);
            a->xyR[i][i] = v;
            bufr->xR[i] = 0.0;
            for (j = i + 1; j < n; j++) {
               v = a->xyR[i][j];
               vv = bufr->xR[j];
               a->xyR[i][j] = cs * v + sn * vv;
               bufr->xR[j] = -sn * v + cs * vv;
            }
         }
      }
   } else {
   // Calculate rows of modified Cholesky factor, row-by-row
   // (updates performed during variable fixing are applied
   // simultaneously to each row)
      vectorsetlengthatleast(bufr, 3 * n);
      for (j = nz; j < n; j++) {
         bufr->xR[j] = u->xR[j];
      }
      for (i = nz; i < n; i++) {
      // Update all previous updates [Idx+1...I-1] to I-th row
         vv = bufr->xR[i];
         for (j = nz; j < i; j++) {
            cs = bufr->xR[n + 2 * j];
            sn = bufr->xR[n + 2 * j + 1];
            v = a->xyR[i][j];
            a->xyR[i][j] = cs * v + sn * vv;
            vv = -sn * v + cs * vv;
         }
      // generate rotation applied to I-th element of update vector
         generaterotation(a->xyR[i][i], vv, &cs, &sn, &v);
         a->xyR[i][i] = v;
         bufr->xR[n + 2 * i] = cs;
         bufr->xR[n + 2 * i + 1] = sn;
      }
   }
}

// Update of Cholesky  decomposition:  "fixing"  some  variables.  "Buffered"
// version which uses preallocated buffer which is saved  between  subsequent
// function calls.
//
// See comments for SPDMatrixCholeskyUpdateFix() for more information.
//
// Inputs:
//     A       -   upper or lower Cholesky factor.
//                 array with elements [0..N-1, 0..N-1].
//                 Exception is thrown if array size is too small.
//     N       -   size of matrix A, N > 0
//     IsUpper -   if IsUpper=True, then A contains  upper  Cholesky  factor;
//                 otherwise A contains a lower one.
//     Fix     -   array[N], I-th element is True if I-th  variable  must  be
//                 fixed. Exception is thrown if array size is too small.
//     BufR    -   possibly preallocated  buffer;  automatically  resized  if
//                 needed. It is recommended to  reuse  this  buffer  if  you
//                 perform a lot of subsequent decompositions.
//
// Outputs:
//     A       -   updated factorization.  If  IsUpper=True,  then  the  upper
//                 triangle contains matrix U, and the elements below the main
//                 diagonal are not modified. Similarly, if IsUpper = False.
// ALGLIB: Copyright 03.02.2014 by Sergey Bochkanov
// API: void spdmatrixcholeskyupdatefixbuf(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix, real_1d_array &bufr);
void spdmatrixcholeskyupdatefixbuf(RMatrix *a, ae_int_t n, bool isupper, BVector *fix, RVector *bufr) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t nfix;
   ae_int_t idx;
   double cs;
   double sn;
   double v;
   double vv;
   ae_assert(n > 0, "SPDMatrixCholeskyUpdateFixBuf: N <= 0");
   ae_assert(a->rows >= n, "SPDMatrixCholeskyUpdateFixBuf: Rows(A)<N");
   ae_assert(a->cols >= n, "SPDMatrixCholeskyUpdateFixBuf: Cols(A)<N");
   ae_assert(fix->cnt >= n, "SPDMatrixCholeskyUpdateFixBuf: Length(Fix)<N");
// Count number of variables to fix.
// Quick exit if NFix=0 or NFix=N
   nfix = 0;
   for (i = 0; i < n; i++) {
      if (fix->xB[i]) {
         nfix++;
      }
   }
   if (nfix == 0) {
   // Nothing to fix
      return;
   }
   if (nfix == n) {
   // All variables are fixed.
   // Set A to identity and exit.
      if (isupper) {
         for (i = 0; i < n; i++) {
            a->xyR[i][i] = 1.0;
            for (j = i + 1; j < n; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
      } else {
         for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
               a->xyR[i][j] = 0.0;
            }
            a->xyR[i][i] = 1.0;
         }
      }
      return;
   }
// If working with upper triangular matrix
   if (isupper) {
   // Perform a sequence of updates which fix variables one by one.
   // This approach is different from one which is used when we work
   // with lower triangular matrix.
      vectorsetlengthatleast(bufr, n);
      for (k = 0; k < n; k++) {
         if (fix->xB[k]) {
            idx = k;
         // Quick exit if it is last variable
            if (idx == n - 1) {
               for (i = 0; i < idx; i++) {
                  a->xyR[i][idx] = 0.0;
               }
               a->xyR[idx][idx] = 1.0;
               continue;
            }
         // We have Cholesky decomposition of quadratic term in A,
         // with upper triangle being stored as given below:
         //
         //         ( U00 u01 U02 )
         //     U = (     u11 u12 )
         //         (         U22 )
         //
         // Here u11 is diagonal element corresponding to variable K. We
         // want to fix this variable, and we do so by modifying U as follows:
         //
         //             ( U00  0  U02 )
         //     U_mod = (      1   0  )
         //             (         U_m )
         //
         // with U_m = CHOLESKY [ (U22^T)*U22 + (u12^T)*u12 ]
         //
         // Of course, we can calculate U_m by calculating (U22^T)*U22 explicitly,
         // modifying it and performing Cholesky decomposition of modified matrix.
         // However, we can treat it as follows:
         // * we already have CHOLESKY[(U22^T)*U22], which is equal to U22
         // * we have rank-1 update (u12^T)*u12 applied to (U22^T)*U22
         // * thus, we can calculate updated Cholesky with O(N^2) algorithm
         //   instead of O(N^3) one
            for (j = idx + 1; j < n; j++) {
               bufr->xR[j] = a->xyR[idx][j];
            }
            for (i = 0; i < idx; i++) {
               a->xyR[i][idx] = 0.0;
            }
            a->xyR[idx][idx] = 1.0;
            for (i = idx + 1; i < n; i++) {
               a->xyR[idx][i] = 0.0;
            }
            for (i = idx + 1; i < n; i++) {
               if (bufr->xR[i] != 0.0) {
                  generaterotation(a->xyR[i][i], bufr->xR[i], &cs, &sn, &v);
                  a->xyR[i][i] = v;
                  bufr->xR[i] = 0.0;
                  for (j = i + 1; j < n; j++) {
                     v = a->xyR[i][j];
                     vv = bufr->xR[j];
                     a->xyR[i][j] = cs * v + sn * vv;
                     bufr->xR[j] = -sn * v + cs * vv;
                  }
               }
            }
         }
      }
   } else {
   // Calculate rows of modified Cholesky factor, row-by-row
   // (updates performed during variable fixing are applied
   // simultaneously to each row)
      vectorsetlengthatleast(bufr, 3 * n);
      for (k = 0; k < n; k++) {
         if (fix->xB[k]) {
            idx = k;
         // Quick exit if it is last variable
            if (idx == n - 1) {
               for (i = 0; i < idx; i++) {
                  a->xyR[idx][i] = 0.0;
               }
               a->xyR[idx][idx] = 1.0;
               continue;
            }
         // store column to buffer and clear row/column of A
            for (j = idx + 1; j < n; j++) {
               bufr->xR[j] = a->xyR[j][idx];
            }
            for (i = 0; i < idx; i++) {
               a->xyR[idx][i] = 0.0;
            }
            a->xyR[idx][idx] = 1.0;
            for (i = idx + 1; i < n; i++) {
               a->xyR[i][idx] = 0.0;
            }
         // Apply update to rows of A
            for (i = idx + 1; i < n; i++) {
            // Update all previous updates [Idx+1...I-1] to I-th row
               vv = bufr->xR[i];
               for (j = idx + 1; j < i; j++) {
                  cs = bufr->xR[n + 2 * j];
                  sn = bufr->xR[n + 2 * j + 1];
                  v = a->xyR[i][j];
                  a->xyR[i][j] = cs * v + sn * vv;
                  vv = -sn * v + cs * vv;
               }
            // generate rotation applied to I-th element of update vector
               generaterotation(a->xyR[i][i], vv, &cs, &sn, &v);
               a->xyR[i][i] = v;
               bufr->xR[n + 2 * i] = cs;
               bufr->xR[n + 2 * i + 1] = sn;
            }
         }
      }
   }
}

// Sparse LU decomposition with column pivoting for sparsity and row pivoting
// for stability. Input must be square sparse matrix stored in CRS format.
//
// The algorithm  computes  LU  decomposition  of  a  general  square  matrix
// (rectangular ones are not supported). The result  of  an  algorithm  is  a
// representation of A as A = P*L*U*Q, where:
// * L is lower unitriangular matrix
// * U is upper triangular matrix
// * P = P0*P1*...*PK, K=N-1, Pi - permutation matrix for I and P[I]
// * Q = QK*...*Q1*Q0, K=N-1, Qi - permutation matrix for I and Q[I]
//
// This function pivots columns for higher sparsity, and then pivots rows for
// stability (larger element at the diagonal).
//
// Inputs:
//     A       -   sparse NxN matrix in CRS format. An exception is generated
//                 if matrix is non-CRS or non-square.
//     PivotType-  pivoting strategy:
//                 * 0 for best pivoting available (2 in current version)
//                 * 1 for row-only pivoting (NOT RECOMMENDED)
//                 * 2 for complete pivoting which produces most sparse outputs
//
// Outputs:
//     A       -   the result of factorization, matrices L and U stored in
//                 compact form using CRS sparse storage format:
//                 * lower unitriangular L is stored strictly under main diagonal
//                 * upper triangilar U is stored ON and ABOVE main diagonal
//     P       -   row permutation matrix in compact form, array[N]
//     Q       -   col permutation matrix in compact form, array[N]
//
// This function always succeeds, i.e. it ALWAYS returns valid factorization,
// but for your convenience it also returns  boolean  value  which  helps  to
// detect symbolically degenerate matrices:
// * function returns TRUE, if the matrix was factorized AND symbolically
//   non-degenerate
// * function returns FALSE, if the matrix was factorized but U has strictly
//   zero elements at the diagonal (the factorization is returned anyway).
//
//
// ALGLIB Routine: Copyright 03.09.2018 by Sergey Bochkanov
// API: bool sparselu(const sparsematrix &a, const ae_int_t pivottype, integer_1d_array &p, integer_1d_array &q);
bool sparselu(sparsematrix *a, ae_int_t pivottype, ZVector *p, ZVector *q) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   SetVector(p);
   SetVector(q);
   NewObj(sluv2buffer, buf2);
   ae_assert(pivottype == 0 || pivottype == 1 || pivottype == 2, "SparseLU: unexpected pivot type");
   ae_assert(sparseiscrs(a), "SparseLU: A is not stored in CRS format");
   ae_assert(sparsegetnrows(a) == sparsegetncols(a), "SparseLU: non-square A");
   result = sptrflu(a, pivottype, p, q, &buf2);
   ae_frame_leave();
   return result;
}

// Sparse Cholesky decomposition for skyline matrix using in-place algorithm
// without allocating additional storage.
//
// The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
// definite sparse matrix. The result of an algorithm is a representation  of
// A as A=U^T*U or A=L*L^T
//
// This  function  is  a  more  efficient alternative to general, but  slower
// SparseCholeskyX(), because it does not  create  temporary  copies  of  the
// target. It performs factorization in-place, which gives  best  performance
// on low-profile matrices. Its drawback, however, is that it can not perform
// profile-reducing permutation of input matrix.
//
// Inputs:
//     A       -   sparse matrix in skyline storage (SKS) format.
//     N       -   size of matrix A (can be smaller than actual size of A)
//     IsUpper -   if IsUpper=True, then factorization is performed on  upper
//                 triangle. Another triangle is ignored (it may contant some
//                 data, but it is not changed).
//
//
// Outputs:
//     A       -   the result of factorization, stored in SKS. If IsUpper=True,
//                 then the upper  triangle  contains  matrix  U,  such  that
//                 A = U^T*U. Lower triangle is not changed.
//                 Similarly, if IsUpper = False. In this case L is returned,
//                 and we have A = L*(L^T).
//                 Note that THIS function does not  perform  permutation  of
//                 rows to reduce bandwidth.
//
// Result:
//     If  the  matrix  is  positive-definite,  the  function  returns  True.
//     Otherwise, the function returns False. Contents of A is not determined
//     in such case.
//
// NOTE: for  performance  reasons  this  function  does NOT check that input
//       matrix  includes  only  finite  values. It is your responsibility to
//       make sure that there are no infinite or NAN values in the matrix.
//
// ALGLIB Routine: Copyright 16.01.2014 by Sergey Bochkanov
// API: bool sparsecholeskyskyline(const sparsematrix &a, const ae_int_t n, const bool isupper);
bool sparsecholeskyskyline(sparsematrix *a, ae_int_t n, bool isupper) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t jnz;
   ae_int_t jnza;
   ae_int_t jnzl;
   double v;
   double vv;
   double a12;
   ae_int_t nready;
   ae_int_t nadd;
   ae_int_t banda;
   ae_int_t offsa;
   ae_int_t offsl;
   bool result;
   ae_assert(n >= 0, "SparseCholeskySkyline: N<0");
   ae_assert(sparsegetnrows(a) >= n, "SparseCholeskySkyline: rows(A)<N");
   ae_assert(sparsegetncols(a) >= n, "SparseCholeskySkyline: cols(A)<N");
   ae_assert(sparseissks(a), "SparseCholeskySkyline: A is not stored in SKS format");
   result = false;
// transpose if needed
   if (isupper) {
      sparsetransposesks(a);
   }
// Perform Cholesky decomposition:
// * we assume than leading NReady*NReady submatrix is done
// * having Cholesky decomposition of NReady*NReady submatrix we
//   obtain decomposition of larger (NReady+NAdd)*(NReady+NAdd) one.
//
// Here is algorithm. At the start we have
//
//     (      |   )
//     (  L   |   )
// S = (      |   )
//     (----------)
//     (  A   | B )
//
// with L being already computed Cholesky factor, A and B being
// unprocessed parts of the matrix. Of course, L/A/B are stored
// in SKS format.
//
// Then, we calculate A1:=(inv(L)*A')' and replace A with A1.
// Then, we calculate B1:=B-A1*A1'     and replace B with B1
//
// Finally, we calculate small NAdd*NAdd Cholesky of B1 with
// dense solver. Now, L/A1/B1 are Cholesky decomposition of the
// larger (NReady+NAdd)*(NReady+NAdd) matrix.
   nready = 0;
   nadd = 1;
   while (nready < n) {
      ae_assert(nadd == 1, "SkylineCholesky: internal error");
   // Calculate A1:=(inv(L)*A')'
   //
   // Elements are calculated row by row (example below is given
   // for NAdd=1):
   // * first, we solve L[0,0]*A1[0]=A[0]
   // * then, we solve  L[1,0]*A1[0]+L[1,1]*A1[1]=A[1]
   // * then, we move to next row and so on
   // * during calculation of A1 we update A12 - squared norm of A1
   //
   // We extensively use sparsity of both A/A1 and L:
   // * first, equations from 0 to BANDWIDTH(A1)-1 are completely zero
   // * second, for I >= BANDWIDTH(A1), I-th equation is reduced from
   //     L[I,0]*A1[0] + L[I,1]*A1[1] + ... + L[I,I]*A1[I] = A[I]
   //   to
   //     L[I,JNZ]*A1[JNZ] + ... + L[I,I]*A1[I] = A[I]
   //   where JNZ = max(NReady-BANDWIDTH(A1),I-BANDWIDTH(L[i]))
   //   (JNZ is an index of the firts column where both A and L become
   //   nonzero).
   //
   // NOTE: we rely on details of SparseMatrix internal storage format.
   //       This is allowed by SparseMatrix specification.
      a12 = 0.0;
      if (a->didx.xZ[nready] > 0) {
         banda = a->didx.xZ[nready];
         for (i = nready - banda; i < nready; i++) {
         // Elements of A1[0:I-1] were computed:
         // * A1[0:NReady-BandA-1] are zero (sparse)
         // * A1[NReady-BandA:I-1] replaced corresponding elements of A
         //
         // Now it is time to get I-th one.
         //
         // First, we calculate:
         // * JNZA  - index of the first column where A become nonzero
         // * JNZL  - index of the first column where L become nonzero
         // * JNZ   - index of the first column where both A and L become nonzero
         // * OffsA - offset of A[JNZ] in A.Vals
         // * OffsL - offset of L[I,JNZ] in A.Vals
         //
         // Then, we solve SUM(A1[j]*L[I,j],j=JNZ..I-1) + A1[I]*L[I,I] = A[I],
         // with A1[JNZ..I-1] already known, and A1[I] unknown.
            jnza = nready - banda;
            jnzl = i - a->didx.xZ[i];
            jnz = imax2(jnza, jnzl);
            offsa = a->ridx.xZ[nready] + (jnz - jnza);
            offsl = a->ridx.xZ[i] + (jnz - jnzl);
            v = 0.0;
            k = i - 1 - jnz;
            for (j = 0; j <= k; j++) {
               v += a->vals.xR[offsa + j] * a->vals.xR[offsl + j];
            }
            vv = (a->vals.xR[offsa + k + 1] - v) / a->vals.xR[offsl + k + 1];
            a->vals.xR[offsa + k + 1] = vv;
            a12 += vv * vv;
         }
      }
   // Calculate CHOLESKY(B-A1*A1')
      offsa = a->ridx.xZ[nready] + a->didx.xZ[nready];
      v = a->vals.xR[offsa];
      if (v <= a12) {
         result = false;
         return result;
      }
      a->vals.xR[offsa] = sqrt(v - a12);
   // Increase size of the updated matrix
      nready++;
   }
// transpose if needed
   if (isupper) {
      sparsetransposesks(a);
   }
   result = true;
   return result;
}

// Sparse Cholesky decomposition: "expert" function.
//
// The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
// definite sparse matrix. The result is representation of A  as  A=U^T*U  or
// A=L*L^T
//
// Triangular factor L or U is written to separate SparseMatrix structure. If
// output buffer already contrains enough memory to store L/U, this memory is
// reused.
//
// Inputs:
//     A       -   upper or lower triangle of sparse matrix.
//                 Matrix can be in any sparse storage format.
//     N       -   size of matrix A (can be smaller than actual size of A)
//     IsUpper -   if IsUpper=True, then A contains an upper triangle of
//                 a symmetric matrix, otherwise A contains a lower one.
//                 Another triangle is ignored.
//     P0, P1  -   integer arrays:
//                 * for Ordering=-3  -  user-supplied permutation  of  rows/
//                   columns, which complies  to  requirements stated  in the
//                   "Outputs:" section.  Both  P0 and  P1  must  be
//                   initialized by user.
//                 * for other values of  Ordering  -  possibly  preallocated
//                   buffer,  which   is   filled   by  internally  generated
//                   permutation. Automatically resized if its  size  is  too
//                   small to store data.
//     Ordering-   sparse matrix reordering algorithm which is used to reduce
//                 fill-in amount:
//                 * -3    use ordering supplied by user in P0/P1
//                 * -2    use random ordering
//                 * -1    use original order
//                 * 0     use best algorithm implemented so far
//                 If input matrix is  given  in  SKS  format,  factorization
//                 function ignores Ordering and uses original order  of  the
//                 columns. The idea is that if you already store  matrix  in
//                 SKS format, it is better not to perform costly reordering.
//     Algo    -   type of algorithm which is used during factorization:
//                 * 0     use best  algorithm  (for  SKS  input  or   output
//                         matrices Algo=2 is used; otherwise Algo=1 is used)
//                 * 1     use CRS-based algorithm
//                 * 2     use skyline-based factorization algorithm.
//                         This algorithm is a  fastest  one  for low-profile
//                         matrices,  but  requires  too  much of memory  for
//                         matrices with large bandwidth.
//     Fmt     -   desired storage format  of  the  output,  as  returned  by
//                 SparseGetMatrixType() function:
//                 * 0 for hash-based storage
//                 * 1 for CRS
//                 * 2 for SKS
//                 If you do not know what format to choose, use 1 (CRS).
//     Buf     -   SparseBuffers structure which is used to store temporaries.
//                 This function may reuse previously allocated  storage,  so
//                 if you perform repeated factorizations it is beneficial to
//                 reuse Buf.
//     C       -   SparseMatrix structure  which  can  be  just  some  random
//                 garbage. In  case  in  contains  enough  memory  to  store
//                 triangular factors, this memory will be reused. Othwerwise,
//                 algorithm will automatically allocate enough memory.
//
//
// Outputs:
//     C       -   the result of factorization, stored in desired format.  If
//                 IsUpper=True, then the upper triangle  contains  matrix U,
//                 such  that  (P'*A*P) = U^T*U,  where  P  is  a permutation
//                 matrix (see below). The elements below the  main  diagonal
//                 are zero.
//                 Similarly, if IsUpper = False. In this case L is returned,
//                 and we have (P'*A*P) = L*(L^T).
//     P0      -   permutation  (according   to   Ordering  parameter)  which
//                 minimizes amount of fill-in:
//                 * P0 is array[N]
//                 * permutation is applied to A before  factorization  takes
//                   place, i.e. we have U'*U = L*L' = P'*A*P
//                 * P0[k]=j means that column/row j of A  is  moved  to k-th
//                   position before starting factorization.
//     P1      -   permutation P in another format, array[N]:
//                 * P1[k]=j means that k-th column/row of A is moved to j-th
//                   position
//
// Result:
//     If  the  matrix  is  positive-definite,  the  function  returns  True.
//     Otherwise, the function returns False. Contents of C is not determined
//     in such case.
//
// NOTE: for  performance  reasons  this  function  does NOT check that input
//       matrix  includes  only  finite  values. It is your responsibility to
//       make sure that there are no infinite or NAN values in the matrix.
//
// ALGLIB Routine: Copyright 16.01.2014 by Sergey Bochkanov
bool sparsecholeskyx(sparsematrix *a, ae_int_t n, bool isupper, ZVector *p0, ZVector *p1, ae_int_t ordering, ae_int_t algo, ae_int_t fmt, sparsebuffers *buf, sparsematrix *c) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t t0;
   ae_int_t t1;
   double v;
   bool result;
   ae_frame_make(&_frame_block);
   NewObj(hqrndstate, rs);
   ae_assert(n >= 0, "SparseMatrixCholeskyBuf: N<0");
   ae_assert(sparsegetnrows(a) >= n, "SparseMatrixCholeskyBuf: rows(A)<N");
   ae_assert(sparsegetncols(a) >= n, "SparseMatrixCholeskyBuf: cols(A)<N");
   ae_assert(ordering >= -3 && ordering <= 0, "SparseMatrixCholeskyBuf: invalid Ordering parameter");
   ae_assert(algo >= 0 && algo <= 2, "SparseMatrixCholeskyBuf: invalid Algo parameter");
   hqrndrandomize(&rs);
// Perform some quick checks.
// Because sparse matrices are expensive data structures, these
// checks are better to perform during early stages of the factorization.
   result = false;
   if (n < 1) {
      ae_frame_leave();
      return result;
   }
   for (i = 0; i < n; i++) {
      if (sparsegetdiagonal(a, i) <= 0.0) {
         ae_frame_leave();
         return result;
      }
   }
// First, determine appropriate ordering:
// * for SKS inputs, Ordering=-1 is automatically chosen (overrides user settings)
   if (ordering == 0) {
      ordering = -1;
   }
   if (sparseissks(a)) {
      ordering = -1;
   }
   if (ordering == -3) {
   // User-supplied ordering.
   // Check its correctness.
      ae_assert(p0->cnt >= n, "SparseCholeskyX: user-supplied permutation is too short");
      ae_assert(p1->cnt >= n, "SparseCholeskyX: user-supplied permutation is too short");
      for (i = 0; i < n; i++) {
         ae_assert(p0->xZ[i] >= 0 && p0->xZ[i] < n, "SparseCholeskyX: user-supplied permutation includes values outside of [0,N)");
         ae_assert(p1->xZ[i] >= 0 && p1->xZ[i] < n, "SparseCholeskyX: user-supplied permutation includes values outside of [0,N)");
         ae_assert(p1->xZ[p0->xZ[i]] == i, "SparseCholeskyX: user-supplied permutation is inconsistent - P1 is not inverse of P0");
      }
   }
   if (ordering == -2) {
   // Use random ordering
      vectorsetlengthatleast(p0, n);
      vectorsetlengthatleast(p1, n);
      for (i = 0; i < n; i++) {
         p0->xZ[i] = i;
      }
      for (i = 0; i < n; i++) {
         j = i + hqrnduniformi(&rs, n - i);
         if (j != i) {
            swapi(&p0->xZ[i], &p0->xZ[j]);
         }
      }
      for (i = 0; i < n; i++) {
         p1->xZ[p0->xZ[i]] = i;
      }
   }
   if (ordering == -1) {
   // Use initial ordering
      vectorsetlengthatleast(p0, n);
      vectorsetlengthatleast(p1, n);
      for (i = 0; i < n; i++) {
         p0->xZ[i] = i;
         p1->xZ[i] = i;
      }
   }
// Determine algorithm to use:
// * for SKS input or output - use SKS solver (overrides user settings)
// * default is to use Algo=1
   if (algo == 0) {
      algo = 1;
   }
   if (sparseissks(a) || fmt == 2) {
      algo = 2;
   }
   algo = 2;
   if (algo == 2) {
   // Skyline Cholesky with non-skyline output.
   //
   // Call CholeskyX() recursively with Buf.S as output matrix,
   // then perform conversion from SKS to desired format. We can
   // use Buf.S in reccurrent call because SKS-to-SKS CholeskyX()
   // does not uses this field.
      if (fmt != 2) {
         result = sparsecholeskyx(a, n, isupper, p0, p1, -3, algo, 2, buf, &buf->s);
         if (result) {
            sparsecopytobuf(&buf->s, fmt, c);
         }
         ae_frame_leave();
         return result;
      }
   // Skyline Cholesky with skyline output
      if (sparseissks(a) && ordering == -1) {
      // Non-permuted skyline matrix.
      //
      // Quickly copy matrix to output buffer without permutation.
      //
      // NOTE: Buf.D is used as dummy vector filled with zeros.
         vectorsetlengthatleast(&buf->d, n);
         for (i = 0; i < n; i++) {
            buf->d.xZ[i] = 0;
         }
         if (isupper) {
         // Create strictly upper-triangular matrix,
         // copy upper triangle of input.
            sparsecreatesksbuf(n, n, &buf->d, &a->uidx, c);
            for (i = 0; i < n; i++) {
               t0 = a->ridx.xZ[i + 1] - a->uidx.xZ[i] - 1;
               t1 = a->ridx.xZ[i + 1] - 1;
               k = c->ridx.xZ[i + 1] - c->uidx.xZ[i] - 1;
               for (j = t0; j <= t1; j++) {
                  c->vals.xR[k] = a->vals.xR[j];
                  k++;
               }
            }
         } else {
         // Create strictly lower-triangular matrix,
         // copy lower triangle of input.
            sparsecreatesksbuf(n, n, &a->didx, &buf->d, c);
            for (i = 0; i < n; i++) {
               t0 = a->ridx.xZ[i];
               t1 = a->ridx.xZ[i] + a->didx.xZ[i];
               k = c->ridx.xZ[i];
               for (j = t0; j <= t1; j++) {
                  c->vals.xR[k] = a->vals.xR[j];
                  k++;
               }
            }
         }
      } else {
      // Non-identity permutations OR non-skyline input:
      // * investigate profile of permuted A
      // * create skyline matrix in output buffer
      // * copy input with permutation
         vectorsetlengthatleast(&buf->d, n);
         vectorsetlengthatleast(&buf->u, n);
         for (i = 0; i < n; i++) {
            buf->d.xZ[i] = 0;
            buf->u.xZ[i] = 0;
         }
         t0 = 0;
         t1 = 0;
         while (sparseenumerate(a, &t0, &t1, &i, &j, &v)) {
            if (isupper? j >= i: j <= i) {
               i = p1->xZ[i];
               j = p1->xZ[j];
               if (isupper? j < i: j > i) {
                  swapi(&i, &j);
               }
               if (i > j) {
                  buf->d.xZ[i] = imax2(buf->d.xZ[i], i - j);
               } else {
                  buf->u.xZ[j] = imax2(buf->u.xZ[j], j - i);
               }
            }
         }
         sparsecreatesksbuf(n, n, &buf->d, &buf->u, c);
         t0 = 0;
         t1 = 0;
         while (sparseenumerate(a, &t0, &t1, &i, &j, &v)) {
            if (isupper? j >= i: j <= i) {
               i = p1->xZ[i];
               j = p1->xZ[j];
               if (isupper? j < i: j > i) {
                  swapi(&j, &i);
               }
               sparserewriteexisting(c, i, j, v);
            }
         }
      }
      result = sparsecholeskyskyline(c, n, isupper);
      ae_frame_leave();
      return result;
   }
   ae_assert(false, "SparseCholeskyX: internal error - unexpected algorithm");
   ae_frame_leave();
   return result;
}

void rmatrixlup(RMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double mx;
   double v;
   ae_frame_make(&_frame_block);
   SetVector(pivots);
   NewVector(tmp, 0, DT_REAL);
// Internal LU decomposition subroutine.
// Never call it directly.
   ae_assert(m > 0, "RMatrixLUP: incorrect M!");
   ae_assert(n > 0, "RMatrixLUP: incorrect N!");
// Scale matrix to avoid overflows,
// decompose it, then scale back.
   mx = 0.0;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         mx = rmax2(mx, fabs(a->xyR[i][j]));
      }
   }
   if (mx != 0.0) {
      v = 1 / mx;
      for (i = 0; i < m; i++) {
         ae_v_muld(a->xyR[i], 1, n, v);
      }
   }
   ae_vector_set_length(pivots, imin2(m, n));
   ae_vector_set_length(&tmp, 2 * imax2(m, n));
   rmatrixluprec(a, 0, m, n, pivots, &tmp);
   if (mx != 0.0) {
      v = mx;
      for (i = 0; i < m; i++) {
         ae_v_muld(a->xyR[i], 1, imin2(i + 1, n), v);
      }
   }
   ae_frame_leave();
}

void cmatrixlup(CMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double mx;
   double v;
   ae_frame_make(&_frame_block);
   SetVector(pivots);
   NewVector(tmp, 0, DT_COMPLEX);
// Internal LU decomposition subroutine.
// Never call it directly.
   ae_assert(m > 0, "CMatrixLUP: incorrect M!");
   ae_assert(n > 0, "CMatrixLUP: incorrect N!");
// Scale matrix to avoid overflows,
// decompose it, then scale back.
   mx = 0.0;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         mx = rmax2(mx, ae_c_abs(a->xyC[i][j]));
      }
   }
   if (mx != 0.0) {
      v = 1 / mx;
      for (i = 0; i < m; i++) {
         ae_v_cmuld(a->xyC[i], 1, n, v);
      }
   }
   ae_vector_set_length(pivots, imin2(m, n));
   ae_vector_set_length(&tmp, 2 * imax2(m, n));
   cmatrixluprec(a, 0, m, n, pivots, &tmp);
   if (mx != 0.0) {
      v = mx;
      for (i = 0; i < m; i++) {
         ae_v_cmuld(a->xyC[i], 1, imin2(i + 1, n), v);
      }
   }
   ae_frame_leave();
}

void rmatrixplu(RMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double mx;
   double v;
   ae_frame_make(&_frame_block);
   SetVector(pivots);
   NewVector(tmp, 0, DT_REAL);
// Internal LU decomposition subroutine.
// Never call it directly.
   ae_assert(m > 0, "RMatrixPLU: incorrect M!");
   ae_assert(n > 0, "RMatrixPLU: incorrect N!");
   ae_vector_set_length(&tmp, 2 * imax2(m, n));
   ae_vector_set_length(pivots, imin2(m, n));
// Scale matrix to avoid overflows,
// decompose it, then scale back.
   mx = 0.0;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         mx = rmax2(mx, fabs(a->xyR[i][j]));
      }
   }
   if (mx != 0.0) {
      v = 1 / mx;
      for (i = 0; i < m; i++) {
         ae_v_muld(a->xyR[i], 1, n, v);
      }
   }
   rmatrixplurec(a, 0, m, n, pivots, &tmp);
   if (mx != 0.0) {
      v = mx;
      for (i = 0; i < imin2(m, n); i++) {
         ae_v_muld(&a->xyR[i][i], 1, n - i, v);
      }
   }
   ae_frame_leave();
}

void cmatrixplu(CMatrix *a, ae_int_t m, ae_int_t n, ZVector *pivots) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double mx;
   ae_complex v;
   ae_frame_make(&_frame_block);
   SetVector(pivots);
   NewVector(tmp, 0, DT_COMPLEX);
// Internal LU decomposition subroutine.
// Never call it directly.
   ae_assert(m > 0, "CMatrixPLU: incorrect M!");
   ae_assert(n > 0, "CMatrixPLU: incorrect N!");
   ae_vector_set_length(&tmp, 2 * imax2(m, n));
   ae_vector_set_length(pivots, imin2(m, n));
// Scale matrix to avoid overflows,
// decompose it, then scale back.
   mx = 0.0;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         mx = rmax2(mx, ae_c_abs(a->xyC[i][j]));
      }
   }
   if (mx != 0.0) {
      v = ae_complex_from_d(1 / mx);
      for (i = 0; i < m; i++) {
         ae_v_cmulc(a->xyC[i], 1, n, v);
      }
   }
   cmatrixplurec(a, 0, m, n, pivots, &tmp);
   if (mx != 0.0) {
      v = ae_complex_from_d(mx);
      for (i = 0; i < imin2(m, n); i++) {
         ae_v_cmulc(&a->xyC[i][i], 1, n - i, v);
      }
   }
   ae_frame_leave();
}

// Level-2 Cholesky subroutine
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static bool trfac_spdmatrixcholesky2(RMatrix *aaa, ae_int_t offs, ae_int_t n, bool isupper, RVector *tmp) {
   ae_int_t i;
   ae_int_t j;
   double ajj;
   double v;
   double r;
   bool result;
   result = true;
   if (n < 0) {
      result = false;
      return result;
   }
// Quick return if possible
   if (n == 0) {
      return result;
   }
   if (isupper) {
   // Compute the Cholesky factorization A = U'*U.
      for (j = 0; j < n; j++) {
      // Compute U(J,J) and test for non-positive-definiteness.
         v = ae_v_dotproduct(&aaa->xyR[offs][offs + j], aaa->stride, &aaa->xyR[offs][offs + j], aaa->stride, j);
         ajj = aaa->xyR[offs + j][offs + j] - v;
         if (ajj <= 0.0) {
            aaa->xyR[offs + j][offs + j] = ajj;
            result = false;
            return result;
         }
         ajj = sqrt(ajj);
         aaa->xyR[offs + j][offs + j] = ajj;
      // Compute elements J+1:N-1 of row J.
         if (j < n - 1) {
            if (j > 0) {
               ae_v_moveneg(tmp->xR, 1, &aaa->xyR[offs][offs + j], aaa->stride, j);
               rmatrixmv(n - j - 1, j, aaa, offs, offs + j + 1, 1, tmp, 0, tmp, n);
               ae_v_add(&aaa->xyR[offs + j][offs + j + 1], 1, &tmp->xR[n], 1, n - j - 1);
            }
            r = 1 / ajj;
            ae_v_muld(&aaa->xyR[offs + j][offs + j + 1], 1, n - j - 1, r);
         }
      }
   } else {
   // Compute the Cholesky factorization A = L*L'.
      for (j = 0; j < n; j++) {
      // Compute L(J+1,J+1) and test for non-positive-definiteness.
         v = ae_v_dotproduct(&aaa->xyR[offs + j][offs], 1, &aaa->xyR[offs + j][offs], 1, j);
         ajj = aaa->xyR[offs + j][offs + j] - v;
         if (ajj <= 0.0) {
            aaa->xyR[offs + j][offs + j] = ajj;
            result = false;
            return result;
         }
         ajj = sqrt(ajj);
         aaa->xyR[offs + j][offs + j] = ajj;
      // Compute elements J+1:N of column J.
         if (j < n - 1) {
            r = 1 / ajj;
            if (j > 0) {
               ae_v_move(tmp->xR, 1, &aaa->xyR[offs + j][offs], 1, j);
               rmatrixmv(n - j - 1, j, aaa, offs + j + 1, offs, 0, tmp, 0, tmp, n);
               for (i = 0; i < n - j - 1; i++) {
                  aaa->xyR[offs + j + 1 + i][offs + j] = (aaa->xyR[offs + j + 1 + i][offs + j] - tmp->xR[n + i]) * r;
               }
            } else {
               for (i = 0; i < n - j - 1; i++) {
                  aaa->xyR[offs + j + 1 + i][offs + j] *= r;
               }
            }
         }
      }
   }
   return result;
}

// Advanced interface for SPDMatrixCholesky, performs no temporary allocations.
//
// Inputs:
//     A       -   matrix given by upper or lower triangle
//     Offs    -   offset of diagonal block to decompose
//     N       -   diagonal block size
//     IsUpper -   what half is given
//     Tmp     -   temporary array; allocated by function, if its size is too
//                 small; can be reused on subsequent calls.
//
// Outputs:
//     A       -   upper (or lower) triangle contains Cholesky decomposition
//
// Result:
//     True, on success
//     False, on failure
//
// ALGLIB Routine: Copyright 15.12.2009 by Sergey Bochkanov
bool spdmatrixcholeskyrec(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, RVector *tmp) {
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   bool result;
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
// Allocate temporaries
   if (tmp->cnt < 2 * n) {
      ae_vector_set_length(tmp, 2 * n);
   }
// Basecases
   if (n < 1) {
      result = false;
      return result;
   }
   if (n == 1) {
      if (a->xyR[offs][offs] > 0.0) {
         a->xyR[offs][offs] = sqrt(a->xyR[offs][offs]);
         result = true;
      } else {
         result = false;
      }
      return result;
   }
   if (n <= tsb) {
      if (spdmatrixcholeskymkl(a, offs, n, isupper, &result)) {
         return result;
      }
   }
   if (n <= tsa) {
      result = trfac_spdmatrixcholesky2(a, offs, n, isupper, tmp);
      return result;
   }
// Split task into smaller ones
   if (n > tsb) {
   // Split leading B-sized block from the beginning (block-matrix approach)
      n1 = tsb;
      n2 = n - n1;
   } else {
   // Smaller than B-size, perform cache-oblivious split
      n1 = tiledsplit(n, tsa), n2 = n - n1;
   }
   result = spdmatrixcholeskyrec(a, offs, n1, isupper, tmp);
   if (!result) {
      return result;
   }
   if (n2 > 0) {
      if (isupper) {
         rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 1, a, offs, offs + n1);
         rmatrixsyrk(n2, n1, -1.0, a, offs, offs + n1, 1, 1.0, a, offs + n1, offs + n1, isupper);
      } else {
         rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 1, a, offs + n1, offs);
         rmatrixsyrk(n2, n1, -1.0, a, offs + n1, offs, 0, 1.0, a, offs + n1, offs + n1, isupper);
      }
      result = spdmatrixcholeskyrec(a, offs + n1, n2, isupper, tmp);
      if (!result) {
         return result;
      }
   }
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
void rmatrixlu(real_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlu(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, pivots));
   alglib_impl::ae_state_clear();
}

void cmatrixlu(complex_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlu(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, pivots));
   alglib_impl::ae_state_clear();
}

bool spdmatrixcholesky(real_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::spdmatrixcholesky(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return Ok;
}

bool hpdmatrixcholesky(complex_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::hpdmatrixcholesky(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return Ok;
}

void spdmatrixcholeskyupdateadd1(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskyupdateadd1(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, u));
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskyupdatefix(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskyupdatefix(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, fix));
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskyupdateadd1buf(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u, real_1d_array &bufr) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskyupdateadd1buf(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, u), ConstT(ae_vector, bufr));
   alglib_impl::ae_state_clear();
}

void spdmatrixcholeskyupdatefixbuf(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix, real_1d_array &bufr) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskyupdatefixbuf(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, fix), ConstT(ae_vector, bufr));
   alglib_impl::ae_state_clear();
}

bool sparselu(const sparsematrix &a, const ae_int_t pivottype, integer_1d_array &p, integer_1d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparselu(ConstT(sparsematrix, a), pivottype, ConstT(ae_vector, p), ConstT(ae_vector, q));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool sparsecholeskyskyline(const sparsematrix &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::sparsecholeskyskyline(ConstT(sparsematrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return Ok;
}
} // end of namespace alglib

// === RCOND Package ===
// Depends on: (AlgLibInternal) TRLINSOLVE, SAFESOLVE
// Depends on: TRFAC
namespace alglib_impl {
// Internal subroutine for matrix norm estimation
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static bool rcond_rmatrixestimatenorm(ae_int_t n, RVector *v, RVector *x, ZVector *isgn, double *est, ae_int_t *kase) {
   ae_int_t i;
   double t;
   bool flg;
   AutoS ae_int_t iter;
   AutoS ae_int_t j;
   AutoS ae_int_t jlast;
   AutoS ae_int_t jump;
   AutoS double altsgn;
   AutoS double estold;
   AutoS double temp;
   const ae_int_t itmax = 5;
// Manually threaded two-way signalling.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (*kase != 0) switch (jump) {
      case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      case 4: goto Resume4; case 5: goto Resume5;
      default: goto Exit;
   }
Spawn:
   ae_vector_set_length(v, n + 1);
   ae_vector_set_length(x, n + 1);
   ae_vector_set_length(isgn, n + 1);
   t = 1.0 / (double)n;
   for (i = 1; i <= n; i++) {
      x->xR[i] = t;
   }
   *kase = 1; jump = 1; goto Pause; Resume1: // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
   if (n == 1) {
      v->xR[1] = x->xR[1];
      *est = fabs(v->xR[1]);
      goto Exit;
   }
   *est = 0.0;
   for (i = 1; i <= n; i++) {
      *est += fabs(x->xR[i]);
   }
   for (i = 1; i <= n; i++) {
      if (x->xR[i] >= 0.0) {
         x->xR[i] = 1.0;
      } else {
         x->xR[i] = -1.0;
      }
      isgn->xZ[i] = ae_sign(x->xR[i]);
   }
   *kase = 2; jump = 2; goto Pause; Resume2: // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
   j = 1;
   for (i = 2; i <= n; i++) {
      if (fabs(x->xR[i]) > fabs(x->xR[j])) {
         j = i;
      }
   }
   iter = 2;
// MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
   do {
      for (i = 1; i <= n; i++) {
         x->xR[i] = 0.0;
      }
      x->xR[j] = 1.0;
      *kase = 1; jump = 3; goto Pause; Resume3: // X HAS BEEN OVERWRITTEN BY A*X.
      ae_v_move(&v->xR[1], 1, &x->xR[1], 1, n);
      estold = *est;
      *est = 0.0;
      for (i = 1; i <= n; i++) {
         *est += fabs(v->xR[i]);
      }
      flg = false;
      for (i = 1; i <= n; i++) {
         if (x->xR[i] >= 0.0 && isgn->xZ[i] < 0 || x->xR[i] < 0.0 && isgn->xZ[i] >= 0) {
            flg = true;
         }
      }
   // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
   // OR MAY BE CYCLING.
      if (!flg || *est <= estold) break;
      for (i = 1; i <= n; i++) {
         if (x->xR[i] >= 0.0) {
            x->xR[i] = 1.0;
            isgn->xZ[i] = 1;
         } else {
            x->xR[i] = -1.0;
            isgn->xZ[i] = -1;
         }
      }
      *kase = 2; jump = 4; goto Pause; Resume4: // X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
      jlast = j;
      j = 1;
      for (i = 2; i <= n; i++) {
         if (fabs(x->xR[i]) > fabs(x->xR[j])) {
            j = i;
         }
      }
   } while (x->xR[jlast] != fabs(x->xR[j]) && iter++ < itmax);
// ITERATION COMPLETE.  FINAL STAGE.
   altsgn = 1.0;
   for (i = 1; i <= n; i++) {
      x->xR[i] = altsgn * (1 + (double)(i - 1) / (double)(n - 1));
      altsgn = -altsgn;
   }
   *kase = 1; jump = 5; goto Pause; Resume5: // X HAS BEEN OVERWRITTEN BY A*X.
   temp = 0.0;
   for (i = 1; i <= n; i++) {
      temp += fabs(x->xR[i]);
   }
   temp = 2 * temp / (3 * n);
   if (temp > *est) {
      ae_v_move(&v->xR[1], 1, &x->xR[1], 1, n);
      *est = temp;
   }
Exit:
   *kase = 0;
   return false;
Pause:
   return true;
}

static double rcond_internalcomplexrcondscsum1(CVector *x, ae_int_t n) {
   ae_int_t i;
   double result;
   result = 0.0;
   for (i = 1; i <= n; i++) {
      result += ae_c_abs(x->xC[i]);
   }
   return result;
}

static ae_int_t rcond_internalcomplexrcondicmax1(CVector *x, ae_int_t n) {
   ae_int_t i;
   double m;
   ae_int_t result;
   result = 1;
   m = ae_c_abs(x->xC[1]);
   for (i = 2; i <= n; i++) {
      if (ae_c_abs(x->xC[i]) > m) {
         result = i;
         m = ae_c_abs(x->xC[i]);
      }
   }
   return result;
}

static bool rcond_cmatrixestimatenorm(ae_int_t n, CVector *v, CVector *x, double *est, ae_int_t *kase) {
   AutoS ae_int_t i;
   AutoS ae_int_t iter;
   AutoS ae_int_t j;
   AutoS ae_int_t jlast;
   AutoS ae_int_t jump;
   AutoS double absxi;
   AutoS double altsgn;
   AutoS double estold;
   AutoS double temp;
// Executable Statements ..
   const ae_int_t itmax = 5;
   const double safmin = ae_minrealnumber;
// Manually threaded two-way signalling.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (*kase != 0) switch (jump) {
      case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      case 4: goto Resume4; case 5: goto Resume5;
      default: goto Exit;
   }
Spawn:
   ae_vector_set_length(v, n + 1);
   ae_vector_set_length(x, n + 1);
   for (i = 1; i <= n; i++) {
      x->xC[i] = ae_complex_from_d(1.0 / (double)n);
   }
   *kase = 1; jump = 1; goto Pause; Resume1: // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
   if (n == 1) {
      v->xC[1] = x->xC[1];
      *est = ae_c_abs(v->xC[1]);
      goto Exit;
   }
   *est = rcond_internalcomplexrcondscsum1(x, n);
   for (i = 1; i <= n; i++) {
      absxi = ae_c_abs(x->xC[i]);
      if (absxi > safmin) {
         x->xC[i] = ae_c_div_d(x->xC[i], absxi);
      } else {
         x->xC[i] = ae_complex_from_i(1);
      }
   }
   *kase = 2; jump = 2; goto Pause; Resume2: // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
   j = rcond_internalcomplexrcondicmax1(x, n);
   iter = 2;
   do {
   // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
      for (i = 1; i <= n; i++) {
         x->xC[i] = ae_complex_from_i(0);
      }
      x->xC[j] = ae_complex_from_i(1);
      *kase = 1; jump = 3; goto Pause; Resume3: // X HAS BEEN OVERWRITTEN BY A*X.
      ae_v_cmove(&v->xC[1], 1, &x->xC[1], 1, "N", n);
      estold = *est;
      *est = rcond_internalcomplexrcondscsum1(v, n);
   // TEST FOR CYCLING.
      if (*est <= estold) {
         break;
      }
      for (i = 1; i <= n; i++) {
         absxi = ae_c_abs(x->xC[i]);
         if (absxi > safmin) {
            x->xC[i] = ae_c_div_d(x->xC[i], absxi);
         } else {
            x->xC[i] = ae_complex_from_i(1);
         }
      }
      *kase = 2; jump = 4; goto Pause; Resume4: // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
      jlast = j;
      j = rcond_internalcomplexrcondicmax1(x, n);
   } while (ae_c_abs(x->xC[jlast]) != ae_c_abs(x->xC[j]) && iter++ < itmax);
// ITERATION COMPLETE.  FINAL STAGE.
   altsgn = 1.0;
   for (i = 1; i <= n; i++) {
      x->xC[i] = ae_complex_from_d(altsgn * (1 + (double)(i - 1) / (double)(n - 1)));
      altsgn = -altsgn;
   }
   *kase = 1; jump = 5; goto Pause; Resume5: // X HAS BEEN OVERWRITTEN BY A*X.
   temp = 2 * (rcond_internalcomplexrcondscsum1(x, n) / (3 * n));
   if (temp > *est) {
      ae_v_cmove(&v->xC[1], 1, &x->xC[1], 1, "N", n);
      *est = temp;
   }
Exit:
   *kase = 0;
   return false;
Pause:
   return true;
}

// Threshold for rcond: matrices with condition number beyond this  threshold
// are considered singular.
//
// Threshold must be far enough from underflow, at least Sqr(Threshold)  must
// be greater than underflow.
double rcondthreshold() {
   double result;
   result = sqrt(sqrt(ae_minrealnumber));
   return result;
}

// Internal subroutine for condition number estimation
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static void rcond_rmatrixrcondluinternal(RMatrix *lua, ae_int_t n, bool onenorm, bool isanormprovided, double anorm, double *rc) {
   ae_frame _frame_block;
   double v;
   ae_int_t i;
   ae_int_t j;
   ae_int_t kase;
   ae_int_t kase1;
   double ainvnm;
   double maxgrowth;
   double su;
   double sl;
   bool mupper;
   bool munit;
   ae_frame_make(&_frame_block);
   *rc = 0;
   NewVector(ex, 0, DT_REAL);
   NewVector(ev, 0, DT_REAL);
   NewVector(iwork, 0, DT_INT);
   NewVector(tmp, 0, DT_REAL);
// RC=0 if something happens
   *rc = 0.0;
// init
   if (onenorm) {
      kase1 = 1;
   } else {
      kase1 = 2;
   }
   mupper = true;
   munit = true;
   ae_vector_set_length(&iwork, n + 1);
   ae_vector_set_length(&tmp, n);
// prepare parameters for triangular solver
   maxgrowth = 1 / rcondthreshold();
   su = 0.0;
   sl = 1.0;
   for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
         sl = rmax2(sl, fabs(lua->xyR[i][j]));
      }
      for (j = i; j < n; j++) {
         su = rmax2(su, fabs(lua->xyR[i][j]));
      }
   }
   if (su == 0.0) {
      su = 1.0;
   }
   su = 1 / su;
   sl = 1 / sl;
// Estimate the norm of A.
   if (!isanormprovided) {
      kase = 0;
      anorm = 0.0;
      while (rcond_rmatrixestimatenorm(n, &ev, &ex, &iwork, &anorm, &kase)) {
         if (kase == kase1) {
         // Multiply by U
            for (i = 1; i <= n; i++) {
               v = ae_v_dotproduct(&lua->xyR[i - 1][i - 1], 1, &ex.xR[i], 1, n - i + 1);
               ex.xR[i] = v;
            }
         // Multiply by L
            for (i = n; i >= 1; i--) {
               if (i > 1) {
                  v = ae_v_dotproduct(lua->xyR[i - 1], 1, &ex.xR[1], 1, i - 1);
               } else {
                  v = 0.0;
               }
               ex.xR[i] += v;
            }
         } else {
         // Multiply by L'
            for (i = 0; i < n; i++) {
               tmp.xR[i] = 0.0;
            }
            for (i = 0; i < n; i++) {
               v = ex.xR[i + 1];
               if (i >= 1) {
                  ae_v_addd(tmp.xR, 1, lua->xyR[i], 1, i, v);
               }
               tmp.xR[i] += v;
            }
            ae_v_move(&ex.xR[1], 1, tmp.xR, 1, n);
         // Multiply by U'
            for (i = 0; i < n; i++) {
               tmp.xR[i] = 0.0;
            }
            for (i = 0; i < n; i++) {
               v = ex.xR[i + 1];
               ae_v_addd(&tmp.xR[i], 1, &lua->xyR[i][i], 1, n - i, v);
            }
            ae_v_move(&ex.xR[1], 1, tmp.xR, 1, n);
         }
      }
   }
// Scale according to SU/SL
   anorm *= su * sl;
// Quick return if possible
// We assume that ANORM != 0 after this block
   if (anorm == 0.0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      *rc = 1.0;
      ae_frame_leave();
      return;
   }
// Estimate the norm of inv(A).
   ainvnm = 0.0;
   kase = 0;
   while (rcond_rmatrixestimatenorm(n, &ev, &ex, &iwork, &ainvnm, &kase)) {
   // from 1-based array to 0-based
      for (i = 0; i < n; i++) {
         ex.xR[i] = ex.xR[i + 1];
      }
   // multiply by inv(A) or inv(A')
      if (kase == kase1) {
      // Multiply by inv(L).
         if (!rmatrixscaledtrsafesolve(lua, sl, n, &ex, !mupper, 0, munit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      // Multiply by inv(U).
         if (!rmatrixscaledtrsafesolve(lua, su, n, &ex, mupper, 0, !munit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      } else {
      // Multiply by inv(U').
         if (!rmatrixscaledtrsafesolve(lua, su, n, &ex, mupper, 1, !munit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      // Multiply by inv(L').
         if (!rmatrixscaledtrsafesolve(lua, sl, n, &ex, !mupper, 1, munit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      }
   // from 0-based array to 1-based
      for (i = n - 1; i >= 0; i--) {
         ex.xR[i + 1] = ex.xR[i];
      }
   }
// Compute the estimate of the reciprocal condition number.
   if (ainvnm != 0.0) {
      *rc = 1 / ainvnm;
      *rc /= anorm;
      if (*rc < rcondthreshold()) {
         *rc = 0.0;
      }
   }
   ae_frame_leave();
}

// Condition number estimation
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      March 31, 1993
static void rcond_cmatrixrcondluinternal(CMatrix *lua, ae_int_t n, bool onenorm, bool isanormprovided, double anorm, double *rc) {
   ae_frame _frame_block;
   ae_int_t kase;
   ae_int_t kase1;
   double ainvnm;
   ae_complex v;
   ae_int_t i;
   ae_int_t j;
   double su;
   double sl;
   double maxgrowth;
   ae_frame_make(&_frame_block);
   *rc = 0;
   NewVector(ex, 0, DT_COMPLEX);
   NewVector(cwork2, 0, DT_COMPLEX);
   NewVector(cwork3, 0, DT_COMPLEX);
   NewVector(cwork4, 0, DT_COMPLEX);
   if (n <= 0) {
      ae_frame_leave();
      return;
   }
   ae_vector_set_length(&cwork2, n + 1);
   *rc = 0.0;
   if (n == 0) {
      *rc = 1.0;
      ae_frame_leave();
      return;
   }
// prepare parameters for triangular solver
   maxgrowth = 1 / rcondthreshold();
   su = 0.0;
   sl = 1.0;
   for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
         sl = rmax2(sl, ae_c_abs(lua->xyC[i][j]));
      }
      for (j = i; j < n; j++) {
         su = rmax2(su, ae_c_abs(lua->xyC[i][j]));
      }
   }
   if (su == 0.0) {
      su = 1.0;
   }
   su = 1 / su;
   sl = 1 / sl;
// Estimate the norm of SU*SL*A.
   if (!isanormprovided) {
      anorm = 0.0;
      if (onenorm) {
         kase1 = 1;
      } else {
         kase1 = 2;
      }
      kase = 0;
      do {
         if (rcond_cmatrixestimatenorm(n, &cwork4, &ex, &anorm, &kase)) {
            if (kase == kase1) {
            // Multiply by U
               for (i = 1; i <= n; i++) {
                  v = ae_v_cdotproduct(&lua->xyC[i - 1][i - 1], 1, "N", &ex.xC[i], 1, "N", n - i + 1);
                  ex.xC[i] = v;
               }
            // Multiply by L
               for (i = n; i >= 1; i--) {
                  v = ae_complex_from_i(0);
                  if (i > 1) {
                     v = ae_v_cdotproduct(lua->xyC[i - 1], 1, "N", &ex.xC[1], 1, "N", i - 1);
                  }
                  ex.xC[i] = ae_c_add(v, ex.xC[i]);
               }
            } else {
            // Multiply by L'
               for (i = 1; i <= n; i++) {
                  cwork2.xC[i] = ae_complex_from_i(0);
               }
               for (i = 1; i <= n; i++) {
                  v = ex.xC[i];
                  if (i > 1) {
                     ae_v_caddc(&cwork2.xC[1], 1, lua->xyC[i - 1], 1, "Conj", i - 1, v);
                  }
                  cwork2.xC[i] = ae_c_add(cwork2.xC[i], v);
               }
            // Multiply by U'
               for (i = 1; i <= n; i++) {
                  ex.xC[i] = ae_complex_from_i(0);
               }
               for (i = 1; i <= n; i++) {
                  v = cwork2.xC[i];
                  ae_v_caddc(&ex.xC[i], 1, &lua->xyC[i - 1][i - 1], 1, "Conj", n - i + 1, v);
               }
            }
         }
      } while (kase != 0);
   }
// Scale according to SU/SL
   anorm *= su * sl;
// Quick return if possible
   if (anorm == 0.0) {
      ae_frame_leave();
      return;
   }
// Estimate the norm of inv(A).
   ainvnm = 0.0;
   if (onenorm) {
      kase1 = 1;
   } else {
      kase1 = 2;
   }
   kase = 0;
   while (rcond_cmatrixestimatenorm(n, &cwork4, &ex, &ainvnm, &kase)) {
   // From 1-based to 0-based
      for (i = 0; i < n; i++) {
         ex.xC[i] = ex.xC[i + 1];
      }
   // multiply by inv(A) or inv(A')
      if (kase == kase1) {
      // Multiply by inv(L).
         if (!cmatrixscaledtrsafesolve(lua, sl, n, &ex, false, 0, true, maxgrowth)) {
            *rc = 0.0;
            ae_frame_leave();
            return;
         }
      // Multiply by inv(U).
         if (!cmatrixscaledtrsafesolve(lua, su, n, &ex, true, 0, false, maxgrowth)) {
            *rc = 0.0;
            ae_frame_leave();
            return;
         }
      } else {
      // Multiply by inv(U').
         if (!cmatrixscaledtrsafesolve(lua, su, n, &ex, true, 2, false, maxgrowth)) {
            *rc = 0.0;
            ae_frame_leave();
            return;
         }
      // Multiply by inv(L').
         if (!cmatrixscaledtrsafesolve(lua, sl, n, &ex, false, 2, true, maxgrowth)) {
            *rc = 0.0;
            ae_frame_leave();
            return;
         }
      }
   // from 0-based to 1-based
      for (i = n - 1; i >= 0; i--) {
         ex.xC[i + 1] = ex.xC[i];
      }
   }
// Compute the estimate of the reciprocal condition number.
   if (ainvnm != 0.0) {
      *rc = 1 / ainvnm;
      *rc /= anorm;
      if (*rc < rcondthreshold()) {
         *rc = 0.0;
      }
   }
   ae_frame_leave();
}

// Estimate of a matrix condition number (1-norm)
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double rmatrixrcond1(const real_2d_array &a, const ae_int_t n);
double rmatrixrcond1(RMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(pivots, 0, DT_INT);
   NewVector(t, 0, DT_REAL);
   ae_assert(n >= 1, "RMatrixRCond1: N<1!");
   ae_vector_set_length(&t, n);
   for (i = 0; i < n; i++) {
      t.xR[i] = 0.0;
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         t.xR[j] += fabs(a->xyR[i][j]);
      }
   }
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      nrm = rmax2(nrm, t.xR[i]);
   }
   rmatrixlu(a, n, n, &pivots);
   rcond_rmatrixrcondluinternal(a, n, true, true, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Estimate of a matrix condition number (1-norm)
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double cmatrixrcond1(const complex_2d_array &a, const ae_int_t n);
double cmatrixrcond1(CMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(pivots, 0, DT_INT);
   NewVector(t, 0, DT_REAL);
   ae_assert(n >= 1, "CMatrixRCond1: N<1!");
   ae_vector_set_length(&t, n);
   for (i = 0; i < n; i++) {
      t.xR[i] = 0.0;
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         t.xR[j] += ae_c_abs(a->xyC[i][j]);
      }
   }
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      nrm = rmax2(nrm, t.xR[i]);
   }
   cmatrixlu(a, n, n, &pivots);
   rcond_cmatrixrcondluinternal(a, n, true, true, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Estimate of a matrix condition number (infinity-norm).
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double rmatrixrcondinf(const real_2d_array &a, const ae_int_t n);
double rmatrixrcondinf(RMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n >= 1, "RMatrixRCondInf: N<1!");
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      v = 0.0;
      for (j = 0; j < n; j++) {
         v += fabs(a->xyR[i][j]);
      }
      nrm = rmax2(nrm, v);
   }
   rmatrixlu(a, n, n, &pivots);
   rcond_rmatrixrcondluinternal(a, n, false, true, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Estimate of a matrix condition number (infinity-norm).
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double cmatrixrcondinf(const complex_2d_array &a, const ae_int_t n);
double cmatrixrcondinf(CMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n >= 1, "CMatrixRCondInf: N<1!");
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      v = 0.0;
      for (j = 0; j < n; j++) {
         v += ae_c_abs(a->xyC[i][j]);
      }
      nrm = rmax2(nrm, v);
   }
   cmatrixlu(a, n, n, &pivots);
   rcond_cmatrixrcondluinternal(a, n, false, true, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Internal subroutine for condition number estimation
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static void rcond_spdmatrixrcondcholeskyinternal(RMatrix *cha, ae_int_t n, bool isupper, bool isnormprovided, double anorm, double *rc) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t kase;
   double ainvnm;
   double sa;
   double v;
   double maxgrowth;
   ae_frame_make(&_frame_block);
   *rc = 0;
   NewVector(ex, 0, DT_REAL);
   NewVector(ev, 0, DT_REAL);
   NewVector(tmp, 0, DT_REAL);
   NewVector(iwork, 0, DT_INT);
   ae_assert(n >= 1, "Assertion failed");
   ae_vector_set_length(&tmp, n);
// RC=0 if something happens
   *rc = 0.0;
// prepare parameters for triangular solver
   maxgrowth = 1 / rcondthreshold();
   sa = 0.0;
   if (isupper) {
      for (i = 0; i < n; i++) {
         for (j = i; j < n; j++) {
            sa = rmax2(sa, ae_c_abs(ae_complex_from_d(cha->xyR[i][j])));
         }
      }
   } else {
      for (i = 0; i < n; i++) {
         for (j = 0; j <= i; j++) {
            sa = rmax2(sa, ae_c_abs(ae_complex_from_d(cha->xyR[i][j])));
         }
      }
   }
   if (sa == 0.0) {
      sa = 1.0;
   }
   sa = 1 / sa;
// Estimate the norm of A.
   if (!isnormprovided) {
      kase = 0;
      anorm = 0.0;
      while (rcond_rmatrixestimatenorm(n, &ev, &ex, &iwork, &anorm, &kase)) {
         if (isupper) {
         // Multiply by U
            for (i = 1; i <= n; i++) {
               v = ae_v_dotproduct(&cha->xyR[i - 1][i - 1], 1, &ex.xR[i], 1, n - i + 1);
               ex.xR[i] = v;
            }
            ae_v_muld(&ex.xR[1], 1, n, sa);
         // Multiply by U'
            for (i = 0; i < n; i++) {
               tmp.xR[i] = 0.0;
            }
            for (i = 0; i < n; i++) {
               v = ex.xR[i + 1];
               ae_v_addd(&tmp.xR[i], 1, &cha->xyR[i][i], 1, n - i, v);
            }
            ae_v_move(&ex.xR[1], 1, tmp.xR, 1, n);
            ae_v_muld(&ex.xR[1], 1, n, sa);
         } else {
         // Multiply by L'
            for (i = 0; i < n; i++) {
               tmp.xR[i] = 0.0;
            }
            for (i = 0; i < n; i++) {
               v = ex.xR[i + 1];
               ae_v_addd(tmp.xR, 1, cha->xyR[i], 1, i + 1, v);
            }
            ae_v_move(&ex.xR[1], 1, tmp.xR, 1, n);
            ae_v_muld(&ex.xR[1], 1, n, sa);
         // Multiply by L
            for (i = n; i >= 1; i--) {
               v = ae_v_dotproduct(cha->xyR[i - 1], 1, &ex.xR[1], 1, i);
               ex.xR[i] = v;
            }
            ae_v_muld(&ex.xR[1], 1, n, sa);
         }
      }
   }
// Quick return if possible
   if (anorm == 0.0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      *rc = 1.0;
      ae_frame_leave();
      return;
   }
// Estimate the 1-norm of inv(A).
   kase = 0;
   while (rcond_rmatrixestimatenorm(n, &ev, &ex, &iwork, &ainvnm, &kase)) {
      for (i = 0; i < n; i++) {
         ex.xR[i] = ex.xR[i + 1];
      }
      if (isupper) {
      // Multiply by inv(U').
         if (!rmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 1, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      // Multiply by inv(U).
         if (!rmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 0, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      } else {
      // Multiply by inv(L).
         if (!rmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 0, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      // Multiply by inv(L').
         if (!rmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 1, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      }
      for (i = n - 1; i >= 0; i--) {
         ex.xR[i + 1] = ex.xR[i];
      }
   }
// Compute the estimate of the reciprocal condition number.
   if (ainvnm != 0.0) {
      v = 1 / ainvnm;
      *rc = v / anorm;
      if (*rc < rcondthreshold()) {
         *rc = 0.0;
      }
   }
   ae_frame_leave();
}

// Internal subroutine for condition number estimation
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static void rcond_hpdmatrixrcondcholeskyinternal(CMatrix *cha, ae_int_t n, bool isupper, bool isnormprovided, double anorm, double *rc) {
   ae_frame _frame_block;
   ae_int_t kase;
   double ainvnm;
   ae_complex v;
   ae_int_t i;
   ae_int_t j;
   double sa;
   double maxgrowth;
   ae_frame_make(&_frame_block);
   *rc = 0;
   NewVector(ex, 0, DT_COMPLEX);
   NewVector(ev, 0, DT_COMPLEX);
   NewVector(tmp, 0, DT_COMPLEX);
   ae_assert(n >= 1, "Assertion failed");
   ae_vector_set_length(&tmp, n);
// RC=0 if something happens
   *rc = 0.0;
// prepare parameters for triangular solver
   maxgrowth = 1 / rcondthreshold();
   sa = 0.0;
   if (isupper) {
      for (i = 0; i < n; i++) {
         for (j = i; j < n; j++) {
            sa = rmax2(sa, ae_c_abs(cha->xyC[i][j]));
         }
      }
   } else {
      for (i = 0; i < n; i++) {
         for (j = 0; j <= i; j++) {
            sa = rmax2(sa, ae_c_abs(cha->xyC[i][j]));
         }
      }
   }
   if (sa == 0.0) {
      sa = 1.0;
   }
   sa = 1 / sa;
// Estimate the norm of A
   if (!isnormprovided) {
      anorm = 0.0;
      kase = 0;
      while (rcond_cmatrixestimatenorm(n, &ev, &ex, &anorm, &kase)) {
         if (isupper) {
         // Multiply by U
            for (i = 1; i <= n; i++) {
               v = ae_v_cdotproduct(&cha->xyC[i - 1][i - 1], 1, "N", &ex.xC[i], 1, "N", n - i + 1);
               ex.xC[i] = v;
            }
            ae_v_cmuld(&ex.xC[1], 1, n, sa);
         // Multiply by U'
            for (i = 0; i < n; i++) {
               tmp.xC[i] = ae_complex_from_i(0);
            }
            for (i = 0; i < n; i++) {
               v = ex.xC[i + 1];
               ae_v_caddc(&tmp.xC[i], 1, &cha->xyC[i][i], 1, "Conj", n - i, v);
            }
            ae_v_cmove(&ex.xC[1], 1, tmp.xC, 1, "N", n);
            ae_v_cmuld(&ex.xC[1], 1, n, sa);
         } else {
         // Multiply by L'
            for (i = 0; i < n; i++) {
               tmp.xC[i] = ae_complex_from_i(0);
            }
            for (i = 0; i < n; i++) {
               v = ex.xC[i + 1];
               ae_v_caddc(tmp.xC, 1, cha->xyC[i], 1, "Conj", i + 1, v);
            }
            ae_v_cmove(&ex.xC[1], 1, tmp.xC, 1, "N", n);
            ae_v_cmuld(&ex.xC[1], 1, n, sa);
         // Multiply by L
            for (i = n; i >= 1; i--) {
               v = ae_v_cdotproduct(cha->xyC[i - 1], 1, "N", &ex.xC[1], 1, "N", i);
               ex.xC[i] = v;
            }
            ae_v_cmuld(&ex.xC[1], 1, n, sa);
         }
      }
   }
// Quick return if possible
// After this block we assume that ANORM != 0
   if (anorm == 0.0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      *rc = 1.0;
      ae_frame_leave();
      return;
   }
// Estimate the norm of inv(A).
   ainvnm = 0.0;
   kase = 0;
   while (rcond_cmatrixestimatenorm(n, &ev, &ex, &ainvnm, &kase)) {
      for (i = 0; i < n; i++) {
         ex.xC[i] = ex.xC[i + 1];
      }
      if (isupper) {
      // Multiply by inv(U').
         if (!cmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 2, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      // Multiply by inv(U).
         if (!cmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 0, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      } else {
      // Multiply by inv(L).
         if (!cmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 0, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      // Multiply by inv(L').
         if (!cmatrixscaledtrsafesolve(cha, sa, n, &ex, isupper, 2, false, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      }
      for (i = n - 1; i >= 0; i--) {
         ex.xC[i + 1] = ex.xC[i];
      }
   }
// Compute the estimate of the reciprocal condition number.
   if (ainvnm != 0.0) {
      *rc = 1 / ainvnm;
      *rc /= anorm;
      if (*rc < rcondthreshold()) {
         *rc = 0.0;
      }
   }
   ae_frame_leave();
}

// Condition number estimate of a symmetric positive definite matrix.
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// It should be noted that 1-norm and inf-norm of condition numbers of symmetric
// matrices are equal, so the algorithm doesn't take into account the
// differences between these types of norms.
//
// Inputs:
//     A       -   symmetric positive definite matrix which is given by its
//                 upper or lower triangle depending on the value of
//                 IsUpper. Array with elements [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   storage format.
//
// Result:
//     1/LowerBound(cond(A)), if matrix A is positive definite,
//    -1, if matrix A is not positive definite, and its condition number
//     could not be found by this algorithm.
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double spdmatrixrcond(const real_2d_array &a, const ae_int_t n, const bool isupper);
double spdmatrixrcond(RMatrix *a, ae_int_t n, bool isupper) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   double v;
   double nrm;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(t, 0, DT_REAL);
   ae_vector_set_length(&t, n);
   for (i = 0; i < n; i++) {
      t.xR[i] = 0.0;
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
         if (i == j) {
            t.xR[i] += fabs(a->xyR[i][i]);
         } else {
            t.xR[i] += fabs(a->xyR[i][j]);
            t.xR[j] += fabs(a->xyR[i][j]);
         }
      }
   }
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      nrm = rmax2(nrm, t.xR[i]);
   }
   if (spdmatrixcholesky(a, n, isupper)) {
      rcond_spdmatrixrcondcholeskyinternal(a, n, isupper, true, nrm, &v);
      result = v;
   } else {
      result = -1.0;
   }
   ae_frame_leave();
   return result;
}

// Condition number estimate of a Hermitian positive definite matrix.
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// It should be noted that 1-norm and inf-norm of condition numbers of symmetric
// matrices are equal, so the algorithm doesn't take into account the
// differences between these types of norms.
//
// Inputs:
//     A       -   Hermitian positive definite matrix which is given by its
//                 upper or lower triangle depending on the value of
//                 IsUpper. Array with elements [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   storage format.
//
// Result:
//     1/LowerBound(cond(A)), if matrix A is positive definite,
//    -1, if matrix A is not positive definite, and its condition number
//     could not be found by this algorithm.
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double hpdmatrixrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper);
double hpdmatrixrcond(CMatrix *a, ae_int_t n, bool isupper) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   double v;
   double nrm;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(t, 0, DT_REAL);
   ae_vector_set_length(&t, n);
   for (i = 0; i < n; i++) {
      t.xR[i] = 0.0;
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
         if (i == j) {
            t.xR[i] += ae_c_abs(a->xyC[i][i]);
         } else {
            t.xR[i] += ae_c_abs(a->xyC[i][j]);
            t.xR[j] += ae_c_abs(a->xyC[i][j]);
         }
      }
   }
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      nrm = rmax2(nrm, t.xR[i]);
   }
   if (hpdmatrixcholesky(a, n, isupper)) {
      rcond_hpdmatrixrcondcholeskyinternal(a, n, isupper, true, nrm, &v);
      result = v;
   } else {
      result = -1.0;
   }
   ae_frame_leave();
   return result;
}

// Internal subroutine for condition number estimation
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992
static void rcond_rmatrixrcondtrinternal(RMatrix *a, ae_int_t n, bool isupper, bool isunit, bool onenorm, double anorm, double *rc) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t kase;
   ae_int_t kase1;
   ae_int_t j1;
   ae_int_t j2;
   double ainvnm;
   double maxgrowth;
   double s;
   ae_frame_make(&_frame_block);
   *rc = 0;
   NewVector(ex, 0, DT_REAL);
   NewVector(ev, 0, DT_REAL);
   NewVector(iwork, 0, DT_INT);
   NewVector(tmp, 0, DT_REAL);
// RC=0 if something happens
   *rc = 0.0;
// init
   if (onenorm) {
      kase1 = 1;
   } else {
      kase1 = 2;
   }
   ae_vector_set_length(&iwork, n + 1);
   ae_vector_set_length(&tmp, n);
// prepare parameters for triangular solver
   maxgrowth = 1 / rcondthreshold();
   s = 0.0;
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i + 1;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i - 1;
      }
      for (j = j1; j <= j2; j++) {
         s = rmax2(s, fabs(a->xyR[i][j]));
      }
      if (isunit) {
         s = rmax2(s, 1.0);
      } else {
         s = rmax2(s, fabs(a->xyR[i][i]));
      }
   }
   if (s == 0.0) {
      s = 1.0;
   }
   s = 1 / s;
// Scale according to S
   anorm *= s;
// Quick return if possible
// We assume that ANORM != 0 after this block
   if (anorm == 0.0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      *rc = 1.0;
      ae_frame_leave();
      return;
   }
// Estimate the norm of inv(A).
   ainvnm = 0.0;
   kase = 0;
   while (rcond_rmatrixestimatenorm(n, &ev, &ex, &iwork, &ainvnm, &kase)) {
   // from 1-based array to 0-based
      for (i = 0; i < n; i++) {
         ex.xR[i] = ex.xR[i + 1];
      }
   // multiply by inv(A) or inv(A')
      if (kase == kase1) {
      // multiply by inv(A)
         if (!rmatrixscaledtrsafesolve(a, s, n, &ex, isupper, 0, isunit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      } else {
      // multiply by inv(A')
         if (!rmatrixscaledtrsafesolve(a, s, n, &ex, isupper, 1, isunit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      }
   // from 0-based array to 1-based
      for (i = n - 1; i >= 0; i--) {
         ex.xR[i + 1] = ex.xR[i];
      }
   }
// Compute the estimate of the reciprocal condition number.
   if (ainvnm != 0.0) {
      *rc = 1 / ainvnm;
      *rc /= anorm;
      if (*rc < rcondthreshold()) {
         *rc = 0.0;
      }
   }
   ae_frame_leave();
}

// Condition number estimation
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      March 31, 1993
static void rcond_cmatrixrcondtrinternal(CMatrix *a, ae_int_t n, bool isupper, bool isunit, bool onenorm, double anorm, double *rc) {
   ae_frame _frame_block;
   ae_int_t kase;
   ae_int_t kase1;
   double ainvnm;
   ae_int_t i;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   double s;
   double maxgrowth;
   ae_frame_make(&_frame_block);
   *rc = 0;
   NewVector(ex, 0, DT_COMPLEX);
   NewVector(cwork2, 0, DT_COMPLEX);
   NewVector(cwork3, 0, DT_COMPLEX);
   NewVector(cwork4, 0, DT_COMPLEX);
// RC=0 if something happens
   *rc = 0.0;
// init
   if (n <= 0) {
      ae_frame_leave();
      return;
   }
   if (n == 0) {
      *rc = 1.0;
      ae_frame_leave();
      return;
   }
   ae_vector_set_length(&cwork2, n + 1);
// prepare parameters for triangular solver
   maxgrowth = 1 / rcondthreshold();
   s = 0.0;
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i + 1;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i - 1;
      }
      for (j = j1; j <= j2; j++) {
         s = rmax2(s, ae_c_abs(a->xyC[i][j]));
      }
      if (isunit) {
         s = rmax2(s, 1.0);
      } else {
         s = rmax2(s, ae_c_abs(a->xyC[i][i]));
      }
   }
   if (s == 0.0) {
      s = 1.0;
   }
   s = 1 / s;
// Scale according to S
   anorm *= s;
// Quick return if possible
   if (anorm == 0.0) {
      ae_frame_leave();
      return;
   }
// Estimate the norm of inv(A).
   ainvnm = 0.0;
   if (onenorm) {
      kase1 = 1;
   } else {
      kase1 = 2;
   }
   kase = 0;
   while (rcond_cmatrixestimatenorm(n, &cwork4, &ex, &ainvnm, &kase)) {
   // From 1-based to 0-based
      for (i = 0; i < n; i++) {
         ex.xC[i] = ex.xC[i + 1];
      }
   // multiply by inv(A) or inv(A')
      if (kase == kase1) {
      // multiply by inv(A)
         if (!cmatrixscaledtrsafesolve(a, s, n, &ex, isupper, 0, isunit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      } else {
      // multiply by inv(A')
         if (!cmatrixscaledtrsafesolve(a, s, n, &ex, isupper, 2, isunit, maxgrowth)) {
            ae_frame_leave();
            return;
         }
      }
   // from 0-based to 1-based
      for (i = n - 1; i >= 0; i--) {
         ex.xC[i + 1] = ex.xC[i];
      }
   }
// Compute the estimate of the reciprocal condition number.
   if (ainvnm != 0.0) {
      *rc = 1 / ainvnm;
      *rc /= anorm;
      if (*rc < rcondthreshold()) {
         *rc = 0.0;
      }
   }
   ae_frame_leave();
}

// Triangular matrix: estimate of a condition number (1-norm)
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A       -   matrix. Array[0..N-1, 0..N-1].
//     N       -   size of A.
//     IsUpper -   True, if the matrix is upper triangular.
//     IsUnit  -   True, if the matrix has a unit diagonal.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double rmatrixtrrcond1(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit);
double rmatrixtrrcond1(RMatrix *a, ae_int_t n, bool isupper, bool isunit) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   ae_int_t j1;
   ae_int_t j2;
   double result;
   ae_frame_make(&_frame_block);
   NewVector(pivots, 0, DT_INT);
   NewVector(t, 0, DT_REAL);
   ae_assert(n >= 1, "RMatrixTRRCond1: N<1!");
   ae_vector_set_length(&t, n);
   for (i = 0; i < n; i++) {
      t.xR[i] = 0.0;
   }
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i + 1;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i - 1;
      }
      for (j = j1; j <= j2; j++) {
         t.xR[j] += fabs(a->xyR[i][j]);
      }
      if (isunit) {
         t.xR[i]++;
      } else {
         t.xR[i] += fabs(a->xyR[i][i]);
      }
   }
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      nrm = rmax2(nrm, t.xR[i]);
   }
   rcond_rmatrixrcondtrinternal(a, n, isupper, isunit, true, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Triangular matrix: estimate of a condition number (1-norm)
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A       -   matrix. Array[0..N-1, 0..N-1].
//     N       -   size of A.
//     IsUpper -   True, if the matrix is upper triangular.
//     IsUnit  -   True, if the matrix has a unit diagonal.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double cmatrixtrrcond1(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit);
double cmatrixtrrcond1(CMatrix *a, ae_int_t n, bool isupper, bool isunit) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   ae_int_t j1;
   ae_int_t j2;
   double result;
   ae_frame_make(&_frame_block);
   NewVector(pivots, 0, DT_INT);
   NewVector(t, 0, DT_REAL);
   ae_assert(n >= 1, "RMatrixTRRCond1: N<1!");
   ae_vector_set_length(&t, n);
   for (i = 0; i < n; i++) {
      t.xR[i] = 0.0;
   }
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i + 1;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i - 1;
      }
      for (j = j1; j <= j2; j++) {
         t.xR[j] += ae_c_abs(a->xyC[i][j]);
      }
      if (isunit) {
         t.xR[i]++;
      } else {
         t.xR[i] += ae_c_abs(a->xyC[i][i]);
      }
   }
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      nrm = rmax2(nrm, t.xR[i]);
   }
   rcond_cmatrixrcondtrinternal(a, n, isupper, isunit, true, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Triangular matrix: estimate of a matrix condition number (infinity-norm).
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of matrix A.
//     IsUpper -   True, if the matrix is upper triangular.
//     IsUnit  -   True, if the matrix has a unit diagonal.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double rmatrixtrrcondinf(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit);
double rmatrixtrrcondinf(RMatrix *a, ae_int_t n, bool isupper, bool isunit) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   ae_int_t j1;
   ae_int_t j2;
   double result;
   ae_frame_make(&_frame_block);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n >= 1, "RMatrixTRRCondInf: N<1!");
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i + 1;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i - 1;
      }
      v = 0.0;
      for (j = j1; j <= j2; j++) {
         v += fabs(a->xyR[i][j]);
      }
      if (isunit) {
         v++;
      } else {
         v += fabs(a->xyR[i][i]);
      }
      nrm = rmax2(nrm, v);
   }
   rcond_rmatrixrcondtrinternal(a, n, isupper, isunit, false, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Triangular matrix: estimate of a matrix condition number (infinity-norm).
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of matrix A.
//     IsUpper -   True, if the matrix is upper triangular.
//     IsUnit  -   True, if the matrix has a unit diagonal.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double cmatrixtrrcondinf(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit);
double cmatrixtrrcondinf(CMatrix *a, ae_int_t n, bool isupper, bool isunit) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   double nrm;
   ae_int_t j1;
   ae_int_t j2;
   double result;
   ae_frame_make(&_frame_block);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n >= 1, "RMatrixTRRCondInf: N<1!");
   nrm = 0.0;
   for (i = 0; i < n; i++) {
      if (isupper) {
         j1 = i + 1;
         j2 = n - 1;
      } else {
         j1 = 0;
         j2 = i - 1;
      }
      v = 0.0;
      for (j = j1; j <= j2; j++) {
         v += ae_c_abs(a->xyC[i][j]);
      }
      if (isunit) {
         v++;
      } else {
         v += ae_c_abs(a->xyC[i][i]);
      }
      nrm = rmax2(nrm, v);
   }
   rcond_cmatrixrcondtrinternal(a, n, isupper, isunit, false, nrm, &v);
   result = v;
   ae_frame_leave();
   return result;
}

// Estimate of the condition number of a matrix given by its LU decomposition (1-norm)
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     LUA         -   LU decomposition of a matrix in compact form. Output of
//                     the RMatrixLU subroutine.
//     N           -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double rmatrixlurcond1(const real_2d_array &lua, const ae_int_t n);
double rmatrixlurcond1(RMatrix *lua, ae_int_t n) {
   double v;
   double result;
   rcond_rmatrixrcondluinternal(lua, n, true, false, 0.0, &v);
   result = v;
   return result;
}

// Estimate of the condition number of a matrix given by its LU decomposition (1-norm)
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     LUA         -   LU decomposition of a matrix in compact form. Output of
//                     the CMatrixLU subroutine.
//     N           -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double cmatrixlurcond1(const complex_2d_array &lua, const ae_int_t n);
double cmatrixlurcond1(CMatrix *lua, ae_int_t n) {
   double v;
   double result;
   ae_assert(n >= 1, "CMatrixLURCond1: N<1!");
   rcond_cmatrixrcondluinternal(lua, n, true, false, 0.0, &v);
   result = v;
   return result;
}

// Estimate of the condition number of a matrix given by its LU decomposition
// (infinity norm).
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     LUA     -   LU decomposition of a matrix in compact form. Output of
//                 the RMatrixLU subroutine.
//     N       -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double rmatrixlurcondinf(const real_2d_array &lua, const ae_int_t n);
double rmatrixlurcondinf(RMatrix *lua, ae_int_t n) {
   double v;
   double result;
   rcond_rmatrixrcondluinternal(lua, n, false, false, 0.0, &v);
   result = v;
   return result;
}

// Estimate of the condition number of a matrix given by its LU decomposition
// (infinity norm).
//
// The algorithm calculates a lower bound of the condition number. In this case,
// the algorithm does not return a lower bound of the condition number, but an
// inverse number (to avoid an overflow in case of a singular matrix).
//
// Inputs:
//     LUA     -   LU decomposition of a matrix in compact form. Output of
//                 the CMatrixLU subroutine.
//     N       -   size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double cmatrixlurcondinf(const complex_2d_array &lua, const ae_int_t n);
double cmatrixlurcondinf(CMatrix *lua, ae_int_t n) {
   double v;
   double result;
   ae_assert(n >= 1, "CMatrixLURCondInf: N<1!");
   rcond_cmatrixrcondluinternal(lua, n, false, false, 0.0, &v);
   result = v;
   return result;
}

// Condition number estimate of a symmetric positive definite matrix given by
// Cholesky decomposition.
//
// The algorithm calculates a lower bound of the condition number. In this
// case, the algorithm does not return a lower bound of the condition number,
// but an inverse number (to avoid an overflow in case of a singular matrix).
//
// It should be noted that 1-norm and inf-norm condition numbers of symmetric
// matrices are equal, so the algorithm doesn't take into account the
// differences between these types of norms.
//
// Inputs:
//     CD  - Cholesky decomposition of matrix A,
//           output of SMatrixCholesky subroutine.
//     N   - size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double spdmatrixcholeskyrcond(const real_2d_array &a, const ae_int_t n, const bool isupper);
double spdmatrixcholeskyrcond(RMatrix *a, ae_int_t n, bool isupper) {
   double v;
   double result;
   rcond_spdmatrixrcondcholeskyinternal(a, n, isupper, false, 0.0, &v);
   result = v;
   return result;
}

// Condition number estimate of a Hermitian positive definite matrix given by
// Cholesky decomposition.
//
// The algorithm calculates a lower bound of the condition number. In this
// case, the algorithm does not return a lower bound of the condition number,
// but an inverse number (to avoid an overflow in case of a singular matrix).
//
// It should be noted that 1-norm and inf-norm condition numbers of symmetric
// matrices are equal, so the algorithm doesn't take into account the
// differences between these types of norms.
//
// Inputs:
//     CD  - Cholesky decomposition of matrix A,
//           output of SMatrixCholesky subroutine.
//     N   - size of matrix A.
//
// Result: 1/LowerBound(cond(A))
//
// NOTE:
//     if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
//     0.0 is returned in such cases.
// API: double hpdmatrixcholeskyrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper);
double hpdmatrixcholeskyrcond(CMatrix *a, ae_int_t n, bool isupper) {
   double v;
   double result;
   rcond_hpdmatrixrcondcholeskyinternal(a, n, isupper, false, 0.0, &v);
   result = v;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double rmatrixrcond1(const real_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixrcond1(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}

double cmatrixrcond1(const complex_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cmatrixrcond1(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}

double rmatrixrcondinf(const real_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixrcondinf(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}

double cmatrixrcondinf(const complex_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cmatrixrcondinf(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}

double spdmatrixrcond(const real_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spdmatrixrcond(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return D;
}

double hpdmatrixrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hpdmatrixrcond(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return D;
}

double rmatrixtrrcond1(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixtrrcond1(ConstT(ae_matrix, a), n, isupper, isunit);
   alglib_impl::ae_state_clear();
   return D;
}

double cmatrixtrrcond1(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cmatrixtrrcond1(ConstT(ae_matrix, a), n, isupper, isunit);
   alglib_impl::ae_state_clear();
   return D;
}

double rmatrixtrrcondinf(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixtrrcondinf(ConstT(ae_matrix, a), n, isupper, isunit);
   alglib_impl::ae_state_clear();
   return D;
}

double cmatrixtrrcondinf(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cmatrixtrrcondinf(ConstT(ae_matrix, a), n, isupper, isunit);
   alglib_impl::ae_state_clear();
   return D;
}

double rmatrixlurcond1(const real_2d_array &lua, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixlurcond1(ConstT(ae_matrix, lua), n);
   alglib_impl::ae_state_clear();
   return D;
}

double cmatrixlurcond1(const complex_2d_array &lua, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cmatrixlurcond1(ConstT(ae_matrix, lua), n);
   alglib_impl::ae_state_clear();
   return D;
}

double rmatrixlurcondinf(const real_2d_array &lua, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixlurcondinf(ConstT(ae_matrix, lua), n);
   alglib_impl::ae_state_clear();
   return D;
}

double cmatrixlurcondinf(const complex_2d_array &lua, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::cmatrixlurcondinf(ConstT(ae_matrix, lua), n);
   alglib_impl::ae_state_clear();
   return D;
}

double spdmatrixcholeskyrcond(const real_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spdmatrixcholeskyrcond(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return D;
}

double hpdmatrixcholeskyrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hpdmatrixcholeskyrcond(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === MATINV Package ===
// Depends on: RCOND
namespace alglib_impl {
// Triangular matrix inversion, recursive subroutine
//
// NOTE: this function sets Info on failure, leaves it unchanged on success.
//
// NOTE: only Tmp[Offs:Offs+N-1] is modified, other entries of the temporary array are not modified
// ALGLIB: Copyright 05.02.2010 by Sergey Bochkanov
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992.
static void matinv_rmatrixtrinverserec(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, bool isunit, RVector *tmp, ae_int_t *info) {
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t i;
   ae_int_t j;
   double v;
   double ajj;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   if (n < 1) {
      *info = -1;
      return;
   }
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (n <= tsb) {
      tscur = tsa;
   }
// Parallelism was activated if: n >= 2 * tsb && (double)n * n * n / 3.0 >= smpactivationlevel()
// Base case
   if (n <= tsa) {
      if (isupper) {
      // Compute inverse of upper triangular matrix.
         for (j = 0; j < n; j++) {
            if (!isunit) {
               if (a->xyR[offs + j][offs + j] == 0.0) {
                  *info = -3;
                  return;
               }
               a->xyR[offs + j][offs + j] = 1 / a->xyR[offs + j][offs + j];
               ajj = -a->xyR[offs + j][offs + j];
            } else {
               ajj = -1.0;
            }
         // Compute elements 1:j-1 of j-th column.
            if (j > 0) {
               ae_v_move(&tmp->xR[offs], 1, &a->xyR[offs][offs + j], a->stride, j);
               for (i = 0; i < j; i++) {
                  if (i < j - 1) {
                     v = ae_v_dotproduct(&a->xyR[offs + i][offs + i + 1], 1, &tmp->xR[offs + i + 1], 1, j - i - 1);
                  } else {
                     v = 0.0;
                  }
                  if (!isunit) {
                     a->xyR[offs + i][offs + j] = v + a->xyR[offs + i][offs + i] * tmp->xR[offs + i];
                  } else {
                     a->xyR[offs + i][offs + j] = v + tmp->xR[offs + i];
                  }
               }
               ae_v_muld(&a->xyR[offs][offs + j], a->stride, j, ajj);
            }
         }
      } else {
      // Compute inverse of lower triangular matrix.
         for (j = n - 1; j >= 0; j--) {
            if (!isunit) {
               if (a->xyR[offs + j][offs + j] == 0.0) {
                  *info = -3;
                  return;
               }
               a->xyR[offs + j][offs + j] = 1 / a->xyR[offs + j][offs + j];
               ajj = -a->xyR[offs + j][offs + j];
            } else {
               ajj = -1.0;
            }
            if (j < n - 1) {
            // Compute elements j+1:n of j-th column.
               ae_v_move(&tmp->xR[offs + j + 1], 1, &a->xyR[offs + j + 1][offs + j], a->stride, n - j - 1);
               for (i = j + 1; i < n; i++) {
                  if (i > j + 1) {
                     v = ae_v_dotproduct(&a->xyR[offs + i][offs + j + 1], 1, &tmp->xR[offs + j + 1], 1, i - j - 1);
                  } else {
                     v = 0.0;
                  }
                  if (!isunit) {
                     a->xyR[offs + i][offs + j] = v + a->xyR[offs + i][offs + i] * tmp->xR[offs + i];
                  } else {
                     a->xyR[offs + i][offs + j] = v + tmp->xR[offs + i];
                  }
               }
               ae_v_muld(&a->xyR[offs + j + 1][offs + j], a->stride, n - j - 1, ajj);
            }
         }
      }
      return;
   }
// Recursive case
   n1 = tiledsplit(n, tscur), n2 = n - n1;
// ae_int_t mn = imin2(n1, n2); //(@) Unused.
   if (n2 > 0) {
      if (isupper) {
         for (i = 0; i < n1; i++) {
            ae_v_muld(&a->xyR[offs + i][offs + n1], 1, n - n1, -1);
         }
         rmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, isupper, isunit, 0, a, offs, offs + n1);
         matinv_rmatrixtrinverserec(a, offs + n1, n2, isupper, isunit, tmp, info);
         rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, isunit, 0, a, offs, offs + n1);
      } else {
         for (i = 0; i < n2; i++) {
            ae_v_muld(&a->xyR[offs + n1 + i][offs], 1, n1, -1);
         }
         rmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, isupper, isunit, 0, a, offs + n1, offs);
         matinv_rmatrixtrinverserec(a, offs + n1, n2, isupper, isunit, tmp, info);
         rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, isunit, 0, a, offs + n1, offs);
      }
   }
   matinv_rmatrixtrinverserec(a, offs, n1, isupper, isunit, tmp, info);
}

// Triangular matrix inversion, recursive subroutine.
//
// Info is modified on failure, left unchanged on success.
// ALGLIB: Copyright 05.02.2010 by Sergey Bochkanov
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      February 29, 1992.
static void matinv_cmatrixtrinverserec(CMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, bool isunit, CVector *tmp, ae_int_t *info) {
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t i;
   ae_int_t j;
   ae_complex v;
   ae_complex ajj;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   if (n < 1) {
      *info = -1;
      return;
   }
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (n <= tsb) {
      tscur = tsa;
   }
// Parallelism was activated if: n >= 2 * tsb && 4.0 * n * n * n / 3.0 >= smpactivationlevel()
// Base case
   if (n <= tsa) {
      if (isupper) {
      // Compute inverse of upper triangular matrix.
         for (j = 0; j < n; j++) {
            if (!isunit) {
               if (ae_c_eq_d(a->xyC[offs + j][offs + j], 0.0)) {
                  *info = -3;
                  return;
               }
               a->xyC[offs + j][offs + j] = ae_c_d_div(1, a->xyC[offs + j][offs + j]);
               ajj = ae_c_neg(a->xyC[offs + j][offs + j]);
            } else {
               ajj = ae_complex_from_i(-1);
            }
         // Compute elements 1:j-1 of j-th column.
            if (j > 0) {
               ae_v_cmove(&tmp->xC[offs], 1, &a->xyC[offs][offs + j], a->stride, "N", j);
               for (i = 0; i < j; i++) {
                  if (i < j - 1) {
                     v = ae_v_cdotproduct(&a->xyC[offs + i][offs + i + 1], 1, "N", &tmp->xC[offs + i + 1], 1, "N", j - i - 1);
                  } else {
                     v = ae_complex_from_i(0);
                  }
                  if (!isunit) {
                     a->xyC[offs + i][offs + j] = ae_c_add(v, ae_c_mul(a->xyC[offs + i][offs + i], tmp->xC[offs + i]));
                  } else {
                     a->xyC[offs + i][offs + j] = ae_c_add(v, tmp->xC[offs + i]);
                  }
               }
               ae_v_cmulc(&a->xyC[offs][offs + j], a->stride, j, ajj);
            }
         }
      } else {
      // Compute inverse of lower triangular matrix.
         for (j = n - 1; j >= 0; j--) {
            if (!isunit) {
               if (ae_c_eq_d(a->xyC[offs + j][offs + j], 0.0)) {
                  *info = -3;
                  return;
               }
               a->xyC[offs + j][offs + j] = ae_c_d_div(1, a->xyC[offs + j][offs + j]);
               ajj = ae_c_neg(a->xyC[offs + j][offs + j]);
            } else {
               ajj = ae_complex_from_i(-1);
            }
            if (j < n - 1) {
            // Compute elements j+1:n of j-th column.
               ae_v_cmove(&tmp->xC[offs + j + 1], 1, &a->xyC[offs + j + 1][offs + j], a->stride, "N", n - j - 1);
               for (i = j + 1; i < n; i++) {
                  if (i > j + 1) {
                     v = ae_v_cdotproduct(&a->xyC[offs + i][offs + j + 1], 1, "N", &tmp->xC[offs + j + 1], 1, "N", i - j - 1);
                  } else {
                     v = ae_complex_from_i(0);
                  }
                  if (!isunit) {
                     a->xyC[offs + i][offs + j] = ae_c_add(v, ae_c_mul(a->xyC[offs + i][offs + i], tmp->xC[offs + i]));
                  } else {
                     a->xyC[offs + i][offs + j] = ae_c_add(v, tmp->xC[offs + i]);
                  }
               }
               ae_v_cmulc(&a->xyC[offs + j + 1][offs + j], a->stride, n - j - 1, ajj);
            }
         }
      }
      return;
   }
// Recursive case
   n1 = tiledsplit(n, tscur), n2 = n - n1;
// ae_int_t mn = imin2(n1, n2); //(@) Unused.
   if (n2 > 0) {
      if (isupper) {
         for (i = 0; i < n1; i++) {
            ae_v_cmuld(&a->xyC[offs + i][offs + n1], 1, n - n1, -1);
         }
         cmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, isupper, isunit, 0, a, offs, offs + n1);
         matinv_cmatrixtrinverserec(a, offs + n1, n2, isupper, isunit, tmp, info);
         cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, isunit, 0, a, offs, offs + n1);
      } else {
         for (i = 0; i < n2; i++) {
            ae_v_cmuld(&a->xyC[offs + n1 + i][offs], 1, n1, -1);
         }
         cmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, isupper, isunit, 0, a, offs + n1, offs);
         matinv_cmatrixtrinverserec(a, offs + n1, n2, isupper, isunit, tmp, info);
         cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, isunit, 0, a, offs + n1, offs);
      }
   }
   matinv_cmatrixtrinverserec(a, offs, n1, isupper, isunit, tmp, info);
}

static void matinv_rmatrixluinverserec(RMatrix *a, ae_int_t offs, ae_int_t n, RVector *work, ae_int_t *info, matinvreport *rep) {
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   if (n < 1) {
      *info = -1;
      return;
   }
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (n <= tsb) {
      tscur = tsa;
   }
// Parallelism was activated if: n >= 2 * tsb && 8.0 / 6.0 * n * n * n >= smpactivationlevel()
// Base case
   if (n <= tsa) {
   // Form inv(U)
      matinv_rmatrixtrinverserec(a, offs, n, true, false, work, info);
      if (*info <= 0) {
         return;
      }
   // Solve the equation inv(A)*L = inv(U) for inv(A).
      for (j = n - 1; j >= 0; j--) {
      // Copy current column of L to WORK and replace with zeros.
         for (i = j + 1; i < n; i++) {
            work->xR[i] = a->xyR[offs + i][offs + j];
            a->xyR[offs + i][offs + j] = 0.0;
         }
      // Compute current column of inv(A).
         if (j < n - 1) {
            for (i = 0; i < n; i++) {
               v = ae_v_dotproduct(&a->xyR[offs + i][offs + j + 1], 1, &work->xR[j + 1], 1, n - j - 1);
               a->xyR[offs + i][offs + j] -= v;
            }
         }
      }
      return;
   }
// Recursive code:
//
//         ( L1      )   ( U1  U12 )
// A    =  (         ) * (         )
//         ( L12  L2 )   (     U2  )
//
//         ( W   X )
// A^-1 =  (       )
//         ( Y   Z )
//
// In-place calculation can be done as follows:
// * X := inv(U1)*U12*inv(U2)
// * Y := inv(L2)*L12*inv(L1)
// * W := inv(L1*U1)+X*Y
// * X := -X*inv(L2)
// * Y := -inv(U2)*Y
// * Z := inv(L2*U2)
//
// Reordering w.r.t. interdependencies gives us:
//
// * X := inv(U1)*U12      \ suitable for parallel execution
// * Y := L12*inv(L1)      /
//
// * X := X*inv(U2)        \
// * Y := inv(L2)*Y        | suitable for parallel execution
// * W := inv(L1*U1)       /
//
// * W := W+X*Y
//
// * X := -X*inv(L2)       \ suitable for parallel execution
// * Y := -inv(U2)*Y       /
//
// * Z := inv(L2*U2)
   n1 = tiledsplit(n, tscur), n2 = n - n1;
// ae_int_t mn = imin2(n1, n2); //(@) Unused.
   ae_assert(n2 > 0, "LUInverseRec: internal error!");
// X := inv(U1)*U12
// Y := L12*inv(L1)
   rmatrixlefttrsm(n1, n2, a, offs, offs, true, false, 0, a, offs, offs + n1);
   rmatrixrighttrsm(n2, n1, a, offs, offs, false, true, 0, a, offs + n1, offs);
// X := X*inv(U2)
// Y := inv(L2)*Y
// W := inv(L1*U1)
   rmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, true, false, 0, a, offs, offs + n1);
   rmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, false, true, 0, a, offs + n1, offs);
   matinv_rmatrixluinverserec(a, offs, n1, work, info, rep);
   if (*info <= 0) {
      return;
   }
// W := W+X*Y
   rmatrixgemm(n1, n1, n2, 1.0, a, offs, offs + n1, 0, a, offs + n1, offs, 0, 1.0, a, offs, offs);
// X := -X*inv(L2)
// Y := -inv(U2)*Y
   rmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, false, true, 0, a, offs, offs + n1);
   rmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, true, false, 0, a, offs + n1, offs);
   for (i = 0; i < n1; i++) {
      ae_v_muld(&a->xyR[offs + i][offs + n1], 1, n - n1, -1);
   }
   for (i = 0; i < n2; i++) {
      ae_v_muld(&a->xyR[offs + n1 + i][offs], 1, n1, -1);
   }
// Z := inv(L2*U2)
   matinv_rmatrixluinverserec(a, offs + n1, n2, work, info, rep);
}

static void matinv_cmatrixluinverserec(CMatrix *a, ae_int_t offs, ae_int_t n, CVector *work, ae_int_t *ssinfo, matinvreport *rep) {
   ae_int_t i;
   ae_int_t j;
   ae_complex v;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   if (n < 1) {
      *ssinfo = -1;
      return;
   }
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (n <= tsb) {
      tscur = tsa;
   }
// Parallelism was activated if: n >= 2 * tsb && 32.0 * n * n * n / 6.0 >= smpactivationlevel()
// Base case
   if (n <= tsa) {
   // Form inv(U)
      matinv_cmatrixtrinverserec(a, offs, n, true, false, work, ssinfo);
      if (*ssinfo <= 0) {
         return;
      }
   // Solve the equation inv(A)*L = inv(U) for inv(A).
      for (j = n - 1; j >= 0; j--) {
      // Copy current column of L to WORK and replace with zeros.
         for (i = j + 1; i < n; i++) {
            work->xC[i] = a->xyC[offs + i][offs + j];
            a->xyC[offs + i][offs + j] = ae_complex_from_i(0);
         }
      // Compute current column of inv(A).
         if (j < n - 1) {
            for (i = 0; i < n; i++) {
               v = ae_v_cdotproduct(&a->xyC[offs + i][offs + j + 1], 1, "N", &work->xC[j + 1], 1, "N", n - j - 1);
               a->xyC[offs + i][offs + j] = ae_c_sub(a->xyC[offs + i][offs + j], v);
            }
         }
      }
      return;
   }
// Recursive code:
//
//         ( L1      )   ( U1  U12 )
// A    =  (         ) * (         )
//         ( L12  L2 )   (     U2  )
//
//         ( W   X )
// A^-1 =  (       )
//         ( Y   Z )
//
// In-place calculation can be done as follows:
// * X := inv(U1)*U12*inv(U2)
// * Y := inv(L2)*L12*inv(L1)
// * W := inv(L1*U1)+X*Y
// * X := -X*inv(L2)
// * Y := -inv(U2)*Y
// * Z := inv(L2*U2)
//
// Reordering w.r.t. interdependencies gives us:
//
// * X := inv(U1)*U12      \ suitable for parallel execution
// * Y := L12*inv(L1)      /
//
// * X := X*inv(U2)        \
// * Y := inv(L2)*Y        | suitable for parallel execution
// * W := inv(L1*U1)       /
//
// * W := W+X*Y
//
// * X := -X*inv(L2)       \ suitable for parallel execution
// * Y := -inv(U2)*Y       /
//
// * Z := inv(L2*U2)
   n1 = tiledsplit(n, tscur), n2 = n - n1;
// ae_int_t mn = imin2(n1, n2); //(@) Unused.
   ae_assert(n2 > 0, "LUInverseRec: internal error!");
// X := inv(U1)*U12
// Y := L12*inv(L1)
   cmatrixlefttrsm(n1, n2, a, offs, offs, true, false, 0, a, offs, offs + n1);
   cmatrixrighttrsm(n2, n1, a, offs, offs, false, true, 0, a, offs + n1, offs);
// X := X*inv(U2)
// Y := inv(L2)*Y
// W := inv(L1*U1)
   cmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, true, false, 0, a, offs, offs + n1);
   cmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, false, true, 0, a, offs + n1, offs);
   matinv_cmatrixluinverserec(a, offs, n1, work, ssinfo, rep);
   if (*ssinfo <= 0) {
      return;
   }
// W := W+X*Y
   cmatrixgemm(n1, n1, n2, ae_complex_from_d(1.0), a, offs, offs + n1, 0, a, offs + n1, offs, 0, ae_complex_from_d(1.0), a, offs, offs);
// X := -X*inv(L2)
// Y := -inv(U2)*Y
   cmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, false, true, 0, a, offs, offs + n1);
   cmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, true, false, 0, a, offs + n1, offs);
   for (i = 0; i < n1; i++) {
      ae_v_cmuld(&a->xyC[offs + i][offs + n1], 1, n - n1, -1);
   }
   for (i = 0; i < n2; i++) {
      ae_v_cmuld(&a->xyC[offs + n1 + i][offs], 1, n1, -1);
   }
// Z := inv(L2*U2)
   matinv_cmatrixluinverserec(a, offs + n1, n2, work, ssinfo, rep);
}

// Inversion of a matrix given by its LU decomposition.
//
// Inputs:
//     A       -   LU decomposition of the matrix
//                 (output of RMatrixLU subroutine).
//     Pivots  -   table of permutations
//                 (the output of RMatrixLU subroutine).
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//
// Outputs:
//     Info    -   return code:
//                 * -3    A is singular, or VERY close to singular.
//                         it is filled by zeros in such cases.
//                 *  1    task is solved (but matrix A may be ill-conditioned,
//                         check R1/RInf parameters for condition numbers).
//     Rep     -   solver report, see below for more info
//     A       -   inverse of matrix A.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//
// SOLVER REPORT
//
// Subroutine sets following fields of the Rep structure:
// * R1        reciprocal of condition number: 1/cond(A), 1-norm.
// * RInf      reciprocal of condition number: 1/cond(A), inf-norm.
//
// ALGLIB Routine: Copyright 05.02.2010 by Sergey Bochkanov
// API: void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep);
// API: void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep);
void rmatrixluinverse(RMatrix *a, ZVector *pivots, ae_int_t n, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t sinfo;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(work, 0, DT_REAL);
   ae_assert(n > 0, "RMatrixLUInverse: N <= 0!");
   ae_assert(a->cols >= n, "RMatrixLUInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "RMatrixLUInverse: rows(A)<N!");
   ae_assert(pivots->cnt >= n, "RMatrixLUInverse: len(Pivots)<N!");
   ae_assert(apservisfinitematrix(a, n, n), "RMatrixLUInverse: A contains infinite or NaN values!");
   *info = 1;
   for (i = 0; i < n; i++) {
      if (pivots->xZ[i] > n - 1 || pivots->xZ[i] < i) {
         *info = -1;
      }
   }
   ae_assert(*info > 0, "RMatrixLUInverse: incorrect Pivots array!");
// calculate condition numbers
   rep->r1 = rmatrixlurcond1(a, n);
   rep->rinf = rmatrixlurcondinf(a, n);
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            a->xyR[i][j] = 0.0;
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
// Call cache-oblivious code
   ae_vector_set_length(&work, n);
   sinfo = 1;
   matinv_rmatrixluinverserec(a, 0, n, &work, &sinfo, rep);
   *info = sinfo;
// apply permutations
   for (i = 0; i < n; i++) {
      for (j = n - 2; j >= 0; j--) {
         swapr(&a->xyR[i][j], &a->xyR[i][pivots->xZ[j]]);
      }
   }
   ae_frame_leave();
}

// Inversion of a matrix given by its LU decomposition.
//
// Inputs:
//     A       -   LU decomposition of the matrix
//                 (output of CMatrixLU subroutine).
//     Pivots  -   table of permutations
//                 (the output of CMatrixLU subroutine).
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
//
// ALGLIB Routine: Copyright 05.02.2010 by Sergey Bochkanov
// API: void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep);
// API: void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep);
void cmatrixluinverse(CMatrix *a, ZVector *pivots, ae_int_t n, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t sinfo;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(work, 0, DT_COMPLEX);
   ae_assert(n > 0, "CMatrixLUInverse: N <= 0!");
   ae_assert(a->cols >= n, "CMatrixLUInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "CMatrixLUInverse: rows(A)<N!");
   ae_assert(pivots->cnt >= n, "CMatrixLUInverse: len(Pivots)<N!");
   ae_assert(apservisfinitecmatrix(a, n, n), "CMatrixLUInverse: A contains infinite or NaN values!");
   *info = 1;
   for (i = 0; i < n; i++) {
      if (pivots->xZ[i] > n - 1 || pivots->xZ[i] < i) {
         *info = -1;
      }
   }
   ae_assert(*info > 0, "CMatrixLUInverse: incorrect Pivots array!");
// calculate condition numbers
   rep->r1 = cmatrixlurcond1(a, n);
   rep->rinf = cmatrixlurcondinf(a, n);
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            a->xyC[i][j] = ae_complex_from_i(0);
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
// Call cache-oblivious code
   ae_vector_set_length(&work, n);
   sinfo = 1;
   matinv_cmatrixluinverserec(a, 0, n, &work, &sinfo, rep);
   *info = sinfo;
// apply permutations
   for (i = 0; i < n; i++) {
      for (j = n - 2; j >= 0; j--) {
         swapc(&a->xyC[i][j], &a->xyC[i][pivots->xZ[j]]);
      }
   }
   ae_frame_leave();
}

// Inversion of a general matrix.
//
// Inputs:
//     A       -   matrix.
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
//
// Result:
//     True, if the matrix is not singular.
//     False, if the matrix is singular.
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixinverse(real_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep);
// API: void rmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
void rmatrixinverse(RMatrix *a, ae_int_t n, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n > 0, "RMatrixInverse: N <= 0!");
   ae_assert(a->cols >= n, "RMatrixInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "RMatrixInverse: rows(A)<N!");
   ae_assert(apservisfinitematrix(a, n, n), "RMatrixInverse: A contains infinite or NaN values!");
   rmatrixlu(a, n, n, &pivots);
   rmatrixluinverse(a, &pivots, n, info, rep);
   ae_frame_leave();
}

// Inversion of a general matrix.
//
// Inputs:
//     A       -   matrix
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: void cmatrixinverse(complex_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep);
// API: void cmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
void cmatrixinverse(CMatrix *a, ae_int_t n, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n > 0, "CRMatrixInverse: N <= 0!");
   ae_assert(a->cols >= n, "CRMatrixInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "CRMatrixInverse: rows(A)<N!");
   ae_assert(apservisfinitecmatrix(a, n, n), "CMatrixInverse: A contains infinite or NaN values!");
   cmatrixlu(a, n, n, &pivots);
   cmatrixluinverse(a, &pivots, n, info, rep);
   ae_frame_leave();
}

// Recursive subroutine for SPD inversion.
//
// NOTE: this function expects that matris is strictly positive-definite.
//
// ALGLIB Routine: Copyright 10.02.2010 by Sergey Bochkanov
static void matinv_spdmatrixcholeskyinverserec(RMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, RVector *tmp) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   ae_int_t sinfo2;
   ae_frame_make(&_frame_block);
   if (n < 1) {
      ae_frame_leave();
      return;
   }
   tsa = matrixtilesizea();
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (n <= tsb) {
      tscur = tsa;
   }
// Base case
   if (n <= tsa) {
      sinfo2 = 1;
      matinv_rmatrixtrinverserec(a, offs, n, isupper, false, tmp, &sinfo2);
      ae_assert(sinfo2 > 0, "SPDMatrixCholeskyInverseRec: integrity check failed");
      if (isupper) {
      // Compute the product U * U'.
      // NOTE: we never assume that diagonal of U is real
         for (i = 0; i < n; i++) {
            if (i == 0) {
            // 1x1 matrix
               a->xyR[offs + i][offs + i] = ae_sqr(a->xyR[offs + i][offs + i]);
            } else {
            // (I+1)x(I+1) matrix,
            //
            // ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
            // (          ) * (              ) = (                                )
            // (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
            //
            // A11 is IxI, A22 is 1x1.
               ae_v_move(tmp->xR, 1, &a->xyR[offs][offs + i], a->stride, i);
               for (j = 0; j < i; j++) {
                  v = a->xyR[offs + j][offs + i];
                  ae_v_addd(&a->xyR[offs + j][offs + j], 1, &tmp->xR[j], 1, i - j, v);
               }
               v = a->xyR[offs + i][offs + i];
               ae_v_muld(&a->xyR[offs][offs + i], a->stride, i, v);
               a->xyR[offs + i][offs + i] = ae_sqr(a->xyR[offs + i][offs + i]);
            }
         }
      } else {
      // Compute the product L' * L
      // NOTE: we never assume that diagonal of L is real
         for (i = 0; i < n; i++) {
            if (i == 0) {
            // 1x1 matrix
               a->xyR[offs + i][offs + i] = ae_sqr(a->xyR[offs + i][offs + i]);
            } else {
            // (I+1)x(I+1) matrix,
            //
            // ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
            // (              ) * (          ) = (                                )
            // (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
            //
            // A11 is IxI, A22 is 1x1.
               ae_v_move(tmp->xR, 1, &a->xyR[offs + i][offs], 1, i);
               for (j = 0; j < i; j++) {
                  v = a->xyR[offs + i][offs + j];
                  ae_v_addd(&a->xyR[offs + j][offs], 1, tmp->xR, 1, j + 1, v);
               }
               v = a->xyR[offs + i][offs + i];
               ae_v_muld(&a->xyR[offs + i][offs], 1, i, v);
               a->xyR[offs + i][offs + i] = ae_sqr(a->xyR[offs + i][offs + i]);
            }
         }
      }
      ae_frame_leave();
      return;
   }
// Recursive code: triangular factor inversion merged with
// UU' or L'L multiplication
   n1 = tiledsplit(n, tscur), n2 = n - n1;
// form off-diagonal block of trangular inverse
   if (isupper) {
      for (i = 0; i < n1; i++) {
         ae_v_muld(&a->xyR[offs + i][offs + n1], 1, n - n1, -1);
      }
      rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 0, a, offs, offs + n1);
      rmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, isupper, false, 0, a, offs, offs + n1);
   } else {
      for (i = 0; i < n2; i++) {
         ae_v_muld(&a->xyR[offs + n1 + i][offs], 1, n1, -1);
      }
      rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 0, a, offs + n1, offs);
      rmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, isupper, false, 0, a, offs + n1, offs);
   }
// invert first diagonal block
   matinv_spdmatrixcholeskyinverserec(a, offs, n1, isupper, tmp);
// update first diagonal block with off-diagonal block,
// update off-diagonal block
   if (isupper) {
      rmatrixsyrk(n1, n2, 1.0, a, offs, offs + n1, 0, 1.0, a, offs, offs, isupper);
      rmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, isupper, false, 1, a, offs, offs + n1);
   } else {
      rmatrixsyrk(n1, n2, 1.0, a, offs + n1, offs, 1, 1.0, a, offs, offs, isupper);
      rmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, isupper, false, 1, a, offs + n1, offs);
   }
// invert second diagonal block
   matinv_spdmatrixcholeskyinverserec(a, offs + n1, n2, isupper, tmp);
   ae_frame_leave();
}

// Recursive subroutine for HPD inversion.
//
// ALGLIB Routine: Copyright 10.02.2010 by Sergey Bochkanov
static void matinv_hpdmatrixcholeskyinverserec(CMatrix *a, ae_int_t offs, ae_int_t n, bool isupper, CVector *tmp) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_complex v;
   ae_int_t n1;
   ae_int_t n2;
   ae_int_t tsa;
   ae_int_t tsb;
   ae_int_t tscur;
   ae_int_t sinfo;
   ae_frame_make(&_frame_block);
   if (n < 1) {
      ae_frame_leave();
      return;
   }
   tsa = matrixtilesizea() / 2;
   tsb = matrixtilesizeb();
   tscur = tsb;
   if (n <= tsb) {
      tscur = tsa;
   }
// Base case
   if (n <= tsa) {
      sinfo = 1;
      matinv_cmatrixtrinverserec(a, offs, n, isupper, false, tmp, &sinfo);
      ae_assert(sinfo > 0, "HPDMatrixCholeskyInverseRec: integrity check failed");
      if (isupper) {
      // Compute the product U * U'.
      // NOTE: we never assume that diagonal of U is real
         for (i = 0; i < n; i++) {
            if (i == 0) {
            // 1x1 matrix
               a->xyC[offs + i][offs + i] = ae_complex_from_d(ae_sqr(a->xyC[offs + i][offs + i].x) + ae_sqr(a->xyC[offs + i][offs + i].y));
            } else {
            // (I+1)x(I+1) matrix,
            //
            // ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
            // (          ) * (              ) = (                                )
            // (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
            //
            // A11 is IxI, A22 is 1x1.
               ae_v_cmove(tmp->xC, 1, &a->xyC[offs][offs + i], a->stride, "Conj", i);
               for (j = 0; j < i; j++) {
                  v = a->xyC[offs + j][offs + i];
                  ae_v_caddc(&a->xyC[offs + j][offs + j], 1, &tmp->xC[j], 1, "N", i - j, v);
               }
               v = ae_c_conj(a->xyC[offs + i][offs + i]);
               ae_v_cmulc(&a->xyC[offs][offs + i], a->stride, i, v);
               a->xyC[offs + i][offs + i] = ae_complex_from_d(ae_sqr(a->xyC[offs + i][offs + i].x) + ae_sqr(a->xyC[offs + i][offs + i].y));
            }
         }
      } else {
      // Compute the product L' * L
      // NOTE: we never assume that diagonal of L is real
         for (i = 0; i < n; i++) {
            if (i == 0) {
            // 1x1 matrix
               a->xyC[offs + i][offs + i] = ae_complex_from_d(ae_sqr(a->xyC[offs + i][offs + i].x) + ae_sqr(a->xyC[offs + i][offs + i].y));
            } else {
            // (I+1)x(I+1) matrix,
            //
            // ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
            // (              ) * (          ) = (                                )
            // (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
            //
            // A11 is IxI, A22 is 1x1.
               ae_v_cmove(tmp->xC, 1, &a->xyC[offs + i][offs], 1, "N", i);
               for (j = 0; j < i; j++) {
                  v = ae_c_conj(a->xyC[offs + i][offs + j]);
                  ae_v_caddc(&a->xyC[offs + j][offs], 1, tmp->xC, 1, "N", j + 1, v);
               }
               v = ae_c_conj(a->xyC[offs + i][offs + i]);
               ae_v_cmulc(&a->xyC[offs + i][offs], 1, i, v);
               a->xyC[offs + i][offs + i] = ae_complex_from_d(ae_sqr(a->xyC[offs + i][offs + i].x) + ae_sqr(a->xyC[offs + i][offs + i].y));
            }
         }
      }
      ae_frame_leave();
      return;
   }
// Recursive code: triangular factor inversion merged with
// UU' or L'L multiplication
   n1 = tiledsplit(n, tscur), n2 = n - n1;
// form off-diagonal block of trangular inverse
   if (isupper) {
      for (i = 0; i < n1; i++) {
         ae_v_cmuld(&a->xyC[offs + i][offs + n1], 1, n - n1, -1);
      }
      cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 0, a, offs, offs + n1);
      cmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, isupper, false, 0, a, offs, offs + n1);
   } else {
      for (i = 0; i < n2; i++) {
         ae_v_cmuld(&a->xyC[offs + n1 + i][offs], 1, n1, -1);
      }
      cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 0, a, offs + n1, offs);
      cmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, isupper, false, 0, a, offs + n1, offs);
   }
// invert first diagonal block
   matinv_hpdmatrixcholeskyinverserec(a, offs, n1, isupper, tmp);
// update first diagonal block with off-diagonal block,
// update off-diagonal block
   if (isupper) {
      cmatrixherk(n1, n2, 1.0, a, offs, offs + n1, 0, 1.0, a, offs, offs, isupper);
      cmatrixrighttrsm(n1, n2, a, offs + n1, offs + n1, isupper, false, 2, a, offs, offs + n1);
   } else {
      cmatrixherk(n1, n2, 1.0, a, offs + n1, offs, 2, 1.0, a, offs, offs, isupper);
      cmatrixlefttrsm(n2, n1, a, offs + n1, offs + n1, isupper, false, 2, a, offs + n1, offs);
   }
// invert second diagonal block
   matinv_hpdmatrixcholeskyinverserec(a, offs + n1, n2, isupper, tmp);
   ae_frame_leave();
}

// Inversion of a symmetric positive definite matrix which is given
// by Cholesky decomposition.
//
// Inputs:
//     A       -   Cholesky decomposition of the matrix to be inverted:
//                 A=U'*U or A = L*L'.
//                 Output of  SPDMatrixCholesky subroutine.
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//     IsUpper -   storage type (optional):
//                 * if True, symmetric  matrix  A  is  given  by  its  upper
//                   triangle, and the lower triangle isn't  used/changed  by
//                   function
//                 * if False,  symmetric matrix  A  is  given  by  its lower
//                   triangle, and the  upper triangle isn't used/changed  by
//                   function
//                 * if not given, lower half is used.
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
//
// ALGLIB Routine: Copyright 10.02.2010 by Sergey Bochkanov
// API: void spdmatrixcholeskyinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
// API: void spdmatrixcholeskyinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
void spdmatrixcholeskyinverse(RMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   bool f;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(tmp, 0, DT_REAL);
   NewObj(matinvreport, rep2);
   ae_assert(n > 0, "SPDMatrixCholeskyInverse: N <= 0!");
   ae_assert(a->cols >= n, "SPDMatrixCholeskyInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "SPDMatrixCholeskyInverse: rows(A)<N!");
   *info = 1;
   f = true;
   for (i = 0; i < n; i++) {
      f = f && isfinite(a->xyR[i][i]);
   }
   ae_assert(f, "SPDMatrixCholeskyInverse: A contains infinite or NaN values!");
// calculate condition numbers
   rep->r1 = spdmatrixcholeskyrcond(a, n, isupper);
   rep->rinf = rep->r1;
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      if (isupper) {
         for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
      } else {
         for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
// Inverse
   ae_vector_set_length(&tmp, n);
   matinv_spdmatrixcholeskyinverserec(a, 0, n, isupper, &tmp);
   ae_frame_leave();
}

// Inversion of a Hermitian positive definite matrix which is given
// by Cholesky decomposition.
//
// Inputs:
//     A       -   Cholesky decomposition of the matrix to be inverted:
//                 A=U'*U or A = L*L'.
//                 Output of  HPDMatrixCholesky subroutine.
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//     IsUpper -   storage type (optional):
//                 * if True, symmetric  matrix  A  is  given  by  its  upper
//                   triangle, and the lower triangle isn't  used/changed  by
//                   function
//                 * if False,  symmetric matrix  A  is  given  by  its lower
//                   triangle, and the  upper triangle isn't used/changed  by
//                   function
//                 * if not given, lower half is used.
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
//
// ALGLIB Routine: Copyright 10.02.2010 by Sergey Bochkanov
// API: void hpdmatrixcholeskyinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
// API: void hpdmatrixcholeskyinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
void hpdmatrixcholeskyinverse(CMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   bool f;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewObj(matinvreport, rep2);
   NewVector(tmp, 0, DT_COMPLEX);
   ae_assert(n > 0, "HPDMatrixCholeskyInverse: N <= 0!");
   ae_assert(a->cols >= n, "HPDMatrixCholeskyInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "HPDMatrixCholeskyInverse: rows(A)<N!");
   f = true;
   for (i = 0; i < n; i++) {
      f = f && isfinite(a->xyC[i][i].x) && isfinite(a->xyC[i][i].y);
   }
   ae_assert(f, "HPDMatrixCholeskyInverse: A contains infinite or NaN values!");
   *info = 1;
// calculate condition numbers
   rep->r1 = hpdmatrixcholeskyrcond(a, n, isupper);
   rep->rinf = rep->r1;
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      if (isupper) {
         for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
               a->xyC[i][j] = ae_complex_from_i(0);
            }
         }
      } else {
         for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
               a->xyC[i][j] = ae_complex_from_i(0);
            }
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
// Inverse
   ae_vector_set_length(&tmp, n);
   matinv_hpdmatrixcholeskyinverserec(a, 0, n, isupper, &tmp);
   ae_frame_leave();
}

// Inversion of a symmetric positive definite matrix.
//
// Given an upper or lower triangle of a symmetric positive definite matrix,
// the algorithm generates matrix A^-1 and saves the upper or lower triangle
// depending on the input.
//
// Inputs:
//     A       -   matrix to be inverted (upper or lower triangle).
//                 Array with elements [0..N-1,0..N-1].
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//     IsUpper -   storage type (optional):
//                 * if True, symmetric  matrix  A  is  given  by  its  upper
//                   triangle, and the lower triangle isn't  used/changed  by
//                   function
//                 * if False,  symmetric matrix  A  is  given  by  its lower
//                   triangle, and the  upper triangle isn't used/changed  by
//                   function
//                 * if not given,  both lower and upper  triangles  must  be
//                   filled.
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
//
// ALGLIB Routine: Copyright 10.02.2010 by Sergey Bochkanov
// API: void spdmatrixinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
// API: void spdmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
void spdmatrixinverse(RMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep) {
   *info = 0;
   SetObj(matinvreport, rep);
   ae_assert(n > 0, "SPDMatrixInverse: N <= 0!");
   ae_assert(a->cols >= n, "SPDMatrixInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "SPDMatrixInverse: rows(A)<N!");
   ae_assert(isfinitertrmatrix(a, n, isupper), "SPDMatrixInverse: A contains infinite or NaN values!");
   *info = 1;
   if (spdmatrixcholesky(a, n, isupper)) {
      spdmatrixcholeskyinverse(a, n, isupper, info, rep);
   } else {
      *info = -3;
   }
}

// Inversion of a Hermitian positive definite matrix.
//
// Given an upper or lower triangle of a Hermitian positive definite matrix,
// the algorithm generates matrix A^-1 and saves the upper or lower triangle
// depending on the input.
//
// Inputs:
//     A       -   matrix to be inverted (upper or lower triangle).
//                 Array with elements [0..N-1,0..N-1].
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//     IsUpper -   storage type (optional):
//                 * if True, symmetric  matrix  A  is  given  by  its  upper
//                   triangle, and the lower triangle isn't  used/changed  by
//                   function
//                 * if False,  symmetric matrix  A  is  given  by  its lower
//                   triangle, and the  upper triangle isn't used/changed  by
//                   function
//                 * if not given,  both lower and upper  triangles  must  be
//                   filled.
//
// Outputs:
//     Info    -   return code, same as in RMatrixLUInverse
//     Rep     -   solver report, same as in RMatrixLUInverse
//     A       -   inverse of matrix A, same as in RMatrixLUInverse
//
// ALGLIB Routine: Copyright 10.02.2010 by Sergey Bochkanov
// API: void hpdmatrixinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
// API: void hpdmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
void hpdmatrixinverse(CMatrix *a, ae_int_t n, bool isupper, ae_int_t *info, matinvreport *rep) {
   *info = 0;
   SetObj(matinvreport, rep);
   ae_assert(n > 0, "HPDMatrixInverse: N <= 0!");
   ae_assert(a->cols >= n, "HPDMatrixInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "HPDMatrixInverse: rows(A)<N!");
   ae_assert(apservisfinitectrmatrix(a, n, isupper), "HPDMatrixInverse: A contains infinite or NaN values!");
   *info = 1;
   if (hpdmatrixcholesky(a, n, isupper)) {
      hpdmatrixcholeskyinverse(a, n, isupper, info, rep);
   } else {
      *info = -3;
   }
}

// Triangular matrix inverse (real)
//
// The subroutine inverts the following types of matrices:
//     * upper triangular
//     * upper triangular with unit diagonal
//     * lower triangular
//     * lower triangular with unit diagonal
//
// In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
// also be upper (lower) triangular, and after the end of the algorithm,  the
// inverse matrix replaces the source matrix. The elements  below (above) the
// main diagonal are not changed by the algorithm.
//
// If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
// diagonal, and the diagonal elements are not passed to the algorithm.
//
// Inputs:
//     A       -   matrix, array[0..N-1, 0..N-1].
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//     IsUpper -   True, if the matrix is upper triangular.
//     IsUnit  -   diagonal type (optional):
//                 * if True, matrix has unit diagonal (a[i,i] are NOT used)
//                 * if False, matrix diagonal is arbitrary
//                 * if not given, False is assumed
//
// Outputs:
//     Info    -   same as for RMatrixLUInverse
//     Rep     -   same as for RMatrixLUInverse
//     A       -   same as for RMatrixLUInverse.
// ALGLIB: Copyright 05.02.2010 by Sergey Bochkanov
// API: void rmatrixtrinverse(real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep);
// API: void rmatrixtrinverse(real_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep);
void rmatrixtrinverse(RMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t sinfo;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(tmp, 0, DT_REAL);
   ae_assert(n > 0, "RMatrixTRInverse: N <= 0!");
   ae_assert(a->cols >= n, "RMatrixTRInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "RMatrixTRInverse: rows(A)<N!");
   ae_assert(isfinitertrmatrix(a, n, isupper), "RMatrixTRInverse: A contains infinite or NaN values!");
// calculate condition numbers
   rep->r1 = rmatrixtrrcond1(a, n, isupper, isunit);
   rep->rinf = rmatrixtrrcondinf(a, n, isupper, isunit);
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            a->xyR[i][j] = 0.0;
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
// Invert
   ae_vector_set_length(&tmp, n);
   sinfo = 1;
   matinv_rmatrixtrinverserec(a, 0, n, isupper, isunit, &tmp, &sinfo);
   *info = sinfo;
   ae_frame_leave();
}

// Triangular matrix inverse (complex)
//
// The subroutine inverts the following types of matrices:
//     * upper triangular
//     * upper triangular with unit diagonal
//     * lower triangular
//     * lower triangular with unit diagonal
//
// In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
// also be upper (lower) triangular, and after the end of the algorithm,  the
// inverse matrix replaces the source matrix. The elements  below (above) the
// main diagonal are not changed by the algorithm.
//
// If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
// diagonal, and the diagonal elements are not passed to the algorithm.
//
// Inputs:
//     A       -   matrix, array[0..N-1, 0..N-1].
//     N       -   size of matrix A (optional) :
//                 * if given, only principal NxN submatrix is processed  and
//                   overwritten. other elements are unchanged.
//                 * if not given,  size  is  automatically  determined  from
//                   matrix size (A must be square matrix)
//     IsUpper -   True, if the matrix is upper triangular.
//     IsUnit  -   diagonal type (optional):
//                 * if True, matrix has unit diagonal (a[i,i] are NOT used)
//                 * if False, matrix diagonal is arbitrary
//                 * if not given, False is assumed
//
// Outputs:
//     Info    -   same as for RMatrixLUInverse
//     Rep     -   same as for RMatrixLUInverse
//     A       -   same as for RMatrixLUInverse.
// ALGLIB: Copyright 05.02.2010 by Sergey Bochkanov
// API: void cmatrixtrinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep);
// API: void cmatrixtrinverse(complex_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep);
void cmatrixtrinverse(CMatrix *a, ae_int_t n, bool isupper, bool isunit, ae_int_t *info, matinvreport *rep) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_int_t sinfo;
   ae_frame_make(&_frame_block);
   *info = 0;
   SetObj(matinvreport, rep);
   NewVector(tmp, 0, DT_COMPLEX);
   ae_assert(n > 0, "CMatrixTRInverse: N <= 0!");
   ae_assert(a->cols >= n, "CMatrixTRInverse: cols(A)<N!");
   ae_assert(a->rows >= n, "CMatrixTRInverse: rows(A)<N!");
   ae_assert(apservisfinitectrmatrix(a, n, isupper), "CMatrixTRInverse: A contains infinite or NaN values!");
// calculate condition numbers
   rep->r1 = cmatrixtrrcond1(a, n, isupper, isunit);
   rep->rinf = cmatrixtrrcondinf(a, n, isupper, isunit);
   if (rep->r1 < rcondthreshold() || rep->rinf < rcondthreshold()) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            a->xyC[i][j] = ae_complex_from_i(0);
         }
      }
      rep->r1 = 0.0;
      rep->rinf = 0.0;
      *info = -3;
      ae_frame_leave();
      return;
   }
// Invert
   ae_vector_set_length(&tmp, n);
   sinfo = 1;
   matinv_cmatrixtrinverserec(a, 0, n, isupper, isunit, &tmp, &sinfo);
   *info = sinfo;
   ae_frame_leave();
}

void matinvreport_init(void *_p, bool make_automatic) {
}

void matinvreport_copy(void *_dst, void *_src, bool make_automatic) {
   matinvreport *dst = (matinvreport *)_dst;
   matinvreport *src = (matinvreport *)_src;
   dst->r1 = src->r1;
   dst->rinf = src->rinf;
}

void matinvreport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
// Matrix inverse report:
// * R1    reciprocal of condition number in 1-norm
// * RInf  reciprocal of condition number in inf-norm
DefClass(matinvreport, AndD DecVal(r1) AndD DecVal(rinf))

void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixluinverse(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows() || a.cols() != pivots.length()) ThrowError("Error while calling 'rmatrixluinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixluinverse(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixluinverse(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows() || a.cols() != pivots.length()) ThrowError("Error while calling 'cmatrixluinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixluinverse(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void rmatrixinverse(real_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixinverse(ConstT(ae_matrix, a), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void rmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'rmatrixinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixinverse(ConstT(ae_matrix, a), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void cmatrixinverse(complex_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixinverse(ConstT(ae_matrix, a), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void cmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'cmatrixinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixinverse(ConstT(ae_matrix, a), n, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void spdmatrixcholeskyinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskyinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void spdmatrixcholeskyinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'spdmatrixcholeskyinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   bool isupper = false;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixcholeskyinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void hpdmatrixcholeskyinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixcholeskyinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void hpdmatrixcholeskyinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'hpdmatrixcholeskyinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   bool isupper = false;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixcholeskyinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void spdmatrixinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void spdmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'spdmatrixinverse': looks like one of arguments has wrong size");
   if (!alglib_impl::ae_is_symmetric(ConstT(ae_matrix, a))) ThrowError("'a' parameter is not symmetric matrix");
   ae_int_t n = a.cols();
   bool isupper = false;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::spdmatrixinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   if (!alglib_impl::ae_force_symmetric(ConstT(ae_matrix, a))) ThrowError("Internal error while forcing symmetricity of 'a' parameter");
   alglib_impl::ae_state_clear();
}
#endif

void hpdmatrixinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void hpdmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'hpdmatrixinverse': looks like one of arguments has wrong size");
   if (!alglib_impl::ae_is_hermitian(ConstT(ae_matrix, a))) ThrowError("'a' parameter is not Hermitian matrix");
   ae_int_t n = a.cols();
   bool isupper = false;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hpdmatrixinverse(ConstT(ae_matrix, a), n, isupper, &info, ConstT(matinvreport, rep));
   if (!alglib_impl::ae_force_hermitian(ConstT(ae_matrix, a))) ThrowError("Internal error while forcing Hermitian properties of 'a' parameter");
   alglib_impl::ae_state_clear();
}
#endif

void rmatrixtrinverse(real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixtrinverse(ConstT(ae_matrix, a), n, isupper, isunit, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void rmatrixtrinverse(real_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'rmatrixtrinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   bool isunit = false;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixtrinverse(ConstT(ae_matrix, a), n, isupper, isunit, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif

void cmatrixtrinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixtrinverse(ConstT(ae_matrix, a), n, isupper, isunit, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#if !defined AE_NO_EXCEPTIONS
void cmatrixtrinverse(complex_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep) {
   if (a.cols() != a.rows()) ThrowError("Error while calling 'cmatrixtrinverse': looks like one of arguments has wrong size");
   ae_int_t n = a.cols();
   bool isunit = false;
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixtrinverse(ConstT(ae_matrix, a), n, isupper, isunit, &info, ConstT(matinvreport, rep));
   alglib_impl::ae_state_clear();
}
#endif
} // end of namespace alglib

// === ORTFAC Package ===
// Depends on: (AlgLibInternal) CREFLECTIONS, HBLAS, SBLAS
// Depends on: (AlgLibMisc) HQRND
// Depends on: ABLAS
namespace alglib_impl {
// Generate block reflector:
// * fill unused parts of reflectors matrix by zeros
// * fill diagonal of reflectors matrix by ones
// * generate triangular factor T
//
// PARAMETERS:
//     A           -   either LengthA*BlockSize (if ColumnwiseA) or
//                     BlockSize*LengthA (if not ColumnwiseA) matrix of
//                     elementary reflectors.
//                     Modified on exit.
//     Tau         -   scalar factors
//     ColumnwiseA -   reflectors are stored in rows or in columns
//     LengthA     -   length of largest reflector
//     BlockSize   -   number of reflectors
//     T           -   array[BlockSize,2*BlockSize]. Left BlockSize*BlockSize
//                     submatrix stores triangular factor on exit.
//     WORK        -   array[BlockSize]
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
static void ortfac_rmatrixblockreflector(RMatrix *a, RVector *tau, bool columnwisea, ae_int_t lengtha, ae_int_t blocksize, RMatrix *t, RVector *work) {
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double v;
// fill beginning of new column with zeros,
// load 1.0 in the first non-zero element
   for (k = 0; k < blocksize; k++) {
      if (columnwisea) {
         for (i = 0; i < k; i++) {
            a->xyR[i][k] = 0.0;
         }
      } else {
         for (i = 0; i < k; i++) {
            a->xyR[k][i] = 0.0;
         }
      }
      a->xyR[k][k] = 1.0;
   }
// Calculate Gram matrix of A
   for (i = 0; i < blocksize; i++) {
      for (j = 0; j < blocksize; j++) {
         t->xyR[i][blocksize + j] = 0.0;
      }
   }
   for (k = 0; k < lengtha; k++) {
      for (j = 1; j < blocksize; j++) {
         if (columnwisea) {
            v = a->xyR[k][j];
            if (v != 0.0) {
               ae_v_addd(&t->xyR[j][blocksize], 1, a->xyR[k], 1, j, v);
            }
         } else {
            v = a->xyR[j][k];
            if (v != 0.0) {
               ae_v_addd(&t->xyR[j][blocksize], 1, &a->xyR[0][k], a->stride, j, v);
            }
         }
      }
   }
// Prepare Y (stored in TmpA) and T (stored in TmpT)
   for (k = 0; k < blocksize; k++) {
   // fill non-zero part of T, use pre-calculated Gram matrix
      ae_v_move(work->xR, 1, &t->xyR[k][blocksize], 1, k);
      for (i = 0; i < k; i++) {
         v = ae_v_dotproduct(&t->xyR[i][i], 1, &work->xR[i], 1, k - i);
         t->xyR[i][k] = -tau->xR[k] * v;
      }
      t->xyR[k][k] = -tau->xR[k];
   // Rest of T is filled by zeros
      for (i = k + 1; i < blocksize; i++) {
         t->xyR[i][k] = 0.0;
      }
   }
}

// Generate block reflector (complex):
// * fill unused parts of reflectors matrix by zeros
// * fill diagonal of reflectors matrix by ones
// * generate triangular factor T
//
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
static void ortfac_cmatrixblockreflector(CMatrix *a, CVector *tau, bool columnwisea, ae_int_t lengtha, ae_int_t blocksize, CMatrix *t, CVector *work) {
   ae_int_t i;
   ae_int_t k;
   ae_complex v;
// Prepare Y (stored in TmpA) and T (stored in TmpT)
   for (k = 0; k < blocksize; k++) {
   // fill beginning of new column with zeros,
   // load 1.0 in the first non-zero element
      if (columnwisea) {
         for (i = 0; i < k; i++) {
            a->xyC[i][k] = ae_complex_from_i(0);
         }
      } else {
         for (i = 0; i < k; i++) {
            a->xyC[k][i] = ae_complex_from_i(0);
         }
      }
      a->xyC[k][k] = ae_complex_from_i(1);
   // fill non-zero part of T,
      for (i = 0; i < k; i++) {
         if (columnwisea) {
            v = ae_v_cdotproduct(&a->xyC[k][i], a->stride, "Conj", &a->xyC[k][k], a->stride, "N", lengtha - k);
         } else {
            v = ae_v_cdotproduct(&a->xyC[i][k], 1, "N", &a->xyC[k][k], 1, "Conj", lengtha - k);
         }
         work->xC[i] = v;
      }
      for (i = 0; i < k; i++) {
         v = ae_v_cdotproduct(&t->xyC[i][i], 1, "N", &work->xC[i], 1, "N", k - i);
         t->xyC[i][k] = ae_c_neg(ae_c_mul(tau->xC[k], v));
      }
      t->xyC[k][k] = ae_c_neg(tau->xC[k]);
   // Rest of T is filled by zeros
      for (i = k + 1; i < blocksize; i++) {
         t->xyC[i][k] = ae_complex_from_i(0);
      }
   }
}

// Base case for real QR
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994.
//      Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
//      pseudocode, 2007-2010.
static void ortfac_rmatrixqrbasecase(RMatrix *a, ae_int_t m, ae_int_t n, RVector *work, RVector *t, RVector *tau) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t minmn;
   double tmp;
   minmn = imin2(m, n);
// Test the input arguments
   k = minmn;
   for (i = 0; i < k; i++) {
   // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
      ae_v_move(&t->xR[1], 1, &a->xyR[i][i], a->stride, m - i);
      generatereflection(t, m - i, &tmp);
      tau->xR[i] = tmp;
      ae_v_move(&a->xyR[i][i], a->stride, &t->xR[1], 1, m - i);
      t->xR[1] = 1.0;
      if (i < n) {
      // Apply H(i) to A(i:m-1,i+1:n-1) from the left
         applyreflectionfromtheleft(a, tau->xR[i], t, i, m - 1, i + 1, n - 1, work);
      }
   }
}

// Base case for complex QR
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994.
//      Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
//      pseudocode, 2007-2010.
static void ortfac_cmatrixqrbasecase(CMatrix *a, ae_int_t m, ae_int_t n, CVector *work, CVector *t, CVector *tau) {
   ae_int_t i;
   ae_int_t k;
   ae_int_t mmi;
   ae_int_t minmn;
   ae_complex tmp;
   minmn = imin2(m, n);
   if (minmn <= 0) {
      return;
   }
// Test the input arguments
   k = imin2(m, n);
   for (i = 0; i < k; i++) {
   // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
      mmi = m - i;
      ae_v_cmove(&t->xC[1], 1, &a->xyC[i][i], a->stride, "N", mmi);
      complexgeneratereflection(t, mmi, &tmp);
      tau->xC[i] = tmp;
      ae_v_cmove(&a->xyC[i][i], a->stride, &t->xC[1], 1, "N", m - i);
      t->xC[1] = ae_complex_from_i(1);
      if (i < n - 1) {
      // Apply H'(i) to A(i:m,i+1:n) from the left
         complexapplyreflectionfromtheleft(a, ae_c_conj(tau->xC[i]), t, i, m - 1, i + 1, n - 1, work);
      }
   }
}

// QR decomposition of a rectangular matrix of size MxN
//
// Inputs:
//     A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
//     M   -   number of rows in matrix A.
//     N   -   number of columns in matrix A.
//
// Outputs:
//     A   -   matrices Q and R in compact form (see below).
//     Tau -   array of scalar factors which are used to form
//             matrix Q. Array whose index ranges within [0.. Min(M-1,N-1)].
//
// Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
// MxM, R - upper triangular (or upper trapezoid) matrix of size M x N.
//
// The elements of matrix R are located on and above the main diagonal of
// matrix A. The elements which are located in Tau array and below the main
// diagonal of matrix A are used to form matrix Q as follows:
//
// Matrix Q is represented as a product of elementary reflections
//
// Q = H(0)*H(2)*...*H(k-1),
//
// where k = min(m,n), and each H(i) is in the form
//
// H(i) = 1 - tau * v * (v^T)
//
// where tau is a scalar stored in Tau[I]; v - real vector,
// so that v(0:i-1) = 0, v(i) = 1, v(i+1:m-1) stored in A(i+1:m-1,i).
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void rmatrixqr(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau);
void rmatrixqr(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t rowscount;
   ae_int_t i;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   NewVector(work, 0, DT_REAL);
   NewVector(t, 0, DT_REAL);
   NewVector(taubuf, 0, DT_REAL);
   NewMatrix(tmpa, 0, 0, DT_REAL);
   NewMatrix(tmpt, 0, 0, DT_REAL);
   NewMatrix(tmpr, 0, 0, DT_REAL);
   if (m <= 0 || n <= 0) {
      ae_frame_leave();
      return;
   }
   minmn = imin2(m, n);
   ts = matrixtilesizeb();
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(tau, minmn);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, m, ts);
   ae_matrix_set_length(&tmpt, ts, 2 * ts);
   ae_matrix_set_length(&tmpr, 2 * ts, n);
// Blocked code
   blockstart = 0;
   while (blockstart != minmn) {
   // Determine block size
      blocksize = minmn - blockstart;
      if (blocksize > ts) {
         blocksize = ts;
      }
      rowscount = m - blockstart;
   // QR decomposition of submatrix.
   // Matrix is copied to temporary storage to solve
   // some TLB issues arising from non-contiguous memory
   // access pattern.
      rmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, &tmpa, 0, 0);
      ortfac_rmatrixqrbasecase(&tmpa, rowscount, blocksize, &work, &t, &taubuf);
      rmatrixcopy(rowscount, blocksize, &tmpa, 0, 0, a, blockstart, blockstart);
      ae_v_move(&tau->xR[blockstart], 1, taubuf.xR, 1, blocksize);
   // Update the rest, choose between:
   // a) Level 2 algorithm (when the rest of the matrix is small enough)
   // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
   //    representation for products of Householder transformations',
   //    by R. Schreiber and C. Van Loan.
      if (blockstart + blocksize < n) {
         if (n - blockstart - blocksize >= 2 * ts || rowscount >= 4 * ts) {
         // Prepare block reflector
            ortfac_rmatrixblockreflector(&tmpa, &taubuf, true, rowscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q'.
         //
         // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
         // Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
            rmatrixgemm(blocksize, n - blockstart - blocksize, rowscount, 1.0, &tmpa, 0, 0, 1, a, blockstart, blockstart + blocksize, 0, 0.0, &tmpr, 0, 0);
            rmatrixgemm(blocksize, n - blockstart - blocksize, blocksize, 1.0, &tmpt, 0, 0, 1, &tmpr, 0, 0, 0, 0.0, &tmpr, blocksize, 0);
            rmatrixgemm(rowscount, n - blockstart - blocksize, blocksize, 1.0, &tmpa, 0, 0, 0, &tmpr, blocksize, 0, 0, 1.0, a, blockstart, blockstart + blocksize);
         } else {
         // Level 2 algorithm
            for (i = 0; i < blocksize; i++) {
               ae_v_move(&t.xR[1], 1, &tmpa.xyR[i][i], tmpa.stride, rowscount - i);
               t.xR[1] = 1.0;
               applyreflectionfromtheleft(a, taubuf.xR[i], &t, blockstart + i, m - 1, blockstart + blocksize, n - 1, &work);
            }
         }
      }
   // Advance
      blockstart += blocksize;
   }
   ae_frame_leave();
}

// QR decomposition of a rectangular complex matrix of size MxN
//
// Inputs:
//     A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
//     M   -   number of rows in matrix A.
//     N   -   number of columns in matrix A.
//
// Outputs:
//     A   -   matrices Q and R in compact form
//     Tau -   array of scalar factors which are used to form matrix Q. Array
//             whose indexes range within [0.. Min(M,N)-1]
//
// Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
// MxM, R - upper triangular (or upper trapezoid) matrix of size MxN.
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
// API: void cmatrixqr(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau);
void cmatrixqr(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t rowscount;
   ae_int_t i;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   NewVector(work, 0, DT_COMPLEX);
   NewVector(t, 0, DT_COMPLEX);
   NewVector(taubuf, 0, DT_COMPLEX);
   NewMatrix(tmpa, 0, 0, DT_COMPLEX);
   NewMatrix(tmpt, 0, 0, DT_COMPLEX);
   NewMatrix(tmpr, 0, 0, DT_COMPLEX);
   if (m <= 0 || n <= 0) {
      ae_frame_leave();
      return;
   }
   ts = matrixtilesizeb() / 2;
   minmn = imin2(m, n);
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(tau, minmn);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, m, ts);
   ae_matrix_set_length(&tmpt, ts, ts);
   ae_matrix_set_length(&tmpr, 2 * ts, n);
// Blocked code
   blockstart = 0;
   while (blockstart != minmn) {
   // Determine block size
      blocksize = minmn - blockstart;
      if (blocksize > ts) {
         blocksize = ts;
      }
      rowscount = m - blockstart;
   // QR decomposition of submatrix.
   // Matrix is copied to temporary storage to solve
   // some TLB issues arising from non-contiguous memory
   // access pattern.
      cmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, &tmpa, 0, 0);
      ortfac_cmatrixqrbasecase(&tmpa, rowscount, blocksize, &work, &t, &taubuf);
      cmatrixcopy(rowscount, blocksize, &tmpa, 0, 0, a, blockstart, blockstart);
      ae_v_cmove(&tau->xC[blockstart], 1, taubuf.xC, 1, "N", blocksize);
   // Update the rest, choose between:
   // a) Level 2 algorithm (when the rest of the matrix is small enough)
   // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
   //    representation for products of Householder transformations',
   //    by R. Schreiber and C. Van Loan.
      if (blockstart + blocksize < n) {
         if (n - blockstart - blocksize >= 2 * ts) {
         // Prepare block reflector
            ortfac_cmatrixblockreflector(&tmpa, &taubuf, true, rowscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q'.
         //
         // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
         // Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
            cmatrixgemm(blocksize, n - blockstart - blocksize, rowscount, ae_complex_from_d(1.0), &tmpa, 0, 0, 2, a, blockstart, blockstart + blocksize, 0, ae_complex_from_d(0.0), &tmpr, 0, 0);
            cmatrixgemm(blocksize, n - blockstart - blocksize, blocksize, ae_complex_from_d(1.0), &tmpt, 0, 0, 2, &tmpr, 0, 0, 0, ae_complex_from_d(0.0), &tmpr, blocksize, 0);
            cmatrixgemm(rowscount, n - blockstart - blocksize, blocksize, ae_complex_from_d(1.0), &tmpa, 0, 0, 0, &tmpr, blocksize, 0, 0, ae_complex_from_d(1.0), a, blockstart, blockstart + blocksize);
         } else {
         // Level 2 algorithm
            for (i = 0; i < blocksize; i++) {
               ae_v_cmove(&t.xC[1], 1, &tmpa.xyC[i][i], tmpa.stride, "N", rowscount - i);
               t.xC[1] = ae_complex_from_i(1);
               complexapplyreflectionfromtheleft(a, ae_c_conj(taubuf.xC[i]), &t, blockstart + i, m - 1, blockstart + blocksize, n - 1, &work);
            }
         }
      }
   // Advance
      blockstart += blocksize;
   }
   ae_frame_leave();
}

// Base case for real LQ
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994.
//      Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
//      pseudocode, 2007-2010.
void rmatrixlqbasecase(RMatrix *a, ae_int_t m, ae_int_t n, RVector *work, RVector *t, RVector *tau) {
   ae_int_t i;
   ae_int_t k;
   double tmp;
   k = imin2(m, n);
   for (i = 0; i < k; i++) {
   // Generate elementary reflector H(i) to annihilate A(i,i+1:n-1)
      ae_v_move(&t->xR[1], 1, &a->xyR[i][i], 1, n - i);
      generatereflection(t, n - i, &tmp);
      tau->xR[i] = tmp;
      ae_v_move(&a->xyR[i][i], 1, &t->xR[1], 1, n - i);
      t->xR[1] = 1.0;
      if (i < n) {
      // Apply H(i) to A(i+1:m,i:n) from the right
         applyreflectionfromtheright(a, tau->xR[i], t, i + 1, m - 1, i, n - 1, work);
      }
   }
}

// Base case for complex LQ
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994.
//      Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
//      pseudocode, 2007-2010.
static void ortfac_cmatrixlqbasecase(CMatrix *a, ae_int_t m, ae_int_t n, CVector *work, CVector *t, CVector *tau) {
   ae_int_t i;
   ae_int_t minmn;
   ae_complex tmp;
   minmn = imin2(m, n);
   if (minmn <= 0) {
      return;
   }
// Test the input arguments
   for (i = 0; i < minmn; i++) {
   // Generate elementary reflector H(i)
   //
   // NOTE: ComplexGenerateReflection() generates left reflector,
   // i.e. H which reduces x by applyiong from the left, but we
   // need RIGHT reflector. So we replace H=E-tau*v*v' by H^H,
   // which changes v to conj(v).
      ae_v_cmove(&t->xC[1], 1, &a->xyC[i][i], 1, "Conj", n - i);
      complexgeneratereflection(t, n - i, &tmp);
      tau->xC[i] = tmp;
      ae_v_cmove(&a->xyC[i][i], 1, &t->xC[1], 1, "Conj", n - i);
      t->xC[1] = ae_complex_from_i(1);
      if (i < m - 1) {
      // Apply H'(i)
         complexapplyreflectionfromtheright(a, tau->xC[i], t, i + 1, m - 1, i, n - 1, work);
      }
   }
}

// LQ decomposition of a rectangular matrix of size MxN
//
// Inputs:
//     A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
//     M   -   number of rows in matrix A.
//     N   -   number of columns in matrix A.
//
// Outputs:
//     A   -   matrices L and Q in compact form (see below)
//     Tau -   array of scalar factors which are used to form
//             matrix Q. Array whose index ranges within [0..Min(M,N)-1].
//
// Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
// MxM, L - lower triangular (or lower trapezoid) matrix of size M x N.
//
// The elements of matrix L are located on and below  the  main  diagonal  of
// matrix A. The elements which are located in Tau array and above  the  main
// diagonal of matrix A are used to form matrix Q as follows:
//
// Matrix Q is represented as a product of elementary reflections
//
// Q = H(k-1)*H(k-2)*...*H(1)*H(0),
//
// where k = min(m,n), and each H(i) is of the form
//
// H(i) = 1 - tau * v * (v^T)
//
// where tau is a scalar stored in Tau[I]; v - real vector, so that v(0:i-1)=0,
// v(i) = 1, v(i+1:n-1) stored in A(i,i+1:n-1).
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void rmatrixlq(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau);
void rmatrixlq(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t columnscount;
   ae_int_t i;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   NewVector(work, 0, DT_REAL);
   NewVector(t, 0, DT_REAL);
   NewVector(taubuf, 0, DT_REAL);
   NewMatrix(tmpa, 0, 0, DT_REAL);
   NewMatrix(tmpt, 0, 0, DT_REAL);
   NewMatrix(tmpr, 0, 0, DT_REAL);
   if (m <= 0 || n <= 0) {
      ae_frame_leave();
      return;
   }
   minmn = imin2(m, n);
   ts = matrixtilesizeb();
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(tau, minmn);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, ts, n);
   ae_matrix_set_length(&tmpt, ts, 2 * ts);
   ae_matrix_set_length(&tmpr, m, 2 * ts);
// Blocked code
   blockstart = 0;
   while (blockstart != minmn) {
   // Determine block size
      blocksize = minmn - blockstart;
      if (blocksize > ts) {
         blocksize = ts;
      }
      columnscount = n - blockstart;
   // LQ decomposition of submatrix.
   // Matrix is copied to temporary storage to solve
   // some TLB issues arising from non-contiguous memory
   // access pattern.
      rmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, &tmpa, 0, 0);
      rmatrixlqbasecase(&tmpa, blocksize, columnscount, &work, &t, &taubuf);
      rmatrixcopy(blocksize, columnscount, &tmpa, 0, 0, a, blockstart, blockstart);
      ae_v_move(&tau->xR[blockstart], 1, taubuf.xR, 1, blocksize);
   // Update the rest, choose between:
   // a) Level 2 algorithm (when the rest of the matrix is small enough)
   // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
   //    representation for products of Householder transformations',
   //    by R. Schreiber and C. Van Loan.
      if (blockstart + blocksize < m) {
         if (m - blockstart - blocksize >= 2 * ts) {
         // Prepare block reflector
            ortfac_rmatrixblockreflector(&tmpa, &taubuf, false, columnscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q.
         //
         // Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
            rmatrixgemm(m - blockstart - blocksize, blocksize, columnscount, 1.0, a, blockstart + blocksize, blockstart, 0, &tmpa, 0, 0, 1, 0.0, &tmpr, 0, 0);
            rmatrixgemm(m - blockstart - blocksize, blocksize, blocksize, 1.0, &tmpr, 0, 0, 0, &tmpt, 0, 0, 0, 0.0, &tmpr, 0, blocksize);
            rmatrixgemm(m - blockstart - blocksize, columnscount, blocksize, 1.0, &tmpr, 0, blocksize, 0, &tmpa, 0, 0, 0, 1.0, a, blockstart + blocksize, blockstart);
         } else {
         // Level 2 algorithm
            for (i = 0; i < blocksize; i++) {
               ae_v_move(&t.xR[1], 1, &tmpa.xyR[i][i], 1, columnscount - i);
               t.xR[1] = 1.0;
               applyreflectionfromtheright(a, taubuf.xR[i], &t, blockstart + blocksize, m - 1, blockstart + i, n - 1, &work);
            }
         }
      }
   // Advance
      blockstart += blocksize;
   }
   ae_frame_leave();
}

// LQ decomposition of a rectangular complex matrix of size MxN
//
// Inputs:
//     A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
//     M   -   number of rows in matrix A.
//     N   -   number of columns in matrix A.
//
// Outputs:
//     A   -   matrices Q and L in compact form
//     Tau -   array of scalar factors which are used to form matrix Q. Array
//             whose indexes range within [0.. Min(M,N)-1]
//
// Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
// MxM, L - lower triangular (or lower trapezoid) matrix of size MxN.
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
// API: void cmatrixlq(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau);
void cmatrixlq(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t columnscount;
   ae_int_t i;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   NewVector(work, 0, DT_COMPLEX);
   NewVector(t, 0, DT_COMPLEX);
   NewVector(taubuf, 0, DT_COMPLEX);
   NewMatrix(tmpa, 0, 0, DT_COMPLEX);
   NewMatrix(tmpt, 0, 0, DT_COMPLEX);
   NewMatrix(tmpr, 0, 0, DT_COMPLEX);
   if (m <= 0 || n <= 0) {
      ae_frame_leave();
      return;
   }
   ts = matrixtilesizeb() / 2;
   minmn = imin2(m, n);
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(tau, minmn);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, ts, n);
   ae_matrix_set_length(&tmpt, ts, ts);
   ae_matrix_set_length(&tmpr, m, 2 * ts);
// Blocked code
   blockstart = 0;
   while (blockstart != minmn) {
   // Determine block size
      blocksize = minmn - blockstart;
      if (blocksize > ts) {
         blocksize = ts;
      }
      columnscount = n - blockstart;
   // LQ decomposition of submatrix.
   // Matrix is copied to temporary storage to solve
   // some TLB issues arising from non-contiguous memory
   // access pattern.
      cmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, &tmpa, 0, 0);
      ortfac_cmatrixlqbasecase(&tmpa, blocksize, columnscount, &work, &t, &taubuf);
      cmatrixcopy(blocksize, columnscount, &tmpa, 0, 0, a, blockstart, blockstart);
      ae_v_cmove(&tau->xC[blockstart], 1, taubuf.xC, 1, "N", blocksize);
   // Update the rest, choose between:
   // a) Level 2 algorithm (when the rest of the matrix is small enough)
   // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
   //    representation for products of Householder transformations',
   //    by R. Schreiber and C. Van Loan.
      if (blockstart + blocksize < m) {
         if (m - blockstart - blocksize >= 2 * ts) {
         // Prepare block reflector
            ortfac_cmatrixblockreflector(&tmpa, &taubuf, false, columnscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q.
         //
         // Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
            cmatrixgemm(m - blockstart - blocksize, blocksize, columnscount, ae_complex_from_d(1.0), a, blockstart + blocksize, blockstart, 0, &tmpa, 0, 0, 2, ae_complex_from_d(0.0), &tmpr, 0, 0);
            cmatrixgemm(m - blockstart - blocksize, blocksize, blocksize, ae_complex_from_d(1.0), &tmpr, 0, 0, 0, &tmpt, 0, 0, 0, ae_complex_from_d(0.0), &tmpr, 0, blocksize);
            cmatrixgemm(m - blockstart - blocksize, columnscount, blocksize, ae_complex_from_d(1.0), &tmpr, 0, blocksize, 0, &tmpa, 0, 0, 0, ae_complex_from_d(1.0), a, blockstart + blocksize, blockstart);
         } else {
         // Level 2 algorithm
            for (i = 0; i < blocksize; i++) {
               ae_v_cmove(&t.xC[1], 1, &tmpa.xyC[i][i], 1, "Conj", columnscount - i);
               t.xC[1] = ae_complex_from_i(1);
               complexapplyreflectionfromtheright(a, taubuf.xC[i], &t, blockstart + blocksize, m - 1, blockstart + i, n - 1, &work);
            }
         }
      }
   // Advance
      blockstart += blocksize;
   }
   ae_frame_leave();
}

// Partial unpacking of matrix Q from the QR decomposition of a matrix A
//
// Inputs:
//     A       -   matrices Q and R in compact form.
//                 Output of RMatrixQR subroutine.
//     M       -   number of rows in given matrix A. M >= 0.
//     N       -   number of columns in given matrix A. N >= 0.
//     Tau     -   scalar factors which are used to form Q.
//                 Output of the RMatrixQR subroutine.
//     QColumns -  required number of columns of matrix Q. M >= QColumns >= 0.
//
// Outputs:
//     Q       -   first QColumns columns of matrix Q.
//                 Array whose indexes range within [0..M-1, 0..QColumns-1].
//                 If QColumns=0, the array remains unchanged.
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void rmatrixqrunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qcolumns, real_2d_array &q);
void rmatrixqrunpackq(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau, ae_int_t qcolumns, RMatrix *q) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t refcnt;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t rowscount;
   ae_int_t i;
   ae_int_t j;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(work, 0, DT_REAL);
   NewVector(t, 0, DT_REAL);
   NewVector(taubuf, 0, DT_REAL);
   NewMatrix(tmpa, 0, 0, DT_REAL);
   NewMatrix(tmpt, 0, 0, DT_REAL);
   NewMatrix(tmpr, 0, 0, DT_REAL);
   ae_assert(qcolumns <= m, "UnpackQFromQR: QColumns>M!");
   if (m <= 0 || n <= 0 || qcolumns <= 0) {
      ae_frame_leave();
      return;
   }
// init
   ts = matrixtilesizeb();
   minmn = imin2(m, n);
   refcnt = imin2(minmn, qcolumns);
   ae_matrix_set_length(q, m, qcolumns);
   for (i = 0; i < m; i++) {
      for (j = 0; j < qcolumns; j++) {
         if (i == j) {
            q->xyR[i][j] = 1.0;
         } else {
            q->xyR[i][j] = 0.0;
         }
      }
   }
   ae_vector_set_length(&work, imax2(m, qcolumns) + 1);
   ae_vector_set_length(&t, imax2(m, qcolumns) + 1);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, m, ts);
   ae_matrix_set_length(&tmpt, ts, 2 * ts);
   ae_matrix_set_length(&tmpr, 2 * ts, qcolumns);
// Blocked code
   blockstart = ts * (refcnt / ts);
   blocksize = refcnt - blockstart;
   while (blockstart >= 0) {
      rowscount = m - blockstart;
      if (blocksize > 0) {
      // Copy current block
         rmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, &tmpa, 0, 0);
         ae_v_move(taubuf.xR, 1, &tau->xR[blockstart], 1, blocksize);
      // Update, choose between:
      // a) Level 2 algorithm (when the rest of the matrix is small enough)
      // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //    representation for products of Householder transformations',
      //    by R. Schreiber and C. Van Loan.
         if (qcolumns >= 2 * ts) {
         // Prepare block reflector
            ortfac_rmatrixblockreflector(&tmpa, &taubuf, true, rowscount, blocksize, &tmpt, &work);
         // Multiply matrix by Q.
         //
         // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            rmatrixgemm(blocksize, qcolumns, rowscount, 1.0, &tmpa, 0, 0, 1, q, blockstart, 0, 0, 0.0, &tmpr, 0, 0);
            rmatrixgemm(blocksize, qcolumns, blocksize, 1.0, &tmpt, 0, 0, 0, &tmpr, 0, 0, 0, 0.0, &tmpr, blocksize, 0);
            rmatrixgemm(rowscount, qcolumns, blocksize, 1.0, &tmpa, 0, 0, 0, &tmpr, blocksize, 0, 0, 1.0, q, blockstart, 0);
         } else {
         // Level 2 algorithm
            for (i = blocksize - 1; i >= 0; i--) {
               ae_v_move(&t.xR[1], 1, &tmpa.xyR[i][i], tmpa.stride, rowscount - i);
               t.xR[1] = 1.0;
               applyreflectionfromtheleft(q, taubuf.xR[i], &t, blockstart + i, m - 1, 0, qcolumns - 1, &work);
            }
         }
      }
   // Advance
      blockstart -= ts;
      blocksize = ts;
   }
   ae_frame_leave();
}

// Partial unpacking of matrix Q from QR decomposition of a complex matrix A.
//
// Inputs:
//     A           -   matrices Q and R in compact form.
//                     Output of CMatrixQR subroutine .
//     M           -   number of rows in matrix A. M >= 0.
//     N           -   number of columns in matrix A. N >= 0.
//     Tau         -   scalar factors which are used to form Q.
//                     Output of CMatrixQR subroutine .
//     QColumns    -   required number of columns in matrix Q. M >= QColumns >= 0.
//
// Outputs:
//     Q           -   first QColumns columns of matrix Q.
//                     Array whose index ranges within [0..M-1, 0..QColumns-1].
//                     If QColumns=0, array isn't changed.
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void cmatrixqrunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qcolumns, complex_2d_array &q);
void cmatrixqrunpackq(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau, ae_int_t qcolumns, CMatrix *q) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t refcnt;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t rowscount;
   ae_int_t i;
   ae_int_t j;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(work, 0, DT_COMPLEX);
   NewVector(t, 0, DT_COMPLEX);
   NewVector(taubuf, 0, DT_COMPLEX);
   NewMatrix(tmpa, 0, 0, DT_COMPLEX);
   NewMatrix(tmpt, 0, 0, DT_COMPLEX);
   NewMatrix(tmpr, 0, 0, DT_COMPLEX);
   ae_assert(qcolumns <= m, "UnpackQFromQR: QColumns>M!");
   if (m <= 0 || n <= 0) {
      ae_frame_leave();
      return;
   }
// init
   ts = matrixtilesizeb() / 2;
   minmn = imin2(m, n);
   refcnt = imin2(minmn, qcolumns);
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, m, ts);
   ae_matrix_set_length(&tmpt, ts, ts);
   ae_matrix_set_length(&tmpr, 2 * ts, qcolumns);
   ae_matrix_set_length(q, m, qcolumns);
   for (i = 0; i < m; i++) {
      for (j = 0; j < qcolumns; j++) {
         if (i == j) {
            q->xyC[i][j] = ae_complex_from_i(1);
         } else {
            q->xyC[i][j] = ae_complex_from_i(0);
         }
      }
   }
// Blocked code
   blockstart = ts * (refcnt / ts);
   blocksize = refcnt - blockstart;
   while (blockstart >= 0) {
      rowscount = m - blockstart;
      if (blocksize > 0) {
      // QR decomposition of submatrix.
      // Matrix is copied to temporary storage to solve
      // some TLB issues arising from non-contiguous memory
      // access pattern.
         cmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, &tmpa, 0, 0);
         ae_v_cmove(taubuf.xC, 1, &tau->xC[blockstart], 1, "N", blocksize);
      // Update matrix, choose between:
      // a) Level 2 algorithm (when the rest of the matrix is small enough)
      // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //    representation for products of Householder transformations',
      //    by R. Schreiber and C. Van Loan.
         if (qcolumns >= 2 * ts) {
         // Prepare block reflector
            ortfac_cmatrixblockreflector(&tmpa, &taubuf, true, rowscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q.
         //
         // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            cmatrixgemm(blocksize, qcolumns, rowscount, ae_complex_from_d(1.0), &tmpa, 0, 0, 2, q, blockstart, 0, 0, ae_complex_from_d(0.0), &tmpr, 0, 0);
            cmatrixgemm(blocksize, qcolumns, blocksize, ae_complex_from_d(1.0), &tmpt, 0, 0, 0, &tmpr, 0, 0, 0, ae_complex_from_d(0.0), &tmpr, blocksize, 0);
            cmatrixgemm(rowscount, qcolumns, blocksize, ae_complex_from_d(1.0), &tmpa, 0, 0, 0, &tmpr, blocksize, 0, 0, ae_complex_from_d(1.0), q, blockstart, 0);
         } else {
         // Level 2 algorithm
            for (i = blocksize - 1; i >= 0; i--) {
               ae_v_cmove(&t.xC[1], 1, &tmpa.xyC[i][i], tmpa.stride, "N", rowscount - i);
               t.xC[1] = ae_complex_from_i(1);
               complexapplyreflectionfromtheleft(q, taubuf.xC[i], &t, blockstart + i, m - 1, 0, qcolumns - 1, &work);
            }
         }
      }
   // Advance
      blockstart -= ts;
      blocksize = ts;
   }
   ae_frame_leave();
}

// Unpacking of matrix R from the QR decomposition of a matrix A
//
// Inputs:
//     A       -   matrices Q and R in compact form.
//                 Output of RMatrixQR subroutine.
//     M       -   number of rows in given matrix A. M >= 0.
//     N       -   number of columns in given matrix A. N >= 0.
//
// Outputs:
//     R       -   matrix R, array[0..M-1, 0..N-1].
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void rmatrixqrunpackr(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &r);
void rmatrixqrunpackr(RMatrix *a, ae_int_t m, ae_int_t n, RMatrix *r) {
   ae_int_t i;
   ae_int_t k;
   SetMatrix(r);
   if (m <= 0 || n <= 0) {
      return;
   }
   k = imin2(m, n);
   ae_matrix_set_length(r, m, n);
   for (i = 0; i < n; i++) {
      r->xyR[0][i] = 0.0;
   }
   for (i = 1; i < m; i++) {
      ae_v_move(r->xyR[i], 1, r->xyR[0], 1, n);
   }
   for (i = 0; i < k; i++) {
      ae_v_move(&r->xyR[i][i], 1, &a->xyR[i][i], 1, n - i);
   }
}

// Unpacking of matrix R from the QR decomposition of a matrix A
//
// Inputs:
//     A       -   matrices Q and R in compact form.
//                 Output of CMatrixQR subroutine.
//     M       -   number of rows in given matrix A. M >= 0.
//     N       -   number of columns in given matrix A. N >= 0.
//
// Outputs:
//     R       -   matrix R, array[0..M-1, 0..N-1].
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void cmatrixqrunpackr(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &r);
void cmatrixqrunpackr(CMatrix *a, ae_int_t m, ae_int_t n, CMatrix *r) {
   ae_int_t i;
   ae_int_t k;
   SetMatrix(r);
   if (m <= 0 || n <= 0) {
      return;
   }
   k = imin2(m, n);
   ae_matrix_set_length(r, m, n);
   for (i = 0; i < n; i++) {
      r->xyC[0][i] = ae_complex_from_i(0);
   }
   for (i = 1; i < m; i++) {
      ae_v_cmove(r->xyC[i], 1, r->xyC[0], 1, "N", n);
   }
   for (i = 0; i < k; i++) {
      ae_v_cmove(&r->xyC[i][i], 1, &a->xyC[i][i], 1, "N", n - i);
   }
}

// Partial unpacking of matrix Q from the LQ decomposition of a matrix A
//
// Inputs:
//     A       -   matrices L and Q in compact form.
//                 Output of RMatrixLQ subroutine.
//     M       -   number of rows in given matrix A. M >= 0.
//     N       -   number of columns in given matrix A. N >= 0.
//     Tau     -   scalar factors which are used to form Q.
//                 Output of the RMatrixLQ subroutine.
//     QRows   -   required number of rows in matrix Q. N >= QRows >= 0.
//
// Outputs:
//     Q       -   first QRows rows of matrix Q. Array whose indexes range
//                 within [0..QRows-1, 0..N-1]. If QRows=0, the array remains
//                 unchanged.
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void rmatrixlqunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qrows, real_2d_array &q);
void rmatrixlqunpackq(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tau, ae_int_t qrows, RMatrix *q) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t refcnt;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t columnscount;
   ae_int_t i;
   ae_int_t j;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(work, 0, DT_REAL);
   NewVector(t, 0, DT_REAL);
   NewVector(taubuf, 0, DT_REAL);
   NewMatrix(tmpa, 0, 0, DT_REAL);
   NewMatrix(tmpt, 0, 0, DT_REAL);
   NewMatrix(tmpr, 0, 0, DT_REAL);
   ae_assert(qrows <= n, "RMatrixLQUnpackQ: QRows>N!");
   if (m <= 0 || n <= 0 || qrows <= 0) {
      ae_frame_leave();
      return;
   }
// init
   ts = matrixtilesizeb();
   minmn = imin2(m, n);
   refcnt = imin2(minmn, qrows);
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, ts, n);
   ae_matrix_set_length(&tmpt, ts, 2 * ts);
   ae_matrix_set_length(&tmpr, qrows, 2 * ts);
   ae_matrix_set_length(q, qrows, n);
   for (i = 0; i < qrows; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            q->xyR[i][j] = 1.0;
         } else {
            q->xyR[i][j] = 0.0;
         }
      }
   }
// Blocked code
   blockstart = ts * (refcnt / ts);
   blocksize = refcnt - blockstart;
   while (blockstart >= 0) {
      columnscount = n - blockstart;
      if (blocksize > 0) {
      // Copy submatrix
         rmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, &tmpa, 0, 0);
         ae_v_move(taubuf.xR, 1, &tau->xR[blockstart], 1, blocksize);
      // Update matrix, choose between:
      // a) Level 2 algorithm (when the rest of the matrix is small enough)
      // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //    representation for products of Householder transformations',
      //    by R. Schreiber and C. Van Loan.
         if (qrows >= 2 * ts) {
         // Prepare block reflector
            ortfac_rmatrixblockreflector(&tmpa, &taubuf, false, columnscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q'.
         //
         // Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
            rmatrixgemm(qrows, blocksize, columnscount, 1.0, q, 0, blockstart, 0, &tmpa, 0, 0, 1, 0.0, &tmpr, 0, 0);
            rmatrixgemm(qrows, blocksize, blocksize, 1.0, &tmpr, 0, 0, 0, &tmpt, 0, 0, 1, 0.0, &tmpr, 0, blocksize);
            rmatrixgemm(qrows, columnscount, blocksize, 1.0, &tmpr, 0, blocksize, 0, &tmpa, 0, 0, 0, 1.0, q, 0, blockstart);
         } else {
         // Level 2 algorithm
            for (i = blocksize - 1; i >= 0; i--) {
               ae_v_move(&t.xR[1], 1, &tmpa.xyR[i][i], 1, columnscount - i);
               t.xR[1] = 1.0;
               applyreflectionfromtheright(q, taubuf.xR[i], &t, 0, qrows - 1, blockstart + i, n - 1, &work);
            }
         }
      }
   // Advance
      blockstart -= ts;
      blocksize = ts;
   }
   ae_frame_leave();
}

// Partial unpacking of matrix Q from LQ decomposition of a complex matrix A.
//
// Inputs:
//     A           -   matrices Q and R in compact form.
//                     Output of CMatrixLQ subroutine .
//     M           -   number of rows in matrix A. M >= 0.
//     N           -   number of columns in matrix A. N >= 0.
//     Tau         -   scalar factors which are used to form Q.
//                     Output of CMatrixLQ subroutine .
//     QRows       -   required number of rows in matrix Q. N >= QColumns >= 0.
//
// Outputs:
//     Q           -   first QRows rows of matrix Q.
//                     Array whose index ranges within [0..QRows-1, 0..N-1].
//                     If QRows=0, array isn't changed.
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void cmatrixlqunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qrows, complex_2d_array &q);
void cmatrixlqunpackq(CMatrix *a, ae_int_t m, ae_int_t n, CVector *tau, ae_int_t qrows, CMatrix *q) {
   ae_frame _frame_block;
   ae_int_t minmn;
   ae_int_t refcnt;
   ae_int_t blockstart;
   ae_int_t blocksize;
   ae_int_t columnscount;
   ae_int_t i;
   ae_int_t j;
   ae_int_t ts;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(work, 0, DT_COMPLEX);
   NewVector(t, 0, DT_COMPLEX);
   NewVector(taubuf, 0, DT_COMPLEX);
   NewMatrix(tmpa, 0, 0, DT_COMPLEX);
   NewMatrix(tmpt, 0, 0, DT_COMPLEX);
   NewMatrix(tmpr, 0, 0, DT_COMPLEX);
   if (m <= 0 || n <= 0) {
      ae_frame_leave();
      return;
   }
// Init
   ts = matrixtilesizeb() / 2;
   minmn = imin2(m, n);
   refcnt = imin2(minmn, qrows);
   ae_vector_set_length(&work, imax2(m, n) + 1);
   ae_vector_set_length(&t, imax2(m, n) + 1);
   ae_vector_set_length(&taubuf, minmn);
   ae_matrix_set_length(&tmpa, ts, n);
   ae_matrix_set_length(&tmpt, ts, ts);
   ae_matrix_set_length(&tmpr, qrows, 2 * ts);
   ae_matrix_set_length(q, qrows, n);
   for (i = 0; i < qrows; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            q->xyC[i][j] = ae_complex_from_i(1);
         } else {
            q->xyC[i][j] = ae_complex_from_i(0);
         }
      }
   }
// Blocked code
   blockstart = ts * (refcnt / ts);
   blocksize = refcnt - blockstart;
   while (blockstart >= 0) {
      columnscount = n - blockstart;
      if (blocksize > 0) {
      // LQ decomposition of submatrix.
      // Matrix is copied to temporary storage to solve
      // some TLB issues arising from non-contiguous memory
      // access pattern.
         cmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, &tmpa, 0, 0);
         ae_v_cmove(taubuf.xC, 1, &tau->xC[blockstart], 1, "N", blocksize);
      // Update matrix, choose between:
      // a) Level 2 algorithm (when the rest of the matrix is small enough)
      // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //    representation for products of Householder transformations',
      //    by R. Schreiber and C. Van Loan.
         if (qrows >= 2 * ts) {
         // Prepare block reflector
            ortfac_cmatrixblockreflector(&tmpa, &taubuf, false, columnscount, blocksize, &tmpt, &work);
         // Multiply the rest of A by Q'.
         //
         // Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
            cmatrixgemm(qrows, blocksize, columnscount, ae_complex_from_d(1.0), q, 0, blockstart, 0, &tmpa, 0, 0, 2, ae_complex_from_d(0.0), &tmpr, 0, 0);
            cmatrixgemm(qrows, blocksize, blocksize, ae_complex_from_d(1.0), &tmpr, 0, 0, 0, &tmpt, 0, 0, 2, ae_complex_from_d(0.0), &tmpr, 0, blocksize);
            cmatrixgemm(qrows, columnscount, blocksize, ae_complex_from_d(1.0), &tmpr, 0, blocksize, 0, &tmpa, 0, 0, 0, ae_complex_from_d(1.0), q, 0, blockstart);
         } else {
         // Level 2 algorithm
            for (i = blocksize - 1; i >= 0; i--) {
               ae_v_cmove(&t.xC[1], 1, &tmpa.xyC[i][i], 1, "Conj", columnscount - i);
               t.xC[1] = ae_complex_from_i(1);
               complexapplyreflectionfromtheright(q, ae_c_conj(taubuf.xC[i]), &t, 0, qrows - 1, blockstart + i, n - 1, &work);
            }
         }
      }
   // Advance
      blockstart -= ts;
      blocksize = ts;
   }
   ae_frame_leave();
}

// Unpacking of matrix L from the LQ decomposition of a matrix A
//
// Inputs:
//     A       -   matrices Q and L in compact form.
//                 Output of RMatrixLQ subroutine.
//     M       -   number of rows in given matrix A. M >= 0.
//     N       -   number of columns in given matrix A. N >= 0.
//
// Outputs:
//     L       -   matrix L, array[0..M-1, 0..N-1].
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void rmatrixlqunpackl(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &l);
void rmatrixlqunpackl(RMatrix *a, ae_int_t m, ae_int_t n, RMatrix *l) {
   ae_int_t i;
   ae_int_t k;
   SetMatrix(l);
   if (m <= 0 || n <= 0) {
      return;
   }
   ae_matrix_set_length(l, m, n);
   for (i = 0; i < n; i++) {
      l->xyR[0][i] = 0.0;
   }
   for (i = 1; i < m; i++) {
      ae_v_move(l->xyR[i], 1, l->xyR[0], 1, n);
   }
   for (i = 0; i < m; i++) {
      k = imin2(i, n - 1);
      ae_v_move(l->xyR[i], 1, a->xyR[i], 1, k + 1);
   }
}

// Unpacking of matrix L from the LQ decomposition of a matrix A
//
// Inputs:
//     A       -   matrices Q and L in compact form.
//                 Output of CMatrixLQ subroutine.
//     M       -   number of rows in given matrix A. M >= 0.
//     N       -   number of columns in given matrix A. N >= 0.
//
// Outputs:
//     L       -   matrix L, array[0..M-1, 0..N-1].
//
// ALGLIB Routine: Copyright 17.02.2010 by Sergey Bochkanov
// API: void cmatrixlqunpackl(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &l);
void cmatrixlqunpackl(CMatrix *a, ae_int_t m, ae_int_t n, CMatrix *l) {
   ae_int_t i;
   ae_int_t k;
   SetMatrix(l);
   if (m <= 0 || n <= 0) {
      return;
   }
   ae_matrix_set_length(l, m, n);
   for (i = 0; i < n; i++) {
      l->xyC[0][i] = ae_complex_from_i(0);
   }
   for (i = 1; i < m; i++) {
      ae_v_cmove(l->xyC[i], 1, l->xyC[0], 1, "N", n);
   }
   for (i = 0; i < m; i++) {
      k = imin2(i, n - 1);
      ae_v_cmove(l->xyC[i], 1, a->xyC[i], 1, "N", k + 1);
   }
}

// Reduction of a rectangular matrix to  bidiagonal form
//
// The algorithm reduces the rectangular matrix A to  bidiagonal form by
// orthogonal transformations P and Q: A = Q*B*(P^T).
//
// Inputs:
//     A       -   source matrix. array[0..M-1, 0..N-1]
//     M       -   number of rows in matrix A.
//     N       -   number of columns in matrix A.
//
// Outputs:
//     A       -   matrices Q, B, P in compact form (see below).
//     TauQ    -   scalar factors which are used to form matrix Q.
//     TauP    -   scalar factors which are used to form matrix P.
//
// The main diagonal and one of the  secondary  diagonals  of  matrix  A  are
// replaced with bidiagonal  matrix  B.  Other  elements  contain  elementary
// reflections which form MxM matrix Q and NxN matrix P, respectively.
//
// If M >= N, B is the upper  bidiagonal  MxN  matrix  and  is  stored  in  the
// corresponding  elements  of  matrix  A.  Matrix  Q  is  represented  as  a
// product   of   elementary   reflections   Q = H(0)*H(1)*...*H(n-1),  where
// H(i) = 1-tau*v*v'. Here tau is a scalar which is stored  in  TauQ[i],  and
// vector v has the following  structure:  v(0:i-1)=0, v(i)=1, v(i+1:m-1)  is
// stored   in   elements   A(i+1:m-1,i).   Matrix   P  is  as  follows:  P =
// G(0)*G(1)*...*G(n-2), where G(i) = 1 - tau*u*u'. Tau is stored in TauP[i],
// u(0:i)=0, u(i+1)=1, u(i+2:n-1) is stored in elements A(i,i+2:n-1).
//
// If M < N, B is the  lower  bidiagonal  MxN  matrix  and  is  stored  in  the
// corresponding   elements  of  matrix  A.  Q = H(0)*H(1)*...*H(m-2),  where
// H(i) = 1 - tau*v*v', tau is stored in TauQ, v(0:i)=0, v(i+1)=1, v(i+2:m-1)
// is    stored    in   elements   A(i+2:m-1,i).    P = G(0)*G(1)*...*G(m-1),
// G(i) = 1-tau*u*u', tau is stored in  TauP,  u(0:i-1)=0, u(i)=1, u(i+1:n-1)
// is stored in A(i,i+1:n-1).
//
// EXAMPLE:
//
// m=6, n=5 (m > n):               m=5, n=6 (m < n):
//
// (  d   e   u1  u1  u1 )         (  d   u1  u1  u1  u1  u1 )
// (  v1  d   e   u2  u2 )         (  e   d   u2  u2  u2  u2 )
// (  v1  v2  d   e   u3 )         (  v1  e   d   u3  u3  u3 )
// (  v1  v2  v3  d   e  )         (  v1  v2  e   d   u4  u4 )
// (  v1  v2  v3  v4  d  )         (  v1  v2  v3  e   d   u5 )
// (  v1  v2  v3  v4  v5 )
//
// Here vi and ui are vectors which form H(i) and G(i), and d and e -
// are the diagonal and off-diagonal elements of matrix B.
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994.
//      Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
//      pseudocode, 2007-2010.
// API: void rmatrixbd(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tauq, real_1d_array &taup);
void rmatrixbd(RMatrix *a, ae_int_t m, ae_int_t n, RVector *tauq, RVector *taup) {
   ae_frame _frame_block;
   ae_int_t maxmn;
   ae_int_t i;
   double ltau;
   ae_frame_make(&_frame_block);
   SetVector(tauq);
   SetVector(taup);
   NewVector(work, 0, DT_REAL);
   NewVector(t, 0, DT_REAL);
// Prepare
   if (n <= 0 || m <= 0) {
      ae_frame_leave();
      return;
   }
   maxmn = imax2(m, n);
   ae_vector_set_length(&work, maxmn + 1);
   ae_vector_set_length(&t, maxmn + 1);
   if (m >= n) {
      ae_vector_set_length(tauq, n);
      ae_vector_set_length(taup, n);
      for (i = 0; i < n; i++) {
         tauq->xR[i] = 0.0;
         taup->xR[i] = 0.0;
      }
   } else {
      ae_vector_set_length(tauq, m);
      ae_vector_set_length(taup, m);
      for (i = 0; i < m; i++) {
         tauq->xR[i] = 0.0;
         taup->xR[i] = 0.0;
      }
   }
// Try to use MKL code
//
// NOTE: buffers Work[] and T[] are used for temporary storage of diagonals;
// because they are present in A[], we do not use them.
   if (rmatrixbdmkl(a, m, n, &work, &t, tauq, taup)) {
      ae_frame_leave();
      return;
   }
// ALGLIB code
   if (m >= n) {
   // Reduce to upper bidiagonal form
      for (i = 0; i < n; i++) {
      // Generate elementary reflector H(i) to annihilate A(i+1:m-1,i)
         ae_v_move(&t.xR[1], 1, &a->xyR[i][i], a->stride, m - i);
         generatereflection(&t, m - i, &ltau);
         tauq->xR[i] = ltau;
         ae_v_move(&a->xyR[i][i], a->stride, &t.xR[1], 1, m - i);
         t.xR[1] = 1.0;
      // Apply H(i) to A(i:m-1,i+1:n-1) from the left
         applyreflectionfromtheleft(a, ltau, &t, i, m - 1, i + 1, n - 1, &work);
         if (i < n - 1) {
         // Generate elementary reflector G(i) to annihilate
         // A(i,i+2:n-1)
            ae_v_move(&t.xR[1], 1, &a->xyR[i][i + 1], 1, n - i - 1);
            generatereflection(&t, n - 1 - i, &ltau);
            taup->xR[i] = ltau;
            ae_v_move(&a->xyR[i][i + 1], 1, &t.xR[1], 1, n - i - 1);
            t.xR[1] = 1.0;
         // Apply G(i) to A(i+1:m-1,i+1:n-1) from the right
            applyreflectionfromtheright(a, ltau, &t, i + 1, m - 1, i + 1, n - 1, &work);
         } else {
            taup->xR[i] = 0.0;
         }
      }
   } else {
   // Reduce to lower bidiagonal form
      for (i = 0; i < m; i++) {
      // Generate elementary reflector G(i) to annihilate A(i,i+1:n-1)
         ae_v_move(&t.xR[1], 1, &a->xyR[i][i], 1, n - i);
         generatereflection(&t, n - i, &ltau);
         taup->xR[i] = ltau;
         ae_v_move(&a->xyR[i][i], 1, &t.xR[1], 1, n - i);
         t.xR[1] = 1.0;
      // Apply G(i) to A(i+1:m-1,i:n-1) from the right
         applyreflectionfromtheright(a, ltau, &t, i + 1, m - 1, i, n - 1, &work);
         if (i < m - 1) {
         // Generate elementary reflector H(i) to annihilate
         // A(i+2:m-1,i)
            ae_v_move(&t.xR[1], 1, &a->xyR[i + 1][i], a->stride, m - 1 - i);
            generatereflection(&t, m - 1 - i, &ltau);
            tauq->xR[i] = ltau;
            ae_v_move(&a->xyR[i + 1][i], a->stride, &t.xR[1], 1, m - i - 1);
            t.xR[1] = 1.0;
         // Apply H(i) to A(i+1:m-1,i+1:n-1) from the left
            applyreflectionfromtheleft(a, ltau, &t, i + 1, m - 1, i + 1, n - 1, &work);
         } else {
            tauq->xR[i] = 0.0;
         }
      }
   }
   ae_frame_leave();
}

// Unpacking matrix Q which reduces a matrix to bidiagonal form.
//
// Inputs:
//     QP          -   matrices Q and P in compact form.
//                     Output of ToBidiagonal subroutine.
//     M           -   number of rows in matrix A.
//     N           -   number of columns in matrix A.
//     TAUQ        -   scalar factors which are used to form Q.
//                     Output of ToBidiagonal subroutine.
//     QColumns    -   required number of columns in matrix Q.
//                     M >= QColumns >= 0.
//
// Outputs:
//     Q           -   first QColumns columns of matrix Q.
//                     Array[0..M-1, 0..QColumns-1]
//                     If QColumns=0, the array is not modified.
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixbdunpackq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, const ae_int_t qcolumns, real_2d_array &q);
void rmatrixbdunpackq(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, ae_int_t qcolumns, RMatrix *q) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(q);
   ae_assert(qcolumns <= m, "RMatrixBDUnpackQ: QColumns>M!");
   ae_assert(qcolumns >= 0, "RMatrixBDUnpackQ: QColumns<0!");
   if (m == 0 || n == 0 || qcolumns == 0) {
      return;
   }
// prepare Q
   ae_matrix_set_length(q, m, qcolumns);
   for (i = 0; i < m; i++) {
      for (j = 0; j < qcolumns; j++) {
         if (i == j) {
            q->xyR[i][j] = 1.0;
         } else {
            q->xyR[i][j] = 0.0;
         }
      }
   }
// Calculate
   rmatrixbdmultiplybyq(qp, m, n, tauq, q, m, qcolumns, false, false);
}

// Multiplication by matrix Q which reduces matrix A to  bidiagonal form.
//
// The algorithm allows pre- or post-multiply by Q or Q'.
//
// Inputs:
//     QP          -   matrices Q and P in compact form.
//                     Output of ToBidiagonal subroutine.
//     M           -   number of rows in matrix A.
//     N           -   number of columns in matrix A.
//     TAUQ        -   scalar factors which are used to form Q.
//                     Output of ToBidiagonal subroutine.
//     Z           -   multiplied matrix.
//                     array[0..ZRows-1,0..ZColumns-1]
//     ZRows       -   number of rows in matrix Z. If FromTheRight=False,
//                     ZRows=M, otherwise ZRows can be arbitrary.
//     ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
//                     ZColumns=M, otherwise ZColumns can be arbitrary.
//     FromTheRight -  pre- or post-multiply.
//     DoTranspose -   multiply by Q or Q'.
//
// Outputs:
//     Z           -   product of Z and Q.
//                     Array[0..ZRows-1,0..ZColumns-1]
//                     If ZRows=0 or ZColumns=0, the array is not modified.
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixbdmultiplybyq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose);
void rmatrixbdmultiplybyq(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *tauq, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool fromtheright, bool dotranspose) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t istep;
   ae_int_t mx;
   ae_frame_make(&_frame_block);
   NewVector(v, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   NewVector(dummy, 0, DT_REAL);
   if (m <= 0 || n <= 0 || zrows <= 0 || zcolumns <= 0) {
      ae_frame_leave();
      return;
   }
   ae_assert((fromtheright? zcolumns: zrows) == m, "RMatrixBDMultiplyByQ: incorrect Z size!");
// Try to use MKL code
   if (rmatrixbdmultiplybymkl(qp, m, n, tauq, &dummy, z, zrows, zcolumns, true, fromtheright, dotranspose)) {
      ae_frame_leave();
      return;
   }
// init
   mx = imax2(m, n);
   mx = imax2(mx, zrows);
   mx = imax2(mx, zcolumns);
   ae_vector_set_length(&v, mx + 1);
   ae_vector_set_length(&work, mx + 1);
   if (m >= n) {
   // setup
      if (fromtheright) {
         i1 = 0;
         i2 = n - 1;
         istep = 1;
      } else {
         i1 = n - 1;
         i2 = 0;
         istep = -1;
      }
      if (dotranspose) {
         i = i1;
         i1 = i2;
         i2 = i;
         istep = -istep;
      }
   // Process
      i = i1;
      do {
         ae_v_move(&v.xR[1], 1, &qp->xyR[i][i], qp->stride, m - i);
         v.xR[1] = 1.0;
         if (fromtheright) {
            applyreflectionfromtheright(z, tauq->xR[i], &v, 0, zrows - 1, i, m - 1, &work);
         } else {
            applyreflectionfromtheleft(z, tauq->xR[i], &v, i, m - 1, 0, zcolumns - 1, &work);
         }
         i += istep;
      } while (i != i2 + istep);
   } else {
   // setup
      if (fromtheright) {
         i1 = 0;
         i2 = m - 2;
         istep = 1;
      } else {
         i1 = m - 2;
         i2 = 0;
         istep = -1;
      }
      if (dotranspose) {
         i = i1;
         i1 = i2;
         i2 = i;
         istep = -istep;
      }
   // Process
      if (m - 1 > 0) {
         i = i1;
         do {
            ae_v_move(&v.xR[1], 1, &qp->xyR[i + 1][i], qp->stride, m - i - 1);
            v.xR[1] = 1.0;
            if (fromtheright) {
               applyreflectionfromtheright(z, tauq->xR[i], &v, 0, zrows - 1, i + 1, m - 1, &work);
            } else {
               applyreflectionfromtheleft(z, tauq->xR[i], &v, i + 1, m - 1, 0, zcolumns - 1, &work);
            }
            i += istep;
         } while (i != i2 + istep);
      }
   }
   ae_frame_leave();
}

// Unpacking matrix P which reduces matrix A to bidiagonal form.
// The subroutine returns transposed matrix P.
//
// Inputs:
//     QP      -   matrices Q and P in compact form.
//                 Output of ToBidiagonal subroutine.
//     M       -   number of rows in matrix A.
//     N       -   number of columns in matrix A.
//     TAUP    -   scalar factors which are used to form P.
//                 Output of ToBidiagonal subroutine.
//     PTRows  -   required number of rows of matrix P^T. N >= PTRows >= 0.
//
// Outputs:
//     PT      -   first PTRows columns of matrix P^T
//                 Array[0..PTRows-1, 0..N-1]
//                 If PTRows=0, the array is not modified.
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixbdunpackpt(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, const ae_int_t ptrows, real_2d_array &pt);
void rmatrixbdunpackpt(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *taup, ae_int_t ptrows, RMatrix *pt) {
   ae_int_t i;
   ae_int_t j;
   SetMatrix(pt);
   ae_assert(ptrows <= n, "RMatrixBDUnpackPT: PTRows>N!");
   ae_assert(ptrows >= 0, "RMatrixBDUnpackPT: PTRows<0!");
   if (m == 0 || n == 0 || ptrows == 0) {
      return;
   }
// prepare PT
   ae_matrix_set_length(pt, ptrows, n);
   for (i = 0; i < ptrows; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            pt->xyR[i][j] = 1.0;
         } else {
            pt->xyR[i][j] = 0.0;
         }
      }
   }
// Calculate
   rmatrixbdmultiplybyp(qp, m, n, taup, pt, ptrows, n, true, true);
}

// Multiplication by matrix P which reduces matrix A to  bidiagonal form.
//
// The algorithm allows pre- or post-multiply by P or P'.
//
// Inputs:
//     QP          -   matrices Q and P in compact form.
//                     Output of RMatrixBD subroutine.
//     M           -   number of rows in matrix A.
//     N           -   number of columns in matrix A.
//     TAUP        -   scalar factors which are used to form P.
//                     Output of RMatrixBD subroutine.
//     Z           -   multiplied matrix.
//                     Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
//     ZRows       -   number of rows in matrix Z. If FromTheRight=False,
//                     ZRows=N, otherwise ZRows can be arbitrary.
//     ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
//                     ZColumns=N, otherwise ZColumns can be arbitrary.
//     FromTheRight -  pre- or post-multiply.
//     DoTranspose -   multiply by P or P'.
//
// Outputs:
//     Z - product of Z and P.
//                 Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
//                 If ZRows=0 or ZColumns=0, the array is not modified.
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixbdmultiplybyp(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose);
void rmatrixbdmultiplybyp(RMatrix *qp, ae_int_t m, ae_int_t n, RVector *taup, RMatrix *z, ae_int_t zrows, ae_int_t zcolumns, bool fromtheright, bool dotranspose) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t mx;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t istep;
   ae_frame_make(&_frame_block);
   NewVector(v, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   NewVector(dummy, 0, DT_REAL);
   if (m <= 0 || n <= 0 || zrows <= 0 || zcolumns <= 0) {
      ae_frame_leave();
      return;
   }
   ae_assert((fromtheright? zcolumns: zrows) == n, "RMatrixBDMultiplyByP: incorrect Z size!");
// init
   mx = imax2(m, n);
   mx = imax2(mx, zrows);
   mx = imax2(mx, zcolumns);
   ae_vector_set_length(&v, mx + 1);
   ae_vector_set_length(&work, mx + 1);
   if (m >= n) {
   // setup
      if (fromtheright) {
         i1 = n - 2;
         i2 = 0;
         istep = -1;
      } else {
         i1 = 0;
         i2 = n - 2;
         istep = 1;
      }
      if (!dotranspose) {
         i = i1;
         i1 = i2;
         i2 = i;
         istep = -istep;
      }
   // Process
      if (n - 1 > 0) {
         i = i1;
         do {
            ae_v_move(&v.xR[1], 1, &qp->xyR[i][i + 1], 1, n - 1 - i);
            v.xR[1] = 1.0;
            if (fromtheright) {
               applyreflectionfromtheright(z, taup->xR[i], &v, 0, zrows - 1, i + 1, n - 1, &work);
            } else {
               applyreflectionfromtheleft(z, taup->xR[i], &v, i + 1, n - 1, 0, zcolumns - 1, &work);
            }
            i += istep;
         } while (i != i2 + istep);
      }
   } else {
   // setup
      if (fromtheright) {
         i1 = m - 1;
         i2 = 0;
         istep = -1;
      } else {
         i1 = 0;
         i2 = m - 1;
         istep = 1;
      }
      if (!dotranspose) {
         i = i1;
         i1 = i2;
         i2 = i;
         istep = -istep;
      }
   // Process
      i = i1;
      do {
         ae_v_move(&v.xR[1], 1, &qp->xyR[i][i], 1, n - i);
         v.xR[1] = 1.0;
         if (fromtheright) {
            applyreflectionfromtheright(z, taup->xR[i], &v, 0, zrows - 1, i, n - 1, &work);
         } else {
            applyreflectionfromtheleft(z, taup->xR[i], &v, i, n - 1, 0, zcolumns - 1, &work);
         }
         i += istep;
      } while (i != i2 + istep);
   }
   ae_frame_leave();
}

// Unpacking of the main and secondary diagonals of bidiagonal decomposition
// of matrix A.
//
// Inputs:
//     B   -   output of RMatrixBD subroutine.
//     M   -   number of rows in matrix B.
//     N   -   number of columns in matrix B.
//
// Outputs:
//     IsUpper -   True, if the matrix is upper bidiagonal.
//                 otherwise IsUpper is False.
//     D       -   the main diagonal.
//                 Array whose index ranges within [0..Min(M,N)-1].
//     E       -   the secondary diagonal (upper or lower, depending on
//                 the value of IsUpper).
//                 Array index ranges within [0..Min(M,N)-1], the last
//                 element is not used.
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixbdunpackdiagonals(const real_2d_array &b, const ae_int_t m, const ae_int_t n, bool &isupper, real_1d_array &d, real_1d_array &e);
void rmatrixbdunpackdiagonals(RMatrix *b, ae_int_t m, ae_int_t n, bool *isupper, RVector *d, RVector *e) {
   ae_int_t i;
   *isupper = false;
   SetVector(d);
   SetVector(e);
   *isupper = m >= n;
   if (m <= 0 || n <= 0) {
      return;
   }
   if (*isupper) {
      ae_vector_set_length(d, n);
      ae_vector_set_length(e, n);
      for (i = 0; i < n - 1; i++) {
         d->xR[i] = b->xyR[i][i];
         e->xR[i] = b->xyR[i][i + 1];
      }
      d->xR[n - 1] = b->xyR[n - 1][n - 1];
   } else {
      ae_vector_set_length(d, m);
      ae_vector_set_length(e, m);
      for (i = 0; i < m - 1; i++) {
         d->xR[i] = b->xyR[i][i];
         e->xR[i] = b->xyR[i + 1][i];
      }
      d->xR[m - 1] = b->xyR[m - 1][m - 1];
   }
}

// Reduction of a square matrix to  upper Hessenberg form: Q'*A*Q = H,
// where Q is an orthogonal matrix, H - Hessenberg matrix.
//
// Inputs:
//     A       -   matrix A with elements [0..N-1, 0..N-1]
//     N       -   size of matrix A.
//
// Outputs:
//     A       -   matrices Q and P in  compact form (see below).
//     Tau     -   array of scalar factors which are used to form matrix Q.
//                 Array whose index ranges within [0..N-2]
//
// Matrix H is located on the main diagonal, on the lower secondary  diagonal
// and above the main diagonal of matrix A. The elements which are used to
// form matrix Q are situated in array Tau and below the lower secondary
// diagonal of matrix A as follows:
//
// Matrix Q is represented as a product of elementary reflections
//
// Q = H(0)*H(2)*...*H(n-2),
//
// where each H(i) is given by
//
// H(i) = 1 - tau * v * (v^T)
//
// where tau is a scalar stored in Tau[I]; v - is a real vector,
// so that v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) stored in A(i+2:n-1,i).
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
// API: void rmatrixhessenberg(real_2d_array &a, const ae_int_t n, real_1d_array &tau);
void rmatrixhessenberg(RMatrix *a, ae_int_t n, RVector *tau) {
   ae_frame _frame_block;
   ae_int_t i;
   double v;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   NewVector(t, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   ae_assert(n >= 0, "RMatrixHessenberg: incorrect N!");
// Quick return if possible
   if (n <= 1) {
      ae_frame_leave();
      return;
   }
// Allocate place
   ae_vector_set_length(tau, n - 1);
   ae_vector_set_length(&t, n + 1);
   ae_vector_set_length(&work, n);
// MKL version
   if (rmatrixhessenbergmkl(a, n, tau)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version
   for (i = 0; i < n - 1; i++) {
   // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
      ae_v_move(&t.xR[1], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
      generatereflection(&t, n - i - 1, &v);
      ae_v_move(&a->xyR[i + 1][i], a->stride, &t.xR[1], 1, n - i - 1);
      tau->xR[i] = v;
      t.xR[1] = 1.0;
   // Apply H(i) to A(1:ihi,i+1:ihi) from the right
      applyreflectionfromtheright(a, v, &t, 0, n - 1, i + 1, n - 1, &work);
   // Apply H(i) to A(i+1:ihi,i+1:n) from the left
      applyreflectionfromtheleft(a, v, &t, i + 1, n - 1, i + 1, n - 1, &work);
   }
   ae_frame_leave();
}

// Unpacking matrix Q which reduces matrix A to upper Hessenberg form
//
// Inputs:
//     A   -   output of RMatrixHessenberg subroutine.
//     N   -   size of matrix A.
//     Tau -   scalar factors which are used to form Q.
//             Output of RMatrixHessenberg subroutine.
//
// Outputs:
//     Q   -   matrix Q.
//             Array whose indexes range within [0..N-1, 0..N-1].
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixhessenbergunpackq(const real_2d_array &a, const ae_int_t n, const real_1d_array &tau, real_2d_array &q);
void rmatrixhessenbergunpackq(RMatrix *a, ae_int_t n, RVector *tau, RMatrix *q) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(v, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   if (n == 0) {
      ae_frame_leave();
      return;
   }
// init
   ae_matrix_set_length(q, n, n);
   ae_vector_set_length(&v, n);
   ae_vector_set_length(&work, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            q->xyR[i][j] = 1.0;
         } else {
            q->xyR[i][j] = 0.0;
         }
      }
   }
// MKL version
   if (rmatrixhessenbergunpackqmkl(a, n, tau, q)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version: unpack Q
   for (i = 0; i < n - 1; i++) {
   // Apply H(i)
      ae_v_move(&v.xR[1], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
      v.xR[1] = 1.0;
      applyreflectionfromtheright(q, tau->xR[i], &v, 0, n - 1, i + 1, n - 1, &work);
   }
   ae_frame_leave();
}

// Unpacking matrix H (the result of matrix A reduction to upper Hessenberg form)
//
// Inputs:
//     A   -   output of RMatrixHessenberg subroutine.
//     N   -   size of matrix A.
//
// Outputs:
//     H   -   matrix H. Array whose indexes range within [0..N-1, 0..N-1].
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void rmatrixhessenbergunpackh(const real_2d_array &a, const ae_int_t n, real_2d_array &h);
void rmatrixhessenbergunpackh(RMatrix *a, ae_int_t n, RMatrix *h) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   SetMatrix(h);
   NewVector(v, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   if (n == 0) {
      ae_frame_leave();
      return;
   }
   ae_matrix_set_length(h, n, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < i - 1; j++) {
         h->xyR[i][j] = 0.0;
      }
      j = imax2(0, i - 1);
      ae_v_move(&h->xyR[i][j], 1, &a->xyR[i][j], 1, n - j);
   }
   ae_frame_leave();
}

// Reduction of a symmetric matrix which is given by its higher or lower
// triangular part to a tridiagonal matrix using orthogonal similarity
// transformation: Q'*A*Q=T.
//
// Inputs:
//     A       -   matrix to be transformed
//                 array with elements [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   storage format. If IsUpper = True, then matrix A is given
//                 by its upper triangle, and the lower triangle is not used
//                 and not modified by the algorithm, and vice versa
//                 if IsUpper = False.
//
// Outputs:
//     A       -   matrices T and Q in  compact form (see lower)
//     Tau     -   array of factors which are forming matrices H(i)
//                 array with elements [0..N-2].
//     D       -   main diagonal of symmetric matrix T.
//                 array with elements [0..N-1].
//     E       -   secondary diagonal of symmetric matrix T.
//                 array with elements [0..N-2].
//
//
//   If IsUpper=True, the matrix Q is represented as a product of elementary
//   reflectors
//
//      Q = H(n-2) . . . H(2) H(0).
//
//   Each H(i) has the form
//
//      H(i) = I - tau * v * v'
//
//   where tau is a real scalar, and v is a real vector with
//   v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
//   A(0:i-1,i+1), and tau in TAU(i).
//
//   If IsUpper=False, the matrix Q is represented as a product of elementary
//   reflectors
//
//      Q = H(0) H(2) . . . H(n-2).
//
//   Each H(i) has the form
//
//      H(i) = I - tau * v * v'
//
//   where tau is a real scalar, and v is a real vector with
//   v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
//   and tau in TAU(i).
//
//   The contents of A on exit are illustrated by the following examples
//   with n = 5:
//
//   if UPLO = 'U':                       if UPLO = 'L':
//
//     (  d   e   v1  v2  v3 )              (  d                  )
//     (      d   e   v2  v3 )              (  e   d              )
//     (          d   e   v3 )              (  v0  e   d          )
//     (              d   e  )              (  v0  v1  e   d      )
//     (                  d  )              (  v0  v1  v2  e   d  )
//
//   where d and e denote diagonal and off-diagonal elements of T, and vi
//   denotes an element of the vector defining H(i).
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
// API: void smatrixtd(real_2d_array &a, const ae_int_t n, const bool isupper, real_1d_array &tau, real_1d_array &d, real_1d_array &e);
void smatrixtd(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RVector *d, RVector *e) {
   ae_frame _frame_block;
   ae_int_t i;
   double alpha;
   double taui;
   double v;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   SetVector(d);
   SetVector(e);
   NewVector(t, 0, DT_REAL);
   NewVector(t2, 0, DT_REAL);
   NewVector(t3, 0, DT_REAL);
   if (n <= 0) {
      ae_frame_leave();
      return;
   }
   ae_vector_set_length(&t, n + 1);
   ae_vector_set_length(&t2, n + 1);
   ae_vector_set_length(&t3, n + 1);
   if (n > 1) {
      ae_vector_set_length(tau, n - 1);
   }
   ae_vector_set_length(d, n);
   if (n > 1) {
      ae_vector_set_length(e, n - 1);
   }
// Try to use MKL
   if (smatrixtdmkl(a, n, isupper, tau, d, e)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version
   if (isupper) {
   // Reduce the upper triangle of A
      for (i = n - 2; i >= 0; i--) {
      // Generate elementary reflector H() = E - tau * v * v'
         if (i >= 1) {
            ae_v_move(&t.xR[2], 1, &a->xyR[0][i + 1], a->stride, i);
         }
         t.xR[1] = a->xyR[i][i + 1];
         generatereflection(&t, i + 1, &taui);
         if (i >= 1) {
            ae_v_move(&a->xyR[0][i + 1], a->stride, &t.xR[2], 1, i);
         }
         a->xyR[i][i + 1] = t.xR[1];
         e->xR[i] = a->xyR[i][i + 1];
         if (taui != 0.0) {
         // Apply H from both sides to A
            a->xyR[i][i + 1] = 1.0;
         // Compute  x := tau * A * v  storing x in TAU
            ae_v_move(&t.xR[1], 1, &a->xyR[0][i + 1], a->stride, i + 1);
            symmetricmatrixvectormultiply(a, isupper, 0, i, &t, taui, &t3);
            ae_v_move(tau->xR, 1, &t3.xR[1], 1, i + 1);
         // Compute  w := x - 1/2 * tau * (x'*v) * v
            v = ae_v_dotproduct(tau->xR, 1, &a->xyR[0][i + 1], a->stride, i + 1);
            alpha = -0.5 * taui * v;
            ae_v_addd(tau->xR, 1, &a->xyR[0][i + 1], a->stride, i + 1, alpha);
         // Apply the transformation as a rank-2 update:
         //    A := A - v * w' - w * v'
            ae_v_move(&t.xR[1], 1, &a->xyR[0][i + 1], a->stride, i + 1);
            ae_v_move(&t3.xR[1], 1, tau->xR, 1, i + 1);
            symmetricrank2update(a, isupper, 0, i, &t, &t3, &t2, -1.0);
            a->xyR[i][i + 1] = e->xR[i];
         }
         d->xR[i + 1] = a->xyR[i + 1][i + 1];
         tau->xR[i] = taui;
      }
      d->xR[0] = a->xyR[0][0];
   } else {
   // Reduce the lower triangle of A
      for (i = 0; i < n - 1; i++) {
      // Generate elementary reflector H = E - tau * v * v'
         ae_v_move(&t.xR[1], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
         generatereflection(&t, n - i - 1, &taui);
         ae_v_move(&a->xyR[i + 1][i], a->stride, &t.xR[1], 1, n - i - 1);
         e->xR[i] = a->xyR[i + 1][i];
         if (taui != 0.0) {
         // Apply H from both sides to A
            a->xyR[i + 1][i] = 1.0;
         // Compute  x := tau * A * v  storing y in TAU
            ae_v_move(&t.xR[1], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
            symmetricmatrixvectormultiply(a, isupper, i + 1, n - 1, &t, taui, &t2);
            ae_v_move(&tau->xR[i], 1, &t2.xR[1], 1, n - i - 1);
         // Compute  w := x - 1/2 * tau * (x'*v) * v
            v = ae_v_dotproduct(&tau->xR[i], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
            alpha = -0.5 * taui * v;
            ae_v_addd(&tau->xR[i], 1, &a->xyR[i + 1][i], a->stride, n - i - 1, alpha);
         // Apply the transformation as a rank-2 update:
         //     A := A - v * w' - w * v'
         //
            ae_v_move(&t.xR[1], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
            ae_v_move(&t2.xR[1], 1, &tau->xR[i], 1, n - i - 1);
            symmetricrank2update(a, isupper, i + 1, n - 1, &t, &t2, &t3, -1.0);
            a->xyR[i + 1][i] = e->xR[i];
         }
         d->xR[i] = a->xyR[i][i];
         tau->xR[i] = taui;
      }
      d->xR[n - 1] = a->xyR[n - 1][n - 1];
   }
   ae_frame_leave();
}

// Reduction of a Hermitian matrix which is given  by  its  higher  or  lower
// triangular part to a real  tridiagonal  matrix  using  unitary  similarity
// transformation: Q'*A*Q = T.
//
// Inputs:
//     A       -   matrix to be transformed
//                 array with elements [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   storage format. If IsUpper = True, then matrix A is  given
//                 by its upper triangle, and the lower triangle is not  used
//                 and not modified by the algorithm, and vice versa
//                 if IsUpper = False.
//
// Outputs:
//     A       -   matrices T and Q in  compact form (see lower)
//     Tau     -   array of factors which are forming matrices H(i)
//                 array with elements [0..N-2].
//     D       -   main diagonal of real symmetric matrix T.
//                 array with elements [0..N-1].
//     E       -   secondary diagonal of real symmetric matrix T.
//                 array with elements [0..N-2].
//
//
//   If IsUpper=True, the matrix Q is represented as a product of elementary
//   reflectors
//
//      Q = H(n-2) . . . H(2) H(0).
//
//   Each H(i) has the form
//
//      H(i) = I - tau * v * v'
//
//   where tau is a complex scalar, and v is a complex vector with
//   v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
//   A(0:i-1,i+1), and tau in TAU(i).
//
//   If IsUpper=False, the matrix Q is represented as a product of elementary
//   reflectors
//
//      Q = H(0) H(2) . . . H(n-2).
//
//   Each H(i) has the form
//
//      H(i) = I - tau * v * v'
//
//   where tau is a complex scalar, and v is a complex vector with
//   v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
//   and tau in TAU(i).
//
//   The contents of A on exit are illustrated by the following examples
//   with n = 5:
//
//   if UPLO = 'U':                       if UPLO = 'L':
//
//     (  d   e   v1  v2  v3 )              (  d                  )
//     (      d   e   v2  v3 )              (  e   d              )
//     (          d   e   v3 )              (  v0  e   d          )
//     (              d   e  )              (  v0  v1  e   d      )
//     (                  d  )              (  v0  v1  v2  e   d  )
//
// where d and e denote diagonal and off-diagonal elements of T, and vi
// denotes an element of the vector defining H(i).
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
// API: void hmatrixtd(complex_2d_array &a, const ae_int_t n, const bool isupper, complex_1d_array &tau, real_1d_array &d, real_1d_array &e);
void hmatrixtd(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, RVector *d, RVector *e) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_complex alpha;
   ae_complex taui;
   ae_complex v;
   ae_frame_make(&_frame_block);
   SetVector(tau);
   SetVector(d);
   SetVector(e);
   NewVector(t, 0, DT_COMPLEX);
   NewVector(t2, 0, DT_COMPLEX);
   NewVector(t3, 0, DT_COMPLEX);
// Init and test
   if (n <= 0) {
      ae_frame_leave();
      return;
   }
   for (i = 0; i < n; i++) {
      ae_assert(a->xyC[i][i].y == 0.0, "Assertion failed");
   }
   if (n > 1) {
      ae_vector_set_length(tau, n - 1);
      ae_vector_set_length(e, n - 1);
   }
   ae_vector_set_length(d, n);
   ae_vector_set_length(&t, n);
   ae_vector_set_length(&t2, n);
   ae_vector_set_length(&t3, n);
// MKL version
   if (hmatrixtdmkl(a, n, isupper, tau, d, e)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version
   if (isupper) {
   // Reduce the upper triangle of A
      a->xyC[n - 1][n - 1] = ae_complex_from_d(a->xyC[n - 1][n - 1].x);
      for (i = n - 2; i >= 0; i--) {
      // Generate elementary reflector H = I+1 - tau * v * v'
         alpha = a->xyC[i][i + 1];
         t.xC[1] = alpha;
         if (i >= 1) {
            ae_v_cmove(&t.xC[2], 1, &a->xyC[0][i + 1], a->stride, "N", i);
         }
         complexgeneratereflection(&t, i + 1, &taui);
         if (i >= 1) {
            ae_v_cmove(&a->xyC[0][i + 1], a->stride, &t.xC[2], 1, "N", i);
         }
         alpha = t.xC[1];
         e->xR[i] = alpha.x;
         if (ae_c_neq_d(taui, 0.0)) {
         // Apply H(I+1) from both sides to A
            a->xyC[i][i + 1] = ae_complex_from_i(1);
         // Compute  x := tau * A * v  storing x in TAU
            ae_v_cmove(&t.xC[1], 1, &a->xyC[0][i + 1], a->stride, "N", i + 1);
            hermitianmatrixvectormultiply(a, isupper, 0, i, &t, taui, &t2);
            ae_v_cmove(tau->xC, 1, &t2.xC[1], 1, "N", i + 1);
         // Compute  w := x - 1/2 * tau * (x'*v) * v
            v = ae_v_cdotproduct(tau->xC, 1, "Conj", &a->xyC[0][i + 1], a->stride, "N", i + 1);
            alpha = ae_c_neg(ae_c_mul(ae_c_mul_d(taui, 0.5), v));
            ae_v_caddc(tau->xC, 1, &a->xyC[0][i + 1], a->stride, "N", i + 1, alpha);
         // Apply the transformation as a rank-2 update:
         //    A := A - v * w' - w * v'
            ae_v_cmove(&t.xC[1], 1, &a->xyC[0][i + 1], a->stride, "N", i + 1);
            ae_v_cmove(&t3.xC[1], 1, tau->xC, 1, "N", i + 1);
            hermitianrank2update(a, isupper, 0, i, &t, &t3, &t2, ae_complex_from_i(-1));
         } else {
            a->xyC[i][i] = ae_complex_from_d(a->xyC[i][i].x);
         }
         a->xyC[i][i + 1] = ae_complex_from_d(e->xR[i]);
         d->xR[i + 1] = a->xyC[i + 1][i + 1].x;
         tau->xC[i] = taui;
      }
      d->xR[0] = a->xyC[0][0].x;
   } else {
   // Reduce the lower triangle of A
      a->xyC[0][0] = ae_complex_from_d(a->xyC[0][0].x);
      for (i = 0; i < n - 1; i++) {
      // Generate elementary reflector H = I - tau * v * v'
         ae_v_cmove(&t.xC[1], 1, &a->xyC[i + 1][i], a->stride, "N", n - i - 1);
         complexgeneratereflection(&t, n - i - 1, &taui);
         ae_v_cmove(&a->xyC[i + 1][i], a->stride, &t.xC[1], 1, "N", n - i - 1);
         e->xR[i] = a->xyC[i + 1][i].x;
         if (ae_c_neq_d(taui, 0.0)) {
         // Apply H(i) from both sides to A(i+1:n,i+1:n)
            a->xyC[i + 1][i] = ae_complex_from_i(1);
         // Compute  x := tau * A * v  storing y in TAU
            ae_v_cmove(&t.xC[1], 1, &a->xyC[i + 1][i], a->stride, "N", n - i - 1);
            hermitianmatrixvectormultiply(a, isupper, i + 1, n - 1, &t, taui, &t2);
            ae_v_cmove(&tau->xC[i], 1, &t2.xC[1], 1, "N", n - i - 1);
         // Compute  w := x - 1/2 * tau * (x'*v) * v
            v = ae_v_cdotproduct(&tau->xC[i], 1, "Conj", &a->xyC[i + 1][i], a->stride, "N", n - i - 1);
            alpha = ae_c_neg(ae_c_mul(ae_c_mul_d(taui, 0.5), v));
            ae_v_caddc(&tau->xC[i], 1, &a->xyC[i + 1][i], a->stride, "N", n - i - 1, alpha);
         // Apply the transformation as a rank-2 update:
         // A := A - v * w' - w * v'
            ae_v_cmove(&t.xC[1], 1, &a->xyC[i + 1][i], a->stride, "N", n - i - 1);
            ae_v_cmove(&t2.xC[1], 1, &tau->xC[i], 1, "N", n - i - 1);
            hermitianrank2update(a, isupper, i + 1, n - 1, &t, &t2, &t3, ae_complex_from_i(-1));
         } else {
            a->xyC[i + 1][i + 1] = ae_complex_from_d(a->xyC[i + 1][i + 1].x);
         }
         a->xyC[i + 1][i] = ae_complex_from_d(e->xR[i]);
         d->xR[i] = a->xyC[i][i].x;
         tau->xC[i] = taui;
      }
      d->xR[n - 1] = a->xyC[n - 1][n - 1].x;
   }
   ae_frame_leave();
}

// Unpacking matrix Q which reduces symmetric matrix to a tridiagonal
// form.
//
// Inputs:
//     A       -   the result of a SMatrixTD subroutine
//     N       -   size of matrix A.
//     IsUpper -   storage format (a parameter of SMatrixTD subroutine)
//     Tau     -   the result of a SMatrixTD subroutine
//
// Outputs:
//     Q       -   transformation matrix.
//                 array with elements [0..N-1, 0..N-1].
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void smatrixtdunpackq(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &tau, real_2d_array &q);
void smatrixtdunpackq(RMatrix *a, ae_int_t n, bool isupper, RVector *tau, RMatrix *q) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(v, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   if (n == 0) {
      ae_frame_leave();
      return;
   }
// init
   ae_matrix_set_length(q, n, n);
   ae_vector_set_length(&v, n + 1);
   ae_vector_set_length(&work, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            q->xyR[i][j] = 1.0;
         } else {
            q->xyR[i][j] = 0.0;
         }
      }
   }
// MKL version
   if (smatrixtdunpackqmkl(a, n, isupper, tau, q)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version: unpack Q
   if (isupper) {
      for (i = 0; i < n - 1; i++) {
      // Apply H(i)
         ae_v_move(&v.xR[1], 1, &a->xyR[0][i + 1], a->stride, i + 1);
         v.xR[i + 1] = 1.0;
         applyreflectionfromtheleft(q, tau->xR[i], &v, 0, i, 0, n - 1, &work);
      }
   } else {
      for (i = n - 2; i >= 0; i--) {
      // Apply H(i)
         ae_v_move(&v.xR[1], 1, &a->xyR[i + 1][i], a->stride, n - i - 1);
         v.xR[1] = 1.0;
         applyreflectionfromtheleft(q, tau->xR[i], &v, i + 1, n - 1, 0, n - 1, &work);
      }
   }
   ae_frame_leave();
}

// Unpacking matrix Q which reduces a Hermitian matrix to a real  tridiagonal
// form.
//
// Inputs:
//     A       -   the result of a HMatrixTD subroutine
//     N       -   size of matrix A.
//     IsUpper -   storage format (a parameter of HMatrixTD subroutine)
//     Tau     -   the result of a HMatrixTD subroutine
//
// Outputs:
//     Q       -   transformation matrix.
//                 array with elements [0..N-1, 0..N-1].
// ALGLIB: Copyright 2005-2010 by Sergey Bochkanov
// API: void hmatrixtdunpackq(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &tau, complex_2d_array &q);
void hmatrixtdunpackq(CMatrix *a, ae_int_t n, bool isupper, CVector *tau, CMatrix *q) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   SetMatrix(q);
   NewVector(v, 0, DT_COMPLEX);
   NewVector(work, 0, DT_COMPLEX);
   if (n == 0) {
      ae_frame_leave();
      return;
   }
// init
   ae_matrix_set_length(q, n, n);
   ae_vector_set_length(&v, n + 1);
   ae_vector_set_length(&work, n);
// MKL version
   if (hmatrixtdunpackqmkl(a, n, isupper, tau, q)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (i == j) {
            q->xyC[i][j] = ae_complex_from_i(1);
         } else {
            q->xyC[i][j] = ae_complex_from_i(0);
         }
      }
   }
   if (isupper) {
      for (i = 0; i < n - 1; i++) {
      // Apply H(i)
         ae_v_cmove(&v.xC[1], 1, &a->xyC[0][i + 1], a->stride, "N", i + 1);
         v.xC[i + 1] = ae_complex_from_i(1);
         complexapplyreflectionfromtheleft(q, tau->xC[i], &v, 0, i, 0, n - 1, &work);
      }
   } else {
      for (i = n - 2; i >= 0; i--) {
      // Apply H(i)
         ae_v_cmove(&v.xC[1], 1, &a->xyC[i + 1][i], a->stride, "N", n - i - 1);
         v.xC[1] = ae_complex_from_i(1);
         complexapplyreflectionfromtheleft(q, tau->xC[i], &v, i + 1, n - 1, 0, n - 1, &work);
      }
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void rmatrixqr(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixqr(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau));
   alglib_impl::ae_state_clear();
}

void cmatrixqr(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixqr(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau));
   alglib_impl::ae_state_clear();
}

void rmatrixlq(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlq(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau));
   alglib_impl::ae_state_clear();
}

void cmatrixlq(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlq(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau));
   alglib_impl::ae_state_clear();
}

void rmatrixqrunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qcolumns, real_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixqrunpackq(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau), qcolumns, ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void cmatrixqrunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qcolumns, complex_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixqrunpackq(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau), qcolumns, ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void rmatrixqrunpackr(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixqrunpackr(ConstT(ae_matrix, a), m, n, ConstT(ae_matrix, r));
   alglib_impl::ae_state_clear();
}

void cmatrixqrunpackr(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &r) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixqrunpackr(ConstT(ae_matrix, a), m, n, ConstT(ae_matrix, r));
   alglib_impl::ae_state_clear();
}

void rmatrixlqunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qrows, real_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlqunpackq(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau), qrows, ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void cmatrixlqunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qrows, complex_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlqunpackq(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tau), qrows, ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void rmatrixlqunpackl(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &l) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixlqunpackl(ConstT(ae_matrix, a), m, n, ConstT(ae_matrix, l));
   alglib_impl::ae_state_clear();
}

void cmatrixlqunpackl(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &l) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::cmatrixlqunpackl(ConstT(ae_matrix, a), m, n, ConstT(ae_matrix, l));
   alglib_impl::ae_state_clear();
}

void rmatrixbd(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tauq, real_1d_array &taup) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixbd(ConstT(ae_matrix, a), m, n, ConstT(ae_vector, tauq), ConstT(ae_vector, taup));
   alglib_impl::ae_state_clear();
}

void rmatrixbdunpackq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, const ae_int_t qcolumns, real_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixbdunpackq(ConstT(ae_matrix, qp), m, n, ConstT(ae_vector, tauq), qcolumns, ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void rmatrixbdmultiplybyq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixbdmultiplybyq(ConstT(ae_matrix, qp), m, n, ConstT(ae_vector, tauq), ConstT(ae_matrix, z), zrows, zcolumns, fromtheright, dotranspose);
   alglib_impl::ae_state_clear();
}

void rmatrixbdunpackpt(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, const ae_int_t ptrows, real_2d_array &pt) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixbdunpackpt(ConstT(ae_matrix, qp), m, n, ConstT(ae_vector, taup), ptrows, ConstT(ae_matrix, pt));
   alglib_impl::ae_state_clear();
}

void rmatrixbdmultiplybyp(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixbdmultiplybyp(ConstT(ae_matrix, qp), m, n, ConstT(ae_vector, taup), ConstT(ae_matrix, z), zrows, zcolumns, fromtheright, dotranspose);
   alglib_impl::ae_state_clear();
}

void rmatrixbdunpackdiagonals(const real_2d_array &b, const ae_int_t m, const ae_int_t n, bool &isupper, real_1d_array &d, real_1d_array &e) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixbdunpackdiagonals(ConstT(ae_matrix, b), m, n, &isupper, ConstT(ae_vector, d), ConstT(ae_vector, e));
   alglib_impl::ae_state_clear();
}

void rmatrixhessenberg(real_2d_array &a, const ae_int_t n, real_1d_array &tau) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixhessenberg(ConstT(ae_matrix, a), n, ConstT(ae_vector, tau));
   alglib_impl::ae_state_clear();
}

void rmatrixhessenbergunpackq(const real_2d_array &a, const ae_int_t n, const real_1d_array &tau, real_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixhessenbergunpackq(ConstT(ae_matrix, a), n, ConstT(ae_vector, tau), ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void rmatrixhessenbergunpackh(const real_2d_array &a, const ae_int_t n, real_2d_array &h) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixhessenbergunpackh(ConstT(ae_matrix, a), n, ConstT(ae_matrix, h));
   alglib_impl::ae_state_clear();
}

void smatrixtd(real_2d_array &a, const ae_int_t n, const bool isupper, real_1d_array &tau, real_1d_array &d, real_1d_array &e) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::smatrixtd(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, tau), ConstT(ae_vector, d), ConstT(ae_vector, e));
   alglib_impl::ae_state_clear();
}

void hmatrixtd(complex_2d_array &a, const ae_int_t n, const bool isupper, complex_1d_array &tau, real_1d_array &d, real_1d_array &e) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hmatrixtd(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, tau), ConstT(ae_vector, d), ConstT(ae_vector, e));
   alglib_impl::ae_state_clear();
}

void smatrixtdunpackq(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &tau, real_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::smatrixtdunpackq(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, tau), ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}

void hmatrixtdunpackq(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &tau, complex_2d_array &q) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hmatrixtdunpackq(ConstT(ae_matrix, a), n, isupper, ConstT(ae_vector, tau), ConstT(ae_matrix, q));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === FBLS Package ===
// Depends on: ORTFAC
namespace alglib_impl {
// Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.
//
// This subroutine assumes that:
// * A*ScaleA is well scaled
// * A is well-conditioned, so no zero divisions or overflow may occur
//
// Inputs:
//     CHA     -   Cholesky decomposition of A
//     SqrtScaleA- square root of scale factor ScaleA
//     N       -   matrix size, N >= 0.
//     IsUpper -   storage type
//     XB      -   right part
//     Tmp     -   buffer; function automatically allocates it, if it is  too
//                 small.  It  can  be  reused  if function is called several
//                 times.
//
// Outputs:
//     XB      -   solution
//
// NOTE 1: no assertion or tests are done during algorithm operation
// NOTE 2: N=0 will force algorithm to silently return
// ALGLIB: Copyright 13.10.2010 by Sergey Bochkanov
void fblscholeskysolve(RMatrix *cha, double sqrtscalea, ae_int_t n, bool isupper, RVector *xb, RVector *tmp) {
   double v;
   if (n <= 0) {
      return;
   }
   if (tmp->cnt < n) {
      ae_vector_set_length(tmp, n);
   }
// Scale right part
   v = 1 / ae_sqr(sqrtscalea);
   ae_v_muld(xb->xR, 1, n, v);
// Solve A = L*L' or A=U'*U
   if (isupper) {
   // Solve U'*y=b first.
      rmatrixtrsv(n, cha, 0, 0, true, false, 1, xb, 0);
   // Solve U*x=y then.
      rmatrixtrsv(n, cha, 0, 0, true, false, 0, xb, 0);
   } else {
   // Solve L*y=b first
      rmatrixtrsv(n, cha, 0, 0, false, false, 0, xb, 0);
   // Solve L'*x=y then.
      rmatrixtrsv(n, cha, 0, 0, false, false, 1, xb, 0);
   }
}

// Fast basic linear solver: linear SPD CG
//
// Solves (A^T*A + alpha*I)*x = b where:
// * A is MxN matrix
// * alpha>0 is a scalar
// * I is NxN identity matrix
// * b is Nx1 vector
// * X is Nx1 unknown vector.
//
// N iterations of linear conjugate gradient are used to solve problem.
//
// Inputs:
//     A   -   array[M,N], matrix
//     M   -   number of rows
//     N   -   number of unknowns
//     B   -   array[N], right part
//     X   -   initial approxumation, array[N]
//     Buf -   buffer; function automatically allocates it, if it is too
//             small. It can be reused if function is called several times
//             with same M and N.
//
// Outputs:
//     X   -   improved solution
//
// NOTES:
// *   solver checks quality of improved solution. If (because of problem
//     condition number, numerical noise, etc.) new solution is WORSE than
//     original approximation, then original approximation is returned.
// *   solver assumes that both A, B, Alpha are well scaled (i.e. they are
//     less than sqrt(overflow) and greater than sqrt(underflow)).
// ALGLIB: Copyright 20.08.2009 by Sergey Bochkanov
void fblssolvecgx(RMatrix *a, ae_int_t m, ae_int_t n, double alpha, RVector *b, RVector *x, RVector *buf) {
   ae_int_t k;
   ae_int_t offsrk;
   ae_int_t offsrk1;
   ae_int_t offsxk;
   ae_int_t offsxk1;
   ae_int_t offspk;
   ae_int_t offspk1;
   ae_int_t offstmp1;
   ae_int_t offstmp2;
   ae_int_t bs;
   double e1;
   double e2;
   double rk2;
   double rk12;
   double pap;
   double s;
   double betak;
   double v1;
   double v2;
// Test for special case: B=0
   v1 = ae_v_dotproduct(b->xR, 1, b->xR, 1, n);
   if (v1 == 0.0) {
      for (k = 0; k < n; k++) {
         x->xR[k] = 0.0;
      }
      return;
   }
// Offsets inside Buf for:
// * R[K], R[K+1]
// * X[K], X[K+1]
// * P[K], P[K+1]
// * Tmp1 - array[M], Tmp2 - array[N]
   offsrk = 0;
   offsrk1 = offsrk + n;
   offsxk = offsrk1 + n;
   offsxk1 = offsxk + n;
   offspk = offsxk1 + n;
   offspk1 = offspk + n;
   offstmp1 = offspk1 + n;
   offstmp2 = offstmp1 + m;
   bs = offstmp2 + n;
   if (buf->cnt < bs) {
      ae_vector_set_length(buf, bs);
   }
// x(0) = x
   ae_v_move(&buf->xR[offsxk], 1, x->xR, 1, n);
// r(0) = b-A*x(0)
// RK2 = r(0)'*r(0)
   rmatrixmv(m, n, a, 0, 0, 0, buf, offsxk, buf, offstmp1);
   rmatrixmv(n, m, a, 0, 0, 1, buf, offstmp1, buf, offstmp2);
   ae_v_addd(&buf->xR[offstmp2], 1, &buf->xR[offsxk], 1, n, alpha);
   ae_v_move(&buf->xR[offsrk], 1, b->xR, 1, n);
   ae_v_sub(&buf->xR[offsrk], 1, &buf->xR[offstmp2], 1, n);
   rk2 = ae_v_dotproduct(&buf->xR[offsrk], 1, &buf->xR[offsrk], 1, n);
   ae_v_move(&buf->xR[offspk], 1, &buf->xR[offsrk], 1, n);
   e1 = sqrt(rk2);
// Cycle
   for (k = 0; k < n; k++) {
   // Calculate A*p(k) - store in Buf[OffsTmp2:OffsTmp2+N-1]
   // and p(k)'*A*p(k)  - store in PAP
   //
   // If PAP=0, break (iteration is over)
      rmatrixmv(m, n, a, 0, 0, 0, buf, offspk, buf, offstmp1);
      v1 = ae_v_dotproduct(&buf->xR[offstmp1], 1, &buf->xR[offstmp1], 1, m);
      v2 = ae_v_dotproduct(&buf->xR[offspk], 1, &buf->xR[offspk], 1, n);
      pap = v1 + alpha * v2;
      rmatrixmv(n, m, a, 0, 0, 1, buf, offstmp1, buf, offstmp2);
      ae_v_addd(&buf->xR[offstmp2], 1, &buf->xR[offspk], 1, n, alpha);
      if (pap == 0.0) {
         break;
      }
   // S = (r(k)'*r(k))/(p(k)'*A*p(k))
      s = rk2 / pap;
   // x(k+1) = x(k) + S*p(k)
      ae_v_move(&buf->xR[offsxk1], 1, &buf->xR[offsxk], 1, n);
      ae_v_addd(&buf->xR[offsxk1], 1, &buf->xR[offspk], 1, n, s);
   // r(k+1) = r(k) - S*A*p(k)
   // RK12 = r(k+1)'*r(k+1)
   //
   // Break if r(k+1) small enough (when compared to r(k))
      ae_v_move(&buf->xR[offsrk1], 1, &buf->xR[offsrk], 1, n);
      ae_v_subd(&buf->xR[offsrk1], 1, &buf->xR[offstmp2], 1, n, s);
      rk12 = ae_v_dotproduct(&buf->xR[offsrk1], 1, &buf->xR[offsrk1], 1, n);
      if (sqrt(rk12) <= 100 * ae_machineepsilon * sqrt(rk2)) {
      // X(k) = x(k+1) before exit -
      // - because we expect to find solution at x(k)
         ae_v_move(&buf->xR[offsxk], 1, &buf->xR[offsxk1], 1, n);
         break;
      }
   // BetaK = RK12/RK2
   // p(k+1) = r(k+1)+betak*p(k)
      betak = rk12 / rk2;
      ae_v_move(&buf->xR[offspk1], 1, &buf->xR[offsrk1], 1, n);
      ae_v_addd(&buf->xR[offspk1], 1, &buf->xR[offspk], 1, n, betak);
   // r(k) := r(k+1)
   // x(k) := x(k+1)
   // p(k) := p(k+1)
      ae_v_move(&buf->xR[offsrk], 1, &buf->xR[offsrk1], 1, n);
      ae_v_move(&buf->xR[offsxk], 1, &buf->xR[offsxk1], 1, n);
      ae_v_move(&buf->xR[offspk], 1, &buf->xR[offspk1], 1, n);
      rk2 = rk12;
   }
// Calculate E2
   rmatrixmv(m, n, a, 0, 0, 0, buf, offsxk, buf, offstmp1);
   rmatrixmv(n, m, a, 0, 0, 1, buf, offstmp1, buf, offstmp2);
   ae_v_addd(&buf->xR[offstmp2], 1, &buf->xR[offsxk], 1, n, alpha);
   ae_v_move(&buf->xR[offsrk], 1, b->xR, 1, n);
   ae_v_sub(&buf->xR[offsrk], 1, &buf->xR[offstmp2], 1, n);
   v1 = ae_v_dotproduct(&buf->xR[offsrk], 1, &buf->xR[offsrk], 1, n);
   e2 = sqrt(v1);
// Output result (if it was improved)
   if (e2 < e1) {
      ae_v_move(x->xR, 1, &buf->xR[offsxk], 1, n);
   }
}

// Construction of linear conjugate gradient solver.
//
// State parameter passed using "var" semantics (i.e. previous state  is  NOT
// erased). When it is already initialized, we can reause prevously allocated
// memory.
//
// Inputs:
//     X       -   initial solution
//     B       -   right part
//     N       -   system size
//     State   -   structure; may be preallocated, if we want to reuse memory
//
// Outputs:
//     State   -   structure which is used by FBLSCGIteration() to store
//                 algorithm state between subsequent calls.
//
// NOTE: no error checking is done; caller must check all parameters, prevent
//       overflows, and so on.
// ALGLIB: Copyright 22.10.2009 by Sergey Bochkanov
void fblscgcreate(RVector *x, RVector *b, ae_int_t n, fblslincgstate *state) {
   if (state->b.cnt < n) {
      ae_vector_set_length(&state->b, n);
   }
   if (state->rk.cnt < n) {
      ae_vector_set_length(&state->rk, n);
   }
   if (state->rk1.cnt < n) {
      ae_vector_set_length(&state->rk1, n);
   }
   if (state->xk.cnt < n) {
      ae_vector_set_length(&state->xk, n);
   }
   if (state->xk1.cnt < n) {
      ae_vector_set_length(&state->xk1, n);
   }
   if (state->pk.cnt < n) {
      ae_vector_set_length(&state->pk, n);
   }
   if (state->pk1.cnt < n) {
      ae_vector_set_length(&state->pk1, n);
   }
   if (state->tmp2.cnt < n) {
      ae_vector_set_length(&state->tmp2, n);
   }
   if (state->x.cnt < n) {
      ae_vector_set_length(&state->x, n);
   }
   if (state->ax.cnt < n) {
      ae_vector_set_length(&state->ax, n);
   }
   state->n = n;
   ae_v_move(state->xk.xR, 1, x->xR, 1, n);
   ae_v_move(state->b.xR, 1, b->xR, 1, n);
   state->PQ = -1;
}

// Linear CG solver, function relying on reverse communication to calculate
// matrix-vector products.
//
// See comments for FBLSLinCGState structure for more info.
// ALGLIB: Copyright 22.10.2009 by Sergey Bochkanov
bool fblscgiteration(fblslincgstate *state) {
   AutoS ae_int_t n;
   AutoS ae_int_t k;
   AutoS double rk2;
   AutoS double rk12;
   AutoS double pap;
   AutoS double s;
   AutoS double betak;
   AutoS double v1;
   AutoS double v2;
// Manually threaded two-way signalling.
// Locals are set arbitrarily the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0; case 1: goto Resume1; case 2: goto Resume2;
      default: goto Exit;
   }
Spawn:
   k = -58;
   rk2 = -919;
   rk12 = -909;
   pap = 81;
   s = 255;
   betak = 74;
   v2 = 809;
// prepare locals
   n = state->n;
// Test for special case: B=0
   v1 = ae_v_dotproduct(state->b.xR, 1, state->b.xR, 1, n);
   if (v1 == 0.0) {
      for (k = 0; k < n; k++) {
         state->xk.xR[k] = 0.0;
      }
      goto Exit;
   }
// r(0) = b-A*x(0)
// RK2 = r(0)'*r(0)
   ae_v_move(state->x.xR, 1, state->xk.xR, 1, n);
   state->PQ = 0; goto Pause; Resume0:
   ae_v_move(state->rk.xR, 1, state->b.xR, 1, n);
   ae_v_sub(state->rk.xR, 1, state->ax.xR, 1, n);
   rk2 = ae_v_dotproduct(state->rk.xR, 1, state->rk.xR, 1, n);
   ae_v_move(state->pk.xR, 1, state->rk.xR, 1, n);
   state->e1 = sqrt(rk2);
// Cycle
   for (k = 0; k < n; k++) {
   // Calculate A*p(k) - store in State.Tmp2
   // and p(k)'*A*p(k)  - store in PAP
   //
   // If PAP=0, break (iteration is over)
      ae_v_move(state->x.xR, 1, state->pk.xR, 1, n);
      state->PQ = 1; goto Pause; Resume1:
      ae_v_move(state->tmp2.xR, 1, state->ax.xR, 1, n);
      pap = state->xax;
      if (!isfinite(pap)) {
         break;
      }
      if (pap <= 0.0) {
         break;
      }
   // S = (r(k)'*r(k))/(p(k)'*A*p(k))
      s = rk2 / pap;
   // x(k+1) = x(k) + S*p(k)
      ae_v_move(state->xk1.xR, 1, state->xk.xR, 1, n);
      ae_v_addd(state->xk1.xR, 1, state->pk.xR, 1, n, s);
   // r(k+1) = r(k) - S*A*p(k)
   // RK12 = r(k+1)'*r(k+1)
   //
   // Break if r(k+1) small enough (when compared to r(k))
      ae_v_move(state->rk1.xR, 1, state->rk.xR, 1, n);
      ae_v_subd(state->rk1.xR, 1, state->tmp2.xR, 1, n, s);
      rk12 = ae_v_dotproduct(state->rk1.xR, 1, state->rk1.xR, 1, n);
      if (sqrt(rk12) <= 100 * ae_machineepsilon * state->e1) {
      // X(k) = x(k+1) before exit -
      // - because we expect to find solution at x(k)
         ae_v_move(state->xk.xR, 1, state->xk1.xR, 1, n);
         break;
      }
   // BetaK = RK12/RK2
   // p(k+1) = r(k+1)+betak*p(k)
   //
   // NOTE: we expect that BetaK won't overflow because of
   // "sqrt(RK12) <= 100*MachineEpsilon*E1" test above.
      betak = rk12 / rk2;
      ae_v_move(state->pk1.xR, 1, state->rk1.xR, 1, n);
      ae_v_addd(state->pk1.xR, 1, state->pk.xR, 1, n, betak);
   // r(k) := r(k+1)
   // x(k) := x(k+1)
   // p(k) := p(k+1)
      ae_v_move(state->rk.xR, 1, state->rk1.xR, 1, n);
      ae_v_move(state->xk.xR, 1, state->xk1.xR, 1, n);
      ae_v_move(state->pk.xR, 1, state->pk1.xR, 1, n);
      rk2 = rk12;
   }
// Calculate E2
   ae_v_move(state->x.xR, 1, state->xk.xR, 1, n);
   state->PQ = 2; goto Pause; Resume2:
   ae_v_move(state->rk.xR, 1, state->b.xR, 1, n);
   ae_v_sub(state->rk.xR, 1, state->ax.xR, 1, n);
   v1 = ae_v_dotproduct(state->rk.xR, 1, state->rk.xR, 1, n);
   state->e2 = sqrt(v1);
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// Fast  least  squares  solver,  solves  well  conditioned  system   without
// performing  any  checks  for  degeneracy,  and using user-provided buffers
// (which are automatically reallocated if too small).
//
// This  function  is  intended  for solution of moderately sized systems. It
// uses factorization algorithms based on Level 2 BLAS  operations,  thus  it
// won't work efficiently on large scale systems.
//
// Inputs:
//     A       -   array[M,N], system matrix.
//                 Contents of A is destroyed during solution.
//     B       -   array[M], right part
//     M       -   number of equations
//     N       -   number of variables, N <= M
//     Tmp0, Tmp1, Tmp2-
//                 buffers; function automatically allocates them, if they are
//                 too  small. They can  be  reused  if  function  is   called
//                 several times.
//
// Outputs:
//     B       -   solution (first N components, next M-N are zero)
// ALGLIB: Copyright 20.01.2012 by Sergey Bochkanov
void fblssolvels(RMatrix *a, RVector *b, ae_int_t m, ae_int_t n, RVector *tmp0, RVector *tmp1, RVector *tmp2) {
   ae_int_t i;
   ae_int_t k;
   double v;
   ae_assert(n > 0, "FBLSSolveLS: N <= 0");
   ae_assert(m >= n, "FBLSSolveLS: M < N");
   ae_assert(a->rows >= m, "FBLSSolveLS: Rows(A)<M");
   ae_assert(a->cols >= n, "FBLSSolveLS: Cols(A)<N");
   ae_assert(b->cnt >= m, "FBLSSolveLS: Length(B)<M");
// Allocate temporaries
   vectorsetlengthatleast(tmp0, imax2(m, n) + 1);
   vectorsetlengthatleast(tmp1, imax2(m, n) + 1);
   vectorsetlengthatleast(tmp2, imin2(m, n));
// Call basecase QR
   ortfac_rmatrixqrbasecase(a, m, n, tmp0, tmp1, tmp2);
// Multiply B by Q'
   for (k = 0; k < n; k++) {
      for (i = 0; i < k; i++) {
         tmp0->xR[i] = 0.0;
      }
      ae_v_move(&tmp0->xR[k], 1, &a->xyR[k][k], a->stride, m - k);
      tmp0->xR[k] = 1.0;
      v = ae_v_dotproduct(&tmp0->xR[k], 1, &b->xR[k], 1, m - k);
      v *= tmp2->xR[k];
      ae_v_subd(&b->xR[k], 1, &tmp0->xR[k], 1, m - k, v);
   }
// Solve triangular system
   b->xR[n - 1] /= a->xyR[n - 1][n - 1];
   for (i = n - 2; i >= 0; i--) {
      v = ae_v_dotproduct(&a->xyR[i][i + 1], 1, &b->xR[i + 1], 1, n - i - 1);
      b->xR[i] = (b->xR[i] - v) / a->xyR[i][i];
   }
   for (i = n; i < m; i++) {
      b->xR[i] = 0.0;
   }
}

void fblslincgstate_init(void *_p, bool make_automatic) {
   fblslincgstate *p = (fblslincgstate *)_p;
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->ax, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rk, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rk1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xk, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xk1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->pk, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->pk1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->b, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tmp2, 0, DT_REAL, make_automatic);
}

void fblslincgstate_copy(void *_dst, void *_src, bool make_automatic) {
   fblslincgstate *dst = (fblslincgstate *)_dst;
   fblslincgstate *src = (fblslincgstate *)_src;
   dst->e1 = src->e1;
   dst->e2 = src->e2;
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   ae_vector_copy(&dst->ax, &src->ax, make_automatic);
   dst->xax = src->xax;
   dst->n = src->n;
   ae_vector_copy(&dst->rk, &src->rk, make_automatic);
   ae_vector_copy(&dst->rk1, &src->rk1, make_automatic);
   ae_vector_copy(&dst->xk, &src->xk, make_automatic);
   ae_vector_copy(&dst->xk1, &src->xk1, make_automatic);
   ae_vector_copy(&dst->pk, &src->pk, make_automatic);
   ae_vector_copy(&dst->pk1, &src->pk1, make_automatic);
   ae_vector_copy(&dst->b, &src->b, make_automatic);
   dst->PQ = src->PQ;
   ae_vector_copy(&dst->tmp2, &src->tmp2, make_automatic);
}

void fblslincgstate_free(void *_p, bool make_automatic) {
   fblslincgstate *p = (fblslincgstate *)_p;
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->ax, make_automatic);
   ae_vector_free(&p->rk, make_automatic);
   ae_vector_free(&p->rk1, make_automatic);
   ae_vector_free(&p->xk, make_automatic);
   ae_vector_free(&p->xk1, make_automatic);
   ae_vector_free(&p->pk, make_automatic);
   ae_vector_free(&p->pk1, make_automatic);
   ae_vector_free(&p->b, make_automatic);
   ae_vector_free(&p->tmp2, make_automatic);
}
} // end of namespace alglib_impl

// === BDSVD Package ===
// Depends on: (AlgLibInternal) ROTATIONS
// Depends on: (AlgLibMisc) HQRND
// Depends on: ABLAS
namespace alglib_impl {
static double bdsvd_extsignbdsqr(double a, double b) {
   double result;
   if (b >= 0.0) {
      result = fabs(a);
   } else {
      result = -fabs(a);
   }
   return result;
}

static void bdsvd_svdv2x2(double f, double g, double h, double *ssmin, double *ssmax, double *snr, double *csr, double *snl, double *csl) {
   bool gasmal;
   bool swp;
   ae_int_t pmax;
   double a;
   double clt;
   double crt;
   double d;
   double fa;
   double ft;
   double ga;
   double gt;
   double ha;
   double ht;
   double l;
   double m;
   double mm;
   double r;
   double s;
   double slt;
   double srt;
   double t;
   double tsign;
   double tt;
   double v;
   *ssmin = 0;
   *ssmax = 0;
   *snr = 0;
   *csr = 0;
   *snl = 0;
   *csl = 0;
   ft = f;
   fa = fabs(ft);
   ht = h;
   ha = fabs(h);
// these initializers are not really necessary,
// but without them compiler complains about uninitialized locals
   clt = 0.0;
   crt = 0.0;
   slt = 0.0;
   srt = 0.0;
   tsign = 0.0;
// PMAX points to the maximum absolute element of matrix
//  PMAX = 1 if F largest in absolute values
//  PMAX = 2 if G largest in absolute values
//  PMAX = 3 if H largest in absolute values
   pmax = 1;
   swp = ha > fa;
   if (swp) {
   // Now FA .ge. HA
      pmax = 3;
      swapr(&ft, &ht);
      swapr(&fa, &ha);
   }
   gt = g;
   ga = fabs(gt);
   if (ga == 0.0) {
   // Diagonal matrix
      *ssmin = ha;
      *ssmax = fa;
      clt = 1.0;
      crt = 1.0;
      slt = 0.0;
      srt = 0.0;
   } else {
      gasmal = true;
      if (ga > fa) {
         pmax = 2;
         if (fa / ga < ae_machineepsilon) {
         // Case of very large GA
            gasmal = false;
            *ssmax = ga;
            if (ha > 1.0) {
               v = ga / ha;
               *ssmin = fa / v;
            } else {
               v = fa / ga;
               *ssmin = v * ha;
            }
            clt = 1.0;
            slt = ht / gt;
            srt = 1.0;
            crt = ft / gt;
         }
      }
      if (gasmal) {
      // Normal case
         d = fa - ha;
         if (d == fa) {
            l = 1.0;
         } else {
            l = d / fa;
         }
         m = gt / ft;
         t = 2 - l;
         mm = m * m;
         tt = t * t;
         s = sqrt(tt + mm);
         if (l == 0.0) {
            r = fabs(m);
         } else {
            r = sqrt(l * l + mm);
         }
         a = 0.5 * (s + r);
         *ssmin = ha / a;
         *ssmax = fa * a;
         if (mm == 0.0) {
         // Note that M is very tiny
            if (l == 0.0) {
               t = bdsvd_extsignbdsqr(2.0, ft) * bdsvd_extsignbdsqr(1.0, gt);
            } else {
               t = gt / bdsvd_extsignbdsqr(d, ft) + m / t;
            }
         } else {
            t = (m / (s + t) + m / (r + l)) * (1 + a);
         }
         l = sqrt(t * t + 4);
         crt = 2 / l;
         srt = t / l;
         clt = (crt + srt * m) / a;
         v = ht / ft;
         slt = v * srt / a;
      }
   }
   if (swp) {
      *csl = srt;
      *snl = crt;
      *csr = slt;
      *snr = clt;
   } else {
      *csl = clt;
      *snl = slt;
      *csr = crt;
      *snr = srt;
   }
// Correct signs of SSMAX and SSMIN
   if (pmax == 1) {
      tsign = bdsvd_extsignbdsqr(1.0, *csr) * bdsvd_extsignbdsqr(1.0, *csl) * bdsvd_extsignbdsqr(1.0, f);
   }
   if (pmax == 2) {
      tsign = bdsvd_extsignbdsqr(1.0, *snr) * bdsvd_extsignbdsqr(1.0, *csl) * bdsvd_extsignbdsqr(1.0, g);
   }
   if (pmax == 3) {
      tsign = bdsvd_extsignbdsqr(1.0, *snr) * bdsvd_extsignbdsqr(1.0, *snl) * bdsvd_extsignbdsqr(1.0, h);
   }
   *ssmax = bdsvd_extsignbdsqr(*ssmax, tsign);
   *ssmin = bdsvd_extsignbdsqr(*ssmin, tsign * bdsvd_extsignbdsqr(1.0, f) * bdsvd_extsignbdsqr(1.0, h));
}

static void bdsvd_svd2x2(double f, double g, double h, double *ssmin, double *ssmax) {
   double aas;
   double at;
   double au;
   double c;
   double fa;
   double fhmn;
   double fhmx;
   double ga;
   double ha;
   *ssmin = 0;
   *ssmax = 0;
   fa = fabs(f);
   ga = fabs(g);
   ha = fabs(h);
   fhmn = rmin2(fa, ha);
   fhmx = rmax2(fa, ha);
   if (fhmn == 0.0) {
      *ssmin = 0.0;
      if (fhmx == 0.0) {
         *ssmax = ga;
      } else {
         *ssmax = rmax2(fhmx, ga) * sqrt(1 + ae_sqr(rmin2(fhmx, ga) / rmax2(fhmx, ga)));
      }
   } else {
      if (ga < fhmx) {
         aas = 1 + fhmn / fhmx;
         at = (fhmx - fhmn) / fhmx;
         au = ae_sqr(ga / fhmx);
         c = 2 / (sqrt(aas * aas + au) + sqrt(at * at + au));
         *ssmin = fhmn * c;
         *ssmax = fhmx / c;
      } else {
         au = fhmx / ga;
         if (au == 0.0) {
         // Avoid possible harmful underflow if exponent range
         // asymmetric (true SSMIN may not underflow even if
         // AU underflows)
            *ssmin = fhmn * fhmx / ga;
            *ssmax = ga;
         } else {
            aas = 1 + fhmn / fhmx;
            at = (fhmx - fhmn) / fhmx;
            c = 1 / (sqrt(1 + ae_sqr(aas * au)) + sqrt(1 + ae_sqr(at * au)));
            *ssmin = fhmn * c * au;
            *ssmin += *ssmin;
            *ssmax = ga / (c + c);
         }
      }
   }
}

// Internal working subroutine for bidiagonal decomposition
static bool bdsvd_bidiagonalsvddecompositioninternal(RVector *d, RVector *e, ae_int_t n, bool isupper, bool isfractionalaccuracyrequired, RMatrix *uu, ae_int_t ustart, ae_int_t nru, RMatrix *c, ae_int_t cstart, ae_int_t ncc, RMatrix *vt, ae_int_t vstart, ae_int_t ncvt) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t idir;
   ae_int_t isub;
   ae_int_t iter;
   ae_int_t j;
   ae_int_t ll;
   ae_int_t lll;
   ae_int_t m;
   ae_int_t maxit;
   ae_int_t oldll;
   ae_int_t oldm;
   double abse;
   double abss;
   double cosl;
   double cosr;
   double cs;
   double eps;
   double f;
   double g;
   double h;
   double mu;
   double oldcs;
   double oldsn;
   double r;
   double shift;
   double sigmn;
   double sigmx;
   double sinl;
   double sinr;
   double sll;
   double smax;
   double smin;
   double sminl;
   double sminoa;
   double sn;
   double thresh;
   double tol;
   double tolmul;
   double unfl;
   ae_int_t maxitr;
   bool matrixsplitflag;
   bool iterflag;
   bool fwddir;
   double tmp;
   ae_int_t mm1;
   ae_int_t mm0;
   bool bchangedir;
   ae_int_t uend;
   ae_int_t cend;
   ae_int_t vend;
   bool result;
   ae_frame_make(&_frame_block);
   DupVector(e);
   NewVector(work0, 0, DT_REAL);
   NewVector(work1, 0, DT_REAL);
   NewVector(work2, 0, DT_REAL);
   NewVector(work3, 0, DT_REAL);
   NewVector(utemp, 0, DT_REAL);
   NewVector(vttemp, 0, DT_REAL);
   NewVector(ctemp, 0, DT_REAL);
   NewVector(etemp, 0, DT_REAL);
   NewMatrix(ut, 0, 0, DT_REAL);
   result = true;
   if (n == 0) {
      ae_frame_leave();
      return result;
   }
   if (n == 1) {
      if (d->xR[1] < 0.0) {
         d->xR[1] = -d->xR[1];
         if (ncvt > 0) {
            ae_v_muld(&vt->xyR[vstart][vstart], 1, ncvt, -1);
         }
      }
      ae_frame_leave();
      return result;
   }
// these initializers are not really necessary,
// but without them compiler complains about uninitialized locals
   ll = 0;
   oldsn = 0.0;
// init
   ae_vector_set_length(&work0, n);
   ae_vector_set_length(&work1, n);
   ae_vector_set_length(&work2, n);
   ae_vector_set_length(&work3, n);
   uend = ustart + imax2(nru - 1, 0);
   vend = vstart + imax2(ncvt - 1, 0);
   cend = cstart + imax2(ncc - 1, 0);
   ae_vector_set_length(&utemp, uend + 1);
   ae_vector_set_length(&vttemp, vend + 1);
   ae_vector_set_length(&ctemp, cend + 1);
   maxitr = 12;
   fwddir = true;
   if (nru > 0) {
      ae_matrix_set_length(&ut, ustart + n, ustart + nru);
      rmatrixtranspose(nru, n, uu, ustart, ustart, &ut, ustart, ustart);
   }
// resize E from N-1 to N
   ae_vector_set_length(&etemp, n + 1);
   for (i = 1; i < n; i++) {
      etemp.xR[i] = e->xR[i];
   }
   ae_vector_set_length(e, n + 1);
   for (i = 1; i < n; i++) {
      e->xR[i] = etemp.xR[i];
   }
   e->xR[n] = 0.0;
   idir = 0;
// Get machine constants
   eps = ae_machineepsilon;
   unfl = ae_minrealnumber;
// If matrix lower bidiagonal, rotate to be upper bidiagonal
// by applying Givens rotations on the left
   if (!isupper) {
      for (i = 1; i < n; i++) {
         generaterotation(d->xR[i], e->xR[i], &cs, &sn, &r);
         d->xR[i] = r;
         e->xR[i] = sn * d->xR[i + 1];
         d->xR[i + 1] *= cs;
         work0.xR[i] = cs;
         work1.xR[i] = sn;
      }
   // Update singular vectors if desired
      if (nru > 0) {
         applyrotationsfromtheleft(fwddir, 1 + ustart - 1, n + ustart - 1, ustart, uend, &work0, &work1, &ut, &utemp);
      }
      if (ncc > 0) {
         applyrotationsfromtheleft(fwddir, 1 + cstart - 1, n + cstart - 1, cstart, cend, &work0, &work1, c, &ctemp);
      }
   }
// Compute singular values to relative accuracy TOL
// (By setting TOL to be negative, algorithm will compute
// singular values to absolute accuracy ABS(TOL)*norm(input matrix))
   tolmul = rmax2(10.0, rmin2(100.0, pow(eps, -0.125)));
   tol = tolmul * eps;
// Compute approximate maximum, minimum singular values
   smax = 0.0;
   for (i = 1; i <= n; i++) {
      smax = rmax2(smax, fabs(d->xR[i]));
   }
   for (i = 1; i < n; i++) {
      smax = rmax2(smax, fabs(e->xR[i]));
   }
   sminl = 0.0;
   if (tol >= 0.0) {
   // Relative accuracy desired
      sminoa = fabs(d->xR[1]);
      if (sminoa != 0.0) {
         mu = sminoa;
         for (i = 2; i <= n; i++) {
            mu = fabs(d->xR[i]) * (mu / (mu + fabs(e->xR[i - 1])));
            sminoa = rmin2(sminoa, mu);
            if (sminoa == 0.0) {
               break;
            }
         }
      }
      sminoa /= sqrt((double)n);
      thresh = rmax2(tol * sminoa, maxitr * n * n * unfl);
   } else {
   // Absolute accuracy desired
      thresh = rmax2(fabs(tol) * smax, maxitr * n * n * unfl);
   }
// Prepare for main iteration loop for the singular values
// (MAXIT is the maximum number of passes through the inner
// loop permitted before nonconvergence signalled.)
   maxit = maxitr * n * n;
   iter = 0;
   oldll = -1;
   oldm = -1;
// M points to last element of unconverged part of matrix
   m = n;
// Begin main iteration loop
   while (true) {
   // Check for convergence or exceeding iteration count
      if (m <= 1) {
         break;
      }
      if (iter > maxit) {
         result = false;
         ae_frame_leave();
         return result;
      }
   // Find diagonal block of matrix to work on
      if (tol < 0.0 && SmallAtR(d->xR[m], thresh)) {
         d->xR[m] = 0.0;
      }
      smax = fabs(d->xR[m]);
      smin = smax;
      matrixsplitflag = false;
      for (lll = 1; lll < m; lll++) {
         ll = m - lll;
         abss = fabs(d->xR[ll]);
         abse = fabs(e->xR[ll]);
         if (tol < 0.0 && abss <= thresh) {
            d->xR[ll] = 0.0;
         }
         if (abse <= thresh) {
            matrixsplitflag = true;
            break;
         }
         smin = rmin2(smin, abss);
         smax = rmax2(smax, rmax2(abss, abse));
      }
      if (!matrixsplitflag) {
         ll = 0;
      } else {
      // Matrix splits since E(LL) = 0
         e->xR[ll] = 0.0;
         if (ll == m - 1) {
         // Convergence of bottom singular value, return to top of loop
            m--;
            continue;
         }
      }
      ll++;
   // E(LL) through E(M-1) are nonzero, E(LL-1) is zero
      if (ll == m - 1) {
      // 2 by 2 block, handle separately
         bdsvd_svdv2x2(d->xR[m - 1], e->xR[m - 1], d->xR[m], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
         d->xR[m - 1] = sigmx;
         e->xR[m - 1] = 0.0;
         d->xR[m] = sigmn;
      // Compute singular vectors, if desired
         if (ncvt > 0) {
            mm0 = m + (vstart - 1);
            mm1 = m - 1 + (vstart - 1);
            ae_v_moved(&vttemp.xR[vstart], 1, &vt->xyR[mm1][vstart], 1, vend - vstart + 1, cosr);
            ae_v_addd(&vttemp.xR[vstart], 1, &vt->xyR[mm0][vstart], 1, vend - vstart + 1, sinr);
            ae_v_muld(&vt->xyR[mm0][vstart], 1, vend - vstart + 1, cosr);
            ae_v_subd(&vt->xyR[mm0][vstart], 1, &vt->xyR[mm1][vstart], 1, vend - vstart + 1, sinr);
            ae_v_move(&vt->xyR[mm1][vstart], 1, &vttemp.xR[vstart], 1, vend - vstart + 1);
         }
         if (nru > 0) {
            mm0 = m + ustart - 1;
            mm1 = m - 1 + ustart - 1;
            ae_v_moved(&utemp.xR[ustart], 1, &ut.xyR[mm1][ustart], 1, uend - ustart + 1, cosl);
            ae_v_addd(&utemp.xR[ustart], 1, &ut.xyR[mm0][ustart], 1, uend - ustart + 1, sinl);
            ae_v_muld(&ut.xyR[mm0][ustart], 1, uend - ustart + 1, cosl);
            ae_v_subd(&ut.xyR[mm0][ustart], 1, &ut.xyR[mm1][ustart], 1, uend - ustart + 1, sinl);
            ae_v_move(&ut.xyR[mm1][ustart], 1, &utemp.xR[ustart], 1, uend - ustart + 1);
         }
         if (ncc > 0) {
            mm0 = m + cstart - 1;
            mm1 = m - 1 + cstart - 1;
            ae_v_moved(&ctemp.xR[cstart], 1, &c->xyR[mm1][cstart], 1, cend - cstart + 1, cosl);
            ae_v_addd(&ctemp.xR[cstart], 1, &c->xyR[mm0][cstart], 1, cend - cstart + 1, sinl);
            ae_v_muld(&c->xyR[mm0][cstart], 1, cend - cstart + 1, cosl);
            ae_v_subd(&c->xyR[mm0][cstart], 1, &c->xyR[mm1][cstart], 1, cend - cstart + 1, sinl);
            ae_v_move(&c->xyR[mm1][cstart], 1, &ctemp.xR[cstart], 1, cend - cstart + 1);
         }
         m -= 2;
         continue;
      }
   // If working on new submatrix, choose shift direction
   // (from larger end diagonal element towards smaller)
   //
   // Previously was
   //     "if (LL>OLDM) or (M<OLDLL) then"
   // fixed thanks to Michael Rolle < m@rolle.name >
   // Very strange that LAPACK still contains it.
      bchangedir = false;
      if (idir == 1 && SmallR(d->xR[ll], 1.0E-3 * fabs(d->xR[m]))) {
         bchangedir = true;
      }
      if (idir == 2 && SmallR(d->xR[m], 1.0E-3 * fabs(d->xR[ll]))) {
         bchangedir = true;
      }
      if (ll != oldll || m != oldm || bchangedir) {
         if (fabs(d->xR[ll]) >= fabs(d->xR[m])) {
         // Chase bulge from top (big end) to bottom (small end)
            idir = 1;
         } else {
         // Chase bulge from bottom (big end) to top (small end)
            idir = 2;
         }
      }
   // Apply convergence tests
      if (idir == 1) {
      // Run convergence test in forward direction
      // First apply standard test to bottom of matrix
         if (SmallAtR(e->xR[m - 1], fabs(tol) * fabs(d->xR[m])) || tol < 0.0 && SmallAtR(e->xR[m - 1], thresh)) {
            e->xR[m - 1] = 0.0;
            continue;
         }
         if (tol >= 0.0) {
         // If relative accuracy desired,
         // apply convergence criterion forward
            mu = fabs(d->xR[ll]);
            sminl = mu;
            iterflag = false;
            for (lll = ll; lll < m; lll++) {
               if (SmallAtR(e->xR[lll], tol * mu)) {
                  e->xR[lll] = 0.0;
                  iterflag = true;
                  break;
               }
               mu = fabs(d->xR[lll + 1]) * (mu / (mu + fabs(e->xR[lll])));
               sminl = rmin2(sminl, mu);
            }
            if (iterflag) {
               continue;
            }
         }
      } else {
      // Run convergence test in backward direction
      // First apply standard test to top of matrix
         if (SmallAtR(e->xR[ll], fabs(tol) * fabs(d->xR[ll])) || tol < 0.0 && SmallAtR(e->xR[ll], thresh)) {
            e->xR[ll] = 0.0;
            continue;
         }
         if (tol >= 0.0) {
         // If relative accuracy desired,
         // apply convergence criterion backward
            mu = fabs(d->xR[m]);
            sminl = mu;
            iterflag = false;
            for (lll = m - 1; lll >= ll; lll--) {
               if (SmallAtR(e->xR[lll], tol * mu)) {
                  e->xR[lll] = 0.0;
                  iterflag = true;
                  break;
               }
               mu = fabs(d->xR[lll]) * (mu / (mu + fabs(e->xR[lll])));
               sminl = rmin2(sminl, mu);
            }
            if (iterflag) {
               continue;
            }
         }
      }
      oldll = ll;
      oldm = m;
   // Compute shift.  First, test if shifting would ruin relative
   // accuracy, and if so set the shift to zero.
      if (tol >= 0.0 && n * tol * (sminl / smax) <= rmax2(eps, 0.01 * tol)) {
      // Use a zero shift to avoid loss of relative accuracy
         shift = 0.0;
      } else {
      // Compute the shift from 2-by-2 block at end of matrix
         if (idir == 1) {
            sll = fabs(d->xR[ll]);
            bdsvd_svd2x2(d->xR[m - 1], e->xR[m - 1], d->xR[m], &shift, &r);
         } else {
            sll = fabs(d->xR[m]);
            bdsvd_svd2x2(d->xR[ll], e->xR[ll], d->xR[ll + 1], &shift, &r);
         }
      // Test if shift negligible, and if so set to zero
         if (sll > 0.0) {
            if (ae_sqr(shift / sll) < eps) {
               shift = 0.0;
            }
         }
      }
   // Increment iteration count
      iter += m - ll;
   // If SHIFT = 0, do simplified QR iteration
      if (shift == 0.0) {
         if (idir == 1) {
         // Chase bulge from top to bottom
         // Save cosines and sines for later singular vector updates
            cs = 1.0;
            oldcs = 1.0;
            for (i = ll; i < m; i++) {
               generaterotation(d->xR[i] * cs, e->xR[i], &cs, &sn, &r);
               if (i > ll) {
                  e->xR[i - 1] = oldsn * r;
               }
               generaterotation(oldcs * r, d->xR[i + 1] * sn, &oldcs, &oldsn, &tmp);
               d->xR[i] = tmp;
               work0.xR[i - ll + 1] = cs;
               work1.xR[i - ll + 1] = sn;
               work2.xR[i - ll + 1] = oldcs;
               work3.xR[i - ll + 1] = oldsn;
            }
            h = d->xR[m] * cs;
            d->xR[m] = h * oldcs;
            e->xR[m - 1] = h * oldsn;
         // Update singular vectors
            if (ncvt > 0) {
               applyrotationsfromtheleft(fwddir, ll + vstart - 1, m + vstart - 1, vstart, vend, &work0, &work1, vt, &vttemp);
            }
            if (nru > 0) {
               applyrotationsfromtheleft(fwddir, ll + ustart - 1, m + ustart - 1, ustart, uend, &work2, &work3, &ut, &utemp);
            }
            if (ncc > 0) {
               applyrotationsfromtheleft(fwddir, ll + cstart - 1, m + cstart - 1, cstart, cend, &work2, &work3, c, &ctemp);
            }
         // Test convergence
            if (SmallAtR(e->xR[m - 1], thresh)) {
               e->xR[m - 1] = 0.0;
            }
         } else {
         // Chase bulge from bottom to top
         // Save cosines and sines for later singular vector updates
            cs = 1.0;
            oldcs = 1.0;
            for (i = m; i >= ll + 1; i--) {
               generaterotation(d->xR[i] * cs, e->xR[i - 1], &cs, &sn, &r);
               if (i < m) {
                  e->xR[i] = oldsn * r;
               }
               generaterotation(oldcs * r, d->xR[i - 1] * sn, &oldcs, &oldsn, &tmp);
               d->xR[i] = tmp;
               work0.xR[i - ll] = cs;
               work1.xR[i - ll] = -sn;
               work2.xR[i - ll] = oldcs;
               work3.xR[i - ll] = -oldsn;
            }
            h = d->xR[ll] * cs;
            d->xR[ll] = h * oldcs;
            e->xR[ll] = h * oldsn;
         // Update singular vectors
            if (ncvt > 0) {
               applyrotationsfromtheleft(!fwddir, ll + vstart - 1, m + vstart - 1, vstart, vend, &work2, &work3, vt, &vttemp);
            }
            if (nru > 0) {
               applyrotationsfromtheleft(!fwddir, ll + ustart - 1, m + ustart - 1, ustart, uend, &work0, &work1, &ut, &utemp);
            }
            if (ncc > 0) {
               applyrotationsfromtheleft(!fwddir, ll + cstart - 1, m + cstart - 1, cstart, cend, &work0, &work1, c, &ctemp);
            }
         // Test convergence
            if (SmallAtR(e->xR[ll], thresh)) {
               e->xR[ll] = 0.0;
            }
         }
      } else {
      // Use nonzero shift
         if (idir == 1) {
         // Chase bulge from top to bottom
         // Save cosines and sines for later singular vector updates
            f = (fabs(d->xR[ll]) - shift) * (bdsvd_extsignbdsqr(1.0, d->xR[ll]) + shift / d->xR[ll]);
            g = e->xR[ll];
            for (i = ll; i < m; i++) {
               generaterotation(f, g, &cosr, &sinr, &r);
               if (i > ll) {
                  e->xR[i - 1] = r;
               }
               f = cosr * d->xR[i] + sinr * e->xR[i];
               e->xR[i] = cosr * e->xR[i] - sinr * d->xR[i];
               g = sinr * d->xR[i + 1];
               d->xR[i + 1] *= cosr;
               generaterotation(f, g, &cosl, &sinl, &r);
               d->xR[i] = r;
               f = cosl * e->xR[i] + sinl * d->xR[i + 1];
               d->xR[i + 1] = cosl * d->xR[i + 1] - sinl * e->xR[i];
               if (i < m - 1) {
                  g = sinl * e->xR[i + 1];
                  e->xR[i + 1] *= cosl;
               }
               work0.xR[i - ll + 1] = cosr;
               work1.xR[i - ll + 1] = sinr;
               work2.xR[i - ll + 1] = cosl;
               work3.xR[i - ll + 1] = sinl;
            }
            e->xR[m - 1] = f;
         // Update singular vectors
            if (ncvt > 0) {
               applyrotationsfromtheleft(fwddir, ll + vstart - 1, m + vstart - 1, vstart, vend, &work0, &work1, vt, &vttemp);
            }
            if (nru > 0) {
               applyrotationsfromtheleft(fwddir, ll + ustart - 1, m + ustart - 1, ustart, uend, &work2, &work3, &ut, &utemp);
            }
            if (ncc > 0) {
               applyrotationsfromtheleft(fwddir, ll + cstart - 1, m + cstart - 1, cstart, cend, &work2, &work3, c, &ctemp);
            }
         // Test convergence
            if (SmallAtR(e->xR[m - 1], thresh)) {
               e->xR[m - 1] = 0.0;
            }
         } else {
         // Chase bulge from bottom to top
         // Save cosines and sines for later singular vector updates
            f = (fabs(d->xR[m]) - shift) * (bdsvd_extsignbdsqr(1.0, d->xR[m]) + shift / d->xR[m]);
            g = e->xR[m - 1];
            for (i = m; i >= ll + 1; i--) {
               generaterotation(f, g, &cosr, &sinr, &r);
               if (i < m) {
                  e->xR[i] = r;
               }
               f = cosr * d->xR[i] + sinr * e->xR[i - 1];
               e->xR[i - 1] = cosr * e->xR[i - 1] - sinr * d->xR[i];
               g = sinr * d->xR[i - 1];
               d->xR[i - 1] *= cosr;
               generaterotation(f, g, &cosl, &sinl, &r);
               d->xR[i] = r;
               f = cosl * e->xR[i - 1] + sinl * d->xR[i - 1];
               d->xR[i - 1] = cosl * d->xR[i - 1] - sinl * e->xR[i - 1];
               if (i > ll + 1) {
                  g = sinl * e->xR[i - 2];
                  e->xR[i - 2] *= cosl;
               }
               work0.xR[i - ll] = cosr;
               work1.xR[i - ll] = -sinr;
               work2.xR[i - ll] = cosl;
               work3.xR[i - ll] = -sinl;
            }
            e->xR[ll] = f;
         // Test convergence
            if (SmallAtR(e->xR[ll], thresh)) {
               e->xR[ll] = 0.0;
            }
         // Update singular vectors if desired
            if (ncvt > 0) {
               applyrotationsfromtheleft(!fwddir, ll + vstart - 1, m + vstart - 1, vstart, vend, &work2, &work3, vt, &vttemp);
            }
            if (nru > 0) {
               applyrotationsfromtheleft(!fwddir, ll + ustart - 1, m + ustart - 1, ustart, uend, &work0, &work1, &ut, &utemp);
            }
            if (ncc > 0) {
               applyrotationsfromtheleft(!fwddir, ll + cstart - 1, m + cstart - 1, cstart, cend, &work0, &work1, c, &ctemp);
            }
         }
      }
   // QR iteration finished, go back and check convergence
      continue;
   }
// All singular values converged, so make them positive
   for (i = 1; i <= n; i++) {
      if (d->xR[i] < 0.0) {
         d->xR[i] = -d->xR[i];
      // Change sign of singular vectors, if desired
         if (ncvt > 0) {
            ae_v_muld(&vt->xyR[i + vstart - 1][vstart], 1, vend - vstart + 1, -1);
         }
      }
   }
// Sort the singular values into decreasing order (insertion sort on
// singular values, but only one transposition per singular vector)
   for (i = 1; i < n; i++) {
   // Scan for smallest D(I)
      isub = 1;
      smin = d->xR[1];
      for (j = 2; j <= n + 1 - i; j++) {
         if (d->xR[j] <= smin) {
            isub = j;
            smin = d->xR[j];
         }
      }
      if (isub != n + 1 - i) {
      // Swap singular values and vectors
         d->xR[isub] = d->xR[n + 1 - i];
         d->xR[n + 1 - i] = smin;
         if (ncvt > 0) {
            j = n + 1 - i;
            ae_v_move(&vttemp.xR[vstart], 1, &vt->xyR[isub + vstart - 1][vstart], 1, vend - vstart + 1);
            ae_v_move(&vt->xyR[isub + vstart - 1][vstart], 1, &vt->xyR[j + vstart - 1][vstart], 1, vend - vstart + 1);
            ae_v_move(&vt->xyR[j + vstart - 1][vstart], 1, &vttemp.xR[vstart], 1, vend - vstart + 1);
         }
         if (nru > 0) {
            j = n + 1 - i;
            ae_v_move(&utemp.xR[ustart], 1, &ut.xyR[isub + ustart - 1][ustart], 1, uend - ustart + 1);
            ae_v_move(&ut.xyR[isub + ustart - 1][ustart], 1, &ut.xyR[j + ustart - 1][ustart], 1, uend - ustart + 1);
            ae_v_move(&ut.xyR[j + ustart - 1][ustart], 1, &utemp.xR[ustart], 1, uend - ustart + 1);
         }
         if (ncc > 0) {
            j = n + 1 - i;
            ae_v_move(&ctemp.xR[cstart], 1, &c->xyR[isub + cstart - 1][cstart], 1, cend - cstart + 1);
            ae_v_move(&c->xyR[isub + cstart - 1][cstart], 1, &c->xyR[j + cstart - 1][cstart], 1, cend - cstart + 1);
            ae_v_move(&c->xyR[j + cstart - 1][cstart], 1, &ctemp.xR[cstart], 1, cend - cstart + 1);
         }
      }
   }
// Copy U back from temporary storage
   if (nru > 0) {
      rmatrixtranspose(n, nru, &ut, ustart, ustart, uu, ustart, ustart);
   }
   ae_frame_leave();
   return result;
}

// Singular value decomposition of a bidiagonal matrix (extended algorithm)
//
// The algorithm performs the singular value decomposition  of  a  bidiagonal
// matrix B (upper or lower) representing it as B = Q*S*P^T, where Q and  P -
// orthogonal matrices, S - diagonal matrix with non-negative elements on the
// main diagonal, in descending order.
//
// The  algorithm  finds  singular  values.  In  addition,  the algorithm can
// calculate  matrices  Q  and P (more precisely, not the matrices, but their
// product  with  given  matrices U and VT - U*Q and (P^T)*VT)).  Of  course,
// matrices U and VT can be of any type, including identity. Furthermore, the
// algorithm can calculate Q'*C (this product is calculated more  effectively
// than U*Q,  because  this calculation operates with rows instead  of matrix
// columns).
//
// The feature of the algorithm is its ability to find  all  singular  values
// including those which are arbitrarily close to 0  with  relative  accuracy
// close to  machine precision. If the parameter IsFractionalAccuracyRequired
// is set to True, all singular values will have high relative accuracy close
// to machine precision. If the parameter is set to False, only  the  biggest
// singular value will have relative accuracy  close  to  machine  precision.
// The absolute error of other singular values is equal to the absolute error
// of the biggest singular value.
//
// Inputs:
//     D       -   main diagonal of matrix B.
//                 Array whose index ranges within [0..N-1].
//     E       -   superdiagonal (or subdiagonal) of matrix B.
//                 Array whose index ranges within [0..N-2].
//     N       -   size of matrix B.
//     IsUpper -   True, if the matrix is upper bidiagonal.
//     IsFractionalAccuracyRequired -
//                 THIS PARAMETER IS IGNORED SINCE ALGLIB 3.5.0
//                 SINGULAR VALUES ARE ALWAYS SEARCHED WITH HIGH ACCURACY.
//     U       -   matrix to be multiplied by Q.
//                 Array whose indexes range within [0..NRU-1, 0..N-1].
//                 The matrix can be bigger, in that case only the  submatrix
//                 [0..NRU-1, 0..N-1] will be multiplied by Q.
//     NRU     -   number of rows in matrix U.
//     C       -   matrix to be multiplied by Q'.
//                 Array whose indexes range within [0..N-1, 0..NCC-1].
//                 The matrix can be bigger, in that case only the  submatrix
//                 [0..N-1, 0..NCC-1] will be multiplied by Q'.
//     NCC     -   number of columns in matrix C.
//     VT      -   matrix to be multiplied by P^T.
//                 Array whose indexes range within [0..N-1, 0..NCVT-1].
//                 The matrix can be bigger, in that case only the  submatrix
//                 [0..N-1, 0..NCVT-1] will be multiplied by P^T.
//     NCVT    -   number of columns in matrix VT.
//
// Outputs:
//     D       -   singular values of matrix B in descending order.
//     U       -   if NRU > 0, contains matrix U*Q.
//     VT      -   if NCVT > 0, contains matrix (P^T)*VT.
//     C       -   if NCC > 0, contains matrix Q'*C.
//
// Result:
//     True, if the algorithm has converged.
//     False, if the algorithm hasn't converged (rare case).
//
// NOTE: multiplication U*Q is performed by means of transposition to internal
//       buffer, multiplication and backward transposition. It helps to avoid
//       costly columnwise operations and speed-up algorithm.
//
// Additional information:
//     The type of convergence is controlled by the internal  parameter  TOL.
//     If the parameter is greater than 0, the singular values will have
//     relative accuracy TOL. If TOL < 0, the singular values will have
//     absolute accuracy ABS(TOL)*norm(B).
//     By default, |TOL| falls within the range of 10*Epsilon and 100*Epsilon,
//     where Epsilon is the machine precision. It is not  recommended  to  use
//     TOL less than 10*Epsilon since this will  considerably  slow  down  the
//     algorithm and may not lead to error decreasing.
//
// History:
//     * 31 March, 2007.
//         changed MAXITR from 6 to 12.
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1999.
// API: bool rmatrixbdsvd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const bool isupper, const bool isfractionalaccuracyrequired, real_2d_array &u, const ae_int_t nru, real_2d_array &c, const ae_int_t ncc, real_2d_array &vt, const ae_int_t ncvt);
bool rmatrixbdsvd(RVector *d, RVector *e, ae_int_t n, bool isupper, bool isfractionalaccuracyrequired, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt) {
   ae_frame _frame_block;
   ae_int_t i;
   bool result;
   ae_frame_make(&_frame_block);
   DupVector(e);
   NewVector(en, 0, DT_REAL);
   NewVector(d1, 0, DT_REAL);
   NewVector(e1, 0, DT_REAL);
   result = false;
// Try to use MKL
   ae_vector_set_length(&en, n);
   for (i = 0; i < n - 1; i++) {
      en.xR[i] = e->xR[i];
   }
   en.xR[n - 1] = 0.0;
   if (rmatrixbdsvdmkl(d, &en, n, isupper, u, nru, c, ncc, vt, ncvt, &result)) {
      ae_frame_leave();
      return result;
   }
// Use ALGLIB code
   ae_vector_set_length(&d1, n + 1);
   ae_v_move(&d1.xR[1], 1, d->xR, 1, n);
   if (n > 1) {
      ae_vector_set_length(&e1, n);
      ae_v_move(&e1.xR[1], 1, e->xR, 1, n - 1);
   }
   result = bdsvd_bidiagonalsvddecompositioninternal(&d1, &e1, n, isupper, isfractionalaccuracyrequired, u, 0, nru, c, 0, ncc, vt, 0, ncvt);
   ae_v_move(d->xR, 1, &d1.xR[1], 1, n);
   ae_frame_leave();
   return result;
}

bool bidiagonalsvddecomposition(RVector *d, RVector *e, ae_int_t n, bool isupper, bool isfractionalaccuracyrequired, RMatrix *u, ae_int_t nru, RMatrix *c, ae_int_t ncc, RMatrix *vt, ae_int_t ncvt) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   DupVector(e);
   result = bdsvd_bidiagonalsvddecompositioninternal(d, e, n, isupper, isfractionalaccuracyrequired, u, 1, nru, c, 1, ncc, vt, 1, ncvt);
   ae_frame_leave();
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
bool rmatrixbdsvd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const bool isupper, const bool isfractionalaccuracyrequired, real_2d_array &u, const ae_int_t nru, real_2d_array &c, const ae_int_t ncc, real_2d_array &vt, const ae_int_t ncvt) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::rmatrixbdsvd(ConstT(ae_vector, d), ConstT(ae_vector, e), n, isupper, isfractionalaccuracyrequired, ConstT(ae_matrix, u), nru, ConstT(ae_matrix, c), ncc, ConstT(ae_matrix, vt), ncvt);
   alglib_impl::ae_state_clear();
   return Ok;
}
} // end of namespace alglib

// === SVD Package ===
// Depends on: (AlgLibInternal) BLAS
// Depends on: ORTFAC, BDSVD
namespace alglib_impl {
// Singular value decomposition of a rectangular matrix.
//
// The algorithm calculates the singular value decomposition of a matrix of
// size MxN: A = U * S * V^T
//
// The algorithm finds the singular values and, optionally, matrices U and V^T.
// The algorithm can find both first min(M,N) columns of matrix U and rows of
// matrix V^T (singular vectors), and matrices U and V^T wholly (of sizes MxM
// and NxN respectively).
//
// Take into account that the subroutine does not return matrix V but V^T.
//
// Inputs:
//     A           -   matrix to be decomposed.
//                     Array whose indexes range within [0..M-1, 0..N-1].
//     M           -   number of rows in matrix A.
//     N           -   number of columns in matrix A.
//     UNeeded     -   0, 1 or 2. See the description of the parameter U.
//     VTNeeded    -   0, 1 or 2. See the description of the parameter VT.
//     AdditionalMemory -
//                     If the parameter:
//                      * equals 0, the algorithm doesn't use additional
//                        memory (lower requirements, lower performance).
//                      * equals 1, the algorithm uses additional
//                        memory of size min(M,N)*min(M,N) of real numbers.
//                        It often speeds up the algorithm.
//                      * equals 2, the algorithm uses additional
//                        memory of size M*min(M,N) of real numbers.
//                        It allows to get a maximum performance.
//                     The recommended value of the parameter is 2.
//
// Outputs:
//     W           -   contains singular values in descending order.
//     U           -   if UNeeded=0, U isn't changed, the left singular vectors
//                     are not calculated.
//                     if Uneeded=1, U contains left singular vectors (first
//                     min(M,N) columns of matrix U). Array whose indexes range
//                     within [0..M-1, 0..Min(M,N)-1].
//                     if UNeeded=2, U contains matrix U wholly. Array whose
//                     indexes range within [0..M-1, 0..M-1].
//     VT          -   if VTNeeded=0, VT isn't changed, the right singular vectors
//                     are not calculated.
//                     if VTNeeded=1, VT contains right singular vectors (first
//                     min(M,N) rows of matrix V^T). Array whose indexes range
//                     within [0..min(M,N)-1, 0..N-1].
//                     if VTNeeded=2, VT contains matrix V^T wholly. Array whose
//                     indexes range within [0..N-1, 0..N-1].
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: bool rmatrixsvd(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const ae_int_t uneeded, const ae_int_t vtneeded, const ae_int_t additionalmemory, real_1d_array &w, real_2d_array &u, real_2d_array &vt);
bool rmatrixsvd(RMatrix *a, ae_int_t m, ae_int_t n, ae_int_t uneeded, ae_int_t vtneeded, ae_int_t additionalmemory, RVector *w, RMatrix *u, RMatrix *vt) {
   ae_frame _frame_block;
   bool isupper;
   ae_int_t minmn;
   ae_int_t ncu;
   ae_int_t nrvt;
   ae_int_t nru;
   ae_int_t ncvt;
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(w);
   SetMatrix(u);
   SetMatrix(vt);
   NewVector(tauq, 0, DT_REAL);
   NewVector(taup, 0, DT_REAL);
   NewVector(tau, 0, DT_REAL);
   NewVector(e, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   NewMatrix(t2, 0, 0, DT_REAL);
   result = true;
   if (m == 0 || n == 0) {
      ae_frame_leave();
      return result;
   }
   ae_assert(uneeded >= 0 && uneeded <= 2, "SVDDecomposition: wrong parameters!");
   ae_assert(vtneeded >= 0 && vtneeded <= 2, "SVDDecomposition: wrong parameters!");
   ae_assert(additionalmemory >= 0 && additionalmemory <= 2, "SVDDecomposition: wrong parameters!");
// initialize
   minmn = imin2(m, n);
   ae_vector_set_length(w, minmn + 1);
   ncu = 0;
   nru = 0;
   if (uneeded == 1) {
      nru = m;
      ncu = minmn;
      ae_matrix_set_length(u, nru, ncu);
   }
   if (uneeded == 2) {
      nru = m;
      ncu = m;
      ae_matrix_set_length(u, nru, ncu);
   }
   nrvt = 0;
   ncvt = 0;
   if (vtneeded == 1) {
      nrvt = minmn;
      ncvt = n;
      ae_matrix_set_length(vt, nrvt, ncvt);
   }
   if (vtneeded == 2) {
      nrvt = n;
      ncvt = n;
      ae_matrix_set_length(vt, nrvt, ncvt);
   }
// M much larger than N
// Use bidiagonal reduction with QR-decomposition
   if ((double)m > 1.6 * n) {
      if (uneeded == 0) {
      // No left singular vectors to be computed
         rmatrixqr(a, m, n, &tau);
         for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
         rmatrixbd(a, n, n, &tauq, &taup);
         rmatrixbdunpackpt(a, n, n, &taup, nrvt, vt);
         rmatrixbdunpackdiagonals(a, n, n, &isupper, w, &e);
         result = rmatrixbdsvd(w, &e, n, isupper, false, u, 0, a, 0, vt, ncvt);
         ae_frame_leave();
         return result;
      } else {
      // Left singular vectors (may be full matrix U) to be computed
         rmatrixqr(a, m, n, &tau);
         rmatrixqrunpackq(a, m, n, &tau, ncu, u);
         for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
         rmatrixbd(a, n, n, &tauq, &taup);
         rmatrixbdunpackpt(a, n, n, &taup, nrvt, vt);
         rmatrixbdunpackdiagonals(a, n, n, &isupper, w, &e);
         if (additionalmemory < 1) {
         // No additional memory can be used
            rmatrixbdmultiplybyq(a, n, n, &tauq, u, m, n, true, false);
            result = rmatrixbdsvd(w, &e, n, isupper, false, u, m, a, 0, vt, ncvt);
         } else {
         // Large U. Transforming intermediate matrix T2
            ae_vector_set_length(&work, imax2(m, n) + 1);
            rmatrixbdunpackq(a, n, n, &tauq, n, &t2);
            copymatrix(u, 0, m - 1, 0, n - 1, a, 0, m - 1, 0, n - 1);
            inplacetranspose(&t2, 0, n - 1, 0, n - 1, &work);
            result = rmatrixbdsvd(w, &e, n, isupper, false, u, 0, &t2, n, vt, ncvt);
            rmatrixgemm(m, n, n, 1.0, a, 0, 0, 0, &t2, 0, 0, 1, 0.0, u, 0, 0);
         }
         ae_frame_leave();
         return result;
      }
   }
// N much larger than M
// Use bidiagonal reduction with LQ-decomposition
   if ((double)n > 1.6 * m) {
      if (vtneeded == 0) {
      // No right singular vectors to be computed
         rmatrixlq(a, m, n, &tau);
         for (i = 0; i < m; i++) {
            for (j = i + 1; j < m; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
         rmatrixbd(a, m, m, &tauq, &taup);
         rmatrixbdunpackq(a, m, m, &tauq, ncu, u);
         rmatrixbdunpackdiagonals(a, m, m, &isupper, w, &e);
         ae_vector_set_length(&work, m + 1);
         inplacetranspose(u, 0, nru - 1, 0, ncu - 1, &work);
         result = rmatrixbdsvd(w, &e, m, isupper, false, a, 0, u, nru, vt, 0);
         inplacetranspose(u, 0, nru - 1, 0, ncu - 1, &work);
         ae_frame_leave();
         return result;
      } else {
      // Right singular vectors (may be full matrix VT) to be computed
         rmatrixlq(a, m, n, &tau);
         rmatrixlqunpackq(a, m, n, &tau, nrvt, vt);
         for (i = 0; i < m; i++) {
            for (j = i + 1; j < m; j++) {
               a->xyR[i][j] = 0.0;
            }
         }
         rmatrixbd(a, m, m, &tauq, &taup);
         rmatrixbdunpackq(a, m, m, &tauq, ncu, u);
         rmatrixbdunpackdiagonals(a, m, m, &isupper, w, &e);
         ae_vector_set_length(&work, imax2(m, n) + 1);
         inplacetranspose(u, 0, nru - 1, 0, ncu - 1, &work);
         if (additionalmemory < 1) {
         // No additional memory available
            rmatrixbdmultiplybyp(a, m, m, &taup, vt, m, n, false, true);
            result = rmatrixbdsvd(w, &e, m, isupper, false, a, 0, u, nru, vt, n);
         } else {
         // Large VT. Transforming intermediate matrix T2
            rmatrixbdunpackpt(a, m, m, &taup, m, &t2);
            result = rmatrixbdsvd(w, &e, m, isupper, false, a, 0, u, nru, &t2, m);
            copymatrix(vt, 0, m - 1, 0, n - 1, a, 0, m - 1, 0, n - 1);
            rmatrixgemm(m, n, m, 1.0, &t2, 0, 0, 0, a, 0, 0, 0, 0.0, vt, 0, 0);
         }
         inplacetranspose(u, 0, nru - 1, 0, ncu - 1, &work);
         ae_frame_leave();
         return result;
      }
   }
// M <= N
// We can use inplace transposition of U to get rid of columnwise operations
   if (m <= n) {
      rmatrixbd(a, m, n, &tauq, &taup);
      rmatrixbdunpackq(a, m, n, &tauq, ncu, u);
      rmatrixbdunpackpt(a, m, n, &taup, nrvt, vt);
      rmatrixbdunpackdiagonals(a, m, n, &isupper, w, &e);
      ae_vector_set_length(&work, m + 1);
      inplacetranspose(u, 0, nru - 1, 0, ncu - 1, &work);
      result = rmatrixbdsvd(w, &e, minmn, isupper, false, a, 0, u, nru, vt, ncvt);
      inplacetranspose(u, 0, nru - 1, 0, ncu - 1, &work);
      ae_frame_leave();
      return result;
   }
// Simple bidiagonal reduction
   rmatrixbd(a, m, n, &tauq, &taup);
   rmatrixbdunpackq(a, m, n, &tauq, ncu, u);
   rmatrixbdunpackpt(a, m, n, &taup, nrvt, vt);
   rmatrixbdunpackdiagonals(a, m, n, &isupper, w, &e);
   if (additionalmemory < 2 || uneeded == 0) {
   // We cant use additional memory or there is no need in such operations
      result = rmatrixbdsvd(w, &e, minmn, isupper, false, u, nru, a, 0, vt, ncvt);
   } else {
   // We can use additional memory
      ae_matrix_set_length(&t2, minmn, m);
      copyandtranspose(u, 0, m - 1, 0, minmn - 1, &t2, 0, minmn - 1, 0, m - 1);
      result = rmatrixbdsvd(w, &e, minmn, isupper, false, u, 0, &t2, m, vt, ncvt);
      copyandtranspose(&t2, 0, minmn - 1, 0, m - 1, u, 0, m - 1, 0, minmn - 1);
   }
   ae_frame_leave();
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
bool rmatrixsvd(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const ae_int_t uneeded, const ae_int_t vtneeded, const ae_int_t additionalmemory, real_1d_array &w, real_2d_array &u, real_2d_array &vt) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::rmatrixsvd(ConstT(ae_matrix, a), m, n, uneeded, vtneeded, additionalmemory, ConstT(ae_vector, w), ConstT(ae_matrix, u), ConstT(ae_matrix, vt));
   alglib_impl::ae_state_clear();
   return Ok;
}
} // end of namespace alglib

// === NORMESTIMATOR Package ===
// Depends on: SPARSE, MATGEN
namespace alglib_impl {
// This procedure initializes matrix norm estimator.
//
// USAGE:
// 1. User initializes algorithm state with NormEstimatorCreate() call
// 2. User calls NormEstimatorEstimateSparse() (or NormEstimatorIteration())
// 3. User calls NormEstimatorResults() to get solution.
//
// Inputs:
//     M       -   number of rows in the matrix being estimated, M > 0
//     N       -   number of columns in the matrix being estimated, N > 0
//     NStart  -   number of random starting vectors
//                 recommended value - at least 5.
//     NIts    -   number of iterations to do with best starting vector
//                 recommended value - at least 5.
//
// Outputs:
//     State   -   structure which stores algorithm state
//
//
// NOTE: this algorithm is effectively deterministic, i.e. it always  returns
// same result when repeatedly called for the same matrix. In fact, algorithm
// uses randomized starting vectors, but internal  random  numbers  generator
// always generates same sequence of the random values (it is a  feature, not
// bug).
//
// Algorithm can be made non-deterministic with NormEstimatorSetSeed(0) call.
// ALGLIB: Copyright 06.12.2011 by Sergey Bochkanov
// API: void normestimatorcreate(const ae_int_t m, const ae_int_t n, const ae_int_t nstart, const ae_int_t nits, normestimatorstate &state);
void normestimatorcreate(ae_int_t m, ae_int_t n, ae_int_t nstart, ae_int_t nits, normestimatorstate *state) {
   SetObj(normestimatorstate, state);
   ae_assert(m > 0, "NormEstimatorCreate: M <= 0");
   ae_assert(n > 0, "NormEstimatorCreate: N <= 0");
   ae_assert(nstart > 0, "NormEstimatorCreate: NStart <= 0");
   ae_assert(nits > 0, "NormEstimatorCreate: NIts <= 0");
   state->m = m;
   state->n = n;
   state->nstart = nstart;
   state->nits = nits;
   state->seedval = 11;
   hqrndrandomize(&state->r);
   ae_vector_set_length(&state->x0, state->n);
   ae_vector_set_length(&state->t, state->m);
   ae_vector_set_length(&state->x1, state->n);
   ae_vector_set_length(&state->xbest, state->n);
   ae_vector_set_length(&state->x, imax2(state->n, state->m));
   ae_vector_set_length(&state->mv, state->m);
   ae_vector_set_length(&state->mtv, state->n);
   state->PQ = -1;
}

// This function changes seed value used by algorithm. In some cases we  need
// deterministic processing, i.e. subsequent calls must return equal results,
// in other cases we need non-deterministic algorithm which returns different
// results for the same matrix on every pass.
//
// Setting zero seed will lead to non-deterministic algorithm, while non-zero
// value will make our algorithm deterministic.
//
// Inputs:
//     State       -   norm estimator state, must be initialized with a  call
//                     to NormEstimatorCreate()
//     SeedVal     -   seed value, >= 0. Zero value = non-deterministic algo.
// ALGLIB: Copyright 06.12.2011 by Sergey Bochkanov
// API: void normestimatorsetseed(const normestimatorstate &state, const ae_int_t seedval);
void normestimatorsetseed(normestimatorstate *state, ae_int_t seedval) {
   ae_assert(seedval >= 0, "NormEstimatorSetSeed: SeedVal<0");
   state->seedval = seedval;
}

// ALGLIB: Copyright 06.12.2011 by Sergey Bochkanov
bool normestimatoriteration(normestimatorstate *state) {
   AutoS ae_int_t n;
   AutoS ae_int_t m;
   AutoS ae_int_t i;
   AutoS ae_int_t itcnt;
   AutoS double v;
   AutoS double growth;
   AutoS double bestgrowth;
// Manually threaded two-way signalling.
// Locals are set arbitrarily the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0; case 1: goto Resume1; case 2: goto Resume2; case 3: goto Resume3;
      default: goto Exit;
   }
Spawn:
   v = 81;
   growth = 255;
   n = state->n;
   m = state->m;
   if (state->seedval > 0) {
      hqrndseed(state->seedval, state->seedval + 2, &state->r);
   }
   bestgrowth = 0.0;
   state->needmtv = state->needmv = false;
   state->xbest.xR[0] = 1.0;
   for (i = 1; i < n; i++) {
      state->xbest.xR[i] = 0.0;
   }
   for (itcnt = 0; itcnt < state->nstart; itcnt++) {
      do {
         v = 0.0;
         for (i = 0; i < n; i++) {
            state->x0.xR[i] = hqrndnormal(&state->r);
            v += ae_sqr(state->x0.xR[i]);
         }
      } while (v == 0.0);
      v = 1 / sqrt(v);
      ae_v_muld(state->x0.xR, 1, n, v);
      ae_v_move(state->x.xR, 1, state->x0.xR, 1, n);
      state->needmv = true, state->PQ = 0; goto Pause; Resume0: state->needmv = false;
      ae_v_move(state->x.xR, 1, state->mv.xR, 1, m);
      state->needmtv = true, state->PQ = 1; goto Pause; Resume1: state->needmtv = false;
      ae_v_move(state->x1.xR, 1, state->mtv.xR, 1, n);
      v = 0.0;
      for (i = 0; i < n; i++) {
         v += ae_sqr(state->x1.xR[i]);
      }
      growth = sqrt(sqrt(v));
      if (growth > bestgrowth) {
         v = 1 / sqrt(v);
         ae_v_moved(state->xbest.xR, 1, state->x1.xR, 1, n, v);
         bestgrowth = growth;
      }
   }
   ae_v_move(state->x0.xR, 1, state->xbest.xR, 1, n);
   for (itcnt = 0; itcnt < state->nits; itcnt++) {
      ae_v_move(state->x.xR, 1, state->x0.xR, 1, n);
      state->needmv = true, state->PQ = 2; goto Pause; Resume2: state->needmv = false;
      ae_v_move(state->x.xR, 1, state->mv.xR, 1, m);
      state->needmtv = true, state->PQ = 3; goto Pause; Resume3: state->needmtv = false;
      ae_v_move(state->x1.xR, 1, state->mtv.xR, 1, n);
      v = 0.0;
      for (i = 0; i < n; i++) {
         v += ae_sqr(state->x1.xR[i]);
      }
      state->repnorm = sqrt(sqrt(v));
      if (v != 0.0) {
         v = 1 / sqrt(v);
         ae_v_moved(state->x0.xR, 1, state->x1.xR, 1, n, v);
      }
   }
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// This function estimates norm of the sparse M*N matrix A.
//
// Inputs:
//     State       -   norm estimator state, must be initialized with a  call
//                     to NormEstimatorCreate()
//     A           -   sparse M*N matrix, must be converted to CRS format
//                     prior to calling this function.
//
// After this function  is  over  you can call NormEstimatorResults() to get
// estimate of the norm(A).
// ALGLIB: Copyright 06.12.2011 by Sergey Bochkanov
// API: void normestimatorestimatesparse(const normestimatorstate &state, const sparsematrix &a);
void normestimatorestimatesparse(normestimatorstate *state, sparsematrix *a) {
   for (normestimatorrestart(state); normestimatoriteration(state); ) {
      if (state->needmv) sparsemv(a, &state->x, &state->mv);
      else if (state->needmtv) sparsemtv(a, &state->x, &state->mtv);
   }
}

// Matrix norm estimation results
//
// Inputs:
//     State   -   algorithm state
//
// Outputs:
//     Nrm     -   estimate of the matrix norm, Nrm >= 0
// ALGLIB: Copyright 06.12.2011 by Sergey Bochkanov
// API: void normestimatorresults(const normestimatorstate &state, double &nrm);
void normestimatorresults(normestimatorstate *state, double *nrm) {
   *nrm = 0;
   *nrm = state->repnorm;
}

// This  function  restarts estimator and prepares it for the next estimation
// round.
//
// Inputs:
//     State   -   algorithm state
// ALGLIB: Copyright 06.12.2011 by Sergey Bochkanov
void normestimatorrestart(normestimatorstate *state) {
   state->PQ = -1;
}

void normestimatorstate_init(void *_p, bool make_automatic) {
   normestimatorstate *p = (normestimatorstate *)_p;
   ae_vector_init(&p->x0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->x1, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->t, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->xbest, 0, DT_REAL, make_automatic);
   hqrndstate_init(&p->r, make_automatic);
   ae_vector_init(&p->x, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->mv, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->mtv, 0, DT_REAL, make_automatic);
}

void normestimatorstate_copy(void *_dst, void *_src, bool make_automatic) {
   normestimatorstate *dst = (normestimatorstate *)_dst;
   normestimatorstate *src = (normestimatorstate *)_src;
   dst->n = src->n;
   dst->m = src->m;
   dst->nstart = src->nstart;
   dst->nits = src->nits;
   dst->seedval = src->seedval;
   ae_vector_copy(&dst->x0, &src->x0, make_automatic);
   ae_vector_copy(&dst->x1, &src->x1, make_automatic);
   ae_vector_copy(&dst->t, &src->t, make_automatic);
   ae_vector_copy(&dst->xbest, &src->xbest, make_automatic);
   hqrndstate_copy(&dst->r, &src->r, make_automatic);
   ae_vector_copy(&dst->x, &src->x, make_automatic);
   ae_vector_copy(&dst->mv, &src->mv, make_automatic);
   ae_vector_copy(&dst->mtv, &src->mtv, make_automatic);
   dst->needmv = src->needmv;
   dst->needmtv = src->needmtv;
   dst->repnorm = src->repnorm;
   dst->PQ = src->PQ;
}

void normestimatorstate_free(void *_p, bool make_automatic) {
   normestimatorstate *p = (normestimatorstate *)_p;
   ae_vector_free(&p->x0, make_automatic);
   ae_vector_free(&p->x1, make_automatic);
   ae_vector_free(&p->t, make_automatic);
   ae_vector_free(&p->xbest, make_automatic);
   hqrndstate_free(&p->r, make_automatic);
   ae_vector_free(&p->x, make_automatic);
   ae_vector_free(&p->mv, make_automatic);
   ae_vector_free(&p->mtv, make_automatic);
}
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the iterative norm estimation algorithm.
// You should use ALGLIB functions to work with this object.
DefClass(normestimatorstate, EndD)

void normestimatorcreate(const ae_int_t m, const ae_int_t n, const ae_int_t nstart, const ae_int_t nits, normestimatorstate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::normestimatorcreate(m, n, nstart, nits, ConstT(normestimatorstate, state));
   alglib_impl::ae_state_clear();
}

void normestimatorsetseed(const normestimatorstate &state, const ae_int_t seedval) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::normestimatorsetseed(ConstT(normestimatorstate, state), seedval);
   alglib_impl::ae_state_clear();
}

void normestimatorestimatesparse(const normestimatorstate &state, const sparsematrix &a) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::normestimatorestimatesparse(ConstT(normestimatorstate, state), ConstT(sparsematrix, a));
   alglib_impl::ae_state_clear();
}

void normestimatorresults(const normestimatorstate &state, double &nrm) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::normestimatorresults(ConstT(normestimatorstate, state), &nrm);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === HSSCHUR Package ===
// Depends on: (AlgLibInternal) ROTATIONS, BLAS
// Depends on: ABLAS
namespace alglib_impl {
void rmatrixinternalschurdecomposition(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   SetVector(wr);
   SetVector(wi);
   *info = 0;
   NewMatrix(h1, 0, 0, DT_REAL);
   NewMatrix(z1, 0, 0, DT_REAL);
   NewVector(wr1, 0, DT_REAL);
   NewVector(wi1, 0, DT_REAL);
// Allocate space
   ae_vector_set_length(wr, n);
   ae_vector_set_length(wi, n);
   if (zneeded == 2) {
      matrixsetlengthatleast(z, n, n);
   }
// MKL version
   if (rmatrixinternalschurdecompositionmkl(h, n, tneeded, zneeded, wr, wi, z, info)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version
   ae_matrix_set_length(&h1, n + 1, n + 1);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         h1.xyR[1 + i][1 + j] = h->xyR[i][j];
      }
   }
   if (zneeded == 1) {
      ae_matrix_set_length(&z1, n + 1, n + 1);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            z1.xyR[1 + i][1 + j] = z->xyR[i][j];
         }
      }
   }
   internalschurdecomposition(&h1, n, tneeded, zneeded, &wr1, &wi1, &z1, info);
   for (i = 0; i < n; i++) {
      wr->xR[i] = wr1.xR[i + 1];
      wi->xR[i] = wi1.xR[i + 1];
   }
   if (tneeded != 0) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            h->xyR[i][j] = h1.xyR[1 + i][1 + j];
         }
      }
   }
   if (zneeded != 0) {
      matrixsetlengthatleast(z, n, n);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            z->xyR[i][j] = z1.xyR[1 + i][1 + j];
         }
      }
   }
   ae_frame_leave();
}

// Subroutine performing  the  Schur  decomposition  of  a  matrix  in  upper
// Hessenberg form using the QR algorithm with multiple shifts.
//
// The  source matrix  H  is  represented as  S'*H*S = T, where H - matrix in
// upper Hessenberg form,  S - orthogonal matrix (Schur vectors),   T - upper
// quasi-triangular matrix (with blocks of sizes  1x1  and  2x2  on  the main
// diagonal).
//
// Inputs:
//     H   -   matrix to be decomposed.
//             Array whose indexes range within [1..N, 1..N].
//     N   -   size of H, N >= 0.
//
//
// Outputs:
//     H   -   contains the matrix T.
//             Array whose indexes range within [1..N, 1..N].
//             All elements below the blocks on the main diagonal are equal
//             to 0.
//     S   -   contains Schur vectors.
//             Array whose indexes range within [1..N, 1..N].
//
// Note 1:
//     The block structure of matrix T could be easily recognized: since  all
//     the elements  below  the blocks are zeros, the elements a[i+1,i] which
//     are equal to 0 show the block border.
//
// Note 2:
//     the algorithm  performance  depends  on  the  value  of  the  internal
//     parameter NS of InternalSchurDecomposition  subroutine  which  defines
//     the number of shifts in the QR algorithm (analog of  the  block  width
//     in block matrix algorithms in linear algebra). If you require  maximum
//     performance  on  your  machine,  it  is  recommended  to  adjust  this
//     parameter manually.
//
// Result:
//     True, if the algorithm has converged and the parameters H and S contain
//         the result.
//     False, if the algorithm has not converged.
//
// Algorithm implemented on the basis of subroutine DHSEQR (LAPACK 3.0 library).
bool upperhessenbergschurdecomposition(RMatrix *h, ae_int_t n, RMatrix *s) {
   ae_frame _frame_block;
   ae_int_t info;
   bool result;
   ae_frame_make(&_frame_block);
   SetMatrix(s);
   NewVector(wi, 0, DT_REAL);
   NewVector(wr, 0, DT_REAL);
   internalschurdecomposition(h, n, 1, 2, &wr, &wi, s, &info);
   result = info == 0;
   ae_frame_leave();
   return result;
}

static double hsschur_extschursign(double a, double b) {
   double result;
   if (b >= 0.0) {
      result = fabs(a);
   } else {
      result = -fabs(a);
   }
   return result;
}

static ae_int_t hsschur_extschursigntoone(double b) {
   ae_int_t result;
   if (b >= 0.0) {
      result = 1;
   } else {
      result = -1;
   }
   return result;
}

static void hsschur_aux2x2schur(double *a, double *b, double *c, double *d, double *rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn) {
   double multpl;
   double aa;
   double bb;
   double bcmax;
   double bcmis;
   double cc;
   double cs1;
   double dd;
   double eps;
   double p;
   double sab;
   double sac;
   double scl;
   double sigma;
   double sn1;
   double tau;
   double temp;
   double z;
   *rt1r = 0;
   *rt1i = 0;
   *rt2r = 0;
   *rt2i = 0;
   *cs = 0;
   *sn = 0;
   multpl = 4.0;
   eps = ae_machineepsilon;
   if (*c == 0.0) {
      *cs = 1.0;
      *sn = 0.0;
   } else {
      if (*b == 0.0) {
      // Swap rows and columns
         *cs = 0.0;
         *sn = 1.0;
         swapr(a, d);
         *b = -*c;
         *c = 0.0;
      } else {
         if (*a - (*d) == 0.0 && hsschur_extschursigntoone(*b) != hsschur_extschursigntoone(*c)) {
            *cs = 1.0;
            *sn = 0.0;
         } else {
            temp = *a - (*d);
            p = 0.5 * temp;
            bcmax = rmax2(fabs(*b), fabs(*c));
            bcmis = rmin2(fabs(*b), fabs(*c)) * hsschur_extschursigntoone(*b) * hsschur_extschursigntoone(*c);
            scl = rmax2(fabs(p), bcmax);
            z = p / scl * p + bcmax / scl * bcmis;
         // If Z is of the order of the machine accuracy, postpone the
         // decision on the nature of eigenvalues
            if (z >= multpl * eps) {
            // Real eigenvalues. Compute A and D.
               z = p + hsschur_extschursign(sqrt(scl) * sqrt(z), p);
               *a = *d + z;
               *d -= bcmax / z * bcmis;
            // Compute B and the rotation matrix
               tau = safepythag2(*c, z);
               *cs = z / tau;
               *sn = *c / tau;
               *b -= *c;
               *c = 0.0;
            } else {
            // Complex eigenvalues, or real (almost) equal eigenvalues.
            // Make diagonal elements equal.
               sigma = *b + (*c);
               tau = safepythag2(sigma, temp);
               *cs = sqrt(0.5 * (1 + fabs(sigma) / tau));
               *sn = -p / (tau * (*cs)) * hsschur_extschursign(1.0, sigma);
            // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
            //         [ CC  DD ]   [ C  D ] [ SN  CS ]
               aa = *a * (*cs) + *b * (*sn);
               bb = -*a * (*sn) + *b * (*cs);
               cc = *c * (*cs) + *d * (*sn);
               dd = -*c * (*sn) + *d * (*cs);
            // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
            //         [ C  D ]   [-SN  CS ] [ CC  DD ]
               *a = aa * (*cs) + cc * (*sn);
               *b = bb * (*cs) + dd * (*sn);
               *c = -aa * (*sn) + cc * (*cs);
               *d = -bb * (*sn) + dd * (*cs);
               temp = 0.5 * (*a + (*d));
               *a = temp;
               *d = temp;
               if (*c != 0.0) {
                  if (*b != 0.0) {
                     if (hsschur_extschursigntoone(*b) == hsschur_extschursigntoone(*c)) {
                     // Real eigenvalues: reduce to upper triangular form
                        sab = sqrt(fabs(*b));
                        sac = sqrt(fabs(*c));
                        p = hsschur_extschursign(sab * sac, *c);
                        tau = 1 / sqrt(fabs(*b + (*c)));
                        *a = temp + p;
                        *d = temp - p;
                        *b -= *c;
                        *c = 0.0;
                        cs1 = sab * tau;
                        sn1 = sac * tau;
                        temp = *cs * cs1 - *sn * sn1;
                        *sn = *cs * sn1 + *sn * cs1;
                        *cs = temp;
                     }
                  } else {
                     *b = -*c;
                     *c = 0.0;
                     temp = *cs;
                     *cs = -*sn;
                     *sn = temp;
                  }
               }
            }
         }
      }
   }
// Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
   *rt1r = *a;
   *rt2r = *d;
   if (*c == 0.0) {
      *rt1i = 0.0;
      *rt2i = 0.0;
   } else {
      *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
      *rt2i = -*rt1i;
   }
}

// Translation of DLAHQR from LAPACK.
static void hsschur_internalauxschur(bool wantt, bool wantz, ae_int_t n, ae_int_t ilo, ae_int_t ihi, RMatrix *h, RVector *wr, RVector *wi, ae_int_t iloz, ae_int_t ihiz, RMatrix *z, RVector *work, RVector *workv3, RVector *workc1, RVector *works1, ae_int_t *info) {
   double safmin;
   double tst;
   double ab;
   double ba;
   double aa;
   double bb;
   double rt1r;
   double rt1i;
   double rt2r;
   double rt2i;
   double tr;
   double det;
   double rtdisc;
   double h21s;
   ae_int_t i;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t itmax;
   ae_int_t its;
   ae_int_t j;
   ae_int_t k;
   ae_int_t l;
   ae_int_t m;
   ae_int_t nh;
   ae_int_t nr;
   ae_int_t nz;
   double cs;
   double h11;
   double h12;
   double h21;
   double h22;
   double s;
   double smlnum;
   double sn;
   double sum;
   double t1;
   double t2;
   double t3;
   double v2;
   double v3;
   bool failflag;
   double dat1;
   double dat2;
   ae_int_t p1;
   double him1im1;
   double him1i;
   double hiim1;
   double hii;
   double wrim1;
   double wri;
   double wiim1;
   double wii;
   double ulp;
   *info = 0;
   *info = 0;
   dat1 = 0.75;
   dat2 = -0.4375;
// Quick return if possible
   if (n == 0) {
      return;
   }
   if (ilo == ihi) {
      wr->xR[ilo] = h->xyR[ilo][ilo];
      wi->xR[ilo] = 0.0;
      return;
   }
// ==== clear out the trash ====
   for (j = ilo; j < ihi - 2; j++) {
      h->xyR[j + 2][j] = 0.0;
      h->xyR[j + 3][j] = 0.0;
   }
   if (ilo < ihi - 1) {
      h->xyR[ihi][ihi - 2] = 0.0;
   }
   nh = ihi - ilo + 1;
   nz = ihiz - iloz + 1;
// Set machine-dependent constants for the stopping criterion.
   safmin = ae_minrealnumber;
   ulp = ae_machineepsilon;
   smlnum = safmin * (nh / ulp);
// I1 and I2 are the indices of the first row and last column of H
// to which transformations must be applied. If eigenvalues only are
// being computed, I1 and I2 are set inside the main loop.
//
// Setting them to large negative value helps to debug possible errors
// due to uninitialized variables; also it helps to avoid compiler
// warnings.
   i1 = -99999;
   i2 = -99999;
   if (wantt) {
      i1 = 1;
      i2 = n;
   }
// ITMAX is the total number of QR iterations allowed.
   itmax = 30 * imax2(10, nh);
// The main loop begins here. I is the loop index and decreases from
// IHI to ILO in steps of 1 or 2. Each iteration of the loop works
// with the active submatrix in rows and columns L to I.
// Eigenvalues I+1 to IHI have already converged. Either L = ILO or
// H(L,L-1) is negligible so that the matrix splits.
   i = ihi;
   while (true) {
      l = ilo;
      if (i < ilo) {
         return;
      }
   // Perform QR iterations on rows and columns ILO to I until a
   // submatrix of order 1 or 2 splits off at the bottom because a
   // subdiagonal element has become negligible.
      failflag = true;
      for (its = 0; its <= itmax; its++) {
      // Look for a single small subdiagonal element.
         for (k = i; k >= l + 1; k--) {
            if (SmallAtR(h->xyR[k][k - 1], smlnum)) {
               break;
            }
            tst = fabs(h->xyR[k - 1][k - 1]) + fabs(h->xyR[k][k]);
            if (tst == 0.0) {
               if (k - 2 >= ilo) {
                  tst += fabs(h->xyR[k - 1][k - 2]);
               }
               if (k + 1 <= ihi) {
                  tst += fabs(h->xyR[k + 1][k]);
               }
            }
         // ==== The following is a conservative small subdiagonal
         // .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
         // .    1997). It has better mathematical foundation and
         // .    improves accuracy in some cases.  ====
            if (SmallAtR(h->xyR[k][k - 1], ulp * tst)) {
               ab = rmax2(fabs(h->xyR[k][k - 1]), fabs(h->xyR[k - 1][k]));
               ba = rmin2(fabs(h->xyR[k][k - 1]), fabs(h->xyR[k - 1][k]));
               aa = rmax2(fabs(h->xyR[k][k]), fabs(h->xyR[k - 1][k - 1] - h->xyR[k][k]));
               bb = rmin2(fabs(h->xyR[k][k]), fabs(h->xyR[k - 1][k - 1] - h->xyR[k][k]));
               s = aa + ab;
               if (ba * (ab / s) <= rmax2(smlnum, ulp * (bb * (aa / s)))) {
                  break;
               }
            }
         }
         l = k;
         if (l > ilo) {
         // H(L,L-1) is negligible
            h->xyR[l][l - 1] = 0.0;
         }
      // Exit from loop if a submatrix of order 1 or 2 has split off.
         if (l >= i - 1) {
            failflag = false;
            break;
         }
      // Now the active submatrix is in rows and columns L to I. If
      // eigenvalues only are being computed, only the active submatrix
      // need be transformed.
         if (!wantt) {
            i1 = l;
            i2 = i;
         }
      // Shifts
         if (its == 10) {
         // Exceptional shift.
            s = fabs(h->xyR[l + 1][l]) + fabs(h->xyR[l + 2][l + 1]);
            h11 = dat1 * s + h->xyR[l][l];
            h12 = dat2 * s;
            h21 = s;
            h22 = h11;
         } else {
            if (its == 20) {
            // Exceptional shift.
               s = fabs(h->xyR[i][i - 1]) + fabs(h->xyR[i - 1][i - 2]);
               h11 = dat1 * s + h->xyR[i][i];
               h12 = dat2 * s;
               h21 = s;
               h22 = h11;
            } else {
            // Prepare to use Francis' double shift
            // (i.e. 2nd degree generalized Rayleigh quotient)
               h11 = h->xyR[i - 1][i - 1];
               h21 = h->xyR[i][i - 1];
               h12 = h->xyR[i - 1][i];
               h22 = h->xyR[i][i];
            }
         }
         s = fabs(h11) + fabs(h12) + fabs(h21) + fabs(h22);
         if (s == 0.0) {
            rt1r = 0.0;
            rt1i = 0.0;
            rt2r = 0.0;
            rt2i = 0.0;
         } else {
            h11 /= s;
            h21 /= s;
            h12 /= s;
            h22 /= s;
            tr = (h11 + h22) / 2;
            det = (h11 - tr) * (h22 - tr) - h12 * h21;
            rtdisc = sqrt(fabs(det));
            if (det >= 0.0) {
            // ==== complex conjugate shifts ====
               rt1r = tr * s;
               rt2r = rt1r;
               rt1i = rtdisc * s;
               rt2i = -rt1i;
            } else {
            // ==== real shifts (use only one of them)  ====
               rt1r = tr + rtdisc;
               rt2r = tr - rtdisc;
               if (fabs(rt1r - h22) <= fabs(rt2r - h22)) {
                  rt1r *= s;
                  rt2r = rt1r;
               } else {
                  rt2r *= s;
                  rt1r = rt2r;
               }
               rt1i = 0.0;
               rt2i = 0.0;
            }
         }
      // Look for two consecutive small subdiagonal elements.
         for (m = i - 2; m >= l; m--) {
         // Determine the effect of starting the double-shift QR
         // iteration at row M, and see if this would make H(M,M-1)
         // negligible.  (The following uses scaling to avoid
         // overflows and most underflows.)
            h21s = h->xyR[m + 1][m];
            s = fabs(h->xyR[m][m] - rt2r) + fabs(rt2i) + fabs(h21s);
            h21s = h->xyR[m + 1][m] / s;
            workv3->xR[1] = h21s * h->xyR[m][m + 1] + (h->xyR[m][m] - rt1r) * ((h->xyR[m][m] - rt2r) / s) - rt1i * (rt2i / s);
            workv3->xR[2] = h21s * (h->xyR[m][m] + h->xyR[m + 1][m + 1] - rt1r - rt2r);
            workv3->xR[3] = h21s * h->xyR[m + 2][m + 1];
            s = fabs(workv3->xR[1]) + fabs(workv3->xR[2]) + fabs(workv3->xR[3]);
            workv3->xR[1] /= s;
            workv3->xR[2] /= s;
            workv3->xR[3] /= s;
            if (m == l) {
               break;
            }
            if (fabs(h->xyR[m][m - 1]) * (fabs(workv3->xR[2]) + fabs(workv3->xR[3])) <= ulp * fabs(workv3->xR[1]) * (fabs(h->xyR[m - 1][m - 1]) + fabs(h->xyR[m][m]) + fabs(h->xyR[m + 1][m + 1]))) {
               break;
            }
         }
      // Double-shift QR step
         for (k = m; k < i; k++) {
         // The first iteration of this loop determines a reflection G
         // from the vector V and applies it from left and right to H,
         // thus creating a nonzero bulge below the subdiagonal.
         //
         // Each subsequent iteration determines a reflection G to
         // restore the Hessenberg form in the (K-1)th column, and thus
         // chases the bulge one step toward the bottom of the active
         // submatrix. NR is the order of G.
            nr = imin2(3, i - k + 1);
            if (k > m) {
               for (p1 = 1; p1 <= nr; p1++) {
                  workv3->xR[p1] = h->xyR[k + p1 - 1][k - 1];
               }
            }
            generatereflection(workv3, nr, &t1);
            if (k > m) {
               h->xyR[k][k - 1] = workv3->xR[1];
               h->xyR[k + 1][k - 1] = 0.0;
               if (k < i - 1) {
                  h->xyR[k + 2][k - 1] = 0.0;
               }
            } else {
               if (m > l) {
               // ==== Use the following instead of
               // H( K, K-1 ) = -H( K, K-1 ) to
               // avoid a bug when v(2) and v(3)
               // underflow. ====
                  h->xyR[k][k - 1] *= 1 - t1;
               }
            }
            v2 = workv3->xR[2];
            t2 = t1 * v2;
            if (nr == 3) {
               v3 = workv3->xR[3];
               t3 = t1 * v3;
            // Apply G from the left to transform the rows of the matrix
            // in columns K to I2.
               for (j = k; j <= i2; j++) {
                  sum = h->xyR[k][j] + v2 * h->xyR[k + 1][j] + v3 * h->xyR[k + 2][j];
                  h->xyR[k][j] -= sum * t1;
                  h->xyR[k + 1][j] -= sum * t2;
                  h->xyR[k + 2][j] -= sum * t3;
               }
            // Apply G from the right to transform the columns of the
            // matrix in rows I1 to min(K+3,I).
               for (j = i1; j <= imin2(k + 3, i); j++) {
                  sum = h->xyR[j][k] + v2 * h->xyR[j][k + 1] + v3 * h->xyR[j][k + 2];
                  h->xyR[j][k] -= sum * t1;
                  h->xyR[j][k + 1] -= sum * t2;
                  h->xyR[j][k + 2] -= sum * t3;
               }
               if (wantz) {
               // Accumulate transformations in the matrix Z
                  for (j = iloz; j <= ihiz; j++) {
                     sum = z->xyR[j][k] + v2 * z->xyR[j][k + 1] + v3 * z->xyR[j][k + 2];
                     z->xyR[j][k] -= sum * t1;
                     z->xyR[j][k + 1] -= sum * t2;
                     z->xyR[j][k + 2] -= sum * t3;
                  }
               }
            } else {
               if (nr == 2) {
               // Apply G from the left to transform the rows of the matrix
               // in columns K to I2.
                  for (j = k; j <= i2; j++) {
                     sum = h->xyR[k][j] + v2 * h->xyR[k + 1][j];
                     h->xyR[k][j] -= sum * t1;
                     h->xyR[k + 1][j] -= sum * t2;
                  }
               // Apply G from the right to transform the columns of the
               // matrix in rows I1 to min(K+3,I).
                  for (j = i1; j <= i; j++) {
                     sum = h->xyR[j][k] + v2 * h->xyR[j][k + 1];
                     h->xyR[j][k] -= sum * t1;
                     h->xyR[j][k + 1] -= sum * t2;
                  }
                  if (wantz) {
                  // Accumulate transformations in the matrix Z
                     for (j = iloz; j <= ihiz; j++) {
                        sum = z->xyR[j][k] + v2 * z->xyR[j][k + 1];
                        z->xyR[j][k] -= sum * t1;
                        z->xyR[j][k + 1] -= sum * t2;
                     }
                  }
               }
            }
         }
      }
   // Failure to converge in remaining number of iterations
      if (failflag) {
         *info = i;
         return;
      }
   // Convergence
      if (l == i) {
      // H(I,I-1) is negligible: one eigenvalue has converged.
         wr->xR[i] = h->xyR[i][i];
         wi->xR[i] = 0.0;
      } else {
         if (l == i - 1) {
         // H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
         //
         // Transform the 2-by-2 submatrix to standard Schur form,
         // and compute and store the eigenvalues.
            him1im1 = h->xyR[i - 1][i - 1];
            him1i = h->xyR[i - 1][i];
            hiim1 = h->xyR[i][i - 1];
            hii = h->xyR[i][i];
            hsschur_aux2x2schur(&him1im1, &him1i, &hiim1, &hii, &wrim1, &wiim1, &wri, &wii, &cs, &sn);
            wr->xR[i - 1] = wrim1;
            wi->xR[i - 1] = wiim1;
            wr->xR[i] = wri;
            wi->xR[i] = wii;
            h->xyR[i - 1][i - 1] = him1im1;
            h->xyR[i - 1][i] = him1i;
            h->xyR[i][i - 1] = hiim1;
            h->xyR[i][i] = hii;
            if (wantt) {
            // Apply the transformation to the rest of H.
               if (i2 > i) {
                  workc1->xR[1] = cs;
                  works1->xR[1] = sn;
                  applyrotationsfromtheleft(true, i - 1, i, i + 1, i2, workc1, works1, h, work);
               }
               workc1->xR[1] = cs;
               works1->xR[1] = sn;
               applyrotationsfromtheright(true, i1, i - 2, i - 1, i, workc1, works1, h, work);
            }
            if (wantz) {
            // Apply the transformation to Z.
               workc1->xR[1] = cs;
               works1->xR[1] = sn;
               applyrotationsfromtheright(true, iloz, iloz + nz - 1, i - 1, i, workc1, works1, z, work);
            }
         }
      }
   // return to start of the main loop with new value of I.
      i = l - 1;
   }
}

void internalschurdecomposition(RMatrix *h, ae_int_t n, ae_int_t tneeded, ae_int_t zneeded, RVector *wr, RVector *wi, RMatrix *z, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t i1;
   ae_int_t i2;
   ae_int_t ierr;
   ae_int_t ii;
   ae_int_t itemp;
   ae_int_t itn;
   ae_int_t its;
   ae_int_t j;
   ae_int_t k;
   ae_int_t l;
   ae_int_t maxb;
   ae_int_t nr;
   ae_int_t ns;
   ae_int_t nv;
   double absw;
   double smlnum;
   double tau;
   double temp;
   double tst1;
   double ulp;
   double unfl;
   bool initz;
   bool wantt;
   bool wantz;
   double cnst;
   bool failflag;
   ae_int_t p1;
   double vt;
   ae_frame_make(&_frame_block);
   SetVector(wr);
   SetVector(wi);
   *info = 0;
   NewVector(work, 0, DT_REAL);
   NewMatrix(s, 0, 0, DT_REAL);
   NewVector(v, 0, DT_REAL);
   NewVector(vv, 0, DT_REAL);
   NewVector(workc1, 0, DT_REAL);
   NewVector(works1, 0, DT_REAL);
   NewVector(workv3, 0, DT_REAL);
   NewVector(tmpwr, 0, DT_REAL);
   NewVector(tmpwi, 0, DT_REAL);
// Set the order of the multi-shift QR algorithm to be used.
// If you want to tune algorithm, change this values
   ns = 12;
   maxb = 50;
// Now 2 < NS <= MAXB < NH.
   maxb = imax2(3, maxb);
   ns = imin2(maxb, ns);
// Initialize
   cnst = 1.5;
   ae_vector_set_length(&work, imax2(n, 1) + 1);
   ae_matrix_set_length(&s, ns + 1, ns + 1);
   ae_vector_set_length(&v, ns + 1 + 1);
   ae_vector_set_length(&vv, ns + 1 + 1);
   ae_vector_set_length(wr, imax2(n, 1) + 1);
   ae_vector_set_length(wi, imax2(n, 1) + 1);
   ae_vector_set_length(&workc1, 1 + 1);
   ae_vector_set_length(&works1, 1 + 1);
   ae_vector_set_length(&workv3, 3 + 1);
   ae_vector_set_length(&tmpwr, imax2(n, 1) + 1);
   ae_vector_set_length(&tmpwi, imax2(n, 1) + 1);
   ae_assert(n >= 0, "InternalSchurDecomposition: incorrect N!");
   ae_assert(tneeded == 0 || tneeded == 1, "InternalSchurDecomposition: incorrect TNeeded!");
   ae_assert(zneeded == 0 || zneeded == 1 || zneeded == 2, "InternalSchurDecomposition: incorrect ZNeeded!");
   wantt = tneeded == 1;
   initz = zneeded == 2;
   wantz = zneeded != 0;
   *info = 0;
// Initialize Z, if necessary
   if (initz) {
      matrixsetlengthatleast(z, n + 1, n + 1);
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= n; j++) {
            if (i == j) {
               z->xyR[i][j] = 1.0;
            } else {
               z->xyR[i][j] = 0.0;
            }
         }
      }
   }
// Quick return if possible
   if (n == 0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      wr->xR[1] = h->xyR[1][1];
      wi->xR[1] = 0.0;
      ae_frame_leave();
      return;
   }
// Set rows and columns 1 to N to zero below the first
// subdiagonal.
   for (j = 1; j < n - 1; j++) {
      for (i = j + 2; i <= n; i++) {
         h->xyR[i][j] = 0.0;
      }
   }
// Test if N is sufficiently small
   if (ns <= 2 || ns > n || maxb >= n) {
   // Use the standard double-shift algorithm
      hsschur_internalauxschur(wantt, wantz, n, 1, n, h, wr, wi, 1, n, z, &work, &workv3, &workc1, &works1, info);
   // fill entries under diagonal blocks of T with zeros
      if (wantt) {
         j = 1;
         while (j <= n) {
            if (wi->xR[j] == 0.0) {
               for (i = j + 1; i <= n; i++) {
                  h->xyR[i][j] = 0.0;
               }
               j++;
            } else {
               for (i = j + 2; i <= n; i++) {
                  h->xyR[i][j] = 0.0;
                  h->xyR[i][j + 1] = 0.0;
               }
               j += 2;
            }
         }
      }
      ae_frame_leave();
      return;
   }
   unfl = ae_minrealnumber;
   ulp = 2 * ae_machineepsilon;
   smlnum = unfl * (n / ulp);
// I1 and I2 are the indices of the first row and last column of H
// to which transformations must be applied. If eigenvalues only are
// being computed, I1 and I2 are set inside the main loop.
   i1 = 1;
   i2 = n;
// ITN is the total number of multiple-shift QR iterations allowed.
   itn = 30 * n;
// The main loop begins here. I is the loop index and decreases from
// IHI to ILO in steps of at most MAXB. Each iteration of the loop
// works with the active submatrix in rows and columns L to I.
// Eigenvalues I+1 to IHI have already converged. Either L = ILO or
// H(L,L-1) is negligible so that the matrix splits.
   i = n;
   while (true) {
      l = 1;
      if (i < 1) {
      // fill entries under diagonal blocks of T with zeros
         if (wantt) {
            j = 1;
            while (j <= n) {
               if (wi->xR[j] == 0.0) {
                  for (i = j + 1; i <= n; i++) {
                     h->xyR[i][j] = 0.0;
                  }
                  j++;
               } else {
                  for (i = j + 2; i <= n; i++) {
                     h->xyR[i][j] = 0.0;
                     h->xyR[i][j + 1] = 0.0;
                  }
                  j += 2;
               }
            }
         }
         break;
      }
   // Perform multiple-shift QR iterations on rows and columns ILO to I
   // until a submatrix of order at most MAXB splits off at the bottom
   // because a subdiagonal element has become negligible.
      failflag = true;
      for (its = 0; its <= itn; its++) {
      // Look for a single small subdiagonal element.
         for (k = i; k >= l + 1; k--) {
            tst1 = fabs(h->xyR[k - 1][k - 1]) + fabs(h->xyR[k][k]);
            if (tst1 == 0.0) {
               tst1 = upperhessenberg1norm(h, l, i, l, i, &work);
            }
            if (SmallAtR(h->xyR[k][k - 1], rmax2(ulp * tst1, smlnum))) {
               break;
            }
         }
         l = k;
         if (l > 1) {
         // H(L,L-1) is negligible.
            h->xyR[l][l - 1] = 0.0;
         }
      // Exit from loop if a submatrix of order <= MAXB has split off.
         if (l >= i - maxb + 1) {
            failflag = false;
            break;
         }
      // Now the active submatrix is in rows and columns L to I. If
      // eigenvalues only are being computed, only the active submatrix
      // need be transformed.
         if (its == 20 || its == 30) {
         // Exceptional shifts.
            for (ii = i - ns + 1; ii <= i; ii++) {
               wr->xR[ii] = cnst * (fabs(h->xyR[ii][ii - 1]) + fabs(h->xyR[ii][ii]));
               wi->xR[ii] = 0.0;
            }
         } else {
         // Use eigenvalues of trailing submatrix of order NS as shifts.
            copymatrix(h, i - ns + 1, i, i - ns + 1, i, &s, 1, ns, 1, ns);
            hsschur_internalauxschur(false, false, ns, 1, ns, &s, &tmpwr, &tmpwi, 1, ns, z, &work, &workv3, &workc1, &works1, &ierr);
            for (p1 = 1; p1 <= ns; p1++) {
               wr->xR[i - ns + p1] = tmpwr.xR[p1];
               wi->xR[i - ns + p1] = tmpwi.xR[p1];
            }
            if (ierr > 0) {
            // If DLAHQR failed to compute all NS eigenvalues, use the
            // unconverged diagonal elements as the remaining shifts.
               for (ii = 1; ii <= ierr; ii++) {
                  wr->xR[i - ns + ii] = s.xyR[ii][ii];
                  wi->xR[i - ns + ii] = 0.0;
               }
            }
         }
      // Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
      // where G is the Hessenberg submatrix H(L:I,L:I) and w is
      // the vector of shifts (stored in WR and WI). The result is
      // stored in the local array V.
         v.xR[1] = 1.0;
         for (ii = 2; ii <= ns + 1; ii++) {
            v.xR[ii] = 0.0;
         }
         nv = 1;
         for (j = i - ns + 1; j <= i; j++) {
            if (wi->xR[j] >= 0.0) {
               if (wi->xR[j] == 0.0) {
               // real shift
                  p1 = nv + 1;
                  ae_v_move(&vv.xR[1], 1, &v.xR[1], 1, p1);
                  matrixvectormultiply(h, l, l + nv, l, l + nv - 1, false, &vv, 1, nv, 1.0, &v, 1, nv + 1, -wr->xR[j]);
                  nv++;
               } else {
                  if (wi->xR[j] > 0.0) {
                  // complex conjugate pair of shifts
                     p1 = nv + 1;
                     ae_v_move(&vv.xR[1], 1, &v.xR[1], 1, p1);
                     matrixvectormultiply(h, l, l + nv, l, l + nv - 1, false, &v, 1, nv, 1.0, &vv, 1, nv + 1, -2 * wr->xR[j]);
                     itemp = vectoridxabsmax(&vv, 1, nv + 1);
                     temp = 1 / rmax2(fabs(vv.xR[itemp]), smlnum);
                     p1 = nv + 1;
                     ae_v_muld(&vv.xR[1], 1, p1, temp);
                     absw = safepythag2(wr->xR[j], wi->xR[j]);
                     temp *= absw * absw;
                     matrixvectormultiply(h, l, l + nv + 1, l, l + nv, false, &vv, 1, nv + 1, 1.0, &v, 1, nv + 2, temp);
                     nv += 2;
                  }
               }
            // Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
            // reset it to the unit vector.
               itemp = vectoridxabsmax(&v, 1, nv);
               temp = fabs(v.xR[itemp]);
               if (temp == 0.0) {
                  v.xR[1] = 1.0;
                  for (ii = 2; ii <= nv; ii++) {
                     v.xR[ii] = 0.0;
                  }
               } else {
                  temp = rmax2(temp, smlnum);
                  vt = 1 / temp;
                  ae_v_muld(&v.xR[1], 1, nv, vt);
               }
            }
         }
      // Multiple-shift QR step
         for (k = l; k < i; k++) {
         // The first iteration of this loop determines a reflection G
         // from the vector V and applies it from left and right to H,
         // thus creating a nonzero bulge below the subdiagonal.
         //
         // Each subsequent iteration determines a reflection G to
         // restore the Hessenberg form in the (K-1)th column, and thus
         // chases the bulge one step toward the bottom of the active
         // submatrix. NR is the order of G.
            nr = imin2(ns + 1, i - k + 1);
            if (k > l) {
               p1 = k - 1;
               ae_v_move(&v.xR[1], 1, &h->xyR[k][p1], h->stride, nr);
            }
            generatereflection(&v, nr, &tau);
            if (k > l) {
               h->xyR[k][k - 1] = v.xR[1];
               for (ii = k + 1; ii <= i; ii++) {
                  h->xyR[ii][k - 1] = 0.0;
               }
            }
            v.xR[1] = 1.0;
         // Apply G from the left to transform the rows of the matrix in
         // columns K to I2.
            applyreflectionfromtheleft(h, tau, &v, k, k + nr - 1, k, i2, &work);
         // Apply G from the right to transform the columns of the
         // matrix in rows I1 to min(K+NR,I).
            applyreflectionfromtheright(h, tau, &v, i1, imin2(k + nr, i), k, k + nr - 1, &work);
            if (wantz) {
            // Accumulate transformations in the matrix Z
               applyreflectionfromtheright(z, tau, &v, 1, n, k, k + nr - 1, &work);
            }
         }
      }
   // Failure to converge in remaining number of iterations
      if (failflag) {
         *info = i;
         break;
      }
   // A submatrix of order <= MAXB in rows and columns L to I has split
   // off. Use the double-shift QR algorithm to handle it.
      hsschur_internalauxschur(wantt, wantz, n, l, i, h, wr, wi, 1, n, z, &work, &workv3, &workc1, &works1, info);
      if (*info > 0) {
         break;
      }
   // Decrement number of remaining iterations, and return to start of
   // the main loop with a new value of I.
      itn -= its;
      i = l - 1;
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

// === EVD Package ===
// Depends on: (AlgLibInternal) BASICSTATOPS
// Depends on: SPARSE, MATGEN, ORTFAC, HSSCHUR
namespace alglib_impl {
static const ae_int_t evd_stepswithintol = 2;

// This function initializes subspace iteration solver. This solver  is  used
// to solve symmetric real eigenproblems where just a few (top K) eigenvalues
// and corresponding eigenvectors is required.
//
// This solver can be significantly faster than  complete  EVD  decomposition
// in the following case:
// * when only just a small fraction  of  top  eigenpairs  of dense matrix is
//   required. When K approaches N, this solver is slower than complete dense
//   EVD
// * when problem matrix is sparse (and/or is not known explicitly, i.e. only
//   matrix-matrix product can be performed)
//
// USAGE (explicit dense/sparse matrix):
// 1. User initializes algorithm state with eigsubspacecreate() call
// 2. [optional] User tunes solver parameters by calling eigsubspacesetcond()
//    or other functions
// 3. User  calls  eigsubspacesolvedense() or eigsubspacesolvesparse() methods,
//    which take algorithm state and 2D array or alglib.sparsematrix object.
//
// USAGE (out-of-core mode):
// 1. User initializes algorithm state with eigsubspacecreate() call
// 2. [optional] User tunes solver parameters by calling eigsubspacesetcond()
//    or other functions
// 3. User activates out-of-core mode of  the  solver  and  repeatedly  calls
//    communication functions in a loop like below:
//    > alglib.eigsubspaceoocstart(state)
//    > while alglib.eigsubspaceooccontinue(state) do
//    >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
//    >     alglib.eigsubspaceoocgetrequestdata(state, out X)
//    >     [calculate  Y=A*X, with X=R^NxM]
//    >     alglib.eigsubspaceoocsendresult(state, in Y)
//    > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
//
// Inputs:
//     N       -   problem dimensionality, N > 0
//     K       -   number of top eigenvector to calculate, 0 < K <= N.
//
// Outputs:
//     State   -   structure which stores algorithm state
//
// NOTE: if you solve many similar EVD problems you may  find  it  useful  to
//       reuse previous subspace as warm-start point for new EVD problem.  It
//       can be done with eigsubspacesetwarmstart() function.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspacecreate(const ae_int_t n, const ae_int_t k, eigsubspacestate &state);
void eigsubspacecreate(ae_int_t n, ae_int_t k, eigsubspacestate *state) {
   SetObj(eigsubspacestate, state);
   ae_assert(n > 0, "EigSubspaceCreate: N <= 0");
   ae_assert(k > 0, "EigSubspaceCreate: K <= 0");
   ae_assert(k <= n, "EigSubspaceCreate: K > N");
   eigsubspacecreatebuf(n, k, state);
}

// Buffered version of constructor which aims to reuse  previously  allocated
// memory as much as possible.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspacecreatebuf(const ae_int_t n, const ae_int_t k, const eigsubspacestate &state);
void eigsubspacecreatebuf(ae_int_t n, ae_int_t k, eigsubspacestate *state) {
   ae_assert(n > 0, "EigSubspaceCreate: N <= 0");
   ae_assert(k > 0, "EigSubspaceCreate: K <= 0");
   ae_assert(k <= n, "EigSubspaceCreate: K > N");
// Initialize algorithm parameters
   state->running = false;
   state->n = n;
   state->k = k;
   state->nwork = imin2(imax2(2 * k, 8), n);
   state->eigenvectorsneeded = 1;
   state->usewarmstart = false;
   state->firstcall = true;
   eigsubspacesetcond(state, 0.0, 0);
// Allocate temporaries
   matrixsetlengthatleast(&state->x, state->n, state->nwork);
   matrixsetlengthatleast(&state->ax, state->n, state->nwork);
}

// This function sets stopping critera for the solver:
// * error in eigenvector/value allowed by solver
// * maximum number of iterations to perform
//
// Inputs:
//     State       -   solver structure
//     Eps         -   eps >= 0,  with non-zero value used to tell solver  that
//                     it can  stop  after  all  eigenvalues  converged  with
//                     error  roughly  proportional  to  eps*MAX(LAMBDA_MAX),
//                     where LAMBDA_MAX is a maximum eigenvalue.
//                     Zero  value  means  that  no  check  for  precision is
//                     performed.
//     MaxIts      -   maxits >= 0,  with non-zero value used  to  tell  solver
//                     that it can stop after maxits  steps  (no  matter  how
//                     precise current estimate is)
//
// NOTE: passing  eps=0  and  maxits=0  results  in  automatic  selection  of
//       moderate eps as stopping criteria (1.0E-6 in current implementation,
//       but it may change without notice).
//
// NOTE: very small values of eps are possible (say, 1.0E-12),  although  the
//       larger problem you solve (N and/or K), the  harder  it  is  to  find
//       precise eigenvectors because rounding errors tend to accumulate.
//
// NOTE: passing non-zero eps results in  some performance  penalty,  roughly
//       equal to 2N*(2K)^2 FLOPs per iteration. These additional computations
//       are required in order to estimate current error in  eigenvalues  via
//       Rayleigh-Ritz process.
//       Most of this additional time is  spent  in  construction  of  ~2Kx2K
//       symmetric  subproblem  whose  eigenvalues  are  checked  with  exact
//       eigensolver.
//       This additional time is negligible if you search for eigenvalues  of
//       the large dense matrix, but may become noticeable on  highly  sparse
//       EVD problems, where cost of matrix-matrix product is low.
//       If you set eps to exactly zero,  Rayleigh-Ritz  phase  is completely
//       turned off.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspacesetcond(const eigsubspacestate &state, const double eps, const ae_int_t maxits);
void eigsubspacesetcond(eigsubspacestate *state, double eps, ae_int_t maxits) {
   ae_assert(!state->running, "EigSubspaceSetCond: solver is already running");
   ae_assert(isfinite(eps) && eps >= 0.0, "EigSubspaceSetCond: Eps < 0 or NAN/INF");
   ae_assert(maxits >= 0, "EigSubspaceSetCond: MaxIts<0");
   if (eps == 0.0 && maxits == 0) {
      eps = 1.0E-6;
   }
   state->eps = eps;
   state->maxits = maxits;
}

// This function sets warm-start mode of the solver: next call to the  solver
// will reuse previous subspace as warm-start  point.  It  can  significantly
// speed-up convergence when you solve many similar eigenproblems.
//
// Inputs:
//     State       -   solver structure
//     UseWarmStart-   either True or False
// ALGLIB: Copyright 12.11.2017 by Sergey Bochkanov
// API: void eigsubspacesetwarmstart(const eigsubspacestate &state, const bool usewarmstart);
void eigsubspacesetwarmstart(eigsubspacestate *state, bool usewarmstart) {
   ae_assert(!state->running, "EigSubspaceSetWarmStart: solver is already running");
   state->usewarmstart = usewarmstart;
}

// Clears request fileds (to be sure that we don't forgot to clear something)
static void evd_clearrfields(eigsubspacestate *state) {
   state->requesttype = -1;
   state->requestsize = -1;
}

// This  function  initiates  out-of-core  mode  of  subspace eigensolver. It
// should be used in conjunction with other out-of-core-related functions  of
// this subspackage in a loop like below:
//
// > alglib.eigsubspaceoocstart(state)
// > while alglib.eigsubspaceooccontinue(state) do
// >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
// >     alglib.eigsubspaceoocgetrequestdata(state, out X)
// >     [calculate  Y=A*X, with X=R^NxM]
// >     alglib.eigsubspaceoocsendresult(state, in Y)
// > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
//
// Inputs:
//     State       -   solver object
//     MType       -   matrix type:
//                     * 0 for real  symmetric  matrix  (solver  assumes that
//                       matrix  being   processed  is  symmetric;  symmetric
//                       direct eigensolver is used for  smaller  subproblems
//                       arising during solution of larger "full" task)
//                     Future versions of ALGLIB may  introduce  support  for
//                     other  matrix   types;   for   now,   only   symmetric
//                     eigenproblems are supported.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspaceoocstart(const eigsubspacestate &state, const ae_int_t mtype);
void eigsubspaceoocstart(eigsubspacestate *state, ae_int_t mtype) {
   ae_assert(!state->running, "EigSubspaceStart: solver is already running");
   ae_assert(mtype == 0, "EigSubspaceStart: incorrect mtype parameter");
   state->PQ = -1;
   evd_clearrfields(state);
   state->running = true;
   state->matrixtype = mtype;
}

// This function performs subspace iteration  in  the  out-of-core  mode.  It
// should be used in conjunction with other out-of-core-related functions  of
// this subspackage in a loop like below:
//
// > alglib.eigsubspaceoocstart(state)
// > while alglib.eigsubspaceooccontinue(state) do
// >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
// >     alglib.eigsubspaceoocgetrequestdata(state, out X)
// >     [calculate  Y=A*X, with X=R^NxM]
// >     alglib.eigsubspaceoocsendresult(state, in Y)
// > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: bool eigsubspaceooccontinue(const eigsubspacestate &state);
bool eigsubspaceooccontinue(eigsubspacestate *state) {
   bool result;
   ae_assert(state->running, "EigSubspaceContinue: solver is not running");
   result = eigsubspaceiteration(state);
   state->running = result;
   return result;
}

// This function is used to retrieve information  about  out-of-core  request
// sent by solver to user code: request type (current version  of  the solver
// sends only requests for matrix-matrix products) and request size (size  of
// the matrices being multiplied).
//
// This function returns just request metrics; in order  to  get contents  of
// the matrices being multiplied, use eigsubspaceoocgetrequestdata().
//
// It should be used in conjunction with other out-of-core-related  functions
// of this subspackage in a loop like below:
//
// > alglib.eigsubspaceoocstart(state)
// > while alglib.eigsubspaceooccontinue(state) do
// >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
// >     alglib.eigsubspaceoocgetrequestdata(state, out X)
// >     [calculate  Y=A*X, with X=R^NxM]
// >     alglib.eigsubspaceoocsendresult(state, in Y)
// > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
//
// Inputs:
//     State           -   solver running in out-of-core mode
//
// Outputs:
//     RequestType     -   type of the request to process:
//                         * 0 - for matrix-matrix product A*X, with A  being
//                           NxN matrix whose eigenvalues/vectors are needed,
//                           and X being NxREQUESTSIZE one which is  returned
//                           by the eigsubspaceoocgetrequestdata().
//     RequestSize     -   size of the X matrix (number of columns),  usually
//                         it is several times larger than number of  vectors
//                         K requested by user.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspaceoocgetrequestinfo(const eigsubspacestate &state, ae_int_t &requesttype, ae_int_t &requestsize);
void eigsubspaceoocgetrequestinfo(eigsubspacestate *state, ae_int_t *requesttype, ae_int_t *requestsize) {
   *requesttype = 0;
   *requestsize = 0;
   ae_assert(state->running, "EigSubspaceOOCGetRequestInfo: solver is not running");
   *requesttype = state->requesttype;
   *requestsize = state->requestsize;
}

// This function is used to retrieve information  about  out-of-core  request
// sent by solver to user code: matrix X (array[N,RequestSize) which have  to
// be multiplied by out-of-core matrix A in a product A*X.
//
// This function returns just request data; in order to get size of  the data
// prior to processing requestm, use eigsubspaceoocgetrequestinfo().
//
// It should be used in conjunction with other out-of-core-related  functions
// of this subspackage in a loop like below:
//
// > alglib.eigsubspaceoocstart(state)
// > while alglib.eigsubspaceooccontinue(state) do
// >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
// >     alglib.eigsubspaceoocgetrequestdata(state, out X)
// >     [calculate  Y=A*X, with X=R^NxM]
// >     alglib.eigsubspaceoocsendresult(state, in Y)
// > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
//
// Inputs:
//     State           -   solver running in out-of-core mode
//     X               -   possibly  preallocated   storage;  reallocated  if
//                         needed, left unchanged, if large enough  to  store
//                         request data.
//
// Outputs:
//     X               -   array[N,RequestSize] or larger, leading  rectangle
//                         is filled with dense matrix X.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspaceoocgetrequestdata(const eigsubspacestate &state, real_2d_array &x);
void eigsubspaceoocgetrequestdata(eigsubspacestate *state, RMatrix *x) {
   ae_int_t i;
   ae_int_t j;
   ae_assert(state->running, "EigSubspaceOOCGetRequestInfo: solver is not running");
   matrixsetlengthatleast(x, state->n, state->requestsize);
   for (i = 0; i < state->n; i++) {
      for (j = 0; j < state->requestsize; j++) {
         x->xyR[i][j] = state->x.xyR[i][j];
      }
   }
}

// This function is used to send user reply to out-of-core  request  sent  by
// solver. Usually it is product A*X for returned by solver matrix X.
//
// It should be used in conjunction with other out-of-core-related  functions
// of this subspackage in a loop like below:
//
// > alglib.eigsubspaceoocstart(state)
// > while alglib.eigsubspaceooccontinue(state) do
// >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
// >     alglib.eigsubspaceoocgetrequestdata(state, out X)
// >     [calculate  Y=A*X, with X=R^NxM]
// >     alglib.eigsubspaceoocsendresult(state, in Y)
// > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
//
// Inputs:
//     State           -   solver running in out-of-core mode
//     AX              -   array[N,RequestSize] or larger, leading  rectangle
//                         is filled with product A*X.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspaceoocsendresult(const eigsubspacestate &state, const real_2d_array &ax);
void eigsubspaceoocsendresult(eigsubspacestate *state, RMatrix *ax) {
   ae_int_t i;
   ae_int_t j;
   ae_assert(state->running, "EigSubspaceOOCGetRequestInfo: solver is not running");
   for (i = 0; i < state->n; i++) {
      for (j = 0; j < state->requestsize; j++) {
         state->ax.xyR[i][j] = ax->xyR[i][j];
      }
   }
}

// This  function  finalizes out-of-core  mode  of  subspace eigensolver.  It
// should be used in conjunction with other out-of-core-related functions  of
// this subspackage in a loop like below:
//
// > alglib.eigsubspaceoocstart(state)
// > while alglib.eigsubspaceooccontinue(state) do
// >     alglib.eigsubspaceoocgetrequestinfo(state, out RequestType, out M)
// >     alglib.eigsubspaceoocgetrequestdata(state, out X)
// >     [calculate  Y=A*X, with X=R^NxM]
// >     alglib.eigsubspaceoocsendresult(state, in Y)
// > alglib.eigsubspaceoocstop(state, out W, out Z, out Report)
//
// Inputs:
//     State       -   solver state
//
// Outputs:
//     W           -   array[K], depending on solver settings:
//                     * top  K  eigenvalues ordered  by  descending   -   if
//                       eigenvectors are returned in Z
//                     * zeros - if invariant subspace is returned in Z
//     Z           -   array[N,K], depending on solver settings either:
//                     * matrix of eigenvectors found
//                     * orthogonal basis of K-dimensional invariant subspace
//     Rep         -   report with additional parameters
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspaceoocstop(const eigsubspacestate &state, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep);
void eigsubspaceoocstop(eigsubspacestate *state, RVector *w, RMatrix *z, eigsubspacereport *rep) {
   ae_int_t n;
   ae_int_t k;
   ae_int_t i;
   ae_int_t j;
   SetVector(w);
   SetMatrix(z);
   SetObj(eigsubspacereport, rep);
   ae_assert(!state->running, "EigSubspaceStop: solver is still running");
   n = state->n;
   k = state->k;
   ae_vector_set_length(w, k);
   ae_matrix_set_length(z, n, k);
   for (i = 0; i < k; i++) {
      w->xR[i] = state->rw.xR[i];
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
         z->xyR[i][j] = state->rq.xyR[i][j];
      }
   }
   rep->iterationscount = state->repiterationscount;
}

// This  function runs eigensolver for dense NxN symmetric matrix A, given by
// upper or lower triangle.
//
// This function can not process nonsymmetric matrices.
//
// Inputs:
//     State       -   solver state
//     A           -   array[N,N], symmetric NxN matrix given by one  of  its
//                     triangles
//     IsUpper     -   whether upper or lower triangle of  A  is  given  (the
//                     other one is not referenced at all).
//
// Outputs:
//     W           -   array[K], top  K  eigenvalues ordered  by   descending
//                     of their absolute values
//     Z           -   array[N,K], matrix of eigenvectors found
//     Rep         -   report with additional parameters
//
// NOTE: internally this function allocates a copy of NxN dense A. You should
//       take it into account when working with very large matrices occupying
//       almost all RAM.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspacesolvedenses(const eigsubspacestate &state, const real_2d_array &a, const bool isupper, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep);
void eigsubspacesolvedenses(eigsubspacestate *state, RMatrix *a, bool isupper, RVector *w, RMatrix *z, eigsubspacereport *rep) {
   ae_frame _frame_block;
   ae_int_t n;
   ae_int_t m;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   double v;
   ae_frame_make(&_frame_block);
   SetVector(w);
   SetMatrix(z);
   SetObj(eigsubspacereport, rep);
   NewMatrix(acopy, 0, 0, DT_REAL);
   ae_assert(!state->running, "EigSubspaceSolveDenseS: solver is still running");
   n = state->n;
// Allocate copy of A, copy one triangle to another
   ae_matrix_set_length(&acopy, n, n);
   for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
         if (isupper) {
            v = a->xyR[i][j];
         } else {
            v = a->xyR[j][i];
         }
         acopy.xyR[i][j] = v;
         acopy.xyR[j][i] = v;
      }
   }
// Start iterations
   state->matrixtype = 0;
   state->PQ = -1;
   evd_clearrfields(state);
   while (eigsubspaceiteration(state)) {
   // Calculate A*X with RMatrixGEMM
      ae_assert(state->requesttype == 0, "EigSubspaceSolveDense: integrity check failed");
      ae_assert(state->requestsize > 0, "EigSubspaceSolveDense: integrity check failed");
      m = state->requestsize;
      rmatrixgemm(n, m, n, 1.0, &acopy, 0, 0, 0, &state->x, 0, 0, 0, 0.0, &state->ax, 0, 0);
   }
   k = state->k;
   ae_vector_set_length(w, k);
   ae_matrix_set_length(z, n, k);
   for (i = 0; i < k; i++) {
      w->xR[i] = state->rw.xR[i];
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
         z->xyR[i][j] = state->rq.xyR[i][j];
      }
   }
   rep->iterationscount = state->repiterationscount;
   ae_frame_leave();
}

// This  function runs eigensolver for dense NxN symmetric matrix A, given by
// upper or lower triangle.
//
// This function can not process nonsymmetric matrices.
//
// Inputs:
//     State       -   solver state
//     A           -   NxN symmetric matrix given by one of its triangles
//     IsUpper     -   whether upper or lower triangle of  A  is  given  (the
//                     other one is not referenced at all).
//
// Outputs:
//     W           -   array[K], top  K  eigenvalues ordered  by   descending
//                     of their absolute values
//     Z           -   array[N,K], matrix of eigenvectors found
//     Rep         -   report with additional parameters
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
// API: void eigsubspacesolvesparses(const eigsubspacestate &state, const sparsematrix &a, const bool isupper, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep);
void eigsubspacesolvesparses(eigsubspacestate *state, sparsematrix *a, bool isupper, RVector *w, RMatrix *z, eigsubspacereport *rep) {
   ae_int_t n;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   SetVector(w);
   SetMatrix(z);
   SetObj(eigsubspacereport, rep);
   ae_assert(!state->running, "EigSubspaceSolveSparseS: solver is still running");
   n = state->n;
   state->matrixtype = 0;
   state->PQ = -1;
   evd_clearrfields(state);
   while (eigsubspaceiteration(state)) {
      ae_assert(state->requesttype == 0, "EigSubspaceSolveDense: integrity check failed");
      ae_assert(state->requestsize > 0, "EigSubspaceSolveDense: integrity check failed");
      sparsesmm(a, isupper, &state->x, state->requestsize, &state->ax);
   }
   k = state->k;
   ae_vector_set_length(w, k);
   ae_matrix_set_length(z, n, k);
   for (i = 0; i < k; i++) {
      w->xR[i] = state->rw.xR[i];
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
         z->xyR[i][j] = state->rq.xyR[i][j];
      }
   }
   rep->iterationscount = state->repiterationscount;
}

// Internal r-comm function.
// ALGLIB: Copyright 16.01.2017 by Sergey Bochkanov
bool eigsubspaceiteration(eigsubspacestate *state) {
   AutoS ae_int_t n;
   AutoS ae_int_t nwork;
   AutoS ae_int_t k;
   AutoS ae_int_t cnt;
   AutoS ae_int_t i;
   AutoS ae_int_t i1;
   AutoS ae_int_t j;
   AutoS double vv;
   AutoS double v;
   AutoS ae_int_t convcnt;
// Manually threaded two-way signalling.
// Locals are set arbitrarily the first time around and are retained between pauses and subsequent resumes.
// A Spawn occurs when the routine is (re-)started.
// A Pause sends an event signal and waits for a response with data before carrying out the matching Resume.
// An Exit sends an exit signal indicating the end of the process.
   if (state->PQ >= 0) switch (state->PQ) {
      case 0: goto Resume0;
      default: goto Exit;
   }
Spawn:
   cnt = -909;
   i1 = 255;
   j = 74;
   convcnt = -788;
   vv = 809;
   v = 205;
   n = state->n;
   k = state->k;
   nwork = state->nwork;
// Initialize RNG. Deterministic initialization (with fixed
// seed) is required because we need deterministic behavior
// of the entire solver.
   hqrndseed(453, 463664, &state->rs);
// Prepare iteration
// Initialize QNew with random orthogonal matrix (or reuse its previous value).
   state->repiterationscount = 0;
   matrixsetlengthatleast(&state->qcur, nwork, n);
   matrixsetlengthatleast(&state->qnew, nwork, n);
   matrixsetlengthatleast(&state->znew, nwork, n);
   vectorsetlengthatleast(&state->wcur, nwork);
   vectorsetlengthatleast(&state->wprev, nwork);
   vectorsetlengthatleast(&state->wrank, nwork);
   matrixsetlengthatleast(&state->x, n, nwork);
   matrixsetlengthatleast(&state->ax, n, nwork);
   matrixsetlengthatleast(&state->rq, n, k);
   vectorsetlengthatleast(&state->rw, k);
   matrixsetlengthatleast(&state->rz, nwork, k);
   matrixsetlengthatleast(&state->r, nwork, nwork);
   for (i = 0; i < nwork; i++) {
      state->wprev.xR[i] = -1.0;
   }
   if (!state->usewarmstart || state->firstcall) {
   // Use Q0 (either no warm start request, or warm start was
   // requested by user - but it is first call).
   //
      if (state->firstcall) {
      // First call, generate Q0
         for (i = 0; i < nwork; i++) {
            for (j = 0; j < n; j++) {
               state->znew.xyR[i][j] = hqrnduniformr(&state->rs) - 0.5;
            }
         }
         rmatrixlq(&state->znew, nwork, n, &state->tau);
         rmatrixlqunpackq(&state->znew, nwork, n, &state->tau, nwork, &state->q0);
         state->firstcall = false;
      }
      rmatrixcopy(nwork, n, &state->q0, 0, 0, &state->qnew, 0, 0);
   }
// Start iteration
   state->repiterationscount = 0;
   convcnt = 0;
   while ((state->maxits == 0 || state->repiterationscount < state->maxits) && convcnt < evd_stepswithintol) {
   // Update QCur := QNew
   //
   // Calculate A*Q'
      rmatrixcopy(nwork, n, &state->qnew, 0, 0, &state->qcur, 0, 0);
      rmatrixtranspose(nwork, n, &state->qcur, 0, 0, &state->x, 0, 0);
      evd_clearrfields(state);
      state->requesttype = 0;
      state->requestsize = nwork;
      state->PQ = 0; goto Pause; Resume0:
   // Perform Rayleigh-Ritz step to estimate convergence of diagonal eigenvalues
      if (state->eps > 0.0) {
         ae_assert(state->matrixtype == 0, "integrity check failed");
         matrixsetlengthatleast(&state->r, nwork, nwork);
         rmatrixgemm(nwork, nwork, n, 1.0, &state->qcur, 0, 0, 0, &state->ax, 0, 0, 0, 0.0, &state->r, 0, 0);
         if (!smatrixevd(&state->r, nwork, 0, true, &state->wcur, &state->dummy)) {
            ae_assert(false, "EigSubspace: direct eigensolver failed to converge");
         }
         for (j = 0; j < nwork; j++) {
            state->wrank.xR[j] = fabs(state->wcur.xR[j]);
         }
         rankxuntied(&state->wrank, nwork, &state->buf);
         v = 0.0;
         vv = 0.0;
         for (j = 0; j < nwork; j++) {
            if (state->wrank.xR[j] >= (double)(nwork - k)) {
               v = rmax2(v, fabs(state->wcur.xR[j] - state->wprev.xR[j]));
               vv = rmax2(vv, fabs(state->wcur.xR[j]));
            }
         }
         if (vv == 0.0) {
            vv = 1.0;
         }
         if (v <= state->eps * vv) {
            convcnt++;
         } else {
            convcnt = 0;
         }
         for (j = 0; j < nwork; j++) {
            state->wprev.xR[j] = state->wcur.xR[j];
         }
      }
   // QR renormalization and update of QNew
      rmatrixtranspose(n, nwork, &state->ax, 0, 0, &state->znew, 0, 0);
      rmatrixlq(&state->znew, nwork, n, &state->tau);
      rmatrixlqunpackq(&state->znew, nwork, n, &state->tau, nwork, &state->qnew);
   // Update iteration index
      state->repiterationscount++;
   }
// Perform Rayleigh-Ritz step: find true eigenpairs in NWork-dimensional
// subspace.
   ae_assert(state->matrixtype == 0, "integrity check failed");
   ae_assert(state->eigenvectorsneeded == 1, "Assertion failed");
   rmatrixgemm(nwork, nwork, n, 1.0, &state->qcur, 0, 0, 0, &state->ax, 0, 0, 0, 0.0, &state->r, 0, 0);
   if (!smatrixevd(&state->r, nwork, 1, true, &state->tw, &state->tz)) {
      ae_assert(false, "EigSubspace: direct eigensolver failed to converge");
   }
// Reorder eigenpairs according to their absolute magnitude, select
// K top ones. This reordering algorithm is very inefficient and has
// O(NWork*K) running time, but it is still faster than other parts
// of the solver, so we may use it.
//
// Then, we transform RZ to RQ (full N-dimensional representation).
// After this part is done, RW and RQ contain solution.
   for (j = 0; j < nwork; j++) {
      state->wrank.xR[j] = fabs(state->tw.xR[j]);
   }
   rankxuntied(&state->wrank, nwork, &state->buf);
   cnt = 0;
   for (i = nwork - 1; i >= nwork - k; i--) {
      for (i1 = 0; i1 < nwork; i1++) {
         if (state->wrank.xR[i1] == (double)i) {
            ae_assert(cnt < k, "EigSubspace: integrity check failed");
            state->rw.xR[cnt] = state->tw.xR[i1];
            for (j = 0; j < nwork; j++) {
               state->rz.xyR[j][cnt] = state->tz.xyR[j][i1];
            }
            cnt++;
         }
      }
   }
   ae_assert(cnt == k, "EigSubspace: integrity check failed");
   rmatrixgemm(n, k, nwork, 1.0, &state->qcur, 0, 0, 1, &state->rz, 0, 0, 0, 0.0, &state->rq, 0, 0);
Exit:
   state->PQ = -1;
   return false;
Pause:
   return true;
}

// Finding the eigenvalues and eigenvectors of a symmetric matrix
//
// The algorithm finds eigen pairs of a symmetric matrix by reducing it to
// tridiagonal form and using the QL/QR algorithm.
//
// Inputs:
//     A       -   symmetric matrix which is given by its upper or lower
//                 triangular part.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     ZNeeded -   flag controlling whether the eigenvectors are needed or not.
//                 If ZNeeded is equal to:
//                  * 0, the eigenvectors are not returned;
//                  * 1, the eigenvectors are returned.
//     IsUpper -   storage format.
//
// Outputs:
//     D       -   eigenvalues in ascending order.
//                 Array whose index ranges within [0..N-1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains the eigenvectors.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//                 The eigenvectors are stored in the matrix columns.
//
// Result:
//     True, if the algorithm has converged.
//     False, if the algorithm hasn't converged (rare case).
// ALGLIB: Copyright 2005-2008 by Sergey Bochkanov
// API: bool smatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, real_2d_array &z);
bool smatrixevd(RMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, RVector *d, RMatrix *z) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(d);
   SetMatrix(z);
   NewVector(tau, 0, DT_REAL);
   NewVector(e, 0, DT_REAL);
   ae_assert(zneeded == 0 || zneeded == 1, "SMatrixEVD: incorrect ZNeeded");
   smatrixtd(a, n, isupper, &tau, d, &e);
   if (zneeded == 1) {
      smatrixtdunpackq(a, n, isupper, &tau, z);
   }
   result = smatrixtdevd(d, &e, n, zneeded, z);
   ae_frame_leave();
   return result;
}

// Subroutine for finding the eigenvalues (and eigenvectors) of  a  symmetric
// matrix  in  a  given half open interval (A, B] by using  a  bisection  and
// inverse iteration
//
// Inputs:
//     A       -   symmetric matrix which is given by its upper or lower
//                 triangular part. Array [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     ZNeeded -   flag controlling whether the eigenvectors are needed or not.
//                 If ZNeeded is equal to:
//                  * 0, the eigenvectors are not returned;
//                  * 1, the eigenvectors are returned.
//     IsUpperA -  storage format of matrix A.
//     B1, B2 -    half open interval (B1, B2] to search eigenvalues in.
//
// Outputs:
//     M       -   number of eigenvalues found in a given half-interval (M >= 0).
//     W       -   array of the eigenvalues found.
//                 Array whose index ranges within [0..M-1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains eigenvectors.
//                 Array whose indexes range within [0..N-1, 0..M-1].
//                 The eigenvectors are stored in the matrix columns.
//
// Result:
//     True, if successful. M contains the number of eigenvalues in the given
//     half-interval (could be equal to 0), W contains the eigenvalues,
//     Z contains the eigenvectors (if needed).
//
//     False, if the bisection method subroutine wasn't able to find the
//     eigenvalues in the given interval or if the inverse iteration subroutine
//     wasn't able to find all the corresponding eigenvectors.
//     In that case, the eigenvalues and eigenvectors are not returned,
//     M is equal to 0.
// ALGLIB: Copyright 07.01.2006 by Sergey Bochkanov
// API: bool smatrixevdr(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, real_2d_array &z);
bool smatrixevdr(RMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, double b1, double b2, ae_int_t *m, RVector *w, RMatrix *z) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *m = 0;
   SetVector(w);
   SetMatrix(z);
   NewVector(tau, 0, DT_REAL);
   NewVector(e, 0, DT_REAL);
   ae_assert(zneeded == 0 || zneeded == 1, "SMatrixTDEVDR: incorrect ZNeeded");
   smatrixtd(a, n, isupper, &tau, w, &e);
   if (zneeded == 1) {
      smatrixtdunpackq(a, n, isupper, &tau, z);
   }
   result = smatrixtdevdr(w, &e, n, zneeded, b1, b2, m, z);
   ae_frame_leave();
   return result;
}

// Subroutine for finding the eigenvalues and  eigenvectors  of  a  symmetric
// matrix with given indexes by using bisection and inverse iteration methods.
//
// Inputs:
//     A       -   symmetric matrix which is given by its upper or lower
//                 triangular part. Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     ZNeeded -   flag controlling whether the eigenvectors are needed or not.
//                 If ZNeeded is equal to:
//                  * 0, the eigenvectors are not returned;
//                  * 1, the eigenvectors are returned.
//     IsUpperA -  storage format of matrix A.
//     I1, I2 -    index interval for searching (from I1 to I2).
//                 0 <= I1 <= I2 <= N-1.
//
// Outputs:
//     W       -   array of the eigenvalues found.
//                 Array whose index ranges within [0..I2-I1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains eigenvectors.
//                 Array whose indexes range within [0..N-1, 0..I2-I1].
//                 In that case, the eigenvectors are stored in the matrix columns.
//
// Result:
//     True, if successful. W contains the eigenvalues, Z contains the
//     eigenvectors (if needed).
//
//     False, if the bisection method subroutine wasn't able to find the
//     eigenvalues in the given interval or if the inverse iteration subroutine
//     wasn't able to find all the corresponding eigenvectors.
//     In that case, the eigenvalues and eigenvectors are not returned.
// ALGLIB: Copyright 07.01.2006 by Sergey Bochkanov
// API: bool smatrixevdi(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, real_2d_array &z);
bool smatrixevdi(RMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, ae_int_t i1, ae_int_t i2, RVector *w, RMatrix *z) {
   ae_frame _frame_block;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(w);
   SetMatrix(z);
   NewVector(tau, 0, DT_REAL);
   NewVector(e, 0, DT_REAL);
   ae_assert(zneeded == 0 || zneeded == 1, "SMatrixEVDI: incorrect ZNeeded");
   smatrixtd(a, n, isupper, &tau, w, &e);
   if (zneeded == 1) {
      smatrixtdunpackq(a, n, isupper, &tau, z);
   }
   result = smatrixtdevdi(w, &e, n, zneeded, i1, i2, z);
   ae_frame_leave();
   return result;
}

// Finding the eigenvalues and eigenvectors of a Hermitian matrix
//
// The algorithm finds eigen pairs of a Hermitian matrix by  reducing  it  to
// real tridiagonal form and using the QL/QR algorithm.
//
// Inputs:
//     A       -   Hermitian matrix which is given  by  its  upper  or  lower
//                 triangular part.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     IsUpper -   storage format.
//     ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
//                 not. If ZNeeded is equal to:
//                  * 0, the eigenvectors are not returned;
//                  * 1, the eigenvectors are returned.
//
// Outputs:
//     D       -   eigenvalues in ascending order.
//                 Array whose index ranges within [0..N-1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains the eigenvectors.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//                 The eigenvectors are stored in the matrix columns.
//
// Result:
//     True, if the algorithm has converged.
//     False, if the algorithm hasn't converged (rare case).
//
// Note:
//     eigenvectors of Hermitian matrix are defined up to  multiplication  by
//     a complex number L, such that |L|=1.
// ALGLIB: Copyright 2005, 2007 March 23 by Sergey Bochkanov
// API: bool hmatrixevd(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, complex_2d_array &z);
bool hmatrixevd(CMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, RVector *d, CMatrix *z) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(d);
   SetMatrix(z);
   NewVector(tau, 0, DT_COMPLEX);
   NewVector(e, 0, DT_REAL);
   NewMatrix(t, 0, 0, DT_REAL);
   NewMatrix(qz, 0, 0, DT_REAL);
   NewMatrix(q, 0, 0, DT_COMPLEX);
   ae_assert(zneeded == 0 || zneeded == 1, "HermitianEVD: incorrect ZNeeded");
// Reduce to tridiagonal form
   hmatrixtd(a, n, isupper, &tau, d, &e);
   if (zneeded == 1) {
      hmatrixtdunpackq(a, n, isupper, &tau, &q);
      zneeded = 2;
   }
// TDEVD
   result = smatrixtdevd(d, &e, n, zneeded, &t);
// Eigenvectors are needed
// Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
   if (result && zneeded != 0) {
      ae_matrix_set_length(z, n, n);
      ae_matrix_set_length(&qz, n, 2 * n);
   // Calculate Re(Q)*T
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            qz.xyR[i][j] = q.xyC[i][j].x;
         }
      }
      rmatrixgemm(n, n, n, 1.0, &qz, 0, 0, 0, &t, 0, 0, 0, 0.0, &qz, 0, n);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            z->xyC[i][j].x = qz.xyR[i][n + j];
         }
      }
   // Calculate Im(Q)*T
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            qz.xyR[i][j] = q.xyC[i][j].y;
         }
      }
      rmatrixgemm(n, n, n, 1.0, &qz, 0, 0, 0, &t, 0, 0, 0, 0.0, &qz, 0, n);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            z->xyC[i][j].y = qz.xyR[i][n + j];
         }
      }
   }
   ae_frame_leave();
   return result;
}

// Subroutine for finding the eigenvalues (and eigenvectors) of  a  Hermitian
// matrix  in  a  given half-interval (A, B] by using a bisection and inverse
// iteration
//
// Inputs:
//     A       -   Hermitian matrix which is given  by  its  upper  or  lower
//                 triangular  part.  Array  whose   indexes   range   within
//                 [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
//                 not. If ZNeeded is equal to:
//                  * 0, the eigenvectors are not returned;
//                  * 1, the eigenvectors are returned.
//     IsUpperA -  storage format of matrix A.
//     B1, B2 -    half-interval (B1, B2] to search eigenvalues in.
//
// Outputs:
//     M       -   number of eigenvalues found in a given half-interval, M >= 0
//     W       -   array of the eigenvalues found.
//                 Array whose index ranges within [0..M-1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains eigenvectors.
//                 Array whose indexes range within [0..N-1, 0..M-1].
//                 The eigenvectors are stored in the matrix columns.
//
// Result:
//     True, if successful. M contains the number of eigenvalues in the given
//     half-interval (could be equal to 0), W contains the eigenvalues,
//     Z contains the eigenvectors (if needed).
//
//     False, if the bisection method subroutine  wasn't  able  to  find  the
//     eigenvalues  in  the  given  interval  or  if  the  inverse  iteration
//     subroutine  wasn't  able  to  find all the corresponding eigenvectors.
//     In that case, the eigenvalues and eigenvectors are not returned, M  is
//     equal to 0.
//
// Note:
//     eigen vectors of Hermitian matrix are defined up to multiplication  by
//     a complex number L, such as |L|=1.
// ALGLIB: Copyright 07.01.2006, 24.03.2007 by Sergey Bochkanov
// API: bool hmatrixevdr(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, complex_2d_array &z);
bool hmatrixevdr(CMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, double b1, double b2, ae_int_t *m, RVector *w, CMatrix *z) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t k;
   double v;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   *m = 0;
   SetVector(w);
   SetMatrix(z);
   NewMatrix(q, 0, 0, DT_COMPLEX);
   NewMatrix(t, 0, 0, DT_REAL);
   NewVector(tau, 0, DT_COMPLEX);
   NewVector(e, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   ae_assert(zneeded == 0 || zneeded == 1, "HermitianEigenValuesAndVectorsInInterval: incorrect ZNeeded");
// Reduce to tridiagonal form
   hmatrixtd(a, n, isupper, &tau, w, &e);
   if (zneeded == 1) {
      hmatrixtdunpackq(a, n, isupper, &tau, &q);
      zneeded = 2;
   }
// Bisection and inverse iteration
   result = smatrixtdevdr(w, &e, n, zneeded, b1, b2, m, &t);
// Eigenvectors are needed
// Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
   if (result && zneeded != 0 && *m != 0) {
      ae_vector_set_length(&work, *m);
      ae_matrix_set_length(z, n, *m);
      for (i = 0; i < n; i++) {
      // Calculate real part
         for (k = 0; k < *m; k++) {
            work.xR[k] = 0.0;
         }
         for (k = 0; k < n; k++) {
            v = q.xyC[i][k].x;
            ae_v_addd(work.xR, 1, t.xyR[k], 1, *m, v);
         }
         for (k = 0; k < *m; k++) {
            z->xyC[i][k].x = work.xR[k];
         }
      // Calculate imaginary part
         for (k = 0; k < *m; k++) {
            work.xR[k] = 0.0;
         }
         for (k = 0; k < n; k++) {
            v = q.xyC[i][k].y;
            ae_v_addd(work.xR, 1, t.xyR[k], 1, *m, v);
         }
         for (k = 0; k < *m; k++) {
            z->xyC[i][k].y = work.xR[k];
         }
      }
   }
   ae_frame_leave();
   return result;
}

// Subroutine for finding the eigenvalues and  eigenvectors  of  a  Hermitian
// matrix with given indexes by using bisection and inverse iteration methods
//
// Inputs:
//     A       -   Hermitian matrix which is given  by  its  upper  or  lower
//                 triangular part.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
//                 not. If ZNeeded is equal to:
//                  * 0, the eigenvectors are not returned;
//                  * 1, the eigenvectors are returned.
//     IsUpperA -  storage format of matrix A.
//     I1, I2 -    index interval for searching (from I1 to I2).
//                 0 <= I1 <= I2 <= N-1.
//
// Outputs:
//     W       -   array of the eigenvalues found.
//                 Array whose index ranges within [0..I2-I1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains eigenvectors.
//                 Array whose indexes range within [0..N-1, 0..I2-I1].
//                 In  that  case,  the eigenvectors are stored in the matrix
//                 columns.
//
// Result:
//     True, if successful. W contains the eigenvalues, Z contains the
//     eigenvectors (if needed).
//
//     False, if the bisection method subroutine  wasn't  able  to  find  the
//     eigenvalues  in  the  given  interval  or  if  the  inverse  iteration
//     subroutine wasn't able to find  all  the  corresponding  eigenvectors.
//     In that case, the eigenvalues and eigenvectors are not returned.
//
// Note:
//     eigen vectors of Hermitian matrix are defined up to multiplication  by
//     a complex number L, such as |L|=1.
// ALGLIB: Copyright 07.01.2006, 24.03.2007 by Sergey Bochkanov
// API: bool hmatrixevdi(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, complex_2d_array &z);
bool hmatrixevdi(CMatrix *a, ae_int_t n, ae_int_t zneeded, bool isupper, ae_int_t i1, ae_int_t i2, RVector *w, CMatrix *z) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t k;
   double v;
   ae_int_t m;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(w);
   SetMatrix(z);
   NewMatrix(q, 0, 0, DT_COMPLEX);
   NewMatrix(t, 0, 0, DT_REAL);
   NewVector(tau, 0, DT_COMPLEX);
   NewVector(e, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   ae_assert(zneeded == 0 || zneeded == 1, "HermitianEigenValuesAndVectorsByIndexes: incorrect ZNeeded");
// Reduce to tridiagonal form
   hmatrixtd(a, n, isupper, &tau, w, &e);
   if (zneeded == 1) {
      hmatrixtdunpackq(a, n, isupper, &tau, &q);
      zneeded = 2;
   }
// Bisection and inverse iteration
   result = smatrixtdevdi(w, &e, n, zneeded, i1, i2, &t);
// Eigenvectors are needed
// Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
   m = i2 - i1 + 1;
   if (result && zneeded != 0) {
      ae_vector_set_length(&work, m);
      ae_matrix_set_length(z, n, m);
      for (i = 0; i < n; i++) {
      // Calculate real part
         for (k = 0; k < m; k++) {
            work.xR[k] = 0.0;
         }
         for (k = 0; k < n; k++) {
            v = q.xyC[i][k].x;
            ae_v_addd(work.xR, 1, t.xyR[k], 1, m, v);
         }
         for (k = 0; k < m; k++) {
            z->xyC[i][k].x = work.xR[k];
         }
      // Calculate imaginary part
         for (k = 0; k < m; k++) {
            work.xR[k] = 0.0;
         }
         for (k = 0; k < n; k++) {
            v = q.xyC[i][k].y;
            ae_v_addd(work.xR, 1, t.xyR[k], 1, m, v);
         }
         for (k = 0; k < m; k++) {
            z->xyC[i][k].y = work.xR[k];
         }
      }
   }
   ae_frame_leave();
   return result;
}

// DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
//    [  A   B  ]
//    [  B   C  ].
// On return, RT1 is the eigenvalue of larger absolute value, and RT2
// is the eigenvalue of smaller absolute value.
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
static void evd_tdevde2(double a, double b, double c, double *rt1, double *rt2) {
   double ab;
   double acmn;
   double acmx;
   double adf;
   double df;
   double rt;
   double sm;
   double tb;
   *rt1 = 0;
   *rt2 = 0;
   sm = a + c;
   df = a - c;
   adf = fabs(df);
   tb = b + b;
   ab = fabs(tb);
   if (fabs(a) > fabs(c)) {
      acmx = a;
      acmn = c;
   } else {
      acmx = c;
      acmn = a;
   }
   if (adf > ab) {
      rt = adf * sqrt(1 + ae_sqr(ab / adf));
   } else {
      if (adf < ab) {
         rt = ab * sqrt(1 + ae_sqr(adf / ab));
      } else {
      // Includes case AB=ADF=0
         rt = ab * sqrt(2.0);
      }
   }
   if (sm < 0.0) {
      *rt1 = 0.5 * (sm - rt);
   // Order of execution important.
   // To get fully accurate smaller eigenvalue,
   // next line needs to be executed in higher precision.
      *rt2 = acmx / (*rt1) * acmn - b / (*rt1) * b;
   } else {
      if (sm > 0.0) {
         *rt1 = 0.5 * (sm + rt);
      // Order of execution important.
      // To get fully accurate smaller eigenvalue,
      // next line needs to be executed in higher precision.
         *rt2 = acmx / (*rt1) * acmn - b / (*rt1) * b;
      } else {
      // Includes case RT1 = RT2 = 0
         *rt1 = 0.5 * rt;
         *rt2 = -0.5 * rt;
      }
   }
}

// DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
//
//    [  A   B  ]
//    [  B   C  ].
//
// On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
// eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
// eigenvector for RT1, giving the decomposition
//
//    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
//    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
//
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
static void evd_tdevdev2(double a, double b, double c, double *rt1, double *rt2, double *cs1, double *sn1) {
   ae_int_t sgn1;
   ae_int_t sgn2;
   double ab;
   double acmn;
   double acmx;
   double acs;
   double adf;
   double cs;
   double ct;
   double df;
   double rt;
   double sm;
   double tb;
   double tn;
   *rt1 = 0;
   *rt2 = 0;
   *cs1 = 0;
   *sn1 = 0;
// Compute the eigenvalues
   sm = a + c;
   df = a - c;
   adf = fabs(df);
   tb = b + b;
   ab = fabs(tb);
   if (fabs(a) > fabs(c)) {
      acmx = a;
      acmn = c;
   } else {
      acmx = c;
      acmn = a;
   }
   if (adf > ab) {
      rt = adf * sqrt(1 + ae_sqr(ab / adf));
   } else {
      if (adf < ab) {
         rt = ab * sqrt(1 + ae_sqr(adf / ab));
      } else {
      // Includes case AB=ADF=0
         rt = ab * sqrt(2.0);
      }
   }
   if (sm < 0.0) {
      *rt1 = 0.5 * (sm - rt);
      sgn1 = -1;
   // Order of execution important.
   // To get fully accurate smaller eigenvalue,
   // next line needs to be executed in higher precision.
      *rt2 = acmx / (*rt1) * acmn - b / (*rt1) * b;
   } else {
      if (sm > 0.0) {
         *rt1 = 0.5 * (sm + rt);
         sgn1 = 1;
      // Order of execution important.
      // To get fully accurate smaller eigenvalue,
      // next line needs to be executed in higher precision.
         *rt2 = acmx / (*rt1) * acmn - b / (*rt1) * b;
      } else {
      // Includes case RT1 = RT2 = 0
         *rt1 = 0.5 * rt;
         *rt2 = -0.5 * rt;
         sgn1 = 1;
      }
   }
// Compute the eigenvector
   if (df >= 0.0) {
      cs = df + rt;
      sgn2 = 1;
   } else {
      cs = df - rt;
      sgn2 = -1;
   }
   acs = fabs(cs);
   if (acs > ab) {
      ct = -tb / cs;
      *sn1 = 1 / sqrt(1 + ct * ct);
      *cs1 = ct * (*sn1);
   } else {
      if (ab == 0.0) {
         *cs1 = 1.0;
         *sn1 = 0.0;
      } else {
         tn = -cs / tb;
         *cs1 = 1 / sqrt(1 + tn * tn);
         *sn1 = tn * (*cs1);
      }
   }
   if (sgn1 == sgn2) {
      tn = *cs1;
      *cs1 = -*sn1;
      *sn1 = tn;
   }
}

// Internal routine
static double evd_tdevdpythag(double a, double b) {
   double result;
   if (fabs(a) < fabs(b)) {
      result = fabs(b) * sqrt(1 + ae_sqr(a / b));
   } else {
      result = fabs(a) * sqrt(1 + ae_sqr(b / a));
   }
   return result;
}

// Internal routine
static double evd_tdevdextsign(double a, double b) {
   double result;
   if (b >= 0.0) {
      result = fabs(a);
   } else {
      result = -fabs(a);
   }
   return result;
}

static bool evd_tridiagonalevd(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z) {
   ae_frame _frame_block;
   ae_int_t maxit;
   ae_int_t i;
   ae_int_t ii;
   ae_int_t iscale;
   ae_int_t j;
   ae_int_t jtot;
   ae_int_t k;
   ae_int_t t;
   ae_int_t l;
   ae_int_t l1;
   ae_int_t lend;
   ae_int_t lendm1;
   ae_int_t lendp1;
   ae_int_t lendsv;
   ae_int_t lm1;
   ae_int_t lsv;
   ae_int_t m;
   ae_int_t mm1;
   ae_int_t nm1;
   ae_int_t nmaxit;
   ae_int_t tmpint;
   double anorm;
   double b;
   double c;
   double eps;
   double eps2;
   double f;
   double g;
   double p;
   double r;
   double rt1;
   double rt2;
   double s;
   double safmax;
   double safmin;
   double ssfmax;
   double ssfmin;
   double tst;
   double tmp;
   bool gotoflag;
   ae_int_t zrows;
   bool wastranspose;
   bool result;
   ae_frame_make(&_frame_block);
   DupVector(e);
   NewVector(work1, 0, DT_REAL);
   NewVector(work2, 0, DT_REAL);
   NewVector(workc, 0, DT_REAL);
   NewVector(works, 0, DT_REAL);
   NewVector(wtemp, 0, DT_REAL);
   ae_assert(zneeded >= 0 && zneeded <= 3, "TridiagonalEVD: Incorrent ZNeeded");
// Quick return if possible
   if (zneeded < 0 || zneeded > 3) {
      result = false;
      ae_frame_leave();
      return result;
   }
   result = true;
   if (n == 0) {
      ae_frame_leave();
      return result;
   }
   if (n == 1) {
      if (zneeded == 2 || zneeded == 3) {
         ae_matrix_set_length(z, 1 + 1, 1 + 1);
         z->xyR[1][1] = 1.0;
      }
      ae_frame_leave();
      return result;
   }
   maxit = 30;
// Initialize arrays
   ae_vector_set_length(&wtemp, n + 1);
   ae_vector_set_length(&work1, n);
   ae_vector_set_length(&work2, n);
   ae_vector_set_length(&workc, n + 1);
   ae_vector_set_length(&works, n + 1);
// Determine the unit roundoff and over/underflow thresholds.
   eps = ae_machineepsilon;
   eps2 = ae_sqr(eps);
   safmin = ae_minrealnumber;
   safmax = ae_maxrealnumber;
   ssfmax = sqrt(safmax) / 3;
   ssfmin = sqrt(safmin) / eps2;
// Prepare Z
//
// Here we are using transposition to get rid of column operations
//
   wastranspose = false;
   zrows = 0;
   if (zneeded == 1) {
      zrows = n;
   }
   if (zneeded == 2) {
      zrows = n;
   }
   if (zneeded == 3) {
      zrows = 1;
   }
   if (zneeded == 1) {
      wastranspose = true;
      inplacetranspose(z, 1, n, 1, n, &wtemp);
   }
   if (zneeded == 2) {
      wastranspose = true;
      ae_matrix_set_length(z, n + 1, n + 1);
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= n; j++) {
            if (i == j) {
               z->xyR[i][j] = 1.0;
            } else {
               z->xyR[i][j] = 0.0;
            }
         }
      }
   }
   if (zneeded == 3) {
      wastranspose = false;
      ae_matrix_set_length(z, 1 + 1, n + 1);
      for (j = 1; j <= n; j++) {
         if (j == 1) {
            z->xyR[1][j] = 1.0;
         } else {
            z->xyR[1][j] = 0.0;
         }
      }
   }
   nmaxit = n * maxit;
   jtot = 0;
// Determine where the matrix splits and choose QL or QR iteration
// for each block, according to whether top or bottom diagonal
// element is smaller.
   l1 = 1;
   nm1 = n - 1;
   while (true) {
      if (l1 > n) {
         break;
      }
      if (l1 > 1) {
         e->xR[l1 - 1] = 0.0;
      }
      gotoflag = false;
      m = l1;
      if (l1 <= nm1) {
         for (m = l1; m <= nm1; m++) {
            tst = fabs(e->xR[m]);
            if (tst == 0.0) {
               gotoflag = true;
               break;
            }
            if (tst <= sqrt(fabs(d->xR[m])) * sqrt(fabs(d->xR[m + 1])) * eps) {
               e->xR[m] = 0.0;
               gotoflag = true;
               break;
            }
         }
      }
      if (!gotoflag) {
         m = n;
      }
   // label 30:
      l = l1;
      lsv = l;
      lend = m;
      lendsv = lend;
      l1 = m + 1;
      if (lend == l) {
         continue;
      }
   // Scale submatrix in rows and columns L to LEND
      if (l == lend) {
         anorm = fabs(d->xR[l]);
      } else {
         anorm = rmax2(fabs(d->xR[l]) + fabs(e->xR[l]), fabs(e->xR[lend - 1]) + fabs(d->xR[lend]));
         for (i = l + 1; i < lend; i++) {
            anorm = rmax2(anorm, fabs(d->xR[i]) + fabs(e->xR[i]) + fabs(e->xR[i - 1]));
         }
      }
      iscale = 0;
      if (anorm == 0.0) {
         continue;
      }
      if (anorm > ssfmax) {
         iscale = 1;
         tmp = ssfmax / anorm;
         tmpint = lend - 1;
         ae_v_muld(&d->xR[l], 1, lend - l + 1, tmp);
         ae_v_muld(&e->xR[l], 1, tmpint - l + 1, tmp);
      }
      if (anorm < ssfmin) {
         iscale = 2;
         tmp = ssfmin / anorm;
         tmpint = lend - 1;
         ae_v_muld(&d->xR[l], 1, lend - l + 1, tmp);
         ae_v_muld(&e->xR[l], 1, tmpint - l + 1, tmp);
      }
   // Choose between QL and QR iteration
      if (fabs(d->xR[lend]) < fabs(d->xR[l])) {
         lend = lsv;
         l = lendsv;
      }
      if (lend > l) {
      // QL Iteration
      //
      // Look for small subdiagonal element.
         while (true) {
            gotoflag = false;
            if (l != lend) {
               lendm1 = lend - 1;
               for (m = l; m <= lendm1; m++) {
                  tst = ae_sqr(fabs(e->xR[m]));
                  if (tst <= eps2 * fabs(d->xR[m]) * fabs(d->xR[m + 1]) + safmin) {
                     gotoflag = true;
                     break;
                  }
               }
            }
            if (!gotoflag) {
               m = lend;
            }
            if (m < lend) {
               e->xR[m] = 0.0;
            }
            p = d->xR[l];
            if (m != l) {
            // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
            // to compute its eigensystem.
               if (m == l + 1) {
                  if (zneeded > 0) {
                     evd_tdevdev2(d->xR[l], e->xR[l], d->xR[l + 1], &rt1, &rt2, &c, &s);
                     work1.xR[l] = c;
                     work2.xR[l] = s;
                     workc.xR[1] = work1.xR[l];
                     works.xR[1] = work2.xR[l];
                     if (!wastranspose) {
                        applyrotationsfromtheright(false, 1, zrows, l, l + 1, &workc, &works, z, &wtemp);
                     } else {
                        applyrotationsfromtheleft(false, l, l + 1, 1, zrows, &workc, &works, z, &wtemp);
                     }
                  } else {
                     evd_tdevde2(d->xR[l], e->xR[l], d->xR[l + 1], &rt1, &rt2);
                  }
                  d->xR[l] = rt1;
                  d->xR[l + 1] = rt2;
                  e->xR[l] = 0.0;
                  l += 2;
                  if (l <= lend) {
                     continue;
                  }
               // GOTO 140
                  break;
               }
               if (jtot == nmaxit) {
               // GOTO 140
                  break;
               }
               jtot++;
            // Form shift.
               g = (d->xR[l + 1] - p) / (2 * e->xR[l]);
               r = evd_tdevdpythag(g, 1.0);
               g = d->xR[m] - p + e->xR[l] / (g + evd_tdevdextsign(r, g));
               s = 1.0;
               c = 1.0;
               p = 0.0;
            // Inner loop
               mm1 = m - 1;
               for (i = mm1; i >= l; i--) {
                  f = s * e->xR[i];
                  b = c * e->xR[i];
                  generaterotation(g, f, &c, &s, &r);
                  if (i != m - 1) {
                     e->xR[i + 1] = r;
                  }
                  g = d->xR[i + 1] - p;
                  r = (d->xR[i] - g) * s + 2 * c * b;
                  p = s * r;
                  d->xR[i + 1] = g + p;
                  g = c * r - b;
               // If eigenvectors are desired, then save rotations.
                  if (zneeded > 0) {
                     work1.xR[i] = c;
                     work2.xR[i] = -s;
                  }
               }
            // If eigenvectors are desired, then apply saved rotations.
               if (zneeded > 0) {
                  for (i = l; i < m; i++) {
                     workc.xR[i - l + 1] = work1.xR[i];
                     works.xR[i - l + 1] = work2.xR[i];
                  }
                  if (!wastranspose) {
                     applyrotationsfromtheright(false, 1, zrows, l, m, &workc, &works, z, &wtemp);
                  } else {
                     applyrotationsfromtheleft(false, l, m, 1, zrows, &workc, &works, z, &wtemp);
                  }
               }
               d->xR[l] -= p;
               e->xR[l] = g;
               continue;
            }
         // Eigenvalue found.
            d->xR[l] = p;
            l++;
            if (l <= lend) {
               continue;
            }
            break;
         }
      } else {
      // QR Iteration
      //
      // Look for small superdiagonal element.
         while (true) {
            gotoflag = false;
            if (l != lend) {
               lendp1 = lend + 1;
               for (m = l; m >= lendp1; m--) {
                  tst = ae_sqr(fabs(e->xR[m - 1]));
                  if (tst <= eps2 * fabs(d->xR[m]) * fabs(d->xR[m - 1]) + safmin) {
                     gotoflag = true;
                     break;
                  }
               }
            }
            if (!gotoflag) {
               m = lend;
            }
            if (m > lend) {
               e->xR[m - 1] = 0.0;
            }
            p = d->xR[l];
            if (m != l) {
            // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
            // to compute its eigensystem.
               if (m == l - 1) {
                  if (zneeded > 0) {
                     evd_tdevdev2(d->xR[l - 1], e->xR[l - 1], d->xR[l], &rt1, &rt2, &c, &s);
                     work1.xR[m] = c;
                     work2.xR[m] = s;
                     workc.xR[1] = c;
                     works.xR[1] = s;
                     if (!wastranspose) {
                        applyrotationsfromtheright(true, 1, zrows, l - 1, l, &workc, &works, z, &wtemp);
                     } else {
                        applyrotationsfromtheleft(true, l - 1, l, 1, zrows, &workc, &works, z, &wtemp);
                     }
                  } else {
                     evd_tdevde2(d->xR[l - 1], e->xR[l - 1], d->xR[l], &rt1, &rt2);
                  }
                  d->xR[l - 1] = rt1;
                  d->xR[l] = rt2;
                  e->xR[l - 1] = 0.0;
                  l -= 2;
                  if (l >= lend) {
                     continue;
                  }
                  break;
               }
               if (jtot == nmaxit) {
                  break;
               }
               jtot++;
            // Form shift.
               g = (d->xR[l - 1] - p) / (2 * e->xR[l - 1]);
               r = evd_tdevdpythag(g, 1.0);
               g = d->xR[m] - p + e->xR[l - 1] / (g + evd_tdevdextsign(r, g));
               s = 1.0;
               c = 1.0;
               p = 0.0;
            // Inner loop
               lm1 = l - 1;
               for (i = m; i <= lm1; i++) {
                  f = s * e->xR[i];
                  b = c * e->xR[i];
                  generaterotation(g, f, &c, &s, &r);
                  if (i != m) {
                     e->xR[i - 1] = r;
                  }
                  g = d->xR[i] - p;
                  r = (d->xR[i + 1] - g) * s + 2 * c * b;
                  p = s * r;
                  d->xR[i] = g + p;
                  g = c * r - b;
               // If eigenvectors are desired, then save rotations.
                  if (zneeded > 0) {
                     work1.xR[i] = c;
                     work2.xR[i] = s;
                  }
               }
            // If eigenvectors are desired, then apply saved rotations.
               if (zneeded > 0) {
                  for (i = m; i < l; i++) {
                     workc.xR[i - m + 1] = work1.xR[i];
                     works.xR[i - m + 1] = work2.xR[i];
                  }
                  if (!wastranspose) {
                     applyrotationsfromtheright(true, 1, zrows, m, l, &workc, &works, z, &wtemp);
                  } else {
                     applyrotationsfromtheleft(true, m, l, 1, zrows, &workc, &works, z, &wtemp);
                  }
               }
               d->xR[l] -= p;
               e->xR[lm1] = g;
               continue;
            }
         // Eigenvalue found.
            d->xR[l] = p;
            l--;
            if (l >= lend) {
               continue;
            }
            break;
         }
      }
   // Undo scaling if necessary
      if (iscale == 1) {
         tmp = anorm / ssfmax;
         tmpint = lendsv - 1;
         ae_v_muld(&d->xR[lsv], 1, lendsv - lsv + 1, tmp);
         ae_v_muld(&e->xR[lsv], 1, tmpint - lsv + 1, tmp);
      }
      if (iscale == 2) {
         tmp = anorm / ssfmin;
         tmpint = lendsv - 1;
         ae_v_muld(&d->xR[lsv], 1, lendsv - lsv + 1, tmp);
         ae_v_muld(&e->xR[lsv], 1, tmpint - lsv + 1, tmp);
      }
   // Check for no convergence to an eigenvalue after a total
   // of N*MAXIT iterations.
      if (jtot >= nmaxit) {
         result = false;
         if (wastranspose) {
            inplacetranspose(z, 1, n, 1, n, &wtemp);
         }
         ae_frame_leave();
         return result;
      }
   }
// Order eigenvalues and eigenvectors.
   if (zneeded == 0) {
   // Sort
      if (n == 1) {
         ae_frame_leave();
         return result;
      }
      if (n == 2) {
         if (d->xR[1] > d->xR[2]) {
            swapr(&d->xR[1], &d->xR[2]);
         }
         ae_frame_leave();
         return result;
      }
      i = 2;
      do {
         t = i;
         while (t != 1) {
            k = t / 2;
            if (d->xR[k] >= d->xR[t]) {
               t = 1;
            } else {
               swapr(&d->xR[k], &d->xR[t]);
               t = k;
            }
         }
         i++;
      } while (i <= n);
      i = n - 1;
      do {
         swapr(&d->xR[i + 1], &d->xR[1]);
         t = 1;
         while (t != 0) {
            k = 2 * t;
            if (k > i) {
               t = 0;
            } else {
               if (k < i) {
                  if (d->xR[k + 1] > d->xR[k]) {
                     k++;
                  }
               }
               if (d->xR[t] >= d->xR[k]) {
                  t = 0;
               } else {
                  swapr(&d->xR[k], &d->xR[t]);
                  t = k;
               }
            }
         }
         i--;
      } while (i >= 1);
   } else {
   // Use Selection Sort to minimize swaps of eigenvectors
      for (ii = 2; ii <= n; ii++) {
         i = ii - 1;
         k = i;
         p = d->xR[i];
         for (j = ii; j <= n; j++) {
            if (d->xR[j] < p) {
               k = j;
               p = d->xR[j];
            }
         }
         if (k != i) {
            d->xR[k] = d->xR[i];
            d->xR[i] = p;
            if (wastranspose) {
               ae_v_move(&wtemp.xR[1], 1, &z->xyR[i][1], 1, n);
               ae_v_move(&z->xyR[i][1], 1, &z->xyR[k][1], 1, n);
               ae_v_move(&z->xyR[k][1], 1, &wtemp.xR[1], 1, n);
            } else {
               ae_v_move(&wtemp.xR[1], 1, &z->xyR[1][i], z->stride, zrows);
               ae_v_move(&z->xyR[1][i], z->stride, &z->xyR[1][k], z->stride, zrows);
               ae_v_move(&z->xyR[1][k], z->stride, &wtemp.xR[1], 1, zrows);
            }
         }
      }
      if (wastranspose) {
         inplacetranspose(z, 1, n, 1, n, &wtemp);
      }
   }
   ae_frame_leave();
   return result;
}

static void evd_internaldlaebz(ae_int_t ijob, ae_int_t nitmax, ae_int_t n, ae_int_t mmax, ae_int_t minp, double abstol, double reltol, double pivmin, RVector *d, RVector *e, RVector *e2, ZVector *nval, RMatrix *ab, RVector *c, ae_int_t *mout, ZMatrix *nab, RVector *work, ZVector *iwork, ae_int_t *info) {
   ae_int_t itmp1;
   ae_int_t itmp2;
   ae_int_t j;
   ae_int_t ji;
   ae_int_t jit;
   ae_int_t jp;
   ae_int_t kf;
   ae_int_t kfnew;
   ae_int_t kl;
   ae_int_t klnew;
   double tmp1;
   double tmp2;
   *mout = 0;
   *info = 0;
   *info = 0;
   if (ijob < 1 || ijob > 3) {
      *info = -1;
      return;
   }
// Initialize NAB
   if (ijob == 1) {
   // Compute the number of eigenvalues in the initial intervals.
      *mout = 0;
   // DIR$ NOVECTOR
      for (ji = 1; ji <= minp; ji++) {
         for (jp = 1; jp <= 2; jp++) {
            tmp1 = d->xR[1] - ab->xyR[ji][jp];
            if (SmallR(tmp1, pivmin)) {
               tmp1 = -pivmin;
            }
            nab->xyZ[ji][jp] = 0;
            if (tmp1 <= 0.0) {
               nab->xyZ[ji][jp] = 1;
            }
            for (j = 2; j <= n; j++) {
               tmp1 = d->xR[j] - e2->xR[j - 1] / tmp1 - ab->xyR[ji][jp];
               if (SmallR(tmp1, pivmin)) {
                  tmp1 = -pivmin;
               }
               if (tmp1 <= 0.0) {
                  nab->xyZ[ji][jp]++;
               }
            }
         }
         *mout += nab->xyZ[ji][2] - nab->xyZ[ji][1];
      }
      return;
   }
// Initialize for loop
//
// KF and KL have the following meaning:
//   Intervals 1,...,KF-1 have converged.
//   Intervals KF,...,KL  still need to be refined.
   kf = 1;
   kl = minp;
// If IJOB=2, initialize C.
// If IJOB=3, use the user-supplied starting point.
   if (ijob == 2) {
      for (ji = 1; ji <= minp; ji++) {
         c->xR[ji] = 0.5 * (ab->xyR[ji][1] + ab->xyR[ji][2]);
      }
   }
// Iteration loop
   for (jit = 1; jit <= nitmax; jit++) {
   // Loop over intervals
   //
   //
   // Serial Version of the loop
      klnew = kl;
      for (ji = kf; ji <= kl; ji++) {
      // Compute N(w), the number of eigenvalues less than w
         tmp1 = c->xR[ji];
         tmp2 = d->xR[1] - tmp1;
         itmp1 = 0;
         if (tmp2 <= pivmin) {
            itmp1 = 1;
            tmp2 = rmin2(tmp2, -pivmin);
         }
      // A series of compiler directives to defeat vectorization
      // for the next loop
      //
      // *$PL$ CMCHAR=' '
      // CDIR$          NEXTSCALAR
      // C$DIR          SCALAR
      // CDIR$          NEXT SCALAR
      // CVD$L          NOVECTOR
      // CDEC$          NOVECTOR
      // CVD$           NOVECTOR
      // *VDIR          NOVECTOR
      // *VOCL          LOOP,SCALAR
      // CIBM           PREFER SCALAR
      // *$PL$ CMCHAR='*'
         for (j = 2; j <= n; j++) {
            tmp2 = d->xR[j] - e2->xR[j - 1] / tmp2 - tmp1;
            if (tmp2 <= pivmin) {
               itmp1++;
               tmp2 = rmin2(tmp2, -pivmin);
            }
         }
         if (ijob <= 2) {
         // IJOB=2: Choose all intervals containing eigenvalues.
         //
         // Insure that N(w) is monotone
            itmp1 = imin2(nab->xyZ[ji][2], imax2(nab->xyZ[ji][1], itmp1));
         // Update the Queue -- add intervals if both halves
         // contain eigenvalues.
            if (itmp1 == nab->xyZ[ji][2]) {
            // No eigenvalue in the upper interval:
            // just use the lower interval.
               ab->xyR[ji][2] = tmp1;
            } else {
               if (itmp1 == nab->xyZ[ji][1]) {
               // No eigenvalue in the lower interval:
               // just use the upper interval.
                  ab->xyR[ji][1] = tmp1;
               } else {
                  if (klnew < mmax) {
                  // Eigenvalue in both intervals -- add upper to queue.
                     klnew++;
                     ab->xyR[klnew][2] = ab->xyR[ji][2];
                     nab->xyZ[klnew][2] = nab->xyZ[ji][2];
                     ab->xyR[klnew][1] = tmp1;
                     nab->xyZ[klnew][1] = itmp1;
                     ab->xyR[ji][2] = tmp1;
                     nab->xyZ[ji][2] = itmp1;
                  } else {
                     *info = mmax + 1;
                     return;
                  }
               }
            }
         } else {
         // IJOB=3: Binary search.  Keep only the interval
         // containing  w  s.t. N(w) = NVAL
            if (itmp1 <= nval->xZ[ji]) {
               ab->xyR[ji][1] = tmp1;
               nab->xyZ[ji][1] = itmp1;
            }
            if (itmp1 >= nval->xZ[ji]) {
               ab->xyR[ji][2] = tmp1;
               nab->xyZ[ji][2] = itmp1;
            }
         }
      }
      kl = klnew;
   // Check for convergence
      kfnew = kf;
      for (ji = kf; ji <= kl; ji++) {
         tmp1 = fabs(ab->xyR[ji][2] - ab->xyR[ji][1]);
         tmp2 = rmax2(fabs(ab->xyR[ji][2]), fabs(ab->xyR[ji][1]));
         if (tmp1 < rmax2(abstol, rmax2(pivmin, reltol * tmp2)) || nab->xyZ[ji][1] >= nab->xyZ[ji][2]) {
         // Converged -- Swap with position KFNEW,
         // then increment KFNEW
            if (ji > kfnew) {
               swapr(&ab->xyR[ji][1], &ab->xyR[kfnew][1]);
               swapr(&ab->xyR[ji][2], &ab->xyR[kfnew][2]);
               swapi(&nab->xyZ[ji][1], &nab->xyZ[kfnew][1]);
               swapi(&nab->xyZ[ji][2], &nab->xyZ[kfnew][2]);
               if (ijob == 3) {
                  swapi(&nval->xZ[ji], &nval->xZ[kfnew]);
               }
            }
            kfnew++;
         }
      }
      kf = kfnew;
   // Choose Midpoints
      for (ji = kf; ji <= kl; ji++) {
         c->xR[ji] = 0.5 * (ab->xyR[ji][1] + ab->xyR[ji][2]);
      }
   // If no more intervals to refine, quit.
      if (kf > kl) {
         break;
      }
   }
// Converged
   *info = imax2(kl + 1 - kf, 0);
   *mout = kl;
}

static bool evd_internalbisectioneigenvalues(RVector *d, RVector *e, ae_int_t n, ae_int_t irange, ae_int_t iorder, double vl, double vu, ae_int_t il, ae_int_t iu, double abstol, RVector *w, ae_int_t *m, ae_int_t *nsplit, ZVector *iblock, ZVector *isplit, ae_int_t *errorcode) {
   ae_frame _frame_block;
   double fudge;
   double relfac;
   bool ncnvrg;
   bool toofew;
   ae_int_t ib;
   ae_int_t ibegin;
   ae_int_t idiscl;
   ae_int_t idiscu;
   ae_int_t ie;
   ae_int_t iend;
   ae_int_t iinfo;
   ae_int_t im;
   ae_int_t iin;
   ae_int_t ioff;
   ae_int_t iout;
   ae_int_t itmax;
   ae_int_t iw;
   ae_int_t iwoff;
   ae_int_t j;
   ae_int_t jb;
   ae_int_t jdisc;
   ae_int_t je;
   ae_int_t nwl;
   ae_int_t nwu;
   double atoli;
   double bnorm;
   double gl;
   double gu;
   double pivmin;
   double rtoli;
   double safemn;
   double tmp1;
   double tmp2;
   double tnorm;
   double ulp;
   double wkill;
   double wl;
   double wlu;
   double wu;
   double wul;
   double scalefactor;
   double t;
   ae_int_t tmpi;
   bool result;
   ae_frame_make(&_frame_block);
   DupVector(d);
   DupVector(e);
   SetVector(w);
   *m = 0;
   *nsplit = 0;
   SetVector(iblock);
   SetVector(isplit);
   *errorcode = 0;
   NewVector(idumma, 0, DT_INT);
   NewVector(work, 0, DT_REAL);
   NewVector(iwork, 0, DT_INT);
   NewVector(ia1s2, 0, DT_INT);
   NewVector(ra1s2, 0, DT_REAL);
   NewMatrix(ra1s2x2, 0, 0, DT_REAL);
   NewMatrix(ia1s2x2, 0, 0, DT_INT);
   NewVector(ra1siin, 0, DT_REAL);
   NewVector(ra2siin, 0, DT_REAL);
   NewVector(ra3siin, 0, DT_REAL);
   NewVector(ra4siin, 0, DT_REAL);
   NewMatrix(ra1siinx2, 0, 0, DT_REAL);
   NewMatrix(ia1siinx2, 0, 0, DT_INT);
   NewVector(iworkspace, 0, DT_INT);
   NewVector(rworkspace, 0, DT_REAL);
// Quick return if possible
   *m = 0;
   if (n == 0) {
      result = true;
      ae_frame_leave();
      return result;
   }
// Get machine constants
// NB is the minimum vector length for vector bisection, or 0
// if only scalar is to be done.
   fudge = 2.0;
   relfac = 2.0;
   safemn = ae_minrealnumber;
   ulp = 2 * ae_machineepsilon;
   rtoli = ulp * relfac;
   ae_vector_set_length(&idumma, 1 + 1);
   ae_vector_set_length(&work, 4 * n + 1);
   ae_vector_set_length(&iwork, 3 * n + 1);
   ae_vector_set_length(w, n + 1);
   ae_vector_set_length(iblock, n + 1);
   ae_vector_set_length(isplit, n + 1);
   ae_vector_set_length(&ia1s2, 2 + 1);
   ae_vector_set_length(&ra1s2, 2 + 1);
   ae_matrix_set_length(&ra1s2x2, 2 + 1, 2 + 1);
   ae_matrix_set_length(&ia1s2x2, 2 + 1, 2 + 1);
   ae_vector_set_length(&ra1siin, n + 1);
   ae_vector_set_length(&ra2siin, n + 1);
   ae_vector_set_length(&ra3siin, n + 1);
   ae_vector_set_length(&ra4siin, n + 1);
   ae_matrix_set_length(&ra1siinx2, n + 1, 2 + 1);
   ae_matrix_set_length(&ia1siinx2, n + 1, 2 + 1);
   ae_vector_set_length(&iworkspace, n + 1);
   ae_vector_set_length(&rworkspace, n + 1);
// these initializers are not really necessary,
// but without them compiler complains about uninitialized locals
   wlu = 0.0;
   wul = 0.0;
// Check for Errors
   result = false;
   *errorcode = 0;
   if (irange <= 0 || irange >= 4) {
      *errorcode = -4;
   }
   if (iorder <= 0 || iorder >= 3) {
      *errorcode = -5;
   }
   if (n < 0) {
      *errorcode = -3;
   }
   if (irange == 2 && vl >= vu) {
      *errorcode = -6;
   }
   if (irange == 3 && (il < 1 || il > imax2(1, n))) {
      *errorcode = -8;
   }
   if (irange == 3 && (iu < imin2(n, il) || iu > n)) {
      *errorcode = -9;
   }
   if (*errorcode != 0) {
      ae_frame_leave();
      return result;
   }
// Initialize error flags
   ncnvrg = false;
   toofew = false;
// Simplifications:
   if (irange == 3 && il == 1 && iu == n) {
      irange = 1;
   }
// Special Case when N=1
   if (n == 1) {
      *nsplit = 1;
      isplit->xZ[1] = 1;
      if (irange == 2 && (vl >= d->xR[1] || vu < d->xR[1])) {
         *m = 0;
      } else {
         w->xR[1] = d->xR[1];
         iblock->xZ[1] = 1;
         *m = 1;
      }
      result = true;
      ae_frame_leave();
      return result;
   }
// Scaling
   t = fabs(d->xR[n]);
   for (j = 1; j < n; j++) {
      t = rmax2(t, fabs(d->xR[j]));
      t = rmax2(t, fabs(e->xR[j]));
   }
   scalefactor = 1.0;
   if (t != 0.0) {
      if (t > sqrt(sqrt(ae_minrealnumber)) * sqrt(ae_maxrealnumber)) {
         scalefactor = t;
      }
      if (t < sqrt(sqrt(ae_maxrealnumber)) * sqrt(ae_minrealnumber)) {
         scalefactor = t;
      }
      for (j = 1; j < n; j++) {
         d->xR[j] /= scalefactor;
         e->xR[j] /= scalefactor;
      }
      d->xR[n] /= scalefactor;
   }
// Compute Splitting Points
   *nsplit = 1;
   work.xR[n] = 0.0;
   pivmin = 1.0;
   for (j = 2; j <= n; j++) {
      tmp1 = ae_sqr(e->xR[j - 1]);
      if (fabs(d->xR[j] * d->xR[j - 1]) * ae_sqr(ulp) + safemn > tmp1) {
         isplit->xZ[*nsplit] = j - 1;
         ++*nsplit;
         work.xR[j - 1] = 0.0;
      } else {
         work.xR[j - 1] = tmp1;
         pivmin = rmax2(pivmin, tmp1);
      }
   }
   isplit->xZ[*nsplit] = n;
   pivmin *= safemn;
// Compute Interval and ATOLI
   if (irange == 3) {
   // RANGE='I': Compute the interval containing eigenvalues
   //     IL through IU.
   //
   // Compute Gershgorin interval for entire (split) matrix
   // and use it as the initial interval
      gu = d->xR[1];
      gl = d->xR[1];
      tmp1 = 0.0;
      for (j = 1; j < n; j++) {
         tmp2 = sqrt(work.xR[j]);
         gu = rmax2(gu, d->xR[j] + tmp1 + tmp2);
         gl = rmin2(gl, d->xR[j] - tmp1 - tmp2);
         tmp1 = tmp2;
      }
      gu = rmax2(gu, d->xR[n] + tmp1);
      gl = rmin2(gl, d->xR[n] - tmp1);
      tnorm = rmax2(fabs(gl), fabs(gu));
      gl -= fudge * tnorm * ulp * n + fudge * 2 * pivmin;
      gu += fudge * tnorm * ulp * n + fudge * pivmin;
   // Compute Iteration parameters
      itmax = CeilZ((log(tnorm + pivmin) - log(pivmin)) / log(2.0)) + 2;
      if (abstol <= 0.0) {
         atoli = ulp * tnorm;
      } else {
         atoli = abstol;
      }
      work.xR[n + 1] = gl;
      work.xR[n + 2] = gl;
      work.xR[n + 3] = gu;
      work.xR[n + 4] = gu;
      work.xR[n + 5] = gl;
      work.xR[n + 6] = gu;
      iwork.xZ[1] = -1;
      iwork.xZ[2] = -1;
      iwork.xZ[3] = n + 1;
      iwork.xZ[4] = n + 1;
      iwork.xZ[5] = il - 1;
      iwork.xZ[6] = iu;
   // Calling DLAEBZ
   //
   // DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,
   //    WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
   //    IWORK, W, IBLOCK, IINFO )
      ia1s2.xZ[1] = iwork.xZ[5];
      ia1s2.xZ[2] = iwork.xZ[6];
      ra1s2.xR[1] = work.xR[n + 5];
      ra1s2.xR[2] = work.xR[n + 6];
      ra1s2x2.xyR[1][1] = work.xR[n + 1];
      ra1s2x2.xyR[2][1] = work.xR[n + 2];
      ra1s2x2.xyR[1][2] = work.xR[n + 3];
      ra1s2x2.xyR[2][2] = work.xR[n + 4];
      ia1s2x2.xyZ[1][1] = iwork.xZ[1];
      ia1s2x2.xyZ[2][1] = iwork.xZ[2];
      ia1s2x2.xyZ[1][2] = iwork.xZ[3];
      ia1s2x2.xyZ[2][2] = iwork.xZ[4];
      evd_internaldlaebz(3, itmax, n, 2, 2, atoli, rtoli, pivmin, d, e, &work, &ia1s2, &ra1s2x2, &ra1s2, &iout, &ia1s2x2, w, iblock, &iinfo);
      iwork.xZ[5] = ia1s2.xZ[1];
      iwork.xZ[6] = ia1s2.xZ[2];
      work.xR[n + 5] = ra1s2.xR[1];
      work.xR[n + 6] = ra1s2.xR[2];
      work.xR[n + 1] = ra1s2x2.xyR[1][1];
      work.xR[n + 2] = ra1s2x2.xyR[2][1];
      work.xR[n + 3] = ra1s2x2.xyR[1][2];
      work.xR[n + 4] = ra1s2x2.xyR[2][2];
      iwork.xZ[1] = ia1s2x2.xyZ[1][1];
      iwork.xZ[2] = ia1s2x2.xyZ[2][1];
      iwork.xZ[3] = ia1s2x2.xyZ[1][2];
      iwork.xZ[4] = ia1s2x2.xyZ[2][2];
      if (iwork.xZ[6] == iu) {
         wl = work.xR[n + 1];
         wlu = work.xR[n + 3];
         nwl = iwork.xZ[1];
         wu = work.xR[n + 4];
         wul = work.xR[n + 2];
         nwu = iwork.xZ[4];
      } else {
         wl = work.xR[n + 2];
         wlu = work.xR[n + 4];
         nwl = iwork.xZ[2];
         wu = work.xR[n + 3];
         wul = work.xR[n + 1];
         nwu = iwork.xZ[3];
      }
      if (nwl < 0 || nwl >= n || nwu < 1 || nwu > n) {
         *errorcode = 4;
         result = false;
         ae_frame_leave();
         return result;
      }
   } else {
   // RANGE='A' or 'V' -- Set ATOLI
      tnorm = rmax2(fabs(d->xR[1]) + fabs(e->xR[1]), fabs(d->xR[n]) + fabs(e->xR[n - 1]));
      for (j = 2; j < n; j++) {
         tnorm = rmax2(tnorm, fabs(d->xR[j]) + fabs(e->xR[j - 1]) + fabs(e->xR[j]));
      }
      if (abstol <= 0.0) {
         atoli = ulp * tnorm;
      } else {
         atoli = abstol;
      }
      if (irange == 2) {
         wl = vl;
         wu = vu;
      } else {
         wl = 0.0;
         wu = 0.0;
      }
   }
// Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
// NWL accumulates the number of eigenvalues .le. WL,
// NWU accumulates the number of eigenvalues .le. WU
   *m = 0;
   iend = 0;
   *errorcode = 0;
   nwl = 0;
   nwu = 0;
   for (jb = 1; jb <= *nsplit; jb++) {
      ioff = iend;
      ibegin = ioff + 1;
      iend = isplit->xZ[jb];
      iin = iend - ioff;
      if (iin == 1) {
      // Special Case -- IIN=1
         if (irange == 1 || wl >= d->xR[ibegin] - pivmin) {
            nwl++;
         }
         if (irange == 1 || wu >= d->xR[ibegin] - pivmin) {
            nwu++;
         }
         if (irange == 1 || wl < d->xR[ibegin] - pivmin && wu >= d->xR[ibegin] - pivmin) {
            ++*m;
            w->xR[*m] = d->xR[ibegin];
            iblock->xZ[*m] = jb;
         }
      } else {
      // General Case -- IIN > 1
      //
      // Compute Gershgorin Interval
      // and use it as the initial interval
         gu = d->xR[ibegin];
         gl = d->xR[ibegin];
         tmp1 = 0.0;
         for (j = ibegin; j < iend; j++) {
            tmp2 = fabs(e->xR[j]);
            gu = rmax2(gu, d->xR[j] + tmp1 + tmp2);
            gl = rmin2(gl, d->xR[j] - tmp1 - tmp2);
            tmp1 = tmp2;
         }
         gu = rmax2(gu, d->xR[iend] + tmp1);
         gl = rmin2(gl, d->xR[iend] - tmp1);
         bnorm = rmax2(fabs(gl), fabs(gu));
         gl -= fudge * bnorm * ulp * iin + fudge * pivmin;
         gu += fudge * bnorm * ulp * iin + fudge * pivmin;
      // Compute ATOLI for the current submatrix
         if (abstol <= 0.0) {
            atoli = ulp * rmax2(fabs(gl), fabs(gu));
         } else {
            atoli = abstol;
         }
         if (irange > 1) {
            if (gu < wl) {
               nwl += iin;
               nwu += iin;
               continue;
            }
            gl = rmax2(gl, wl);
            gu = rmin2(gu, wu);
            if (gl >= gu) {
               continue;
            }
         }
      // Set Up Initial Interval
         work.xR[n + 1] = gl;
         work.xR[n + iin + 1] = gu;
      // Calling DLAEBZ
      //
      // CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
      //    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
      //    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
      //    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
         for (tmpi = 1; tmpi <= iin; tmpi++) {
            ra1siin.xR[tmpi] = d->xR[ibegin - 1 + tmpi];
            if (ibegin - 1 + tmpi < n) {
               ra2siin.xR[tmpi] = e->xR[ibegin - 1 + tmpi];
            }
            ra3siin.xR[tmpi] = work.xR[ibegin - 1 + tmpi];
            ra1siinx2.xyR[tmpi][1] = work.xR[n + tmpi];
            ra1siinx2.xyR[tmpi][2] = work.xR[n + tmpi + iin];
            ra4siin.xR[tmpi] = work.xR[n + 2 * iin + tmpi];
            rworkspace.xR[tmpi] = w->xR[*m + tmpi];
            iworkspace.xZ[tmpi] = iblock->xZ[*m + tmpi];
            ia1siinx2.xyZ[tmpi][1] = iwork.xZ[tmpi];
            ia1siinx2.xyZ[tmpi][2] = iwork.xZ[tmpi + iin];
         }
         evd_internaldlaebz(1, 0, iin, iin, 1, atoli, rtoli, pivmin, &ra1siin, &ra2siin, &ra3siin, &idumma, &ra1siinx2, &ra4siin, &im, &ia1siinx2, &rworkspace, &iworkspace, &iinfo);
         for (tmpi = 1; tmpi <= iin; tmpi++) {
            work.xR[n + tmpi] = ra1siinx2.xyR[tmpi][1];
            work.xR[n + tmpi + iin] = ra1siinx2.xyR[tmpi][2];
            work.xR[n + 2 * iin + tmpi] = ra4siin.xR[tmpi];
            w->xR[*m + tmpi] = rworkspace.xR[tmpi];
            iblock->xZ[*m + tmpi] = iworkspace.xZ[tmpi];
            iwork.xZ[tmpi] = ia1siinx2.xyZ[tmpi][1];
            iwork.xZ[tmpi + iin] = ia1siinx2.xyZ[tmpi][2];
         }
         nwl += iwork.xZ[1];
         nwu += iwork.xZ[iin + 1];
         iwoff = *m - iwork.xZ[1];
      // Compute Eigenvalues
         itmax = CeilZ((log(gu - gl + pivmin) - log(pivmin)) / log(2.0)) + 2;
      // Calling DLAEBZ
      //
      // CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
      //    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
      //    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
      //    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
         for (tmpi = 1; tmpi <= iin; tmpi++) {
            ra1siin.xR[tmpi] = d->xR[ibegin - 1 + tmpi];
            if (ibegin - 1 + tmpi < n) {
               ra2siin.xR[tmpi] = e->xR[ibegin - 1 + tmpi];
            }
            ra3siin.xR[tmpi] = work.xR[ibegin - 1 + tmpi];
            ra1siinx2.xyR[tmpi][1] = work.xR[n + tmpi];
            ra1siinx2.xyR[tmpi][2] = work.xR[n + tmpi + iin];
            ra4siin.xR[tmpi] = work.xR[n + 2 * iin + tmpi];
            rworkspace.xR[tmpi] = w->xR[*m + tmpi];
            iworkspace.xZ[tmpi] = iblock->xZ[*m + tmpi];
            ia1siinx2.xyZ[tmpi][1] = iwork.xZ[tmpi];
            ia1siinx2.xyZ[tmpi][2] = iwork.xZ[tmpi + iin];
         }
         evd_internaldlaebz(2, itmax, iin, iin, 1, atoli, rtoli, pivmin, &ra1siin, &ra2siin, &ra3siin, &idumma, &ra1siinx2, &ra4siin, &iout, &ia1siinx2, &rworkspace, &iworkspace, &iinfo);
         for (tmpi = 1; tmpi <= iin; tmpi++) {
            work.xR[n + tmpi] = ra1siinx2.xyR[tmpi][1];
            work.xR[n + tmpi + iin] = ra1siinx2.xyR[tmpi][2];
            work.xR[n + 2 * iin + tmpi] = ra4siin.xR[tmpi];
            w->xR[*m + tmpi] = rworkspace.xR[tmpi];
            iblock->xZ[*m + tmpi] = iworkspace.xZ[tmpi];
            iwork.xZ[tmpi] = ia1siinx2.xyZ[tmpi][1];
            iwork.xZ[tmpi + iin] = ia1siinx2.xyZ[tmpi][2];
         }
      // Copy Eigenvalues Into W and IBLOCK
      // Use -JB for block number for unconverged eigenvalues.
         for (j = 1; j <= iout; j++) {
            tmp1 = 0.5 * (work.xR[j + n] + work.xR[j + iin + n]);
         // Flag non-convergence.
            if (j > iout - iinfo) {
               ncnvrg = true;
               ib = -jb;
            } else {
               ib = jb;
            }
            for (je = iwork.xZ[j] + 1 + iwoff; je <= iwork.xZ[j + iin] + iwoff; je++) {
               w->xR[je] = tmp1;
               iblock->xZ[je] = ib;
            }
         }
         *m += im;
      }
   }
// If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
// If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
   if (irange == 3) {
      im = 0;
      idiscl = il - 1 - nwl;
      idiscu = nwu - iu;
      if (idiscl > 0 || idiscu > 0) {
         for (je = 1; je <= *m; je++) {
            if (w->xR[je] <= wlu && idiscl > 0) {
               idiscl--;
            } else {
               if (w->xR[je] >= wul && idiscu > 0) {
                  idiscu--;
               } else {
                  im++;
                  w->xR[im] = w->xR[je];
                  iblock->xZ[im] = iblock->xZ[je];
               }
            }
         }
         *m = im;
      }
      if (idiscl > 0 || idiscu > 0) {
      // Code to deal with effects of bad arithmetic:
      // Some low eigenvalues to be discarded are not in (WL,WLU],
      // or high eigenvalues to be discarded are not in (WUL,WU]
      // so just kill off the smallest IDISCL/largest IDISCU
      // eigenvalues, by simply finding the smallest/largest
      // eigenvalue(s).
      //
      // (If N(w) is monotone non-decreasing, this should never
      // happen.)
         if (idiscl > 0) {
            wkill = wu;
            for (jdisc = 1; jdisc <= idiscl; jdisc++) {
               iw = 0;
               for (je = 1; je <= *m; je++) {
                  if (iblock->xZ[je] != 0 && (w->xR[je] < wkill || iw == 0)) {
                     iw = je;
                     wkill = w->xR[je];
                  }
               }
               iblock->xZ[iw] = 0;
            }
         }
         if (idiscu > 0) {
            wkill = wl;
            for (jdisc = 1; jdisc <= idiscu; jdisc++) {
               iw = 0;
               for (je = 1; je <= *m; je++) {
                  if (iblock->xZ[je] != 0 && (w->xR[je] > wkill || iw == 0)) {
                     iw = je;
                     wkill = w->xR[je];
                  }
               }
               iblock->xZ[iw] = 0;
            }
         }
         im = 0;
         for (je = 1; je <= *m; je++) {
            if (iblock->xZ[je] != 0) {
               im++;
               w->xR[im] = w->xR[je];
               iblock->xZ[im] = iblock->xZ[je];
            }
         }
         *m = im;
      }
      if (idiscl < 0 || idiscu < 0) {
         toofew = true;
      }
   }
// If ORDER='B', do nothing -- the eigenvalues are already sorted
//    by block.
// If ORDER='E', sort the eigenvalues from smallest to largest
   if (iorder == 1 && *nsplit > 1) {
      for (je = 1; je < *m; je++) {
         ie = 0;
         tmp1 = w->xR[je];
         for (j = je + 1; j <= *m; j++) {
            if (w->xR[j] < tmp1) {
               ie = j;
               tmp1 = w->xR[j];
            }
         }
         if (ie != 0) {
            w->xR[ie] = w->xR[je];
            w->xR[je] = tmp1;
            swapi(&iblock->xZ[ie], &iblock->xZ[je]);
         }
      }
   }
   for (j = 1; j <= *m; j++) {
      w->xR[j] *= scalefactor;
   }
   *errorcode = 0;
   if (ncnvrg) {
      ++*errorcode;
   }
   if (toofew) {
      *errorcode += 2;
   }
   result = *errorcode == 0;
   ae_frame_leave();
   return result;
}

static void evd_tdininternaldlagtf(ae_int_t n, RVector *a, double lambdav, RVector *b, RVector *c, double tol, RVector *d, ZVector *iin, ae_int_t *info) {
   ae_int_t k;
   double eps;
   double mult;
   double piv1;
   double piv2;
   double scale1;
   double scale2;
   double temp;
   double tl;
   *info = 0;
   *info = 0;
   if (n < 0) {
      *info = -1;
      return;
   }
   if (n == 0) {
      return;
   }
   a->xR[1] -= lambdav;
   iin->xZ[n] = 0;
   if (n == 1) {
      if (a->xR[1] == 0.0) {
         iin->xZ[1] = 1;
      }
      return;
   }
   eps = ae_machineepsilon;
   tl = rmax2(tol, eps);
   scale1 = fabs(a->xR[1]) + fabs(b->xR[1]);
   for (k = 1; k < n; k++) {
      a->xR[k + 1] -= lambdav;
      scale2 = fabs(c->xR[k]) + fabs(a->xR[k + 1]);
      if (k < n - 1) {
         scale2 += fabs(b->xR[k + 1]);
      }
      if (a->xR[k] == 0.0) {
         piv1 = 0.0;
      } else {
         piv1 = fabs(a->xR[k]) / scale1;
      }
      if (c->xR[k] == 0.0) {
         iin->xZ[k] = 0;
         piv2 = 0.0;
         scale1 = scale2;
         if (k < n - 1) {
            d->xR[k] = 0.0;
         }
      } else {
         piv2 = fabs(c->xR[k]) / scale2;
         if (piv2 <= piv1) {
            iin->xZ[k] = 0;
            scale1 = scale2;
            c->xR[k] /= a->xR[k];
            a->xR[k + 1] -= c->xR[k] * b->xR[k];
            if (k < n - 1) {
               d->xR[k] = 0.0;
            }
         } else {
            iin->xZ[k] = 1;
            mult = a->xR[k] / c->xR[k];
            a->xR[k] = c->xR[k];
            temp = a->xR[k + 1];
            a->xR[k + 1] = b->xR[k] - mult * temp;
            if (k < n - 1) {
               d->xR[k] = b->xR[k + 1];
               b->xR[k + 1] = -mult * d->xR[k];
            }
            b->xR[k] = temp;
            c->xR[k] = mult;
         }
      }
      if (rmax2(piv1, piv2) <= tl && iin->xZ[n] == 0) {
         iin->xZ[n] = k;
      }
   }
   if (SmallAtR(a->xR[n], scale1 * tl) && iin->xZ[n] == 0) {
      iin->xZ[n] = n;
   }
}

static void evd_tdininternaldlagts(ae_int_t n, RVector *a, RVector *b, RVector *c, RVector *d, ZVector *iin, RVector *y, double *tol, ae_int_t *info) {
   ae_int_t k;
   double absak;
   double ak;
   double bignum;
   double eps;
   double pert;
   double sfmin;
   double temp;
   *info = 0;
   *info = 0;
   if (n < 0) {
      *info = -1;
      return;
   }
   if (n == 0) {
      return;
   }
   eps = ae_machineepsilon;
   sfmin = ae_minrealnumber;
   bignum = 1 / sfmin;
   if (*tol <= 0.0) {
      *tol = fabs(a->xR[1]);
      if (n > 1) {
         *tol = rmax2(*tol, rmax2(fabs(a->xR[2]), fabs(b->xR[1])));
      }
      for (k = 3; k <= n; k++) {
         *tol = rmax2(*tol, rmax2(fabs(a->xR[k]), rmax2(fabs(b->xR[k - 1]), fabs(d->xR[k - 2]))));
      }
      *tol *= eps;
      if (*tol == 0.0) {
         *tol = eps;
      }
   }
   for (k = 2; k <= n; k++) {
      if (iin->xZ[k - 1] == 0) {
         y->xR[k] -= c->xR[k - 1] * y->xR[k - 1];
      } else {
         temp = y->xR[k - 1];
         y->xR[k - 1] = y->xR[k];
         y->xR[k] = temp - c->xR[k - 1] * y->xR[k];
      }
   }
   for (k = n; k >= 1; k--) {
      if (k < n - 1) {
         temp = y->xR[k] - b->xR[k] * y->xR[k + 1] - d->xR[k] * y->xR[k + 2];
      } else {
         if (k == n - 1) {
            temp = y->xR[k] - b->xR[k] * y->xR[k + 1];
         } else {
            temp = y->xR[k];
         }
      }
      ak = a->xR[k];
      pert = fabs(*tol);
      if (ak < 0.0) {
         pert = -pert;
      }
      while (true) {
         absak = fabs(ak);
         if (absak < 1.0) {
            if (absak < sfmin) {
               if (absak == 0.0 || fabs(temp) * sfmin > absak) {
                  ak += pert;
                  pert *= 2;
                  continue;
               } else {
                  temp *= bignum;
                  ak *= bignum;
               }
            } else {
               if (!SmallAtR(temp, absak * bignum)) {
                  ak += pert;
                  pert *= 2;
                  continue;
               }
            }
         }
         break;
      }
      y->xR[k] = temp / ak;
   }
}

static void evd_internaldstein(ae_int_t n, RVector *d, RVector *e, ae_int_t m, RVector *w, ZVector *iblock, ZVector *isplit, RMatrix *z, ZVector *ifail, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t maxits;
   ae_int_t extra;
   ae_int_t b1;
   ae_int_t blksiz;
   ae_int_t bn;
   ae_int_t gpind;
   ae_int_t i;
   ae_int_t iinfo;
   ae_int_t its;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t jblk;
   ae_int_t jmax;
   ae_int_t nblk;
   ae_int_t nrmchk;
   double dtpcrt;
   double eps;
   double eps1;
   double nrm;
   double onenrm;
   double ortol;
   double pertol;
   double scl;
   double sep;
   double tol;
   double xj;
   double xjm;
   double ztr;
   bool tmpcriterion;
   ae_int_t ti;
   ae_int_t i1;
   double v;
   ae_frame_make(&_frame_block);
   DupVector(e);
   DupVector(w);
   SetMatrix(z);
   SetVector(ifail);
   *info = 0;
   NewVector(work1, 0, DT_REAL);
   NewVector(work2, 0, DT_REAL);
   NewVector(work3, 0, DT_REAL);
   NewVector(work4, 0, DT_REAL);
   NewVector(work5, 0, DT_REAL);
   NewVector(iwork, 0, DT_INT);
   NewObj(hqrndstate, rs);
   hqrndseed(346436, 2434, &rs);
   maxits = 5;
   extra = 2;
   ae_vector_set_length(&work1, imax2(n, 1) + 1);
   ae_vector_set_length(&work2, imax2(n - 1, 1) + 1);
   ae_vector_set_length(&work3, imax2(n, 1) + 1);
   ae_vector_set_length(&work4, imax2(n, 1) + 1);
   ae_vector_set_length(&work5, imax2(n, 1) + 1);
   ae_vector_set_length(&iwork, imax2(n, 1) + 1);
   ae_vector_set_length(ifail, imax2(m, 1) + 1);
   ae_matrix_set_length(z, imax2(n, 1) + 1, imax2(m, 1) + 1);
// these initializers are not really necessary,
// but without them compiler complains about uninitialized locals
   gpind = 0;
   onenrm = 0.0;
   ortol = 0.0;
   dtpcrt = 0.0;
   xjm = 0.0;
// Test the input parameters.
   *info = 0;
   for (i = 1; i <= m; i++) {
      ifail->xZ[i] = 0;
   }
   if (n < 0) {
      *info = -1;
      ae_frame_leave();
      return;
   }
   if (m < 0 || m > n) {
      *info = -4;
      ae_frame_leave();
      return;
   }
   for (j = 2; j <= m; j++) {
      if (iblock->xZ[j] < iblock->xZ[j - 1]) {
         *info = -6;
         break;
      }
      if (iblock->xZ[j] == iblock->xZ[j - 1] && w->xR[j] < w->xR[j - 1]) {
         *info = -5;
         break;
      }
   }
   if (*info != 0) {
      ae_frame_leave();
      return;
   }
// Quick return if possible
   if (n == 0 || m == 0) {
      ae_frame_leave();
      return;
   }
   if (n == 1) {
      z->xyR[1][1] = 1.0;
      ae_frame_leave();
      return;
   }
// Some preparations
   ti = n - 1;
   ae_v_move(&work1.xR[1], 1, &e->xR[1], 1, ti);
   ae_vector_set_length(e, n + 1);
   ae_v_move(&e->xR[1], 1, &work1.xR[1], 1, ti);
   ae_v_move(&work1.xR[1], 1, &w->xR[1], 1, m);
   ae_vector_set_length(w, n + 1);
   ae_v_move(&w->xR[1], 1, &work1.xR[1], 1, m);
// Get machine constants.
   eps = ae_machineepsilon;
// Compute eigenvectors of matrix blocks.
   j1 = 1;
   for (nblk = 1; nblk <= iblock->xZ[m]; nblk++) {
   // Find starting and ending indices of block nblk.
      if (nblk == 1) {
         b1 = 1;
      } else {
         b1 = isplit->xZ[nblk - 1] + 1;
      }
      bn = isplit->xZ[nblk];
      blksiz = bn - b1 + 1;
      if (blksiz != 1) {
      // Compute reorthogonalization criterion and stopping criterion.
         gpind = b1;
         onenrm = fabs(d->xR[b1]) + fabs(e->xR[b1]);
         onenrm = rmax2(onenrm, fabs(d->xR[bn]) + fabs(e->xR[bn - 1]));
         for (i = b1 + 1; i < bn; i++) {
            onenrm = rmax2(onenrm, fabs(d->xR[i]) + fabs(e->xR[i - 1]) + fabs(e->xR[i]));
         }
         ortol = 0.001 * onenrm;
         dtpcrt = sqrt(0.1 / blksiz);
      }
   // Loop through eigenvalues of block nblk.
      jblk = 0;
      for (j = j1; j <= m; j++) {
         if (iblock->xZ[j] != nblk) {
            j1 = j;
            break;
         }
         jblk++;
         xj = w->xR[j];
         if (blksiz == 1) {
         // Skip all the work if the block size is one.
            work1.xR[1] = 1.0;
         } else {
         // If eigenvalues j and j-1 are too close, add a relatively
         // small perturbation.
            if (jblk > 1) {
               eps1 = fabs(eps * xj);
               pertol = 10 * eps1;
               sep = xj - xjm;
               if (sep < pertol) {
                  xj = xjm + pertol;
               }
            }
            its = 0;
            nrmchk = 0;
         // Get random starting vector.
            for (ti = 1; ti <= blksiz; ti++) {
               work1.xR[ti] = hqrndmiduniformr(&rs);
            }
         // Copy the matrix T so it won't be destroyed in factorization.
            for (ti = 1; ti < blksiz; ti++) {
               work2.xR[ti] = e->xR[b1 + ti - 1];
               work3.xR[ti] = e->xR[b1 + ti - 1];
               work4.xR[ti] = d->xR[b1 + ti - 1];
            }
            work4.xR[blksiz] = d->xR[b1 + blksiz - 1];
         // Compute LU factors with partial pivoting  ( PT = LU )
            tol = 0.0;
            evd_tdininternaldlagtf(blksiz, &work4, xj, &work2, &work3, tol, &work5, &iwork, &iinfo);
         // Update iteration count.
            do {
               its++;
               if (its > maxits) {
               // If stopping criterion was not satisfied, update info and
               // store eigenvector number in array ifail.
                  ++*info;
                  ifail->xZ[*info] = j;
                  break;
               }
            // Normalize and scale the righthand side vector Pb.
               v = 0.0;
               for (ti = 1; ti <= blksiz; ti++) {
                  v += fabs(work1.xR[ti]);
               }
               scl = blksiz * onenrm * rmax2(eps, fabs(work4.xR[blksiz])) / v;
               ae_v_muld(&work1.xR[1], 1, blksiz, scl);
            // Solve the system LU = Pb.
               evd_tdininternaldlagts(blksiz, &work4, &work2, &work3, &work5, &iwork, &work1, &tol, &iinfo);
            // Reorthogonalize by modified Gram-Schmidt if eigenvalues are
            // close enough.
               if (jblk != 1) {
                  if (!NearAtR(xj, xjm, ortol)) {
                     gpind = j;
                  }
                  if (gpind != j) {
                     for (i = gpind; i < j; i++) {
                        i1 = b1;
                        ztr = ae_v_dotproduct(&work1.xR[1], 1, &z->xyR[i1][i], z->stride, blksiz);
                        ae_v_subd(&work1.xR[1], 1, &z->xyR[i1][i], z->stride, blksiz, ztr);
                     }
                  }
               }
            // Check the infinity norm of the iterate.
               jmax = vectoridxabsmax(&work1, 1, blksiz);
               nrm = fabs(work1.xR[jmax]);
            // Continue for additional iterations after norm reaches
            // stopping criterion.
               tmpcriterion = false;
               if (nrm < dtpcrt) {
                  tmpcriterion = true;
               } else {
                  nrmchk++;
                  if (nrmchk < extra + 1) {
                     tmpcriterion = true;
                  }
               }
            } while (tmpcriterion);
         // Accept iterate as jth eigenvector.
            scl = 1 / vectornorm2(&work1, 1, blksiz);
            jmax = vectoridxabsmax(&work1, 1, blksiz);
            if (work1.xR[jmax] < 0.0) {
               scl = -scl;
            }
            ae_v_muld(&work1.xR[1], 1, blksiz, scl);
         }
         for (i = 1; i <= n; i++) {
            z->xyR[i][j] = 0.0;
         }
         for (i = 1; i <= blksiz; i++) {
            z->xyR[b1 + i - 1][j] = work1.xR[i];
         }
      // Save the shift to check eigenvalue spacing at next
      // iteration.
         xjm = xj;
      }
   }
   ae_frame_leave();
}

// performs complex division in  real arithmetic
//
//                         a + i*b
//              p + i*q = ---------
//                         c + i*d
//
// The algorithm is due to Robert L. Smith and can be found
// in D. Knuth, The art of Computer Programming, Vol.2, p.195
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
static void evd_internalhsevdladiv(double a, double b, double c, double d, double *p, double *q) {
   double e;
   double f;
   *p = 0;
   *q = 0;
   if (fabs(d) < fabs(c)) {
      e = d / c;
      f = c + d * e;
      *p = (a + b * e) / f;
      *q = (b - a * e) / f;
   } else {
      e = c / d;
      f = d + c * e;
      *p = (b + a * e) / f;
      *q = (-a + b * e) / f;
   }
}

// DLALN2 solves a system of the form  (ca A - w D ) X = s B
// or (ca A' - w D) X = s B   with possible scaling ("s") and
// perturbation of A.  (A' means A-transpose.)
//
// A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
// real diagonal matrix, w is a real or complex value, and X and B are
// NA x 1 matrices -- real if w is real, complex if w is complex.  NA
// may be 1 or 2.
//
// If w is complex, X and B are represented as NA x 2 matrices,
// the first column of each being the real part and the second
// being the imaginary part.
//
// "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
// so chosen that X can be computed without overflow.  X is further
// scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
// than overflow.
//
// If both singular values of (ca A - w D) are less than SMIN,
// SMIN*identity will be used instead of (ca A - w D).  If only one
// singular value is less than SMIN, one element of (ca A - w D) will be
// perturbed enough to make the smallest singular value roughly SMIN.
// If both singular values are at least SMIN, (ca A - w D) will not be
// perturbed.  In any case, the perturbation will be at most some small
// multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
// are computed by infinity-norm approximations, and thus will only be
// correct to a factor of 2 or so.
//
// Note: all input quantities are assumed to be smaller than overflow
// by a reasonable factor.  (See BIGNUM.)
//
//   -- LAPACK auxiliary routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      October 31, 1992
static void evd_internalhsevdlaln2(bool ltrans, ae_int_t na, ae_int_t nw, double smin, double ca, RMatrix *a, double d1, double d2, RMatrix *b, double wr, double wi, BVector *rswap4, BVector *zswap4, ZMatrix *ipivot44, RVector *civ4, RVector *crv4, RMatrix *x, double *scl, double *xnorm, ae_int_t *info) {
   ae_int_t icmax;
   ae_int_t j;
   double bbnd;
   double bi1;
   double bi2;
   double bignum;
   double bnorm;
   double br1;
   double br2;
   double ci21;
   double ci22;
   double cmax;
   double cnorm;
   double cr21;
   double cr22;
   double csi;
   double csr;
   double li21;
   double lr21;
   double smini;
   double smlnum;
   double temp;
   double u22abs;
   double ui11;
   double ui11r;
   double ui12;
   double ui12s;
   double ui22;
   double ur11;
   double ur11r;
   double ur12;
   double ur12s;
   double ur22;
   double xi1;
   double xi2;
   double xr1;
   double xr2;
   double tmp1;
   double tmp2;
   *scl = 0;
   *xnorm = 0;
   *info = 0;
   zswap4->xB[1] = false;
   zswap4->xB[2] = false;
   zswap4->xB[3] = true;
   zswap4->xB[4] = true;
   rswap4->xB[1] = false;
   rswap4->xB[2] = true;
   rswap4->xB[3] = false;
   rswap4->xB[4] = true;
   ipivot44->xyZ[1][1] = 1;
   ipivot44->xyZ[2][1] = 2;
   ipivot44->xyZ[3][1] = 3;
   ipivot44->xyZ[4][1] = 4;
   ipivot44->xyZ[1][2] = 2;
   ipivot44->xyZ[2][2] = 1;
   ipivot44->xyZ[3][2] = 4;
   ipivot44->xyZ[4][2] = 3;
   ipivot44->xyZ[1][3] = 3;
   ipivot44->xyZ[2][3] = 4;
   ipivot44->xyZ[3][3] = 1;
   ipivot44->xyZ[4][3] = 2;
   ipivot44->xyZ[1][4] = 4;
   ipivot44->xyZ[2][4] = 3;
   ipivot44->xyZ[3][4] = 2;
   ipivot44->xyZ[4][4] = 1;
   smlnum = 2 * ae_minrealnumber;
   bignum = 1 / smlnum;
   smini = rmax2(smin, smlnum);
// Don't check for input errors
   *info = 0;
// Standard Initializations
   *scl = 1.0;
   if (na == 1) {
   // 1 x 1  (i.e., scalar) system   C X = B
      if (nw == 1) {
      // Real 1x1 system.
      //
      // C = ca A - w D
         csr = ca * a->xyR[1][1] - wr * d1;
         cnorm = fabs(csr);
      // If | C | < SMINI, use C = SMINI
         if (cnorm < smini) {
            csr = smini;
            cnorm = smini;
            *info = 1;
         }
      // Check scaling for  X = B / C
         bnorm = fabs(b->xyR[1][1]);
         if (cnorm < 1.0 && bnorm > 1.0) {
            if (bnorm > bignum * cnorm) {
               *scl = 1 / bnorm;
            }
         }
      // Compute X
         x->xyR[1][1] = b->xyR[1][1] * (*scl) / csr;
         *xnorm = fabs(x->xyR[1][1]);
      } else {
      // Complex 1x1 system (w is complex)
      //
      // C = ca A - w D
         csr = ca * a->xyR[1][1] - wr * d1;
         csi = -wi * d1;
         cnorm = fabs(csr) + fabs(csi);
      // If | C | < SMINI, use C = SMINI
         if (cnorm < smini) {
            csr = smini;
            csi = 0.0;
            cnorm = smini;
            *info = 1;
         }
      // Check scaling for  X = B / C
         bnorm = fabs(b->xyR[1][1]) + fabs(b->xyR[1][2]);
         if (cnorm < 1.0 && bnorm > 1.0) {
            if (bnorm > bignum * cnorm) {
               *scl = 1 / bnorm;
            }
         }
      // Compute X
         evd_internalhsevdladiv(*scl * b->xyR[1][1], *scl * b->xyR[1][2], csr, csi, &tmp1, &tmp2);
         x->xyR[1][1] = tmp1;
         x->xyR[1][2] = tmp2;
         *xnorm = fabs(x->xyR[1][1]) + fabs(x->xyR[1][2]);
      }
   } else {
   // 2x2 System
   //
   // Compute the real part of  C = ca A - w D  (or  ca A' - w D )
      crv4->xR[1] = ca * a->xyR[1][1] - wr * d1;
      crv4->xR[2 + 2] = ca * a->xyR[2][2] - wr * d2;
      if (ltrans) {
         crv4->xR[1 + 2] = ca * a->xyR[2][1];
         crv4->xR[2] = ca * a->xyR[1][2];
      } else {
         crv4->xR[2] = ca * a->xyR[2][1];
         crv4->xR[1 + 2] = ca * a->xyR[1][2];
      }
      if (nw == 1) {
      // Real 2x2 system  (w is real)
      //
      // Find the largest element in C
         cmax = 0.0;
         icmax = 0;
         for (j = 1; j <= 4; j++) {
            if (!SmallAtR(crv4->xR[j], cmax)) {
               cmax = fabs(crv4->xR[j]);
               icmax = j;
            }
         }
      // If norm(C) < SMINI, use SMINI*identity.
         if (cmax < smini) {
            bnorm = rmax2(fabs(b->xyR[1][1]), fabs(b->xyR[2][1]));
            if (smini < 1.0 && bnorm > 1.0) {
               if (bnorm > bignum * smini) {
                  *scl = 1 / bnorm;
               }
            }
            temp = *scl / smini;
            x->xyR[1][1] = temp * b->xyR[1][1];
            x->xyR[2][1] = temp * b->xyR[2][1];
            *xnorm = temp * bnorm;
            *info = 1;
            return;
         }
      // Gaussian elimination with complete pivoting.
         ur11 = crv4->xR[icmax];
         cr21 = crv4->xR[ipivot44->xyZ[2][icmax]];
         ur12 = crv4->xR[ipivot44->xyZ[3][icmax]];
         cr22 = crv4->xR[ipivot44->xyZ[4][icmax]];
         ur11r = 1 / ur11;
         lr21 = ur11r * cr21;
         ur22 = cr22 - ur12 * lr21;
      // If smaller pivot < SMINI, use SMINI
         if (SmallR(ur22, smini)) {
            ur22 = smini;
            *info = 1;
         }
         if (rswap4->xB[icmax]) {
            br1 = b->xyR[2][1];
            br2 = b->xyR[1][1];
         } else {
            br1 = b->xyR[1][1];
            br2 = b->xyR[2][1];
         }
         br2 -= lr21 * br1;
         bbnd = rmax2(fabs(br1 * (ur22 * ur11r)), fabs(br2));
         if (bbnd > 1.0 && SmallR(ur22, 1.0)) {
            if (bbnd >= bignum * fabs(ur22)) {
               *scl = 1 / bbnd;
            }
         }
         xr2 = br2 * (*scl) / ur22;
         xr1 = *scl * br1 * ur11r - xr2 * (ur11r * ur12);
         if (zswap4->xB[icmax]) {
            x->xyR[1][1] = xr2;
            x->xyR[2][1] = xr1;
         } else {
            x->xyR[1][1] = xr1;
            x->xyR[2][1] = xr2;
         }
         *xnorm = rmax2(fabs(xr1), fabs(xr2));
      // Further scaling if  norm(A) norm(X) > overflow
         if (*xnorm > 1.0 && cmax > 1.0) {
            if (*xnorm > bignum / cmax) {
               temp = cmax / bignum;
               x->xyR[1][1] *= temp;
               x->xyR[2][1] *= temp;
               *xnorm *= temp;
               *scl *= temp;
            }
         }
      } else {
      // Complex 2x2 system  (w is complex)
      //
      // Find the largest element in C
         civ4->xR[1] = -wi * d1;
         civ4->xR[2] = 0.0;
         civ4->xR[1 + 2] = 0.0;
         civ4->xR[2 + 2] = -wi * d2;
         cmax = 0.0;
         icmax = 0;
         for (j = 1; j <= 4; j++) {
            if (fabs(crv4->xR[j]) + fabs(civ4->xR[j]) > cmax) {
               cmax = fabs(crv4->xR[j]) + fabs(civ4->xR[j]);
               icmax = j;
            }
         }
      // If norm(C) < SMINI, use SMINI*identity.
         if (cmax < smini) {
            bnorm = rmax2(fabs(b->xyR[1][1]) + fabs(b->xyR[1][2]), fabs(b->xyR[2][1]) + fabs(b->xyR[2][2]));
            if (smini < 1.0 && bnorm > 1.0) {
               if (bnorm > bignum * smini) {
                  *scl = 1 / bnorm;
               }
            }
            temp = *scl / smini;
            x->xyR[1][1] = temp * b->xyR[1][1];
            x->xyR[2][1] = temp * b->xyR[2][1];
            x->xyR[1][2] = temp * b->xyR[1][2];
            x->xyR[2][2] = temp * b->xyR[2][2];
            *xnorm = temp * bnorm;
            *info = 1;
            return;
         }
      // Gaussian elimination with complete pivoting.
         ur11 = crv4->xR[icmax];
         ui11 = civ4->xR[icmax];
         cr21 = crv4->xR[ipivot44->xyZ[2][icmax]];
         ci21 = civ4->xR[ipivot44->xyZ[2][icmax]];
         ur12 = crv4->xR[ipivot44->xyZ[3][icmax]];
         ui12 = civ4->xR[ipivot44->xyZ[3][icmax]];
         cr22 = crv4->xR[ipivot44->xyZ[4][icmax]];
         ci22 = civ4->xR[ipivot44->xyZ[4][icmax]];
         if (icmax == 1 || icmax == 4) {
         // Code when off-diagonals of pivoted C are real
            if (fabs(ur11) > fabs(ui11)) {
               temp = ui11 / ur11;
               ur11r = 1 / (ur11 * (1 + ae_sqr(temp)));
               ui11r = -temp * ur11r;
            } else {
               temp = ur11 / ui11;
               ui11r = -1 / (ui11 * (1 + ae_sqr(temp)));
               ur11r = -temp * ui11r;
            }
            lr21 = cr21 * ur11r;
            li21 = cr21 * ui11r;
            ur12s = ur12 * ur11r;
            ui12s = ur12 * ui11r;
            ur22 = cr22 - ur12 * lr21;
            ui22 = ci22 - ur12 * li21;
         } else {
         // Code when diagonals of pivoted C are real
            ur11r = 1 / ur11;
            ui11r = 0.0;
            lr21 = cr21 * ur11r;
            li21 = ci21 * ur11r;
            ur12s = ur12 * ur11r;
            ui12s = ui12 * ur11r;
            ur22 = cr22 - ur12 * lr21 + ui12 * li21;
            ui22 = -ur12 * li21 - ui12 * lr21;
         }
         u22abs = fabs(ur22) + fabs(ui22);
      // If smaller pivot < SMINI, use SMINI
         if (u22abs < smini) {
            ur22 = smini;
            ui22 = 0.0;
            *info = 1;
         }
         if (rswap4->xB[icmax]) {
            br2 = b->xyR[1][1];
            br1 = b->xyR[2][1];
            bi2 = b->xyR[1][2];
            bi1 = b->xyR[2][2];
         } else {
            br1 = b->xyR[1][1];
            br2 = b->xyR[2][1];
            bi1 = b->xyR[1][2];
            bi2 = b->xyR[2][2];
         }
         br2 -= lr21 * br1 - li21 * bi1;
         bi2 -= li21 * br1 + lr21 * bi1;
         bbnd = rmax2((fabs(br1) + fabs(bi1)) * (u22abs * (fabs(ur11r) + fabs(ui11r))), fabs(br2) + fabs(bi2));
         if (bbnd > 1.0 && u22abs < 1.0) {
            if (bbnd >= bignum * u22abs) {
               *scl = 1 / bbnd;
               br1 *= *scl;
               bi1 *= *scl;
               br2 *= *scl;
               bi2 *= *scl;
            }
         }
         evd_internalhsevdladiv(br2, bi2, ur22, ui22, &xr2, &xi2);
         xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
         xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
         if (zswap4->xB[icmax]) {
            x->xyR[1][1] = xr2;
            x->xyR[2][1] = xr1;
            x->xyR[1][2] = xi2;
            x->xyR[2][2] = xi1;
         } else {
            x->xyR[1][1] = xr1;
            x->xyR[2][1] = xr2;
            x->xyR[1][2] = xi1;
            x->xyR[2][2] = xi2;
         }
         *xnorm = rmax2(fabs(xr1) + fabs(xi1), fabs(xr2) + fabs(xi2));
      // Further scaling if  norm(A) norm(X) > overflow
         if (*xnorm > 1.0 && cmax > 1.0) {
            if (*xnorm > bignum / cmax) {
               temp = cmax / bignum;
               x->xyR[1][1] *= temp;
               x->xyR[2][1] *= temp;
               x->xyR[1][2] *= temp;
               x->xyR[2][2] *= temp;
               *xnorm *= temp;
               *scl *= temp;
            }
         }
      }
   }
}

// Internal subroutine
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      June 30, 1999
static void evd_internaltrevc(RMatrix *t, ae_int_t n, ae_int_t side, ae_int_t howmny, BVector *vselect, RMatrix *vl, RMatrix *vr, ae_int_t *m, ae_int_t *info) {
   ae_frame _frame_block;
   bool allv;
   bool bothv;
   bool leftv;
   bool over;
   bool pair;
   bool rightv;
   bool somev;
   ae_int_t i;
   ae_int_t ierr;
   ae_int_t ii;
   ae_int_t ip;
   ae_int_t iis;
   ae_int_t j;
   ae_int_t j1;
   ae_int_t j2;
   ae_int_t jnxt;
   ae_int_t k;
   ae_int_t ki;
   ae_int_t n2;
   double beta;
   double bignum;
   double emax;
   double rec;
   double remax;
   double scl;
   double smin;
   double smlnum;
   double ulp;
   double unfl;
   double vcrit;
   double vmax;
   double wi;
   double wr;
   double xnorm;
   bool skipflag;
   ae_int_t k1;
   ae_int_t k2;
   ae_int_t k3;
   ae_int_t k4;
   double vt;
   ae_frame_make(&_frame_block);
   DupVector(vselect);
   *m = 0;
   *info = 0;
   NewMatrix(x, 0, 0, DT_REAL);
   NewVector(work, 0, DT_REAL);
   NewVector(temp, 0, DT_REAL);
   NewMatrix(temp11, 0, 0, DT_REAL);
   NewMatrix(temp22, 0, 0, DT_REAL);
   NewMatrix(temp11b, 0, 0, DT_REAL);
   NewMatrix(temp21b, 0, 0, DT_REAL);
   NewMatrix(temp12b, 0, 0, DT_REAL);
   NewMatrix(temp22b, 0, 0, DT_REAL);
   NewVector(rswap4, 0, DT_BOOL);
   NewVector(zswap4, 0, DT_BOOL);
   NewMatrix(ipivot44, 0, 0, DT_INT);
   NewVector(civ4, 0, DT_REAL);
   NewVector(crv4, 0, DT_REAL);
   ae_matrix_set_length(&x, 2 + 1, 2 + 1);
   ae_matrix_set_length(&temp11, 1 + 1, 1 + 1);
   ae_matrix_set_length(&temp11b, 1 + 1, 1 + 1);
   ae_matrix_set_length(&temp21b, 2 + 1, 1 + 1);
   ae_matrix_set_length(&temp12b, 1 + 1, 2 + 1);
   ae_matrix_set_length(&temp22b, 2 + 1, 2 + 1);
   ae_matrix_set_length(&temp22, 2 + 1, 2 + 1);
   ae_vector_set_length(&work, 3 * n + 1);
   ae_vector_set_length(&temp, n + 1);
   ae_vector_set_length(&rswap4, 4 + 1);
   ae_vector_set_length(&zswap4, 4 + 1);
   ae_matrix_set_length(&ipivot44, 4 + 1, 4 + 1);
   ae_vector_set_length(&civ4, 4 + 1);
   ae_vector_set_length(&crv4, 4 + 1);
   if (howmny != 1) {
      if (side == 1 || side == 3) {
         ae_matrix_set_length(vr, n + 1, n + 1);
      }
      if (side == 2 || side == 3) {
         ae_matrix_set_length(vl, n + 1, n + 1);
      }
   }
// Decode and test the input parameters
   bothv = side == 3;
   rightv = side == 1 || bothv;
   leftv = side == 2 || bothv;
   allv = howmny == 2;
   over = howmny == 1;
   somev = howmny == 3;
   *info = 0;
   if (n < 0) {
      *info = -2;
      ae_frame_leave();
      return;
   }
   if (!rightv && !leftv) {
      *info = -3;
      ae_frame_leave();
      return;
   }
   if (!allv && !over && !somev) {
      *info = -4;
      ae_frame_leave();
      return;
   }
// Set M to the number of columns required to store the selected
// eigenvectors, standardize the array SELECT if necessary, and
// test MM.
   if (somev) {
      *m = 0;
      pair = false;
      for (j = 1; j <= n; j++) {
         if (pair) {
            pair = false;
            vselect->xB[j] = false;
         } else {
            if (j < n) {
               if (t->xyR[j + 1][j] == 0.0) {
                  if (vselect->xB[j]) {
                     ++*m;
                  }
               } else {
                  pair = true;
                  if (vselect->xB[j] || vselect->xB[j + 1]) {
                     vselect->xB[j] = true;
                     *m += 2;
                  }
               }
            } else {
               if (vselect->xB[n]) {
                  ++*m;
               }
            }
         }
      }
   } else {
      *m = n;
   }
// Quick return if possible.
   if (n == 0) {
      ae_frame_leave();
      return;
   }
// Set the constants to control overflow.
   unfl = ae_minrealnumber;
   ulp = ae_machineepsilon;
   smlnum = unfl * (n / ulp);
   bignum = (1 - ulp) / smlnum;
// Compute 1-norm of each column of strictly upper triangular
// part of T to control overflow in triangular solver.
   work.xR[1] = 0.0;
   for (j = 2; j <= n; j++) {
      work.xR[j] = 0.0;
      for (i = 1; i < j; i++) {
         work.xR[j] += fabs(t->xyR[i][j]);
      }
   }
// Index IP is used to specify the real or complex eigenvalue:
// IP = 0, real eigenvalue,
//      1, first of conjugate complex pair: (wr,wi)
//     -1, second of conjugate complex pair: (wr,wi)
   n2 = 2 * n;
   if (rightv) {
   // Compute right eigenvectors.
      ip = 0;
      iis = *m;
      for (ki = n; ki >= 1; ki--) {
         skipflag = false;
         if (ip == 1) {
            skipflag = true;
         } else {
            if (ki != 1) {
               if (t->xyR[ki][ki - 1] != 0.0) {
                  ip = -1;
               }
            }
            if (somev) {
               if (ip == 0) {
                  if (!vselect->xB[ki]) {
                     skipflag = true;
                  }
               } else {
                  if (!vselect->xB[ki - 1]) {
                     skipflag = true;
                  }
               }
            }
         }
         if (!skipflag) {
         // Compute the KI-th eigenvalue (WR,WI).
            wr = t->xyR[ki][ki];
            wi = 0.0;
            if (ip != 0) {
               wi = sqrt(fabs(t->xyR[ki][ki - 1])) * sqrt(fabs(t->xyR[ki - 1][ki]));
            }
            smin = rmax2(ulp * (fabs(wr) + fabs(wi)), smlnum);
            if (ip == 0) {
            // Real right eigenvector
               work.xR[ki + n] = 1.0;
            // Form right-hand side
               for (k = 1; k < ki; k++) {
                  work.xR[k + n] = -t->xyR[k][ki];
               }
            // Solve the upper quasi-triangular system:
            //   (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
               jnxt = ki - 1;
               for (j = ki - 1; j >= 1; j--) {
                  if (j > jnxt) {
                     continue;
                  }
                  j1 = j;
                  j2 = j;
                  jnxt = j - 1;
                  if (j > 1) {
                     if (t->xyR[j][j - 1] != 0.0) {
                        j1 = j - 1;
                        jnxt = j - 2;
                     }
                  }
                  if (j1 == j2) {
                  // 1-by-1 diagonal block
                     temp11.xyR[1][1] = t->xyR[j][j];
                     temp11b.xyR[1][1] = work.xR[j + n];
                     evd_internalhsevdlaln2(false, 1, 1, smin, 1.0, &temp11, 1.0, 1.0, &temp11b, wr, 0.0, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale X(1,1) to avoid overflow when updating
                  // the right-hand side.
                     if (xnorm > 1.0) {
                        if (work.xR[j] > bignum / xnorm) {
                           x.xyR[1][1] /= xnorm;
                           scl /= xnorm;
                        }
                     }
                  // Scale if necessary
                     if (scl != 1.0) {
                        k1 = n + 1;
                        k2 = n + ki;
                        ae_v_muld(&work.xR[k1], 1, k2 - k1 + 1, scl);
                     }
                     work.xR[j + n] = x.xyR[1][1];
                  // Update right-hand side
                     k1 = 1 + n;
                     k2 = j - 1 + n;
                     k3 = j - 1;
                     vt = -x.xyR[1][1];
                     ae_v_addd(&work.xR[k1], 1, &t->xyR[1][j], t->stride, k2 - k1 + 1, vt);
                  } else {
                  // 2-by-2 diagonal block
                     temp22.xyR[1][1] = t->xyR[j - 1][j - 1];
                     temp22.xyR[1][2] = t->xyR[j - 1][j];
                     temp22.xyR[2][1] = t->xyR[j][j - 1];
                     temp22.xyR[2][2] = t->xyR[j][j];
                     temp21b.xyR[1][1] = work.xR[j - 1 + n];
                     temp21b.xyR[2][1] = work.xR[j + n];
                     evd_internalhsevdlaln2(false, 2, 1, smin, 1.0, &temp22, 1.0, 1.0, &temp21b, wr, 0.0, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale X(1,1) and X(2,1) to avoid overflow when
                  // updating the right-hand side.
                     if (xnorm > 1.0) {
                        beta = rmax2(work.xR[j - 1], work.xR[j]);
                        if (beta > bignum / xnorm) {
                           x.xyR[1][1] /= xnorm;
                           x.xyR[2][1] /= xnorm;
                           scl /= xnorm;
                        }
                     }
                  // Scale if necessary
                     if (scl != 1.0) {
                        k1 = 1 + n;
                        k2 = ki + n;
                        ae_v_muld(&work.xR[k1], 1, k2 - k1 + 1, scl);
                     }
                     work.xR[j - 1 + n] = x.xyR[1][1];
                     work.xR[j + n] = x.xyR[2][1];
                  // Update right-hand side
                     k1 = 1 + n;
                     k2 = j - 2 + n;
                     k3 = j - 2;
                     k4 = j - 1;
                     vt = -x.xyR[1][1];
                     ae_v_addd(&work.xR[k1], 1, &t->xyR[1][k4], t->stride, k2 - k1 + 1, vt);
                     vt = -x.xyR[2][1];
                     ae_v_addd(&work.xR[k1], 1, &t->xyR[1][j], t->stride, k2 - k1 + 1, vt);
                  }
               }
            // Copy the vector x or Q*x to VR and normalize.
               if (!over) {
                  k1 = 1 + n;
                  k2 = ki + n;
                  ae_v_move(&vr->xyR[1][iis], vr->stride, &work.xR[k1], 1, ki);
                  ii = columnidxabsmax(vr, 1, ki, iis);
                  remax = 1 / fabs(vr->xyR[ii][iis]);
                  ae_v_muld(&vr->xyR[1][iis], vr->stride, ki, remax);
                  for (k = ki + 1; k <= n; k++) {
                     vr->xyR[k][iis] = 0.0;
                  }
               } else {
                  if (ki > 1) {
                     ae_v_move(&temp.xR[1], 1, &vr->xyR[1][ki], vr->stride, n);
                     matrixvectormultiply(vr, 1, n, 1, ki - 1, false, &work, 1 + n, ki - 1 + n, 1.0, &temp, 1, n, work.xR[ki + n]);
                     ae_v_move(&vr->xyR[1][ki], vr->stride, &temp.xR[1], 1, n);
                  }
                  ii = columnidxabsmax(vr, 1, n, ki);
                  remax = 1 / fabs(vr->xyR[ii][ki]);
                  ae_v_muld(&vr->xyR[1][ki], vr->stride, n, remax);
               }
            } else {
            // Complex right eigenvector.
            //
            // Initial solve
            //     [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
            //     [ (T(KI,KI-1)   T(KI,KI)   )               ]
               if (fabs(t->xyR[ki - 1][ki]) >= fabs(t->xyR[ki][ki - 1])) {
                  work.xR[ki - 1 + n] = 1.0;
                  work.xR[ki + n2] = wi / t->xyR[ki - 1][ki];
               } else {
                  work.xR[ki - 1 + n] = -wi / t->xyR[ki][ki - 1];
                  work.xR[ki + n2] = 1.0;
               }
               work.xR[ki + n] = 0.0;
               work.xR[ki - 1 + n2] = 0.0;
            // Form right-hand side
               for (k = 1; k < ki - 1; k++) {
                  work.xR[k + n] = -work.xR[ki - 1 + n] * t->xyR[k][ki - 1];
                  work.xR[k + n2] = -work.xR[ki + n2] * t->xyR[k][ki];
               }
            // Solve upper quasi-triangular system:
            // (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
               jnxt = ki - 2;
               for (j = ki - 2; j >= 1; j--) {
                  if (j > jnxt) {
                     continue;
                  }
                  j1 = j;
                  j2 = j;
                  jnxt = j - 1;
                  if (j > 1) {
                     if (t->xyR[j][j - 1] != 0.0) {
                        j1 = j - 1;
                        jnxt = j - 2;
                     }
                  }
                  if (j1 == j2) {
                  // 1-by-1 diagonal block
                     temp11.xyR[1][1] = t->xyR[j][j];
                     temp12b.xyR[1][1] = work.xR[j + n];
                     temp12b.xyR[1][2] = work.xR[j + n + n];
                     evd_internalhsevdlaln2(false, 1, 2, smin, 1.0, &temp11, 1.0, 1.0, &temp12b, wr, wi, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale X(1,1) and X(1,2) to avoid overflow when
                  // updating the right-hand side.
                     if (xnorm > 1.0) {
                        if (work.xR[j] > bignum / xnorm) {
                           x.xyR[1][1] /= xnorm;
                           x.xyR[1][2] /= xnorm;
                           scl /= xnorm;
                        }
                     }
                  // Scale if necessary
                     if (scl != 1.0) {
                        k1 = 1 + n;
                        k2 = ki + n;
                        ae_v_muld(&work.xR[k1], 1, k2 - k1 + 1, scl);
                        k1 = 1 + n2;
                        k2 = ki + n2;
                        ae_v_muld(&work.xR[k1], 1, k2 - k1 + 1, scl);
                     }
                     work.xR[j + n] = x.xyR[1][1];
                     work.xR[j + n2] = x.xyR[1][2];
                  // Update the right-hand side
                     k1 = 1 + n;
                     k2 = j - 1 + n;
                     k3 = 1;
                     k4 = j - 1;
                     vt = -x.xyR[1][1];
                     ae_v_addd(&work.xR[k1], 1, &t->xyR[k3][j], t->stride, k2 - k1 + 1, vt);
                     k1 = 1 + n2;
                     k2 = j - 1 + n2;
                     k3 = 1;
                     k4 = j - 1;
                     vt = -x.xyR[1][2];
                     ae_v_addd(&work.xR[k1], 1, &t->xyR[k3][j], t->stride, k2 - k1 + 1, vt);
                  } else {
                  // 2-by-2 diagonal block
                     temp22.xyR[1][1] = t->xyR[j - 1][j - 1];
                     temp22.xyR[1][2] = t->xyR[j - 1][j];
                     temp22.xyR[2][1] = t->xyR[j][j - 1];
                     temp22.xyR[2][2] = t->xyR[j][j];
                     temp22b.xyR[1][1] = work.xR[j - 1 + n];
                     temp22b.xyR[1][2] = work.xR[j - 1 + n + n];
                     temp22b.xyR[2][1] = work.xR[j + n];
                     temp22b.xyR[2][2] = work.xR[j + n + n];
                     evd_internalhsevdlaln2(false, 2, 2, smin, 1.0, &temp22, 1.0, 1.0, &temp22b, wr, wi, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale X to avoid overflow when updating
                  // the right-hand side.
                     if (xnorm > 1.0) {
                        beta = rmax2(work.xR[j - 1], work.xR[j]);
                        if (beta > bignum / xnorm) {
                           rec = 1 / xnorm;
                           x.xyR[1][1] *= rec;
                           x.xyR[1][2] *= rec;
                           x.xyR[2][1] *= rec;
                           x.xyR[2][2] *= rec;
                           scl *= rec;
                        }
                     }
                  // Scale if necessary
                     if (scl != 1.0) {
                        ae_v_muld(&work.xR[1 + n], 1, ki, scl);
                        ae_v_muld(&work.xR[1 + n2], 1, ki, scl);
                     }
                     work.xR[j - 1 + n] = x.xyR[1][1];
                     work.xR[j + n] = x.xyR[2][1];
                     work.xR[j - 1 + n2] = x.xyR[1][2];
                     work.xR[j + n2] = x.xyR[2][2];
                  // Update the right-hand side
                     vt = -x.xyR[1][1];
                     ae_v_addd(&work.xR[n + 1], 1, &t->xyR[1][j - 1], t->stride, j - 2, vt);
                     vt = -x.xyR[2][1];
                     ae_v_addd(&work.xR[n + 1], 1, &t->xyR[1][j], t->stride, j - 2, vt);
                     vt = -x.xyR[1][2];
                     ae_v_addd(&work.xR[n2 + 1], 1, &t->xyR[1][j - 1], t->stride, j - 2, vt);
                     vt = -x.xyR[2][2];
                     ae_v_addd(&work.xR[n2 + 1], 1, &t->xyR[1][j], t->stride, j - 2, vt);
                  }
               }
            // Copy the vector x or Q*x to VR and normalize.
               if (!over) {
                  ae_v_move(&vr->xyR[1][iis - 1], vr->stride, &work.xR[n + 1], 1, ki);
                  ae_v_move(&vr->xyR[1][iis], vr->stride, &work.xR[n2 + 1], 1, ki);
                  emax = 0.0;
                  for (k = 1; k <= ki; k++) {
                     emax = rmax2(emax, fabs(vr->xyR[k][iis - 1]) + fabs(vr->xyR[k][iis]));
                  }
                  remax = 1 / emax;
                  ae_v_muld(&vr->xyR[1][iis - 1], vr->stride, ki, remax);
                  ae_v_muld(&vr->xyR[1][iis], vr->stride, ki, remax);
                  for (k = ki + 1; k <= n; k++) {
                     vr->xyR[k][iis - 1] = 0.0;
                     vr->xyR[k][iis] = 0.0;
                  }
               } else {
                  if (ki > 2) {
                     ae_v_move(&temp.xR[1], 1, &vr->xyR[1][ki - 1], vr->stride, n);
                     matrixvectormultiply(vr, 1, n, 1, ki - 2, false, &work, 1 + n, ki - 2 + n, 1.0, &temp, 1, n, work.xR[ki - 1 + n]);
                     ae_v_move(&vr->xyR[1][ki - 1], vr->stride, &temp.xR[1], 1, n);
                     ae_v_move(&temp.xR[1], 1, &vr->xyR[1][ki], vr->stride, n);
                     matrixvectormultiply(vr, 1, n, 1, ki - 2, false, &work, 1 + n2, ki - 2 + n2, 1.0, &temp, 1, n, work.xR[ki + n2]);
                     ae_v_move(&vr->xyR[1][ki], vr->stride, &temp.xR[1], 1, n);
                  } else {
                     vt = work.xR[ki - 1 + n];
                     ae_v_muld(&vr->xyR[1][ki - 1], vr->stride, n, vt);
                     vt = work.xR[ki + n2];
                     ae_v_muld(&vr->xyR[1][ki], vr->stride, n, vt);
                  }
                  emax = 0.0;
                  for (k = 1; k <= n; k++) {
                     emax = rmax2(emax, fabs(vr->xyR[k][ki - 1]) + fabs(vr->xyR[k][ki]));
                  }
                  remax = 1 / emax;
                  ae_v_muld(&vr->xyR[1][ki - 1], vr->stride, n, remax);
                  ae_v_muld(&vr->xyR[1][ki], vr->stride, n, remax);
               }
            }
            iis--;
            if (ip != 0) {
               iis--;
            }
         }
         if (ip == 1) {
            ip = 0;
         }
         if (ip == -1) {
            ip = 1;
         }
      }
   }
   if (leftv) {
   // Compute left eigenvectors.
      ip = 0;
      iis = 1;
      for (ki = 1; ki <= n; ki++) {
         skipflag = false;
         if (ip == -1) {
            skipflag = true;
         } else {
            if (ki != n) {
               if (t->xyR[ki + 1][ki] != 0.0) {
                  ip = 1;
               }
            }
            if (somev) {
               if (!vselect->xB[ki]) {
                  skipflag = true;
               }
            }
         }
         if (!skipflag) {
         // Compute the KI-th eigenvalue (WR,WI).
            wr = t->xyR[ki][ki];
            wi = 0.0;
            if (ip != 0) {
               wi = sqrt(fabs(t->xyR[ki][ki + 1])) * sqrt(fabs(t->xyR[ki + 1][ki]));
            }
            smin = rmax2(ulp * (fabs(wr) + fabs(wi)), smlnum);
            if (ip == 0) {
            // Real left eigenvector.
               work.xR[ki + n] = 1.0;
            // Form right-hand side
               for (k = ki + 1; k <= n; k++) {
                  work.xR[k + n] = -t->xyR[ki][k];
               }
            // Solve the quasi-triangular system:
            // (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
               vmax = 1.0;
               vcrit = bignum;
               jnxt = ki + 1;
               for (j = ki + 1; j <= n; j++) {
                  if (j < jnxt) {
                     continue;
                  }
                  j1 = j;
                  j2 = j;
                  jnxt = j + 1;
                  if (j < n) {
                     if (t->xyR[j + 1][j] != 0.0) {
                        j2 = j + 1;
                        jnxt = j + 2;
                     }
                  }
                  if (j1 == j2) {
                  // 1-by-1 diagonal block
                  //
                  // Scale if necessary to avoid overflow when forming
                  // the right-hand side.
                     if (work.xR[j] > vcrit) {
                        rec = 1 / vmax;
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, rec);
                        vmax = 1.0;
                        vcrit = bignum;
                     }
                     vt = ae_v_dotproduct(&t->xyR[ki + 1][j], t->stride, &work.xR[ki + 1 + n], 1, j - ki - 1);
                     work.xR[j + n] -= vt;
                  // Solve (T(J,J)-WR)'*X = WORK
                     temp11.xyR[1][1] = t->xyR[j][j];
                     temp11b.xyR[1][1] = work.xR[j + n];
                     evd_internalhsevdlaln2(false, 1, 1, smin, 1.0, &temp11, 1.0, 1.0, &temp11b, wr, 0.0, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale if necessary
                     if (scl != 1.0) {
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, scl);
                     }
                     work.xR[j + n] = x.xyR[1][1];
                     vmax = rmax2(fabs(work.xR[j + n]), vmax);
                     vcrit = bignum / vmax;
                  } else {
                  // 2-by-2 diagonal block
                  //
                  // Scale if necessary to avoid overflow when forming
                  // the right-hand side.
                     beta = rmax2(work.xR[j], work.xR[j + 1]);
                     if (beta > vcrit) {
                        rec = 1 / vmax;
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, rec);
                        vmax = 1.0;
                        vcrit = bignum;
                     }
                     vt = ae_v_dotproduct(&t->xyR[ki + 1][j], t->stride, &work.xR[ki + 1 + n], 1, j - ki - 1);
                     work.xR[j + n] -= vt;
                     vt = ae_v_dotproduct(&t->xyR[ki + 1][j + 1], t->stride, &work.xR[ki + 1 + n], 1, j - ki - 1);
                     work.xR[j + 1 + n] -= vt;
                  // Solve
                  //    [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
                  //    [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
                     temp22.xyR[1][1] = t->xyR[j][j];
                     temp22.xyR[1][2] = t->xyR[j][j + 1];
                     temp22.xyR[2][1] = t->xyR[j + 1][j];
                     temp22.xyR[2][2] = t->xyR[j + 1][j + 1];
                     temp21b.xyR[1][1] = work.xR[j + n];
                     temp21b.xyR[2][1] = work.xR[j + 1 + n];
                     evd_internalhsevdlaln2(true, 2, 1, smin, 1.0, &temp22, 1.0, 1.0, &temp21b, wr, 0.0, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale if necessary
                     if (scl != 1.0) {
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, scl);
                     }
                     work.xR[j + n] = x.xyR[1][1];
                     work.xR[j + 1 + n] = x.xyR[2][1];
                     vmax = rmax2(fabs(work.xR[j + n]), rmax2(fabs(work.xR[j + 1 + n]), vmax));
                     vcrit = bignum / vmax;
                  }
               }
            // Copy the vector x or Q*x to VL and normalize.
               if (!over) {
                  ae_v_move(&vl->xyR[ki][iis], vl->stride, &work.xR[ki + n], 1, n - ki + 1);
                  ii = columnidxabsmax(vl, ki, n, iis);
                  remax = 1 / fabs(vl->xyR[ii][iis]);
                  ae_v_muld(&vl->xyR[ki][iis], vl->stride, n - ki + 1, remax);
                  for (k = 1; k < ki; k++) {
                     vl->xyR[k][iis] = 0.0;
                  }
               } else {
                  if (ki < n) {
                     ae_v_move(&temp.xR[1], 1, &vl->xyR[1][ki], vl->stride, n);
                     matrixvectormultiply(vl, 1, n, ki + 1, n, false, &work, ki + 1 + n, n + n, 1.0, &temp, 1, n, work.xR[ki + n]);
                     ae_v_move(&vl->xyR[1][ki], vl->stride, &temp.xR[1], 1, n);
                  }
                  ii = columnidxabsmax(vl, 1, n, ki);
                  remax = 1 / fabs(vl->xyR[ii][ki]);
                  ae_v_muld(&vl->xyR[1][ki], vl->stride, n, remax);
               }
            } else {
            // Complex left eigenvector.
            //
            // Initial solve:
            //   ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
            //   ((T(KI+1,KI) T(KI+1,KI+1))                )
               if (fabs(t->xyR[ki][ki + 1]) >= fabs(t->xyR[ki + 1][ki])) {
                  work.xR[ki + n] = wi / t->xyR[ki][ki + 1];
                  work.xR[ki + 1 + n2] = 1.0;
               } else {
                  work.xR[ki + n] = 1.0;
                  work.xR[ki + 1 + n2] = -wi / t->xyR[ki + 1][ki];
               }
               work.xR[ki + 1 + n] = 0.0;
               work.xR[ki + n2] = 0.0;
            // Form right-hand side
               for (k = ki + 2; k <= n; k++) {
                  work.xR[k + n] = -work.xR[ki + n] * t->xyR[ki][k];
                  work.xR[k + n2] = -work.xR[ki + 1 + n2] * t->xyR[ki + 1][k];
               }
            // Solve complex quasi-triangular system:
            // ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
               vmax = 1.0;
               vcrit = bignum;
               jnxt = ki + 2;
               for (j = ki + 2; j <= n; j++) {
                  if (j < jnxt) {
                     continue;
                  }
                  j1 = j;
                  j2 = j;
                  jnxt = j + 1;
                  if (j < n) {
                     if (t->xyR[j + 1][j] != 0.0) {
                        j2 = j + 1;
                        jnxt = j + 2;
                     }
                  }
                  if (j1 == j2) {
                  // 1-by-1 diagonal block
                  //
                  // Scale if necessary to avoid overflow when
                  // forming the right-hand side elements.
                     if (work.xR[j] > vcrit) {
                        rec = 1 / vmax;
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, rec);
                        ae_v_muld(&work.xR[ki + n2], 1, n - ki + 1, rec);
                        vmax = 1.0;
                        vcrit = bignum;
                     }
                     vt = ae_v_dotproduct(&t->xyR[ki + 2][j], t->stride, &work.xR[ki + 2 + n], 1, j - ki - 2);
                     work.xR[j + n] -= vt;
                     vt = ae_v_dotproduct(&t->xyR[ki + 2][j], t->stride, &work.xR[ki + 2 + n2], 1, j - ki - 2);
                     work.xR[j + n2] -= vt;
                  // Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
                     temp11.xyR[1][1] = t->xyR[j][j];
                     temp12b.xyR[1][1] = work.xR[j + n];
                     temp12b.xyR[1][2] = work.xR[j + n + n];
                     evd_internalhsevdlaln2(false, 1, 2, smin, 1.0, &temp11, 1.0, 1.0, &temp12b, wr, -wi, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale if necessary
                     if (scl != 1.0) {
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, scl);
                        ae_v_muld(&work.xR[ki + n2], 1, n - ki + 1, scl);
                     }
                     work.xR[j + n] = x.xyR[1][1];
                     work.xR[j + n2] = x.xyR[1][2];
                     vmax = rmax2(fabs(work.xR[j + n]), rmax2(fabs(work.xR[j + n2]), vmax));
                     vcrit = bignum / vmax;
                  } else {
                  // 2-by-2 diagonal block
                  //
                  // Scale if necessary to avoid overflow when forming
                  // the right-hand side elements.
                     beta = rmax2(work.xR[j], work.xR[j + 1]);
                     if (beta > vcrit) {
                        rec = 1 / vmax;
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, rec);
                        ae_v_muld(&work.xR[ki + n2], 1, n - ki + 1, rec);
                        vmax = 1.0;
                        vcrit = bignum;
                     }
                     vt = ae_v_dotproduct(&t->xyR[ki + 2][j], t->stride, &work.xR[ki + 2 + n], 1, j - ki - 2);
                     work.xR[j + n] -= vt;
                     vt = ae_v_dotproduct(&t->xyR[ki + 2][j], t->stride, &work.xR[ki + 2 + n2], 1, j - ki - 2);
                     work.xR[j + n2] -= vt;
                     vt = ae_v_dotproduct(&t->xyR[ki + 2][j + 1], t->stride, &work.xR[ki + 2 + n], 1, j - ki - 2);
                     work.xR[j + 1 + n] -= vt;
                     vt = ae_v_dotproduct(&t->xyR[ki + 2][j + 1], t->stride, &work.xR[ki + 2 + n2], 1, j - ki - 2);
                     work.xR[j + 1 + n2] -= vt;
                  // Solve 2-by-2 complex linear equation
                  //   ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
                  //   ([T(j+1,j) T(j+1,j+1)]             )
                     temp22.xyR[1][1] = t->xyR[j][j];
                     temp22.xyR[1][2] = t->xyR[j][j + 1];
                     temp22.xyR[2][1] = t->xyR[j + 1][j];
                     temp22.xyR[2][2] = t->xyR[j + 1][j + 1];
                     temp22b.xyR[1][1] = work.xR[j + n];
                     temp22b.xyR[1][2] = work.xR[j + n + n];
                     temp22b.xyR[2][1] = work.xR[j + 1 + n];
                     temp22b.xyR[2][2] = work.xR[j + 1 + n + n];
                     evd_internalhsevdlaln2(true, 2, 2, smin, 1.0, &temp22, 1.0, 1.0, &temp22b, wr, -wi, &rswap4, &zswap4, &ipivot44, &civ4, &crv4, &x, &scl, &xnorm, &ierr);
                  // Scale if necessary
                     if (scl != 1.0) {
                        ae_v_muld(&work.xR[ki + n], 1, n - ki + 1, scl);
                        ae_v_muld(&work.xR[ki + n2], 1, n - ki + 1, scl);
                     }
                     work.xR[j + n] = x.xyR[1][1];
                     work.xR[j + n2] = x.xyR[1][2];
                     work.xR[j + 1 + n] = x.xyR[2][1];
                     work.xR[j + 1 + n2] = x.xyR[2][2];
                     vmax = rmax2(fabs(x.xyR[1][1]), vmax);
                     vmax = rmax2(fabs(x.xyR[1][2]), vmax);
                     vmax = rmax2(fabs(x.xyR[2][1]), vmax);
                     vmax = rmax2(fabs(x.xyR[2][2]), vmax);
                     vcrit = bignum / vmax;
                  }
               }
            // Copy the vector x or Q*x to VL and normalize.
               if (!over) {
                  ae_v_move(&vl->xyR[ki][iis], vl->stride, &work.xR[ki + n], 1, n - ki + 1);
                  ae_v_move(&vl->xyR[ki][iis + 1], vl->stride, &work.xR[ki + n2], 1, n - ki + 1);
                  emax = 0.0;
                  for (k = ki; k <= n; k++) {
                     emax = rmax2(emax, fabs(vl->xyR[k][iis]) + fabs(vl->xyR[k][iis + 1]));
                  }
                  remax = 1 / emax;
                  ae_v_muld(&vl->xyR[ki][iis], vl->stride, n - ki + 1, remax);
                  ae_v_muld(&vl->xyR[ki][iis + 1], vl->stride, n - ki + 1, remax);
                  for (k = 1; k < ki; k++) {
                     vl->xyR[k][iis] = 0.0;
                     vl->xyR[k][iis + 1] = 0.0;
                  }
               } else {
                  if (ki < n - 1) {
                     ae_v_move(&temp.xR[1], 1, &vl->xyR[1][ki], vl->stride, n);
                     matrixvectormultiply(vl, 1, n, ki + 2, n, false, &work, ki + 2 + n, n + n, 1.0, &temp, 1, n, work.xR[ki + n]);
                     ae_v_move(&vl->xyR[1][ki], vl->stride, &temp.xR[1], 1, n);
                     ae_v_move(&temp.xR[1], 1, &vl->xyR[1][ki + 1], vl->stride, n);
                     matrixvectormultiply(vl, 1, n, ki + 2, n, false, &work, ki + 2 + n2, n + n2, 1.0, &temp, 1, n, work.xR[ki + 1 + n2]);
                     ae_v_move(&vl->xyR[1][ki + 1], vl->stride, &temp.xR[1], 1, n);
                  } else {
                     vt = work.xR[ki + n];
                     ae_v_muld(&vl->xyR[1][ki], vl->stride, n, vt);
                     vt = work.xR[ki + 1 + n2];
                     ae_v_muld(&vl->xyR[1][ki + 1], vl->stride, n, vt);
                  }
                  emax = 0.0;
                  for (k = 1; k <= n; k++) {
                     emax = rmax2(emax, fabs(vl->xyR[k][ki]) + fabs(vl->xyR[k][ki + 1]));
                  }
                  remax = 1 / emax;
                  ae_v_muld(&vl->xyR[1][ki], vl->stride, n, remax);
                  ae_v_muld(&vl->xyR[1][ki + 1], vl->stride, n, remax);
               }
            }
            iis++;
            if (ip != 0) {
               iis++;
            }
         }
         if (ip == -1) {
            ip = 0;
         }
         if (ip == 1) {
            ip = -1;
         }
      }
   }
   ae_frame_leave();
}

// Internal subroutine
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      June 30, 1999
static void evd_rmatrixinternaltrevc(RMatrix *t, ae_int_t n, ae_int_t side, ae_int_t howmny, BVector *vselect, RMatrix *vl, RMatrix *vr, ae_int_t *m, ae_int_t *info) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   ae_frame_make(&_frame_block);
   DupVector(vselect);
   *m = 0;
   *info = 0;
   NewMatrix(t1, 0, 0, DT_REAL);
   NewMatrix(vl1, 0, 0, DT_REAL);
   NewMatrix(vr1, 0, 0, DT_REAL);
   NewVector(vselect1, 0, DT_BOOL);
// Allocate VL/VR, if needed
   if (howmny == 2 || howmny == 3) {
      if (side == 1 || side == 3) {
         matrixsetlengthatleast(vr, n, n);
      }
      if (side == 2 || side == 3) {
         matrixsetlengthatleast(vl, n, n);
      }
   }
// Try to use MKL kernel
   if (rmatrixinternaltrevcmkl(t, n, side, howmny, vl, vr, m, info)) {
      ae_frame_leave();
      return;
   }
// ALGLIB version
   ae_matrix_set_length(&t1, n + 1, n + 1);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         t1.xyR[i + 1][j + 1] = t->xyR[i][j];
      }
   }
   if (howmny == 3) {
      ae_vector_set_length(&vselect1, n + 1);
      for (i = 0; i < n; i++) {
         vselect1.xB[1 + i] = vselect->xB[i];
      }
   }
   if ((side == 2 || side == 3) && howmny == 1) {
      ae_matrix_set_length(&vl1, n + 1, n + 1);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            vl1.xyR[i + 1][j + 1] = vl->xyR[i][j];
         }
      }
   }
   if ((side == 1 || side == 3) && howmny == 1) {
      ae_matrix_set_length(&vr1, n + 1, n + 1);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            vr1.xyR[i + 1][j + 1] = vr->xyR[i][j];
         }
      }
   }
   evd_internaltrevc(&t1, n, side, howmny, &vselect1, &vl1, &vr1, m, info);
   if (side != 1) {
      matrixsetlengthatleast(vl, n, n);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            vl->xyR[i][j] = vl1.xyR[i + 1][j + 1];
         }
      }
   }
   if (side != 2) {
      matrixsetlengthatleast(vr, n, n);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            vr->xyR[i][j] = vr1.xyR[i + 1][j + 1];
         }
      }
   }
   ae_frame_leave();
}

// Finding the eigenvalues and eigenvectors of a tridiagonal symmetric matrix
//
// The algorithm finds the eigen pairs of a tridiagonal symmetric matrix by
// using an QL/QR algorithm with implicit shifts.
//
// Inputs:
//     D       -   the main diagonal of a tridiagonal matrix.
//                 Array whose index ranges within [0..N-1].
//     E       -   the secondary diagonal of a tridiagonal matrix.
//                 Array whose index ranges within [0..N-2].
//     N       -   size of matrix A.
//     ZNeeded -   flag controlling whether the eigenvectors are needed or not.
//                 If ZNeeded is equal to:
//                  * 0, the eigenvectors are not needed;
//                  * 1, the eigenvectors of a tridiagonal matrix
//                    are multiplied by the square matrix Z. It is used if the
//                    tridiagonal matrix is obtained by the similarity
//                    transformation of a symmetric matrix;
//                  * 2, the eigenvectors of a tridiagonal matrix replace the
//                    square matrix Z;
//                  * 3, matrix Z contains the first row of the eigenvectors
//                    matrix.
//     Z       -   if ZNeeded=1, Z contains the square matrix by which the
//                 eigenvectors are multiplied.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//
// Outputs:
//     D       -   eigenvalues in ascending order.
//                 Array whose index ranges within [0..N-1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z hasn't changed;
//                  * 1, Z contains the product of a given matrix (from the left)
//                    and the eigenvectors matrix (from the right);
//                  * 2, Z contains the eigenvectors.
//                  * 3, Z contains the first row of the eigenvectors matrix.
//                 If ZNeeded < 3, Z is the array whose indexes range within [0..N-1, 0..N-1].
//                 In that case, the eigenvectors are stored in the matrix columns.
//                 If ZNeeded=3, Z is the array whose indexes range within [0..0, 0..N-1].
//
// Result:
//     True, if the algorithm has converged.
//     False, if the algorithm hasn't converged.
//
//   -- LAPACK routine (version 3.0) --
//      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//      Courant Institute, Argonne National Lab, and Rice University
//      September 30, 1994
// API: bool smatrixtdevd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, real_2d_array &z);
bool smatrixtdevd(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, RMatrix *z) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   bool result;
   ae_frame_make(&_frame_block);
   DupVector(e);
   NewVector(d1, 0, DT_REAL);
   NewVector(e1, 0, DT_REAL);
   NewVector(ex, 0, DT_REAL);
   NewMatrix(z1, 0, 0, DT_REAL);
   ae_assert(n >= 1, "SMatrixTDEVD: N <= 0");
   ae_assert(zneeded >= 0 && zneeded <= 3, "SMatrixTDEVD: incorrect ZNeeded");
   result = false;
// Preprocess Z: make ZNeeded equal to 0, 1 or 3.
// Ensure that memory for Z is allocated.
   if (zneeded == 2) {
   // Load identity to Z
      matrixsetlengthatleast(z, n, n);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            z->xyR[i][j] = 0.0;
         }
         z->xyR[i][i] = 1.0;
      }
      zneeded = 1;
   }
   if (zneeded == 3) {
   // Allocate memory
      matrixsetlengthatleast(z, 1, n);
   }
// Try to solve problem with MKL
   ae_vector_set_length(&ex, n);
   for (i = 0; i < n - 1; i++) {
      ex.xR[i] = e->xR[i];
   }
   if (smatrixtdevdmkl(d, &ex, n, zneeded, z, &result)) {
      ae_frame_leave();
      return result;
   }
// Prepare 1-based task
   ae_vector_set_length(&d1, n + 1);
   ae_vector_set_length(&e1, n + 1);
   ae_v_move(&d1.xR[1], 1, d->xR, 1, n);
   if (n > 1) {
      ae_v_move(&e1.xR[1], 1, e->xR, 1, n - 1);
   }
   if (zneeded == 1) {
      ae_matrix_set_length(&z1, n + 1, n + 1);
      for (i = 1; i <= n; i++) {
         ae_v_move(&z1.xyR[i][1], 1, z->xyR[i - 1], 1, n);
      }
   }
// Solve 1-based task
   result = evd_tridiagonalevd(&d1, &e1, n, zneeded, &z1);
   if (!result) {
      ae_frame_leave();
      return result;
   }
// Convert back to 0-based result
   ae_v_move(d->xR, 1, &d1.xR[1], 1, n);
   if (zneeded != 0) {
      if (zneeded == 1) {
         for (i = 1; i <= n; i++) {
            ae_v_move(z->xyR[i - 1], 1, &z1.xyR[i][1], 1, n);
         }
         ae_frame_leave();
         return result;
      }
      if (zneeded == 2) {
         ae_matrix_set_length(z, n, n);
         for (i = 1; i <= n; i++) {
            ae_v_move(z->xyR[i - 1], 1, &z1.xyR[i][1], 1, n);
         }
         ae_frame_leave();
         return result;
      }
      if (zneeded == 3) {
         ae_matrix_set_length(z, 0 + 1, n);
         ae_v_move(z->xyR[0], 1, &z1.xyR[1][1], 1, n);
         ae_frame_leave();
         return result;
      }
      ae_assert(false, "SMatrixTDEVD: Incorrect ZNeeded!");
   }
   ae_frame_leave();
   return result;
}

// Subroutine for finding the tridiagonal matrix eigenvalues/vectors in a
// given half-interval (A, B] by using bisection and inverse iteration.
//
// Inputs:
//     D       -   the main diagonal of a tridiagonal matrix.
//                 Array whose index ranges within [0..N-1].
//     E       -   the secondary diagonal of a tridiagonal matrix.
//                 Array whose index ranges within [0..N-2].
//     N       -   size of matrix, N >= 0.
//     ZNeeded -   flag controlling whether the eigenvectors are needed or not.
//                 If ZNeeded is equal to:
//                  * 0, the eigenvectors are not needed;
//                  * 1, the eigenvectors of a tridiagonal matrix are multiplied
//                    by the square matrix Z. It is used if the tridiagonal
//                    matrix is obtained by the similarity transformation
//                    of a symmetric matrix.
//                  * 2, the eigenvectors of a tridiagonal matrix replace matrix Z.
//     A, B    -   half-interval (A, B] to search eigenvalues in.
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z isn't used and remains unchanged;
//                  * 1, Z contains the square matrix (array whose indexes range
//                    within [0..N-1, 0..N-1]) which reduces the given symmetric
//                    matrix to tridiagonal form;
//                  * 2, Z isn't used (but changed on the exit).
//
// Outputs:
//     D       -   array of the eigenvalues found.
//                 Array whose index ranges within [0..M-1].
//     M       -   number of eigenvalues found in the given half-interval (M >= 0).
//     Z       -   if ZNeeded is equal to:
//                  * 0, doesn't contain any information;
//                  * 1, contains the product of a given NxN matrix Z (from the
//                    left) and NxM matrix of the eigenvectors found (from the
//                    right). Array whose indexes range within [0..N-1, 0..M-1].
//                  * 2, contains the matrix of the eigenvectors found.
//                    Array whose indexes range within [0..N-1, 0..M-1].
//
// Result:
//
//     True, if successful. In that case, M contains the number of eigenvalues
//     in the given half-interval (could be equal to 0), D contains the eigenvalues,
//     Z contains the eigenvectors (if needed).
//     It should be noted that the subroutine changes the size of arrays D and Z.
//
//     False, if the bisection method subroutine wasn't able to find the
//     eigenvalues in the given interval or if the inverse iteration subroutine
//     wasn't able to find all the corresponding eigenvectors. In that case,
//     the eigenvalues and eigenvectors are not returned, M is equal to 0.
// ALGLIB: Copyright 31.03.2008 by Sergey Bochkanov
// API: bool smatrixtdevdr(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const double a, const double b, ae_int_t &m, real_2d_array &z);
bool smatrixtdevdr(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, double a, double b, ae_int_t *m, RMatrix *z) {
   ae_frame _frame_block;
   ae_int_t errorcode;
   ae_int_t nsplit;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t cr;
   double v;
   bool result;
   ae_frame_make(&_frame_block);
   *m = 0;
   NewVector(iblock, 0, DT_INT);
   NewVector(isplit, 0, DT_INT);
   NewVector(ifail, 0, DT_INT);
   NewVector(d1, 0, DT_REAL);
   NewVector(e1, 0, DT_REAL);
   NewVector(w, 0, DT_REAL);
   NewMatrix(z2, 0, 0, DT_REAL);
   NewMatrix(z3, 0, 0, DT_REAL);
   ae_assert(zneeded >= 0 && zneeded <= 2, "SMatrixTDEVDR: incorrect ZNeeded!");
// Special cases
   if (b <= a) {
      *m = 0;
      result = true;
      ae_frame_leave();
      return result;
   }
   if (n <= 0) {
      *m = 0;
      result = true;
      ae_frame_leave();
      return result;
   }
// Copy D,E to D1, E1
   ae_vector_set_length(&d1, n + 1);
   ae_v_move(&d1.xR[1], 1, d->xR, 1, n);
   if (n > 1) {
      ae_vector_set_length(&e1, n);
      ae_v_move(&e1.xR[1], 1, e->xR, 1, n - 1);
   }
// No eigen vectors
   if (zneeded == 0) {
      result = evd_internalbisectioneigenvalues(&d1, &e1, n, 2, 1, a, b, 0, 0, -1.0, &w, m, &nsplit, &iblock, &isplit, &errorcode);
      if (!result || *m == 0) {
         *m = 0;
         ae_frame_leave();
         return result;
      }
      ae_vector_set_length(d, *m);
      ae_v_move(d->xR, 1, &w.xR[1], 1, *m);
      ae_frame_leave();
      return result;
   }
// Eigen vectors are multiplied by Z
   if (zneeded == 1) {
   // Find eigen pairs
      result = evd_internalbisectioneigenvalues(&d1, &e1, n, 2, 2, a, b, 0, 0, -1.0, &w, m, &nsplit, &iblock, &isplit, &errorcode);
      if (!result || *m == 0) {
         *m = 0;
         ae_frame_leave();
         return result;
      }
      evd_internaldstein(n, &d1, &e1, *m, &w, &iblock, &isplit, &z2, &ifail, &cr);
      if (cr != 0) {
         *m = 0;
         result = false;
         ae_frame_leave();
         return result;
      }
   // Sort eigen values and vectors
      for (i = 1; i <= *m; i++) {
         k = i;
         for (j = i; j <= *m; j++) {
            if (w.xR[j] < w.xR[k]) {
               k = j;
            }
         }
         swapr(&w.xR[i], &w.xR[k]);
         for (j = 1; j <= n; j++) {
            swapr(&z2.xyR[j][i], &z2.xyR[j][k]);
         }
      }
   // Transform Z2 and overwrite Z
      ae_matrix_set_length(&z3, *m + 1, n + 1);
      for (i = 1; i <= *m; i++) {
         ae_v_move(&z3.xyR[i][1], 1, &z2.xyR[1][i], z2.stride, n);
      }
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= *m; j++) {
            v = ae_v_dotproduct(z->xyR[i - 1], 1, &z3.xyR[j][1], 1, n);
            z2.xyR[i][j] = v;
         }
      }
      ae_matrix_set_length(z, n, *m);
      for (i = 1; i <= *m; i++) {
         ae_v_move(&z->xyR[0][i - 1], z->stride, &z2.xyR[1][i], z2.stride, n);
      }
   // Store W
      ae_vector_set_length(d, *m);
      for (i = 1; i <= *m; i++) {
         d->xR[i - 1] = w.xR[i];
      }
      ae_frame_leave();
      return result;
   }
// Eigen vectors are stored in Z
   if (zneeded == 2) {
   // Find eigen pairs
      result = evd_internalbisectioneigenvalues(&d1, &e1, n, 2, 2, a, b, 0, 0, -1.0, &w, m, &nsplit, &iblock, &isplit, &errorcode);
      if (!result || *m == 0) {
         *m = 0;
         ae_frame_leave();
         return result;
      }
      evd_internaldstein(n, &d1, &e1, *m, &w, &iblock, &isplit, &z2, &ifail, &cr);
      if (cr != 0) {
         *m = 0;
         result = false;
         ae_frame_leave();
         return result;
      }
   // Sort eigen values and vectors
      for (i = 1; i <= *m; i++) {
         k = i;
         for (j = i; j <= *m; j++) {
            if (w.xR[j] < w.xR[k]) {
               k = j;
            }
         }
         swapr(&w.xR[i], &w.xR[k]);
         for (j = 1; j <= n; j++) {
            swapr(&z2.xyR[j][i], &z2.xyR[j][k]);
         }
      }
   // Store W
      ae_vector_set_length(d, *m);
      for (i = 1; i <= *m; i++) {
         d->xR[i - 1] = w.xR[i];
      }
      ae_matrix_set_length(z, n, *m);
      for (i = 1; i <= *m; i++) {
         ae_v_move(&z->xyR[0][i - 1], z->stride, &z2.xyR[1][i], z2.stride, n);
      }
      ae_frame_leave();
      return result;
   }
   result = false;
   ae_frame_leave();
   return result;
}

// Subroutine for finding tridiagonal matrix eigenvalues/vectors with given
// indexes (in ascending order) by using the bisection and inverse iteraion.
//
// Inputs:
//     D       -   the main diagonal of a tridiagonal matrix.
//                 Array whose index ranges within [0..N-1].
//     E       -   the secondary diagonal of a tridiagonal matrix.
//                 Array whose index ranges within [0..N-2].
//     N       -   size of matrix. N >= 0.
//     ZNeeded -   flag controlling whether the eigenvectors are needed or not.
//                 If ZNeeded is equal to:
//                  * 0, the eigenvectors are not needed;
//                  * 1, the eigenvectors of a tridiagonal matrix are multiplied
//                    by the square matrix Z. It is used if the
//                    tridiagonal matrix is obtained by the similarity transformation
//                    of a symmetric matrix.
//                  * 2, the eigenvectors of a tridiagonal matrix replace
//                    matrix Z.
//     I1, I2  -   index interval for searching (from I1 to I2).
//                 0 <= I1 <= I2 <= N-1.
//     Z       -   if ZNeeded is equal to:
//                  * 0, Z isn't used and remains unchanged;
//                  * 1, Z contains the square matrix (array whose indexes range within [0..N-1, 0..N-1])
//                    which reduces the given symmetric matrix to  tridiagonal form;
//                  * 2, Z isn't used (but changed on the exit).
//
// Outputs:
//     D       -   array of the eigenvalues found.
//                 Array whose index ranges within [0..I2-I1].
//     Z       -   if ZNeeded is equal to:
//                  * 0, doesn't contain any information;
//                  * 1, contains the product of a given NxN matrix Z (from the left) and
//                    Nx(I2-I1) matrix of the eigenvectors found (from the right).
//                    Array whose indexes range within [0..N-1, 0..I2-I1].
//                  * 2, contains the matrix of the eigenvalues found.
//                    Array whose indexes range within [0..N-1, 0..I2-I1].
//
//
// Result:
//
//     True, if successful. In that case, D contains the eigenvalues,
//     Z contains the eigenvectors (if needed).
//     It should be noted that the subroutine changes the size of arrays D and Z.
//
//     False, if the bisection method subroutine wasn't able to find the eigenvalues
//     in the given interval or if the inverse iteration subroutine wasn't able
//     to find all the corresponding eigenvectors. In that case, the eigenvalues
//     and eigenvectors are not returned.
// ALGLIB: Copyright 25.12.2005 by Sergey Bochkanov
// API: bool smatrixtdevdi(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const ae_int_t i1, const ae_int_t i2, real_2d_array &z);
bool smatrixtdevdi(RVector *d, RVector *e, ae_int_t n, ae_int_t zneeded, ae_int_t i1, ae_int_t i2, RMatrix *z) {
   ae_frame _frame_block;
   ae_int_t errorcode;
   ae_int_t nsplit;
   ae_int_t i;
   ae_int_t j;
   ae_int_t k;
   ae_int_t m;
   ae_int_t cr;
   double v;
   bool result;
   ae_frame_make(&_frame_block);
   NewVector(iblock, 0, DT_INT);
   NewVector(isplit, 0, DT_INT);
   NewVector(ifail, 0, DT_INT);
   NewVector(w, 0, DT_REAL);
   NewVector(d1, 0, DT_REAL);
   NewVector(e1, 0, DT_REAL);
   NewMatrix(z2, 0, 0, DT_REAL);
   NewMatrix(z3, 0, 0, DT_REAL);
   ae_assert(0 <= i1 && i1 <= i2 && i2 < n, "SMatrixTDEVDI: incorrect I1/I2!");
// Copy D,E to D1, E1
   ae_vector_set_length(&d1, n + 1);
   ae_v_move(&d1.xR[1], 1, d->xR, 1, n);
   if (n > 1) {
      ae_vector_set_length(&e1, n);
      ae_v_move(&e1.xR[1], 1, e->xR, 1, n - 1);
   }
// No eigen vectors
   if (zneeded == 0) {
      result = evd_internalbisectioneigenvalues(&d1, &e1, n, 3, 1, 0.0, 0.0, i1 + 1, i2 + 1, -1.0, &w, &m, &nsplit, &iblock, &isplit, &errorcode);
      if (!result) {
         ae_frame_leave();
         return result;
      }
      if (m != i2 - i1 + 1) {
         result = false;
         ae_frame_leave();
         return result;
      }
      ae_vector_set_length(d, m);
      for (i = 1; i <= m; i++) {
         d->xR[i - 1] = w.xR[i];
      }
      ae_frame_leave();
      return result;
   }
// Eigen vectors are multiplied by Z
   if (zneeded == 1) {
   // Find eigen pairs
      result = evd_internalbisectioneigenvalues(&d1, &e1, n, 3, 2, 0.0, 0.0, i1 + 1, i2 + 1, -1.0, &w, &m, &nsplit, &iblock, &isplit, &errorcode);
      if (!result) {
         ae_frame_leave();
         return result;
      }
      if (m != i2 - i1 + 1) {
         result = false;
         ae_frame_leave();
         return result;
      }
      evd_internaldstein(n, &d1, &e1, m, &w, &iblock, &isplit, &z2, &ifail, &cr);
      if (cr != 0) {
         result = false;
         ae_frame_leave();
         return result;
      }
   // Sort eigen values and vectors
      for (i = 1; i <= m; i++) {
         k = i;
         for (j = i; j <= m; j++) {
            if (w.xR[j] < w.xR[k]) {
               k = j;
            }
         }
         swapr(&w.xR[i], &w.xR[k]);
         for (j = 1; j <= n; j++) {
            swapr(&z2.xyR[j][i], &z2.xyR[j][k]);
         }
      }
   // Transform Z2 and overwrite Z
      ae_matrix_set_length(&z3, m + 1, n + 1);
      for (i = 1; i <= m; i++) {
         ae_v_move(&z3.xyR[i][1], 1, &z2.xyR[1][i], z2.stride, n);
      }
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= m; j++) {
            v = ae_v_dotproduct(z->xyR[i - 1], 1, &z3.xyR[j][1], 1, n);
            z2.xyR[i][j] = v;
         }
      }
      ae_matrix_set_length(z, n, m);
      for (i = 1; i <= m; i++) {
         ae_v_move(&z->xyR[0][i - 1], z->stride, &z2.xyR[1][i], z2.stride, n);
      }
   // Store W
      ae_vector_set_length(d, m);
      for (i = 1; i <= m; i++) {
         d->xR[i - 1] = w.xR[i];
      }
      ae_frame_leave();
      return result;
   }
// Eigen vectors are stored in Z
   if (zneeded == 2) {
   // Find eigen pairs
      result = evd_internalbisectioneigenvalues(&d1, &e1, n, 3, 2, 0.0, 0.0, i1 + 1, i2 + 1, -1.0, &w, &m, &nsplit, &iblock, &isplit, &errorcode);
      if (!result) {
         ae_frame_leave();
         return result;
      }
      if (m != i2 - i1 + 1) {
         result = false;
         ae_frame_leave();
         return result;
      }
      evd_internaldstein(n, &d1, &e1, m, &w, &iblock, &isplit, &z2, &ifail, &cr);
      if (cr != 0) {
         result = false;
         ae_frame_leave();
         return result;
      }
   // Sort eigen values and vectors
      for (i = 1; i <= m; i++) {
         k = i;
         for (j = i; j <= m; j++) {
            if (w.xR[j] < w.xR[k]) {
               k = j;
            }
         }
         swapr(&w.xR[i], &w.xR[k]);
         for (j = 1; j <= n; j++) {
            swapr(&z2.xyR[j][i], &z2.xyR[j][k]);
         }
      }
   // Store Z
      ae_matrix_set_length(z, n, m);
      for (i = 1; i <= m; i++) {
         ae_v_move(&z->xyR[0][i - 1], z->stride, &z2.xyR[1][i], z2.stride, n);
      }
   // Store W
      ae_vector_set_length(d, m);
      for (i = 1; i <= m; i++) {
         d->xR[i - 1] = w.xR[i];
      }
      ae_frame_leave();
      return result;
   }
   result = false;
   ae_frame_leave();
   return result;
}

// Finding eigenvalues and eigenvectors of a general (unsymmetric) matrix
//
// The algorithm finds eigenvalues and eigenvectors of a general matrix by
// using the QR algorithm with multiple shifts. The algorithm can find
// eigenvalues and both left and right eigenvectors.
//
// The right eigenvector is a vector x such that A*x = w*x, and the left
// eigenvector is a vector y such that y'*A = w*y' (here y' implies a complex
// conjugate transposition of vector y).
//
// Inputs:
//     A       -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     VNeeded -   flag controlling whether eigenvectors are needed or not.
//                 If VNeeded is equal to:
//                  * 0, eigenvectors are not returned;
//                  * 1, right eigenvectors are returned;
//                  * 2, left eigenvectors are returned;
//                  * 3, both left and right eigenvectors are returned.
//
// Outputs:
//     WR      -   real parts of eigenvalues.
//                 Array whose index ranges within [0..N-1].
//     WR      -   imaginary parts of eigenvalues.
//                 Array whose index ranges within [0..N-1].
//     VL, VR  -   arrays of left and right eigenvectors (if they are needed).
//                 If WI[i]=0, the respective eigenvalue is a real number,
//                 and it corresponds to the column number I of matrices VL/VR.
//                 If WI[i] > 0, we have a pair of complex conjugate numbers with
//                 positive and negative imaginary parts:
//                     the first eigenvalue WR[i] + sqrt(-1)*WI[i];
//                     the second eigenvalue WR[i+1] + sqrt(-1)*WI[i+1];
//                     WI[i] > 0
//                     WI[i+1] = -WI[i] < 0
//                 In that case, the eigenvector  corresponding to the first
//                 eigenvalue is located in i and i+1 columns of matrices
//                 VL/VR (the column number i contains the real part, and the
//                 column number i+1 contains the imaginary part), and the vector
//                 corresponding to the second eigenvalue is a complex conjugate to
//                 the first vector.
//                 Arrays whose indexes range within [0..N-1, 0..N-1].
//
// Result:
//     True, if the algorithm has converged.
//     False, if the algorithm has not converged.
//
// Note 1:
//     Some users may ask the following question: what if WI[N-1] > 0?
//     WI[N] must contain an eigenvalue which is complex conjugate to the
//     N-th eigenvalue, but the array has only size N?
//     The answer is as follows: such a situation cannot occur because the
//     algorithm finds a pairs of eigenvalues, therefore, if WI[i] > 0, I is
//     strictly less than N-1.
//
// Note 2:
//     The algorithm performance depends on the value of the internal parameter
//     NS of the InternalSchurDecomposition subroutine which defines the number
//     of shifts in the QR algorithm (similarly to the block width in block-matrix
//     algorithms of linear algebra). If you require maximum performance
//     on your machine, it is recommended to adjust this parameter manually.
//
//
// See also the InternalTREVC subroutine.
//
// The algorithm is based on the LAPACK 3.0 library.
// API: bool rmatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t vneeded, real_1d_array &wr, real_1d_array &wi, real_2d_array &vl, real_2d_array &vr);
bool rmatrixevd(RMatrix *a, ae_int_t n, ae_int_t vneeded, RVector *wr, RVector *wi, RMatrix *vl, RMatrix *vr) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t info;
   ae_int_t m1;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(wr);
   SetVector(wi);
   SetMatrix(vl);
   SetMatrix(vr);
   NewMatrix(a1, 0, 0, DT_REAL);
   NewMatrix(vl1, 0, 0, DT_REAL);
   NewMatrix(vr1, 0, 0, DT_REAL);
   NewMatrix(s1, 0, 0, DT_REAL);
   NewMatrix(s, 0, 0, DT_REAL);
   NewMatrix(dummy, 0, 0, DT_REAL);
   NewVector(wr1, 0, DT_REAL);
   NewVector(wi1, 0, DT_REAL);
   NewVector(tau, 0, DT_REAL);
   NewVector(sel1, 0, DT_BOOL);
   ae_assert(vneeded >= 0 && vneeded <= 3, "RMatrixEVD: incorrect VNeeded!");
   if (vneeded == 0) {
   // Eigen values only
      rmatrixhessenberg(a, n, &tau);
      rmatrixinternalschurdecomposition(a, n, 0, 0, wr, wi, &dummy, &info);
      result = info == 0;
      ae_frame_leave();
      return result;
   }
// Eigen values and vectors
   rmatrixhessenberg(a, n, &tau);
   rmatrixhessenbergunpackq(a, n, &tau, &s);
   rmatrixinternalschurdecomposition(a, n, 1, 1, wr, wi, &s, &info);
   result = info == 0;
   if (!result) {
      ae_frame_leave();
      return result;
   }
   if (vneeded == 1 || vneeded == 3) {
      ae_matrix_set_length(vr, n, n);
      for (i = 0; i < n; i++) {
         ae_v_move(vr->xyR[i], 1, s.xyR[i], 1, n);
      }
   }
   if (vneeded == 2 || vneeded == 3) {
      ae_matrix_set_length(vl, n, n);
      for (i = 0; i < n; i++) {
         ae_v_move(vl->xyR[i], 1, s.xyR[i], 1, n);
      }
   }
   evd_rmatrixinternaltrevc(a, n, vneeded, 1, &sel1, vl, vr, &m1, &info);
   result = info == 0;
   ae_frame_leave();
   return result;
}

void eigsubspacestate_init(void *_p, bool make_automatic) {
   eigsubspacestate *p = (eigsubspacestate *)_p;
   hqrndstate_init(&p->rs, make_automatic);
   ae_vector_init(&p->tau, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->q0, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->qcur, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->qnew, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->znew, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->r, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->rz, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->tz, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->rq, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->dummy, 0, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->rw, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->tw, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->wcur, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->wprev, 0, DT_REAL, make_automatic);
   ae_vector_init(&p->wrank, 0, DT_REAL, make_automatic);
   apbuffers_init(&p->buf, make_automatic);
   ae_matrix_init(&p->x, 0, 0, DT_REAL, make_automatic);
   ae_matrix_init(&p->ax, 0, 0, DT_REAL, make_automatic);
}

void eigsubspacestate_copy(void *_dst, void *_src, bool make_automatic) {
   eigsubspacestate *dst = (eigsubspacestate *)_dst;
   eigsubspacestate *src = (eigsubspacestate *)_src;
   dst->n = src->n;
   dst->k = src->k;
   dst->nwork = src->nwork;
   dst->maxits = src->maxits;
   dst->eps = src->eps;
   dst->eigenvectorsneeded = src->eigenvectorsneeded;
   dst->matrixtype = src->matrixtype;
   dst->usewarmstart = src->usewarmstart;
   dst->firstcall = src->firstcall;
   hqrndstate_copy(&dst->rs, &src->rs, make_automatic);
   dst->running = src->running;
   ae_vector_copy(&dst->tau, &src->tau, make_automatic);
   ae_matrix_copy(&dst->q0, &src->q0, make_automatic);
   ae_matrix_copy(&dst->qcur, &src->qcur, make_automatic);
   ae_matrix_copy(&dst->qnew, &src->qnew, make_automatic);
   ae_matrix_copy(&dst->znew, &src->znew, make_automatic);
   ae_matrix_copy(&dst->r, &src->r, make_automatic);
   ae_matrix_copy(&dst->rz, &src->rz, make_automatic);
   ae_matrix_copy(&dst->tz, &src->tz, make_automatic);
   ae_matrix_copy(&dst->rq, &src->rq, make_automatic);
   ae_matrix_copy(&dst->dummy, &src->dummy, make_automatic);
   ae_vector_copy(&dst->rw, &src->rw, make_automatic);
   ae_vector_copy(&dst->tw, &src->tw, make_automatic);
   ae_vector_copy(&dst->wcur, &src->wcur, make_automatic);
   ae_vector_copy(&dst->wprev, &src->wprev, make_automatic);
   ae_vector_copy(&dst->wrank, &src->wrank, make_automatic);
   apbuffers_copy(&dst->buf, &src->buf, make_automatic);
   ae_matrix_copy(&dst->x, &src->x, make_automatic);
   ae_matrix_copy(&dst->ax, &src->ax, make_automatic);
   dst->requesttype = src->requesttype;
   dst->requestsize = src->requestsize;
   dst->repiterationscount = src->repiterationscount;
   dst->PQ = src->PQ;
}

void eigsubspacestate_free(void *_p, bool make_automatic) {
   eigsubspacestate *p = (eigsubspacestate *)_p;
   hqrndstate_free(&p->rs, make_automatic);
   ae_vector_free(&p->tau, make_automatic);
   ae_matrix_free(&p->q0, make_automatic);
   ae_matrix_free(&p->qcur, make_automatic);
   ae_matrix_free(&p->qnew, make_automatic);
   ae_matrix_free(&p->znew, make_automatic);
   ae_matrix_free(&p->r, make_automatic);
   ae_matrix_free(&p->rz, make_automatic);
   ae_matrix_free(&p->tz, make_automatic);
   ae_matrix_free(&p->rq, make_automatic);
   ae_matrix_free(&p->dummy, make_automatic);
   ae_vector_free(&p->rw, make_automatic);
   ae_vector_free(&p->tw, make_automatic);
   ae_vector_free(&p->wcur, make_automatic);
   ae_vector_free(&p->wprev, make_automatic);
   ae_vector_free(&p->wrank, make_automatic);
   apbuffers_free(&p->buf, make_automatic);
   ae_matrix_free(&p->x, make_automatic);
   ae_matrix_free(&p->ax, make_automatic);
}

void eigsubspacereport_init(void *_p, bool make_automatic) {
}

void eigsubspacereport_copy(void *_dst, void *_src, bool make_automatic) {
   eigsubspacereport *dst = (eigsubspacereport *)_dst;
   eigsubspacereport *src = (eigsubspacereport *)_src;
   dst->iterationscount = src->iterationscount;
}

void eigsubspacereport_free(void *_p, bool make_automatic) {
}
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the subspace iteration algorithm.
// You should use ALGLIB functions to work with this object.
DefClass(eigsubspacestate, EndD)

// This object stores state of the subspace iteration algorithm.
// You should use ALGLIB functions to work with this object.
DefClass(eigsubspacereport, AndD DecVal(iterationscount))

void eigsubspacecreate(const ae_int_t n, const ae_int_t k, eigsubspacestate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspacecreate(n, k, ConstT(eigsubspacestate, state));
   alglib_impl::ae_state_clear();
}

void eigsubspacecreatebuf(const ae_int_t n, const ae_int_t k, const eigsubspacestate &state) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspacecreatebuf(n, k, ConstT(eigsubspacestate, state));
   alglib_impl::ae_state_clear();
}

void eigsubspacesetcond(const eigsubspacestate &state, const double eps, const ae_int_t maxits) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspacesetcond(ConstT(eigsubspacestate, state), eps, maxits);
   alglib_impl::ae_state_clear();
}

void eigsubspacesetwarmstart(const eigsubspacestate &state, const bool usewarmstart) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspacesetwarmstart(ConstT(eigsubspacestate, state), usewarmstart);
   alglib_impl::ae_state_clear();
}

void eigsubspaceoocstart(const eigsubspacestate &state, const ae_int_t mtype) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspaceoocstart(ConstT(eigsubspacestate, state), mtype);
   alglib_impl::ae_state_clear();
}

bool eigsubspaceooccontinue(const eigsubspacestate &state) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::eigsubspaceooccontinue(ConstT(eigsubspacestate, state));
   alglib_impl::ae_state_clear();
   return Ok;
}

void eigsubspaceoocgetrequestinfo(const eigsubspacestate &state, ae_int_t &requesttype, ae_int_t &requestsize) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspaceoocgetrequestinfo(ConstT(eigsubspacestate, state), &requesttype, &requestsize);
   alglib_impl::ae_state_clear();
}

void eigsubspaceoocgetrequestdata(const eigsubspacestate &state, real_2d_array &x) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspaceoocgetrequestdata(ConstT(eigsubspacestate, state), ConstT(ae_matrix, x));
   alglib_impl::ae_state_clear();
}

void eigsubspaceoocsendresult(const eigsubspacestate &state, const real_2d_array &ax) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspaceoocsendresult(ConstT(eigsubspacestate, state), ConstT(ae_matrix, ax));
   alglib_impl::ae_state_clear();
}

void eigsubspaceoocstop(const eigsubspacestate &state, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspaceoocstop(ConstT(eigsubspacestate, state), ConstT(ae_vector, w), ConstT(ae_matrix, z), ConstT(eigsubspacereport, rep));
   alglib_impl::ae_state_clear();
}

void eigsubspacesolvedenses(const eigsubspacestate &state, const real_2d_array &a, const bool isupper, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspacesolvedenses(ConstT(eigsubspacestate, state), ConstT(ae_matrix, a), isupper, ConstT(ae_vector, w), ConstT(ae_matrix, z), ConstT(eigsubspacereport, rep));
   alglib_impl::ae_state_clear();
}

void eigsubspacesolvesparses(const eigsubspacestate &state, const sparsematrix &a, const bool isupper, real_1d_array &w, real_2d_array &z, eigsubspacereport &rep) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::eigsubspacesolvesparses(ConstT(eigsubspacestate, state), ConstT(sparsematrix, a), isupper, ConstT(ae_vector, w), ConstT(ae_matrix, z), ConstT(eigsubspacereport, rep));
   alglib_impl::ae_state_clear();
}

bool smatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixevd(ConstT(ae_matrix, a), n, zneeded, isupper, ConstT(ae_vector, d), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool smatrixevdr(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixevdr(ConstT(ae_matrix, a), n, zneeded, isupper, b1, b2, &m, ConstT(ae_vector, w), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool smatrixevdi(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixevdi(ConstT(ae_matrix, a), n, zneeded, isupper, i1, i2, ConstT(ae_vector, w), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool hmatrixevd(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, complex_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::hmatrixevd(ConstT(ae_matrix, a), n, zneeded, isupper, ConstT(ae_vector, d), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool hmatrixevdr(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, complex_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::hmatrixevdr(ConstT(ae_matrix, a), n, zneeded, isupper, b1, b2, &m, ConstT(ae_vector, w), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool hmatrixevdi(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, complex_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::hmatrixevdi(ConstT(ae_matrix, a), n, zneeded, isupper, i1, i2, ConstT(ae_vector, w), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool smatrixtdevd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixtdevd(ConstT(ae_vector, d), ConstT(ae_vector, e), n, zneeded, ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool smatrixtdevdr(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const double a, const double b, ae_int_t &m, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixtdevdr(ConstT(ae_vector, d), ConstT(ae_vector, e), n, zneeded, a, b, &m, ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool smatrixtdevdi(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const ae_int_t i1, const ae_int_t i2, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixtdevdi(ConstT(ae_vector, d), ConstT(ae_vector, e), n, zneeded, i1, i2, ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool rmatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t vneeded, real_1d_array &wr, real_1d_array &wi, real_2d_array &vl, real_2d_array &vr) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::rmatrixevd(ConstT(ae_matrix, a), n, vneeded, ConstT(ae_vector, wr), ConstT(ae_vector, wi), ConstT(ae_matrix, vl), ConstT(ae_matrix, vr));
   alglib_impl::ae_state_clear();
   return Ok;
}
} // end of namespace alglib

// === SCHUR Package ===
// Depends on: ORTFAC, HSSCHUR
namespace alglib_impl {
// Subroutine performing the Schur decomposition of a general matrix by using
// the QR algorithm with multiple shifts.
//
// The source matrix A is represented as S'*A*S = T, where S is an orthogonal
// matrix (Schur vectors), T - upper quasi-triangular matrix (with blocks of
// sizes 1x1 and 2x2 on the main diagonal).
//
// Inputs:
//     A   -   matrix to be decomposed.
//             Array whose indexes range within [0..N-1, 0..N-1].
//     N   -   size of A, N >= 0.
//
//
// Outputs:
//     A   -   contains matrix T.
//             Array whose indexes range within [0..N-1, 0..N-1].
//     S   -   contains Schur vectors.
//             Array whose indexes range within [0..N-1, 0..N-1].
//
// Note 1:
//     The block structure of matrix T can be easily recognized: since all
//     the elements below the blocks are zeros, the elements a[i+1,i] which
//     are equal to 0 show the block border.
//
// Note 2:
//     The algorithm performance depends on the value of the internal parameter
//     NS of the InternalSchurDecomposition subroutine which defines the number
//     of shifts in the QR algorithm (similarly to the block width in block-matrix
//     algorithms in linear algebra). If you require maximum performance on
//     your machine, it is recommended to adjust this parameter manually.
//
// Result:
//     True,
//         if the algorithm has converged and parameters A and S contain the result.
//     False,
//         if the algorithm has not converged.
//
// Algorithm implemented on the basis of the DHSEQR subroutine (LAPACK 3.0 library).
// API: bool rmatrixschur(real_2d_array &a, const ae_int_t n, real_2d_array &s);
bool rmatrixschur(RMatrix *a, ae_int_t n, RMatrix *s) {
   ae_frame _frame_block;
   ae_int_t info;
   bool result;
   ae_frame_make(&_frame_block);
   SetMatrix(s);
   NewVector(tau, 0, DT_REAL);
   NewVector(wi, 0, DT_REAL);
   NewVector(wr, 0, DT_REAL);
// Upper Hessenberg form of the 0-based matrix
   rmatrixhessenberg(a, n, &tau);
   rmatrixhessenbergunpackq(a, n, &tau, s);
// Schur decomposition
   rmatrixinternalschurdecomposition(a, n, 1, 1, &wr, &wi, s, &info);
   result = info == 0;
   ae_frame_leave();
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
bool rmatrixschur(real_2d_array &a, const ae_int_t n, real_2d_array &s) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::rmatrixschur(ConstT(ae_matrix, a), n, ConstT(ae_matrix, s));
   alglib_impl::ae_state_clear();
   return Ok;
}
} // end of namespace alglib

// === SPDGEVD Package ===
// Depends on: MATINV, EVD
namespace alglib_impl {
// Algorithm for solving the following generalized symmetric positive-definite
// eigenproblem:
//     A*x = lambda*B*x (1) or
//     A*B*x = lambda*x (2) or
//     B*A*x = lambda*x (3).
// where A is a symmetric matrix, B - symmetric positive-definite matrix.
// The problem is solved by reducing it to an ordinary  symmetric  eigenvalue
// problem.
//
// Inputs:
//     A           -   symmetric matrix which is given by its upper or lower
//                     triangular part.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     N           -   size of matrices A and B.
//     IsUpperA    -   storage format of matrix A.
//     B           -   symmetric positive-definite matrix which is given by
//                     its upper or lower triangular part.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     IsUpperB    -   storage format of matrix B.
//     ZNeeded     -   if ZNeeded is equal to:
//                      * 0, the eigenvectors are not returned;
//                      * 1, the eigenvectors are returned.
//     ProblemType -   if ProblemType is equal to:
//                      * 1, the following problem is solved: A*x = lambda*B*x;
//                      * 2, the following problem is solved: A*B*x = lambda*x;
//                      * 3, the following problem is solved: B*A*x = lambda*x.
//
// Outputs:
//     D           -   eigenvalues in ascending order.
//                     Array whose index ranges within [0..N-1].
//     Z           -   if ZNeeded is equal to:
//                      * 0, Z hasn't changed;
//                      * 1, Z contains eigenvectors.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//                     The eigenvectors are stored in matrix columns. It should
//                     be noted that the eigenvectors in such problems do not
//                     form an orthogonal system.
//
// Result:
//     True, if the problem was solved successfully.
//     False, if the error occurred during the Cholesky decomposition of matrix
//     B (the matrix isn't positive-definite) or during the work of the iterative
//     algorithm for solving the symmetric eigenproblem.
//
// See also the GeneralizedSymmetricDefiniteEVDReduce subroutine.
// ALGLIB: Copyright 01.28.2006 by Sergey Bochkanov
// API: bool smatrixgevd(const real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t zneeded, const ae_int_t problemtype, real_1d_array &d, real_2d_array &z);
bool smatrixgevd(RMatrix *a, ae_int_t n, bool isuppera, RMatrix *b, bool isupperb, ae_int_t zneeded, ae_int_t problemtype, RVector *d, RMatrix *z) {
   ae_frame _frame_block;
   bool isupperr;
   ae_int_t j1;
   ae_int_t j2;
   ae_int_t j1inc;
   ae_int_t j2inc;
   ae_int_t i;
   ae_int_t j;
   double v;
   bool result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   SetVector(d);
   SetMatrix(z);
   NewMatrix(r, 0, 0, DT_REAL);
   NewMatrix(t, 0, 0, DT_REAL);
// Reduce and solve
   result = smatrixgevdreduce(a, n, isuppera, b, isupperb, problemtype, &r, &isupperr);
   if (!result) {
      ae_frame_leave();
      return result;
   }
   result = smatrixevd(a, n, zneeded, isuppera, d, &t);
   if (!result) {
      ae_frame_leave();
      return result;
   }
// Transform eigenvectors if needed
   if (zneeded != 0) {
   // fill Z with zeros
      ae_matrix_set_length(z, n, n);
      for (j = 0; j < n; j++) {
         z->xyR[0][j] = 0.0;
      }
      for (i = 1; i < n; i++) {
         ae_v_move(z->xyR[i], 1, z->xyR[0], 1, n);
      }
   // Setup R properties
      if (isupperr) {
         j1 = 0;
         j2 = n - 1;
         j1inc = 1;
         j2inc = 0;
      } else {
         j1 = 0;
         j2 = 0;
         j1inc = 0;
         j2inc = 1;
      }
   // Calculate R*Z
      for (i = 0; i < n; i++) {
         for (j = j1; j <= j2; j++) {
            v = r.xyR[i][j];
            ae_v_addd(z->xyR[i], 1, t.xyR[j], 1, n, v);
         }
         j1 += j1inc;
         j2 += j2inc;
      }
   }
   ae_frame_leave();
   return result;
}

// Algorithm for reduction of the following generalized symmetric positive-
// definite eigenvalue problem:
//     A*x = lambda*B*x (1) or
//     A*B*x = lambda*x (2) or
//     B*A*x = lambda*x (3)
// to the symmetric eigenvalues problem C*y = lambda*y (eigenvalues of this and
// the given problems are the same, and the eigenvectors of the given problem
// could be obtained by multiplying the obtained eigenvectors by the
// transformation matrix x = R*y).
//
// Here A is a symmetric matrix, B - symmetric positive-definite matrix.
//
// Inputs:
//     A           -   symmetric matrix which is given by its upper or lower
//                     triangular part.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     N           -   size of matrices A and B.
//     IsUpperA    -   storage format of matrix A.
//     B           -   symmetric positive-definite matrix which is given by
//                     its upper or lower triangular part.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     IsUpperB    -   storage format of matrix B.
//     ProblemType -   if ProblemType is equal to:
//                      * 1, the following problem is solved: A*x = lambda*B*x;
//                      * 2, the following problem is solved: A*B*x = lambda*x;
//                      * 3, the following problem is solved: B*A*x = lambda*x.
//
// Outputs:
//     A           -   symmetric matrix which is given by its upper or lower
//                     triangle depending on IsUpperA. Contains matrix C.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     R           -   upper triangular or low triangular transformation matrix
//                     which is used to obtain the eigenvectors of a given problem
//                     as the product of eigenvectors of C (from the right) and
//                     matrix R (from the left). If the matrix is upper
//                     triangular, the elements below the main diagonal
//                     are equal to 0 (and vice versa). Thus, we can perform
//                     the multiplication without taking into account the
//                     internal structure (which is an easier though less
//                     effective way).
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     IsUpperR    -   type of matrix R (upper or lower triangular).
//
// Result:
//     True, if the problem was reduced successfully.
//     False, if the error occurred during the Cholesky decomposition of
//         matrix B (the matrix is not positive-definite).
// ALGLIB: Copyright 01.28.2006 by Sergey Bochkanov
// API: bool smatrixgevdreduce(real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t problemtype, real_2d_array &r, bool &isupperr);
bool smatrixgevdreduce(RMatrix *a, ae_int_t n, bool isuppera, RMatrix *b, bool isupperb, ae_int_t problemtype, RMatrix *r, bool *isupperr) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double v;
   ae_int_t info;
   bool result;
   ae_frame_make(&_frame_block);
   SetMatrix(r);
   *isupperr = false;
   NewMatrix(t, 0, 0, DT_REAL);
   NewVector(w1, 0, DT_REAL);
   NewVector(w2, 0, DT_REAL);
   NewVector(w3, 0, DT_REAL);
   NewObj(matinvreport, rep);
   ae_assert(n > 0, "SMatrixGEVDReduce: N <= 0!");
   ae_assert(problemtype == 1 || problemtype == 2 || problemtype == 3, "SMatrixGEVDReduce: incorrect ProblemType!");
   result = true;
// Problem 1:  A*x = lambda*B*x
//
// Reducing to:
//     C*y = lambda*y
//     C = L^(-1) * A * L^(-T)
//     x = L^(-T) * y
   if (problemtype == 1) {
   // Factorize B in T: B = LL'
      ae_matrix_set_length(&t, n, n);
      if (isupperb) {
         for (i = 0; i < n; i++) {
            ae_v_move(&t.xyR[i][i], t.stride, &b->xyR[i][i], 1, n - i);
         }
      } else {
         for (i = 0; i < n; i++) {
            ae_v_move(t.xyR[i], 1, b->xyR[i], 1, i + 1);
         }
      }
      if (!spdmatrixcholesky(&t, n, false)) {
         result = false;
         ae_frame_leave();
         return result;
      }
   // Invert L in T
      rmatrixtrinverse(&t, n, false, false, &info, &rep);
      if (info <= 0) {
         result = false;
         ae_frame_leave();
         return result;
      }
   // Build L^(-1) * A * L^(-T) in R
      ae_vector_set_length(&w1, n + 1);
      ae_vector_set_length(&w2, n + 1);
      ae_matrix_set_length(r, n, n);
      for (j = 1; j <= n; j++) {
      // Form w2 = A * l'(j) (here l'(j) is j-th column of L^(-T))
         ae_v_move(&w1.xR[1], 1, t.xyR[j - 1], 1, j);
         symmetricmatrixvectormultiply(a, isuppera, 0, j - 1, &w1, 1.0, &w2);
         if (isuppera) {
            matrixvectormultiply(a, 0, j - 1, j, n - 1, true, &w1, 1, j, 1.0, &w2, j + 1, n, 0.0);
         } else {
            matrixvectormultiply(a, j, n - 1, 0, j - 1, false, &w1, 1, j, 1.0, &w2, j + 1, n, 0.0);
         }
      // Form l(i)*w2 (here l(i) is i-th row of L^(-1))
         for (i = 1; i <= n; i++) {
            v = ae_v_dotproduct(t.xyR[i - 1], 1, &w2.xR[1], 1, i);
            r->xyR[i - 1][j - 1] = v;
         }
      }
   // Copy R to A
      for (i = 0; i < n; i++) {
         ae_v_move(a->xyR[i], 1, r->xyR[i], 1, n);
      }
   // Copy L^(-1) from T to R and transpose
      *isupperr = true;
      for (i = 0; i < n; i++) {
         for (j = 0; j < i; j++) {
            r->xyR[i][j] = 0.0;
         }
      }
      for (i = 0; i < n; i++) {
         ae_v_move(&r->xyR[i][i], 1, &t.xyR[i][i], t.stride, n - i);
      }
      ae_frame_leave();
      return result;
   }
// Problem 2:  A*B*x = lambda*x
// or
// problem 3:  B*A*x = lambda*x
//
// Reducing to:
//     C*y = lambda*y
//     C = U * A * U'
//     B = U'* U
   if (problemtype == 2 || problemtype == 3) {
   // Factorize B in T: B = U'*U
      ae_matrix_set_length(&t, n, n);
      if (isupperb) {
         for (i = 0; i < n; i++) {
            ae_v_move(&t.xyR[i][i], 1, &b->xyR[i][i], 1, n - i);
         }
      } else {
         for (i = 0; i < n; i++) {
            ae_v_move(&t.xyR[i][i], 1, &b->xyR[i][i], b->stride, n - i);
         }
      }
      if (!spdmatrixcholesky(&t, n, true)) {
         result = false;
         ae_frame_leave();
         return result;
      }
   // Build U * A * U' in R
      ae_vector_set_length(&w1, n + 1);
      ae_vector_set_length(&w2, n + 1);
      ae_vector_set_length(&w3, n + 1);
      ae_matrix_set_length(r, n, n);
      for (j = 1; j <= n; j++) {
      // Form w2 = A * u'(j) (here u'(j) is j-th column of U')
         ae_v_move(&w1.xR[1], 1, &t.xyR[j - 1][j - 1], 1, n - j + 1);
         symmetricmatrixvectormultiply(a, isuppera, j - 1, n - 1, &w1, 1.0, &w3);
         ae_v_move(&w2.xR[j], 1, &w3.xR[1], 1, n - j + 1);
         ae_v_move(&w1.xR[j], 1, &t.xyR[j - 1][j - 1], 1, n - j + 1);
         if (isuppera) {
            matrixvectormultiply(a, 0, j - 2, j - 1, n - 1, false, &w1, j, n, 1.0, &w2, 1, j - 1, 0.0);
         } else {
            matrixvectormultiply(a, j - 1, n - 1, 0, j - 2, true, &w1, j, n, 1.0, &w2, 1, j - 1, 0.0);
         }
      // Form u(i)*w2 (here u(i) is i-th row of U)
         for (i = 1; i <= n; i++) {
            v = ae_v_dotproduct(&t.xyR[i - 1][i - 1], 1, &w2.xR[i], 1, n - i + 1);
            r->xyR[i - 1][j - 1] = v;
         }
      }
   // Copy R to A
      for (i = 0; i < n; i++) {
         ae_v_move(a->xyR[i], 1, r->xyR[i], 1, n);
      }
      if (problemtype == 2) {
      // Invert U in T
         rmatrixtrinverse(&t, n, true, false, &info, &rep);
         if (info <= 0) {
            result = false;
            ae_frame_leave();
            return result;
         }
      // Copy U^-1 from T to R
         *isupperr = true;
         for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
               r->xyR[i][j] = 0.0;
            }
         }
         for (i = 0; i < n; i++) {
            ae_v_move(&r->xyR[i][i], 1, &t.xyR[i][i], 1, n - i);
         }
      } else {
      // Copy U from T to R and transpose
         *isupperr = false;
         for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
               r->xyR[i][j] = 0.0;
            }
         }
         for (i = 0; i < n; i++) {
            ae_v_move(&r->xyR[i][i], r->stride, &t.xyR[i][i], 1, n - i);
         }
      }
   }
   ae_frame_leave();
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
bool smatrixgevd(const real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t zneeded, const ae_int_t problemtype, real_1d_array &d, real_2d_array &z) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixgevd(ConstT(ae_matrix, a), n, isuppera, ConstT(ae_matrix, b), isupperb, zneeded, problemtype, ConstT(ae_vector, d), ConstT(ae_matrix, z));
   alglib_impl::ae_state_clear();
   return Ok;
}

bool smatrixgevdreduce(real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t problemtype, real_2d_array &r, bool &isupperr) {
   alglib_impl::ae_state_init();
   TryCatch(false)
   bool Ok = alglib_impl::smatrixgevdreduce(ConstT(ae_matrix, a), n, isuppera, ConstT(ae_matrix, b), isupperb, problemtype, ConstT(ae_matrix, r), &isupperr);
   alglib_impl::ae_state_clear();
   return Ok;
}
} // end of namespace alglib

// === INVERSEUPDATE Package ===
namespace alglib_impl {
// Inverse matrix update by the Sherman-Morrison formula
//
// The algorithm updates matrix A^-1 when adding a number to an element
// of matrix A.
//
// Inputs:
//     InvA    -   inverse of matrix A.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     UpdRow  -   row where the element to be updated is stored.
//     UpdColumn - column where the element to be updated is stored.
//     UpdVal  -   a number to be added to the element.
//
//
// Outputs:
//     InvA    -   inverse of modified matrix A.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: void rmatrixinvupdatesimple(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const ae_int_t updcolumn, const double updval);
void rmatrixinvupdatesimple(RMatrix *inva, ae_int_t n, ae_int_t updrow, ae_int_t updcolumn, double updval) {
   ae_frame _frame_block;
   ae_int_t i;
   double lambdav;
   double vt;
   ae_frame_make(&_frame_block);
   NewVector(t1, 0, DT_REAL);
   NewVector(t2, 0, DT_REAL);
   ae_assert(updrow >= 0 && updrow < n, "RMatrixInvUpdateSimple: incorrect UpdRow!");
   ae_assert(updcolumn >= 0 && updcolumn < n, "RMatrixInvUpdateSimple: incorrect UpdColumn!");
   ae_vector_set_length(&t1, n);
   ae_vector_set_length(&t2, n);
// T1 = InvA * U
   ae_v_move(t1.xR, 1, &inva->xyR[0][updrow], inva->stride, n);
// T2 = v*InvA
   ae_v_move(t2.xR, 1, inva->xyR[updcolumn], 1, n);
// Lambda = v * InvA * U
   lambdav = updval * inva->xyR[updcolumn][updrow];
// InvA = InvA - correction
   for (i = 0; i < n; i++) {
      vt = updval * t1.xR[i];
      vt /= 1 + lambdav;
      ae_v_subd(inva->xyR[i], 1, t2.xR, 1, n, vt);
   }
   ae_frame_leave();
}

// Inverse matrix update by the Sherman-Morrison formula
//
// The algorithm updates matrix A^-1 when adding a vector to a row
// of matrix A.
//
// Inputs:
//     InvA    -   inverse of matrix A.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     UpdRow  -   the row of A whose vector V was added.
//                 0 <= Row <= N-1
//     V       -   the vector to be added to a row.
//                 Array whose index ranges within [0..N-1].
//
// Outputs:
//     InvA    -   inverse of modified matrix A.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: void rmatrixinvupdaterow(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const real_1d_array &v);
void rmatrixinvupdaterow(RMatrix *inva, ae_int_t n, ae_int_t updrow, RVector *v) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double lambdav;
   double vt;
   ae_frame_make(&_frame_block);
   NewVector(t1, 0, DT_REAL);
   NewVector(t2, 0, DT_REAL);
   ae_vector_set_length(&t1, n);
   ae_vector_set_length(&t2, n);
// T1 = InvA * U
   ae_v_move(t1.xR, 1, &inva->xyR[0][updrow], inva->stride, n);
// T2 = v*InvA
// Lambda = v * InvA * U
   for (j = 0; j < n; j++) {
      vt = ae_v_dotproduct(v->xR, 1, &inva->xyR[0][j], inva->stride, n);
      t2.xR[j] = vt;
   }
   lambdav = t2.xR[updrow];
// InvA = InvA - correction
   for (i = 0; i < n; i++) {
      vt = t1.xR[i] / (1 + lambdav);
      ae_v_subd(inva->xyR[i], 1, t2.xR, 1, n, vt);
   }
   ae_frame_leave();
}

// Inverse matrix update by the Sherman-Morrison formula
//
// The algorithm updates matrix A^-1 when adding a vector to a column
// of matrix A.
//
// Inputs:
//     InvA        -   inverse of matrix A.
//                     Array whose indexes range within [0..N-1, 0..N-1].
//     N           -   size of matrix A.
//     UpdColumn   -   the column of A whose vector U was added.
//                     0 <= UpdColumn <= N-1
//     U           -   the vector to be added to a column.
//                     Array whose index ranges within [0..N-1].
//
// Outputs:
//     InvA        -   inverse of modified matrix A.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: void rmatrixinvupdatecolumn(real_2d_array &inva, const ae_int_t n, const ae_int_t updcolumn, const real_1d_array &u);
void rmatrixinvupdatecolumn(RMatrix *inva, ae_int_t n, ae_int_t updcolumn, RVector *u) {
   ae_frame _frame_block;
   ae_int_t i;
   double lambdav;
   double vt;
   ae_frame_make(&_frame_block);
   NewVector(t1, 0, DT_REAL);
   NewVector(t2, 0, DT_REAL);
   ae_vector_set_length(&t1, n);
   ae_vector_set_length(&t2, n);
// T1 = InvA * U
// Lambda = v * InvA * U
   for (i = 0; i < n; i++) {
      vt = ae_v_dotproduct(inva->xyR[i], 1, u->xR, 1, n);
      t1.xR[i] = vt;
   }
   lambdav = t1.xR[updcolumn];
// T2 = v*InvA
   ae_v_move(t2.xR, 1, inva->xyR[updcolumn], 1, n);
// InvA = InvA - correction
   for (i = 0; i < n; i++) {
      vt = t1.xR[i] / (1 + lambdav);
      ae_v_subd(inva->xyR[i], 1, t2.xR, 1, n, vt);
   }
   ae_frame_leave();
}

// Inverse matrix update by the Sherman-Morrison formula
//
// The algorithm computes the inverse of matrix A+u*v' by using the given matrix
// A^-1 and the vectors u and v.
//
// Inputs:
//     InvA    -   inverse of matrix A.
//                 Array whose indexes range within [0..N-1, 0..N-1].
//     N       -   size of matrix A.
//     U       -   the vector modifying the matrix.
//                 Array whose index ranges within [0..N-1].
//     V       -   the vector modifying the matrix.
//                 Array whose index ranges within [0..N-1].
//
// Outputs:
//     InvA - inverse of matrix A + u*v'.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: void rmatrixinvupdateuv(real_2d_array &inva, const ae_int_t n, const real_1d_array &u, const real_1d_array &v);
void rmatrixinvupdateuv(RMatrix *inva, ae_int_t n, RVector *u, RVector *v) {
   ae_frame _frame_block;
   ae_int_t i;
   ae_int_t j;
   double lambdav;
   double vt;
   ae_frame_make(&_frame_block);
   NewVector(t1, 0, DT_REAL);
   NewVector(t2, 0, DT_REAL);
   ae_vector_set_length(&t1, n);
   ae_vector_set_length(&t2, n);
// T1 = InvA * U
// Lambda = v * T1
   for (i = 0; i < n; i++) {
      vt = ae_v_dotproduct(inva->xyR[i], 1, u->xR, 1, n);
      t1.xR[i] = vt;
   }
   lambdav = ae_v_dotproduct(v->xR, 1, t1.xR, 1, n);
// T2 = v*InvA
   for (j = 0; j < n; j++) {
      vt = ae_v_dotproduct(v->xR, 1, &inva->xyR[0][j], inva->stride, n);
      t2.xR[j] = vt;
   }
// InvA = InvA - correction
   for (i = 0; i < n; i++) {
      vt = t1.xR[i] / (1 + lambdav);
      ae_v_subd(inva->xyR[i], 1, t2.xR, 1, n, vt);
   }
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void rmatrixinvupdatesimple(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const ae_int_t updcolumn, const double updval) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixinvupdatesimple(ConstT(ae_matrix, inva), n, updrow, updcolumn, updval);
   alglib_impl::ae_state_clear();
}

void rmatrixinvupdaterow(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const real_1d_array &v) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixinvupdaterow(ConstT(ae_matrix, inva), n, updrow, ConstT(ae_vector, v));
   alglib_impl::ae_state_clear();
}

void rmatrixinvupdatecolumn(real_2d_array &inva, const ae_int_t n, const ae_int_t updcolumn, const real_1d_array &u) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixinvupdatecolumn(ConstT(ae_matrix, inva), n, updcolumn, ConstT(ae_vector, u));
   alglib_impl::ae_state_clear();
}

void rmatrixinvupdateuv(real_2d_array &inva, const ae_int_t n, const real_1d_array &u, const real_1d_array &v) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::rmatrixinvupdateuv(ConstT(ae_matrix, inva), n, ConstT(ae_vector, u), ConstT(ae_vector, v));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === MATDET Package ===
// Depends on: TRFAC
namespace alglib_impl {
// Determinant calculation of the matrix given by its LU decomposition.
//
// Inputs:
//     A       -   LU decomposition of the matrix (output of
//                 RMatrixLU subroutine).
//     Pivots  -   table of permutations which were made during
//                 the LU decomposition.
//                 Output of RMatrixLU subroutine.
//     N       -   (optional) size of matrix A:
//                 * if given, only principal NxN submatrix is processed and
//                   overwritten. other elements are unchanged.
//                 * if not given, automatically determined from matrix size
//                   (A must be square matrix)
//
// Result: matrix determinant.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n);
// API: double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots);
double rmatrixludet(RMatrix *a, ZVector *pivots, ae_int_t n) {
   ae_int_t i;
   ae_int_t s;
   double result;
   ae_assert(n >= 1, "RMatrixLUDet: N<1!");
   ae_assert(pivots->cnt >= n, "RMatrixLUDet: Pivots array is too short!");
   ae_assert(a->rows >= n, "RMatrixLUDet: rows(A)<N!");
   ae_assert(a->cols >= n, "RMatrixLUDet: cols(A)<N!");
   ae_assert(apservisfinitematrix(a, n, n), "RMatrixLUDet: A contains infinite or NaN values!");
   result = 1.0;
   s = 1;
   for (i = 0; i < n; i++) {
      result *= a->xyR[i][i];
      if (pivots->xZ[i] != i) {
         s = -s;
      }
   }
   result *= s;
   return result;
}

// Calculation of the determinant of a general matrix
//
// Inputs:
//     A       -   matrix, array[0..N-1, 0..N-1]
//     N       -   (optional) size of matrix A:
//                 * if given, only principal NxN submatrix is processed and
//                   overwritten. other elements are unchanged.
//                 * if not given, automatically determined from matrix size
//                   (A must be square matrix)
//
// Result: determinant of matrix A.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: double rmatrixdet(const real_2d_array &a, const ae_int_t n);
// API: double rmatrixdet(const real_2d_array &a);
double rmatrixdet(RMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n >= 1, "RMatrixDet: N<1!");
   ae_assert(a->rows >= n, "RMatrixDet: rows(A)<N!");
   ae_assert(a->cols >= n, "RMatrixDet: cols(A)<N!");
   ae_assert(apservisfinitematrix(a, n, n), "RMatrixDet: A contains infinite or NaN values!");
   rmatrixlu(a, n, n, &pivots);
   result = rmatrixludet(a, &pivots, n);
   ae_frame_leave();
   return result;
}

// Determinant calculation of the matrix given by its LU decomposition.
//
// Inputs:
//     A       -   LU decomposition of the matrix (output of
//                 RMatrixLU subroutine).
//     Pivots  -   table of permutations which were made during
//                 the LU decomposition.
//                 Output of RMatrixLU subroutine.
//     N       -   (optional) size of matrix A:
//                 * if given, only principal NxN submatrix is processed and
//                   overwritten. other elements are unchanged.
//                 * if not given, automatically determined from matrix size
//                   (A must be square matrix)
//
// Result: matrix determinant.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n);
// API: complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots);
ae_complex cmatrixludet(CMatrix *a, ZVector *pivots, ae_int_t n) {
   ae_int_t i;
   ae_int_t s;
   ae_complex result;
   ae_assert(n >= 1, "CMatrixLUDet: N<1!");
   ae_assert(pivots->cnt >= n, "CMatrixLUDet: Pivots array is too short!");
   ae_assert(a->rows >= n, "CMatrixLUDet: rows(A)<N!");
   ae_assert(a->cols >= n, "CMatrixLUDet: cols(A)<N!");
   ae_assert(apservisfinitecmatrix(a, n, n), "CMatrixLUDet: A contains infinite or NaN values!");
   result = ae_complex_from_i(1);
   s = 1;
   for (i = 0; i < n; i++) {
      result = ae_c_mul(result, a->xyC[i][i]);
      if (pivots->xZ[i] != i) {
         s = -s;
      }
   }
   result = ae_c_mul_d(result, (double)s);
   return result;
}

// Calculation of the determinant of a general matrix
//
// Inputs:
//     A       -   matrix, array[0..N-1, 0..N-1]
//     N       -   (optional) size of matrix A:
//                 * if given, only principal NxN submatrix is processed and
//                   overwritten. other elements are unchanged.
//                 * if not given, automatically determined from matrix size
//                   (A must be square matrix)
//
// Result: determinant of matrix A.
// ALGLIB: Copyright 2005 by Sergey Bochkanov
// API: complex cmatrixdet(const complex_2d_array &a, const ae_int_t n);
// API: complex cmatrixdet(const complex_2d_array &a);
ae_complex cmatrixdet(CMatrix *a, ae_int_t n) {
   ae_frame _frame_block;
   ae_complex result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   NewVector(pivots, 0, DT_INT);
   ae_assert(n >= 1, "CMatrixDet: N<1!");
   ae_assert(a->rows >= n, "CMatrixDet: rows(A)<N!");
   ae_assert(a->cols >= n, "CMatrixDet: cols(A)<N!");
   ae_assert(apservisfinitecmatrix(a, n, n), "CMatrixDet: A contains infinite or NaN values!");
   cmatrixlu(a, n, n, &pivots);
   result = cmatrixludet(a, &pivots, n);
   ae_frame_leave();
   return result;
}

// Determinant calculation of the matrix given by the Cholesky decomposition.
//
// Inputs:
//     A       -   Cholesky decomposition,
//                 output of SMatrixCholesky subroutine.
//     N       -   (optional) size of matrix A:
//                 * if given, only principal NxN submatrix is processed and
//                   overwritten. other elements are unchanged.
//                 * if not given, automatically determined from matrix size
//                   (A must be square matrix)
//
// As the determinant is equal to the product of squares of diagonal elements,
// it's not necessary to specify which triangle - lower or upper - the matrix
// is stored in.
//
// Result:
//     matrix determinant.
// ALGLIB: Copyright 2005-2008 by Sergey Bochkanov
// API: double spdmatrixcholeskydet(const real_2d_array &a, const ae_int_t n);
// API: double spdmatrixcholeskydet(const real_2d_array &a);
double spdmatrixcholeskydet(RMatrix *a, ae_int_t n) {
   ae_int_t i;
   bool f;
   double result;
   ae_assert(n >= 1, "SPDMatrixCholeskyDet: N<1!");
   ae_assert(a->rows >= n, "SPDMatrixCholeskyDet: rows(A)<N!");
   ae_assert(a->cols >= n, "SPDMatrixCholeskyDet: cols(A)<N!");
   f = true;
   for (i = 0; i < n; i++) {
      f = f && isfinite(a->xyR[i][i]);
   }
   ae_assert(f, "SPDMatrixCholeskyDet: A contains infinite or NaN values!");
   result = 1.0;
   for (i = 0; i < n; i++) {
      result *= ae_sqr(a->xyR[i][i]);
   }
   return result;
}

// Determinant calculation of the symmetric positive definite matrix.
//
// Inputs:
//     A       -   matrix. Array with elements [0..N-1, 0..N-1].
//     N       -   (optional) size of matrix A:
//                 * if given, only principal NxN submatrix is processed and
//                   overwritten. other elements are unchanged.
//                 * if not given, automatically determined from matrix size
//                   (A must be square matrix)
//     IsUpper -   (optional) storage type:
//                 * if True, symmetric matrix  A  is  given  by  its  upper
//                   triangle, and the lower triangle isn't used/changed  by
//                   function
//                 * if False, symmetric matrix  A  is  given  by  its lower
//                   triangle, and the upper triangle isn't used/changed  by
//                   function
//                 * if not given, both lower and upper  triangles  must  be
//                   filled.
//
// Result:
//     determinant of matrix A.
//     If matrix A is not positive definite, exception is thrown.
// ALGLIB: Copyright 2005-2008 by Sergey Bochkanov
// API: double spdmatrixdet(const real_2d_array &a, const ae_int_t n, const bool isupper);
// API: double spdmatrixdet(const real_2d_array &a);
double spdmatrixdet(RMatrix *a, ae_int_t n, bool isupper) {
   ae_frame _frame_block;
   bool b;
   double result;
   ae_frame_make(&_frame_block);
   DupMatrix(a);
   ae_assert(n >= 1, "SPDMatrixDet: N<1!");
   ae_assert(a->rows >= n, "SPDMatrixDet: rows(A)<N!");
   ae_assert(a->cols >= n, "SPDMatrixDet: cols(A)<N!");
   ae_assert(isfinitertrmatrix(a, n, isupper), "SPDMatrixDet: A contains infinite or NaN values!");
   b = spdmatrixcholesky(a, n, isupper);
   ae_assert(b, "SPDMatrixDet: A is not SPD!");
   result = spdmatrixcholeskydet(a, n);
   ae_frame_leave();
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixludet(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots) {
   if (a.rows() != a.cols() || a.rows() != pivots.length()) ThrowError("Error while calling 'rmatrixludet': looks like one of arguments has wrong size");
   ae_int_t n = a.rows();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixludet(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double rmatrixdet(const real_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixdet(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double rmatrixdet(const real_2d_array &a) {
   if (a.rows() != a.cols()) ThrowError("Error while calling 'rmatrixdet': looks like one of arguments has wrong size");
   ae_int_t n = a.rows();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::rmatrixdet(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(complex(0.0))
   alglib_impl::ae_complex C = alglib_impl::cmatrixludet(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n);
   alglib_impl::ae_state_clear();
   return ComplexOf(C);
}
#if !defined AE_NO_EXCEPTIONS
complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots) {
   if (a.rows() != a.cols() || a.rows() != pivots.length()) ThrowError("Error while calling 'cmatrixludet': looks like one of arguments has wrong size");
   ae_int_t n = a.rows();
   alglib_impl::ae_state_init();
   TryCatch(complex(0.0))
   alglib_impl::ae_complex C = alglib_impl::cmatrixludet(ConstT(ae_matrix, a), ConstT(ae_vector, pivots), n);
   alglib_impl::ae_state_clear();
   return ComplexOf(C);
}
#endif

complex cmatrixdet(const complex_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(complex(0.0))
   alglib_impl::ae_complex C = alglib_impl::cmatrixdet(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return ComplexOf(C);
}
#if !defined AE_NO_EXCEPTIONS
complex cmatrixdet(const complex_2d_array &a) {
   if (a.rows() != a.cols()) ThrowError("Error while calling 'cmatrixdet': looks like one of arguments has wrong size");
   ae_int_t n = a.rows();
   alglib_impl::ae_state_init();
   TryCatch(complex(0.0))
   alglib_impl::ae_complex C = alglib_impl::cmatrixdet(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return ComplexOf(C);
}
#endif

double spdmatrixcholeskydet(const real_2d_array &a, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spdmatrixcholeskydet(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double spdmatrixcholeskydet(const real_2d_array &a) {
   if (a.rows() != a.cols()) ThrowError("Error while calling 'spdmatrixcholeskydet': looks like one of arguments has wrong size");
   ae_int_t n = a.rows();
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spdmatrixcholeskydet(ConstT(ae_matrix, a), n);
   alglib_impl::ae_state_clear();
   return D;
}
#endif

double spdmatrixdet(const real_2d_array &a, const ae_int_t n, const bool isupper) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spdmatrixdet(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return D;
}
#if !defined AE_NO_EXCEPTIONS
double spdmatrixdet(const real_2d_array &a) {
   if (a.rows() != a.cols()) ThrowError("Error while calling 'spdmatrixdet': looks like one of arguments has wrong size");
   if (!alglib_impl::ae_is_symmetric(ConstT(ae_matrix, a))) ThrowError("'a' parameter is not symmetric matrix");
   ae_int_t n = a.rows();
   bool isupper = false;
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::spdmatrixdet(ConstT(ae_matrix, a), n, isupper);
   alglib_impl::ae_state_clear();
   return D;
}
#endif
} // end of namespace alglib
