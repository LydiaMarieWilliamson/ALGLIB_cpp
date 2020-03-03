// ALGLIB++
// Copyright (c) Lydia Marie Williamson, Mark Hopkins Consulting
// Based on ALGLIB: Copyright (c) Sergey Bochkanov (ALGLIB project).
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
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <sys/time.h> // For Ticks()
#include <string.h>
#include <stdint.h>
#include <ctype.h>

// === Ap Infrastructure ===
// Basic shared C/C++ functionality: like memory management for vectors/matrices, etc.
// Define Integer, StatusT and TypeT.
// A boolean type was also originally (and unnecessarily) defined.
// For C (as of 2011) one needs only to include <stdbool.h> to define "bool", "false" and "true",
// For C++: it is already a part of the language.
typedef ptrdiff_t Integer;
typedef enum { OkST, NoMemoryST, BadAssertST } StatusT;

#define MaxPosReal 1E300

typedef void (*FreeFnT)(void *);

// Dynamic block which may be automatically deallocated during stack unwinding.
typedef void *BlockS, **Block, *Frame;

// Abnormally abort the program, using one of several ways:
// *	else abort().
// Clear the frame stack, if any.
static void StateBreak(StatusT Status, const char *Msg) {
   const char *StatS = "";
   switch (Status) {
      default: case OkST: break;
      case NoMemoryST: StatS = "Out of memory"; break;
      case BadAssertST: StatS = "Assertion failure"; break;
   }
   fprintf(stderr, "Fatal (%s): %s\n", StatS, Msg), abort();
}

// Assertion
// Upon failure abort().
// IMPORTANT:
// *	These functions always evaluate Cond, and cannot be replaced by macros which do nothing.
//	In particular, when invoked, a function call may be used as the Cond argument, and it will be carried out.
void Assert(bool Cond, const char *Msg) { if (!Cond) StateBreak(BadAssertST, Msg); }
static void Impose(bool Cond, const char *Msg) { if (!Cond) fprintf(stderr, Msg), abort(); }
#define Demand(Cond) if (!(Cond)) abort()

// Allocate memory with automatic alignment.
// Return NULL when size N == 0 is specified.
// Upon failure return NULL.
void *Allocate(size_t N) { return N == 0? NULL: malloc((size_t)N); }
void *ReAllocate(void *X, size_t N) {
   if (X != NULL) free(X), X = NULL;
   return N == 0? NULL: malloc((size_t)N);
}
void Free(void *X) { if (X != NULL) free(X); }

template<typename T> struct Vector {
   Vector<T>(Integer N = 0) { InitVector(N); }
   Vector<T>(const Vector<T> &A) { CopyVector(A); }
   ~Vector<T>() { FreeVector(); }
// A new vector of size N >= 0.
// Upon allocation failure or N < 0 call StateBreak().
// NOTE:
// *	no memory allocation is performed for initialization with N == 0.
   void InitVector(Integer N = 0) { _N = N, _X = (T *)Allocate((size_t)(N*sizeof _X[0])); }
// A copy of vector A.
// Upon allocation failure call StateBreak().
   void CopyVector(const Vector<T> &A) { InitVector(A._N); if (A._N > 0) memmove(_X, A._X, (size_t)(A._N*sizeof(T))); }
// Resize a previously-initialized vector to size N >= 0.
   void ReSizeVector(Integer N) { if (_N != N) _N = N, _X = (T *)ReAllocate(_X, N*sizeof _X[0]); }
// The Free functionality for vectors (clearing contents and freeing all internal structures).
// Clear the vector (releasing all dynamically allocated memory), but leave the structure intact in state.
// The vector may be on the frame - in which case it will NOT be removed from list.
// IMPORTANT:
// *	This function does NOT invalidate the vector; it just releases all dynamically allocated storage, but it still may be used afterwards.
   void FreeVector() { _N = 0, Free(_X), _X = NULL; }
// === ApServ Package ===
// Make the vector have length at least N. (One routine for each vector type).
// AlgLib: Copyright 20.03.2009 by Sergey Bochkanov
   void GrowVector(Integer N) { if (_N < N) ReSizeVector(N); }
// Resize the vector, preserving its contents and zero-padding the area expanded into.
// AlgLib: Copyright 20.03.2009 by Sergey Bochkanov
   void ExpandVector(Integer N) {
      Vector<T> ExA(N);
      Integer ExN = _N; _N = ExA._N, ExA._N = ExN;
      T *ExX = _X; _X = ExA._X, ExA._X = ExX;
      T *X = _X; for (Integer n = 0; n < N; n++) X[n] = n < ExN? ExX[n]: 0;
   }
// private:
   Integer _N; T *_X;
};

typedef Vector<Integer> ZVector;
typedef Vector<double> RVector;

template<typename T> struct Matrix {
   Matrix<T>(Integer Ys = 0, Integer Xs = 0) { InitMatrix(Ys, Xs); }
   Matrix<T>(const Matrix<T> &A) { CopyMatrix(A); }
   ~Matrix<T>() { FreeMatrix(); }
// A new Xs x Ys type matrix.
// Xs and Ys may be 0, but must both be, if so.
// Upon allocation failure call StateBreak().
// NOTE:
// *	No memory allocation is performed for initialization with Ys = Xs = 0.
   void InitMatrix(Integer Ys = 0, Integer Xs = 0) {
      Assert(Xs >= 0 && Ys >= 0, "Matrix::InitMatrix: negative size");
   // If either Xs or Ys is 0, then they both must be made so.
      if (Xs == 0) Ys = 0; else if (Ys == 0) Xs = 0;
   // Initialize.
      _Xs = Xs, _Ys = Ys;
   // If Xs and Ys are 0; perform a quick exit.
      if (Xs == 0 || Ys == 0) { _XY = NULL; return; }
      _XY = (T **)Allocate((size_t)(Ys*sizeof _XY[0]));
      T *Map = (T *)Allocate((size_t)(Xs*Ys*sizeof Map[0]));
   // Set pointers to the matrix rows.
      for (Integer Y = 0; Y < _Ys; Map += _Xs, Y++) _XY[Y] = Map;
   }
// A copy of matrix A.
// Upon allocation failure call StateBreak().
   void CopyMatrix(const Matrix<T> &A) {
      Integer Xs = A._Xs, Ys = A._Ys; T *Map = A._XY[0];
      InitMatrix(Ys, Xs);
      if (Xs > 0 && Ys > 0) memmove(_XY[0], Map, (size_t)(Xs*Ys*sizeof Map[0]));
   }
// Resize matrix to size Xs x Ys, where either Xs, Ys or both may be 0.
// Upon allocation failure return an indication of success or failure.
   void ReSizeMatrix(Integer Ys, Integer Xs) {
      Assert(Xs >= 0 && Ys >= 0, "Matrix::ReSizeMatrix: negative size");
      if (_Xs == Xs && _Ys == Ys) return;
   // Re-initialize.
      _Xs = Xs, _Ys = Ys;
      T *Map = _XY == NULL? NULL: _XY[0];
      _XY = (T **)ReAllocate(_XY, (size_t)(Ys*sizeof _XY[0]));
      Map = (T *)ReAllocate(Map, (size_t)(Xs*Ys*sizeof Map[0]));
   // Set pointers to the matrix rows.
      for (Integer Y = 0; Y < _Ys; Map += _Xs, Y++) _XY[Y] = Map;
   }
// The Free functionality for matrices (clearing contents but leaving the structure in a valid state).
// Clear the contents of the matrix, releasing all dynamically allocated memory, but leaving the structure intact in state.
// The matrix may be in automatic management list - in this case it will NOT be removed from list.
// It still may be used after call to ReSizeMatrix().
// IMPORTANT:
// *	This function does NOT invalidate the matrix; it just releases all dynamically allocated storage, but it still may be used afterwards.
   void FreeMatrix() {
      if (_XY != NULL) Free(_XY[0]);
      _Ys = _Xs = 0, Free(_XY), _XY = NULL;
   }
// === ApServ Package ===
// Make the matrix have at least Xs x Ys size.
// AlgLib: Copyright 20.03.2009 by Sergey Bochkanov
   void GrowMatrix(Integer Ys, Integer Xs) { if (Xs > 0 && Ys > 0 && (_Xs < Xs || _Ys < Ys)) ReSizeMatrix(Ys, Xs); }
// private:
   Integer _Xs, _Ys; T **_XY;
};

typedef Matrix<double> RMatrix;

// Real math functions.
template<typename T> T Max(T A, T B) { return A > B? A: B; }
template<typename T> T Min(T A, T B) { return A < B? A: B; }
template<typename T> T Max4(T A, T B, T C, T D) { T AB = A > B? A: B, CD = C > D? C: D; return AB > CD? AB: CD; }
template<typename T> T Squ(T X) { return X*X; }
template<typename T> void CopyV(T A[], const T B[], Integer N) {
   Integer Nq = N/2, Nr = N%2;
   for (Integer q = 0; q < Nq; A += 2, B += 2, q++) A[0] = B[0], A[1] = B[1];
   if (Nr) A[0] = B[0];
}
template<typename T> void CopyV(T A[], Integer dA, const T B[], Integer dB, Integer N) {
   if (dA == 1 && dB == 1) CopyV(A, B, N);
   else for (Integer n = 0; n < N; A += dA, B += dB, n++) *A = *B;
}

int Ticks() {
   struct timeval Now; gettimeofday(&Now, NULL);
   return 1000*Now.tv_sec + Now.tv_usec/1000;
}

// === ApServ Package ===
// Check that A has at least length N and that its first N values are finite.
// AlgLib: Copyright 18.06.2010 by Sergey Bochkanov
bool IsFiniteRVector(const RVector &A, const Integer N) {
   Assert(N >= 0, "IsFiniteRVector: internal error (N < 0)");
   if (N == 0) return true; else if (A._N < N) return false;
   double *AP = A._X;
// (#) This is used, instead, from Version 3.15.0 on:
// double V;
// for (Integer n = 0; n < N; n++) V = 0.01*V + AP[n];
// return isfinite(V);
   for (Integer n = 0; n < N; n++) if (!isfinite(AP[n])) return false;
   return true;
}

// Check that A has at least Xs x Ys and its values over [0, Xs) x [0, Ys) are all finite.
// AlgLib: Copyright 18.06.2010 by Sergey Bochkanov
bool IsFiniteRMatrix(const RMatrix &A, const Integer Ys, const Integer Xs) {
   Assert(Xs >= 0, "IsFiniteRMatrix: internal error (Xs < 0)");
   Assert(Ys >= 0, "IsFiniteRMatrix: internal error (Ys < 0)");
   if (Xs == 0 || Ys == 0) return true; else if (A._Xs < Xs || A._Ys < Ys) return false;
   double **AP = A._XY;
   for (Integer Y = 0; Y < Ys; Y++) {
      double *Row = AP[Y];
      for (Integer X = 0; X < Xs; X++) if (!isfinite(Row[X])) return false;
   }
   return true;
}

// Swap the values A and B.
template<typename T> void Swap(T &A, T &B) { T ExA = A; A = B, B = ExA; }

// === TSort Package ===
// Heap operations: add an element to the heap.
// Inputs:
//	A, B:	The heap: a vector of size at least N + 1 and its permutation tags.
//	N:	The size of the heap, without the new element.
//	Av, Bv:	The value and tag of the element being added.
// Output:
//	N:	The size of the updated heap.
// AlgLib: Copyright 28.02.2010 by Sergey Bochkanov
void TagHeapPush(RVector &A, ZVector &B, Integer &N, double Av, Integer Bv) {
   if (N < 0) return;
   double *AP = A._X; Integer *BP = B._X;
// Add the current point to the heap (add to the bottom, then move everything else up).
// We dont' write the point to the heap until its final position is determined.
// (so as to reduce the number of array access operations).
   Integer n = N++;
   while (n > 0) {
      Integer nL = (n - 1)/2;
      double a = AP[nL];
      if (a < Av) AP[n] = a, BP[n] = BP[nL], n = nL; // Swap with the higher element.
      else break; // The element is in its place. terminate.
   }
   AP[n] = Av, BP[n] = Bv;
}

// Heap operations: replace the top element with a new element (which is moved down)
// Parameters:
//	A:	The heap: a vector of size at least N.
//	B:	The vector of integer tags, which are updated according to permutations in the heap.
//	N:	The size of the heap.
//	Av:	The value of the element which replaces top element
//	Bv:	The value of the tag
// AlgLib: Copyright 28.02.2010 by Sergey Bochkanov
void TagHeapReTop(RVector &A, ZVector &B, Integer N, double Av, Integer Bv) {
   if (N < 1) return;
   double *AP = A._X; Integer *BP = B._X;
// Move down through heap:
// *	n:	The current element
// *	n0:	The firstborn.
// We don't write the point to the heap until its final position is determined
// (so as to reduce the number of array access operations).
   Integer n = 0;
   for (Integer n0 = 1; n0 < N; n0 = 2*n + 1) {
      if (n0 + 1 >= N) { // Only one child: swap and terminate (because this child has no siblings due to heap structure)
         double a = AP[n0];
         if (a > Av) AP[n] = a, BP[n] = BP[n0], n = n0;
         break;
      }
   // Two children.
      Integer nL = AP[n0] > AP[n0 + 1]? n0: n0 + 1;
      if (Av >= AP[nL]) break;
      AP[n] = AP[nL], BP[n] = BP[nL], n = nL;
   }
   AP[n] = Av, BP[n] = Bv;
}

// Heap operations: pops top element from the heap
// Inputs:
//	A:	The heap itself, a vector of size at least N.
//	B:	The vector of integer tags, which are updated according to permutations in the heap.
//	N:	The size of the heap, N > 0.
// On output the top element is moved to A[N-1], B[N-1], heap is reordered, N is decreased by 1.
// AlgLib: Copyright 28.02.2010 by Sergey Bochkanov
void TagHeapPop(RVector &A, ZVector &B, Integer &N) {
   if (N < 1) return;
   else if (--N > 0) {
      double *AP = A._X; Integer *BP = B._X;
   // Swap the top and last elements, then reorder the heap.
      double Av = AP[N]; Integer Bv = BP[N]; AP[N] = AP[0], BP[N] = BP[0], TagHeapReTop(A, B, N, Av, Bv);
   }
}

// === KdTree Package ===
struct KdTree;

struct KdTreeReq {
   KdTreeReq() { InitKdTreeReq(); }
   KdTreeReq(const KdTreeReq &A) { CopyKdTreeReq(A); }
   ~KdTreeReq() { FreeKdTreeReq(); }
   void MakeKdBuf(const KdTree &Kd);
// private:
   void InitKdTreeReq();
   void CopyKdTreeReq(const KdTreeReq &A);
   void FreeKdTreeReq();
   RVector _X, _Lo, _Hi;
   Integer _InKs; double _InR; bool _Selfie; double _EstF;
   Integer _CurK; ZVector _Ix; RVector _R, _Buf, _CurLo, _CurHi; double _Dist;
};
void KdTreeReq::InitKdTreeReq() {
   _Lo.InitVector(), _Hi.InitVector();
   _X.InitVector(), _Ix.InitVector();
   _R.InitVector(), _Buf.InitVector();
   _CurLo.InitVector(), _CurHi.InitVector();
}
void KdTreeReq::CopyKdTreeReq(const KdTreeReq &A) {
   _Lo.CopyVector(A._Lo), _Hi.CopyVector(A._Hi);
   _InKs = A._InKs, _InR = A._InR, _Selfie = A._Selfie, _EstF = A._EstF, _CurK = A._CurK;
   _X.CopyVector(A._X), _Ix.CopyVector(A._Ix);
   _R.CopyVector(A._R), _Buf.CopyVector(A._Buf);
   _CurLo.CopyVector(A._CurLo), _CurHi.CopyVector(A._CurHi);
   _Dist = A._Dist;
}
void KdTreeReq::FreeKdTreeReq() {
   _Lo.FreeVector(), _Hi.FreeVector();
   _X.FreeVector(), _Ix.FreeVector();
   _R.FreeVector(), _Buf.FreeVector();
   _CurLo.FreeVector(), _CurHi.FreeVector();
}

struct KdTree {
   KdTree() { InitKdTree(); }
   KdTree(const KdTree &A) { CopyKdTree(A); }
   ~KdTree() { FreeKdTree(); }
   void TagKd(const RMatrix XY, const ZVector &Tags, const Integer N, const Integer Xs, const Integer Ys, const Integer Metric);
   void MakeKd(const RMatrix XY, const Integer N, const Integer Xs, const Integer Ys, const Integer Metric);
// private:
   void InitKdTree();
   void CopyKdTree(const KdTree &A);
   void FreeKdTree();
   Integer SplitKd(Integer Lo, Integer Hi, const Integer D, const double S) const;
   void GenerateKd(Integer &NodeN, Integer &SplitN, const Integer Lo, const Integer Hi, const Integer HiLeafN);
   void NewIndependentsKd(const Integer Xs, const Integer Ys);
   void NewDependentsKd(const Integer N, const Integer Xs, const Integer Ys);
   Integer _N, _Xs, _Ys, _Metric;
   RMatrix _XY;
   ZVector _Tags;
   RVector _Lo, _Hi;
   ZVector _Node; RVector _Split;
   KdTreeReq _Kb;
   Integer _DebugN;
};
void KdTree::InitKdTree() {
   _XY.InitMatrix(), _Tags.InitVector();
   _Lo.InitVector(), _Hi.InitVector();
   _Node.InitVector(), _Split.InitVector();
   _Kb.InitKdTreeReq();
}
void KdTree::CopyKdTree(const KdTree &A) {
   _N = A._N, _Xs = A._Xs, _Ys = A._Ys, _Metric = A._Metric;
   _XY.CopyMatrix(A._XY), _Tags.CopyVector(A._Tags);
   _Node.CopyVector(A._Node), _Split.CopyVector(A._Split);
   _Lo.CopyVector(A._Lo), _Hi.CopyVector(A._Hi);
   _Kb.CopyKdTreeReq(A._Kb);
   _DebugN = A._DebugN;
}
void KdTree::FreeKdTree() {
   _XY.FreeMatrix(), _Tags.FreeVector();
   _Node.FreeVector(), _Split.FreeVector();
   _Lo.FreeVector(), _Hi.FreeVector();
   _Kb.FreeKdTreeReq();
}

const Integer SplitNodeSizeKd = 6, FirstVersionKd = 0;

// Segregate nodes [Lo, Hi) in dimension D by threshold S into [Lo, Mid) + [Mid, Hi).
// This subroutine doesn't create tree structures, but only rearranges nodes.
Integer KdTree::SplitKd(Integer Lo, Integer Hi, const Integer D, const double S) const {
   Assert(_N > 0, "KdTree::SplitKd: internal error");
// Split XY/Tags in two parts: [Lo, Hi) is the non-processed part of XY/Tags.
// After the cycle is done, we have [Lo, Hi) = [Mid, Mid].
// We deal with this element separately.
// This establishes the splitting point for the original interval [Lo, Hi).
   Hi--;
   while (Lo < Hi) {
      double *LoV = _XY._XY[Lo];
      if (LoV[D] <= S) Lo++; // XY[Lo] is where it belongs. Advance Lo.
      else { // XY[Lo, ...] must be at Hi. Swap and advance Hi.
         double *HiV = _XY._XY[Hi]; Integer XYs = 2*_Xs + _Ys;
         for (Integer xy = 0; xy < XYs; xy++) Swap(LoV[xy], HiV[xy]);
         Swap(_Tags._X[Lo], _Tags._X[Hi--]);
      }
   }
   if (_XY._XY[Lo][D] <= S) Lo++; else Hi--;
   return Lo;
}

// Recursive KD-tree generation subroutine.
// Parameters:
//	NodeN:		The unused part of Nodes[] to be filled by the tree.
//	SplitN:		The unused part of Splits[]
//	[Lo, Hi)	The interval that points are to be processed from.
// NodesOffs[] and SplitsOffs[] must be large enough.
// AlgLib: Copyright 28.02.2010 by Sergey Bochkanov
void KdTree::GenerateKd(Integer &NodeN, Integer &SplitN, const Integer Lo, const Integer Hi, const Integer HiLeafN) {
   Assert(_N > 0, "KdTree::GenerateKd: internal error (empty tree)");
   Assert(Hi > Lo, "KdTree::GenerateKd: internal error (Hi <= Lo)");
// Generate leaf if needed.
   if (Hi - Lo <= HiLeafN) { _Node._X[NodeN++] = Hi - Lo, _Node._X[NodeN++] = Lo; return; }
// Select which dimension D to split.
// Create a leaf node if the bounding box has size zero.
   Integer D = 0; double Ds = _Kb._CurHi._X[0] - _Kb._CurLo._X[0];
   for (Integer X = 1; X < _Xs; X++) {
      double V = _Kb._CurHi._X[X] - _Kb._CurLo._X[X];
      if (V > Ds) Ds = V, D = X;
   }
   if (Ds == 0.0) { _Node._X[NodeN++] = Hi - Lo, _Node._X[NodeN++] = Lo; return; }
// Select split position S using the sliding midpoint rule, rearrange points into [Lo, Mid) + [Mid, Hi).
// If all the points have the same D component LoV == HiV,
// we zero dimension D of the bounding box and repeat the tree construction.
   double S = _Kb._CurLo._X[D] + 0.5*Ds;
   CopyV(_Kb._Buf._X, 1, &_XY._XY[Lo][D], _XY._Xs, Hi - Lo);
   Integer N = Hi - Lo;
   Integer LoNs = 0, LoIx = Lo; double LoV = _Kb._Buf._X[0];
   Integer HiNs = 0, HiIx = Lo; double HiV = _Kb._Buf._X[0];
   for (Integer n = 0; n < N; n++) {
      double V = _Kb._Buf._X[n];
      if (V < LoV) LoV = V, LoIx = Lo + n; else if (V > HiV) HiV = V, HiIx = Lo + n;
      if (V < S) LoNs++; else if (V > S) HiNs++;
   }
   if (LoV == HiV) {
   // If all the points have the same D component LoV == HiV,
   // we zero dimension D of the bounding box and repeat the tree construction.
      double CurLo = _Kb._CurLo._X[D]; _Kb._CurLo._X[D] = LoV;
      double CurHi = _Kb._CurHi._X[D]; _Kb._CurHi._X[D] = HiV;
      GenerateKd(NodeN, SplitN, Lo, Hi, HiLeafN);
      _Kb._CurLo._X[D] = CurLo, _Kb._CurHi._X[D] = CurHi;
      return;
   }
   Integer Mid;
   if (LoNs > 0 && HiNs > 0) Mid = SplitKd(Lo, Hi, D, S); // Normal midpoint split.
   else if (LoNs == 0) { // Sliding down midpoint.
   // Move split to LoV, place one point to the left bin (move to Lo), others - to the right bin
      S = LoV;
      if (LoIx != Lo) {
         for (Integer xy = 0; xy < 2*_Xs + _Ys; xy++) Swap(_XY._XY[LoIx][xy], _XY._XY[Lo][xy]);
         Swap(_Tags._X[LoIx], _Tags._X[Lo]);
      }
      Mid = Lo + 1;
   } else {
   // Move split to HiV, place one point to the right bin (move to Hi - 1), others - to the left bin
      S = HiV;
      if (HiIx != Hi - 1) {
         for (Integer xy = 0; xy < 2*_Xs + _Ys; xy++) Swap(_XY._XY[HiIx][xy], _XY._XY[Hi - 1][xy]);
         Swap(_Tags._X[HiIx], _Tags._X[Hi - 1]);
      }
      Mid = Hi - 1;
   }
// Generate the 'split' node.
   _Node._X[NodeN] = 0, _Node._X[NodeN + 1] = D, _Node._X[NodeN + 2] = SplitN;
   _Split._X[SplitN] = S;
   Integer ExOff = NodeN; NodeN += SplitNodeSizeKd, SplitN++;
// Update CurBox, do recursive calls on the split subtrees, restore CurBox.
   double V;
   _Node._X[ExOff + 3] = NodeN, V = _Kb._CurHi._X[D], _Kb._CurHi._X[D] = S;
   GenerateKd(NodeN, SplitN, Lo, Mid, HiLeafN), _Kb._CurHi._X[D] = V;
   _Node._X[ExOff + 4] = NodeN, V = _Kb._CurLo._X[D], _Kb._CurLo._X[D] = S;
   GenerateKd(NodeN, SplitN, Mid, Hi, HiLeafN), _Kb._CurLo._X[D] = V;
// Zero-fill unused portions of the node (avoid false warnings by Valgrind about attempt to serialize uninitialized values)
   Assert(SplitNodeSizeKd == 6, "KdTree::GenerateKd: node size has unexpectedly changed");
   _Node._X[ExOff + 5] = 0;
}

// Recursive subroutine for NN queries.
// AlgLib: Copyright 28.02.2010 by Sergey Bochkanov
static void NnQueryKd(const KdTree &Kd, KdTreeReq &Kb, const Integer Off) {
   Assert(Kd._N > 0, "NnQueryKd: internal error");
   if (Kd._Node._X[Off] > 0) { // Leaf node. Process points.
      Integer LoY = Kd._Node._X[Off + 1], HiY = LoY + Kd._Node._X[Off];
      for (Integer y = LoY; y < HiY; y++) {
      // Calculate the distance.
         Integer Xs = Kd._Xs; double Dist = 0.0;
         switch (Kd._Metric) {
            case 0:
               for (Integer x = 0; x < Xs; x++) Dist = Max(Dist, fabs(Kd._XY._XY[y][x] - Kb._X._X[x]));
            break;
            case 1:
               for (Integer x = 0; x < Xs; x++) Dist += fabs(Kd._XY._XY[y][x] - Kb._X._X[x]);
            break;
            case 2:
               for (Integer x = 0; x < Xs; x++) Dist += Squ(Kd._XY._XY[y][x] - Kb._X._X[x]);
            break;
         }
      // Skip points with zero distance if self-matches are turned off.
         if (Dist == 0.0 && !Kb._Selfie) continue;
      // We can NOT process the point if the R-criterion fails to hold, i.e. if InR != 0.0 && Dist > R.
         if (Kb._InR == 0.0 || Dist <= Kb._InR) {
         // The R-criterion is satisfied, we must either:
         // * replace the worst point (or skip, if it is better), if InKs != 0 && CurK == InKs,
         // * add the current point without replacement otherwise.
            if (Kb._CurK < Kb._InKs || Kb._InKs == 0) // Add the current point to heap without replacement.
               TagHeapPush(Kb._R, Kb._Ix, Kb._CurK, Dist, y);
            else if (Dist < Kb._R._X[0])
            // New points are added or not, depending on their distance.
            // If added, they replace element at the top of the heap
               if (Kb._InKs == 1) Kb._Ix._X[0] = y, Kb._R._X[0] = Dist;
               else TagHeapReTop(Kb._R, Kb._Ix, Kb._InKs, Dist, y);
         }
      }
   } else if (Kd._Node._X[Off] == 0) { // Simple split.
   // D = the dimension to split, S = the split position.
      Integer D = Kd._Node._X[Off + 1]; double S = Kd._Split._X[Kd._Node._X[Off + 2]];
   // Determine the range [LoSubOff, HiSubOff] of chances for subboxes.
      Integer LoSubOff, HiSubOff; bool BestLeft = Kb._X._X[D] <= S;
      if (BestLeft)
         HiSubOff = Kd._Node._X[Off + 3], LoSubOff = Kd._Node._X[Off + 4];
      else
         LoSubOff = Kd._Node._X[Off + 3], HiSubOff = Kd._Node._X[Off + 4];
   // Navigate through subboxes.
      for (Integer Sub = 0; Sub < 2; Sub++) {
      // Select which subbox to process:
      // * SubOff = current subbox offset in Nodes[],
      // * IsLo indicates whether the low or high value of the bounding box is to be changed on update.
         Integer SubOff = Sub != 0? LoSubOff: HiSubOff; bool IsLo = (Sub != 0) == BestLeft;
      // Update the bounding box and current distance.
         double Dist, T1, V;
         if (IsLo) {
            Dist = Kb._Dist, T1 = Kb._X._X[D], V = Kb._CurLo._X[D];
            if (T1 <= S) switch (Kd._Metric) {
               case 0: Kb._Dist = Max(Kb._Dist, S - T1); break;
               case 1: Kb._Dist += Min(S - V, S - T1); break;
               case 2: Kb._Dist += Squ(S - T1) - Squ(Max(V - T1, 0.0)); break;
            }
            Kb._CurLo._X[D] = S;
         } else {
            Dist = Kb._Dist, T1 = Kb._X._X[D], V = Kb._CurHi._X[D];
            if (T1 >= S) switch (Kd._Metric) {
               case 0: Kb._Dist = Max(Kb._Dist, T1 - S); break;
               case 1: Kb._Dist += Min(V - S, T1 - S); break;
               case 2: Kb._Dist += Squ(T1 - S) - Squ(Max(T1 - V, 0.0)); break;
            }
            Kb._CurHi._X[D] = S;
         }
      // Decide: to dive into cell or not to dive.
      // CurK < InKs means not all points are found, and
      // CurK == InKs means decide to dive or not to dive using point position relative to bounding box.
         if (
            !(Kb._InR != 0.0 && Kb._Dist > Kb._InR) &&
            (Kb._CurK < Kb._InKs || Kb._InKs == 0 || Kb._Dist <= Kb._R._X[0]*Kb._EstF)
         ) NnQueryKd(Kd, Kb, SubOff);
      // Restore the bounding box and distance.
         if (IsLo) Kb._CurLo._X[D] = V; else Kb._CurHi._X[D] = V;
         Kb._Dist = Dist;
      }
   }
}

// Copy X[] to Kb._X[].
// Load the distance from X[] to the bounding box.
// Initialize Kb._CurBox[].
// AlgLib: Copyright 28.02.2010 by Sergey Bochkanov
static void InitBoxKd(const KdTree &Kd, const RVector &X, KdTreeReq &Kb) {
   Assert(Kd._N > 0, "InitBoxKd: internal error");
// Calculate the distance from the point to the current bounding box.
   Kb._Dist = 0.0;
   switch (Kd._Metric) {
      case 0: for (Integer x = 0; x < Kd._Xs; x++) {
         double Vx = X._X[x], LoV = Kd._Lo._X[x], HiV = Kd._Hi._X[x];
         Kb._X._X[x] = Vx, Kb._CurLo._X[x] = LoV, Kb._CurHi._X[x] = HiV;
         if (Vx < LoV) Kb._Dist = Max(Kb._Dist, LoV - Vx); else if (Vx > HiV) Kb._Dist = Max(Kb._Dist, Vx - HiV);
      }
      break;
      case 1: for (Integer x = 0; x < Kd._Xs; x++) {
         double Vx = X._X[x], LoV = Kd._Lo._X[x], HiV = Kd._Hi._X[x];
         Kb._X._X[x] = Vx, Kb._CurLo._X[x] = LoV, Kb._CurHi._X[x] = HiV;
         if (Vx < LoV) Kb._Dist += LoV - Vx; else if (Vx > HiV) Kb._Dist += Vx - HiV;
      }
      break;
      case 2: for (Integer x = 0; x < Kd._Xs; x++) {
         double Vx = X._X[x], LoV = Kd._Lo._X[x], HiV = Kd._Hi._X[x];
         Kb._X._X[x] = Vx, Kb._CurLo._X[x] = LoV, Kb._CurHi._X[x] = HiV;
         if (Vx < LoV) Kb._Dist += Squ(LoV - Vx); else if (Vx > HiV) Kb._Dist += Squ(Vx - HiV);
      }
      break;
   }
}

// Allocate all dataset-independent array fields of KdTree, i.e. array fields whose dimensions do not depend on the dataset size.
// Do not set Kd._Xs or Kd._Ys; just allocate the arrays.
// AlgLib: Copyright 14.03.2011 by Sergey Bochkanov
void KdTree::NewIndependentsKd(const Integer Xs, const Integer Ys) {
   Assert(_N > 0, "KdTree::NewIndependentsKd: internal error");
   _Lo.InitVector(Xs), _Hi.InitVector(Xs);
}

// Allocate all dataset-dependent array fields of KdTree, i.e. array fields whose dimensions depend on the dataset size.
// Do not set Kd._N, Kd._Xs or Kd._Ys; just allocate the arrays.
// AlgLib: Copyright 14.03.2011 by Sergey Bochkanov
void KdTree::NewDependentsKd(const Integer N, const Integer Xs, const Integer Ys) {
   Assert(N > 0, "KdTree::NewDependentsKd: internal error");
   _XY.InitMatrix(N, 2*Xs + Ys), _Tags.InitVector(N);
   _Node.InitVector(2*N*SplitNodeSizeKd), _Split.InitVector(2*N);
}

// This function checks consistency of request buffer structure with dimensions of kd-tree object.
// AlgLib: Copyright 02.04.2016 by Sergey Bochkanov
static void CheckBufKd(const KdTree &Kd, const KdTreeReq &Kb) {
   Assert(Kb._X._N >= Kd._Xs, "CheckBufKd: dimensions of Kb are inconsistent with Kd");
   Assert(Kb._Ix._N >= Kd._N, "CheckBufKd: dimensions of Kb are inconsistent with Kd");
   Assert(Kb._R._N >= Kd._N, "CheckBufKd: dimensions of Kb are inconsistent with Kd");
   Assert(Kb._Buf._N >= Max(Kd._N, Kd._Xs), "CheckBufKd: dimensions of Kb are inconsistent with Kd");
   Assert(Kb._CurLo._N >= Kd._Xs, "CheckBufKd: dimensions of Kb are inconsistent with Kd");
   Assert(Kb._CurHi._N >= Kd._Xs, "CheckBufKd: dimensions of Kb are inconsistent with Kd");
}

// R-NN query: All points within an R-sphere centered at X, using an external thread-local buffer, in ascending order of distance from X.
// You can call this function from multiple threads for same Kd-Tree, if each thread uses its own buffer object.
// NOTE:
// *	It is also possible to perform unordered queries performed using QueryDistUKdBuf().
//	Such queries are faster because we do not have to use heap structure for sorting.
// Inputs:
//	Kd:	A KD-Tree.
//	Kb:	The request buffer: created for this particular KD-tree structure with MakeKdBuf().
//	X, Xs:	A point, represented as an Xs-vector.
//	R:	The radius of the sphere R > 0 (in the corresponding norm).
//	Selfie:	true (by default) iff self-matches are to be allowed.
// Result
// *	The number of neighbors found; >= 0
// The result are stored internally in the buffer object and can be queried by these routines:
// *	ResultsXsKdBuf()	the X-values.
// *	ResultsXYsKdBuf()	the X- and Y-values
// *	ResultsTagsKdBuf()	the tags.
// *	ResultsDistKdBuf()	the distances.
// IMPORTANT:
// *	The KD-tree buffer should be used only with the KD-tree it was initialized for.
//	Any attempt to use it with different objects may produce an integrity check failure exception
//	because the sizes of the internal arrays may not fit the dimensions of KD-tree.
// AlgLib: Copyright 18.03.2016 by Sergey Bochkanov
static Integer QueryRnnKd(const KdTree &Kd, KdTreeReq &Kb, const RVector &X, const double R, const bool Selfie, const bool Sorted) {
// Handle the special case: Kd._N == 0.
   if (Kd._N == 0) return Kb._CurK = 0;
// Check consistency of request buffer
   CheckBufKd(Kd, Kb);
// Prepare parameters
   Kb._InKs = 0, Kb._InR = Kd._Metric != 2? R: R*R, Kb._Selfie = Selfie, Kb._EstF = 1.0, Kb._CurK = 0;
// Calculate distance from the point to the current bounding box.
   InitBoxKd(Kd, X, Kb);
// Call recursive searchand return the results as a heap.
   NnQueryKd(Kd, Kb, 0);
   Integer Dist = Kb._CurK;
// Pop from heap to sort the representation.
// The last element is left alone, because it is already where it belongs.
   if (Sorted) for (Integer K = Dist, k = Dist; k > 1; k--) TagHeapPop(Kb._R, Kb._Ix, K);
   return Dist;
}

// KD-Tree creation.
// Create a KD-Tree from a set of X-values and optional Y-values (and optionally, integer tags).
// Described in further detail in AlgLibMisc.h.

// A new buffer structure which can be used to perform parallel KD-tree requests.
void KdTreeReq::MakeKdBuf(const KdTree &Kd) {
   _X.InitVector(Kd._Xs), _Ix.InitVector(Kd._N);
   _R.InitVector(Kd._N), _Buf.InitVector(Max(Kd._N, Kd._Xs));
   _Lo.InitVector(Kd._Xs), _Hi.InitVector(Kd._Xs);
   _CurLo.InitVector(Kd._Xs), _CurHi.InitVector(Kd._Xs);
   _CurK = 0;
}

// Create a KD-tree from set of X-values, integer tags and optional Y-values
void KdTree::TagKd(const RMatrix XY, const ZVector &Tags, const Integer N, const Integer Xs, const Integer Ys, const Integer Metric) {
   Assert(N >= 0, "KdTree::TagKd: N < 0");
   Assert(Xs >= 1, "KdTree::TagKd: Xs < 1");
   Assert(Ys >= 0, "KdTree::TagKd: Ys < 0");
   Assert(Metric >= 0 && Metric <= 2, "KdTree::TagKd: incorrect Metric");
   Assert(XY._Ys >= N, "KdTree::TagKd: Ys(X) < N");
   Assert(XY._Xs >= Xs + Ys || N == 0, "KdTree::TagKd: Xs(X) < Xs + Ys");
   Assert(IsFiniteRMatrix(XY, N, Xs + Ys), "KdTree::TagKd: XY contains infinite or NaN values");
// Initialize.
   _N = N, _Xs = Xs, _Ys = Ys, _Metric = Metric, _Kb._CurK = 0;
   if (N == 0) { InitKdTree(); return; } // Nothing to do: quick exit.
// Allocate and initialize.
   NewIndependentsKd(Xs, Ys), NewDependentsKd(N, Xs, Ys);
   _Kb.MakeKdBuf(*this);
// Initial fill.
   for (Integer n = 0; n < N; n++) CopyV(_XY._XY[n], XY._XY[n], Xs), CopyV(&_XY._XY[n][Xs], XY._XY[n], Xs + Ys), _Tags._X[n] = Tags._X[n];
// Determine the bounding box.
   CopyV(_Lo._X, _XY._XY[0], Xs), CopyV(_Hi._X, _XY._XY[0], Xs);
   for (Integer n = 1; n < N; n++) for (Integer x = 0; x < Xs; x++)
      _Lo._X[x] = Min(_Lo._X[x], _XY._XY[n][x]), _Hi._X[x] = Max(_Hi._X[x], _XY._XY[n][x]);
   CopyV(_Kb._CurLo._X, _Lo._X, Xs), CopyV(_Kb._CurHi._X, _Hi._X, Xs);
// Generate the tree.
   Integer NodeN = 0, SplitN = 0; GenerateKd(NodeN, SplitN, 0, N, 8);
   _Node.ExpandVector(NodeN), _Split.ExpandVector(SplitN);
}

// Create a KD-Tree from a set of X-values and optional Y-values.
void KdTree::MakeKd(const RMatrix XY, const Integer N, const Integer Xs, const Integer Ys, const Integer Metric) {
   Assert(N >= 0, "KdTree::MakeKd: N < 0");
   Assert(Xs >= 1, "KdTree::MakeKd: Xs < 1");
   Assert(Ys >= 0, "KdTree::MakeKd: Ys < 0");
   Assert(Metric >= 0 && Metric <= 2, "KdTree::MakeKd: incorrect Metric");
   Assert(XY._Ys >= N, "KdTree::MakeKd: Ys(X) < N");
   Assert(XY._Xs >= Xs + Ys || N == 0, "KdTree::MakeKd: Xs(X) < Xs + Ys");
   Assert(IsFiniteRMatrix(XY, N, Xs + Ys), "KdTree::MakeKd: XY contains infinite or NaN values");
   ZVector Tags(N); for (Integer n = 0; n < N; n++) Tags._X[n] = 0;
   TagKd(XY, Tags, N, Xs, Ys, Metric);
}

// KD-Tree Query: [K-NN] (approximate) K nearest neighbors and [R-NN] all points within an R-sphere centered at X.
// Described in further detail in AlgLibMisc.h.

// R-NN query: All points within an R-sphere centered at X, using an external thread-local buffer, in ascending order of distance from X.
Integer QueryDistKdBuf(const KdTree &Kd, KdTreeReq &Kb, const RVector &X, const double R, const bool Selfie = true) {
   Assert(isfinite(R) && R > 0.0, "QueryDistKdBuf: incorrect R!");
   Assert(X._N >= Kd._Xs, "QueryDistKdBuf: N(X) < Xs!");
   Assert(IsFiniteRVector(X, Kd._Xs), "QueryDistKdBuf: X contains infinite or NaN values!");
   return QueryRnnKd(Kd, Kb, X, R, Selfie, true);
}

// R-NN query: All points within an R-sphere centered at X, in ascending order of distance from X.
Integer QueryDistKd(KdTree &Kd, const RVector &X, const double R, const bool Selfie = true) {
   Assert(R > 0.0, "QueryDistKd: incorrect R!");
   Assert(X._N >= Kd._Xs, "QueryDistKd: N(X) < Xs!");
   Assert(IsFiniteRVector(X, Kd._Xs), "QueryDistKd: X contains infinite or NaN values!");
   return QueryDistKdBuf(Kd, Kd._Kb, X, R, Selfie);
}

// R-NN query: All points within an R-sphere centered at X, using an external thread-local buffer, unsorted (to significantly speed up large queries).
Integer QueryDistUKdBuf(const KdTree &Kd, KdTreeReq &Kb, const RVector &X, const double R, const bool Selfie = true) {
   Assert(isfinite(R) && R > 0.0, "QueryDistUKdBuf: incorrect R!");
   Assert(X._N >= Kd._Xs, "QueryDistUKdBuf: N(X) < Xs!");
   Assert(IsFiniteRVector(X, Kd._Xs), "QueryDistUKdBuf: X contains infinite or NaN values!");
   return QueryRnnKd(Kd, Kb, X, R, Selfie, false);
}

// Component values (X, XY, Tags, R) of the most recent query; with/without the use of external KdTreeReq buffer.
// Described in further detail in AlgLibMisc.h.

// X- and Y-values from the last query associated with the KdTreeReq Kb.
void ResultXYsKdBuf(const KdTree &Kd, const KdTreeReq &Kb, RMatrix &XY) {
   Integer K = Kb._CurK; if (K == 0) return;
   Integer Xs = Kd._Xs, Ys = Kd._Ys, *IxP = Kb._Ix._X; double **KdXY = Kd._XY._XY;
   XY.GrowMatrix(K, Xs + Ys); double **XP = XY._XY;
   for (Integer k = 0; k < K; k++) CopyV(XP[k], &KdXY[IxP[k]][Xs], Xs + Ys);
}
// Tags from last query associated with the KdTreeReq Kb.
void ResultTagsKdBuf(const KdTree &Kd, const KdTreeReq &Kb, ZVector &Tags) {
   Integer K = Kb._CurK; if (K == 0) return;
   Integer *IxP = Kb._Ix._X, *KdT = Kd._Tags._X;
   Tags.GrowVector(K); Integer *TP = Tags._X;
   for (Integer k = 0; k < K; k++) TP[k] = KdT[IxP[k]];
}
// Tags from the last query.
void ResultTagsKd(const KdTree &Kd, ZVector &Tags) { ResultTagsKdBuf(Kd, Kd._Kb, Tags); }
// Distances from last query associated with the KdTreeReq Kb.
void ResultDistKdBuf(const KdTree &Kd, const KdTreeReq &Kb, RVector &R) {
   Integer K = Kb._CurK; if (K == 0) return;
   double *KbR = Kb._R._X;
   R.GrowVector(K); double *RP = R._X;
// Unload the norms.
// Absolute values are used to handle the cases of negative norms, which KFN requests may produce.
   switch (Kd._Metric) {
      case 0:
         for (Integer k = 0; k < K; k++) RP[k] = fabs(KbR[k]);
      break;
      case 1:
         for (Integer k = 0; k < K; k++) RP[k] = fabs(KbR[k]);
      break;
      case 2:
         for (Integer k = 0; k < K; k++) RP[k] = sqrt(fabs(KbR[k]));
      break;
   }
}
// Distances from the last query.
void ResultDistKd(const KdTree &Kd, RVector &R) { ResultDistKdBuf(Kd, Kd._Kb, R); }

// === IdwInt Package ===
struct IdwQ;

struct IdwBuf {
   IdwBuf() { InitIdwBuf(); }
   IdwBuf(const IdwBuf &A) { CopyIdwBuf(A); }
   ~IdwBuf() { FreeIdwBuf(); }
   IdwBuf(const IdwQ &Q) { IdwMakeBuf(Q); }
// private:
   void InitIdwBuf();
   void CopyIdwBuf(const IdwBuf &A);
   void FreeIdwBuf();
   void IdwMakeBuf(const IdwQ &Q);
   RVector _X, _Y, _YWb, _Wb, _Rb;
   RMatrix _XYb; KdTreeReq _Kb;
};
void IdwBuf::InitIdwBuf() {
   _X.InitVector(), _Y.InitVector(), _YWb.InitVector(), _Wb.InitVector(), _XYb.InitMatrix();
   _Rb.InitVector(); _Kb.InitKdTreeReq();
}
void IdwBuf::CopyIdwBuf(const IdwBuf &A) {
   _X.CopyVector(A._X), _Y.CopyVector(A._Y);
   _YWb.CopyVector(A._YWb), _Wb.CopyVector(A._Wb), _XYb.CopyMatrix(A._XYb);
   _Rb.CopyVector(A._Rb), _Kb.CopyKdTreeReq(A._Kb);
}
void IdwBuf::FreeIdwBuf() {
   _X.FreeVector(), _Y.FreeVector(), _YWb.FreeVector(), _Wb.FreeVector(), _XYb.FreeMatrix();
   _Rb.FreeVector(), _Kb.FreeKdTreeReq();
}

struct IdwQ {
   IdwQ() { InitIdwQ(); }
   IdwQ(const IdwQ &A) { CopyIdwQ(A); }
   ~IdwQ() { FreeIdwQ(); }
   void IdwCalcTs(IdwBuf &Buf, const RVector &X, RVector &Y) const;
   double IdwCalc1(const double X);
   double IdwCalc2(const double X, const double Y);
   double IdwCalc3(const double X, const double Y, const double Z);
   void IdwCalc(const RVector &X, RVector &Y);
   void IdwCalcBuf(const RVector &X, RVector &Y);
// private:
   void InitIdwQ();
   void CopyIdwQ(const IdwQ &A);
   void FreeIdwQ();
   Integer _Xs, _Ys, _Op, _Ls, _Ps;
   double _R0, _dR, _Lam0, _Lam1, _dLam, _ShepP;
   RVector _GlobP, _ShepXY;
   KdTree _Kd; IdwBuf _Buf;
};
void IdwQ::InitIdwQ() {
   _GlobP.InitVector(), _ShepXY.InitVector();
   _Kd.InitKdTree(), _Buf.InitIdwBuf();
}
void IdwQ::CopyIdwQ(const IdwQ &A) {
   _Xs = A._Xs, _Ys = A._Ys, _Op = A._Op, _Ls = A._Ls, _Ps = A._Ps;
   _R0 = A._R0, _dR = A._dR, _Lam0 = A._Lam0, _Lam1 = A._Lam1, _dLam = A._dLam, _ShepP = A._ShepP;
   _GlobP.CopyVector(A._GlobP), _ShepXY.CopyVector(A._ShepXY);
   _Kd.CopyKdTree(A._Kd), _Buf.CopyIdwBuf(A._Buf);
}
void IdwQ::FreeIdwQ() {
   _GlobP.FreeVector(), _ShepXY.FreeVector();
   _Kd.FreeKdTree(), _Buf.FreeIdwBuf();
}

struct IdwB {
   IdwB() { InitIdwB(); }
   IdwB(const IdwB &A) { CopyIdwB(A); }
   ~IdwB() { FreeIdwB(); }
   IdwB(const Integer Xs, const Integer Ys);
   void IdbSetLayers(const Integer Ls);
   void IdbSetPoints(const RMatrix &XY, const Integer Ps);
   void IdbSetMStab(const double RadS);
   void IdbSetShepard(const double P);
   void IdbSetModShepard(const double RadS);
   void IdbSetUserTerm(const double V);
   void IdbSetConstTerm();
   void IdbSetZeroTerm();
// private:
   void InitIdwB();
   void CopyIdwB(const IdwB &A);
   void FreeIdwB();
   Integer _PriorType; RVector _PriorVal;
   Integer _Op, _Ls, _Ps, _Xs, _Ys;
   double _R0, _dR, _Lam0, _Lam1, _dLam, _ShepP;
   RVector _XY;
   RMatrix _TmpXY, _TmpLs;
   RVector _TmpR, _TmpX, _TmpWY, _TmpW, _TmpMu;
   ZVector _TmpTags;
   KdTree _TmpKd;
};
void IdwB::InitIdwB() {
   _PriorVal.InitVector(), _XY.InitVector();
   _TmpXY.InitMatrix(), _TmpLs.InitMatrix(), _TmpTags.InitVector();
   _TmpR.InitVector(), _TmpX.InitVector();
   _TmpWY.InitVector(), _TmpW.InitVector();
   _TmpKd.InitKdTree(), _TmpMu.InitVector();
}
void IdwB::CopyIdwB(const IdwB &A) {
   _PriorType = A._PriorType, _XY.CopyVector(A._XY);
   _PriorVal.CopyVector(A._PriorVal);
   _Op = A._Op, _Ls = A._Ls, _Ps = A._Ps, _Xs = A._Xs, _Ys = A._Ys;
   _R0 = A._R0, _dR = A._dR, _Lam0 = A._Lam0, _Lam1 = A._Lam1, _dLam = A._dLam, _ShepP = A._ShepP;
   _TmpXY.CopyMatrix(A._TmpXY), _TmpLs.CopyMatrix(A._TmpLs), _TmpTags.CopyVector(A._TmpTags);
   _TmpR.CopyVector(A._TmpR), _TmpX.CopyVector(A._TmpX);
   _TmpWY.CopyVector(A._TmpWY), _TmpW.CopyVector(A._TmpW);
   _TmpKd.CopyKdTree(A._TmpKd), _TmpMu.CopyVector(A._TmpMu);
}
void IdwB::FreeIdwB() {
   _PriorVal.FreeVector(), _XY.FreeVector();
   _TmpXY.FreeMatrix(), _TmpLs.FreeMatrix(), _TmpTags.FreeVector();
   _TmpR.FreeVector(), _TmpX.FreeVector();
   _TmpWY.FreeVector(), _TmpW.FreeVector();
   _TmpKd.FreeKdTree(), _TmpMu.FreeVector();
}

struct IdwR {
   IdwR() { InitIdwR(); }
   IdwR(const IdwR &A) { CopyIdwR(A); }
   ~IdwR() { FreeIdwR(); }
   IdwR(IdwB &Qb);
// public:
   double _RmsErr, _AveErr, _MaxErr, _R2;
// private:
   void InitIdwR();
   void CopyIdwR(const IdwR &A);
   void FreeIdwR();
   void IdwErrorMetricsViaCalc(const IdwB &Qb);
   IdwQ _Q;
};
void IdwR::InitIdwR() {
   _Q.InitIdwQ();
}
void IdwR::CopyIdwR(const IdwR &A) {
   _Q.CopyIdwQ(A._Q);
   _RmsErr = A._RmsErr, _AveErr = A._AveErr, _MaxErr = A._MaxErr, _R2 = A._R2;
}
void IdwR::FreeIdwR() {
   _Q.FreeIdwQ();
}

static double IdwW0 = 1.0, IdwLam0 = 0.3333;

// A new buffer which can be used for parallel IDW model evaluations from multiple threads, with each thread accessing the IDW model from its own buffer.
void IdwBuf::IdwMakeBuf(const IdwQ &Q) {
   Assert(Q._Xs >= 1, "IdwBuf::IdwMakeBuf: integrity check failed");
   Assert(Q._Ys >= 1, "IdwBuf::IdwMakeBuf: integrity check failed");
   Assert(Q._Ls >= 0, "IdwBuf::IdwMakeBuf: integrity check failed");
   Assert(Q._Op >= 0, "IdwBuf::IdwMakeBuf: integrity check failed");
   _X.InitVector(Q._Xs), _Y.InitVector(Q._Ys);
   _YWb.InitVector(Q._Ys*Max(Q._Ls, 1)), _Wb.InitVector(Max(Q._Ls, 1)), _XYb.InitMatrix();
   _Rb.InitVector();
   if (Q._Ls >= 1 && Q._Op != 0) _Kb.MakeKdBuf(Q._Kd); else _Kb.InitKdTreeReq();
}

// A new builder used to generate IDW model from irregularly sampled (scattered) dataset.
IdwB::IdwB(const Integer Xs, const Integer Ys) {
   const Integer IdwLs = 16;
   InitIdwB();
   Assert(Xs >= 1, "IdwB::IdwB: Xs <= 0");
   Assert(Ys >= 1, "IdwB::IdwB: Ys <= 0");
// We choose reasonable defaults for the algorithm:
// *	MSTAB algorithm
// *	12 layers
// *	default radius
// *	default Lambda0
   _Op = 2, _PriorType = 2;
   _PriorVal.GrowVector(Ys);
   _Ls = IdwLs, _R0 = 0.0, _dR = 0.5;
   _Lam0 = IdwLam0, _Lam1 = 0.0, _dLam = 1.0;
// Other parameters, not used but initialized.
   _ShepP = 0.0;
// Initial dataset is empty.
   _Ps = 0, _Xs = Xs, _Ys = Ys;
}

// Change the number of layers used by the IDW-MSTAB algorithm.
void IdwB::IdbSetLayers(const Integer Ls) {
   Assert(Ls >= 1, "IdwB::IdbSetLayers: Ls < 1");
   _Ls = Ls;
}

// Add a data set to the Idw builder, ovewriting any data set already there.
void IdwB::IdbSetPoints(const RMatrix &XY, const Integer Ps) {
   Assert(Ps >= 0, "IdwB::IdbSetPoints: Ps < 0");
   Assert(XY._Ys >= Ps, "IdwB::IdbSetPoints: Ys(XY) < Ps");
   Integer XYs = _Xs + _Ys;
   Assert(Ps == 0 || XY._Xs >= XYs, "IdwB::IdbSetPoints: Xs(XY) < Xs + Ys");
   Assert(IsFiniteRMatrix(XY, Ps, XYs), "IdwB::IdbSetPoints: XY contains infinite or NaN values!");
   _Ps = Ps;
   _XY.GrowVector(Ps*XYs);
   for (Integer xyp = 0, p = 0; p < Ps; p++) for (Integer xy = 0; xy < XYs; xyp++, xy++) _XY._X[xyp] = XY._XY[p][xy];
}

// Configure the IDW model to use, as its construction algorithm, the Multilayer Stabilized IDW method (IDW-MSTAB).
void IdwB::IdbSetMStab(const double RadS) {
   Assert(isfinite(RadS), "IdwB::IdbSetMStab: RadS is not finite");
   Assert(RadS > 0.0, "IdwB::IdbSetMStab: RadS <= 0");
// Set the algorithm.
   _Op = 2;
// Set the options.
   _R0 = RadS, _dR = 0.5;
   _Lam0 = IdwLam0, _Lam1 = 0.0, _dLam = 1.0;
}

// Configure the IDW model to use, as its construction algorithm, the textbook Shepard's algorithm with custom (user-specified) power parameter.
void IdwB::IdbSetShepard(const double P) {
   Assert(isfinite(P), "IdwB::IdbSetShepard: P is not finite");
   Assert(P > 0.0, "IdwB::IdbSetShepard: P <= 0");
// Set the algorithm and options.
   _Op = 0, _ShepP = P;
}

// Configure the IDW model to use, as its construction algorithm, the 'textbook' modified Shepard's algorithm with user-specified search radius.
void IdwB::IdbSetModShepard(const double RadS) {
   Assert(isfinite(RadS), "IdwB::IdbSetModShephard: RadS is not finite");
   Assert(RadS > 0.0, "IdwB::IdbSetModShephard: RadS <= 0");
// Set the algorithm and options.
   _Op = 1, _R0 = RadS;
}

// Set the prior term (the model value at infinity) as user-specified value.
void IdwB::IdbSetUserTerm(const double V) {
   Assert(isfinite(V), "IdwB::IdbSetUserTerm: infinite/NAN value passed");
   _PriorType = 0;
   for (Integer Y = 0; Y < _Ys; Y++) _PriorVal._X[Y] = V;
}

// Set the constant prior term (model value at infinity) to the mean value over dataset.
void IdwB::IdbSetConstTerm() { _PriorType = 2; }

// Set the zero prior term (model value at infinity).
void IdwB::IdbSetZeroTerm() { _PriorType = 3; }

// The values of the IDW model at the given point, using an external buffer object (the internal temporaries of IDW model are not modified).
void IdwQ::IdwCalcTs(IdwBuf &Buf, const RVector &X, RVector &Y) const {
   Assert(X._N >= _Xs, "IdwQ::IdwCalcTs: N(X) < _Xs");
   Assert(IsFiniteRVector(X, _Xs), "IdwQ::IdwCalcTs: X contains infinite or NaN values");
// Allocate the output.
   Y.GrowVector(_Ys);
// Quick exit for _Ls == 0 (no data set).
   if (_Ls == 0) {
      for (Integer y = 0; y < _Ys; y++) Y._X[y] = _GlobP._X[y];
      return;
   }
   if (_Op == 0) { // Textbook Shepard's method.
      Assert(_Ps > 0, "IdwQ::IdwCalcTs: integrity check failed");
      double Eps = 1.0E-50;
      Integer XYs = _Xs + _Ys;
      double ShepP = _ShepP;
      for (Integer y = 0; y < _Ys; y++) Y._X[y] = 0.0, Buf._YWb._X[y] = Eps;
      for (Integer p = 0; p < _Ps; p++) {
      // Compute squared distance.
         double V = 0.0;
         for (Integer x = 0; x < _Xs; x++) V += Squ(_ShepXY._X[p*XYs + x] - X._X[x]);
      // Compute weight (with small regularizing addition)
         V = 1.0/(Eps + pow(V, ShepP*0.5));
      // Accumulate
         for (Integer y = 0; y < _Ys; y++) Y._X[y] += V*_ShepXY._X[p*XYs + _Xs + y], Buf._YWb._X[y] += V;
      }
      for (Integer y = 0; y < _Ys; y++) Y._X[y] = Y._X[y]/Buf._YWb._X[y] + _GlobP._X[y];
   } else if (_Op == 1) { // Textbook modified Shepard's method.
      double Eps = 1.0E-50, R = _R0;
      for (Integer y = 0; y < _Ys; y++) Y._X[y] = 0.0, Buf._YWb._X[y] = Eps;
      Integer K = QueryDistKdBuf(_Kd, Buf._Kb, X, R, true);
      ResultXYsKdBuf(_Kd, Buf._Kb, Buf._XYb);
      ResultDistKdBuf(_Kd, Buf._Kb, Buf._Rb);
      for (Integer k = 0; k < K; k++) {
         double V = Buf._Rb._X[k];
         V = Squ((R - V)/(R*V + Eps));
         for (Integer y = 0; y < _Ys; y++) Y._X[y] += V*Buf._XYb._XY[k][_Xs + y], Buf._YWb._X[y] += V;
      }
      for (Integer y = 0; y < _Ys; y++) Y._X[y] = Y._X[y]/Buf._YWb._X[y] + _GlobP._X[y];
   } else if (_Op == 2) { // MSTAB.
      Assert(IdwW0 == 1.0, "IdwQ::IdwCalcTs: unexpected W0, integrity check failed");
      double UDecay = 1.0/_dR, U = 1.0/_R0, LamDecay = _dLam;
      bool SpeedUp = _Ys == 1 && _Ls >= 3 && LamDecay == 1.0;
      double Wf0 = 0.0, Ws0 = SpeedUp? IdwW0: 0.0;
      double Wf1 = 0.0, Ws1 = SpeedUp? IdwW0: 0.0;
      if (SpeedUp) // Important special case, _Ys == 1, no lambda-decay, we can perform optimized fast evaluation
         for (Integer L = 0; L < _Ls; L++) Buf._YWb._X[L] = 0.0, Buf._Wb._X[L] = IdwW0;
      else {
      // Setup variables for generic evaluation path
         for (Integer yl = 0; yl < _Ys*_Ls; yl++) Buf._YWb._X[yl] = 0.0;
         for (Integer L = 0; L < _Ls; L++) Buf._Wb._X[L] = IdwW0;
      }
      Integer K = QueryDistUKdBuf(_Kd, Buf._Kb, X, _R0, true);
      ResultXYsKdBuf(_Kd, Buf._Kb, Buf._XYb);
      ResultDistKdBuf(_Kd, Buf._Kb, Buf._Rb);
      for (Integer k = 0; k < K; k++) {
         double Lam = _Lam0, V0 = Buf._Rb._X[k]*U;
         if (SpeedUp) { // An important special case, fast evaluation possible.
            double V = V0*V0; V = (1.0 - V)*(1.0 - V)/(V + Lam);
            Wf0 += V*Buf._XYb._XY[k][_Xs], Ws0 += V;
            V0 *= UDecay; if (V0 >= 1.0) continue;
            V = V0*V0, V = (1.0 - V)*(1.0 - V)/(V + Lam);
            Wf1 += V*Buf._XYb._XY[k][_Xs + 1], Ws1 += V;
            V0 *= UDecay; if (V0 >= 1.0) continue;
            for (Integer L = 2; L < _Ls; L++) {
               if (L == _Ls - 1) Lam = _Lam1;
               double V = V0*V0; V = (1.0 - V)*(1.0 - V)/(V + Lam);
               Buf._YWb._X[L] += V*Buf._XYb._XY[k][_Xs + L], Buf._Wb._X[L] += V;
               V0 *= UDecay; if (V0 >= 1.0) break;
            }
         } else { // The general case.
            for (Integer L = 0; L < _Ls; L++) {
               if (L == _Ls - 1) Lam = _Lam1;
               if (V0 >= 1.0) break;
               double V = V0*V0; V = (1.0 - V)*(1.0 - V)/(V + Lam);
               for (Integer y = 0; y < _Ys; y++) Buf._YWb._X[L*_Ys + y] += V*Buf._XYb._XY[k][_Xs + L*_Ys + y];
               Buf._Wb._X[L] += V, Lam *= LamDecay, V0 *= UDecay;
            }
         }
      }
      if (SpeedUp) // Important special case, finish up the evaluations.
         Buf._YWb._X[0] = Wf0, Buf._Wb._X[0] = Ws0,
         Buf._YWb._X[1] = Wf1, Buf._Wb._X[1] = Ws1;
      for (Integer y = 0; y < _Ys; y++) Y._X[y] = _GlobP._X[y];
      for (Integer L = 0; L < _Ls; L++) for (Integer y = 0; y < _Ys; y++) Y._X[y] += Buf._YWb._X[L*_Ys + y]/Buf._Wb._X[L];
   } else Assert(false, "IdwQ::IdwCalcTs: unexpected Op");
}

// IDW interpolation: scalar target, 1-dimensional argument.
double IdwQ::IdwCalc1(const double X) {
   Assert(_Xs == 1, "IdwQ::IdwCalc1: _Xs != 1");
   Assert(_Ys == 1, "IdwQ::IdwCalc1: _Ys != 1");
   Assert(isfinite(X), "IdwQ::IdwCalc1: X is INFINITY or NAN");
   _Buf._X._X[0] = X, IdwCalcTs(_Buf, _Buf._X, _Buf._Y);
   return _Buf._Y._X[0];
}

// IDW interpolation: scalar target, 2-dimensional argument.
double IdwQ::IdwCalc2(const double X, const double Y) {
   Assert(_Xs == 2, "IdwQ::IdwCalc2: _Xs != 2");
   Assert(_Ys == 1, "IdwQ::IdwCalc2: _Ys != 1");
   Assert(isfinite(X), "IdwQ::IdwCalc2: X is INFINITY or NAN");
   Assert(isfinite(Y), "IdwQ::IdwCalc2: Y is INFINITY or NAN");
   _Buf._X._X[0] = X, _Buf._X._X[1] = Y, IdwCalcTs(_Buf, _Buf._X, _Buf._Y);
   return _Buf._Y._X[0];
}

// IDW interpolation: scalar target, 3-dimensional argument.
double IdwQ::IdwCalc3(const double X, const double Y, const double Z) {
   Assert(_Xs == 3, "IdwQ::IdwCalc3: _Xs != 3");
   Assert(_Ys == 1, "IdwQ::IdwCalc3: _Ys != 1");
   Assert(isfinite(X), "IdwQ::IdwCalc3: X is INFINITY or NAN");
   Assert(isfinite(Y), "IdwQ::IdwCalc3: Y is INFINITY or NAN");
   Assert(isfinite(Z), "IdwQ::IdwCalc3: Z is INFINITY or NAN");
   _Buf._X._X[0] = X, _Buf._X._X[1] = Y, _Buf._X._X[2] = Z, IdwCalcTs(_Buf, _Buf._X, _Buf._Y);
   return _Buf._Y._X[0];
}

// The values of the IDW model at the given point.
void IdwQ::IdwCalc(const RVector &X, RVector &Y) { Y.FreeVector(), IdwCalcTs(_Buf, X, Y); }

// The values of the IDW model at the given point, with a pre-allocated output buffer.
void IdwQ::IdwCalcBuf(const RVector &X, RVector &Y) { IdwCalcTs(_Buf, X, Y); }

// Evaluate error metrics for the model using IDbCalcBuf() to calculate model at each point.
// NOTE:
// *	The current IDW algorithms (MSTAB, MSMOOTH) can generate residuals during model construction, so they do not need this function in order to evaluate error metrics.
// The following fields are filled:
// *	_RmsErr
// *	_AveErr
// *	_MaxErr
// *	_R2
// AlgLib: Copyright 22.10.2018 by Sergey Bochkanov
void IdwR::IdwErrorMetricsViaCalc(const IdwB &Qb) {
   Integer Ps = Qb._Ps, Xs = Qb._Xs, Ys = Qb._Ys;
   if (Ps == 0) { _MaxErr = _AveErr = _RmsErr = 0.0, _R2 = 1.0; return; }
   _MaxErr = _AveErr = _RmsErr = 0.0;
   double Rss = 0.0, Tss = 0.0;
   for (Integer p = 0; p < Ps; p++) {
      for (Integer x = 0; x < Xs; x++) _Q._Buf._X._X[x] = Qb._XY._X[p*(Xs + Ys) + x];
      _Q.IdwCalcTs(_Q._Buf, _Q._Buf._X, _Q._Buf._Y);
      for (Integer y = 0; y < Ys; y++) {
         double V = Qb._XY._X[p*(Xs + Ys) + Xs + y], dV = fabs(V - _Q._Buf._Y._X[y]);
         _RmsErr += dV*dV, _AveErr += dV, _MaxErr = Max(_MaxErr, dV);
         Rss += dV*dV, Tss += Squ(V - Qb._TmpMu._X[y]);
      }
   }
   _RmsErr = sqrt(_RmsErr/(Ps*Ys)), _AveErr /= Ps*Ys, _R2 = 1.0 - Rss/(Tss == 0.0? 1.0: Tss);
}

// Fit the IDW model to the data set configured with the current IDW construction algorithm to produce the model being built and fitting report.
IdwR::IdwR(IdwB &Qb) {
   const double IdwEpsM = 1.0E-50;
   InitIdwR();
   Integer Xs = Qb._Xs, Ys = Qb._Ys, Ps = Qb._Ps;
// Clear the report fields.
   _MaxErr = _AveErr = _RmsErr = 0.0, _R2 = 1.0;
// Quick exit for empty dataset.
   if (Qb._Ps == 0) {
      _Q._Xs = Xs, _Q._Ys = Ys;
      _Q._GlobP.ReSizeVector(Ys); for (Integer y = 0; y < Ys; y++) _Q._GlobP._X[y] = 0.0;
      _Q._Op = 0, _Q._Ls = 0;
      _Q._R0 = 1.0, _Q._dR = 0.5;
      _Q._Lam0 = 0.0, _Q._Lam1 = 0.0, _Q._dLam = 1.0;
      _Q._ShepP = 2.0, _Q._Ps = 0;
      _Q._Buf.IdwMakeBuf(_Q);
      return;
   }
// Compute temporaries which will be required later:
// *	global mean
   Assert(Qb._Ps > 0, "IdwR::IdwR: integrity check failed");
   Qb._TmpMu.GrowVector(Ys);
   for (Integer y = 0; y < Ys; y++) Qb._TmpMu._X[y] = 0.0;
   for (Integer p = 0; p < Ps; p++) for (Integer y = 0; y < Ys; y++) Qb._TmpMu._X[y] += Qb._XY._X[p*(Xs + Ys) + Xs + y];
   for (Integer y = 0; y < Ys; y++) Qb._TmpMu._X[y] /= Ps;
// Compute the global prior.
// NOTE:
// *	For the original Shepard's method it is always the mean value.
   _Q._GlobP.GrowVector(Ys); for (Integer y = 0; y < Ys; y++) _Q._GlobP._X[y] = Qb._TmpMu._X[y];
   if (Qb._Op != 0) {
   // The algorithm is set to one of the "advanced" versions with search radius which can handle non-mean prior term.
      if (Qb._PriorType == 0) // User-specified prior.
         for (Integer y = 0; y < Ys; y++) _Q._GlobP._X[y] = Qb._PriorVal._X[y];
      else if (Qb._PriorType == 3) // Zero prior.
         for (Integer y = 0; y < Ys; y++) _Q._GlobP._X[y] = 0.0;
   }
   if (Qb._Op == 0) { // Textbook Shepard.
   // Initialize the model.
      _Q._Op = 0;
      _Q._Xs = Xs, _Q._Ys = Ys, _Q._Ls = 1;
      _Q._R0 = 1.0, _Q._dR = 0.5;
      _Q._Lam0 = 0.0, _Q._Lam1 = 0.0, _Q._dLam = 1.0;
      _Q._ShepP = Qb._ShepP;
   // Copy the data set.
      _Q._ShepXY.GrowVector(Ps*(Xs + Ys));
      for (Integer p = 0; p < Ps; p++) {
         for (Integer x = 0; x < Xs; x++) _Q._ShepXY._X[p*(Xs + Ys) + x] = Qb._XY._X[p*(Xs + Ys) + x];
         for (Integer y = 0; y < Ys; y++) _Q._ShepXY._X[p*(Xs + Ys) + Xs + y] = Qb._XY._X[p*(Xs + Ys) + Xs + y] - _Q._GlobP._X[y];
      }
      _Q._Ps = Ps;
   // Prepare internal buffer
   // Evaluate report fields
      _Q._Buf.IdwMakeBuf(_Q), IdwErrorMetricsViaCalc(Qb);
   } else if (Qb._Op == 1) { // Textbook modified Shepard's method.
   // Initialize the model.
      _Q._Op = 1;
      _Q._Xs = Xs, _Q._Ys = Ys, _Q._Ls = 1;
      _Q._R0 = Qb._R0, _Q._dR = 1.0;
      _Q._Lam0 = 0.0, _Q._Lam1 = 0.0, _Q._dLam = 1.0;
      _Q._ShepP = 0.0;
   // Build a KD search tree.
      Qb._TmpXY.GrowMatrix(Ps, Xs + Ys);
      for (Integer p = 0; p < Ps; p++) {
         for (Integer x = 0; x < Xs; x++) Qb._TmpXY._XY[p][x] = Qb._XY._X[p*(Xs + Ys) + x];
         for (Integer y = 0; y < Ys; y++) Qb._TmpXY._XY[p][Xs + y] = Qb._XY._X[p*(Xs + Ys) + Xs + y] - _Q._GlobP._X[y];
      }
      _Q._Kd.MakeKd(Qb._TmpXY, Ps, Xs, Ys, 2);
   // Prepare the internal buffer.
   // Evaluate the report fields.
      _Q._Buf.IdwMakeBuf(_Q),
      IdwErrorMetricsViaCalc(Qb);
   } else if (Qb._Op == 2) { // MSTAB algorithm.
      Assert(Qb._Ls >= 1, "IdwR::IdwR: integrity check failed");
   // Initialize the model.
      _Q._Op = 2;
      _Q._Xs = Xs, _Q._Ys = Ys, _Q._Ls = Qb._Ls;
      _Q._R0 = Qb._R0, _Q._dR = 0.5;
      _Q._Lam0 = Qb._Lam0, _Q._dLam = 1.0, _Q._Lam1 = IdwEpsM;
      _Q._ShepP = 0.0;
   // Build a KD search tree, prepare input residuals for the first layer of the model.
      Qb._TmpXY.GrowMatrix(Ps, Xs), Qb._TmpLs.GrowMatrix(Ps, Xs + Ys*(Qb._Ls + 1)), Qb._TmpTags.GrowVector(Ps);
      for (Integer p = 0; p < Ps; p++) {
         for (Integer x = 0; x < Xs; x++) Qb._TmpLs._XY[p][x] = Qb._TmpXY._XY[p][x] = Qb._XY._X[p*(Xs + Ys) + x];
         Qb._TmpTags._X[p] = p;
         for (Integer y = 0; y < Ys; y++) Qb._TmpLs._XY[p][Xs + y] = Qb._XY._X[p*(Xs + Ys) + Xs + y] - _Q._GlobP._X[y];
      }
      Qb._TmpKd.TagKd(Qb._TmpXY, Qb._TmpTags, Ps, Xs, 0, 2);
   // Iteratively build layer by layer.
      Qb._TmpX.GrowVector(Xs), Qb._TmpWY.GrowVector(Ys), Qb._TmpW.GrowVector(Ys);
      for (Integer L = 0; L < Qb._Ls; L++) {
      // Determine the layer metrics.
         double CurR = _Q._R0*pow(_Q._dR, L), Lam = _Q._Lam0*pow(_Q._dLam, L);
         if (L == Qb._Ls - 1) Lam = _Q._Lam1;
      // For each point, compute the residual from fitting with the current layer.
         for (Integer p = 0; p < Ps; p++) {
            for (Integer x = 0; x < Xs; x++) Qb._TmpX._X[x] = Qb._TmpLs._XY[p][x];
            Integer K = QueryDistKd(Qb._TmpKd, Qb._TmpX, CurR, true);
            ResultTagsKd(Qb._TmpKd, Qb._TmpTags);
            ResultDistKd(Qb._TmpKd, Qb._TmpR);
            for (Integer y = 0; y < Ys; y++) Qb._TmpWY._X[y] = 0.0, Qb._TmpW._X[y] = IdwW0;
            for (Integer k = 0; k < K; k++) {
               double V = Squ(Qb._TmpR._X[k]/CurR); V = (1.0 - V)*(1.0 - V)/(V + Lam);
               Integer AdIx = Qb._TmpTags._X[k];
               for (Integer y = 0; y < Ys; y++) Qb._TmpWY._X[y] += V*Qb._TmpLs._XY[AdIx][Xs + L*Ys + y], Qb._TmpW._X[y] += V;
            }
            for (Integer y = 0; y < Ys; y++)
               Qb._TmpLs._XY[p][Xs + (L + 1)*Ys + y] = Qb._TmpLs._XY[p][Xs + L*Ys + y] - Qb._TmpWY._X[y]/Qb._TmpW._X[y];
         }
      }
      _Q._Kd.MakeKd(Qb._TmpLs, Ps, Xs, Ys*Qb._Ls, 2);
   // Evaluate the report fields.
      _MaxErr = _AveErr = _RmsErr = 0.0;
      double Rss = 0.0, Tss = 0.0;
      for (Integer p = 0; p < Ps; p++) for (Integer y = 0; y < Ys; y++) {
         double V = fabs(Qb._TmpLs._XY[p][Xs + Qb._Ls*Ys + y]);
         _RmsErr += V*V, _AveErr += V, _MaxErr = Max(_MaxErr, fabs(V));
         Rss += V*V, Tss += Squ(Qb._XY._X[p*(Xs + Ys) + Xs + y] - Qb._TmpMu._X[y]);
      }
      _RmsErr = sqrt(_RmsErr/(Ps*Ys)), _AveErr /= Ps*Ys, _R2 = 1.0 - Rss/(Tss == 0.0? 1.0: Tss);
   // Prepare the internal buffer.
      _Q._Buf.IdwMakeBuf(_Q);
   } else Assert(false, "IdwR::IdwR: integrity check failed, unexpected algorithm"); // Unknown algorithm.
}

// === Random Package ===
unsigned CurSeed = 1;
void SRand(unsigned Seed) { srand(CurSeed = Seed); }
Integer RandomZ(Integer HiZ) { return ((Integer)rand_r(&CurSeed))%HiZ; }

struct RandomQ {
   RandomQ();
   RandomQ(RandomQ &A): _S0(A._S0), _S1(A._S1), _MagicV(A._MagicV) { }
   ~RandomQ() { }
   double RandomUniformR();
   double RandMidUniformR();
   Integer RandomUniformZ(Integer N);
   double RandomNormal1();
private:
   Integer RandBase();
   Integer _S0, _S1, _MagicV;
};

const Integer RandMax = 0x7fffffa9, RandM0 = 0x7fffffab, RandM1 = 0x7fffff07, RandMagic = 1634357784;

// A random integer in [0, RandMax].
// L'Ecuyer, Efficient and portable combined random number generators
Integer RandomQ::RandBase() {
   Assert(_MagicV == RandMagic, "RandomQ::RandBase: this is not correctly initialized!");
   Integer K0 = _S0/53668; _S0 = 40014*(_S0 - K0*53668) - K0*12211; if (_S0 < 0) _S0 += RandM0;
   Integer K1 = _S1/52774; _S1 = 40692*(_S1 - K1*52774) - K1*3791; if (_S1 < 0) _S1 += RandM1;
   Integer Base = _S0 - _S1; return Base < 1? Base += RandMax: --Base;
}

// RandomQ initialization with random values which come from a standard random number generator.
// AlgLib: Copyright 02.12.2009 by Sergey Bochkanov
RandomQ::RandomQ() {
   _S0 = RandomZ(RandM0), _S1 = RandomZ(RandM1);
// Protection against negative seeds: SEED = -(SEED + 1).
// We can use just "-SEED" because there exists such integer number N that N < 0, -N = N < 0 too.
// (This number is equal to 0x800...000).
// Need to handle such seed correctly forces us to use a bit more complicated formula.
   if (_S0 < 0) _S0 = -(_S0 + 1); _S0 = _S0%(RandM0 - 1) + 1;
   if (_S1 < 0) _S1 = -(_S1 + 1); _S1 = _S1%(RandM1 - 1) + 1;
   _MagicV = RandMagic;
}

// A random real number in (0, 1), not including interval boundaries.
// AlgLib: Copyright 02.12.2009 by Sergey Bochkanov
double RandomQ::RandomUniformR() { return (double)(RandBase() + 1)/(double)(RandMax + 2); }
double RandomQ::RandMidUniformR() { return (double)(2*RandBase() - RandMax)/(double)(RandMax + 2); }

// A random integer number in [0, N).
// N can be any positive number except for very large numbers:
// *	close to 2^31 on 32-bit systems
// *	close to 2^62 on 64-bit systems
// An exception will be generated if N is too large.
// AlgLib: Copyright 02.12.2009 by Sergey Bochkanov
Integer RandomQ::RandomUniformZ(Integer N) {
   Assert(N > 0, "RandomQ::RandomUniformZ: N <= 0!");
   Integer HiN = RandMax + 1;
   if (N > HiN) { // Reduce the problem on the interval spanning [0, N) to several subproblems on intervals spanning [0, HiN).
      Integer Nq = N/HiN; if (N%HiN > 0) Nq++; Assert(Nq <= HiN, "RandomQ::RandomUniformZ: N is too large");
   // Produce a random B in Nq, a random offset A within bin B
   // and filter out any results at or beyond N to avoid bias in the result.
      Integer A, B, C;
      do A = RandomUniformZ(HiN), B = RandomUniformZ(Nq); while ((C = A + HiN*B) >= N);
      return C;
   } else { // We can not simply return "RandBase() mod N", but need to skew it for large N's in [0.1*RandMax...RandMax].
      Integer Hi = HiN - HiN%N, C;
      do C = RandBase(); while (C >= Hi);
      return C%N;
   }
}

// A random number from the normal distribution.
// AlgLib: Copyright 02.12.2009 by Sergey Bochkanov
double RandomQ::RandomNormal1() {
   double X, Y, R;
// Two independently distributed random numbers from the normal distribution.
   do X = RandMidUniformR(), Y = RandMidUniformR(), R = X*X + Y*Y; while (R <= 0.0 || R >= 1.0);
// Do the square roots separately to avoid overflow for small R.
   return sqrt(-2.0*log(R))/sqrt(R)*X;
}

// === Idw (Inverse Distance Weighted Interpolation) Testing Unit ===
#define QUIET true

// Test number A for equality to B (or 0) within the range Tiny.
// AlgLib: Copyright 02.12.2009 by Sergey Bochkanov
static inline bool NearR(double A, double B, double Tiny) { return fabs(A - B) < Tiny; }
static inline bool SmallR(double A, double Tiny) { return fabs(A) < Tiny; }

// Testing continuity properties: C0 (Deg == 0) or C1 (Deg == 1) continuity.
static bool TestIdwContinuity(IdwQ &Q, Integer Xs, Integer Ys, RVector &X0, RVector &X1, Integer Steps, Integer Deg) {
   bool Ok = true;
   Assert(Steps >= 10, "TestIdwContinuity: Steps is too small");
   Assert(Deg == 0 || Deg == 1, "TestIdwContinuity: incorrect Deg");
// Compute sequence of function values.
   RVector Xc(Xs); RMatrix Yv(Steps, Ys); RVector Yc;
   for (Integer s = 0; s < Steps; s++) {
      double T = (double)s/(double)(Steps - 1);
      for (Integer x = 0; x < Xs; x++) Xc._X[x] = X0._X[x]*T + X1._X[x]*(1 - T);
      Q.IdwCalcBuf(Xc, Yc); for (Integer y = 0; y < Ys; y++) Yv._XY[s][y] = Yc._X[y];
   }
// Evaluate all differentiability levels (C0, C1) requested by user.
   for (Integer IxC = 0; IxC <= Deg; IxC++) {
   // Compute Lipschitz constant for original and increased steps.
      double Lc0 = 0.0, Lc1 = 0.0;
      for (Integer s = 0; s < Steps - 2; s++) for (Integer y = 0; y < Ys; y++)
         Lc0 = Max(Lc0, fabs(Yv._XY[s][y] - Yv._XY[s + 1][y])), Lc1 = Max(Lc1, fabs(Yv._XY[s][y] - Yv._XY[s + 2][y])/2);
      Ok = Ok && (Lc1 <= 0.0001 || Lc0 <= 1.750*Lc1);
   // Differentiate the function, repeat one more time.
      for (Integer s = 0; s < Steps - 1; s++) for (Integer y = 0; y < Ys; y++) Yv._XY[s][y] = Yv._XY[s + 1][y] - Yv._XY[s][y];
      Steps--;
   }
   return Ok;
}

// Test MSTAB.
static bool TestIdwCommon() {
   RMatrix XY; RVector X, Xx, MuY, Y;
   RandomQ RQ;
   bool Ok = true;
   const double Tiny = 1.0E-10, LoR = 0.05;
// Try all algorithms
   for (Integer Op = 0; Op < 3; Op++) {
   // Test empty dataset
      for (Integer Xs = 1; Xs < 6; Xs++) for (Integer Ys = 1; Ys < 6; Ys++) {
         IdwB Qb(Xs, Ys);
         switch (Op) {
            case 0: Qb.IdbSetShepard(1 + (Xs + 1)*RQ.RandomUniformR()); break;
            case 1: Qb.IdbSetModShepard(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            case 2: Qb.IdbSetMStab(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            default: Assert(false, "TestIdwCommon: unexpected Op"); break;
         }
      // Fit and store result directly into the variable
         IdwR R(Qb); IdwBuf Buf(R._Q);
      // The test report.
         Ok = Ok && R._RmsErr == 0.0 && R._AveErr == 0.0 && R._MaxErr == 0.0 && R._R2 == 1.0;
      // Test simplified evaluation.
         double X0 = RQ.RandomNormal1(), X1 = RQ.RandomNormal1(), X2 = RQ.RandomNormal1();
         if (Xs == 1 && Ys == 1) Ok = Ok && R._Q.IdwCalc1(X0) == 0.0;
         if (Xs == 2 && Ys == 1) Ok = Ok && R._Q.IdwCalc2(X0, X1) == 0.0;
         if (Xs == 3 && Ys == 1) Ok = Ok && R._Q.IdwCalc3(X0, X1, X2) == 0.0;
      // Test generic evaluation.
         X.ReSizeVector(Xs); for (Integer x = 0; x < Xs; x++) X._X[x] = RQ.RandomNormal1();
         Y.ReSizeVector(0), R._Q.IdwCalc(X, Y), Ok = Ok && Y._N == Ys; if (!Ok) return Ok;
         for (Integer y = 0; y < Ys; y++) Ok = Ok && Y._X[y] == 0.0;
         Y.ReSizeVector(0), R._Q.IdwCalcBuf(X, Y), Ok = Ok && Y._N == Ys; if (!Ok) return Ok;
         for (Integer y = 0; y < Ys; y++) Ok = Ok && Y._X[y] == 0.0;
         Y.ReSizeVector(0), R._Q.IdwCalcTs(Buf, X, Y), Ok = Ok && Y._N == Ys; if (!Ok) return Ok;
         for (Integer y = 0; y < Ys; y++) Ok = Ok && Y._X[y] == 0.0;
      }
   // Generate a random data set with distinct points, test interpolation properties (the target function is reproduced almost exactly, the model is continuous).
      for (Integer Pass = 0; Pass < 20; Pass++) {
         Integer Ps = 1 + RQ.RandomUniformZ(25);
         Integer Xs = 1 + RQ.RandomUniformZ(4), Ys = 1 + RQ.RandomUniformZ(4);
         X.ReSizeVector(Xs), Xx.ReSizeVector(Xs), Y.ReSizeVector(Ys);
      // Generate a data set with distinct points.
         XY.ReSizeMatrix(Ps, Xs + Ys);
         MuY.ReSizeVector(Ys); for (Integer y = 0; y < Ys; y++) MuY._X[y] = 0.0;
         for (Integer p = 0; p < Ps; ) {
         // Generate a random point.
            for (Integer x = 0; x < Xs; x++) XY._XY[p][x] = RQ.RandomNormal1();
         // Test the distance between the newly-generated point and other ones.
         // Repeat point generation if it is too close to some other point.
            double V = MaxPosReal;
            for (Integer p0 = 0; p0 < p; p0++) {
               double V2 = 0.0;
               for (Integer x = 0; x < Xs; x++) V2 = Max(V2, fabs(XY._XY[p][x] - XY._XY[p0][x]));
               V = Min(V, V2);
            }
            if (V >= LoR) { // The point is accepted.
               for (Integer y = 0; y < Ys; y++) XY._XY[p][Xs + y] = RQ.RandomNormal1(), MuY._X[y] += XY._XY[p][Xs + y]/Ps;
               p++;
            }
         }
      // Build the IDW model.
         IdwB Qb(Xs, Ys);
         switch (Op) {
            case 0: Qb.IdbSetShepard(1 + (Xs + 1)*RQ.RandomUniformR()); break;
            case 1: Qb.IdbSetModShepard(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            case 2: Qb.IdbSetMStab(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            default: Assert(false, "TestIdwCommon: unexpected Op");
         }
         Qb.IdbSetPoints(XY, Ps);
      // Fit and store result directly into the variable.
         IdwR R(Qb); IdwBuf Buf(R._Q);
      // Test error metrics.
      // NOTE:
      // *	We expect that dataset is reproduced exactly.
         Ok = Ok && SmallR(R._RmsErr, Tiny) && SmallR(R._AveErr, Tiny) && SmallR(R._MaxErr, Tiny) && !SmallR(R._R2, 1.0 - Tiny);
      // Test that dataset is actually exactly reproduced.
         for (Integer p = 0; p < Ps; p++) {
         // Test generic evaluation.
            for (Integer x = 0; x < Xs; x++) X._X[x] = XY._XY[p][x];
            Y.ReSizeVector(0), R._Q.IdwCalc(X, Y), Ok = Ok && Y._N == Ys;
            for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], XY._XY[p][Xs + y], Tiny);
            Integer K = RQ.RandomUniformZ(2*Ys + 1);
            Y.ReSizeVector(K); for (Integer n = 0; n < Y._N; n++) Y._X[n] = 0.0;
            R._Q.IdwCalcBuf(X, Y), Ok = Ok && Y._N == Max(Ys, K);
            for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], XY._XY[p][Xs + y], Tiny);
            K = RQ.RandomUniformZ(2*Ys + 1);
            Y.ReSizeVector(K); for (Integer n = 0; n < Y._N; n++) Y._X[n] = 0.0;
            R._Q.IdwCalcTs(Buf, X, Y), Ok = Ok && Y._N == Max(Ys, K);
            for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], XY._XY[p][Xs + y], Tiny);
         // Specialized 1, 2, 3-dimensional cases
            if (Ys == 1) switch (Xs) {
               case 1: Ok = Ok && NearR(R._Q.IdwCalc1(X._X[0]), XY._XY[p][Xs], Tiny); break;
               case 2: Ok = Ok && NearR(R._Q.IdwCalc2(X._X[0], X._X[1]), XY._XY[p][Xs], Tiny); break;
               case 3: Ok = Ok && NearR(R._Q.IdwCalc3(X._X[0], X._X[1], X._X[2]), XY._XY[p][Xs], Tiny); break;
            }
         }
      // Test continuity properties:
      // *	continuity is guaranteed for the original Shepard's method, MSTAB and MSMOOTH.
      // *	the modified Shepard method does not guarantee continuity of the model,
      //	but we can be sure that the model is continuous along the line connecting the two nearest points.
         for (Integer K = 0; K < 2; K++) {
            Integer p0 = RQ.RandomUniformZ(Ps);
            for (Integer x = 0; x < Xs; x++) X._X[x] = XY._XY[p0][x];
            Integer p1 = -1;
            double V = MaxPosReal;
            for (Integer p = 0; p < Ps; p++) {
               double V2 = 0.0;
               for (Integer x = 0; x < Xs; x++) V2 += Squ(X._X[x] - XY._XY[p][x]);
               if (V2 < V && V2 > 0.0) {
                  p1 = p;
                  for (Integer x = 0; x < Xs; x++) Xx._X[x] = XY._XY[p][x];
                  V = V2;
               }
            }
            if (p1 < 0) {
               p1 = RQ.RandomUniformZ(Ps);
               for (Integer x = 0; x < Xs; x++) Xx._X[x] = XY._XY[p1][x];
            }
            Integer Steps = Op == 0? 0: Op == 1? -1: +1;
            if (Steps >= 0) Ok = Ok && TestIdwContinuity(R._Q, Xs, Ys, X, Xx, 10000, Steps);
         }
      // Test evaluation at remote points.
         X.ReSizeVector(Xs); for (Integer x = 0; x < Xs; x++) X._X[x] = 1.0E20*(2*RQ.RandomUniformZ(2) - 1);
         R._Q.IdwCalc(X, Y);
         for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], MuY._X[y], Tiny);
      }
   // Generate a random data set with NONdistinct points, test approximation properties and error reports.
      for (Integer Pass = 0; Pass < 20; Pass++) {
         Integer Ps = 2*(1 + RQ.RandomUniformZ(10)), Xs = 1 + RQ.RandomUniformZ(4), Ys = 1 + RQ.RandomUniformZ(4);
         X.ReSizeVector(Xs), Xx.ReSizeVector(Xs);
         Y.ReSizeVector(Ys), MuY.ReSizeVector(Ys);
      // Generate a data set with non distinct points, each point is repeated; compute reference values of the error metrics.
         XY.ReSizeMatrix(Ps, Xs + Ys);
         double RmsErr = 0.0, AveErr = 0.0, MaxErr = 0.0;
         double Rss = 0.0, Tss = 0.0;
         for (Integer y = 0; y < Ys; y++) MuY._X[y] = 0.0;
         for (Integer p = 0; p < Ps/2; ) {
         // Generate two copies of the same point.
            for (Integer xy = 0; xy < Xs + Ys; xy++) XY._XY[2*p + 1][xy] = XY._XY[2*p][xy] = RQ.RandomNormal1();
         // Test the distance between the newly-generated point and other ones.
         // Repeat point generation if it is too close to some other point.
            double V = MaxPosReal;
            for (Integer p0 = 0; p0 < 2*p; p0++) {
               double V2 = 0.0;
               for (Integer x = 0; x < Xs; x++) V2 = Max(V2, fabs(XY._XY[2*p][x] - XY._XY[p0][x]));
               V = Min(V, V2);
            }
            if (V >= LoR) { // Update MuY.
               for (Integer y = 0; y < Ys; y++) MuY._X[y] += (XY._XY[2*p][Xs + y] + XY._XY[2*p + 1][Xs + y])/Ps;
            // Perturbation the target value.
               for (Integer y = 0; y < Ys; y++) {
                  double V = pow(2.0, RQ.RandomNormal1());
                  XY._XY[2*p][Xs + y] += V, XY._XY[2*p + 1][Xs + y] -= V;
                  V = fabs(V);
                  RmsErr += 2.0*V*V, AveErr += 2.0*V, MaxErr = Max(MaxErr, V), Rss += 2.0*V*V;
               }
            // Next point.
               p++;
            }
         }
         for (Integer p = 0; p < Ps; p++) for (Integer y = 0; y < Ys; y++) Tss += Squ(XY._XY[p][Xs + y] - MuY._X[y]);
         RmsErr = sqrt(RmsErr/(Ps*Ys)), AveErr /= Ps*Ys;
         double R2 = 1.0 - Rss/(Tss == 0.0? 1.0: Tss);
      // Build IDW model
         IdwB Qb(Xs, Ys);
         switch (Op) {
            case 0: Qb.IdbSetShepard(Xs*(1 + RQ.RandomUniformR())); break;
            case 1: Qb.IdbSetModShepard(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            case 2: Qb.IdbSetMStab(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            default: Assert(false, "TestIdwCommon: unexpected Op"); break;
         }
         Qb.IdbSetPoints(XY, Ps);
      // Fit and store the result directly into the variable.
         IdwR R(Qb);
      // Test error metrics
         Ok = Ok && NearR(R._RmsErr, RmsErr, Tiny) && NearR(R._AveErr, AveErr, Tiny) && NearR(R._MaxErr, MaxErr, Tiny) && NearR(R._R2, R2, Tiny);
      // Test the ability to reproduce the mean over non-distinct points.
      // NOTE:
      // *	We do not test all evaluation functions, just IdwCalc().
         for (Integer p = 0; p < Ps/2; p++) {
            for (Integer x = 0; x < Xs; x++) X._X[x] = XY._XY[2*p][x];
            R._Q.IdwCalc(X, Y); for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], 0.5*(XY._XY[2*p][Xs + y] + XY._XY[2*p + 1][Xs + y]), Tiny);
         }
      // Test continuity properties:
      // *	continuity is guaranteed for original Shepard's method, MSTAB and MSMOOTH,
      // *	the modified Shepard method does not guarantee continuity of the model,
      //	but we can be sure that the model is continuous along the line connecting the two nearest points.
         for (Integer K = 0; K < 2; K++) {
            Integer p0 = RQ.RandomUniformZ(Ps);
            for (Integer x = 0; x < Xs; x++) X._X[x] = XY._XY[p0][x];
            Integer p1 = -1;
            double V = MaxPosReal;
            for (Integer p = 0; p < Ps; p++) {
               double V2 = 0.0;
               for (Integer x = 0; x < Xs; x++) V2 += Squ(X._X[x] - XY._XY[p][x]);
               if (V2 < V && V2 > 0.0) {
                  p1 = p;
                  for (Integer x = 0; x < Xs; x++) Xx._X[x] = XY._XY[p][x];
                  V = V2;
               }
            }
            if (p1 < 0) {
               p1 = RQ.RandomUniformZ(Ps);
               for (Integer x = 0; x < Xs; x++) Xx._X[x] = XY._XY[p1][x];
            }
            Integer Steps = Op == 0? 0: Op == 1? -1: +1;
            if (Steps >= 0) Ok = Ok && TestIdwContinuity(R._Q, Xs, Ys, X, Xx, 10000, Steps);
         }
      }
   // Test correct handling of the prior term.
      Integer Ps = 10;
      for (Integer Pass = 0; Pass < 20; Pass++) {
         Integer Xs = 1 + RQ.RandomUniformZ(4), Ys = 1 + RQ.RandomUniformZ(4);
         X.ReSizeVector(Xs), Y.ReSizeVector(Ys);
         MuY.ReSizeVector(Ys); for (Integer y = 0; y < Ys; y++) MuY._X[y] = 0.0;
         XY.ReSizeMatrix(Ps, Xs + Ys);
         for (Integer p = 0; p < Ps; p++) {
            for (Integer x = 0; x < Xs; x++) XY._XY[p][x] = RQ.RandomNormal1();
            for (Integer y = 0; y < Ys; y++) XY._XY[p][Xs + y] = RQ.RandomNormal1(), MuY._X[y] += XY._XY[p][Xs + y]/Ps;
         }
         IdwB Qb(Xs, Ys);
         switch (Op) {
            case 0: Qb.IdbSetShepard(1 + (Xs + 1)*RQ.RandomUniformR()); break;
            case 1: Qb.IdbSetModShepard(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            case 2: Qb.IdbSetMStab(pow(2.0, 2.0*RQ.RandMidUniformR())); break;
            default: Assert(false, "TestIdwCommon: unexpected Op (prior test)"); break;
         }
         Qb.IdbSetPoints(XY, Ps);
      // Zero prior (not tested with the textbook Shepard method).
         if (Op != 0) {
            Qb.IdbSetZeroTerm(); IdwR R(Qb); for (Integer x = 0; x < Xs; x++) X._X[x] = 1.0E20*(2*RQ.RandomUniformZ(2) - 1);
            R._Q.IdwCalc(X, Y); for (Integer y = 0; y < Ys; y++) Ok = Ok && SmallR(Y._X[y], Tiny);
         }
      {
      // Mean prior.
         Qb.IdbSetConstTerm(); IdwR R(Qb);
         for (Integer x = 0; x < Xs; x++) X._X[x] = 1.0E20*(2*RQ.RandomUniformZ(2) - 1);
         R._Q.IdwCalc(X, Y); for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], MuY._X[y], Tiny);
      }
      // User-specified prior (not tested with the textbook Shepard method).
         if (Op != 0) {
            double V = RQ.RandomNormal1();
            Qb.IdbSetUserTerm(V); IdwR R(Qb); for (Integer x = 0; x < Xs; x++) X._X[x] = 1.0E20*(2*RQ.RandomUniformZ(2) - 1);
            R._Q.IdwCalc(X, Y); for (Integer y = 0; y < Ys; y++) Ok = Ok && NearR(Y._X[y], V, Tiny);
         }
      }
   }
   return Ok;
}

// Test MSTAB.
static bool TestIdwMStab() {
   RandomQ RQ;
   bool Ok = true;
// Basic test #1: non-zero derivative.
// *	XY = [[-1,-1],[0,0,],[1,1]]
// *	Rad0 >= 2
// *	derivative at x == 0 must be positive and bigger than 0.1
   RMatrix XY(3, 2);
   for (Integer I = 0; I < 3; I++) XY._XY[I][1] = XY._XY[I][0] = (double)(I - 1);
{
   IdwB Qb(1, 1); Qb.IdbSetMStab(pow(2.0, 1.0 + RQ.RandomUniformR())), Qb.IdbSetPoints(XY, 3);
   IdwR R(Qb);
   double V = 0.01;
   Ok = Ok && (R._Q.IdwCalc1(V) - R._Q.IdwCalc1(-V))/(2*V) >= 0.1;
}
// Basic test #2: good smoothness
// *	2D task, the data set is composed from 3 parallel lines along Y == -0.1, Y == 0 and Y == +0.1,
//	with the outer lines having constant zero target value, the inner line having constant target equal to 1.
// *	Rad0 == 1 is used
// *	we test that the function value does not change significantly along the line.
   Integer Ps = 100;
   XY.ReSizeMatrix(3*Ps, 3);
   for (Integer n = 0; n < Ps; n++) {
      XY._XY[3*n][0] = (double)n/(double)(Ps - 1), XY._XY[3*n][1] = -0.1, XY._XY[3*n][2] = 0.0;
      XY._XY[3*n + 1][0] = (double)n/(double)(Ps - 1), XY._XY[3*n + 1][1] = 0.0, XY._XY[3*n + 1][2] = 1.0;
      XY._XY[3*n + 2][0] = (double)n/(double)(Ps - 1), XY._XY[3*n + 2][1] = 0.1, XY._XY[3*n + 2][2] = 0.0;
   }
{
   IdwB Qb(2, 1); Qb.IdbSetMStab(1.0), Qb.IdbSetPoints(XY, 3*Ps);
   IdwR R(Qb);
   double MaxV = 0.0;
   for (Integer I = 0; I <= 1000; I++)
      MaxV = Max4(MaxV, fabs(R._Q.IdwCalc2(RQ.RandomUniformR(), -0.1)), fabs(R._Q.IdwCalc2(RQ.RandomUniformR(), 0.1)), fabs(R._Q.IdwCalc2(RQ.RandomUniformR(), 0.0) - 1));
   Ok = Ok && MaxV <= 0.001;
}
// Continuity when moving away from the data set.
   XY.ReSizeMatrix(1, 2), XY._XY[0][0] = 0.0, XY._XY[0][1] = 1.0;
   double Rad0 = 1.0;
{
   IdwB Qb(1, 1); Qb.IdbSetMStab(Rad0), Qb.IdbSetPoints(XY, 1), Qb.IdbSetZeroTerm();
   IdwR R(Qb);
   Ok = Ok && R._Q.IdwCalc1(100000.0) == 0.0;
   double MaxDv = 0.0;
   for (Integer I = 0; I <= 500; I++) {
      double X0 = 1.2*Rad0*(I/500.0), X1 = 1.2*Rad0*((I + 1)/500.0);
      MaxDv = Max(MaxDv, fabs((R._Q.IdwCalc1(X1) - R._Q.IdwCalc1(X0))/(X1 - X0)));
   }
   double V2 = 0.0;
   for (Integer I = 0; I <= 1000; I++) {
      double X0 = 1.2*Rad0*(I/1000.0), X1 = 1.2*Rad0*((I + 1)/1000.0);
      V2 = Max(V2, fabs((R._Q.IdwCalc1(X1) - R._Q.IdwCalc1(X0))/(X1 - X0)));
   }
   Ok = Ok && V2/MaxDv <= 1.333;
}
   return Ok;
}

// Testing IDW interpolation.
int main(int AC, char **AV) {
   unsigned Seed;
   if (AC == 2) Seed = (unsigned)atoi(AV[1]); else { time_t Now; Seed = (unsigned)time(&Now); }
   printf("Seed: %u\n", (unsigned)Seed), SRand(Seed), fflush(stdout);
// Now we are ready to test!
   time_t BegT; time(&BegT);
   Seed = CurSeed;
// HPC AlgLib uses multithreaded code for the second TestFn, all of which are added to the list above.
   bool CommonOk = TestIdwCommon(), MStabOk = TestIdwMStab();
   bool Ok = CommonOk && MStabOk;
   if (!Ok)
      printf("Inverse Distance Weighting Tests\n"),
      printf("* Common Properties                      %s\n", CommonOk? "Ok": "Failed"),
      printf("* MSTAB-Specific Tests                   %s\n", MStabOk? "Ok": "Failed"),
      printf("Test %s\n", Ok? "Passed": "Failed");
   printf("[%08x]: %-32s %s\n", Seed, "Idw", Ok? "Ok": "Failed"), fflush(stdout);
   time_t EndT; time(&EndT), printf("Done in %ld seconds\n", (long)difftime(EndT, BegT));
// Return the result.
   return Ok? 0: 1;
}
