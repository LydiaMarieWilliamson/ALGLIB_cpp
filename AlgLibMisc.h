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
#ifndef OnceOnlyAlgLibMisc_h
#define OnceOnlyAlgLibMisc_h

#include "AlgLibInternal.h"

// === NEARESTNEIGHBOR Package ===
// Depends on: (AlgLibInternal) SCODES, TSORT
namespace alglib_impl {
struct kdtreerequestbuffer {
   ae_vector x;
   ae_vector boxmin;
   ae_vector boxmax;
   ae_int_t kneeded;
   double rneeded;
   bool selfmatch;
   double approxf;
   ae_int_t kcur;
   ae_vector idx;
   ae_vector r;
   ae_vector buf;
   ae_vector curboxmin;
   ae_vector curboxmax;
   double curdist;
};
void kdtreerequestbuffer_init(void *_p, bool make_automatic);
void kdtreerequestbuffer_copy(void *_dst, void *_src, bool make_automatic);
void kdtreerequestbuffer_free(void *_p, bool make_automatic);

struct kdtree {
   ae_int_t n;
   ae_int_t nx;
   ae_int_t ny;
   ae_int_t normtype;
   ae_matrix xy;
   ae_vector tags;
   ae_vector boxmin;
   ae_vector boxmax;
   ae_vector nodes;
   ae_vector splits;
   kdtreerequestbuffer innerbuf;
   ae_int_t debugcounter;
};
void kdtree_init(void *_p, bool make_automatic);
void kdtree_copy(void *_dst, void *_src, bool make_automatic);
void kdtree_free(void *_p, bool make_automatic);
void kdtreealloc(ae_serializer *s, kdtree *tree);
void kdtreeserialize(ae_serializer *s, kdtree *tree);
void kdtreeunserialize(ae_serializer *s, kdtree *tree);

void kdtreebuild(RMatrix *xy, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree *kdt);
void kdtreebuildtagged(RMatrix *xy, ZVector *tags, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree *kdt);
void kdtreecreaterequestbuffer(kdtree *kdt, kdtreerequestbuffer *buf);
ae_int_t kdtreequeryknn(kdtree *kdt, RVector *x, ae_int_t k, bool selfmatch);
ae_int_t kdtreetsqueryknn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, ae_int_t k, bool selfmatch);
ae_int_t kdtreequeryrnn(kdtree *kdt, RVector *x, double r, bool selfmatch);
ae_int_t kdtreequeryrnnu(kdtree *kdt, RVector *x, double r, bool selfmatch);
ae_int_t kdtreetsqueryrnn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch);
ae_int_t kdtreetsqueryrnnu(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, double r, bool selfmatch);
ae_int_t kdtreequeryaknn(kdtree *kdt, RVector *x, ae_int_t k, bool selfmatch, double eps);
ae_int_t kdtreetsqueryaknn(kdtree *kdt, kdtreerequestbuffer *buf, RVector *x, ae_int_t k, bool selfmatch, double eps);
ae_int_t kdtreequerybox(kdtree *kdt, RVector *boxmin, RVector *boxmax);
ae_int_t kdtreetsquerybox(kdtree *kdt, kdtreerequestbuffer *buf, RVector *boxmin, RVector *boxmax);
void kdtreequeryresultsx(kdtree *kdt, RMatrix *x);
void kdtreequeryresultsxy(kdtree *kdt, RMatrix *xy);
void kdtreequeryresultstags(kdtree *kdt, ZVector *tags);
void kdtreequeryresultsdistances(kdtree *kdt, RVector *r);
void kdtreetsqueryresultsx(kdtree *kdt, kdtreerequestbuffer *buf, RMatrix *x);
void kdtreetsqueryresultsxy(kdtree *kdt, kdtreerequestbuffer *buf, RMatrix *xy);
void kdtreetsqueryresultstags(kdtree *kdt, kdtreerequestbuffer *buf, ZVector *tags);
void kdtreetsqueryresultsdistances(kdtree *kdt, kdtreerequestbuffer *buf, RVector *r);
void kdtreequeryresultsxi(kdtree *kdt, RMatrix *x);
void kdtreequeryresultsxyi(kdtree *kdt, RMatrix *xy);
void kdtreequeryresultstagsi(kdtree *kdt, ZVector *tags);
void kdtreequeryresultsdistancesi(kdtree *kdt, RVector *r);
void kdtreeexplorebox(kdtree *kdt, RVector *boxmin, RVector *boxmax);
void kdtreeexplorenodetype(kdtree *kdt, ae_int_t node, ae_int_t *nodetype);
void kdtreeexploreleaf(kdtree *kdt, ae_int_t node, RMatrix *xy, ae_int_t *k);
void kdtreeexploresplit(kdtree *kdt, ae_int_t node, ae_int_t *d, double *s, ae_int_t *nodele, ae_int_t *nodege);
} // end of namespace alglib_impl

namespace alglib {
// Buffer object which is used to perform nearest neighbor  requests  in  the
// multithreaded mode (multiple threads working with same KD-tree object).
//
// This object should be created with KDTreeCreateRequestBuffer().
DecClass(kdtreerequestbuffer, EndD);

// KD-tree object.
DecClass(kdtree, EndD);

// These functions serialize a data structure to a C++ string or stream.
// * serialization can be freely moved across 32-bit and 64-bit systems,
//   and different byte orders. For example, you can serialize a string
//   on a SPARC and unserialize it on an x86.
// * ALGLIB++ serialization is compatible with serialization in ALGLIB,
//   in both directions.
// Important properties of s_out:
// * it contains alphanumeric characters, dots, underscores, minus signs
// * these symbols are grouped into words, which are separated by spaces
//   and Windows-style (CR+LF) newlines
void kdtreeserialize(kdtree &obj, std::string &s_out);
void kdtreeserialize(kdtree &obj, std::ostream &s_out);

// These functions unserialize a data structure from a C++ string or stream.
// Important properties of s_in:
// * any combination of spaces, tabs, Windows or Unix stype newlines can
//   be used as separators, so as to allow flexible reformatting of the
//   stream or string from text or XML files.
// * But you should not insert separators into the middle of the "words"
//   nor you should change case of letters.
void kdtreeunserialize(const std::string &s_in, kdtree &obj);
void kdtreeunserialize(const std::istream &s_in, kdtree &obj);

// KD-tree creation
//
// This subroutine creates KD-tree from set of X-values and optional Y-values
//
// Inputs:
//     XY      -   dataset, array[0..N-1,0..NX+NY-1].
//                 one row corresponds to one point.
//                 first NX columns contain X-values, next NY (NY may be zero)
//                 columns may contain associated Y-values
//     N       -   number of points, N >= 0.
//     NX      -   space dimension, NX >= 1.
//     NY      -   number of optional Y-values, NY >= 0.
//     NormType-   norm type:
//                 * 0 denotes infinity-norm
//                 * 1 denotes 1-norm
//                 * 2 denotes 2-norm (Euclidean norm)
//
// Outputs:
//     KDT     -   KD-tree
//
//
// NOTES
//
// 1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
//    requirements.
// 2. Although KD-trees may be used with any combination of N  and  NX,  they
//    are more efficient than brute-force search only when N >> 4^NX. So they
//    are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
//    inefficient case, because  simple  binary  search  (without  additional
//    structures) is much more efficient in such tasks than KD-trees.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreebuild(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuild(const real_2d_array &xy, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);

// KD-tree creation
//
// This  subroutine  creates  KD-tree  from set of X-values, integer tags and
// optional Y-values
//
// Inputs:
//     XY      -   dataset, array[0..N-1,0..NX+NY-1].
//                 one row corresponds to one point.
//                 first NX columns contain X-values, next NY (NY may be zero)
//                 columns may contain associated Y-values
//     Tags    -   tags, array[0..N-1], contains integer tags associated
//                 with points.
//     N       -   number of points, N >= 0
//     NX      -   space dimension, NX >= 1.
//     NY      -   number of optional Y-values, NY >= 0.
//     NormType-   norm type:
//                 * 0 denotes infinity-norm
//                 * 1 denotes 1-norm
//                 * 2 denotes 2-norm (Euclidean norm)
//
// Outputs:
//     KDT     -   KD-tree
//
// NOTES
//
// 1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
//    requirements.
// 2. Although KD-trees may be used with any combination of N  and  NX,  they
//    are more efficient than brute-force search only when N >> 4^NX. So they
//    are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
//    inefficient case, because  simple  binary  search  (without  additional
//    structures) is much more efficient in such tasks than KD-trees.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);

// This function creates buffer  structure  which  can  be  used  to  perform
// parallel KD-tree requests.
//
// KD-tree subpackage provides two sets of request functions - ones which use
// internal buffer of KD-tree object  (these  functions  are  single-threaded
// because they use same buffer, which can not shared between  threads),  and
// ones which use external buffer.
//
// This function is used to initialize external buffer.
//
// Inputs:
//     KDT         -   KD-tree which is associated with newly created buffer
//
// Outputs:
//     Buf         -   external buffer.
//
//
// IMPORTANT: KD-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use buffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
void kdtreecreaterequestbuffer(const kdtree &kdt, kdtreerequestbuffer &buf);

// K-NN query: K nearest neighbors
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: KDTreeTsQueryKNN() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// these results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k);

// K-NN query: K nearest neighbors, using external thread-local buffer.
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   kd-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
ae_int_t kdtreetsqueryknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k);

// R-NN query: all points within R-sphere centered at X, ordered by  distance
// between point and X (by ascending).
//
// NOTE: it is also possible to perform undordered queries performed by means
//       of kdtreequeryrnnu() and kdtreetsqueryrnnu() functions. Such queries
//       are faster because we do not have to use heap structure for sorting.
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: kdtreetsqueryrnn() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// actual results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r);

// R-NN query: all points within R-sphere  centered  at  X,  no  ordering  by
// distance as undicated by "U" suffix (faster that ordered query, for  large
// queries - significantly faster).
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: kdtreetsqueryrnn() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// actual results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
//
// As indicated by "U" suffix, this function returns unordered results.
// ALGLIB: Copyright 01.11.2018 by Sergey Bochkanov
ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreequeryrnnu(const kdtree &kdt, const real_1d_array &x, const double r);

// R-NN query: all points within  R-sphere  centered  at  X,  using  external
// thread-local buffer, sorted by distance between point and X (by ascending)
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// NOTE: it is also possible to perform undordered queries performed by means
//       of kdtreequeryrnnu() and kdtreetsqueryrnnu() functions. Such queries
//       are faster because we do not have to use heap structure for sorting.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreetsqueryrnn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r);

// R-NN query: all points within  R-sphere  centered  at  X,  using  external
// thread-local buffer, no ordering by distance as undicated  by  "U"  suffix
// (faster that ordered query, for large queries - significantly faster).
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     R           -   radius of sphere (in corresponding norm), R > 0
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//
// Result:
//     number of neighbors found, >= 0
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// As indicated by "U" suffix, this function returns unordered results.
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r, const bool selfmatch);
ae_int_t kdtreetsqueryrnnu(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const double r);

// K-NN query: approximate K nearest neighbors
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: KDTreeTsQueryAKNN() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//     Eps         -   approximation factor, Eps >= 0. eps-approximate  nearest
//                     neighbor  is  a  neighbor  whose distance from X is at
//                     most (1+eps) times distance of true nearest neighbor.
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// NOTES
//     significant performance gain may be achieved only when Eps  is  is  on
//     the order of magnitude of 1 or larger.
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// these results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() to get distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const double eps);

// K-NN query: approximate K nearest neighbors, using thread-local buffer.
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     X           -   point, array[0..NX-1].
//     K           -   number of neighbors to return, K >= 1
//     SelfMatch   -   whether self-matches are allowed:
//                     * if True, nearest neighbor may be the point itself
//                       (if it exists in original dataset)
//                     * if False, then only points with non-zero distance
//                       are returned
//                     * if not given, considered True
//     Eps         -   approximation factor, Eps >= 0. eps-approximate  nearest
//                     neighbor  is  a  neighbor  whose distance from X is at
//                     most (1+eps) times distance of true nearest neighbor.
//
// Result:
//     number of actual neighbors found (either K or N, if K > N).
//
// NOTES
//     significant performance gain may be achieved only when Eps  is  is  on
//     the order of magnitude of 1 or larger.
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "buf" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() to get distances
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 18.03.2016 by Sergey Bochkanov
ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
ae_int_t kdtreetsqueryaknn(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &x, const ae_int_t k, const double eps);

// Box query: all points within user-specified box.
//
// IMPORTANT: this function can not be used in multithreaded code because  it
//            uses internal temporary buffer of kd-tree object, which can not
//            be shared between multiple threads.  If  you  want  to  perform
//            parallel requests, use function  which  uses  external  request
//            buffer: KDTreeTsQueryBox() ("Ts" stands for "thread-safe").
//
// Inputs:
//     KDT         -   KD-tree
//     BoxMin      -   lower bounds, array[0..NX-1].
//     BoxMax      -   upper bounds, array[0..NX-1].
//
//
// Result:
//     number of actual neighbors found (in [0,N]).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures of the KD-tree. You can use  following  subroutines  to  obtain
// these results:
// * KDTreeQueryResultsX() to get X-values
// * KDTreeQueryResultsXY() to get X- and Y-values
// * KDTreeQueryResultsTags() to get tag values
// * KDTreeQueryResultsDistances() returns zeros for this request
//
// NOTE: this particular query returns unordered results, because there is no
//       meaningful way of  ordering  points.  Furthermore,  no 'distance' is
//       associated with points - it is either INSIDE  or OUTSIDE (so request
//       for distances will return zeros).
// ALGLIB: Copyright 14.05.2016 by Sergey Bochkanov
ae_int_t kdtreequerybox(const kdtree &kdt, const real_1d_array &boxmin, const real_1d_array &boxmax);

// Box query: all points within user-specified box, using thread-local buffer.
//
// You can call this function from multiple threads for same kd-tree instance,
// assuming that different instances of buffer object are passed to different
// threads.
//
// Inputs:
//     KDT         -   KD-tree
//     Buf         -   request buffer  object  created  for  this  particular
//                     instance of kd-tree structure with kdtreecreaterequestbuffer()
//                     function.
//     BoxMin      -   lower bounds, array[0..NX-1].
//     BoxMax      -   upper bounds, array[0..NX-1].
//
// Result:
//     number of actual neighbors found (in [0,N]).
//
// This  subroutine  performs  query  and  stores  its result in the internal
// structures  of  the  buffer object. You can use following  subroutines  to
// obtain these results (pay attention to "ts" in their names):
// * KDTreeTsQueryResultsX() to get X-values
// * KDTreeTsQueryResultsXY() to get X- and Y-values
// * KDTreeTsQueryResultsTags() to get tag values
// * KDTreeTsQueryResultsDistances() returns zeros for this query
//
// NOTE: this particular query returns unordered results, because there is no
//       meaningful way of  ordering  points.  Furthermore,  no 'distance' is
//       associated with points - it is either INSIDE  or OUTSIDE (so request
//       for distances will return zeros).
//
// IMPORTANT: kd-tree buffer should be used only with  KD-tree  object  which
//            was used to initialize buffer. Any attempt to use biffer   with
//            different object is dangerous - you  may  get  integrity  check
//            failure (exception) because sizes of internal arrays do not fit
//            to dimensions of KD-tree structure.
// ALGLIB: Copyright 14.05.2016 by Sergey Bochkanov
ae_int_t kdtreetsquerybox(const kdtree &kdt, const kdtreerequestbuffer &buf, const real_1d_array &boxmin, const real_1d_array &boxmax);

// X-values from last query.
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultsx().
//
// Inputs:
//     KDT     -   KD-tree
//     X       -   possibly pre-allocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     X       -   rows are filled with X-values
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultsx(const kdtree &kdt, real_2d_array &x);

// X- and Y-values from last query
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultsxy().
//
// Inputs:
//     KDT     -   KD-tree
//     XY      -   possibly pre-allocated buffer. If XY is too small to store
//                 result, it is resized. If size(XY) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     XY      -   rows are filled with points: first NX columns with
//                 X-values, next NY columns - with Y-values.
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultsxy(const kdtree &kdt, real_2d_array &xy);

// Tags from last query
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultstags().
//
// Inputs:
//     KDT     -   KD-tree
//     Tags    -   possibly pre-allocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     Tags    -   filled with tags associated with points,
//                 or, when no tags were supplied, with zeros
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultstags(const kdtree &kdt, integer_1d_array &tags);

// Distances from last query
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - kdtreetsqueryresultsdistances().
//
// Inputs:
//     KDT     -   KD-tree
//     R       -   possibly pre-allocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     R       -   filled with distances (in corresponding norm)
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultsdistances(const kdtree &kdt, real_1d_array &r);

// X-values from last query associated with kdtreerequestbuffer object.
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     X       -   possibly pre-allocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     X       -   rows are filled with X-values
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreetsqueryresultsx(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &x);

// X- and Y-values from last query associated with kdtreerequestbuffer object.
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     XY      -   possibly pre-allocated buffer. If XY is too small to store
//                 result, it is resized. If size(XY) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     XY      -   rows are filled with points: first NX columns with
//                 X-values, next NY columns - with Y-values.
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsTags()          tag values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreetsqueryresultsxy(const kdtree &kdt, const kdtreerequestbuffer &buf, real_2d_array &xy);

// Tags from last query associated with kdtreerequestbuffer object.
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - KDTreeTsqueryresultstags().
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     Tags    -   possibly pre-allocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     Tags    -   filled with tags associated with points,
//                 or, when no tags were supplied, with zeros
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsDistances()     distances
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreetsqueryresultstags(const kdtree &kdt, const kdtreerequestbuffer &buf, integer_1d_array &tags);

// Distances from last query associated with kdtreerequestbuffer object.
//
// This function retuns results stored in  the  internal  buffer  of  kd-tree
// object. If you performed buffered requests (ones which  use  instances  of
// kdtreerequestbuffer class), you  should  call  buffered  version  of  this
// function - KDTreeTsqueryresultsdistances().
//
// Inputs:
//     KDT     -   KD-tree
//     Buf     -   request  buffer  object  created   for   this   particular
//                 instance of kd-tree structure.
//     R       -   possibly pre-allocated buffer. If X is too small to store
//                 result, it is resized. If size(X) is enough to store
//                 result, it is left unchanged.
//
// Outputs:
//     R       -   filled with distances (in corresponding norm)
//
// NOTES
// 1. points are ordered by distance from the query point (first = closest)
// 2. if  XY is larger than required to store result, only leading part  will
//    be overwritten; trailing part will be left unchanged. So  if  on  input
//    XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
//    XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
//    you want function  to  resize  array  according  to  result  size,  use
//    function with same name and suffix 'I'.
//
// SEE ALSO
// * KDTreeQueryResultsX()             X-values
// * KDTreeQueryResultsXY()            X- and Y-values
// * KDTreeQueryResultsTags()          tag values
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreetsqueryresultsdistances(const kdtree &kdt, const kdtreerequestbuffer &buf, real_1d_array &r);

// X-values from last query; 'interactive' variant for languages like  Python
// which   support    constructs   like  "X = KDTreeQueryResultsXI(KDT)"  and
// interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultsxi(const kdtree &kdt, real_2d_array &x);

// XY-values from last query; 'interactive' variant for languages like Python
// which   support    constructs   like "XY = KDTreeQueryResultsXYI(KDT)" and
// interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultsxyi(const kdtree &kdt, real_2d_array &xy);

// Tags  from  last  query;  'interactive' variant for languages like  Python
// which  support  constructs  like "Tags = KDTreeQueryResultsTagsI(KDT)" and
// interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultstagsi(const kdtree &kdt, integer_1d_array &tags);

// Distances from last query; 'interactive' variant for languages like Python
// which  support  constructs   like  "R = KDTreeQueryResultsDistancesI(KDT)"
// and interactive mode of interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void kdtreequeryresultsdistancesi(const kdtree &kdt, real_1d_array &r);
} // end of namespace alglib

// === HQRND Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
struct hqrndstate {
   ae_int_t s1;
   ae_int_t s2;
   ae_int_t magicv;
};
void hqrndstate_init(void *_p, bool make_automatic);
void hqrndstate_copy(void *_dst, void *_src, bool make_automatic);
void hqrndstate_free(void *_p, bool make_automatic);

void hqrndrandomize(hqrndstate *state);
void hqrndseed(ae_int_t s1, ae_int_t s2, hqrndstate *state);
double hqrnduniformr(hqrndstate *state);
ae_int_t hqrnduniformi(hqrndstate *state, ae_int_t n);
double hqrndnormal(hqrndstate *state);
void hqrndunit2(hqrndstate *state, double *x, double *y);
void hqrndnormal2(hqrndstate *state, double *x1, double *x2);
double hqrndexponential(hqrndstate *state, double lambdav);
double hqrnddiscrete(hqrndstate *state, RVector *x, ae_int_t n);
double hqrndcontinuous(hqrndstate *state, RVector *x, ae_int_t n);
} // end of namespace alglib_impl

namespace alglib {
// Portable high quality random number generator state.
// Initialized with HQRNDRandomize() or HQRNDSeed().
//
// Fields:
//     S1, S2      -   seed values
//     V           -   precomputed value
//     MagicV      -   'magic' value used to determine whether State structure
//                     was correctly initialized.
DecClass(hqrndstate, EndD);

// HQRNDState  initialization  with  random  values  which come from standard
// RNG.
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void hqrndrandomize(hqrndstate &state);

// HQRNDState initialization with seed values
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void hqrndseed(const ae_int_t s1, const ae_int_t s2, hqrndstate &state);

// This function generates random real number in (0,1),
// not including interval boundaries
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
double hqrnduniformr(const hqrndstate &state);

// This function generates random integer number in [0, N)
//
// 1. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()
// 2. N can be any positive number except for very large numbers:
//    * close to 2^31 on 32-bit systems
//    * close to 2^62 on 64-bit systems
//    An exception will be generated if N is too large.
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
ae_int_t hqrnduniformi(const hqrndstate &state, const ae_int_t n);

// Random number generator: normal numbers
//
// This function generates one random number from normal distribution.
// Its performance is equal to that of HQRNDNormal2()
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
double hqrndnormal(const hqrndstate &state);

// Random number generator: random X and Y such that X^2+Y^2=1
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void hqrndunit2(const hqrndstate &state, double &x, double &y);

// Random number generator: normal numbers
//
// This function generates two independent random numbers from normal
// distribution. Its performance is equal to that of HQRNDNormal()
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 02.12.2009 by Sergey Bochkanov
void hqrndnormal2(const hqrndstate &state, double &x1, double &x2);

// Random number generator: exponential distribution
//
// State structure must be initialized with HQRNDRandomize() or HQRNDSeed().
// ALGLIB: Copyright 11.08.2007 by Sergey Bochkanov
double hqrndexponential(const hqrndstate &state, const double lambdav);

// This function generates  random number from discrete distribution given by
// finite sample X.
//
// Inputs:
//     State   -   high quality random number generator, must be
//                 initialized with HQRNDRandomize() or HQRNDSeed().
//         X   -   finite sample
//         N   -   number of elements to use, N >= 1
//
// Result:
//     this function returns one of the X[i] for random i=0..N-1
// ALGLIB: Copyright 08.11.2011 by Sergey Bochkanov
double hqrnddiscrete(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);

// This function generates random number from continuous  distribution  given
// by finite sample X.
//
// Inputs:
//     State   -   high quality random number generator, must be
//                 initialized with HQRNDRandomize() or HQRNDSeed().
//         X   -   finite sample, array[N] (can be larger, in this  case only
//                 leading N elements are used). THIS ARRAY MUST BE SORTED BY
//                 ASCENDING.
//         N   -   number of elements to use, N >= 1
//
// Result:
//     this function returns random number from continuous distribution which
//     tries to approximate X as mush as possible. min(X) <= Result <= max(X).
// ALGLIB: Copyright 08.11.2011 by Sergey Bochkanov
double hqrndcontinuous(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
} // end of namespace alglib

// === XDEBUG Package ===
namespace alglib_impl {
struct xdebugrecord1 {
   ae_int_t i;
   ae_complex c;
   ae_vector a;
};
void xdebugrecord1_init(void *_p, bool make_automatic);
void xdebugrecord1_copy(void *_dst, void *_src, bool make_automatic);
void xdebugrecord1_free(void *_p, bool make_automatic);

void xdebuginitrecord1(xdebugrecord1 *rec1);
ae_int_t xdebugb1count(BVector *a);
void xdebugb1not(BVector *a);
void xdebugb1appendcopy(BVector *a);
void xdebugb1outeven(ae_int_t n, BVector *a);
ae_int_t xdebugi1sum(ZVector *a);
void xdebugi1neg(ZVector *a);
void xdebugi1appendcopy(ZVector *a);
void xdebugi1outeven(ae_int_t n, ZVector *a);
double xdebugr1sum(RVector *a);
void xdebugr1neg(RVector *a);
void xdebugr1appendcopy(RVector *a);
void xdebugr1outeven(ae_int_t n, RVector *a);
ae_complex xdebugc1sum(CVector *a);
void xdebugc1neg(CVector *a);
void xdebugc1appendcopy(CVector *a);
void xdebugc1outeven(ae_int_t n, CVector *a);
ae_int_t xdebugb2count(BMatrix *a);
void xdebugb2not(BMatrix *a);
void xdebugb2transpose(BMatrix *a);
void xdebugb2outsin(ae_int_t m, ae_int_t n, BMatrix *a);
ae_int_t xdebugi2sum(ZMatrix *a);
void xdebugi2neg(ZMatrix *a);
void xdebugi2transpose(ZMatrix *a);
void xdebugi2outsin(ae_int_t m, ae_int_t n, ZMatrix *a);
double xdebugr2sum(RMatrix *a);
void xdebugr2neg(RMatrix *a);
void xdebugr2transpose(RMatrix *a);
void xdebugr2outsin(ae_int_t m, ae_int_t n, RMatrix *a);
ae_complex xdebugc2sum(CMatrix *a);
void xdebugc2neg(CMatrix *a);
void xdebugc2transpose(CMatrix *a);
void xdebugc2outsincos(ae_int_t m, ae_int_t n, CMatrix *a);
double xdebugmaskedbiasedproductsum(ae_int_t m, ae_int_t n, RMatrix *a, RMatrix *b, BMatrix *c);
} // end of namespace alglib_impl

namespace alglib {
DecClass(xdebugrecord1, ae_int_t &i; complex &c; real_1d_array a;);

// Debug functions to test the ALGLIB interface generator: not meant for use in production code.

// Creates and returns XDebugRecord1 structure:
// * integer and complex fields of Rec1 are set to 1 and 1+i correspondingly
// * array field of Rec1 is set to [2,3]
// ALGLIB: Copyright 27.05.2014 by Sergey Bochkanov
void xdebuginitrecord1(xdebugrecord1 &rec1);

// Counts number of True values in the boolean 1D array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
ae_int_t xdebugb1count(const boolean_1d_array &a);

// Replace all values in array by NOT(a[i]).
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugb1not(const boolean_1d_array &a);

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugb1appendcopy(boolean_1d_array &a);

// Generate N-element array with even-numbered elements set to True.
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugb1outeven(const ae_int_t n, boolean_1d_array &a);

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
ae_int_t xdebugi1sum(const integer_1d_array &a);

// Replace all values in array by -A[I]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugi1neg(const integer_1d_array &a);

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugi1appendcopy(integer_1d_array &a);

// Generate N-element array with even-numbered A[I] set to I, and odd-numbered
// ones set to 0.
//
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugi1outeven(const ae_int_t n, integer_1d_array &a);

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
double xdebugr1sum(const real_1d_array &a);

// Replace all values in array by -A[I]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugr1neg(const real_1d_array &a);

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugr1appendcopy(real_1d_array &a);

// Generate N-element array with even-numbered A[I] set to I*0.25,
// and odd-numbered ones are set to 0.
//
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugr1outeven(const ae_int_t n, real_1d_array &a);

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
complex xdebugc1sum(const complex_1d_array &a);

// Replace all values in array by -A[I]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugc1neg(const complex_1d_array &a);

// Appends copy of array to itself.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugc1appendcopy(complex_1d_array &a);

// Generate N-element array with even-numbered A[K] set to (x,y) = (K*0.25, K*0.125)
// and odd-numbered ones are set to 0.
//
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugc1outeven(const ae_int_t n, complex_1d_array &a);

// Counts number of True values in the boolean 2D array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
ae_int_t xdebugb2count(const boolean_2d_array &a);

// Replace all values in array by NOT(a[i]).
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugb2not(const boolean_2d_array &a);

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugb2transpose(boolean_2d_array &a);

// Generate MxN matrix with elements set to "sin(3*I+5*J) > 0"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugb2outsin(const ae_int_t m, const ae_int_t n, boolean_2d_array &a);

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
ae_int_t xdebugi2sum(const integer_2d_array &a);

// Replace all values in array by -a[i,j]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugi2neg(const integer_2d_array &a);

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugi2transpose(integer_2d_array &a);

// Generate MxN matrix with elements set to "Sign(sin(3*I+5*J))"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugi2outsin(const ae_int_t m, const ae_int_t n, integer_2d_array &a);

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
double xdebugr2sum(const real_2d_array &a);

// Replace all values in array by -a[i,j]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugr2neg(const real_2d_array &a);

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugr2transpose(real_2d_array &a);

// Generate MxN matrix with elements set to "sin(3*I+5*J)"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugr2outsin(const ae_int_t m, const ae_int_t n, real_2d_array &a);

// Returns sum of elements in the array.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
complex xdebugc2sum(const complex_2d_array &a);

// Replace all values in array by -a[i,j]
// Array is passed using "shared" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugc2neg(const complex_2d_array &a);

// Transposes array.
// Array is passed using "var" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugc2transpose(complex_2d_array &a);

// Generate MxN matrix with elements set to "sin(3*I+5*J),cos(3*I+5*J)"
// Array is passed using "out" convention.
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
void xdebugc2outsincos(const ae_int_t m, const ae_int_t n, complex_2d_array &a);

// Returns sum of a[i,j]*(1+b[i,j]) such that c[i,j] is True
// ALGLIB: Copyright 11.10.2013 by Sergey Bochkanov
double xdebugmaskedbiasedproductsum(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const real_2d_array &b, const boolean_2d_array &c);
} // end of namespace alglib

#endif // OnceOnly
