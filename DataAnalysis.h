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
#ifndef OnceOnlyDataAnalysis_h
#define OnceOnlyDataAnalysis_h

#include "Statistics.h"
#include "Optimization.h"

// === PCA Package ===
// Depends on: (LinAlg) SVD, EVD
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
void pcabuildbasis(RMatrix x, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector s2, RMatrix v);
void pcatruncatedsubspace(RMatrix x, ae_int_t npoints, ae_int_t nvars, ae_int_t nneeded, double eps, ae_int_t maxits, RVector s2, RMatrix v);
void pcatruncatedsubspacesparse(sparsematrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t nneeded, double eps, ae_int_t maxits, RVector s2, RMatrix v);
} // end of namespace alglib_impl

namespace alglib {
// Principal components analysis
//
// This function builds orthogonal basis  where  first  axis  corresponds  to
// direction with maximum variance, second axis  maximizes  variance  in  the
// subspace orthogonal to first axis and so on.
//
// This function builds FULL basis, i.e. returns N vectors  corresponding  to
// ALL directions, no matter how informative. If you need  just a  few  (say,
// 10 or 50) of the most important directions, you may find it faster to  use
// one of the reduced versions:
// * pcatruncatedsubspace() - for subspace iteration based method
//
// It should be noted that, unlike LDA, PCA does not use class labels.
//
// Inputs:
//     X           -   dataset, array[0..NPoints-1,0..NVars-1].
//                     matrix contains ONLY INDEPENDENT VARIABLES.
//     NPoints     -   dataset size, NPoints >= 0
//     NVars       -   number of independent variables, NVars >= 1
//
// Outputs:
//     Info        -   return code:
//                     * -4, if SVD subroutine haven't converged
//                     * -1, if wrong parameters has been passed (NPoints<0,
//                           NVars<1)
//                     *  1, if task is solved
//     S2          -   array[0..NVars-1]. variance values corresponding
//                     to basis vectors.
//     V           -   array[0..NVars-1,0..NVars-1]
//                     matrix, whose columns store basis vectors.
//
// ALGLIB: Copyright 25.08.2008 by Sergey Bochkanov
void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v);

// Principal components analysis
//
// This function performs truncated PCA, i.e. returns just a few most important
// directions.
//
// Internally it uses iterative eigensolver which is very efficient when only
// a minor fraction of full basis is required. Thus, if you need full  basis,
// it is better to use pcabuildbasis() function.
//
// It should be noted that, unlike LDA, PCA does not use class labels.
//
// Inputs:
//     X           -   dataset, array[0..NPoints-1,0..NVars-1].
//                     matrix contains ONLY INDEPENDENT VARIABLES.
//     NPoints     -   dataset size, NPoints >= 0
//     NVars       -   number of independent variables, NVars >= 1
//     NNeeded     -   number of requested components, in [1,NVars] range;
//                     this function is efficient only for NNeeded<<NVars.
//     Eps         -   desired  precision  of  vectors  returned;  underlying
//                     solver will stop iterations as soon as absolute  error
//                     in corresponding singular values  reduces  to  roughly
//                     eps*MAX(lambda[]), with lambda[] being array of  eigen
//                     values.
//                     Zero value means that  algorithm  performs  number  of
//                     iterations  specified  by  maxits  parameter,  without
//                     paying attention to precision.
//     MaxIts      -   number of iterations performed by  subspace  iteration
//                     method. Zero value means that no  limit  on  iteration
//                     count is placed (eps-based stopping condition is used).
//
//
// Outputs:
//     S2          -   array[NNeeded]. Variance values corresponding
//                     to basis vectors.
//     V           -   array[NVars,NNeeded]
//                     matrix, whose columns store basis vectors.
//
// NOTE: passing eps=0 and maxits=0 results in small eps  being  selected  as
// stopping condition. Exact value of automatically selected eps is  version-
// -dependent.
//
// ALGLIB: Copyright 10.01.2017 by Sergey Bochkanov
void pcatruncatedsubspace(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nneeded, const double eps, const ae_int_t maxits, real_1d_array &s2, real_2d_array &v);

// Sparse truncated principal components analysis
//
// This function performs sparse truncated PCA, i.e. returns just a few  most
// important principal components for a sparse input X.
//
// Internally it uses iterative eigensolver which is very efficient when only
// a minor fraction of full basis is required.
//
// It should be noted that, unlike LDA, PCA does not use class labels.
//
// Inputs:
//     X           -   sparse dataset, sparse  npoints*nvars  matrix.  It  is
//                     recommended to use CRS sparse storage format;  non-CRS
//                     input will be internally converted to CRS.
//                     Matrix contains ONLY INDEPENDENT VARIABLES,  and  must
//                     be EXACTLY npoints*nvars.
//     NPoints     -   dataset size, NPoints >= 0
//     NVars       -   number of independent variables, NVars >= 1
//     NNeeded     -   number of requested components, in [1,NVars] range;
//                     this function is efficient only for NNeeded<<NVars.
//     Eps         -   desired  precision  of  vectors  returned;  underlying
//                     solver will stop iterations as soon as absolute  error
//                     in corresponding singular values  reduces  to  roughly
//                     eps*MAX(lambda[]), with lambda[] being array of  eigen
//                     values.
//                     Zero value means that  algorithm  performs  number  of
//                     iterations  specified  by  maxits  parameter,  without
//                     paying attention to precision.
//     MaxIts      -   number of iterations performed by  subspace  iteration
//                     method. Zero value means that no  limit  on  iteration
//                     count is placed (eps-based stopping condition is used).
//
//
// Outputs:
//     S2          -   array[NNeeded]. Variance values corresponding
//                     to basis vectors.
//     V           -   array[NVars,NNeeded]
//                     matrix, whose columns store basis vectors.
//
// NOTE: passing eps=0 and maxits=0 results in small eps  being  selected  as
//       a stopping condition. Exact value of automatically selected  eps  is
//       version-dependent.
//
// NOTE: zero  MaxIts  is  silently  replaced  by some reasonable value which
//       prevents eternal loops (possible when inputs are degenerate and  too
//       stringent stopping criteria are specified). In  current  version  it
//       is 50+2*NVars.
//
// ALGLIB: Copyright 10.01.2017 by Sergey Bochkanov
void pcatruncatedsubspacesparse(const sparsematrix &x, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nneeded, const double eps, const ae_int_t maxits, real_1d_array &s2, real_2d_array &v);
} // end of namespace alglib

// === BDSS Package ===
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
} cvreport;
void cvreport_init(void *_p, bool make_automatic);
void cvreport_copy(void *_dst, void *_src, bool make_automatic);
void cvreport_free(void *_p, bool make_automatic);

void dserrallocate(ae_int_t nclasses, RVector buf);
void dserraccumulate(RVector buf, RVector y, RVector desiredy);
void dserrfinish(RVector buf);
void dsnormalize(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector means, RVector sigmas);
void dsnormalizec(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector means, RVector sigmas);
double dsgetmeanmindistance(RMatrix xy, ae_int_t npoints, ae_int_t nvars);
void dstie(RVector a, ae_int_t n, ZVector ties, ae_int_t *tiecount, ZVector p1, ZVector p2);
void dstiefasti(RVector a, ZVector b, ae_int_t n, ZVector ties, ae_int_t *tiecount, RVector bufr, ZVector bufi);
void dsoptimalsplit2(RVector a, ZVector c, ae_int_t n, ae_int_t *info, double *threshold, double *pal, double *pbl, double *par, double *pbr, double *cve);
void dsoptimalsplit2fast(RVector a, ZVector c, ZVector tiesbuf, ZVector cntbuf, RVector bufr, ZVector bufi, ae_int_t n, ae_int_t nc, double alpha, ae_int_t *info, double *threshold, double *rms, double *cvrms);
void dssplitk(RVector a, ZVector c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t *info, RVector thresholds, ae_int_t *ni, double *cve);
void dsoptimalsplitk(RVector a, ZVector c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t *info, RVector thresholds, ae_int_t *ni, double *cve);
} // end of namespace alglib_impl

namespace alglib {
// Optimal binary classification
//
// Algorithms finds optimal (=with minimal cross-entropy) binary partition.
// Internal subroutine.
//
// Inputs:
//     A       -   array[0..N-1], variable
//     C       -   array[0..N-1], class numbers (0 or 1).
//     N       -   array size
//
// Outputs:
//     Info    -   completetion code:
//                 * -3, all values of A[] are same (partition is impossible)
//                 * -2, one of C[] is incorrect (<0, >1)
//                 * -1, incorrect pararemets were passed (N <= 0).
//                 *  1, OK
//     Threshold-  partiton boundary. Left part contains values which are
//                 strictly less than Threshold. Right part contains values
//                 which are greater than or equal to Threshold.
//     PAL, PBL-   probabilities P(0|v<Threshold) and P(1|v<Threshold)
//     PAR, PBR-   probabilities P(0|v >= Threshold) and P(1|v >= Threshold)
//     CVE     -   cross-validation estimate of cross-entropy
//
// ALGLIB: Copyright 22.05.2008 by Sergey Bochkanov
void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve);

// Optimal partition, internal subroutine. Fast version.
//
// Accepts:
//     A       array[0..N-1]       array of attributes     array[0..N-1]
//     C       array[0..N-1]       array of class labels
//     TiesBuf array[0..N]         temporaries (ties)
//     CntBuf  array[0..2*NC-1]    temporaries (counts)
//     Alpha                       centering factor (0 <= alpha <= 1, recommended value - 0.05)
//     BufR    array[0..N-1]       temporaries
//     BufI    array[0..N-1]       temporaries
//
// Output:
//     Info    error code (">0"=OK, "<0"=bad)
//     RMS     training set RMS error
//     CVRMS   leave-one-out RMS error
//
// Note:
//     content of all arrays is changed by subroutine;
//     it doesn't allocate temporaries.
//
// ALGLIB: Copyright 11.12.2008 by Sergey Bochkanov
void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms);
} // end of namespace alglib

// === MLPBASE Package ===
// Depends on: (AlgLibInternal) SCODES, HPCCORES
// Depends on: (LinAlg) SPARSE
// Depends on: BDSS
namespace alglib_impl {
typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
} modelerrors;
void modelerrors_init(void *_p, bool make_automatic);
void modelerrors_copy(void *_dst, void *_src, bool make_automatic);
void modelerrors_free(void *_p, bool make_automatic);

typedef struct {
   double f;
   ae_vector g;
} smlpgrad;
void smlpgrad_init(void *_p, bool make_automatic);
void smlpgrad_copy(void *_dst, void *_src, bool make_automatic);
void smlpgrad_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t hlnetworktype;
   ae_int_t hlnormtype;
   ae_vector hllayersizes;
   ae_vector hlconnections;
   ae_vector hlneurons;
   ae_vector structinfo;
   ae_vector weights;
   ae_vector columnmeans;
   ae_vector columnsigmas;
   ae_vector neurons;
   ae_vector dfdnet;
   ae_vector derror;
   ae_vector x;
   ae_vector y;
   ae_matrix xy;
   ae_vector xyrow;
   ae_vector nwbuf;
   ae_vector integerbuf;
   modelerrors err;
   ae_vector rndbuf;
   ae_shared_pool buf;
   ae_shared_pool gradbuf;
   ae_matrix dummydxy;
   sparsematrix dummysxy;
   ae_vector dummyidx;
   ae_shared_pool dummypool;
} multilayerperceptron;
void multilayerperceptron_init(void *_p, bool make_automatic);
void multilayerperceptron_copy(void *_dst, void *_src, bool make_automatic);
void multilayerperceptron_free(void *_p, bool make_automatic);

ae_int_t mlpgradsplitcost();
ae_int_t mlpgradsplitsize();
void mlpcreate0(ae_int_t nin, ae_int_t nout, multilayerperceptron *network);
void mlpcreate1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, multilayerperceptron *network);
void mlpcreate2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, multilayerperceptron *network);
void mlpcreateb0(ae_int_t nin, ae_int_t nout, double b, double d, multilayerperceptron *network);
void mlpcreateb1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double b, double d, multilayerperceptron *network);
void mlpcreateb2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double b, double d, multilayerperceptron *network);
void mlpcreater0(ae_int_t nin, ae_int_t nout, double a, double b, multilayerperceptron *network);
void mlpcreater1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double a, double b, multilayerperceptron *network);
void mlpcreater2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double a, double b, multilayerperceptron *network);
void mlpcreatec0(ae_int_t nin, ae_int_t nout, multilayerperceptron *network);
void mlpcreatec1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, multilayerperceptron *network);
void mlpcreatec2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, multilayerperceptron *network);
void mlpcopy(multilayerperceptron *network1, multilayerperceptron *network2);
void mlpcopyshared(multilayerperceptron *network1, multilayerperceptron *network2);
bool mlpsamearchitecture(multilayerperceptron *network1, multilayerperceptron *network2);
void mlpcopytunableparameters(multilayerperceptron *network1, multilayerperceptron *network2);
void mlpexporttunableparameters(multilayerperceptron *network, RVector p, ae_int_t *pcount);
void mlpimporttunableparameters(multilayerperceptron *network, RVector p);
void mlpserializeold(multilayerperceptron *network, RVector ra, ae_int_t *rlen);
void mlpunserializeold(RVector ra, multilayerperceptron *network);
void mlprandomize(multilayerperceptron *network);
void mlprandomizefull(multilayerperceptron *network);
void mlpinitpreprocessor(multilayerperceptron *network, RMatrix xy, ae_int_t ssize);
void mlpinitpreprocessorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t ssize);
void mlpinitpreprocessorsubset(multilayerperceptron *network, RMatrix xy, ae_int_t setsize, ZVector idx, ae_int_t subsetsize);
void mlpinitpreprocessorsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector idx, ae_int_t subsetsize);
void mlpproperties(multilayerperceptron *network, ae_int_t *nin, ae_int_t *nout, ae_int_t *wcount);
ae_int_t mlpntotal(multilayerperceptron *network);
ae_int_t mlpgetinputscount(multilayerperceptron *network);
ae_int_t mlpgetoutputscount(multilayerperceptron *network);
ae_int_t mlpgetweightscount(multilayerperceptron *network);
bool mlpissoftmax(multilayerperceptron *network);
ae_int_t mlpgetlayerscount(multilayerperceptron *network);
ae_int_t mlpgetlayersize(multilayerperceptron *network, ae_int_t k);
void mlpgetinputscaling(multilayerperceptron *network, ae_int_t i, double *mean, double *sigma);
void mlpgetoutputscaling(multilayerperceptron *network, ae_int_t i, double *mean, double *sigma);
void mlpgetneuroninfo(multilayerperceptron *network, ae_int_t k, ae_int_t i, ae_int_t *fkind, double *threshold);
double mlpgetweight(multilayerperceptron *network, ae_int_t k0, ae_int_t i0, ae_int_t k1, ae_int_t i1);
void mlpsetinputscaling(multilayerperceptron *network, ae_int_t i, double mean, double sigma);
void mlpsetoutputscaling(multilayerperceptron *network, ae_int_t i, double mean, double sigma);
void mlpsetneuroninfo(multilayerperceptron *network, ae_int_t k, ae_int_t i, ae_int_t fkind, double threshold);
void mlpsetweight(multilayerperceptron *network, ae_int_t k0, ae_int_t i0, ae_int_t k1, ae_int_t i1, double w);
void mlpactivationfunction(double net, ae_int_t k, double *f, double *df, double *d2f);
void mlpprocess(multilayerperceptron *network, RVector x, RVector y);
void mlpprocessi(multilayerperceptron *network, RVector x, RVector y);
double mlperror(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlperrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlperrorn(multilayerperceptron *network, RMatrix xy, ae_int_t ssize);
ae_int_t mlpclserror(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlprelclserror(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlprelclserrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlpavgce(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlpavgcesparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlprmserror(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlprmserrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlpavgerror(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlpavgerrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlpavgrelerror(multilayerperceptron *network, RMatrix xy, ae_int_t npoints);
double mlpavgrelerrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
void mlpgrad(multilayerperceptron *network, RVector x, RVector desiredy, double *e, RVector grad);
void mlpgradn(multilayerperceptron *network, RVector x, RVector desiredy, double *e, RVector grad);
void mlpgradbatch(multilayerperceptron *network, RMatrix xy, ae_int_t ssize, double *e, RVector grad);
void mlpgradbatchsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t ssize, double *e, RVector grad);
void mlpgradbatchsubset(multilayerperceptron *network, RMatrix xy, ae_int_t setsize, ZVector idx, ae_int_t subsetsize, double *e, RVector grad);
void mlpgradbatchsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector idx, ae_int_t subsetsize, double *e, RVector grad);
void mlpgradbatchx(multilayerperceptron *network, RMatrix densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, ae_shared_pool *gradbuf);
void mlpgradnbatch(multilayerperceptron *network, RMatrix xy, ae_int_t ssize, double *e, RVector grad);
void mlphessiannbatch(multilayerperceptron *network, RMatrix xy, ae_int_t ssize, double *e, RVector grad, RMatrix h);
void mlphessianbatch(multilayerperceptron *network, RMatrix xy, ae_int_t ssize, double *e, RVector grad, RMatrix h);
void mlpinternalprocessvector(ZVector structinfo, RVector weights, RVector columnmeans, RVector columnsigmas, RVector neurons, RVector dfdnet, RVector x, RVector y);
void mlpalloc(ae_serializer *s, multilayerperceptron *network);
void mlpserialize(ae_serializer *s, multilayerperceptron *network);
void mlpunserialize(ae_serializer *s, multilayerperceptron *network);
void mlpallerrorssubset(multilayerperceptron *network, RMatrix xy, ae_int_t setsize, ZVector subset, ae_int_t subsetsize, modelerrors *rep);
void mlpallerrorssparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector subset, ae_int_t subsetsize, modelerrors *rep);
double mlperrorsubset(multilayerperceptron *network, RMatrix xy, ae_int_t setsize, ZVector subset, ae_int_t subsetsize);
double mlperrorsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector subset, ae_int_t subsetsize);
void mlpallerrorsx(multilayerperceptron *network, RMatrix densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep);
} // end of namespace alglib_impl

namespace alglib {
// Model's errors:
//     * RelCLSError   -   fraction of misclassified cases.
//     * AvgCE         -   acerage cross-entropy
//     * RMSError      -   root-mean-square error
//     * AvgError      -   average error
//     * AvgRelError   -   average relative error
//
// NOTE 1: RelCLSError/AvgCE are zero on regression problems.
//
// NOTE 2: on classification problems  RMSError/AvgError/AvgRelError  contain
//         errors in prediction of posterior probabilities
DecClass(modelerrors, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);

//
DecClass(multilayerperceptron, EndD);

// This function serializes data structure to string.
//
// Important properties of s_out:
// * it contains alphanumeric characters, dots, underscores, minus signs
// * these symbols are grouped into words, which are separated by spaces
//   and Windows-style (CR+LF) newlines
// * although  serializer  uses  spaces and CR+LF as separators, you can
//   replace any separator character by arbitrary combination of spaces,
//   tabs, Windows or Unix newlines. It allows flexible reformatting  of
//   the  string  in  case you want to include it into text or XML file.
//   But you should not insert separators into the middle of the "words"
//   nor you should change case of letters.
// * s_out can be freely moved between 32-bit and 64-bit systems, little
//   and big endian machines, and so on. You can serialize structure  on
//   32-bit machine and unserialize it on 64-bit one (or vice versa), or
//   serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also
//   serialize  it  in  C++ version of ALGLIB and unserialize in C# one,
//   and vice versa.
void mlpserialize(multilayerperceptron &obj, std::string &s_out);

// This function unserializes data structure from string.
void mlpunserialize(const std::string &s_in, multilayerperceptron &obj);

// This function serializes data structure to C++ stream.
//
// Data stream generated by this function is same as  string  representation
// generated  by  string  version  of  serializer - alphanumeric characters,
// dots, underscores, minus signs, which are grouped into words separated by
// spaces and CR+LF.
//
// We recommend you to read comments on string version of serializer to find
// out more about serialization of AlGLIB objects.
void mlpserialize(multilayerperceptron &obj, std::ostream &s_out);

// This function unserializes data structure from stream.
void mlpunserialize(const std::istream &s_in, multilayerperceptron &obj);

// Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
// layers, with linear output layer. Network weights are  filled  with  small
// random values.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);

// Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
// non-linear activation function. Output layer is linear.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);

// Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
// with non-linear activation function. Output layer is linear.
//  $ALL
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);

// Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
// layers with non-linear output layer. Network weights are filled with small
// random values.
//
// Activation function of the output layer takes values:
//
//     (B, +INF), if D >= 0
//
// or
//
//     (-INF, B), if D<0.
//
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);

// Same as MLPCreateB0 but with non-linear hidden layer.
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);

// Same as MLPCreateB0 but with two non-linear hidden layers.
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);

// Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
// layers with non-linear output layer. Network weights are filled with small
// random values. Activation function of the output layer takes values [A,B].
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);

// Same as MLPCreateR0, but with non-linear hidden layer.
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);

// Same as MLPCreateR0, but with two non-linear hidden layers.
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);

// Creates classifier network with NIn  inputs  and  NOut  possible  classes.
// Network contains no hidden layers and linear output  layer  with  SOFTMAX-
// normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
// probabilities).
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);

// Same as MLPCreateC0, but with one non-linear hidden layer.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);

// Same as MLPCreateC0, but with two non-linear hidden layers.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);

// Copying of neural network
//
// Inputs:
//     Network1 -   original
//
// Outputs:
//     Network2 -   copy
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpcopy(const multilayerperceptron &network1, multilayerperceptron &network2);

// This function copies tunable  parameters (weights/means/sigmas)  from  one
// network to another with same architecture. It  performs  some  rudimentary
// checks that architectures are same, and throws exception if check fails.
//
// It is intended for fast copying of states between two  network  which  are
// known to have same geometry.
//
// Inputs:
//     Network1 -   source, must be correctly initialized
//     Network2 -   target, must have same architecture
//
// Outputs:
//     Network2 -   network state is copied from source to target
//
// ALGLIB: Copyright 20.06.2013 by Sergey Bochkanov
void mlpcopytunableparameters(const multilayerperceptron &network1, const multilayerperceptron &network2);

// Randomization of neural network weights
//
// ALGLIB: Copyright 06.11.2007 by Sergey Bochkanov
void mlprandomize(const multilayerperceptron &network);

// Randomization of neural network weights and standartisator
//
// ALGLIB: Copyright 10.03.2008 by Sergey Bochkanov
void mlprandomizefull(const multilayerperceptron &network);

// Internal subroutine.
//
// ALGLIB: Copyright 30.03.2008 by Sergey Bochkanov
void mlpinitpreprocessor(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);

// Returns information about initialized network: number of inputs, outputs,
// weights.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount);

// Returns number of inputs.
//
// ALGLIB: Copyright 19.10.2011 by Sergey Bochkanov
ae_int_t mlpgetinputscount(const multilayerperceptron &network);

// Returns number of outputs.
//
// ALGLIB: Copyright 19.10.2011 by Sergey Bochkanov
ae_int_t mlpgetoutputscount(const multilayerperceptron &network);

// Returns number of weights.
//
// ALGLIB: Copyright 19.10.2011 by Sergey Bochkanov
ae_int_t mlpgetweightscount(const multilayerperceptron &network);

// Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
bool mlpissoftmax(const multilayerperceptron &network);

// This function returns total number of layers (including input, hidden and
// output layers).
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
ae_int_t mlpgetlayerscount(const multilayerperceptron &network);

// This function returns size of K-th layer.
//
// K=0 corresponds to input layer, K=CNT-1 corresponds to output layer.
//
// Size of the output layer is always equal to the number of outputs, although
// when we have softmax-normalized network, last neuron doesn't have any
// connections - it is just zero.
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k);

// This function returns offset/scaling coefficients for I-th input of the
// network.
//
// Inputs:
//     Network     -   network
//     I           -   input index
//
// Outputs:
//     Mean        -   mean term
//     Sigma       -   sigma term, guaranteed to be nonzero.
//
// I-th input is passed through linear transformation
//     IN[i] = (IN[i]-Mean)/Sigma
// before feeding to the network
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);

// This function returns offset/scaling coefficients for I-th output of the
// network.
//
// Inputs:
//     Network     -   network
//     I           -   input index
//
// Outputs:
//     Mean        -   mean term
//     Sigma       -   sigma term, guaranteed to be nonzero.
//
// I-th output is passed through linear transformation
//     OUT[i] = OUT[i]*Sigma+Mean
// before returning it to user. In case we have SOFTMAX-normalized network,
// we return (Mean,Sigma)=(0.0,1.0).
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);

// This function returns information about Ith neuron of Kth layer
//
// Inputs:
//     Network     -   network
//     K           -   layer index
//     I           -   neuron index (within layer)
//
// Outputs:
//     FKind       -   activation function type (used by MLPActivationFunction())
//                     this value is zero for input or linear neurons
//     Threshold   -   also called offset, bias
//                     zero for input neurons
//
// NOTE: this function throws exception if layer or neuron with  given  index
// do not exists.
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold);

// This function returns information about connection from I0-th neuron of
// K0-th layer to I1-th neuron of K1-th layer.
//
// Inputs:
//     Network     -   network
//     K0          -   layer index
//     I0          -   neuron index (within layer)
//     K1          -   layer index
//     I1          -   neuron index (within layer)
//
// Result:
//     connection weight (zero for non-existent connections)
//
// This function:
// 1. throws exception if layer or neuron with given index do not exists.
// 2. returns zero if neurons exist, but there is no connection between them
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1);

// This function sets offset/scaling coefficients for I-th input of the
// network.
//
// Inputs:
//     Network     -   network
//     I           -   input index
//     Mean        -   mean term
//     Sigma       -   sigma term (if zero, will be replaced by 1.0)
//
// NTE: I-th input is passed through linear transformation
//     IN[i] = (IN[i]-Mean)/Sigma
// before feeding to the network. This function sets Mean and Sigma.
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);

// This function sets offset/scaling coefficients for I-th output of the
// network.
//
// Inputs:
//     Network     -   network
//     I           -   input index
//     Mean        -   mean term
//     Sigma       -   sigma term (if zero, will be replaced by 1.0)
//
// Outputs:
//
// NOTE: I-th output is passed through linear transformation
//     OUT[i] = OUT[i]*Sigma+Mean
// before returning it to user. This function sets Sigma/Mean. In case we
// have SOFTMAX-normalized network, you can not set (Sigma,Mean) to anything
// other than(0.0,1.0) - this function will throw exception.
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);

// This function modifies information about Ith neuron of Kth layer
//
// Inputs:
//     Network     -   network
//     K           -   layer index
//     I           -   neuron index (within layer)
//     FKind       -   activation function type (used by MLPActivationFunction())
//                     this value must be zero for input neurons
//                     (you can not set activation function for input neurons)
//     Threshold   -   also called offset, bias
//                     this value must be zero for input neurons
//                     (you can not set threshold for input neurons)
//
// NOTES:
// 1. this function throws exception if layer or neuron with given index do
//    not exists.
// 2. this function also throws exception when you try to set non-linear
//    activation function for input neurons (any kind of network) or for output
//    neurons of classifier network.
// 3. this function throws exception when you try to set non-zero threshold for
//    input neurons (any kind of network).
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold);

// This function modifies information about connection from I0-th neuron of
// K0-th layer to I1-th neuron of K1-th layer.
//
// Inputs:
//     Network     -   network
//     K0          -   layer index
//     I0          -   neuron index (within layer)
//     K1          -   layer index
//     I1          -   neuron index (within layer)
//     W           -   connection weight (must be zero for non-existent
//                     connections)
//
// This function:
// 1. throws exception if layer or neuron with given index do not exists.
// 2. throws exception if you try to set non-zero weight for non-existent
//    connection
//
// ALGLIB: Copyright 25.03.2011 by Sergey Bochkanov
void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w);

// Neural network activation function
//
// Inputs:
//     NET         -   neuron input
//     K           -   function index (zero for linear function)
//
// Outputs:
//     F           -   function
//     DF          -   its derivative
//     D2F         -   its second derivative
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f);

// Procesing
//
// Inputs:
//     Network -   neural network
//     X       -   input vector,  array[0..NIn-1].
//
// Outputs:
//     Y       -   result. Regression estimate when solving regression  task,
//                 vector of posterior probabilities for classification task.
//
// See also MLPProcessI
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);

// 'interactive'  variant  of  MLPProcess  for  languages  like  Python which
// support constructs like "Y = MLPProcess(NN,X)" and interactive mode of the
// interpreter
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
//
// ALGLIB: Copyright 21.09.2010 by Sergey Bochkanov
void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);

// Error of the neural network on dataset.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
//     sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// Error of the neural network on dataset given by sparse matrix.
//
// Inputs:
//     Network     -   neural network
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.  Sparse  matrix  must  use  CRS  format for
//                     storage.
//     NPoints     -   points count, >= 0
//
// Result:
//     sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
double mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);

// Natural error function for neural network, internal subroutine.
//
// NOTE: this function is single-threaded. Unlike other  error  function,  it
// receives no speed-up from being executed in SMP mode.
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);

// Classification error of the neural network on dataset.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
//     classification error (number of misclassified cases)
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// Relative classification error on the test set.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
// Percent   of incorrectly   classified  cases.  Works  both  for classifier
// networks and general purpose networks used as classifiers.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 25.12.2008 by Sergey Bochkanov
double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// Relative classification error on the test set given by sparse matrix.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. Sparse matrix must use CRS format
//                     for storage.
//     NPoints     -   points count, >= 0.
//
// Result:
// Percent   of incorrectly   classified  cases.  Works  both  for classifier
// networks and general purpose networks used as classifiers.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 09.08.2012 by Sergey Bochkanov
double mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);

// Average cross-entropy  (in bits  per element) on the test set.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
// CrossEntropy/(NPoints*LN(2)).
// Zero if network solves regression task.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 08.01.2009 by Sergey Bochkanov
double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// Average  cross-entropy  (in bits  per element)  on the  test set  given by
// sparse matrix.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.  Sparse  matrix  must  use  CRS  format for
//                     storage.
//     NPoints     -   points count, >= 0.
//
// Result:
// CrossEntropy/(NPoints*LN(2)).
// Zero if network solves regression task.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 09.08.2012 by Sergey Bochkanov
double mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);

// RMS error on the test set given.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
// Root mean  square error. Its meaning for regression task is obvious. As for
// classification  task,  RMS  error  means  error  when estimating  posterior
// probabilities.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// RMS error on the test set given by sparse matrix.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.  Sparse  matrix  must  use  CRS  format for
//                     storage.
//     NPoints     -   points count, >= 0.
//
// Result:
// Root mean  square error. Its meaning for regression task is obvious. As for
// classification  task,  RMS  error  means  error  when estimating  posterior
// probabilities.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 09.08.2012 by Sergey Bochkanov
double mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);

// Average absolute error on the test set.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
// Its meaning for regression task is obvious. As for classification task, it
// means average error when estimating posterior probabilities.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 11.03.2008 by Sergey Bochkanov
double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// Average absolute error on the test set given by sparse matrix.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.  Sparse  matrix  must  use  CRS  format for
//                     storage.
//     NPoints     -   points count, >= 0.
//
// Result:
// Its meaning for regression task is obvious. As for classification task, it
// means average error when estimating posterior probabilities.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 09.08.2012 by Sergey Bochkanov
double mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);

// Average relative error on the test set.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format;
//     NPoints     -   points count.
//
// Result:
// Its meaning for regression task is obvious. As for classification task, it
// means  average  relative  error  when  estimating posterior probability of
// belonging to the correct class.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 11.03.2008 by Sergey Bochkanov
double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);

// Average relative error on the test set given by sparse matrix.
//
// Inputs:
//     Network     -   neural network;
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.  Sparse  matrix  must  use  CRS  format for
//                     storage.
//     NPoints     -   points count, >= 0.
//
// Result:
// Its meaning for regression task is obvious. As for classification task, it
// means  average  relative  error  when  estimating posterior probability of
// belonging to the correct class.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 09.08.2012 by Sergey Bochkanov
double mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);

// Gradient calculation
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     X       -   input vector, length of array must be at least NIn
//     DesiredY-   desired outputs, length of array must be at least NOut
//     Grad    -   possibly preallocated array. If size of array is smaller
//                 than WCount, it will be reallocated. It is recommended to
//                 reuse previously allocated array to reduce allocation
//                 overhead.
//
// Outputs:
//     E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
//     Grad    -   gradient of E with respect to weights of network, array[WCount]
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);

// Gradient calculation (natural error function is used)
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     X       -   input vector, length of array must be at least NIn
//     DesiredY-   desired outputs, length of array must be at least NOut
//     Grad    -   possibly preallocated array. If size of array is smaller
//                 than WCount, it will be reallocated. It is recommended to
//                 reuse previously allocated array to reduce allocation
//                 overhead.
//
// Outputs:
//     E       -   error function, sum-of-squares for regression networks,
//                 cross-entropy for classification networks.
//     Grad    -   gradient of E with respect to weights of network, array[WCount]
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);

// Batch gradient calculation for a set of inputs/outputs
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   original dataset in dense format; one sample = one row:
//                 * first NIn columns contain inputs,
//                 * for regression problem, next NOut columns store
//                   desired outputs.
//                 * for classification problem, next column (just one!)
//                   stores class number.
//     SSize   -   number of elements in XY
//     Grad    -   possibly preallocated array. If size of array is smaller
//                 than WCount, it will be reallocated. It is recommended to
//                 reuse previously allocated array to reduce allocation
//                 overhead.
//
// Outputs:
//     E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
//     Grad    -   gradient of E with respect to weights of network, array[WCount]
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);

// Batch gradient calculation for a set  of inputs/outputs  given  by  sparse
// matrices
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   original dataset in sparse format; one sample = one row:
//                 * MATRIX MUST BE STORED IN CRS FORMAT
//                 * first NIn columns contain inputs.
//                 * for regression problem, next NOut columns store
//                   desired outputs.
//                 * for classification problem, next column (just one!)
//                   stores class number.
//     SSize   -   number of elements in XY
//     Grad    -   possibly preallocated array. If size of array is smaller
//                 than WCount, it will be reallocated. It is recommended to
//                 reuse previously allocated array to reduce allocation
//                 overhead.
//
// Outputs:
//     E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
//     Grad    -   gradient of E with respect to weights of network, array[WCount]
//
// ALGLIB: Copyright 26.07.2012 by Sergey Bochkanov
void mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad);

// Batch gradient calculation for a subset of dataset
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   original dataset in dense format; one sample = one row:
//                 * first NIn columns contain inputs,
//                 * for regression problem, next NOut columns store
//                   desired outputs.
//                 * for classification problem, next column (just one!)
//                   stores class number.
//     SetSize -   real size of XY, SetSize >= 0;
//     Idx     -   subset of SubsetSize elements, array[SubsetSize]:
//                 * Idx[I] stores row index in the original dataset which is
//                   given by XY. Gradient is calculated with respect to rows
//                   whose indexes are stored in Idx[].
//                 * Idx[]  must store correct indexes; this function  throws
//                   an  exception  in  case  incorrect index (less than 0 or
//                   larger than rows(XY)) is given
//                 * Idx[]  may  store  indexes  in  any  order and even with
//                   repetitions.
//     SubsetSize- number of elements in Idx[] array:
//                 * positive value means that subset given by Idx[] is processed
//                 * zero value results in zero gradient
//                 * negative value means that full dataset is processed
//     Grad      - possibly  preallocated array. If size of array is  smaller
//                 than WCount, it will be reallocated. It is  recommended to
//                 reuse  previously  allocated  array  to  reduce allocation
//                 overhead.
//
// Outputs:
//     E         - error function, SUM(sqr(y[i]-desiredy[i])/2,i)
//     Grad      - gradient  of  E  with  respect   to  weights  of  network,
//                 array[WCount]
//
// ALGLIB: Copyright 26.07.2012 by Sergey Bochkanov
void mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);

// Batch gradient calculation for a set of inputs/outputs  for  a  subset  of
// dataset given by set of indexes.
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   original dataset in sparse format; one sample = one row:
//                 * MATRIX MUST BE STORED IN CRS FORMAT
//                 * first NIn columns contain inputs,
//                 * for regression problem, next NOut columns store
//                   desired outputs.
//                 * for classification problem, next column (just one!)
//                   stores class number.
//     SetSize -   real size of XY, SetSize >= 0;
//     Idx     -   subset of SubsetSize elements, array[SubsetSize]:
//                 * Idx[I] stores row index in the original dataset which is
//                   given by XY. Gradient is calculated with respect to rows
//                   whose indexes are stored in Idx[].
//                 * Idx[]  must store correct indexes; this function  throws
//                   an  exception  in  case  incorrect index (less than 0 or
//                   larger than rows(XY)) is given
//                 * Idx[]  may  store  indexes  in  any  order and even with
//                   repetitions.
//     SubsetSize- number of elements in Idx[] array:
//                 * positive value means that subset given by Idx[] is processed
//                 * zero value results in zero gradient
//                 * negative value means that full dataset is processed
//     Grad      - possibly  preallocated array. If size of array is  smaller
//                 than WCount, it will be reallocated. It is  recommended to
//                 reuse  previously  allocated  array  to  reduce allocation
//                 overhead.
//
// Outputs:
//     E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
//     Grad    -   gradient  of  E  with  respect   to  weights  of  network,
//                 array[WCount]
//
// NOTE: when  SubsetSize<0 is used full dataset by call MLPGradBatchSparse
//       function.
//
// ALGLIB: Copyright 26.07.2012 by Sergey Bochkanov
void mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);

// Batch gradient calculation for a set of inputs/outputs
// (natural error function is used)
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   set of inputs/outputs; one sample = one row;
//                 first NIn columns contain inputs,
//                 next NOut columns - desired outputs.
//     SSize   -   number of elements in XY
//     Grad    -   possibly preallocated array. If size of array is smaller
//                 than WCount, it will be reallocated. It is recommended to
//                 reuse previously allocated array to reduce allocation
//                 overhead.
//
// Outputs:
//     E       -   error function, sum-of-squares for regression networks,
//                 cross-entropy for classification networks.
//     Grad    -   gradient of E with respect to weights of network, array[WCount]
//
// ALGLIB: Copyright 04.11.2007 by Sergey Bochkanov
void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);

// Batch Hessian calculation (natural error function) using R-algorithm.
// Internal subroutine.
//
// ALGLIB: Copyright 26.01.2008 by Sergey Bochkanov
//
//      Hessian calculation based on R-algorithm described in
//      "Fast Exact Multiplication by the Hessian",
//      B. A. Pearlmutter,
//      Neural Computation, 1994.
void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);

// Batch Hessian calculation using R-algorithm.
// Internal subroutine.
//
// ALGLIB: Copyright 26.01.2008 by Sergey Bochkanov
//
//      Hessian calculation based on R-algorithm described in
//      "Fast Exact Multiplication by the Hessian",
//      B. A. Pearlmutter,
//      Neural Computation, 1994.
void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);

// Calculation of all types of errors on subset of dataset.
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   original dataset; one sample = one row;
//                 first NIn columns contain inputs,
//                 next NOut columns - desired outputs.
//     SetSize -   real size of XY, SetSize >= 0;
//     Subset  -   subset of SubsetSize elements, array[SubsetSize];
//     SubsetSize- number of elements in Subset[] array:
//                 * if SubsetSize>0, rows of XY with indices Subset[0]...
//                   ...Subset[SubsetSize-1] are processed
//                 * if SubsetSize=0, zeros are returned
//                 * if SubsetSize<0, entire dataset is  processed;  Subset[]
//                   array is ignored in this case.
//
// Outputs:
//     Rep     -   it contains all type of errors.
//
// ALGLIB: Copyright 04.09.2012 by Sergey Bochkanov
void mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);

// Calculation of all types of errors on subset of dataset.
//
// Inputs:
//     Network -   network initialized with one of the network creation funcs
//     XY      -   original dataset given by sparse matrix;
//                 one sample = one row;
//                 first NIn columns contain inputs,
//                 next NOut columns - desired outputs.
//     SetSize -   real size of XY, SetSize >= 0;
//     Subset  -   subset of SubsetSize elements, array[SubsetSize];
//     SubsetSize- number of elements in Subset[] array:
//                 * if SubsetSize>0, rows of XY with indices Subset[0]...
//                   ...Subset[SubsetSize-1] are processed
//                 * if SubsetSize=0, zeros are returned
//                 * if SubsetSize<0, entire dataset is  processed;  Subset[]
//                   array is ignored in this case.
//
// Outputs:
//     Rep     -   it contains all type of errors.
//
//
// ALGLIB: Copyright 04.09.2012 by Sergey Bochkanov
void mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);

// Error of the neural network on subset of dataset.
//
// Inputs:
//     Network   -     neural network;
//     XY        -     training  set,  see  below  for  information  on   the
//                     training set format;
//     SetSize   -     real size of XY, SetSize >= 0;
//     Subset    -     subset of SubsetSize elements, array[SubsetSize];
//     SubsetSize-     number of elements in Subset[] array:
//                     * if SubsetSize>0, rows of XY with indices Subset[0]...
//                       ...Subset[SubsetSize-1] are processed
//                     * if SubsetSize=0, zeros are returned
//                     * if SubsetSize<0, entire dataset is  processed;  Subset[]
//                       array is ignored in this case.
//
// Result:
//     sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 04.09.2012 by Sergey Bochkanov
double mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);

// Error of the neural network on subset of sparse dataset.
//
// Inputs:
//     Network   -     neural network;
//     XY        -     training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.  Sparse  matrix  must  use  CRS  format for
//                     storage.
//     SetSize   -     real size of XY, SetSize >= 0;
//                     it is used when SubsetSize<0;
//     Subset    -     subset of SubsetSize elements, array[SubsetSize];
//     SubsetSize-     number of elements in Subset[] array:
//                     * if SubsetSize>0, rows of XY with indices Subset[0]...
//                       ...Subset[SubsetSize-1] are processed
//                     * if SubsetSize=0, zeros are returned
//                     * if SubsetSize<0, entire dataset is  processed;  Subset[]
//                       array is ignored in this case.
//
// Result:
//     sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// dataset format is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 04.09.2012 by Sergey Bochkanov
double mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
} // end of namespace alglib

// === LDA Package ===
// Depends on: (LinAlg) MATINV, EVD
namespace alglib_impl {
void fisherlda(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, RVector w);
void fisherldan(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, RMatrix w);
} // end of namespace alglib_impl

namespace alglib {
// Multiclass Fisher LDA
//
// Subroutine finds coefficients of linear combination which optimally separates
// training set on classes.
//
// Inputs:
//     XY          -   training set, array[0..NPoints-1,0..NVars].
//                     First NVars columns store values of independent
//                     variables, next column stores number of class (from 0
//                     to NClasses-1) which dataset element belongs to. Fractional
//                     values are rounded to nearest integer.
//     NPoints     -   training set size, NPoints >= 0
//     NVars       -   number of independent variables, NVars >= 1
//     NClasses    -   number of classes, NClasses >= 2
//
//
// Outputs:
//     Info        -   return code:
//                     * -4, if internal EVD subroutine hasn't converged
//                     * -2, if there is a point with class number
//                           outside of [0..NClasses-1].
//                     * -1, if incorrect parameters was passed (NPoints<0,
//                           NVars<1, NClasses<2)
//                     *  1, if task has been solved
//                     *  2, if there was a multicollinearity in training set,
//                           but task has been solved.
//     W           -   linear combination coefficients, array[0..NVars-1]
//
// ALGLIB: Copyright 31.05.2008 by Sergey Bochkanov
void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w);

// N-dimensional multiclass Fisher LDA
//
// Subroutine finds coefficients of linear combinations which optimally separates
// training set on classes. It returns N-dimensional basis whose vector are sorted
// by quality of training set separation (in descending order).
//
// Inputs:
//     XY          -   training set, array[0..NPoints-1,0..NVars].
//                     First NVars columns store values of independent
//                     variables, next column stores number of class (from 0
//                     to NClasses-1) which dataset element belongs to. Fractional
//                     values are rounded to nearest integer.
//     NPoints     -   training set size, NPoints >= 0
//     NVars       -   number of independent variables, NVars >= 1
//     NClasses    -   number of classes, NClasses >= 2
//
//
// Outputs:
//     Info        -   return code:
//                     * -4, if internal EVD subroutine hasn't converged
//                     * -2, if there is a point with class number
//                           outside of [0..NClasses-1].
//                     * -1, if incorrect parameters was passed (NPoints<0,
//                           NVars<1, NClasses<2)
//                     *  1, if task has been solved
//                     *  2, if there was a multicollinearity in training set,
//                           but task has been solved.
//     W           -   basis, array[0..NVars-1,0..NVars-1]
//                     columns of matrix stores basis vectors, sorted by
//                     quality of training set separation (in descending order)
//
// ALGLIB: Copyright 31.05.2008 by Sergey Bochkanov
void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);
} // end of namespace alglib

// === SSA Package ===
// Depends on: (LinAlg) SVD, EVD
namespace alglib_impl {
typedef struct {
   ae_int_t nsequences;
   ae_vector sequenceidx;
   ae_vector sequencedata;
   ae_int_t algotype;
   ae_int_t windowwidth;
   ae_int_t rtpowerup;
   ae_int_t topk;
   ae_int_t precomputedwidth;
   ae_int_t precomputednbasis;
   ae_matrix precomputedbasis;
   ae_int_t defaultsubspaceits;
   ae_int_t memorylimit;
   bool arebasisandsolvervalid;
   ae_matrix basis;
   ae_matrix basist;
   ae_vector sv;
   ae_vector forecasta;
   ae_int_t nbasis;
   eigsubspacestate solver;
   ae_matrix xxt;
   hqrndstate rs;
   ae_int_t rngseed;
   ae_vector rtqueue;
   ae_int_t rtqueuecnt;
   ae_int_t rtqueuechunk;
   ae_int_t dbgcntevd;
   ae_vector tmp0;
   ae_vector tmp1;
   eigsubspacereport solverrep;
   ae_vector alongtrend;
   ae_vector alongnoise;
   ae_matrix aseqtrajectory;
   ae_matrix aseqtbproduct;
   ae_vector aseqcounts;
   ae_vector fctrend;
   ae_vector fcnoise;
   ae_matrix fctrendm;
   ae_matrix uxbatch;
   ae_int_t uxbatchwidth;
   ae_int_t uxbatchsize;
   ae_int_t uxbatchlimit;
} ssamodel;
void ssamodel_init(void *_p, bool make_automatic);
void ssamodel_copy(void *_dst, void *_src, bool make_automatic);
void ssamodel_free(void *_p, bool make_automatic);

void ssacreate(ssamodel *s);
void ssasetwindow(ssamodel *s, ae_int_t windowwidth);
void ssasetseed(ssamodel *s, ae_int_t seed);
void ssasetpoweruplength(ssamodel *s, ae_int_t pwlen);
void ssasetmemorylimit(ssamodel *s, ae_int_t memlimit);
void ssaaddsequence(ssamodel *s, RVector x, ae_int_t n);
void ssaappendpointandupdate(ssamodel *s, double x, double updateits);
void ssaappendsequenceandupdate(ssamodel *s, RVector x, ae_int_t nticks, double updateits);
void ssasetalgoprecomputed(ssamodel *s, RMatrix a, ae_int_t windowwidth, ae_int_t nbasis);
void ssasetalgotopkdirect(ssamodel *s, ae_int_t topk);
void ssasetalgotopkrealtime(ssamodel *s, ae_int_t topk);
void ssacleardata(ssamodel *s);
void ssagetbasis(ssamodel *s, RMatrix a, RVector sv, ae_int_t *windowwidth, ae_int_t *nbasis);
void ssagetlrr(ssamodel *s, RVector a, ae_int_t *windowwidth);
void ssaanalyzelastwindow(ssamodel *s, RVector trend, RVector noise, ae_int_t *nticks);
void ssaanalyzelast(ssamodel *s, ae_int_t nticks, RVector trend, RVector noise);
void ssaanalyzesequence(ssamodel *s, RVector data, ae_int_t nticks, RVector trend, RVector noise);
void ssaforecastlast(ssamodel *s, ae_int_t nticks, RVector trend);
void ssaforecastsequence(ssamodel *s, RVector data, ae_int_t datalen, ae_int_t forecastlen, bool applysmoothing, RVector trend);
void ssaforecastavglast(ssamodel *s, ae_int_t m, ae_int_t nticks, RVector trend);
void ssaforecastavgsequence(ssamodel *s, RVector data, ae_int_t datalen, ae_int_t m, ae_int_t forecastlen, bool applysmoothing, RVector trend);
} // end of namespace alglib_impl

namespace alglib {
// This object stores state of the SSA model.
//
// You should use ALGLIB functions to work with this object.
DecClass(ssamodel, EndD);

// This function creates SSA model object.  Right after creation model is  in
// "dummy" mode - you can add data,  but   analyzing/prediction  will  return
// just zeros (it assumes that basis is empty).
//
// HOW TO USE SSA MODEL:
//
// 1. create model with ssacreate()
// 2. add data with one/many ssaaddsequence() calls
// 3. choose SSA algorithm with one of ssasetalgo...() functions:
//    * ssasetalgotopkdirect() for direct one-run analysis
//    * ssasetalgotopkrealtime() for algorithm optimized for many  subsequent
//      runs with warm-start capabilities
//    * ssasetalgoprecomputed() for user-supplied basis
// 4. set window width with ssasetwindow()
// 5. perform one of the analysis-related activities:
//    a) call ssagetbasis() to get basis
//    b) call ssaanalyzelast() ssaanalyzesequence() or ssaanalyzelastwindow()
//       to perform analysis (trend/noise separation)
//    c) call  one  of   the   forecasting   functions  (ssaforecastlast() or
//       ssaforecastsequence()) to perform prediction; alternatively, you can
//       extract linear recurrence coefficients with ssagetlrr().
//    SSA analysis will be performed during first  call  to  analysis-related
//    function. SSA model is smart enough to track all changes in the dataset
//    and  model  settings,  to  cache  previously  computed  basis  and   to
//    re-evaluate basis only when necessary.
//
// Additionally, if your setting involves constant stream  of  incoming data,
// you can perform quick update already calculated  model  with  one  of  the
// incremental   append-and-update   functions:  ssaappendpointandupdate() or
// ssaappendsequenceandupdate().
//
// NOTE: steps (2), (3), (4) can be performed in arbitrary order.
//
// Inputs:
//     none
//
// Outputs:
//     S               -   structure which stores model state
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssacreate(ssamodel &s);

// This function sets window width for SSA model. You should call  it  before
// analysis phase. Default window width is 1 (not for real use).
//
// Special notes:
// * this function call can be performed at any moment before  first call  to
//   analysis-related functions
// * changing window width invalidates internally stored basis; if you change
//   window width AFTER you call analysis-related  function,  next  analysis
//   phase will require re-calculation of  the  basis  according  to  current
//   algorithm.
// * calling this function with exactly  same window width as current one has
//   no effect
// * if you specify window width larger  than any data sequence stored in the
//   model, analysis will return zero basis.
//
// Inputs:
//     S               -   SSA model created with ssacreate()
//     WindowWidth     - >= 1, new window width
//
// Outputs:
//     S               -   SSA model, updated
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssasetwindow(const ssamodel &s, const ae_int_t windowwidth);

// This  function  sets  seed  which  is used to initialize internal RNG when
// we make pseudorandom decisions on model updates.
//
// By default, deterministic seed is used - which results in same sequence of
// pseudorandom decisions every time you run SSA model. If you  specify  non-
// deterministic seed value, then SSA  model  may  return  slightly different
// results after each run.
//
// This function can be useful when you have several SSA models updated  with
// sseappendpointandupdate() called with 0<UpdateIts<1 (fractional value) and
// due to performance limitations want them to perform updates  at  different
// moments.
//
// Inputs:
//     S       -   SSA model
//     Seed    -   seed:
//                 * positive values = use deterministic seed for each run of
//                   algorithms which depend on random initialization
//                 * zero or negative values = use non-deterministic seed
//
// ALGLIB: Copyright 03.11.2017 by Sergey Bochkanov
void ssasetseed(const ssamodel &s, const ae_int_t seed);

// This function sets length of power-up cycle for real-time algorithm.
//
// By default, this algorithm performs costly O(N*WindowWidth^2)  init  phase
// followed by full run of truncated  EVD.  However,  if  you  are  ready  to
// live with a bit lower-quality basis during first few iterations,  you  can
// split this O(N*WindowWidth^2) initialization  between  several  subsequent
// append-and-update rounds. It results in better latency of the algorithm.
//
// This function invalidates basis/solver, next analysis call will result  in
// full recalculation of everything.
//
// Inputs:
//     S       -   SSA model
//     PWLen   -   length of the power-up stage:
//                 * 0 means that no power-up is requested
//                 * 1 is the same as 0
//                 * >1 means that delayed power-up is performed
//
// ALGLIB: Copyright 03.11.2017 by Sergey Bochkanov
void ssasetpoweruplength(const ssamodel &s, const ae_int_t pwlen);

// This function sets memory limit of SSA analysis.
//
// Straightforward SSA with sequence length T and window width W needs O(T*W)
// memory. It is possible to reduce memory consumption by splitting task into
// smaller chunks.
//
// Thus function allows you to specify approximate memory limit (measured  in
// double precision numbers used for buffers). Actual memory consumption will
// be comparable to the number specified by you.
//
// Default memory limit is 50.000.000 (400Mbytes) in current version.
//
// Inputs:
//     S       -   SSA model
//     MemLimit-   memory limit, >= 0. Zero value means no limit.
//
// ALGLIB: Copyright 20.12.2017 by Sergey Bochkanov
void ssasetmemorylimit(const ssamodel &s, const ae_int_t memlimit);

// This function adds data sequence to SSA  model.  Only   single-dimensional
// sequences are supported.
//
// What is a sequences? Following definitions/requirements apply:
// * a sequence  is  an  array of  values  measured  in  subsequent,  equally
//   separated time moments (ticks).
// * you may have many sequences  in your  dataset;  say,  one  sequence  may
//   correspond to one trading session.
// * sequence length should be larger  than current  window  length  (shorter
//   sequences will be ignored during analysis).
// * analysis is performed within a  sequence; different  sequences  are  NOT
//   stacked together to produce one large contiguous stream of data.
// * analysis is performed for all  sequences at once, i.e. same set of basis
//   vectors is computed for all sequences
//
// INCREMENTAL ANALYSIS
//
// This function is non intended for  incremental updates of previously found
// SSA basis. Calling it invalidates  all previous analysis results (basis is
// reset and will be recalculated from zero during next analysis).
//
// If  you  want  to  perform   incremental/real-time  SSA,  consider   using
// following functions:
// * ssaappendpointandupdate() for appending one point
// * ssaappendsequenceandupdate() for appending new sequence
//
// Inputs:
//     S               -   SSA model created with ssacreate()
//     X               -   array[N], data, can be larger (additional values
//                         are ignored)
//     N               -   data length, can be automatically determined from
//                         the array length. N >= 0.
//
// Outputs:
//     S               -   SSA model, updated
//
// NOTE: you can clear dataset with ssacleardata()
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaaddsequence(const ssamodel &s, const real_1d_array &x, const ae_int_t n);
void ssaaddsequence(const ssamodel &s, const real_1d_array &x);

// This function appends single point to last data sequence stored in the SSA
// model and tries to update model in the  incremental  manner  (if  possible
// with current algorithm).
//
// If you want to add more than one point at once:
// * if you want to add M points to the same sequence, perform M-1 calls with
//   UpdateIts parameter set to 0.0, and last call with non-zero UpdateIts.
// * if you want to add new sequence, use ssaappendsequenceandupdate()
//
// Running time of this function does NOT depend on  dataset  size,  only  on
// window width and number of singular vectors. Depending on algorithm  being
// used, incremental update has complexity:
// * for top-K real time   - O(UpdateIts*K*Width^2), with fractional UpdateIts
// * for top-K direct      - O(Width^3) for any non-zero UpdateIts
// * for precomputed basis - O(1), no update is performed
//
// Inputs:
//     S               -   SSA model created with ssacreate()
//     X               -   new point
//     UpdateIts       - >= 0,  floating  point (!)  value,  desired  update
//                         frequency:
//                         * zero value means that point is  stored,  but  no
//                           update is performed
//                         * integer part of the value means  that  specified
//                           number of iterations is always performed
//                         * fractional part of  the  value  means  that  one
//                           iteration is performed with this probability.
//
//                         Recommended value: 0<UpdateIts <= 1.  Values  larger
//                         than 1 are VERY seldom  needed.  If  your  dataset
//                         changes slowly, you can set it  to  0.1  and  skip
//                         90% of updates.
//
//                         In any case, no information is lost even with zero
//                         value of UpdateIts! It will be  incorporated  into
//                         model, sooner or later.
//
// Outputs:
//     S               -   SSA model, updated
//
// NOTE: this function uses internal  RNG  to  handle  fractional  values  of
//       UpdateIts. By default it  is  initialized  with  fixed  seed  during
//       initial calculation of basis. Thus subsequent calls to this function
//       will result in the same sequence of pseudorandom decisions.
//
//       However, if  you  have  several  SSA  models  which  are  calculated
//       simultaneously, and if you want to reduce computational  bottlenecks
//       by performing random updates at random moments, then fixed  seed  is
//       not an option - all updates will fire at same moments.
//
//       You may change it with ssasetseed() function.
//
// NOTE: this function throws an exception if called for empty dataset (there
//       is no "last" sequence to modify).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaappendpointandupdate(const ssamodel &s, const double x, const double updateits);

// This function appends new sequence to dataset stored in the SSA  model and
// tries to update model in the incremental manner (if possible  with current
// algorithm).
//
// Notes:
// * if you want to add M sequences at once, perform M-1 calls with UpdateIts
//   parameter set to 0.0, and last call with non-zero UpdateIts.
// * if you want to add just one point, use ssaappendpointandupdate()
//
// Running time of this function does NOT depend on  dataset  size,  only  on
// sequence length, window width and number of singular vectors. Depending on
// algorithm being used, incremental update has complexity:
// * for top-K real time   - O(UpdateIts*K*Width^2+(NTicks-Width)*Width^2)
// * for top-K direct      - O(Width^3+(NTicks-Width)*Width^2)
// * for precomputed basis - O(1), no update is performed
//
// Inputs:
//     S               -   SSA model created with ssacreate()
//     X               -   new sequence, array[NTicks] or larget
//     NTicks          - >= 1, number of ticks in the sequence
//     UpdateIts       - >= 0,  floating  point (!)  value,  desired  update
//                         frequency:
//                         * zero value means that point is  stored,  but  no
//                           update is performed
//                         * integer part of the value means  that  specified
//                           number of iterations is always performed
//                         * fractional part of  the  value  means  that  one
//                           iteration is performed with this probability.
//
//                         Recommended value: 0<UpdateIts <= 1.  Values  larger
//                         than 1 are VERY seldom  needed.  If  your  dataset
//                         changes slowly, you can set it  to  0.1  and  skip
//                         90% of updates.
//
//                         In any case, no information is lost even with zero
//                         value of UpdateIts! It will be  incorporated  into
//                         model, sooner or later.
//
// Outputs:
//     S               -   SSA model, updated
//
// NOTE: this function uses internal  RNG  to  handle  fractional  values  of
//       UpdateIts. By default it  is  initialized  with  fixed  seed  during
//       initial calculation of basis. Thus subsequent calls to this function
//       will result in the same sequence of pseudorandom decisions.
//
//       However, if  you  have  several  SSA  models  which  are  calculated
//       simultaneously, and if you want to reduce computational  bottlenecks
//       by performing random updates at random moments, then fixed  seed  is
//       not an option - all updates will fire at same moments.
//
//       You may change it with ssasetseed() function.
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaappendsequenceandupdate(const ssamodel &s, const real_1d_array &x, const ae_int_t nticks, const double updateits);
void ssaappendsequenceandupdate(const ssamodel &s, const real_1d_array &x, const double updateits);

// This  function sets SSA algorithm to "precomputed vectors" algorithm.
//
// This  algorithm  uses  precomputed  set  of  orthonormal  (orthogonal  AND
// normalized) basis vectors supplied by user. Thus, basis calculation  phase
// is not performed -  we  already  have  our  basis  -  and  only  analysis/
// forecasting phase requires actual calculations.
//
// This algorithm may handle "append" requests which add just  one/few  ticks
// to the end of the last sequence in O(1) time.
//
// NOTE: this algorithm accepts both basis and window  width,  because  these
//       two parameters are naturally aligned.  Calling  this  function  sets
//       window width; if you call ssasetwindow() with  other  window  width,
//       then during analysis stage algorithm will detect conflict and  reset
//       to zero basis.
//
// Inputs:
//     S               -   SSA model
//     A               -   array[WindowWidth,NBasis], orthonormalized  basis;
//                         this function does NOT control  orthogonality  and
//                         does NOT perform any kind of  renormalization.  It
//                         is your responsibility to provide it with  correct
//                         basis.
//     WindowWidth     -   window width, >= 1
//     NBasis          -   number of basis vectors, 1 <= NBasis <= WindowWidth
//
// Outputs:
//     S               -   updated model
//
// NOTE: calling this function invalidates basis in all cases.
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssasetalgoprecomputed(const ssamodel &s, const real_2d_array &a, const ae_int_t windowwidth, const ae_int_t nbasis);
void ssasetalgoprecomputed(const ssamodel &s, const real_2d_array &a);

// This  function sets SSA algorithm to "direct top-K" algorithm.
//
// "Direct top-K" algorithm performs full  SVD  of  the  N*WINDOW  trajectory
// matrix (hence its name - direct solver  is  used),  then  extracts  top  K
// components. Overall running time is O(N*WINDOW^2), where N is a number  of
// ticks in the dataset, WINDOW is window width.
//
// This algorithm may handle "append" requests which add just  one/few  ticks
// to the end of the last sequence in O(WINDOW^3) time,  which  is  ~N/WINDOW
// times faster than re-computing everything from scratch.
//
// Inputs:
//     S               -   SSA model
//     TopK            -   number of components to analyze; TopK >= 1.
//
// Outputs:
//     S               -   updated model
//
//
// NOTE: TopK>WindowWidth is silently decreased to WindowWidth during analysis
//       phase
//
// NOTE: calling this function invalidates basis, except  for  the  situation
//       when this algorithm was already set with same parameters.
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssasetalgotopkdirect(const ssamodel &s, const ae_int_t topk);

// This function sets SSA algorithm to "top-K real time algorithm". This algo
// extracts K components with largest singular values.
//
// It  is  real-time  version  of  top-K  algorithm  which  is  optimized for
// incremental processing and  fast  start-up. Internally  it  uses  subspace
// eigensolver for truncated SVD. It results  in  ability  to  perform  quick
// updates of the basis when only a few points/sequences is added to dataset.
//
// Performance profile of the algorithm is given below:
// * O(K*WindowWidth^2) running time for incremental update  of  the  dataset
//   with one of the "append-and-update" functions (ssaappendpointandupdate()
//   or ssaappendsequenceandupdate()).
// * O(N*WindowWidth^2) running time for initial basis evaluation (N=size  of
//   dataset)
// * ability  to  split  costly  initialization  across  several  incremental
//   updates of the basis (so called "Power-Up" functionality,  activated  by
//   ssasetpoweruplength() function)
//
// Inputs:
//     S               -   SSA model
//     TopK            -   number of components to analyze; TopK >= 1.
//
// Outputs:
//     S               -   updated model
//
// NOTE: this  algorithm  is  optimized  for  large-scale  tasks  with  large
//       datasets. On toy problems with just  5-10 points it can return basis
//       which is slightly different from that returned by  direct  algorithm
//       (ssasetalgotopkdirect() function). However, the  difference  becomes
//       negligible as dataset grows.
//
// NOTE: TopK>WindowWidth is silently decreased to WindowWidth during analysis
//       phase
//
// NOTE: calling this function invalidates basis, except  for  the  situation
//       when this algorithm was already set with same parameters.
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssasetalgotopkrealtime(const ssamodel &s, const ae_int_t topk);

// This function clears all data stored in the  model  and  invalidates  all
// basis components found so far.
//
// Inputs:
//     S               -   SSA model created with ssacreate()
//
// Outputs:
//     S               -   SSA model, updated
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssacleardata(const ssamodel &s);

// This function executes SSA on internally stored dataset and returns  basis
// found by current method.
//
// Inputs:
//     S               -   SSA model
//
// Outputs:
//     A               -   array[WindowWidth,NBasis],   basis;  vectors  are
//                         stored in matrix columns, by descreasing variance
//     SV              -   array[NBasis]:
//                         * zeros - for model initialized with SSASetAlgoPrecomputed()
//                         * singular values - for other algorithms
//     WindowWidth     -   current window
//     NBasis          -   basis size
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
//
// HANDLING OF DEGENERATE CASES
//
// Calling  this  function  in  degenerate  cases  (no  data  or all data are
// shorter than window size; no algorithm is specified)  returns  basis  with
// just one zero vector.
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssagetbasis(const ssamodel &s, real_2d_array &a, real_1d_array &sv, ae_int_t &windowwidth, ae_int_t &nbasis);

// This function returns linear recurrence relation (LRR) coefficients  found
// by current SSA algorithm.
//
// Inputs:
//     S               -   SSA model
//
// Outputs:
//     A               -   array[WindowWidth-1]. Coefficients  of  the
//                         linear recurrence of the form:
//                         X[W-1] = X[W-2]*A[W-2] + X[W-3]*A[W-3] + ... + X[0]*A[0].
//                         Empty array for WindowWidth=1.
//     WindowWidth     -   current window width
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
//
// HANDLING OF DEGENERATE CASES
//
// Calling  this  function  in  degenerate  cases  (no  data  or all data are
// shorter than window size; no algorithm is specified) returns zeros.
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssagetlrr(const ssamodel &s, real_1d_array &a, ae_int_t &windowwidth);

// This  function  executes  SSA  on  internally  stored  dataset and returns
// analysis  for  the  last  window  of  the  last sequence. Such analysis is
// an lightweight alternative for full scale reconstruction (see below).
//
// Typical use case for this function is  real-time  setting,  when  you  are
// interested in quick-and-dirty (very quick and very  dirty)  processing  of
// just a few last ticks of the trend.
//
// IMPORTANT: full  scale  SSA  involves  analysis  of  the  ENTIRE  dataset,
//            with reconstruction being done for  all  positions  of  sliding
//            window with subsequent hankelization  (diagonal  averaging)  of
//            the resulting matrix.
//
//            Such analysis requires O((DataLen-Window)*Window*NBasis)  FLOPs
//            and can be quite costly. However, it has  nice  noise-canceling
//            effects due to averaging.
//
//            This function performs REDUCED analysis of the last window.  It
//            is much faster - just O(Window*NBasis),  but  its  results  are
//            DIFFERENT from that of ssaanalyzelast(). In  particular,  first
//            few points of the trend are much more prone to noise.
//
// Inputs:
//     S               -   SSA model
//
// Outputs:
//     Trend           -   array[WindowSize], reconstructed trend line
//     Noise           -   array[WindowSize], the rest of the signal;
//                         it holds that ActualData = Trend+Noise.
//     NTicks          -   current WindowSize
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
// In  any  case,  only  basis  is  reused. Reconstruction is performed  from
// scratch every time you call this function.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * last sequence is shorter than the window length (analysis can  be  done,
//   but we can not perform reconstruction on the last sequence)
//
// Calling this function in degenerate cases returns following result:
// * in any case, WindowWidth ticks is returned
// * trend is assumed to be zero
// * noise is initialized by the last sequence; if last sequence  is  shorter
//   than the window size, it is moved to  the  end  of  the  array, and  the
//   beginning of the noise array is filled by zeros
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaanalyzelastwindow(const ssamodel &s, real_1d_array &trend, real_1d_array &noise, ae_int_t &nticks);

// This function:
// * builds SSA basis using internally stored (entire) dataset
// * returns reconstruction for the last NTicks of the last sequence
//
// If you want to analyze some other sequence, use ssaanalyzesequence().
//
// Reconstruction phase involves  generation  of  NTicks-WindowWidth  sliding
// windows, their decomposition using empirical orthogonal functions found by
// SSA, followed by averaging of each data point across  several  overlapping
// windows. Thus, every point in the output trend is reconstructed  using  up
// to WindowWidth overlapping  windows  (WindowWidth windows exactly  in  the
// inner points, just one window at the extremal points).
//
// IMPORTANT: due to averaging this function returns  different  results  for
//            different values of NTicks. It is expected and not a bug.
//
//            For example:
//            * Trend[NTicks-1] is always same because it is not averaged  in
//              any case (same applies to Trend[0]).
//            * Trend[NTicks-2] has different values  for  NTicks=WindowWidth
//              and NTicks=WindowWidth+1 because former  case  means that  no
//              averaging is performed, and latter  case means that averaging
//              using two sliding windows  is  performed.  Larger  values  of
//              NTicks produce same results as NTicks=WindowWidth+1.
//            * ...and so on...
//
// PERFORMANCE: this  function has O((NTicks-WindowWidth)*WindowWidth*NBasis)
//              running time. If you work  in  time-constrained  setting  and
//              have to analyze just a few last ticks, choosing NTicks  equal
//              to WindowWidth+SmoothingLen, with SmoothingLen=1...WindowWidth
//              will result in good compromise between noise cancellation and
//              analysis speed.
//
// Inputs:
//     S               -   SSA model
//     NTicks          -   number of ticks to analyze, Nticks >= 1.
//                         * special case of NTicks <= WindowWidth  is  handled
//                           by analyzing last window and  returning   NTicks
//                           last ticks.
//                         * special case NTicks>LastSequenceLen  is  handled
//                           by prepending result with NTicks-LastSequenceLen
//                           zeros.
//
// Outputs:
//     Trend           -   array[NTicks], reconstructed trend line
//     Noise           -   array[NTicks], the rest of the signal;
//                         it holds that ActualData = Trend+Noise.
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
// In  any  case,  only  basis  is  reused. Reconstruction is performed  from
// scratch every time you call this function.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * last sequence is shorter than the window length (analysis  can  be done,
//   but we can not perform reconstruction on the last sequence)
//
// Calling this function in degenerate cases returns following result:
// * in any case, NTicks ticks is returned
// * trend is assumed to be zero
// * noise is initialized by the last sequence; if last sequence  is  shorter
//   than the window size, it is moved to  the  end  of  the  array, and  the
//   beginning of the noise array is filled by zeros
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaanalyzelast(const ssamodel &s, const ae_int_t nticks, real_1d_array &trend, real_1d_array &noise);

// This function:
// * builds SSA basis using internally stored (entire) dataset
// * returns reconstruction for the sequence being passed to this function
//
// If  you  want  to  analyze  last  sequence  stored  in   the   model,  use
// ssaanalyzelast().
//
// Reconstruction phase involves  generation  of  NTicks-WindowWidth  sliding
// windows, their decomposition using empirical orthogonal functions found by
// SSA, followed by averaging of each data point across  several  overlapping
// windows. Thus, every point in the output trend is reconstructed  using  up
// to WindowWidth overlapping  windows  (WindowWidth windows exactly  in  the
// inner points, just one window at the extremal points).
//
// PERFORMANCE: this  function has O((NTicks-WindowWidth)*WindowWidth*NBasis)
//              running time. If you work  in  time-constrained  setting  and
//              have to analyze just a few last ticks, choosing NTicks  equal
//              to WindowWidth+SmoothingLen, with SmoothingLen=1...WindowWidth
//              will result in good compromise between noise cancellation and
//              analysis speed.
//
// Inputs:
//     S               -   SSA model
//     Data            -   array[NTicks], can be larger (only NTicks  leading
//                         elements will be used)
//     NTicks          -   number of ticks to analyze, Nticks >= 1.
//                         * special case of NTicks<WindowWidth  is   handled
//                           by returning zeros as trend, and signal as noise
//
// Outputs:
//     Trend           -   array[NTicks], reconstructed trend line
//     Noise           -   array[NTicks], the rest of the signal;
//                         it holds that ActualData = Trend+Noise.
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
// In  any  case,  only  basis  is  reused. Reconstruction is performed  from
// scratch every time you call this function.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * sequence being passed is shorter than the window length
//
// Calling this function in degenerate cases returns following result:
// * in any case, NTicks ticks is returned
// * trend is assumed to be zero
// * noise is initialized by the sequence.
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaanalyzesequence(const ssamodel &s, const real_1d_array &data, const ae_int_t nticks, real_1d_array &trend, real_1d_array &noise);
void ssaanalyzesequence(const ssamodel &s, const real_1d_array &data, real_1d_array &trend, real_1d_array &noise);

// This function builds SSA basis and performs forecasting  for  a  specified
// number of ticks, returning value of trend.
//
// Forecast is performed as follows:
// * SSA  trend  extraction  is  applied  to last WindowWidth elements of the
//   internally stored dataset; this step is basically a noise reduction.
// * linear recurrence relation is applied to extracted trend
//
// This function has following running time:
// * O(NBasis*WindowWidth) for trend extraction phase (always performed)
// * O(WindowWidth*NTicks) for forecast phase
//
// NOTE: noise reduction is ALWAYS applied by this algorithm; if you want  to
//       apply recurrence relation  to  raw  unprocessed  data,  use  another
//       function - ssaforecastsequence() which allows to  turn  on  and  off
//       noise reduction phase.
//
// NOTE: this algorithm performs prediction using only one - last  -  sliding
//       window.  Predictions  produced   by   such   approach   are   smooth
//       continuations of the reconstructed  trend  line,  but  they  can  be
//       easily corrupted by noise. If you need  noise-resistant  prediction,
//       use ssaforecastavglast() function, which averages predictions  built
//       using several sliding windows.
//
// Inputs:
//     S               -   SSA model
//     NTicks          -   number of ticks to forecast, NTicks >= 1
//
// Outputs:
//     Trend           -   array[NTicks], predicted trend line
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * last sequence is shorter than the WindowWidth   (analysis  can  be done,
//   but we can not perform forecasting on the last sequence)
// * window lentgh is 1 (impossible to use for forecasting)
// * SSA analysis algorithm is  configured  to  extract  basis  whose size is
//   equal to window length (impossible to use for  forecasting;  only  basis
//   whose size is less than window length can be used).
//
// Calling this function in degenerate cases returns following result:
// * NTicks  copies  of  the  last  value is returned for non-empty task with
//   large enough dataset, but with overcomplete  basis  (window  width=1  or
//   basis size is equal to window width)
// * zero trend with length=NTicks is returned for empty task
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is ever constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaforecastlast(const ssamodel &s, const ae_int_t nticks, real_1d_array &trend);

// This function builds SSA  basis  and  performs  forecasting  for  a  user-
// specified sequence, returning value of trend.
//
// Forecasting is done in two stages:
// * first,  we  extract  trend  from the WindowWidth  last  elements of  the
//   sequence. This stage is optional, you  can  turn  it  off  if  you  pass
//   data which are already processed with SSA. Of course, you  can  turn  it
//   off even for raw data, but it is not recommended - noise suppression  is
//   very important for correct prediction.
// * then, we apply LRR for last  WindowWidth-1  elements  of  the  extracted
//   trend.
//
// This function has following running time:
// * O(NBasis*WindowWidth) for trend extraction phase
// * O(WindowWidth*NTicks) for forecast phase
//
// NOTE: this algorithm performs prediction using only one - last  -  sliding
//       window.  Predictions  produced   by   such   approach   are   smooth
//       continuations of the reconstructed  trend  line,  but  they  can  be
//       easily corrupted by noise. If you need  noise-resistant  prediction,
//       use ssaforecastavgsequence() function,  which  averages  predictions
//       built using several sliding windows.
//
// Inputs:
//     S               -   SSA model
//     Data            -   array[NTicks], data to forecast
//     DataLen         -   number of ticks in the data, DataLen >= 1
//     ForecastLen     -   number of ticks to predict, ForecastLen >= 1
//     ApplySmoothing  -   whether to apply smoothing trend extraction or not;
//                         if you do not know what to specify, pass True.
//
// Outputs:
//     Trend           -   array[ForecastLen], forecasted trend
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * data sequence is shorter than the WindowWidth   (analysis  can  be done,
//   but we can not perform forecasting on the last sequence)
// * window lentgh is 1 (impossible to use for forecasting)
// * SSA analysis algorithm is  configured  to  extract  basis  whose size is
//   equal to window length (impossible to use for  forecasting;  only  basis
//   whose size is less than window length can be used).
//
// Calling this function in degenerate cases returns following result:
// * ForecastLen copies of the last value is returned for non-empty task with
//   large enough dataset, but with overcomplete  basis  (window  width=1  or
//   basis size is equal to window width)
// * zero trend with length=ForecastLen is returned for empty task
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is ever constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaforecastsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t datalen, const ae_int_t forecastlen, const bool applysmoothing, real_1d_array &trend);
void ssaforecastsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t forecastlen, real_1d_array &trend);

// This function builds SSA basis and performs forecasting  for  a  specified
// number of ticks, returning value of trend.
//
// Forecast is performed as follows:
// * SSA  trend  extraction  is  applied to last  M  sliding windows  of  the
//   internally stored dataset
// * for each of M sliding windows, M predictions are built
// * average value of M predictions is returned
//
// This function has following running time:
// * O(NBasis*WindowWidth*M) for trend extraction phase (always performed)
// * O(WindowWidth*NTicks*M) for forecast phase
//
// NOTE: noise reduction is ALWAYS applied by this algorithm; if you want  to
//       apply recurrence relation  to  raw  unprocessed  data,  use  another
//       function - ssaforecastsequence() which allows to  turn  on  and  off
//       noise reduction phase.
//
// NOTE: combination of several predictions results in lesser sensitivity  to
//       noise, but it may produce undesirable discontinuities  between  last
//       point of the trend and first point of the prediction. The reason  is
//       that  last  point  of  the  trend is usually corrupted by noise, but
//       average  value of  several  predictions  is less sensitive to noise,
//       thus discontinuity appears. It is not a bug.
//
// Inputs:
//     S               -   SSA model
//     M               -   number  of  sliding  windows  to combine, M >= 1. If
//                         your dataset has less than M sliding windows, this
//                         parameter will be silently reduced.
//     NTicks          -   number of ticks to forecast, NTicks >= 1
//
// Outputs:
//     Trend           -   array[NTicks], predicted trend line
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * last sequence is shorter than the WindowWidth   (analysis  can  be done,
//   but we can not perform forecasting on the last sequence)
// * window lentgh is 1 (impossible to use for forecasting)
// * SSA analysis algorithm is  configured  to  extract  basis  whose size is
//   equal to window length (impossible to use for  forecasting;  only  basis
//   whose size is less than window length can be used).
//
// Calling this function in degenerate cases returns following result:
// * NTicks  copies  of  the  last  value is returned for non-empty task with
//   large enough dataset, but with overcomplete  basis  (window  width=1  or
//   basis size is equal to window width)
// * zero trend with length=NTicks is returned for empty task
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is ever constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaforecastavglast(const ssamodel &s, const ae_int_t m, const ae_int_t nticks, real_1d_array &trend);

// This function builds SSA  basis  and  performs  forecasting  for  a  user-
// specified sequence, returning value of trend.
//
// Forecasting is done in two stages:
// * first,  we  extract  trend  from M last sliding windows of the sequence.
//   This stage is optional, you can  turn  it  off  if  you  pass data which
//   are already processed with SSA. Of course, you  can  turn  it  off  even
//   for raw data, but it is not recommended  -  noise  suppression  is  very
//   important for correct prediction.
// * then, we apply LRR independently for M sliding windows
// * average of M predictions is returned
//
// This function has following running time:
// * O(NBasis*WindowWidth*M) for trend extraction phase
// * O(WindowWidth*NTicks*M) for forecast phase
//
// NOTE: combination of several predictions results in lesser sensitivity  to
//       noise, but it may produce undesirable discontinuities  between  last
//       point of the trend and first point of the prediction. The reason  is
//       that  last  point  of  the  trend is usually corrupted by noise, but
//       average  value of  several  predictions  is less sensitive to noise,
//       thus discontinuity appears. It is not a bug.
//
// Inputs:
//     S               -   SSA model
//     Data            -   array[NTicks], data to forecast
//     DataLen         -   number of ticks in the data, DataLen >= 1
//     M               -   number  of  sliding  windows  to combine, M >= 1. If
//                         your dataset has less than M sliding windows, this
//                         parameter will be silently reduced.
//     ForecastLen     -   number of ticks to predict, ForecastLen >= 1
//     ApplySmoothing  -   whether to apply smoothing trend extraction or not.
//                         if you do not know what to specify, pass true.
//
// Outputs:
//     Trend           -   array[ForecastLen], forecasted trend
//
//
// CACHING/REUSE OF THE BASIS
//
// Caching/reuse of previous results is performed:
// * first call performs full run of SSA; basis is stored in the cache
// * subsequent calls reuse previously cached basis
// * if you call any function which changes model properties (window  length,
//   algorithm, dataset), internal basis will be invalidated.
// * the only calls which do NOT invalidate basis are listed below:
//   a) ssasetwindow() with same window length
//   b) ssaappendpointandupdate()
//   c) ssaappendsequenceandupdate()
//   d) ssasetalgotopk...() with exactly same K
//   Calling these functions will result in reuse of previously found basis.
//
//
// HANDLING OF DEGENERATE CASES
//
// Following degenerate cases may happen:
// * dataset is empty (no analysis can be done)
// * all sequences are shorter than the window length,no analysis can be done
// * no algorithm is specified (no analysis can be done)
// * data sequence is shorter than the WindowWidth   (analysis  can  be done,
//   but we can not perform forecasting on the last sequence)
// * window lentgh is 1 (impossible to use for forecasting)
// * SSA analysis algorithm is  configured  to  extract  basis  whose size is
//   equal to window length (impossible to use for  forecasting;  only  basis
//   whose size is less than window length can be used).
//
// Calling this function in degenerate cases returns following result:
// * ForecastLen copies of the last value is returned for non-empty task with
//   large enough dataset, but with overcomplete  basis  (window  width=1  or
//   basis size is equal to window width)
// * zero trend with length=ForecastLen is returned for empty task
//
// No analysis is performed in degenerate cases (we immediately return  dummy
// values, no basis is ever constructed).
//
// ALGLIB: Copyright 30.10.2017 by Sergey Bochkanov
void ssaforecastavgsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t datalen, const ae_int_t m, const ae_int_t forecastlen, const bool applysmoothing, real_1d_array &trend);
void ssaforecastavgsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t m, const ae_int_t forecastlen, real_1d_array &trend);
} // end of namespace alglib

// === LINREG Package ===
// Depends on: (SpecialFunctions) IGAMMAF
// Depends on: (LinAlg) SVD
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
typedef struct {
   ae_vector w;
} linearmodel;
void linearmodel_init(void *_p, bool make_automatic);
void linearmodel_copy(void *_dst, void *_src, bool make_automatic);
void linearmodel_free(void *_p, bool make_automatic);

typedef struct {
   ae_matrix c;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double cvrmserror;
   double cvavgerror;
   double cvavgrelerror;
   ae_int_t ncvdefects;
   ae_vector cvdefects;
} lrreport;
void lrreport_init(void *_p, bool make_automatic);
void lrreport_copy(void *_dst, void *_src, bool make_automatic);
void lrreport_free(void *_p, bool make_automatic);

void lrbuild(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrbuilds(RMatrix xy, RVector s, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrbuildzs(RMatrix xy, RVector s, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrbuildz(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrunpack(linearmodel *lm, RVector v, ae_int_t *nvars);
void lrpack(RVector v, ae_int_t nvars, linearmodel *lm);
double lrprocess(linearmodel *lm, RVector x);
double lrrmserror(linearmodel *lm, RMatrix xy, ae_int_t npoints);
double lravgerror(linearmodel *lm, RMatrix xy, ae_int_t npoints);
double lravgrelerror(linearmodel *lm, RMatrix xy, ae_int_t npoints);
void lrcopy(linearmodel *lm1, linearmodel *lm2);
void lrlines(RMatrix xy, RVector s, ae_int_t n, ae_int_t *info, double *a, double *b, double *vara, double *varb, double *covab, double *corrab, double *p);
void lrline(RMatrix xy, ae_int_t n, ae_int_t *info, double *a, double *b);
} // end of namespace alglib_impl

namespace alglib {
DecClass(linearmodel, EndD);

// LRReport structure contains additional information about linear model:
// * C             -   covariation matrix,  array[0..NVars,0..NVars].
//                     C[i,j] = Cov(A[i],A[j])
// * RMSError      -   root mean square error on a training set
// * AvgError      -   average error on a training set
// * AvgRelError   -   average relative error on a training set (excluding
//                     observations with zero function value).
// * CVRMSError    -   leave-one-out cross-validation estimate of
//                     generalization error. Calculated using fast algorithm
//                     with O(NVars*NPoints) complexity.
// * CVAvgError    -   cross-validation estimate of average error
// * CVAvgRelError -   cross-validation estimate of average relative error
//
// All other fields of the structure are intended for internal use and should
// not be used outside ALGLIB.
DecClass(lrreport, real_2d_array c; double &rmserror; double &avgerror; double &avgrelerror; double &cvrmserror; double &cvavgerror; double &cvavgrelerror; ae_int_t &ncvdefects; integer_1d_array cvdefects;);

// Linear regression
//
// Subroutine builds model:
//
//     Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)
//
// and model found in ALGLIB format, covariation matrix, training set  errors
// (rms,  average,  average  relative)   and  leave-one-out  cross-validation
// estimate of the generalization error. CV  estimate calculated  using  fast
// algorithm with O(NPoints*NVars) complexity.
//
// When  covariation  matrix  is  calculated  standard deviations of function
// values are assumed to be equal to RMS error on the training set.
//
// Inputs:
//     XY          -   training set, array [0..NPoints-1,0..NVars]:
//                     * NVars columns - independent variables
//                     * last column - dependent variable
//     NPoints     -   training set size, NPoints>NVars+1
//     NVars       -   number of independent variables
//
// Outputs:
//     Info        -   return code:
//                     * -255, in case of unknown internal error
//                     * -4, if internal SVD subroutine haven't converged
//                     * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
//                     *  1, if subroutine successfully finished
//     LM          -   linear model in the ALGLIB format. Use subroutines of
//                     this unit to work with the model.
//     AR          -   additional results
//
//
// ALGLIB: Copyright 02.08.2008 by Sergey Bochkanov
void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);

// Linear regression
//
// Variant of LRBuild which uses vector of standatd deviations (errors in
// function values).
//
// Inputs:
//     XY          -   training set, array [0..NPoints-1,0..NVars]:
//                     * NVars columns - independent variables
//                     * last column - dependent variable
//     S           -   standard deviations (errors in function values)
//                     array[0..NPoints-1], S[i]>0.
//     NPoints     -   training set size, NPoints>NVars+1
//     NVars       -   number of independent variables
//
// Outputs:
//     Info        -   return code:
//                     * -255, in case of unknown internal error
//                     * -4, if internal SVD subroutine haven't converged
//                     * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
//                     * -2, if S[I] <= 0
//                     *  1, if subroutine successfully finished
//     LM          -   linear model in the ALGLIB format. Use subroutines of
//                     this unit to work with the model.
//     AR          -   additional results
//
//
// ALGLIB: Copyright 02.08.2008 by Sergey Bochkanov
void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);

// Like LRBuildS, but builds model
//
//     Y = A(0)*X[0] + ... + A(N-1)*X[N-1]
//
// i.e. with zero constant term.
//
// ALGLIB: Copyright 30.10.2008 by Sergey Bochkanov
void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);

// Like LRBuild but builds model
//
//     Y = A(0)*X[0] + ... + A(N-1)*X[N-1]
//
// i.e. with zero constant term.
//
// ALGLIB: Copyright 30.10.2008 by Sergey Bochkanov
void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);

// Unpacks coefficients of linear model.
//
// Inputs:
//     LM          -   linear model in ALGLIB format
//
// Outputs:
//     V           -   coefficients, array[0..NVars]
//                     constant term (intercept) is stored in the V[NVars].
//     NVars       -   number of independent variables (one less than number
//                     of coefficients)
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars);

// "Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
// reversed).
//
// Inputs:
//     V           -   coefficients, array[0..NVars]
//     NVars       -   number of independent variables
//
// OUTPUT PAREMETERS:
//     LM          -   linear model.
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm);

// Procesing
//
// Inputs:
//     LM      -   linear model
//     X       -   input vector,  array[0..NVars-1].
//
// Result:
//     value of linear model regression estimate
//
// ALGLIB: Copyright 03.09.2008 by Sergey Bochkanov
double lrprocess(const linearmodel &lm, const real_1d_array &x);

// RMS error on the test set
//
// Inputs:
//     LM      -   linear model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     root mean square error.
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

// Average error on the test set
//
// Inputs:
//     LM      -   linear model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     average error.
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

// RMS error on the test set
//
// Inputs:
//     LM      -   linear model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     average relative error.
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === FILTERS Package ===
// Depends on: LINREG
namespace alglib_impl {
void filtersma(RVector x, ae_int_t n, ae_int_t k);
void filterema(RVector x, ae_int_t n, double alpha);
void filterlrma(RVector x, ae_int_t n, ae_int_t k);
} // end of namespace alglib_impl

namespace alglib {
// Filters: simple moving averages (unsymmetric).
//
// This filter replaces array by results of SMA(K) filter. SMA(K) is defined
// as filter which averages at most K previous points (previous - not points
// AROUND central point) - or less, in case of the first K-1 points.
//
// Inputs:
//     X           -   array[N], array to process. It can be larger than N,
//                     in this case only first N points are processed.
//     N           -   points count, N >= 0
//     K           -   K >= 1 (K can be larger than N ,  such  cases  will  be
//                     correctly handled). Window width. K=1 corresponds  to
//                     identity transformation (nothing changes).
//
// Outputs:
//     X           -   array, whose first N elements were processed with SMA(K)
//
// NOTE 1: this function uses efficient in-place  algorithm  which  does not
//         allocate temporary arrays.
//
// NOTE 2: this algorithm makes only one pass through array and uses running
//         sum  to speed-up calculation of the averages. Additional measures
//         are taken to ensure that running sum on a long sequence  of  zero
//         elements will be correctly reset to zero even in the presence  of
//         round-off error.
//
// NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
//         averages points after the current one. Only X[i], X[i-1], ... are
//         used when calculating new value of X[i]. We should also note that
//         this algorithm uses BOTH previous points and  current  one,  i.e.
//         new value of X[i] depends on BOTH previous point and X[i] itself.
//
// ALGLIB: Copyright 25.10.2011 by Sergey Bochkanov
void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filtersma(real_1d_array &x, const ae_int_t k);

// Filters: exponential moving averages.
//
// This filter replaces array by results of EMA(alpha) filter. EMA(alpha) is
// defined as filter which replaces X[] by S[]:
//     S[0] = X[0]
//     S[t] = alpha*X[t] + (1-alpha)*S[t-1]
//
// Inputs:
//     X           -   array[N], array to process. It can be larger than N,
//                     in this case only first N points are processed.
//     N           -   points count, N >= 0
//     alpha       -   0<alpha <= 1, smoothing parameter.
//
// Outputs:
//     X           -   array, whose first N elements were processed
//                     with EMA(alpha)
//
// NOTE 1: this function uses efficient in-place  algorithm  which  does not
//         allocate temporary arrays.
//
// NOTE 2: this algorithm uses BOTH previous points and  current  one,  i.e.
//         new value of X[i] depends on BOTH previous point and X[i] itself.
//
// NOTE 3: technical analytis users quite often work  with  EMA  coefficient
//         expressed in DAYS instead of fractions. If you want to  calculate
//         EMA(N), where N is a number of days, you can use alpha=2/(N+1).
//
// ALGLIB: Copyright 25.10.2011 by Sergey Bochkanov
void filterema(real_1d_array &x, const ae_int_t n, const double alpha);
void filterema(real_1d_array &x, const double alpha);

// Filters: linear regression moving averages.
//
// This filter replaces array by results of LRMA(K) filter.
//
// LRMA(K) is defined as filter which, for each data  point,  builds  linear
// regression  model  using  K  prevous  points (point itself is included in
// these K points) and calculates value of this linear model at the point in
// question.
//
// Inputs:
//     X           -   array[N], array to process. It can be larger than N,
//                     in this case only first N points are processed.
//     N           -   points count, N >= 0
//     K           -   K >= 1 (K can be larger than N ,  such  cases  will  be
//                     correctly handled). Window width. K=1 corresponds  to
//                     identity transformation (nothing changes).
//
// Outputs:
//     X           -   array, whose first N elements were processed with SMA(K)
//
// NOTE 1: this function uses efficient in-place  algorithm  which  does not
//         allocate temporary arrays.
//
// NOTE 2: this algorithm makes only one pass through array and uses running
//         sum  to speed-up calculation of the averages. Additional measures
//         are taken to ensure that running sum on a long sequence  of  zero
//         elements will be correctly reset to zero even in the presence  of
//         round-off error.
//
// NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
//         averages points after the current one. Only X[i], X[i-1], ... are
//         used when calculating new value of X[i]. We should also note that
//         this algorithm uses BOTH previous points and  current  one,  i.e.
//         new value of X[i] depends on BOTH previous point and X[i] itself.
//
// ALGLIB: Copyright 25.10.2011 by Sergey Bochkanov
void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filterlrma(real_1d_array &x, const ae_int_t k);
} // end of namespace alglib

// === LOGIT Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS
// Depends on: MLPBASE
namespace alglib_impl {
typedef struct {
   ae_vector w;
} logitmodel;
void logitmodel_init(void *_p, bool make_automatic);
void logitmodel_copy(void *_dst, void *_src, bool make_automatic);
void logitmodel_free(void *_p, bool make_automatic);

typedef struct {
   bool brackt;
   bool stage1;
   ae_int_t infoc;
   double dg;
   double dgm;
   double dginit;
   double dgtest;
   double dgx;
   double dgxm;
   double dgy;
   double dgym;
   double finit;
   double ftest1;
   double fm;
   double fx;
   double fxm;
   double fy;
   double fym;
   double stx;
   double sty;
   double stmin;
   double stmax;
   double width;
   double width1;
   double xtrapf;
} logitmcstate;
void logitmcstate_init(void *_p, bool make_automatic);
void logitmcstate_copy(void *_dst, void *_src, bool make_automatic);
void logitmcstate_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t ngrad;
   ae_int_t nhess;
} mnlreport;
void mnlreport_init(void *_p, bool make_automatic);
void mnlreport_copy(void *_dst, void *_src, bool make_automatic);
void mnlreport_free(void *_p, bool make_automatic);

void mnltrainh(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, logitmodel *lm, mnlreport *rep);
void mnlprocess(logitmodel *lm, RVector x, RVector y);
void mnlprocessi(logitmodel *lm, RVector x, RVector y);
void mnlunpack(logitmodel *lm, RMatrix a, ae_int_t *nvars, ae_int_t *nclasses);
void mnlpack(RMatrix a, ae_int_t nvars, ae_int_t nclasses, logitmodel *lm);
void mnlcopy(logitmodel *lm1, logitmodel *lm2);
double mnlavgce(logitmodel *lm, RMatrix xy, ae_int_t npoints);
double mnlrelclserror(logitmodel *lm, RMatrix xy, ae_int_t npoints);
double mnlrmserror(logitmodel *lm, RMatrix xy, ae_int_t npoints);
double mnlavgerror(logitmodel *lm, RMatrix xy, ae_int_t npoints);
double mnlavgrelerror(logitmodel *lm, RMatrix xy, ae_int_t ssize);
ae_int_t mnlclserror(logitmodel *lm, RMatrix xy, ae_int_t npoints);
} // end of namespace alglib_impl

namespace alglib {
DecClass(logitmodel, EndD);

// MNLReport structure contains information about training process:
// * NGrad     -   number of gradient calculations
// * NHess     -   number of Hessian calculations
DecClass(mnlreport, ae_int_t &ngrad; ae_int_t &nhess;);

// This subroutine trains logit model.
//
// Inputs:
//     XY          -   training set, array[0..NPoints-1,0..NVars]
//                     First NVars columns store values of independent
//                     variables, next column stores number of class (from 0
//                     to NClasses-1) which dataset element belongs to. Fractional
//                     values are rounded to nearest integer.
//     NPoints     -   training set size, NPoints >= 1
//     NVars       -   number of independent variables, NVars >= 1
//     NClasses    -   number of classes, NClasses >= 2
//
// Outputs:
//     Info        -   return code:
//                     * -2, if there is a point with class number
//                           outside of [0..NClasses-1].
//                     * -1, if incorrect parameters was passed
//                           (NPoints<NVars+2, NVars<1, NClasses<2).
//                     *  1, if task has been solved
//     LM          -   model built
//     Rep         -   training report
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep);

// Procesing
//
// Inputs:
//     LM      -   logit model, passed by non-constant reference
//                 (some fields of structure are used as temporaries
//                 when calculating model output).
//     X       -   input vector,  array[0..NVars-1].
//     Y       -   (possibly) preallocated buffer; if size of Y is less than
//                 NClasses, it will be reallocated.If it is large enough, it
//                 is NOT reallocated, so we can save some time on reallocation.
//
// Outputs:
//     Y       -   result, array[0..NClasses-1]
//                 Vector of posterior probabilities for classification task.
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);

// 'interactive'  variant  of  MNLProcess  for  languages  like  Python which
// support constructs like "Y = MNLProcess(LM,X)" and interactive mode of the
// interpreter
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);

// Unpacks coefficients of logit model. Logit model have form:
//
//     P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))
//           S(i) = exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]), when i<M-1
//         S(M-1) = 1
//
// Inputs:
//     LM          -   logit model in ALGLIB format
//
// Outputs:
//     V           -   coefficients, array[0..NClasses-2,0..NVars]
//     NVars       -   number of independent variables
//     NClasses    -   number of classes
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses);

// "Packs" coefficients and creates logit model in ALGLIB format (MNLUnpack
// reversed).
//
// Inputs:
//     A           -   model (see MNLUnpack)
//     NVars       -   number of independent variables
//     NClasses    -   number of classes
//
// Outputs:
//     LM          -   logit model.
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm);

// Average cross-entropy (in bits per element) on the test set
//
// Inputs:
//     LM      -   logit model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     CrossEntropy/(NPoints*ln(2)).
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

// Relative classification error on the test set
//
// Inputs:
//     LM      -   logit model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     percent of incorrectly classified cases.
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

// RMS error on the test set
//
// Inputs:
//     LM      -   logit model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     root mean square error (error when estimating posterior probabilities).
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

// Average error on the test set
//
// Inputs:
//     LM      -   logit model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     average error (error when estimating posterior probabilities).
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

// Average relative error on the test set
//
// Inputs:
//     LM      -   logit model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     average relative error (error when estimating posterior probabilities).
//
// ALGLIB: Copyright 30.08.2008 by Sergey Bochkanov
double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize);

// Classification error on test set = MNLRelClsError*NPoints
//
// ALGLIB: Copyright 10.09.2008 by Sergey Bochkanov
ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === MCPD Package ===
// Depends on: (Optimization) MINBLEIC
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_vector states;
   ae_int_t npairs;
   ae_matrix data;
   ae_matrix ec;
   ae_matrix bndl;
   ae_matrix bndu;
   ae_matrix c;
   ae_vector ct;
   ae_int_t ccnt;
   ae_vector pw;
   ae_matrix priorp;
   double regterm;
   minbleicstate bs;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repnfev;
   ae_int_t repterminationtype;
   minbleicreport br;
   ae_vector tmpp;
   ae_vector effectivew;
   ae_vector effectivebndl;
   ae_vector effectivebndu;
   ae_matrix effectivec;
   ae_vector effectivect;
   ae_vector h;
   ae_matrix p;
} mcpdstate;
void mcpdstate_init(void *_p, bool make_automatic);
void mcpdstate_copy(void *_dst, void *_src, bool make_automatic);
void mcpdstate_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t inneriterationscount;
   ae_int_t outeriterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
} mcpdreport;
void mcpdreport_init(void *_p, bool make_automatic);
void mcpdreport_copy(void *_dst, void *_src, bool make_automatic);
void mcpdreport_free(void *_p, bool make_automatic);

void mcpdcreate(ae_int_t n, mcpdstate *s);
void mcpdcreateentry(ae_int_t n, ae_int_t entrystate, mcpdstate *s);
void mcpdcreateexit(ae_int_t n, ae_int_t exitstate, mcpdstate *s);
void mcpdcreateentryexit(ae_int_t n, ae_int_t entrystate, ae_int_t exitstate, mcpdstate *s);
void mcpdaddtrack(mcpdstate *s, RMatrix xy, ae_int_t k);
void mcpdsetec(mcpdstate *s, RMatrix ec);
void mcpdaddec(mcpdstate *s, ae_int_t i, ae_int_t j, double c);
void mcpdsetbc(mcpdstate *s, RMatrix bndl, RMatrix bndu);
void mcpdaddbc(mcpdstate *s, ae_int_t i, ae_int_t j, double bndl, double bndu);
void mcpdsetlc(mcpdstate *s, RMatrix c, ZVector ct, ae_int_t k);
void mcpdsettikhonovregularizer(mcpdstate *s, double v);
void mcpdsetprior(mcpdstate *s, RMatrix pp);
void mcpdsetpredictionweights(mcpdstate *s, RVector pw);
void mcpdsolve(mcpdstate *s);
void mcpdresults(mcpdstate *s, RMatrix p, mcpdreport *rep);
} // end of namespace alglib_impl

namespace alglib {
// This structure is a MCPD (Markov Chains for Population Data) solver.
//
// You should use ALGLIB functions in order to work with this object.
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
DecClass(mcpdstate, EndD);

// This structure is a MCPD training report:
//     InnerIterationsCount    -   number of inner iterations of the
//                                 underlying optimization algorithm
//     OuterIterationsCount    -   number of outer iterations of the
//                                 underlying optimization algorithm
//     NFEV                    -   number of merit function evaluations
//     TerminationType         -   termination type
//                                 (same as for MinBLEIC optimizer, positive
//                                 values denote success, negative ones -
//                                 failure)
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
DecClass(mcpdreport, ae_int_t &inneriterationscount; ae_int_t &outeriterationscount; ae_int_t &nfev; ae_int_t &terminationtype;);

// DESCRIPTION:
//
// This function creates MCPD (Markov Chains for Population Data) solver.
//
// This  solver  can  be  used  to find transition matrix P for N-dimensional
// prediction  problem  where transition from X[i] to X[i+1] is  modelled  as
//     X[i+1] = P*X[i]
// where X[i] and X[i+1] are N-dimensional population vectors (components  of
// each X are non-negative), and P is a N*N transition matrix (elements of  P
// are non-negative, each column sums to 1.0).
//
// Such models arise when when:
// * there is some population of individuals
// * individuals can have different states
// * individuals can transit from one state to another
// * population size is constant, i.e. there is no new individuals and no one
//   leaves population
// * you want to model transitions of individuals from one state into another
//
// USAGE:
//
// Here we give very brief outline of the MCPD. We strongly recommend you  to
// read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
// on data analysis which is available at http://www.alglib.net/dataanalysis/
//
// 1. User initializes algorithm state with MCPDCreate() call
//
// 2. User  adds  one  or  more  tracks -  sequences of states which describe
//    evolution of a system being modelled from different starting conditions
//
// 3. User may add optional boundary, equality  and/or  linear constraints on
//    the coefficients of P by calling one of the following functions:
//    * MCPDSetEC() to set equality constraints
//    * MCPDSetBC() to set bound constraints
//    * MCPDSetLC() to set linear constraints
//
// 4. Optionally,  user  may  set  custom  weights  for prediction errors (by
//    default, algorithm assigns non-equal, automatically chosen weights  for
//    errors in the prediction of different components of X). It can be  done
//    with a call of MCPDSetPredictionWeights() function.
//
// 5. User calls MCPDSolve() function which takes algorithm  state and
//    pointer (delegate, etc.) to callback function which calculates F/G.
//
// 6. User calls MCPDResults() to get solution
//
// Inputs:
//     N       -   problem dimension, N >= 1
//
// Outputs:
//     State   -   structure stores algorithm state
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdcreate(const ae_int_t n, mcpdstate &s);

// DESCRIPTION:
//
// This function is a specialized version of MCPDCreate()  function,  and  we
// recommend  you  to read comments for this function for general information
// about MCPD solver.
//
// This  function  creates  MCPD (Markov Chains for Population  Data)  solver
// for "Entry-state" model,  i.e. model  where transition from X[i] to X[i+1]
// is modelled as
//     X[i+1] = P*X[i]
// where
//     X[i] and X[i+1] are N-dimensional state vectors
//     P is a N*N transition matrix
// and  one  selected component of X[] is called "entry" state and is treated
// in a special way:
//     system state always transits from "entry" state to some another state
//     system state can not transit from any state into "entry" state
// Such conditions basically mean that row of P which corresponds to  "entry"
// state is zero.
//
// Such models arise when:
// * there is some population of individuals
// * individuals can have different states
// * individuals can transit from one state to another
// * population size is NOT constant -  at every moment of time there is some
//   (unpredictable) amount of "new" individuals, which can transit into  one
//   of the states at the next turn, but still no one leaves population
// * you want to model transitions of individuals from one state into another
// * but you do NOT want to predict amount of "new"  individuals  because  it
//   does not depends on individuals already present (hence  system  can  not
//   transit INTO entry state - it can only transit FROM it).
//
// This model is discussed  in  more  details  in  the ALGLIB User Guide (see
// http://www.alglib.net/dataanalysis/ for more data).
//
// Inputs:
//     N       -   problem dimension, N >= 2
//     EntryState- index of entry state, in 0..N-1
//
// Outputs:
//     State   -   structure stores algorithm state
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s);

// DESCRIPTION:
//
// This function is a specialized version of MCPDCreate()  function,  and  we
// recommend  you  to read comments for this function for general information
// about MCPD solver.
//
// This  function  creates  MCPD (Markov Chains for Population  Data)  solver
// for "Exit-state" model,  i.e. model  where  transition from X[i] to X[i+1]
// is modelled as
//     X[i+1] = P*X[i]
// where
//     X[i] and X[i+1] are N-dimensional state vectors
//     P is a N*N transition matrix
// and  one  selected component of X[] is called "exit"  state and is treated
// in a special way:
//     system state can transit from any state into "exit" state
//     system state can not transit from "exit" state into any other state
//     transition operator discards "exit" state (makes it zero at each turn)
// Such  conditions  basically  mean  that  column  of P which corresponds to
// "exit" state is zero. Multiplication by such P may decrease sum of  vector
// components.
//
// Such models arise when:
// * there is some population of individuals
// * individuals can have different states
// * individuals can transit from one state to another
// * population size is NOT constant - individuals can move into "exit" state
//   and leave population at the next turn, but there are no new individuals
// * amount of individuals which leave population can be predicted
// * you want to model transitions of individuals from one state into another
//   (including transitions into the "exit" state)
//
// This model is discussed  in  more  details  in  the ALGLIB User Guide (see
// http://www.alglib.net/dataanalysis/ for more data).
//
// Inputs:
//     N       -   problem dimension, N >= 2
//     ExitState-  index of exit state, in 0..N-1
//
// Outputs:
//     State   -   structure stores algorithm state
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s);

// DESCRIPTION:
//
// This function is a specialized version of MCPDCreate()  function,  and  we
// recommend  you  to read comments for this function for general information
// about MCPD solver.
//
// This  function  creates  MCPD (Markov Chains for Population  Data)  solver
// for "Entry-Exit-states" model, i.e. model where  transition  from  X[i] to
// X[i+1] is modelled as
//     X[i+1] = P*X[i]
// where
//     X[i] and X[i+1] are N-dimensional state vectors
//     P is a N*N transition matrix
// one selected component of X[] is called "entry" state and is treated in  a
// special way:
//     system state always transits from "entry" state to some another state
//     system state can not transit from any state into "entry" state
// and another one component of X[] is called "exit" state and is treated  in
// a special way too:
//     system state can transit from any state into "exit" state
//     system state can not transit from "exit" state into any other state
//     transition operator discards "exit" state (makes it zero at each turn)
// Such conditions basically mean that:
//     row of P which corresponds to "entry" state is zero
//     column of P which corresponds to "exit" state is zero
// Multiplication by such P may decrease sum of vector components.
//
// Such models arise when:
// * there is some population of individuals
// * individuals can have different states
// * individuals can transit from one state to another
// * population size is NOT constant
// * at every moment of time there is some (unpredictable)  amount  of  "new"
//   individuals, which can transit into one of the states at the next turn
// * some  individuals  can  move  (predictably)  into "exit" state and leave
//   population at the next turn
// * you want to model transitions of individuals from one state into another,
//   including transitions from the "entry" state and into the "exit" state.
// * but you do NOT want to predict amount of "new"  individuals  because  it
//   does not depends on individuals already present (hence  system  can  not
//   transit INTO entry state - it can only transit FROM it).
//
// This model is discussed  in  more  details  in  the ALGLIB User Guide (see
// http://www.alglib.net/dataanalysis/ for more data).
//
// Inputs:
//     N       -   problem dimension, N >= 2
//     EntryState- index of entry state, in 0..N-1
//     ExitState-  index of exit state, in 0..N-1
//
// Outputs:
//     State   -   structure stores algorithm state
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s);

// This  function  is  used to add a track - sequence of system states at the
// different moments of its evolution.
//
// You  may  add  one  or several tracks to the MCPD solver. In case you have
// several tracks, they won't overwrite each other. For example,  if you pass
// two tracks, A1-A2-A3 (system at t=A+1, t=A+2 and t=A+3) and B1-B2-B3, then
// solver will try to model transitions from t=A+1 to t=A+2, t=A+2 to  t=A+3,
// t=B+1 to t=B+2, t=B+2 to t=B+3. But it WONT mix these two tracks - i.e. it
// wont try to model transition from t=A+3 to t=B+1.
//
// Inputs:
//     S       -   solver
//     XY      -   track, array[K,N]:
//                 * I-th row is a state at t=I
//                 * elements of XY must be non-negative (exception will be
//                   thrown on negative elements)
//     K       -   number of points in a track
//                 * if given, only leading K rows of XY are used
//                 * if not given, automatically determined from size of XY
//
// NOTES:
//
// 1. Track may contain either proportional or population data:
//    * with proportional data all rows of XY must sum to 1.0, i.e. we have
//      proportions instead of absolute population values
//    * with population data rows of XY contain population counts and generally
//      do not sum to 1.0 (although they still must be non-negative)
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy);

// This function is used to add equality constraints on the elements  of  the
// transition matrix P.
//
// MCPD solver has four types of constraints which can be placed on P:
// * user-specified equality constraints (optional)
// * user-specified bound constraints (optional)
// * user-specified general linear constraints (optional)
// * basic constraints (always present):
//   * non-negativity: P[i,j] >= 0
//   * consistency: every column of P sums to 1.0
//
// Final  constraints  which  are  passed  to  the  underlying  optimizer are
// calculated  as  intersection  of all present constraints. For example, you
// may specify boundary constraint on P[0,0] and equality one:
//     0.1 <= P[0,0] <= 0.9
//     P[0,0]=0.5
// Such  combination  of  constraints  will  be  silently  reduced  to  their
// intersection, which is P[0,0]=0.5.
//
// This  function  can  be  used  to  place equality constraints on arbitrary
// subset of elements of P. Set of constraints is specified by EC, which  may
// contain either NAN's or finite numbers from [0,1]. NAN denotes absence  of
// constraint, finite number denotes equality constraint on specific  element
// of P.
//
// You can also  use  MCPDAddEC()  function  which  allows  to  ADD  equality
// constraint  for  one  element  of P without changing constraints for other
// elements.
//
// These functions (MCPDSetEC and MCPDAddEC) interact as follows:
// * there is internal matrix of equality constraints which is stored in  the
//   MCPD solver
// * MCPDSetEC() replaces this matrix by another one (SET)
// * MCPDAddEC() modifies one element of this matrix and  leaves  other  ones
//   unchanged (ADD)
// * thus  MCPDAddEC()  call  preserves  all  modifications  done by previous
//   calls,  while  MCPDSetEC()  completely discards all changes  done to the
//   equality constraints.
//
// Inputs:
//     S       -   solver
//     EC      -   equality constraints, array[N,N]. Elements of  EC  can  be
//                 either NAN's or finite  numbers from  [0,1].  NAN  denotes
//                 absence  of  constraints,  while  finite  value    denotes
//                 equality constraint on the corresponding element of P.
//
// NOTES:
//
// 1. infinite values of EC will lead to exception being thrown. Values  less
// than 0.0 or greater than 1.0 will lead to error code being returned  after
// call to MCPDSolve().
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsetec(const mcpdstate &s, const real_2d_array &ec);

// This function is used to add equality constraints on the elements  of  the
// transition matrix P.
//
// MCPD solver has four types of constraints which can be placed on P:
// * user-specified equality constraints (optional)
// * user-specified bound constraints (optional)
// * user-specified general linear constraints (optional)
// * basic constraints (always present):
//   * non-negativity: P[i,j] >= 0
//   * consistency: every column of P sums to 1.0
//
// Final  constraints  which  are  passed  to  the  underlying  optimizer are
// calculated  as  intersection  of all present constraints. For example, you
// may specify boundary constraint on P[0,0] and equality one:
//     0.1 <= P[0,0] <= 0.9
//     P[0,0]=0.5
// Such  combination  of  constraints  will  be  silently  reduced  to  their
// intersection, which is P[0,0]=0.5.
//
// This function can be used to ADD equality constraint for one element of  P
// without changing constraints for other elements.
//
// You  can  also  use  MCPDSetEC()  function  which  allows  you  to specify
// arbitrary set of equality constraints in one call.
//
// These functions (MCPDSetEC and MCPDAddEC) interact as follows:
// * there is internal matrix of equality constraints which is stored in the
//   MCPD solver
// * MCPDSetEC() replaces this matrix by another one (SET)
// * MCPDAddEC() modifies one element of this matrix and leaves  other  ones
//   unchanged (ADD)
// * thus  MCPDAddEC()  call  preserves  all  modifications done by previous
//   calls,  while  MCPDSetEC()  completely discards all changes done to the
//   equality constraints.
//
// Inputs:
//     S       -   solver
//     I       -   row index of element being constrained
//     J       -   column index of element being constrained
//     C       -   value (constraint for P[I,J]).  Can  be  either  NAN  (no
//                 constraint) or finite value from [0,1].
//
// NOTES:
//
// 1. infinite values of C  will lead to exception being thrown. Values  less
// than 0.0 or greater than 1.0 will lead to error code being returned  after
// call to MCPDSolve().
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c);

// This function is used to add bound constraints  on  the  elements  of  the
// transition matrix P.
//
// MCPD solver has four types of constraints which can be placed on P:
// * user-specified equality constraints (optional)
// * user-specified bound constraints (optional)
// * user-specified general linear constraints (optional)
// * basic constraints (always present):
//   * non-negativity: P[i,j] >= 0
//   * consistency: every column of P sums to 1.0
//
// Final  constraints  which  are  passed  to  the  underlying  optimizer are
// calculated  as  intersection  of all present constraints. For example, you
// may specify boundary constraint on P[0,0] and equality one:
//     0.1 <= P[0,0] <= 0.9
//     P[0,0]=0.5
// Such  combination  of  constraints  will  be  silently  reduced  to  their
// intersection, which is P[0,0]=0.5.
//
// This  function  can  be  used  to  place bound   constraints  on arbitrary
// subset  of  elements  of  P.  Set of constraints is specified by BndL/BndU
// matrices, which may contain arbitrary combination  of  finite  numbers  or
// infinities (like -INF<x <= 0.5 or 0.1 <= x<+INF).
//
// You can also use MCPDAddBC() function which allows to ADD bound constraint
// for one element of P without changing constraints for other elements.
//
// These functions (MCPDSetBC and MCPDAddBC) interact as follows:
// * there is internal matrix of bound constraints which is stored in the
//   MCPD solver
// * MCPDSetBC() replaces this matrix by another one (SET)
// * MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
//   unchanged (ADD)
// * thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
//   calls,  while  MCPDSetBC()  completely discards all changes  done to the
//   equality constraints.
//
// Inputs:
//     S       -   solver
//     BndL    -   lower bounds constraints, array[N,N]. Elements of BndL can
//                 be finite numbers or -INF.
//     BndU    -   upper bounds constraints, array[N,N]. Elements of BndU can
//                 be finite numbers or +INF.
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu);

// This function is used to add bound constraints  on  the  elements  of  the
// transition matrix P.
//
// MCPD solver has four types of constraints which can be placed on P:
// * user-specified equality constraints (optional)
// * user-specified bound constraints (optional)
// * user-specified general linear constraints (optional)
// * basic constraints (always present):
//   * non-negativity: P[i,j] >= 0
//   * consistency: every column of P sums to 1.0
//
// Final  constraints  which  are  passed  to  the  underlying  optimizer are
// calculated  as  intersection  of all present constraints. For example, you
// may specify boundary constraint on P[0,0] and equality one:
//     0.1 <= P[0,0] <= 0.9
//     P[0,0]=0.5
// Such  combination  of  constraints  will  be  silently  reduced  to  their
// intersection, which is P[0,0]=0.5.
//
// This  function  can  be  used to ADD bound constraint for one element of P
// without changing constraints for other elements.
//
// You  can  also  use  MCPDSetBC()  function  which  allows to  place  bound
// constraints  on arbitrary subset of elements of P.   Set of constraints is
// specified  by  BndL/BndU matrices, which may contain arbitrary combination
// of finite numbers or infinities (like -INF<x <= 0.5 or 0.1 <= x<+INF).
//
// These functions (MCPDSetBC and MCPDAddBC) interact as follows:
// * there is internal matrix of bound constraints which is stored in the
//   MCPD solver
// * MCPDSetBC() replaces this matrix by another one (SET)
// * MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
//   unchanged (ADD)
// * thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
//   calls,  while  MCPDSetBC()  completely discards all changes  done to the
//   equality constraints.
//
// Inputs:
//     S       -   solver
//     I       -   row index of element being constrained
//     J       -   column index of element being constrained
//     BndL    -   lower bound
//     BndU    -   upper bound
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu);

// This function is used to set linear equality/inequality constraints on the
// elements of the transition matrix P.
//
// This function can be used to set one or several general linear constraints
// on the elements of P. Two types of constraints are supported:
// * equality constraints
// * inequality constraints (both less-or-equal and greater-or-equal)
//
// Coefficients  of  constraints  are  specified  by  matrix  C (one  of  the
// parameters).  One  row  of  C  corresponds  to  one  constraint.   Because
// transition  matrix P has N*N elements,  we  need  N*N columns to store all
// coefficients  (they  are  stored row by row), and one more column to store
// right part - hence C has N*N+1 columns.  Constraint  kind is stored in the
// CT array.
//
// Thus, I-th linear constraint is
//     P[0,0]*C[I,0] + P[0,1]*C[I,1] + .. + P[0,N-1]*C[I,N-1] +
//         + P[1,0]*C[I,N] + P[1,1]*C[I,N+1] + ... +
//         + P[N-1,N-1]*C[I,N*N-1]  ?=?  C[I,N*N]
// where ?=? can be either "=" (CT[i]=0), "<=" (CT[i]<0) or ">=" (CT[i]>0).
//
// Your constraint may involve only some subset of P (less than N*N elements).
// For example it can be something like
//     P[0,0] + P[0,1] = 0.5
// In this case you still should pass matrix  with N*N+1 columns, but all its
// elements (except for C[0,0], C[0,1] and C[0,N*N-1]) will be zero.
//
// Inputs:
//     S       -   solver
//     C       -   array[K,N*N+1] - coefficients of constraints
//                 (see above for complete description)
//     CT      -   array[K] - constraint types
//                 (see above for complete description)
//     K       -   number of equality/inequality constraints, K >= 0:
//                 * if given, only leading K elements of C/CT are used
//                 * if not given, automatically determined from sizes of C/CT
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct);

// This function allows to  tune  amount  of  Tikhonov  regularization  being
// applied to your problem.
//
// By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
// small non-zero value,  P is transition matrix, prior_P is identity matrix,
// ||X||^2 is a sum of squared elements of X.
//
// This  function  allows  you to change coefficient r. You can  also  change
// prior values with MCPDSetPrior() function.
//
// Inputs:
//     S       -   solver
//     V       -   regularization  coefficient, finite non-negative value. It
//                 is  not  recommended  to specify zero value unless you are
//                 pretty sure that you want it.
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsettikhonovregularizer(const mcpdstate &s, const double v);

// This  function  allows to set prior values used for regularization of your
// problem.
//
// By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
// small non-zero value,  P is transition matrix, prior_P is identity matrix,
// ||X||^2 is a sum of squared elements of X.
//
// This  function  allows  you to change prior values prior_P. You  can  also
// change r with MCPDSetTikhonovRegularizer() function.
//
// Inputs:
//     S       -   solver
//     PP      -   array[N,N], matrix of prior values:
//                 1. elements must be real numbers from [0,1]
//                 2. columns must sum to 1.0.
//                 First property is checked (exception is thrown otherwise),
//                 while second one is not checked/enforced.
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp);

// This function is used to change prediction weights
//
// MCPD solver scales prediction errors as follows
//     Error(P) = ||W*(y-P*x)||^2
// where
//     x is a system state at time t
//     y is a system state at time t+1
//     P is a transition matrix
//     W is a diagonal scaling matrix
//
// By default, weights are chosen in order  to  minimize  relative prediction
// error instead of absolute one. For example, if one component of  state  is
// about 0.5 in magnitude and another one is about 0.05, then algorithm  will
// make corresponding weights equal to 2.0 and 20.0.
//
// Inputs:
//     S       -   solver
//     PW      -   array[N], weights:
//                 * must be non-negative values (exception will be thrown otherwise)
//                 * zero values will be replaced by automatically chosen values
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw);

// This function is used to start solution of the MCPD problem.
//
// After return from this function, you can use MCPDResults() to get solution
// and completion code.
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdsolve(const mcpdstate &s);

// MCPD results
//
// Inputs:
//     State   -   algorithm state
//
// Outputs:
//     P       -   array[N,N], transition matrix
//     Rep     -   optimization report. You should check Rep.TerminationType
//                 in  order  to  distinguish  successful  termination  from
//                 unsuccessful one. Speaking short, positive values  denote
//                 success, negative ones are failures.
//                 More information about fields of this  structure  can  be
//                 found in the comments on MCPDReport datatype.
//
//
// ALGLIB: Copyright 23.05.2010 by Sergey Bochkanov
void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep);
} // end of namespace alglib

// === MLPE Package ===
// Depends on: MLPBASE
namespace alglib_impl {
typedef struct {
   ae_int_t ensemblesize;
   ae_vector weights;
   ae_vector columnmeans;
   ae_vector columnsigmas;
   multilayerperceptron network;
   ae_vector y;
} mlpensemble;
void mlpensemble_init(void *_p, bool make_automatic);
void mlpensemble_copy(void *_dst, void *_src, bool make_automatic);
void mlpensemble_free(void *_p, bool make_automatic);

void mlpecreate0(ae_int_t nin, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreate1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreate2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreateb0(ae_int_t nin, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreateb1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreateb2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreater0(ae_int_t nin, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreater1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreater2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreatec0(ae_int_t nin, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreatec1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreatec2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecreatefromnetwork(multilayerperceptron *network, ae_int_t ensemblesize, mlpensemble *ensemble);
void mlpecopy(mlpensemble *ensemble1, mlpensemble *ensemble2);
void mlperandomize(mlpensemble *ensemble);
void mlpeproperties(mlpensemble *ensemble, ae_int_t *nin, ae_int_t *nout);
bool mlpeissoftmax(mlpensemble *ensemble);
void mlpeprocess(mlpensemble *ensemble, RVector x, RVector y);
void mlpeprocessi(mlpensemble *ensemble, RVector x, RVector y);
void mlpeallerrorsx(mlpensemble *ensemble, RMatrix densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep);
void mlpeallerrorssparse(mlpensemble *ensemble, sparsematrix *xy, ae_int_t npoints, double *relcls, double *avgce, double *rms, double *avg, double *avgrel);
double mlperelclserror(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints);
double mlpeavgce(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints);
double mlpermserror(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints);
double mlpeavgerror(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints);
double mlpeavgrelerror(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints);
void mlpealloc(ae_serializer *s, mlpensemble *ensemble);
void mlpeserialize(ae_serializer *s, mlpensemble *ensemble);
void mlpeunserialize(ae_serializer *s, mlpensemble *ensemble);
} // end of namespace alglib_impl

namespace alglib {
// Neural networks ensemble
DecClass(mlpensemble, EndD);

// This function serializes data structure to string.
//
// Important properties of s_out:
// * it contains alphanumeric characters, dots, underscores, minus signs
// * these symbols are grouped into words, which are separated by spaces
//   and Windows-style (CR+LF) newlines
// * although  serializer  uses  spaces and CR+LF as separators, you can
//   replace any separator character by arbitrary combination of spaces,
//   tabs, Windows or Unix newlines. It allows flexible reformatting  of
//   the  string  in  case you want to include it into text or XML file.
//   But you should not insert separators into the middle of the "words"
//   nor you should change case of letters.
// * s_out can be freely moved between 32-bit and 64-bit systems, little
//   and big endian machines, and so on. You can serialize structure  on
//   32-bit machine and unserialize it on 64-bit one (or vice versa), or
//   serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also
//   serialize  it  in  C++ version of ALGLIB and unserialize in C# one,
//   and vice versa.
void mlpeserialize(mlpensemble &obj, std::string &s_out);

// This function unserializes data structure from string.
void mlpeunserialize(const std::string &s_in, mlpensemble &obj);

// This function serializes data structure to C++ stream.
//
// Data stream generated by this function is same as  string  representation
// generated  by  string  version  of  serializer - alphanumeric characters,
// dots, underscores, minus signs, which are grouped into words separated by
// spaces and CR+LF.
//
// We recommend you to read comments on string version of serializer to find
// out more about serialization of AlGLIB objects.
void mlpeserialize(mlpensemble &obj, std::ostream &s_out);

// This function unserializes data structure from stream.
void mlpeunserialize(const std::istream &s_in, mlpensemble &obj);

// Like MLPCreate0, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreate1, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreate2, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateB0, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateB1, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateB2, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateR0, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateR1, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateR2, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateC0, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateC1, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Like MLPCreateC2, but for ensembles.
//
// ALGLIB: Copyright 18.02.2009 by Sergey Bochkanov
void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Creates ensemble from network. Only network geometry is copied.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble);

// Randomization of MLP ensemble
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlperandomize(const mlpensemble &ensemble);

// Return ensemble properties (number of inputs and outputs).
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout);

// Return normalization type (whether ensemble is SOFTMAX-normalized or not).
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
bool mlpeissoftmax(const mlpensemble &ensemble);

// Procesing
//
// Inputs:
//     Ensemble-   neural networks ensemble
//     X       -   input vector,  array[0..NIn-1].
//     Y       -   (possibly) preallocated buffer; if size of Y is less than
//                 NOut, it will be reallocated. If it is large enough, it
//                 is NOT reallocated, so we can save some time on reallocation.
//
//
// Outputs:
//     Y       -   result. Regression estimate when solving regression  task,
//                 vector of posterior probabilities for classification task.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);

// 'interactive'  variant  of  MLPEProcess  for  languages  like Python which
// support constructs like "Y = MLPEProcess(LM,X)" and interactive mode of the
// interpreter
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);

// Relative classification error on the test set
//
// Inputs:
//     Ensemble-   ensemble
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     percent of incorrectly classified cases.
//     Works both for classifier betwork and for regression networks which
// are used as classifiers.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);

// Average cross-entropy (in bits per element) on the test set
//
// Inputs:
//     Ensemble-   ensemble
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     CrossEntropy/(NPoints*LN(2)).
//     Zero if ensemble solves regression task.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);

// RMS error on the test set
//
// Inputs:
//     Ensemble-   ensemble
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     root mean square error.
//     Its meaning for regression task is obvious. As for classification task
// RMS error means error when estimating posterior probabilities.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);

// Average error on the test set
//
// Inputs:
//     Ensemble-   ensemble
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     Its meaning for regression task is obvious. As for classification task
// it means average error when estimating posterior probabilities.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);

// Average relative error on the test set
//
// Inputs:
//     Ensemble-   ensemble
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     Its meaning for regression task is obvious. As for classification task
// it means average relative error when estimating posterior probabilities.
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === MLPTRAIN Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS
// Depends on: (Optimization) MINLBFGS
// Depends on: MLPE
namespace alglib_impl {
typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
   ae_int_t ngrad;
   ae_int_t nhess;
   ae_int_t ncholesky;
} mlpreport;
void mlpreport_init(void *_p, bool make_automatic);
void mlpreport_copy(void *_dst, void *_src, bool make_automatic);
void mlpreport_free(void *_p, bool make_automatic);

typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
} mlpcvreport;
void mlpcvreport_init(void *_p, bool make_automatic);
void mlpcvreport_copy(void *_dst, void *_src, bool make_automatic);
void mlpcvreport_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector bestparameters;
   double bestrmserror;
   bool randomizenetwork;
   multilayerperceptron network;
   minlbfgsstate optimizer;
   minlbfgsreport optimizerrep;
   ae_vector wbuf0;
   ae_vector wbuf1;
   ae_vector allminibatches;
   ae_vector currentminibatch;
   ae_int_t PQ;
   ae_int_t algoused;
   ae_int_t minibatchsize;
   hqrndstate generator;
} smlptrnsession;
void smlptrnsession_init(void *_p, bool make_automatic);
void smlptrnsession_copy(void *_dst, void *_src, bool make_automatic);
void smlptrnsession_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector trnsubset;
   ae_vector valsubset;
   ae_shared_pool mlpsessions;
   mlpreport mlprep;
   multilayerperceptron network;
} mlpetrnsession;
void mlpetrnsession_init(void *_p, bool make_automatic);
void mlpetrnsession_copy(void *_dst, void *_src, bool make_automatic);
void mlpetrnsession_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t nin;
   ae_int_t nout;
   bool rcpar;
   ae_int_t lbfgsfactor;
   double decay;
   double wstep;
   ae_int_t maxits;
   ae_int_t datatype;
   ae_int_t npoints;
   ae_matrix densexy;
   sparsematrix sparsexy;
   smlptrnsession session;
   ae_int_t ngradbatch;
   ae_vector subset;
   ae_int_t subsetsize;
   ae_vector valsubset;
   ae_int_t valsubsetsize;
   ae_int_t algokind;
   ae_int_t minibatchsize;
} mlptrainer;
void mlptrainer_init(void *_p, bool make_automatic);
void mlptrainer_copy(void *_dst, void *_src, bool make_automatic);
void mlptrainer_free(void *_p, bool make_automatic);

typedef struct {
   multilayerperceptron network;
   mlpreport rep;
   ae_vector subset;
   ae_int_t subsetsize;
   ae_vector xyrow;
   ae_vector y;
   ae_int_t ngrad;
   ae_shared_pool trnpool;
} mlpparallelizationcv;
void mlpparallelizationcv_init(void *_p, bool make_automatic);
void mlpparallelizationcv_copy(void *_dst, void *_src, bool make_automatic);
void mlpparallelizationcv_free(void *_p, bool make_automatic);

void mlptrainlm(multilayerperceptron *network, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep);
void mlptrainlbfgs(multilayerperceptron *network, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t *info, mlpreport *rep);
void mlptraines(multilayerperceptron *network, RMatrix trnxy, ae_int_t trnsize, RMatrix valxy, ae_int_t valsize, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep);
void mlpkfoldcvlbfgs(multilayerperceptron *network, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t foldscount, ae_int_t *info, mlpreport *rep, mlpcvreport *cvrep);
void mlpkfoldcvlm(multilayerperceptron *network, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t foldscount, ae_int_t *info, mlpreport *rep, mlpcvreport *cvrep);
void mlpkfoldcv(mlptrainer *s, multilayerperceptron *network, ae_int_t nrestarts, ae_int_t foldscount, mlpreport *rep);
void mlpcreatetrainer(ae_int_t nin, ae_int_t nout, mlptrainer *s);
void mlpcreatetrainercls(ae_int_t nin, ae_int_t nclasses, mlptrainer *s);
void mlpsetdataset(mlptrainer *s, RMatrix xy, ae_int_t npoints);
void mlpsetsparsedataset(mlptrainer *s, sparsematrix *xy, ae_int_t npoints);
void mlpsetdecay(mlptrainer *s, double decay);
void mlpsetcond(mlptrainer *s, double wstep, ae_int_t maxits);
void mlpsetalgobatch(mlptrainer *s);
void mlptrainnetwork(mlptrainer *s, multilayerperceptron *network, ae_int_t nrestarts, mlpreport *rep);
void mlpstarttraining(mlptrainer *s, multilayerperceptron *network, bool randomstart);
bool mlpcontinuetraining(mlptrainer *s, multilayerperceptron *network);
void mlpebagginglm(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep, mlpcvreport *ooberrors);
void mlpebagginglbfgs(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t *info, mlpreport *rep, mlpcvreport *ooberrors);
void mlpetraines(mlpensemble *ensemble, RMatrix xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep);
void mlptrainensemblees(mlptrainer *s, mlpensemble *ensemble, ae_int_t nrestarts, mlpreport *rep);
} // end of namespace alglib_impl

namespace alglib {
// Training report:
//     * RelCLSError   -   fraction of misclassified cases.
//     * AvgCE         -   acerage cross-entropy
//     * RMSError      -   root-mean-square error
//     * AvgError      -   average error
//     * AvgRelError   -   average relative error
//     * NGrad         -   number of gradient calculations
//     * NHess         -   number of Hessian calculations
//     * NCholesky     -   number of Cholesky decompositions
//
// NOTE 1: RelCLSError/AvgCE are zero on regression problems.
//
// NOTE 2: on classification problems  RMSError/AvgError/AvgRelError  contain
//         errors in prediction of posterior probabilities
DecClass(mlpreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror; ae_int_t &ngrad; ae_int_t &nhess; ae_int_t &ncholesky;);

// Cross-validation estimates of generalization error
DecClass(mlpcvreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);

// Trainer object for neural network.
//
// You should not try to access fields of this object directly -  use  ALGLIB
// functions to work with this object.
DecClass(mlptrainer, EndD);

// Neural network training  using  modified  Levenberg-Marquardt  with  exact
// Hessian calculation and regularization. Subroutine trains  neural  network
// with restarts from random positions. Algorithm is well  suited  for  small
// and medium scale problems (hundreds of weights).
//
// Inputs:
//     Network     -   neural network with initialized geometry
//     XY          -   training set
//     NPoints     -   training set size
//     Decay       -   weight decay constant, >= 0.001
//                     Decay term 'Decay*||Weights||^2' is added to error
//                     function.
//                     If you don't know what Decay to choose, use 0.001.
//     Restarts    -   number of restarts from random position, >0.
//                     If you don't know what Restarts to choose, use 2.
//
// Outputs:
//     Network     -   trained neural network.
//     Info        -   return code:
//                     * -9, if internal matrix inverse subroutine failed
//                     * -2, if there is a point with class number
//                           outside of [0..NOut-1].
//                     * -1, if wrong parameters specified
//                           (NPoints<0, Restarts<1).
//                     *  2, if task has been solved.
//     Rep         -   training report
//
// ALGLIB: Copyright 10.03.2009 by Sergey Bochkanov
void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);

// Neural  network  training  using  L-BFGS  algorithm  with  regularization.
// Subroutine  trains  neural  network  with  restarts from random positions.
// Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
// requirements and step complexity are linear by weights number).
//
// Inputs:
//     Network     -   neural network with initialized geometry
//     XY          -   training set
//     NPoints     -   training set size
//     Decay       -   weight decay constant, >= 0.001
//                     Decay term 'Decay*||Weights||^2' is added to error
//                     function.
//                     If you don't know what Decay to choose, use 0.001.
//     Restarts    -   number of restarts from random position, >0.
//                     If you don't know what Restarts to choose, use 2.
//     WStep       -   stopping criterion. Algorithm stops if  step  size  is
//                     less than WStep. Recommended value - 0.01.  Zero  step
//                     size means stopping after MaxIts iterations.
//     MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
//                     iterations (NOT gradient  calculations).  Zero  MaxIts
//                     means stopping when step is sufficiently small.
//
// Outputs:
//     Network     -   trained neural network.
//     Info        -   return code:
//                     * -8, if both WStep=0 and MaxIts=0
//                     * -2, if there is a point with class number
//                           outside of [0..NOut-1].
//                     * -1, if wrong parameters specified
//                           (NPoints<0, Restarts<1).
//                     *  2, if task has been solved.
//     Rep         -   training report
//
// ALGLIB: Copyright 09.12.2007 by Sergey Bochkanov
void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep);

// Neural network training using early stopping (base algorithm - L-BFGS with
// regularization).
//
// Inputs:
//     Network     -   neural network with initialized geometry
//     TrnXY       -   training set
//     TrnSize     -   training set size, TrnSize>0
//     ValXY       -   validation set
//     ValSize     -   validation set size, ValSize>0
//     Decay       -   weight decay constant, >= 0.001
//                     Decay term 'Decay*||Weights||^2' is added to error
//                     function.
//                     If you don't know what Decay to choose, use 0.001.
//     Restarts    -   number of restarts, either:
//                     * strictly positive number - algorithm make specified
//                       number of restarts from random position.
//                     * -1, in which case algorithm makes exactly one run
//                       from the initial state of the network (no randomization).
//                     If you don't know what Restarts to choose, choose one
//                     one the following:
//                     * -1 (deterministic start)
//                     * +1 (one random restart)
//                     * +5 (moderate amount of random restarts)
//
// Outputs:
//     Network     -   trained neural network.
//     Info        -   return code:
//                     * -2, if there is a point with class number
//                           outside of [0..NOut-1].
//                     * -1, if wrong parameters specified
//                           (NPoints<0, Restarts<1, ...).
//                     *  2, task has been solved, stopping  criterion  met -
//                           sufficiently small step size.  Not expected  (we
//                           use  EARLY  stopping)  but  possible  and not an
//                           error.
//                     *  6, task has been solved, stopping  criterion  met -
//                           increasing of validation set error.
//     Rep         -   training report
//
// NOTE:
//
// Algorithm stops if validation set error increases for  a  long  enough  or
// step size is small enought  (there  are  task  where  validation  set  may
// decrease for eternity). In any case solution returned corresponds  to  the
// minimum of validation set error.
//
// ALGLIB: Copyright 10.03.2009 by Sergey Bochkanov
void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);

// Cross-validation estimate of generalization error.
//
// Base algorithm - L-BFGS.
//
// Inputs:
//     Network     -   neural network with initialized geometry.   Network is
//                     not changed during cross-validation -  it is used only
//                     as a representative of its architecture.
//     XY          -   training set.
//     SSize       -   training set size
//     Decay       -   weight  decay, same as in MLPTrainLBFGS
//     Restarts    -   number of restarts, >0.
//                     restarts are counted for each partition separately, so
//                     total number of restarts will be Restarts*FoldsCount.
//     WStep       -   stopping criterion, same as in MLPTrainLBFGS
//     MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
//     FoldsCount  -   number of folds in k-fold cross-validation,
//                     2 <= FoldsCount <= SSize.
//                     recommended value: 10.
//
// Outputs:
//     Info        -   return code, same as in MLPTrainLBFGS
//     Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
//     CVRep       -   generalization error estimates
//
// ALGLIB: Copyright 09.12.2007 by Sergey Bochkanov
void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);

// Cross-validation estimate of generalization error.
//
// Base algorithm - Levenberg-Marquardt.
//
// Inputs:
//     Network     -   neural network with initialized geometry.   Network is
//                     not changed during cross-validation -  it is used only
//                     as a representative of its architecture.
//     XY          -   training set.
//     SSize       -   training set size
//     Decay       -   weight  decay, same as in MLPTrainLBFGS
//     Restarts    -   number of restarts, >0.
//                     restarts are counted for each partition separately, so
//                     total number of restarts will be Restarts*FoldsCount.
//     FoldsCount  -   number of folds in k-fold cross-validation,
//                     2 <= FoldsCount <= SSize.
//                     recommended value: 10.
//
// Outputs:
//     Info        -   return code, same as in MLPTrainLBFGS
//     Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
//     CVRep       -   generalization error estimates
//
// ALGLIB: Copyright 09.12.2007 by Sergey Bochkanov
void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);

// This function estimates generalization error using cross-validation on the
// current dataset with current training settings.
//
// Inputs:
//     S           -   trainer object
//     Network     -   neural network. It must have same number of inputs and
//                     output/classes as was specified during creation of the
//                     trainer object. Network is not changed  during  cross-
//                     validation and is not trained - it  is  used  only  as
//                     representative of its architecture. I.e., we  estimate
//                     generalization properties of  ARCHITECTURE,  not  some
//                     specific network.
//     NRestarts   -   number of restarts, >= 0:
//                     * NRestarts>0  means  that  for  each cross-validation
//                       round   specified  number   of  random  restarts  is
//                       performed,  with  best  network  being  chosen after
//                       training.
//                     * NRestarts=0 is same as NRestarts=1
//     FoldsCount  -   number of folds in k-fold cross-validation:
//                     * 2 <= FoldsCount <= size of dataset
//                     * recommended value: 10.
//                     * values larger than dataset size will be silently
//                       truncated down to dataset size
//
// Outputs:
//     Rep         -   structure which contains cross-validation estimates:
//                     * Rep.RelCLSError - fraction of misclassified cases.
//                     * Rep.AvgCE - acerage cross-entropy
//                     * Rep.RMSError - root-mean-square error
//                     * Rep.AvgError - average error
//                     * Rep.AvgRelError - average relative error
//
// NOTE: when no dataset was specified with MLPSetDataset/SetSparseDataset(),
//       or subset with only one point  was  given,  zeros  are  returned  as
//       estimates.
//
// NOTE: this method performs FoldsCount cross-validation  rounds,  each  one
//       with NRestarts random starts.  Thus,  FoldsCount*NRestarts  networks
//       are trained in total.
//
// NOTE: Rep.RelCLSError/Rep.AvgCE are zero on regression problems.
//
// NOTE: on classification problems Rep.RMSError/Rep.AvgError/Rep.AvgRelError
//       contain errors in prediction of posterior probabilities.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep);

// Creation of the network trainer object for regression networks
//
// Inputs:
//     NIn         -   number of inputs, NIn >= 1
//     NOut        -   number of outputs, NOut >= 1
//
// Outputs:
//     S           -   neural network trainer object.
//                     This structure can be used to train any regression
//                     network with NIn inputs and NOut outputs.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpcreatetrainer(const ae_int_t nin, const ae_int_t nout, mlptrainer &s);

// Creation of the network trainer object for classification networks
//
// Inputs:
//     NIn         -   number of inputs, NIn >= 1
//     NClasses    -   number of classes, NClasses >= 2
//
// Outputs:
//     S           -   neural network trainer object.
//                     This structure can be used to train any classification
//                     network with NIn inputs and NOut outputs.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpcreatetrainercls(const ae_int_t nin, const ae_int_t nclasses, mlptrainer &s);

// This function sets "current dataset" of the trainer object to  one  passed
// by user.
//
// Inputs:
//     S           -   trainer object
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed.
//     NPoints     -   points count, >= 0.
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// datasetformat is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpsetdataset(const mlptrainer &s, const real_2d_array &xy, const ae_int_t npoints);

// This function sets "current dataset" of the trainer object to  one  passed
// by user (sparse matrix is used to store dataset).
//
// Inputs:
//     S           -   trainer object
//     XY          -   training  set,  see  below  for  information  on   the
//                     training set format. This function checks  correctness
//                     of  the  dataset  (no  NANs/INFs,  class  numbers  are
//                     correct) and throws exception when  incorrect  dataset
//                     is passed. Any  sparse  storage  format  can be  used:
//                     Hash-table, CRS...
//     NPoints     -   points count, >= 0
//
// DATASET FORMAT:
//
// This  function  uses  two  different  dataset formats - one for regression
// networks, another one for classification networks.
//
// For regression networks with NIn inputs and NOut outputs following dataset
// format is used:
// * dataset is given by NPoints*(NIn+NOut) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, next NOut columns are outputs
//
// For classification networks with NIn inputs and NClasses clases  following
// datasetformat is used:
// * dataset is given by NPoints*(NIn+1) matrix
// * each row corresponds to one example
// * first NIn columns are inputs, last column stores class number (from 0 to
//   NClasses-1).
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpsetsparsedataset(const mlptrainer &s, const sparsematrix &xy, const ae_int_t npoints);

// This function sets weight decay coefficient which is used for training.
//
// Inputs:
//     S           -   trainer object
//     Decay       -   weight  decay  coefficient, >= 0.  Weight  decay  term
//                     'Decay*||Weights||^2' is added to error  function.  If
//                     you don't know what Decay to choose, use 1.0E-3.
//                     Weight decay can be set to zero,  in this case network
//                     is trained without weight decay.
//
// NOTE: by default network uses some small nonzero value for weight decay.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpsetdecay(const mlptrainer &s, const double decay);

// This function sets stopping criteria for the optimizer.
//
// Inputs:
//     S           -   trainer object
//     WStep       -   stopping criterion. Algorithm stops if  step  size  is
//                     less than WStep. Recommended value - 0.01.  Zero  step
//                     size means stopping after MaxIts iterations.
//                     WStep >= 0.
//     MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
//                     epochs (full passes over entire dataset).  Zero MaxIts
//                     means stopping when step is sufficiently small.
//                     MaxIts >= 0.
//
// NOTE: by default, WStep=0.005 and MaxIts=0 are used. These values are also
//       used when MLPSetCond() is called with WStep=0 and MaxIts=0.
//
// NOTE: these stopping criteria are used for all kinds of neural training  -
//       from "conventional" networks to early stopping ensembles. When  used
//       for "conventional" networks, they are  used  as  the  only  stopping
//       criteria. When combined with early stopping, they used as ADDITIONAL
//       stopping criteria which can terminate early stopping algorithm.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpsetcond(const mlptrainer &s, const double wstep, const ae_int_t maxits);

// This function sets training algorithm: batch training using L-BFGS will be
// used.
//
// This algorithm:
// * the most robust for small-scale problems, but may be too slow for  large
//   scale ones.
// * perfoms full pass through the dataset before performing step
// * uses conditions specified by MLPSetCond() for stopping
// * is default one used by trainer object
//
// Inputs:
//     S           -   trainer object
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpsetalgobatch(const mlptrainer &s);

// This function trains neural network passed to this function, using current
// dataset (one which was passed to MLPSetDataset() or MLPSetSparseDataset())
// and current training settings. Training  from  NRestarts  random  starting
// positions is performed, best network is chosen.
//
// Training is performed using current training algorithm.
//
// Inputs:
//     S           -   trainer object
//     Network     -   neural network. It must have same number of inputs and
//                     output/classes as was specified during creation of the
//                     trainer object.
//     NRestarts   -   number of restarts, >= 0:
//                     * NRestarts>0 means that specified  number  of  random
//                       restarts are performed, best network is chosen after
//                       training
//                     * NRestarts=0 means that current state of the  network
//                       is used for training.
//
// Outputs:
//     Network     -   trained network
//
// NOTE: when no dataset was specified with MLPSetDataset/SetSparseDataset(),
//       network  is  filled  by zero  values.  Same  behavior  for functions
//       MLPStartTraining and MLPContinueTraining.
//
// NOTE: this method uses sum-of-squares error function for training.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep);

// IMPORTANT: this is an "expert" version of the MLPTrain() function.  We  do
//            not recommend you to use it unless you are pretty sure that you
//            need ability to monitor training progress.
//
// This function performs step-by-step training of the neural  network.  Here
// "step-by-step" means that training  starts  with  MLPStartTraining() call,
// and then user subsequently calls MLPContinueTraining() to perform one more
// iteration of the training.
//
// After call to this function trainer object remembers network and  is ready
// to  train  it.  However,  no  training  is  performed  until first call to
// MLPContinueTraining() function. Subsequent calls  to MLPContinueTraining()
// will advance training progress one iteration further.
//
// EXAMPLE:
//     >
//     > ...initialize network and trainer object....
//     >
//     > MLPStartTraining(Trainer, Network, True)
//     > while MLPContinueTraining(Trainer, Network) do
//     >     ...visualize training progress...
//     >
//
// Inputs:
//     S           -   trainer object
//     Network     -   neural network. It must have same number of inputs and
//                     output/classes as was specified during creation of the
//                     trainer object.
//     RandomStart -   randomize network before training or not:
//                     * True  means  that  network  is  randomized  and  its
//                       initial state (one which was passed to  the  trainer
//                       object) is lost.
//                     * False  means  that  training  is  started  from  the
//                       current state of the network
//
// Outputs:
//     Network     -   neural network which is ready to training (weights are
//                     initialized, preprocessor is initialized using current
//                     training set)
//
// NOTE: this method uses sum-of-squares error function for training.
//
// NOTE: it is expected that trainer object settings are NOT  changed  during
//       step-by-step training, i.e. no  one  changes  stopping  criteria  or
//       training set during training. It is possible and there is no defense
//       against  such  actions,  but  algorithm  behavior  in  such cases is
//       undefined and can be unpredictable.
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
void mlpstarttraining(const mlptrainer &s, const multilayerperceptron &network, const bool randomstart);

// IMPORTANT: this is an "expert" version of the MLPTrain() function.  We  do
//            not recommend you to use it unless you are pretty sure that you
//            need ability to monitor training progress.
//
// This function performs step-by-step training of the neural  network.  Here
// "step-by-step" means that training starts  with  MLPStartTraining()  call,
// and then user subsequently calls MLPContinueTraining() to perform one more
// iteration of the training.
//
// This  function  performs  one  more  iteration of the training and returns
// either True (training continues) or False (training stopped). In case True
// was returned, Network weights are updated according to the  current  state
// of the optimization progress. In case False was  returned,  no  additional
// updates is performed (previous update of  the  network weights moved us to
// the final point, and no additional updates is needed).
//
// EXAMPLE:
//     >
//     > [initialize network and trainer object]
//     >
//     > MLPStartTraining(Trainer, Network, True)
//     > while MLPContinueTraining(Trainer, Network) do
//     >     [visualize training progress]
//     >
//
// Inputs:
//     S           -   trainer object
//     Network     -   neural  network  structure,  which  is  used to  store
//                     current state of the training process.
//
// Outputs:
//     Network     -   weights of the neural network  are  rewritten  by  the
//                     current approximation.
//
// NOTE: this method uses sum-of-squares error function for training.
//
// NOTE: it is expected that trainer object settings are NOT  changed  during
//       step-by-step training, i.e. no  one  changes  stopping  criteria  or
//       training set during training. It is possible and there is no defense
//       against  such  actions,  but  algorithm  behavior  in  such cases is
//       undefined and can be unpredictable.
//
// NOTE: It  is  expected that Network is the same one which  was  passed  to
//       MLPStartTraining() function.  However,  THIS  function  checks  only
//       following:
//       * that number of network inputs is consistent with trainer object
//         settings
//       * that number of network outputs/classes is consistent with  trainer
//         object settings
//       * that number of network weights is the same as number of weights in
//         the network passed to MLPStartTraining() function
//       Exception is thrown when these conditions are violated.
//
//       It is also expected that you do not change state of the  network  on
//       your own - the only party who has right to change network during its
//       training is a trainer object. Any attempt to interfere with  trainer
//       may lead to unpredictable results.
//
//
// ALGLIB: Copyright 23.07.2012 by Sergey Bochkanov
bool mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network);

// Training neural networks ensemble using  bootstrap  aggregating (bagging).
// Modified Levenberg-Marquardt algorithm is used as base training method.
//
// Inputs:
//     Ensemble    -   model with initialized geometry
//     XY          -   training set
//     NPoints     -   training set size
//     Decay       -   weight decay coefficient, >= 0.001
//     Restarts    -   restarts, >0.
//
// Outputs:
//     Ensemble    -   trained model
//     Info        -   return code:
//                     * -2, if there is a point with class number
//                           outside of [0..NClasses-1].
//                     * -1, if incorrect parameters was passed
//                           (NPoints<0, Restarts<1).
//                     *  2, if task has been solved.
//     Rep         -   training report.
//     OOBErrors   -   out-of-bag generalization error estimate
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);

// Training neural networks ensemble using  bootstrap  aggregating (bagging).
// L-BFGS algorithm is used as base training method.
//
// Inputs:
//     Ensemble    -   model with initialized geometry
//     XY          -   training set
//     NPoints     -   training set size
//     Decay       -   weight decay coefficient, >= 0.001
//     Restarts    -   restarts, >0.
//     WStep       -   stopping criterion, same as in MLPTrainLBFGS
//     MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
//
// Outputs:
//     Ensemble    -   trained model
//     Info        -   return code:
//                     * -8, if both WStep=0 and MaxIts=0
//                     * -2, if there is a point with class number
//                           outside of [0..NClasses-1].
//                     * -1, if incorrect parameters was passed
//                           (NPoints<0, Restarts<1).
//                     *  2, if task has been solved.
//     Rep         -   training report.
//     OOBErrors   -   out-of-bag generalization error estimate
//
// ALGLIB: Copyright 17.02.2009 by Sergey Bochkanov
void mlpebagginglbfgs(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);

// Training neural networks ensemble using early stopping.
//
// Inputs:
//     Ensemble    -   model with initialized geometry
//     XY          -   training set
//     NPoints     -   training set size
//     Decay       -   weight decay coefficient, >= 0.001
//     Restarts    -   restarts, >0.
//
// Outputs:
//     Ensemble    -   trained model
//     Info        -   return code:
//                     * -2, if there is a point with class number
//                           outside of [0..NClasses-1].
//                     * -1, if incorrect parameters was passed
//                           (NPoints<0, Restarts<1).
//                     *  6, if task has been solved.
//     Rep         -   training report.
//     OOBErrors   -   out-of-bag generalization error estimate
//
// ALGLIB: Copyright 10.03.2009 by Sergey Bochkanov
void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);

// This function trains neural network ensemble passed to this function using
// current dataset and early stopping training algorithm. Each early stopping
// round performs NRestarts  random  restarts  (thus,  EnsembleSize*NRestarts
// training rounds is performed in total).
//
// Inputs:
//     S           -   trainer object;
//     Ensemble    -   neural network ensemble. It must have same  number  of
//                     inputs and outputs/classes  as  was  specified  during
//                     creation of the trainer object.
//     NRestarts   -   number of restarts, >= 0:
//                     * NRestarts>0 means that specified  number  of  random
//                       restarts are performed during each ES round;
//                     * NRestarts=0 is silently replaced by 1.
//
// Outputs:
//     Ensemble    -   trained ensemble;
//     Rep         -   it contains all type of errors.
//
// NOTE: this training method uses BOTH early stopping and weight decay!  So,
//       you should select weight decay before starting training just as  you
//       select it before training "conventional" networks.
//
// NOTE: when no dataset was specified with MLPSetDataset/SetSparseDataset(),
//       or  single-point  dataset  was  passed,  ensemble  is filled by zero
//       values.
//
// NOTE: this method uses sum-of-squares error function for training.
//
// ALGLIB: Copyright 22.08.2012 by Sergey Bochkanov
void mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep);
} // end of namespace alglib

// === CLUSTERING Package ===
// Depends on: (AlgLibInternal) BLAS
// Depends on: (AlgLibMisc) HQRND
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
typedef struct {
   ae_matrix ct;
   ae_matrix ctbest;
   ae_vector xycbest;
   ae_vector xycprev;
   ae_vector d2;
   ae_vector csizes;
   apbuffers initbuf;
   ae_shared_pool updatepool;
} kmeansbuffers;
void kmeansbuffers_init(void *_p, bool make_automatic);
void kmeansbuffers_copy(void *_dst, void *_src, bool make_automatic);
void kmeansbuffers_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t npoints;
   ae_int_t nfeatures;
   ae_int_t disttype;
   ae_matrix xy;
   ae_matrix d;
   ae_int_t ahcalgo;
   ae_int_t kmeansrestarts;
   ae_int_t kmeansmaxits;
   ae_int_t kmeansinitalgo;
   bool kmeansdbgnoits;
   ae_int_t seed;
   ae_matrix tmpd;
   apbuffers distbuf;
   kmeansbuffers kmeanstmp;
} clusterizerstate;
void clusterizerstate_init(void *_p, bool make_automatic);
void clusterizerstate_copy(void *_dst, void *_src, bool make_automatic);
void clusterizerstate_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t terminationtype;
   ae_int_t npoints;
   ae_vector p;
   ae_matrix z;
   ae_matrix pz;
   ae_matrix pm;
   ae_vector mergedist;
} ahcreport;
void ahcreport_init(void *_p, bool make_automatic);
void ahcreport_copy(void *_dst, void *_src, bool make_automatic);
void ahcreport_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t npoints;
   ae_int_t nfeatures;
   ae_int_t terminationtype;
   ae_int_t iterationscount;
   double energy;
   ae_int_t k;
   ae_matrix c;
   ae_vector cidx;
} kmeansreport;
void kmeansreport_init(void *_p, bool make_automatic);
void kmeansreport_copy(void *_dst, void *_src, bool make_automatic);
void kmeansreport_free(void *_p, bool make_automatic);

void clusterizercreate(clusterizerstate *s);
void clusterizersetpoints(clusterizerstate *s, RMatrix xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype);
void clusterizersetdistances(clusterizerstate *s, RMatrix d, ae_int_t npoints, bool isupper);
void clusterizersetahcalgo(clusterizerstate *s, ae_int_t algo);
void clusterizersetkmeanslimits(clusterizerstate *s, ae_int_t restarts, ae_int_t maxits);
void clusterizersetkmeansinit(clusterizerstate *s, ae_int_t initalgo);
void clusterizersetseed(clusterizerstate *s, ae_int_t seed);
void clusterizerrunahc(clusterizerstate *s, ahcreport *rep);
void clusterizerrunkmeans(clusterizerstate *s, ae_int_t k, kmeansreport *rep);
void clusterizergetdistances(RMatrix xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, RMatrix d);
void clusterizergetdistancesbuf(apbuffers *buf, RMatrix xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, RMatrix d);
void clusterizergetkclusters(ahcreport *rep, ae_int_t k, ZVector cidx, ZVector cz);
void clusterizerseparatedbydist(ahcreport *rep, double r, ae_int_t *k, ZVector cidx, ZVector cz);
void clusterizerseparatedbycorr(ahcreport *rep, double r, ae_int_t *k, ZVector cidx, ZVector cz);
void kmeansinitbuf(kmeansbuffers *buf);
void kmeansgenerateinternal(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t initalgo, ae_int_t seed, ae_int_t maxits, ae_int_t restarts, bool kmeansdbgnoits, ae_int_t *info, ae_int_t *iterationscount, RMatrix ccol, bool needccol, RMatrix crow, bool needcrow, ZVector xyc, double *energy, kmeansbuffers *buf);
void kmeansupdatedistances(RMatrix xy, ae_int_t idx0, ae_int_t idx1, ae_int_t nvars, RMatrix ct, ae_int_t cidx0, ae_int_t cidx1, ZVector xyc, RVector xydist2, ae_shared_pool *bufferpool);
} // end of namespace alglib_impl

namespace alglib {
// This structure is a clusterization engine.
//
// You should not try to access its fields directly.
// Use ALGLIB functions in order to work with this object.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
DecClass(clusterizerstate, EndD);

// This structure  is used to store results of the agglomerative hierarchical
// clustering (AHC).
//
// Following information is returned:
//
// * TerminationType - completion code:
//   * 1   for successful completion of algorithm
//   * -5  inappropriate combination of  clustering  algorithm  and  distance
//         function was used. As for now, it  is  possible  only when  Ward's
//         method is called for dataset with non-Euclidean distance function.
//   In case negative completion code is returned,  other  fields  of  report
//   structure are invalid and should not be used.
//
// * NPoints contains number of points in the original dataset
//
// * Z contains information about merges performed  (see below).  Z  contains
//   indexes from the original (unsorted) dataset and it can be used when you
//   need to know what points were merged. However, it is not convenient when
//   you want to build a dendrograd (see below).
//
// * if  you  want  to  build  dendrogram, you  can use Z, but it is not good
//   option, because Z contains  indexes from  unsorted  dataset.  Dendrogram
//   built from such dataset is likely to have intersections. So, you have to
//   reorder you points before building dendrogram.
//   Permutation which reorders point is returned in P. Another representation
//   of  merges,  which  is  more  convenient for dendorgram construction, is
//   returned in PM.
//
// * more information on format of Z, P and PM can be found below and in the
//   examples from ALGLIB Reference Manual.
//
// FORMAL DESCRIPTION OF FIELDS:
//     NPoints         number of points
//     Z               array[NPoints-1,2],  contains   indexes   of  clusters
//                     linked in pairs to  form  clustering  tree.  I-th  row
//                     corresponds to I-th merge:
//                     * Z[I,0] - index of the first cluster to merge
//                     * Z[I,1] - index of the second cluster to merge
//                     * Z[I,0]<Z[I,1]
//                     * clusters are  numbered  from 0 to 2*NPoints-2,  with
//                       indexes from 0 to NPoints-1 corresponding to  points
//                       of the original dataset, and indexes from NPoints to
//                       2*NPoints-2  correspond  to  clusters  generated  by
//                       subsequent  merges  (I-th  row  of Z creates cluster
//                       with index NPoints+I).
//
//                     IMPORTANT: indexes in Z[] are indexes in the ORIGINAL,
//                     unsorted dataset. In addition to  Z algorithm  outputs
//                     permutation which rearranges points in such  way  that
//                     subsequent merges are  performed  on  adjacent  points
//                     (such order is needed if you want to build dendrogram).
//                     However,  indexes  in  Z  are  related  to   original,
//                     unrearranged sequence of points.
//
//     P               array[NPoints], permutation which reorders points  for
//                     dendrogram  construction.  P[i] contains  index of the
//                     position  where  we  should  move  I-th  point  of the
//                     original dataset in order to apply merges PZ/PM.
//
//     PZ              same as Z, but for permutation of points given  by  P.
//                     The  only  thing  which  changed  are  indexes  of the
//                     original points; indexes of clusters remained same.
//
//     MergeDist       array[NPoints-1], contains distances between  clusters
//                     being merged (MergeDist[i] correspond to merge  stored
//                     in Z[i,...]):
//                     * CLINK, SLINK and  average  linkage algorithms report
//                       "raw", unmodified distance metric.
//                     * Ward's   method   reports   weighted   intra-cluster
//                       variance, which is equal to ||Ca-Cb||^2 * Sa*Sb/(Sa+Sb).
//                       Here  A  and  B  are  clusters being merged, Ca is a
//                       center of A, Cb is a center of B, Sa is a size of A,
//                       Sb is a size of B.
//
//     PM              array[NPoints-1,6], another representation of  merges,
//                     which is suited for dendrogram construction. It  deals
//                     with rearranged points (permutation P is applied)  and
//                     represents merges in a form which different  from  one
//                     used by Z.
//                     For each I from 0 to NPoints-2, I-th row of PM represents
//                     merge performed on two clusters C0 and C1. Here:
//                     * C0 contains points with indexes PM[I,0]...PM[I,1]
//                     * C1 contains points with indexes PM[I,2]...PM[I,3]
//                     * indexes stored in PM are given for dataset sorted
//                       according to permutation P
//                     * PM[I,1]=PM[I,2]-1 (only adjacent clusters are merged)
//                     * PM[I,0] <= PM[I,1], PM[I,2] <= PM[I,3], i.e. both
//                       clusters contain at least one point
//                     * heights of "subdendrograms" corresponding  to  C0/C1
//                       are stored in PM[I,4]  and  PM[I,5].  Subdendrograms
//                       corresponding   to   single-point   clusters    have
//                       height=0. Dendrogram of the merge result has  height
//                       H=max(H0,H1)+1.
//
// NOTE: there is one-to-one correspondence between merges described by Z and
//       PM. I-th row of Z describes same merge of clusters as I-th row of PM,
//       with "left" cluster from Z corresponding to the "left" one from PM.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
DecClass(ahcreport, ae_int_t &terminationtype; ae_int_t &npoints; integer_1d_array p; integer_2d_array z; integer_2d_array pz; integer_2d_array pm; real_1d_array mergedist;);

// This  structure   is  used  to  store  results of the  k-means  clustering
// algorithm.
//
// Following information is always returned:
// * NPoints contains number of points in the original dataset
// * TerminationType contains completion code, negative on failure, positive
//   on success
// * K contains number of clusters
//
// For positive TerminationType we return:
// * NFeatures contains number of variables in the original dataset
// * C, which contains centers found by algorithm
// * CIdx, which maps points of the original dataset to clusters
//
// FORMAL DESCRIPTION OF FIELDS:
//     NPoints         number of points, >= 0
//     NFeatures       number of variables, >= 1
//     TerminationType completion code:
//                     * -5 if  distance  type  is  anything  different  from
//                          Euclidean metric
//                     * -3 for degenerate dataset: a) less  than  K  distinct
//                          points, b) K=0 for non-empty dataset.
//                     * +1 for successful completion
//     K               number of clusters
//     C               array[K,NFeatures], rows of the array store centers
//     CIdx            array[NPoints], which contains cluster indexes
//     IterationsCount actual number of iterations performed by clusterizer.
//                     If algorithm performed more than one random restart,
//                     total number of iterations is returned.
//     Energy          merit function, "energy", sum  of  squared  deviations
//                     from cluster centers
//
// ALGLIB: Copyright 27.11.2012 by Sergey Bochkanov
DecClass(kmeansreport, ae_int_t &npoints; ae_int_t &nfeatures; ae_int_t &terminationtype; ae_int_t &iterationscount; double &energy; ae_int_t &k; real_2d_array c; integer_1d_array cidx;);

// This function initializes clusterizer object. Newly initialized object  is
// empty, i.e. it does not contain dataset. You should use it as follows:
// 1. creation
// 2. dataset is added with ClusterizerSetPoints()
// 3. additional parameters are set
// 3. clusterization is performed with one of the clustering functions
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizercreate(clusterizerstate &s);

// This function adds dataset to the clusterizer structure.
//
// This function overrides all previous calls  of  ClusterizerSetPoints()  or
// ClusterizerSetDistances().
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     XY      -   array[NPoints,NFeatures], dataset
//     NPoints -   number of points, >= 0
//     NFeatures-  number of features, >= 1
//     DistType-   distance function:
//                 *  0    Chebyshev distance  (L-inf norm)
//                 *  1    city block distance (L1 norm)
//                 *  2    Euclidean distance  (L2 norm), non-squared
//                 * 10    Pearson correlation:
//                         dist(a,b) = 1-corr(a,b)
//                 * 11    Absolute Pearson correlation:
//                         dist(a,b) = 1-|corr(a,b)|
//                 * 12    Uncentered Pearson correlation (cosine of the angle):
//                         dist(a,b) = a'*b/(|a|*|b|)
//                 * 13    Absolute uncentered Pearson correlation
//                         dist(a,b) = |a'*b|/(|a|*|b|)
//                 * 20    Spearman rank correlation:
//                         dist(a,b) = 1-rankcorr(a,b)
//                 * 21    Absolute Spearman rank correlation
//                         dist(a,b) = 1-|rankcorr(a,b)|
//
// NOTE 1: different distance functions have different performance penalty:
//         * Euclidean or Pearson correlation distances are the fastest ones
//         * Spearman correlation distance function is a bit slower
//         * city block and Chebyshev distances are order of magnitude slower
//
//         The reason behing difference in performance is that correlation-based
//         distance functions are computed using optimized linear algebra kernels,
//         while Chebyshev and city block distance functions are computed using
//         simple nested loops with two branches at each iteration.
//
// NOTE 2: different clustering algorithms have different limitations:
//         * agglomerative hierarchical clustering algorithms may be used with
//           any kind of distance metric
//         * k-means++ clustering algorithm may be used only  with  Euclidean
//           distance function
//         Thus, list of specific clustering algorithms you may  use  depends
//         on distance function you specify when you set your dataset.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype);
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t disttype);

// This function adds dataset given by distance  matrix  to  the  clusterizer
// structure. It is important that dataset is not  given  explicitly  -  only
// distance matrix is given.
//
// This function overrides all previous calls  of  ClusterizerSetPoints()  or
// ClusterizerSetDistances().
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     D       -   array[NPoints,NPoints], distance matrix given by its upper
//                 or lower triangle (main diagonal is  ignored  because  its
//                 entries are expected to be zero).
//     NPoints -   number of points
//     IsUpper -   whether upper or lower triangle of D is given.
//
// NOTE 1: different clustering algorithms have different limitations:
//         * agglomerative hierarchical clustering algorithms may be used with
//           any kind of distance metric, including one  which  is  given  by
//           distance matrix
//         * k-means++ clustering algorithm may be used only  with  Euclidean
//           distance function and explicitly given points - it  can  not  be
//           used with dataset given by distance matrix
//         Thus, if you call this function, you will be unable to use k-means
//         clustering algorithm to process your problem.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const ae_int_t npoints, const bool isupper);
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const bool isupper);

// This function sets agglomerative hierarchical clustering algorithm
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     Algo    -   algorithm type:
//                 * 0     complete linkage (default algorithm)
//                 * 1     single linkage
//                 * 2     unweighted average linkage
//                 * 3     weighted average linkage
//                 * 4     Ward's method
//
// NOTE: Ward's method works correctly only with Euclidean  distance,  that's
//       why algorithm will return negative termination  code  (failure)  for
//       any other distance type.
//
//       It is possible, however,  to  use  this  method  with  user-supplied
//       distance matrix. It  is  your  responsibility  to pass one which was
//       calculated with Euclidean distance function.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizersetahcalgo(const clusterizerstate &s, const ae_int_t algo);

// This  function  sets k-means properties:  number  of  restarts and maximum
// number of iterations per one run.
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     Restarts-   restarts count, >= 1.
//                 k-means++ algorithm performs several restarts and  chooses
//                 best set of centers (one with minimum squared distance).
//     MaxIts  -   maximum number of k-means iterations performed during  one
//                 run. >= 0, zero value means that algorithm performs unlimited
//                 number of iterations.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizersetkmeanslimits(const clusterizerstate &s, const ae_int_t restarts, const ae_int_t maxits);

// This function sets k-means  initialization  algorithm.  Several  different
// algorithms can be chosen, including k-means++.
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     InitAlgo-   initialization algorithm:
//                 * 0  automatic selection ( different  versions  of  ALGLIB
//                      may select different algorithms)
//                 * 1  random initialization
//                 * 2  k-means++ initialization  (best  quality  of  initial
//                      centers, but long  non-parallelizable  initialization
//                      phase with bad cache locality)
//                 * 3  "fast-greedy"  algorithm  with  efficient,  easy   to
//                      parallelize initialization. Quality of initial centers
//                      is  somewhat  worse  than  that  of  k-means++.  This
//                      algorithm is a default one in the current version  of
//                      ALGLIB.
//                 *-1  "debug" algorithm which always selects first  K  rows
//                      of dataset; this algorithm is used for debug purposes
//                      only. Do not use it in the industrial code!
//
// ALGLIB: Copyright 21.01.2015 by Sergey Bochkanov
void clusterizersetkmeansinit(const clusterizerstate &s, const ae_int_t initalgo);

// This  function  sets  seed  which  is  used to initialize internal RNG. By
// default, deterministic seed is used - same for each run of clusterizer. If
// you specify non-deterministic  seed  value,  then  some  algorithms  which
// depend on random initialization (in current version: k-means)  may  return
// slightly different results after each run.
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     Seed    -   seed:
//                 * positive values = use deterministic seed for each run of
//                   algorithms which depend on random initialization
//                 * zero or negative values = use non-deterministic seed
//
// ALGLIB: Copyright 08.06.2017 by Sergey Bochkanov
void clusterizersetseed(const clusterizerstate &s, const ae_int_t seed);

// This function performs agglomerative hierarchical clustering
//
// NOTE: Agglomerative  hierarchical  clustering  algorithm  has two  phases:
//       distance matrix calculation and clustering  itself. Only first phase
//       (distance matrix  calculation)  is  accelerated  by  Intel  MKL  and
//       multithreading. Thus, acceleration is significant only for medium or
//       high-dimensional problems.
//
//       Although activating multithreading gives some speedup  over  single-
//       threaded execution, you  should  not  expect  nearly-linear  scaling
//       with respect to cores count.
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//
// Outputs:
//     Rep     -   clustering results; see description of AHCReport
//                 structure for more information.
//
// NOTE 1: hierarchical clustering algorithms require large amounts of memory.
//         In particular, this implementation needs  sizeof(double)*NPoints^2
//         bytes, which are used to store distance matrix. In  case  we  work
//         with user-supplied matrix, this amount is multiplied by 2 (we have
//         to store original matrix and to work with its copy).
//
//         For example, problem with 10000 points  would require 800M of RAM,
//         even when working in a 1-dimensional space.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizerrunahc(const clusterizerstate &s, ahcreport &rep);

// This function performs clustering by k-means++ algorithm.
//
// You may change algorithm properties by calling:
// * ClusterizerSetKMeansLimits() to change number of restarts or iterations
// * ClusterizerSetKMeansInit() to change initialization algorithm
//
// By  default,  one  restart  and  unlimited number of iterations are  used.
// Initialization algorithm is chosen automatically.
//
// NOTE: k-means clustering  algorithm has two  phases:  selection of initial
//       centers and clustering  itself.  ALGLIB  parallelizes  both  phases.
//       Parallel version is optimized for the following  scenario: medium or
//       high-dimensional problem (8 or more dimensions) with large number of
//       points and clusters. However, some speed-up  can  be  obtained  even
//       when assumptions above are violated.
//
// Inputs:
//     S       -   clusterizer state, initialized by ClusterizerCreate()
//     K       -   number of clusters, K >= 0.
//                 K  can  be  zero only when algorithm is called  for  empty
//                 dataset,  in   this   case   completion  code  is  set  to
//                 success (+1).
//                 If  K=0  and  dataset  size  is  non-zero,  we   can   not
//                 meaningfully assign points to some center  (there  are  no
//                 centers because K=0) and  return  -3  as  completion  code
//                 (failure).
//
// Outputs:
//     Rep     -   clustering results; see description of KMeansReport
//                 structure for more information.
//
// NOTE 1: k-means  clustering  can  be  performed  only  for  datasets  with
//         Euclidean  distance  function.  Algorithm  will  return   negative
//         completion code in Rep.TerminationType in case dataset  was  added
//         to clusterizer with DistType other than Euclidean (or dataset  was
//         specified by distance matrix instead of explicitly given points).
//
// NOTE 2: by default, k-means uses non-deterministic seed to initialize  RNG
//         which is used to select initial centers. As  result,  each  run of
//         algorithm may return different values. If you  need  deterministic
//         behavior, use ClusterizerSetSeed() function.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep);

// This function returns distance matrix for dataset
//
// Inputs:
//     XY      -   array[NPoints,NFeatures], dataset
//     NPoints -   number of points, >= 0
//     NFeatures-  number of features, >= 1
//     DistType-   distance function:
//                 *  0    Chebyshev distance  (L-inf norm)
//                 *  1    city block distance (L1 norm)
//                 *  2    Euclidean distance  (L2 norm, non-squared)
//                 * 10    Pearson correlation:
//                         dist(a,b) = 1-corr(a,b)
//                 * 11    Absolute Pearson correlation:
//                         dist(a,b) = 1-|corr(a,b)|
//                 * 12    Uncentered Pearson correlation (cosine of the angle):
//                         dist(a,b) = a'*b/(|a|*|b|)
//                 * 13    Absolute uncentered Pearson correlation
//                         dist(a,b) = |a'*b|/(|a|*|b|)
//                 * 20    Spearman rank correlation:
//                         dist(a,b) = 1-rankcorr(a,b)
//                 * 21    Absolute Spearman rank correlation
//                         dist(a,b) = 1-|rankcorr(a,b)|
//
// Outputs:
//     D       -   array[NPoints,NPoints], distance matrix
//                 (full matrix is returned, with lower and upper triangles)
//
// NOTE:  different distance functions have different performance penalty:
//        * Euclidean or Pearson correlation distances are the fastest ones
//        * Spearman correlation distance function is a bit slower
//        * city block and Chebyshev distances are order of magnitude slower
//
//        The reason behing difference in performance is that correlation-based
//        distance functions are computed using optimized linear algebra kernels,
//        while Chebyshev and city block distance functions are computed using
//        simple nested loops with two branches at each iteration.
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d);

// This function takes as input clusterization report Rep,  desired  clusters
// count K, and builds top K clusters from hierarchical clusterization  tree.
// It returns assignment of points to clusters (array of cluster indexes).
//
// Inputs:
//     Rep     -   report from ClusterizerRunAHC() performed on XY
//     K       -   desired number of clusters, 1 <= K <= NPoints.
//                 K can be zero only when NPoints=0.
//
// Outputs:
//     CIdx    -   array[NPoints], I-th element contains cluster index  (from
//                 0 to K-1) for I-th point of the dataset.
//     CZ      -   array[K]. This array allows  to  convert  cluster  indexes
//                 returned by this function to indexes used by  Rep.Z.  J-th
//                 cluster returned by this function corresponds to  CZ[J]-th
//                 cluster stored in Rep.Z/PZ/PM.
//                 It is guaranteed that CZ[I]<CZ[I+1].
//
// NOTE: K clusters built by this subroutine are assumed to have no hierarchy.
//       Although  they  were  obtained  by  manipulation with top K nodes of
//       dendrogram  (i.e.  hierarchical  decomposition  of  dataset),   this
//       function does not return information about hierarchy.  Each  of  the
//       clusters stand on its own.
//
// NOTE: Cluster indexes returned by this function  does  not  correspond  to
//       indexes returned in Rep.Z/PZ/PM. Either you work  with  hierarchical
//       representation of the dataset (dendrogram), or you work with  "flat"
//       representation returned by this function.  Each  of  representations
//       has its own clusters indexing system (former uses [0, 2*NPoints-2]),
//       while latter uses [0..K-1]), although  it  is  possible  to  perform
//       conversion from one system to another by means of CZ array, returned
//       by this function, which allows you to convert indexes stored in CIdx
//       to the numeration system used by Rep.Z.
//
// NOTE: this subroutine is optimized for moderate values of K. Say, for  K=5
//       it will perform many times faster than  for  K=100.  Its  worst-case
//       performance is O(N*K), although in average case  it  perform  better
//       (up to O(N*log(K))).
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizergetkclusters(const ahcreport &rep, const ae_int_t k, integer_1d_array &cidx, integer_1d_array &cz);

// This  function  accepts  AHC  report  Rep,  desired  minimum  intercluster
// distance and returns top clusters from  hierarchical  clusterization  tree
// which are separated by distance R or HIGHER.
//
// It returns assignment of points to clusters (array of cluster indexes).
//
// There is one more function with similar name - ClusterizerSeparatedByCorr,
// which returns clusters with intercluster correlation equal to R  or  LOWER
// (note: higher for distance, lower for correlation).
//
// Inputs:
//     Rep     -   report from ClusterizerRunAHC() performed on XY
//     R       -   desired minimum intercluster distance, R >= 0
//
// Outputs:
//     K       -   number of clusters, 1 <= K <= NPoints
//     CIdx    -   array[NPoints], I-th element contains cluster index  (from
//                 0 to K-1) for I-th point of the dataset.
//     CZ      -   array[K]. This array allows  to  convert  cluster  indexes
//                 returned by this function to indexes used by  Rep.Z.  J-th
//                 cluster returned by this function corresponds to  CZ[J]-th
//                 cluster stored in Rep.Z/PZ/PM.
//                 It is guaranteed that CZ[I]<CZ[I+1].
//
// NOTE: K clusters built by this subroutine are assumed to have no hierarchy.
//       Although  they  were  obtained  by  manipulation with top K nodes of
//       dendrogram  (i.e.  hierarchical  decomposition  of  dataset),   this
//       function does not return information about hierarchy.  Each  of  the
//       clusters stand on its own.
//
// NOTE: Cluster indexes returned by this function  does  not  correspond  to
//       indexes returned in Rep.Z/PZ/PM. Either you work  with  hierarchical
//       representation of the dataset (dendrogram), or you work with  "flat"
//       representation returned by this function.  Each  of  representations
//       has its own clusters indexing system (former uses [0, 2*NPoints-2]),
//       while latter uses [0..K-1]), although  it  is  possible  to  perform
//       conversion from one system to another by means of CZ array, returned
//       by this function, which allows you to convert indexes stored in CIdx
//       to the numeration system used by Rep.Z.
//
// NOTE: this subroutine is optimized for moderate values of K. Say, for  K=5
//       it will perform many times faster than  for  K=100.  Its  worst-case
//       performance is O(N*K), although in average case  it  perform  better
//       (up to O(N*log(K))).
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizerseparatedbydist(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);

// This  function  accepts  AHC  report  Rep,  desired  maximum  intercluster
// correlation and returns top clusters from hierarchical clusterization tree
// which are separated by correlation R or LOWER.
//
// It returns assignment of points to clusters (array of cluster indexes).
//
// There is one more function with similar name - ClusterizerSeparatedByDist,
// which returns clusters with intercluster distance equal  to  R  or  HIGHER
// (note: higher for distance, lower for correlation).
//
// Inputs:
//     Rep     -   report from ClusterizerRunAHC() performed on XY
//     R       -   desired maximum intercluster correlation, -1 <= R <= +1
//
// Outputs:
//     K       -   number of clusters, 1 <= K <= NPoints
//     CIdx    -   array[NPoints], I-th element contains cluster index  (from
//                 0 to K-1) for I-th point of the dataset.
//     CZ      -   array[K]. This array allows  to  convert  cluster  indexes
//                 returned by this function to indexes used by  Rep.Z.  J-th
//                 cluster returned by this function corresponds to  CZ[J]-th
//                 cluster stored in Rep.Z/PZ/PM.
//                 It is guaranteed that CZ[I]<CZ[I+1].
//
// NOTE: K clusters built by this subroutine are assumed to have no hierarchy.
//       Although  they  were  obtained  by  manipulation with top K nodes of
//       dendrogram  (i.e.  hierarchical  decomposition  of  dataset),   this
//       function does not return information about hierarchy.  Each  of  the
//       clusters stand on its own.
//
// NOTE: Cluster indexes returned by this function  does  not  correspond  to
//       indexes returned in Rep.Z/PZ/PM. Either you work  with  hierarchical
//       representation of the dataset (dendrogram), or you work with  "flat"
//       representation returned by this function.  Each  of  representations
//       has its own clusters indexing system (former uses [0, 2*NPoints-2]),
//       while latter uses [0..K-1]), although  it  is  possible  to  perform
//       conversion from one system to another by means of CZ array, returned
//       by this function, which allows you to convert indexes stored in CIdx
//       to the numeration system used by Rep.Z.
//
// NOTE: this subroutine is optimized for moderate values of K. Say, for  K=5
//       it will perform many times faster than  for  K=100.  Its  worst-case
//       performance is O(N*K), although in average case  it  perform  better
//       (up to O(N*log(K))).
//
// ALGLIB: Copyright 10.07.2012 by Sergey Bochkanov
void clusterizerseparatedbycorr(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);
} // end of namespace alglib

// === DFOREST Package ===
// Depends on: (AlgLibInternal) SCODES
// Depends on: (AlgLibMisc) HQRND
// Depends on: BDSS
namespace alglib_impl {
typedef struct {
   ae_int_t dstype;
   ae_int_t npoints;
   ae_int_t nvars;
   ae_int_t nclasses;
   ae_vector dsdata;
   ae_vector dsrval;
   ae_vector dsival;
   ae_int_t rdfalgo;
   double rdfratio;
   double rdfvars;
   ae_int_t rdfglobalseed;
   ae_int_t rdfsplitstrength;
   ae_int_t rdfimportance;
   ae_vector dsmin;
   ae_vector dsmax;
   ae_vector dsbinary;
   double dsravg;
   ae_vector dsctotals;
   ae_int_t rdfprogress;
   ae_int_t rdftotal;
   ae_shared_pool workpool;
   ae_shared_pool votepool;
   ae_shared_pool treepool;
   ae_shared_pool treefactory;
   bool neediobmatrix;
   ae_matrix iobmatrix;
   ae_vector varimpshuffle2;
} decisionforestbuilder;
void decisionforestbuilder_init(void *_p, bool make_automatic);
void decisionforestbuilder_copy(void *_dst, void *_src, bool make_automatic);
void decisionforestbuilder_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector classpriors;
   ae_vector varpool;
   ae_int_t varpoolsize;
   ae_vector trnset;
   ae_int_t trnsize;
   ae_vector trnlabelsr;
   ae_vector trnlabelsi;
   ae_vector oobset;
   ae_int_t oobsize;
   ae_vector ooblabelsr;
   ae_vector ooblabelsi;
   ae_vector treebuf;
   ae_vector curvals;
   ae_vector bestvals;
   ae_vector tmp0i;
   ae_vector tmp1i;
   ae_vector tmp0r;
   ae_vector tmp1r;
   ae_vector tmp2r;
   ae_vector tmp3r;
   ae_vector tmpnrms2;
   ae_vector classtotals0;
   ae_vector classtotals1;
   ae_vector classtotals01;
} dfworkbuf;
void dfworkbuf_init(void *_p, bool make_automatic);
void dfworkbuf_copy(void *_dst, void *_src, bool make_automatic);
void dfworkbuf_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector trntotals;
   ae_vector oobtotals;
   ae_vector trncounts;
   ae_vector oobcounts;
   ae_vector giniimportances;
} dfvotebuf;
void dfvotebuf_init(void *_p, bool make_automatic);
void dfvotebuf_copy(void *_dst, void *_src, bool make_automatic);
void dfvotebuf_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector losses;
   ae_vector xraw;
   ae_vector xdist;
   ae_vector xcur;
   ae_vector y;
   ae_vector yv;
   ae_vector targety;
   ae_vector startnodes;
} dfpermimpbuf;
void dfpermimpbuf_init(void *_p, bool make_automatic);
void dfpermimpbuf_copy(void *_dst, void *_src, bool make_automatic);
void dfpermimpbuf_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector treebuf;
   ae_int_t treeidx;
} dftreebuf;
void dftreebuf_init(void *_p, bool make_automatic);
void dftreebuf_copy(void *_dst, void *_src, bool make_automatic);
void dftreebuf_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector x;
   ae_vector y;
} decisionforestbuffer;
void decisionforestbuffer_init(void *_p, bool make_automatic);
void decisionforestbuffer_copy(void *_dst, void *_src, bool make_automatic);
void decisionforestbuffer_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t forestformat;
   bool usemantissa8;
   ae_int_t nvars;
   ae_int_t nclasses;
   ae_int_t ntrees;
   ae_int_t bufsize;
   ae_vector trees;
   decisionforestbuffer buffer;
   ae_vector trees8;
} decisionforest;
void decisionforest_init(void *_p, bool make_automatic);
void decisionforest_copy(void *_dst, void *_src, bool make_automatic);
void decisionforest_free(void *_p, bool make_automatic);

typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double oobrelclserror;
   double oobavgce;
   double oobrmserror;
   double oobavgerror;
   double oobavgrelerror;
   ae_vector topvars;
   ae_vector varimportances;
} dfreport;
void dfreport_init(void *_p, bool make_automatic);
void dfreport_copy(void *_dst, void *_src, bool make_automatic);
void dfreport_free(void *_p, bool make_automatic);

typedef struct {
   ae_vector treebuf;
   ae_vector idxbuf;
   ae_vector tmpbufr;
   ae_vector tmpbufr2;
   ae_vector tmpbufi;
   ae_vector classibuf;
   ae_vector sortrbuf;
   ae_vector sortrbuf2;
   ae_vector sortibuf;
   ae_vector varpool;
   ae_vector evsbin;
   ae_vector evssplits;
} dfinternalbuffers;
void dfinternalbuffers_init(void *_p, bool make_automatic);
void dfinternalbuffers_copy(void *_dst, void *_src, bool make_automatic);
void dfinternalbuffers_free(void *_p, bool make_automatic);

void dfcreatebuffer(decisionforest *model, decisionforestbuffer *buf);
void dfbuildercreate(decisionforestbuilder *s);
void dfbuildersetdataset(decisionforestbuilder *s, RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses);
void dfbuildersetrndvars(decisionforestbuilder *s, ae_int_t rndvars);
void dfbuildersetrndvarsratio(decisionforestbuilder *s, double f);
void dfbuildersetrndvarsauto(decisionforestbuilder *s);
void dfbuildersetsubsampleratio(decisionforestbuilder *s, double f);
void dfbuildersetseed(decisionforestbuilder *s, ae_int_t seedval);
void dfbuildersetrdfalgo(decisionforestbuilder *s, ae_int_t algotype);
void dfbuildersetrdfsplitstrength(decisionforestbuilder *s, ae_int_t splitstrength);
void dfbuildersetimportancetrngini(decisionforestbuilder *s);
void dfbuildersetimportanceoobgini(decisionforestbuilder *s);
void dfbuildersetimportancepermutation(decisionforestbuilder *s);
void dfbuildersetimportancenone(decisionforestbuilder *s);
double dfbuildergetprogress(decisionforestbuilder *s);
double dfbuilderpeekprogress(decisionforestbuilder *s);
void dfbuilderbuildrandomforest(decisionforestbuilder *s, ae_int_t ntrees, decisionforest *df, dfreport *rep);
double dfbinarycompression(decisionforest *df);
double dfbinarycompression8(decisionforest *df);
void dfprocess(decisionforest *df, RVector x, RVector y);
void dfprocessi(decisionforest *df, RVector x, RVector y);
double dfprocess0(decisionforest *model, RVector x);
ae_int_t dfclassify(decisionforest *model, RVector x);
void dftsprocess(decisionforest *df, decisionforestbuffer *buf, RVector x, RVector y);
double dfrelclserror(decisionforest *df, RMatrix xy, ae_int_t npoints);
double dfavgce(decisionforest *df, RMatrix xy, ae_int_t npoints);
double dfrmserror(decisionforest *df, RMatrix xy, ae_int_t npoints);
double dfavgerror(decisionforest *df, RMatrix xy, ae_int_t npoints);
double dfavgrelerror(decisionforest *df, RMatrix xy, ae_int_t npoints);
void dfcopy(decisionforest *df1, decisionforest *df2);
void dfalloc(ae_serializer *s, decisionforest *forest);
void dfserialize(ae_serializer *s, decisionforest *forest);
void dfunserialize(ae_serializer *s, decisionforest *forest);
void dfbuildrandomdecisionforest(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, double r, ae_int_t *info, decisionforest *df, dfreport *rep);
void dfbuildrandomdecisionforestx1(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t nrndvars, double r, ae_int_t *info, decisionforest *df, dfreport *rep);
void dfbuildinternal(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t samplesize, ae_int_t nfeatures, ae_int_t flags, ae_int_t *info, decisionforest *df, dfreport *rep);
} // end of namespace alglib_impl

namespace alglib {
// A random forest (decision forest) builder object.
//
// Used to store dataset and specify decision forest training algorithm settings.
DecClass(decisionforestbuilder, EndD);

// Buffer object which is used to perform  various  requests  (usually  model
// inference) in the multithreaded mode (multiple threads working  with  same
// DF object).
//
// This object should be created with DFCreateBuffer().
DecClass(decisionforestbuffer, EndD);

// Decision forest (random forest) model.
DecClass(decisionforest, EndD);

// Decision forest training report.
//
// ==== training/oob errors ====
//
// Following fields store training set errors:
// * relclserror           -   fraction of misclassified cases, [0,1]
// * avgce                 -   average cross-entropy in bits per symbol
// * rmserror              -   root-mean-square error
// * avgerror              -   average error
// * avgrelerror           -   average relative error
//
// Out-of-bag estimates are stored in fields with same names, but "oob" prefix.
//
// For classification problems:
// * RMS, AVG and AVGREL errors are calculated for posterior probabilities
//
// For regression problems:
// * RELCLS and AVGCE errors are zero
//
// ==== variable importance ====
//
// Following fields are used to store variable importance information:
//
// * topvars               -   variables ordered from the most  important  to
//                             less  important  ones  (according  to  current
//                             choice of importance raiting).
//                             For example, topvars[0] contains index of  the
//                             most important variable, and topvars[0:2]  are
//                             indexes of 3 most important ones and so on.
//
// * varimportances        -   array[nvars], ratings (the  larger,  the  more
//                             important the variable  is,  always  in  [0,1]
//                             range).
//                             By default, filled  by  zeros  (no  importance
//                             ratings are  provided  unless  you  explicitly
//                             request them).
//                             Zero rating means that variable is not important,
//                             however you will rarely encounter such a thing,
//                             in many cases  unimportant  variables  produce
//                             nearly-zero (but nonzero) ratings.
//
// Variable importance report must be EXPLICITLY requested by calling:
// * dfbuildersetimportancegini() function, if you need out-of-bag Gini-based
//   importance rating also known as MDI  (fast to  calculate,  resistant  to
//   overfitting  issues,   but   has   some   bias  towards  continuous  and
//   high-cardinality categorical variables)
// * dfbuildersetimportancetrngini() function, if you need training set Gini-
//   -based importance rating (what other packages typically report).
// * dfbuildersetimportancepermutation() function, if you  need  permutation-
//   based importance rating also known as MDA (slower to calculate, but less
//   biased)
// * dfbuildersetimportancenone() function,  if  you  do  not  need  importance
//   ratings - ratings will be zero, topvars[] will be [0,1,2,...]
//
// Different importance ratings (Gini or permutation) produce  non-comparable
// values. Although in all cases rating values lie in [0,1] range, there  are
// exist differences:
// * informally speaking, Gini importance rating tends to divide "unit amount
//   of importance"  between  several  important  variables, i.e. it produces
//   estimates which roughly sum to 1.0 (or less than 1.0, if your  task  can
//   not be solved exactly). If all variables  are  equally  important,  they
//   will have same rating,  roughly  1/NVars,  even  if  every  variable  is
//   critically important.
// * from the other side, permutation importance tells us what percentage  of
//   the model predictive power will be ruined  by  permuting  this  specific
//   variable. It does not produce estimates which  sum  to  one.  Critically
//   important variable will have rating close  to  1.0,  and  you  may  have
//   multiple variables with such a rating.
//
// More information on variable importance ratings can be found  in  comments
// on the dfbuildersetimportancegini() and dfbuildersetimportancepermutation()
// functions.
DecClass(dfreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror; double &oobrelclserror; double &oobavgce; double &oobrmserror; double &oobavgerror; double &oobavgrelerror; integer_1d_array topvars; real_1d_array varimportances;);

// This function serializes data structure to string.
//
// Important properties of s_out:
// * it contains alphanumeric characters, dots, underscores, minus signs
// * these symbols are grouped into words, which are separated by spaces
//   and Windows-style (CR+LF) newlines
// * although  serializer  uses  spaces and CR+LF as separators, you can
//   replace any separator character by arbitrary combination of spaces,
//   tabs, Windows or Unix newlines. It allows flexible reformatting  of
//   the  string  in  case you want to include it into text or XML file.
//   But you should not insert separators into the middle of the "words"
//   nor you should change case of letters.
// * s_out can be freely moved between 32-bit and 64-bit systems, little
//   and big endian machines, and so on. You can serialize structure  on
//   32-bit machine and unserialize it on 64-bit one (or vice versa), or
//   serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also
//   serialize  it  in  C++ version of ALGLIB and unserialize in C# one,
//   and vice versa.
void dfserialize(decisionforest &obj, std::string &s_out);

// This function unserializes data structure from string.
void dfunserialize(const std::string &s_in, decisionforest &obj);

// This function serializes data structure to C++ stream.
//
// Data stream generated by this function is same as  string  representation
// generated  by  string  version  of  serializer - alphanumeric characters,
// dots, underscores, minus signs, which are grouped into words separated by
// spaces and CR+LF.
//
// We recommend you to read comments on string version of serializer to find
// out more about serialization of AlGLIB objects.
void dfserialize(decisionforest &obj, std::ostream &s_out);

// This function unserializes data structure from stream.
void dfunserialize(const std::istream &s_in, decisionforest &obj);

// This function creates buffer  structure  which  can  be  used  to  perform
// parallel inference requests.
//
// DF subpackage  provides two sets of computing functions - ones  which  use
// internal buffer of DF model  (these  functions are single-threaded because
// they use same buffer, which can not  shared  between  threads),  and  ones
// which use external buffer.
//
// This function is used to initialize external buffer.
//
// Inputs:
//     Model       -   DF model which is associated with newly created buffer
//
// Outputs:
//     Buf         -   external buffer.
//
//
// IMPORTANT: buffer object should be used only with model which was used  to
//            initialize buffer. Any attempt to  use  buffer  with  different
//            object is dangerous - you  may   get  integrity  check  failure
//            (exception) because sizes of internal  arrays  do  not  fit  to
//            dimensions of the model structure.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void dfcreatebuffer(const decisionforest &model, decisionforestbuffer &buf);

// This subroutine creates DecisionForestBuilder  object  which  is  used  to
// train decision forests.
//
// By default, new builder stores empty dataset and some  reasonable  default
// settings. At the very least, you should specify dataset prior to  building
// decision forest. You can also tweak settings of  the  forest  construction
// algorithm (recommended, although default setting should work well).
//
// Following actions are mandatory:
// * calling dfbuildersetdataset() to specify dataset
// * calling dfbuilderbuildrandomforest()  to  build  decision  forest  using
//   current dataset and default settings
//
// Additionally, you may call:
// * dfbuildersetrndvars() or dfbuildersetrndvarsratio() to specify number of
//   variables randomly chosen for each split
// * dfbuildersetsubsampleratio() to specify fraction of the dataset randomly
//   subsampled to build each tree
// * dfbuildersetseed() to control random seed chosen for tree construction
//
// Inputs:
//     none
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildercreate(decisionforestbuilder &s);

// This subroutine adds dense dataset to the internal storage of the  builder
// object. Specifying your dataset in the dense format means that  the  dense
// version of the forest construction algorithm will be invoked.
//
// Inputs:
//     S           -   decision forest builder object
//     XY          -   array[NPoints,NVars+1] (minimum size; actual size  can
//                     be larger, only leading part is used anyway), dataset:
//                     * first NVars elements of each row store values of the
//                       independent variables
//                     * last  column  store class number (in 0...NClasses-1)
//                       or real value of the dependent variable
//     NPoints     -   number of rows in the dataset, NPoints >= 1
//     NVars       -   number of independent variables, NVars >= 1
//     NClasses    -   indicates type of the problem being solved:
//                     * NClasses >= 2 means  that  classification  problem  is
//                       solved  (last  column  of  the  dataset stores class
//                       number)
//                     * NClasses=1 means that regression problem  is  solved
//                       (last column of the dataset stores variable value)
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetdataset(const decisionforestbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses);

// This function sets number  of  variables  (in  [1,NVars]  range)  used  by
// decision forest construction algorithm.
//
// The default option is to use roughly sqrt(NVars) variables.
//
// Inputs:
//     S           -   decision forest builder object
//     RndVars     -   number of randomly selected variables; values  outside
//                     of [1,NVars] range are silently clipped.
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetrndvars(const decisionforestbuilder &s, const ae_int_t rndvars);

// This function sets number of variables used by decision forest construction
// algorithm as a fraction of total variable count (0,1) range.
//
// The default option is to use roughly sqrt(NVars) variables.
//
// Inputs:
//     S           -   decision forest builder object
//     F           -   round(NVars*F) variables are selected
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetrndvarsratio(const decisionforestbuilder &s, const double f);

// This function tells decision forest builder to automatically choose number
// of  variables  used  by  decision forest construction  algorithm.  Roughly
// sqrt(NVars) variables will be used.
//
// Inputs:
//     S           -   decision forest builder object
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetrndvarsauto(const decisionforestbuilder &s);

// This function sets size of dataset subsample generated the decision forest
// construction algorithm. Size is specified as a fraction of  total  dataset
// size.
//
// The default option is to use 50% of the dataset for training, 50% for  the
// OOB estimates. You can decrease fraction F down to 10%, 1% or  even  below
// in order to reduce overfitting.
//
// Inputs:
//     S           -   decision forest builder object
//     F           -   fraction of the dataset to use, in (0,1] range. Values
//                     outside of this range will  be  silently  clipped.  At
//                     least one element is always selected for the  training
//                     set.
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetsubsampleratio(const decisionforestbuilder &s, const double f);

// This function sets seed used by internal RNG for  random  subsampling  and
// random selection of variable subsets.
//
// By default random seed is used, i.e. every time you build decision forest,
// we seed generator with new value  obtained  from  system-wide  RNG.  Thus,
// decision forest builder returns non-deterministic results. You can  change
// such behavior by specyfing fixed positive seed value.
//
// Inputs:
//     S           -   decision forest builder object
//     SeedVal     -   seed value:
//                     * positive values are used for seeding RNG with fixed
//                       seed, i.e. subsequent runs on same data will return
//                       same decision forests
//                     * non-positive seed means that random seed is used
//                       for every run of builder, i.e. subsequent  runs  on
//                       same  datasets  will  return   slightly   different
//                       decision forests
//
// Outputs:
//     S           -   decision forest builder, see
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetseed(const decisionforestbuilder &s, const ae_int_t seedval);

// This function sets random decision forest construction algorithm.
//
// As for now, only one decision forest construction algorithm is supported -
// a dense "baseline" RDF algorithm.
//
// Inputs:
//     S           -   decision forest builder object
//     AlgoType    -   algorithm type:
//                     * 0 = baseline dense RDF
//
// Outputs:
//     S           -   decision forest builder, see
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetrdfalgo(const decisionforestbuilder &s, const ae_int_t algotype);

// This  function  sets  split  selection  algorithm used by decision  forest
// classifier. You may choose several algorithms, with  different  speed  and
// quality of the results.
//
// Inputs:
//     S           -   decision forest builder object
//     SplitStrength-  split type:
//                     * 0 = split at the random position, fastest one
//                     * 1 = split at the middle of the range
//                     * 2 = strong split at the best point of the range (default)
//
// Outputs:
//     S           -   decision forest builder, see
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuildersetrdfsplitstrength(const decisionforestbuilder &s, const ae_int_t splitstrength);

// This  function  tells  decision  forest  construction  algorithm  to   use
// Gini impurity based variable importance estimation (also known as MDI).
//
// This version of importance estimation algorithm analyzes mean decrease  in
// impurity (MDI) on training sample during  splits.  The result  is  divided
// by impurity at the root node in order to produce estimate in [0,1] range.
//
// Such estimates are fast to calculate and beautifully  normalized  (sum  to
// one) but have following downsides:
// * They ALWAYS sum to 1.0, even if output is completely unpredictable. I.e.
//   MDI allows to order variables by importance, but does not  tell us about
//   "absolute" importances of variables
// * there exist some bias towards continuous and high-cardinality categorical
//   variables
//
// NOTE: informally speaking, MDA (permutation importance) rating answers the
//       question  "what  part  of  the  model  predictive power is ruined by
//       permuting k-th variable?" while MDI tells us "what part of the model
//       predictive power was achieved due to usage of k-th variable".
//
//       Thus, MDA rates each variable independently at "0 to 1"  scale while
//       MDI (and OOB-MDI too) tends to divide "unit  amount  of  importance"
//       between several important variables.
//
//       If  all  variables  are  equally  important,  they  will  have  same
//       MDI/OOB-MDI rating, equal (for OOB-MDI: roughly equal)  to  1/NVars.
//       However, roughly  same  picture  will  be  produced   for  the  "all
//       variables provide information no one is critical" situation  and for
//       the "all variables are critical, drop any one, everything is ruined"
//       situation.
//
//       Contrary to that, MDA will rate critical variable as ~1.0 important,
//       and important but non-critical variable will  have  less  than  unit
//       rating.
//
// NOTE: quite an often MDA and MDI return same results. It generally happens
//       on problems with low test set error (a few  percents  at  most)  and
//       large enough training set to avoid overfitting.
//
//       The difference between MDA, MDI and OOB-MDI becomes  important  only
//       on "hard" tasks with high test set error and/or small training set.
//
// Inputs:
//     S           -   decision forest builder object
//
// Outputs:
//     S           -   decision forest builder object. Next call to the forest
//                     construction function will produce:
//                     * importance estimates in rep.varimportances field
//                     * variable ranks in rep.topvars field
//
// ALGLIB: Copyright 29.07.2019 by Sergey Bochkanov
void dfbuildersetimportancetrngini(const decisionforestbuilder &s);

// This  function  tells  decision  forest  construction  algorithm  to   use
// out-of-bag version of Gini variable importance estimation (also  known  as
// OOB-MDI).
//
// This version of importance estimation algorithm analyzes mean decrease  in
// impurity (MDI) on out-of-bag sample during splits. The result  is  divided
// by impurity at the root node in order to produce estimate in [0,1] range.
//
// Such estimates are fast to calculate and resistant to  overfitting  issues
// (thanks to the  out-of-bag  estimates  used). However, OOB Gini rating has
// following downsides:
// * there exist some bias towards continuous and high-cardinality categorical
//   variables
// * Gini rating allows us to order variables by importance, but it  is  hard
//   to define importance of the variable by itself.
//
// NOTE: informally speaking, MDA (permutation importance) rating answers the
//       question  "what  part  of  the  model  predictive power is ruined by
//       permuting k-th variable?" while MDI tells us "what part of the model
//       predictive power was achieved due to usage of k-th variable".
//
//       Thus, MDA rates each variable independently at "0 to 1"  scale while
//       MDI (and OOB-MDI too) tends to divide "unit  amount  of  importance"
//       between several important variables.
//
//       If  all  variables  are  equally  important,  they  will  have  same
//       MDI/OOB-MDI rating, equal (for OOB-MDI: roughly equal)  to  1/NVars.
//       However, roughly  same  picture  will  be  produced   for  the  "all
//       variables provide information no one is critical" situation  and for
//       the "all variables are critical, drop any one, everything is ruined"
//       situation.
//
//       Contrary to that, MDA will rate critical variable as ~1.0 important,
//       and important but non-critical variable will  have  less  than  unit
//       rating.
//
// NOTE: quite an often MDA and MDI return same results. It generally happens
//       on problems with low test set error (a few  percents  at  most)  and
//       large enough training set to avoid overfitting.
//
//       The difference between MDA, MDI and OOB-MDI becomes  important  only
//       on "hard" tasks with high test set error and/or small training set.
//
// Inputs:
//     S           -   decision forest builder object
//
// Outputs:
//     S           -   decision forest builder object. Next call to the forest
//                     construction function will produce:
//                     * importance estimates in rep.varimportances field
//                     * variable ranks in rep.topvars field
//
// ALGLIB: Copyright 29.07.2019 by Sergey Bochkanov
void dfbuildersetimportanceoobgini(const decisionforestbuilder &s);

// This  function  tells  decision  forest  construction  algorithm  to   use
// permutation variable importance estimator (also known as MDA).
//
// This version of importance estimation algorithm analyzes mean increase  in
// out-of-bag sum of squared  residuals  after  random  permutation  of  J-th
// variable. The result is divided by error computed with all variables being
// perturbed in order to produce R-squared-like estimate in [0,1] range.
//
// Such estimate  is  slower to calculate than Gini-based rating  because  it
// needs multiple inference runs for each of variables being studied.
//
// ALGLIB uses parallelized and highly  optimized  algorithm  which  analyzes
// path through the decision tree and allows  to  handle  most  perturbations
// in O(1) time; nevertheless, requesting MDA importances may increase forest
// construction time from 10% to 200% (or more,  if  you  have  thousands  of
// variables).
//
// However, MDA rating has following benefits over Gini-based ones:
// * no bias towards specific variable types
// * ability to directly evaluate "absolute" importance of some  variable  at
//   "0 to 1" scale (contrary to Gini-based rating, which returns comparative
//   importances).
//
// NOTE: informally speaking, MDA (permutation importance) rating answers the
//       question  "what  part  of  the  model  predictive power is ruined by
//       permuting k-th variable?" while MDI tells us "what part of the model
//       predictive power was achieved due to usage of k-th variable".
//
//       Thus, MDA rates each variable independently at "0 to 1"  scale while
//       MDI (and OOB-MDI too) tends to divide "unit  amount  of  importance"
//       between several important variables.
//
//       If  all  variables  are  equally  important,  they  will  have  same
//       MDI/OOB-MDI rating, equal (for OOB-MDI: roughly equal)  to  1/NVars.
//       However, roughly  same  picture  will  be  produced   for  the  "all
//       variables provide information no one is critical" situation  and for
//       the "all variables are critical, drop any one, everything is ruined"
//       situation.
//
//       Contrary to that, MDA will rate critical variable as ~1.0 important,
//       and important but non-critical variable will  have  less  than  unit
//       rating.
//
// NOTE: quite an often MDA and MDI return same results. It generally happens
//       on problems with low test set error (a few  percents  at  most)  and
//       large enough training set to avoid overfitting.
//
//       The difference between MDA, MDI and OOB-MDI becomes  important  only
//       on "hard" tasks with high test set error and/or small training set.
//
// Inputs:
//     S           -   decision forest builder object
//
// Outputs:
//     S           -   decision forest builder object. Next call to the forest
//                     construction function will produce:
//                     * importance estimates in rep.varimportances field
//                     * variable ranks in rep.topvars field
//
// ALGLIB: Copyright 29.07.2019 by Sergey Bochkanov
void dfbuildersetimportancepermutation(const decisionforestbuilder &s);

// This  function  tells  decision  forest  construction  algorithm  to  skip
// variable importance estimation.
//
// Inputs:
//     S           -   decision forest builder object
//
// Outputs:
//     S           -   decision forest builder object. Next call to the forest
//                     construction function will result in forest being built
//                     without variable importance estimation.
//
// ALGLIB: Copyright 29.07.2019 by Sergey Bochkanov
void dfbuildersetimportancenone(const decisionforestbuilder &s);

// This function is an alias for dfbuilderpeekprogress(), left in ALGLIB  for
// backward compatibility reasons.
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
double dfbuildergetprogress(const decisionforestbuilder &s);

// This function is used to peek into  decision  forest  construction process
// from some other thread and get current progress indicator.
//
// It returns value in [0,1].
//
// Inputs:
//     S           -   decision forest builder object used  to  build  forest
//                     in some other thread
//
// Result:
//     progress value, in [0,1]
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
double dfbuilderpeekprogress(const decisionforestbuilder &s);

// This subroutine builds decision forest according to current settings using
// dataset internally stored in the builder object. Dense algorithm is used.
//
// NOTE: this   function   uses   dense  algorithm  for  forest  construction
//       independently from the dataset format (dense or sparse).
//
// NOTE: forest built with this function is  stored  in-memory  using  64-bit
//       data structures for offsets/indexes/split values. It is possible  to
//       convert  forest  into  more  memory-efficient   compressed    binary
//       representation.  Depending  on  the  problem  properties,  3.7x-5.7x
//       compression factors are possible.
//
//       The downsides of compression are (a) slight reduction in  the  model
//       accuracy and (b) ~1.5x reduction in  the  inference  speed  (due  to
//       increased complexity of the storage format).
//
//       See comments on dfbinarycompression() for more info.
//
// Default settings are used by the algorithm; you can tweak  them  with  the
// help of the following functions:
// * dfbuildersetrfactor() - to control a fraction of the  dataset  used  for
//   subsampling
// * dfbuildersetrandomvars() - to control number of variables randomly chosen
//   for decision rule creation
//
// Inputs:
//     S           -   decision forest builder object
//     NTrees      -   NTrees >= 1, number of trees to train
//
// Outputs:
//     DF          -   decision forest. You can compress this forest to  more
//                     compact 16-bit representation with dfbinarycompression()
//     Rep         -   report, see below for information on its fields.
//
// ==== report information produced by forest construction function ====
//
// Decision forest training report includes following information:
// * training set errors
// * out-of-bag estimates of errors
// * variable importance ratings
//
// Following fields are used to store information:
// * training set errors are stored in rep.relclserror, rep.avgce, rep.rmserror,
//   rep.avgerror and rep.avgrelerror
// * out-of-bag estimates of errors are stored in rep.oobrelclserror, rep.oobavgce,
//   rep.oobrmserror, rep.oobavgerror and rep.oobavgrelerror
//
// Variable importance reports, if requested by dfbuildersetimportancegini(),
// dfbuildersetimportancetrngini() or dfbuildersetimportancepermutation()
// call, are stored in:
// * rep.varimportances field stores importance ratings
// * rep.topvars stores variable indexes ordered from the most important to
//   less important ones
//
// You can find more information about report fields in:
// * comments on dfreport structure
// * comments on dfbuildersetimportancegini function
// * comments on dfbuildersetimportancetrngini function
// * comments on dfbuildersetimportancepermutation function
//
// ALGLIB: Copyright 21.05.2018 by Sergey Bochkanov
void dfbuilderbuildrandomforest(const decisionforestbuilder &s, const ae_int_t ntrees, decisionforest &df, dfreport &rep);

// This function performs binary compression of the decision forest.
//
// Original decision forest produced by the  forest  builder  is stored using
// 64-bit representation for all numbers - offsets, variable  indexes,  split
// points.
//
// It is possible to significantly reduce model size by means of:
// * using compressed  dynamic encoding for integers  (offsets  and  variable
//   indexes), which uses just 1 byte to store small ints  (less  than  128),
//   just 2 bytes for larger values (less than 128^2) and so on
// * storing floating point numbers using 8-bit exponent and 16-bit mantissa
//
// As  result,  model  needs  significantly  less  memory (compression factor
// depends on  variable and class counts). In particular:
// * NVars<128   and NClasses<128 result in 4.4x-5.7x model size reduction
// * NVars<16384 and NClasses<128 result in 3.7x-4.5x model size reduction
//
// Such storage format performs lossless compression  of  all  integers,  but
// compression of floating point values (split values) is lossy, with roughly
// 0.01% relative error introduced during rounding. Thus, we recommend you to
// re-evaluate model accuracy after compression.
//
// Another downside  of  compression  is  ~1.5x reduction  in  the  inference
// speed due to necessity of dynamic decompression of the compressed model.
//
// Inputs:
//     DF      -   decision forest built by forest builder
//
// Outputs:
//     DF      -   replaced by compressed forest
//
// Result:
//     compression factor (in-RAM size of the compressed model vs than of the
//     uncompressed one), positive number larger than 1.0
//
// ALGLIB: Copyright 22.07.2019 by Sergey Bochkanov
double dfbinarycompression(const decisionforest &df);

// Inference using decision forest
//
// IMPORTANT: this  function  is  thread-unsafe  and  may   modify   internal
//            structures of the model! You can not use same model  object for
//            parallel evaluation from several threads.
//
//            Use dftsprocess()  with  independent  thread-local  buffers  if
//            you need thread-safe evaluation.
//
// Inputs:
//     DF      -   decision forest model
//     X       -   input vector,  array[NVars]
//     Y       -   possibly preallocated buffer, reallocated if too small
//
// Outputs:
//     Y       -   result. Regression estimate when solving regression  task,
//                 vector of posterior probabilities for classification task.
//
// See also DFProcessI.
//
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y);

// 'interactive' variant of DFProcess for languages like Python which support
// constructs like "Y = DFProcessI(DF,X)" and interactive mode of interpreter
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
//
// IMPORTANT: this  function  is  thread-unsafe  and  may   modify   internal
//            structures of the model! You can not use same model  object for
//            parallel evaluation from several threads.
//
//            Use dftsprocess()  with  independent  thread-local  buffers  if
//            you need thread-safe evaluation.
//
// ALGLIB: Copyright 28.02.2010 by Sergey Bochkanov
void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y);

// This function returns first component of the  inferred  vector  (i.e.  one
// with index #0).
//
// It is a convenience wrapper for dfprocess() intended for either:
// * 1-dimensional regression problems
// * 2-class classification problems
//
// In the former case this function returns inference result as scalar, which
// is definitely more convenient that wrapping it as vector.  In  the  latter
// case it returns probability of object belonging to class #0.
//
// If you call it for anything different from two cases above, it  will  work
// as defined, i.e. return y[0], although it is of less use in such cases.
//
// IMPORTANT: this function is thread-unsafe and modifies internal structures
//            of the model! You can not use same model  object  for  parallel
//            evaluation from several threads.
//
//            Use dftsprocess() with  independent  thread-local  buffers,  if
//            you need thread-safe evaluation.
//
// Inputs:
//     Model   -   DF model
//     X       -   input vector,  array[0..NVars-1].
//
// Result:
//     Y[0]
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double dfprocess0(const decisionforest &model, const real_1d_array &x);

// This function returns most probable class number for an  input  X.  It  is
// same as calling  dfprocess(model,x,y), then determining i=argmax(y[i]) and
// returning i.
//
// A class number in [0,NOut) range in returned for classification  problems,
// -1 is returned when this function is called for regression problems.
//
// IMPORTANT: this function is thread-unsafe and modifies internal structures
//            of the model! You can not use same model  object  for  parallel
//            evaluation from several threads.
//
//            Use dftsprocess()  with independent  thread-local  buffers,  if
//            you need thread-safe evaluation.
//
// Inputs:
//     Model   -   decision forest model
//     X       -   input vector,  array[0..NVars-1].
//
// Result:
//     class number, -1 for regression tasks
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
ae_int_t dfclassify(const decisionforest &model, const real_1d_array &x);

// Inference using decision forest
//
// Thread-safe procesing using external buffer for temporaries.
//
// This function is thread-safe (i.e .  you  can  use  same  DF   model  from
// multiple threads) as long as you use different buffer objects for different
// threads.
//
// Inputs:
//     DF      -   decision forest model
//     Buf     -   buffer object, must be  allocated  specifically  for  this
//                 model with dfcreatebuffer().
//     X       -   input vector,  array[NVars]
//     Y       -   possibly preallocated buffer, reallocated if too small
//
// Outputs:
//     Y       -   result. Regression estimate when solving regression  task,
//                 vector of posterior probabilities for classification task.
//
// See also DFProcessI.
//
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
void dftsprocess(const decisionforest &df, const decisionforestbuffer &buf, const real_1d_array &x, real_1d_array &y);

// Relative classification error on the test set
//
// Inputs:
//     DF      -   decision forest model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     percent of incorrectly classified cases.
//     Zero if model solves regression task.
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

// Average cross-entropy (in bits per element) on the test set
//
// Inputs:
//     DF      -   decision forest model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     CrossEntropy/(NPoints*LN(2)).
//     Zero if model solves regression task.
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

// RMS error on the test set
//
// Inputs:
//     DF      -   decision forest model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     root mean square error.
//     Its meaning for regression task is obvious. As for
//     classification task, RMS error means error when estimating posterior
//     probabilities.
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

// Average error on the test set
//
// Inputs:
//     DF      -   decision forest model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     Its meaning for regression task is obvious. As for
//     classification task, it means average error when estimating posterior
//     probabilities.
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

// Average relative error on the test set
//
// Inputs:
//     DF      -   decision forest model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     Its meaning for regression task is obvious. As for
//     classification task, it means average relative error when estimating
//     posterior probability of belonging to the correct class.
//
// ALGLIB: Copyright 16.02.2009 by Sergey Bochkanov
double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

// This subroutine builds random decision forest.
//
// --------- DEPRECATED VERSION! USE DECISION FOREST BUILDER OBJECT ---------
//
// ALGLIB: Copyright 19.02.2009 by Sergey Bochkanov
void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);

// This subroutine builds random decision forest.
//
// --------- DEPRECATED VERSION! USE DECISION FOREST BUILDER OBJECT ---------
//
// ALGLIB: Copyright 19.02.2009 by Sergey Bochkanov
void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);
} // end of namespace alglib

// === KNN Package ===
// Depends on: (AlgLibMisc) HQRND, NEARESTNEIGHBOR
// Depends on: BDSS
namespace alglib_impl {
typedef struct {
   kdtreerequestbuffer treebuf;
   ae_vector x;
   ae_vector y;
   ae_vector tags;
   ae_matrix xy;
} knnbuffer;
void knnbuffer_init(void *_p, bool make_automatic);
void knnbuffer_copy(void *_dst, void *_src, bool make_automatic);
void knnbuffer_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t dstype;
   ae_int_t npoints;
   ae_int_t nvars;
   bool iscls;
   ae_int_t nout;
   ae_matrix dsdata;
   ae_vector dsrval;
   ae_vector dsival;
   ae_int_t knnnrm;
} knnbuilder;
void knnbuilder_init(void *_p, bool make_automatic);
void knnbuilder_copy(void *_dst, void *_src, bool make_automatic);
void knnbuilder_free(void *_p, bool make_automatic);

typedef struct {
   ae_int_t nvars;
   ae_int_t nout;
   ae_int_t k;
   double eps;
   bool iscls;
   bool isdummy;
   kdtree tree;
   knnbuffer buffer;
} knnmodel;
void knnmodel_init(void *_p, bool make_automatic);
void knnmodel_copy(void *_dst, void *_src, bool make_automatic);
void knnmodel_free(void *_p, bool make_automatic);

typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
} knnreport;
void knnreport_init(void *_p, bool make_automatic);
void knnreport_copy(void *_dst, void *_src, bool make_automatic);
void knnreport_free(void *_p, bool make_automatic);

void knncreatebuffer(knnmodel *model, knnbuffer *buf);
void knnbuildercreate(knnbuilder *s);
void knnbuildersetdatasetreg(knnbuilder *s, RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nout);
void knnbuildersetdatasetcls(knnbuilder *s, RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses);
void knnbuildersetnorm(knnbuilder *s, ae_int_t nrmtype);
void knnbuilderbuildknnmodel(knnbuilder *s, ae_int_t k, double eps, knnmodel *model, knnreport *rep);
void knnrewritekeps(knnmodel *model, ae_int_t k, double eps);
void knnprocess(knnmodel *model, RVector x, RVector y);
double knnprocess0(knnmodel *model, RVector x);
ae_int_t knnclassify(knnmodel *model, RVector x);
void knnprocessi(knnmodel *model, RVector x, RVector y);
void knntsprocess(knnmodel *model, knnbuffer *buf, RVector x, RVector y);
double knnrelclserror(knnmodel *model, RMatrix xy, ae_int_t npoints);
double knnavgce(knnmodel *model, RMatrix xy, ae_int_t npoints);
double knnrmserror(knnmodel *model, RMatrix xy, ae_int_t npoints);
double knnavgerror(knnmodel *model, RMatrix xy, ae_int_t npoints);
double knnavgrelerror(knnmodel *model, RMatrix xy, ae_int_t npoints);
void knnallerrors(knnmodel *model, RMatrix xy, ae_int_t npoints, knnreport *rep);
void knnalloc(ae_serializer *s, knnmodel *model);
void knnserialize(ae_serializer *s, knnmodel *model);
void knnunserialize(ae_serializer *s, knnmodel *model);
} // end of namespace alglib_impl

namespace alglib {
// Buffer object which is used to perform  various  requests  (usually  model
// inference) in the multithreaded mode (multiple threads working  with  same
// KNN object).
//
// This object should be created with KNNCreateBuffer().
DecClass(knnbuffer, EndD);

// A KNN builder object; this object encapsulates  dataset  and  all  related
// settings, it is used to create an actual instance of KNN model.
DecClass(knnbuilder, EndD);

// KNN model, can be used for classification or regression
DecClass(knnmodel, EndD);

// KNN training report.
//
// Following fields store training set errors:
// * relclserror       -   fraction of misclassified cases, [0,1]
// * avgce             -   average cross-entropy in bits per symbol
// * rmserror          -   root-mean-square error
// * avgerror          -   average error
// * avgrelerror       -   average relative error
//
// For classification problems:
// * RMS, AVG and AVGREL errors are calculated for posterior probabilities
//
// For regression problems:
// * RELCLS and AVGCE errors are zero
DecClass(knnreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);

// This function serializes data structure to string.
//
// Important properties of s_out:
// * it contains alphanumeric characters, dots, underscores, minus signs
// * these symbols are grouped into words, which are separated by spaces
//   and Windows-style (CR+LF) newlines
// * although  serializer  uses  spaces and CR+LF as separators, you can
//   replace any separator character by arbitrary combination of spaces,
//   tabs, Windows or Unix newlines. It allows flexible reformatting  of
//   the  string  in  case you want to include it into text or XML file.
//   But you should not insert separators into the middle of the "words"
//   nor you should change case of letters.
// * s_out can be freely moved between 32-bit and 64-bit systems, little
//   and big endian machines, and so on. You can serialize structure  on
//   32-bit machine and unserialize it on 64-bit one (or vice versa), or
//   serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also
//   serialize  it  in  C++ version of ALGLIB and unserialize in C# one,
//   and vice versa.
void knnserialize(knnmodel &obj, std::string &s_out);

// This function unserializes data structure from string.
void knnunserialize(const std::string &s_in, knnmodel &obj);

// This function serializes data structure to C++ stream.
//
// Data stream generated by this function is same as  string  representation
// generated  by  string  version  of  serializer - alphanumeric characters,
// dots, underscores, minus signs, which are grouped into words separated by
// spaces and CR+LF.
//
// We recommend you to read comments on string version of serializer to find
// out more about serialization of AlGLIB objects.
void knnserialize(knnmodel &obj, std::ostream &s_out);

// This function unserializes data structure from stream.
void knnunserialize(const std::istream &s_in, knnmodel &obj);

// This function creates buffer  structure  which  can  be  used  to  perform
// parallel KNN requests.
//
// KNN subpackage provides two sets of computing functions - ones  which  use
// internal buffer of KNN model (these  functions are single-threaded because
// they use same buffer, which can not  shared  between  threads),  and  ones
// which use external buffer.
//
// This function is used to initialize external buffer.
//
// Inputs:
//     Model       -   KNN model which is associated with newly created buffer
//
// Outputs:
//     Buf         -   external buffer.
//
//
// IMPORTANT: buffer object should be used only with model which was used  to
//            initialize buffer. Any attempt to  use  buffer  with  different
//            object is dangerous - you  may   get  integrity  check  failure
//            (exception) because sizes of internal  arrays  do  not  fit  to
//            dimensions of the model structure.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knncreatebuffer(const knnmodel &model, knnbuffer &buf);

// This subroutine creates KNNBuilder object which is used to train KNN models.
//
// By default, new builder stores empty dataset and some  reasonable  default
// settings. At the very least, you should specify dataset prior to  building
// KNN model. You can also tweak settings of the model construction algorithm
// (recommended, although default settings should work well).
//
// Following actions are mandatory:
// * calling knnbuildersetdataset() to specify dataset
// * calling knnbuilderbuildknnmodel() to build KNN model using current
//   dataset and default settings
//
// Additionally, you may call:
// * knnbuildersetnorm() to change norm being used
//
// Inputs:
//     none
//
// Outputs:
//     S           -   KNN builder
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnbuildercreate(knnbuilder &s);

// Specifies regression problem (one or more continuous  output variables are
// predicted). There also exists "classification" version of this function.
//
// This subroutine adds dense dataset to the internal storage of the  builder
// object. Specifying your dataset in the dense format means that  the  dense
// version of the KNN construction algorithm will be invoked.
//
// Inputs:
//     S           -   KNN builder object
//     XY          -   array[NPoints,NVars+NOut] (note: actual  size  can  be
//                     larger, only leading part is used anyway), dataset:
//                     * first NVars elements of each row store values of the
//                       independent variables
//                     * next NOut elements store  values  of  the  dependent
//                       variables
//     NPoints     -   number of rows in the dataset, NPoints >= 1
//     NVars       -   number of independent variables, NVars >= 1
//     NOut        -   number of dependent variables, NOut >= 1
//
// Outputs:
//     S           -   KNN builder
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnbuildersetdatasetreg(const knnbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nout);

// Specifies classification problem (two  or  more  classes  are  predicted).
// There also exists "regression" version of this function.
//
// This subroutine adds dense dataset to the internal storage of the  builder
// object. Specifying your dataset in the dense format means that  the  dense
// version of the KNN construction algorithm will be invoked.
//
// Inputs:
//     S           -   KNN builder object
//     XY          -   array[NPoints,NVars+1] (note:   actual   size  can  be
//                     larger, only leading part is used anyway), dataset:
//                     * first NVars elements of each row store values of the
//                       independent variables
//                     * next element stores class index, in [0,NClasses)
//     NPoints     -   number of rows in the dataset, NPoints >= 1
//     NVars       -   number of independent variables, NVars >= 1
//     NClasses    -   number of classes, NClasses >= 2
//
// Outputs:
//     S           -   KNN builder
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnbuildersetdatasetcls(const knnbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses);

// This function sets norm type used for neighbor search.
//
// Inputs:
//     S           -   decision forest builder object
//     NormType    -   norm type:
//                     * 0      inf-norm
//                     * 1      1-norm
//                     * 2      Euclidean norm (default)
//
// Outputs:
//     S           -   decision forest builder
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnbuildersetnorm(const knnbuilder &s, const ae_int_t nrmtype);

// This subroutine builds KNN model  according  to  current  settings,  using
// dataset internally stored in the builder object.
//
// The model being built performs inference using Eps-approximate  K  nearest
// neighbors search algorithm, with:
// * K=1,  Eps=0 corresponding to the "nearest neighbor algorithm"
// * K>1,  Eps=0 corresponding to the "K nearest neighbors algorithm"
// * K >= 1, Eps>0 corresponding to "approximate nearest neighbors algorithm"
//
// An approximate KNN is a good option for high-dimensional  datasets  (exact
// KNN works slowly when dimensions count grows).
//
// An ALGLIB implementation of kd-trees is used to perform k-nn searches.
//
// Inputs:
//     S       -   KNN builder object
//     K       -   number of neighbors to search for, K >= 1
//     Eps     -   approximation factor:
//                 * Eps=0 means that exact kNN search is performed
//                 * Eps>0 means that (1+Eps)-approximate search is performed
//
// Outputs:
//     Model       -   KNN model
//     Rep         -   report
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnbuilderbuildknnmodel(const knnbuilder &s, const ae_int_t k, const double eps, knnmodel &model, knnreport &rep);

// Changing search settings of KNN model.
//
// K and EPS parameters of KNN  (AKNN)  search  are  specified  during  model
// construction. However, plain KNN algorithm with Euclidean distance  allows
// you to change them at any moment.
//
// NOTE: future versions of KNN model may support advanced versions  of  KNN,
//       such as NCA or LMNN. It is possible that such algorithms won't allow
//       you to change search settings on the fly. If you call this  function
//       for an algorithm which does not support on-the-fly changes, it  will
//       throw an exception.
//
// Inputs:
//     Model   -   KNN model
//     K       -   K >= 1, neighbors count
//     EPS     -   accuracy of the EPS-approximate NN search. Set to 0.0,  if
//                 you want to perform "classic" KNN search.  Specify  larger
//                 values  if  you  need  to  speed-up  high-dimensional  KNN
//                 queries.
//
// Outputs:
//     nothing on success, exception on failure
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnrewritekeps(const knnmodel &model, const ae_int_t k, const double eps);

// Inference using KNN model.
//
// See also knnprocess0(), knnprocessi() and knnclassify() for options with a
// bit more convenient interface.
//
// IMPORTANT: this function is thread-unsafe and modifies internal structures
//            of the model! You can not use same model  object  for  parallel
//            evaluation from several threads.
//
//            Use knntsprocess() with independent  thread-local  buffers,  if
//            you need thread-safe evaluation.
//
// Inputs:
//     Model   -   KNN model
//     X       -   input vector,  array[0..NVars-1].
//     Y       -   possible preallocated buffer. Reused if long enough.
//
// Outputs:
//     Y       -   result. Regression estimate when solving regression  task,
//                 vector of posterior probabilities for classification task.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnprocess(const knnmodel &model, const real_1d_array &x, real_1d_array &y);

// This function returns first component of the  inferred  vector  (i.e.  one
// with index #0).
//
// It is a convenience wrapper for knnprocess() intended for either:
// * 1-dimensional regression problems
// * 2-class classification problems
//
// In the former case this function returns inference result as scalar, which
// is definitely more convenient that wrapping it as vector.  In  the  latter
// case it returns probability of object belonging to class #0.
//
// If you call it for anything different from two cases above, it  will  work
// as defined, i.e. return y[0], although it is of less use in such cases.
//
// IMPORTANT: this function is thread-unsafe and modifies internal structures
//            of the model! You can not use same model  object  for  parallel
//            evaluation from several threads.
//
//            Use knntsprocess() with independent  thread-local  buffers,  if
//            you need thread-safe evaluation.
//
// Inputs:
//     Model   -   KNN model
//     X       -   input vector,  array[0..NVars-1].
//
// Result:
//     Y[0]
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double knnprocess0(const knnmodel &model, const real_1d_array &x);

// This function returns most probable class number for an  input  X.  It  is
// same as calling knnprocess(model,x,y), then determining i=argmax(y[i]) and
// returning i.
//
// A class number in [0,NOut) range in returned for classification  problems,
// -1 is returned when this function is called for regression problems.
//
// IMPORTANT: this function is thread-unsafe and modifies internal structures
//            of the model! You can not use same model  object  for  parallel
//            evaluation from several threads.
//
//            Use knntsprocess() with independent  thread-local  buffers,  if
//            you need thread-safe evaluation.
//
// Inputs:
//     Model   -   KNN model
//     X       -   input vector,  array[0..NVars-1].
//
// Result:
//     class number, -1 for regression tasks
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
ae_int_t knnclassify(const knnmodel &model, const real_1d_array &x);

// 'interactive' variant of knnprocess()  for  languages  like  Python  which
// support constructs like "y = knnprocessi(model,x)" and interactive mode of
// the interpreter.
//
// This function allocates new array on each call,  so  it  is  significantly
// slower than its 'non-interactive' counterpart, but it is  more  convenient
// when you call it from command line.
//
// IMPORTANT: this  function  is  thread-unsafe  and  may   modify   internal
//            structures of the model! You can not use same model  object for
//            parallel evaluation from several threads.
//
//            Use knntsprocess()  with  independent  thread-local  buffers if
//            you need thread-safe evaluation.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnprocessi(const knnmodel &model, const real_1d_array &x, real_1d_array &y);

// Thread-safe procesing using external buffer for temporaries.
//
// This function is thread-safe (i.e .  you  can  use  same  KNN  model  from
// multiple threads) as long as you use different buffer objects for different
// threads.
//
// Inputs:
//     Model   -   KNN model
//     Buf     -   buffer object, must be  allocated  specifically  for  this
//                 model with knncreatebuffer().
//     X       -   input vector,  array[NVars]
//
// Outputs:
//     Y       -   result, array[NOut].   Regression  estimate  when  solving
//                 regression task,  vector  of  posterior  probabilities for
//                 a classification task.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knntsprocess(const knnmodel &model, const knnbuffer &buf, const real_1d_array &x, real_1d_array &y);

// Relative classification error on the test set
//
// Inputs:
//     Model   -   KNN model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     percent of incorrectly classified cases.
//     Zero if model solves regression task.
//
// NOTE: if  you  need several different kinds of error metrics, it is better
//       to use knnallerrors() which computes all error metric  with just one
//       pass over dataset.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double knnrelclserror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);

// Average cross-entropy (in bits per element) on the test set
//
// Inputs:
//     Model   -   KNN model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     CrossEntropy/NPoints.
//     Zero if model solves regression task.
//
// NOTE: the cross-entropy metric is too unstable when used to  evaluate  KNN
//       models (such models can report exactly  zero probabilities),  so  we
//       do not recommend using it.
//
// NOTE: if  you  need several different kinds of error metrics, it is better
//       to use knnallerrors() which computes all error metric  with just one
//       pass over dataset.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double knnavgce(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);

// RMS error on the test set.
//
// Its meaning for regression task is obvious. As for classification problems,
// RMS error means error when estimating posterior probabilities.
//
// Inputs:
//     Model   -   KNN model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     root mean square error.
//
// NOTE: if  you  need several different kinds of error metrics, it is better
//       to use knnallerrors() which computes all error metric  with just one
//       pass over dataset.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double knnrmserror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);

// Average error on the test set
//
// Its meaning for regression task is obvious. As for classification problems,
// average error means error when estimating posterior probabilities.
//
// Inputs:
//     Model   -   KNN model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     average error
//
// NOTE: if  you  need several different kinds of error metrics, it is better
//       to use knnallerrors() which computes all error metric  with just one
//       pass over dataset.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double knnavgerror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);

// Average relative error on the test set
//
// Its meaning for regression task is obvious. As for classification problems,
// average relative error means error when estimating posterior probabilities.
//
// Inputs:
//     Model   -   KNN model
//     XY      -   test set
//     NPoints -   test set size
//
// Result:
//     average relative error
//
// NOTE: if  you  need several different kinds of error metrics, it is better
//       to use knnallerrors() which computes all error metric  with just one
//       pass over dataset.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
double knnavgrelerror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);

// Calculates all kinds of errors for the model in one call.
//
// Inputs:
//     Model   -   KNN model
//     XY      -   test set:
//                 * one row per point
//                 * first NVars columns store independent variables
//                 * depending on problem type:
//                   * next column stores class number in [0,NClasses) -  for
//                     classification problems
//                   * next NOut columns  store  dependent  variables  -  for
//                     regression problems
//     NPoints -   test set size, NPoints >= 0
//
// Outputs:
//     Rep     -   following fields are loaded with errors for both regression
//                 and classification models:
//                 * rep.rmserror - RMS error for the output
//                 * rep.avgerror - average error
//                 * rep.avgrelerror - average relative error
//                 following fields are set only  for classification  models,
//                 zero for regression ones:
//                 * relclserror   - relative classification error, in [0,1]
//                 * avgce - average cross-entropy in bits per dataset entry
//
// NOTE: the cross-entropy metric is too unstable when used to  evaluate  KNN
//       models (such models can report exactly  zero probabilities),  so  we
//       do not recommend using it.
//
// ALGLIB: Copyright 15.02.2019 by Sergey Bochkanov
void knnallerrors(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, knnreport &rep);
} // end of namespace alglib

// === DATACOMP Package ===
// Depends on: CLUSTERING
namespace alglib_impl {
void kmeansgenerate(RMatrix xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t restarts, ae_int_t *info, RMatrix c, ZVector xyc);
} // end of namespace alglib_impl

namespace alglib {
// k-means++ clusterization.
// Backward compatibility function, we recommend to use CLUSTERING subpackage
// as better replacement.
//
// ALGLIB: Copyright 21.03.2009 by Sergey Bochkanov
void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc);
} // end of namespace alglib

#endif // OnceOnly