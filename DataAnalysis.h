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
// Depends on: (LinAlg) EVD, SVD
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
void pcabuildbasis(RMatrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector *s2, RMatrix *v, ae_state *_state);
void pcatruncatedsubspace(RMatrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t nneeded, double eps, ae_int_t maxits, RVector *s2, RMatrix *v, ae_state *_state);
void pcatruncatedsubspacesparse(sparsematrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t nneeded, double eps, ae_int_t maxits, RVector *s2, RMatrix *v, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v, const xparams _xparams = xdefault);
void pcatruncatedsubspace(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nneeded, const double eps, const ae_int_t maxits, real_1d_array &s2, real_2d_array &v, const xparams _xparams = xdefault);
void pcatruncatedsubspacesparse(const sparsematrix &x, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nneeded, const double eps, const ae_int_t maxits, real_1d_array &s2, real_2d_array &v, const xparams _xparams = xdefault);
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

void dserrallocate(ae_int_t nclasses, RVector *buf, ae_state *_state);
void dserraccumulate(RVector *buf, RVector *y, RVector *desiredy, ae_state *_state);
void dserrfinish(RVector *buf, ae_state *_state);
void dsnormalize(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector *means, RVector *sigmas, ae_state *_state);
void dsnormalizec(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector *means, RVector *sigmas, ae_state *_state);
double dsgetmeanmindistance(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_state *_state);
void dstie(RVector *a, ae_int_t n, ZVector *ties, ae_int_t *tiecount, ZVector *p1, ZVector *p2, ae_state *_state);
void dstiefasti(RVector *a, ZVector *b, ae_int_t n, ZVector *ties, ae_int_t *tiecount, RVector *bufr, ZVector *bufi, ae_state *_state);
void dsoptimalsplit2(RVector *a, ZVector *c, ae_int_t n, ae_int_t *info, double *threshold, double *pal, double *pbl, double *par, double *pbr, double *cve, ae_state *_state);
void dsoptimalsplit2fast(RVector *a, ZVector *c, ZVector *tiesbuf, ZVector *cntbuf, RVector *bufr, ZVector *bufi, ae_int_t n, ae_int_t nc, double alpha, ae_int_t *info, double *threshold, double *rms, double *cvrms, ae_state *_state);
void dssplitk(RVector *a, ZVector *c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t *info, RVector *thresholds, ae_int_t *ni, double *cve, ae_state *_state);
void dsoptimalsplitk(RVector *a, ZVector *c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t *info, RVector *thresholds, ae_int_t *ni, double *cve, ae_state *_state);
void cvreport_init(void *_p, ae_state *_state, bool make_automatic);
void cvreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void cvreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve, const xparams _xparams = xdefault);
void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MLPBASE Package ===
// Depends on: (AlgLibInternal) HPCCORES
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
typedef struct {
   double f;
   ae_vector g;
} smlpgrad;
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

ae_int_t mlpgradsplitcost(ae_state *_state);
ae_int_t mlpgradsplitsize(ae_state *_state);
void mlpcreate0(ae_int_t nin, ae_int_t nout, multilayerperceptron *network, ae_state *_state);
void mlpcreate1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, multilayerperceptron *network, ae_state *_state);
void mlpcreate2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, multilayerperceptron *network, ae_state *_state);
void mlpcreateb0(ae_int_t nin, ae_int_t nout, double b, double d, multilayerperceptron *network, ae_state *_state);
void mlpcreateb1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double b, double d, multilayerperceptron *network, ae_state *_state);
void mlpcreateb2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double b, double d, multilayerperceptron *network, ae_state *_state);
void mlpcreater0(ae_int_t nin, ae_int_t nout, double a, double b, multilayerperceptron *network, ae_state *_state);
void mlpcreater1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double a, double b, multilayerperceptron *network, ae_state *_state);
void mlpcreater2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double a, double b, multilayerperceptron *network, ae_state *_state);
void mlpcreatec0(ae_int_t nin, ae_int_t nout, multilayerperceptron *network, ae_state *_state);
void mlpcreatec1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, multilayerperceptron *network, ae_state *_state);
void mlpcreatec2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, multilayerperceptron *network, ae_state *_state);
void mlpcopy(multilayerperceptron *network1, multilayerperceptron *network2, ae_state *_state);
void mlpcopyshared(multilayerperceptron *network1, multilayerperceptron *network2, ae_state *_state);
bool mlpsamearchitecture(multilayerperceptron *network1, multilayerperceptron *network2, ae_state *_state);
void mlpcopytunableparameters(multilayerperceptron *network1, multilayerperceptron *network2, ae_state *_state);
void mlpexporttunableparameters(multilayerperceptron *network, RVector *p, ae_int_t *pcount, ae_state *_state);
void mlpimporttunableparameters(multilayerperceptron *network, RVector *p, ae_state *_state);
void mlpserializeold(multilayerperceptron *network, RVector *ra, ae_int_t *rlen, ae_state *_state);
void mlpunserializeold(RVector *ra, multilayerperceptron *network, ae_state *_state);
void mlprandomize(multilayerperceptron *network, ae_state *_state);
void mlprandomizefull(multilayerperceptron *network, ae_state *_state);
void mlpinitpreprocessor(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, ae_state *_state);
void mlpinitpreprocessorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t ssize, ae_state *_state);
void mlpinitpreprocessorsubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize, ae_state *_state);
void mlpinitpreprocessorsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize, ae_state *_state);
void mlpproperties(multilayerperceptron *network, ae_int_t *nin, ae_int_t *nout, ae_int_t *wcount, ae_state *_state);
ae_int_t mlpntotal(multilayerperceptron *network, ae_state *_state);
ae_int_t mlpgetinputscount(multilayerperceptron *network, ae_state *_state);
ae_int_t mlpgetoutputscount(multilayerperceptron *network, ae_state *_state);
ae_int_t mlpgetweightscount(multilayerperceptron *network, ae_state *_state);
bool mlpissoftmax(multilayerperceptron *network, ae_state *_state);
ae_int_t mlpgetlayerscount(multilayerperceptron *network, ae_state *_state);
ae_int_t mlpgetlayersize(multilayerperceptron *network, ae_int_t k, ae_state *_state);
void mlpgetinputscaling(multilayerperceptron *network, ae_int_t i, double *mean, double *sigma, ae_state *_state);
void mlpgetoutputscaling(multilayerperceptron *network, ae_int_t i, double *mean, double *sigma, ae_state *_state);
void mlpgetneuroninfo(multilayerperceptron *network, ae_int_t k, ae_int_t i, ae_int_t *fkind, double *threshold, ae_state *_state);
double mlpgetweight(multilayerperceptron *network, ae_int_t k0, ae_int_t i0, ae_int_t k1, ae_int_t i1, ae_state *_state);
void mlpsetinputscaling(multilayerperceptron *network, ae_int_t i, double mean, double sigma, ae_state *_state);
void mlpsetoutputscaling(multilayerperceptron *network, ae_int_t i, double mean, double sigma, ae_state *_state);
void mlpsetneuroninfo(multilayerperceptron *network, ae_int_t k, ae_int_t i, ae_int_t fkind, double threshold, ae_state *_state);
void mlpsetweight(multilayerperceptron *network, ae_int_t k0, ae_int_t i0, ae_int_t k1, ae_int_t i1, double w, ae_state *_state);
void mlpactivationfunction(double net, ae_int_t k, double *f, double *df, double *d2f, ae_state *_state);
void mlpprocess(multilayerperceptron *network, RVector *x, RVector *y, ae_state *_state);
void mlpprocessi(multilayerperceptron *network, RVector *x, RVector *y, ae_state *_state);
double mlperror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlperrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
double mlperrorn(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, ae_state *_state);
ae_int_t mlpclserror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlprelclserror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlprelclserrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
double mlpavgce(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpavgcesparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
double mlprmserror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlprmserrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
double mlpavgerror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpavgerrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
double mlpavgrelerror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpavgrelerrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
void mlpgrad(multilayerperceptron *network, RVector *x, RVector *desiredy, double *e, RVector *grad, ae_state *_state);
void mlpgradn(multilayerperceptron *network, RVector *x, RVector *desiredy, double *e, RVector *grad, ae_state *_state);
void mlpgradbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad, ae_state *_state);
void mlpgradbatchsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t ssize, double *e, RVector *grad, ae_state *_state);
void mlpgradbatchsubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize, double *e, RVector *grad, ae_state *_state);
void mlpgradbatchsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize, double *e, RVector *grad, ae_state *_state);
void mlpgradbatchx(multilayerperceptron *network, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, ae_shared_pool *gradbuf, ae_state *_state);
bool _trypexec_mlpgradbatchx(multilayerperceptron *network, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, ae_shared_pool *gradbuf, ae_state *_state);
void mlpgradnbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad, ae_state *_state);
void mlphessiannbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad, RMatrix *h, ae_state *_state);
void mlphessianbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad, RMatrix *h, ae_state *_state);
void mlpinternalprocessvector(ZVector *structinfo, RVector *weights, RVector *columnmeans, RVector *columnsigmas, RVector *neurons, RVector *dfdnet, RVector *x, RVector *y, ae_state *_state);
void mlpalloc(ae_serializer *s, multilayerperceptron *network, ae_state *_state);
void mlpserialize(ae_serializer *s, multilayerperceptron *network, ae_state *_state);
void mlpunserialize(ae_serializer *s, multilayerperceptron *network, ae_state *_state);
void mlpallerrorssubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize, modelerrors *rep, ae_state *_state);
void mlpallerrorssparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize, modelerrors *rep, ae_state *_state);
double mlperrorsubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize, ae_state *_state);
double mlperrorsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize, ae_state *_state);
void mlpallerrorsx(multilayerperceptron *network, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep, ae_state *_state);
bool _trypexec_mlpallerrorsx(multilayerperceptron *network, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep, ae_state *_state);
void modelerrors_init(void *_p, ae_state *_state, bool make_automatic);
void modelerrors_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void modelerrors_free(void *_p, bool make_automatic);
void smlpgrad_init(void *_p, ae_state *_state, bool make_automatic);
void smlpgrad_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void smlpgrad_free(void *_p, bool make_automatic);
void multilayerperceptron_init(void *_p, ae_state *_state, bool make_automatic);
void multilayerperceptron_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void multilayerperceptron_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(modelerrors, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);
DecClass(multilayerperceptron, );

void mlpserialize(multilayerperceptron &obj, std::string &s_out);
void mlpunserialize(const std::string &s_in, multilayerperceptron &obj);
void mlpserialize(multilayerperceptron &obj, std::ostream &s_out);
void mlpunserialize(const std::istream &s_in, multilayerperceptron &obj);
void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpcopy(const multilayerperceptron &network1, multilayerperceptron &network2, const xparams _xparams = xdefault);
void mlpcopytunableparameters(const multilayerperceptron &network1, const multilayerperceptron &network2, const xparams _xparams = xdefault);
void mlprandomize(const multilayerperceptron &network, const xparams _xparams = xdefault);
void mlprandomizefull(const multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpinitpreprocessor(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, const xparams _xparams = xdefault);
void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount, const xparams _xparams = xdefault);
ae_int_t mlpgetinputscount(const multilayerperceptron &network, const xparams _xparams = xdefault);
ae_int_t mlpgetoutputscount(const multilayerperceptron &network, const xparams _xparams = xdefault);
ae_int_t mlpgetweightscount(const multilayerperceptron &network, const xparams _xparams = xdefault);
bool mlpissoftmax(const multilayerperceptron &network, const xparams _xparams = xdefault);
ae_int_t mlpgetlayerscount(const multilayerperceptron &network, const xparams _xparams = xdefault);
ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k, const xparams _xparams = xdefault);
void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma, const xparams _xparams = xdefault);
void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma, const xparams _xparams = xdefault);
void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold, const xparams _xparams = xdefault);
double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const xparams _xparams = xdefault);
void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma, const xparams _xparams = xdefault);
void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma, const xparams _xparams = xdefault);
void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold, const xparams _xparams = xdefault);
void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w, const xparams _xparams = xdefault);
void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f, const xparams _xparams = xdefault);
void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, const xparams _xparams = xdefault);
ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, const xparams _xparams = xdefault);
void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h, const xparams _xparams = xdefault);
void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h, const xparams _xparams = xdefault);
void mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep, const xparams _xparams = xdefault);
void mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep, const xparams _xparams = xdefault);
double mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, const xparams _xparams = xdefault);
double mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, const xparams _xparams = xdefault);
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

void mlpecreate0(ae_int_t nin, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreate1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreate2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreateb0(ae_int_t nin, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreateb1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreateb2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreater0(ae_int_t nin, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreater1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreater2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreatec0(ae_int_t nin, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreatec1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreatec2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecreatefromnetwork(multilayerperceptron *network, ae_int_t ensemblesize, mlpensemble *ensemble, ae_state *_state);
void mlpecopy(mlpensemble *ensemble1, mlpensemble *ensemble2, ae_state *_state);
void mlperandomize(mlpensemble *ensemble, ae_state *_state);
void mlpeproperties(mlpensemble *ensemble, ae_int_t *nin, ae_int_t *nout, ae_state *_state);
bool mlpeissoftmax(mlpensemble *ensemble, ae_state *_state);
void mlpeprocess(mlpensemble *ensemble, RVector *x, RVector *y, ae_state *_state);
void mlpeprocessi(mlpensemble *ensemble, RVector *x, RVector *y, ae_state *_state);
void mlpeallerrorsx(mlpensemble *ensemble, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep, ae_state *_state);
void mlpeallerrorssparse(mlpensemble *ensemble, sparsematrix *xy, ae_int_t npoints, double *relcls, double *avgce, double *rms, double *avg, double *avgrel, ae_state *_state);
double mlperelclserror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpeavgce(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpermserror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpeavgerror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mlpeavgrelerror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, ae_state *_state);
void mlpealloc(ae_serializer *s, mlpensemble *ensemble, ae_state *_state);
void mlpeserialize(ae_serializer *s, mlpensemble *ensemble, ae_state *_state);
void mlpeunserialize(ae_serializer *s, mlpensemble *ensemble, ae_state *_state);
void mlpensemble_init(void *_p, ae_state *_state, bool make_automatic);
void mlpensemble_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlpensemble_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mlpensemble, );

void mlpeserialize(mlpensemble &obj, std::string &s_out);
void mlpeunserialize(const std::string &s_in, mlpensemble &obj);
void mlpeserialize(mlpensemble &obj, std::ostream &s_out);
void mlpeunserialize(const std::istream &s_in, mlpensemble &obj);
void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlperandomize(const mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout, const xparams _xparams = xdefault);
bool mlpeissoftmax(const mlpensemble &ensemble, const xparams _xparams = xdefault);
void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
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
typedef struct {
   ae_int_t terminationtype;
   ae_int_t npoints;
   ae_vector p;
   ae_matrix z;
   ae_matrix pz;
   ae_matrix pm;
   ae_vector mergedist;
} ahcreport;
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

void clusterizercreate(clusterizerstate *s, ae_state *_state);
void clusterizersetpoints(clusterizerstate *s, RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, ae_state *_state);
void clusterizersetdistances(clusterizerstate *s, RMatrix *d, ae_int_t npoints, bool isupper, ae_state *_state);
void clusterizersetahcalgo(clusterizerstate *s, ae_int_t algo, ae_state *_state);
void clusterizersetkmeanslimits(clusterizerstate *s, ae_int_t restarts, ae_int_t maxits, ae_state *_state);
void clusterizersetkmeansinit(clusterizerstate *s, ae_int_t initalgo, ae_state *_state);
void clusterizersetseed(clusterizerstate *s, ae_int_t seed, ae_state *_state);
void clusterizerrunahc(clusterizerstate *s, ahcreport *rep, ae_state *_state);
void clusterizerrunkmeans(clusterizerstate *s, ae_int_t k, kmeansreport *rep, ae_state *_state);
void clusterizergetdistances(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, RMatrix *d, ae_state *_state);
void clusterizergetdistancesbuf(apbuffers *buf, RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, RMatrix *d, ae_state *_state);
void clusterizergetkclusters(ahcreport *rep, ae_int_t k, ZVector *cidx, ZVector *cz, ae_state *_state);
void clusterizerseparatedbydist(ahcreport *rep, double r, ae_int_t *k, ZVector *cidx, ZVector *cz, ae_state *_state);
void clusterizerseparatedbycorr(ahcreport *rep, double r, ae_int_t *k, ZVector *cidx, ZVector *cz, ae_state *_state);
void kmeansinitbuf(kmeansbuffers *buf, ae_state *_state);
void kmeansgenerateinternal(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t initalgo, ae_int_t seed, ae_int_t maxits, ae_int_t restarts, bool kmeansdbgnoits, ae_int_t *info, ae_int_t *iterationscount, RMatrix *ccol, bool needccol, RMatrix *crow, bool needcrow, ZVector *xyc, double *energy, kmeansbuffers *buf, ae_state *_state);
void kmeansupdatedistances(RMatrix *xy, ae_int_t idx0, ae_int_t idx1, ae_int_t nvars, RMatrix *ct, ae_int_t cidx0, ae_int_t cidx1, ZVector *xyc, RVector *xydist2, ae_shared_pool *bufferpool, ae_state *_state);
bool _trypexec_kmeansupdatedistances(RMatrix *xy, ae_int_t idx0, ae_int_t idx1, ae_int_t nvars, RMatrix *ct, ae_int_t cidx0, ae_int_t cidx1, ZVector *xyc, RVector *xydist2, ae_shared_pool *bufferpool, ae_state *_state);
void kmeansbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void kmeansbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void kmeansbuffers_free(void *_p, bool make_automatic);
void clusterizerstate_init(void *_p, ae_state *_state, bool make_automatic);
void clusterizerstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void clusterizerstate_free(void *_p, bool make_automatic);
void ahcreport_init(void *_p, ae_state *_state, bool make_automatic);
void ahcreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void ahcreport_free(void *_p, bool make_automatic);
void kmeansreport_init(void *_p, ae_state *_state, bool make_automatic);
void kmeansreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void kmeansreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(clusterizerstate, );
DecClass(ahcreport, ae_int_t &terminationtype; ae_int_t &npoints; integer_1d_array p; integer_2d_array z; integer_2d_array pz; integer_2d_array pm; real_1d_array mergedist;);
DecClass(kmeansreport, ae_int_t &npoints; ae_int_t &nfeatures; ae_int_t &terminationtype; ae_int_t &iterationscount; double &energy; ae_int_t &k; real_2d_array c; integer_1d_array cidx;);

void clusterizercreate(clusterizerstate &s, const xparams _xparams = xdefault);
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, const xparams _xparams = xdefault);
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t disttype, const xparams _xparams = xdefault);
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const ae_int_t npoints, const bool isupper, const xparams _xparams = xdefault);
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const bool isupper, const xparams _xparams = xdefault);
void clusterizersetahcalgo(const clusterizerstate &s, const ae_int_t algo, const xparams _xparams = xdefault);
void clusterizersetkmeanslimits(const clusterizerstate &s, const ae_int_t restarts, const ae_int_t maxits, const xparams _xparams = xdefault);
void clusterizersetkmeansinit(const clusterizerstate &s, const ae_int_t initalgo, const xparams _xparams = xdefault);
void clusterizersetseed(const clusterizerstate &s, const ae_int_t seed, const xparams _xparams = xdefault);
void clusterizerrunahc(const clusterizerstate &s, ahcreport &rep, const xparams _xparams = xdefault);
void clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep, const xparams _xparams = xdefault);
void clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d, const xparams _xparams = xdefault);
void clusterizergetkclusters(const ahcreport &rep, const ae_int_t k, integer_1d_array &cidx, integer_1d_array &cz, const xparams _xparams = xdefault);
void clusterizerseparatedbydist(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz, const xparams _xparams = xdefault);
void clusterizerseparatedbycorr(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz, const xparams _xparams = xdefault);
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
typedef struct {
   ae_vector trntotals;
   ae_vector oobtotals;
   ae_vector trncounts;
   ae_vector oobcounts;
   ae_vector giniimportances;
} dfvotebuf;
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
typedef struct {
   ae_vector treebuf;
   ae_int_t treeidx;
} dftreebuf;
typedef struct {
   ae_vector x;
   ae_vector y;
} decisionforestbuffer;
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

void dfcreatebuffer(decisionforest *model, decisionforestbuffer *buf, ae_state *_state);
void dfbuildercreate(decisionforestbuilder *s, ae_state *_state);
void dfbuildersetdataset(decisionforestbuilder *s, RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_state *_state);
void dfbuildersetrndvars(decisionforestbuilder *s, ae_int_t rndvars, ae_state *_state);
void dfbuildersetrndvarsratio(decisionforestbuilder *s, double f, ae_state *_state);
void dfbuildersetrndvarsauto(decisionforestbuilder *s, ae_state *_state);
void dfbuildersetsubsampleratio(decisionforestbuilder *s, double f, ae_state *_state);
void dfbuildersetseed(decisionforestbuilder *s, ae_int_t seedval, ae_state *_state);
void dfbuildersetrdfalgo(decisionforestbuilder *s, ae_int_t algotype, ae_state *_state);
void dfbuildersetrdfsplitstrength(decisionforestbuilder *s, ae_int_t splitstrength, ae_state *_state);
void dfbuildersetimportancetrngini(decisionforestbuilder *s, ae_state *_state);
void dfbuildersetimportanceoobgini(decisionforestbuilder *s, ae_state *_state);
void dfbuildersetimportancepermutation(decisionforestbuilder *s, ae_state *_state);
void dfbuildersetimportancenone(decisionforestbuilder *s, ae_state *_state);
double dfbuildergetprogress(decisionforestbuilder *s, ae_state *_state);
double dfbuilderpeekprogress(decisionforestbuilder *s, ae_state *_state);
void dfbuilderbuildrandomforest(decisionforestbuilder *s, ae_int_t ntrees, decisionforest *df, dfreport *rep, ae_state *_state);
double dfbinarycompression(decisionforest *df, ae_state *_state);
double dfbinarycompression8(decisionforest *df, ae_state *_state);
void dfprocess(decisionforest *df, RVector *x, RVector *y, ae_state *_state);
void dfprocessi(decisionforest *df, RVector *x, RVector *y, ae_state *_state);
double dfprocess0(decisionforest *model, RVector *x, ae_state *_state);
ae_int_t dfclassify(decisionforest *model, RVector *x, ae_state *_state);
void dftsprocess(decisionforest *df, decisionforestbuffer *buf, RVector *x, RVector *y, ae_state *_state);
double dfrelclserror(decisionforest *df, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double dfavgce(decisionforest *df, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double dfrmserror(decisionforest *df, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double dfavgerror(decisionforest *df, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double dfavgrelerror(decisionforest *df, RMatrix *xy, ae_int_t npoints, ae_state *_state);
void dfcopy(decisionforest *df1, decisionforest *df2, ae_state *_state);
void dfalloc(ae_serializer *s, decisionforest *forest, ae_state *_state);
void dfserialize(ae_serializer *s, decisionforest *forest, ae_state *_state);
void dfunserialize(ae_serializer *s, decisionforest *forest, ae_state *_state);
void dfbuildrandomdecisionforest(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, double r, ae_int_t *info, decisionforest *df, dfreport *rep, ae_state *_state);
void dfbuildrandomdecisionforestx1(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t nrndvars, double r, ae_int_t *info, decisionforest *df, dfreport *rep, ae_state *_state);
void dfbuildinternal(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t samplesize, ae_int_t nfeatures, ae_int_t flags, ae_int_t *info, decisionforest *df, dfreport *rep, ae_state *_state);
void decisionforestbuilder_init(void *_p, ae_state *_state, bool make_automatic);
void decisionforestbuilder_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void decisionforestbuilder_free(void *_p, bool make_automatic);
void dfworkbuf_init(void *_p, ae_state *_state, bool make_automatic);
void dfworkbuf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dfworkbuf_free(void *_p, bool make_automatic);
void dfvotebuf_init(void *_p, ae_state *_state, bool make_automatic);
void dfvotebuf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dfvotebuf_free(void *_p, bool make_automatic);
void dfpermimpbuf_init(void *_p, ae_state *_state, bool make_automatic);
void dfpermimpbuf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dfpermimpbuf_free(void *_p, bool make_automatic);
void dftreebuf_init(void *_p, ae_state *_state, bool make_automatic);
void dftreebuf_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dftreebuf_free(void *_p, bool make_automatic);
void decisionforestbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void decisionforestbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void decisionforestbuffer_free(void *_p, bool make_automatic);
void decisionforest_init(void *_p, ae_state *_state, bool make_automatic);
void decisionforest_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void decisionforest_free(void *_p, bool make_automatic);
void dfreport_init(void *_p, ae_state *_state, bool make_automatic);
void dfreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dfreport_free(void *_p, bool make_automatic);
void dfinternalbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void dfinternalbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dfinternalbuffers_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(decisionforestbuilder, );
DecClass(decisionforestbuffer, );
DecClass(decisionforest, );
DecClass(dfreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror; double &oobrelclserror; double &oobavgce; double &oobrmserror; double &oobavgerror; double &oobavgrelerror; integer_1d_array topvars; real_1d_array varimportances;);

void dfserialize(decisionforest &obj, std::string &s_out);
void dfunserialize(const std::string &s_in, decisionforest &obj);
void dfserialize(decisionforest &obj, std::ostream &s_out);
void dfunserialize(const std::istream &s_in, decisionforest &obj);
void dfcreatebuffer(const decisionforest &model, decisionforestbuffer &buf, const xparams _xparams = xdefault);
void dfbuildercreate(decisionforestbuilder &s, const xparams _xparams = xdefault);
void dfbuildersetdataset(const decisionforestbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const xparams _xparams = xdefault);
void dfbuildersetrndvars(const decisionforestbuilder &s, const ae_int_t rndvars, const xparams _xparams = xdefault);
void dfbuildersetrndvarsratio(const decisionforestbuilder &s, const double f, const xparams _xparams = xdefault);
void dfbuildersetrndvarsauto(const decisionforestbuilder &s, const xparams _xparams = xdefault);
void dfbuildersetsubsampleratio(const decisionforestbuilder &s, const double f, const xparams _xparams = xdefault);
void dfbuildersetseed(const decisionforestbuilder &s, const ae_int_t seedval, const xparams _xparams = xdefault);
void dfbuildersetrdfalgo(const decisionforestbuilder &s, const ae_int_t algotype, const xparams _xparams = xdefault);
void dfbuildersetrdfsplitstrength(const decisionforestbuilder &s, const ae_int_t splitstrength, const xparams _xparams = xdefault);
void dfbuildersetimportancetrngini(const decisionforestbuilder &s, const xparams _xparams = xdefault);
void dfbuildersetimportanceoobgini(const decisionforestbuilder &s, const xparams _xparams = xdefault);
void dfbuildersetimportancepermutation(const decisionforestbuilder &s, const xparams _xparams = xdefault);
void dfbuildersetimportancenone(const decisionforestbuilder &s, const xparams _xparams = xdefault);
double dfbuildergetprogress(const decisionforestbuilder &s, const xparams _xparams = xdefault);
double dfbuilderpeekprogress(const decisionforestbuilder &s, const xparams _xparams = xdefault);
void dfbuilderbuildrandomforest(const decisionforestbuilder &s, const ae_int_t ntrees, decisionforest &df, dfreport &rep, const xparams _xparams = xdefault);
double dfbinarycompression(const decisionforest &df, const xparams _xparams = xdefault);
void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
double dfprocess0(const decisionforest &model, const real_1d_array &x, const xparams _xparams = xdefault);
ae_int_t dfclassify(const decisionforest &model, const real_1d_array &x, const xparams _xparams = xdefault);
void dftsprocess(const decisionforest &df, const decisionforestbuffer &buf, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep, const xparams _xparams = xdefault);
void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LINREG Package ===
// Depends on: (SpecialFunctions) IGAMMAF
// Depends on: (LinAlg) SVD
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
typedef struct {
   ae_vector w;
} linearmodel;
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

void lrbuild(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar, ae_state *_state);
void lrbuilds(RMatrix *xy, RVector *s, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar, ae_state *_state);
void lrbuildzs(RMatrix *xy, RVector *s, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar, ae_state *_state);
void lrbuildz(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar, ae_state *_state);
void lrunpack(linearmodel *lm, RVector *v, ae_int_t *nvars, ae_state *_state);
void lrpack(RVector *v, ae_int_t nvars, linearmodel *lm, ae_state *_state);
double lrprocess(linearmodel *lm, RVector *x, ae_state *_state);
double lrrmserror(linearmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double lravgerror(linearmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double lravgrelerror(linearmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
void lrcopy(linearmodel *lm1, linearmodel *lm2, ae_state *_state);
void lrlines(RMatrix *xy, RVector *s, ae_int_t n, ae_int_t *info, double *a, double *b, double *vara, double *varb, double *covab, double *corrab, double *p, ae_state *_state);
void lrline(RMatrix *xy, ae_int_t n, ae_int_t *info, double *a, double *b, ae_state *_state);
void linearmodel_init(void *_p, ae_state *_state, bool make_automatic);
void linearmodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void linearmodel_free(void *_p, bool make_automatic);
void lrreport_init(void *_p, ae_state *_state, bool make_automatic);
void lrreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void lrreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(linearmodel, );
DecClass(lrreport, real_2d_array c; double &rmserror; double &avgerror; double &avgrelerror; double &cvrmserror; double &cvavgerror; double &cvavgrelerror; ae_int_t &ncvdefects; integer_1d_array cvdefects;);

void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar, const xparams _xparams = xdefault);
void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar, const xparams _xparams = xdefault);
void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar, const xparams _xparams = xdefault);
void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar, const xparams _xparams = xdefault);
void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars, const xparams _xparams = xdefault);
void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm, const xparams _xparams = xdefault);
double lrprocess(const linearmodel &lm, const real_1d_array &x, const xparams _xparams = xdefault);
double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
} // end of namespace alglib

// === FILTERS Package ===
// Depends on: LINREG
namespace alglib_impl {
void filtersma(RVector *x, ae_int_t n, ae_int_t k, ae_state *_state);
void filterema(RVector *x, ae_int_t n, double alpha, ae_state *_state);
void filterlrma(RVector *x, ae_int_t n, ae_int_t k, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k, const xparams _xparams = xdefault);
void filtersma(real_1d_array &x, const ae_int_t k, const xparams _xparams = xdefault);
void filterema(real_1d_array &x, const ae_int_t n, const double alpha, const xparams _xparams = xdefault);
void filterema(real_1d_array &x, const double alpha, const xparams _xparams = xdefault);
void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k, const xparams _xparams = xdefault);
void filterlrma(real_1d_array &x, const ae_int_t k, const xparams _xparams = xdefault);
} // end of namespace alglib

// === SSA Package ===
// Depends on: (LinAlg) EVD, SVD
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

void ssacreate(ssamodel *s, ae_state *_state);
void ssasetwindow(ssamodel *s, ae_int_t windowwidth, ae_state *_state);
void ssasetseed(ssamodel *s, ae_int_t seed, ae_state *_state);
void ssasetpoweruplength(ssamodel *s, ae_int_t pwlen, ae_state *_state);
void ssasetmemorylimit(ssamodel *s, ae_int_t memlimit, ae_state *_state);
void ssaaddsequence(ssamodel *s, RVector *x, ae_int_t n, ae_state *_state);
void ssaappendpointandupdate(ssamodel *s, double x, double updateits, ae_state *_state);
void ssaappendsequenceandupdate(ssamodel *s, RVector *x, ae_int_t nticks, double updateits, ae_state *_state);
void ssasetalgoprecomputed(ssamodel *s, RMatrix *a, ae_int_t windowwidth, ae_int_t nbasis, ae_state *_state);
void ssasetalgotopkdirect(ssamodel *s, ae_int_t topk, ae_state *_state);
void ssasetalgotopkrealtime(ssamodel *s, ae_int_t topk, ae_state *_state);
void ssacleardata(ssamodel *s, ae_state *_state);
void ssagetbasis(ssamodel *s, RMatrix *a, RVector *sv, ae_int_t *windowwidth, ae_int_t *nbasis, ae_state *_state);
void ssagetlrr(ssamodel *s, RVector *a, ae_int_t *windowwidth, ae_state *_state);
void ssaanalyzelastwindow(ssamodel *s, RVector *trend, RVector *noise, ae_int_t *nticks, ae_state *_state);
void ssaanalyzelast(ssamodel *s, ae_int_t nticks, RVector *trend, RVector *noise, ae_state *_state);
void ssaanalyzesequence(ssamodel *s, RVector *data, ae_int_t nticks, RVector *trend, RVector *noise, ae_state *_state);
void ssaforecastlast(ssamodel *s, ae_int_t nticks, RVector *trend, ae_state *_state);
void ssaforecastsequence(ssamodel *s, RVector *data, ae_int_t datalen, ae_int_t forecastlen, bool applysmoothing, RVector *trend, ae_state *_state);
void ssaforecastavglast(ssamodel *s, ae_int_t m, ae_int_t nticks, RVector *trend, ae_state *_state);
void ssaforecastavgsequence(ssamodel *s, RVector *data, ae_int_t datalen, ae_int_t m, ae_int_t forecastlen, bool applysmoothing, RVector *trend, ae_state *_state);
void ssamodel_init(void *_p, ae_state *_state, bool make_automatic);
void ssamodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void ssamodel_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(ssamodel, );

void ssacreate(ssamodel &s, const xparams _xparams = xdefault);
void ssasetwindow(const ssamodel &s, const ae_int_t windowwidth, const xparams _xparams = xdefault);
void ssasetseed(const ssamodel &s, const ae_int_t seed, const xparams _xparams = xdefault);
void ssasetpoweruplength(const ssamodel &s, const ae_int_t pwlen, const xparams _xparams = xdefault);
void ssasetmemorylimit(const ssamodel &s, const ae_int_t memlimit, const xparams _xparams = xdefault);
void ssaaddsequence(const ssamodel &s, const real_1d_array &x, const ae_int_t n, const xparams _xparams = xdefault);
void ssaaddsequence(const ssamodel &s, const real_1d_array &x, const xparams _xparams = xdefault);
void ssaappendpointandupdate(const ssamodel &s, const double x, const double updateits, const xparams _xparams = xdefault);
void ssaappendsequenceandupdate(const ssamodel &s, const real_1d_array &x, const ae_int_t nticks, const double updateits, const xparams _xparams = xdefault);
void ssaappendsequenceandupdate(const ssamodel &s, const real_1d_array &x, const double updateits, const xparams _xparams = xdefault);
void ssasetalgoprecomputed(const ssamodel &s, const real_2d_array &a, const ae_int_t windowwidth, const ae_int_t nbasis, const xparams _xparams = xdefault);
void ssasetalgoprecomputed(const ssamodel &s, const real_2d_array &a, const xparams _xparams = xdefault);
void ssasetalgotopkdirect(const ssamodel &s, const ae_int_t topk, const xparams _xparams = xdefault);
void ssasetalgotopkrealtime(const ssamodel &s, const ae_int_t topk, const xparams _xparams = xdefault);
void ssacleardata(const ssamodel &s, const xparams _xparams = xdefault);
void ssagetbasis(const ssamodel &s, real_2d_array &a, real_1d_array &sv, ae_int_t &windowwidth, ae_int_t &nbasis, const xparams _xparams = xdefault);
void ssagetlrr(const ssamodel &s, real_1d_array &a, ae_int_t &windowwidth, const xparams _xparams = xdefault);
void ssaanalyzelastwindow(const ssamodel &s, real_1d_array &trend, real_1d_array &noise, ae_int_t &nticks, const xparams _xparams = xdefault);
void ssaanalyzelast(const ssamodel &s, const ae_int_t nticks, real_1d_array &trend, real_1d_array &noise, const xparams _xparams = xdefault);
void ssaanalyzesequence(const ssamodel &s, const real_1d_array &data, const ae_int_t nticks, real_1d_array &trend, real_1d_array &noise, const xparams _xparams = xdefault);
void ssaanalyzesequence(const ssamodel &s, const real_1d_array &data, real_1d_array &trend, real_1d_array &noise, const xparams _xparams = xdefault);
void ssaforecastlast(const ssamodel &s, const ae_int_t nticks, real_1d_array &trend, const xparams _xparams = xdefault);
void ssaforecastsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t datalen, const ae_int_t forecastlen, const bool applysmoothing, real_1d_array &trend, const xparams _xparams = xdefault);
void ssaforecastsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t forecastlen, real_1d_array &trend, const xparams _xparams = xdefault);
void ssaforecastavglast(const ssamodel &s, const ae_int_t m, const ae_int_t nticks, real_1d_array &trend, const xparams _xparams = xdefault);
void ssaforecastavgsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t datalen, const ae_int_t m, const ae_int_t forecastlen, const bool applysmoothing, real_1d_array &trend, const xparams _xparams = xdefault);
void ssaforecastavgsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t m, const ae_int_t forecastlen, real_1d_array &trend, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LDA Package ===
// Depends on: (LinAlg) EVD, MATINV
namespace alglib_impl {
void fisherlda(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, RVector *w, ae_state *_state);
void fisherldan(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, RMatrix *w, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w, const xparams _xparams = xdefault);
void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w, const xparams _xparams = xdefault);
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
typedef struct {
   ae_int_t inneriterationscount;
   ae_int_t outeriterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
} mcpdreport;

void mcpdcreate(ae_int_t n, mcpdstate *s, ae_state *_state);
void mcpdcreateentry(ae_int_t n, ae_int_t entrystate, mcpdstate *s, ae_state *_state);
void mcpdcreateexit(ae_int_t n, ae_int_t exitstate, mcpdstate *s, ae_state *_state);
void mcpdcreateentryexit(ae_int_t n, ae_int_t entrystate, ae_int_t exitstate, mcpdstate *s, ae_state *_state);
void mcpdaddtrack(mcpdstate *s, RMatrix *xy, ae_int_t k, ae_state *_state);
void mcpdsetec(mcpdstate *s, RMatrix *ec, ae_state *_state);
void mcpdaddec(mcpdstate *s, ae_int_t i, ae_int_t j, double c, ae_state *_state);
void mcpdsetbc(mcpdstate *s, RMatrix *bndl, RMatrix *bndu, ae_state *_state);
void mcpdaddbc(mcpdstate *s, ae_int_t i, ae_int_t j, double bndl, double bndu, ae_state *_state);
void mcpdsetlc(mcpdstate *s, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void mcpdsettikhonovregularizer(mcpdstate *s, double v, ae_state *_state);
void mcpdsetprior(mcpdstate *s, RMatrix *pp, ae_state *_state);
void mcpdsetpredictionweights(mcpdstate *s, RVector *pw, ae_state *_state);
void mcpdsolve(mcpdstate *s, ae_state *_state);
void mcpdresults(mcpdstate *s, RMatrix *p, mcpdreport *rep, ae_state *_state);
void mcpdstate_init(void *_p, ae_state *_state, bool make_automatic);
void mcpdstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mcpdstate_free(void *_p, bool make_automatic);
void mcpdreport_init(void *_p, ae_state *_state, bool make_automatic);
void mcpdreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mcpdreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mcpdstate, );
DecClass(mcpdreport, ae_int_t &inneriterationscount; ae_int_t &outeriterationscount; ae_int_t &nfev; ae_int_t &terminationtype;);

void mcpdcreate(const ae_int_t n, mcpdstate &s, const xparams _xparams = xdefault);
void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s, const xparams _xparams = xdefault);
void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s, const xparams _xparams = xdefault);
void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s, const xparams _xparams = xdefault);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k, const xparams _xparams = xdefault);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const xparams _xparams = xdefault);
void mcpdsetec(const mcpdstate &s, const real_2d_array &ec, const xparams _xparams = xdefault);
void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c, const xparams _xparams = xdefault);
void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu, const xparams _xparams = xdefault);
void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu, const xparams _xparams = xdefault);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = xdefault);
void mcpdsettikhonovregularizer(const mcpdstate &s, const double v, const xparams _xparams = xdefault);
void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp, const xparams _xparams = xdefault);
void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw, const xparams _xparams = xdefault);
void mcpdsolve(const mcpdstate &s, const xparams _xparams = xdefault);
void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === LOGIT Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS
// Depends on: MLPBASE
namespace alglib_impl {
typedef struct {
   ae_vector w;
} logitmodel;
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
typedef struct {
   ae_int_t ngrad;
   ae_int_t nhess;
} mnlreport;

void mnltrainh(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, logitmodel *lm, mnlreport *rep, ae_state *_state);
void mnlprocess(logitmodel *lm, RVector *x, RVector *y, ae_state *_state);
void mnlprocessi(logitmodel *lm, RVector *x, RVector *y, ae_state *_state);
void mnlunpack(logitmodel *lm, RMatrix *a, ae_int_t *nvars, ae_int_t *nclasses, ae_state *_state);
void mnlpack(RMatrix *a, ae_int_t nvars, ae_int_t nclasses, logitmodel *lm, ae_state *_state);
void mnlcopy(logitmodel *lm1, logitmodel *lm2, ae_state *_state);
double mnlavgce(logitmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mnlrelclserror(logitmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mnlrmserror(logitmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mnlavgerror(logitmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double mnlavgrelerror(logitmodel *lm, RMatrix *xy, ae_int_t ssize, ae_state *_state);
ae_int_t mnlclserror(logitmodel *lm, RMatrix *xy, ae_int_t npoints, ae_state *_state);
void logitmodel_init(void *_p, ae_state *_state, bool make_automatic);
void logitmodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void logitmodel_free(void *_p, bool make_automatic);
void logitmcstate_init(void *_p, ae_state *_state, bool make_automatic);
void logitmcstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void logitmcstate_free(void *_p, bool make_automatic);
void mnlreport_init(void *_p, ae_state *_state, bool make_automatic);
void mnlreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mnlreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(logitmodel, );
DecClass(mnlreport, ae_int_t &ngrad; ae_int_t &nhess;);

void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep, const xparams _xparams = xdefault);
void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses, const xparams _xparams = xdefault);
void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm, const xparams _xparams = xdefault);
double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize, const xparams _xparams = xdefault);
ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
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
typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
} knnreport;

void knncreatebuffer(knnmodel *model, knnbuffer *buf, ae_state *_state);
void knnbuildercreate(knnbuilder *s, ae_state *_state);
void knnbuildersetdatasetreg(knnbuilder *s, RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nout, ae_state *_state);
void knnbuildersetdatasetcls(knnbuilder *s, RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_state *_state);
void knnbuildersetnorm(knnbuilder *s, ae_int_t nrmtype, ae_state *_state);
void knnbuilderbuildknnmodel(knnbuilder *s, ae_int_t k, double eps, knnmodel *model, knnreport *rep, ae_state *_state);
void knnrewritekeps(knnmodel *model, ae_int_t k, double eps, ae_state *_state);
void knnprocess(knnmodel *model, RVector *x, RVector *y, ae_state *_state);
double knnprocess0(knnmodel *model, RVector *x, ae_state *_state);
ae_int_t knnclassify(knnmodel *model, RVector *x, ae_state *_state);
void knnprocessi(knnmodel *model, RVector *x, RVector *y, ae_state *_state);
void knntsprocess(knnmodel *model, knnbuffer *buf, RVector *x, RVector *y, ae_state *_state);
double knnrelclserror(knnmodel *model, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double knnavgce(knnmodel *model, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double knnrmserror(knnmodel *model, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double knnavgerror(knnmodel *model, RMatrix *xy, ae_int_t npoints, ae_state *_state);
double knnavgrelerror(knnmodel *model, RMatrix *xy, ae_int_t npoints, ae_state *_state);
void knnallerrors(knnmodel *model, RMatrix *xy, ae_int_t npoints, knnreport *rep, ae_state *_state);
void knnalloc(ae_serializer *s, knnmodel *model, ae_state *_state);
void knnserialize(ae_serializer *s, knnmodel *model, ae_state *_state);
void knnunserialize(ae_serializer *s, knnmodel *model, ae_state *_state);
void knnbuffer_init(void *_p, ae_state *_state, bool make_automatic);
void knnbuffer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void knnbuffer_free(void *_p, bool make_automatic);
void knnbuilder_init(void *_p, ae_state *_state, bool make_automatic);
void knnbuilder_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void knnbuilder_free(void *_p, bool make_automatic);
void knnmodel_init(void *_p, ae_state *_state, bool make_automatic);
void knnmodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void knnmodel_free(void *_p, bool make_automatic);
void knnreport_init(void *_p, ae_state *_state, bool make_automatic);
void knnreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void knnreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(knnbuffer, );
DecClass(knnbuilder, );
DecClass(knnmodel, );
DecClass(knnreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);

void knnserialize(knnmodel &obj, std::string &s_out);
void knnunserialize(const std::string &s_in, knnmodel &obj);
void knnserialize(knnmodel &obj, std::ostream &s_out);
void knnunserialize(const std::istream &s_in, knnmodel &obj);
void knncreatebuffer(const knnmodel &model, knnbuffer &buf, const xparams _xparams = xdefault);
void knnbuildercreate(knnbuilder &s, const xparams _xparams = xdefault);
void knnbuildersetdatasetreg(const knnbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nout, const xparams _xparams = xdefault);
void knnbuildersetdatasetcls(const knnbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const xparams _xparams = xdefault);
void knnbuildersetnorm(const knnbuilder &s, const ae_int_t nrmtype, const xparams _xparams = xdefault);
void knnbuilderbuildknnmodel(const knnbuilder &s, const ae_int_t k, const double eps, knnmodel &model, knnreport &rep, const xparams _xparams = xdefault);
void knnrewritekeps(const knnmodel &model, const ae_int_t k, const double eps, const xparams _xparams = xdefault);
void knnprocess(const knnmodel &model, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
double knnprocess0(const knnmodel &model, const real_1d_array &x, const xparams _xparams = xdefault);
ae_int_t knnclassify(const knnmodel &model, const real_1d_array &x, const xparams _xparams = xdefault);
void knnprocessi(const knnmodel &model, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
void knntsprocess(const knnmodel &model, const knnbuffer &buf, const real_1d_array &x, real_1d_array &y, const xparams _xparams = xdefault);
double knnrelclserror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double knnavgce(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double knnrmserror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double knnavgerror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
double knnavgrelerror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
void knnallerrors(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, knnreport &rep, const xparams _xparams = xdefault);
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
typedef struct {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
} mlpcvreport;
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
   rcommstate rstate;
   ae_int_t algoused;
   ae_int_t minibatchsize;
   hqrndstate generator;
} smlptrnsession;
typedef struct {
   ae_vector trnsubset;
   ae_vector valsubset;
   ae_shared_pool mlpsessions;
   mlpreport mlprep;
   multilayerperceptron network;
} mlpetrnsession;
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

void mlptrainlm(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep, ae_state *_state);
void mlptrainlbfgs(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t *info, mlpreport *rep, ae_state *_state);
void mlptraines(multilayerperceptron *network, RMatrix *trnxy, ae_int_t trnsize, RMatrix *valxy, ae_int_t valsize, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep, ae_state *_state);
void mlpkfoldcvlbfgs(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t foldscount, ae_int_t *info, mlpreport *rep, mlpcvreport *cvrep, ae_state *_state);
void mlpkfoldcvlm(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t foldscount, ae_int_t *info, mlpreport *rep, mlpcvreport *cvrep, ae_state *_state);
void mlpkfoldcv(mlptrainer *s, multilayerperceptron *network, ae_int_t nrestarts, ae_int_t foldscount, mlpreport *rep, ae_state *_state);
void mlpcreatetrainer(ae_int_t nin, ae_int_t nout, mlptrainer *s, ae_state *_state);
void mlpcreatetrainercls(ae_int_t nin, ae_int_t nclasses, mlptrainer *s, ae_state *_state);
void mlpsetdataset(mlptrainer *s, RMatrix *xy, ae_int_t npoints, ae_state *_state);
void mlpsetsparsedataset(mlptrainer *s, sparsematrix *xy, ae_int_t npoints, ae_state *_state);
void mlpsetdecay(mlptrainer *s, double decay, ae_state *_state);
void mlpsetcond(mlptrainer *s, double wstep, ae_int_t maxits, ae_state *_state);
void mlpsetalgobatch(mlptrainer *s, ae_state *_state);
void mlptrainnetwork(mlptrainer *s, multilayerperceptron *network, ae_int_t nrestarts, mlpreport *rep, ae_state *_state);
void mlpstarttraining(mlptrainer *s, multilayerperceptron *network, bool randomstart, ae_state *_state);
bool mlpcontinuetraining(mlptrainer *s, multilayerperceptron *network, ae_state *_state);
void mlpebagginglm(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep, mlpcvreport *ooberrors, ae_state *_state);
void mlpebagginglbfgs(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t *info, mlpreport *rep, mlpcvreport *ooberrors, ae_state *_state);
void mlpetraines(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep, ae_state *_state);
void mlptrainensemblees(mlptrainer *s, mlpensemble *ensemble, ae_int_t nrestarts, mlpreport *rep, ae_state *_state);
void mlpreport_init(void *_p, ae_state *_state, bool make_automatic);
void mlpreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlpreport_free(void *_p, bool make_automatic);
void mlpcvreport_init(void *_p, ae_state *_state, bool make_automatic);
void mlpcvreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlpcvreport_free(void *_p, bool make_automatic);
void smlptrnsession_init(void *_p, ae_state *_state, bool make_automatic);
void smlptrnsession_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void smlptrnsession_free(void *_p, bool make_automatic);
void mlpetrnsession_init(void *_p, ae_state *_state, bool make_automatic);
void mlpetrnsession_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlpetrnsession_free(void *_p, bool make_automatic);
void mlptrainer_init(void *_p, ae_state *_state, bool make_automatic);
void mlptrainer_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlptrainer_free(void *_p, bool make_automatic);
void mlpparallelizationcv_init(void *_p, ae_state *_state, bool make_automatic);
void mlpparallelizationcv_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mlpparallelizationcv_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mlpreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror; ae_int_t &ngrad; ae_int_t &nhess; ae_int_t &ncholesky;);
DecClass(mlpcvreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);
DecClass(mlptrainer, );

void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, const xparams _xparams = xdefault);
void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, const xparams _xparams = xdefault);
void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, const xparams _xparams = xdefault);
void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep, const xparams _xparams = xdefault);
void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep, const xparams _xparams = xdefault);
void mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep, const xparams _xparams = xdefault);
void mlpcreatetrainer(const ae_int_t nin, const ae_int_t nout, mlptrainer &s, const xparams _xparams = xdefault);
void mlpcreatetrainercls(const ae_int_t nin, const ae_int_t nclasses, mlptrainer &s, const xparams _xparams = xdefault);
void mlpsetdataset(const mlptrainer &s, const real_2d_array &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
void mlpsetsparsedataset(const mlptrainer &s, const sparsematrix &xy, const ae_int_t npoints, const xparams _xparams = xdefault);
void mlpsetdecay(const mlptrainer &s, const double decay, const xparams _xparams = xdefault);
void mlpsetcond(const mlptrainer &s, const double wstep, const ae_int_t maxits, const xparams _xparams = xdefault);
void mlpsetalgobatch(const mlptrainer &s, const xparams _xparams = xdefault);
void mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep, const xparams _xparams = xdefault);
void mlpstarttraining(const mlptrainer &s, const multilayerperceptron &network, const bool randomstart, const xparams _xparams = xdefault);
bool mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network, const xparams _xparams = xdefault);
void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors, const xparams _xparams = xdefault);
void mlpebagginglbfgs(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors, const xparams _xparams = xdefault);
void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, const xparams _xparams = xdefault);
void mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === DATACOMP Package ===
// Depends on: CLUSTERING
namespace alglib_impl {
void kmeansgenerate(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t restarts, ae_int_t *info, RMatrix *c, ZVector *xyc, ae_state *_state);
} // end of namespace alglib_impl

namespace alglib {
void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc, const xparams _xparams = xdefault);
} // end of namespace alglib

#endif // OnceOnly
