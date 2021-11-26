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
void pcabuildbasis(RMatrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector *s2, RMatrix *v);
void pcatruncatedsubspace(RMatrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t nneeded, double eps, ae_int_t maxits, RVector *s2, RMatrix *v);
void pcatruncatedsubspacesparse(sparsematrix *x, ae_int_t npoints, ae_int_t nvars, ae_int_t nneeded, double eps, ae_int_t maxits, RVector *s2, RMatrix *v);
} // end of namespace alglib_impl

namespace alglib {
void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v);
void pcatruncatedsubspace(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nneeded, const double eps, const ae_int_t maxits, real_1d_array &s2, real_2d_array &v);
void pcatruncatedsubspacesparse(const sparsematrix &x, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nneeded, const double eps, const ae_int_t maxits, real_1d_array &s2, real_2d_array &v);
} // end of namespace alglib

// === BDSS Package ===
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
struct cvreport {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
};
void cvreport_init(void *_p, bool make_automatic);
void cvreport_copy(void *_dst, void *_src, bool make_automatic);
void cvreport_free(void *_p, bool make_automatic);

void dserrallocate(ae_int_t nclasses, RVector *buf);
void dserraccumulate(RVector *buf, RVector *y, RVector *desiredy);
void dserrfinish(RVector *buf);
void dsnormalize(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector *means, RVector *sigmas);
void dsnormalizec(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, RVector *means, RVector *sigmas);
double dsgetmeanmindistance(RMatrix *xy, ae_int_t npoints, ae_int_t nvars);
void dstie(RVector *a, ae_int_t n, ZVector *ties, ae_int_t *tiecount, ZVector *p1, ZVector *p2);
void dstiefasti(RVector *a, ZVector *b, ae_int_t n, ZVector *ties, ae_int_t *tiecount, RVector *bufr, ZVector *bufi);
void dsoptimalsplit2(RVector *a, ZVector *c, ae_int_t n, ae_int_t *info, double *threshold, double *pal, double *pbl, double *par, double *pbr, double *cve);
void dsoptimalsplit2fast(RVector *a, ZVector *c, ZVector *tiesbuf, ZVector *cntbuf, RVector *bufr, ZVector *bufi, ae_int_t n, ae_int_t nc, double alpha, ae_int_t *info, double *threshold, double *rms, double *cvrms);
void dssplitk(RVector *a, ZVector *c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t *info, RVector *thresholds, ae_int_t *ni, double *cve);
void dsoptimalsplitk(RVector *a, ZVector *c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t *info, RVector *thresholds, ae_int_t *ni, double *cve);
} // end of namespace alglib_impl

namespace alglib {
void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve);
void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms);
} // end of namespace alglib

// === MLPBASE Package ===
// Depends on: (AlgLibInternal) SCODES, HPCCORES
// Depends on: (LinAlg) SPARSE
// Depends on: BDSS
namespace alglib_impl {
struct modelerrors {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
};
void modelerrors_init(void *_p, bool make_automatic);
void modelerrors_copy(void *_dst, void *_src, bool make_automatic);
void modelerrors_free(void *_p, bool make_automatic);

struct smlpgrad {
   double f;
   ae_vector g;
};
void smlpgrad_init(void *_p, bool make_automatic);
void smlpgrad_copy(void *_dst, void *_src, bool make_automatic);
void smlpgrad_free(void *_p, bool make_automatic);

struct multilayerperceptron {
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
};
void multilayerperceptron_init(void *_p, bool make_automatic);
void multilayerperceptron_copy(void *_dst, void *_src, bool make_automatic);
void multilayerperceptron_free(void *_p, bool make_automatic);
void mlpalloc(ae_serializer *s, multilayerperceptron *network);
void mlpserialize(ae_serializer *s, multilayerperceptron *network);
void mlpunserialize(ae_serializer *s, multilayerperceptron *network);

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
void mlpexporttunableparameters(multilayerperceptron *network, RVector *p, ae_int_t *pcount);
void mlpimporttunableparameters(multilayerperceptron *network, RVector *p);
void mlpserializeold(multilayerperceptron *network, RVector *ra, ae_int_t *rlen);
void mlpunserializeold(RVector *ra, multilayerperceptron *network);
void mlprandomize(multilayerperceptron *network);
void mlprandomizefull(multilayerperceptron *network);
void mlpinitpreprocessor(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize);
void mlpinitpreprocessorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t ssize);
void mlpinitpreprocessorsubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize);
void mlpinitpreprocessorsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize);
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
void mlpprocess(multilayerperceptron *network, RVector *x, RVector *y);
void mlpprocessi(multilayerperceptron *network, RVector *x, RVector *y);
double mlperror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlperrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlperrorn(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize);
ae_int_t mlpclserror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlprelclserror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlprelclserrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlpavgce(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlpavgcesparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlprmserror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlprmserrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlpavgerror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlpavgerrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
double mlpavgrelerror(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints);
double mlpavgrelerrorsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t npoints);
void mlpgrad(multilayerperceptron *network, RVector *x, RVector *desiredy, double *e, RVector *grad);
void mlpgradn(multilayerperceptron *network, RVector *x, RVector *desiredy, double *e, RVector *grad);
void mlpgradbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad);
void mlpgradbatchsparse(multilayerperceptron *network, sparsematrix *xy, ae_int_t ssize, double *e, RVector *grad);
void mlpgradbatchsubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize, double *e, RVector *grad);
void mlpgradbatchsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *idx, ae_int_t subsetsize, double *e, RVector *grad);
void mlpgradbatchx(multilayerperceptron *network, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, ae_shared_pool *gradbuf);
void mlpgradnbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad);
void mlphessiannbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad, RMatrix *h);
void mlphessianbatch(multilayerperceptron *network, RMatrix *xy, ae_int_t ssize, double *e, RVector *grad, RMatrix *h);
void mlpinternalprocessvector(ZVector *structinfo, RVector *weights, RVector *columnmeans, RVector *columnsigmas, RVector *neurons, RVector *dfdnet, RVector *x, RVector *y);
void mlpallerrorssubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize, modelerrors *rep);
void mlpallerrorssparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize, modelerrors *rep);
double mlperrorsubset(multilayerperceptron *network, RMatrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize);
double mlperrorsparsesubset(multilayerperceptron *network, sparsematrix *xy, ae_int_t setsize, ZVector *subset, ae_int_t subsetsize);
void mlpallerrorsx(multilayerperceptron *network, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(modelerrors, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);
DecClass(multilayerperceptron, );
void mlpserialize(multilayerperceptron &obj, std::string &s_out);
void mlpserialize(multilayerperceptron &obj, std::ostream &s_out);
void mlpunserialize(const std::string &s_in, multilayerperceptron &obj);
void mlpunserialize(const std::istream &s_in, multilayerperceptron &obj);

void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);
void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);
void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);
void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);
void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);
void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);
void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);
void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);
void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);
void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);
void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);
void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);
void mlpcopy(const multilayerperceptron &network1, multilayerperceptron &network2);
void mlpcopytunableparameters(const multilayerperceptron &network1, const multilayerperceptron &network2);
void mlprandomize(const multilayerperceptron &network);
void mlprandomizefull(const multilayerperceptron &network);
void mlpinitpreprocessor(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);
void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount);
ae_int_t mlpgetinputscount(const multilayerperceptron &network);
ae_int_t mlpgetoutputscount(const multilayerperceptron &network);
ae_int_t mlpgetweightscount(const multilayerperceptron &network);
bool mlpissoftmax(const multilayerperceptron &network);
ae_int_t mlpgetlayerscount(const multilayerperceptron &network);
ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k);
void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);
void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);
void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold);
double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1);
void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);
void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);
void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold);
void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w);
void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f);
void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);
void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);
double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);
ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);
void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);
void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
void mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
void mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
void mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);
void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);
void mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
void mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
double mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
double mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
} // end of namespace alglib

// === MLPE Package ===
// Depends on: MLPBASE
namespace alglib_impl {
struct mlpensemble {
   ae_int_t ensemblesize;
   ae_vector weights;
   ae_vector columnmeans;
   ae_vector columnsigmas;
   multilayerperceptron network;
   ae_vector y;
};
void mlpensemble_init(void *_p, bool make_automatic);
void mlpensemble_copy(void *_dst, void *_src, bool make_automatic);
void mlpensemble_free(void *_p, bool make_automatic);
void mlpealloc(ae_serializer *s, mlpensemble *ensemble);
void mlpeserialize(ae_serializer *s, mlpensemble *ensemble);
void mlpeunserialize(ae_serializer *s, mlpensemble *ensemble);

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
void mlpeprocess(mlpensemble *ensemble, RVector *x, RVector *y);
void mlpeprocessi(mlpensemble *ensemble, RVector *x, RVector *y);
void mlpeallerrorsx(mlpensemble *ensemble, RMatrix *densexy, sparsematrix *sparsexy, ae_int_t datasetsize, ae_int_t datasettype, ZVector *idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool *buf, modelerrors *rep);
void mlpeallerrorssparse(mlpensemble *ensemble, sparsematrix *xy, ae_int_t npoints, double *relcls, double *avgce, double *rms, double *avg, double *avgrel);
double mlperelclserror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints);
double mlpeavgce(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints);
double mlpermserror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints);
double mlpeavgerror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints);
double mlpeavgrelerror(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mlpensemble, );
void mlpeserialize(mlpensemble &obj, std::string &s_out);
void mlpeserialize(mlpensemble &obj, std::ostream &s_out);
void mlpeunserialize(const std::string &s_in, mlpensemble &obj);
void mlpeunserialize(const std::istream &s_in, mlpensemble &obj);

void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble);
void mlperandomize(const mlpensemble &ensemble);
void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout);
bool mlpeissoftmax(const mlpensemble &ensemble);
void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);
void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);
double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === CLUSTERING Package ===
// Depends on: (AlgLibInternal) BLAS
// Depends on: (AlgLibMisc) HQRND
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
struct kmeansbuffers {
   ae_matrix ct;
   ae_matrix ctbest;
   ae_vector xycbest;
   ae_vector xycprev;
   ae_vector d2;
   ae_vector csizes;
   apbuffers initbuf;
   ae_shared_pool updatepool;
};
void kmeansbuffers_init(void *_p, bool make_automatic);
void kmeansbuffers_copy(void *_dst, void *_src, bool make_automatic);
void kmeansbuffers_free(void *_p, bool make_automatic);

struct clusterizerstate {
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
};
void clusterizerstate_init(void *_p, bool make_automatic);
void clusterizerstate_copy(void *_dst, void *_src, bool make_automatic);
void clusterizerstate_free(void *_p, bool make_automatic);

struct ahcreport {
   ae_int_t terminationtype;
   ae_int_t npoints;
   ae_vector p;
   ae_matrix z;
   ae_matrix pz;
   ae_matrix pm;
   ae_vector mergedist;
};
void ahcreport_init(void *_p, bool make_automatic);
void ahcreport_copy(void *_dst, void *_src, bool make_automatic);
void ahcreport_free(void *_p, bool make_automatic);

struct kmeansreport {
   ae_int_t npoints;
   ae_int_t nfeatures;
   ae_int_t terminationtype;
   ae_int_t iterationscount;
   double energy;
   ae_int_t k;
   ae_matrix c;
   ae_vector cidx;
};
void kmeansreport_init(void *_p, bool make_automatic);
void kmeansreport_copy(void *_dst, void *_src, bool make_automatic);
void kmeansreport_free(void *_p, bool make_automatic);

void clusterizercreate(clusterizerstate *s);
void clusterizersetpoints(clusterizerstate *s, RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype);
void clusterizersetdistances(clusterizerstate *s, RMatrix *d, ae_int_t npoints, bool isupper);
void clusterizersetahcalgo(clusterizerstate *s, ae_int_t algo);
void clusterizersetkmeanslimits(clusterizerstate *s, ae_int_t restarts, ae_int_t maxits);
void clusterizersetkmeansinit(clusterizerstate *s, ae_int_t initalgo);
void clusterizersetseed(clusterizerstate *s, ae_int_t seed);
void clusterizerrunahc(clusterizerstate *s, ahcreport *rep);
void clusterizerrunkmeans(clusterizerstate *s, ae_int_t k, kmeansreport *rep);
void clusterizergetdistances(RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, RMatrix *d);
void clusterizergetdistancesbuf(apbuffers *buf, RMatrix *xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, RMatrix *d);
void clusterizergetkclusters(ahcreport *rep, ae_int_t k, ZVector *cidx, ZVector *cz);
void clusterizerseparatedbydist(ahcreport *rep, double r, ae_int_t *k, ZVector *cidx, ZVector *cz);
void clusterizerseparatedbycorr(ahcreport *rep, double r, ae_int_t *k, ZVector *cidx, ZVector *cz);
void kmeansinitbuf(kmeansbuffers *buf);
void kmeansgenerateinternal(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t initalgo, ae_int_t seed, ae_int_t maxits, ae_int_t restarts, bool kmeansdbgnoits, ae_int_t *info, ae_int_t *iterationscount, RMatrix *ccol, bool needccol, RMatrix *crow, bool needcrow, ZVector *xyc, double *energy, kmeansbuffers *buf);
void kmeansupdatedistances(RMatrix *xy, ae_int_t idx0, ae_int_t idx1, ae_int_t nvars, RMatrix *ct, ae_int_t cidx0, ae_int_t cidx1, ZVector *xyc, RVector *xydist2, ae_shared_pool *bufferpool);
} // end of namespace alglib_impl

namespace alglib {
DecClass(clusterizerstate, );
DecClass(ahcreport, ae_int_t &terminationtype; ae_int_t &npoints; integer_1d_array p; integer_2d_array z; integer_2d_array pz; integer_2d_array pm; real_1d_array mergedist;);
DecClass(kmeansreport, ae_int_t &npoints; ae_int_t &nfeatures; ae_int_t &terminationtype; ae_int_t &iterationscount; double &energy; ae_int_t &k; real_2d_array c; integer_1d_array cidx;);

void clusterizercreate(clusterizerstate &s);
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype);
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t disttype);
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const ae_int_t npoints, const bool isupper);
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const bool isupper);
void clusterizersetahcalgo(const clusterizerstate &s, const ae_int_t algo);
void clusterizersetkmeanslimits(const clusterizerstate &s, const ae_int_t restarts, const ae_int_t maxits);
void clusterizersetkmeansinit(const clusterizerstate &s, const ae_int_t initalgo);
void clusterizersetseed(const clusterizerstate &s, const ae_int_t seed);
void clusterizerrunahc(const clusterizerstate &s, ahcreport &rep);
void clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep);
void clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d);
void clusterizergetkclusters(const ahcreport &rep, const ae_int_t k, integer_1d_array &cidx, integer_1d_array &cz);
void clusterizerseparatedbydist(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);
void clusterizerseparatedbycorr(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);
} // end of namespace alglib

// === DFOREST Package ===
// Depends on: (AlgLibInternal) SCODES
// Depends on: (AlgLibMisc) HQRND
// Depends on: BDSS
namespace alglib_impl {
struct decisionforestbuilder {
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
};
void decisionforestbuilder_init(void *_p, bool make_automatic);
void decisionforestbuilder_copy(void *_dst, void *_src, bool make_automatic);
void decisionforestbuilder_free(void *_p, bool make_automatic);

struct dfworkbuf {
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
};
void dfworkbuf_init(void *_p, bool make_automatic);
void dfworkbuf_copy(void *_dst, void *_src, bool make_automatic);
void dfworkbuf_free(void *_p, bool make_automatic);

struct dfvotebuf {
   ae_vector trntotals;
   ae_vector oobtotals;
   ae_vector trncounts;
   ae_vector oobcounts;
   ae_vector giniimportances;
};
void dfvotebuf_init(void *_p, bool make_automatic);
void dfvotebuf_copy(void *_dst, void *_src, bool make_automatic);
void dfvotebuf_free(void *_p, bool make_automatic);

struct dfpermimpbuf {
   ae_vector losses;
   ae_vector xraw;
   ae_vector xdist;
   ae_vector xcur;
   ae_vector y;
   ae_vector yv;
   ae_vector targety;
   ae_vector startnodes;
};
void dfpermimpbuf_init(void *_p, bool make_automatic);
void dfpermimpbuf_copy(void *_dst, void *_src, bool make_automatic);
void dfpermimpbuf_free(void *_p, bool make_automatic);

struct dftreebuf {
   ae_vector treebuf;
   ae_int_t treeidx;
};
void dftreebuf_init(void *_p, bool make_automatic);
void dftreebuf_copy(void *_dst, void *_src, bool make_automatic);
void dftreebuf_free(void *_p, bool make_automatic);

struct decisionforestbuffer {
   ae_vector x;
   ae_vector y;
};
void decisionforestbuffer_init(void *_p, bool make_automatic);
void decisionforestbuffer_copy(void *_dst, void *_src, bool make_automatic);
void decisionforestbuffer_free(void *_p, bool make_automatic);

struct decisionforest {
   ae_int_t forestformat;
   bool usemantissa8;
   ae_int_t nvars;
   ae_int_t nclasses;
   ae_int_t ntrees;
   ae_int_t bufsize;
   ae_vector trees;
   decisionforestbuffer buffer;
   ae_vector trees8;
};
void decisionforest_init(void *_p, bool make_automatic);
void decisionforest_copy(void *_dst, void *_src, bool make_automatic);
void decisionforest_free(void *_p, bool make_automatic);

struct dfreport {
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
};
void dfreport_init(void *_p, bool make_automatic);
void dfreport_copy(void *_dst, void *_src, bool make_automatic);
void dfreport_free(void *_p, bool make_automatic);

struct dfinternalbuffers {
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
};
void dfinternalbuffers_init(void *_p, bool make_automatic);
void dfinternalbuffers_copy(void *_dst, void *_src, bool make_automatic);
void dfinternalbuffers_free(void *_p, bool make_automatic);
void dfalloc(ae_serializer *s, decisionforest *forest);
void dfserialize(ae_serializer *s, decisionforest *forest);
void dfunserialize(ae_serializer *s, decisionforest *forest);

void dfcreatebuffer(decisionforest *model, decisionforestbuffer *buf);
void dfbuildrandomdecisionforest(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, double r, ae_int_t *info, decisionforest *df, dfreport *rep);
void dfbuildrandomdecisionforestx1(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t nrndvars, double r, ae_int_t *info, decisionforest *df, dfreport *rep);
void dfbuildercreate(decisionforestbuilder *s);
void dfbuildersetdataset(decisionforestbuilder *s, RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses);
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
double dfbuilderpeekprogress(decisionforestbuilder *s);
double dfbuildergetprogress(decisionforestbuilder *s);
void dfbuilderbuildrandomforest(decisionforestbuilder *s, ae_int_t ntrees, decisionforest *df, dfreport *rep);
double dfbinarycompression(decisionforest *df);
double dfbinarycompression8(decisionforest *df);
void dfprocess(decisionforest *df, RVector *x, RVector *y);
void dfprocessi(decisionforest *df, RVector *x, RVector *y);
double dfprocess0(decisionforest *model, RVector *x);
ae_int_t dfclassify(decisionforest *model, RVector *x);
void dftsprocess(decisionforest *df, decisionforestbuffer *buf, RVector *x, RVector *y);
double dfrelclserror(decisionforest *df, RMatrix *xy, ae_int_t npoints);
double dfavgce(decisionforest *df, RMatrix *xy, ae_int_t npoints);
double dfrmserror(decisionforest *df, RMatrix *xy, ae_int_t npoints);
double dfavgerror(decisionforest *df, RMatrix *xy, ae_int_t npoints);
double dfavgrelerror(decisionforest *df, RMatrix *xy, ae_int_t npoints);
void dfcopy(decisionforest *df1, decisionforest *df2);
void dfbuildinternal(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t samplesize, ae_int_t nfeatures, ae_int_t flags, ae_int_t *info, decisionforest *df, dfreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(decisionforestbuilder, );
DecClass(decisionforestbuffer, );
DecClass(decisionforest, );
DecClass(dfreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror; double &oobrelclserror; double &oobavgce; double &oobrmserror; double &oobavgerror; double &oobavgrelerror; integer_1d_array topvars; real_1d_array varimportances;);
void dfserialize(decisionforest &obj, std::string &s_out);
void dfserialize(decisionforest &obj, std::ostream &s_out);
void dfunserialize(const std::string &s_in, decisionforest &obj);
void dfunserialize(const std::istream &s_in, decisionforest &obj);

void dfcreatebuffer(const decisionforest &model, decisionforestbuffer &buf);
void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);
void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);
void dfbuildercreate(decisionforestbuilder &s);
void dfbuildersetdataset(const decisionforestbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses);
void dfbuildersetrndvars(const decisionforestbuilder &s, const ae_int_t rndvars);
void dfbuildersetrndvarsratio(const decisionforestbuilder &s, const double f);
void dfbuildersetrndvarsauto(const decisionforestbuilder &s);
void dfbuildersetsubsampleratio(const decisionforestbuilder &s, const double f);
void dfbuildersetseed(const decisionforestbuilder &s, const ae_int_t seedval);
void dfbuildersetrdfalgo(const decisionforestbuilder &s, const ae_int_t algotype);
void dfbuildersetrdfsplitstrength(const decisionforestbuilder &s, const ae_int_t splitstrength);
void dfbuildersetimportancetrngini(const decisionforestbuilder &s);
void dfbuildersetimportanceoobgini(const decisionforestbuilder &s);
void dfbuildersetimportancepermutation(const decisionforestbuilder &s);
void dfbuildersetimportancenone(const decisionforestbuilder &s);
double dfbuilderpeekprogress(const decisionforestbuilder &s);
double dfbuildergetprogress(const decisionforestbuilder &s);
void dfbuilderbuildrandomforest(const decisionforestbuilder &s, const ae_int_t ntrees, decisionforest &df, dfreport &rep);
double dfbinarycompression(const decisionforest &df);
void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y);
void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y);
double dfprocess0(const decisionforest &model, const real_1d_array &x);
ae_int_t dfclassify(const decisionforest &model, const real_1d_array &x);
void dftsprocess(const decisionforest &df, const decisionforestbuffer &buf, const real_1d_array &x, real_1d_array &y);
double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === LINREG Package ===
// Depends on: (SpecialFunctions) IGAMMAF
// Depends on: (LinAlg) SVD
// Depends on: (Statistics) BASESTAT
namespace alglib_impl {
struct linearmodel {
   ae_vector w;
};
void linearmodel_init(void *_p, bool make_automatic);
void linearmodel_copy(void *_dst, void *_src, bool make_automatic);
void linearmodel_free(void *_p, bool make_automatic);

struct lrreport {
   ae_matrix c;
   double rmserror;
   double avgerror;
   double avgrelerror;
   double cvrmserror;
   double cvavgerror;
   double cvavgrelerror;
   ae_int_t ncvdefects;
   ae_vector cvdefects;
};
void lrreport_init(void *_p, bool make_automatic);
void lrreport_copy(void *_dst, void *_src, bool make_automatic);
void lrreport_free(void *_p, bool make_automatic);

void lrbuild(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrbuilds(RMatrix *xy, RVector *s, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrbuildzs(RMatrix *xy, RVector *s, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrbuildz(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t *info, linearmodel *lm, lrreport *ar);
void lrunpack(linearmodel *lm, RVector *v, ae_int_t *nvars);
void lrpack(RVector *v, ae_int_t nvars, linearmodel *lm);
double lrprocess(linearmodel *lm, RVector *x);
double lrrmserror(linearmodel *lm, RMatrix *xy, ae_int_t npoints);
double lravgerror(linearmodel *lm, RMatrix *xy, ae_int_t npoints);
double lravgrelerror(linearmodel *lm, RMatrix *xy, ae_int_t npoints);
void lrcopy(linearmodel *lm1, linearmodel *lm2);
void lrlines(RMatrix *xy, RVector *s, ae_int_t n, ae_int_t *info, double *a, double *b, double *vara, double *varb, double *covab, double *corrab, double *p);
void lrline(RMatrix *xy, ae_int_t n, ae_int_t *info, double *a, double *b);
} // end of namespace alglib_impl

namespace alglib {
DecClass(linearmodel, );
DecClass(lrreport, real_2d_array c; double &rmserror; double &avgerror; double &avgrelerror; double &cvrmserror; double &cvavgerror; double &cvavgrelerror; ae_int_t &ncvdefects; integer_1d_array cvdefects;);

void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars);
void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm);
double lrprocess(const linearmodel &lm, const real_1d_array &x);
double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === FILTERS Package ===
// Depends on: LINREG
namespace alglib_impl {
void filtersma(RVector *x, ae_int_t n, ae_int_t k);
void filterema(RVector *x, ae_int_t n, double alpha);
void filterlrma(RVector *x, ae_int_t n, ae_int_t k);
} // end of namespace alglib_impl

namespace alglib {
void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filtersma(real_1d_array &x, const ae_int_t k);
void filterema(real_1d_array &x, const ae_int_t n, const double alpha);
void filterema(real_1d_array &x, const double alpha);
void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filterlrma(real_1d_array &x, const ae_int_t k);
} // end of namespace alglib

// === SSA Package ===
// Depends on: (LinAlg) SVD, EVD
namespace alglib_impl {
struct ssamodel {
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
};
void ssamodel_init(void *_p, bool make_automatic);
void ssamodel_copy(void *_dst, void *_src, bool make_automatic);
void ssamodel_free(void *_p, bool make_automatic);

void ssacreate(ssamodel *s);
void ssasetwindow(ssamodel *s, ae_int_t windowwidth);
void ssasetseed(ssamodel *s, ae_int_t seed);
void ssasetpoweruplength(ssamodel *s, ae_int_t pwlen);
void ssasetmemorylimit(ssamodel *s, ae_int_t memlimit);
void ssaaddsequence(ssamodel *s, RVector *x, ae_int_t n);
void ssaappendpointandupdate(ssamodel *s, double x, double updateits);
void ssaappendsequenceandupdate(ssamodel *s, RVector *x, ae_int_t nticks, double updateits);
void ssasetalgoprecomputed(ssamodel *s, RMatrix *a, ae_int_t windowwidth, ae_int_t nbasis);
void ssasetalgotopkdirect(ssamodel *s, ae_int_t topk);
void ssasetalgotopkrealtime(ssamodel *s, ae_int_t topk);
void ssacleardata(ssamodel *s);
void ssagetbasis(ssamodel *s, RMatrix *a, RVector *sv, ae_int_t *windowwidth, ae_int_t *nbasis);
void ssagetlrr(ssamodel *s, RVector *a, ae_int_t *windowwidth);
void ssaanalyzelastwindow(ssamodel *s, RVector *trend, RVector *noise, ae_int_t *nticks);
void ssaanalyzelast(ssamodel *s, ae_int_t nticks, RVector *trend, RVector *noise);
void ssaanalyzesequence(ssamodel *s, RVector *data, ae_int_t nticks, RVector *trend, RVector *noise);
void ssaforecastlast(ssamodel *s, ae_int_t nticks, RVector *trend);
void ssaforecastsequence(ssamodel *s, RVector *data, ae_int_t datalen, ae_int_t forecastlen, bool applysmoothing, RVector *trend);
void ssaforecastavglast(ssamodel *s, ae_int_t m, ae_int_t nticks, RVector *trend);
void ssaforecastavgsequence(ssamodel *s, RVector *data, ae_int_t datalen, ae_int_t m, ae_int_t forecastlen, bool applysmoothing, RVector *trend);
} // end of namespace alglib_impl

namespace alglib {
DecClass(ssamodel, );

void ssacreate(ssamodel &s);
void ssasetwindow(const ssamodel &s, const ae_int_t windowwidth);
void ssasetseed(const ssamodel &s, const ae_int_t seed);
void ssasetpoweruplength(const ssamodel &s, const ae_int_t pwlen);
void ssasetmemorylimit(const ssamodel &s, const ae_int_t memlimit);
void ssaaddsequence(const ssamodel &s, const real_1d_array &x, const ae_int_t n);
void ssaaddsequence(const ssamodel &s, const real_1d_array &x);
void ssaappendpointandupdate(const ssamodel &s, const double x, const double updateits);
void ssaappendsequenceandupdate(const ssamodel &s, const real_1d_array &x, const ae_int_t nticks, const double updateits);
void ssaappendsequenceandupdate(const ssamodel &s, const real_1d_array &x, const double updateits);
void ssasetalgoprecomputed(const ssamodel &s, const real_2d_array &a, const ae_int_t windowwidth, const ae_int_t nbasis);
void ssasetalgoprecomputed(const ssamodel &s, const real_2d_array &a);
void ssasetalgotopkdirect(const ssamodel &s, const ae_int_t topk);
void ssasetalgotopkrealtime(const ssamodel &s, const ae_int_t topk);
void ssacleardata(const ssamodel &s);
void ssagetbasis(const ssamodel &s, real_2d_array &a, real_1d_array &sv, ae_int_t &windowwidth, ae_int_t &nbasis);
void ssagetlrr(const ssamodel &s, real_1d_array &a, ae_int_t &windowwidth);
void ssaanalyzelastwindow(const ssamodel &s, real_1d_array &trend, real_1d_array &noise, ae_int_t &nticks);
void ssaanalyzelast(const ssamodel &s, const ae_int_t nticks, real_1d_array &trend, real_1d_array &noise);
void ssaanalyzesequence(const ssamodel &s, const real_1d_array &data, const ae_int_t nticks, real_1d_array &trend, real_1d_array &noise);
void ssaanalyzesequence(const ssamodel &s, const real_1d_array &data, real_1d_array &trend, real_1d_array &noise);
void ssaforecastlast(const ssamodel &s, const ae_int_t nticks, real_1d_array &trend);
void ssaforecastsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t datalen, const ae_int_t forecastlen, const bool applysmoothing, real_1d_array &trend);
void ssaforecastsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t forecastlen, real_1d_array &trend);
void ssaforecastavglast(const ssamodel &s, const ae_int_t m, const ae_int_t nticks, real_1d_array &trend);
void ssaforecastavgsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t datalen, const ae_int_t m, const ae_int_t forecastlen, const bool applysmoothing, real_1d_array &trend);
void ssaforecastavgsequence(const ssamodel &s, const real_1d_array &data, const ae_int_t m, const ae_int_t forecastlen, real_1d_array &trend);
} // end of namespace alglib

// === LDA Package ===
// Depends on: (LinAlg) MATINV, EVD
namespace alglib_impl {
void fisherlda(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, RVector *w);
void fisherldan(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, RMatrix *w);
} // end of namespace alglib_impl

namespace alglib {
void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w);
void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);
} // end of namespace alglib

// === MCPD Package ===
// Depends on: (Optimization) MINBLEIC
namespace alglib_impl {
struct mcpdstate {
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
};
void mcpdstate_init(void *_p, bool make_automatic);
void mcpdstate_copy(void *_dst, void *_src, bool make_automatic);
void mcpdstate_free(void *_p, bool make_automatic);

struct mcpdreport {
   ae_int_t inneriterationscount;
   ae_int_t outeriterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
};
void mcpdreport_init(void *_p, bool make_automatic);
void mcpdreport_copy(void *_dst, void *_src, bool make_automatic);
void mcpdreport_free(void *_p, bool make_automatic);

void mcpdcreate(ae_int_t n, mcpdstate *s);
void mcpdcreateentry(ae_int_t n, ae_int_t entrystate, mcpdstate *s);
void mcpdcreateexit(ae_int_t n, ae_int_t exitstate, mcpdstate *s);
void mcpdcreateentryexit(ae_int_t n, ae_int_t entrystate, ae_int_t exitstate, mcpdstate *s);
void mcpdaddtrack(mcpdstate *s, RMatrix *xy, ae_int_t k);
void mcpdsetec(mcpdstate *s, RMatrix *ec);
void mcpdaddec(mcpdstate *s, ae_int_t i, ae_int_t j, double c);
void mcpdsetbc(mcpdstate *s, RMatrix *bndl, RMatrix *bndu);
void mcpdaddbc(mcpdstate *s, ae_int_t i, ae_int_t j, double bndl, double bndu);
void mcpdsetlc(mcpdstate *s, RMatrix *c, ZVector *ct, ae_int_t k);
void mcpdsettikhonovregularizer(mcpdstate *s, double v);
void mcpdsetprior(mcpdstate *s, RMatrix *pp);
void mcpdsetpredictionweights(mcpdstate *s, RVector *pw);
void mcpdsolve(mcpdstate *s);
void mcpdresults(mcpdstate *s, RMatrix *p, mcpdreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mcpdstate, );
DecClass(mcpdreport, ae_int_t &inneriterationscount; ae_int_t &outeriterationscount; ae_int_t &nfev; ae_int_t &terminationtype;);

void mcpdcreate(const ae_int_t n, mcpdstate &s);
void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s);
void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s);
void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy);
void mcpdsetec(const mcpdstate &s, const real_2d_array &ec);
void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c);
void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu);
void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct);
void mcpdsettikhonovregularizer(const mcpdstate &s, const double v);
void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp);
void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw);
void mcpdsolve(const mcpdstate &s);
void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep);
} // end of namespace alglib

// === LOGIT Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS
// Depends on: MLPBASE
namespace alglib_impl {
struct logitmodel {
   ae_vector w;
};
void logitmodel_init(void *_p, bool make_automatic);
void logitmodel_copy(void *_dst, void *_src, bool make_automatic);
void logitmodel_free(void *_p, bool make_automatic);

struct logitmcstate {
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
};
void logitmcstate_init(void *_p, bool make_automatic);
void logitmcstate_copy(void *_dst, void *_src, bool make_automatic);
void logitmcstate_free(void *_p, bool make_automatic);

struct mnlreport {
   ae_int_t ngrad;
   ae_int_t nhess;
};
void mnlreport_init(void *_p, bool make_automatic);
void mnlreport_copy(void *_dst, void *_src, bool make_automatic);
void mnlreport_free(void *_p, bool make_automatic);

void mnltrainh(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t *info, logitmodel *lm, mnlreport *rep);
void mnlprocess(logitmodel *lm, RVector *x, RVector *y);
void mnlprocessi(logitmodel *lm, RVector *x, RVector *y);
void mnlunpack(logitmodel *lm, RMatrix *a, ae_int_t *nvars, ae_int_t *nclasses);
void mnlpack(RMatrix *a, ae_int_t nvars, ae_int_t nclasses, logitmodel *lm);
void mnlcopy(logitmodel *lm1, logitmodel *lm2);
double mnlavgce(logitmodel *lm, RMatrix *xy, ae_int_t npoints);
double mnlrelclserror(logitmodel *lm, RMatrix *xy, ae_int_t npoints);
double mnlrmserror(logitmodel *lm, RMatrix *xy, ae_int_t npoints);
double mnlavgerror(logitmodel *lm, RMatrix *xy, ae_int_t npoints);
double mnlavgrelerror(logitmodel *lm, RMatrix *xy, ae_int_t ssize);
ae_int_t mnlclserror(logitmodel *lm, RMatrix *xy, ae_int_t npoints);
} // end of namespace alglib_impl

namespace alglib {
DecClass(logitmodel, );
DecClass(mnlreport, ae_int_t &ngrad; ae_int_t &nhess;);

void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep);
void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);
void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);
void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses);
void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm);
double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize);
ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
} // end of namespace alglib

// === KNN Package ===
// Depends on: (AlgLibMisc) HQRND, NEARESTNEIGHBOR
// Depends on: BDSS
namespace alglib_impl {
struct knnbuffer {
   kdtreerequestbuffer treebuf;
   ae_vector x;
   ae_vector y;
   ae_vector tags;
   ae_matrix xy;
};
void knnbuffer_init(void *_p, bool make_automatic);
void knnbuffer_copy(void *_dst, void *_src, bool make_automatic);
void knnbuffer_free(void *_p, bool make_automatic);

struct knnbuilder {
   ae_int_t dstype;
   ae_int_t npoints;
   ae_int_t nvars;
   bool iscls;
   ae_int_t nout;
   ae_matrix dsdata;
   ae_vector dsrval;
   ae_vector dsival;
   ae_int_t knnnrm;
};
void knnbuilder_init(void *_p, bool make_automatic);
void knnbuilder_copy(void *_dst, void *_src, bool make_automatic);
void knnbuilder_free(void *_p, bool make_automatic);

struct knnmodel {
   ae_int_t nvars;
   ae_int_t nout;
   ae_int_t k;
   double eps;
   bool iscls;
   bool isdummy;
   kdtree tree;
   knnbuffer buffer;
};
void knnmodel_init(void *_p, bool make_automatic);
void knnmodel_copy(void *_dst, void *_src, bool make_automatic);
void knnmodel_free(void *_p, bool make_automatic);
void knnalloc(ae_serializer *s, knnmodel *model);
void knnserialize(ae_serializer *s, knnmodel *model);
void knnunserialize(ae_serializer *s, knnmodel *model);

struct knnreport {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
};
void knnreport_init(void *_p, bool make_automatic);
void knnreport_copy(void *_dst, void *_src, bool make_automatic);
void knnreport_free(void *_p, bool make_automatic);

void knncreatebuffer(knnmodel *model, knnbuffer *buf);
void knnbuildercreate(knnbuilder *s);
void knnbuildersetdatasetreg(knnbuilder *s, RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nout);
void knnbuildersetdatasetcls(knnbuilder *s, RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses);
void knnbuildersetnorm(knnbuilder *s, ae_int_t nrmtype);
void knnbuilderbuildknnmodel(knnbuilder *s, ae_int_t k, double eps, knnmodel *model, knnreport *rep);
void knnrewritekeps(knnmodel *model, ae_int_t k, double eps);
void knnprocess(knnmodel *model, RVector *x, RVector *y);
double knnprocess0(knnmodel *model, RVector *x);
ae_int_t knnclassify(knnmodel *model, RVector *x);
void knnprocessi(knnmodel *model, RVector *x, RVector *y);
void knntsprocess(knnmodel *model, knnbuffer *buf, RVector *x, RVector *y);
double knnrelclserror(knnmodel *model, RMatrix *xy, ae_int_t npoints);
double knnavgce(knnmodel *model, RMatrix *xy, ae_int_t npoints);
double knnrmserror(knnmodel *model, RMatrix *xy, ae_int_t npoints);
double knnavgerror(knnmodel *model, RMatrix *xy, ae_int_t npoints);
double knnavgrelerror(knnmodel *model, RMatrix *xy, ae_int_t npoints);
void knnallerrors(knnmodel *model, RMatrix *xy, ae_int_t npoints, knnreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(knnbuffer, );
DecClass(knnbuilder, );
DecClass(knnmodel, );
DecClass(knnreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);
void knnserialize(knnmodel &obj, std::string &s_out);
void knnserialize(knnmodel &obj, std::ostream &s_out);
void knnunserialize(const std::string &s_in, knnmodel &obj);
void knnunserialize(const std::istream &s_in, knnmodel &obj);

void knncreatebuffer(const knnmodel &model, knnbuffer &buf);
void knnbuildercreate(knnbuilder &s);
void knnbuildersetdatasetreg(const knnbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nout);
void knnbuildersetdatasetcls(const knnbuilder &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses);
void knnbuildersetnorm(const knnbuilder &s, const ae_int_t nrmtype);
void knnbuilderbuildknnmodel(const knnbuilder &s, const ae_int_t k, const double eps, knnmodel &model, knnreport &rep);
void knnrewritekeps(const knnmodel &model, const ae_int_t k, const double eps);
void knnprocess(const knnmodel &model, const real_1d_array &x, real_1d_array &y);
double knnprocess0(const knnmodel &model, const real_1d_array &x);
ae_int_t knnclassify(const knnmodel &model, const real_1d_array &x);
void knnprocessi(const knnmodel &model, const real_1d_array &x, real_1d_array &y);
void knntsprocess(const knnmodel &model, const knnbuffer &buf, const real_1d_array &x, real_1d_array &y);
double knnrelclserror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);
double knnavgce(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);
double knnrmserror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);
double knnavgerror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);
double knnavgrelerror(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints);
void knnallerrors(const knnmodel &model, const real_2d_array &xy, const ae_int_t npoints, knnreport &rep);
} // end of namespace alglib

// === MLPTRAIN Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS
// Depends on: (Optimization) MINLBFGS
// Depends on: MLPE
namespace alglib_impl {
struct mlpreport {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
   ae_int_t ngrad;
   ae_int_t nhess;
   ae_int_t ncholesky;
};
void mlpreport_init(void *_p, bool make_automatic);
void mlpreport_copy(void *_dst, void *_src, bool make_automatic);
void mlpreport_free(void *_p, bool make_automatic);

struct mlpcvreport {
   double relclserror;
   double avgce;
   double rmserror;
   double avgerror;
   double avgrelerror;
};
void mlpcvreport_init(void *_p, bool make_automatic);
void mlpcvreport_copy(void *_dst, void *_src, bool make_automatic);
void mlpcvreport_free(void *_p, bool make_automatic);

struct smlptrnsession {
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
};
void smlptrnsession_init(void *_p, bool make_automatic);
void smlptrnsession_copy(void *_dst, void *_src, bool make_automatic);
void smlptrnsession_free(void *_p, bool make_automatic);

struct mlpetrnsession {
   ae_vector trnsubset;
   ae_vector valsubset;
   ae_shared_pool mlpsessions;
   mlpreport mlprep;
   multilayerperceptron network;
};
void mlpetrnsession_init(void *_p, bool make_automatic);
void mlpetrnsession_copy(void *_dst, void *_src, bool make_automatic);
void mlpetrnsession_free(void *_p, bool make_automatic);

struct mlptrainer {
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
};
void mlptrainer_init(void *_p, bool make_automatic);
void mlptrainer_copy(void *_dst, void *_src, bool make_automatic);
void mlptrainer_free(void *_p, bool make_automatic);

struct mlpparallelizationcv {
   multilayerperceptron network;
   mlpreport rep;
   ae_vector subset;
   ae_int_t subsetsize;
   ae_vector xyrow;
   ae_vector y;
   ae_int_t ngrad;
   ae_shared_pool trnpool;
};
void mlpparallelizationcv_init(void *_p, bool make_automatic);
void mlpparallelizationcv_copy(void *_dst, void *_src, bool make_automatic);
void mlpparallelizationcv_free(void *_p, bool make_automatic);

void mlptrainlm(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep);
void mlptrainlbfgs(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t *info, mlpreport *rep);
void mlptraines(multilayerperceptron *network, RMatrix *trnxy, ae_int_t trnsize, RMatrix *valxy, ae_int_t valsize, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep);
void mlpkfoldcvlbfgs(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t foldscount, ae_int_t *info, mlpreport *rep, mlpcvreport *cvrep);
void mlpkfoldcvlm(multilayerperceptron *network, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t foldscount, ae_int_t *info, mlpreport *rep, mlpcvreport *cvrep);
void mlpkfoldcv(mlptrainer *s, multilayerperceptron *network, ae_int_t nrestarts, ae_int_t foldscount, mlpreport *rep);
void mlpcreatetrainer(ae_int_t nin, ae_int_t nout, mlptrainer *s);
void mlpcreatetrainercls(ae_int_t nin, ae_int_t nclasses, mlptrainer *s);
void mlpsetdataset(mlptrainer *s, RMatrix *xy, ae_int_t npoints);
void mlpsetsparsedataset(mlptrainer *s, sparsematrix *xy, ae_int_t npoints);
void mlpsetdecay(mlptrainer *s, double decay);
void mlpsetcond(mlptrainer *s, double wstep, ae_int_t maxits);
void mlpsetalgobatch(mlptrainer *s);
void mlptrainnetwork(mlptrainer *s, multilayerperceptron *network, ae_int_t nrestarts, mlpreport *rep);
void mlpstarttraining(mlptrainer *s, multilayerperceptron *network, bool randomstart);
bool mlpcontinuetraining(mlptrainer *s, multilayerperceptron *network);
void mlpebagginglm(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep, mlpcvreport *ooberrors);
void mlpebagginglbfgs(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t *info, mlpreport *rep, mlpcvreport *ooberrors);
void mlpetraines(mlpensemble *ensemble, RMatrix *xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t *info, mlpreport *rep);
void mlptrainensemblees(mlptrainer *s, mlpensemble *ensemble, ae_int_t nrestarts, mlpreport *rep);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mlpreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror; ae_int_t &ngrad; ae_int_t &nhess; ae_int_t &ncholesky;);
DecClass(mlpcvreport, double &relclserror; double &avgce; double &rmserror; double &avgerror; double &avgrelerror;);
DecClass(mlptrainer, );

void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);
void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep);
void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);
void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);
void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);
void mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep);
void mlpcreatetrainer(const ae_int_t nin, const ae_int_t nout, mlptrainer &s);
void mlpcreatetrainercls(const ae_int_t nin, const ae_int_t nclasses, mlptrainer &s);
void mlpsetdataset(const mlptrainer &s, const real_2d_array &xy, const ae_int_t npoints);
void mlpsetsparsedataset(const mlptrainer &s, const sparsematrix &xy, const ae_int_t npoints);
void mlpsetdecay(const mlptrainer &s, const double decay);
void mlpsetcond(const mlptrainer &s, const double wstep, const ae_int_t maxits);
void mlpsetalgobatch(const mlptrainer &s);
void mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep);
void mlpstarttraining(const mlptrainer &s, const multilayerperceptron &network, const bool randomstart);
bool mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network);
void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);
void mlpebagginglbfgs(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);
void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);
void mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep);
} // end of namespace alglib

// === DATACOMP Package ===
// Depends on: CLUSTERING
namespace alglib_impl {
void kmeansgenerate(RMatrix *xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t restarts, ae_int_t *info, RMatrix *c, ZVector *xyc);
} // end of namespace alglib_impl

namespace alglib {
void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc);
} // end of namespace alglib

#endif // OnceOnly
