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
#ifndef OnceOnlyOptimization_h
#define OnceOnlyOptimization_h

#include "Solvers.h"

// === OPTGUARDAPI Package ===
// Depends on: (AlgLibInternal) APSERV
namespace alglib_impl {
typedef struct {
   bool nonc0suspected;
   bool nonc0test0positive;
   ae_int_t nonc0fidx;
   double nonc0lipschitzc;
   bool nonc1suspected;
   bool nonc1test0positive;
   bool nonc1test1positive;
   ae_int_t nonc1fidx;
   double nonc1lipschitzc;
   bool badgradsuspected;
   ae_int_t badgradfidx;
   ae_int_t badgradvidx;
   ae_vector badgradxbase;
   ae_matrix badgraduser;
   ae_matrix badgradnum;
} optguardreport;
typedef struct {
   bool positive;
   ae_int_t fidx;
   ae_vector x0;
   ae_vector d;
   ae_int_t n;
   ae_vector stp;
   ae_vector f;
   ae_int_t cnt;
   ae_int_t stpidxa;
   ae_int_t stpidxb;
} optguardnonc0report;
typedef struct {
   bool positive;
   ae_int_t fidx;
   ae_vector x0;
   ae_vector d;
   ae_int_t n;
   ae_vector stp;
   ae_vector f;
   ae_int_t cnt;
   ae_int_t stpidxa;
   ae_int_t stpidxb;
} optguardnonc1test0report;
typedef struct {
   bool positive;
   ae_int_t fidx;
   ae_int_t vidx;
   ae_vector x0;
   ae_vector d;
   ae_int_t n;
   ae_vector stp;
   ae_vector g;
   ae_int_t cnt;
   ae_int_t stpidxa;
   ae_int_t stpidxb;
} optguardnonc1test1report;

void optguardinitinternal(optguardreport *rep, ae_int_t n, ae_int_t k, ae_state *_state);
void optguardexportreport(optguardreport *srcrep, ae_int_t n, ae_int_t k, bool badgradhasxj, optguardreport *dstrep, ae_state *_state);
void smoothnessmonitorexportc1test0report(optguardnonc1test0report *srcrep, RVector *s, optguardnonc1test0report *dstrep, ae_state *_state);
void smoothnessmonitorexportc1test1report(optguardnonc1test1report *srcrep, RVector *s, optguardnonc1test1report *dstrep, ae_state *_state);
bool optguardallclear(optguardreport *rep, ae_state *_state);
void optguardreport_init(void *_p, ae_state *_state, bool make_automatic);
void optguardreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void optguardreport_free(void *_p, bool make_automatic);
void optguardnonc0report_init(void *_p, ae_state *_state, bool make_automatic);
void optguardnonc0report_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void optguardnonc0report_free(void *_p, bool make_automatic);
void optguardnonc1test0report_init(void *_p, ae_state *_state, bool make_automatic);
void optguardnonc1test0report_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void optguardnonc1test0report_free(void *_p, bool make_automatic);
void optguardnonc1test1report_init(void *_p, ae_state *_state, bool make_automatic);
void optguardnonc1test1report_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void optguardnonc1test1report_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(optguardreport, bool &nonc0suspected; bool &nonc0test0positive; ae_int_t &nonc0fidx; double &nonc0lipschitzc; bool &nonc1suspected; bool &nonc1test0positive; bool &nonc1test1positive; ae_int_t &nonc1fidx; double &nonc1lipschitzc; bool &badgradsuspected; ae_int_t &badgradfidx; ae_int_t &badgradvidx; real_1d_array badgradxbase; real_2d_array badgraduser; real_2d_array badgradnum;);
DecClass(optguardnonc0report, bool &positive; ae_int_t &fidx; real_1d_array x0; real_1d_array d; ae_int_t &n; real_1d_array stp; real_1d_array f; ae_int_t &cnt; ae_int_t &stpidxa; ae_int_t &stpidxb;);
DecClass(optguardnonc1test0report, bool &positive; ae_int_t &fidx; real_1d_array x0; real_1d_array d; ae_int_t &n; real_1d_array stp; real_1d_array f; ae_int_t &cnt; ae_int_t &stpidxa; ae_int_t &stpidxb;);
DecClass(optguardnonc1test1report, bool &positive; ae_int_t &fidx; ae_int_t &vidx; real_1d_array x0; real_1d_array d; ae_int_t &n; real_1d_array stp; real_1d_array g; ae_int_t &cnt; ae_int_t &stpidxa; ae_int_t &stpidxb;);
} // end of namespace alglib

// === OPTSERV Package ===
// Depends on: (LinAlg) SVD, MATINV
// Depends on: OPTGUARDAPI
namespace alglib_impl {
typedef struct {
   ae_vector norms;
   ae_vector alpha;
   ae_vector rho;
   ae_matrix yk;
   ae_vector idx;
   ae_vector bufa;
   ae_vector bufb;
} precbuflbfgs;
typedef struct {
   ae_int_t n;
   ae_int_t k;
   ae_vector d;
   ae_matrix v;
   ae_vector bufc;
   ae_matrix bufz;
   ae_matrix bufw;
   ae_vector tmp;
} precbuflowrank;
typedef struct {
   ae_int_t n;
   ae_int_t k;
   bool checksmoothness;
   ae_vector s;
   ae_vector dcur;
   ae_int_t enqueuedcnt;
   ae_vector enqueuedstp;
   ae_vector enqueuedx;
   ae_vector enqueuedfunc;
   ae_matrix enqueuedjac;
   ae_vector sortedstp;
   ae_vector sortedidx;
   ae_int_t sortedcnt;
   double probingstp;
   ae_vector probingf;
   ae_int_t probingnvalues;
   double probingstepmax;
   double probingstepscale;
   ae_int_t probingnstepsstored;
   ae_vector probingsteps;
   ae_matrix probingvalues;
   ae_matrix probingslopes;
   rcommstate probingrcomm;
   bool linesearchspoiled;
   bool linesearchstarted;
   double nonc0currentrating;
   double nonc1currentrating;
   bool badgradhasxj;
   optguardreport rep;
   double nonc0strrating;
   double nonc0lngrating;
   optguardnonc0report nonc0strrep;
   optguardnonc0report nonc0lngrep;
   double nonc1test0strrating;
   double nonc1test0lngrating;
   optguardnonc1test0report nonc1test0strrep;
   optguardnonc1test0report nonc1test0lngrep;
   double nonc1test1strrating;
   double nonc1test1lngrating;
   optguardnonc1test1report nonc1test1strrep;
   optguardnonc1test1report nonc1test1lngrep;
   bool needfij;
   ae_vector x;
   ae_vector fi;
   ae_matrix j;
   rcommstate rstateg0;
   ae_vector xbase;
   ae_vector fbase;
   ae_vector fm;
   ae_vector fc;
   ae_vector fp;
   ae_vector jm;
   ae_vector jc;
   ae_vector jp;
   ae_matrix jbaseusr;
   ae_matrix jbasenum;
   ae_vector stp;
   ae_vector bufr;
   ae_vector f;
   ae_vector g;
   ae_vector deltax;
   ae_vector tmpidx;
   ae_vector bufi;
   ae_vector xu;
   ae_vector du;
   ae_vector f0;
   ae_matrix j0;
} smoothnessmonitor;

void checkbcviolation(BVector *hasbndl, RVector *bndl, BVector *hasbndu, RVector *bndu, RVector *x, ae_int_t n, RVector *s, bool nonunits, double *bcerr, ae_int_t *bcidx, ae_state *_state);
void checklcviolation(RMatrix *cleic, ZVector *lcsrcidx, ae_int_t nec, ae_int_t nic, RVector *x, ae_int_t n, double *lcerr, ae_int_t *lcidx, ae_state *_state);
void checknlcviolation(RVector *fi, ae_int_t ng, ae_int_t nh, double *nlcerr, ae_int_t *nlcidx, ae_state *_state);
void unscaleandchecknlcviolation(RVector *fi, RVector *fscales, ae_int_t ng, ae_int_t nh, double *nlcerr, ae_int_t *nlcidx, ae_state *_state);
void trimprepare(double f, double *threshold, ae_state *_state);
void trimfunction(double *f, RVector *g, ae_int_t n, double threshold, ae_state *_state);
bool enforceboundaryconstraints(RVector *x, RVector *bl, BVector *havebl, RVector *bu, BVector *havebu, ae_int_t nmain, ae_int_t nslack, ae_state *_state);
void projectgradientintobc(RVector *x, RVector *g, RVector *bl, BVector *havebl, RVector *bu, BVector *havebu, ae_int_t nmain, ae_int_t nslack, ae_state *_state);
void calculatestepbound(RVector *x, RVector *d, double alpha, RVector *bndl, BVector *havebndl, RVector *bndu, BVector *havebndu, ae_int_t nmain, ae_int_t nslack, ae_int_t *variabletofreeze, double *valuetofreeze, double *maxsteplen, ae_state *_state);
ae_int_t postprocessboundedstep(RVector *x, RVector *xprev, RVector *bndl, BVector *havebndl, RVector *bndu, BVector *havebndu, ae_int_t nmain, ae_int_t nslack, ae_int_t variabletofreeze, double valuetofreeze, double steptaken, double maxsteplen, ae_state *_state);
void filterdirection(RVector *d, RVector *x, RVector *bndl, BVector *havebndl, RVector *bndu, BVector *havebndu, RVector *s, ae_int_t nmain, ae_int_t nslack, double droptol, ae_state *_state);
ae_int_t numberofchangedconstraints(RVector *x, RVector *xprev, RVector *bndl, BVector *havebndl, RVector *bndu, BVector *havebndu, ae_int_t nmain, ae_int_t nslack, ae_state *_state);
bool findfeasiblepoint(RVector *x, RVector *bndl, BVector *havebndl, RVector *bndu, BVector *havebndu, ae_int_t nmain, ae_int_t nslack, RMatrix *ce, ae_int_t k, double epsi, ae_int_t *qpits, ae_int_t *gpaits, ae_state *_state);
bool derivativecheck(double f0, double df0, double f1, double df1, double f, double df, double width, ae_state *_state);
void estimateparabolicmodel(double absasum, double absasum2, double mx, double mb, double md, double d1, double d2, ae_int_t *d1est, ae_int_t *d2est, ae_state *_state);
void inexactlbfgspreconditioner(RVector *s, ae_int_t n, RVector *d, RVector *c, RMatrix *w, ae_int_t k, precbuflbfgs *buf, ae_state *_state);
void preparelowrankpreconditioner(RVector *d, RVector *c, RMatrix *w, ae_int_t n, ae_int_t k, precbuflowrank *buf, ae_state *_state);
void applylowrankpreconditioner(RVector *s, precbuflowrank *buf, ae_state *_state);
void smoothnessmonitorinit(smoothnessmonitor *monitor, RVector *s, ae_int_t n, ae_int_t k, bool checksmoothness, ae_state *_state);
void smoothnessmonitorstartlinesearch(smoothnessmonitor *monitor, RVector *x, RVector *fi, RMatrix *jac, ae_state *_state);
void smoothnessmonitorstartlinesearch1u(smoothnessmonitor *monitor, RVector *s, RVector *invs, RVector *x, double f0, RVector *j0, ae_state *_state);
void smoothnessmonitorenqueuepoint(smoothnessmonitor *monitor, RVector *d, double stp, RVector *x, RVector *fi, RMatrix *jac, ae_state *_state);
void smoothnessmonitorenqueuepoint1u(smoothnessmonitor *monitor, RVector *s, RVector *invs, RVector *d, double stp, RVector *x, double f0, RVector *j0, ae_state *_state);
void smoothnessmonitorfinalizelinesearch(smoothnessmonitor *monitor, ae_state *_state);
void smoothnessmonitorstartprobing(smoothnessmonitor *monitor, double stpmax, ae_int_t nvalues, double stepscale, ae_state *_state);
bool smoothnessmonitorprobe(smoothnessmonitor *monitor, ae_state *_state);
void smoothnessmonitortraceprobingresults(smoothnessmonitor *monitor, ae_state *_state);
void smoothnessmonitortracestatus(smoothnessmonitor *monitor, bool callersuggeststrace, ae_state *_state);
void smoothnessmonitorexportreport(smoothnessmonitor *monitor, optguardreport *rep, ae_state *_state);
bool smoothnessmonitorcheckgradientatx0(smoothnessmonitor *monitor, RVector *unscaledx0, RVector *s, RVector *bndl, RVector *bndu, bool hasboxconstraints, double teststep, ae_state *_state);
void precbuflbfgs_init(void *_p, ae_state *_state, bool make_automatic);
void precbuflbfgs_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void precbuflbfgs_free(void *_p, bool make_automatic);
void precbuflowrank_init(void *_p, ae_state *_state, bool make_automatic);
void precbuflowrank_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void precbuflowrank_free(void *_p, bool make_automatic);
void smoothnessmonitor_init(void *_p, ae_state *_state, bool make_automatic);
void smoothnessmonitor_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void smoothnessmonitor_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === MINLBFGS Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: (LinAlg) FBLS
// Depends on: OPTSERV
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t m;
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxits;
   bool xrep;
   double stpmax;
   ae_vector s;
   double diffstep;
   ae_int_t nfev;
   ae_int_t mcstage;
   ae_int_t k;
   ae_int_t q;
   ae_int_t p;
   ae_vector rho;
   ae_matrix yk;
   ae_matrix sk;
   ae_vector xp;
   ae_vector theta;
   ae_vector d;
   double stp;
   ae_vector work;
   double fold;
   double trimthreshold;
   ae_vector xbase;
   ae_int_t prectype;
   double gammak;
   ae_matrix denseh;
   ae_vector diagh;
   ae_vector precc;
   ae_vector precd;
   ae_matrix precw;
   ae_int_t preck;
   precbuflbfgs precbuf;
   precbuflowrank lowrankbuf;
   double fbase;
   double fm2;
   double fm1;
   double fp1;
   double fp2;
   ae_vector autobuf;
   ae_vector invs;
   ae_vector x;
   double f;
   ae_vector g;
   bool needf;
   bool needfg;
   bool xupdated;
   bool userterminationneeded;
   double teststep;
   rcommstate rstate;
   ae_int_t repiterationscount;
   ae_int_t repnfev;
   ae_int_t repterminationtype;
   linminstate lstate;
   ae_int_t smoothnessguardlevel;
   smoothnessmonitor smonitor;
   ae_vector lastscaleused;
} minlbfgsstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
} minlbfgsreport;

void minlbfgscreate(ae_int_t n, ae_int_t m, RVector *x, minlbfgsstate *state, ae_state *_state);
void minlbfgscreatef(ae_int_t n, ae_int_t m, RVector *x, double diffstep, minlbfgsstate *state, ae_state *_state);
void minlbfgssetcond(minlbfgsstate *state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
void minlbfgssetxrep(minlbfgsstate *state, bool needxrep, ae_state *_state);
void minlbfgssetstpmax(minlbfgsstate *state, double stpmax, ae_state *_state);
void minlbfgssetscale(minlbfgsstate *state, RVector *s, ae_state *_state);
void minlbfgscreatex(ae_int_t n, ae_int_t m, RVector *x, ae_int_t flags, double diffstep, minlbfgsstate *state, ae_state *_state);
void minlbfgssetprecdefault(minlbfgsstate *state, ae_state *_state);
void minlbfgssetpreccholesky(minlbfgsstate *state, RMatrix *p, bool isupper, ae_state *_state);
void minlbfgssetprecdiag(minlbfgsstate *state, RVector *d, ae_state *_state);
void minlbfgssetprecscale(minlbfgsstate *state, ae_state *_state);
void minlbfgssetprecrankklbfgsfast(minlbfgsstate *state, RVector *d, RVector *c, RMatrix *w, ae_int_t cnt, ae_state *_state);
void minlbfgssetpreclowrankexact(minlbfgsstate *state, RVector *d, RVector *c, RMatrix *w, ae_int_t cnt, ae_state *_state);
bool minlbfgsiteration(minlbfgsstate *state, ae_state *_state);
void minlbfgsoptguardgradient(minlbfgsstate *state, double teststep, ae_state *_state);
void minlbfgsoptguardsmoothness(minlbfgsstate *state, ae_int_t level, ae_state *_state);
void minlbfgsoptguardresults(minlbfgsstate *state, optguardreport *rep, ae_state *_state);
void minlbfgsoptguardnonc1test0results(minlbfgsstate *state, optguardnonc1test0report *strrep, optguardnonc1test0report *lngrep, ae_state *_state);
void minlbfgsoptguardnonc1test1results(minlbfgsstate *state, optguardnonc1test1report *strrep, optguardnonc1test1report *lngrep, ae_state *_state);
void minlbfgsresults(minlbfgsstate *state, RVector *x, minlbfgsreport *rep, ae_state *_state);
void minlbfgsresultsbuf(minlbfgsstate *state, RVector *x, minlbfgsreport *rep, ae_state *_state);
void minlbfgsrestartfrom(minlbfgsstate *state, RVector *x, ae_state *_state);
void minlbfgsrequesttermination(minlbfgsstate *state, ae_state *_state);
void minlbfgsstate_init(void *_p, ae_state *_state, bool make_automatic);
void minlbfgsstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlbfgsstate_free(void *_p, bool make_automatic);
void minlbfgsreport_init(void *_p, ae_state *_state, bool make_automatic);
void minlbfgsreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlbfgsreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minlbfgsstate, bool &needf; bool &needfg; bool &xupdated; double &f; real_1d_array g; real_1d_array x;);
DecClass(minlbfgsreport, ae_int_t &iterationscount; ae_int_t &nfev; ae_int_t &terminationtype;);

void minlbfgscreate(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgscreate(const ae_int_t m, const real_1d_array &x, minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgscreatef(const ae_int_t n, const ae_int_t m, const real_1d_array &x, const double diffstep, minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgscreatef(const ae_int_t m, const real_1d_array &x, const double diffstep, minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgssetcond(const minlbfgsstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minlbfgssetxrep(const minlbfgsstate &state, const bool needxrep, const xparams _xparams = xdefault);
void minlbfgssetstpmax(const minlbfgsstate &state, const double stpmax, const xparams _xparams = xdefault);
void minlbfgssetscale(const minlbfgsstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minlbfgssetprecdefault(const minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgssetpreccholesky(const minlbfgsstate &state, const real_2d_array &p, const bool isupper, const xparams _xparams = xdefault);
void minlbfgssetprecdiag(const minlbfgsstate &state, const real_1d_array &d, const xparams _xparams = xdefault);
void minlbfgssetprecscale(const minlbfgsstate &state, const xparams _xparams = xdefault);
bool minlbfgsiteration(const minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgsoptimize(minlbfgsstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlbfgsoptimize(minlbfgsstate &state, void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlbfgsoptguardgradient(const minlbfgsstate &state, const double teststep, const xparams _xparams = xdefault);
void minlbfgsoptguardsmoothness(const minlbfgsstate &state, const ae_int_t level, const xparams _xparams = xdefault);
void minlbfgsoptguardsmoothness(const minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgsoptguardresults(const minlbfgsstate &state, optguardreport &rep, const xparams _xparams = xdefault);
void minlbfgsoptguardnonc1test0results(const minlbfgsstate &state, optguardnonc1test0report &strrep, optguardnonc1test0report &lngrep, const xparams _xparams = xdefault);
void minlbfgsoptguardnonc1test1results(const minlbfgsstate &state, optguardnonc1test1report &strrep, optguardnonc1test1report &lngrep, const xparams _xparams = xdefault);
void minlbfgsresults(const minlbfgsstate &state, real_1d_array &x, minlbfgsreport &rep, const xparams _xparams = xdefault);
void minlbfgsresultsbuf(const minlbfgsstate &state, real_1d_array &x, minlbfgsreport &rep, const xparams _xparams = xdefault);
void minlbfgsrestartfrom(const minlbfgsstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void minlbfgsrequesttermination(const minlbfgsstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === CQMODELS Package ===
// Depends on: (LinAlg) TRFAC, FBLS
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t k;
   double alpha;
   double tau;
   double theta;
   ae_matrix a;
   ae_matrix q;
   ae_vector b;
   ae_vector r;
   ae_vector xc;
   ae_vector d;
   ae_vector activeset;
   ae_matrix tq2dense;
   ae_matrix tk2;
   ae_vector tq2diag;
   ae_vector tq1;
   ae_vector tk1;
   double tq0;
   double tk0;
   ae_vector txc;
   ae_vector tb;
   ae_int_t nfree;
   ae_int_t ecakind;
   ae_matrix ecadense;
   ae_matrix eq;
   ae_matrix eccm;
   ae_vector ecadiag;
   ae_vector eb;
   double ec;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmpg;
   ae_matrix tmp2;
   bool ismaintermchanged;
   bool issecondarytermchanged;
   bool islineartermchanged;
   bool isactivesetchanged;
} convexquadraticmodel;

void cqminit(ae_int_t n, convexquadraticmodel *s, ae_state *_state);
void cqmseta(convexquadraticmodel *s, RMatrix *a, bool isupper, double alpha, ae_state *_state);
void cqmgeta(convexquadraticmodel *s, RMatrix *a, ae_state *_state);
void cqmrewritedensediagonal(convexquadraticmodel *s, RVector *z, ae_state *_state);
void cqmsetd(convexquadraticmodel *s, RVector *d, double tau, ae_state *_state);
void cqmdropa(convexquadraticmodel *s, ae_state *_state);
void cqmsetb(convexquadraticmodel *s, RVector *b, ae_state *_state);
void cqmsetq(convexquadraticmodel *s, RMatrix *q, RVector *r, ae_int_t k, double theta, ae_state *_state);
void cqmsetactiveset(convexquadraticmodel *s, RVector *x, BVector *activeset, ae_state *_state);
double cqmeval(convexquadraticmodel *s, RVector *x, ae_state *_state);
void cqmevalx(convexquadraticmodel *s, RVector *x, double *r, double *noise, ae_state *_state);
void cqmgradunconstrained(convexquadraticmodel *s, RVector *x, RVector *g, ae_state *_state);
double cqmxtadx2(convexquadraticmodel *s, RVector *x, RVector *tmp, ae_state *_state);
void cqmadx(convexquadraticmodel *s, RVector *x, RVector *y, ae_state *_state);
bool cqmconstrainedoptimum(convexquadraticmodel *s, RVector *x, ae_state *_state);
void cqmscalevector(convexquadraticmodel *s, RVector *x, ae_state *_state);
void cqmgetdiaga(convexquadraticmodel *s, RVector *x, ae_state *_state);
double cqmdebugconstrainedevalt(convexquadraticmodel *s, RVector *x, ae_state *_state);
double cqmdebugconstrainedevale(convexquadraticmodel *s, RVector *x, ae_state *_state);
void convexquadraticmodel_init(void *_p, ae_state *_state, bool make_automatic);
void convexquadraticmodel_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void convexquadraticmodel_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === LPQPSERV Package ===
// Depends on: (LinAlg) SPARSE
namespace alglib_impl {
void scaleshiftbcinplace(RVector *s, RVector *xorigin, RVector *bndl, RVector *bndu, ae_int_t n, ae_state *_state);
void scaleshiftdensebrlcinplace(RVector *s, RVector *xorigin, ae_int_t n, RMatrix *densea, RVector *ab, RVector *ar, ae_int_t m, ae_state *_state);
void scaleshiftmixedbrlcinplace(RVector *s, RVector *xorigin, ae_int_t n, sparsematrix *sparsea, ae_int_t msparse, RMatrix *densea, ae_int_t mdense, RVector *ab, RVector *ar, ae_state *_state);
void scaledenseqpinplace(RMatrix *densea, bool isupper, ae_int_t nmain, RVector *denseb, ae_int_t ntotal, RVector *s, ae_state *_state);
void scalesparseqpinplace(RVector *s, ae_int_t n, sparsematrix *sparsea, RVector *denseb, ae_state *_state);
void normalizedensebrlcinplace(RMatrix *densea, RVector *ab, RVector *ar, ae_int_t n, ae_int_t m, RVector *rownorms, bool neednorms, ae_state *_state);
void normalizemixedbrlcinplace(sparsematrix *sparsea, ae_int_t msparse, RMatrix *densea, ae_int_t mdense, RVector *ab, RVector *ar, ae_int_t n, bool limitedamplification, RVector *rownorms, bool neednorms, ae_state *_state);
double normalizedenseqpinplace(RMatrix *densea, bool isupper, ae_int_t nmain, RVector *denseb, ae_int_t ntotal, ae_state *_state);
double normalizesparseqpinplace(sparsematrix *sparsea, bool isupper, RVector *denseb, ae_int_t n, ae_state *_state);
void unscaleunshiftpointbc(RVector *s, RVector *xorigin, RVector *rawbndl, RVector *rawbndu, RVector *sclsftbndl, RVector *sclsftbndu, BVector *hasbndl, BVector *hasbndu, RVector *x, ae_int_t n, ae_state *_state);
} // end of namespace alglib_impl

// === SNNLS Package ===
// Depends on: (LinAlg) TRFAC, FBLS
namespace alglib_impl {
typedef struct {
   ae_int_t ns;
   ae_int_t nd;
   ae_int_t nr;
   ae_matrix densea;
   ae_vector b;
   ae_vector nnc;
   double debugflops;
   ae_int_t debugmaxinnerits;
   ae_vector xn;
   ae_vector xp;
   ae_matrix tmpca;
   ae_matrix tmplq;
   ae_matrix trda;
   ae_vector trdd;
   ae_vector crb;
   ae_vector g;
   ae_vector d;
   ae_vector dx;
   ae_vector diagaa;
   ae_vector cb;
   ae_vector cx;
   ae_vector cborg;
   ae_vector tmpcholesky;
   ae_vector r;
   ae_vector regdiag;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmp2;
   ae_vector rdtmprowmap;
} snnlssolver;

void snnlsinit(ae_int_t nsmax, ae_int_t ndmax, ae_int_t nrmax, snnlssolver *s, ae_state *_state);
void snnlssetproblem(snnlssolver *s, RMatrix *a, RVector *b, ae_int_t ns, ae_int_t nd, ae_int_t nr, ae_state *_state);
void snnlsdropnnc(snnlssolver *s, ae_int_t idx, ae_state *_state);
void snnlssolve(snnlssolver *s, RVector *x, ae_state *_state);
void snnlssolver_init(void *_p, ae_state *_state, bool make_automatic);
void snnlssolver_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void snnlssolver_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === SACTIVESETS Package ===
// Depends on: OPTSERV, SNNLS
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t algostate;
   ae_vector xc;
   bool hasxc;
   ae_vector s;
   ae_vector h;
   ae_vector cstatus;
   bool basisisready;
   ae_matrix sdensebatch;
   ae_matrix pdensebatch;
   ae_matrix idensebatch;
   ae_int_t densebatchsize;
   ae_vector sparsebatch;
   ae_int_t sparsebatchsize;
   ae_int_t basisage;
   bool feasinitpt;
   bool constraintschanged;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_vector bndl;
   ae_vector bndu;
   ae_matrix cleic;
   ae_int_t nec;
   ae_int_t nic;
   ae_vector mtnew;
   ae_vector mtx;
   ae_vector mtas;
   ae_vector cdtmp;
   ae_vector corrtmp;
   ae_vector unitdiagonal;
   snnlssolver solver;
   ae_vector scntmp;
   ae_vector tmp0;
   ae_vector tmpfeas;
   ae_matrix tmpm0;
   ae_vector rctmps;
   ae_vector rctmpg;
   ae_vector rctmprightpart;
   ae_matrix rctmpdense0;
   ae_matrix rctmpdense1;
   ae_vector rctmpisequality;
   ae_vector rctmpconstraintidx;
   ae_vector rctmplambdas;
   ae_matrix tmpbasis;
   ae_vector tmpnormestimates;
   ae_vector tmpreciph;
   ae_vector tmpprodp;
   ae_vector tmpprods;
   ae_vector tmpcp;
   ae_vector tmpcs;
   ae_vector tmpci;
} sactiveset;

void sasinit(ae_int_t n, sactiveset *s, ae_state *_state);
void sassetscale(sactiveset *state, RVector *s, ae_state *_state);
void sassetprecdiag(sactiveset *state, RVector *d, ae_state *_state);
void sassetbc(sactiveset *state, RVector *bndl, RVector *bndu, ae_state *_state);
void sassetlc(sactiveset *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void sassetlcx(sactiveset *state, RMatrix *cleic, ae_int_t nec, ae_int_t nic, ae_state *_state);
bool sasstartoptimization(sactiveset *state, RVector *x, ae_state *_state);
void sasexploredirection(sactiveset *state, RVector *d, double *stpmax, ae_int_t *cidx, double *vval, ae_state *_state);
ae_int_t sasmoveto(sactiveset *state, RVector *xn, bool needact, ae_int_t cidx, double cval, ae_state *_state);
void sasimmediateactivation(sactiveset *state, ae_int_t cidx, double cval, ae_state *_state);
void sasconstraineddescent(sactiveset *state, RVector *g, RVector *d, ae_state *_state);
void sasconstraineddescentprec(sactiveset *state, RVector *g, RVector *d, ae_state *_state);
void sasconstraineddirection(sactiveset *state, RVector *d, ae_state *_state);
void sasconstraineddirectionprec(sactiveset *state, RVector *d, ae_state *_state);
void sascorrection(sactiveset *state, RVector *x, double *penalty, ae_state *_state);
double sasactivelcpenalty1(sactiveset *state, RVector *x, ae_state *_state);
double sasscaledconstrainednorm(sactiveset *state, RVector *d, ae_state *_state);
void sasstopoptimization(sactiveset *state, ae_state *_state);
void sasreactivateconstraints(sactiveset *state, RVector *gc, ae_state *_state);
void sasreactivateconstraintsprec(sactiveset *state, RVector *gc, ae_state *_state);
void sasrebuildbasis(sactiveset *state, ae_state *_state);
void sasappendtobasis(sactiveset *state, BVector *newentries, ae_state *_state);
void sactiveset_init(void *_p, ae_state *_state, bool make_automatic);
void sactiveset_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void sactiveset_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === QQPSOLVER Package ===
// Depends on: CQMODELS, SACTIVESETS
namespace alglib_impl {
typedef struct {
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxouterits;
   bool cgphase;
   bool cnphase;
   ae_int_t cgminits;
   ae_int_t cgmaxits;
   ae_int_t cnmaxupdates;
   ae_int_t sparsesolver;
} qqpsettings;
typedef struct {
   ae_int_t n;
   ae_int_t akind;
   ae_matrix densea;
   sparsematrix sparsea;
   bool sparseupper;
   double absamax;
   double absasum;
   double absasum2;
   ae_vector b;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector havebndl;
   ae_vector havebndu;
   ae_vector xs;
   ae_vector xf;
   ae_vector gc;
   ae_vector xp;
   ae_vector dc;
   ae_vector dp;
   ae_vector cgc;
   ae_vector cgp;
   sactiveset sas;
   ae_vector activated;
   ae_int_t nfree;
   ae_int_t cnmodelage;
   ae_matrix densez;
   sparsematrix sparsecca;
   ae_vector yidx;
   ae_vector regdiag;
   ae_vector regx0;
   ae_vector tmpcn;
   ae_vector tmpcni;
   ae_vector tmpcnb;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector stpbuf;
   sparsebuffers sbuf;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repncholesky;
   ae_int_t repncupdates;
} qqpbuffers;

void qqploaddefaults(ae_int_t n, qqpsettings *s, ae_state *_state);
void qqpcopysettings(qqpsettings *src, qqpsettings *dst, ae_state *_state);
void qqppreallocategrowdense(qqpbuffers *sstate, ae_int_t nexpected, ae_int_t ngrowto, ae_state *_state);
void qqpoptimize(convexquadraticmodel *cqmac, sparsematrix *sparseac, RMatrix *denseac, ae_int_t akind, bool isupper, RVector *bc, RVector *bndlc, RVector *bnduc, RVector *sc, RVector *xoriginc, ae_int_t nc, qqpsettings *settings, qqpbuffers *sstate, RVector *xs, ae_int_t *terminationtype, ae_state *_state);
void qqpsettings_init(void *_p, ae_state *_state, bool make_automatic);
void qqpsettings_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void qqpsettings_free(void *_p, bool make_automatic);
void qqpbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void qqpbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void qqpbuffers_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === QPDENSEAULSOLVER Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS, LINLSQR
// Depends on: MINLBFGS, LPQPSERV, QQPSOLVER
namespace alglib_impl {
typedef struct {
   double epsx;
   ae_int_t outerits;
   double rho;
} qpdenseaulsettings;
typedef struct {
   ae_vector nulc;
   ae_matrix sclsfta;
   ae_vector sclsftb;
   ae_vector sclsfthasbndl;
   ae_vector sclsfthasbndu;
   ae_vector sclsftbndl;
   ae_vector sclsftbndu;
   ae_vector sclsftxc;
   ae_matrix sclsftcleic;
   ae_vector cidx;
   ae_vector cscales;
   ae_matrix exa;
   ae_vector exb;
   ae_vector exxc;
   ae_vector exbndl;
   ae_vector exbndu;
   ae_vector exscale;
   ae_vector exxorigin;
   qqpsettings qqpsettingsuser;
   qqpbuffers qqpbuf;
   ae_vector nulcest;
   ae_vector tmpg;
   ae_vector tmp0;
   ae_matrix tmp2;
   ae_vector modelg;
   ae_vector d;
   ae_vector deltax;
   convexquadraticmodel dummycqm;
   sparsematrix dummysparse;
   ae_matrix qrkkt;
   ae_vector qrrightpart;
   ae_vector qrtau;
   ae_vector qrsv0;
   ae_vector qrsvx1;
   ae_vector nicerr;
   ae_vector nicnact;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repncholesky;
   ae_int_t repnwrkchanges;
   ae_int_t repnwrk0;
   ae_int_t repnwrk1;
   ae_int_t repnwrkf;
   ae_int_t repnmv;
} qpdenseaulbuffers;

void qpdenseaulloaddefaults(ae_int_t nmain, qpdenseaulsettings *s, ae_state *_state);
void qpdenseauloptimize(convexquadraticmodel *a, sparsematrix *sparsea, ae_int_t akind, bool sparseaupper, RVector *b, RVector *bndl, RVector *bndu, RVector *s, RVector *xorigin, ae_int_t nn, RMatrix *cleic, ae_int_t dnec, ae_int_t dnic, sparsematrix *scleic, ae_int_t snec, ae_int_t snic, bool renormlc, qpdenseaulsettings *settings, qpdenseaulbuffers *state, RVector *xs, RVector *lagbc, RVector *laglc, ae_int_t *terminationtype, ae_state *_state);
void qpdenseaulsettings_init(void *_p, ae_state *_state, bool make_automatic);
void qpdenseaulsettings_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void qpdenseaulsettings_free(void *_p, bool make_automatic);
void qpdenseaulbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void qpdenseaulbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void qpdenseaulbuffers_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === MINBLEIC Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: CQMODELS, SACTIVESETS
namespace alglib_impl {
typedef struct {
   ae_int_t nmain;
   ae_int_t nslack;
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxits;
   bool xrep;
   bool drep;
   double stpmax;
   double diffstep;
   sactiveset sas;
   ae_vector s;
   ae_int_t prectype;
   ae_vector diagh;
   ae_vector x;
   double f;
   ae_vector g;
   bool needf;
   bool needfg;
   bool xupdated;
   bool lsstart;
   bool steepestdescentstep;
   bool boundedstep;
   bool userterminationneeded;
   rcommstate rstate;
   ae_vector ugc;
   ae_vector cgc;
   ae_vector xn;
   ae_vector ugn;
   ae_vector cgn;
   ae_vector xp;
   double fc;
   double fn;
   double fp;
   ae_vector d;
   ae_matrix cleic;
   ae_int_t nec;
   ae_int_t nic;
   double lastgoodstep;
   double lastscaledgoodstep;
   double maxscaledgrad;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_vector bndl;
   ae_vector bndu;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repnfev;
   ae_int_t repvaridx;
   ae_int_t repterminationtype;
   double repdebugeqerr;
   double repdebugfs;
   double repdebugff;
   double repdebugdx;
   ae_int_t repdebugfeasqpits;
   ae_int_t repdebugfeasgpaits;
   ae_vector xstart;
   snnlssolver solver;
   double fbase;
   double fm2;
   double fm1;
   double fp1;
   double fp2;
   double xm1;
   double xp1;
   double gm1;
   double gp1;
   ae_int_t cidx;
   double cval;
   ae_vector tmpprec;
   ae_vector tmp0;
   ae_int_t nfev;
   ae_int_t mcstage;
   double stp;
   double curstpmax;
   double activationstep;
   ae_vector work;
   linminstate lstate;
   double trimthreshold;
   ae_int_t nonmonotoniccnt;
   ae_matrix bufyk;
   ae_matrix bufsk;
   ae_vector bufrho;
   ae_vector buftheta;
   ae_int_t bufsize;
   double teststep;
   ae_int_t smoothnessguardlevel;
   smoothnessmonitor smonitor;
   ae_vector lastscaleused;
   ae_vector invs;
} minbleicstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   ae_int_t varidx;
   ae_int_t terminationtype;
   double debugeqerr;
   double debugfs;
   double debugff;
   double debugdx;
   ae_int_t debugfeasqpits;
   ae_int_t debugfeasgpaits;
   ae_int_t inneriterationscount;
   ae_int_t outeriterationscount;
} minbleicreport;

void minbleiccreate(ae_int_t n, RVector *x, minbleicstate *state, ae_state *_state);
void minbleiccreatef(ae_int_t n, RVector *x, double diffstep, minbleicstate *state, ae_state *_state);
void minbleicsetbc(minbleicstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minbleicsetlc(minbleicstate *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void minbleicsetcond(minbleicstate *state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
void minbleicsetscale(minbleicstate *state, RVector *s, ae_state *_state);
void minbleicsetprecdefault(minbleicstate *state, ae_state *_state);
void minbleicsetprecdiag(minbleicstate *state, RVector *d, ae_state *_state);
void minbleicsetprecscale(minbleicstate *state, ae_state *_state);
void minbleicsetxrep(minbleicstate *state, bool needxrep, ae_state *_state);
void minbleicsetdrep(minbleicstate *state, bool needdrep, ae_state *_state);
void minbleicsetstpmax(minbleicstate *state, double stpmax, ae_state *_state);
bool minbleiciteration(minbleicstate *state, ae_state *_state);
void minbleicoptguardgradient(minbleicstate *state, double teststep, ae_state *_state);
void minbleicoptguardsmoothness(minbleicstate *state, ae_int_t level, ae_state *_state);
void minbleicoptguardresults(minbleicstate *state, optguardreport *rep, ae_state *_state);
void minbleicoptguardnonc1test0results(minbleicstate *state, optguardnonc1test0report *strrep, optguardnonc1test0report *lngrep, ae_state *_state);
void minbleicoptguardnonc1test1results(minbleicstate *state, optguardnonc1test1report *strrep, optguardnonc1test1report *lngrep, ae_state *_state);
void minbleicresults(minbleicstate *state, RVector *x, minbleicreport *rep, ae_state *_state);
void minbleicresultsbuf(minbleicstate *state, RVector *x, minbleicreport *rep, ae_state *_state);
void minbleicrestartfrom(minbleicstate *state, RVector *x, ae_state *_state);
void minbleicrequesttermination(minbleicstate *state, ae_state *_state);
void minbleicemergencytermination(minbleicstate *state, ae_state *_state);
void minbleicstate_init(void *_p, ae_state *_state, bool make_automatic);
void minbleicstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minbleicstate_free(void *_p, bool make_automatic);
void minbleicreport_init(void *_p, ae_state *_state, bool make_automatic);
void minbleicreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minbleicreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minbleicstate, bool &needf; bool &needfg; bool &xupdated; double &f; real_1d_array g; real_1d_array x;);
DecClass(minbleicreport, ae_int_t &iterationscount; ae_int_t &nfev; ae_int_t &varidx; ae_int_t &terminationtype; double &debugeqerr; double &debugfs; double &debugff; double &debugdx; ae_int_t &debugfeasqpits; ae_int_t &debugfeasgpaits; ae_int_t &inneriterationscount; ae_int_t &outeriterationscount;);

void minbleiccreate(const ae_int_t n, const real_1d_array &x, minbleicstate &state, const xparams _xparams = xdefault);
void minbleiccreate(const real_1d_array &x, minbleicstate &state, const xparams _xparams = xdefault);
void minbleiccreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minbleicstate &state, const xparams _xparams = xdefault);
void minbleiccreatef(const real_1d_array &x, const double diffstep, minbleicstate &state, const xparams _xparams = xdefault);
void minbleicsetbc(const minbleicstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minbleicsetlc(const minbleicstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minbleicsetlc(const minbleicstate &state, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = xdefault);
void minbleicsetcond(const minbleicstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minbleicsetscale(const minbleicstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minbleicsetprecdefault(const minbleicstate &state, const xparams _xparams = xdefault);
void minbleicsetprecdiag(const minbleicstate &state, const real_1d_array &d, const xparams _xparams = xdefault);
void minbleicsetprecscale(const minbleicstate &state, const xparams _xparams = xdefault);
void minbleicsetxrep(const minbleicstate &state, const bool needxrep, const xparams _xparams = xdefault);
void minbleicsetstpmax(const minbleicstate &state, const double stpmax, const xparams _xparams = xdefault);
bool minbleiciteration(const minbleicstate &state, const xparams _xparams = xdefault);
void minbleicoptimize(minbleicstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minbleicoptimize(minbleicstate &state, void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minbleicoptguardgradient(const minbleicstate &state, const double teststep, const xparams _xparams = xdefault);
void minbleicoptguardsmoothness(const minbleicstate &state, const ae_int_t level, const xparams _xparams = xdefault);
void minbleicoptguardsmoothness(const minbleicstate &state, const xparams _xparams = xdefault);
void minbleicoptguardresults(const minbleicstate &state, optguardreport &rep, const xparams _xparams = xdefault);
void minbleicoptguardnonc1test0results(const minbleicstate &state, optguardnonc1test0report &strrep, optguardnonc1test0report &lngrep, const xparams _xparams = xdefault);
void minbleicoptguardnonc1test1results(const minbleicstate &state, optguardnonc1test1report &strrep, optguardnonc1test1report &lngrep, const xparams _xparams = xdefault);
void minbleicresults(const minbleicstate &state, real_1d_array &x, minbleicreport &rep, const xparams _xparams = xdefault);
void minbleicresultsbuf(const minbleicstate &state, real_1d_array &x, minbleicreport &rep, const xparams _xparams = xdefault);
void minbleicrestartfrom(const minbleicstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void minbleicrequesttermination(const minbleicstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === QPBLEICSOLVER Package ===
// Depends on: MINBLEIC
namespace alglib_impl {
typedef struct {
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxits;
} qpbleicsettings;
typedef struct {
   minbleicstate solver;
   minbleicreport solverrep;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmpi;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
} qpbleicbuffers;

void qpbleicloaddefaults(ae_int_t nmain, qpbleicsettings *s, ae_state *_state);
void qpbleiccopysettings(qpbleicsettings *src, qpbleicsettings *dst, ae_state *_state);
void qpbleicoptimize(convexquadraticmodel *a, sparsematrix *sparsea, ae_int_t akind, bool sparseaupper, double absasum, double absasum2, RVector *b, RVector *bndl, RVector *bndu, RVector *s, RVector *xorigin, ae_int_t n, RMatrix *cleic, ae_int_t nec, ae_int_t nic, qpbleicsettings *settings, qpbleicbuffers *sstate, bool *firstcall, RVector *xs, ae_int_t *terminationtype, ae_state *_state);
void qpbleicsettings_init(void *_p, ae_state *_state, bool make_automatic);
void qpbleicsettings_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void qpbleicsettings_free(void *_p, bool make_automatic);
void qpbleicbuffers_init(void *_p, ae_state *_state, bool make_automatic);
void qpbleicbuffers_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void qpbleicbuffers_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === VIPMSOLVER Package ===
// Depends on: (Solvers) DIRECTDENSESOLVERS
// Depends on: MINLBFGS, CQMODELS, LPQPSERV
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t m;
   ae_vector x;
   ae_vector g;
   ae_vector w;
   ae_vector t;
   ae_vector p;
   ae_vector y;
   ae_vector z;
   ae_vector v;
   ae_vector s;
   ae_vector q;
} vipmvars;
typedef struct {
   ae_vector sigma;
   ae_vector beta;
   ae_vector rho;
   ae_vector nu;
   ae_vector tau;
   ae_vector alpha;
   ae_vector gammaz;
   ae_vector gammas;
   ae_vector gammaw;
   ae_vector gammaq;
} vipmrighthandside;
typedef struct {
   bool slacksforequalityconstraints;
   ae_int_t n;
   ae_int_t nmain;
   double epsp;
   double epsd;
   double epsgap;
   bool islinear;
   ae_vector scl;
   ae_vector invscl;
   ae_vector xorigin;
   double targetscale;
   ae_vector c;
   ae_matrix denseh;
   sparsematrix sparseh;
   ae_vector diagr;
   ae_int_t hkind;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector rawbndl;
   ae_vector rawbndu;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_matrix denseafull;
   ae_matrix denseamain;
   sparsematrix sparseafull;
   sparsematrix sparseamain;
   sparsematrix combinedaslack;
   ae_vector ascales;
   ae_vector aflips;
   ae_vector b;
   ae_vector r;
   ae_vector hasr;
   ae_int_t mdense;
   ae_int_t msparse;
   vipmvars x0;
   vipmvars current;
   vipmvars best;
   vipmvars trial;
   vipmvars deltaaff;
   vipmvars deltacorr;
   ae_vector isfrozen;
   ae_vector hasgz;
   ae_vector hasts;
   ae_vector haswv;
   ae_vector haspq;
   ae_int_t repiterationscount;
   ae_int_t repncholesky;
   bool dotrace;
   bool dodetailedtrace;
   ae_int_t factorizationtype;
   bool factorizationpoweredup;
   bool factorizationpresent;
   ae_vector diagdz;
   ae_vector diagdzi;
   ae_vector diagdziri;
   ae_vector diagds;
   ae_vector diagdsi;
   ae_vector diagdsiri;
   ae_vector diagdw;
   ae_vector diagdwi;
   ae_vector diagdwir;
   ae_vector diagdq;
   ae_vector diagdqi;
   ae_vector diagdqiri;
   ae_vector diagddr;
   ae_vector diagde;
   ae_vector diagder;
   ae_matrix factdensehaug;
   ae_vector factregdhrh;
   ae_vector factinvregdzrz;
   ae_vector factregewave;
   sparsematrix factsparsekkttmpl;
   sparsematrix factsparsekkt;
   ae_vector factsparsekktpivp;
   ae_vector facttmpdiag;
   spcholanalysis ldltanalysis;
   ae_vector factsparsediagd;
   vipmrighthandside rhs;
   ae_vector rhsalphacap;
   ae_vector rhsbetacap;
   ae_vector rhsnucap;
   ae_vector rhstaucap;
   ae_vector deltaxy;
   ae_vector tmphx;
   ae_vector tmpax;
   ae_vector tmpaty;
   vipmvars zerovars;
   ae_vector dummyr;
   ae_vector tmpy;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmp2;
   ae_matrix tmpr2;
   ae_vector tmplaggrad;
   ae_vector tmpi;
   sparsematrix tmpsparse0;
} vipmstate;

void vipminitdense(vipmstate *state, RVector *s, RVector *xorigin, ae_int_t n, ae_state *_state);
void vipminitdensewithslacks(vipmstate *state, RVector *s, RVector *xorigin, ae_int_t nmain, ae_int_t n, ae_state *_state);
void vipminitsparse(vipmstate *state, RVector *s, RVector *xorigin, ae_int_t n, ae_state *_state);
void vipmsetquadraticlinear(vipmstate *state, RMatrix *denseh, sparsematrix *sparseh, ae_int_t hkind, bool isupper, RVector *c, ae_state *_state);
void vipmsetconstraints(vipmstate *state, RVector *bndl, RVector *bndu, sparsematrix *sparsea, ae_int_t msparse, RMatrix *densea, ae_int_t mdense, RVector *cl, RVector *cu, ae_state *_state);
void vipmsetcond(vipmstate *state, double epsp, double epsd, double epsgap, ae_state *_state);
void vipmoptimize(vipmstate *state, bool dropbigbounds, RVector *xs, RVector *lagbc, RVector *laglc, ae_int_t *terminationtype, ae_state *_state);
void vipmvars_init(void *_p, ae_state *_state, bool make_automatic);
void vipmvars_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void vipmvars_free(void *_p, bool make_automatic);
void vipmrighthandside_init(void *_p, ae_state *_state, bool make_automatic);
void vipmrighthandside_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void vipmrighthandside_free(void *_p, bool make_automatic);
void vipmstate_init(void *_p, ae_state *_state, bool make_automatic);
void vipmstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void vipmstate_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === MINQP Package ===
// Depends on: QPDENSEAULSOLVER, QPBLEICSOLVER, VIPMSOLVER
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   qqpsettings qqpsettingsuser;
   qpbleicsettings qpbleicsettingsuser;
   qpdenseaulsettings qpdenseaulsettingsuser;
   double veps;
   bool dbgskipconstraintnormalization;
   ae_int_t algokind;
   ae_int_t akind;
   convexquadraticmodel a;
   sparsematrix sparsea;
   bool sparseaupper;
   double absamax;
   double absasum;
   double absasum2;
   ae_vector b;
   ae_vector bndl;
   ae_vector bndu;
   ae_int_t stype;
   ae_vector s;
   ae_vector havebndl;
   ae_vector havebndu;
   ae_vector xorigin;
   ae_vector startx;
   bool havex;
   ae_matrix densec;
   sparsematrix sparsec;
   ae_vector cl;
   ae_vector cu;
   ae_int_t mdense;
   ae_int_t msparse;
   ae_vector xs;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repncholesky;
   ae_int_t repnmv;
   ae_int_t repterminationtype;
   ae_vector replagbc;
   ae_vector replaglc;
   ae_vector effectives;
   ae_vector tmp0;
   ae_matrix ecleic;
   ae_vector elaglc;
   ae_vector elagmlt;
   ae_vector elagidx;
   ae_matrix dummyr2;
   sparsematrix dummysparse;
   ae_matrix tmpr2;
   ae_vector wrkbndl;
   ae_vector wrkbndu;
   ae_vector wrkcl;
   ae_vector wrkcu;
   ae_matrix wrkdensec;
   sparsematrix wrksparsec;
   bool qpbleicfirstcall;
   qpbleicbuffers qpbleicbuf;
   qqpbuffers qqpbuf;
   qpdenseaulbuffers qpdenseaulbuf;
   vipmstate vsolver;
} minqpstate;
typedef struct {
   ae_int_t inneriterationscount;
   ae_int_t outeriterationscount;
   ae_int_t nmv;
   ae_int_t ncholesky;
   ae_int_t terminationtype;
   ae_vector lagbc;
   ae_vector laglc;
} minqpreport;

void minqpcreate(ae_int_t n, minqpstate *state, ae_state *_state);
void minqpsetlinearterm(minqpstate *state, RVector *b, ae_state *_state);
void minqpsetquadraticterm(minqpstate *state, RMatrix *a, bool isupper, ae_state *_state);
void minqpsetquadratictermsparse(minqpstate *state, sparsematrix *a, bool isupper, ae_state *_state);
void minqpsetstartingpoint(minqpstate *state, RVector *x, ae_state *_state);
void minqpsetorigin(minqpstate *state, RVector *xorigin, ae_state *_state);
void minqpsetscale(minqpstate *state, RVector *s, ae_state *_state);
void minqpsetscaleautodiag(minqpstate *state, ae_state *_state);
void minqpsetalgobleic(minqpstate *state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
void minqpsetalgodenseaul(minqpstate *state, double epsx, double rho, ae_int_t itscnt, ae_state *_state);
void minqpsetalgodenseipm(minqpstate *state, double eps, ae_state *_state);
void minqpsetalgosparseipm(minqpstate *state, double eps, ae_state *_state);
void minqpsetalgoquickqp(minqpstate *state, double epsg, double epsf, double epsx, ae_int_t maxouterits, bool usenewton, ae_state *_state);
void minqpsetbc(minqpstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minqpsetbcall(minqpstate *state, double bndl, double bndu, ae_state *_state);
void minqpsetbci(minqpstate *state, ae_int_t i, double bndl, double bndu, ae_state *_state);
void minqpsetlc(minqpstate *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void minqpsetlcsparse(minqpstate *state, sparsematrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void minqpsetlcmixed(minqpstate *state, sparsematrix *sparsec, ZVector *sparsect, ae_int_t sparsek, RMatrix *densec, ZVector *densect, ae_int_t densek, ae_state *_state);
void minqpsetlcmixedlegacy(minqpstate *state, RMatrix *densec, ZVector *densect, ae_int_t densek, sparsematrix *sparsec, ZVector *sparsect, ae_int_t sparsek, ae_state *_state);
void minqpsetlc2dense(minqpstate *state, RMatrix *a, RVector *al, RVector *au, ae_int_t k, ae_state *_state);
void minqpsetlc2(minqpstate *state, sparsematrix *a, RVector *al, RVector *au, ae_int_t k, ae_state *_state);
void minqpsetlc2mixed(minqpstate *state, sparsematrix *sparsea, ae_int_t ksparse, RMatrix *densea, ae_int_t kdense, RVector *al, RVector *au, ae_state *_state);
void minqpaddlc2dense(minqpstate *state, RVector *a, double al, double au, ae_state *_state);
void minqpaddlc2(minqpstate *state, ZVector *idxa, RVector *vala, ae_int_t nnz, double al, double au, ae_state *_state);
void minqpaddlc2sparsefromdense(minqpstate *state, RVector *da, double al, double au, ae_state *_state);
void minqpoptimize(minqpstate *state, ae_state *_state);
void minqpresults(minqpstate *state, RVector *x, minqpreport *rep, ae_state *_state);
void minqpresultsbuf(minqpstate *state, RVector *x, minqpreport *rep, ae_state *_state);
void minqpsetlineartermfast(minqpstate *state, RVector *b, ae_state *_state);
void minqpsetquadratictermfast(minqpstate *state, RMatrix *a, bool isupper, double s, ae_state *_state);
void minqprewritediagonal(minqpstate *state, RVector *s, ae_state *_state);
void minqpsetstartingpointfast(minqpstate *state, RVector *x, ae_state *_state);
void minqpsetoriginfast(minqpstate *state, RVector *xorigin, ae_state *_state);
void minqpstate_init(void *_p, ae_state *_state, bool make_automatic);
void minqpstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minqpstate_free(void *_p, bool make_automatic);
void minqpreport_init(void *_p, ae_state *_state, bool make_automatic);
void minqpreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minqpreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minqpstate, );
DecClass(minqpreport, ae_int_t &inneriterationscount; ae_int_t &outeriterationscount; ae_int_t &nmv; ae_int_t &ncholesky; ae_int_t &terminationtype; real_1d_array lagbc; real_1d_array laglc;);

void minqpcreate(const ae_int_t n, minqpstate &state, const xparams _xparams = xdefault);
void minqpsetlinearterm(const minqpstate &state, const real_1d_array &b, const xparams _xparams = xdefault);
void minqpsetquadraticterm(const minqpstate &state, const real_2d_array &a, const bool isupper, const xparams _xparams = xdefault);
void minqpsetquadraticterm(const minqpstate &state, const real_2d_array &a, const xparams _xparams = xdefault);
void minqpsetquadratictermsparse(const minqpstate &state, const sparsematrix &a, const bool isupper, const xparams _xparams = xdefault);
void minqpsetstartingpoint(const minqpstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void minqpsetorigin(const minqpstate &state, const real_1d_array &xorigin, const xparams _xparams = xdefault);
void minqpsetscale(const minqpstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minqpsetscaleautodiag(const minqpstate &state, const xparams _xparams = xdefault);
void minqpsetalgobleic(const minqpstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minqpsetalgodenseaul(const minqpstate &state, const double epsx, const double rho, const ae_int_t itscnt, const xparams _xparams = xdefault);
void minqpsetalgodenseipm(const minqpstate &state, const double eps, const xparams _xparams = xdefault);
void minqpsetalgosparseipm(const minqpstate &state, const double eps, const xparams _xparams = xdefault);
void minqpsetalgoquickqp(const minqpstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxouterits, const bool usenewton, const xparams _xparams = xdefault);
void minqpsetbc(const minqpstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minqpsetbcall(const minqpstate &state, const double bndl, const double bndu, const xparams _xparams = xdefault);
void minqpsetbci(const minqpstate &state, const ae_int_t i, const double bndl, const double bndu, const xparams _xparams = xdefault);
void minqpsetlc(const minqpstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minqpsetlc(const minqpstate &state, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = xdefault);
void minqpsetlcsparse(const minqpstate &state, const sparsematrix &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minqpsetlcmixed(const minqpstate &state, const sparsematrix &sparsec, const integer_1d_array &sparsect, const ae_int_t sparsek, const real_2d_array &densec, const integer_1d_array &densect, const ae_int_t densek, const xparams _xparams = xdefault);
void minqpsetlcmixedlegacy(const minqpstate &state, const real_2d_array &densec, const integer_1d_array &densect, const ae_int_t densek, const sparsematrix &sparsec, const integer_1d_array &sparsect, const ae_int_t sparsek, const xparams _xparams = xdefault);
void minqpsetlc2dense(const minqpstate &state, const real_2d_array &a, const real_1d_array &al, const real_1d_array &au, const ae_int_t k, const xparams _xparams = xdefault);
void minqpsetlc2dense(const minqpstate &state, const real_2d_array &a, const real_1d_array &al, const real_1d_array &au, const xparams _xparams = xdefault);
void minqpsetlc2(const minqpstate &state, const sparsematrix &a, const real_1d_array &al, const real_1d_array &au, const ae_int_t k, const xparams _xparams = xdefault);
void minqpsetlc2mixed(const minqpstate &state, const sparsematrix &sparsea, const ae_int_t ksparse, const real_2d_array &densea, const ae_int_t kdense, const real_1d_array &al, const real_1d_array &au, const xparams _xparams = xdefault);
void minqpaddlc2dense(const minqpstate &state, const real_1d_array &a, const double al, const double au, const xparams _xparams = xdefault);
void minqpaddlc2(const minqpstate &state, const integer_1d_array &idxa, const real_1d_array &vala, const ae_int_t nnz, const double al, const double au, const xparams _xparams = xdefault);
void minqpaddlc2sparsefromdense(const minqpstate &state, const real_1d_array &da, const double al, const double au, const xparams _xparams = xdefault);
void minqpoptimize(const minqpstate &state, const xparams _xparams = xdefault);
void minqpresults(const minqpstate &state, real_1d_array &x, minqpreport &rep, const xparams _xparams = xdefault);
void minqpresultsbuf(const minqpstate &state, real_1d_array &x, minqpreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MINLM Package ===
// Depends on: MINQP
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t m;
   double stpmax;
   ae_int_t modelage;
   ae_int_t maxmodelage;
   bool hasfi;
   double epsx;
   ae_vector x;
   double f;
   ae_vector fi;
   bool needf;
   bool needfi;
   double fbase;
   ae_vector modeldiag;
   ae_vector xbase;
   ae_vector fibase;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector havebndl;
   ae_vector havebndu;
   ae_vector s;
   rcommstate rstate;
   ae_vector xdir;
   ae_vector choleskybuf;
   ae_vector tmp0;
   ae_vector tmpct;
   double actualdecrease;
   double predicteddecrease;
   minqpstate qpstate;
   minqpreport qprep;
   sparsematrix tmpsp;
} minlmstepfinder;
typedef struct {
   ae_int_t n;
   ae_int_t m;
   double diffstep;
   double epsx;
   ae_int_t maxits;
   bool xrep;
   double stpmax;
   ae_int_t maxmodelage;
   bool makeadditers;
   ae_vector x;
   double f;
   ae_vector fi;
   ae_matrix j;
   ae_matrix h;
   ae_vector g;
   bool needf;
   bool needfg;
   bool needfgh;
   bool needfij;
   bool needfi;
   bool xupdated;
   bool userterminationneeded;
   ae_int_t algomode;
   bool hasf;
   bool hasfi;
   bool hasg;
   ae_vector xbase;
   double fbase;
   ae_vector fibase;
   ae_vector gbase;
   ae_matrix quadraticmodel;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector havebndl;
   ae_vector havebndu;
   ae_vector s;
   ae_matrix cleic;
   ae_int_t nec;
   ae_int_t nic;
   double lambdav;
   double nu;
   ae_int_t modelage;
   ae_vector xnew;
   ae_vector xdir;
   ae_vector deltax;
   ae_vector deltaf;
   bool deltaxready;
   bool deltafready;
   smoothnessmonitor smonitor;
   double teststep;
   ae_vector lastscaleused;
   ae_int_t repiterationscount;
   ae_int_t repterminationtype;
   ae_int_t repnfunc;
   ae_int_t repnjac;
   ae_int_t repngrad;
   ae_int_t repnhess;
   ae_int_t repncholesky;
   rcommstate rstate;
   ae_vector choleskybuf;
   ae_vector tmp0;
   double actualdecrease;
   double predicteddecrease;
   double xm1;
   double xp1;
   ae_vector fm1;
   ae_vector fp1;
   ae_vector fc1;
   ae_vector gm1;
   ae_vector gp1;
   ae_vector gc1;
   minlbfgsstate internalstate;
   minlbfgsreport internalrep;
   minqpstate qpstate;
   minqpreport qprep;
   minlmstepfinder finderstate;
} minlmstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t terminationtype;
   ae_int_t nfunc;
   ae_int_t njac;
   ae_int_t ngrad;
   ae_int_t nhess;
   ae_int_t ncholesky;
} minlmreport;

void minlmcreatevj(ae_int_t n, ae_int_t m, RVector *x, minlmstate *state, ae_state *_state);
void minlmcreatev(ae_int_t n, ae_int_t m, RVector *x, double diffstep, minlmstate *state, ae_state *_state);
void minlmcreatefgh(ae_int_t n, RVector *x, minlmstate *state, ae_state *_state);
void minlmsetcond(minlmstate *state, double epsx, ae_int_t maxits, ae_state *_state);
void minlmsetxrep(minlmstate *state, bool needxrep, ae_state *_state);
void minlmsetstpmax(minlmstate *state, double stpmax, ae_state *_state);
void minlmsetscale(minlmstate *state, RVector *s, ae_state *_state);
void minlmsetbc(minlmstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minlmsetlc(minlmstate *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void minlmsetacctype(minlmstate *state, ae_int_t acctype, ae_state *_state);
bool minlmiteration(minlmstate *state, ae_state *_state);
void minlmoptguardgradient(minlmstate *state, double teststep, ae_state *_state);
void minlmoptguardresults(minlmstate *state, optguardreport *rep, ae_state *_state);
void minlmresults(minlmstate *state, RVector *x, minlmreport *rep, ae_state *_state);
void minlmresultsbuf(minlmstate *state, RVector *x, minlmreport *rep, ae_state *_state);
void minlmrestartfrom(minlmstate *state, RVector *x, ae_state *_state);
void minlmrequesttermination(minlmstate *state, ae_state *_state);
void minlmcreatevgj(ae_int_t n, ae_int_t m, RVector *x, minlmstate *state, ae_state *_state);
void minlmcreatefgj(ae_int_t n, ae_int_t m, RVector *x, minlmstate *state, ae_state *_state);
void minlmcreatefj(ae_int_t n, ae_int_t m, RVector *x, minlmstate *state, ae_state *_state);
void minlmstepfinder_init(void *_p, ae_state *_state, bool make_automatic);
void minlmstepfinder_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlmstepfinder_free(void *_p, bool make_automatic);
void minlmstate_init(void *_p, ae_state *_state, bool make_automatic);
void minlmstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlmstate_free(void *_p, bool make_automatic);
void minlmreport_init(void *_p, ae_state *_state, bool make_automatic);
void minlmreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlmreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minlmstate, bool &needf; bool &needfg; bool &needfgh; bool &needfi; bool &needfij; bool &xupdated; double &f; real_1d_array fi; real_1d_array g; real_2d_array h; real_2d_array j; real_1d_array x;);
DecClass(minlmreport, ae_int_t &iterationscount; ae_int_t &terminationtype; ae_int_t &nfunc; ae_int_t &njac; ae_int_t &ngrad; ae_int_t &nhess; ae_int_t &ncholesky;);

void minlmcreatevj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatevj(const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatev(const ae_int_t n, const ae_int_t m, const real_1d_array &x, const double diffstep, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatev(const ae_int_t m, const real_1d_array &x, const double diffstep, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatefgh(const ae_int_t n, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatefgh(const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmsetcond(const minlmstate &state, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minlmsetxrep(const minlmstate &state, const bool needxrep, const xparams _xparams = xdefault);
void minlmsetstpmax(const minlmstate &state, const double stpmax, const xparams _xparams = xdefault);
void minlmsetscale(const minlmstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minlmsetbc(const minlmstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minlmsetlc(const minlmstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minlmsetlc(const minlmstate &state, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = xdefault);
void minlmsetacctype(const minlmstate &state, const ae_int_t acctype, const xparams _xparams = xdefault);
bool minlmiteration(const minlmstate &state, const xparams _xparams = xdefault);
void minlmoptimize(minlmstate &state, void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlmoptimize(minlmstate &state, void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlmoptimize(minlmstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*hess)(const real_1d_array &x, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlmoptimize(minlmstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlmoptimize(minlmstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minlmoptguardgradient(const minlmstate &state, const double teststep, const xparams _xparams = xdefault);
void minlmoptguardresults(const minlmstate &state, optguardreport &rep, const xparams _xparams = xdefault);
void minlmresults(const minlmstate &state, real_1d_array &x, minlmreport &rep, const xparams _xparams = xdefault);
void minlmresultsbuf(const minlmstate &state, real_1d_array &x, minlmreport &rep, const xparams _xparams = xdefault);
void minlmrestartfrom(const minlmstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void minlmrequesttermination(const minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatevgj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatevgj(const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatefgj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatefgj(const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatefj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
void minlmcreatefj(const ae_int_t m, const real_1d_array &x, minlmstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MINCG Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: OPTSERV
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxits;
   double stpmax;
   double suggestedstep;
   bool xrep;
   bool drep;
   ae_int_t cgtype;
   ae_int_t prectype;
   ae_vector diagh;
   ae_vector diaghl2;
   ae_matrix vcorr;
   ae_int_t vcnt;
   ae_vector s;
   double diffstep;
   ae_int_t nfev;
   ae_int_t mcstage;
   ae_int_t k;
   ae_vector xk;
   ae_vector dk;
   ae_vector xn;
   ae_vector dn;
   ae_vector d;
   double fold;
   double stp;
   double curstpmax;
   ae_vector yk;
   double lastgoodstep;
   double lastscaledstep;
   ae_int_t mcinfo;
   bool innerresetneeded;
   bool terminationneeded;
   double trimthreshold;
   ae_vector xbase;
   ae_int_t rstimer;
   ae_vector x;
   double f;
   ae_vector g;
   bool needf;
   bool needfg;
   bool xupdated;
   bool algpowerup;
   bool lsstart;
   bool lsend;
   bool userterminationneeded;
   rcommstate rstate;
   ae_int_t repiterationscount;
   ae_int_t repnfev;
   ae_int_t repterminationtype;
   ae_int_t debugrestartscount;
   linminstate lstate;
   double fbase;
   double fm2;
   double fm1;
   double fp1;
   double fp2;
   double betahs;
   double betady;
   ae_vector work0;
   ae_vector work1;
   ae_vector invs;
   double teststep;
   ae_int_t smoothnessguardlevel;
   smoothnessmonitor smonitor;
   ae_vector lastscaleused;
} mincgstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
} mincgreport;

void mincgcreate(ae_int_t n, RVector *x, mincgstate *state, ae_state *_state);
void mincgcreatef(ae_int_t n, RVector *x, double diffstep, mincgstate *state, ae_state *_state);
void mincgsetcond(mincgstate *state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
void mincgsetscale(mincgstate *state, RVector *s, ae_state *_state);
void mincgsetxrep(mincgstate *state, bool needxrep, ae_state *_state);
void mincgsetdrep(mincgstate *state, bool needdrep, ae_state *_state);
void mincgsetcgtype(mincgstate *state, ae_int_t cgtype, ae_state *_state);
void mincgsetstpmax(mincgstate *state, double stpmax, ae_state *_state);
void mincgsuggeststep(mincgstate *state, double stp, ae_state *_state);
double mincglastgoodstep(mincgstate *state, ae_state *_state);
void mincgsetprecdefault(mincgstate *state, ae_state *_state);
void mincgsetprecdiag(mincgstate *state, RVector *d, ae_state *_state);
void mincgsetprecscale(mincgstate *state, ae_state *_state);
bool mincgiteration(mincgstate *state, ae_state *_state);
void mincgoptguardgradient(mincgstate *state, double teststep, ae_state *_state);
void mincgoptguardsmoothness(mincgstate *state, ae_int_t level, ae_state *_state);
void mincgoptguardresults(mincgstate *state, optguardreport *rep, ae_state *_state);
void mincgoptguardnonc1test0results(mincgstate *state, optguardnonc1test0report *strrep, optguardnonc1test0report *lngrep, ae_state *_state);
void mincgoptguardnonc1test1results(mincgstate *state, optguardnonc1test1report *strrep, optguardnonc1test1report *lngrep, ae_state *_state);
void mincgresults(mincgstate *state, RVector *x, mincgreport *rep, ae_state *_state);
void mincgresultsbuf(mincgstate *state, RVector *x, mincgreport *rep, ae_state *_state);
void mincgrestartfrom(mincgstate *state, RVector *x, ae_state *_state);
void mincgrequesttermination(mincgstate *state, ae_state *_state);
void mincgsetprecdiagfast(mincgstate *state, RVector *d, ae_state *_state);
void mincgsetpreclowrankfast(mincgstate *state, RVector *d1, RVector *c, RMatrix *v, ae_int_t vcnt, ae_state *_state);
void mincgsetprecvarpart(mincgstate *state, RVector *d2, ae_state *_state);
void mincgstate_init(void *_p, ae_state *_state, bool make_automatic);
void mincgstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mincgstate_free(void *_p, bool make_automatic);
void mincgreport_init(void *_p, ae_state *_state, bool make_automatic);
void mincgreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void mincgreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(mincgstate, bool &needf; bool &needfg; bool &xupdated; double &f; real_1d_array g; real_1d_array x;);
DecClass(mincgreport, ae_int_t &iterationscount; ae_int_t &nfev; ae_int_t &terminationtype;);

void mincgcreate(const ae_int_t n, const real_1d_array &x, mincgstate &state, const xparams _xparams = xdefault);
void mincgcreate(const real_1d_array &x, mincgstate &state, const xparams _xparams = xdefault);
void mincgcreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, mincgstate &state, const xparams _xparams = xdefault);
void mincgcreatef(const real_1d_array &x, const double diffstep, mincgstate &state, const xparams _xparams = xdefault);
void mincgsetcond(const mincgstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void mincgsetscale(const mincgstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void mincgsetxrep(const mincgstate &state, const bool needxrep, const xparams _xparams = xdefault);
void mincgsetcgtype(const mincgstate &state, const ae_int_t cgtype, const xparams _xparams = xdefault);
void mincgsetstpmax(const mincgstate &state, const double stpmax, const xparams _xparams = xdefault);
void mincgsuggeststep(const mincgstate &state, const double stp, const xparams _xparams = xdefault);
void mincgsetprecdefault(const mincgstate &state, const xparams _xparams = xdefault);
void mincgsetprecdiag(const mincgstate &state, const real_1d_array &d, const xparams _xparams = xdefault);
void mincgsetprecscale(const mincgstate &state, const xparams _xparams = xdefault);
bool mincgiteration(const mincgstate &state, const xparams _xparams = xdefault);
void mincgoptimize(mincgstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void mincgoptimize(mincgstate &state, void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void mincgoptguardgradient(const mincgstate &state, const double teststep, const xparams _xparams = xdefault);
void mincgoptguardsmoothness(const mincgstate &state, const ae_int_t level, const xparams _xparams = xdefault);
void mincgoptguardsmoothness(const mincgstate &state, const xparams _xparams = xdefault);
void mincgoptguardresults(const mincgstate &state, optguardreport &rep, const xparams _xparams = xdefault);
void mincgoptguardnonc1test0results(const mincgstate &state, optguardnonc1test0report &strrep, optguardnonc1test0report &lngrep, const xparams _xparams = xdefault);
void mincgoptguardnonc1test1results(const mincgstate &state, optguardnonc1test1report &strrep, optguardnonc1test1report &lngrep, const xparams _xparams = xdefault);
void mincgresults(const mincgstate &state, real_1d_array &x, mincgreport &rep, const xparams _xparams = xdefault);
void mincgresultsbuf(const mincgstate &state, real_1d_array &x, mincgreport &rep, const xparams _xparams = xdefault);
void mincgrestartfrom(const mincgstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void mincgrequesttermination(const mincgstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === NLCSQP Package ===
// Depends on: VIPMSOLVER
namespace alglib_impl {
typedef struct {
   ae_int_t algokind;
   vipmstate ipmsolver;
   ae_vector curb;
   ae_vector curbndl;
   ae_vector curbndu;
   ae_vector cural;
   ae_vector curau;
   sparsematrix sparserawlc;
   sparsematrix sparseefflc;
   ae_vector d0;
   ae_matrix h;
   ae_matrix densedummy;
   sparsematrix sparsedummy;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmp2;
   ae_vector sk;
   ae_vector yk;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_vector hasal;
   ae_vector hasau;
   ae_matrix activea;
   ae_vector activerhs;
   ae_vector activeidx;
   ae_int_t activesetsize;
} minsqpsubsolver;
typedef struct {
   ae_vector sclagtmp0;
   ae_vector sclagtmp1;
} minsqptmplagrangian;
typedef struct {
   ae_vector mftmp0;
} minsqptmpmerit;
typedef struct {
   ae_int_t n;
   ae_int_t nec;
   ae_int_t nic;
   ae_int_t nlec;
   ae_int_t nlic;
   ae_vector d;
   ae_vector dx;
   ae_vector stepkx;
   ae_vector stepkxc;
   ae_vector stepkxn;
   ae_vector stepkfi;
   ae_vector stepkfic;
   ae_vector stepkfin;
   ae_matrix stepkj;
   ae_matrix stepkjc;
   ae_matrix stepkjn;
   ae_vector lagmult;
   ae_vector dummylagmult;
   ae_vector penalties;
   minsqptmpmerit tmpmerit;
   minsqptmplagrangian tmplagrangianfg;
   ae_vector stepklaggrad;
   ae_vector stepknlaggrad;
   ae_int_t status;
   bool increasebigc;
   rcommstate rmeritphasestate;
} minsqpmeritphasestate;
typedef struct {
   ae_int_t n;
   ae_int_t nec;
   ae_int_t nic;
   ae_int_t nlec;
   ae_int_t nlic;
   ae_vector s;
   ae_matrix scaledcleic;
   ae_vector lcsrcidx;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_vector scaledbndl;
   ae_vector scaledbndu;
   double epsx;
   ae_int_t maxits;
   ae_vector x;
   ae_vector fi;
   ae_matrix j;
   double f;
   bool needfij;
   bool xupdated;
   minsqpmeritphasestate meritstate;
   double bigc;
   double trustrad;
   ae_int_t trustradstagnationcnt;
   ae_int_t fstagnationcnt;
   ae_vector step0x;
   ae_vector stepkx;
   ae_vector backupx;
   ae_vector step0fi;
   ae_vector stepkfi;
   ae_vector backupfi;
   ae_matrix step0j;
   ae_matrix stepkj;
   bool haslagmult;
   ae_vector meritlagmult;
   ae_vector dummylagmult;
   ae_matrix abslagmemory;
   ae_vector fscales;
   ae_vector tracegamma;
   minsqpsubsolver subsolver;
   minsqptmpmerit tmpmerit;
   ae_int_t repsimplexiterations;
   ae_int_t repsimplexiterations1;
   ae_int_t repsimplexiterations2;
   ae_int_t repsimplexiterations3;
   ae_int_t repiterationscount;
   ae_int_t repterminationtype;
   double repbcerr;
   ae_int_t repbcidx;
   double replcerr;
   ae_int_t replcidx;
   double repnlcerr;
   ae_int_t repnlcidx;
   rcommstate rstate;
} minsqpstate;

void minsqpinitbuf(RVector *bndl, RVector *bndu, RVector *s, RVector *x0, ae_int_t n, RMatrix *cleic, ZVector *lcsrcidx, ae_int_t nec, ae_int_t nic, ae_int_t nlec, ae_int_t nlic, double epsx, ae_int_t maxits, minsqpstate *state, ae_state *_state);
bool minsqpiteration(minsqpstate *state, smoothnessmonitor *smonitor, bool userterminationneeded, ae_state *_state);
void minsqpsubsolver_init(void *_p, ae_state *_state, bool make_automatic);
void minsqpsubsolver_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minsqpsubsolver_free(void *_p, bool make_automatic);
void minsqptmplagrangian_init(void *_p, ae_state *_state, bool make_automatic);
void minsqptmplagrangian_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minsqptmplagrangian_free(void *_p, bool make_automatic);
void minsqptmpmerit_init(void *_p, ae_state *_state, bool make_automatic);
void minsqptmpmerit_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minsqptmpmerit_free(void *_p, bool make_automatic);
void minsqpmeritphasestate_init(void *_p, ae_state *_state, bool make_automatic);
void minsqpmeritphasestate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minsqpmeritphasestate_free(void *_p, bool make_automatic);
void minsqpstate_init(void *_p, ae_state *_state, bool make_automatic);
void minsqpstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minsqpstate_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === LPQPPRESOLVE Package ===
// Depends on: (LinAlg) SPARSE
namespace alglib_impl {
typedef struct {
   ae_int_t newn;
   ae_int_t oldn;
   ae_int_t newm;
   ae_int_t oldm;
   ae_vector rawbndl;
   ae_vector rawbndu;
   ae_vector colscales;
   ae_vector rowscales;
   double costscale;
   ae_vector c;
   ae_vector bndl;
   ae_vector bndu;
   sparsematrix sparsea;
   ae_vector al;
   ae_vector au;
} presolveinfo;

void presolvenonescaleuser(RVector *s, RVector *c, RVector *bndl, RVector *bndu, ae_int_t n, sparsematrix *sparsea, RVector *al, RVector *au, ae_int_t k, presolveinfo *info, ae_state *_state);
void presolvebwd(presolveinfo *info, RVector *x, ZVector *stats, RVector *lagbc, RVector *laglc, ae_state *_state);
void presolveinfo_init(void *_p, ae_state *_state, bool make_automatic);
void presolveinfo_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void presolveinfo_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === REVISEDDUALSIMPLEX Package ===
// Depends on: (LinAlg) TRFAC
// Depends on: LPQPPRESOLVE
namespace alglib_impl {
typedef struct {
   double pivottol;
   double perturbmag;
   ae_int_t maxtrfage;
   ae_int_t trftype;
   ae_int_t ratiotest;
   ae_int_t pricing;
   ae_int_t shifting;
   double xtolabs;
   double xtolrelabs;
   double dtolabs;
} dualsimplexsettings;
typedef struct {
   ae_int_t n;
   ae_int_t k;
   ae_vector idx;
   ae_vector vals;
   ae_vector dense;
} dssvector;
typedef struct {
   ae_int_t ns;
   ae_int_t m;
   ae_vector idx;
   ae_vector nidx;
   ae_vector isbasic;
   ae_int_t trftype;
   bool isvalidtrf;
   ae_int_t trfage;
   ae_matrix denselu;
   sparsematrix sparsel;
   sparsematrix sparseu;
   sparsematrix sparseut;
   ae_vector rowpermbwd;
   ae_vector colpermbwd;
   ae_vector densepfieta;
   ae_vector densemu;
   ae_vector rk;
   ae_vector dk;
   ae_vector dseweights;
   bool dsevalid;
   double eminu;
   ae_int_t statfact;
   ae_int_t statupdt;
   double statoffdiag;
   ae_vector wtmp0;
   ae_vector wtmp1;
   ae_vector wtmp2;
   ae_vector nrs;
   ae_vector tcinvidx;
   ae_matrix denselu2;
   ae_vector densep2;
   ae_vector densep2c;
   sparsematrix sparselu1;
   sparsematrix sparselu2;
   sluv2buffer lubuf2;
   ae_vector tmpi;
   ae_vector utmp0;
   ae_vector utmpi;
   sparsematrix sparseludbg;
} dualsimplexbasis;
typedef struct {
   ae_int_t ns;
   ae_int_t m;
   ae_vector rawc;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector bndt;
   ae_vector xa;
   ae_vector d;
   ae_int_t state;
   ae_vector xb;
   ae_vector bndlb;
   ae_vector bndub;
   ae_vector bndtb;
   ae_vector bndtollb;
   ae_vector bndtolub;
   ae_vector effc;
} dualsimplexsubproblem;
typedef struct {
   ae_vector rowscales;
   ae_vector rawbndl;
   ae_vector rawbndu;
   ae_int_t ns;
   ae_int_t m;
   sparsematrix a;
   sparsematrix at;
   dualsimplexbasis basis;
   dualsimplexsubproblem primary;
   dualsimplexsubproblem phase1;
   dualsimplexsubproblem phase3;
   ae_vector repx;
   ae_vector replagbc;
   ae_vector replaglc;
   ae_vector repstats;
   ae_int_t repterminationtype;
   ae_int_t repiterationscount;
   ae_int_t repiterationscount1;
   ae_int_t repiterationscount2;
   ae_int_t repiterationscount3;
   ae_int_t repphase1time;
   ae_int_t repphase2time;
   ae_int_t repphase3time;
   ae_int_t repdualpricingtime;
   ae_int_t repdualbtrantime;
   ae_int_t repdualpivotrowtime;
   ae_int_t repdualratiotesttime;
   ae_int_t repdualftrantime;
   ae_int_t repdualupdatesteptime;
   double repfillpivotrow;
   ae_int_t repfillpivotrowcnt;
   double repfillrhor;
   ae_int_t repfillrhorcnt;
   double repfilldensemu;
   ae_int_t repfilldensemucnt;
   bool dotrace;
   bool dodetailedtrace;
   bool dotimers;
   ae_vector btrantmp0;
   ae_vector btrantmp1;
   ae_vector btrantmp2;
   ae_vector ftrantmp0;
   ae_vector ftrantmp1;
   ae_vector possibleflips;
   ae_int_t possibleflipscnt;
   ae_vector dfctmp0;
   ae_vector dfctmp1;
   ae_vector dfctmp2;
   ae_vector ustmpi;
   apbuffers xydsbuf;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector tmp2;
   dssvector alphar;
   dssvector rhor;
   ae_vector tau;
   ae_vector alphaq;
   ae_vector alphaqim;
   ae_vector eligiblealphar;
   ae_vector harrisset;
} dualsimplexstate;

void dsssettingsinit(dualsimplexsettings *settings, ae_state *_state);
void dssinit(ae_int_t n, dualsimplexstate *s, ae_state *_state);
void dsssetproblem(dualsimplexstate *state, RVector *c, RVector *bndl, RVector *bndu, RMatrix *densea, sparsematrix *sparsea, ae_int_t akind, RVector *al, RVector *au, ae_int_t k, dualsimplexbasis *proposedbasis, ae_int_t basisinittype, dualsimplexsettings *settings, ae_state *_state);
void dssexportbasis(dualsimplexstate *state, dualsimplexbasis *basis, ae_state *_state);
void dssoptimize(dualsimplexstate *state, dualsimplexsettings *settings, ae_state *_state);
void dualsimplexsettings_init(void *_p, ae_state *_state, bool make_automatic);
void dualsimplexsettings_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dualsimplexsettings_free(void *_p, bool make_automatic);
void dssvector_init(void *_p, ae_state *_state, bool make_automatic);
void dssvector_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dssvector_free(void *_p, bool make_automatic);
void dualsimplexbasis_init(void *_p, ae_state *_state, bool make_automatic);
void dualsimplexbasis_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dualsimplexbasis_free(void *_p, bool make_automatic);
void dualsimplexsubproblem_init(void *_p, ae_state *_state, bool make_automatic);
void dualsimplexsubproblem_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dualsimplexsubproblem_free(void *_p, bool make_automatic);
void dualsimplexstate_init(void *_p, ae_state *_state, bool make_automatic);
void dualsimplexstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void dualsimplexstate_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === MINLP Package ===
// Depends on: VIPMSOLVER, REVISEDDUALSIMPLEX
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   ae_int_t algokind;
   double ipmlambda;
   ae_vector s;
   ae_vector c;
   ae_vector bndl;
   ae_vector bndu;
   ae_int_t m;
   sparsematrix a;
   ae_vector al;
   ae_vector au;
   ae_vector xs;
   ae_vector lagbc;
   ae_vector laglc;
   ae_vector cs;
   double repf;
   double repprimalerror;
   double repdualerror;
   double repslackerror;
   ae_int_t repiterationscount;
   ae_int_t repterminationtype;
   ae_int_t repn;
   ae_int_t repm;
   double dsseps;
   double ipmeps;
   dualsimplexstate dss;
   vipmstate ipm;
   ae_vector adddtmpi;
   ae_vector adddtmpr;
   ae_vector tmpax;
   ae_vector tmpg;
   presolveinfo presolver;
   ae_vector zeroorigin;
   ae_vector units;
   sparsematrix ipmquadratic;
} minlpstate;
typedef struct {
   double f;
   ae_vector lagbc;
   ae_vector laglc;
   ae_vector y;
   ae_vector stats;
   double primalerror;
   double dualerror;
   double slackerror;
   ae_int_t iterationscount;
   ae_int_t terminationtype;
} minlpreport;

void minlpcreate(ae_int_t n, minlpstate *state, ae_state *_state);
void minlpsetalgodss(minlpstate *state, double eps, ae_state *_state);
void minlpsetalgoipm(minlpstate *state, double eps, ae_state *_state);
void minlpsetcost(minlpstate *state, RVector *c, ae_state *_state);
void minlpsetscale(minlpstate *state, RVector *s, ae_state *_state);
void minlpsetbc(minlpstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minlpsetbcall(minlpstate *state, double bndl, double bndu, ae_state *_state);
void minlpsetbci(minlpstate *state, ae_int_t i, double bndl, double bndu, ae_state *_state);
void minlpsetlc(minlpstate *state, RMatrix *a, ZVector *ct, ae_int_t k, ae_state *_state);
void minlpsetlc2dense(minlpstate *state, RMatrix *a, RVector *al, RVector *au, ae_int_t k, ae_state *_state);
void minlpsetlc2(minlpstate *state, sparsematrix *a, RVector *al, RVector *au, ae_int_t k, ae_state *_state);
void minlpaddlc2dense(minlpstate *state, RVector *a, double al, double au, ae_state *_state);
void minlpaddlc2(minlpstate *state, ZVector *idxa, RVector *vala, ae_int_t nnz, double al, double au, ae_state *_state);
void minlpoptimize(minlpstate *state, ae_state *_state);
void minlpresults(minlpstate *state, RVector *x, minlpreport *rep, ae_state *_state);
void minlpresultsbuf(minlpstate *state, RVector *x, minlpreport *rep, ae_state *_state);
void minlpstate_init(void *_p, ae_state *_state, bool make_automatic);
void minlpstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlpstate_free(void *_p, bool make_automatic);
void minlpreport_init(void *_p, ae_state *_state, bool make_automatic);
void minlpreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minlpreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minlpstate, );
DecClass(minlpreport, double &f; real_1d_array lagbc; real_1d_array laglc; real_1d_array y; integer_1d_array stats; double &primalerror; double &dualerror; double &slackerror; ae_int_t &iterationscount; ae_int_t &terminationtype;);

void minlpcreate(const ae_int_t n, minlpstate &state, const xparams _xparams = xdefault);
void minlpsetalgodss(const minlpstate &state, const double eps, const xparams _xparams = xdefault);
void minlpsetalgoipm(const minlpstate &state, const double eps, const xparams _xparams = xdefault);
void minlpsetalgoipm(const minlpstate &state, const xparams _xparams = xdefault);
void minlpsetcost(const minlpstate &state, const real_1d_array &c, const xparams _xparams = xdefault);
void minlpsetscale(const minlpstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minlpsetbc(const minlpstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minlpsetbcall(const minlpstate &state, const double bndl, const double bndu, const xparams _xparams = xdefault);
void minlpsetbci(const minlpstate &state, const ae_int_t i, const double bndl, const double bndu, const xparams _xparams = xdefault);
void minlpsetlc(const minlpstate &state, const real_2d_array &a, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minlpsetlc(const minlpstate &state, const real_2d_array &a, const integer_1d_array &ct, const xparams _xparams = xdefault);
void minlpsetlc2dense(const minlpstate &state, const real_2d_array &a, const real_1d_array &al, const real_1d_array &au, const ae_int_t k, const xparams _xparams = xdefault);
void minlpsetlc2dense(const minlpstate &state, const real_2d_array &a, const real_1d_array &al, const real_1d_array &au, const xparams _xparams = xdefault);
void minlpsetlc2(const minlpstate &state, const sparsematrix &a, const real_1d_array &al, const real_1d_array &au, const ae_int_t k, const xparams _xparams = xdefault);
void minlpaddlc2dense(const minlpstate &state, const real_1d_array &a, const double al, const double au, const xparams _xparams = xdefault);
void minlpaddlc2(const minlpstate &state, const integer_1d_array &idxa, const real_1d_array &vala, const ae_int_t nnz, const double al, const double au, const xparams _xparams = xdefault);
void minlpoptimize(const minlpstate &state, const xparams _xparams = xdefault);
void minlpresults(const minlpstate &state, real_1d_array &x, minlpreport &rep, const xparams _xparams = xdefault);
void minlpresultsbuf(const minlpstate &state, real_1d_array &x, minlpreport &rep, const xparams _xparams = xdefault);
} // end of namespace alglib

// === NLCSLP Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: OPTSERV, REVISEDDUALSIMPLEX
namespace alglib_impl {
typedef struct {
   presolveinfo presolver;
   dualsimplexstate dss;
   dualsimplexsettings dsssettings;
   dualsimplexbasis lastbasis;
   bool basispresent;
   ae_matrix curd;
   ae_int_t curdcnt;
   ae_vector curb;
   ae_vector curbndl;
   ae_vector curbndu;
   ae_vector cural;
   ae_vector curau;
   sparsematrix sparserawlc;
   sparsematrix sparseefflc;
   ae_int_t hessiantype;
   ae_matrix h;
   ae_matrix curhd;
   ae_matrix densedummy;
   sparsematrix sparsedummy;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_vector sk;
   ae_vector yk;
   ae_vector xs;
   ae_vector laglc;
   ae_vector lagbc;
   ae_vector cs;
} minslpsubsolver;
typedef struct {
   ae_vector sclagtmp0;
   ae_vector sclagtmp1;
} minslptmplagrangian;
typedef struct {
   ae_vector mftmp0;
} minslptmpmerit;
typedef struct {
   bool usecorrection;
   ae_vector d;
   ae_vector dx;
   ae_vector stepkxc;
   ae_vector stepkxn;
   ae_vector stepkfic;
   ae_vector stepkfin;
   ae_matrix stepkjc;
   ae_matrix stepkjn;
   ae_vector dummylagmult;
   minslptmpmerit tmpmerit;
   rcommstate rphase13state;
} minslpphase13state;
typedef struct {
   ae_vector stepkxn;
   ae_vector stepkxc;
   ae_vector stepkfin;
   ae_vector stepkfic;
   ae_matrix stepkjn;
   ae_matrix stepkjc;
   ae_vector stepklaggrad;
   ae_vector stepknlaggrad;
   ae_vector stepknlagmult;
   ae_vector meritlagmult;
   minslptmplagrangian tmplagrangianfg;
   double lastlcerr;
   ae_int_t lastlcidx;
   double lastnlcerr;
   ae_int_t lastnlcidx;
   ae_vector tmp0;
   ae_vector d;
   linminstate mcstate;
   minslptmpmerit tmpmerit;
   rcommstate rphase2state;
} minslpphase2state;
typedef struct {
   ae_int_t n;
   ae_int_t nec;
   ae_int_t nic;
   ae_int_t nlec;
   ae_int_t nlic;
   ae_vector s;
   ae_matrix scaledcleic;
   ae_vector lcsrcidx;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_vector scaledbndl;
   ae_vector scaledbndu;
   double epsx;
   ae_int_t maxits;
   ae_int_t hessiantype;
   ae_vector x;
   ae_vector fi;
   ae_matrix j;
   double f;
   bool needfij;
   bool xupdated;
   minslpphase13state state13;
   minslpphase2state state2;
   double trustrad;
   ae_int_t lpfailurecnt;
   ae_int_t fstagnationcnt;
   ae_vector step0x;
   ae_vector stepkx;
   ae_vector backupx;
   ae_vector step0fi;
   ae_vector stepkfi;
   ae_vector backupfi;
   ae_matrix step0j;
   ae_matrix stepkj;
   ae_matrix backupj;
   ae_vector meritlagmult;
   ae_vector dummylagmult;
   ae_vector fscales;
   ae_vector meritfunctionhistory;
   ae_vector maxlaghistory;
   ae_int_t historylen;
   minslpsubsolver subsolver;
   minslptmpmerit tmpmerit;
   ae_int_t repsimplexiterations;
   ae_int_t repsimplexiterations1;
   ae_int_t repsimplexiterations2;
   ae_int_t repsimplexiterations3;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repterminationtype;
   double repbcerr;
   ae_int_t repbcidx;
   double replcerr;
   ae_int_t replcidx;
   double repnlcerr;
   ae_int_t repnlcidx;
   rcommstate rstate;
} minslpstate;

void minslpinitbuf(RVector *bndl, RVector *bndu, RVector *s, RVector *x0, ae_int_t n, RMatrix *cleic, ZVector *lcsrcidx, ae_int_t nec, ae_int_t nic, ae_int_t nlec, ae_int_t nlic, double epsx, ae_int_t maxits, minslpstate *state, ae_state *_state);
bool minslpiteration(minslpstate *state, smoothnessmonitor *smonitor, bool userterminationneeded, ae_state *_state);
void minslpsubsolver_init(void *_p, ae_state *_state, bool make_automatic);
void minslpsubsolver_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minslpsubsolver_free(void *_p, bool make_automatic);
void minslptmplagrangian_init(void *_p, ae_state *_state, bool make_automatic);
void minslptmplagrangian_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minslptmplagrangian_free(void *_p, bool make_automatic);
void minslptmpmerit_init(void *_p, ae_state *_state, bool make_automatic);
void minslptmpmerit_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minslptmpmerit_free(void *_p, bool make_automatic);
void minslpphase13state_init(void *_p, ae_state *_state, bool make_automatic);
void minslpphase13state_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minslpphase13state_free(void *_p, bool make_automatic);
void minslpphase2state_init(void *_p, ae_state *_state, bool make_automatic);
void minslpphase2state_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minslpphase2state_free(void *_p, bool make_automatic);
void minslpstate_init(void *_p, ae_state *_state, bool make_automatic);
void minslpstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minslpstate_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

// === MINNLC Package ===
// Depends on: MINBLEIC, NLCSQP, NLCSLP
namespace alglib_impl {
typedef struct {
   double stabilizingpoint;
   double initialinequalitymultiplier;
   ae_int_t solvertype;
   ae_int_t prectype;
   ae_int_t updatefreq;
   double rho;
   ae_int_t n;
   double epsx;
   ae_int_t maxits;
   ae_int_t aulitscnt;
   bool xrep;
   double stpmax;
   double diffstep;
   double teststep;
   ae_vector s;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_int_t nec;
   ae_int_t nic;
   ae_matrix cleic;
   ae_vector lcsrcidx;
   ae_int_t ng;
   ae_int_t nh;
   ae_vector x;
   double f;
   ae_vector fi;
   ae_matrix j;
   bool needfij;
   bool needfi;
   bool xupdated;
   rcommstate rstate;
   rcommstate rstateaul;
   rcommstate rstateslp;
   ae_vector scaledbndl;
   ae_vector scaledbndu;
   ae_matrix scaledcleic;
   ae_vector xc;
   ae_vector xstart;
   ae_vector xbase;
   ae_vector fbase;
   ae_vector dfbase;
   ae_vector fm2;
   ae_vector fm1;
   ae_vector fp1;
   ae_vector fp2;
   ae_vector dfm1;
   ae_vector dfp1;
   ae_vector bufd;
   ae_vector bufc;
   ae_vector tmp0;
   ae_matrix bufw;
   ae_matrix bufz;
   ae_vector xk;
   ae_vector xk1;
   ae_vector gk;
   ae_vector gk1;
   double gammak;
   bool xkpresent;
   minlbfgsstate auloptimizer;
   minlbfgsreport aulreport;
   ae_vector nubc;
   ae_vector nulc;
   ae_vector nunlc;
   bool userterminationneeded;
   minslpstate slpsolverstate;
   minsqpstate sqpsolverstate;
   ae_int_t smoothnessguardlevel;
   smoothnessmonitor smonitor;
   ae_vector lastscaleused;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repnfev;
   ae_int_t repterminationtype;
   double repbcerr;
   ae_int_t repbcidx;
   double replcerr;
   ae_int_t replcidx;
   double repnlcerr;
   ae_int_t repnlcidx;
   ae_int_t repdbgphase0its;
} minnlcstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
   double bcerr;
   ae_int_t bcidx;
   double lcerr;
   ae_int_t lcidx;
   double nlcerr;
   ae_int_t nlcidx;
   ae_int_t dbgphase0its;
} minnlcreport;

void minnlccreate(ae_int_t n, RVector *x, minnlcstate *state, ae_state *_state);
void minnlccreatef(ae_int_t n, RVector *x, double diffstep, minnlcstate *state, ae_state *_state);
void minnlcsetbc(minnlcstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minnlcsetlc(minnlcstate *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void minnlcsetnlc(minnlcstate *state, ae_int_t nlec, ae_int_t nlic, ae_state *_state);
void minnlcsetcond(minnlcstate *state, double epsx, ae_int_t maxits, ae_state *_state);
void minnlcsetscale(minnlcstate *state, RVector *s, ae_state *_state);
void minnlcsetprecinexact(minnlcstate *state, ae_state *_state);
void minnlcsetprecexactlowrank(minnlcstate *state, ae_int_t updatefreq, ae_state *_state);
void minnlcsetprecexactrobust(minnlcstate *state, ae_int_t updatefreq, ae_state *_state);
void minnlcsetprecnone(minnlcstate *state, ae_state *_state);
void minnlcsetstpmax(minnlcstate *state, double stpmax, ae_state *_state);
void minnlcsetalgoaul(minnlcstate *state, double rho, ae_int_t itscnt, ae_state *_state);
void minnlcsetalgoslp(minnlcstate *state, ae_state *_state);
void minnlcsetalgosqp(minnlcstate *state, ae_state *_state);
void minnlcsetxrep(minnlcstate *state, bool needxrep, ae_state *_state);
bool minnlciteration(minnlcstate *state, ae_state *_state);
void minnlcoptguardgradient(minnlcstate *state, double teststep, ae_state *_state);
void minnlcoptguardsmoothness(minnlcstate *state, ae_int_t level, ae_state *_state);
void minnlcoptguardresults(minnlcstate *state, optguardreport *rep, ae_state *_state);
void minnlcoptguardnonc1test0results(minnlcstate *state, optguardnonc1test0report *strrep, optguardnonc1test0report *lngrep, ae_state *_state);
void minnlcoptguardnonc1test1results(minnlcstate *state, optguardnonc1test1report *strrep, optguardnonc1test1report *lngrep, ae_state *_state);
void minnlcresults(minnlcstate *state, RVector *x, minnlcreport *rep, ae_state *_state);
void minnlcresultsbuf(minnlcstate *state, RVector *x, minnlcreport *rep, ae_state *_state);
void minnlcrequesttermination(minnlcstate *state, ae_state *_state);
void minnlcrestartfrom(minnlcstate *state, RVector *x, ae_state *_state);
void minnlcequalitypenaltyfunction(double alpha, double *f, double *df, double *d2f, ae_state *_state);
void minnlcinequalitypenaltyfunction(double alpha, double stabilizingpoint, double *f, double *df, double *d2f, ae_state *_state);
void minnlcinequalityshiftfunction(double alpha, double *f, double *df, double *d2f, ae_state *_state);
void minnlcstate_init(void *_p, ae_state *_state, bool make_automatic);
void minnlcstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minnlcstate_free(void *_p, bool make_automatic);
void minnlcreport_init(void *_p, ae_state *_state, bool make_automatic);
void minnlcreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minnlcreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minnlcstate, bool &needfi; bool &needfij; bool &xupdated; double &f; real_1d_array fi; real_2d_array j; real_1d_array x;);
DecClass(minnlcreport, ae_int_t &iterationscount; ae_int_t &nfev; ae_int_t &terminationtype; double &bcerr; ae_int_t &bcidx; double &lcerr; ae_int_t &lcidx; double &nlcerr; ae_int_t &nlcidx; ae_int_t &dbgphase0its;);

void minnlccreate(const ae_int_t n, const real_1d_array &x, minnlcstate &state, const xparams _xparams = xdefault);
void minnlccreate(const real_1d_array &x, minnlcstate &state, const xparams _xparams = xdefault);
void minnlccreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minnlcstate &state, const xparams _xparams = xdefault);
void minnlccreatef(const real_1d_array &x, const double diffstep, minnlcstate &state, const xparams _xparams = xdefault);
void minnlcsetbc(const minnlcstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minnlcsetlc(const minnlcstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minnlcsetlc(const minnlcstate &state, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = xdefault);
void minnlcsetnlc(const minnlcstate &state, const ae_int_t nlec, const ae_int_t nlic, const xparams _xparams = xdefault);
void minnlcsetcond(const minnlcstate &state, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minnlcsetscale(const minnlcstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minnlcsetprecinexact(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcsetprecexactlowrank(const minnlcstate &state, const ae_int_t updatefreq, const xparams _xparams = xdefault);
void minnlcsetprecexactrobust(const minnlcstate &state, const ae_int_t updatefreq, const xparams _xparams = xdefault);
void minnlcsetprecnone(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcsetstpmax(const minnlcstate &state, const double stpmax, const xparams _xparams = xdefault);
void minnlcsetalgoaul(const minnlcstate &state, const double rho, const ae_int_t itscnt, const xparams _xparams = xdefault);
void minnlcsetalgoslp(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcsetalgosqp(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcsetxrep(const minnlcstate &state, const bool needxrep, const xparams _xparams = xdefault);
bool minnlciteration(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcoptimize(minnlcstate &state, void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minnlcoptimize(minnlcstate &state, void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minnlcoptguardgradient(const minnlcstate &state, const double teststep, const xparams _xparams = xdefault);
void minnlcoptguardsmoothness(const minnlcstate &state, const ae_int_t level, const xparams _xparams = xdefault);
void minnlcoptguardsmoothness(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcoptguardresults(const minnlcstate &state, optguardreport &rep, const xparams _xparams = xdefault);
void minnlcoptguardnonc1test0results(const minnlcstate &state, optguardnonc1test0report &strrep, optguardnonc1test0report &lngrep, const xparams _xparams = xdefault);
void minnlcoptguardnonc1test1results(const minnlcstate &state, optguardnonc1test1report &strrep, optguardnonc1test1report &lngrep, const xparams _xparams = xdefault);
void minnlcresults(const minnlcstate &state, real_1d_array &x, minnlcreport &rep, const xparams _xparams = xdefault);
void minnlcresultsbuf(const minnlcstate &state, real_1d_array &x, minnlcreport &rep, const xparams _xparams = xdefault);
void minnlcrequesttermination(const minnlcstate &state, const xparams _xparams = xdefault);
void minnlcrestartfrom(const minnlcstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MINNS Package ===
// Depends on: MINBLEIC
namespace alglib_impl {
typedef struct {
   double fc;
   double fn;
   ae_vector xc;
   ae_vector xn;
   ae_vector x0;
   ae_vector gc;
   ae_vector d;
   ae_matrix uh;
   ae_matrix ch;
   ae_matrix rk;
   ae_vector invutc;
   ae_vector tmp0;
   ae_vector tmpidx;
   ae_vector tmpd;
   ae_vector tmpc;
   ae_vector tmplambdas;
   ae_matrix tmpc2;
   ae_vector tmpb;
   snnlssolver nnls;
} minnsqp;
typedef struct {
   ae_int_t solvertype;
   ae_int_t n;
   double epsx;
   ae_int_t maxits;
   bool xrep;
   double diffstep;
   ae_vector s;
   ae_vector bndl;
   ae_vector bndu;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_int_t nec;
   ae_int_t nic;
   ae_matrix cleic;
   ae_int_t ng;
   ae_int_t nh;
   ae_vector x;
   double f;
   ae_vector fi;
   ae_matrix j;
   bool needfij;
   bool needfi;
   bool xupdated;
   rcommstate rstate;
   rcommstate rstateags;
   hqrndstate agsrs;
   double agsradius;
   ae_int_t agssamplesize;
   double agsraddecay;
   double agsalphadecay;
   double agsdecrease;
   double agsinitstp;
   double agsstattold;
   double agsshortstpabs;
   double agsshortstprel;
   double agsshortf;
   ae_int_t agsshortlimit;
   double agsrhononlinear;
   ae_int_t agsminupdate;
   ae_int_t agsmaxraddecays;
   ae_int_t agsmaxbacktrack;
   ae_int_t agsmaxbacktracknonfull;
   double agspenaltylevel;
   double agspenaltyincrease;
   ae_vector xstart;
   ae_vector xc;
   ae_vector xn;
   ae_vector rawg;
   ae_vector meritg;
   double rawf;
   double meritf;
   ae_vector d;
   ae_vector colmax;
   ae_vector diagh;
   ae_vector signmin;
   ae_vector signmax;
   bool userterminationneeded;
   ae_vector scaledbndl;
   ae_vector scaledbndu;
   ae_matrix scaledcleic;
   double rholinear;
   ae_matrix samplex;
   ae_matrix samplegm;
   ae_matrix samplegmbc;
   ae_vector samplef;
   minnsqp nsqp;
   ae_vector tmp0;
   ae_vector tmp1;
   ae_matrix tmp2;
   ae_vector tmp3;
   ae_vector xbase;
   ae_vector fbase;
   ae_vector fp;
   ae_vector fm;
   ae_vector xscaled;
   ae_int_t repinneriterationscount;
   ae_int_t repouteriterationscount;
   ae_int_t repnfev;
   ae_int_t repvaridx;
   ae_int_t repfuncidx;
   ae_int_t repterminationtype;
   double replcerr;
   double repnlcerr;
   ae_int_t dbgncholesky;
} minnsstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   double cerr;
   double lcerr;
   double nlcerr;
   ae_int_t terminationtype;
   ae_int_t varidx;
   ae_int_t funcidx;
} minnsreport;

void minnscreate(ae_int_t n, RVector *x, minnsstate *state, ae_state *_state);
void minnscreatef(ae_int_t n, RVector *x, double diffstep, minnsstate *state, ae_state *_state);
void minnssetbc(minnsstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minnssetlc(minnsstate *state, RMatrix *c, ZVector *ct, ae_int_t k, ae_state *_state);
void minnssetnlc(minnsstate *state, ae_int_t nlec, ae_int_t nlic, ae_state *_state);
void minnssetcond(minnsstate *state, double epsx, ae_int_t maxits, ae_state *_state);
void minnssetscale(minnsstate *state, RVector *s, ae_state *_state);
void minnssetalgoags(minnsstate *state, double radius, double penalty, ae_state *_state);
void minnssetxrep(minnsstate *state, bool needxrep, ae_state *_state);
void minnsrequesttermination(minnsstate *state, ae_state *_state);
bool minnsiteration(minnsstate *state, ae_state *_state);
void minnsresults(minnsstate *state, RVector *x, minnsreport *rep, ae_state *_state);
void minnsresultsbuf(minnsstate *state, RVector *x, minnsreport *rep, ae_state *_state);
void minnsrestartfrom(minnsstate *state, RVector *x, ae_state *_state);
void minnsqp_init(void *_p, ae_state *_state, bool make_automatic);
void minnsqp_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minnsqp_free(void *_p, bool make_automatic);
void minnsstate_init(void *_p, ae_state *_state, bool make_automatic);
void minnsstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minnsstate_free(void *_p, bool make_automatic);
void minnsreport_init(void *_p, ae_state *_state, bool make_automatic);
void minnsreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minnsreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minnsstate, bool &needfi; bool &needfij; bool &xupdated; double &f; real_1d_array fi; real_2d_array j; real_1d_array x;);
DecClass(minnsreport, ae_int_t &iterationscount; ae_int_t &nfev; double &cerr; double &lcerr; double &nlcerr; ae_int_t &terminationtype; ae_int_t &varidx; ae_int_t &funcidx;);

void minnscreate(const ae_int_t n, const real_1d_array &x, minnsstate &state, const xparams _xparams = xdefault);
void minnscreate(const real_1d_array &x, minnsstate &state, const xparams _xparams = xdefault);
void minnscreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minnsstate &state, const xparams _xparams = xdefault);
void minnscreatef(const real_1d_array &x, const double diffstep, minnsstate &state, const xparams _xparams = xdefault);
void minnssetbc(const minnsstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minnssetlc(const minnsstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k, const xparams _xparams = xdefault);
void minnssetlc(const minnsstate &state, const real_2d_array &c, const integer_1d_array &ct, const xparams _xparams = xdefault);
void minnssetnlc(const minnsstate &state, const ae_int_t nlec, const ae_int_t nlic, const xparams _xparams = xdefault);
void minnssetcond(const minnsstate &state, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minnssetscale(const minnsstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minnssetalgoags(const minnsstate &state, const double radius, const double penalty, const xparams _xparams = xdefault);
void minnssetxrep(const minnsstate &state, const bool needxrep, const xparams _xparams = xdefault);
void minnsrequesttermination(const minnsstate &state, const xparams _xparams = xdefault);
bool minnsiteration(const minnsstate &state, const xparams _xparams = xdefault);
void minnsoptimize(minnsstate &state, void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minnsoptimize(minnsstate &state, void (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minnsresults(const minnsstate &state, real_1d_array &x, minnsreport &rep, const xparams _xparams = xdefault);
void minnsresultsbuf(const minnsstate &state, real_1d_array &x, minnsreport &rep, const xparams _xparams = xdefault);
void minnsrestartfrom(const minnsstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MINCOMP Package ===
// Depends on: MINLBFGS, MINBLEIC
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxits;
   bool xrep;
   double stpmax;
   ae_int_t cgtype;
   ae_int_t k;
   ae_int_t nfev;
   ae_int_t mcstage;
   ae_vector bndl;
   ae_vector bndu;
   ae_int_t curalgo;
   ae_int_t acount;
   double mu;
   double finit;
   double dginit;
   ae_vector ak;
   ae_vector xk;
   ae_vector dk;
   ae_vector an;
   ae_vector xn;
   ae_vector dn;
   ae_vector d;
   double fold;
   double stp;
   ae_vector work;
   ae_vector yk;
   ae_vector gc;
   double laststep;
   ae_vector x;
   double f;
   ae_vector g;
   bool needfg;
   bool xupdated;
   rcommstate rstate;
   ae_int_t repiterationscount;
   ae_int_t repnfev;
   ae_int_t repterminationtype;
   ae_int_t debugrestartscount;
   linminstate lstate;
   double betahs;
   double betady;
} minasastate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   ae_int_t terminationtype;
   ae_int_t activeconstraints;
} minasareport;

void minlbfgssetdefaultpreconditioner(minlbfgsstate *state, ae_state *_state);
void minlbfgssetcholeskypreconditioner(minlbfgsstate *state, RMatrix *p, bool isupper, ae_state *_state);
void minbleicsetbarrierwidth(minbleicstate *state, double mu, ae_state *_state);
void minbleicsetbarrierdecay(minbleicstate *state, double mudecay, ae_state *_state);
void minasacreate(ae_int_t n, RVector *x, RVector *bndl, RVector *bndu, minasastate *state, ae_state *_state);
void minasasetcond(minasastate *state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
void minasasetxrep(minasastate *state, bool needxrep, ae_state *_state);
void minasasetalgorithm(minasastate *state, ae_int_t algotype, ae_state *_state);
void minasasetstpmax(minasastate *state, double stpmax, ae_state *_state);
bool minasaiteration(minasastate *state, ae_state *_state);
void minasaresults(minasastate *state, RVector *x, minasareport *rep, ae_state *_state);
void minasaresultsbuf(minasastate *state, RVector *x, minasareport *rep, ae_state *_state);
void minasarestartfrom(minasastate *state, RVector *x, RVector *bndl, RVector *bndu, ae_state *_state);
void minasastate_init(void *_p, ae_state *_state, bool make_automatic);
void minasastate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minasastate_free(void *_p, bool make_automatic);
void minasareport_init(void *_p, ae_state *_state, bool make_automatic);
void minasareport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minasareport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minasastate, bool &needfg; bool &xupdated; double &f; real_1d_array g; real_1d_array x;);
DecClass(minasareport, ae_int_t &iterationscount; ae_int_t &nfev; ae_int_t &terminationtype; ae_int_t &activeconstraints;);

void minlbfgssetdefaultpreconditioner(const minlbfgsstate &state, const xparams _xparams = xdefault);
void minlbfgssetcholeskypreconditioner(const minlbfgsstate &state, const real_2d_array &p, const bool isupper, const xparams _xparams = xdefault);
void minbleicsetbarrierwidth(const minbleicstate &state, const double mu, const xparams _xparams = xdefault);
void minbleicsetbarrierdecay(const minbleicstate &state, const double mudecay, const xparams _xparams = xdefault);
void minasacreate(const ae_int_t n, const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, minasastate &state, const xparams _xparams = xdefault);
void minasacreate(const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, minasastate &state, const xparams _xparams = xdefault);
void minasasetcond(const minasastate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minasasetxrep(const minasastate &state, const bool needxrep, const xparams _xparams = xdefault);
void minasasetalgorithm(const minasastate &state, const ae_int_t algotype, const xparams _xparams = xdefault);
void minasasetstpmax(const minasastate &state, const double stpmax, const xparams _xparams = xdefault);
bool minasaiteration(const minasastate &state, const xparams _xparams = xdefault);
void minasaoptimize(minasastate &state, void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minasaresults(const minasastate &state, real_1d_array &x, minasareport &rep, const xparams _xparams = xdefault);
void minasaresultsbuf(const minasastate &state, real_1d_array &x, minasareport &rep, const xparams _xparams = xdefault);
void minasarestartfrom(const minasastate &state, const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
} // end of namespace alglib

// === MINBC Package ===
// Depends on: (AlgLibInternal) LINMIN
// Depends on: OPTSERV
namespace alglib_impl {
typedef struct {
   ae_int_t nmain;
   double epsg;
   double epsf;
   double epsx;
   ae_int_t maxits;
   bool xrep;
   double stpmax;
   double diffstep;
   ae_vector s;
   ae_int_t prectype;
   ae_vector diagh;
   ae_vector x;
   double f;
   ae_vector g;
   bool needf;
   bool needfg;
   bool xupdated;
   bool userterminationneeded;
   rcommstate rstate;
   ae_vector xc;
   ae_vector ugc;
   ae_vector cgc;
   ae_vector xn;
   ae_vector ugn;
   ae_vector cgn;
   ae_vector xp;
   double fc;
   double fn;
   double fp;
   ae_vector d;
   double lastscaledgoodstep;
   ae_vector hasbndl;
   ae_vector hasbndu;
   ae_vector bndl;
   ae_vector bndu;
   ae_int_t repiterationscount;
   ae_int_t repnfev;
   ae_int_t repvaridx;
   ae_int_t repterminationtype;
   ae_vector xstart;
   double fbase;
   double fm2;
   double fm1;
   double fp1;
   double fp2;
   double xm1;
   double xp1;
   double gm1;
   double gp1;
   ae_vector tmpprec;
   ae_vector tmp0;
   ae_int_t nfev;
   ae_int_t mcstage;
   double stp;
   double curstpmax;
   ae_vector work;
   linminstate lstate;
   double trimthreshold;
   ae_int_t nonmonotoniccnt;
   ae_matrix bufyk;
   ae_matrix bufsk;
   ae_vector bufrho;
   ae_vector buftheta;
   ae_int_t bufsize;
   double teststep;
   ae_int_t smoothnessguardlevel;
   smoothnessmonitor smonitor;
   ae_vector lastscaleused;
   ae_vector invs;
} minbcstate;
typedef struct {
   ae_int_t iterationscount;
   ae_int_t nfev;
   ae_int_t varidx;
   ae_int_t terminationtype;
} minbcreport;

void minbccreate(ae_int_t n, RVector *x, minbcstate *state, ae_state *_state);
void minbccreatef(ae_int_t n, RVector *x, double diffstep, minbcstate *state, ae_state *_state);
void minbcsetbc(minbcstate *state, RVector *bndl, RVector *bndu, ae_state *_state);
void minbcsetcond(minbcstate *state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
void minbcsetscale(minbcstate *state, RVector *s, ae_state *_state);
void minbcsetprecdefault(minbcstate *state, ae_state *_state);
void minbcsetprecdiag(minbcstate *state, RVector *d, ae_state *_state);
void minbcsetprecscale(minbcstate *state, ae_state *_state);
void minbcsetxrep(minbcstate *state, bool needxrep, ae_state *_state);
void minbcsetstpmax(minbcstate *state, double stpmax, ae_state *_state);
bool minbciteration(minbcstate *state, ae_state *_state);
void minbcoptguardgradient(minbcstate *state, double teststep, ae_state *_state);
void minbcoptguardsmoothness(minbcstate *state, ae_int_t level, ae_state *_state);
void minbcoptguardresults(minbcstate *state, optguardreport *rep, ae_state *_state);
void minbcoptguardnonc1test0results(minbcstate *state, optguardnonc1test0report *strrep, optguardnonc1test0report *lngrep, ae_state *_state);
void minbcoptguardnonc1test1results(minbcstate *state, optguardnonc1test1report *strrep, optguardnonc1test1report *lngrep, ae_state *_state);
void minbcresults(minbcstate *state, RVector *x, minbcreport *rep, ae_state *_state);
void minbcresultsbuf(minbcstate *state, RVector *x, minbcreport *rep, ae_state *_state);
void minbcrestartfrom(minbcstate *state, RVector *x, ae_state *_state);
void minbcrequesttermination(minbcstate *state, ae_state *_state);
void minbcstate_init(void *_p, ae_state *_state, bool make_automatic);
void minbcstate_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minbcstate_free(void *_p, bool make_automatic);
void minbcreport_init(void *_p, ae_state *_state, bool make_automatic);
void minbcreport_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void minbcreport_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(minbcstate, bool &needf; bool &needfg; bool &xupdated; double &f; real_1d_array g; real_1d_array x;);
DecClass(minbcreport, ae_int_t &iterationscount; ae_int_t &nfev; ae_int_t &varidx; ae_int_t &terminationtype;);

void minbccreate(const ae_int_t n, const real_1d_array &x, minbcstate &state, const xparams _xparams = xdefault);
void minbccreate(const real_1d_array &x, minbcstate &state, const xparams _xparams = xdefault);
void minbccreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minbcstate &state, const xparams _xparams = xdefault);
void minbccreatef(const real_1d_array &x, const double diffstep, minbcstate &state, const xparams _xparams = xdefault);
void minbcsetbc(const minbcstate &state, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void minbcsetcond(const minbcstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits, const xparams _xparams = xdefault);
void minbcsetscale(const minbcstate &state, const real_1d_array &s, const xparams _xparams = xdefault);
void minbcsetprecdefault(const minbcstate &state, const xparams _xparams = xdefault);
void minbcsetprecdiag(const minbcstate &state, const real_1d_array &d, const xparams _xparams = xdefault);
void minbcsetprecscale(const minbcstate &state, const xparams _xparams = xdefault);
void minbcsetxrep(const minbcstate &state, const bool needxrep, const xparams _xparams = xdefault);
void minbcsetstpmax(const minbcstate &state, const double stpmax, const xparams _xparams = xdefault);
bool minbciteration(const minbcstate &state, const xparams _xparams = xdefault);
void minbcoptimize(minbcstate &state, void (*func)(const real_1d_array &x, double &func, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minbcoptimize(minbcstate &state, void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr), void (*rep)(const real_1d_array &x, double func, void *ptr) = NULL, void *ptr = NULL, const xparams _xparams = xdefault);
void minbcoptguardgradient(const minbcstate &state, const double teststep, const xparams _xparams = xdefault);
void minbcoptguardsmoothness(const minbcstate &state, const ae_int_t level, const xparams _xparams = xdefault);
void minbcoptguardsmoothness(const minbcstate &state, const xparams _xparams = xdefault);
void minbcoptguardresults(const minbcstate &state, optguardreport &rep, const xparams _xparams = xdefault);
void minbcoptguardnonc1test0results(const minbcstate &state, optguardnonc1test0report &strrep, optguardnonc1test0report &lngrep, const xparams _xparams = xdefault);
void minbcoptguardnonc1test1results(const minbcstate &state, optguardnonc1test1report &strrep, optguardnonc1test1report &lngrep, const xparams _xparams = xdefault);
void minbcresults(const minbcstate &state, real_1d_array &x, minbcreport &rep, const xparams _xparams = xdefault);
void minbcresultsbuf(const minbcstate &state, real_1d_array &x, minbcreport &rep, const xparams _xparams = xdefault);
void minbcrestartfrom(const minbcstate &state, const real_1d_array &x, const xparams _xparams = xdefault);
void minbcrequesttermination(const minbcstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

// === OPTS Package ===
// Depends on: MINLP
namespace alglib_impl {
typedef struct {
   ae_int_t n;
   bool hasknowntarget;
   double targetf;
   ae_vector s;
   ae_vector c;
   ae_vector bndl;
   ae_vector bndu;
   ae_int_t m;
   sparsematrix a;
   ae_vector al;
   ae_vector au;
} lptestproblem;

void lptestproblemcreate(ae_int_t n, bool hasknowntarget, double targetf, lptestproblem *p, ae_state *_state);
void lptestproblemsetscale(lptestproblem *p, RVector *s, ae_state *_state);
void lptestproblemsetcost(lptestproblem *p, RVector *c, ae_state *_state);
void lptestproblemsetbc(lptestproblem *p, RVector *bndl, RVector *bndu, ae_state *_state);
void lptestproblemsetlc2(lptestproblem *p, sparsematrix *a, RVector *al, RVector *au, ae_int_t m, ae_state *_state);
void lptestproblemalloc(ae_serializer *s, lptestproblem *p, ae_state *_state);
void lptestproblemserialize(ae_serializer *s, lptestproblem *p, ae_state *_state);
void lptestproblemunserialize(ae_serializer *s, lptestproblem *p, ae_state *_state);
void xdbgminlpcreatefromtestproblem(lptestproblem *p, minlpstate *state, ae_state *_state);
void lptestproblem_init(void *_p, ae_state *_state, bool make_automatic);
void lptestproblem_copy(void *_dst, void *_src, ae_state *_state, bool make_automatic);
void lptestproblem_free(void *_p, bool make_automatic);
} // end of namespace alglib_impl

namespace alglib {
DecClass(lptestproblem, );

void lptestproblemserialize(lptestproblem &obj, std::string &s_out);
void lptestproblemunserialize(const std::string &s_in, lptestproblem &obj);
void lptestproblemserialize(lptestproblem &obj, std::ostream &s_out);
void lptestproblemunserialize(const std::istream &s_in, lptestproblem &obj);
void lptestproblemcreate(const ae_int_t n, const bool hasknowntarget, const double targetf, lptestproblem &p, const xparams _xparams = xdefault);
void lptestproblemsetscale(const lptestproblem &p, const real_1d_array &s, const xparams _xparams = xdefault);
void lptestproblemsetcost(const lptestproblem &p, const real_1d_array &c, const xparams _xparams = xdefault);
void lptestproblemsetbc(const lptestproblem &p, const real_1d_array &bndl, const real_1d_array &bndu, const xparams _xparams = xdefault);
void lptestproblemsetlc2(const lptestproblem &p, const sparsematrix &a, const real_1d_array &al, const real_1d_array &au, const ae_int_t m, const xparams _xparams = xdefault);
void xdbgminlpcreatefromtestproblem(const lptestproblem &p, minlpstate &state, const xparams _xparams = xdefault);
} // end of namespace alglib

#endif // OnceOnly
