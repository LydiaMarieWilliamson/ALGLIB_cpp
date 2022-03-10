#ifndef OnceOnlyKernelsSse2_h
#define OnceOnlyKernelsSse2_h

#include "Ap.h"
namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_SSE2_INTRINSICS
double sse2_rdotv(const ae_int_t n, const Real *x, const Real *y);
double sse2_rdotv2(const ae_int_t n, const Real *x);
void sse2_rcopyv(const ae_int_t n, const Real *__restrict x, Real *__restrict y);
void sse2_rcopymulv(const ae_int_t n, const double v, const Real *__restrict x, Real *__restrict y);
void sse2_icopyv(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y);
void sse2_bcopyv(const ae_int_t n, const bool *__restrict x, bool *__restrict y);
void sse2_rsetv(const ae_int_t n, const double v, Real *__restrict x);
void sse2_rsetvx(const ae_int_t n, const double v, Real *__restrict x);
void sse2_isetv(const ae_int_t n, const ae_int_t v, ae_int_t *__restrict x);
void sse2_bsetv(const ae_int_t n, const bool v, bool *__restrict x);
void sse2_rmulv(const ae_int_t n, const double v, double *__restrict x);
void sse2_rmulvx(const ae_int_t n, const double v, double *__restrict x);
void sse2_raddv(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x);
void sse2_raddvx(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x);
void sse2_rmergemulv(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
void sse2_rmergemaxv(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
void sse2_rmergeminv(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
double sse2_rmaxv(ae_int_t n, const Real *__restrict x);
double sse2_rmaxabsv(ae_int_t n, const Real *__restrict x);
void sse2_rcopyvx(const ae_int_t n, const double *__restrict x, double *__restrict y);
void sse2_icopyvx(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y);
#endif // ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_SSE2_INTRINSICS
} // end of namespace alglib_impl

#endif // OnceOnly
