#ifndef OnceOnlyKernelsSse2_h
#define OnceOnlyKernelsSse2_h

#include "Ap.h"
namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_SSE2_INTRINSICS
double rdotv_sse2(ae_int_t n, const Real *x, const Real *y);
double rdotv2_sse2(ae_int_t n, const Real *x);
void rcopyv_sse2(const ae_int_t n, const Real *__restrict x, Real *__restrict y);
void rcopymulv_sse2(const ae_int_t n, const double v, const Real *__restrict x, Real *__restrict y);
void icopyv_sse2(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y);
void bcopyv_sse2(const ae_int_t n, const bool *__restrict x, bool *__restrict y);
void rsetv_sse2(const ae_int_t n, const double v, Real *__restrict x);
void rsetvx_sse2(const ae_int_t n, const double v, Real *__restrict x);
void isetv_sse2(const ae_int_t n, const ae_int_t v, ae_int_t *__restrict x);
void bsetv_sse2(const ae_int_t n, const bool v, bool *__restrict x);
void rmulv_sse2(const ae_int_t n, const double v, double *__restrict x);
void rmulvx_sse2(const ae_int_t n, const double v, double *__restrict x);
void raddv_sse2(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x);
void raddvx_sse2(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x);
void rmergemulv_sse2(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
void rmergemaxv_sse2(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
void rmergeminv_sse2(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
double rmaxv_sse2(ae_int_t n, const Real *__restrict x);
double rmaxabsv_sse2(ae_int_t n, const Real *__restrict x);
void rcopyvx_sse2(const ae_int_t n, const double *__restrict x, double *__restrict y);
void icopyvx_sse2(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y);
// ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_SSE2_INTRINSICS
#endif
} // end of namespace alglib_impl

#endif // OnceOnly
