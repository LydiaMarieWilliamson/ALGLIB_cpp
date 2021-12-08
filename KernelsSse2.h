#ifndef OnceOnlyKernelsSse2_h
#define OnceOnlyKernelsSse2_h

#include "Ap.h"

namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_SSE2_INTRINSICS

double rdotv_sse2(ae_int_t n, const Real *x, const Real *y, ae_state *_state);
double rdotv2_sse2(ae_int_t n, const Real *x, ae_state *_state);
void rcopyv_sse2(const ae_int_t n, const Real *__restrict x, Real *__restrict y, ae_state *__restrict _state);
void rcopymulv_sse2(const ae_int_t n, const double v, const Real *__restrict x, Real *__restrict y, const ae_state *__restrict _state);
void icopyv_sse2(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y, ae_state *__restrict _state);
void bcopyv_sse2(const ae_int_t n, const bool *__restrict x, bool *__restrict y, ae_state *__restrict _state);
void rsetv_sse2(const ae_int_t n, const double v, Real *__restrict x, const ae_state *__restrict _state);
void rsetvx_sse2(const ae_int_t n, const double v, Real *__restrict x, const ae_state *__restrict _state);
void isetv_sse2(const ae_int_t n, const ae_int_t v, ae_int_t *__restrict x, ae_state *__restrict _state);
void bsetv_sse2(const ae_int_t n, const bool v, bool *__restrict x, ae_state *__restrict _state);
void rmulv_sse2(const ae_int_t n, const double v, double *__restrict x, const ae_state *__restrict _state);
void rmulvx_sse2(const ae_int_t n, const double v, double *__restrict x, const ae_state *__restrict _state);
void raddv_sse2(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x, const ae_state *__restrict _state);
void raddvx_sse2(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x, ae_state *_state);
void rmergemulv_sse2(const ae_int_t n, const Real *__restrict y, Real *__restrict x, const ae_state *__restrict _state);
void rmergemaxv_sse2(const ae_int_t n, const Real *__restrict y, Real *__restrict x, ae_state *__restrict _state);
void rmergeminv_sse2(const ae_int_t n, const Real *__restrict y, Real *__restrict x, ae_state *__restrict _state);
double rmaxv_sse2(ae_int_t n, const Real *__restrict x, ae_state *__restrict _state);
double rmaxabsv_sse2(ae_int_t n, const Real *__restrict x, ae_state *__restrict _state);
void rcopyvx_sse2(const ae_int_t n, const double *__restrict x, double *__restrict y, ae_state *_state);
void icopyvx_sse2(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y, ae_state *__restrict _state);
// ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_SSE2_INTRINSICS
#endif
} // end of namespace alglib_impl

#endif // OnceOnly
