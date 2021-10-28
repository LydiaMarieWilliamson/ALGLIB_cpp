#ifndef OnceOnlyKernelsFma_h
#define OnceOnlyKernelsFma_h

#include "Ap.h"

#define AE_USE_CPP

namespace alglib_impl {
#if !defined(ALGLIB_NO_FAST_KERNELS) && defined(_ALGLIB_HAS_AVX2_INTRINSICS)

double rdotv_fma(const ae_int_t n, const Real *__restrict x, const Real *__restrict y, const ae_state *__restrict _state);
double rdotv2_fma(const ae_int_t n, const Real *__restrict x, const ae_state *__restrict _state);
void raddv_fma(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x, const ae_state *__restrict _state);
void raddvx_fma(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x, ae_state *_state);
void rgemv_straight_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y, ae_state *_state);
void rgemv_transposed_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y, ae_state *_state);
void rgemvx_straight_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y, ae_state *_state);
void rgemvx_transposed_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y, ae_state *_state);

void ablasf_dotblkh_fma(const double *src_a, const double *src_b, ae_int_t round_length, ae_int_t block_size, ae_int_t micro_size, double *dst, ae_int_t dst_stride);
void spchol_propagatefwd_fma(RVector *x, ae_int_t cols0, ae_int_t blocksize, ZVector *superrowidx, ae_int_t rbase, ae_int_t offdiagsize, RVector *rowstorage, ae_int_t offss, ae_int_t sstride, RVector *simdbuf, ae_int_t simdwidth, ae_state *_state);
ae_bool spchol_updatekernelabc4_fma(double *rowstorage, ae_int_t offss, ae_int_t twidth, ae_int_t offsu, ae_int_t uheight, ae_int_t urank, ae_int_t urowstride, ae_int_t uwidth, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase, ae_state *_state);
ae_bool spchol_updatekernel4444_fma(double *rowstorage, ae_int_t offss, ae_int_t sheight, ae_int_t offsu, ae_int_t uheight, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase, ae_state *_state);

// ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_AVX2_INTRINSICS
#endif
} // end of namespace alglib_impl

#endif // OnceOnly
