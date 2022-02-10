#ifndef OnceOnlyKernelsFma_h
#define OnceOnlyKernelsFma_h

#include "Ap.h"
namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_AVX2_INTRINSICS
double rdotv_fma(const ae_int_t n, const Real *__restrict x, const Real *__restrict y);
double rdotv2_fma(const ae_int_t n, const Real *__restrict x);
void raddv_fma(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x);
void raddvx_fma(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x);
void rgemv_straight_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y);
void rgemv_transposed_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y);
void rgemvx_straight_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y);
void rgemvx_transposed_fma(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y);
void ablasf_dotblkh_fma(const double *src_a, const double *src_b, ae_int_t round_length, ae_int_t block_size, ae_int_t micro_size, double *dst, ae_int_t dst_stride);
void spchol_propagatefwd_fma(RVector *x, ae_int_t cols0, ae_int_t blocksize, ZVector *superrowidx, ae_int_t rbase, ae_int_t offdiagsize, RVector *rowstorage, ae_int_t offss, ae_int_t sstride, RVector *simdbuf, ae_int_t simdwidth);
bool spchol_updatekernelabc4_fma(double *rowstorage, ae_int_t offss, ae_int_t twidth, ae_int_t offsu, ae_int_t uheight, ae_int_t urank, ae_int_t urowstride, ae_int_t uwidth, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase);
bool spchol_updatekernel4444_fma(double *rowstorage, ae_int_t offss, ae_int_t sheight, ae_int_t offsu, ae_int_t uheight, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase);
// ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_AVX2_INTRINSICS
#endif
} // end of namespace alglib_impl

#endif // OnceOnly
