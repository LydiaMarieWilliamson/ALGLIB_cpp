#ifndef OnceOnlyKernelsFma_h
#define OnceOnlyKernelsFma_h

#include "Ap.h"
namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_FMA_INTRINSICS
double fma_rdotv(const ae_int_t n, const Real *__restrict x, const Real *__restrict y);
double fma_rdotv2(const ae_int_t n, const Real *__restrict x);
void fma_raddv(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x);
void fma_raddvx(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x);
void fma_rgemv_straight(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y);
void fma_rgemv_transposed(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y);
void fma_rgemvx_straight(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y);
void fma_rgemvx_transposed(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y);
void fma_ablasf_dotblkh(const double *src_a, const double *src_b, ae_int_t round_length, ae_int_t block_size, ae_int_t micro_size, double *dst, ae_int_t dst_stride);
void fma_spchol_propagatefwd(RVector *x, ae_int_t cols0, ae_int_t blocksize, ZVector *superrowidx, ae_int_t rbase, ae_int_t offdiagsize, RVector *rowstorage, ae_int_t offss, ae_int_t sstride, RVector *simdbuf, ae_int_t simdwidth);
bool fma_spchol_updatekernelabc4(double *rowstorage, ae_int_t offss, ae_int_t twidth, ae_int_t offsu, ae_int_t uheight, ae_int_t urank, ae_int_t urowstride, ae_int_t uwidth, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase);
bool fma_spchol_updatekernel4444(double *rowstorage, ae_int_t offss, ae_int_t sheight, ae_int_t offsu, ae_int_t uheight, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase);
#endif // ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_FMA_INTRINSICS
} // end of namespace alglib_impl

#endif // OnceOnly
