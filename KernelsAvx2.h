#ifndef OnceOnlyKernelsAvx2_h
#define OnceOnlyKernelsAvx2_h

#include "Ap.h"
namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_AVX2_INTRINSICS
double avx2_rdotv(const ae_int_t n, const Real *__restrict x, const Real *__restrict y);
double avx2_rdotv2(const ae_int_t n, const Real *__restrict x);
void avx2_rcopyv(const ae_int_t n, const Real *__restrict x, Real *__restrict y);
void avx2_rcopymulv(const ae_int_t n, const double v, const Real *__restrict x, Real *__restrict y);
void avx2_icopyv(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y);
void avx2_bcopyv(const ae_int_t n, const bool *__restrict x, bool *__restrict y);
void avx2_rsetv(const ae_int_t n, const double v, Real *__restrict x);
void avx2_rsetvx(const ae_int_t n, const double v, Real *__restrict x);
void avx2_isetv(const ae_int_t n, const ae_int_t v, ae_int_t *__restrict x);
void avx2_bsetv(const ae_int_t n, const bool v, bool *__restrict x);
void avx2_rmulv(const ae_int_t n, const double v, double *__restrict x);
void avx2_rmulvx(const ae_int_t n, const double v, double *__restrict x);
void avx2_raddv(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x);
void avx2_raddvx(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x);
void avx2_rmergemulv(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
void avx2_rmergemaxv(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
void avx2_rmergeminv(const ae_int_t n, const Real *__restrict y, Real *__restrict x);
double avx2_rmaxv(ae_int_t n, const Real *__restrict x);
double avx2_rmaxabsv(ae_int_t n, const Real *__restrict x);
void avx2_rcopyvx(const ae_int_t n, const double *__restrict x, double *__restrict y);
void avx2_icopyvx(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y);
void avx2_rgemv_straight(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y);
void avx2_rgemv_transposed(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y);
void avx2_rgemvx_straight(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y);
void avx2_rgemvx_transposed(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y);
ae_int_t avx2_ablasf_packblkh(const double *src, ae_int_t src_stride, ae_int_t op, ae_int_t opsrc_length, ae_int_t opsrc_width, double *dst, ae_int_t block_size, ae_int_t micro_size);
ae_int_t avx2_ablasf_packblkh32(const double *src, ae_int_t src_stride, ae_int_t op, ae_int_t ignore_opsrc_length, ae_int_t opsrc_width, double *dst, ae_int_t ignore_block_size, ae_int_t micro_size);
void avx2_ablasf_dotblkh(const double *src_a, const double *src_b, ae_int_t round_length, ae_int_t block_size, ae_int_t micro_size, double *dst, ae_int_t dst_stride);
void avx2_ablasf_daxpby(ae_int_t n, double alpha, const double *src, double beta, double *dst);
bool avx2_spchol_updatekernelabc4(double *rowstorage, ae_int_t offss, ae_int_t twidth, ae_int_t offsu, ae_int_t uheight, ae_int_t urank, ae_int_t urowstride, ae_int_t uwidth, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase);
bool avx2_spchol_updatekernel4444(double *rowstorage, ae_int_t offss, ae_int_t sheight, ae_int_t offsu, ae_int_t uheight, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase);
#endif // ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_AVX2_INTRINSICS
} // end of namespace alglib_impl

#endif // OnceOnly
