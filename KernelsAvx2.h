#ifndef OnceOnlyKernelsAvx2_h
#define OnceOnlyKernelsAvx2_h

#include "Ap.h"
namespace alglib_impl {
#if !defined ALGLIB_NO_FAST_KERNELS && defined _ALGLIB_HAS_AVX2_INTRINSICS
double rdotv_avx2(const ae_int_t n, const Real *__restrict x, const Real *__restrict y, const ae_state *__restrict _state);
double rdotv2_avx2(const ae_int_t n, const Real *__restrict x, const ae_state *__restrict _state);
void rcopyv_avx2(ae_int_t n, const Real *__restrict x, Real *__restrict y, ae_state *__restrict _state);
void rcopymulv_avx2(const ae_int_t n, const double v, const Real *__restrict x, Real *__restrict y, const ae_state *__restrict _state);
void icopyv_avx2(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y, ae_state *__restrict _state);
void bcopyv_avx2(const ae_int_t n, const bool *__restrict x, bool *__restrict y, ae_state *__restrict _state);
void rsetv_avx2(const ae_int_t n, const double v, Real *__restrict x, const ae_state *__restrict _state);
void rsetvx_avx2(const ae_int_t n, const double v, double *__restrict x, const ae_state *__restrict _state);
void isetv_avx2(const ae_int_t n, const ae_int_t v, ae_int_t *__restrict x, ae_state *__restrict _state);
void bsetv_avx2(const ae_int_t n, const bool v, bool *__restrict x, ae_state *__restrict _state);
void rmulv_avx2(const ae_int_t n, const double v, double *__restrict x, const ae_state *__restrict _state);
void rmulvx_avx2(const ae_int_t n, const double v, double *__restrict x, const ae_state *__restrict _state);
void raddv_avx2(const ae_int_t n, const double alpha, const Real *__restrict y, Real *__restrict x, const ae_state *__restrict _state);
void raddvx_avx2(const ae_int_t n, const double alpha, const double *__restrict y, double *__restrict x, ae_state *_state);
void rmergemulv_avx2(ae_int_t n, const Real *__restrict y, Real *__restrict x, const ae_state *__restrict _state);
void rmergemaxv_avx2(ae_int_t n, const Real *__restrict y, Real *__restrict x, ae_state *__restrict _state);
void rmergeminv_avx2(ae_int_t n, const Real *__restrict y, Real *__restrict x, ae_state *__restrict _state);
double rmaxv_avx2(ae_int_t n, const Real *__restrict x, ae_state *__restrict _state);
double rmaxabsv_avx2(ae_int_t n, const Real *__restrict x, ae_state *__restrict _state);
void rcopyvx_avx2(const ae_int_t n, const double *__restrict x, double *__restrict y, ae_state *_state);
void icopyvx_avx2(const ae_int_t n, const ae_int_t *__restrict x, ae_int_t *__restrict y, ae_state *__restrict _state);
void rgemv_straight_avx2(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y, ae_state *_state);
void rgemv_transposed_avx2(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const double *__restrict x, double *__restrict y, ae_state *_state);
void rgemvx_straight_avx2(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y, ae_state *_state);
void rgemvx_transposed_avx2(const ae_int_t m, const ae_int_t n, const double alpha, ae_matrix *__restrict a, const ae_int_t ia, const ae_int_t ja, const double *__restrict x, double *__restrict y, ae_state *_state);
ae_int_t ablasf_packblkh_avx2(const double *src, ae_int_t src_stride, ae_int_t op, ae_int_t opsrc_length, ae_int_t opsrc_width, double *dst, ae_int_t block_size, ae_int_t micro_size);
ae_int_t ablasf_packblkh32_avx2(const double *src, ae_int_t src_stride, ae_int_t op, ae_int_t ignore_opsrc_length, ae_int_t opsrc_width, double *dst, ae_int_t ignore_block_size, ae_int_t micro_size);
void ablasf_dotblkh_avx2(const double *src_a, const double *src_b, ae_int_t round_length, ae_int_t block_size, ae_int_t micro_size, double *dst, ae_int_t dst_stride);
void ablasf_daxpby_avx2(ae_int_t n, double alpha, const double *src, double beta, double *dst);
bool spchol_updatekernelabc4_avx2(double *rowstorage, ae_int_t offss, ae_int_t twidth, ae_int_t offsu, ae_int_t uheight, ae_int_t urank, ae_int_t urowstride, ae_int_t uwidth, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase, ae_state *_state);
bool spchol_updatekernel4444_avx2(double *rowstorage, ae_int_t offss, ae_int_t sheight, ae_int_t offsu, ae_int_t uheight, double *diagd, ae_int_t offsd, ae_int_t *raw2smap, ae_int_t *superrowidx, ae_int_t urbase, ae_state *_state);
// ALGLIB_NO_FAST_KERNELS, _ALGLIB_HAS_AVX2_INTRINSICS
#endif
} // end of namespace alglib_impl

#endif // OnceOnly
