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
#ifndef OnceOnlyFastTransforms_h
#define OnceOnlyFastTransforms_h

#include "AlgLibInternal.h"

// === FFT Package ===
// Depends on: (AlgLibInternal) FTBASE
namespace alglib_impl {
void fftc1d(CVector *a, ae_int_t n);
void fftc1dinv(CVector *a, ae_int_t n);
void fftr1d(RVector *a, ae_int_t n, CVector *f);
void fftr1dinv(CVector *f, ae_int_t n, RVector *a);
void fftr1dinternaleven(RVector *a, ae_int_t n, RVector *buf, fasttransformplan *plan);
void fftr1dinvinternaleven(RVector *a, ae_int_t n, RVector *buf, fasttransformplan *plan);
} // end of namespace alglib_impl

namespace alglib {
void fftc1d(complex_1d_array &a, const ae_int_t n);
void fftc1d(complex_1d_array &a);
void fftc1dinv(complex_1d_array &a, const ae_int_t n);
void fftc1dinv(complex_1d_array &a);
void fftr1d(const real_1d_array &a, const ae_int_t n, complex_1d_array &f);
void fftr1d(const real_1d_array &a, complex_1d_array &f);
void fftr1dinv(const complex_1d_array &f, const ae_int_t n, real_1d_array &a);
void fftr1dinv(const complex_1d_array &f, real_1d_array &a);
} // end of namespace alglib

// === FHT Package ===
// Depends on: FFT
namespace alglib_impl {
void fhtr1d(RVector *a, ae_int_t n);
void fhtr1dinv(RVector *a, ae_int_t n);
} // end of namespace alglib_impl

namespace alglib {
void fhtr1d(real_1d_array &a, const ae_int_t n);
void fhtr1dinv(real_1d_array &a, const ae_int_t n);
} // end of namespace alglib

// === CONV Package ===
// Depends on: FFT
namespace alglib_impl {
void convc1d(CVector *a, ae_int_t m, CVector *b, ae_int_t n, CVector *r);
void convc1dinv(CVector *a, ae_int_t m, CVector *b, ae_int_t n, CVector *r);
void convc1dcircular(CVector *s, ae_int_t m, CVector *r, ae_int_t n, CVector *c);
void convc1dcircularinv(CVector *a, ae_int_t m, CVector *b, ae_int_t n, CVector *r);
void convr1d(RVector *a, ae_int_t m, RVector *b, ae_int_t n, RVector *r);
void convr1dinv(RVector *a, ae_int_t m, RVector *b, ae_int_t n, RVector *r);
void convr1dcircular(RVector *s, ae_int_t m, RVector *r, ae_int_t n, RVector *c);
void convr1dcircularinv(RVector *a, ae_int_t m, RVector *b, ae_int_t n, RVector *r);
void convc1dx(CVector *a, ae_int_t m, CVector *b, ae_int_t n, bool circular, ae_int_t alg, ae_int_t q, CVector *r);
void convr1dx(RVector *a, ae_int_t m, RVector *b, ae_int_t n, bool circular, ae_int_t alg, ae_int_t q, RVector *r);
} // end of namespace alglib_impl

namespace alglib {
void convc1d(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
void convc1dinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
void convc1dcircular(const complex_1d_array &s, const ae_int_t m, const complex_1d_array &r, const ae_int_t n, complex_1d_array &c);
void convc1dcircularinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
void convr1d(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
void convr1dinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
void convr1dcircular(const real_1d_array &s, const ae_int_t m, const real_1d_array &r, const ae_int_t n, real_1d_array &c);
void convr1dcircularinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
} // end of namespace alglib

// === CORR Package ===
// Depends on: CONV
namespace alglib_impl {
void corrc1d(CVector *signal, ae_int_t n, CVector *pattern, ae_int_t m, CVector *r);
void corrc1dcircular(CVector *signal, ae_int_t m, CVector *pattern, ae_int_t n, CVector *c);
void corrr1d(RVector *signal, ae_int_t n, RVector *pattern, ae_int_t m, RVector *r);
void corrr1dcircular(RVector *signal, ae_int_t m, RVector *pattern, ae_int_t n, RVector *c);
} // end of namespace alglib_impl

namespace alglib {
void corrc1d(const complex_1d_array &signal, const ae_int_t n, const complex_1d_array &pattern, const ae_int_t m, complex_1d_array &r);
void corrc1dcircular(const complex_1d_array &signal, const ae_int_t m, const complex_1d_array &pattern, const ae_int_t n, complex_1d_array &c);
void corrr1d(const real_1d_array &signal, const ae_int_t n, const real_1d_array &pattern, const ae_int_t m, real_1d_array &r);
void corrr1dcircular(const real_1d_array &signal, const ae_int_t m, const real_1d_array &pattern, const ae_int_t n, real_1d_array &c);
} // end of namespace alglib

#endif // OnceOnly
