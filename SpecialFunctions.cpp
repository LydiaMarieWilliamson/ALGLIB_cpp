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
#define InAlgLib
#include "SpecialFunctions.h"

static const double Eul = 0.5772156649015328606065, HalfPi = 1.57079632679489661923, Root2 = 1.41421356237309504880;

// === GAMMAFUNC Package ===
namespace alglib_impl {
#if !defined ALGLIB_INTERCEPTS_SPECFUNCS
static double gammafunc_gammastirf(double x) {
   double y;
   double w;
   double v;
   double stir;
   double result;
   w = 1.0 / x;
   stir = 7.87311395793093628397E-4;
   stir = -2.29549961613378126380E-4 + w * stir;
   stir = -2.68132617805781232825E-3 + w * stir;
   stir = 3.47222221605458667310E-3 + w * stir;
   stir = 8.33333333333482257126E-2 + w * stir;
   w = 1.0 + w * stir;
   y = exp(x);
   if (x > 143.01608) {
      v = pow(x, 0.5 * x - 0.25);
      y = v * (v / y);
   } else {
      y = pow(x, x - 0.5) / y;
   }
   result = 2.50662827463100050242 * y * w;
   return result;
}
#endif

// Gamma function
//
// Inputs:
//     X   -   argument
//
// Domain:
//     0 < X < 171.6
//     -170 < X < 0, X is not an integer.
//
// Relative error:
//  arithmetic   domain     # trials      peak         rms
//     IEEE    -170,-33      20000       2.3e-15     3.3e-16
//     IEEE     -33,  33     20000       9.4e-16     2.2e-16
//     IEEE      33, 171.6   20000       2.3e-15     3.2e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Original copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
// Translated to AlgoPascal by Sergey Bochkanov (2005, 2006, 2007).
// API: double gammafunction(const double x);
double gammafunction(double x) {
#if defined ALGLIB_INTERCEPTS_SPECFUNCS
   return _ialglib_i_gammafunction(x);
#else
   double p;
   double pp;
   double q;
   double qq;
   double z;
   ae_int_t i;
   double sgngam;
   double result;
   sgngam = 1.0;
   q = fabs(x);
   if (q > 33.0) {
      if (x < 0.0) {
         p = floor(q);
         i = iround(p);
         if (i % 2 == 0) {
            sgngam = -1.0;
         }
         z = q - p;
         if (z > 0.5) {
            p++;
            z = q - p;
         }
         z = q * sin(pi * z);
         z = fabs(z);
         z = pi / (z * gammafunc_gammastirf(q));
      } else {
         z = gammafunc_gammastirf(x);
      }
      result = sgngam * z;
      return result;
   }
   z = 1.0;
   while (x >= 3.0) {
      x--;
      z *= x;
   }
   while (x < 0.0) {
      if (x > -0.000000001) {
         result = z / ((1.0 + Eul * x) * x);
         return result;
      }
      z /= x;
      x++;
   }
   while (x < 2.0) {
      if (x < 0.000000001) {
         result = z / ((1.0 + Eul * x) * x);
         return result;
      }
      z /= x;
      x++;
   }
   if (x == 2.0) {
      result = z;
      return result;
   }
   x -= 2.0;
   pp = 1.60119522476751861407E-4;
   pp = 1.19135147006586384913E-3 + x * pp;
   pp = 1.04213797561761569935E-2 + x * pp;
   pp = 4.76367800457137231464E-2 + x * pp;
   pp = 2.07448227648435975150E-1 + x * pp;
   pp = 4.94214826801497100753E-1 + x * pp;
   pp = 9.99999999999999996796E-1 + x * pp;
   qq = -2.31581873324120129819E-5;
   qq = 5.39605580493303397842E-4 + x * qq;
   qq = -4.45641913851797240494E-3 + x * qq;
   qq = 1.18139785222060435552E-2 + x * qq;
   qq = 3.58236398605498653373E-2 + x * qq;
   qq = -2.34591795718243348568E-1 + x * qq;
   qq = 7.14304917030273074085E-2 + x * qq;
   qq = 1.00000000000000000320 + x * qq;
   result = z * pp / qq;
   return result;
#endif
}

// Natural logarithm of gamma function
//
// Inputs:
//     X       -   argument
//
// Result:
//     logarithm of the absolute value of the Gamma(X).
//
// Outputs:
//     SgnGam  -   sign(Gamma(X))
//
// Domain:
//     0 < X < 2.55e305
//     -2.55e305 < X < 0, X is not an integer.
//
// ACCURACY:
// arithmetic      domain        # trials     peak         rms
//    IEEE    0, 3                 28000     5.4e-16     1.1e-16
//    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
// The error criterion was relative when the function magnitude
// was greater than one but absolute when it was less than one.
//
// The following test used the relative error criterion, though
// at certain points the relative error could be much higher than
// indicated.
//    IEEE    -200, -4             10000     4.8e-16     1.3e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
// Translated to AlgoPascal by Sergey Bochkanov (2005, 2006, 2007).
// API: double lngamma(const double x, double &sgngam);
double lngamma(double x, double *sgngam) {
#if defined ALGLIB_INTERCEPTS_SPECFUNCS
   return _ialglib_i_lngamma(x, sgngam);
#else
   double a;
   double b;
   double c;
   double p;
   double q;
   double u;
   double w;
   double z;
   ae_int_t i;
   double logpi;
   double ls2pi;
   double tmp;
   double result;
   *sgngam = 0.0;
   *sgngam = 1.0;
   logpi = 1.14472988584940017414;
   ls2pi = 0.91893853320467274178;
   if (x < -34.0) {
      q = -x;
      w = lngamma(q, &tmp);
      p = floor(q);
      i = iround(p);
      if (i % 2 == 0) {
         *sgngam = -1.0;
      } else {
         *sgngam = 1.0;
      }
      z = q - p;
      if (z > 0.5) {
         p++;
         z = p - q;
      }
      z = q * sin(pi * z);
      result = logpi - log(z) - w;
      return result;
   }
   if (x < 13.0) {
      z = 1.0;
      p = 0.0;
      u = x;
      while (u >= 3.0) {
         p--;
         u = x + p;
         z *= u;
      }
      while (u < 2.0) {
         z /= u;
         p++;
         u = x + p;
      }
      if (z < 0.0) {
         *sgngam = -1.0;
         z = -z;
      } else {
         *sgngam = 1.0;
      }
      if (u == 2.0) {
         result = log(z);
         return result;
      }
      p -= 2.0;
      x += p;
      b = -1378.25152569120859100;
      b = -38801.6315134637840924 + x * b;
      b = -331612.992738871184744 + x * b;
      b = -1162370.97492762307383 + x * b;
      b = -1721737.00820839662146 + x * b;
      b = -853555.664245765465627 + x * b;
      c = 1.0;
      c = -351.815701436523470549 + x * c;
      c = -17064.2106651881159223 + x * c;
      c = -220528.590553854454839 + x * c;
      c = -1139334.44367982507207 + x * c;
      c = -2532523.07177582951285 + x * c;
      c = -2018891.41433532773231 + x * c;
      p = x * b / c;
      result = log(z) + p;
      return result;
   }
   q = (x - 0.5) * log(x) - x + ls2pi;
   if (x > 100000000.0) {
      result = q;
      return result;
   }
   p = 1.0 / (x * x);
   if (x >= 1000.0) {
      q += ((7.9365079365079365079365 * 0.0001 * p - 2.7777777777777777777778 * 0.001) * p + 0.0833333333333333333333) / x;
   } else {
      a = 8.11614167470508450300 * 0.0001;
      a = -5.95061904284301438324 * 0.0001 + p * a;
      a = 7.93650340457716943945 * 0.0001 + p * a;
      a = -2.77777777730099687205 * 0.001 + p * a;
      a = 8.33333333333331927722 * 0.01 + p * a;
      q += a / x;
   }
   result = q;
   return result;
#endif
}
} // end of namespace alglib_impl

namespace alglib {
double gammafunction(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::gammafunction(x);
   alglib_impl::ae_state_clear();
   return D;
}

double lngamma(const double x, double &sgngam) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::lngamma(x, &sgngam);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === NORMALDISTR Package ===
// Depends on: (AlgLibMisc) HQRND
namespace alglib_impl {
// Rationally-based approximation of erf(x) over x in (0.0, 0.5).
static double erfr0(double x) {
   double xx = x * x;
   const double p13 = +0.007547728033418631287834, p11 = -0.288805137207594084924010;
   const double p09 = +14.3383842191748205576712, p07 = +38.0140318123903008244444;
   const double p05 = +3017.82788536507577809226, p03 = +7404.07142710151470082064, p01 = +80437.3630960840172832162;
   double p = x * (p01 + xx * (p03 + xx * (p05 + xx * (p07 + xx * (p09 + xx * (p11 + xx * p13))))));
   const double q8 = 38.0190713951939403753468, q6 = 658.070155459240506326937;
   const double q4 = 6379.60017324428279487120, q2 = 34216.5257924628539769006, q0 = 80437.3630960840172826266;
   double q = q0 + xx * (q2 + xx * (q4 + xx * (q6 + xx * (q8 + xx))));
   return 1.1283791670955125738961589031 * p / q;
}

// Rationally-based approximation of erfc(x) over x in [0.5, 10.0).
static double erfr1(double x) {
   const double p7 = 0.5641877825507397413087057563, p6 = 9.675807882987265400604202961;
   const double p5 = 77.08161730368428609781633646, p4 = 368.5196154710010637133875746, p3 = 1143.262070703886173606073338;
   const double p2 = 2320.439590251635247384768711, p1 = 2898.0293292167655611275846, p0 = 1826.3348842295112592168999;
   double p = p0 + x * (p1 + x * (p2 + x * (p3 + x * (p4 + x * (p5 + x * (p6 + x * p7))))));
   const double q7 = 17.14980943627607849376131193, q6 = 137.1255960500622202878443578;
   const double q5 = 661.7361207107653469211984771, q4 = 2094.384367789539593790281779, q3 = 4429.612803883682726711528526;
   const double q2 = 6089.5424232724435504633068, q1 = 4958.82756472114071495438422, q0 = 1826.3348842295112595576438;
   double q = q0 + x * (q1 + x * (q2 + x * (q3 + x * (q4 + x * (q5 + x * (q6 + x * (q7 + x)))))));
   return exp(-(x * x)) * p / q;
}

// Error function
//
// The integral is
//
//                           x
//                            -
//                 2         | |          2
//   erf(x)  =  --------     |    exp( - t  ) dt.
//              sqrt(pi)   | |
//                          -
//                           0
//
// For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
// erf(x) = 1 - erfc(x).
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,1         30000       3.7e-16     1.0e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// API: double errorfunction(const double x);
double errorfunction(double x) {
   if (x == 0.0) return 0.0;
   else if (-0.5 < x && x < +0.5) return erfr0(x);
   else if (-10.0 < x && x < +10.0) return x < 0.0 ? erfr1(-x) - 1.0 : 1.0 - erfr1(x);
   else return x < 0.0 ? -1.0 : +1.0;
}

// Complementary error function
//
//  1 - erf(x) =
//
//                           inf.
//                             -
//                  2         | |          2
//   erfc(x)  =  --------     |    exp( - t  ) dt
//               sqrt(pi)   | |
//                           -
//                            x
//
// For small x, erfc(x) = 1 - erf(x); otherwise rational
// approximations are computed.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,26.6417   30000       5.7e-14     1.5e-14
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// API: double errorfunctionc(const double x);
double errorfunctionc(double x) {
   if (x == 0.0) return 1.0;
   else if (-0.5 < x && x < +0.5) return 1.0 - erfr0(x);
   else if (-10.0 < x && x < +10.0) return x < 0.0 ? 2.0 - erfr1(-x) : erfr1(x);
   else return x < 0.0 ? 2.0 : 0.0;
}

// Same as normalcdf(), obsolete name.
// API: double normaldistribution(const double x);
double normaldistribution(double x) {
   return normalcdf(x);
}

// Normal distribution PDF
//
// Returns Gaussian probability density function:
//
//                1
//    f(x)  = --------- * exp(-x^2/2)
//            sqrt(2pi)
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// API: double normalpdf(const double x);
double normalpdf(double x) {
   ae_assert(isfinite(x), "NormalPDF: X is infinite");
   return exp(-x * x / 2.0) / sqrt(2.0 * pi);
}

// Normal distribution CDF
//
// Returns the area under the Gaussian probability density
// function, integrated from minus infinity to x:
//
//                            x
//                             -
//                   1        | |          2
//    ndtr(x)  = ---------    |    exp( - t /2 ) dt
//               sqrt(2pi)  | |
//                           -
//                          -inf.
//
//             =  ( 1 + erf(z) ) / 2
//             =  erfc(z) / 2
//
// where z = x/sqrt(2). Computation is via the functions
// erf and erfc.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE     -13,0        30000       3.4e-14     6.7e-15
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// API: double normalcdf(const double x);
double normalcdf(double x) {
// return 0.5 * (errorfunction(x / Root2) + 1.0);
   x /= Root2;
   if (-0.5 < x && x < +0.5) return 0.5 * (erfr0(x) + 1.0);
   else if (-10.0 < x && x < +10.0) return x < 0.0 ? 0.5 * erfr1(-x) : 1.0 - 0.5 * erfr1(x);
   else return x < 0.0 ? 0.0 : +1.0;
}

static double inverfr0(double y) {
   const double s2pi = 2.50662827463100050242;
   double yy = y * y;
   const double p11 = -59.9633501014107895267, p09 = +98.0010754185999661536;
   const double p07 = -56.6762857469070293439, p05 = +13.9312609387279679503, p03 = -1.23916583867381258016;
   double p = y * yy * (p03 + yy * (p05 + yy * (p07 + yy * (p09 + yy * p11))));
   const double q14 = +1.95448858338141759834, q12 = +4.67627912898881538453;
   const double q10 = +86.3602421390890590575, q08 = -225.462687854119370527, q06 = +200.260212380060660359;
   const double q04 = -82.0372256168333339912, q02 = +15.9056225126211695515, q00 = -1.18331621121330003142;
   double q = q00 + yy * (q02 + yy * (q04 + yy * (q06 + yy * (q08 + yy * (q10 + yy * (q12 + yy * (q14 + yy)))))));
   return (y + p / q) * s2pi;
}

static double inverfr1(double z) {
   const double p0 = +4.05544892305962419923, p1 = +31.5251094599893866154, p2 = +57.1628192246421288162;
   const double p3 = +44.0805073893200834700, p4 = +14.6849561928858024014, p5 = +2.18663306850790267539;
   const double p6 = -1.40256079171354495875 * 0.1, p7 = -3.50424626827848203418 * 0.01, p8 = -8.57456785154685413611 * 0.0001;
   double p = z * (p8 + z * (p7 + z * (p6 + z * (p5 + z * (p4 + z * (p3 + z * (p2 + z * (p1 + z * p0))))))));
   const double q2 = +15.7799883256466749731, q3 = +45.3907635128879210584, q4 = +41.3172038254672030440;
   const double q5 = +15.0425385692907503408, q6 = +2.50464946208309415979, q7 = -1.42182922854787788574 * 0.1;
   const double q8 = -3.80806407691578277194 * 0.01, q9 = -9.33259480895457427372 * 0.0001;
   double q = q9 + z * (q8 + z * (q7 + z * (q6 + z * (q5 + z * (q4 + z * (q3 + z * (q2 + z)))))));
   return p / q;
}

static double inverfr2(double z) {
   const double p0 = 3.23774891776946035970, p1 = 6.91522889068984211695, p2 = 3.93881025292474443415;
   const double p3 = 1.33303460815807542389, p4 = 2.01485389549179081538 * 0.1, p5 = 1.23716634817820021358 * 0.01;
   const double p6 = 3.01581553508235416007 * 0.0001, p7 = 2.65806974686737550832 * 0.000001, p8 = 6.23974539184983293730 * 0.000000001;
   double p = z * (p8 + z * (p7 + z * (p6 + z * (p5 + z * (p4 + z * (p3 + z * (p2 + z * (p1 + z * p0))))))));
   const double q2 = 6.02427039364742014255, q3 = 3.67983563856160859403, q4 = 1.37702099489081330271;
   const double q5 = 2.16236993594496635890 * 0.1, q6 = 1.34204006088543189037 * 0.01, q7 = 3.28014464682127739104 * 0.0001;
   const double q8 = 2.89247864745380683936 * 0.000001, q9 = 6.79019408009981274425 * 0.000000001;
   double q = q9 + z * (q8 + z * (q7 + z * (q6 + z * (q5 + z * (q4 + z * (q3 + z * (q2 + z)))))));
   return p / q;
}

// Inverse of Normal CDF
//
// Returns the argument, x, for which the area under the
// Gaussian probability density function (integrated from
// minus infinity to x) is equal to y.
//
// For small arguments 0 < y < exp(-2), the program computes
// z = sqrt( -2.0 * log(y) );  then the approximation is
// x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
// There are two rational functions P/Q, one for 0 < y < exp(-32)
// and the other for y up to exp(-2).  For larger arguments,
// w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain        # trials      peak         rms
//    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
//    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// API: double invnormalcdf(const double y0);
double invnormalcdf(double y0) {
   const double expm2 = 0.13533528323661269189;
   if (y0 <= 0.0) return -maxrealnumber;
   else if (y0 >= 1.0) return +maxrealnumber;
   ae_int_t sign = y0 <= expm2 ? -1 : y0 > 1.0 - expm2 ? +1 : 0;
   if (sign == 0) return inverfr0(y0 - 0.5);
   else if (sign > 0) y0 = 1.0 - y0;
   double x = sqrt(-2.0 * log(y0));
   return sign * (x - log(x) / x - (x < 8.0 ? inverfr1(1.0 / x) : inverfr2(1.0 / x)));
}

// Same as invnormalcdf(), deprecated name
// API: double invnormaldistribution(const double y0);
double invnormaldistribution(double y0) {
   return invnormalcdf(y0);
}

// Inverse of the error function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// API: double inverf(const double e);
double inverf(double e) {
   return invnormaldistribution(0.5 * (e + 1.0)) / Root2;
}

// Bivariate normal PDF
//
// Returns probability density function of the bivariate  Gaussian  with
// correlation parameter equal to Rho:
//
//                          1              (    x^2 - 2*rho*x*y + y^2  )
//     f(x,y,rho) = ----------------- * exp( - ----------------------- )
//                  2pi*sqrt(1-rho^2)      (        2*(1-rho^2)        )
//
// with -1 < rho < +1 and arbitrary x, y.
//
// This function won't fail as long as Rho is in (-1,+1) range.
// ALGLIB: Copyright 15.11.2019 by Sergey Bochkanov
// API: double bivariatenormalpdf(const double x, const double y, const double rho);
double bivariatenormalpdf(double x, double y, double rho) {
   double onerho2;
   double result;
   ae_assert(isfinite(x), "BivariateNormalCDF: X is infinite");
   ae_assert(isfinite(y), "BivariateNormalCDF: Y is infinite");
   ae_assert(isfinite(rho), "BivariateNormalCDF: Rho is infinite");
   ae_assert(-1.0 < rho && rho < 1.0, "BivariateNormalCDF: Rho is not in (-1,+1) range");
   onerho2 = (1.0 - rho) * (1.0 + rho);
   result = exp(-(x * x + y * y - 2.0 * rho * x * y) / (2.0 * onerho2)) / (2.0 * pi * sqrt(onerho2));
   return result;
}

// Internal function which computes integrand of  formula  (3)  by  Alan
// Genz times Gaussian weights (passed by user).
// ALGLIB: Copyright 15.11.2019 by Sergey Bochkanov
static double normaldistr_bvnintegrate3(double rangea, double rangeb, double x, double y, double gw, double gx) {
   double r;
   double t2;
   double dd;
   double sinr;
   double cosr;
   double result;
   r = (rangeb - rangea) * 0.5 * gx + (rangeb + rangea) * 0.5;
   t2 = tan(0.5 * r);
   dd = 1.0 / (1.0 + t2 * t2);
   sinr = 2.0 * t2 * dd;
   cosr = (1.0 - t2 * t2) * dd;
   result = gw * exp(-(x * x + y * y - 2.0 * x * y * sinr) / (2.0 * cosr * cosr));
   return result;
}

// Internal function which computes integrand of  formula  (6)  by  Alan
// Genz times Gaussian weights (passed by user).
// ALGLIB: Copyright 15.11.2019 by Sergey Bochkanov
static double normaldistr_bvnintegrate6(double rangea, double rangeb, double x, double y, double s, double gw, double gx) {
   double r;
   double exphsk22x2;
   double exphsk2;
   double sqrt1x2;
   double exphsk1sqrt1x2;
   double result;
   r = (rangeb - rangea) * 0.5 * gx + (rangeb + rangea) * 0.5;
   exphsk22x2 = exp(-(x - s * y) * (x - s * y) / (2.0 * r * r));
   exphsk2 = exp(-x * s * y / 2.0);
   sqrt1x2 = sqrt((1.0 - r) * (1.0 + r));
   exphsk1sqrt1x2 = exp(-x * s * y / (1.0 + sqrt1x2));
   result = gw * exphsk22x2 * (exphsk1sqrt1x2 / sqrt1x2 - exphsk2 * (1.0 + (4.0 - x * y * s) * r * r / 8.0));
   return result;
}

// Bivariate normal CDF
//
// Returns the area under the bivariate Gaussian  PDF  with  correlation
// parameter equal to Rho, integrated from minus infinity to (x,y):
//
//                                           x      y
//                                           -      -
//                             1            | |    | |
//     bvn(x,y,rho) = -------------------   |      |   f(u,v,rho)*du*dv
//                     2pi*sqrt(1-rho^2)  | |    | |
//                                         -      -
//                                        -INF   -INF
//
// where
//
//                       (    u^2 - 2*rho*u*v + v^2  )
//     f(u,v,rho)   = exp( - ----------------------- )
//                       (        2*(1-rho^2)        )
//
// with -1 < rho < +1 and arbitrary x, y.
//
// This subroutine uses high-precision approximation scheme proposed  by
// Alan Genz in "Numerical  Computation  of  Rectangular  Bivariate  and
// Trivariate Normal and  t  probabilities",  which  computes  CDF  with
// absolute error roughly equal to 1e-14.
//
// This function won't fail as long as Rho is in (-1,+1) range.
// ALGLIB: Copyright 15.11.2019 by Sergey Bochkanov
// API: double bivariatenormalcdf(const double x, const double y, const double rho);
double bivariatenormalcdf(double x, double y, double rho) {
   double rangea;
   double rangeb;
   double s;
   double v;
   double v0;
   double v1;
   double fxys;
   double ta;
   double tb;
   double tc;
   double result;
   ae_assert(isfinite(x), "BivariateNormalCDF: X is infinite");
   ae_assert(isfinite(y), "BivariateNormalCDF: Y is infinite");
   ae_assert(isfinite(rho), "BivariateNormalCDF: Rho is infinite");
   ae_assert(-1.0 < rho && rho < 1.0, "BivariateNormalCDF: Rho is not in (-1,+1) range");
   if (rho == 0.0) {
      result = normalcdf(x) * normalcdf(y);
      return result;
   }
   if (SmallAtR(rho, 0.8)) {
   // Rho is small, compute integral using using formula (3) by Alan Genz, integrated
   // by means of 10-point Gauss-Legendre quadrature
      rangea = 0.0;
      rangeb = asin(rho);
      v = 0.0;
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.2491470458134028, -0.1252334085114689);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.2491470458134028, 0.1252334085114689);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.2334925365383548, -0.3678314989981802);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.2334925365383548, 0.3678314989981802);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.2031674267230659, -0.5873179542866175);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.2031674267230659, 0.5873179542866175);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.1600783285433462, -0.7699026741943047);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.1600783285433462, 0.7699026741943047);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.1069393259953184, -0.9041172563704749);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.1069393259953184, 0.9041172563704749);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.0471753363865118, -0.9815606342467192);
      v += normaldistr_bvnintegrate3(rangea, rangeb, x, y, 0.0471753363865118, 0.9815606342467192);
      v = v * 0.5 * (rangeb - rangea) / (2.0 * pi);
      result = normalcdf(x) * normalcdf(y) + v;
   } else {
   // Rho is large, compute integral using using formula (6) by Alan Genz, integrated
   // by means of 20-point Gauss-Legendre quadrature.
      x = -x;
      y = -y;
      s = sign(rho);
      if (s > 0.0) {
         fxys = normalcdf(-rmax2(x, y));
      } else {
         fxys = rmax2(0.0, normalcdf(-x) - normalcdf(y));
      }
      rangea = 0.0;
      rangeb = sqrt((1.0 - rho) * (1.0 + rho));
   // Compute first term (analytic integral) from formula (6)
      ta = rangeb;
      tb = fabs(x - s * y);
      tc = (4.0 - s * x * y) / 8.0;
      v0 = ta * (1.0 - tc * (tb * tb - ta * ta) / 3) * exp(-tb * tb / (2.0 * ta * ta)) - tb * (1.0 - tc * tb * tb / 3) * sqrt(2.0 * pi) * normalcdf(-tb / ta);
      v0 = v0 * exp(-s * x * y / 2.0) / (2.0 * pi);
   // Compute second term (numerical integral, 20-point Gauss-Legendre rule) from formula (6)
      v1 = 0.0;
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1527533871307258, -0.0765265211334973);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1527533871307258, 0.0765265211334973);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1491729864726037, -0.2277858511416451);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1491729864726037, 0.2277858511416451);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1420961093183820, -0.3737060887154195);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1420961093183820, 0.3737060887154195);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1316886384491766, -0.5108670019508271);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1316886384491766, 0.5108670019508271);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1181945319615184, -0.6360536807265150);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1181945319615184, 0.6360536807265150);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1019301198172404, -0.7463319064601508);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.1019301198172404, 0.7463319064601508);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0832767415767048, -0.8391169718222188);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0832767415767048, 0.8391169718222188);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0626720483341091, -0.9122344282513259);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0626720483341091, 0.9122344282513259);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0406014298003869, -0.9639719272779138);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0406014298003869, 0.9639719272779138);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0176140071391521, -0.9931285991850949);
      v1 += normaldistr_bvnintegrate6(rangea, rangeb, x, y, s, 0.0176140071391521, 0.9931285991850949);
      v1 = v1 * 0.5 * (rangeb - rangea) / (2.0 * pi);
      result = fxys - s * (v0 + v1);
   }
   result = rmax2(result, 0.0);
   result = rmin2(result, 1.0);
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double errorfunction(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::errorfunction(x);
   alglib_impl::ae_state_clear();
   return D;
}

double errorfunctionc(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::errorfunctionc(x);
   alglib_impl::ae_state_clear();
   return D;
}

double normaldistribution(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::normaldistribution(x);
   alglib_impl::ae_state_clear();
   return D;
}

double normalpdf(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::normalpdf(x);
   alglib_impl::ae_state_clear();
   return D;
}

double normalcdf(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::normalcdf(x);
   alglib_impl::ae_state_clear();
   return D;
}

double invnormalcdf(const double y0) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invnormalcdf(y0);
   alglib_impl::ae_state_clear();
   return D;
}

double invnormaldistribution(const double y0) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invnormaldistribution(y0);
   alglib_impl::ae_state_clear();
   return D;
}

double inverf(const double e) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::inverf(e);
   alglib_impl::ae_state_clear();
   return D;
}

double bivariatenormalpdf(const double x, const double y, const double rho) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::bivariatenormalpdf(x, y, rho);
   alglib_impl::ae_state_clear();
   return D;
}

double bivariatenormalcdf(const double x, const double y, const double rho) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::bivariatenormalcdf(x, y, rho);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === IBETAF Package ===
// Depends on: GAMMAFUNC, NORMALDISTR
namespace alglib_impl {
// Continued fraction expansion #1 for incomplete beta integral
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
static double ibetaf_incompletebetafe(double a, double b, double x, double big, double biginv) {
   double xk;
   double pk;
   double pkm1;
   double pkm2;
   double qk;
   double qkm1;
   double qkm2;
   double k1;
   double k2;
   double k3;
   double k4;
   double k5;
   double k6;
   double k7;
   double k8;
   double r;
   double t;
   double ans;
   double thresh;
   ae_int_t n;
   double result;
   k1 = a;
   k2 = a + b;
   k3 = a;
   k4 = a + 1.0;
   k5 = 1.0;
   k6 = b - 1.0;
   k7 = k4;
   k8 = a + 2.0;
   pkm2 = 0.0;
   qkm2 = 1.0;
   pkm1 = 1.0;
   qkm1 = 1.0;
   ans = 1.0;
   r = 1.0;
   n = 0;
   thresh = 3.0 * machineepsilon;
   do {
      xk = -x * k1 * k2 / (k3 * k4);
      pk = pkm1 + pkm2 * xk;
      qk = qkm1 + qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      xk = x * k5 * k6 / (k7 * k8);
      pk = pkm1 + pkm2 * xk;
      qk = qkm1 + qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if (qk != 0.0) {
         r = pk / qk;
      }
      if (r != 0.0) {
         t = fabs((ans - r) / r);
         ans = r;
      } else {
         t = 1.0;
      }
      if (t < thresh) {
         break;
      }
      k1++;
      k2++;
      k3 += 2.0;
      k4 += 2.0;
      k5++;
      k6--;
      k7 += 2.0;
      k8 += 2.0;
      if (fabs(qk) + fabs(pk) > big) {
         pkm2 *= biginv;
         pkm1 *= biginv;
         qkm2 *= biginv;
         qkm1 *= biginv;
      }
      if (SmallR(qk, biginv) || SmallR(pk, biginv)) {
         pkm2 *= big;
         pkm1 *= big;
         qkm2 *= big;
         qkm1 *= big;
      }
      n++;
   } while (n != 300);
   result = ans;
   return result;
}

// Continued fraction expansion #2
// for incomplete beta integral
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
static double ibetaf_incompletebetafe2(double a, double b, double x, double big, double biginv) {
   double xk;
   double pk;
   double pkm1;
   double pkm2;
   double qk;
   double qkm1;
   double qkm2;
   double k1;
   double k2;
   double k3;
   double k4;
   double k5;
   double k6;
   double k7;
   double k8;
   double r;
   double t;
   double ans;
   double z;
   double thresh;
   ae_int_t n;
   double result;
   k1 = a;
   k2 = b - 1.0;
   k3 = a;
   k4 = a + 1.0;
   k5 = 1.0;
   k6 = a + b;
   k7 = a + 1.0;
   k8 = a + 2.0;
   pkm2 = 0.0;
   qkm2 = 1.0;
   pkm1 = 1.0;
   qkm1 = 1.0;
   z = x / (1.0 - x);
   ans = 1.0;
   r = 1.0;
   n = 0;
   thresh = 3.0 * machineepsilon;
   do {
      xk = -z * k1 * k2 / (k3 * k4);
      pk = pkm1 + pkm2 * xk;
      qk = qkm1 + qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      xk = z * k5 * k6 / (k7 * k8);
      pk = pkm1 + pkm2 * xk;
      qk = qkm1 + qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if (qk != 0.0) {
         r = pk / qk;
      }
      if (r != 0.0) {
         t = fabs((ans - r) / r);
         ans = r;
      } else {
         t = 1.0;
      }
      if (t < thresh) {
         break;
      }
      k1++;
      k2--;
      k3 += 2.0;
      k4 += 2.0;
      k5++;
      k6++;
      k7 += 2.0;
      k8 += 2.0;
      if (fabs(qk) + fabs(pk) > big) {
         pkm2 *= biginv;
         pkm1 *= biginv;
         qkm2 *= biginv;
         qkm1 *= biginv;
      }
      if (SmallR(qk, biginv) || SmallR(pk, biginv)) {
         pkm2 *= big;
         pkm1 *= big;
         qkm2 *= big;
         qkm1 *= big;
      }
      n++;
   } while (n != 300);
   result = ans;
   return result;
}

// Power series for incomplete beta integral.
// Use when b*x is small and x not too close to 1.
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
static double ibetaf_incompletebetaps(double a, double b, double x, double maxgam) {
   double s;
   double t;
   double u;
   double v;
   double n;
   double t1;
   double z;
   double ai;
   double sg;
   double result;
   ai = 1.0 / a;
   u = (1.0 - b) * x;
   v = u / (a + 1.0);
   t1 = v;
   t = u;
   n = 2.0;
   s = 0.0;
   z = machineepsilon * ai;
   while (!SmallAtR(v, z)) {
      u = (n - b) * x / n;
      t *= u;
      v = t / (a + n);
      s += v;
      n++;
   }
   s += t1;
   s += ai;
   u = a * log(x);
   if (a + b < maxgam && SmallR(u, log(maxrealnumber))) {
      t = gammafunction(a + b) / (gammafunction(a) * gammafunction(b));
      s *= t * pow(x, a);
   } else {
      t = lngamma(a + b, &sg) - lngamma(a, &sg) - lngamma(b, &sg) + u + log(s);
      if (t < log(minrealnumber)) {
         s = 0.0;
      } else {
         s = exp(t);
      }
   }
   result = s;
   return result;
}

// Incomplete beta integral
//
// Returns incomplete beta integral of the arguments, evaluated
// from zero to x.  The function is defined as
//
//                  x
//     -            -
//    | (a+b)      | |  a-1     b-1
//  -----------    |   t   (1-t)   dt.
//   -     -     | |
//  | (a) | (b)   -
//                 0
//
// The domain of definition is 0 <= x <= 1.  In this
// implementation a and b are restricted to positive values.
// The integral from x to 1 may be obtained by the symmetry
// relation
//
//    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
//
// The integral is evaluated by a continued fraction expansion
// or, when b*x is small, by a power series.
//
// ACCURACY:
//
// Tested at uniformly distributed random points (a,b,x) with a and b
// in "domain" and x between 0 and 1.
//                                        Relative error
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,5         10000       6.9e-15     4.5e-16
//    IEEE      0,85       250000       2.2e-13     1.7e-14
//    IEEE      0,1000      30000       5.3e-12     6.3e-13
//    IEEE      0,10000    250000       9.3e-11     7.1e-12
//    IEEE      0,100000    10000       8.7e-10     4.8e-11
// Outputs smaller than the IEEE gradual underflow threshold
// were excluded from these statistics.
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// API: double incompletebeta(const double a, const double b, const double x);
double incompletebeta(double a, double b, double x) {
   double t;
   double xc;
   double w;
   double y;
   ae_int_t flag;
   double sg;
   double big;
   double biginv;
   double maxgam;
   double minlog;
   double maxlog;
   double result;
   big = 4.503599627370496e15;
   biginv = 2.22044604925031308085e-16;
   maxgam = 171.624376956302725;
   minlog = log(minrealnumber);
   maxlog = log(maxrealnumber);
   ae_assert(a > 0.0 && b > 0.0, "Domain error in IncompleteBeta");
   ae_assert(x >= 0.0 && x <= 1.0, "Domain error in IncompleteBeta");
   if (x == 0.0) {
      result = 0.0;
      return result;
   }
   if (x == 1.0) {
      result = 1.0;
      return result;
   }
   flag = 0;
   if (b * x <= 1.0 && x <= 0.95) {
      result = ibetaf_incompletebetaps(a, b, x, maxgam);
      return result;
   }
   w = 1.0 - x;
   if (x > a / (a + b)) {
      flag = 1;
      swapr(&a, &b);
      xc = x;
      x = w;
   } else {
      xc = w;
   }
   if (flag == 1 && b * x <= 1.0 && x <= 0.95) {
      t = ibetaf_incompletebetaps(a, b, x, maxgam);
      if (t <= machineepsilon) {
         result = 1.0 - machineepsilon;
      } else {
         result = 1.0 - t;
      }
      return result;
   }
   y = x * (a + b - 2.0) - (a - 1.0);
   if (y < 0.0) {
      w = ibetaf_incompletebetafe(a, b, x, big, biginv);
   } else {
      w = ibetaf_incompletebetafe2(a, b, x, big, biginv) / xc;
   }
   y = a * log(x);
   t = b * log(xc);
   if (a + b < maxgam && SmallR(y, maxlog) && SmallR(t, maxlog)) {
      t = pow(xc, b);
      t *= pow(x, a);
      t /= a;
      t *= w;
      t *= gammafunction(a + b) / (gammafunction(a) * gammafunction(b));
      if (flag == 1) {
         if (t <= machineepsilon) {
            result = 1.0 - machineepsilon;
         } else {
            result = 1.0 - t;
         }
      } else {
         result = t;
      }
      return result;
   }
   y += t + lngamma(a + b, &sg) - lngamma(a, &sg) - lngamma(b, &sg);
   y += log(w / a);
   if (y < minlog) {
      t = 0.0;
   } else {
      t = exp(y);
   }
   if (flag == 1) {
      if (t <= machineepsilon) {
         t = 1.0 - machineepsilon;
      } else {
         t = 1.0 - t;
      }
   }
   result = t;
   return result;
}

// Inverse of imcomplete beta integral
//
// Given y, the function finds x such that
//
//  incbet( a, b, x ) = y .
//
// The routine performs interval halving or Newton iterations to find the
// root of incbet(a,b,x) - y = 0.
//
// ACCURACY:
//
//                      Relative error:
//                x     a,b
// arithmetic   domain  domain  # trials    peak       rms
//    IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
//    IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
//    IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
// With a and b constrained to half-integer or integer values:
//    IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
//    IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
// With a = .5, b constrained to half-integer or integer values:
//    IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1996, 2000 by Stephen L. Moshier
// API: double invincompletebeta(const double a, const double b, const double y);
double invincompletebeta(double a, double b, double y) {
   double aaa;
   double bbb;
   double y0;
   double d;
   double yyy;
   double x;
   double x0;
   double x1;
   double lgm;
   double yp;
   double di;
   double dithresh;
   double yl;
   double yh;
   double xt;
   ae_int_t i;
   ae_int_t rflg;
   ae_int_t dir;
   ae_int_t nflg;
   double s;
   ae_int_t mainlooppos;
   ae_int_t ihalve;
   ae_int_t ihalvecycle;
   ae_int_t newt;
   ae_int_t newtcycle;
   ae_int_t breaknewtcycle;
   ae_int_t breakihalvecycle;
   double result;
   i = 0;
   ae_assert(y >= 0.0 && y <= 1.0, "Domain error in InvIncompleteBeta");
// special cases
   if (y == 0.0) {
      result = 0.0;
      return result;
   }
   if (y == 1.0) {
      result = 1.0;
      return result;
   }
// these initializations are not really necessary,
// but without them compiler complains about 'possibly uninitialized variables'.
   dithresh = 0.0;
   rflg = 0;
   aaa = 0.0;
   bbb = 0.0;
   y0 = 0.0;
   x = 0.0;
   yyy = 0.0;
   lgm = 0.0;
   dir = 0;
   di = 0.0;
// normal initializations
   x0 = 0.0;
   yl = 0.0;
   x1 = 1.0;
   yh = 1.0;
   nflg = 0;
   mainlooppos = 0;
   ihalve = 1;
   ihalvecycle = 2;
   newt = 3;
   newtcycle = 4;
   breaknewtcycle = 5;
   breakihalvecycle = 6;
// main loop
   while (true) {
   // start
      if (mainlooppos == 0) {
         if (a <= 1.0 || b <= 1.0) {
            dithresh = 0.000001;
            rflg = 0;
            aaa = a;
            bbb = b;
            y0 = y;
            x = aaa / (aaa + bbb);
            yyy = incompletebeta(aaa, bbb, x);
            mainlooppos = ihalve;
            continue;
         } else {
            dithresh = 0.0001;
         }
         yp = -invnormaldistribution(y);
         if (y > 0.5) {
            rflg = 1;
            aaa = b;
            bbb = a;
            y0 = 1.0 - y;
            yp = -yp;
         } else {
            rflg = 0;
            aaa = a;
            bbb = b;
            y0 = y;
         }
         lgm = (yp * yp - 3.0) / 6.0;
         x = 2.0 / (1.0 / (2.0 * aaa - 1.0) + 1.0 / (2.0 * bbb - 1.0));
         d = yp * sqrt(x + lgm) / x - (1.0 / (2.0 * bbb - 1.0) - 1.0 / (2.0 * aaa - 1.0)) * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * x));
         d *= 2.0;
         if (d < log(minrealnumber)) {
            x = 0.0;
            break;
         }
         x = aaa / (aaa + bbb * exp(d));
         yyy = incompletebeta(aaa, bbb, x);
         yp = (yyy - y0) / y0;
         if (SmallR(yp, 0.2)) {
            mainlooppos = newt;
            continue;
         }
         mainlooppos = ihalve;
         continue;
      }
   // ihalve
      if (mainlooppos == ihalve) {
         dir = 0;
         di = 0.5;
         i = 0;
         mainlooppos = ihalvecycle;
         continue;
      }
   // ihalvecycle
      if (mainlooppos == ihalvecycle) {
         if (i <= 99) {
            if (i != 0) {
               x = x0 + di * (x1 - x0);
               if (x == 1.0) {
                  x = 1.0 - machineepsilon;
               }
               if (x == 0.0) {
                  di = 0.5;
                  x = x0 + di * (x1 - x0);
                  if (x == 0.0) {
                     break;
                  }
               }
               yyy = incompletebeta(aaa, bbb, x);
               yp = (x1 - x0) / (x1 + x0);
               if (SmallR(yp, dithresh)) {
                  mainlooppos = newt;
                  continue;
               }
               yp = (yyy - y0) / y0;
               if (SmallR(yp, dithresh)) {
                  mainlooppos = newt;
                  continue;
               }
            }
            if (yyy < y0) {
               x0 = x;
               yl = yyy;
               if (dir < 0) {
                  dir = 0;
                  di = 0.5;
               } else {
                  if (dir > 3) {
                     di = 1.0 - (1.0 - di) * (1.0 - di);
                  } else {
                     if (dir > 1) {
                        di = 0.5 * di + 0.5;
                     } else {
                        di = (y0 - yyy) / (yh - yl);
                     }
                  }
               }
               dir++;
               if (x0 > 0.75) {
                  if (rflg == 1) {
                     rflg = 0;
                     aaa = a;
                     bbb = b;
                     y0 = y;
                  } else {
                     rflg = 1;
                     aaa = b;
                     bbb = a;
                     y0 = 1.0 - y;
                  }
                  x = 1.0 - x;
                  yyy = incompletebeta(aaa, bbb, x);
                  x0 = 0.0;
                  yl = 0.0;
                  x1 = 1.0;
                  yh = 1.0;
                  mainlooppos = ihalve;
                  continue;
               }
            } else {
               x1 = x;
               if (rflg == 1 && x1 < machineepsilon) {
                  x = 0.0;
                  break;
               }
               yh = yyy;
               if (dir > 0) {
                  dir = 0;
                  di = 0.5;
               } else {
                  if (dir < -3) {
                     di *= di;
                  } else {
                     if (dir < -1) {
                        di *= 0.5;
                     } else {
                        di = (yyy - y0) / (yh - yl);
                     }
                  }
               }
               dir--;
            }
            i++;
            mainlooppos = ihalvecycle;
            continue;
         } else {
            mainlooppos = breakihalvecycle;
            continue;
         }
      }
   // breakihalvecycle
      if (mainlooppos == breakihalvecycle) {
         if (x0 >= 1.0) {
            x = 1.0 - machineepsilon;
            break;
         }
         if (x <= 0.0) {
            x = 0.0;
            break;
         }
         mainlooppos = newt;
         continue;
      }
   // newt
      if (mainlooppos == newt) {
         if (nflg != 0) {
            break;
         }
         nflg = 1;
         lgm = lngamma(aaa + bbb, &s) - lngamma(aaa, &s) - lngamma(bbb, &s);
         i = 0;
         mainlooppos = newtcycle;
         continue;
      }
   // newtcycle
      if (mainlooppos == newtcycle) {
         if (i <= 7) {
            if (i != 0) {
               yyy = incompletebeta(aaa, bbb, x);
            }
            if (yyy < yl) {
               x = x0;
               yyy = yl;
            } else {
               if (yyy > yh) {
                  x = x1;
                  yyy = yh;
               } else {
                  if (yyy < y0) {
                     x0 = x;
                     yl = yyy;
                  } else {
                     x1 = x;
                     yh = yyy;
                  }
               }
            }
            if (x == 1.0 || x == 0.0) {
               mainlooppos = breaknewtcycle;
               continue;
            }
            d = (aaa - 1.0) * log(x) + (bbb - 1.0) * log(1.0 - x) + lgm;
            if (d < log(minrealnumber)) {
               break;
            }
            if (d > log(maxrealnumber)) {
               mainlooppos = breaknewtcycle;
               continue;
            }
            d = exp(d);
            d = (yyy - y0) / d;
            xt = x - d;
            if (xt <= x0) {
               yyy = (x - x0) / (x1 - x0);
               xt = x0 + 0.5 * yyy * (x - x0);
               if (xt <= 0.0) {
                  mainlooppos = breaknewtcycle;
                  continue;
               }
            }
            if (xt >= x1) {
               yyy = (x1 - x) / (x1 - x0);
               xt = x1 - 0.5 * yyy * (x1 - x);
               if (xt >= 1.0) {
                  mainlooppos = breaknewtcycle;
                  continue;
               }
            }
            x = xt;
            if (SmallR(d / x, 128.0 * machineepsilon)) {
               break;
            }
            i++;
            mainlooppos = newtcycle;
            continue;
         } else {
            mainlooppos = breaknewtcycle;
            continue;
         }
      }
   // breaknewtcycle
      if (mainlooppos == breaknewtcycle) {
         dithresh = 256.0 * machineepsilon;
         mainlooppos = ihalve;
         continue;
      }
   }
// done
   if (rflg != 0) {
      if (x <= machineepsilon) {
         x = 1.0 - machineepsilon;
      } else {
         x = 1.0 - x;
      }
   }
   result = x;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double incompletebeta(const double a, const double b, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::incompletebeta(a, b, x);
   alglib_impl::ae_state_clear();
   return D;
}

double invincompletebeta(const double a, const double b, const double y) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invincompletebeta(a, b, y);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === STUDENTTDISTR Package ===
// Depends on: IBETAF
namespace alglib_impl {
// Student's t distribution
//
// Computes the integral from minus infinity to t of the Student
// t distribution with integer k > 0 degrees of freedom:
//
//                                      t
//                                      -
//                                     | |
//              -                      |         2   -(k+1)/2
//             | ( (k+1)/2 )           |  (     x   )
//       ----------------------        |  ( 1 + --- )        dx
//                     -               |  (      k  )
//       sqrt( k pi ) | ( k/2 )        |
//                                   | |
//                                    -
//                                   -inf.
//
// Relation to incomplete beta integral:
//
//        1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
// where
//        z = k/(k + t**2).
//
// For t < -2, this is the method of computation.  For higher t,
// a direct method is derived from integration by parts.
// Since the function is symmetric about t = 0, the area under the
// right tail of the density is found by calling the function
// with -t instead of t.
//
// ACCURACY:
//
// Tested at random 1 <= k <= 25.  The "domain" refers to t.
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE     -100,-2      50000       5.9e-15     1.4e-15
//    IEEE     -2,100      500000       2.7e-15     4.9e-17
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double studenttdistribution(const ae_int_t k, const double t);
double studenttdistribution(ae_int_t k, double t) {
   double x;
   double rk;
   double z;
   double f;
   double tz;
   double p;
   double xsqk;
   ae_int_t j;
   double result;
   ae_assert(k > 0, "Domain error in StudentTDistribution");
   if (t == 0.0) {
      result = 0.5;
      return result;
   }
   if (t < -2.0) {
      rk = k;
      z = rk / (rk + t * t);
      result = 0.5 * incompletebeta(0.5 * rk, 0.5, z);
      return result;
   }
   if (t < 0.0) {
      x = -t;
   } else {
      x = t;
   }
   rk = k;
   z = 1.0 + x * x / rk;
   if (k % 2 != 0) {
      xsqk = x / sqrt(rk);
      p = atan(xsqk);
      if (k > 1) {
         f = 1.0;
         tz = 1.0;
         j = 3;
         while (j < k - 1 && tz / f > machineepsilon) {
            tz *= (j - 1) / (z * j);
            f += tz;
            j += 2;
         }
         p += f * xsqk / z;
      }
      p = p * 2.0 / pi;
   } else {
      f = 1.0;
      tz = 1.0;
      j = 2;
      while (j < k - 1 && tz / f > machineepsilon) {
         tz *= (j - 1) / (z * j);
         f += tz;
         j += 2;
      }
      p = f * x / sqrt(z * rk);
   }
   if (t < 0.0) {
      p = -p;
   }
   result = 0.5 + 0.5 * p;
   return result;
}

// Functional inverse of Student's t distribution
//
// Given probability p, finds the argument t such that stdtr(k,t)
// is equal to p.
//
// ACCURACY:
//
// Tested at random 1 <= k <= 100.  The "domain" refers to p:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE    .001,.999     25000       5.7e-15     8.0e-16
//    IEEE    10^-6,.001    25000       2.0e-12     2.9e-14
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double invstudenttdistribution(const ae_int_t k, const double p);
double invstudenttdistribution(ae_int_t k, double p) {
   double t;
   double rk;
   double z;
   ae_int_t rflg;
   double result;
   ae_assert(k > 0 && p > 0.0 && p < 1.0, "Domain error in InvStudentTDistribution");
   rk = k;
   if (p > 0.25 && p < 0.75) {
      if (p == 0.5) {
         result = 0.0;
         return result;
      }
      z = 1.0 - 2.0 * p;
      z = invincompletebeta(0.5, 0.5 * rk, fabs(z));
      t = sqrt(rk * z / (1.0 - z));
      if (p < 0.5) {
         t = -t;
      }
      result = t;
      return result;
   }
   rflg = -1;
   if (p >= 0.5) {
      p = 1.0 - p;
      rflg = 1;
   }
   z = invincompletebeta(0.5 * rk, 0.5, 2.0 * p);
   if (maxrealnumber * z < rk) {
      result = rflg * maxrealnumber;
      return result;
   }
   t = sqrt(rk / z - rk);
   result = rflg * t;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double studenttdistribution(const ae_int_t k, const double t) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::studenttdistribution(k, t);
   alglib_impl::ae_state_clear();
   return D;
}

double invstudenttdistribution(const ae_int_t k, const double p) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invstudenttdistribution(k, p);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === FDISTR Package ===
// Depends on: IBETAF
namespace alglib_impl {
// F distribution
//
// Returns the area from zero to x under the F density
// function (also known as Snedcor's density or the
// variance ratio density).  This is the density
// of x = (u1/df1)/(u2/df2), where u1 and u2 are random
// variables having Chi square distributions with df1
// and df2 degrees of freedom, respectively.
// The incomplete beta integral is used, according to the
// formula
//
// P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).
//
// The arguments a and b are greater than zero, and x is
// nonnegative.
//
// ACCURACY:
//
// Tested at random points (a,b,x).
//
//                x     a,b                     Relative error:
// arithmetic  domain  domain     # trials      peak         rms
//    IEEE      0,1    0,100       100000      9.8e-15     1.7e-15
//    IEEE      1,5    0,100       100000      6.5e-15     3.5e-16
//    IEEE      0,1    1,10000     100000      2.2e-11     3.3e-12
//    IEEE      1,5    1,10000     100000      1.1e-11     1.7e-13
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double fdistribution(const ae_int_t a, const ae_int_t b, const double x);
double fdistribution(ae_int_t a, ae_int_t b, double x) {
   double w;
   double result;
   ae_assert(a >= 1 && b >= 1 && x >= 0.0, "Domain error in FDistribution");
   w = a * x;
   w /= b + w;
   result = incompletebeta(0.5 * a, 0.5 * b, w);
   return result;
}

// Complemented F distribution
//
// Returns the area from x to infinity under the F density
// function (also known as Snedcor's density or the
// variance ratio density).
//
//                      inf.
//                       -
//              1       | |  a-1      b-1
// 1-P(x)  =  ------    |   t    (1-t)    dt
//            B(a,b)  | |
//                     -
//                      x
//
// The incomplete beta integral is used, according to the
// formula
//
// P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).
//
// ACCURACY:
//
// Tested at random points (a,b,x) in the indicated intervals.
//                x     a,b                     Relative error:
// arithmetic  domain  domain     # trials      peak         rms
//    IEEE      0,1    1,100       100000      3.7e-14     5.9e-16
//    IEEE      1,5    1,100       100000      8.0e-15     1.6e-15
//    IEEE      0,1    1,10000     100000      1.8e-11     3.5e-13
//    IEEE      1,5    1,10000     100000      2.0e-11     3.0e-12
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double fcdistribution(const ae_int_t a, const ae_int_t b, const double x);
double fcdistribution(ae_int_t a, ae_int_t b, double x) {
   double w;
   double result;
   ae_assert(a >= 1 && b >= 1 && x >= 0.0, "Domain error in FCDistribution");
   w = b / (b + a * x);
   result = incompletebeta(0.5 * b, 0.5 * a, w);
   return result;
}

// Inverse of complemented F distribution
//
// Finds the F density argument x such that the integral
// from x to infinity of the F density is equal to the
// given probability p.
//
// This is accomplished using the inverse beta integral
// function and the relations
//
//      z = incbi( df2/2, df1/2, p )
//      x = df2 (1-z) / (df1 z).
//
// Note: the following relations hold for the inverse of
// the uncomplemented F distribution:
//
//      z = incbi( df1/2, df2/2, p )
//      x = df2 z / (df1 (1-z)).
//
// ACCURACY:
//
// Tested at random points (a,b,p).
//
//              a,b                     Relative error:
// arithmetic  domain     # trials      peak         rms
// For p between .001 and 1:
//    IEEE     1,100       100000      8.3e-15     4.7e-16
//    IEEE     1,10000     100000      2.1e-11     1.4e-13
// For p between 10^-6 and 10^-3:
//    IEEE     1,100        50000      1.3e-12     8.4e-15
//    IEEE     1,10000      50000      3.0e-12     4.8e-14
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double invfdistribution(const ae_int_t a, const ae_int_t b, const double y);
double invfdistribution(ae_int_t a, ae_int_t b, double y) {
   double w;
   double result;
   ae_assert(a >= 1 && b >= 1 && y > 0.0 && y <= 1.0, "Domain error in InvFDistribution");
// Compute probability for x = 0.5
   w = incompletebeta(0.5 * b, 0.5 * a, 0.5);
// If that is greater than y, then the solution w < .5
// Otherwise, solve at 1-y to remove cancellation in (b - b*w)
   if (w > y || y < 0.001) {
      w = invincompletebeta(0.5 * b, 0.5 * a, y);
      result = (b - b * w) / (a * w);
   } else {
      w = invincompletebeta(0.5 * a, 0.5 * b, 1.0 - y);
      result = b * w / (a * (1.0 - w));
   }
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double fdistribution(const ae_int_t a, const ae_int_t b, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::fdistribution(a, b, x);
   alglib_impl::ae_state_clear();
   return D;
}

double fcdistribution(const ae_int_t a, const ae_int_t b, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::fcdistribution(a, b, x);
   alglib_impl::ae_state_clear();
   return D;
}

double invfdistribution(const ae_int_t a, const ae_int_t b, const double y) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invfdistribution(a, b, y);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === IGAMMAF Package ===
// Depends on: GAMMAFUNC, NORMALDISTR
namespace alglib_impl {
#if 0
// Forward reference to an indirect recursive call. //(@) Already declared externally.
double incompletegammac(double a, double x);
#endif

// Incomplete gamma integral
//
// The function is defined by
//
//                           x
//                            -
//                   1       | |  -t  a-1
//  igam(a,x)  =   -----     |   e   t   dt.
//                  -      | |
//                 | (a)    -
//                           0
//
// In this implementation both arguments must be positive.
// The integral is evaluated by either a power series or
// continued fraction expansion, depending on the relative
// values of a and x.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,30       200000       3.6e-14     2.9e-15
//    IEEE      0,100      300000       9.9e-14     1.5e-14
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1985, 1987, 2000 by Stephen L. Moshier
// API: double incompletegamma(const double a, const double x);
double incompletegamma(double a, double x) {
   double igammaepsilon;
   double ans;
   double ax;
   double c;
   double r;
   double tmp;
   double result;
   igammaepsilon = 0.000000000000001;
   if (x <= 0.0 || a <= 0.0) {
      result = 0.0;
      return result;
   }
   if (x > 1.0 && x > a) {
      result = 1.0 - incompletegammac(a, x);
      return result;
   }
   ax = a * log(x) - x - lngamma(a, &tmp);
   if (ax < -709.78271289338399) {
      result = 0.0;
      return result;
   }
   ax = exp(ax);
   r = a;
   c = 1.0;
   ans = 1.0;
   do {
      r++;
      c = c * x / r;
      ans += c;
   } while (c / ans > igammaepsilon);
   result = ans * ax / a;
   return result;
}

// Complemented incomplete gamma integral
//
// The function is defined by
//
//  igamc(a,x)   =   1 - igam(a,x)
//
//                            inf.
//                              -
//                     1       | |  -t  a-1
//               =   -----     |   e   t   dt.
//                    -      | |
//                   | (a)    -
//                             x
//
// In this implementation both arguments must be positive.
// The integral is evaluated by either a power series or
// continued fraction expansion, depending on the relative
// values of a and x.
//
// ACCURACY:
//
// Tested at random a, x.
//                a         x                      Relative error:
// arithmetic   domain   domain     # trials      peak         rms
//    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
//    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1985, 1987, 2000 by Stephen L. Moshier
// API: double incompletegammac(const double a, const double x);
double incompletegammac(double a, double x) {
   double igammaepsilon;
   double igammabignumber;
   double igammabignumberinv;
   double ans;
   double ax;
   double c;
   double yc;
   double r;
   double t;
   double y;
   double z;
   double pk;
   double pkm1;
   double pkm2;
   double qk;
   double qkm1;
   double qkm2;
   double tmp;
   double result;
   igammaepsilon = 0.000000000000001;
   igammabignumber = 4503599627370496.0;
   igammabignumberinv = 2.22044604925031308085 * 0.0000000000000001;
   if (x <= 0.0 || a <= 0.0) {
      result = 1.0;
      return result;
   }
   if (x < 1.0 || x < a) {
      result = 1.0 - incompletegamma(a, x);
      return result;
   }
   ax = a * log(x) - x - lngamma(a, &tmp);
   if (ax < -709.78271289338399) {
      result = 0.0;
      return result;
   }
   ax = exp(ax);
   y = 1.0 - a;
   z = x + y + 1.0;
   c = 0.0;
   pkm2 = 1.0;
   qkm2 = x;
   pkm1 = x + 1.0;
   qkm1 = z * x;
   ans = pkm1 / qkm1;
   do {
      c++;
      y++;
      z += 2.0;
      yc = y * c;
      pk = pkm1 * z - pkm2 * yc;
      qk = qkm1 * z - qkm2 * yc;
      if (qk != 0.0) {
         r = pk / qk;
         t = fabs((ans - r) / r);
         ans = r;
      } else {
         t = 1.0;
      }
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if (!SmallAtR(pk, igammabignumber)) {
         pkm2 *= igammabignumberinv;
         pkm1 *= igammabignumberinv;
         qkm2 *= igammabignumberinv;
         qkm1 *= igammabignumberinv;
      }
   } while (t > igammaepsilon);
   result = ans * ax;
   return result;
}

// Inverse of complemented imcomplete gamma integral
//
// Given p, the function finds x such that
//
//  igamc( a, x ) = p.
//
// Starting with the approximate value
//
//         3
//  x = a t
//
//  where
//
//  t = 1 - d - ndtri(p) sqrt(d)
//
// and
//
//  d = 1/9a,
//
// the routine performs up to 10 Newton iterations to find the
// root of igamc(a,x) - p = 0.
//
// ACCURACY:
//
// Tested at random a, p in the intervals indicated.
//
//                a        p                      Relative error:
// arithmetic   domain   domain     # trials      peak         rms
//    IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
//    IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
//    IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double invincompletegammac(const double a, const double y0);
double invincompletegammac(double a, double y0) {
   double igammaepsilon;
   double iinvgammabignumber;
   double x0;
   double x1;
   double x;
   double yl;
   double yh;
   double y;
   double d;
   double lgm;
   double dithresh;
   ae_int_t i;
   ae_int_t dir;
   double tmp;
   double result;
   igammaepsilon = 0.000000000000001;
   iinvgammabignumber = 4503599627370496.0;
   x0 = iinvgammabignumber;
   yl = 0.0;
   x1 = 0.0;
   yh = 1.0;
   dithresh = 5.0 * igammaepsilon;
   d = 1.0 / (9.0 * a);
   y = 1.0 - d - invnormaldistribution(y0) * sqrt(d);
   x = a * y * y * y;
   lgm = lngamma(a, &tmp);
   i = 0;
   while (i < 10) {
      if (x > x0 || x < x1) {
         d = 0.0625;
         break;
      }
      y = incompletegammac(a, x);
      if (y < yl || y > yh) {
         d = 0.0625;
         break;
      }
      if (y < y0) {
         x0 = x;
         yl = y;
      } else {
         x1 = x;
         yh = y;
      }
      d = (a - 1.0) * log(x) - x - lgm;
      if (d < -709.78271289338399) {
         d = 0.0625;
         break;
      }
      d = -exp(d);
      d = (y - y0) / d;
      if (SmallR(d / x, igammaepsilon)) {
         result = x;
         return result;
      }
      x -= d;
      i++;
   }
   if (x0 == iinvgammabignumber) {
      if (x <= 0.0) {
         x = 1.0;
      }
      while (x0 == iinvgammabignumber) {
         x *= 1.0 + d;
         y = incompletegammac(a, x);
         if (y < y0) {
            x0 = x;
            yl = y;
            break;
         }
         d += d;
      }
   }
   d = 0.5;
   dir = 0;
   i = 0;
   while (i < 400) {
      x = x1 + d * (x0 - x1);
      y = incompletegammac(a, x);
      lgm = (x0 - x1) / (x1 + x0);
      if (SmallR(lgm, dithresh)) {
         break;
      }
      lgm = (y - y0) / y0;
      if (SmallR(lgm, dithresh)) {
         break;
      }
      if (x <= 0.0) {
         break;
      }
      if (y >= y0) {
         x1 = x;
         yh = y;
         if (dir < 0) {
            dir = 0;
            d = 0.5;
         } else {
            if (dir > 1) {
               d = 0.5 * d + 0.5;
            } else {
               d = (y0 - yl) / (yh - yl);
            }
         }
         dir++;
      } else {
         x0 = x;
         yl = y;
         if (dir > 0) {
            dir = 0;
            d = 0.5;
         } else {
            if (dir < -1) {
               d *= 0.5;
            } else {
               d = (y0 - yl) / (yh - yl);
            }
         }
         dir--;
      }
      i++;
   }
   result = x;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double incompletegamma(const double a, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::incompletegamma(a, x);
   alglib_impl::ae_state_clear();
   return D;
}

double incompletegammac(const double a, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::incompletegammac(a, x);
   alglib_impl::ae_state_clear();
   return D;
}

double invincompletegammac(const double a, const double y0) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invincompletegammac(a, y0);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === CHISQUAREDISTR Package ===
// Depends on: IGAMMAF
namespace alglib_impl {
// Chi-square distribution
//
// Returns the area under the left hand tail (from 0 to x)
// of the Chi square probability density function with
// v degrees of freedom.
//
//                                   x
//                                    -
//                        1          | |  v/2-1  -t/2
//  P( x | v )   =   -----------     |   t      e     dt
//                    v/2  -       | |
//                   2    | (v/2)   -
//                                   0
//
// where x is the Chi-square variable.
//
// The incomplete gamma integral is used, according to the
// formula
//
// y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
//
// The arguments must both be positive.
//
// ACCURACY:
//
// See incomplete gamma function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double chisquaredistribution(const double v, const double x);
double chisquaredistribution(double v, double x) {
   double result;
   ae_assert(x >= 0.0 && v >= 1.0, "Domain error in ChiSquareDistribution");
   result = incompletegamma(v / 2.0, x / 2.0);
   return result;
}

// Complemented Chi-square distribution
//
// Returns the area under the right hand tail (from x to
// infinity) of the Chi square probability density function
// with v degrees of freedom:
//
//                                  inf.
//                                    -
//                        1          | |  v/2-1  -t/2
//  P( x | v )   =   -----------     |   t      e     dt
//                    v/2  -       | |
//                   2    | (v/2)   -
//                                   x
//
// where x is the Chi-square variable.
//
// The incomplete gamma integral is used, according to the
// formula
//
// y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).
//
// The arguments must both be positive.
//
// ACCURACY:
//
// See incomplete gamma function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double chisquarecdistribution(const double v, const double x);
double chisquarecdistribution(double v, double x) {
   double result;
   ae_assert(x >= 0.0 && v >= 1.0, "Domain error in ChiSquareDistributionC");
   result = incompletegammac(v / 2.0, x / 2.0);
   return result;
}

// Inverse of complemented Chi-square distribution
//
// Finds the Chi-square argument x such that the integral
// from x to infinity of the Chi-square density is equal
// to the given cumulative probability y.
//
// This is accomplished using the inverse gamma integral
// function and the relation
//
//    x/2 = igami( df/2, y );
//
// ACCURACY:
//
// See inverse incomplete gamma function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double invchisquaredistribution(const double v, const double y);
double invchisquaredistribution(double v, double y) {
   double result;
   ae_assert(y >= 0.0 && y <= 1.0 && v >= 1.0, "Domain error in InvChiSquareDistribution");
   result = 2.0 * invincompletegammac(0.5 * v, y);
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double chisquaredistribution(const double v, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::chisquaredistribution(v, x);
   alglib_impl::ae_state_clear();
   return D;
}

double chisquarecdistribution(const double v, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::chisquarecdistribution(v, x);
   alglib_impl::ae_state_clear();
   return D;
}

double invchisquaredistribution(const double v, const double y) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invchisquaredistribution(v, y);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === BINOMIALDISTR Package ===
// Depends on: (AlgLibInternal) NEARUNITYUNIT
// Depends on: IBETAF
namespace alglib_impl {
// Binomial distribution
//
// Returns the sum of the terms 0 through k of the Binomial
// probability density:
//
//   k
//   --  ( n )   j      n-j
//   >   (   )  p  (1-p)
//   --  ( j )
//  j=0
//
// The terms are not summed directly; instead the incomplete
// beta integral is employed, according to the formula
//
// y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
//
// The arguments must be positive, with p ranging from 0 to 1.
//
// ACCURACY:
//
// Tested at random points (a,b,p), with p between 0 and 1.
//
//               a,b                     Relative error:
// arithmetic  domain     # trials      peak         rms
// For p between 0.001 and 1:
//    IEEE     0,100       100000      4.3e-15     2.6e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p);
double binomialdistribution(ae_int_t k, ae_int_t n, double p) {
   double dk;
   double dn;
   double result;
   ae_assert(p >= 0.0 && p <= 1.0, "Domain error in BinomialDistribution");
   ae_assert(k >= -1 && k <= n, "Domain error in BinomialDistribution");
   if (k == -1) {
      result = 0.0;
      return result;
   }
   if (k == n) {
      result = 1.0;
      return result;
   }
   dn = n - k;
   if (k == 0) {
      dk = pow(1.0 - p, dn);
   } else {
      dk = k + 1.0;
      dk = incompletebeta(dn, dk, 1.0 - p);
   }
   result = dk;
   return result;
}

// Complemented binomial distribution
//
// Returns the sum of the terms k+1 through n of the Binomial
// probability density:
//
//   n
//   --  ( n )   j      n-j
//   >   (   )  p  (1-p)
//   --  ( j )
//  j=k+1
//
// The terms are not summed directly; instead the incomplete
// beta integral is employed, according to the formula
//
// y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
//
// The arguments must be positive, with p ranging from 0 to 1.
//
// ACCURACY:
//
// Tested at random points (a,b,p).
//
//               a,b                     Relative error:
// arithmetic  domain     # trials      peak         rms
// For p between 0.001 and 1:
//    IEEE     0,100       100000      6.7e-15     8.2e-16
// For p between 0 and .001:
//    IEEE     0,100       100000      1.5e-13     2.7e-15
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p);
double binomialcdistribution(ae_int_t k, ae_int_t n, double p) {
   double dk;
   double dn;
   double result;
   ae_assert(p >= 0.0 && p <= 1.0, "Domain error in BinomialDistributionC");
   ae_assert(k >= -1 && k <= n, "Domain error in BinomialDistributionC");
   if (k == -1) {
      result = 1.0;
      return result;
   }
   if (k == n) {
      result = 0.0;
      return result;
   }
   dn = n - k;
   if (k == 0) {
      if (p < 0.01) {
         dk = -nuexpm1(dn * nulog1p(-p));
      } else {
         dk = 1.0 - pow(1.0 - p, dn);
      }
   } else {
      dk = k + 1.0;
      dk = incompletebeta(dk, dn, p);
   }
   result = dk;
   return result;
}

// Inverse binomial distribution
//
// Finds the event probability p such that the sum of the
// terms 0 through k of the Binomial probability density
// is equal to the given cumulative probability y.
//
// This is accomplished using the inverse beta integral
// function and the relation
//
// 1 - p = incbi( n-k, k+1, y ).
//
// ACCURACY:
//
// Tested at random points (a,b,p).
//
//               a,b                     Relative error:
// arithmetic  domain     # trials      peak         rms
// For p between 0.001 and 1:
//    IEEE     0,100       100000      2.3e-14     6.4e-16
//    IEEE     0,10000     100000      6.6e-12     1.2e-13
// For p between 10^-6 and 0.001:
//    IEEE     0,100       100000      2.0e-12     1.3e-14
//    IEEE     0,10000     100000      1.5e-12     3.2e-14
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y);
double invbinomialdistribution(ae_int_t k, ae_int_t n, double y) {
   double dk;
   double dn;
   double p;
   double result;
   ae_assert(k >= 0 && k < n, "Domain error in InvBinomialDistribution");
   dn = n - k;
   if (k == 0) {
      if (y > 0.8) {
         p = -nuexpm1(nulog1p(y - 1.0) / dn);
      } else {
         p = 1.0 - pow(y, 1.0 / dn);
      }
   } else {
      dk = k + 1.0;
      p = incompletebeta(dn, dk, 0.5);
      if (p > 0.5) {
         p = invincompletebeta(dk, dn, 1.0 - y);
      } else {
         p = 1.0 - invincompletebeta(dn, dk, y);
      }
   }
   result = p;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::binomialdistribution(k, n, p);
   alglib_impl::ae_state_clear();
   return D;
}

double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::binomialcdistribution(k, n, p);
   alglib_impl::ae_state_clear();
   return D;
}

double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invbinomialdistribution(k, n, y);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === EXPINTEGRALS Package ===
namespace alglib_impl {
// Exponential integral Ei(x)
//
//               x
//                -     t
//               | |   e
//    Ei(x) =   -|-   ---  dt .
//             | |     t
//              -
//             -inf
//
// Not defined for x <= 0.
// See also expn.c.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE       0,100       50000      8.6e-16     1.3e-16
//
// Cephes Math Library Release 2.8:  May, 1999
// Copyright 1999 by Stephen L. Moshier
// API: double exponentialintegralei(const double x);
double exponentialintegralei(double x) {
   double f;
   double f1;
   double f2;
   double w;
   double result;
   if (x <= 0.0) {
      result = 0.0;
      return result;
   }
   if (x < 2.0) {
      f1 = -5.350447357812542947283;
      f1 = f1 * x + 218.5049168816613393830;
      f1 = f1 * x - 4176.572384826693777058;
      f1 = f1 * x + 55411.76756393557601232;
      f1 = f1 * x - 331338.1331178144034309;
      f1 = f1 * x + 1592627.163384945414220;
      f2 = 1.000000000000000000000;
      f2 = f2 * x - 52.50547959112862969197;
      f2 = f2 * x + 1259.616186786790571525;
      f2 = f2 * x - 17565.49581973534652631;
      f2 = f2 * x + 149306.2117002725991967;
      f2 = f2 * x - 729494.9239640527645655;
      f2 = f2 * x + 1592627.163384945429726;
      f = f1 / f2;
      result = Eul + log(x) + x * f;
      return result;
   }
   if (x < 4.0) {
      w = 1.0 / x;
      f1 = 1.981808503259689673238E-2;
      f1 = f1 * w - 1.271645625984917501326;
      f1 = f1 * w - 2.088160335681228318920;
      f1 = f1 * w + 2.755544509187936721172;
      f1 = f1 * w - 4.409507048701600257171E-1;
      f1 = f1 * w + 4.665623805935891391017E-2;
      f1 = f1 * w - 1.545042679673485262580E-3;
      f1 = f1 * w + 7.059980605299617478514E-5;
      f2 = 1.000000000000000000000;
      f2 = f2 * w + 1.476498670914921440652;
      f2 = f2 * w + 5.629177174822436244827E-1;
      f2 = f2 * w + 1.699017897879307263248E-1;
      f2 = f2 * w + 2.291647179034212017463E-2;
      f2 = f2 * w + 4.450150439728752875043E-3;
      f2 = f2 * w + 1.727439612206521482874E-4;
      f2 = f2 * w + 3.953167195549672482304E-5;
      f = f1 / f2;
      result = exp(x) * w * (1.0 + w * f);
      return result;
   }
   if (x < 8.0) {
      w = 1.0 / x;
      f1 = -1.373215375871208729803;
      f1 = f1 * w - 7.084559133740838761406E-1;
      f1 = f1 * w + 1.580806855547941010501;
      f1 = f1 * w - 2.601500427425622944234E-1;
      f1 = f1 * w + 2.994674694113713763365E-2;
      f1 = f1 * w - 1.038086040188744005513E-3;
      f1 = f1 * w + 4.371064420753005429514E-5;
      f1 = f1 * w + 2.141783679522602903795E-6;
      f2 = 1.000000000000000000000;
      f2 = f2 * w + 8.585231423622028380768E-1;
      f2 = f2 * w + 4.483285822873995129957E-1;
      f2 = f2 * w + 7.687932158124475434091E-2;
      f2 = f2 * w + 2.449868241021887685904E-2;
      f2 = f2 * w + 8.832165941927796567926E-4;
      f2 = f2 * w + 4.590952299511353531215E-4;
      f2 = f2 * w + (-4.729848351866523044863E-6);
      f2 = f2 * w + 2.665195537390710170105E-6;
      f = f1 / f2;
      result = exp(x) * w * (1.0 + w * f);
      return result;
   }
   if (x < 16.0) {
      w = 1.0 / x;
      f1 = -2.106934601691916512584;
      f1 = f1 * w + 1.732733869664688041885;
      f1 = f1 * w - 2.423619178935841904839E-1;
      f1 = f1 * w + 2.322724180937565842585E-2;
      f1 = f1 * w + 2.372880440493179832059E-4;
      f1 = f1 * w - 8.343219561192552752335E-5;
      f1 = f1 * w + 1.363408795605250394881E-5;
      f1 = f1 * w - 3.655412321999253963714E-7;
      f1 = f1 * w + 1.464941733975961318456E-8;
      f1 = f1 * w + 6.176407863710360207074E-10;
      f2 = 1.000000000000000000000;
      f2 = f2 * w - 2.298062239901678075778E-1;
      f2 = f2 * w + 1.105077041474037862347E-1;
      f2 = f2 * w - 1.566542966630792353556E-2;
      f2 = f2 * w + 2.761106850817352773874E-3;
      f2 = f2 * w - 2.089148012284048449115E-4;
      f2 = f2 * w + 1.708528938807675304186E-5;
      f2 = f2 * w - 4.459311796356686423199E-7;
      f2 = f2 * w + 1.394634930353847498145E-8;
      f2 = f2 * w + 6.150865933977338354138E-10;
      f = f1 / f2;
      result = exp(x) * w * (1.0 + w * f);
      return result;
   }
   if (x < 32.0) {
      w = 1.0 / x;
      f1 = -2.458119367674020323359E-1;
      f1 = f1 * w - 1.483382253322077687183E-1;
      f1 = f1 * w + 7.248291795735551591813E-2;
      f1 = f1 * w - 1.348315687380940523823E-2;
      f1 = f1 * w + 1.342775069788636972294E-3;
      f1 = f1 * w - 7.942465637159712264564E-5;
      f1 = f1 * w + 2.644179518984235952241E-6;
      f1 = f1 * w - 4.239473659313765177195E-8;
      f2 = 1.000000000000000000000;
      f2 = f2 * w - 1.044225908443871106315E-1;
      f2 = f2 * w - 2.676453128101402655055E-1;
      f2 = f2 * w + 9.695000254621984627876E-2;
      f2 = f2 * w - 1.601745692712991078208E-2;
      f2 = f2 * w + 1.496414899205908021882E-3;
      f2 = f2 * w - 8.462452563778485013756E-5;
      f2 = f2 * w + 2.728938403476726394024E-6;
      f2 = f2 * w - 4.239462431819542051337E-8;
      f = f1 / f2;
      result = exp(x) * w * (1.0 + w * f);
      return result;
   }
   if (x < 64.0) {
      w = 1.0 / x;
      f1 = 1.212561118105456670844E-1;
      f1 = f1 * w - 5.823133179043894485122E-1;
      f1 = f1 * w + 2.348887314557016779211E-1;
      f1 = f1 * w - 3.040034318113248237280E-2;
      f1 = f1 * w + 1.510082146865190661777E-3;
      f1 = f1 * w - 2.523137095499571377122E-5;
      f2 = 1.000000000000000000000;
      f2 = f2 * w - 1.002252150365854016662;
      f2 = f2 * w + 2.928709694872224144953E-1;
      f2 = f2 * w - 3.337004338674007801307E-2;
      f2 = f2 * w + 1.560544881127388842819E-3;
      f2 = f2 * w - 2.523137093603234562648E-5;
      f = f1 / f2;
      result = exp(x) * w * (1.0 + w * f);
      return result;
   }
   w = 1.0 / x;
   f1 = -7.657847078286127362028E-1;
   f1 = f1 * w + 6.886192415566705051750E-1;
   f1 = f1 * w - 2.132598113545206124553E-1;
   f1 = f1 * w + 3.346107552384193813594E-2;
   f1 = f1 * w - 3.076541477344756050249E-3;
   f1 = f1 * w + 1.747119316454907477380E-4;
   f1 = f1 * w - 6.103711682274170530369E-6;
   f1 = f1 * w + 1.218032765428652199087E-7;
   f1 = f1 * w - 1.086076102793290233007E-9;
   f2 = 1.000000000000000000000;
   f2 = f2 * w - 1.888802868662308731041;
   f2 = f2 * w + 1.066691687211408896850;
   f2 = f2 * w - 2.751915982306380647738E-1;
   f2 = f2 * w + 3.930852688233823569726E-2;
   f2 = f2 * w - 3.414684558602365085394E-3;
   f2 = f2 * w + 1.866844370703555398195E-4;
   f2 = f2 * w - 6.345146083130515357861E-6;
   f2 = f2 * w + 1.239754287483206878024E-7;
   f2 = f2 * w - 1.086076102793126632978E-9;
   f = f1 / f2;
   result = exp(x) * w * (1.0 + w * f);
   return result;
}

// Exponential integral En(x)
//
// Evaluates the exponential integral
//
//                 inf.
//                   -
//                  | |   -xt
//                  |    e
//      E (x)  =    |    ----  dt.
//       n          |      n
//                | |     t
//                 -
//                  1
//
// Both n and x must be nonnegative.
//
// The routine employs either a power series, a continued
// fraction, or an asymptotic formula depending on the
// relative values of n and x.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       10000       1.7e-15     3.6e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1985, 2000 by Stephen L. Moshier
// API: double exponentialintegralen(const double x, const ae_int_t n);
double exponentialintegralen(double x, ae_int_t n) {
   double r;
   double t;
   double yk;
   double xk;
   double pk;
   double pkm1;
   double pkm2;
   double qk;
   double qkm1;
   double qkm2;
   double psi;
   double z;
   ae_int_t i;
   ae_int_t k;
   double big;
   double result;
   big = 1.44115188075855872 * pow(10.0, 17.0);
   if (n < 0 || x < 0.0 || x > 170.0 || x == 0.0 && n < 2) {
      result = -1.0;
      return result;
   }
   if (x == 0.0) {
      result = 1.0 / (n - 1);
      return result;
   }
   if (n == 0) {
      result = exp(-x) / x;
      return result;
   }
   if (n > 5000) {
      xk = x + n;
      yk = 1.0 / (xk * xk);
      t = n;
      result = yk * t * (6.0 * x * x - 8.0 * t * x + t * t);
      result = yk * (result + t * (t - 2.0 * x));
      result = yk * (result + t);
      result = (result + 1.0) * exp(-x) / xk;
      return result;
   }
   if (x <= 1.0) {
      psi = -Eul - log(x);
      for (i = 1; i < n; i++) {
         psi += 1.0 / i;
      }
      z = -x;
      xk = 0.0;
      yk = 1.0;
      pk = 1 - n;
      if (n == 1) {
         result = 0.0;
      } else {
         result = 1.0 / pk;
      }
      do {
         xk++;
         yk = yk * z / xk;
         pk++;
         if (pk != 0.0) {
            result += yk / pk;
         }
      } while (result == 0.0 || !SmallR(yk / result, machineepsilon));
      t = 1.0;
      for (i = 1; i < n; i++) {
         t = t * z / i;
      }
      result = psi * t - result;
      return result;
   } else {
      k = 1;
      pkm2 = 1.0;
      qkm2 = x;
      pkm1 = 1.0;
      qkm1 = x + n;
      result = pkm1 / qkm1;
      do {
         k++;
         if (k % 2 == 1) {
            yk = 1.0;
            xk = n + (k - 1) / 2.0;
         } else {
            yk = x;
            xk = k / 2.0;
         }
         pk = pkm1 * yk + pkm2 * xk;
         qk = qkm1 * yk + qkm2 * xk;
         if (qk != 0.0) {
            r = pk / qk;
            t = fabs((result - r) / r);
            result = r;
         } else {
            t = 1.0;
         }
         pkm2 = pkm1;
         pkm1 = pk;
         qkm2 = qkm1;
         qkm1 = qk;
         if (!SmallAtR(pk, big)) {
            pkm2 /= big;
            pkm1 /= big;
            qkm2 /= big;
            qkm1 /= big;
         }
      } while (t >= machineepsilon);
      result *= exp(-x);
   }
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double exponentialintegralei(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::exponentialintegralei(x);
   alglib_impl::ae_state_clear();
   return D;
}

double exponentialintegralen(const double x, const ae_int_t n) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::exponentialintegralen(x, n);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === JACOBIANELLIPTIC Package ===
namespace alglib_impl {
// Jacobian Elliptic Functions
//
// Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
// and dn(u|m) of parameter m between 0 and 1, and real
// argument u.
//
// These functions are periodic, with quarter-period on the
// real axis equal to the complete elliptic integral
// ellpk(1.0-m).
//
// Relation to incomplete elliptic integral:
// If u = ellik(phi,m), then sn(u|m) = sin(phi),
// and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
//
// Computation is by means of the arithmetic-geometric mean
// algorithm, except when m is within 1e-9 of 0 or 1.  In the
// latter case with m close to 1, the approximation applies
// only for phi < pi/2.
//
// ACCURACY:
//
// Tested at random points with u between 0 and 10, m between
// 0 and 1.
//
//            Absolute error (* = relative error):
// arithmetic   function   # trials      peak         rms
//    IEEE      phi         10000       9.2e-16*    1.4e-16*
//    IEEE      sn          50000       4.1e-15     4.6e-16
//    IEEE      cn          40000       3.6e-15     4.4e-16
//    IEEE      dn          10000       1.3e-12     1.8e-14
//
//  Peak error observed in consistency check using addition
// theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
// the above relation to the incomplete elliptic integral.
// Accuracy deteriorates when u is large.
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph);
void jacobianellipticfunctions(double u, double m, double *sn, double *cn, double *dn, double *ph) {
   ae_frame _frame_block;
   double ai;
   double b;
   double phi;
   double t;
   double twon;
   ae_int_t i;
   ae_frame_make(&_frame_block);
   *sn = 0.0;
   *cn = 0.0;
   *dn = 0.0;
   *ph = 0.0;
   NewVector(a, 0, DT_REAL);
   NewVector(c, 0, DT_REAL);
   ae_assert(m >= 0.0 && m <= 1.0, "Domain error in JacobianEllipticFunctions: m < 0 or m > 1");
   ae_vector_set_length(&a, 9);
   ae_vector_set_length(&c, 9);
   if (m < 1.0E-9) {
      t = sin(u);
      b = cos(u);
      ai = 0.25 * m * (u - t * b);
      *sn = t - ai * b;
      *cn = b + ai * t;
      *ph = u - ai;
      *dn = 1.0 - 0.5 * m * t * t;
      ae_frame_leave();
      return;
   }
   if (m >= 0.9999999999) {
      ai = 0.25 * (1.0 - m);
      b = cosh(u);
      t = tanh(u);
      phi = 1.0 / b;
      twon = b * sinh(u);
      *sn = t + ai * (twon - u) / (b * b);
      *ph = 2.0 * atan(exp(u)) - HalfPi + ai * (twon - u) / b;
      ai *= t * phi;
      *cn = phi - ai * (twon - u);
      *dn = phi + ai * (twon + u);
      ae_frame_leave();
      return;
   }
   a.xR[0] = 1.0;
   b = sqrt(1.0 - m);
   c.xR[0] = sqrt(m);
   twon = 1.0;
   i = 0;
   while (!SmallAtR(c.xR[i] / a.xR[i], machineepsilon)) {
      if (i > 7) {
         ae_assert(false, "Overflow in JacobianEllipticFunctions");
         break;
      }
      ai = a.xR[i];
      i++;
      c.xR[i] = 0.5 * (ai - b);
      t = sqrt(ai * b);
      a.xR[i] = 0.5 * (ai + b);
      b = t;
      twon *= 2.0;
   }
   phi = twon * a.xR[i] * u;
   do {
      t = c.xR[i] * sin(phi) / a.xR[i];
      b = phi;
      phi = (asin(t) + phi) / 2.0;
      i--;
   } while (i != 0);
   *sn = sin(phi);
   t = cos(phi);
   *cn = t;
   *dn = t / cos(phi - b);
   *ph = phi;
   ae_frame_leave();
}
} // end of namespace alglib_impl

namespace alglib {
void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::jacobianellipticfunctions(u, m, &sn, &cn, &dn, &ph);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === TRIGINTEGRALS Package ===
namespace alglib_impl {
// Sine and cosine integrals
//
// Evaluates the integrals
//
//                          x
//                          -
//                         |  cos t - 1
//   Ci(x) = eul + ln x +  |  --------- dt,
//                         |      t
//                        -
//                         0
//             x
//             -
//            |  sin t
//   Si(x) =  |  ----- dt
//            |    t
//           -
//            0
//
// where eul = 0.57721566490153286061 is Euler's constant.
// The integrals are approximated by rational functions.
// For x > 8 auxiliary functions f(x) and g(x) are employed
// such that
//
// Ci(x) = f(x) sin(x) - g(x) cos(x)
// Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)
//
// ACCURACY:
//    Test interval = [0,50].
// Absolute error, except relative when > 1:
// arithmetic   function   # trials      peak         rms
//    IEEE        Si        30000       4.4e-16     7.3e-17
//    IEEE        Ci        30000       6.9e-16     5.1e-17
//
// Cephes Math Library Release 2.1:  January, 1989
// Copyright 1984, 1987, 1989 by Stephen L. Moshier
// API: void sinecosineintegrals(const double x, double &si, double &ci);
void sinecosineintegrals(double x, double *si, double *ci) {
   *si = 0.0;
   *ci = 0.0;
   bool neg = x < 0.0;
   if (neg) x = -x;
   if (x == 0.0) {
      *si = 0.0;
      *ci = -maxrealnumber;
      return;
   } else if (x > 1.0E9) {
      *si = HalfPi - cos(x) / x;
      *ci = sin(x) / x;
      return;
   } else if (x <= 4.0) {
      double xx = x * x;
      const double sn11 = -8.39167827910303881427E-11, sn09 = +4.62591714427012837309E-8, sn07 = -9.75759303843632795789E-6;
      const double sn05 = +9.76945438170435310816E-4, sn03 = -4.13470316229406538752E-2, sn01 = +1.00000000000000000302E0;
      double sn = (((((sn11 * xx + sn09) * xx + sn07) * xx + sn05) * xx + sn03) * xx + sn01) * x;
      const double sd10 = +2.03269266195951942049E-12, sd08 = +1.27997891179943299903E-9, sd06 = +4.41827842801218905784E-7;
      const double sd04 = +9.96412122043875552487E-5, sd02 = +1.42085239326149893930E-2, sd00 = +9.99999999999999996984E-1;
      double sd = ((((sd10 * xx + sd08) * xx + sd06) * xx + sd04) * xx + sd02) * xx + sd00;
      double s = sn / sd;
      const double cn12 = +2.02524002389102268789E-11, cn10 = -1.35249504915790756375E-8, cn08 = +3.59325051419993077021E-6;
      const double cn06 = -4.74007206873407909465E-4, cn04 = +2.89159652607555242092E-2, cn02 = -1.00000000000000000080E0;
      double cn = (((((cn12 * xx + cn10) * xx + cn08) * xx + cn06) * xx + cn04) * xx + cn02) * xx;
      const double cd10 = +4.07746040061880559506E-12, cd08 = +3.06780997581887812692E-9, cd06 = +1.23210355685883423679E-6;
      const double cd04 = +3.17442024775032769882E-4, cd02 = +5.10028056236446052392E-2, cd00 = +4.00000000000000000080E0;
      double cd = ((((cd10 * xx + cd08) * xx + cd06) * xx + cd04) * xx + cd02) * xx + cd00;
      double c = cn / cd;
      if (neg) s = -s;
      *si = s;
      *ci = Eul + log(x) + c;
      return;
   }
   double s = sin(x), c = cos(x);
   double w = 1.0 / x, ww = 1.0 / (x * x);
   double f, g;
   if (x < 8.0) {
      const double fn01 = 4.23612862892216586994E0, fn03 = 5.45937717161812843388E0, fn05 = 1.62083287701538329132E0;
      const double fn07 = 1.67006611831323023771E-1, fn09 = 6.81020132472518137426E-3;
      const double fn11 = 1.08936580650328664411E-4, fn13 = 5.48900223421373614008E-7;
      double fn = ((((((fn01 * ww + fn03) * ww + fn05) * ww + fn07) * ww + fn09) * ww + fn11) * ww + fn13) * w;
      const double fd00 = 8.16496634205391016773E0, fd02 = 7.30828822505564552187E0, fd04 = 1.86792257950184183883E0;
      const double fd06 = 1.78792052963149907262E-1, fd08 = 7.01710668322789753610E-3;
      const double fd10 = 1.10034357153915731354E-4, fd12 = 5.48900252756255700982E-7;
      double fd = ((((((ww + fd00) * ww + fd02) * ww + fd04) * ww + fd06) * ww + fd08) * ww + fd10) * ww + fd12;
      f = fn / fd;
      const double gn00 = 8.71001698973114191777E-2, gn02 = 6.11379109952219284151E-1, gn04 = 3.97180296392337498885E-1;
      const double gn06 = 7.48527737628469092119E-2, gn08 = 5.38868681462177273157E-3, gn10 = 1.61999794598934024525E-4;
      const double gn12 = 1.97963874140963632189E-6, gn14 = 7.82579040744090311069E-9, gn16 = 1.00000000000000000000E0;
      double gn = (((((((gn00 * ww + gn02) * ww + gn04) * ww + gn06) * ww + gn08) * ww + gn10) * ww + gn12) * ww + gn14) * ww;
      const double gd04 = 1.64402202413355338886E0, gd06 = 6.66296701268987968381E-1, gd08 = 9.88771761277688796203E-2;
      const double gd10 = 6.22396345441768420760E-3, gd12 = 1.73221081474177119497E-4;
      const double gd14 = 2.02659182086343991969E-6, gd16 = 7.82579218933534490868E-9;
      double gd = ((((((ww + gd04) * ww + gd06) * ww + gd08) * ww + gd10) * ww + gd12) * ww + gd14) * ww + gd16;
      g = gn / gd;
   } else {
      const double fn00 = 4.55880873470465315206E-1, fn02 = 7.13715274100146711374E-1, fn04 = 1.60300158222319456320E-1;
      const double fn06 = 1.16064229408124407915E-2, fn08 = 3.49556442447859055605E-4, fn10 = 4.86215430826454749482E-6;
      const double fn12 = 3.20092790091004902806E-8, fn14 = 9.41779576128512936592E-11, fn16 = 9.70507110881952024631E-14;
      double fn = ((((((((fn00 * ww + fn02) * ww + fn04) * ww + fn06) * ww + fn08) * ww + fn10) * ww + fn12) * ww + fn14) * ww + fn16) * w;
      const double fd03 = 9.17463611873684053703E-1, fd05 = 1.78685545332074536321E-1, fd07 = 1.22253594771971293032E-2;
      const double fd09 = 3.58696481881851580297E-4, fd11 = 4.92435064317881464393E-6, fd13 = 3.21956939101046018377E-8;
      const double fd15 = 9.43720590350276732376E-11, fd17 = 9.70507110881952025725E-14;
      double fd = (((((((ww + fd03) * ww + fd05) * ww + fd07) * ww + fd09) * ww + fd11) * ww + fd13) * ww + fd15) * ww + fd17;
      f = fn / fd;
      const double gn00 = 6.97359953443276214934E-1, gn02 = 3.30410979305632063225E-1, gn04 = 3.84878767649974295920E-2;
      const double gn06 = 1.71718239052347903558E-3, gn08 = 3.48941165502279436777E-5, gn10 = 3.47131167084116673800E-7;
      const double gn12 = 1.70404452782044526189E-9, gn14 = 3.85945925430276600453E-12, gn16 = 3.14040098946363334640E-15;
      double gn = ((((((((gn00 * ww + gn02) * ww + gn04) * ww + gn06) * ww + gn08) * ww + gn10) * ww + gn12) * ww + gn14) * ww + gn16) * ww;
      const double gd02 = 1.68548898811011640017E0, gd04 = 4.87852258695304967486E-1, gd06 = 4.67913194259625806320E-2;
      const double gd08 = 1.90284426674399523638E-3, gd10 = 3.68475504442561108162E-5, gd12 = 3.57043223443740838771E-7;
      const double gd14 = 1.72693748966316146736E-9, gd16 = 3.87830166023954706752E-12, gd18 = 3.14040098946363335242E-15;
      double gd = ((((((((ww + gd02) * ww + gd04) * ww + gd06) * ww + gd08) * ww + gd10) * ww + gd12) * ww + gd14) * ww + gd16) * ww + gd18;
      g = gn / gd;
   }
   *si = HalfPi - f * c - g * s;
   if (neg) *si = -*si;
   *ci = f * s - g * c;
}

// Hyperbolic sine and cosine integrals
//
// Approximates the integrals
//
//                            x
//                            -
//                           | |   cosh t - 1
//   Chi(x) = eul + ln x +   |    -----------  dt,
//                         | |          t
//                          -
//                          0
//
//               x
//               -
//              | |  sinh t
//   Shi(x) =   |    ------  dt
//            | |       t
//             -
//             0
//
// where eul = 0.57721566490153286061 is Euler's constant.
// The integrals are evaluated by power series for x < 8
// and by Chebyshev expansions for x between 8 and 88.
// For large x, both functions approach exp(x)/2x.
// Arguments greater than 88 in magnitude return MAXNUM.
//
// ACCURACY:
//
// Test interval 0 to 88.
//                      Relative error:
// arithmetic   function  # trials      peak         rms
//    IEEE         Shi      30000       6.9e-16     1.6e-16
//        Absolute error, except relative when |Chi| > 1:
//    IEEE         Chi      30000       8.4e-16     1.4e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi);
void hyperbolicsinecosineintegrals(double x, double *shi, double *chi) {
   ae_int_t neg = x < 0.0;
   if (neg) x = -x;
   if (x == 0.0) {
      *shi = 0.0;
      *chi = -maxrealnumber;
      return;
   }
   double s, c;
   if (x < 8.0) {
      double z = x * x, k = 2.0, a = 1.0;
      s = 1.0;
      c = 0.0;
      do {
         a = a * z / k;
         c += a / k;
         k++;
         a /= k;
         s += a / k;
         k++;
      } while (!SmallR(a / s, machineepsilon));
      s *= x;
      c += Eul + log(x);
   } else if (x < 18.0) {
      double k = exp(x) / x, a = (576.0 / x - 52.0) / 10.0;
      const double s21 = +1.83889230173399459482E-17, s20 = -9.55485532279655569575E-17;
      const double s19 = +2.04326105980879882648E-16, s18 = +1.09896949074905343022E-15;
      const double s17 = -1.31313534344092599234E-14, s16 = +5.93976226264314278932E-14, s15 = -3.47197010497749154755E-14;
      const double s14 = -1.40059764613117131000E-12, s13 = +9.49044626224223543299E-12, s12 = -1.61596181145435454033E-11;
      const double s11 = -1.77899784436430310321E-10, s10 = +1.35455469767246947469E-9, s09 = -1.03257121792819495123E-9;
      const double s08 = -3.56699611114982536845E-8, s07 = +1.44818877384267342057E-7, s06 = +7.82018215184051295296E-7;
      const double s05 = -5.39919118403805073710E-6, s04 = -3.12458202168959833422E-5, s03 = +8.90136741950727517826E-5;
      const double s02 = +2.02558474743846862168E-3, s01 = +2.96064440855633256972E-2, s00 = +1.11847751047257036625E0;
      double b0 = s21, b2 = a*b0 + s20, b1 = a*b2 - b0 + s19;
      b0 = a*b1 - b2 + s18;
      b2 = a*b0 - b1 + s17, b1 = a*b2 - b0 + s16, b0 = a*b1 - b2 + s15;
      b2 = a*b0 - b1 + s14, b1 = a*b2 - b0 + s13, b0 = a*b1 - b2 + s12;
      b2 = a*b0 - b1 + s11, b1 = a*b2 - b0 + s10, b0 = a*b1 - b2 + s09;
      b2 = a*b0 - b1 + s08, b1 = a*b2 - b0 + s07, b0 = a*b1 - b2 + s06;
      b2 = a*b0 - b1 + s05, b1 = a*b2 - b0 + s04, b0 = a*b1 - b2 + s03;
      b2 = a*b0 - b1 + s02, b1 = a*b2 - b0 + s01, b0 = a*b1 - b2 + s00;
      s = k * 0.5 * (b0 - b2);
      const double c22 = -8.12435385225864036372E-18, c21 = +2.17586413290339214377E-17;
      const double c20 = +5.22624394924072204667E-17, c19 = -9.48812110591690559363E-16, c18 = +5.35546311647465209166E-15;
      const double c17 = -1.21009970113732918701E-14, c16 = -6.00865178553447437951E-14, c15 = +7.16339649156028587775E-13;
      const double c14 = -2.93496072607599856104E-12, c13 = -1.40359438136491256904E-12, c12 = +8.76302288609054966081E-11;
      const double c11 = -4.40092476213282340617E-10, c10 = -1.87992075640569295479E-10, c09 = +1.31458150989474594064E-8;
      const double c08 = -4.75513930924765465590E-8, c07 = -2.21775018801848880741E-7, c06 = +1.94635531373272490962E-6;
      const double c05 = +4.33505889257316408893E-6, c04 = -6.13387001076494349496E-5, c03 = -3.13085477492997465138E-4;
      const double c02 = +4.97164789823116062801E-4, c01 = +2.64347496031374526641E-2, c00 = +1.11446150876699213025E0;
      b1 = c22, b0 = a*b1 + c21;
      b2 = a*b0 - b1 + c20, b1 = a*b2 - b0 + c19, b0 = a*b1 - b2 + c18;
      b2 = a*b0 - b1 + c17, b1 = a*b2 - b0 + c16, b0 = a*b1 - b2 + c15;
      b2 = a*b0 - b1 + c14, b1 = a*b2 - b0 + c13, b0 = a*b1 - b2 + c12;
      b2 = a*b0 - b1 + c11, b1 = a*b2 - b0 + c10, b0 = a*b1 - b2 + c09;
      b2 = a*b0 - b1 + c08, b1 = a*b2 - b0 + c07, b0 = a*b1 - b2 + c06;
      b2 = a*b0 - b1 + c05, b1 = a*b2 - b0 + c04, b0 = a*b1 - b2 + c03;
      b2 = a*b0 - b1 + c02, b1 = a*b2 - b0 + c01, b0 = a*b1 - b2 + c00;
      c = Eul + log(x) + k * 0.5 * (b0 - b2);
   } else if (x <= 88.0) {
      double k = exp(x) / x, a = (6336.0 / x - 212.0) / 70.0;
      const double s22 = -1.05311574154850938805E-17, s21 = +2.62446095596355225821E-17;
      const double s20 = +8.82090135625368160657E-17, s19 = -3.38459811878103047136E-16, s18 = -8.30608026366935789136E-16;
      const double s17 = +3.93397875437050071776E-15, s16 = +1.01765565969729044505E-14, s15 = -4.21128170307640802703E-14;
      const double s14 = -1.60818204519802480035E-13, s13 = +3.34714954175994481761E-13, s12 = +2.72600352129153073807E-12;
      const double s11 = +1.66894954752839083608E-12, s10 = -3.49278141024730899554E-11, s09 = -1.58580661666482709598E-10;
      const double s08 = -1.79289437183355633342E-10, s07 = +1.76281629144264523277E-9, s06 = +1.69050228879421288846E-8;
      const double s05 = +1.25391771228487041649E-7, s04 = +1.16229947068677338732E-6, s03 = +1.61038260117376323993E-5;
      const double s02 = +3.49810375601053973070E-4, s01 = +1.28478065259647610779E-2, s00 = +1.03665722588798326712E0;
      double b1 = s22, b0 = a*b1 + s21, b2 = a*b0 - b1 + s20;
      b1 = a*b2 - b0 + s19, b0 = a*b1 - b2 + s18;
      b2 = a*b0 - b1 + s17, b1 = a*b2 - b0 + s16, b0 = a*b1 - b2 + s15;
      b2 = a*b0 - b1 + s14, b1 = a*b2 - b0 + s13, b0 = a*b1 - b2 + s12;
      b2 = a*b0 - b1 + s11, b1 = a*b2 - b0 + s10, b0 = a*b1 - b2 + s09;
      b2 = a*b0 - b1 + s08, b1 = a*b2 - b0 + s07, b0 = a*b1 - b2 + s06;
      b2 = a*b0 - b1 + s05, b1 = a*b2 - b0 + s04, b0 = a*b1 - b2 + s03;
      b2 = a*b0 - b1 + s02, b1 = a*b2 - b0 + s01, b0 = a*b1 - b2 + s00;
      s = k * 0.5 * (b0 - b2);
      const double c23 = 8.06913408255155572081E-18, c22 = -2.08074168180148170312E-17, c21 = -5.98111329658272336816E-17;
      const double c20 = +2.68533951085945765591E-16, c19 = +4.52313941698904694774E-16, c18 = -3.10734917335299464535E-15;
      const double c17 = -4.42823207332531972288E-15, c16 = +3.49639695410806959872E-14, c15 = +6.63406731718911586609E-14;
      const double c14 = -3.71902448093119218395E-13, c13 = -1.27135418132338309016E-12, c12 = +2.74851141935315395333E-12;
      const double c11 = +2.33781843985453438400E-11, c10 = +2.71436006377612442764E-11, c09 = -2.56600180000355990529E-10;
      const double c08 = -1.61021375163803438552E-9, c07 = -4.72543064876271773512E-9, c06 = -3.00095178028681682282E-9;
      const double c05 = +7.79387474390914922337E-8, c04 = +1.06942765566401507066E-6, c03 = +1.59503164802313196374E-5;
      const double c02 = +3.49592575153777996871E-4, c01 = +1.28475387530065247392E-2, c00 = +1.03665693917934275131E0;
      b2 = c23, b1 = a*b2 + c22, b0 = a*b1 - b2 + c21;
      b2 = a*b0 - b1 + c20, b1 = a*b2 - b0 + c19, b0 = a*b1 - b2 + c18;
      b2 = a*b0 - b1 + c17, b1 = a*b2 - b0 + c16, b0 = a*b1 - b2 + c15;
      b2 = a*b0 - b1 + c14, b1 = a*b2 - b0 + c13, b0 = a*b1 - b2 + c12;
      b2 = a*b0 - b1 + c11, b1 = a*b2 - b0 + c10, b0 = a*b1 - b2 + c09;
      b2 = a*b0 - b1 + c08, b1 = a*b2 - b0 + c07, b0 = a*b1 - b2 + c06;
      b2 = a*b0 - b1 + c05, b1 = a*b2 - b0 + c04, b0 = a*b1 - b2 + c03;
      b2 = a*b0 - b1 + c02, b1 = a*b2 - b0 + c01, b0 = a*b1 - b2 + c00;
      c = Eul + log(x) + k * 0.5 * (b0 - b2);
   } else {
      c = s = maxrealnumber;
   }
   if (neg) s = -s;
   *shi = s;
   *chi = c;
}
} // end of namespace alglib_impl

namespace alglib {
void sinecosineintegrals(const double x, double &si, double &ci) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::sinecosineintegrals(x, &si, &ci);
   alglib_impl::ae_state_clear();
}

void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hyperbolicsinecosineintegrals(x, &shi, &chi);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === CHEBYSHEV Package ===
namespace alglib_impl {
// Calculation of the value of the Chebyshev polynomials of the
// first and second kinds.
//
// Parameters:
//     r   -   polynomial kind, either 1 or 2.
//     n   -   degree, n >= 0
//     x   -   argument, -1 <= x <= 1
//
// Result:
//     the value of the Chebyshev polynomial at x
// API: double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x);
double chebyshevcalculate(ae_int_t r, ae_int_t n, double x) {
   ae_int_t i;
   double a;
   double b;
   double result;
   result = 0.0;
// Prepare A and B
   if (r == 1) {
      a = 1.0;
      b = x;
   } else {
      a = 1.0;
      b = 2.0 * x;
   }
// Special cases: N = 0 or N = 1
   if (n == 0) {
      result = a;
      return result;
   }
   if (n == 1) {
      result = b;
      return result;
   }
// General case: N >= 2
   for (i = 2; i <= n; i++) {
      result = 2.0 * x * b - a;
      a = b;
      b = result;
   }
   return result;
}

// Summation of Chebyshev polynomials using Clenshaw's recurrence formula.
//
// This routine calculates
//     c[0]*T0(x) + c[1]*T1(x) + ... + c[N]*TN(x)
// or
//     c[0]*U0(x) + c[1]*U1(x) + ... + c[N]*UN(x)
// depending on the R.
//
// Parameters:
//     r   -   polynomial kind, either 1 or 2.
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Chebyshev polynomial at x
// API: double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x);
double chebyshevsum(RVector *c, ae_int_t r, ae_int_t n, double x) {
   double b1;
   double b2;
   ae_int_t i;
   double result;
   b1 = 0.0;
   b2 = 0.0;
   for (i = n; i >= 1; i--) {
      result = 2.0 * x * b1 - b2 + c->xR[i];
      b2 = b1;
      b1 = result;
   }
   if (r == 1) {
      result = -b2 + x * b1 + c->xR[0];
   } else {
      result = -b2 + 2.0 * x * b1 + c->xR[0];
   }
   return result;
}

// Representation of Tn as C[0] + C[1]*X + ... + C[N]*X^N
//
// Inputs:
//     N   -   polynomial degree, n >= 0
//
// Outputs:
//     C   -   coefficients
// API: void chebyshevcoefficients(const ae_int_t n, real_1d_array &c);
void chebyshevcoefficients(ae_int_t n, RVector *c) {
   ae_int_t i;
   SetVector(c);
   ae_vector_set_length(c, n + 1);
   for (i = 0; i <= n; i++) {
      c->xR[i] = 0.0;
   }
   if (n == 0 || n == 1) {
      c->xR[n] = 1.0;
   } else {
      c->xR[n] = exp((n - 1) * log(2.0));
      for (i = 0; i < n / 2; i++) {
         c->xR[n - 2 * (i + 1)] = -c->xR[n - 2 * i] * (n - 2 * i) * (n - 2 * i - 1) / 4.0 / (i + 1) / (n - i - 1);
      }
   }
}

// Conversion of a series of Chebyshev polynomials to a power series.
//
// Represents A[0]*T0(x) + A[1]*T1(x) + ... + A[N]*Tn(x) as
// B[0] + B[1]*X + ... + B[N]*X^N.
//
// Inputs:
//     A   -   Chebyshev series coefficients
//     N   -   degree, N >= 0
//
// Outputs:
//     B   -   power series coefficients
// API: void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b);
void fromchebyshev(RVector *a, ae_int_t n, RVector *b) {
//(#) The original did not have a check on n.
   ae_assert(n >= 0, "Domain error in fromchebyshev: n < 0.");
   ae_int_t i;
   ae_int_t k;
   double e;
   double d;
   SetVector(b);
   ae_vector_set_length(b, n + 1);
   for (i = 0; i <= n; i++) {
      b->xR[i] = 0.0;
   }
   d = 0.0;
   i = 0;
   do {
      k = i;
      do {
         e = b->xR[k];
         b->xR[k] = 0.0;
         if (i <= 1 && k == i) {
            b->xR[k] = 1.0;
         } else {
            if (i != 0) {
               b->xR[k] = 2.0 * d;
            }
            if (k > i + 1) {
               b->xR[k] -= b->xR[k - 2];
            }
         }
         d = e;
         k++;
      } while (k <= n);
      d = b->xR[i];
      e = 0.0;
      k = i;
      while (k <= n) {
         e += b->xR[k] * a->xR[k];
         k += 2;
      }
      b->xR[i] = e;
      i++;
   } while (i <= n);
}
} // end of namespace alglib_impl

namespace alglib {
double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::chebyshevcalculate(r, n, x);
   alglib_impl::ae_state_clear();
   return D;
}

double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::chebyshevsum(ConstT(ae_vector, c), r, n, x);
   alglib_impl::ae_state_clear();
   return D;
}

void chebyshevcoefficients(const ae_int_t n, real_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::chebyshevcoefficients(n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}

void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fromchebyshev(ConstT(ae_vector, a), n, ConstT(ae_vector, b));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === POISSONDISTR Package ===
// Depends on: IGAMMAF
namespace alglib_impl {
// Poisson distribution
//
// Returns the sum of the first k+1 terms of the Poisson
// distribution:
//
//   k         j
//   --   -m  m
//   >   e    --
//   --       j!
//  j=0
//
// The terms are not summed directly; instead the incomplete
// gamma integral is employed, according to the relation
//
// y = pdtr( k, m ) = igamc( k+1, m ).
//
// The arguments must both be positive.
//
// ACCURACY:
//
// See incomplete gamma function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double poissondistribution(const ae_int_t k, const double m);
double poissondistribution(ae_int_t k, double m) {
   double result;
   ae_assert(k >= 0 && m > 0.0, "Domain error in PoissonDistribution");
   result = incompletegammac(k + 1, m);
   return result;
}

// Complemented Poisson distribution
//
// Returns the sum of the terms k+1 to infinity of the Poisson
// distribution:
//
//  inf.       j
//   --   -m  m
//   >   e    --
//   --       j!
//  j=k+1
//
// The terms are not summed directly; instead the incomplete
// gamma integral is employed, according to the formula
//
// y = pdtrc( k, m ) = igam( k+1, m ).
//
// The arguments must both be positive.
//
// ACCURACY:
//
// See incomplete gamma function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double poissoncdistribution(const ae_int_t k, const double m);
double poissoncdistribution(ae_int_t k, double m) {
   double result;
   ae_assert(k >= 0 && m > 0.0, "Domain error in PoissonDistributionC");
   result = incompletegamma(k + 1, m);
   return result;
}

// Inverse Poisson distribution
//
// Finds the Poisson variable x such that the integral
// from 0 to x of the Poisson density is equal to the
// given probability y.
//
// This is accomplished using the inverse gamma integral
// function and the relation
//
//    m = igami( k+1, y ).
//
// ACCURACY:
//
// See inverse incomplete gamma function
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
// API: double invpoissondistribution(const ae_int_t k, const double y);
double invpoissondistribution(ae_int_t k, double y) {
   double result;
   ae_assert(k >= 0 && y >= 0.0 && y < 1.0, "Domain error in InvPoissonDistribution");
   result = invincompletegammac(k + 1, y);
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double poissondistribution(const ae_int_t k, const double m) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::poissondistribution(k, m);
   alglib_impl::ae_state_clear();
   return D;
}

double poissoncdistribution(const ae_int_t k, const double m) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::poissoncdistribution(k, m);
   alglib_impl::ae_state_clear();
   return D;
}

double invpoissondistribution(const ae_int_t k, const double y) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::invpoissondistribution(k, y);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === BETAF Package ===
// Depends on: GAMMAFUNC
namespace alglib_impl {
// Beta function
//
//                   -     -
//                  | (a) | (b)
// beta( a, b )  =  -----------.
//                     -
//                    | (a+b)
//
// For large arguments the logarithm of the function is
// evaluated using lgam(), then exponentiated.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE       0,30       30000       8.1e-14     1.1e-14
//
// Cephes Math Library Release 2.0:  April, 1987
// Copyright 1984, 1987 by Stephen L. Moshier
// API: double beta(const double a, const double b);
double beta(double a, double b) {
   double y;
   double sg;
   double s;
   double result;
   sg = 1.0;
   ae_assert(a > 0.0 || a != floor(a), "Overflow in Beta");
   ae_assert(b > 0.0 || b != floor(b), "Overflow in Beta");
   y = a + b;
   if (!SmallAtR(y, 171.624376956302725)) {
      y = lngamma(y, &s);
      sg *= s;
      y = lngamma(b, &s) - y;
      sg *= s;
      y += lngamma(a, &s);
      sg *= s;
      ae_assert(y <= log(maxrealnumber), "Overflow in Beta");
      result = sg * exp(y);
      return result;
   }
   y = gammafunction(y);
   ae_assert(y != 0.0, "Overflow in Beta");
   if (a > b) {
      y = gammafunction(a) / y;
      y *= gammafunction(b);
   } else {
      y = gammafunction(b) / y;
      y *= gammafunction(a);
   }
   result = y;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double beta(const double a, const double b) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::beta(a, b);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === FRESNEL Package ===
namespace alglib_impl {
// Fresnel integral
//
// Evaluates the Fresnel integrals
//
//           x
//           -
//          | |
// C(x) =   |   cos(pi/2 t**2) dt,
//        | |
//         -
//          0
//
//           x
//           -
//          | |
// S(x) =   |   sin(pi/2 t**2) dt.
//        | |
//         -
//          0
//
// The integrals are evaluated by a power series for x < 1.
// For x >= 1 auxiliary functions f(x) and g(x) are employed
// such that
//
// C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
// S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
//
// ACCURACY:
//
//  Relative error.
//
// Arithmetic  function   domain     # trials      peak         rms
//   IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
//   IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: void fresnelintegral(const double x, double &c, double &s);
void fresnelintegral(double x, double *c, double *s) {
   double xxa;
   double f;
   double g;
   double cc;
   double ss;
   double t;
   double u;
   double x2;
   double sn;
   double sd;
   double cn;
   double cd;
   double fn;
   double fd;
   double gn;
   double gd;
   xxa = x;
   x = fabs(xxa);
   x2 = x * x;
   if (x2 < 2.5625) {
      t = x2 * x2;
      sn = -2.99181919401019853726E3;
      sn = sn * t + 7.08840045257738576863E5;
      sn = sn * t - 6.29741486205862506537E7;
      sn = sn * t + 2.54890880573376359104E9;
      sn = sn * t - 4.42979518059697779103E10;
      sn = sn * t + 3.18016297876567817986E11;
      sd = 1.00000000000000000000E0;
      sd = sd * t + 2.81376268889994315696E2;
      sd = sd * t + 4.55847810806532581675E4;
      sd = sd * t + 5.17343888770096400730E6;
      sd = sd * t + 4.19320245898111231129E8;
      sd = sd * t + 2.24411795645340920940E10;
      sd = sd * t + 6.07366389490084639049E11;
      cn = -4.98843114573573548651E-8;
      cn = cn * t + 9.50428062829859605134E-6;
      cn = cn * t - 6.45191435683965050962E-4;
      cn = cn * t + 1.88843319396703850064E-2;
      cn = cn * t - 2.05525900955013891793E-1;
      cn = cn * t + 9.99999999999999998822E-1;
      cd = 3.99982968972495980367E-12;
      cd = cd * t + 9.15439215774657478799E-10;
      cd = cd * t + 1.25001862479598821474E-7;
      cd = cd * t + 1.22262789024179030997E-5;
      cd = cd * t + 8.68029542941784300606E-4;
      cd = cd * t + 4.12142090722199792936E-2;
      cd = cd * t + 1.00000000000000000118E0;
      *s = sign(xxa) * x * x2 * sn / sd;
      *c = sign(xxa) * x * cn / cd;
      return;
   }
   if (x > 36974.0) {
      *c = sign(xxa) * 0.5;
      *s = sign(xxa) * 0.5;
      return;
   }
   x2 = x * x;
   t = pi * x2;
   u = 1.0 / (t * t);
   t = 1.0 / t;
   fn = 4.21543555043677546506E-1;
   fn = fn * u + 1.43407919780758885261E-1;
   fn = fn * u + 1.15220955073585758835E-2;
   fn = fn * u + 3.45017939782574027900E-4;
   fn = fn * u + 4.63613749287867322088E-6;
   fn = fn * u + 3.05568983790257605827E-8;
   fn = fn * u + 1.02304514164907233465E-10;
   fn = fn * u + 1.72010743268161828879E-13;
   fn = fn * u + 1.34283276233062758925E-16;
   fn = fn * u + 3.76329711269987889006E-20;
   fd = 1.00000000000000000000E0;
   fd = fd * u + 7.51586398353378947175E-1;
   fd = fd * u + 1.16888925859191382142E-1;
   fd = fd * u + 6.44051526508858611005E-3;
   fd = fd * u + 1.55934409164153020873E-4;
   fd = fd * u + 1.84627567348930545870E-6;
   fd = fd * u + 1.12699224763999035261E-8;
   fd = fd * u + 3.60140029589371370404E-11;
   fd = fd * u + 5.88754533621578410010E-14;
   fd = fd * u + 4.52001434074129701496E-17;
   fd = fd * u + 1.25443237090011264384E-20;
   gn = 5.04442073643383265887E-1;
   gn = gn * u + 1.97102833525523411709E-1;
   gn = gn * u + 1.87648584092575249293E-2;
   gn = gn * u + 6.84079380915393090172E-4;
   gn = gn * u + 1.15138826111884280931E-5;
   gn = gn * u + 9.82852443688422223854E-8;
   gn = gn * u + 4.45344415861750144738E-10;
   gn = gn * u + 1.08268041139020870318E-12;
   gn = gn * u + 1.37555460633261799868E-15;
   gn = gn * u + 8.36354435630677421531E-19;
   gn = gn * u + 1.86958710162783235106E-22;
   gd = 1.00000000000000000000E0;
   gd = gd * u + 1.47495759925128324529E0;
   gd = gd * u + 3.37748989120019970451E-1;
   gd = gd * u + 2.53603741420338795122E-2;
   gd = gd * u + 8.14679107184306179049E-4;
   gd = gd * u + 1.27545075667729118702E-5;
   gd = gd * u + 1.04314589657571990585E-7;
   gd = gd * u + 4.60680728146520428211E-10;
   gd = gd * u + 1.10273215066240270757E-12;
   gd = gd * u + 1.38796531259578871258E-15;
   gd = gd * u + 8.39158816283118707363E-19;
   gd = gd * u + 1.86958710162783236342E-22;
   f = 1.0 - u * fn / fd;
   g = t * gn / gd;
   t = HalfPi * x2;
   cc = cos(t);
   ss = sin(t);
   t = pi * x;
   *c = 0.5 + (f * ss - g * cc) / t;
   *s = 0.5 - (f * cc + g * ss) / t;
   *c *= sign(xxa);
   *s *= sign(xxa);
}
} // end of namespace alglib_impl

namespace alglib {
void fresnelintegral(const double x, double &c, double &s) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::fresnelintegral(x, &c, &s);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === PSIF Package ===
namespace alglib_impl {
// Psi (digamma) function
//
//              d      -
//   psi(x)  =  -- ln | (x)
//              dx
//
// is the logarithmic derivative of the gamma function.
// For integer x,
//                   n-1
//                    -
// psi(n) = -EUL  +   >  1/k.
//                    -
//                   k=1
//
// This formula is used for 0 < n <= 10.  If x is negative, it
// is transformed to a positive argument by the reflection
// formula  psi(1-x) = psi(x) + pi cot(pi x).
// For general positive x, the argument is made greater than 10
// using the recurrence  psi(x+1) = psi(x) + 1/x.
// Then the following asymptotic expansion is applied:
//
//                           inf.   B
//                            -      2k
// psi(x) = log(x) - 1/2x -   >   -------
//                            -        2k
//                           k=1   2k x
//
// where the B2k are Bernoulli numbers.
//
// ACCURACY:
//    Relative error (except absolute when |psi| < 1):
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,30        30000       1.3e-15     1.4e-16
//    IEEE      -30,0       40000       1.5e-15     2.2e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
// API: double psi(const double x);
double psi(double x) {
   double p;
   double q;
   double nz;
   double s;
   double w;
   double y;
   double z;
   double polv;
   ae_int_t i;
   ae_int_t n;
   ae_int_t negative;
   double result;
   negative = 0;
   nz = 0.0;
   if (x <= 0.0) {
      negative = 1;
      q = x;
      p = floor(q);
      if (p == q) {
         ae_assert(false, "Singularity in Psi(x)");
         result = maxrealnumber;
         return result;
      }
      nz = q - p;
      if (nz != 0.5) {
         if (nz > 0.5) {
            p++;
            nz = q - p;
         }
         nz = pi / tan(pi * nz);
      } else {
         nz = 0.0;
      }
      x = 1.0 - x;
   }
   if (x <= 10.0 && x == floor(x)) {
      y = 0.0;
      n = ifloor(x);
      for (i = 1; i < n; i++) {
         w = i;
         y += 1.0 / w;
      }
      y -= Eul;
   } else {
      s = x;
      w = 0.0;
      while (s < 10.0) {
         w += 1.0 / s;
         s++;
      }
      if (s < 1.0E17) {
         z = 1.0 / (s * s);
         polv = 8.33333333333333333333E-2;
         polv = polv * z - 2.10927960927960927961E-2;
         polv = polv * z + 7.57575757575757575758E-3;
         polv = polv * z - 4.16666666666666666667E-3;
         polv = polv * z + 3.96825396825396825397E-3;
         polv = polv * z - 8.33333333333333333333E-3;
         polv = polv * z + 8.33333333333333333333E-2;
         y = z * polv;
      } else {
         y = 0.0;
      }
      y = log(s) - 0.5 / s - y - w;
   }
   if (negative != 0) {
      y -= nz;
   }
   result = y;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double psi(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::psi(x);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === AIRYF Package ===
namespace alglib_impl {
// Airy function
//
// Solution of the differential equation
//
// y"(x) = xy.
//
// The function returns the two independent solutions Ai, Bi
// and their first derivatives Ai'(x), Bi'(x).
//
// Evaluation is by power series summation for small x,
// by rational minimax approximations for large x.
//
// ACCURACY:
// Error criterion is absolute when function <= 1, relative
// when function > 1, except * denotes relative error criterion.
// For large negative x, the absolute error increases as x^1.5.
// For large positive x, the relative error increases as x^1.5.
//
// Arithmetic  domain   function  # trials      peak         rms
// IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
// IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
// IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
// IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
// IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
// IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: void airy(const double x, double &ai, double &aip, double &bi, double &bip);
void airy(double x, double *ai, double *aip, double *bi, double *bip) {
   double z;
   double zz;
   double t;
   double f;
   double g;
   double uf;
   double ug;
   double k;
   double zeta;
   double theta;
   ae_int_t domflg;
   double c1;
   double c2;
   double sqrt3;
   double sqpii;
   double afn;
   double afd;
   double agn;
   double agd;
   double apfn;
   double apfd;
   double apgn;
   double apgd;
   double an;
   double ad;
   double apn;
   double apd;
   double bn16;
   double bd16;
   double bppn;
   double bppd;
   *ai = 0.0;
   *aip = 0.0;
   *bi = 0.0;
   *bip = 0.0;
   sqpii = 5.64189583547756286948E-1;
   c1 = 0.35502805388781723926;
   c2 = 0.258819403792806798405;
   sqrt3 = 1.732050807568877293527;
   domflg = 0;
   if (x > 25.77) {
      *ai = 0.0;
      *aip = 0.0;
      *bi = maxrealnumber;
      *bip = maxrealnumber;
      return;
   }
   if (x < -2.09) {
      domflg = 15;
      t = sqrt(-x);
      zeta = -2.0 * x * t / 3.0;
      t = sqrt(t);
      k = sqpii / t;
      z = 1.0 / zeta;
      zz = z * z;
      afn = -1.31696323418331795333E-1;
      afn = afn * zz - 6.26456544431912369773E-1;
      afn = afn * zz - 6.93158036036933542233E-1;
      afn = afn * zz - 2.79779981545119124951E-1;
      afn = afn * zz - 4.91900132609500318020E-2;
      afn = afn * zz - 4.06265923594885404393E-3;
      afn = afn * zz - 1.59276496239262096340E-4;
      afn = afn * zz - 2.77649108155232920844E-6;
      afn = afn * zz - 1.67787698489114633780E-8;
      afd = 1.00000000000000000000E0;
      afd = afd * zz + 1.33560420706553243746E1;
      afd = afd * zz + 3.26825032795224613948E1;
      afd = afd * zz + 2.67367040941499554804E1;
      afd = afd * zz + 9.18707402907259625840E0;
      afd = afd * zz + 1.47529146771666414581E0;
      afd = afd * zz + 1.15687173795188044134E-1;
      afd = afd * zz + 4.40291641615211203805E-3;
      afd = afd * zz + 7.54720348287414296618E-5;
      afd = afd * zz + 4.51850092970580378464E-7;
      uf = 1.0 + zz * afn / afd;
      agn = 1.97339932091685679179E-2;
      agn = agn * zz + 3.91103029615688277255E-1;
      agn = agn * zz + 1.06579897599595591108E0;
      agn = agn * zz + 9.39169229816650230044E-1;
      agn = agn * zz + 3.51465656105547619242E-1;
      agn = agn * zz + 6.33888919628925490927E-2;
      agn = agn * zz + 5.85804113048388458567E-3;
      agn = agn * zz + 2.82851600836737019778E-4;
      agn = agn * zz + 6.98793669997260967291E-6;
      agn = agn * zz + 8.11789239554389293311E-8;
      agn = agn * zz + 3.41551784765923618484E-10;
      agd = 1.00000000000000000000E0;
      agd = agd * zz + 9.30892908077441974853E0;
      agd = agd * zz + 1.98352928718312140417E1;
      agd = agd * zz + 1.55646628932864612953E1;
      agd = agd * zz + 5.47686069422975497931E0;
      agd = agd * zz + 9.54293611618961883998E-1;
      agd = agd * zz + 8.64580826352392193095E-2;
      agd = agd * zz + 4.12656523824222607191E-3;
      agd = agd * zz + 1.01259085116509135510E-4;
      agd = agd * zz + 1.17166733214413521882E-6;
      agd = agd * zz + 4.91834570062930015649E-9;
      ug = z * agn / agd;
      theta = zeta + 0.25 * pi;
      f = sin(theta);
      g = cos(theta);
      *ai = k * (f * uf - g * ug);
      *bi = k * (g * uf + f * ug);
      apfn = 1.85365624022535566142E-1;
      apfn = apfn * zz + 8.86712188052584095637E-1;
      apfn = apfn * zz + 9.87391981747398547272E-1;
      apfn = apfn * zz + 4.01241082318003734092E-1;
      apfn = apfn * zz + 7.10304926289631174579E-2;
      apfn = apfn * zz + 5.90618657995661810071E-3;
      apfn = apfn * zz + 2.33051409401776799569E-4;
      apfn = apfn * zz + 4.08718778289035454598E-6;
      apfn = apfn * zz + 2.48379932900442457853E-8;
      apfd = 1.00000000000000000000E0;
      apfd = apfd * zz + 1.47345854687502542552E1;
      apfd = apfd * zz + 3.75423933435489594466E1;
      apfd = apfd * zz + 3.14657751203046424330E1;
      apfd = apfd * zz + 1.09969125207298778536E1;
      apfd = apfd * zz + 1.78885054766999417817E0;
      apfd = apfd * zz + 1.41733275753662636873E-1;
      apfd = apfd * zz + 5.44066067017226003627E-3;
      apfd = apfd * zz + 9.39421290654511171663E-5;
      apfd = apfd * zz + 5.65978713036027009243E-7;
      uf = 1.0 + zz * apfn / apfd;
      apgn = -3.55615429033082288335E-2;
      apgn = apgn * zz - 6.37311518129435504426E-1;
      apgn = apgn * zz - 1.70856738884312371053E0;
      apgn = apgn * zz - 1.50221872117316635393E0;
      apgn = apgn * zz - 5.63606665822102676611E-1;
      apgn = apgn * zz - 1.02101031120216891789E-1;
      apgn = apgn * zz - 9.48396695961445269093E-3;
      apgn = apgn * zz - 4.60325307486780994357E-4;
      apgn = apgn * zz - 1.14300836484517375919E-5;
      apgn = apgn * zz - 1.33415518685547420648E-7;
      apgn = apgn * zz - 5.63803833958893494476E-10;
      apgd = 1.00000000000000000000E0;
      apgd = apgd * zz + 9.85865801696130355144E0;
      apgd = apgd * zz + 2.16401867356585941885E1;
      apgd = apgd * zz + 1.73130776389749389525E1;
      apgd = apgd * zz + 6.17872175280828766327E0;
      apgd = apgd * zz + 1.08848694396321495475E0;
      apgd = apgd * zz + 9.95005543440888479402E-2;
      apgd = apgd * zz + 4.78468199683886610842E-3;
      apgd = apgd * zz + 1.18159633322838625562E-4;
      apgd = apgd * zz + 1.37480673554219441465E-6;
      apgd = apgd * zz + 5.79912514929147598821E-9;
      ug = z * apgn / apgd;
      k = sqpii * t;
      *aip = -k * (g * uf + f * ug);
      *bip = k * (f * uf - g * ug);
      return;
   }
   if (x >= 2.09) {
      domflg = 5;
      t = sqrt(x);
      zeta = 2.0 * x * t / 3.0;
      g = exp(zeta);
      t = sqrt(t);
      k = 2.0 * t * g;
      z = 1.0 / zeta;
      an = 3.46538101525629032477E-1;
      an = an * z + 1.20075952739645805542E1;
      an = an * z + 7.62796053615234516538E1;
      an = an * z + 1.68089224934630576269E2;
      an = an * z + 1.59756391350164413639E2;
      an = an * z + 7.05360906840444183113E1;
      an = an * z + 1.40264691163389668864E1;
      an = an * z + 9.99999999999999995305E-1;
      ad = 5.67594532638770212846E-1;
      ad = ad * z + 1.47562562584847203173E1;
      ad = ad * z + 8.45138970141474626562E1;
      ad = ad * z + 1.77318088145400459522E2;
      ad = ad * z + 1.64234692871529701831E2;
      ad = ad * z + 7.14778400825575695274E1;
      ad = ad * z + 1.40959135607834029598E1;
      ad = ad * z + 1.00000000000000000470E0;
      f = an / ad;
      *ai = sqpii * f / k;
      k = -0.5 * sqpii * t / g;
      apn = 6.13759184814035759225E-1;
      apn = apn * z + 1.47454670787755323881E1;
      apn = apn * z + 8.20584123476060982430E1;
      apn = apn * z + 1.71184781360976385540E2;
      apn = apn * z + 1.59317847137141783523E2;
      apn = apn * z + 6.99778599330103016170E1;
      apn = apn * z + 1.39470856980481566958E1;
      apn = apn * z + 1.00000000000000000550E0;
      apd = 3.34203677749736953049E-1;
      apd = apd * z + 1.11810297306158156705E1;
      apd = apd * z + 7.11727352147859965283E1;
      apd = apd * z + 1.58778084372838313640E2;
      apd = apd * z + 1.53206427475809220834E2;
      apd = apd * z + 6.86752304592780337944E1;
      apd = apd * z + 1.38498634758259442477E1;
      apd = apd * z + 9.99999999999999994502E-1;
      f = apn / apd;
      *aip = f * k;
      if (x > 8.3203353) {
         bn16 = -2.53240795869364152689E-1;
         bn16 = bn16 * z + 5.75285167332467384228E-1;
         bn16 = bn16 * z - 3.29907036873225371650E-1;
         bn16 = bn16 * z + 6.44404068948199951727E-2;
         bn16 = bn16 * z - 3.82519546641336734394E-3;
         bd16 = 1.00000000000000000000E0;
         bd16 = bd16 * z - 7.15685095054035237902E0;
         bd16 = bd16 * z + 1.06039580715664694291E1;
         bd16 = bd16 * z - 5.23246636471251500874E0;
         bd16 = bd16 * z + 9.57395864378383833152E-1;
         bd16 = bd16 * z - 5.50828147163549611107E-2;
         f = z * bn16 / bd16;
         k = sqpii * g;
         *bi = k * (1.0 + f) / t;
         bppn = 4.65461162774651610328E-1;
         bppn = bppn * z - 1.08992173800493920734E0;
         bppn = bppn * z + 6.38800117371827987759E-1;
         bppn = bppn * z - 1.26844349553102907034E-1;
         bppn = bppn * z + 7.62487844342109852105E-3;
         bppd = 1.00000000000000000000E0;
         bppd = bppd * z - 8.70622787633159124240E0;
         bppd = bppd * z + 1.38993162704553213172E1;
         bppd = bppd * z - 7.14116144616431159572E0;
         bppd = bppd * z + 1.34008595960680518666E0;
         bppd = bppd * z - 7.84273211323341930448E-2;
         f = z * bppn / bppd;
         *bip = k * t * (1.0 + f);
         return;
      }
   }
   f = 1.0;
   g = x;
   t = 1.0;
   uf = 1.0;
   ug = x;
   k = 1.0;
   z = x * x * x;
   while (t > machineepsilon) {
      uf *= z;
      k++;
      uf /= k;
      ug *= z;
      k++;
      ug /= k;
      uf /= k;
      f += uf;
      k++;
      ug /= k;
      g += ug;
      t = fabs(uf / f);
   }
   uf = c1 * f;
   ug = c2 * g;
   if (domflg % 2 == 0) {
      *ai = uf - ug;
   }
   if (domflg / 2 % 2 == 0) {
      *bi = sqrt3 * (uf + ug);
   }
   k = 4.0;
   uf = x * x / 2.0;
   ug = z / 3.0;
   f = uf;
   g = 1.0 + ug;
   uf /= 3.0;
   t = 1.0;
   while (t > machineepsilon) {
      uf *= z;
      ug /= k;
      k++;
      ug *= z;
      uf /= k;
      f += uf;
      k++;
      ug /= k;
      uf /= k;
      g += ug;
      k++;
      t = fabs(ug / g);
   }
   uf = c1 * f;
   ug = c2 * g;
   if (domflg / 4 % 2 == 0) {
      *aip = uf - ug;
   }
   if (domflg / 8 % 2 == 0) {
      *bip = sqrt3 * (uf + ug);
   }
}
} // end of namespace alglib_impl

namespace alglib {
void airy(const double x, double &ai, double &aip, double &bi, double &bip) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::airy(x, &ai, &aip, &bi, &bip);
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === DAWSON Package ===
namespace alglib_impl {
// Dawson's Integral
//
// Approximates the integral
//
//                             x
//                             -
//                      2     | |        2
//  dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
//                          | |
//                           -
//                           0
//
// Three different rational approximations are employed, for
// the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,10        10000       6.9e-16     1.0e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: double dawsonintegral(const double x);
double dawsonintegral(double x) {
   double x2;
   double y;
   ae_int_t sg;
   double an;
   double ad;
   double bn;
   double bd;
   double cn;
   double cd;
   double result;
   sg = 1;
   if (x < 0.0) {
      sg = -1;
      x = -x;
   }
   if (x < 3.25) {
      x2 = x * x;
      an = 1.13681498971755972054E-11;
      an = an * x2 + 8.49262267667473811108E-10;
      an = an * x2 + 1.94434204175553054283E-8;
      an = an * x2 + 9.53151741254484363489E-7;
      an = an * x2 + 3.07828309874913200438E-6;
      an = an * x2 + 3.52513368520288738649E-4;
      an = an * x2 + (-8.50149846724410912031E-4);
      an = an * x2 + 4.22618223005546594270E-2;
      an = an * x2 + (-9.17480371773452345351E-2);
      an = an * x2 + 9.99999999999999994612E-1;
      ad = 2.40372073066762605484E-11;
      ad = ad * x2 + 1.48864681368493396752E-9;
      ad = ad * x2 + 5.21265281010541664570E-8;
      ad = ad * x2 + 1.27258478273186970203E-6;
      ad = ad * x2 + 2.32490249820789513991E-5;
      ad = ad * x2 + 3.25524741826057911661E-4;
      ad = ad * x2 + 3.48805814657162590916E-3;
      ad = ad * x2 + 2.79448531198828973716E-2;
      ad = ad * x2 + 1.58874241960120565368E-1;
      ad = ad * x2 + 5.74918629489320327824E-1;
      ad = ad * x2 + 1.00000000000000000539E0;
      y = x * an / ad;
      result = sg * y;
      return result;
   }
   x2 = 1.0 / (x * x);
   if (x < 6.25) {
      bn = 5.08955156417900903354E-1;
      bn = bn * x2 - 2.44754418142697847934E-1;
      bn = bn * x2 + 9.41512335303534411857E-2;
      bn = bn * x2 - 2.18711255142039025206E-2;
      bn = bn * x2 + 3.66207612329569181322E-3;
      bn = bn * x2 - 4.23209114460388756528E-4;
      bn = bn * x2 + 3.59641304793896631888E-5;
      bn = bn * x2 - 2.14640351719968974225E-6;
      bn = bn * x2 + 9.10010780076391431042E-8;
      bn = bn * x2 - 2.40274520828250956942E-9;
      bn = bn * x2 + 3.59233385440928410398E-11;
      bd = 1.00000000000000000000E0;
      bd = bd * x2 - 6.31839869873368190192E-1;
      bd = bd * x2 + 2.36706788228248691528E-1;
      bd = bd * x2 - 5.31806367003223277662E-2;
      bd = bd * x2 + 8.48041718586295374409E-3;
      bd = bd * x2 - 9.47996768486665330168E-4;
      bd = bd * x2 + 7.81025592944552338085E-5;
      bd = bd * x2 - 4.55875153252442634831E-6;
      bd = bd * x2 + 1.89100358111421846170E-7;
      bd = bd * x2 - 4.91324691331920606875E-9;
      bd = bd * x2 + 7.18466403235734541950E-11;
      y = 1.0 / x + x2 * bn / (bd * x);
      result = sg * 0.5 * y;
      return result;
   }
   if (x > 1.0E9) {
      result = sg * 0.5 / x;
      return result;
   }
   cn = -5.90592860534773254987E-1;
   cn = cn * x2 + 6.29235242724368800674E-1;
   cn = cn * x2 - 1.72858975380388136411E-1;
   cn = cn * x2 + 1.64837047825189632310E-2;
   cn = cn * x2 - 4.86827613020462700845E-4;
   cd = 1.00000000000000000000E0;
   cd = cd * x2 - 2.69820057197544900361E0;
   cd = cd * x2 + 1.73270799045947845857E0;
   cd = cd * x2 - 3.93708582281939493482E-1;
   cd = cd * x2 + 3.44278924041233391079E-2;
   cd = cd * x2 - 9.73655226040941223894E-4;
   y = 1.0 / x + x2 * cn / (cd * x);
   result = sg * 0.5 * y;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double dawsonintegral(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::dawsonintegral(x);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === HERMITE Package ===
namespace alglib_impl {
// Calculation of the value of the Hermite polynomial.
//
// Parameters:
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Hermite polynomial Hn at x
// API: double hermitecalculate(const ae_int_t n, const double x);
double hermitecalculate(ae_int_t n, double x) {
   ae_int_t i;
   double a;
   double b;
   double result;
   result = 0.0;
// Prepare A and B
   a = 1.0;
   b = 2.0 * x;
// Special cases: N = 0 or N = 1
   if (n == 0) {
      result = a;
      return result;
   }
   if (n == 1) {
      result = b;
      return result;
   }
// General case: N >= 2
   for (i = 2; i <= n; i++) {
      result = 2.0 * x * b - 2.0 * (i - 1) * a;
      a = b;
      b = result;
   }
   return result;
}

// Summation of Hermite polynomials using Clenshaw's recurrence formula.
//
// This routine calculates
//     c[0]*H0(x) + c[1]*H1(x) + ... + c[N]*HN(x)
//
// Parameters:
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Hermite polynomial at x
// API: double hermitesum(const real_1d_array &c, const ae_int_t n, const double x);
double hermitesum(RVector *c, ae_int_t n, double x) {
   double b1;
   double b2;
   ae_int_t i;
   double result;
   b1 = 0.0;
   b2 = 0.0;
   result = 0.0;
   for (i = n; i >= 0; i--) {
      result = 2.0 * (x * b1 - (i + 1) * b2) + c->xR[i];
      b2 = b1;
      b1 = result;
   }
   return result;
}

// Representation of Hn as C[0] + C[1]*X + ... + C[N]*X^N
//
// Inputs:
//     N   -   polynomial degree, n >= 0
//
// Outputs:
//     C   -   coefficients
// API: void hermitecoefficients(const ae_int_t n, real_1d_array &c);
void hermitecoefficients(ae_int_t n, RVector *c) {
   ae_int_t i;
   SetVector(c);
   ae_vector_set_length(c, n + 1);
   for (i = 0; i <= n; i++) {
      c->xR[i] = 0.0;
   }
   c->xR[n] = exp(n * log(2.0));
   for (i = 0; i < n / 2; i++) {
      c->xR[n - 2 * (i + 1)] = -c->xR[n - 2 * i] * (n - 2 * i) * (n - 2 * i - 1) / 4.0 / (i + 1);
   }
}
} // end of namespace alglib_impl

namespace alglib {
double hermitecalculate(const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hermitecalculate(n, x);
   alglib_impl::ae_state_clear();
   return D;
}

double hermitesum(const real_1d_array &c, const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::hermitesum(ConstT(ae_vector, c), n, x);
   alglib_impl::ae_state_clear();
   return D;
}

void hermitecoefficients(const ae_int_t n, real_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::hermitecoefficients(n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === LEGENDRE Package ===
namespace alglib_impl {
// Calculation of the value of the Legendre polynomial Pn.
//
// Parameters:
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Legendre polynomial Pn at x
// API: double legendrecalculate(const ae_int_t n, const double x);
double legendrecalculate(ae_int_t n, double x) {
   double a;
   double b;
   ae_int_t i;
   double result;
   result = 1.0;
   a = 1.0;
   b = x;
   if (n == 0) {
      result = a;
      return result;
   }
   if (n == 1) {
      result = b;
      return result;
   }
   for (i = 2; i <= n; i++) {
      result = ((2 * i - 1) * x * b - (i - 1) * a) / i;
      a = b;
      b = result;
   }
   return result;
}

// Summation of Legendre polynomials using Clenshaw's recurrence formula.
//
// This routine calculates
//     c[0]*P0(x) + c[1]*P1(x) + ... + c[N]*PN(x)
//
// Parameters:
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Legendre polynomial at x
// API: double legendresum(const real_1d_array &c, const ae_int_t n, const double x);
double legendresum(RVector *c, ae_int_t n, double x) {
   double b1;
   double b2;
   ae_int_t i;
   double result;
   b1 = 0.0;
   b2 = 0.0;
   result = 0.0;
   for (i = n; i >= 0; i--) {
      result = (2 * i + 1) * x * b1 / (i + 1) - (i + 1) * b2 / (i + 2) + c->xR[i];
      b2 = b1;
      b1 = result;
   }
   return result;
}

// Representation of Pn as C[0] + C[1]*X + ... + C[N]*X^N
//
// Inputs:
//     N   -   polynomial degree, n >= 0
//
// Outputs:
//     C   -   coefficients
// API: void legendrecoefficients(const ae_int_t n, real_1d_array &c);
void legendrecoefficients(ae_int_t n, RVector *c) {
   ae_int_t i;
   SetVector(c);
   ae_vector_set_length(c, n + 1);
   for (i = 0; i <= n; i++) {
      c->xR[i] = 0.0;
   }
   c->xR[n] = 1.0;
   for (i = 1; i <= n; i++) {
      c->xR[n] = c->xR[n] * (n + i) / 2.0 / i;
   }
   for (i = 0; i < n / 2; i++) {
      c->xR[n - 2 * (i + 1)] = -c->xR[n - 2 * i] * (n - 2 * i) * (n - 2 * i - 1) / 2.0 / (i + 1) / (2 * (n - i) - 1);
   }
}
} // end of namespace alglib_impl

namespace alglib {
double legendrecalculate(const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::legendrecalculate(n, x);
   alglib_impl::ae_state_clear();
   return D;
}

double legendresum(const real_1d_array &c, const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::legendresum(ConstT(ae_vector, c), n, x);
   alglib_impl::ae_state_clear();
   return D;
}

void legendrecoefficients(const ae_int_t n, real_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::legendrecoefficients(n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === BESSEL Package ===
namespace alglib_impl {
static void bessel_besselasympt0(double w, double *pzero, double *qzero) {
   double ww = w * w;
   const double pp02 = 2485.271928957404011288128951, pp04 = 153982.6532623911470917825993, pp06 = 2016135.283049983642487182349;
   const double pp08 = 8413041.456550439208464315611, pp10 = 12332384.76817638145232406055, pp12 = 5393485.083869438325262122897;
   double pp = pp12 + ww * (pp10 + ww * (pp08 + ww * (pp06 + ww * (pp04 + ww * pp02))));
   const double pq02 = 2615.700736920839685159081813, pq04 = 156001.7276940030940592769933, pq06 = 2025066.801570134013891035236;
   const double pq08 = 8426449.050629797331554404810, pq10 = 12338310.22786324960844856182, pq12 = 5393485.083869438325560444960;
   double pq = pq12 + ww * (pq10 + ww * (pq08 + ww * (pq06 + ww * (pq04 + ww * (pq02 + ww)))));
   const double qp01 = -4.887199395841261531199129300, qp03 = -226.2630641933704113967255053, qp05 = -2365.956170779108192723612816;
   const double qp07 = -8239.066313485606568803548860, qp09 = -10381.41698748464093880530341, qp11 = -3984.617357595222463506790588;
   double qp = w * (qp11 + ww * (qp09 + ww * (qp07 + ww * (qp05 + ww * (qp03 + ww * (qp01))))));
   const double qq02 = 408.7714673983499223402830260, qq04 = 15704.89191515395519392882766, qq06 = 156021.3206679291652539287109;
   const double qq08 = 533291.3634216897168722255057, qq10 = 666745.4239319826986004038103, qq12 = 255015.5108860942382983170882;
   double qq = qq12 + ww * (qq10 + ww * (qq08 + ww * (qq06 + ww * (qq04 + ww * (qq02 + ww)))));
   *pzero = pp / pq;
   *qzero = qp / qq;
}

static void bessel_besselasympt1(double w, double *pzero, double *qzero) {
   double ww = w * w;
   const double pp02 = -1611.616644324610116477412898, pp04 = -109824.0554345934672737413139, pp06 = -1523529.351181137383255105722;
   const double pp08 = -6603373.248364939109255245434, pp10 = -9942246.505077641195658377899, pp12 = -4435757.816794127857114720794;
   double pp = pp12 + ww * (pp10 + ww * (pp08 + ww * (pp06 + ww * (pp04 + ww * pp02))));
   const double pq02 = -1455.009440190496182453565068, pq04 = -107263.8599110382011903063867, pq06 = -1511809.506634160881644546358;
   const double pq08 = -6585339.479723087072826915069, pq10 = -9934124.389934585658967556309, pq12 = -4435757.816794127856828016962;
   double pq = pq12 + ww * (pq10 + ww * (pq08 + ww * (pq06 + ww * (pq04 + ww * (pq02 + ww)))));
   const double qp01 = 35.26513384663603218592175580, qp03 = 1706.375429020768002061283546, qp05 = 18494.26287322386679652009819;
   const double qp07 = 66178.83658127083517939992166, qp09 = 85145.16067533570196555001171, qp11 = 33220.91340985722351859704442;
   double qp = w * (qp11 + ww * (qp09 + ww * (qp07 + ww * (qp05 + ww * (qp03 + ww * qp01)))));
   const double qq02 = 863.8367769604990967475517183, qq04 = 37890.22974577220264142952256, qq06 = 400294.4358226697511708610813;
   const double qq08 = 1419460.669603720892855755253, qq10 = 1819458.042243997298924553839, qq12 = 708712.8194102874357377502472;
   double qq = qq12 + ww * (qq10 + ww * (qq08 + ww * (qq06 + ww * (qq04 + ww * (qq02 + ww)))));
   *pzero = pp / pq;
   *qzero = qp / qq;
}

// Bessel function of order zero
//
// Returns Bessel function of order zero of the argument.
//
// The domain is divided into the intervals [0, 5] and
// (5, infinity). In the first interval the following rational
// approximation is used:
//
//        2         2
// (w - r  ) (w - r  ) P (w) / Q (w)
//       1         2    3       8
//
//            2
// where w = x  and the two r's are zeros of the function.
//
// In the second interval, the Hankel asymptotic expansion
// is employed with two rational functions of degree 6/6
// and 7/7.
//
// ACCURACY:
//
//                      Absolute error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       60000       4.2e-16     1.1e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: double besselj0(const double x);
double besselj0(double x) {
   if (x < 0.0) x = -x;
   if (x > 8.0) {
      double pzero, qzero;
      bessel_besselasympt0(8.0 / x, &pzero, &qzero);
      double nn = x - pi / 4.0;
      return sqrt(2.0 / pi / x) * (pzero * cos(nn) - qzero * sin(nn));
   }
   double xsq = x * x;
   const double p16 = +26857.86856980014981415848441, p14 = -40504123.71833132706360663322, p12 = +25071582855.36881945555156435;
   const double p10 = -8085222034853.793871199468171, p08 = +1434354939140344.111664316553, p06 = -136762035308817138.6865416609;
   const double p04 = +6382059341072356562.289432465, p02 = -117915762910761053603.8440800, p00 = +493378725179413356181.6813446;
   double p1 = p00 + xsq * (p02 + xsq * (p04 + xsq * (p06 + xsq * (p08 + xsq * (p10 + xsq * (p12 + xsq * (p14 + xsq * p16)))))));
   const double q14 = 1363.063652328970604442810507, q12 = 1114636.098462985378182402543;
   const double q10 = 669998767.2982239671814028660, q08 = 312304311494.1213172572469442, q06 = 112775673967979.8507056031594;
   const double q04 = 30246356167094626.98627330784, q02 = 5428918384092285160.200195092, q00 = 493378725179413356211.3278438;
   double q1 = q00 + xsq * (q02 + xsq * (q04 + xsq * (q06 + xsq * (q08 + xsq * (q10 + xsq * (q12 + xsq * (q14 + xsq)))))));
   return p1 / q1;
}

// Bessel function of order one
//
// Returns Bessel function of order one of the argument.
//
// The domain is divided into the intervals [0, 8] and
// (8, infinity). In the first interval a 24 term Chebyshev
// expansion is used. In the second, the asymptotic
// trigonometric representation is employed using two
// rational functions of degree 5/5.
//
// ACCURACY:
//
//                      Absolute error:
// arithmetic   domain      # trials      peak         rms
//    IEEE      0, 30       30000       2.6e-16     1.1e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: double besselj1(const double x);
double besselj1(double x) {
   bool neg = x < 0.0;
   if (neg) x = -x;
   if (x > 8.0) {
      double pzero, qzero;
      bessel_besselasympt1(8.0 / x, &pzero, &qzero);
      double nn = x - 3.0 * pi / 4.0;
      double result = sqrt(2.0 / pi / x) * (pzero * cos(nn) - qzero * sin(nn));
      return neg ? -result : +result;
   }
   double xsq = x * x;
   const double p17 = +2701.122710892323414856790990, p15 = -4695753.530642995859767162166, p13 = +3413234182.301700539091292655;
   const double p11 = -1322983480332.126453125473247, p09 = +290879526383477.5409737601689, p07 = -35888175699101060.50743641413;
   const double p05 = +2316433580634002297.931815435, p03 = -66721065689249162980.20941484, p01 = +581199354001606143928.050809;
   double p1 = x * (p01 + xsq * (p03 + xsq * (p05 + xsq * (p07 + xsq * (p09 + xsq * (p11 + xsq * (p13 + xsq * (p15 + xsq * p17))))))));
   const double q14 = 1606.931573481487801970916749, q12 = 1501793.594998585505921097578;
   const double q10 = 1013863514.358673989967045588, q08 = 524371026216.7649715406728642, q06 = 208166122130760.7351240184229;
   const double q04 = 60920613989175217.46105196863, q02 = 11857707121903209998.37113348, q00 = 1162398708003212287858.529400;
   double q1 = q00 + xsq * (q02 + xsq * (q04 + xsq * (q06 + xsq * (q08 + xsq * (q10 + xsq * (q12 + xsq * (q14 + xsq)))))));
   return neg ? -p1 / q1 : +p1 / q1;
}

// Bessel function of integer order
//
// Returns Bessel function of order n, where n is a
// (possibly negative) integer.
//
// The ratio of jn(x) to j0(x) is computed by backward
// recurrence.  First the ratio jn/jn-1 is found by a
// continued fraction expansion.  Then the recurrence
// relating successive orders is applied until j0 or j1 is
// reached.
//
// If n = 0 or 1 the routine for j0 or j1 is called
// directly.
//
// ACCURACY:
//
//                      Absolute error:
// arithmetic   range      # trials      peak         rms
//    IEEE      0, 30        5000       4.4e-16     7.9e-17
//
// Not suitable for large n or x. Use jv() (fractional order) instead.
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double besseljn(const ae_int_t n, const double x);
double besseljn(ae_int_t n, double x) {
   ae_int_t sg;
   if (n < 0) {
      n = -n;
      sg = n % 2 == 0 ? +1 : -1;
   } else sg = +1;
   if (x < 0.0) {
      if (n % 2 != 0) sg = -sg;
      x = -x;
   }
   if (n == 0) return sg * besselj0(x);
   else if (n == 1) return sg * besselj1(x);
   else if (n == 2) return x == 0.0 ? 0.0 : sg * (2.0 * besselj1(x) / x - besselj0(x));
   else if (x < machineepsilon) return 0.0;
   ae_int_t k = 53;
   double pk = 2.0 * (n + k);
   double ans = pk;
   double xk = x * x;
   do {
      pk -= 2.0;
      ans = pk - xk / ans;
      k--;
   } while (k != 0);
   ans = x / ans;
   pk = 1.0;
   double pkm1 = 1.0 / ans;
   k = n - 1;
   double r = 2.0 * k;
   do {
      double pkm2 = (pkm1 * r - pk * x) / x;
      pk = pkm1;
      pkm1 = pkm2;
      r -= 2.0;
   } while (--k > 0);
   return sg * (fabs(pk) > fabs(pkm1) ? besselj1(x) / pk : besselj0(x) / pkm1);
}

// Bessel function of the second kind, order zero
//
// Returns Bessel function of the second kind, of order
// zero, of the argument.
//
// The domain is divided into the intervals [0, 5] and
// (5, infinity). In the first interval a rational approximation
// R(x) is employed to compute
//   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
// Thus a call to j0() is required.
//
// In the second interval, the Hankel asymptotic expansion
// is employed with two rational functions of degree 6/6
// and 7/7.
//
// ACCURACY:
//
//  Absolute error, when y0(x) < 1; else relative error:
//
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       30000       1.3e-15     1.6e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: double bessely0(const double x);
double bessely0(double x) {
   if (x > 8.0) {
      double pzero, qzero;
      bessel_besselasympt0(8.0 / x, &pzero, &qzero);
      double nn = x - pi / 4.0;
      return sqrt(2.0 / pi / x) * (pzero * sin(nn) + qzero * cos(nn));
   }
   double xsq = sqr(x);
   const double p16 = -41370.35497933148554125235152, p14 = +59152134.65686889654273830069, p12 = -34363712229.79040378171030138;
   const double p10 = +10255208596863.94284509167421, p08 = -1648605817185729.473122082537, p06 = +137562431639934407.8571335453;
   const double p04 = -5247065581112764941.297350814, p02 = +65874732757195549259.99402049, p00 = -27502866786291095837.01933175;
   double p4 = p00 + xsq * (p02 + xsq * (p04 + xsq * (p06 + xsq * (p08 + xsq * (p10 + xsq * (p12 + xsq * (p14 + xsq * p16)))))));
   const double q16 = 1282.452772478993804176329391, q14 = 1001702.641288906265666651753;
   const double q12 = 579512264.0700729537480087915, q10 = 261306575504.1081249568482092, q08 = 91620380340751.85262489147968;
   const double q06 = 23928830434997818.57439356652, q04 = 4192417043410839973.904769661, q02 = 372645883898616588198.9980;
   double q4 = q02 + xsq * (q04 + xsq * (q06 + xsq * (q08 + xsq * (q10 + xsq * (q12 + xsq * (q14 + xsq * (q16 + xsq)))))));
   return p4 / q4 + 2.0 / pi * besselj0(x) * log(x);
}

// Bessel function of second kind of order one
//
// Returns Bessel function of the second kind of order one
// of the argument.
//
// The domain is divided into the intervals [0, 8] and
// (8, infinity). In the first interval a 25 term Chebyshev
// expansion is used, and a call to j1() is required.
// In the second, the asymptotic trigonometric representation
// is employed using two rational functions of degree 5/5.
//
// ACCURACY:
//
//                      Absolute error:
// arithmetic   domain      # trials      peak         rms
//    IEEE      0, 30       30000       1.0e-15     1.3e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: double bessely1(const double x);
double bessely1(double x) {
   if (x > 8.0) {
      double pzero, qzero;
      bessel_besselasympt1(8.0 / x, &pzero, &qzero);
      double nn = x - 3.0 * pi / 4.0;
      return sqrt(2.0 / pi / x) * (pzero * sin(nn) + qzero * cos(nn));
   }
   double xsq = x * x;
   const double p17 = -2108847.540133123652824139923, p15 = +3639488548.124002058278999428, p13 = -2580681702194.450950541426399;
   const double p11 = +956993023992168.3481121552788, p09 = -196588746272214065.8820322248, p07 = +21931073399177975921.11427556;
   const double p05 = -1212297555414509577913.561535, p03 = +26554738314348543268942.48968, p01 = -99637534243069222259967.44354;
   double p4 = x * (p01 + xsq * (p03 + xsq * (p05 + xsq * (p07 + xsq * (p09 + xsq * (p11 + xsq * (p13 + xsq * (p15 + xsq * p17))))))));
   const double q16 = 1612.361029677000859332072312, q14 = 1563282.754899580604737366452, q12 = 1128686837.169442121732366891;
   const double q10 = 646534088126.5275571961681500, q08 = 297663212564727.6729292742282, q06 = 108225825940881955.2553850180;
   const double q04 = 29549879358971486742.90758119, q02 = 5435310377188854170800.653097, q00 = 508206736694124324531442.4152;
   double q4 = q00 + xsq * (q02 + xsq * (q04 + xsq * (q06 + xsq * (q08 + xsq * (q10 + xsq * (q12 + xsq * (q14 + xsq * (q16 + xsq))))))));
   return p4 / q4 + 2.0 / pi * (besselj1(x) * log(x) - 1.0 / x);
}

// Bessel function of second kind of integer order
//
// Returns Bessel function of order n, where n is a
// (possibly negative) integer.
//
// The function is evaluated by forward recurrence on
// n, starting with values computed by the routines
// y0() and y1().
//
// If n = 0 or 1 the routine for y0 or y1 is called
// directly.
//
// ACCURACY:
//                      Absolute error, except relative
//                      when y > 1:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       30000       3.4e-15     4.3e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double besselyn(const ae_int_t n, const double x);
double besselyn(ae_int_t n, double x) {
   double s = 1.0;
   if (n < 0) {
      n = -n;
      if (n % 2 != 0) s = -1.0;
   }
   if (n == 0) return bessely0(x);
   else if (n == 1) return s * bessely1(x);
   double a = bessely0(x), b = bessely1(x);
   for (ae_int_t i = 1; i < n; i++) {
      double tmp = b;
      b = 2.0 * i / x * b - a;
      a = tmp;
   }
   return s * b;
}

// Modified Bessel function of order zero
//
// Returns modified Bessel function of order zero of the
// argument.
//
// The function is defined as i0(x) = j0( ix ).
//
// The range is partitioned into the two intervals [0,8] and
// (8, infinity).  Chebyshev polynomial expansions are employed
// in each interval.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,30        30000       5.8e-16     1.4e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double besseli0(const double x);
double besseli0(double x) {
   if (x < 0.0) x = -x;
   if (x <= 8.0) {
      double y = x / 2.0 - 2.0;
      const double b29 = -4.41534164647933937950E-18, b28 = +3.33079451882223809783E-17, b27 = -2.43127984654795469359E-16;
      const double b26 = +1.71539128555513303061E-15, b25 = -1.16853328779934516808E-14, b24 = +7.67618549860493561688E-14;
      const double b23 = -4.85644678311192946090E-13, b22 = +2.95505266312963983461E-12, b21 = -1.72682629144155570723E-11;
      const double b20 = +9.67580903537323691224E-11, b19 = -5.18979560163526290666E-10, b18 = +2.65982372468238665035E-9;
      const double b17 = -1.30002500998624804212E-8, b16 = +6.04699502254191894932E-8, b15 = -2.67079385394061173391E-7;
      const double b14 = +1.11738753912010371815E-6, b13 = -4.41673835845875056359E-6, b12 = +1.64484480707288970893E-5;
      const double b11 = -5.75419501008210370398E-5, b10 = +1.88502885095841655729E-4, b09 = -5.76375574538582365885E-4;
      const double b08 = +1.63947561694133579842E-3, b07 = -4.32430999505057594430E-3, b06 = +1.05464603945949983183E-2;
      const double b05 = -2.37374148058994688156E-2, b04 = +4.93052842396707084878E-2, b03 = -9.49010970480476444210E-2;
      const double b02 = +1.71620901522208775349E-1, b01 = -3.04682672343198398683E-1, b00 = +6.76795274409476084995E-1;
      double b2 = b29, b1 = y * b2 + b28, b0 = y * b1 - b2 + b27;
      b2 = y * b0 - b1 + b26, b1 = y * b2 - b0 + b25, b0 = y * b1 - b2 + b24;
      b2 = y * b0 - b1 + b23, b1 = y * b2 - b0 + b22, b0 = y * b1 - b2 + b21;
      b2 = y * b0 - b1 + b20, b1 = y * b2 - b0 + b19, b0 = y * b1 - b2 + b18;
      b2 = y * b0 - b1 + b17, b1 = y * b2 - b0 + b16, b0 = y * b1 - b2 + b15;
      b2 = y * b0 - b1 + b14, b1 = y * b2 - b0 + b13, b0 = y * b1 - b2 + b12;
      b2 = y * b0 - b1 + b11, b1 = y * b2 - b0 + b10, b0 = y * b1 - b2 + b09;
      b2 = y * b0 - b1 + b08, b1 = y * b2 - b0 + b07, b0 = y * b1 - b2 + b06;
      b2 = y * b0 - b1 + b05, b1 = y * b2 - b0 + b04, b0 = y * b1 - b2 + b03;
      b2 = y * b0 - b1 + b02, b1 = y * b2 - b0 + b01, b0 = y * b1 - b2 + b00;
      return 0.5 * (b0 - b2) * exp(x);
   } else {
      double z = 32.0 / x - 2.0;
      const double b24 = -7.23318048787475395456E-18, b23 = -4.83050448594418207126E-18;
      const double b22 = +4.46562142029675999901E-17, b21 = +3.46122286769746109310E-17;
      const double b20 = -2.82762398051658348494E-16, b19 = -3.42548561967721913462E-16, b18 = +1.77256013305652638360E-15;
      const double b17 = +3.81168066935262242075E-15, b16 = -9.55484669882830764870E-15, b15 = -4.15056934728722208663E-14;
      const double b14 = +1.54008621752140982691E-14, b13 = +3.85277838274214270114E-13, b12 = +7.18012445138366623367E-13;
      const double b11 = -1.79417853150680611778E-12, b10 = -1.32158118404477131188E-11, b09 = -3.14991652796324136454E-11;
      const double b08 = +1.18891471078464383424E-11, b07 = +4.94060238822496958910E-10, b06 = +3.39623202570838634515E-9;
      const double b05 = +2.26666899049817806459E-8, b04 = +2.04891858946906374183E-7, b03 = +2.89137052083475648297E-6;
      const double b02 = +6.88975834691682398426E-5, b01 = +3.36911647825569408990E-3, b00 = +8.04490411014108831608E-1;
      double b0 = b24, b2 = z * b0 + b23, b1 = z * b2 - b0 + b22;
      b0 = z * b1 - b2 + b21;
      b2 = z * b0 - b1 + b20, b1 = z * b2 - b0 + b19, b0 = z * b1 - b2 + b18;
      b2 = z * b0 - b1 + b17, b1 = z * b2 - b0 + b16, b0 = z * b1 - b2 + b15;
      b2 = z * b0 - b1 + b14, b1 = z * b2 - b0 + b13, b0 = z * b1 - b2 + b12;
      b2 = z * b0 - b1 + b11, b1 = z * b2 - b0 + b10, b0 = z * b1 - b2 + b09;
      b2 = z * b0 - b1 + b08, b1 = z * b2 - b0 + b07, b0 = z * b1 - b2 + b06;
      b2 = z * b0 - b1 + b05, b1 = z * b2 - b0 + b04, b0 = z * b1 - b2 + b03;
      b2 = z * b0 - b1 + b02, b1 = z * b2 - b0 + b01, b0 = z * b1 - b2 + b00;
      return 0.5 * (b0 - b2) * exp(x) / sqrt(x);
   }
}

// Modified Bessel function of order one
//
// Returns modified Bessel function of order one of the
// argument.
//
// The function is defined as i1(x) = -i j1( ix ).
//
// The range is partitioned into the two intervals [0,8] and
// (8, infinity).  Chebyshev polynomial expansions are employed
// in each interval.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       30000       1.9e-15     2.1e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1985, 1987, 2000 by Stephen L. Moshier
// API: double besseli1(const double x);
double besseli1(double x) {
   double z = fabs(x);
   if (z <= 8.0) {
      double y = z / 2.0 - 2.0;
      const double b28 = +2.77791411276104639959E-18, b27 = -2.11142121435816608115E-17;
      const double b26 = +1.55363195773620046921E-16, b25 = -1.10559694773538630805E-15, b24 = +7.60068429473540693410E-15;
      const double b23 = -5.04218550472791168711E-14, b22 = +3.22379336594557470981E-13, b21 = -1.98397439776494371520E-12;
      const double b20 = +1.17361862988909016308E-11, b19 = -6.66348972350202774223E-11, b18 = +3.62559028155211703701E-10;
      const double b17 = -1.88724975172282928790E-9, b16 = +9.38153738649577178388E-9, b15 = -4.44505912879632808065E-8;
      const double b14 = +2.00329475355213526229E-7, b13 = -8.56872026469545474066E-7, b12 = +3.47025130813767847674E-6;
      const double b11 = -1.32731636560394358279E-5, b10 = +4.78156510755005422638E-5, b09 = -1.61760815825896745588E-4;
      const double b08 = +5.12285956168575772895E-4, b07 = -1.51357245063125314899E-3, b06 = +4.15642294431288815669E-3;
      const double b05 = -1.05640848946261981558E-2, b04 = +2.47264490306265168283E-2, b03 = -5.29459812080949914269E-2;
      const double b02 = +1.02643658689847095384E-1, b01 = -1.76416518357834055153E-1, b00 = +2.52587186443633654823E-1;
      double b1 = b28, b0 = y * b1 + b27, b2 = y * b0 - b1 + b26;
      b1 = y * b2 - b0 + b25, b0 = y * b1 - b2 + b24;
      b2 = y * b0 - b1 + b23, b1 = y * b2 - b0 + b22, b0 = y * b1 - b2 + b21;
      b2 = y * b0 - b1 + b20, b1 = y * b2 - b0 + b19, b0 = y * b1 - b2 + b18;
      b2 = y * b0 - b1 + b17, b1 = y * b2 - b0 + b16, b0 = y * b1 - b2 + b15;
      b2 = y * b0 - b1 + b14, b1 = y * b2 - b0 + b13, b0 = y * b1 - b2 + b12;
      b2 = y * b0 - b1 + b11, b1 = y * b2 - b0 + b10, b0 = y * b1 - b2 + b09;
      b2 = y * b0 - b1 + b08, b1 = y * b2 - b0 + b07, b0 = y * b1 - b2 + b06;
      b2 = y * b0 - b1 + b05, b1 = y * b2 - b0 + b04, b0 = y * b1 - b2 + b03;
      b2 = y * b0 - b1 + b02, b1 = y * b2 - b0 + b01, b0 = y * b1 - b2 + b00;
      z *= 0.5 * (b0 - b2) * exp(z);
   } else {
      double y = 32.0 / z - 2.0;
      const double b24 = +7.51729631084210481353E-18, b23 = +4.41434832307170791151E-18;
      const double b22 = -4.65030536848935832153E-17, b21 = -3.20952592199342395980E-17;
      const double b20 = +2.96262899764595013876E-16, b19 = +3.30820231092092828324E-16, b18 = -1.88035477551078244854E-15;
      const double b17 = -3.81440307243700780478E-15, b16 = +1.04202769841288027642E-14, b15 = +4.27244001671195135429E-14;
      const double b14 = -2.10154184277266431302E-14, b13 = -4.08355111109219731823E-13, b12 = -7.19855177624590851209E-13;
      const double b11 = +2.03562854414708950722E-12, b10 = +1.41258074366137813316E-11, b09 = +3.25260358301548823856E-11;
      const double b08 = -1.89749581235054123450E-11, b07 = -5.58974346219658380687E-10, b06 = -3.83538038596423702205E-9;
      const double b05 = -2.63146884688951950684E-8, b04 = -2.51223623787020892529E-7, b03 = -3.88256480887769039346E-6;
      const double b02 = -1.10588938762623716291E-4, b01 = -9.76109749136146840777E-3, b00 = +7.78576235018280120474E-1;
      double b0 = b24, b2 = y * b0 + b23, b1 = y * b2 - b0 + b22;
      b0 = y * b1 - b2 + b21;
      b2 = y * b0 - b1 + b20, b1 = y * b2 - b0 + b19, b0 = y * b1 - b2 + b18;
      b2 = y * b0 - b1 + b17, b1 = y * b2 - b0 + b16, b0 = y * b1 - b2 + b15;
      b2 = y * b0 - b1 + b14, b1 = y * b2 - b0 + b13, b0 = y * b1 - b2 + b12;
      b2 = y * b0 - b1 + b11, b1 = y * b2 - b0 + b10, b0 = y * b1 - b2 + b09;
      b2 = y * b0 - b1 + b08, b1 = y * b2 - b0 + b07, b0 = y * b1 - b2 + b06;
      b2 = y * b0 - b1 + b05, b1 = y * b2 - b0 + b04, b0 = y * b1 - b2 + b03;
      b2 = y * b0 - b1 + b02, b1 = y * b2 - b0 + b01, b0 = y * b1 - b2 + b00;
      z = 0.5 * (b0 - b2) * exp(z) / sqrt(z);
   }
   return x < 0.0 ? -z : +z;
}

// Modified Bessel function, second kind, order zero
//
// Returns modified Bessel function of the second kind
// of order zero of the argument.
//
// The range is partitioned into the two intervals [0,8] and
// (8, infinity).  Chebyshev polynomial expansions are employed
// in each interval.
//
// ACCURACY:
//
// Tested at 2000 random points between 0 and 8.  Peak absolute
// error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       30000       1.2e-15     1.6e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double besselk0(const double x);
double besselk0(double x) {
   double v;
   ae_assert(x > 0.0, "Domain error in BesselK0: x <= 0");
   if (x <= 2.0) {
      double y = x * x - 2.0;
      const double b09 = +1.37446543561352307156E-16, b08 = +4.25981614279661018399E-14;
      const double b07 = +1.03496952576338420167E-11, b06 = +1.90451637722020886025E-9;
      const double b05 = +2.53479107902614945675E-7, b04 = +2.28621210311945178607E-5, b03 = +1.26461541144692592338E-3;
      const double b02 = +3.59799365153615016266E-2, b01 = +3.44289899924628486886E-1, b00 = -5.35327393233902768720E-1;
      double b0 = b09, b2 = y * b0 + b08, b1 = y * b2 - b0 + b07;
      b0 = y * b1 - b2 + b06;
      b2 = y * b0 - b1 + b05, b1 = y * b2 - b0 + b04, b0 = y * b1 - b2 + b03;
      b2 = y * b0 - b1 + b02, b1 = y * b2 - b0 + b01, b0 = y * b1 - b2 + b00;
      v = 0.5 * (b0 - b2) - log(0.5 * x) * besseli0(x);
   } else {
      double z = 8.0 / x - 2.0;
      const double b24 = +5.30043377268626276149E-18, b23 = -1.64758043015242134646E-17;
      const double b22 = +5.21039150503902756861E-17, b21 = -1.67823109680541210385E-16;
      const double b20 = +5.51205597852431940784E-16, b19 = -1.84859337734377901440E-15, b18 = +6.34007647740507060557E-15;
      const double b17 = -2.22751332699166985548E-14, b16 = +8.03289077536357521100E-14, b15 = -2.98009692317273043925E-13;
      const double b14 = +1.14034058820847496303E-12, b13 = -4.51459788337394416547E-12, b12 = +1.85594911495471785253E-11;
      const double b11 = -7.95748924447710747776E-11, b10 = +3.57739728140030116597E-10, b09 = -1.69753450938905987466E-9;
      const double b08 = +8.57403401741422608519E-9, b07 = -4.66048989768794782956E-8, b06 = +2.76681363944501510342E-7;
      const double b05 = -1.83175552271911948767E-6, b04 = +1.39498137188764993662E-5, b03 = -1.28495495816278026384E-4;
      const double b02 = +1.56988388573005337491E-3, b01 = -3.14481013119645005427E-2, b00 = +2.44030308206595545468E0;
      double b0 = b24, b2 = z * b0 + b23, b1 = z * b2 - b0 + b22;
      b0 = z * b1 - b2 + b21;
      b2 = z * b0 - b1 + b20, b1 = z * b2 - b0 + b19, b0 = z * b1 - b2 + b18;
      b2 = z * b0 - b1 + b17, b1 = z * b2 - b0 + b16, b0 = z * b1 - b2 + b15;
      b2 = z * b0 - b1 + b14, b1 = z * b2 - b0 + b13, b0 = z * b1 - b2 + b12;
      b2 = z * b0 - b1 + b11, b1 = z * b2 - b0 + b10, b0 = z * b1 - b2 + b09;
      b2 = z * b0 - b1 + b08, b1 = z * b2 - b0 + b07, b0 = z * b1 - b2 + b06;
      b2 = z * b0 - b1 + b05, b1 = z * b2 - b0 + b04, b0 = z * b1 - b2 + b03;
      b2 = z * b0 - b1 + b02, b1 = z * b2 - b0 + b01, b0 = z * b1 - b2 + b00;
      v = 0.5 * (b0 - b2) * exp(-x) / sqrt(x);
   }
   return v;
}

// Modified Bessel function, second kind, order one
//
// Computes the modified Bessel function of the second kind
// of order one of the argument.
//
// The range is partitioned into the two intervals [0,2] and
// (2, infinity).  Chebyshev polynomial expansions are employed
// in each interval.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       30000       1.2e-15     1.6e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double besselk1(const double x);
double besselk1(double x) {
   double z = 0.5 * x;
   ae_assert(z > 0.0, "Domain error in K1");
   if (x <= 2.0) {
      double y = x * x - 2.0;
      const double b10 = -7.02386347938628759343E-18, b09 = -2.42744985051936593393E-15;
      const double b08 = -6.66690169419932900609E-13, b07 = -1.41148839263352776110E-10, b06 = -2.21338763073472585583E-8;
      const double b05 = -2.43340614156596823496E-6, b04 = -1.73028895751305206302E-4, b03 = -6.97572385963986435018E-3;
      const double b02 = -1.22611180822657148235E-1, b01 = -3.53155960776544875667E-1, b00 = +1.52530022733894777053E0;
      double b1 = b10, b0 = y * b1 + b09, b2 = y * b0 - b1 + b08;
      b1 = y * b2 - b0 + b07, b0 = y * b1 - b2 + b06;
      b2 = y * b0 - b1 + b05, b1 = y * b2 - b0 + b04, b0 = y * b1 - b2 + b03;
      b2 = y * b0 - b1 + b02, b1 = y * b2 - b0 + b01, b0 = y * b1 - b2 + b00;
      return log(z) * besseli1(x) + 0.5 * (b0 - b2) / x;
   } else {
      double y = 8.0 / x - 2.0;
      const double b24 = -5.75674448366501715755E-18, b23 = +1.79405087314755922667E-17;
      const double b22 = -5.68946255844285935196E-17, b21 = +1.83809354436663880070E-16;
      const double b20 = -6.05704724837331885336E-16, b19 = +2.03870316562433424052E-15, b18 = -7.01983709041831346144E-15;
      const double b17 = +2.47715442448130437068E-14, b16 = -8.97670518232499435011E-14, b15 = +3.34841966607842919884E-13;
      const double b14 = -1.28917396095102890680E-12, b13 = +5.13963967348173025100E-12, b12 = -2.12996783842756842877E-11;
      const double b11 = +9.21831518760500529508E-11, b10 = -4.19035475934189648750E-10, b09 = +2.01504975519703286596E-9;
      const double b08 = -1.03457624656780970260E-8, b07 = +5.74108412545004946722E-8, b06 = -3.50196060308781257119E-7;
      const double b05 = +2.40648494783721712015E-6, b04 = -1.93619797416608296024E-5, b03 = +1.95215518471351631108E-4;
      const double b02 = -2.85781685962277938680E-3, b01 = +1.03923736576817238437E-1, b00 = +2.72062619048444266945E0;
      double b0 = b24, b2 = y * b0 + b23, b1 = y * b2 - b0 + b22;
      b0 = y * b1 - b2 + b21;
      b2 = y * b0 - b1 + b20, b1 = y * b2 - b0 + b19, b0 = y * b1 - b2 + b18;
      b2 = y * b0 - b1 + b17, b1 = y * b2 - b0 + b16, b0 = y * b1 - b2 + b15;
      b2 = y * b0 - b1 + b14, b1 = y * b2 - b0 + b13, b0 = y * b1 - b2 + b12;
      b2 = y * b0 - b1 + b11, b1 = y * b2 - b0 + b10, b0 = y * b1 - b2 + b09;
      b2 = y * b0 - b1 + b08, b1 = y * b2 - b0 + b07, b0 = y * b1 - b2 + b06;
      b2 = y * b0 - b1 + b05, b1 = y * b2 - b0 + b04, b0 = y * b1 - b2 + b03;
      b2 = y * b0 - b1 + b02, b1 = y * b2 - b0 + b01, b0 = y * b1 - b2 + b00;
      return exp(-x) * 0.5 * (b0 - b2) / sqrt(x);
   }
}

// Modified Bessel function, second kind, integer order
//
// Returns modified Bessel function of the second kind
// of order n of the argument.
//
// The range is partitioned into the two intervals [0,9.55] and
// (9.55, infinity).  An ascending power series is used in the
// low range, and an asymptotic expansion in the high range.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,30        90000       1.8e-8      3.0e-10
//
// Error is high only near the crossover point x = 9.55
// between the two expansions used.
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
// API: double besselkn(const ae_int_t nn, const double x);
double besselkn(ae_int_t nn, double x) {
   ae_int_t n = nn < 0 ? -nn : +nn;
   ae_assert(n <= 31, "Overflow in BesselKN");
   ae_assert(x > 0.0, "Domain error in BesselKN");
   if (x <= 9.55) {
      double ans = 0.0, z0 = 0.25 * x * x, tox = 2.0 / x;
      double fn = 1.0, pn = 0.0, zmn = 1.0;
      if (n > 0) {
         pn = -Eul;
         for (ae_int_t i = 1; i < n; i++) {
            pn += 1.0 / i;
            fn *= i;
         }
         zmn = tox;
         if (n == 1) ans = 1.0 / x;
         else {
            double nk1f = fn / n, kf = 1.0, s = nk1f;
            double z = -z0, zn = 1.0;
            for (ae_int_t i = 1; i < n; i++) {
               nk1f /= n - i;
               kf *= i;
               zn *= z;
               double t = nk1f * zn / kf;
               s += t;
               ae_assert(maxrealnumber - fabs(t) > fabs(s), "Overflow in BesselKN");
               ae_assert(!(tox > 1.0 && maxrealnumber / tox < zmn), "Overflow in BesselKN");
               zmn *= tox;
            }
            s *= 0.5;
            double t = fabs(s);
            ae_assert(!(zmn > 1.0 && maxrealnumber / zmn < t), "Overflow in BesselKN");
            ae_assert(!(t > 1.0 && maxrealnumber / t < zmn), "Overflow in BesselKN");
            ans = s * zmn;
         }
      }
      double tlg = 2.0 * log(0.5 * x), pk = -Eul;
      double t;
      if (n == 0) {
         pn = pk;
         t = 1.0;
      } else {
         pn += 1.0 / n;
         t = 1.0 / fn;
      }
      double s = (pk + pn - tlg) * t;
      double k = 1.0;
      do {
         t *= z0 / (k * (k + n));
         pk += 1.0 / k;
         pn += 1.0 / (k + n);
         s += (pk + pn - tlg) * t;
         k++;
      } while (!SmallAtR(t / s, machineepsilon));
      s = 0.5 * s / zmn;
      if (n % 2 != 0) s = -s;
      return ans + s;
   }
   if (x > log(maxrealnumber)) return 0.0;
   double pn = 4.0 * (n * n), pk = 1.0, z0 = 8.0 * x, fn = 1.0, nkf = maxrealnumber;
   double t = 1.0, s = t;
   ae_int_t i = 0;
   do {
      double z = pn - pk * pk;
      t = t * z / (fn * z0);
      double nk1f = fabs(t);
      if (i >= n && nk1f > nkf) break;
      nkf = nk1f;
      s += t;
      fn++;
      pk += 2.0;
      i++;
   } while (!SmallAtR(t / s, machineepsilon));
   return exp(-x) * sqrt(pi / (2.0 * x)) * s;
}
} // end of namespace alglib_impl

namespace alglib {
double besselj0(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besselj0(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besselj1(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besselj1(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besseljn(const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besseljn(n, x);
   alglib_impl::ae_state_clear();
   return D;
}

double bessely0(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::bessely0(x);
   alglib_impl::ae_state_clear();
   return D;
}

double bessely1(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::bessely1(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besselyn(const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besselyn(n, x);
   alglib_impl::ae_state_clear();
   return D;
}

double besseli0(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besseli0(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besseli1(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besseli1(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besselk0(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besselk0(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besselk1(const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besselk1(x);
   alglib_impl::ae_state_clear();
   return D;
}

double besselkn(const ae_int_t nn, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::besselkn(nn, x);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib

// === LAGUERRE Package ===
namespace alglib_impl {
// Calculation of the value of the Laguerre polynomial.
//
// Parameters:
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Laguerre polynomial Ln at x
// API: double laguerrecalculate(const ae_int_t n, const double x);
double laguerrecalculate(ae_int_t n, double x) {
   double a;
   double b;
   double i;
   double result;
   result = 1.0;
   a = 1.0;
   b = 1.0 - x;
   if (n == 1) {
      result = b;
   }
   i = 2.0;
   while (i <= n) {
      result = ((2.0 * i - 1.0 - x) * b - (i - 1.0) * a) / i;
      a = b;
      b = result;
      i++;
   }
   return result;
}

// Summation of Laguerre polynomials using Clenshaw's recurrence formula.
//
// This routine calculates c[0]*L0(x) + c[1]*L1(x) + ... + c[N]*LN(x)
//
// Parameters:
//     n   -   degree, n >= 0
//     x   -   argument
//
// Result:
//     the value of the Laguerre polynomial at x
// API: double laguerresum(const real_1d_array &c, const ae_int_t n, const double x);
double laguerresum(RVector *c, ae_int_t n, double x) {
   double b1;
   double b2;
   ae_int_t i;
   double result;
   b1 = 0.0;
   b2 = 0.0;
   result = 0.0;
   for (i = n; i >= 0; i--) {
      result = (2.0 * i + 1.0 - x) * b1 / (i + 1) - (i + 1) * b2 / (i + 2) + c->xR[i];
      b2 = b1;
      b1 = result;
   }
   return result;
}

// Representation of Ln as C[0] + C[1]*X + ... + C[N]*X^N
//
// Inputs:
//     N   -   polynomial degree, n >= 0
//
// Outputs:
//     C   -   coefficients
// API: void laguerrecoefficients(const ae_int_t n, real_1d_array &c);
void laguerrecoefficients(ae_int_t n, RVector *c) {
   ae_int_t i;
   SetVector(c);
   ae_vector_set_length(c, n + 1);
   c->xR[0] = 1.0;
   for (i = 0; i < n; i++) {
      c->xR[i + 1] = -c->xR[i] * (n - i) / (i + 1) / (i + 1);
   }
}
} // end of namespace alglib_impl

namespace alglib {
double laguerrecalculate(const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::laguerrecalculate(n, x);
   alglib_impl::ae_state_clear();
   return D;
}

double laguerresum(const real_1d_array &c, const ae_int_t n, const double x) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::laguerresum(ConstT(ae_vector, c), n, x);
   alglib_impl::ae_state_clear();
   return D;
}

void laguerrecoefficients(const ae_int_t n, real_1d_array &c) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::laguerrecoefficients(n, ConstT(ae_vector, c));
   alglib_impl::ae_state_clear();
}
} // end of namespace alglib

// === ELLIPTIC Package ===
namespace alglib_impl {
// Complete elliptic integral of the first kind
//
// Approximates the integral
//
//            pi/2
//             -
//            | |
//            |           dt
// K(m)  =    |    ------------------
//            |                   2
//          | |    sqrt( 1 - m sin t )
//           -
//            0
//
// where m = 1 - m1, using the approximation
//
//     P(x)  -  log x Q(x).
//
// The argument m1 is used rather than m so that the logarithmic
// singularity at m = 1 will be shifted to the origin; this
// preserves maximum accuracy.
//
// K(0) = pi/2.
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE       0,1        30000       2.5e-16     6.8e-17
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double ellipticintegralkhighprecision(const double m1);
double ellipticintegralkhighprecision(double m1) {
   double p;
   double q;
   double result;
   if (m1 <= machineepsilon) {
      result = 1.3862943611198906188E0 - 0.5 * log(m1);
   } else {
      p = 1.37982864606273237150E-4;
      p = p * m1 + 2.28025724005875567385E-3;
      p = p * m1 + 7.97404013220415179367E-3;
      p = p * m1 + 9.85821379021226008714E-3;
      p = p * m1 + 6.87489687449949877925E-3;
      p = p * m1 + 6.18901033637687613229E-3;
      p = p * m1 + 8.79078273952743772254E-3;
      p = p * m1 + 1.49380448916805252718E-2;
      p = p * m1 + 3.08851465246711995998E-2;
      p = p * m1 + 9.65735902811690126535E-2;
      p = p * m1 + 1.38629436111989062502E0;
      q = 2.94078955048598507511E-5;
      q = q * m1 + 9.14184723865917226571E-4;
      q = q * m1 + 5.94058303753167793257E-3;
      q = q * m1 + 1.54850516649762399335E-2;
      q = q * m1 + 2.39089602715924892727E-2;
      q = q * m1 + 3.01204715227604046988E-2;
      q = q * m1 + 3.73774314173823228969E-2;
      q = q * m1 + 4.88280347570998239232E-2;
      q = q * m1 + 7.03124996963957469739E-2;
      q = q * m1 + 1.24999999999870820058E-1;
      q = q * m1 + 4.99999999999999999821E-1;
      result = p - q * log(m1);
   }
   return result;
}

// Complete elliptic integral of the first kind
//
// Approximates the integral
//
//            pi/2
//             -
//            | |
//            |           dt
// K(m)  =    |    ------------------
//            |                   2
//          | |    sqrt( 1 - m sin t )
//           -
//            0
//
// using the approximation
//
//     P(x)  -  log x Q(x).
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE       0,1        30000       2.5e-16     6.8e-17
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double ellipticintegralk(const double m);
double ellipticintegralk(double m) {
   double result;
   result = ellipticintegralkhighprecision(1.0 - m);
   return result;
}

// Incomplete elliptic integral of the first kind F(phi|m)
//
// Approximates the integral
//
//                phi
//                 -
//                | |
//                |           dt
// F(phi_\m)  =    |    ------------------
//                |                   2
//              | |    sqrt( 1 - m sin t )
//               -
//                0
//
// of amplitude phi and modulus m, using the arithmetic -
// geometric mean algorithm.
//
// ACCURACY:
//
// Tested at random points with m in [0, 1] and phi as indicated.
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE     -10,10       200000      7.4e-16     1.0e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// API: double incompleteellipticintegralk(const double phi, const double m);
double incompleteellipticintegralk(double phi, double m) {
   double a;
   double b;
   double c;
   double e;
   double temp;
   double t;
   double k;
   ae_int_t d;
   ae_int_t md;
   ae_int_t s;
   ae_int_t npio2;
   double result;
   if (m == 0.0) {
      result = phi;
      return result;
   }
   a = 1.0 - m;
   if (a == 0.0) {
      result = log(tan(0.5 * (HalfPi + phi)));
      return result;
   }
   npio2 = ifloor(phi / HalfPi);
   if (npio2 % 2 != 0) {
      npio2++;
   }
   if (npio2 != 0) {
      k = ellipticintegralk(1.0 - a);
      phi -= npio2 * HalfPi;
   } else {
      k = 0.0;
   }
   if (phi < 0.0) {
      phi = -phi;
      s = -1;
   } else {
      s = 0;
   }
   b = sqrt(a);
   t = tan(phi);
   if (!SmallAtR(t, 10.0)) {
      e = 1.0 / (b * t);
      if (SmallR(e, 10.0)) {
         e = atan(e);
         if (npio2 == 0) {
            k = ellipticintegralk(1.0 - a);
         }
         temp = k - incompleteellipticintegralk(e, m);
         if (s < 0) {
            temp = -temp;
         }
         result = temp + npio2 * k;
         return result;
      }
   }
   a = 1.0;
   c = sqrt(m);
   d = 1;
   md = 0;
   while (!SmallAtR(c / a, machineepsilon)) {
      temp = b / a;
      phi += atan(t * temp) + md * pi;
      md = itrunc((phi + HalfPi) / pi);
      t = t * (1.0 + temp) / (1.0 - temp * t * t);
      c = 0.5 * (a - b);
      temp = sqrt(a * b);
      a = 0.5 * (a + b);
      b = temp;
      d += d;
   }
   temp = (atan(t) + md * pi) / (d * a);
   if (s < 0) {
      temp = -temp;
   }
   result = temp + npio2 * k;
   return result;
}

// Complete elliptic integral of the second kind
//
// Approximates the integral
//
//            pi/2
//             -
//            | |                 2
// E(m)  =    |    sqrt( 1 - m sin t ) dt
//          | |
//           -
//            0
//
// using the approximation
//
//      P(x)  -  x log x Q(x).
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE       0, 1       10000       2.1e-16     7.3e-17
//
// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// API: double ellipticintegrale(const double m);
double ellipticintegrale(double m) {
   double p;
   double q;
   double result;
   ae_assert(m >= 0.0 && m <= 1.0, "Domain error in EllipticIntegralE: m < 0 or m > 1");
   m = 1.0 - m;
   if (m == 0.0) {
      result = 1.0;
      return result;
   }
   p = 1.53552577301013293365E-4;
   p = p * m + 2.50888492163602060990E-3;
   p = p * m + 8.68786816565889628429E-3;
   p = p * m + 1.07350949056076193403E-2;
   p = p * m + 7.77395492516787092951E-3;
   p = p * m + 7.58395289413514708519E-3;
   p = p * m + 1.15688436810574127319E-2;
   p = p * m + 2.18317996015557253103E-2;
   p = p * m + 5.68051945617860553470E-2;
   p = p * m + 4.43147180560990850618E-1;
   p = p * m + 1.00000000000000000299E0;
   q = 3.27954898576485872656E-5;
   q = q * m + 1.00962792679356715133E-3;
   q = q * m + 6.50609489976927491433E-3;
   q = q * m + 1.68862163993311317300E-2;
   q = q * m + 2.61769742454493659583E-2;
   q = q * m + 3.34833904888224918614E-2;
   q = q * m + 4.27180926518931511717E-2;
   q = q * m + 5.85936634471101055642E-2;
   q = q * m + 9.37499997197644278445E-2;
   q = q * m + 2.49999999999888314361E-1;
   result = p - q * m * log(m);
   return result;
}

// Incomplete elliptic integral of the second kind
//
// Approximates the integral
//
//                phi
//                 -
//                | |
//                |                   2
// E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
//                |
//              | |
//               -
//                0
//
// of amplitude phi and modulus m, using the arithmetic -
// geometric mean algorithm.
//
// ACCURACY:
//
// Tested at random arguments with phi in [-10, 10] and m in
// [0, 1].
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE     -10,10      150000       3.3e-15     1.4e-16
//
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1993, 2000 by Stephen L. Moshier
// API: double incompleteellipticintegrale(const double phi, const double m);
double incompleteellipticintegrale(double phi, double m) {
   double a;
   double b;
   double c;
   double e;
   double temp;
   double lphi;
   double t;
   double ebig;
   ae_int_t d;
   ae_int_t md;
   ae_int_t npio2;
   ae_int_t s;
   double result;
   if (m == 0.0) {
      result = phi;
      return result;
   }
   lphi = phi;
   npio2 = ifloor(lphi / HalfPi);
   if (npio2 % 2 != 0) {
      npio2++;
   }
   lphi -= npio2 * HalfPi;
   if (lphi < 0.0) {
      lphi = -lphi;
      s = -1;
   } else {
      s = 1;
   }
   a = 1.0 - m;
   ebig = ellipticintegrale(m);
   if (a == 0.0) {
      temp = sin(lphi);
      if (s < 0) {
         temp = -temp;
      }
      result = temp + npio2 * ebig;
      return result;
   }
   t = tan(lphi);
   b = sqrt(a);
// Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
// for pointing out an instability near odd multiples of pi/2
   if (!SmallAtR(t, 10.0)) {
   // Transform the amplitude
      e = 1.0 / (b * t);
   // ... but avoid multiple recursions.
      if (SmallR(e, 10.0)) {
         e = atan(e);
         temp = ebig + m * sin(lphi) * sin(e) - incompleteellipticintegrale(e, m);
         if (s < 0) {
            temp = -temp;
         }
         result = temp + npio2 * ebig;
         return result;
      }
   }
   c = sqrt(m);
   a = 1.0;
   d = 1;
   e = 0.0;
   md = 0;
   while (!SmallAtR(c / a, machineepsilon)) {
      temp = b / a;
      lphi += atan(t * temp) + md * pi;
      md = itrunc((lphi + HalfPi) / pi);
      t = t * (1.0 + temp) / (1.0 - temp * t * t);
      c = 0.5 * (a - b);
      temp = sqrt(a * b);
      a = 0.5 * (a + b);
      b = temp;
      d += d;
      e += c * sin(lphi);
   }
   temp = ebig / ellipticintegralk(m);
   temp *= (atan(t) + md * pi) / (d * a);
   temp += e;
   if (s < 0) {
      temp = -temp;
   }
   result = temp + npio2 * ebig;
   return result;
}
} // end of namespace alglib_impl

namespace alglib {
double ellipticintegralkhighprecision(const double m1) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::ellipticintegralkhighprecision(m1);
   alglib_impl::ae_state_clear();
   return D;
}

double ellipticintegralk(const double m) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::ellipticintegralk(m);
   alglib_impl::ae_state_clear();
   return D;
}

double incompleteellipticintegralk(const double phi, const double m) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::incompleteellipticintegralk(phi, m);
   alglib_impl::ae_state_clear();
   return D;
}

double ellipticintegrale(const double m) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::ellipticintegrale(m);
   alglib_impl::ae_state_clear();
   return D;
}

double incompleteellipticintegrale(const double phi, const double m) {
   alglib_impl::ae_state_init();
   TryCatch(0.0)
   double D = alglib_impl::incompleteellipticintegrale(phi, m);
   alglib_impl::ae_state_clear();
   return D;
}
} // end of namespace alglib
