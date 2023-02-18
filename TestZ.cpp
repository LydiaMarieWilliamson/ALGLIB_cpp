#include <sstream>
#include "Interpolation.h" // Only SPLINR1D is tested here.

#if !defined AE_NO_EXCEPTIONS
#   error "This test should be compiled with AE_NO_EXCEPTIONS defined"
#endif

using namespace alglib;

static inline bool NearR(double A, double B, double Tiny) { return fabs(A - B) <= Tiny; }

const char *fmt_str = "%-29s %s\n";

int main() {
   printf("Test exception-free error handling:\n");
#ifdef AE_USE_ALLOC_COUNTER
   printf("Allocation counter activated...\n");
   alglib_impl::_use_alloc_counter = true;
   if (alglib_impl::_alloc_counter != 0) {
      printf("Failure: _alloc_counter is non-zero on start!\n");
      return 1;
   }
{{
   real_1d_array x;
   x.setlength(1);
   if (alglib_impl::_alloc_counter == 0) printf(":::: Warning: _alloc_counter is inactive!!! :::::\\n");
}
   if (alglib_impl::_alloc_counter != 0) {
      printf("Failure: _alloc_counter does not decrease!\n");
      return 1;
   }
}
#endif
{ // Test the default state of the error flag.
   bool Ok = true;
   Ok = Ok && !get_error_flag();
   printf(fmt_str, "* default flag value", Ok ? "Ok" : "Failed");
   fflush(stdout);
   if (!Ok) return 1;
} { // Test errors during array creation.
   bool Ok = true;
   clear_error_flag();
// Constructors succeed with working malloc().
   clear_error_flag();
   real_1d_array r1;
   Ok = Ok && !get_error_flag();
// Even constructors with broken malloc() succeed.
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   real_1d_array r2;
   Ok = Ok && !get_error_flag();
// but setlength() fails without malloc().
   clear_error_flag();
   r2.setlength(5);
   Ok = Ok && get_error_flag();
// clear_error_flag() clears the error flag.
   clear_error_flag();
   Ok = Ok && !get_error_flag();
// Without clear_error_flag(), the error flag is not reset by successful calls.
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   r2.setlength(5);
   Ok = Ok && get_error_flag() && r2.length() == 0;
   alglib_impl::_force_malloc_failure = false;
   r2.setlength(6);
   Ok = Ok && get_error_flag() && r2.length() == 6;
   clear_error_flag();
   r2.setlength(7);
   Ok = Ok && !get_error_flag() && r2.length() == 7;
// Assignments to empty arrays require malloc().
   clear_error_flag();
   alglib_impl::_force_malloc_failure = false;
   real_1d_array r3;
   r2.setlength(1);
   r2[0] = 123.25;
   alglib_impl::_force_malloc_failure = true;
   r3 = r2;
   Ok = Ok && get_error_flag() && r3.length() == 0;
   alglib_impl::_force_malloc_failure = false;
   clear_error_flag();
   r3 = r2;
   Ok = Ok && !get_error_flag() && r3.length() == 1 && r3[0] == 123.25;
// Assignments to non-empty arrays do NOT require malloc().
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   r2[0] = 345;
   r3 = r2;
   Ok = Ok && !get_error_flag() && r3.length() == 1 && r3[0] == 345;
   alglib_impl::_force_malloc_failure = false;
   printf(fmt_str, "* 1D arrays", Ok ? "Ok" : "Failed");
   fflush(stdout);
   if (!Ok) return 1;
} {
   bool Ok = true;
   clear_error_flag();
// Constructors succeed with working malloc().
   clear_error_flag();
   real_2d_array r1;
   Ok = Ok && !get_error_flag();
// Constructors even with broken malloc() succeed,
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   real_2d_array r2;
   Ok = Ok && !get_error_flag();
// but setlength() fails without malloc().
   clear_error_flag();
   r2.setlength(5, 6);
   Ok = Ok && get_error_flag();
// clear_error_flag() clears the error flag.
   clear_error_flag();
   Ok = Ok && !get_error_flag();
// Without clear_error_flag(), the error flag is not reset by successful calls.
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   r2.setlength(5, 6);
   Ok = Ok && get_error_flag() && r2.rows() == 0 && r2.cols() == 0;
   alglib_impl::_force_malloc_failure = false;
   r2.setlength(6, 7);
   Ok = Ok && get_error_flag() && r2.rows() == 6 && r2.cols() == 7;
   clear_error_flag();
   r2.setlength(7, 8);
   Ok = Ok && !get_error_flag() && r2.rows() == 7 && r2.cols() == 8;
// Assignments to empty arrays require malloc().
   clear_error_flag();
   alglib_impl::_force_malloc_failure = false;
   real_2d_array r3;
   r2.setlength(1, 1);
   r2[0][0] = 123.25;
   alglib_impl::_force_malloc_failure = true;
   r3 = r2;
   Ok = Ok && get_error_flag() && r3.rows() == 0 && r3.cols() == 0;
   alglib_impl::_force_malloc_failure = false;
   clear_error_flag();
   r3 = r2;
   Ok = Ok && !get_error_flag() && r3.rows() == 1 && r3.cols() == 1 && r3[0][0] == 123.25;
// Assignments to non-empty arrays do NOT require malloc().
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   r2[0][0] = 345;
   r3 = r2;
   Ok = Ok && !get_error_flag() && r3.rows() == 1 && r3.cols() == 1 && r3[0][0] == 345;
   alglib_impl::_force_malloc_failure = false;
   printf(fmt_str, "* 2D arrays", Ok ? "Ok" : "Failed");
   fflush(stdout);
   if (!Ok) return 1;
} { // Test ALGLIB objects.
   bool Ok = true;
   clear_error_flag();
// Prepare the data for the tests.
   real_1d_array x, y;
   x.setlength(2), x[0] = 0, x[1] = 1;
   y.setlength(2), y[0] = 2, y[1] = 3;
//(@) The tests marked (@) were changed, because the constructors and destructors no longer use malloc() and free().
// Constructors with working malloc() succeed.
   clear_error_flag();
   spline1dinterpolant s1;
   Ok = Ok && !get_error_flag();
// Constructors with broken malloc() fail. //(@) Not anymore.
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   spline1dinterpolant s2;
#if 0
   Ok = Ok && get_error_flag(); //(@) Was originally this.
#else
   Ok = Ok && !get_error_flag(); //(@) No more memory management is used.
#endif
   alglib_impl::_force_malloc_failure = false;
// Constructors with correct malloc() succeed.
   clear_error_flag();
   spline1dbuildlinear(x, y, 2, s1);
   Ok = Ok && !get_error_flag() && NearR(spline1dcalc(s1, 0.5), 2.5, 1.0E-12);
// Assignments with broken malloc() fail.
   clear_error_flag();
   spline1dinterpolant s3;
   alglib_impl::_force_malloc_failure = true;
   s3 = s1;
   alglib_impl::_force_malloc_failure = false;
   Ok = Ok && get_error_flag();
// Assignments with broken objects fail, but do not crash.
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   spline1dinterpolant s3b;
#if 0
   Ok = Ok && s3b.c_ptr() == NULL; //(@) Was originally this.
#else
   Ok = Ok && s3b.c_ptr() != NULL; //(@) s3b.c_ptr() != NULL for object types, now.
#endif
   s3b = s1;
   alglib_impl::_force_malloc_failure = false;
   Ok = Ok && get_error_flag();
// Assignments with working malloc() succeed.
   clear_error_flag();
   s3 = s1;
   Ok = Ok && !get_error_flag() && NearR(spline1dcalc(s3, 0.5), 2.5, 1.0E-12);
// Copy constructors with broken malloc() fail.
   clear_error_flag();
   alglib_impl::_force_malloc_failure = true;
   spline1dinterpolant s4(s1);
   alglib_impl::_force_malloc_failure = false;
   Ok = Ok && get_error_flag();
// Copy constructors with working malloc() succeed.
   clear_error_flag();
   spline1dinterpolant s5(s1);
   Ok = Ok && !get_error_flag() && NearR(spline1dcalc(s5, 0.5), 2.5, 1.0E-12);
   printf(fmt_str, "* ALGLIB objects", Ok ? "Ok" : "Failed");
   fflush(stdout);
   if (!Ok) return 1;
} { // Test ALGLIB functions.
   bool Ok = true;
   clear_error_flag();
   spline1dinterpolant s1;
   real_1d_array x, y;
   x.setlength(2), x[0] = 0, x[1] = 1;
   y.setlength(2), y[0] = 2, y[1] = NAN;
   Ok = Ok && !get_error_flag();
   spline1dbuildlinear(x, y, 2, s1);
   Ok = Ok && get_error_flag();
   printf(fmt_str, "* ALGLIB functions", Ok ? "Ok" : "Failed");
   fflush(stdout);
   if (!Ok) return 1;
}
// Allocation counter.
#ifdef AE_USE_ALLOC_COUNTER
   printf("Allocation counter checked... ");
   if (alglib_impl::_alloc_counter != 0) {
      printf("Failure: _alloc_counter is non-zero on end!\n");
      return 1;
   } else printf("Ok\n");
#endif
// Return.
   return 0;
}
