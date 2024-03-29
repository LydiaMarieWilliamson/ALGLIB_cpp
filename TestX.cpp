#include <sstream> // For the stringstream type.
#include "DataAnalysis.h"
#include "Interpolation.h"

#if AE_OS == AE_POSIX
#   include <pthread.h>
typedef const pthread_attr_t ThAttr_t;
typedef void *ThArg_t;
typedef void *ThRet_t;
const ThRet_t ThNoRet = (ThRet_t)NULL;
typedef ThRet_t (*ThOp_t)(ThArg_t);
typedef pthread_t Thread_t;
inline int init_thread(Thread_t *Th, ThAttr_t *Attr, ThOp_t Op, ThArg_t Arg) { return pthread_create(Th, Attr, Op, Arg); }
#elif AE_OS == AE_WINDOWS
#   include <windows.h>
typedef LPSECURITY_ATTRIBUTES ThAttr_t;
typedef LPVOID ThArg_t;
typedef DWORD WINAPI *ThRet_t;
const ThRet_t ThNoRet = (ThRet_t)0;
typedef LPTHREAD_START_ROUTINE ThOp_t;
typedef HANDLE Thread_t;
inline int init_thread(Thread_t *Th, ThAttr_t *Attr, ThOp_t Op, ThArg_t Arg) { *Th = CreateThread(Attr, 0, Op, Arg, 0, NULL); return *Th != NULL; }
#else
// These are totally bogus stubs.
// You need to replace them with whatever specializations you need for your target configuration.
typedef void *ThAttr_t;
typedef void *ThArg_t;
typedef void ThRet_t;
const ThRet_t ThNoRet = (ThRet_t)NULL;
typedef ThRet_t (*ThOp_t)(ThArg_t);
struct Thread_t { ThAttr_t Attr; ThOp_t Op; ThArg_t Arg; };
inline int init_thread(Thread_t *Th, ThAttr_t *Attr, ThOp_t Op, ThArg_t Arg) { Th->Attr = Attr, Th->Op = Op, Th->Arg = Arg; return 0; }
#endif

using namespace alglib;

static inline bool NearR(double A, double B, double Tiny) { return fabs(A - B) <= Tiny; }
static inline bool SmallR(double A, double Tiny) { return fabs(A) <= Tiny; }

const char *fmt_str = "%-29s %s\n";
const char *fmt_speedup = "%-25s %5.1fx\n";

// Flag variables
bool issue505Ok = true;
bool issue478Ok = true;
bool issue528Ok = true;
bool issue591Ok = true;
bool issue594Ok = true;
bool issue764Ok = true;
bool issue813Ok = true;
bool issue824Ok = true;

// Service datatypes
struct innerrec {
   alglib_impl::complex cval;
   double rval;
   ae_int_t ival;
   bool bval;
   alglib_impl::ae_vector i1val;
};

static void innerrec_init(void *_p, bool make_automatic) {
   innerrec *p = (innerrec *)_p;
   alglib_impl::ae_vector_init(&p->i1val, 0, alglib_impl::DT_INT, make_automatic);
}

static void innerrec_copy(void *_dst, const void *_src, bool make_automatic) {
   innerrec *dst = (innerrec *)_dst;
   const innerrec *src = (const innerrec *)_src;
   dst->cval = src->cval;
   dst->rval = src->rval;
   dst->ival = src->ival;
   dst->bval = src->bval;
   alglib_impl::ae_vector_copy(&dst->i1val, &src->i1val, make_automatic);
}

static void innerrec_free(void *_p, bool make_automatic) {
   innerrec *p = (innerrec *)_p;
   alglib_impl::ae_vector_free(&p->i1val, make_automatic);
}

struct seedrec {
   bool bval;
   innerrec recval;
   alglib_impl::ae_shared_pool pool;
};

static void seedrec_init(void *_p, bool make_automatic) {
   seedrec *p = (seedrec *)_p;
   innerrec_init(&p->recval, make_automatic);
   alglib_impl::ae_shared_pool_init(&p->pool, make_automatic);
}

static void seedrec_copy(void *_dst, const void *_src, bool make_automatic) {
   seedrec *dst = (seedrec *)_dst;
   const seedrec *src = (const seedrec *)_src;
   dst->bval = src->bval;
   innerrec_copy(&dst->recval, &src->recval, make_automatic);
   alglib_impl::ae_shared_pool_copy(&dst->pool, &src->pool, make_automatic);
}

static void seedrec_free(void *_p, bool make_automatic) {
   seedrec *p = (seedrec *)_p;
   innerrec_free(&p->recval, make_automatic);
   alglib_impl::ae_shared_pool_free(&p->pool, make_automatic);
}

void func505_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) {
   double x0 = *((double *)ptr);
// This block assigns zero vector to gradient. Because gradient is a proxy vector
// (vector which uses another object as storage), sizes of gradient and vector being
// assigned must be equal. In this case data are copied in the memory linked with
// proxy.
//
// Early versions of ALGLIB failed to handle such assignment (it disrupted link
// between proxy vector and actual gradient stored in the internals of ALGLIB).
   real_1d_array z = "[0]";
   grad = "[0]";
   grad = z;
// This block tries to perform operations which are forbidden for proxy vector:
// * assign vector of non-matching size
// * change length of the vector
// Correct implementation must throw an exception without breaking a link between
// proxy object and actual vector.
   z = "[0,1]";
   try {
      grad = "[0,1]";
      issue505Ok = false;
   } catch(...) { }
   try {
      grad = z;
      issue505Ok = false;
   } catch(...) { }
   try {
      grad.setlength(10);
      issue505Ok = false;
   } catch(...) { }
   try {
      grad.setlength(1);
      issue505Ok = false;
   } catch(...) { }
//
// This block actually calculates function/gradient
//
   func = pow(x[0] - x0, 4);
   grad[0] = 4 * pow(x[0] - x0, 3);
}

void func505_vec(const real_1d_array &x, real_1d_array &fi, void *ptr) {
   double x0 = *((double *)ptr);
   fi[0] = x[0] - x0;
   fi[1] = pow(x[0] - x0, 2);
}

void func505_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr) {
   double x0 = *((double *)ptr);
//
// This block assigns zero matrix to Jacobian. Because Jacobian is a proxy matrix
// (matrix which uses another object as storage), sizes of Jacobian and matrix being
// assigned must be equal. In this case data are copied in the memory linked with
// proxy.
//
// Early versions of ALGLIB failed to handle such assignment (it disrupted link
// between proxy and actual matrix stored in the internals of ALGLIB).
//
   real_2d_array z = "[[0],[0]]";
   jac = "[[0],[0]]";
   jac = z;
//
// This block tries to perform operations which are forbidden for proxy vector:
// * assign vector of non-matching size
// * change length of the vector
// Correct implementation must throw an exception without breaking a link between
// proxy object and actual vector.
//
   try {
      jac = "[[0]]";
      issue505Ok = false;
   } catch(...) { }
   try {
      jac = "[[0,0],[1,1]]";
      issue505Ok = false;
   } catch(...) { }
   try {
      z = "[[0,1]]";
      jac = z;
      issue505Ok = false;
   } catch(...) { }
   try {
      jac.setlength(10, 6);
      issue505Ok = false;
   } catch(...) { }
   try {
      jac.setlength(2, 1);
      issue505Ok = false;
   } catch(...) { }
//
// This block actually calculates function/gradient
//
   fi[0] = x[0] - x0;
   fi[1] = pow(x[0] - x0, 2);
   jac[0][0] = 1.0;
   jac[1][0] = 2 * (x[0] - x0);
}

void issue813_callback(const real_1d_array &, real_1d_array &, void *) {
   throw 0;
}

void issue824_callback_i(const real_1d_array &, double &, void *) {
   throw (int *)NULL;
}

void issue824_callback_d(const real_1d_array &, double &, void *) {
   throw (double *)NULL;
}

void file_put_contents(const char *filename, const char *contents) {
   FILE *f = fopen(filename, "wb");
   if (f == NULL)
      ThrowError("file_put_contents: failed opening file");
   if (fwrite((void *)contents, 1, strlen(contents), f) != strlen(contents))
      ThrowError("file_put_contents: failed writing to file");
   fclose(f);
}

struct async_rbf_record {
   rbfmodel *p_model;
   rbfreport *p_report;
   bool thread_finished;
};

ThRet_t async_build_rbf_model(ThArg_t T) {
   async_rbf_record *p = (async_rbf_record *)T;
   rbfbuildmodel(*p->p_model, *p->p_report);
   p->thread_finished = true;
   return ThNoRet;
}

int main() {
// Report system properties
   printf("System:\n");
#ifdef AE_HPC
   printf("* cores count                %3ld\n", (long)alglib_impl::ae_cores_count());
#else
   printf("* cores count                %3ld\n", (long)1);
#endif
// Check status of allocation counter
#ifdef AE_USE_ALLOC_COUNTER
   printf("Allocation counter activated...\n");
   alglib_impl::_use_alloc_counter = true;
   if (alglib_impl::_alloc_counter != 0) {
      printf("Failed: _alloc_counter is non-zero on start!\n");
      return 1;
   }
   {
      real_1d_array x;
      x.setlength(1);
      if (alglib_impl::_alloc_counter == 0)
         printf("Warning: _alloc_counter is inactive!!!\n");
   }
   if (alglib_impl::_alloc_counter != 0) {
      printf("Failed: _alloc_counter does not decrease!\n");
      return 1;
   }
#else
   printf("No alloc counter.\nSome tests are skipped.\n");
#endif
// Testing basic functionality
   printf("Basic functions:\n");
   {
   // Testing 1D array functionality
      bool Ok = true;
      try {
      //
      // 1D boolean
      //
      // Default constructor, string constructor, copy constructor, assignment constructors:
      // * test that array sizes as reported by length match to what was specified
      // * test item-by-item access
      // * test to_string()
      // * test that modification of the copied array does not change original
      // * test that setlength() changes length
      // * test setcontent/getcontent
      // * test getcontent(), operator() and operator[] on constant arrays
      //   (in this case distinct implementation is used which must be tested separately)
      //
         boolean_1d_array arr_0, arr_1("[]"), arr_2("[true,false,true]"), arr_3(arr_2), arr_4, arr_5;
         arr_4 = arr_2;
         arr_5 = "[true,true,false]";
         Ok = Ok && arr_0.length() == 0;
         Ok = Ok && arr_1.length() == 0;
         Ok = Ok && arr_2.length() == 3;
         Ok = Ok && arr_3.length() == 3;
         Ok = Ok && arr_2[0] == arr_2(0) && arr_2[1] == arr_2(1) && arr_2[2] == arr_2(2);
         Ok = Ok && arr_2[0] && !arr_2[1] && arr_2[2];
         Ok = Ok && arr_3[0] && !arr_3[1] && arr_3[2];
         Ok = Ok && arr_4[0] && !arr_4[1] && arr_4[2];
         Ok = Ok && arr_5[0] && arr_5[1] && !arr_5[2];
         Ok = Ok && arr_2.tostring() == "[true,false,true]";
         Ok = Ok && arr_3.tostring() == "[true,false,true]";
         Ok = Ok && arr_4.tostring() == "[true,false,true]";
         Ok = Ok && arr_5.tostring() == "[true,true,false]";
         arr_2[0] = false;
         Ok = Ok && !arr_2[0] && arr_3[0] && arr_4[0];
         arr_5.setlength(99);
         Ok = Ok && arr_5.length() == 99;
      // setcontent/getcontent
         bool a0[] = { true, false, true, false, false };
         bool a0_mod = false;
         bool a0_orig = true;
         bool *p6;
         boolean_1d_array arr_6;
         arr_6.setcontent(5, a0);
         Ok = Ok && arr_6[0] == a0[0] && arr_6[1] == a0[1] && arr_6[2] == a0[2] && arr_6[3] == a0[3] && arr_6[4] == a0[4];
         p6 = arr_6.getcontent();
         Ok = Ok && p6 != a0;
         Ok = Ok && p6[0] == a0[0] && p6[1] == a0[1] && p6[2] == a0[2] && p6[3] == a0[3] && p6[4] == a0[4];
         a0[0] = a0_mod;
         Ok = Ok && arr_6[0] != a0[0];
         a0[0] = a0_orig;
      // operations on constant arrays
         {
            const boolean_1d_array &ac = arr_6;
            Ok = Ok && ac[0] == a0[0] && ac[1] == a0[1] && ac[2] == a0[2] && ac[3] == a0[3] && ac[4] == a0[4];
            Ok = Ok && ac(0) == a0[0] && ac(1) == a0[1] && ac(2) == a0[2] && ac(3) == a0[3] && ac(4) == a0[4];
            const bool *p = ac.getcontent();
            Ok = Ok && p[0] == a0[0] && p[1] == a0[1] && p[2] == a0[2] && p[3] == a0[3] && p[4] == a0[4];
         }
      //
      // Operations with proxy arrays:
      // * changes in target are propagated to proxy and vice versa
      // * assignments where proxy is source create new independent copy
      // * assignments to proxy are checked (their size must match to that of the target)
      // * incorrect assignments or attempts to change length must generate exception
      // * attempts to call setlength() must fail even when new size match original size
      //   of the array
      //
         boolean_1d_array targt, acopy;
         targt = "[true,false,false,true]";
         boolean_1d_array proxy(targt.c_ptr());
         acopy = proxy;
         Ok = Ok && targt[0] && !targt[1] && !targt[2] && targt[3];
         Ok = Ok && proxy[0] && !proxy[1] && !proxy[2] && proxy[3];
         Ok = Ok && acopy[0] && !acopy[1] && !acopy[2] && acopy[3];
         targt[0] = false;
         Ok = Ok && !targt[0] && !proxy[0] && acopy[0];
         proxy[0] = true;
         Ok = Ok && targt[0] && proxy[0] && acopy[0];
         acopy = "[false,true,true,true]";
         proxy = acopy;
         Ok = Ok && !targt[0] && targt[1] && targt[2] && targt[3];
         Ok = Ok && !proxy[0] && proxy[1] && proxy[2] && proxy[3];
         proxy = "[true,false,true,true]";
         Ok = Ok && targt[0] && !targt[1] && targt[2] && targt[3];
         Ok = Ok && proxy[0] && !proxy[1] && proxy[2] && proxy[3];
         try {
            acopy = "[false,true,true]";
            proxy = acopy;
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy = "[true,true,true]";
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy.setlength(100);
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy.setlength(proxy.length());
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
      } catch(...) {
         Ok = false;
      }
      try {
      //
      // 1D integer
      //
      // Default constructor, string constructor, copy constructor, assignment constructors:
      // * test that array sizes as reported by length match to what was specified
      // * test item-by-item access
      // * test to_string()
      // * test that modification of the copied array does not change original
      // * test that setlength() changes length
      //
         const char *s1 = "[2,3,-1]";
         const char *s2 = "[5,4,3]";
         const char *s3 = "[6,7,3,-4]";
         const char *s4 = "[9,5,-12,-0]";
         const char *s5 = "[1,7,2,1]";
         const char *s6 = "[7,7,7]";
         int v10 = 2, v11 = 3, v12 = -1, v10_mod = 9;
         int v20 = 5, v21 = 4, v22 = 3;
         int v30 = 6, v31 = 7, v32 = 3, v33 = -4, v30_mod = -6;
         int v40 = 9, v41 = 5, v42 = -12, v43 = 0;
         int v50 = 1, v51 = 7, v52 = 2, v53 = 1;
         integer_1d_array arr_0, arr_1("[]"), arr_2(s1), arr_3(arr_2), arr_4, arr_5;
         arr_4 = arr_2;
         arr_5 = s2;
         Ok = Ok && arr_0.length() == 0;
         Ok = Ok && arr_1.length() == 0;
         Ok = Ok && arr_2.length() == 3;
         Ok = Ok && arr_3.length() == 3;
         Ok = Ok && arr_2[0] == arr_2(0) && arr_2[1] == arr_2(1) && arr_2[2] == arr_2(2);
         Ok = Ok && arr_2[0] == v10 && arr_2[1] == v11 && arr_2[2] == v12;
         Ok = Ok && arr_3[0] == v10 && arr_3[1] == v11 && arr_3[2] == v12;
         Ok = Ok && arr_4[0] == v10 && arr_4[1] == v11 && arr_4[2] == v12;
         Ok = Ok && arr_5[0] == v20 && arr_5[1] == v21 && arr_5[2] == v22;
         Ok = Ok && arr_2.tostring() == s1;
         Ok = Ok && arr_3.tostring() == s1;
         Ok = Ok && arr_4.tostring() == s1;
         Ok = Ok && arr_5.tostring() == s2;
         arr_2[0] = v10_mod;
         Ok = Ok && arr_2[0] == v10_mod && arr_3[0] == v10 && arr_4[0] == v10;
         arr_5.setlength(99);
         Ok = Ok && arr_5.length() == 99;
      // setcontent/getcontent
         ae_int_t a0[] = { 2, 3, 1, 9, 2 };
         ae_int_t a0_mod = 7;
         ae_int_t a0_orig = 2;
         ae_int_t *p6;
         integer_1d_array arr_6;
         arr_6.setcontent(5, a0);
         Ok = Ok && arr_6[0] == a0[0] && arr_6[1] == a0[1] && arr_6[2] == a0[2] && arr_6[3] == a0[3] && arr_6[4] == a0[4];
         p6 = arr_6.getcontent();
         Ok = Ok && p6 != a0;
         Ok = Ok && p6[0] == a0[0] && p6[1] == a0[1] && p6[2] == a0[2] && p6[3] == a0[3] && p6[4] == a0[4];
         a0[0] = a0_mod;
         Ok = Ok && arr_6[0] != a0[0];
         a0[0] = a0_orig;
      // operations on constant arrays
         {
            const integer_1d_array &ac = arr_6;
            Ok = Ok && ac[0] == a0[0] && ac[1] == a0[1] && ac[2] == a0[2] && ac[3] == a0[3] && ac[4] == a0[4];
            Ok = Ok && ac(0) == a0[0] && ac(1) == a0[1] && ac(2) == a0[2] && ac(3) == a0[3] && ac(4) == a0[4];
            const ae_int_t *p = ac.getcontent();
            Ok = Ok && p[0] == a0[0] && p[1] == a0[1] && p[2] == a0[2] && p[3] == a0[3] && p[4] == a0[4];
         }
      //
      // Operations with proxy arrays:
      // * changes in target are propagated to proxy and vice versa
      // * assignments where proxy is source create new independent copy
      // * assignments to proxy are checked (their size must match to that of the target)
      // * incorrect assignments or attempts to change length must generate exception
      // * attempts to call setlength() must fail even when new size match original size
      //   of the array
      //
         integer_1d_array targt, acopy;
         targt = s3;
         integer_1d_array proxy(targt.c_ptr());
         acopy = proxy;
         Ok = Ok && targt[0] == v30 && targt[1] == v31 && targt[2] == v32 && targt[3] == v33;
         Ok = Ok && proxy[0] == v30 && proxy[1] == v31 && proxy[2] == v32 && proxy[3] == v33;
         Ok = Ok && acopy[0] == v30 && acopy[1] == v31 && acopy[2] == v32 && acopy[3] == v33;
         targt[0] = v30_mod;
         Ok = Ok && targt[0] == v30_mod && proxy[0] == v30_mod && acopy[0] == v30;
         proxy[0] = v30;
         Ok = Ok && targt[0] == v30 && proxy[0] == v30 && acopy[0] == v30;
         acopy = s4;
         proxy = acopy;
         Ok = Ok && targt[0] == v40 && targt[1] == v41 && targt[2] == v42 && targt[3] == v43;
         Ok = Ok && proxy[0] == v40 && proxy[1] == v41 && proxy[2] == v42 && proxy[3] == v43;
         proxy = s5;
         Ok = Ok && targt[0] == v50 && targt[1] == v51 && targt[2] == v52 && targt[3] == v53;
         Ok = Ok && proxy[0] == v50 && proxy[1] == v51 && proxy[2] == v52 && proxy[3] == v53;
         try {
            acopy = s6;
            proxy = acopy;
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy = s6;
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy.setlength(100);
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy.setlength(proxy.length());
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
      } catch(...) {
         Ok = false;
      }
      try {
      //
      // 1D real
      //
      // Default constructor, string constructor, copy constructor, assignment constructors:
      // * test that array sizes as reported by length match to what was specified
      // * test item-by-item access
      // * test to_string()
      // * test that modification of the copied array does not change original
      // * test that setlength() changes length
      //
         const char *s1 = "[2,3.5,-2.5E-1]";
         const char *s1_fmt = "[2.00,3.50,-0.25]";
         const char *s2 = "[5,4,3.126]";
         const char *s2_fmt = "[5.00,4.00,3.13]";
         const char *s3 = "[6,7,3,-4E2]";
         const char *s4 = "[9,5,-12,-0.01]";
         const char *s5 = "[1,7,2,1]";
         const char *s6 = "[7,7,7]";
         const int dps = 2;
         const double v10 = 2.0, v11 = 3.5, v12 = -0.25, v10_mod = 9.0;
         const double v20 = 5.0, v21 = 4.0, v22 = 3.126;
         const double v30 = 6.0, v31 = 7.0, v32 = 3.0, v33 = -400.0, v30_mod = -6.0;
         const double v40 = 9.0, v41 = 5.0, v42 = -12.0, v43 = -0.01;
         const double v50 = 1.0, v51 = 7.0, v52 = 2.0, v53 = 1.0;
         real_1d_array arr_0, arr_1("[]"), arr_2(s1), arr_3(arr_2), arr_4, arr_5;
         arr_4 = arr_2;
         arr_5 = s2;
         Ok = Ok && arr_0.length() == 0;
         Ok = Ok && arr_1.length() == 0;
         Ok = Ok && arr_2.length() == 3;
         Ok = Ok && arr_3.length() == 3;
         Ok = Ok && arr_2[0] == arr_2(0) && arr_2[1] == arr_2(1) && arr_2[2] == arr_2(2);
         Ok = Ok && arr_2[0] == v10 && arr_2[1] == v11 && arr_2[2] == v12;
         Ok = Ok && arr_3[0] == v10 && arr_3[1] == v11 && arr_3[2] == v12;
         Ok = Ok && arr_4[0] == v10 && arr_4[1] == v11 && arr_4[2] == v12;
         Ok = Ok && arr_5[0] == v20 && arr_5[1] == v21 && arr_5[2] == v22;
         Ok = Ok && arr_2.tostring(dps) == s1_fmt;
         Ok = Ok && arr_3.tostring(dps) == s1_fmt;
         Ok = Ok && arr_4.tostring(dps) == s1_fmt;
         Ok = Ok && arr_5.tostring(dps) == s2_fmt;
         arr_2[0] = v10_mod;
         Ok = Ok && arr_2[0] == v10_mod && arr_3[0] == v10 && arr_4[0] == v10;
         arr_5.setlength(99);
         Ok = Ok && arr_5.length() == 99;
      // setcontent/getcontent
         double a0[] = { 2, 3.5, 1, 9.125, 2 };
         double a0_mod = 7.0;
         double a0_orig = 2.0;
         double *p6;
         real_1d_array arr_6;
         arr_6.setcontent(5, a0);
         Ok = Ok && arr_6[0] == a0[0] && arr_6[1] == a0[1] && arr_6[2] == a0[2] && arr_6[3] == a0[3] && arr_6[4] == a0[4];
         p6 = arr_6.getcontent();
         Ok = Ok && p6 != a0;
         Ok = Ok && p6[0] == a0[0] && p6[1] == a0[1] && p6[2] == a0[2] && p6[3] == a0[3] && p6[4] == a0[4];
         a0[0] = a0_mod;
         Ok = Ok && arr_6[0] != a0[0];
         a0[0] = a0_orig;
      // operations on constant arrays
         {
            const real_1d_array &ac = arr_6;
            Ok = Ok && ac[0] == a0[0] && ac[1] == a0[1] && ac[2] == a0[2] && ac[3] == a0[3] && ac[4] == a0[4];
            Ok = Ok && ac(0) == a0[0] && ac(1) == a0[1] && ac(2) == a0[2] && ac(3) == a0[3] && ac(4) == a0[4];
            const double *p = ac.getcontent();
            Ok = Ok && p[0] == a0[0] && p[1] == a0[1] && p[2] == a0[2] && p[3] == a0[3] && p[4] == a0[4];
         }
      //
      // Operations with proxy arrays attached via attach_to(x_vector *):
      // * changes in target are propagated to proxy and vice versa
      // * assignments where proxy is source create new independent copy
      // * assignments to proxy are checked (their size must match to that of the target)
      // * incorrect assignments or attempts to change length must generate exception
      // * attempts to call setlength() must fail even when new size match original size
      //   of the array
      //
         {
            real_1d_array targt, acopy;
            targt = s3;
            real_1d_array proxy(targt.c_ptr());
            acopy = proxy;
            Ok = Ok && targt[0] == v30 && targt[1] == v31 && targt[2] == v32 && targt[3] == v33;
            Ok = Ok && proxy[0] == v30 && proxy[1] == v31 && proxy[2] == v32 && proxy[3] == v33;
            Ok = Ok && acopy[0] == v30 && acopy[1] == v31 && acopy[2] == v32 && acopy[3] == v33;
            targt[0] = v30_mod;
            Ok = Ok && targt[0] == v30_mod && proxy[0] == v30_mod && acopy[0] == v30;
            proxy[0] = v30;
            Ok = Ok && targt[0] == v30 && proxy[0] == v30 && acopy[0] == v30;
            acopy = s4;
            proxy = acopy;
            Ok = Ok && targt[0] == v40 && targt[1] == v41 && targt[2] == v42 && targt[3] == v43;
            Ok = Ok && proxy[0] == v40 && proxy[1] == v41 && proxy[2] == v42 && proxy[3] == v43;
            proxy = s5;
            Ok = Ok && targt[0] == v50 && targt[1] == v51 && targt[2] == v52 && targt[3] == v53;
            Ok = Ok && proxy[0] == v50 && proxy[1] == v51 && proxy[2] == v52 && proxy[3] == v53;
            try {
               acopy = s6;
               proxy = acopy;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
            try {
               proxy = s6;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
            try {
               proxy.setlength(100);
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
            try {
               proxy.setlength(proxy.length());
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         }
      //
      // >>> Unique for real_1d_array >>>
      //
      // Operations with proxy arrays attached via attach_to_ptr(double *):
      // * changes in target are propagated to proxy and vice versa
      // * assignments where proxy is source create new independent copy
      // * assignments to proxy are checked (their size must match to that of the target)
      // * incorrect assignments or attempts to change length must generate exception
      // * attempts to call setlength() must fail even when new size match original size
      //   of the array
      //
         {
            real_1d_array proxy, acopy;
            double targt[] = { v30, v31, v32, v33 };
            proxy.attach_to_ptr(4, targt);
            acopy = proxy;
            Ok = Ok && targt[0] == v30 && targt[1] == v31 && targt[2] == v32 && targt[3] == v33;
            Ok = Ok && proxy[0] == v30 && proxy[1] == v31 && proxy[2] == v32 && proxy[3] == v33;
            Ok = Ok && acopy[0] == v30 && acopy[1] == v31 && acopy[2] == v32 && acopy[3] == v33;
            targt[0] = v30_mod;
            Ok = Ok && targt[0] == v30_mod && proxy[0] == v30_mod && acopy[0] == v30;
            proxy[0] = v30;
            Ok = Ok && targt[0] == v30 && proxy[0] == v30 && acopy[0] == v30;
            acopy = s4;
            proxy = acopy;
            Ok = Ok && targt[0] == v40 && targt[1] == v41 && targt[2] == v42 && targt[3] == v43;
            Ok = Ok && proxy[0] == v40 && proxy[1] == v41 && proxy[2] == v42 && proxy[3] == v43;
            proxy = s5;
            Ok = Ok && targt[0] == v50 && targt[1] == v51 && targt[2] == v52 && targt[3] == v53;
            Ok = Ok && proxy[0] == v50 && proxy[1] == v51 && proxy[2] == v52 && proxy[3] == v53;
            try {
               acopy = s6;
               proxy = acopy;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
            try {
               proxy = s6;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
            try {
               proxy.setlength(100);
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
            try {
               proxy.setlength(proxy.length());
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         }
      } catch(...) {
         Ok = false;
      }
      try {
      //
      // 1D complex
      //
      // Default constructor, string constructor, copy constructor, assignment constructors:
      // * test that array sizes as reported by length match to what was specified
      // * test item-by-item access
      // * test to_string()
      // * test that modification of the copied array does not change original
      // * test that setlength() changes length
      //
         const char *s1 = "[2,3.5i,1-2.5E-1i]";
         const char *s1_fmt = "[2.00,3.50i,1.00-0.25i]";
         const char *s2 = "[5,-4+1i,3.126]";
         const char *s2_fmt = "[5.00,-4.00+1.00i,3.13]";
         const char *s3 = "[6,7,3,-4E2]";
         const char *s4 = "[9,5,-12,-0.01]";
         const char *s5 = "[1,7,2,1]";
         const char *s6 = "[7,7,7]";
         const int dps = 2;
         complex v10 = 2, v11 = complex(0, 3.5), v12 = complex(1, -0.25), v10_mod = 9;
         complex v20 = 5, v21 = complex(-4, 1), v22 = 3.126;
         complex v30 = 6, v31 = 7, v32 = 3, v33 = -400, v30_mod = -6;
         complex v40 = 9, v41 = 5, v42 = -12, v43 = -0.01;
         complex v50 = 1, v51 = 7, v52 = 2, v53 = 1;
         complex_1d_array arr_0, arr_1("[]"), arr_2(s1), arr_3(arr_2), arr_4, arr_5;
         arr_4 = arr_2;
         arr_5 = s2;
         Ok = Ok && arr_0.length() == 0;
         Ok = Ok && arr_1.length() == 0;
         Ok = Ok && arr_2.length() == 3;
         Ok = Ok && arr_3.length() == 3;
         Ok = Ok && arr_2[0] == arr_2(0) && arr_2[1] == arr_2(1) && arr_2[2] == arr_2(2);
         Ok = Ok && arr_2[0] == v10 && arr_2[1] == v11 && arr_2[2] == v12;
         Ok = Ok && arr_3[0] == v10 && arr_3[1] == v11 && arr_3[2] == v12;
         Ok = Ok && arr_4[0] == v10 && arr_4[1] == v11 && arr_4[2] == v12;
         Ok = Ok && arr_5[0] == v20 && arr_5[1] == v21 && arr_5[2] == v22;
         Ok = Ok && arr_2.tostring(dps) == s1_fmt;
         Ok = Ok && arr_3.tostring(dps) == s1_fmt;
         Ok = Ok && arr_4.tostring(dps) == s1_fmt;
         Ok = Ok && arr_5.tostring(dps) == s2_fmt;
         arr_2[0] = v10_mod;
         Ok = Ok && arr_2[0] == v10_mod && arr_3[0] == v10 && arr_4[0] == v10;
         arr_5.setlength(99);
         Ok = Ok && arr_5.length() == 99;
      // setcontent/getcontent
         complex a0[] = { 2, 3.5, 1, 9.125, 2 };
         complex a0_mod = 7;
         complex a0_orig = 2;
         complex *p6;
         complex_1d_array arr_6;
         arr_6.setcontent(5, a0);
         Ok = Ok && arr_6[0] == a0[0] && arr_6[1] == a0[1] && arr_6[2] == a0[2] && arr_6[3] == a0[3] && arr_6[4] == a0[4];
         p6 = arr_6.getcontent();
         Ok = Ok && p6 != a0;
         Ok = Ok && p6[0] == a0[0] && p6[1] == a0[1] && p6[2] == a0[2] && p6[3] == a0[3] && p6[4] == a0[4];
         a0[0] = a0_mod;
         Ok = Ok && arr_6[0] != a0[0];
         a0[0] = a0_orig;
      // operations on constant arrays
         {
            const complex_1d_array &ac = arr_6;
            Ok = Ok && ac[0] == a0[0] && ac[1] == a0[1] && ac[2] == a0[2] && ac[3] == a0[3] && ac[4] == a0[4];
            Ok = Ok && ac(0) == a0[0] && ac(1) == a0[1] && ac(2) == a0[2] && ac(3) == a0[3] && ac(4) == a0[4];
            const complex *p = ac.getcontent();
            Ok = Ok && p[0] == a0[0] && p[1] == a0[1] && p[2] == a0[2] && p[3] == a0[3] && p[4] == a0[4];
         }
      //
      // Operations with proxy arrays:
      // * changes in target are propagated to proxy and vice versa
      // * assignments where proxy is source create new independent copy
      // * assignments to proxy are checked (their size must match to that of the target)
      // * incorrect assignments or attempts to change length must generate exception
      // * attempts to call setlength() must fail even when new size match original size
      //   of the array
      //
         complex_1d_array targt, acopy;
         targt = s3;
         complex_1d_array proxy(targt.c_ptr());
         acopy = proxy;
         Ok = Ok && targt[0] == v30 && targt[1] == v31 && targt[2] == v32 && targt[3] == v33;
         Ok = Ok && proxy[0] == v30 && proxy[1] == v31 && proxy[2] == v32 && proxy[3] == v33;
         Ok = Ok && acopy[0] == v30 && acopy[1] == v31 && acopy[2] == v32 && acopy[3] == v33;
         targt[0] = v30_mod;
         Ok = Ok && targt[0] == v30_mod && proxy[0] == v30_mod && acopy[0] == v30;
         proxy[0] = v30;
         Ok = Ok && targt[0] == v30 && proxy[0] == v30 && acopy[0] == v30;
         acopy = s4;
         proxy = acopy;
         Ok = Ok && targt[0] == v40 && targt[1] == v41 && targt[2] == v42 && targt[3] == v43;
         Ok = Ok && proxy[0] == v40 && proxy[1] == v41 && proxy[2] == v42 && proxy[3] == v43;
         proxy = s5;
         Ok = Ok && targt[0] == v50 && targt[1] == v51 && targt[2] == v52 && targt[3] == v53;
         Ok = Ok && proxy[0] == v50 && proxy[1] == v51 && proxy[2] == v52 && proxy[3] == v53;
         try {
            acopy = s6;
            proxy = acopy;
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy = s6;
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy.setlength(100);
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
         try {
            proxy.setlength(proxy.length());
            Ok = false;
         } catch(ap_error e) {
         } catch(...) {
            Ok = false;
         }
      } catch(...) {
         Ok = false;
      }
   //
   // Report
   //
      printf(fmt_str, "* 1D arrays", Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!Ok)
         return 1;
   }
   {
   //
   // Testing 2D array functionality
   //
      bool Ok = true;
      try {
      //
      // 2D real
      //
      // Default constructor, string constructor, copy constructor, assignment constructors:
      // * test that array sizes as reported by length match to what was specified
      // * test item-by-item access
      // * test to_string()
      // * test that modification of the copied array does not change original
      // * test that setlength() changes length
      //
         const char *s1 = "[[2,3.5,-2.5E-1],[1,2,3]]";
         const char *s1_fmt = "[[2.00,3.50,-0.25],[1.00,2.00,3.00]]";
         const char *s2 = "[[5],[4],[3.126]]";
         const char *s2_fmt = "[[5.00],[4.00],[3.13]]";
         const char *s3 = "[[6,7],[3,-4E2],[-3,-1]]";
         const char *s4 = "[[9,5],[-12,-0.01],[-1,-2]]";
         const char *s5 = "[[1,7],[2,1],[0,4]]";
         const char *s60 = "[[7,7],[7,7]]";
         const char *s61 = "[[7],[7],[7]]";
         const int dps = 2;
         const double v10 = 2.0, v11 = 3.5, v12 = -0.25, v13 = 1.0, v14 = 2.0, v15 = 3.0, v10_mod = 9.0;
         const double v20 = 5.0, v21 = 4.0, v22 = 3.126;
         const double v30 = 6.0, v31 = 7.0, v32 = 3.0, v33 = -400.0, v34 = -3.0, v35 = -1.0, v30_mod = -6.0;
#if 0
         double v40 = 9.0, v41 = 5.0, v42 = -12.0, v43 = -0.01;
         double v50 = 1.0, v51 = 7.0, v52 = 2.0, v53 = 1.0;
#endif
         double r;
         real_2d_array arr_0, arr_1("[[]]"), arr_2(s1), arr_3(arr_2), arr_4, arr_5;
         arr_4 = arr_2;
         arr_5 = s2;
         Ok = Ok && arr_0.rows() == 0 && arr_0.cols() == 0 && arr_0.getstride() == 0;
         Ok = Ok && arr_1.rows() == 0 && arr_1.cols() == 0 && arr_1.getstride() == 0;
         Ok = Ok && arr_2.rows() == 2 && arr_2.cols() == 3 && arr_2.getstride() >= arr_2.cols();
         Ok = Ok && arr_3.rows() == 2 && arr_3.cols() == 3 && arr_3.getstride() >= arr_3.cols();
         Ok = Ok && arr_4.rows() == 2 && arr_4.cols() == 3 && arr_4.getstride() >= arr_4.cols();
         Ok = Ok && arr_5.rows() == 3 && arr_5.cols() == 1 && arr_5.getstride() >= arr_5.cols();
         Ok = Ok && arr_2[0][0] == arr_2(0, 0) && arr_2[0][1] == arr_2(0, 1) && arr_2[0][2] == arr_2(0, 2);
         Ok = Ok && arr_2[1][0] == arr_2(1, 0) && arr_2[1][1] == arr_2(1, 1) && arr_2[1][2] == arr_2(1, 2);
         Ok = Ok && arr_2[0][0] == v10 && arr_2[0][1] == v11 && arr_2[0][2] == v12;
         Ok = Ok && arr_2[1][0] == v13 && arr_2[1][1] == v14 && arr_2[1][2] == v15;
         Ok = Ok && arr_3[0][0] == v10 && arr_3[0][1] == v11 && arr_3[0][2] == v12;
         Ok = Ok && arr_3[1][0] == v13 && arr_3[1][1] == v14 && arr_3[1][2] == v15;
         Ok = Ok && arr_4[0][0] == v10 && arr_4[0][1] == v11 && arr_4[0][2] == v12;
         Ok = Ok && arr_4[1][0] == v13 && arr_4[1][1] == v14 && arr_4[1][2] == v15;
         Ok = Ok && arr_5[0][0] == v20 && arr_5[1][0] == v21 && arr_5[2][0] == v22;
         Ok = Ok && arr_2.tostring(dps) == s1_fmt;
         Ok = Ok && arr_3.tostring(dps) == s1_fmt;
         Ok = Ok && arr_4.tostring(dps) == s1_fmt;
         Ok = Ok && arr_5.tostring(dps) == s2_fmt;
         arr_2[0][0] = v10_mod;
         Ok = Ok && arr_2[0][0] == v10_mod && arr_3[0][0] == v10 && arr_4[0][0] == v10;
         arr_5.setlength(99, 97);
         Ok = Ok && arr_5.rows() == 99 && arr_5.cols() == 97;
      //
      // setcontent/elementwise access/constant arrays
      //
         ae_int_t n, m, i, j;
         for (n = 1; n <= 10; n++)
            for (m = 1; m <= 10; m++) {
               real_2d_array arr_6;
               double a0[100];
            // fill array by random values, test setcontent(0
               for (i = 0; i < m * n; i++)
                  a0[i] = randomreal();
               arr_6.setcontent(m, n, a0);
               for (i = 0; i < m; i++)
                  for (j = 0; j < n; j++) {
                     Ok = Ok && arr_6[i][j] == a0[i * n + j];
                     Ok = Ok && arr_6(i, j) == a0[i * n + j];
                  }
            // test that setcontent() actually copies data instead of creating just reference
               r = a0[0];
               a0[0]++;
               Ok = Ok && arr_6[0][0] != a0[0];
               a0[0] = r;
            // operations on constant arrays
               {
                  const real_2d_array &ac = arr_6;
                  for (i = 0; i < m; i++)
                     for (j = 0; j < n; j++) {
                        Ok = Ok && ac[i][j] == a0[i * n + j];
                        Ok = Ok && ac(i, j) == a0[i * n + j];
                     }
               }
            }
      //
      // Operations with proxy arrays:
      // * changes in target are propagated to proxy and vice versa
      // * assignments where proxy is source create new independent copy
      // * assignments to proxy are checked (their size must match to that of the target)
      // * incorrect assignments or attempts to change length must generate exception
      // * attempts to call setlength() must fail even when new size match original size
      //   of the array
      //
         { // test attach_to_ptr(x_matrix *)
         // subtest 0
            real_2d_array targt, acopy, acopy2;
            targt = s3;
            real_2d_array proxy(targt.c_ptr());
            acopy = proxy;
            for (i = 0; i < targt.rows(); i++)
               for (j = 0; j < targt.cols(); j++) {
                  Ok = Ok && proxy[i][j] == targt[i][j];
                  Ok = Ok && acopy[i][j] == targt[i][j];
               }
            r = targt[0][0];
            targt[0][0] = r + 1;
            Ok = Ok && targt[0][0] != r && proxy[0][0] != r && acopy[0][0] == r;
            proxy[0][0] = r;
            Ok = Ok && targt[0][0] == r && proxy[0][0] == r && acopy[0][0] == r;
         // subtest 1
            acopy = s4;
            proxy = acopy;
            for (i = 0; i < acopy.rows(); i++)
               for (j = 0; j < acopy.cols(); j++) {
                  Ok = Ok && proxy[i][j] == acopy[i][j];
                  Ok = Ok && targt[i][j] == acopy[i][j];
               }
            r = targt[0][0];
            targt[0][0] = r + 1;
            Ok = Ok && targt[0][0] != r && proxy[0][0] != r && acopy[0][0] == r;
            proxy[0][0] = r;
            Ok = Ok && targt[0][0] == r && proxy[0][0] == r && acopy[0][0] == r;
         // subtest 2
            acopy2 = s5;
            proxy = s5;
            for (i = 0; i < acopy.rows(); i++)
               for (j = 0; j < acopy.cols(); j++) {
                  Ok = Ok && proxy[i][j] == acopy2[i][j];
                  Ok = Ok && targt[i][j] == acopy2[i][j];
               }
         // error handling test 0
            try {
               acopy = s60;
               proxy = acopy;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 1
            try {
               acopy = s61;
               proxy = acopy;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 2
            try {
               proxy = s60;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 3
            try {
               proxy = s61;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 4
            try {
               proxy.setlength(100, 99);
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 5
            try {
               proxy.setlength(proxy.rows(), proxy.cols());
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         }
         { // test attach_to(double *)
         // subtest 0
            real_2d_array proxy, acopy, acopy2;
            double targt[] = { v30, v31, v32, v33, v34, v35 };
            const int NCOLS = 2;
            proxy.attach_to_ptr(3, 2, targt);
            acopy = proxy;
            for (i = 0; i < proxy.rows(); i++)
               for (j = 0; j < proxy.cols(); j++) {
                  Ok = Ok && proxy[i][j] == targt[i * NCOLS + j];
                  Ok = Ok && acopy[i][j] == targt[i * NCOLS + j];
               }
            r = targt[0 * NCOLS];
            targt[0 * NCOLS] = r + 1;
            Ok = Ok && targt[0 * NCOLS] != r && proxy[0][0] != r && acopy[0][0] == r;
            proxy[0][0] = r;
            Ok = Ok && targt[0 * NCOLS] == r && proxy[0][0] == r && acopy[0][0] == r;
         // subtest 1
            acopy = s4;
            proxy = acopy;
            for (i = 0; i < acopy.rows(); i++)
               for (j = 0; j < acopy.cols(); j++) {
                  Ok = Ok && proxy[i][j] == acopy[i][j];
                  Ok = Ok && targt[i * NCOLS + j] == acopy[i][j];
               }
            r = targt[0 * NCOLS];
            targt[0 * NCOLS] = r + 1;
            Ok = Ok && targt[0 * NCOLS] != r && proxy[0][0] != r && acopy[0][0] == r;
            proxy[0][0] = r;
            Ok = Ok && targt[0 * NCOLS] == r && proxy[0][0] == r && acopy[0][0] == r;
         // subtest 2
            acopy2 = s5;
            proxy = s5;
            for (i = 0; i < acopy.rows(); i++)
               for (j = 0; j < acopy.cols(); j++) {
                  Ok = Ok && proxy[i][j] == acopy2[i][j];
                  Ok = Ok && targt[i * NCOLS + j] == acopy2[i][j];
               }
         // error handling test 0
            try {
               acopy = s60;
               proxy = acopy;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 1
            try {
               acopy = s61;
               proxy = acopy;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 2
            try {
               proxy = s60;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 3
            try {
               proxy = s61;
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 4
            try {
               proxy.setlength(100, 99);
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         // error handling test 5
            try {
               proxy.setlength(proxy.rows(), proxy.cols());
               Ok = false;
            } catch(ap_error e) {
            } catch(...) {
               Ok = false;
            }
         }
      } catch(...) {
         Ok = false;
      }
   //
   // Report
   //
      printf(fmt_str, "* 2D arrays", Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!Ok)
         return 1;
   }
   {
   //
   // Testing CSV functionality
   //
      const char *csv_name = "alglib-tst-35252-ndg4sf.csv";
      bool Ok = true;
      try {
      // CSV_DEFAULT must be zero
         Ok = Ok && CSV_DEFAULT == 0;
      // absent file - must fail
         try {
            real_2d_array arr;
            read_csv("nonexistent123foralgtestinglib", '\t', CSV_DEFAULT, arr);
            Ok = false;
         } catch(ap_error) {
         } catch(...) {
            Ok = false;
         }
      // non-rectangular file - must fail
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "a,b,c\r\n1,2");
            read_csv(csv_name, ',', CSV_SKIP_HEADERS, arr);
            remove(csv_name);
            Ok = false;
         } catch(ap_error) {
         } catch(...) {
            Ok = false;
         }
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "a,b,c\r\n1,2,3,4");
            read_csv(csv_name, ',', CSV_SKIP_HEADERS, arr);
            remove(csv_name);
            Ok = false;
         } catch(ap_error) {
         } catch(...) {
            Ok = false;
         }
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "1,2,3,4\n1,2,3\n1,2,3");
            read_csv(csv_name, ',', CSV_DEFAULT, arr);
            remove(csv_name);
            Ok = false;
         } catch(ap_error) {
         } catch(...) {
            Ok = false;
         }
      // empty file
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "");
            read_csv(csv_name, '\t', CSV_DEFAULT, arr);
            remove(csv_name);
            Ok = Ok && arr.rows() == 0 && arr.cols() == 0;
         } catch(...) {
            Ok = false;
         }
      // one row with header, tab separator
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "a\tb\tc\n");
            read_csv(csv_name, '\t', CSV_SKIP_HEADERS, arr);
            remove(csv_name);
            Ok = Ok && arr.rows() == 0 && arr.cols() == 0;
         } catch(...) {
            Ok = false;
         }
      // no header, comma-separated, full stop as decimal point
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "1.5,2,3.25\n4,5,6");
            read_csv(csv_name, ',', CSV_DEFAULT, arr);
            remove(csv_name);
            Ok = Ok && arr.tostring(2) == "[[1.50,2.00,3.25],[4.00,5.00,6.00]]";
         } catch(...) {
            Ok = false;
         }
      // header, tab-separated, mixed use of comma and full stop as decimal points
         try {
            real_2d_array arr;
            file_put_contents(csv_name, "a\tb\tc\n1.5\t2\t3,25\n4\t5.25\t6,1\n");
            read_csv(csv_name, '\t', CSV_SKIP_HEADERS, arr);
            remove(csv_name);
            Ok = Ok && arr.tostring(2) == "[[1.50,2.00,3.25],[4.00,5.25,6.10]]";
         } catch(...) {
            Ok = false;
         }
      // header, tab-separated, fixed/exponential, spaces, mixed use of comma and full stop as decimal points
         try {
            real_2d_array arr;
            file_put_contents(csv_name, " a\t b \tc\n1,1\t 2.9\t -3.5  \n  1.1E1  \t 2.0E-1 \t-3E+1 \n+1  \t -2\t 3.    \n.1\t-.2\t+.3\n");
            read_csv(csv_name, '\t', CSV_SKIP_HEADERS, arr);
            remove(csv_name);
            Ok = Ok && arr.tostring(2) == "[[1.10,2.90,-3.50],[11.00,0.20,-30.00],[1.00,-2.00,3.00],[0.10,-0.20,0.30]]";
         } catch(...) {
            Ok = false;
         }
      } catch(...) {
         Ok = false;
      }
   //
   // Report
   //
      printf(fmt_str, "* CSV support", Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!Ok)
         return 1;
   }
//
// Serialization properties
//
   {
   //
   // Test kd-tree serialization
   //
      bool Ok = true;
      hqrndstate rs;
      kdtree tree0;
      real_2d_array xy, rxy0, rxy1;
      real_1d_array qx;
      const int npts = 50;
      const int nx = 2;
      const int ny = 1;
      int cnt0, cnt1;
      hqrndrandomize(rs);
      xy.setlength(npts, nx + ny);
      for (int i = 0; i < npts; i++)
         for (int j = 0; j < nx + ny; j++)
            xy[i][j] = hqrndnormal(rs);
      kdtreebuild(xy, npts, nx, ny, 2, tree0);
      qx.setlength(nx);
      try {
      // test string serialization/unserialization
         kdtree tree1;
         std::string s;
         kdtreeserialize(tree0, s);
         kdtreeunserialize(s, tree1);
         for (int i = 0; i < 100; i++) {
            for (int j = 0; j < nx; j++)
               qx[j] = hqrndnormal(rs);
            cnt0 = kdtreequeryknn(tree0, qx, 1, true);
            cnt1 = kdtreequeryknn(tree1, qx, 1, true);
            if (cnt0 != 1 || cnt1 != 1) {
               Ok = false;
               break;
            }
            kdtreequeryresultsxy(tree0, rxy0);
            kdtreequeryresultsxy(tree1, rxy1);
            for (int j = 0; j < nx + ny; j++)
               Ok = Ok && rxy0[0][j] == rxy1[0][j];
         }
      } catch(...) {
         Ok = false;
      }
      try {
      // test stream serialization/unserialization
      //
      // NOTE: we add a few symbols at the beginning and after the end of the data
      //       in order to test algorithm ability to work in the middle of the stream
         kdtree tree1;
         std::stringstream s;
         s.put('b');
         s.put('e');
         s.put('g');
         kdtreeserialize(tree0, s);
         s.put('e');
         s.put('n');
         s.put('d');
         s.seekg(0);
         Ok = Ok && s.get() == 'b';
         Ok = Ok && s.get() == 'e';
         Ok = Ok && s.get() == 'g';
         kdtreeunserialize(s, tree1);
         Ok = Ok && s.get() == 'e';
         Ok = Ok && s.get() == 'n';
         Ok = Ok && s.get() == 'd';
         for (int i = 0; i < 100; i++) {
            for (int j = 0; j < nx; j++)
               qx[j] = hqrndnormal(rs);
            cnt0 = kdtreequeryknn(tree0, qx, 1, true);
            cnt1 = kdtreequeryknn(tree1, qx, 1, true);
            if (cnt0 != 1 || cnt1 != 1) {
               Ok = false;
               break;
            }
            kdtreequeryresultsxy(tree0, rxy0);
            kdtreequeryresultsxy(tree1, rxy1);
            for (int j = 0; j < nx + ny; j++)
               Ok = Ok && rxy0[0][j] == rxy1[0][j];
         }
      } catch(...) {
         Ok = false;
      }
      try {
      // test string-to-stream serialization/unserialization
         kdtree tree1;
         std::string s0;
         kdtreeserialize(tree0, s0);
         std::stringstream s1(s0);
         kdtreeunserialize(s1, tree1);
         for (int i = 0; i < 100; i++) {
            for (int j = 0; j < nx; j++)
               qx[j] = hqrndnormal(rs);
            cnt0 = kdtreequeryknn(tree0, qx, 1, true);
            cnt1 = kdtreequeryknn(tree1, qx, 1, true);
            if (cnt0 != 1 || cnt1 != 1) {
               Ok = false;
               break;
            }
            kdtreequeryresultsxy(tree0, rxy0);
            kdtreequeryresultsxy(tree1, rxy1);
            for (int j = 0; j < nx + ny; j++)
               Ok = Ok && rxy0[0][j] == rxy1[0][j];
         }
      } catch(...) {
         Ok = false;
      }
      try {
      // test stream-to-string serialization/unserialization
         kdtree tree1;
         std::stringstream s0;
         kdtreeserialize(tree0, s0);
         std::string s1 = s0.str();
         kdtreeunserialize(s1, tree1);
         for (int i = 0; i < 100; i++) {
            for (int j = 0; j < nx; j++)
               qx[j] = hqrndnormal(rs);
            cnt0 = kdtreequeryknn(tree0, qx, 1, true);
            cnt1 = kdtreequeryknn(tree1, qx, 1, true);
            if (cnt0 != 1 || cnt1 != 1) {
               Ok = false;
               break;
            }
            kdtreequeryresultsxy(tree0, rxy0);
            kdtreequeryresultsxy(tree1, rxy1);
            for (int j = 0; j < nx + ny; j++)
               Ok = Ok && rxy0[0][j] == rxy1[0][j];
         }
      } catch(...) {
         Ok = false;
      }
   //
   // Report
   //
      printf(fmt_str, "* Serialization (kd-tree)", Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!Ok)
         return 1;
   }
   {
   //
   // Test legacy RBF interface
   //
      const char *pc_str =
         "50000000000 00000000000 20000000000 10000000000 A0000000000 30000000000 20000000000 00000000000 A0000000000 30000000000\r"
         "00000000000 20000000000 A0000000000 60000000000 00000000000\n"
         "00000000000 00000000000 00000000000 00000000000 00000000000\r\n"
         "00000000m_3 00000000000 00000000000 00000000m_3 00000000000\n\r"
         "00000000000 00000000004 00000000000 00000000000\t00000000004 00000000000 00000000000 00000000804 "
         "00000000000 00000000000 00000000804 00000000000 00000000000 00000000G04 00000000000 00000000000 "
         "00000000G04 00000000000 00000000000 00000000O04 00000000000 00000000000 00000000O04 00000000000 "
         "00000000000 00000000S04 00000000000 00000000000 00000000S04 00000000000 00000000000 00000000W04 "
         "00000000000 00000000000 00000000W04 00000000000 00000000000 00000000Y04 00000000000 00000000000 "
         "00000000Y04 00000000000 00000000000 00000000K04 00000000000 00000000000 00000000K04 00000000000 "
         "00000000000 A0000000000 00000000000 10000000000 20000000000 30000000000 40000000000 60000000000 "
         "70000000000 80000000000 90000000000 50000000000 30000000000 00000000000 00000000000 00000000000 "
         "30000000000 00000000Y04 00000000000 00000000000 u1000000000 00000000000 00000000000 00000000000 "
         "60000000000 80000000000 00000000000 50000000000 00000000000 50000000000 50000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 K0000000000 00000000I04 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 00000000000 "
         "00000000000 00000000000 A0000000000 30000000000 00000000000 00000000000 00000000000 00000000m_3 "
         "00000000000 00000000000 00000000004 00000000000 00000000000 00000000804 00000000000 00000000000 "
         "00000000G04 00000000000 00000000000 00000000K04 00000000000 00000000000 00000000O04 00000000000 "
         "00000000000 00000000S04 00000000000 00000000000 00000000W04 00000000000 00000000000 00000000Y04 "
         "00000000000 00000000000 A0000000000 40000000000 00000000q04 -pAGQnQBI14 UqUWierJ91C esm8ag6G61C "
         "00000000q04 4wcFMyCtu04 oPDvwHqst04 CExQXp8Ct04 00000000q04 litzPFhRb0C oKJvjcct314 5-fT-X8w614 "
         "00000000q04 3HSOsPVH11C vZWf4dgfv04 GbZg4MTJn04 00000000q04 iv7rMhuR71C hRtixp15r_3 EvCEDtLu-0C "
         "00000000q04 41CXzA_q71C umRYLK2yp0C 1zzY3Zqd91C 00000000q04 JvxJzDeI21C TVbyd7Ygz0C JLywRdR1n0C "
         "00000000q04 KmFarhc4g0C 1ehrn2tUt0C AECfwTIX814 00000000q04 Big__6hwt04 nSPzmAQrh_B 2H3o-KftH14 "
         "00000000q04 n1b9361vI14 mhJhviUE114 54a_qyBrH1C 00000000q04 10000000000 40000000000 StLCgor39-3 "
         "00000000000 00000000000 6qTG7Ae-1_3\n";
      real_1d_array ref_val("[-0.042560546916643, 0.942523544654062, 0.875197036560778, 0.0656948997826632, -0.743065973803404, -0.8903682039297, -0.26994815318748, 0.602248517290195, 0.980011992233124, 0.436594293214176]");
      bool Ok = true;
      try {
      // test unserialization from string without trailing end-of-stream symbol (dot)
      // this test is necessary for backward compatibility
         double eps = 0.0000000001;
         rbfmodel model;
         rbfunserialize(std::string(pc_str), model);
         for (int i = 0; i < ref_val.length(); i++)
            Ok = Ok && NearR(rbfcalc2(model, i, 0), ref_val[i], eps);
      } catch(...) {
         Ok = false;
      }
      try {
      // test unserialization from string with trailing end-of-stream symbol (dot)
      // this test is necessary for forward compatibility
         double eps = 0.0000000001;
         rbfmodel model;
         rbfunserialize(std::string(pc_str) + ".", model);
         for (int i = 0; i < ref_val.length(); i++)
            Ok = Ok && NearR(rbfcalc2(model, i, 0), ref_val[i], eps);
      } catch(...) {
         Ok = false;
      }
      try {
      // test unserialization from stream WITHOUT trailing end-of-stream symbol (dot)
      // this test MUST fail
         double eps = 0.0000000001;
         std::string _s(pc_str);
         std::istringstream stream(_s);
         rbfmodel model;
         rbfunserialize(stream, model);
         Ok = false;
      } catch(...) {
      } // Do nothing, it is expected to fail.
      try {
      // test unserialization from stream WITH trailing end-of-stream symbol (dot)
      // this test must succeed
         double eps = 0.0000000001;
         std::string _s = std::string(pc_str) + ".";
         std::istringstream stream(_s);
         rbfmodel model;
         rbfunserialize(stream, model);
         for (int i = 0; i < ref_val.length(); i++)
            Ok = Ok && NearR(rbfcalc2(model, i, 0), ref_val[i], eps);
      } catch(...) {
         Ok = false;
      }
      try {
      // test that we can read from the stream after unserialization
         double eps = 0.0000000001;
         std::string _s = std::string(pc_str) + ".<az>";
         std::istringstream stream(_s);
         rbfmodel model;
         rbfunserialize(stream, model);
         for (int i = 0; i < ref_val.length(); i++)
            Ok = Ok && NearR(rbfcalc2(model, i, 0), ref_val[i], eps);
         Ok = Ok && stream.get() == '<';
         Ok = Ok && stream.get() == 'a';
         Ok = Ok && stream.get() == 'z';
         Ok = Ok && stream.get() == '>';
      } catch(...) {
         Ok = false;
      }
   //
   // Report
   //
      printf(fmt_str, "* Serialization (RBF)", Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!Ok)
         return 1;
   }
   {
      bool Ok = true;
#if AE_OS == AE_POSIX || AE_OS == AE_WINDOWS
      hqrndstate rs;
      rbfmodel rbf;
      rbfreport rep;
      real_2d_array xy;
      int n = 1000, nx = 2, ny = 1;
      double rbase = 1.0;
      async_rbf_record async_rec;
      hqrndseed(464, 764, rs);
      xy.setlength(n, nx + ny);
      for (int i = 0; i < n; i++)
         for (int j = 0; j <= nx + ny; j++)
            xy[i][j] = hqrndnormal(rs);
      rbfcreate(nx, ny, rbf);
      rbfsetalgohierarchical(rbf, rbase, 1, 0.0);
      rbfsetpoints(rbf, xy);
      rbfsetv2its(rbf, 100000);
      Ok = Ok && rbfpeekprogress(rbf) == 0;
      async_rec.p_model = &rbf;
      async_rec.p_report = &rep;
      async_rec.thread_finished = false;
      Thread_t thread;
      if (init_thread(&thread, NULL, async_build_rbf_model, &async_rec) != 0) {
         printf(fmt_str, "* Progress/termination (RBF)", "Failed");
         printf(">>> unable to create background thread\n");
         fflush(stdout);
         return 1;
      }
//(@) Drop down to alglib_impl here, for now, because the state and frame accesses are no longer thread-safe.
#   define FullMultiThreaded 0
      for (double last_progress = 0.0; last_progress < 0.001; ) {
#   if FullMultiThreaded
         double new_progress = rbfpeekprogress(rbf);
#   else
         double new_progress = alglib_impl::rbfpeekprogress(rbf.c_ptr());
#   endif
         Ok = Ok && new_progress >= last_progress && new_progress <= 0.1; // We expect to terminate well before reaching 10%.
         last_progress = new_progress;
      }
#   if FullMultiThreaded
      rbfrequesttermination(rbf);
#   else
      alglib_impl::rbfrequesttermination(rbf.c_ptr());
#   endif
      while (!async_rec.thread_finished) {
#   if FullMultiThreaded
         double new_progress = rbfpeekprogress(rbf);
#   else
         double new_progress = alglib_impl::rbfpeekprogress(rbf.c_ptr());
#   endif
         Ok = Ok && (new_progress <= 0.1 || new_progress == 1.0); // We expect to terminate well before reaching 10%.
      }
#   undef FullMultiThreaded
      Ok = Ok && rbfpeekprogress(rbf) == 1;
      Ok = Ok && rep.terminationtype == 8;
      Ok = Ok && rbfcalc2(rbf, hqrndnormal(rs), hqrndnormal(rs)) == 0.0;
      printf(fmt_str, "* Progress/termination (RBF)", Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!Ok)
         return 1;
#else
      printf(fmt_str, "* Progress/termination (RBF)", "Multi-threading test disabled");
      fflush(stdout);
      if (!Ok)
         return 1;
#endif
   }
   {
   //
   // Test malloc() exceptions in constructors
   //
#ifdef AE_USE_ALLOC_COUNTER
      bool Ok = true;
      bool were_exceptions = false;
      for (int eidx = 0;; eidx++) // loop is terminated when we survive through all the tests
      {
      //
      // Select moment when we generate exception in the constructor
      //
         alglib_impl::_malloc_failure_after = alglib_impl::_alloc_counter_total + eidx;
      //
      // Perform many activities with ALGLIB, catch exceptions.
      // It is survival test, it checks that we survive exceptions.
      // If it fails, we are likely to end with abort().
      //
         try {
            {
               real_1d_array x0 = "[1,2,3]";
               real_1d_array x1(x0);
               real_1d_array x2;
               x2 = x1;
               real_1d_array *p = new real_1d_array("[1]");
               delete p;
            }
            {
               real_2d_array x0 = "[[1,2,3],[0,0,0]]";
               real_2d_array x1(x0);
               real_2d_array x2;
               x2 = x1;
               real_2d_array *p = new real_2d_array("[[1],[5]]");
               delete p;
            }
            {
               sparsematrix s;
               sparsecreate(2, 2, s);
               sparseset(s, 0, 0, 2.0);
               sparseset(s, 1, 1, 1.0);
               sparseset(s, 0, 1, 1.0);
               sparseadd(s, 1, 1, 4.0);
               sparsematrix s2(s), s3;
               s3 = s;
               double v;
               v = sparseget(s, 0, 0);
               v = sparseget(s, 0, 1);
               v = sparseget(s, 1, 0);
               v = sparseget(s, 1, 1);
               sparseconverttocrs(s);
               real_1d_array x = "[1,-1]";
               real_1d_array y = "[]";
               sparsemv(s, x, y);
            }
            {
               real_1d_array x = "[0,0]";
               double epsg = 0.0000000001;
               double epsf = 0.0;
               double epsx = 0.0;
               double stpmax = 0.1;
               ae_int_t maxits = 0;
               mincgstate state;
               mincgreport rep;
               mincgcreate(x, state);
               mincgsetcond(state, epsg, epsf, epsx, maxits);
               mincgsetstpmax(state, stpmax);
               mincgstate state2;
               state2 = state;
               mincgstate state3(state2);
               mincgrestartfrom(state, x);
            }
            {
               mlptrainer trn;
               multilayerperceptron network;
               mlpreport rep;
               real_2d_array xy = "[[1,1,1],[1,2,2],[2,1,2],[2,2,4]]";
               mlpcreatetrainer(2, 1, trn);
               mlpcreate1(2, 5, 1, network);
               mlpsetdataset(trn, xy, 4);
               mlptrainnetwork(trn, network, 5, rep);
               real_1d_array x = "[2,2]";
               real_1d_array y = "[0]";
               mlpprocess(network, x, y);
            }
            {
               std::string s;
               double v;
               rbfmodel model0;
               rbfmodel model1;
               real_2d_array xy = "[[-1,0,2],[+1,0,3]]";
               rbfreport rep;
               rbfcreate(2, 1, model0);
               rbfsetpoints(model0, xy);
               rbfsetalgohierarchical(model0, 1.0, 3, 0.0);
               rbfbuildmodel(model0, rep);
               rbfserialize(model0, s);
               rbfunserialize(s, model1);
               v = rbfcalc2(model0, 0.0, 0.0);
               v = rbfcalc2(model1, 0.0, 0.0);
               rbfbuildmodel(model0, rep);
               v = rbfcalc2(model0, 0.0, 0.0);
               rbfbuildmodel(model1, rep);
               v = rbfcalc2(model1, 0.0, 0.0);
            }
         //
         // We survived all tests, next iteration will bring no changed, terminate loop!
         //
            break;
         } catch(ap_error) {
            were_exceptions = true;
         }
      }
      alglib_impl::_malloc_failure_after = 0; // turn off artificial malloc failures
      printf(fmt_str, "* Exceptions in constructors", were_exceptions ? (Ok ? "Ok" : "Failed") : "..");
      fflush(stdout);
      if (!Ok)
         return 1;
#else
      printf(fmt_str, "* Exceptions in constructors", "??");
      fflush(stdout);
#endif
   }
   {
   //
   // Test multithreading-related settings
   //
   // For this test we measure performance of large NxNxN GEMMs
   // with different threading settings.
   //
      printf("SMP settings vs GEMM speedup:\n");
      if (ae_cores_count() > 1) {
         bool Ok = true;
         alglib_impl::ae_uint64_t default_global_threading = _ae_get_global_threading();
         ae_int_t default_nworkers = getnworkers();
         double time_default = 0.0, time_glob_ser = 0.0, time_glob_smp = 0.0, time_glob_ser_loc_ser = 0.0, time_glob_ser_loc_smp = 0.0, time_glob_smp_loc_ser = 0.0, time_glob_smp_loc_smp = 0.0, time_glob_smp_nw1 = 0.0;
         ae_int_t n = 800, mintime = 2000, cnt, t0;
         try {
         // allocate temporary matrices
            real_2d_array a, b, c;
            int i, j;
            a.setlength(n, n);
            b.setlength(n, n);
            c.setlength(n, n);
            for (i = 0; i < n; i++)
               for (j = 0; j < n; j++) {
                  a[i][j] = randomreal() - 0.5;
                  b[i][j] = randomreal() - 0.5;
                  c[i][j] = 0.0;
               }
         // measure time; interleave measurements with different settings in order to
         // reduce variance of results
            while (time_default < mintime) {
            // default threading
               t0 = alglib_impl::tickcount();
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_default += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global serial
               t0 = alglib_impl::tickcount();
               setglobalthreading(SerTH);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_ser += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global parallel
               t0 = alglib_impl::tickcount();
               setglobalthreading(ParTH);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_smp += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global serial, local serial
               t0 = alglib_impl::tickcount();
               setglobalthreading(SerTH);
               alglib_impl::ae_state_set_flags(SerTH);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_ser_loc_ser += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global serial, local parallel
               t0 = alglib_impl::tickcount();
               setglobalthreading(SerTH);
               alglib_impl::ae_state_set_flags(ParTH);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_ser_loc_smp += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global parallel, local serial
               t0 = alglib_impl::tickcount();
               setglobalthreading(ParTH);
               alglib_impl::ae_state_set_flags(SerTH);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_smp_loc_ser += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global parallel, local parallel
               t0 = alglib_impl::tickcount();
               setglobalthreading(ParTH);
               alglib_impl::ae_state_set_flags(ParTH);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_smp_loc_smp += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
            // global parallel, nworkers == 1
               t0 = alglib_impl::tickcount();
               setglobalthreading(ParTH);
               alglib_impl::ae_state_set_flags(NonTH);
               setnworkers(1);
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0);
               time_glob_smp_nw1 += alglib_impl::tickcount() - t0;
               _ae_set_global_threading(default_global_threading); // restore
               setnworkers(default_nworkers);
            }
         } catch(ap_error) {
            Ok = false;
         }
         printf(fmt_speedup, "* default speedup", time_glob_ser / time_glob_ser);
         printf(fmt_speedup, "* serial (global)", time_glob_ser / time_default);
         printf(fmt_speedup, "* serial (local)", time_glob_ser / time_glob_ser_loc_ser);
         printf(fmt_speedup, "* serial (nworkers == 1)", time_glob_ser / time_glob_smp_nw1);
         printf(fmt_speedup, "* parallel (global)", time_glob_ser / time_glob_smp);
         printf(fmt_speedup, "* parallel (local) v1", time_glob_ser / time_glob_ser_loc_smp);
         Ok = Ok && time_glob_ser / time_default > 0.85 && time_glob_ser / time_default < 1.15;
         Ok = Ok && time_glob_ser / time_glob_ser > 0.85 && time_glob_ser / time_glob_ser < 1.15;
         Ok = Ok && time_glob_ser / time_glob_ser_loc_ser > 0.85 && time_glob_ser / time_glob_ser_loc_ser < 1.15;
         Ok = Ok && time_glob_ser / time_glob_smp_loc_ser > 0.85 && time_glob_ser / time_glob_smp_loc_ser < 1.15;
         Ok = Ok && time_glob_ser / time_glob_smp_nw1 > 0.85 && time_glob_ser / time_glob_smp_nw1 < 1.15;
         Ok = Ok && time_glob_ser / time_glob_smp > 1.30;
         Ok = Ok && time_glob_ser / time_glob_ser_loc_smp > 1.30;
         Ok = Ok && time_glob_ser / time_glob_smp_loc_smp > 1.30;
         printf(fmt_str, "* test result", Ok ? "Ok" : "Failed (soft failure)");
         fflush(stdout);
      //
      // soft failure:
      // // if (!Ok)
      // //   return 1;
      //
      } else {
         printf(fmt_str, "* test skipped (no SMP)", "??");
         fflush(stdout);
      }
   }
//
// Testing issues which must be fixed
//
   printf("Issues:\n");
   {
   //
   // Testing issue #505 (http://bugs.alglib.net/view.php?id=505) in optimizers.
   // This issue was present in ALL optimizers, but we test it only on two: CG and LM.
   //
      try {
      //
      // Test CG
      // Stopping criteria - EpsX
      //
         mincgstate state;
         mincgreport rep;
         real_1d_array x = "[0.0]";
         double x0 = 10.0 * randommid();
         double epsx = 1.0E-9;
         mincgcreate(1, x, state);
         mincgsetcond(state, 0.0, 0.0, epsx, 0);
         mincgoptimize(state, func505_grad, NULL, &x0);
         mincgresults(state, x, rep);
         issue505Ok = issue505Ok && SmallR(4.0 * pow(x[0] - x0, 3), 0.001);
      } catch(...) {
         issue505Ok = false;
      }
      try {
      // Test LM
      // Stopping criteria - after |grad| < epsG
         minlmstate state;
         minlmreport rep;
         real_1d_array x = "[0.0]";
         double x0 = 10.0 * randommid();
         double epsx = 1.0E-9;
         minlmcreatevj(1, 2, x, state);
         minlmsetcond(state, epsx, 0);
         minlmoptimize(state, func505_vec, func505_jac, NULL, &x0);
         minlmresults(state, x, rep);
         issue505Ok = issue505Ok && NearR(x[0], x0, 0.001);
      } catch(...) {
         issue505Ok = false;
      }
      printf(fmt_str, "* Issue 505", issue505Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!issue505Ok)
         return 1;
   // Testing issue #478 (http://bugs.alglib.net/view.php?id=478)
   // in high-quality RNG. It has to correctly handle random numbers
   // larger than 2^31.
   //
   // This test is performed only in 64-bit mode.
      if (sizeof(ae_int_t) > 4) {
      //
      // 64-bit mode, perform test:
      // * use large NMax > 2^31
      // * generate 1.000.000 random numbers
      // * use two bins - one for numbers less then NMax/2,
      //   another one for the rest of them
      // * bin sizes are equal to n0, n1
      // * both bins should be approximately equal, we use
      //   ad hoc threshold 0.45 < n0,n1 < 0.55.
      //
         try {
            hqrndstate rs;
            ae_int_t nmax[3];
            ae_int_t ncnt = 3, nidx;
            double n0, n1;
            hqrndrandomize(rs);
         //
         // nmax:
         // * first nmax is just large value to test basic uniformity of generator
         //
            nmax[0] = 1000000;
            nmax[0] *= nmax[0];
            nmax[1] = 2147483562;
            nmax[1] *= 1.5;
            nmax[2] = 2147483562;
            nmax[2] *= 3;
            for (nidx = 0; nidx < ncnt; nidx++) {
               n0 = 0;
               n1 = 0;
               for (int i = 0; i < 1000000; i++) {
                  ae_int_t v = hqrnduniformi(rs, nmax[nidx]);
                  if (v < nmax[nidx] / 2)
                     n0++;
                  else
                     n1++;
                  issue478Ok = issue478Ok && v >= 0 && v < nmax[nidx];
               }
               issue478Ok = issue478Ok && n0 / (n0 + n1) > 0.45;
               issue478Ok = issue478Ok && n0 / (n0 + n1) < 0.55;
               issue478Ok = issue478Ok && n1 / (n0 + n1) > 0.45;
               issue478Ok = issue478Ok && n1 / (n0 + n1) < 0.55;
            }
         } catch(...) {
            issue478Ok = false;
         }
         printf(fmt_str, "* Issue 478", issue478Ok ? "Ok" : "Failed");
         fflush(stdout);
         if (!issue478Ok)
            return 1;
      } else {
      //
      // 32-bit mode, skip test
      //
         printf(fmt_str, "* Issue 478", "Ok (skipped in 32-bit mode)");
         fflush(stdout);
      }
   // Testing issue #528 (http://bugs.alglib.net/view.php?id=528)
   // in shared pool and smart pointer which leak memory.
   //
   // In order to test it we create pool, seed it with specially
   // created structure, perform several operations, then clear it.
   // We test allocation counter before and after this operation.
   //
#ifdef AE_USE_ALLOC_COUNTER
      try {
         int alloc_cnt;
         alglib_impl::ae_frame _frame_block;
      // case #0: just seeding the pool
         alloc_cnt = alglib_impl::_alloc_counter;
         alglib_impl::ae_state_init();
         alglib_impl::ae_frame_make(&_frame_block);
         NewObj(alglib_impl::ae_shared_pool, pool);
         NewObj(seedrec, seed);
         alglib_impl::ae_shared_pool_set_seed(&pool, &seed, sizeof seed, seedrec_init, seedrec_copy, seedrec_free);
         alglib_impl::ae_state_clear();
         issue528Ok = issue528Ok && alloc_cnt == alglib_impl::_alloc_counter;
      // case #1: seeding and retrieving, not recycling
         alloc_cnt = alglib_impl::_alloc_counter;
         alglib_impl::ae_state_init();
         alglib_impl::ae_frame_make(&_frame_block);
         RefObj(void, p0);
      // NewObj(alglib_impl::ae_shared_pool, pool); NewObj(seedrec, seed); // ... but without the declarations.
         memset(&pool, 0, sizeof pool), alglib_impl::ae_shared_pool_init(&pool, true);
         memset(&seed, 0, sizeof seed), seedrec_init(&seed, true);
         alglib_impl::ae_shared_pool_set_seed(&pool, &seed, sizeof seed, seedrec_init, seedrec_copy, seedrec_free);
         alglib_impl::ae_shared_pool_retrieve(&pool, &_p0);
         alglib_impl::ae_state_clear();
         issue528Ok = issue528Ok && alloc_cnt == alglib_impl::_alloc_counter;
      // case #2: seeding and retrieving twice to different pointers, recycling both
         alloc_cnt = alglib_impl::_alloc_counter;
         alglib_impl::ae_state_init();
         alglib_impl::ae_frame_make(&_frame_block);
      // RefObj(void, p0); // ... but without the declaration.
         memset(&_p0, 0, sizeof _p0), alglib_impl::ae_smart_ptr_init(&_p0, (void **)&p0, true);
         RefObj(void, p1);
      // NewObj(alglib_impl::ae_shared_pool, pool); NewObj(seedrec, seed); // ... but without the declarations.
         memset(&pool, 0, sizeof pool), alglib_impl::ae_shared_pool_init(&pool, true);
         memset(&seed, 0, sizeof seed), seedrec_init(&seed, true);
         alglib_impl::ae_shared_pool_set_seed(&pool, &seed, sizeof seed, seedrec_init, seedrec_copy, seedrec_free);
         alglib_impl::ae_shared_pool_retrieve(&pool, &_p0);
         alglib_impl::ae_shared_pool_retrieve(&pool, &_p1);
         alglib_impl::ae_shared_pool_recycle(&pool, &_p0);
         alglib_impl::ae_shared_pool_recycle(&pool, &_p1);
         alglib_impl::ae_state_clear();
         issue528Ok = issue528Ok && alloc_cnt == alglib_impl::_alloc_counter;
      } catch(...) {
         issue528Ok = false;
      }
      printf(fmt_str, "* Issue 528", issue528Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!issue528Ok)
         return 1;
#else
      printf(fmt_str, "* Issue 528", "??");
      fflush(stdout);
#endif
   // Testing issue #591 (http://bugs.alglib.net/view.php?id=591)
   // in copying of object containing shared pool as one of its
   // fields.
   //
   // Unfixed ALGLIB crashes because of unneeded assertion in the
   // ae_shared_pool_copy() function.
   //
      try {
         multilayerperceptron net0, net1;
         real_1d_array x("[1,2]"), y0("[0,0]"), y1("[0,0]"), y2("[0,0]");
         mlpcreate0(2, 2, net0);
         mlpprocess(net0, x, y0);
      //
      // Test assignment constructor
      //
         net1 = net0;
         mlpprocess(net1, x, y1);
         issue591Ok = issue591Ok && NearR(y0[0], y1[0], 1.0E-9) && NearR(y0[1], y1[1], 1.0E-9);
      //
      // Test copy constructor
      //
         multilayerperceptron net2(net0);
         mlpprocess(net2, x, y2);
         issue591Ok = issue591Ok && NearR(y0[0], y2[0], 1.0E-9) && NearR(y0[1], y2[1], 1.0E-9);
      } catch(...) {
         issue591Ok = false;
      }
      printf(fmt_str, "* Issue 591", issue591Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!issue591Ok)
         return 1;
   //
   // Task #594 (http://bugs.alglib.net/view.php?id=594) - additional
   // test for correctness of copying of objects. When we copy ALGLIB
   // object, indenendent new copy is created.
   //
   // This test checks both copying with copy constructor and assignment
   // constructor
   //
      try {
         multilayerperceptron net0, net1;
         real_1d_array x("[1,2]"), y0("[0,0]"), y1("[0,0]"), y2("[0,0]");
         mlpcreate0(2, 2, net0);
         mlpprocess(net0, x, y0);
      //
      // Test assignment and copy constructors:
      // * copy object with one of the constructors
      // * process vector with original network
      // * randomize original network
      // * process vector with copied networks and compare
      //
         net1 = net0;
         multilayerperceptron net2(net0);
         mlprandomize(net0);
         mlpprocess(net1, x, y1);
         mlpprocess(net2, x, y2);
         issue594Ok = issue594Ok && NearR(y0[0], y1[0], 1.0E-9) && NearR(y0[1], y1[1], 1.0E-9);
         issue594Ok = issue594Ok && NearR(y0[0], y2[0], 1.0E-9) && NearR(y0[1], y2[1], 1.0E-9);
      } catch(...) {
         issue594Ok = false;
      }
      printf(fmt_str, "* Issue 594", issue594Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!issue594Ok)
         return 1;
   //
   // Issue #764, potential memory leak in the smart pointer
   //
#ifdef AE_USE_ALLOC_COUNTER
      try {
         int alloc_cnt;
         alglib_impl::ae_frame _frame_block;
      // seeding shared pool and retrieving twice to same pointer, no recycling
         alloc_cnt = alglib_impl::_alloc_counter;
         alglib_impl::ae_state_init();
         alglib_impl::ae_frame_make(&_frame_block);
         RefObj(void, p0);
         NewObj(alglib_impl::ae_shared_pool, pool);
         NewObj(seedrec, seed);
         alglib_impl::ae_shared_pool_set_seed(&pool, &seed, sizeof seed, seedrec_init, seedrec_copy, seedrec_free);
         alglib_impl::ae_shared_pool_retrieve(&pool, &_p0);
         alglib_impl::ae_shared_pool_retrieve(&pool, &_p0);
         alglib_impl::ae_state_clear();
         issue764Ok = issue764Ok && alloc_cnt == alglib_impl::_alloc_counter;
      } catch(...) {
         issue764Ok = false;
      }
      printf(fmt_str, "* Issue 764", issue764Ok ? "Ok" : "Failed");
      fflush(stdout);
      if (!issue764Ok)
         return 1;
#else
      printf(fmt_str, "* Issue 764", "??");
      fflush(stdout);
#endif
   //
   // Issue 813: MSVC is unable to handle longjmp() from catch() block; the code
   //            cycles forever.
   //
      {
         minlmstate state;
         real_1d_array x;
         x.setlength(1);
         x[0] = 0;
         minlmcreatev(1, x, 1e-5, state);
         issue813Ok = false;
         try {
            minlmoptimize(state, &issue813_callback);
         } catch(...) {
            issue813Ok = true;
         }
         printf(fmt_str, "* Issue 813", issue813Ok ? "Ok" : "Failed");
         fflush(stdout);
         if (!issue813Ok)
            return 1;
      }
   //
   // Issue 824: pre-3.16 versions of ALGLIB hide exceptions generated in user callbacks
   //
      {
         mincgstate state;
         real_1d_array x;
         x.setlength(1);
         x[0] = 0;
         mincgcreatef(1, x, 1e-5, state);
         issue824Ok = true;
      // throw int*
         try {
            mincgoptimize(state, &issue824_callback_i);
         } catch(int *) {
         } catch(double *) {
            issue824Ok = false;
         } catch(...) {
            issue824Ok = false;
         }
      // throw double*
         try {
            mincgoptimize(state, &issue824_callback_d);
         } catch(int *) {
            issue824Ok = false;
         } catch(double *) {
         } catch(...) {
            issue824Ok = false;
         }
      // done
         printf(fmt_str, "* Issue 824", issue824Ok ? "Ok" : "Failed");
         fflush(stdout);
         if (!issue824Ok)
            return 1;
      }
   }
//
// Performance testing
//
   printf("Performance:\n");
   {
      {
         int _n[] = { 16, 32, 64, 1024, 0 };
         int i, j, k, t, nidx;
         for (nidx = 0; _n[nidx] != 0; nidx++) {
         //
         // Settings:
         // * n - matrix size
         // * nrepeat - number of repeated multiplications, always divisible by 4
         //
            int n = _n[nidx];
            double desiredflops = n > 64 ? 1.0E10 : 1.0E9;
            int nrepeat = (int)(desiredflops / (2 * pow(n, 3.0)));
            nrepeat = 4 * (nrepeat / 4 + 1);
         //
         // Actual processing
         //
            real_2d_array a, b, c;
            double perf0, perf1, perf2;
            a.setlength(n, n);
            b.setlength(n, n);
            c.setlength(n, n);
            for (i = 0; i < n; i++)
               for (j = 0; j < n; j++) {
                  a[i][j] = randomreal() - 0.5;
                  b[i][j] = randomreal() - 0.5;
                  c[i][j] = 0.0;
               }
            t = alglib_impl::tickcount();
            for (k = 0; k < nrepeat; k++)
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, k % 2, b, 0, 0, (k / 2) % 2, 0.0, c, 0, 0);
            t = alglib_impl::tickcount() - t;
            perf0 = 0.000001 * pow(n, 3.0) * 2.0 * nrepeat / (0.001 * t);
            printf("* RGEMM-SEQ-%-4ld (MFLOPS)  %5.0lf\n", (long)n, perf0);
            setnworkers(0);
            t = alglib_impl::tickcount();
            for (k = 0; k < nrepeat; k++)
               alglib_impl::ae_state_set_flags(ParTH),
               rmatrixgemm(n, n, n, 1.0, a, 0, 0, k % 2, b, 0, 0, (k / 2) % 2, 0.0, c, 0, 0),
               alglib_impl::ae_state_set_flags(NonTH);
            t = alglib_impl::tickcount() - t;
            perf2 = 0.000001 * pow(n, 3.0) * 2.0 * nrepeat / (0.001 * t);
            printf("* RGEMM-MTN-%-4ld           %4.1lfx\n", (long)n, perf2 / perf0);
            setnworkers(1);
         }
      }
   }
// Check allocation counter on exit
#ifdef AE_USE_ALLOC_COUNTER
   printf("Allocation counter checked... ");
#   ifdef _ALGLIB_HAS_WORKSTEALING
   alglib_impl::ae_free_disposed_items();
   alglib_impl::ae_complete_finalization_before_exit();
#   endif
   if (alglib_impl::_alloc_counter != 0) {
      printf("Failed: _alloc_counter is non-zero on end!\n");
      return 1;
   } else printf("Ok\n");
#endif
// Return
   return 0;
}
