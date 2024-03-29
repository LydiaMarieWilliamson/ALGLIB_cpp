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

// if AE_OS == AE_LINUX (will be redefined to AE_POSIX in Ap.h),
// set _GNU_SOURCE flag BEFORE any #includes to get affinity management functions.
#if AE_OS == AE_LINUX && !defined _GNU_SOURCE
#   define _GNU_SOURCE
#endif

// Must be defined before we include Ap.h.
#define _ALGLIB_IMPL_DEFINES
#define _ALGLIB_INTEGRITY_CHECKS_ONCE
#include "Ap.h"
#include <limits>
#include <locale.h>
#include <ctype.h>
#if defined AE_CPU && AE_CPU == AE_INTEL
#   if AE_COMPILER == AE_GNUC && 0
#      include <fpu_control.h> //(@) For _FPU_GETCW and _FPU_SETCW, for a fix that's apparently no longer needed.
#   elif AE_COMPILER == AE_MSVC
#      include <intrin.h>
#   endif
#endif

// OS-specific includes.
#if AE_OS == AE_POSIX || defined AE_DEBUG4POSIX
#   include <time.h>
#   include <unistd.h>
#   include <pthread.h>
#   include <sched.h>
#   include <sys/time.h> // For tickcount().
#elif AE_OS == AE_WINDOWS || defined AE_DEBUG4WINDOWS
#   ifndef _WIN32_WINNT
#      define _WIN32_WINNT 0x0501
#   endif
#   include <windows.h>
#   include <process.h>
#endif

// Entropy source
#if ALGLIB_ENTROPY_SRC == ALGLIB_ENTROPY_SRC_OPENSSL
#   include <openssl/rand.h>
#   define ALGLIB_OPENSSL_RAND_MAX             0x7fffffff
#endif

// Debugging helpers for Windows.
#ifdef AE_DEBUG4WINDOWS
#   include <windows.h>
#endif

namespace alglib_impl {
// Core Code (Vectors, Matrices, Memory Management, etc.)
// Local definitions.
#define AE_DATA_ALIGN	0x40
#define AE_PTR_ALIGN	sizeof(void *)

#define AE_LITTLE_ENDIAN	1
#define AE_BIG_ENDIAN		2
#define AE_MIXED_ENDIAN		3

#define AE_SER_ENTRY_LENGTH	11
#define AE_SER_ENTRIES_PER_ROW	5

#if defined ALGLIB_REDZONE
#   define _ALGLIB_REDZONE_VAL	0x3c
#endif

// These declarations are used to ensure, at compile-time, that
//	sizeof(bool) == 1, sizeof(ae_int32_t) == 4, sizeof(ae_int64_t) == 8, sizeof(ae_uint64_t) == 8, sizeof(ae_int_t) == sizeof(void *).
// They are implemented by the following general method to verify the condition Cond:
//	static char DummyArray[Cond ? +1 : -1];
// that would lead to a syntax error if the condition failed (by declaring a negative array size).
// You can remove them, if you want, since they are not used anywhere else.
#define ForceCond(Arr, Cond)	static char Arr[(Cond) ? +1 : -1]
ForceCond(_ae_bool_must_be_8_bits_wide, (int)sizeof(bool) == 1);
ForceCond(_ae_int32_t_must_be_32_bits_wide, (int)sizeof(ae_int32_t) == 4);
ForceCond(_ae_int64_t_must_be_64_bits_wide, (int)sizeof(ae_int64_t) == 8);
ForceCond(_ae_uint64_t_must_be_64_bits_wide, (int)sizeof(ae_uint64_t) == 8);
ForceCond(_ae_int_t_must_be_pointer_sized, (int)sizeof(ae_int_t) == (int)sizeof(void *));
#if defined ALGLIB_REDZONE
ForceCond(_ae_redzone_must_be_multiple_of_64, (ALGLIB_REDZONE) >= (AE_DATA_ALIGN) && (ALGLIB_REDZONE) % (AE_DATA_ALIGN) == 0);
#endif
#undef ForceCond

// Allocation tracking, for debugging.
ae_int_t _alloc_counter = 0;
ae_int_t _alloc_counter_total = 0;
bool _use_alloc_counter = false;

// Allocation debugging.
bool _force_malloc_failure = false;
ae_int_t _malloc_failure_after = 0;

//(@) Originally a part of the global environment structure, these should all be made thread-local.
// A pointer to the jmp_buf for cases when C-style exception handling is used.
// It may be NULL.
AutoS jmp_buf *volatile CurBreakAt;

// Set the jump buffer for error handling.
void ae_state_set_break_jump(jmp_buf *buf) { CurBreakAt = buf; }

// The ae_error_type of the last error and a legible message for it; filled when an exception is thrown.
AutoS ae_error_type volatile CurStatus;
AutoS const char *volatile CurMsg;

// Flags: call-local settings for ALGLIB++.
AutoS ae_uint64_t CurFlags;

// Set CurFlags.
void ae_state_set_flags(ae_uint64_t flags) { CurFlags = flags; }

// A pointer to the top block in a stack of frames which hold dynamically allocated objects.
AutoS ae_frame *volatile TopFr;

#if 0 //(@) Not used, but retained for possible inclusion into the multi-threading core being created anew for ALGLIB++.
// Threading information.
// NOTES:
// *	These are remnants from the Commercial Version, which worked with the (Commercial-only) file smp.h.
// *	They were declared as generic pointers of type (void *) in order to avoid explicit dependency on smp.h.
// *	The current thread pool.
// AutoS void *CurPool = NULL; //(@) Was never included in the Free Version.
// *	The current worker thread.
AutoS void *CurThread = NULL;
// *	The parent task: the one we are solving right now.
AutoS void *CurTask = NULL;
// *	The thread exception handler: the function which must be called by ae_assert() before raising an exception.
AutoS void (*ErrorOp)(void *) = NULL;
#endif

// The stack and frame boundary special blocks.
static unsigned char DynBottom = 1, DynFrame = 2;

// Make a new stack frame for the environment.
// Fr points to the place in the dynamic block that marks where the frame begins.
// The dynamic block is assumed to be initialized by the caller and must be left alone
// (no changes, deallocations, or reuse) until ae_leave_frame() is called.
// It may be a global or (preferrably) a local variable.
void ae_frame_make(ae_frame *Fr) { Fr->p_next = TopFr, Fr->deallocator = NULL, Fr->ptr = &DynFrame, TopFr = Fr; }

// Leave the current stack frame and deallocate all automatic dynamic blocks which were attached to this frame.
void ae_frame_leave() {
   for (; TopFr->ptr != &DynFrame && TopFr->ptr != &DynBottom; TopFr = TopFr->p_next)
      if (TopFr->ptr != NULL && TopFr->deallocator != NULL) TopFr->deallocator(TopFr->ptr);
   if (TopFr->ptr == &DynFrame) TopFr = TopFr->p_next;
}

// Initialize the ALGLIB++ frame stack environment.
// NOTE:
// *	Stacks contain no frames, so ae_make_frame() must be called before attaching dynamic blocks.
//	Without it ae_leave_frame() will cycle forever -- as intended.
void ae_state_init() {
// The base of the current stack of frames.
// p_next points to itself because a correct program should be able to detect the end of the list by looking at the ptr field.
// p_next == NULL may be used to distinguish automatic blocks (in the list) from non-automatic (not in the list).
   static ae_frame BotFr = { &BotFr, NULL, &DynBottom };
// Set the status indicators and clear the frame.
   CurBreakAt = NULL, CurMsg = "", CurFlags = NonTH, TopFr = &BotFr;
#if 0 //(@) Not used, but retained for possible inclusion into the multi-threading core being created anew for ALGLIB++.
// Set the threading information.
   CurThread = NULL, CurTask = NULL, ErrorOp = NULL;
#endif
}

// Clear the ALGLIB++ frame stack environment, freeing all the dynamic data in it that it controls.
void ae_state_clear() {
   if (TopFr == NULL) return;
   for (; TopFr->ptr != &DynBottom; TopFr = TopFr->p_next)
      if (TopFr->ptr != NULL && TopFr->deallocator != NULL) TopFr->deallocator(TopFr->ptr);
   TopFr = NULL;
}

// Clean up automatically managed memory before the caller terminates ALGLIB++.
// For TopFr != NULL call ErrorOp(), if defined.
// For TopFr == NULL do nothing.
void ae_clean_up() {
#if 0 //(@) Not used, but retained for possible inclusion into the multi-threading core being created anew for ALGLIB++.
   if (TopFr != NULL && ErrorOp != NULL) ErrorOp(TopFr);
#endif
}

// Abnormally abort the program, using one of several ways:
// *	if TopFr != NULL and CurBreakAt points to a jmp_buf - longjmp() to the return site,
// *	else abort().
// In all cases, for TopFr != NULL, set CurStatus and CurMsg.
// Clear the frame stack, if any, with ae_state_clear().
#if 0 //(@) Not used, but retained for possible inclusion into the multi-threading core being created anew for ALGLIB++.
// If TopFr != NULL and ErrorOp() != NULL, call ErrorOp() before handling errors and clearing TopFr.
#endif
static void ae_break(ae_error_type error_type, const char *msg) {
   if (TopFr == NULL) abort();
#if 0 //(@) Not used, but retained for possible inclusion into the multi-threading core being created anew for ALGLIB++.
   if (ErrorOp != NULL) ErrorOp(TopFr);
#endif
   ae_state_clear();
   CurStatus = error_type, CurMsg = msg;
   if (CurBreakAt != NULL) longjmp(*CurBreakAt, 1); else abort();
}

// Assertion.
// Upon failure with TopFr != NULL, gracefully exit ALGLIB++, removing all frames and deallocating registered dynamic data structures.
// Otherwise, just abort().
// IMPORTANT:
// *	This function ALWAYS evaluates cond, and cannot be replaced by macros which do nothing.
//	In particular, when invoked, a function call may be used as the cond argument, and it will be carried out.
void ae_assert(bool cond, const char *msg) {
   if (!cond) ae_break(ERR_ASSERTION_FAILED, msg);
}

#define AE_CRITICAL_ASSERT(x) if (!(x)) abort()

// Make flags variables into one or more char-sized variables in order to avoid problems with non-atomic reads/writes
// (single-byte ops are atomic on all contemporary CPUs).
#define _ALGLIB_FLG_THREADING_MASK	0x7
#define _ALGLIB_FLG_THREADING_SHIFT	0

static unsigned char _alglib_global_threading_flags = SerTH >> _ALGLIB_FLG_THREADING_SHIFT;

// Get/Set the default (global) threading model:
// *	serial execution,
// *	multithreading, if cores_to_use allows it.
ae_uint64_t ae_get_global_threading() {
   return (ae_uint64_t)_alglib_global_threading_flags << _ALGLIB_FLG_THREADING_SHIFT;
}

void ae_set_global_threading(ae_uint64_t flg_value) {
   flg_value &= _ALGLIB_FLG_THREADING_MASK;
   AE_CRITICAL_ASSERT(flg_value == SerTH || flg_value == ParTH);
   _alglib_global_threading_flags = (unsigned char)(flg_value >> _ALGLIB_FLG_THREADING_SHIFT);
}

// The recommended number of active workers:
// *	The exact number of cores to use, if AE_NWORKERS > 0,
// *	ALL available cores, if AE_NWORKERS == 0,
// *	max(_alglib_cores_to_use + AE_NWORKERS, 1), if AE_NWORKERS < 0,
// *	Default value == 0: either full parallelism if AE_NWORKERS is not defined,
//	or a manually-set number of cores if AE_NWORKERS is defined.
// PROTECTION:
// *	not needed; runtime modification is possible, but we do not need exact synchronization.
#if defined AE_NWORKERS && AE_NWORKERS <= 0
#   error AE_NWORKERS must be positive number or not defined at all.
#endif

static ae_int_t _alglib_cores_to_use = 0;

// CPUID
// Information about the features the CPU and compiler support.
// You must tell ALGLIB++ what CPU family is used by defining AE_CPU (without this hint zero will be returned).
// NOTE:
// *	The results of this function depend on both the CPU and compiler;
//	if the compiler doesn't support SSE intrinsics, then the function won't set the corresponding flag.
static ae_cpuid_t ae_cpuid() {
// Determine the CPU characteristics and perform CPU-specific initializations.
// Previous calls are cached to speed things up.
// There is no synchronization, but this can be safely done on a per-thread basis,
// provided that simultaneous writes by different cores to the same location will be executed in serial manner,
// which is true of contemporary architectures.
   static volatile bool _ae_cpuid_initialized = false, _ae_cpuid_has_sse2 = false, _ae_cpuid_has_avx2 = false, _ae_cpuid_has_fma = false;
// If not initialized, then determine the system properties.
   if (!_ae_cpuid_initialized) {
#if defined AE_CPU && AE_CPU == AE_INTEL
   { // SSE2
#   if AE_COMPILER == AE_GNUC || AE_COMPILER == AE_SUNC
      ae_int_t a, b, c, d;
#   elif AE_COMPILER == AE_MSVC
      int CPUInfo[4];
#   endif
   // SSE2 support.
#   if defined _ALGLIB_HAS_SSE2_INTRINSICS || AE_COMPILER == AE_SUNC
#      if AE_COMPILER == AE_GNUC || AE_COMPILER == AE_SUNC
      __asm__ __volatile__("cpuid":"=a"(a), "=b"(b), "=c"(c), "=d"(d):"a"(1));
      if (d & 0x04000000) _ae_cpuid_has_sse2 = true;
#      elif AE_COMPILER == AE_MSVC
      __cpuid(CPUInfo, 1);
      if (CPUInfo[3] & 0x04000000) _ae_cpuid_has_sse2 = true;
#      endif
#   endif
#   if defined _ALGLIB_HAS_AVX2_INTRINSICS
   // Check OS support for XSAVE XGETBV.
#      if AE_COMPILER == AE_GNUC
      __asm__ __volatile__("cpuid":"=a"(a), "=b"(b), "=c"(c), "=d"(d):"a"(1));
      if (c & (0x1 << 27)) {
         __asm__ volatile ("xgetbv":"=a" (a), "=d"(d):"c"(0));
         if ((a & 0x6) == 0x6) {
            if (_ae_cpuid_has_sse2) { // AVX2 support.
               __asm__ __volatile__("cpuid":"=a"(a), "=b"(b), "=c"(c), "=d"(d):"a"(7), "c"(0));
               if (b & (0x1 << 5)) _ae_cpuid_has_avx2 = true;
            }
#         if defined _ALGLIB_HAS_FMA_INTRINSICS
            if (_ae_cpuid_has_avx2) { // FMA support.
               __asm__ __volatile__("cpuid":"=a"(a), "=b"(b), "=c"(c), "=d"(d):"a"(1));
               if (c & (0x1 << 12)) _ae_cpuid_has_fma = true;
            }
#         endif
         }
      }
#      elif AE_COMPILER == AE_MSVC && _MSC_VER >= 1600
      __cpuid(CPUInfo, 1);
      if ((CPUInfo[2] & (0x1 << 27)) && (_xgetbv(0) & 0x6) == 0x6) {
         if (_ae_cpuid_has_sse2) { // AVX2 support.
            __cpuidex(CPUInfo, 7, 0);
            if (CPUInfo[1] & (0x1 << 5)) _ae_cpuid_has_avx2 = true;
         }
#         if defined _ALGLIB_HAS_FMA_INTRINSICS
         if (_ae_cpuid_has_avx2) { // FMA support.
            __cpuid(CPUInfo, 1);
            if (CPUInfo[2] & (0x1 << 12)) _ae_cpuid_has_fma = true;
         }
#         endif
      }
#      endif
#   endif
   } { // Perform one more CPUID call to generate memory fence
#   if AE_COMPILER == AE_GNUC || AE_COMPILER == AE_SUNC
      ae_int_t a, b, c, d;
      __asm__ __volatile__("cpuid":"=a"(a), "=b"(b), "=c"(c), "=d"(d):"a"(1));
#   elif AE_COMPILER == AE_MSVC
      int CPUInfo[4];
      __cpuid(CPUInfo, 1);
#   endif
   } { // Perform other CPU-related initializations.
#   if AE_COMPILER == AE_GNUC && 0
   // (@) Legacy code required for earlier versions of GCC, no longer needed here.
   // (@) TODO: If it still needed to be used, then make ModeCPU externally accessible, so that it can be used to reset the CPU.
   // Set rounding for floating-point math to double precision for x86/x64 processors under GCC.
      fp_control_t ModeCPU; _FPU_GETCW(ModeCPU);
      _FPU_SETCW(ModeCPU & ~(_FPU_EXTENDED | _FPU_SINGLE) | FPU_DOUBLE);
#   endif
   }
#endif
   // Set the initialization flag.
      _ae_cpuid_initialized = true;
   }
   return (ae_cpuid_t)(
      (_ae_cpuid_has_sse2 ? (int)CPU_SSE2 : 0) |
      (_ae_cpuid_has_avx2 ? (int)CPU_AVX2 : 0) |
      (_ae_cpuid_has_fma ? (int)CPU_FMA : 0)
   );
}

const/* AutoS */ae_cpuid_t CurCPU = ae_cpuid();

// The number of cores in the system: values < 1 may be returned on error.
ae_int_t ae_cores_count() {
#if AE_OS == AE_POSIX
   return sysconf(_SC_NPROCESSORS_ONLN);
#elif AE_OS == AE_WINDOWS
   SYSTEM_INFO sysInfo; GetSystemInfo(&sysInfo);
   return (ae_int_t)sysInfo.dwNumberOfProcessors;
#else
   return 1;
#endif
}

// The maximum concurrency on the given system, with the given compilation settings.
// Set to 1 on error or if fast kernels are disabled.
ae_int_t maxconcurrency() {
#if defined ALGLIB_NO_FAST_KERNELS
   return 1;
#else
   ae_int_t Cores = ae_cores_count(); return Cores < 0 ? 1 : Cores;
#endif
}

// Map the nworkers number (which can be positive, zero for "all cores" or negative, e.g. -1 meaning "all cores -1")
// to an "effective", strictly positive workers count.
//
// This is meant for use by debugging/testing code which tests different number of worker threads.
// It is NOT aligned in any way with the ALGLIB++ multithreading framework
// (i.e. it can return a non-zero worker count even for single-threaded GPLed ALGLIB++).
ae_int_t ae_get_effective_workers(ae_int_t nworkers) {
// Count the cores.
#if defined AE_NWORKERS
   ae_int_t ncores = AE_NWORKERS;
#else
   ae_int_t ncores = ae_cores_count();
#endif
   AE_CRITICAL_ASSERT(ncores >= 1);
// Map nworkers to its effective value.
   return nworkers >= 1 ? (nworkers > ncores ? ncores : nworkers) : (ncores + nworkers >= 1 ? ncores + nworkers : 1);
}

// Debug counters and flags.
static ae_int_t _dbg_alloc_total = 0;
static bool _use_dbg_counters = false;
static bool _use_vendor_kernels = true;
static bool debug_workstealing = false; // Debug work-stealing environment? false by default.
static ae_int_t dbgws_pushroot_ok = 0;
static ae_int_t dbgws_pushroot_failed = 0;

// Standard function wrappers for better GLIBC portability
#if defined X_FOR_LINUX
__asm__(".symver exp,exp@GLIBC_2.2.5");
__asm__(".symver log,log@GLIBC_2.2.5");
__asm__(".symver pow,pow@GLIBC_2.2.5");
double __wrap_exp(double x) { return exp(x); }
double __wrap_log(double x) { return log(x); }
double __wrap_pow(double x, double y) { return pow(x, y); }
#endif

ae_int64_t ae_get_dbg_value(debug_flag_t id) {
   switch (id) {
      case _ALGLIB_ALLOC_COUNTER: return _alloc_counter;
      case _ALGLIB_TOTAL_ALLOC_SIZE: return _dbg_alloc_total;
      case _ALGLIB_TOTAL_ALLOC_COUNT: return _alloc_counter_total;
#if defined AE_MKL
      case _ALGLIB_VENDOR_MEMSTAT: return ae_mkl_memstat();
#endif
   // Work-stealing counters.
#if defined AE_SMP
      case _ALGLIB_WSDBG_NCORES: return ae_cores_count();
#endif
      case _ALGLIB_WSDBG_PUSHROOT_OK: return dbgws_pushroot_ok;
      case _ALGLIB_WSDBG_PUSHROOT_FAILED: return dbgws_pushroot_failed;
#if defined AE_SMP
      case _ALGLIB_CORES_COUNT: return ae_cores_count();
#endif
      case _ALGLIB_GLOBAL_THREADING: return ae_get_global_threading();
      case _ALGLIB_NWORKERS: return _alglib_cores_to_use;
   // Unknown value.
      default: return 0;
   }
}

void ae_set_dbg_value(debug_flag_t flag_id, ae_int64_t flag_val) {
   switch (flag_id) {
      case _ALGLIB_ALLOC_COUNTER: _use_alloc_counter = flag_val != 0; break;
      case _ALGLIB_TOTAL_ALLOC_SIZE: _use_dbg_counters = flag_val != 0; break;
      case _ALGLIB_USE_VENDOR_KERNELS: _use_vendor_kernels = flag_val != 0; break;
      case _ALGLIB_DEBUG_WORKSTEALING: debug_workstealing = flag_val != 0; break;
      case _ALGLIB_GLOBAL_THREADING: ae_set_global_threading((ae_uint64_t)flag_val); break;
      case _ALGLIB_NWORKERS: _alglib_cores_to_use = (ae_int_t)flag_val; break;
   }
}

// A wrapper around OS-dependent clock routines.
#if AE_OS == AE_POSIX || defined AE_DEBUG4POSIX
ae_int_t tickcount() {
   struct timeval now;
   ae_int64_t r, v;
   gettimeofday(&now, NULL);
   v = now.tv_sec;
   r = v * 1000;
   v = now.tv_usec / (suseconds_t)1000;
   r += v;
   return (ae_int_t)r;
#   if 0
   struct timespec now;
   return clock_gettime(CLOCK_MONOTONIC, &now) ? 0 : now.tv_sec * 1000.0 + now.tv_nsec / 1000000.0;
#   endif
}
#elif AE_OS == AE_WINDOWS || defined AE_DEBUG4WINDOWS
ae_int_t tickcount() {
   return (ae_int_t)GetTickCount();
}
#else
ae_int_t tickcount() {
   return 0;
}
#endif

ae_int_t ae_misalignment(const void *ptr, size_t alignment) {
   union {
      const void *ptr;
      ae_int_t iptr;
   } u;
   u.ptr = ptr;
   return u.iptr % (ae_int_t)alignment;
}

void *ae_align(void *ptr, size_t alignment) {
   char *result = (char *)ptr;
   if ((result - (char *)NULL) % alignment != 0)
      result += alignment - (result - (char *)NULL) % alignment;
   return result;
}

// The "optional atomics" functions:
// functions which either perform atomic changes - or fall back to an ordinary operation,
// if the current compiler settings cannot generate atomic code.
// They are all synchronized, i.e. either all of them work - or none of them do.

// Perform atomic addition on a pointer-sized and pointer-size aligned value.
// NOTE:
// *	This function is not keyed in for high performance, so use it only when necessary.
static void ae_optional_atomic_add_i(ae_int_t *p, ae_int_t v) {
   AE_CRITICAL_ASSERT(ae_misalignment(p, sizeof(void *)) == 0);
#if AE_COMPILER == AE_GNUC && AE_CPU == AE_INTEL && 100 * __GNUC__ + __GNUC__ >= 470
   __atomic_add_fetch(p, v, __ATOMIC_RELAXED);
#elif defined __clang__ && AE_CPU == AE_INTEL
   __atomic_fetch_add(p, v, __ATOMIC_RELAXED);
#elif AE_OS == AE_WINDOWS
   while (true) {
   // Convert between ae_int_t * and void ** without compiler warnings about indirection levels.
      union {
         PVOID volatile *volatile ptr;
         volatile ae_int_t *volatile iptr;
      } u;
      u.iptr = p;
   // Atomic read the initial value, convert it to a 1-byte pointer, then increment and store it.
      PVOID v0 = InterlockedCompareExchangePointer(u.ptr, NULL, NULL);
      if (InterlockedCompareExchangePointer(u.ptr, (PVOID)((char *)v0 + v), v0) == v0) break;
   }
#else
   *p += v; // At least do something for older compilers!
#endif
}

// Perform atomic subtraction on a pointer-sized and pointer-size aligned value.
// NOTE:
// *	This function is not keyed in for high performance, so use it only when necessary.
static void ae_optional_atomic_sub_i(ae_int_t *p, ae_int_t v) {
   AE_CRITICAL_ASSERT(ae_misalignment(p, sizeof(void *)) == 0);
#if AE_COMPILER == AE_GNUC && AE_CPU == AE_INTEL && 100 * __GNUC__ + __GNUC__ >= 470
   __atomic_sub_fetch(p, v, __ATOMIC_RELAXED);
#elif defined __clang__ && AE_CPU == AE_INTEL
   __atomic_fetch_sub(p, v, __ATOMIC_RELAXED);
#elif AE_OS == AE_WINDOWS
   while (true) {
   // Convert between ae_int_t * and void ** without compiler warnings about indirection levels.
      union {
         PVOID volatile *volatile ptr;
         volatile ae_int_t *volatile iptr;
      } u;
      u.iptr = p;
   // Atomic read the initial value, convert it to a 1-byte pointer, then decrement and store it.
      PVOID v0 = InterlockedCompareExchangePointer(u.ptr, NULL, NULL);
      if (InterlockedCompareExchangePointer(u.ptr, (PVOID)((char *)v0 - v), v0) == v0) break;
   }
#else
   *p -= v; // At least do something for older compilers!
#endif
}

#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
// Fields for memory allocation over static array.
#   if AE_THREADING != NonTH
#      error Basic static malloc is thread-unsafe; define AE_THREADING = NonTH to prove that you know it.
#   endif

static ae_int_t sm_page_size = 0;
static ae_int_t sm_page_cnt = 0;
static ae_int_t *sm_page_tbl = NULL;
static unsigned char *sm_mem = NULL;

void set_memory_pool(void *ptr, size_t size) {
// Integrity checks.
   AE_CRITICAL_ASSERT(sm_page_size == 0);
   AE_CRITICAL_ASSERT(sm_page_cnt == 0);
   AE_CRITICAL_ASSERT(sm_page_tbl == NULL);
   AE_CRITICAL_ASSERT(sm_mem == NULL);
   AE_CRITICAL_ASSERT(size > 0);
// Align the pointer.
   size -= ae_misalignment(ptr, sizeof *sm_page_tbl);
   ptr = ae_align(ptr, sizeof *sm_page_tbl);
// Calculate the page size and page count, prepare pointers to the page table and memory.
   sm_page_size = 0x100;
// We expect to have memory for at least one page + table entry + alignment.
   AE_CRITICAL_ASSERT(size >= (sm_page_size + sizeof *sm_page_tbl) + sm_page_size);
   sm_page_cnt = (size - sm_page_size) / (sm_page_size + sizeof *sm_page_tbl);
   AE_CRITICAL_ASSERT(sm_page_cnt > 0);
   sm_page_tbl = (ae_int_t *)ptr;
   sm_mem = (unsigned char *)ae_align(sm_page_tbl + sm_page_cnt, sm_page_size);
// Mark all pages as free.
   memset(sm_page_tbl, 0, sm_page_cnt * sizeof *sm_page_tbl);
}

static void *ae_static_malloc(size_t size, size_t alignment) {
   AE_CRITICAL_ASSERT(size >= 0);
   AE_CRITICAL_ASSERT(sm_page_size > 0);
   AE_CRITICAL_ASSERT(sm_page_cnt > 0);
   AE_CRITICAL_ASSERT(sm_page_tbl != NULL);
   AE_CRITICAL_ASSERT(sm_mem != NULL);
   if (size == 0 || _force_malloc_failure) return NULL;
// Check that the page alignment and requested alignment match each other.
   AE_CRITICAL_ASSERT(alignment <= sm_page_size);
   AE_CRITICAL_ASSERT((sm_page_size % alignment) == 0);
// Search a long enough sequence of pages.
   int rq_pages = size / sm_page_size;
   if (size % sm_page_size) rq_pages++;
   int cur_len = 0;
   for (int i = 0; i < sm_page_cnt; ) {
   // Determine the length of the free page sequence.
      if (sm_page_tbl[i] == 0) cur_len++;
      else {
         AE_CRITICAL_ASSERT(sm_page_tbl[i] > 0);
         cur_len = 0;
         i += sm_page_tbl[i];
         continue;
      }
   // Found it?
      if (cur_len >= rq_pages) {
      // Update whichever counters the use-flags are set for.
         if (_use_alloc_counter) {
            ae_optional_atomic_add_i(&_alloc_counter, 1);
            ae_optional_atomic_add_i(&_alloc_counter_total, 1);
         }
         if (_use_dbg_counters) ae_optional_atomic_add_i(&_dbg_alloc_total, size);
      // Mark pages and return.
         for (int j = 0; j < rq_pages; j++) sm_page_tbl[i - j] = -1;
         sm_page_tbl[i - (rq_pages - 1)] = rq_pages;
         return sm_mem + (i - (rq_pages - 1)) * sm_page_size;
      }
   // The next page.
      i++;
   }
   return NULL;
}

static void ae_static_free(void *block) {
   if (block == NULL) return;
   ae_int_t page_idx = (unsigned char *)block - sm_mem;
   AE_CRITICAL_ASSERT(page_idx >= 0);
   AE_CRITICAL_ASSERT((page_idx % sm_page_size) == 0);
   page_idx /= sm_page_size;
   AE_CRITICAL_ASSERT(page_idx < sm_page_cnt);
   ae_int_t page_cnt = sm_page_tbl[page_idx];
   AE_CRITICAL_ASSERT(page_cnt >= 1);
   for (ae_int_t i = 0; i < page_cnt; i++) sm_page_tbl[page_idx + i] = 0;
// Update the counters (if the use-flag is set).
   if (_use_alloc_counter) ae_optional_atomic_sub_i(&_alloc_counter, 1);
}

void memory_pool_stats(ae_int_t *bytes_used, ae_int_t *bytes_free) {
   AE_CRITICAL_ASSERT(sm_page_size > 0);
   AE_CRITICAL_ASSERT(sm_page_cnt > 0);
   AE_CRITICAL_ASSERT(sm_page_tbl != NULL);
   AE_CRITICAL_ASSERT(sm_mem != NULL);
// Scan the page table.
   *bytes_free = *bytes_used = 0;
   for (int i = 0; i < sm_page_cnt; ) {
      if (sm_page_tbl[i] == 0) {
         ++*bytes_free;
         i++;
      } else {
         AE_CRITICAL_ASSERT(sm_page_tbl[i] > 0);
         *bytes_used += sm_page_tbl[i];
         i += sm_page_tbl[i];
      }
   }
   *bytes_used *= sm_page_size;
   *bytes_free *= sm_page_size;
}
#endif

void *aligned_malloc(size_t size, size_t alignment) {
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
   return ae_static_malloc(size, alignment);
#else
   if (size == 0 || _force_malloc_failure || _malloc_failure_after > 0 && _alloc_counter_total >= _malloc_failure_after) return NULL;
// Allocate, making the appropriate padding adjustments for any alignment > 1.
   size_t alloc_size = 2 * sizeof(void *) + size;
   if (alignment > 1) alloc_size += alignment - 1;
#   if defined ALGLIB_REDZONE
   alloc_size += 2 * (ALGLIB_REDZONE);
#   endif
   void *block = malloc(alloc_size); if (block == NULL) return NULL;
   char *result = (char *)block + 2 * sizeof block;
   result = (char *)ae_align(result, alignment);
   *(void **)(result - sizeof block) = block;
#   if defined ALGLIB_REDZONE
   char *redzone0 = result;
   result = redzone0 + (ALGLIB_REDZONE);
   char *redzone1 = result + size;
   ae_assert(ae_misalignment(result, alignment) == 0, "ALGLIB: improperly configured red zone size - is not multiple of the current alignment");
   *(void **)(redzone0 - 2 * sizeof block) = redzone1;
   memset(redzone0, _ALGLIB_REDZONE_VAL, ALGLIB_REDZONE);
   memset(redzone1, _ALGLIB_REDZONE_VAL, ALGLIB_REDZONE);
#   endif
// Update whichever counters the use-flags are set for.
   if (_use_alloc_counter) {
      ae_optional_atomic_add_i(&_alloc_counter, 1);
      ae_optional_atomic_add_i(&_alloc_counter_total, 1);
   }
   if (_use_dbg_counters) ae_optional_atomic_add_i(&_dbg_alloc_total, (ae_int_t)size);
// Return the result as a generic pointer.
   return (void *)result;
#endif
}

static void *aligned_extract_ptr(void *block) {
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
   return NULL;
#elif defined ALGLIB_REDZONE
   return block == NULL ? NULL : *(void **)((char *)block - (ALGLIB_REDZONE) - sizeof block);
#else
   return block == NULL ? NULL : *(void **)((char *)block - sizeof block);
#endif
}

void aligned_free(void *block) {
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
   ae_static_free(block);
#else
// Handle NULL input.
   if (block == NULL) return;
// If red zone is activated, check it before deallocation.
#   if defined ALGLIB_REDZONE
   char *redzone0 = (char *)block - (ALGLIB_REDZONE);
   char *redzone1 = (char *)*(void **)(redzone0 - 2 * sizeof block);
   for (ae_int_t i = 0; i < (ALGLIB_REDZONE); i++) {
      if (redzone0[i] != _ALGLIB_REDZONE_VAL) {
         const char *msg = "ALGLIB: red zone corruption is detected (write prior to the block beginning?)";
         fprintf(stderr, "%s\n", msg);
         ae_assert(false, msg, NULL);
      }
      if (redzone1[i] != _ALGLIB_REDZONE_VAL) {
         const char *msg = "ALGLIB: red zone corruption is detected (write past the end of the block?)";
         fprintf(stderr, "%s\n", msg);
         ae_assert(false, msg, NULL);
      }
   }
#   endif
// Free the memory and optionally update allocation counters.
   free(aligned_extract_ptr(block));
   if (_use_alloc_counter) ae_optional_atomic_sub_i(&_alloc_counter, 1);
#endif
}

// Allocate memory with automatic alignment.
// Return NULL when size == 0 is specified.
// Upon failure with TopFr != NULL, call ae_break(), otherwise return NULL.
void *ae_malloc(size_t size) {
   if (size == 0) return NULL;
   void *result = aligned_malloc(size, AE_DATA_ALIGN);
   if (result == NULL && TopFr != NULL) ae_break(ERR_OUT_OF_MEMORY, "ae_malloc: out of memory");
   return result;
}

// Allocate memory with automatic alignment and zero-fill.
// Returns NULL when zero size is specified.
// Error handling:
// *	if TopFr == NULL, return NULL on allocation error,
// *	if TopFr != NULL, call ae_break() on allocation error.
static void *ae_malloc_zero(size_t size) {
   void *result = ae_malloc(size);
   if (result != NULL)
      memset(result, 0, size);
   return result;
}

void ae_free(void *p) {
   if (p != NULL) aligned_free(p);
}

// Attach block to the dynamic block list for the ALGLIB++ environment.
// This function does NOT generate exceptions.
// NOTE:
// *	Avoid calling it for the special blocks which mark frame boundaries!
static void ae_db_attach(ae_dyn_block *block) { block->p_next = TopFr, TopFr = block; }

// Allocate and initialize a dynamic block of size >= 0 bytes for the ALGLIB++ environment.
// It is assumed to be uninitialized, its fields are ignored.
// make_automatic indicates that the block is to be added to the dynamic block list.
// Upon allocation failure with TopFr != NULL, call ae_break(), leaving block in a valid (but empty) state.
// NOTES:
// *	Avoid calling it for blocks which are already in the list.
//	Use ae_db_realloc() for already-allocated blocks.
// *	No memory allocation is performed for initialization with size == 0.
void ae_db_init(ae_dyn_block *block, size_t size, bool make_automatic) {
//(@) TopFr != NULL check and zero-check removed.
// NOTE:
// *	These strange dances around block->ptr are necessary in order to correctly handle possible exceptions during memory allocation.
   ae_assert(size >= 0, "ae_db_init: negative size");
   block->ptr = NULL;
   if (make_automatic) ae_db_attach(block); else block->p_next = NULL;
   if (size != 0) block->ptr = ae_malloc(size);
   block->deallocator = ae_free;
}

// Reallocate the dynamic block (assumed to be initialized) to size bytes for the ALGLIB++ environment.
// Delete the old contents but preserve the automatic state.
// Upon allocation failure with TopFr != NULL, call ae_break(), leaving block in a valid (but empty) state.
// NOTE:
// *	Avoid calling it for the special blocks which mark frame boundaries!
void ae_db_realloc(ae_dyn_block *block, ae_int_t size) {
//(@) TopFr != NULL check removed.
// NOTE:
// *	These strange dances around block->ptr are necessary in order to correctly handle possible exceptions during memory allocation.
   ae_assert(size >= 0, "ae_db_realloc: negative size");
   if (block->ptr != NULL) block->deallocator(block->ptr), block->ptr = NULL;
   block->ptr = ae_malloc((size_t)size);
   block->deallocator = ae_free;
}

// Clear the dynamic block (assumed to be initialized), releasing all dynamically allocated memory.
// The dynamic block may be in the automatic management list - in this case it will NOT be removed from the list.
// NOTE:
// *	Avoid calling it for the special blocks which mark frame boundaries!
void ae_db_free(ae_dyn_block *block) {
   if (block->ptr != NULL) block->deallocator(block->ptr), block->ptr = NULL;
   block->deallocator = ae_free;
}

// Swap dynamic blocks block1 and block2 (pointers and deallocators)
// leaving other parameters (automatic management settings, etc.) unchanged.
// NOTE:
// *	Avoid calling it for the special blocks which mark frame boundaries!
void ae_db_swap(ae_dyn_block *block1, ae_dyn_block *block2) {
   void *volatile ptr = block1->ptr;
   ae_deallocator deallocator = block1->deallocator;
   block1->ptr = block2->ptr;
   block1->deallocator = block2->deallocator;
   block2->ptr = ptr;
   block2->deallocator = deallocator;
}

// The size of datatype or zero for dynamic types like strings or multiple precision types.
ae_int_t ae_sizeof(ae_datatype datatype) {
   switch (datatype) {
   // case DT_BYTE: // The same as DT_BOOL.
      case DT_BOOL: return (ae_int_t)sizeof(bool);
      case DT_INT: return (ae_int_t)sizeof(ae_int_t);
      case DT_REAL: return (ae_int_t)sizeof(double);
      case DT_COMPLEX: return 2 * (ae_int_t)sizeof(double);
      default: return 0;
   }
}

// Make dst into a new datatype ae_vector of size >= 0.
// Its contents are assumed to be uninitialized, and its fields are ignored.
// make_automatic indicates whether or not the vector is to be added to the dynamic block list.
// Upon allocation failure or size < 0, call ae_break().
// NOTE:
// *	No memory allocation is performed for initialization with size == 0.
void ae_vector_init(ae_vector *dst, ae_int_t size, ae_datatype datatype, bool make_automatic) {
// Integrity checks.
//(@) TopFr != NULL check and zero-check removed.
   ae_assert(size >= 0, "ae_vector_init: negative size");
// Prepare for possible errors during allocation.
   dst->cnt = 0;
   dst->xX = NULL;
// Initialize.
   ae_db_init(&dst->data, (size_t)(size * ae_sizeof(datatype)), make_automatic);
   dst->cnt = size;
   dst->datatype = datatype;
   dst->xX = dst->data.ptr;
   dst->is_attached = false;
}

// Copy ae_vector src into ae_vector dst.
// dst is assumed to be uninitialized, its fields are ignored.
// The fields copied to dst are to be managed and owned by dst.
// make_automatic indicates whether or not the vector is to be added to the dynamic block list.
// Upon allocation failure, call ae_break().
void ae_vector_copy(ae_vector *dst, const ae_vector *src, bool make_automatic) {
//(@) TopFr != NULL check removed.
   ae_vector_init(dst, src->cnt, src->datatype, make_automatic);
   if (src->cnt != 0) memmove(dst->xX, src->xX, (size_t)(src->cnt * ae_sizeof(src->datatype)));
}

// Resize the ae_vector dst to size newsize >= 0.
// dst must be initialized.
// Its contents are freed by setlength().
// Upon allocation failure with TopFr != NULL, call ae_break(), otherwise return an indication of success or failure.
void ae_vector_set_length(ae_vector *dst, ae_int_t newsize) {
//(@) TopFr != NULL check removed.
   ae_assert(newsize >= 0, "ae_vector_set_length: negative size");
   if (dst->cnt == newsize) return;
// Reallocate, preparing first for possible errors.
   dst->cnt = 0;
   dst->xX = NULL;
   ae_db_realloc(&dst->data, newsize * ae_sizeof(dst->datatype));
   dst->cnt = newsize;
   dst->xX = dst->data.ptr;
}

// Resize the ae_vector dst to size newsize >= 0, preserving previously existing elements.
// dst must be initialized.
// The values of elements added during vector growth are undefined.
// Upon allocation failure, call ae_break().
void ae_vector_resize(ae_vector *dst, ae_int_t newsize) {
   ae_vector tmp; memset(&tmp, 0, sizeof tmp), ae_vector_init(&tmp, newsize, dst->datatype, false);
   ae_int_t bytes_total = (dst->cnt < newsize ? dst->cnt : newsize) * ae_sizeof(dst->datatype);
   if (bytes_total > 0) memmove(tmp.xX, dst->xX, bytes_total);
   ae_swap_vectors(dst, &tmp);
   ae_vector_free(&tmp, true);
}

// The "FREE" functionality for vector dst (cleared contents and freeing all internal structures).
// Clear vector dst (releasing all dynamically allocated memory).
// dst may be on the frame - in which case it will NOT be removed from the list.
// IMPORTANT:
// *	This function does NOT invalidate dst; it just releases all dynamically allocated storage,
//	but dst still may be used after calling ae_vector_set_length().
void ae_vector_free(ae_vector *dst, bool/* make_automatic*/) {
   dst->cnt = 0;
   ae_db_free(&dst->data);
   dst->xX = 0;
   dst->is_attached = false;
}

// Efficiently swap ae_vector vec1 with ae_vector vec2, leaving other pararemeters (automatic management, etc.) intact.
void ae_swap_vectors(ae_vector *vec1, ae_vector *vec2) {
   ae_assert(!vec1->is_attached, "ae_swap_vectors: internal error, attempt to swap vectors attached to X-object");
   ae_assert(!vec2->is_attached, "ae_swap_vectors: internal error, attempt to swap vectors attached to X-object");
   ae_db_swap(&vec1->data, &vec2->data);
   ae_int_t cnt = vec1->cnt;
   ae_datatype datatype = vec1->datatype;
   void *p_ptr = vec1->xX;
   vec1->cnt = vec2->cnt;
   vec1->datatype = vec2->datatype;
   vec1->xX = vec2->xX;
   vec2->cnt = cnt;
   vec2->datatype = datatype;
   vec2->xX = p_ptr;
}

// Lay out the raster for matrix dst from storage.
// *	dst must be a correctly initialized matrix.
// *	dst->data.ptr points to the beginning of memory block allocated for row pointers.
// *	dst->ptr - undefined (initialized during algorithm processing).
// *	storage points to the beginning of actual storage.
static void ae_matrix_update_row_pointers(ae_matrix *dst, void *storage) {
   if (dst->cols > 0 && dst->rows > 0) {
      char *p_base = (char *)storage;
      void **pp_ptr = (void **)dst->data.ptr;
      dst->xyX = pp_ptr;
      for (ae_int_t i = 0; i < dst->rows; i++, p_base += dst->stride * ae_sizeof(dst->datatype))
         pp_ptr[i] = p_base;
   } else dst->xyX = NULL;
}

// Make dst into a new rows x cols datatype ae_matrix.
// The matrix size may be zero, in such cases both cols and rows will be zero.
// Its contents are assumed to be uninitialized, and its fields are ignored.
// make_automatic indicates whether or not the matrix is to be added to the dynamic block list,
// as opposed to being a global object or field of some other object.
// Upon allocation failure or cols < 0 or rows < 0, call ae_break().
// NOTE:
// *	No memory allocation is performed for initialization with cols == 0 or rows == 0.
void ae_matrix_init(ae_matrix *dst, ae_int_t rows, ae_int_t cols, ae_datatype datatype, bool make_automatic) {
//(@) TopFr != NULL check and zero-check removed.
   ae_assert(cols >= 0 && rows >= 0, "ae_matrix_init: negative length");
// If either cols or rows is 0, then they both must be made so.
   if (cols == 0) rows = 0; else if (rows == 0) cols = 0;
// Initialize.
   dst->datatype = datatype;
   dst->stride = cols;
// Prepare for possible errors during allocation.
   dst->rows = dst->cols = 0;
   dst->xyX = NULL;
   dst->is_attached = false;
// If cols and rows are 0; perform a quick exit.
   if (cols == 0 || rows == 0) { ae_db_init(&dst->data, 0, make_automatic); return; }
// Initialize, preparing for possible errors during allocation.
   for (; dst->stride * ae_sizeof(datatype) % AE_DATA_ALIGN != 0; dst->stride++);
   ae_db_init(&dst->data, (size_t)(rows * (sizeof(void *) + dst->stride * ae_sizeof(datatype)) + AE_DATA_ALIGN - 1), make_automatic);
   dst->cols = cols;
   dst->rows = rows;
// Set the pointers to the matrix rows.
   ae_matrix_update_row_pointers(dst, ae_align((char *)dst->data.ptr + rows * sizeof(void *), AE_DATA_ALIGN));
}

// Copy ae_matrix src to ae_matrix dst.
// dst is assumed to be uninitialized, its fields are ignored.
// make_automatic indicates whether or not dst is to be added to the dynamic block list,
// as opposed to being a global object or field of some other object.
// Upon allocation failure, call ae_break().
void ae_matrix_copy(ae_matrix *dst, const ae_matrix *src, bool make_automatic) {
   ae_matrix_init(dst, src->rows, src->cols, src->datatype, make_automatic);
   if (src->cols > 0 && src->rows > 0)
      if (dst->stride == src->stride)
         memmove(dst->xyX[0], src->xyX[0], (size_t)(src->rows * src->stride * ae_sizeof(src->datatype)));
      else for (ae_int_t i = 0; i < dst->rows; i++)
         memmove(dst->xyX[i], src->xyX[i], (size_t)(dst->cols * ae_sizeof(dst->datatype)));
}

// Resize ae_matrix dst to size rows x cols.
// Either cols, rows or both may be 0.
// The matrix dst must be initialized.
// Its contents are freed after setlength().
// Upon allocation failure with TopFr != NULL, call ae_break(), otherwise return an indication of success or failure.
void ae_matrix_set_length(ae_matrix *dst, ae_int_t rows, ae_int_t cols) {
//(@) TopFr != NULL check removed.
   ae_assert(cols >= 0 && rows >= 0, "ae_matrix_set_length: negative length");
   if (dst->cols == cols && dst->rows == rows) return;
// Prepare the stride.
   for (dst->stride = cols; dst->stride * ae_sizeof(dst->datatype) % AE_DATA_ALIGN != 0; dst->stride++);
// Prepare for possible errors during reallocation.
   dst->rows = dst->cols = 0;
   dst->xyX = NULL;
   ae_db_realloc(&dst->data, rows * ((ae_int_t)sizeof(void *) + dst->stride * ae_sizeof(dst->datatype)) + AE_DATA_ALIGN - 1);
   dst->cols = cols;
   dst->rows = rows;
// Set the pointers to the matrix rows.
   ae_matrix_update_row_pointers(dst, ae_align((char *)dst->data.ptr + dst->rows * sizeof(void *), AE_DATA_ALIGN));
}

// The "FREE" functionality for ae_matrix dst.
// Clear the contents of matrix dst, releasing all dynamically allocated memory, but leaving the structure intact in a valid state.
// dst may be on the frame - in which case it will NOT be removed from the list.
// IMPORTANT:
// *	This function does NOT invalidate dst; it just releases all dynamically allocated storage,
//	but dst still may be used after calling ae_matrix_set_length().
void ae_matrix_free(ae_matrix *dst, bool/* make_automatic*/) {
   dst->rows = dst->cols = 0;
   dst->stride = 0;
   ae_db_free(&dst->data);
   dst->xX = 0;
   dst->is_attached = false;
}

// Efficiently swap ae_matrix mat1 with ae_matrix mat2, leaving other pararemeters (automatic management, etc.) intact.
void ae_swap_matrices(ae_matrix *mat1, ae_matrix *mat2) {
   ae_assert(!mat1->is_attached, "ae_swap_matrices: internal error, attempt to swap matrices attached to X-object");
   ae_assert(!mat2->is_attached, "ae_swap_matrices: internal error, attempt to swap matrices attached to X-object");
   ae_db_swap(&mat1->data, &mat2->data);
   ae_int_t cols = mat1->cols;
   ae_int_t rows = mat1->rows;
   ae_int_t stride = mat1->stride;
   ae_datatype datatype = mat1->datatype;
   void *p_ptr = mat1->xX;
   mat1->cols = mat2->cols;
   mat1->rows = mat2->rows;
   mat1->stride = mat2->stride;
   mat1->datatype = mat2->datatype;
   mat1->xX = mat2->xX;
   mat2->cols = cols;
   mat2->rows = rows;
   mat2->stride = stride;
   mat2->datatype = datatype;
   mat2->xX = p_ptr;
}

// Make dst into a new smart pointer.
// dst is assumed to be uninitialized, but preallocated, by aliasing subscriber (which may be NULL) with dst->ptr.
// After initialization, dst stores the NULL pointer.
// make_automatic indicates whether or not dst is to be added to the dynamic block list.
// Upon allocation failure, call ae_break().
void ae_smart_ptr_init(ae_smart_ptr *dst, void **subscriber, bool make_automatic) {
//(@) TopFr != NULL check and zero-check removed.
   dst->subscriber = subscriber;
   dst->ptr = NULL;
   if (dst->subscriber != NULL) *dst->subscriber = dst->ptr;
   dst->is_owner = false;
   dst->is_dynamic = false;
   dst->size_of_object = 0;
   dst->copy = NULL;
   dst->free = NULL;
   dst->frame_entry.deallocator = ae_smart_ptr_free;
   dst->frame_entry.ptr = dst;
   if (make_automatic) ae_db_attach(&dst->frame_entry);
}

// Free the smart pointer _dst so that it contains the NULL reference.
// The change is propagated to its subscriber, if _dst was created with a non-NULL subscriber.
void ae_smart_ptr_free(void *_dst) {
   ae_smart_ptr *dst = (ae_smart_ptr *)_dst;
   if (dst->is_owner && dst->ptr != NULL) {
      dst->free(dst->ptr, false);
      if (dst->is_dynamic) ae_free(dst->ptr);
   }
   dst->is_owner = false;
   dst->is_dynamic = false;
   dst->ptr = NULL;
   dst->size_of_object = 0;
   dst->copy = NULL;
   dst->free = NULL;
   if (dst->subscriber != NULL) *dst->subscriber = NULL;
}

// Assign pointer new_ptr to smart pointer dst.
// Any non-NULL value already contained in and owned by dst is freed beforehand.
// The change is propagated to its subscriber, if dst was created with a non-NULL subscriber.
// is_owner indicates whether dst is to own new_ptr.
// obj_size is the in-memory size of the object. Ignored for is_owner == false.
// copy is the function used to copy it; can not be NULL for new_ptr != NULL.
// free is the function used to free it; can be NULL for is_owner == false.
// is_dynamic indicates whether dst is to be dynamic,
// so that clearing dst would require BOTH calling ae_free() on it AND free() on the memory occupied by dst.
// You can specify NULL new_ptr, in which case is_owner, free() and is_dynamic are all ignored.
void ae_smart_ptr_assign(ae_smart_ptr *dst, void *new_ptr, bool is_owner, bool is_dynamic, size_t obj_size, ae_copy_op copy, ae_free_op free) {
   ae_assert(new_ptr == NULL || !is_owner || copy != NULL, "ae_smart_ptr_assign: new_ptr != NULL, is_owner, but copy constructor is NULL");
   ae_assert(new_ptr == NULL || !is_owner || free != NULL, "ae_smart_ptr_assign: new_ptr != NULL, is_owner, but destructor is NULL");
   ae_assert(new_ptr == NULL || !is_owner || obj_size > 0, "ae_smart_ptr_assign: new_ptr != NULL, is_owner, but object size is zero");
   if (dst->is_owner && dst->ptr != NULL) {
      dst->free(dst->ptr, false);
      if (dst->is_dynamic) ae_free(dst->ptr);
   }
   bool not_null = new_ptr != NULL;
   dst->ptr = new_ptr;
   dst->is_owner = not_null && is_owner;
   dst->is_dynamic = not_null && is_dynamic;
   dst->size_of_object = not_null && is_owner ? obj_size : 0;
   dst->copy = not_null ? copy : NULL;
   dst->free = not_null ? free : NULL;
   if (dst->subscriber != NULL) *dst->subscriber = dst->ptr;
}

// Release the pointer owned by the smart pointer dst by NULLing all internal fields
// and passing any ownership it has to the caller, instead of applying the destructor function to the internal pointer.
// The change is propagated to its subscriber, if the smart pointer was created with subscriber != NULL.
void ae_smart_ptr_release(ae_smart_ptr *dst) {
   dst->is_owner = false;
   dst->is_dynamic = false;
   dst->ptr = NULL;
   dst->size_of_object = 0;
   dst->copy = NULL;
   dst->free = NULL;
   if (dst->subscriber != NULL) *dst->subscriber = NULL;
}

// Copy x_vector src into ae_vector dst.
// dst is assumed to be uninitialized, its fields are ignored.
// The newly-created copy of src is to be owned/managed by dst.
// Both src and dst remain completely independent afterwards.
// make_automatic indicates whether or not the vector will be registered in the ALGLIB++ environment.
void ae_vector_init_from_x(ae_vector *dst, x_vector *src, bool make_automatic) {
//(@) TopFr != NULL check removed.
   ae_vector_init(dst, (ae_int_t)src->cnt, (ae_datatype)src->datatype, make_automatic);
   if (src->cnt > 0) memmove(dst->xX, src->x_ptr, (size_t)((ae_int_t)src->cnt * ae_sizeof((ae_datatype)src->datatype)));
}

// Copy x_vector src into ae_vector dst by attaching dst to src.
// dst is assumed to be uninitialized, its fields are ignored.
// make_automatic indicates whether or not the vector will be registered in the ALGLIB++ environment.
// The new vector is attached to the source:
// *	dst shares memory with src: a write to one changes both.
// *	dst can be reallocated with ae_vector_set_length(), but src remains untouched.
// *	src, however, can NOT be reallocated as long as dst exists.
// *	dst->is_attached is set to true to indicate that dst does not own its memory.
// Upon allocation failure, call ae_break().
void ae_vector_init_attach_to_x(ae_vector *dst, x_vector *src, bool make_automatic) {
//(@) TopFr != NULL check and zero-check removed.
   volatile ae_int_t cnt = (ae_int_t)src->cnt;
// Ensure the correct size.
   ae_assert(cnt == src->cnt, "ae_vector_init_attach_to_x: 32/64 overflow");
   ae_assert(cnt >= 0, "ae_vector_init_attach_to_x: negative length");
// Prepare for possible errors during allocation.
   dst->cnt = 0;
   dst->xX = NULL;
   dst->datatype = (ae_datatype)src->datatype;
// Zero-size initialize in order to correctly register in the frame.
   ae_db_init(&dst->data, 0, make_automatic);
// Initialize.
   dst->cnt = cnt;
   dst->xX = src->x_ptr;
   dst->is_attached = true;
}

// Copy ae_vector src to x_vector dst.
// Not meant for use when dst is attached to src, though it may be used when src is attached to dst.
// One of the following is then applied:
// *	if src is attached to dst, no action is required or done,
// *	for independent vectors of different size: allocate storage in dst and copy src to dst.
//	dst->last_action is set to ACT_NEW_LOCATION, and dst->owner is set to true.
// *	for independent vectors of the same size: no (re)allocation is required or done.
//	Just copy src to the already-existing place.
//	dst->last_action is set to ACT_SAME_LOCATION (unless it was ACT_NEW_LOCATION), dst->owner is unmodified.
// NOTE:
// *	dst is assumed to be initialized.
//	Its contents are freed before copying data from src (if the size/type are different)
//	or overwritten (if possible, given the destination size).
void ae_x_set_vector(x_vector *dst, ae_vector *src) {
   if (src->xX == dst->x_ptr) {
   // src->ptr points to the beginning of dst, attached matrices, no need to copy.
      return;
   }
   if (dst->cnt != src->cnt || dst->datatype != src->datatype) {
      if (dst->owner) ae_free(dst->x_ptr);
      dst->x_ptr = ae_malloc((size_t)(src->cnt * ae_sizeof(src->datatype)));
      if (src->cnt != 0 && dst->x_ptr == NULL) ae_break(ERR_OUT_OF_MEMORY, "ae_x_set_vector: out of memory");
      dst->last_action = ACT_NEW_LOCATION;
      dst->cnt = src->cnt;
      dst->datatype = src->datatype;
      dst->owner = true;
   } else {
      if (dst->last_action == ACT_UNCHANGED) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_SAME_LOCATION) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_NEW_LOCATION) dst->last_action = ACT_NEW_LOCATION;
      else ae_assert(false, "ae_x_set_vector: internal error");
   }
   if (src->cnt != 0) memmove(dst->x_ptr, src->xX, (size_t)(src->cnt * ae_sizeof(src->datatype)));
}

// Attach the x_vector dst to the contents of the ae_vector src.
// Ownership of memory allocated is not changed (it is still managed by the ae_vector).
// NOTES:
// *	dst is assumed to be initialized.
//	Its contents are freed before attaching to src.
// *	Assuming correctly initialized src, this function can't fail, and so doesn't need the global stack frame.
void ae_x_attach_to_vector(x_vector *dst, ae_vector *src) {
   if (dst->owner) ae_free(dst->x_ptr);
   dst->x_ptr = src->xX;
   dst->last_action = ACT_NEW_LOCATION;
   dst->cnt = src->cnt;
   dst->datatype = src->datatype;
   dst->owner = false;
}

// Clear the x_vector dst.
// Do nothing if vector is not owned by the ALGLIB++ environment.
void x_vector_free(x_vector *dst, bool/* make_automatic*/) {
   if (dst->owner) aligned_free(dst->x_ptr);
   dst->x_ptr = NULL;
   dst->cnt = 0;
}

// Copy x_matrix src into ae_matrix dst.
// dst is assumed to be uninitialized, its fields are ignored.
// The newly-created copy of src is to be owned/managed by dst.
// Both src and dst remain completely independent afterwards.
// make_automatic indicates whether or not the matrix will be registered in the ALGLIB++ environment.
void ae_matrix_init_from_x(ae_matrix *dst, x_matrix *src, bool make_automatic) {
//(@) TopFr != NULL check removed.
   ae_matrix_init(dst, (ae_int_t)src->rows, (ae_int_t)src->cols, (ae_datatype)src->datatype, make_automatic);
   if (src->cols > 0 && src->rows > 0) {
      char *p_src_row = (char *)src->x_ptr;
      char *p_dst_row = (char *)(dst->xyX[0]);
      ae_int_t row_size = ae_sizeof((ae_datatype)src->datatype) * (ae_int_t)src->cols;
      for (ae_int_t i = 0; i < src->rows; i++, p_src_row += src->stride * ae_sizeof((ae_datatype)src->datatype), p_dst_row += dst->stride * ae_sizeof((ae_datatype)src->datatype))
         memmove(p_dst_row, p_src_row, (size_t)(row_size));
   }
}

// Copy x_matrix src into ae_matrix dst by attaching dst to src.
// dst is assumed to be uninitialized, its fields are ignored.
// make_automatic indicates whether or not the matrix will be registered in the ALGLIB++ environment.
// The new matrix is attached to the source:
// *	dst shares memory with src: a write to one changes both.
// *	dst can be reallocated with ae_matrix_set_length(), but src remains untouched.
// *	src, however, can NOT be reallocated as long as dst exists.
// *	dst->is_attached is set to true to indicate that dst does not own its memory.
// Upon allocation failure, call ae_break().
void ae_matrix_init_attach_to_x(ae_matrix *dst, x_matrix *src, bool make_automatic) {
//(@) TopFr != NULL check and zero-check removed.
   ae_int_t cols = (ae_int_t)src->cols;
   ae_int_t rows = (ae_int_t)src->rows;
// Check that the X-source is densely packed.
   ae_assert(src->cols == src->stride, "ae_matrix_init_attach_to_x: unsupported stride");
// Ensure the correct size.
   ae_assert(cols == src->cols, "ae_matrix_init_attach_to_x: 32/64 overflow");
   ae_assert(rows == src->rows, "ae_matrix_init_attach_to_x: 32/64 overflow");
   ae_assert(cols >= 0 && rows >= 0, "ae_matrix_init_attach_to_x: negative length");
// If either cols or rows is 0, then they both must be made so.
   if (cols == 0) rows = 0; else if (rows == 0) cols = 0;
// Initialize.
   dst->datatype = (ae_datatype)src->datatype;
   dst->stride = cols;
   dst->is_attached = true;
// Prepare for possible errors during allocation.
   dst->rows = dst->cols = 0;
   dst->xyX = NULL;
   ae_db_init(&dst->data, (size_t)rows * sizeof(void *), make_automatic);
   dst->cols = cols;
   dst->rows = rows;
   if (dst->cols > 0 && dst->rows > 0) {
      char *p_row = (char *)src->x_ptr;
      ae_int_t rowsize = dst->stride * ae_sizeof(dst->datatype);
      void **pp_ptr = (void **)dst->data.ptr;
      dst->xyX = pp_ptr;
      for (ae_int_t i = 0; i < dst->rows; i++, p_row += rowsize) pp_ptr[i] = p_row;
   }
}

// Copy ae_matrix src to x_matrix dst.
// Not meant for use when dst is attached to src, though it may be used when src is attached to dst.
// One of the following is then applied:
// *	If src is attached to dst, no action is required or done.
// *	For independent matrices of different size: allocate storage in dst and copy src to dst.
//	dst->last_action is set to ACT_NEW_LOCATION, and dst->owner is set to true.
// *	For independent matrices of the same size: no (re)allocation is required or done.
//	Just copy src to the already-existing place.
//	dst->last_action is set to ACT_SAME_LOCATION (unless it was ACT_NEW_LOCATION), dst->owner is unmodified.
// NOTE:
// *	dst is assumed to be initialized.
//	Its contents are freed before copying data from src (if the size/type are different)
//	or overwritten (if possible, given the destination size).
void ae_x_set_matrix(x_matrix *dst, ae_matrix *src) {
   if (src->xyX != NULL && src->xyX[0] == dst->x_ptr) {
   // src->ptr points to the beginning of dst, attached matrices, no need to copy.
      return;
   }
   if (dst->cols != src->cols || dst->rows != src->rows || dst->datatype != src->datatype) {
      if (dst->owner) ae_free(dst->x_ptr);
      dst->cols = src->cols;
      dst->rows = src->rows;
      dst->stride = src->cols;
      dst->datatype = src->datatype;
      dst->x_ptr = ae_malloc((size_t)(dst->rows * (ae_int_t)dst->stride * ae_sizeof(src->datatype)));
      if (dst->rows > 0 && dst->stride != 0 && dst->x_ptr == NULL) ae_break(ERR_OUT_OF_MEMORY, "ae_x_set_matrix: out of memory");
      dst->last_action = ACT_NEW_LOCATION;
      dst->owner = true;
   } else {
      if (dst->last_action == ACT_UNCHANGED) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_SAME_LOCATION) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_NEW_LOCATION) dst->last_action = ACT_NEW_LOCATION;
      else ae_assert(false, "ae_x_set_matrix: internal error");
   }
   if (src->cols > 0 && src->rows > 0) {
      char *p_src_row = (char *)(src->xyX[0]);
      char *p_dst_row = (char *)dst->x_ptr;
      ae_int_t row_size = ae_sizeof(src->datatype) * src->cols;
      for (ae_int_t i = 0; i < src->rows; i++, p_src_row += src->stride * ae_sizeof(src->datatype), p_dst_row += dst->stride * ae_sizeof(src->datatype))
         memmove(p_dst_row, p_src_row, (size_t)(row_size));
   }
}

// Attach the x_matrix dst to the contents of the ae_matrix src.
// Ownership of memory allocated is not changed (it is still managed by the ae_matrix).
// NOTES:
// *	dst is assumed to be initialized.
//	Its contents are freed before attaching to src.
// *	Assuming correctly initialized src, this function can't fail, and so doesn't need the global stack frame.
void ae_x_attach_to_matrix(x_matrix *dst, ae_matrix *src) {
   if (dst->owner) ae_free(dst->x_ptr);
   dst->cols = src->cols;
   dst->rows = src->rows;
   dst->stride = src->stride;
   dst->datatype = src->datatype;
   dst->x_ptr = src->xyR[0];
   dst->last_action = ACT_NEW_LOCATION;
   dst->owner = false;
}

static const ae_int_t x_nb = 16; // A cut-off for recursion in the divide-and-conquer x_matrix routines.

// Symmetric/Hermitian properties: check and force
static void x_split_length(ae_int_t n, ae_int_t nb, ae_int_t *n1, ae_int_t *n2) {
   if (n <= nb) {
      *n1 = n;
      *n2 = 0;
   } else if (n % nb != 0) {
      *n2 = n % nb;
      *n1 = n - *n2;
   } else {
      *n2 = n / 2;
      *n1 = n - *n2;
      if (*n1 % nb == 0) return;
      ae_int_t r = nb - *n1 % nb;
      *n1 += r;
      *n2 -= r;
   }
}

static double x_safepythag2(double x, double y) {
   double xabs = fabs(x);
   double yabs = fabs(y);
   double w = xabs > yabs ? xabs : yabs;
   double z = xabs < yabs ? xabs : yabs;
   if (z == 0.0) return w;
   else {
      double t = z / w;
      return w * sqrt(1.0 + t * t);
   }
}

// Verify the symmetry of the len0 x len1 off-diagonal block at (offset0, offset1) with its symmetric counterpart.
// Also, accumulate updates on:
// *	*mx:	the componentwise maximum of a,
// *	*err:	the maximum componentwise difference between the lower block and the transpose of its upper counterpart.
static void is_symmetric_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1, bool *nonfinite, double *mx, double *err) {
   if (len0 <= x_nb && len1 <= x_nb) { // The base case.
      double *p1 = (double *)a->x_ptr + offset0 * a->stride + offset1;
      double *p2 = (double *)a->x_ptr + offset1 * a->stride + offset0;
      for (ae_int_t i = 0; i < len0; i++) {
         double *pcol = p2 + i;
         double *prow = p1 + i * a->stride;
         for (ae_int_t j = 0; j < len1; j++) {
            if (!isfinite(*pcol) || !isfinite(*prow)) {
               *nonfinite = true;
            } else {
               double xa = fabs(*pcol);
               *mx = *mx > xa ? *mx : xa;
               double ya = fabs(*prow);
               *mx = *mx > ya ? *mx : ya;
               double xya = fabs(*pcol - *prow);
               *err = *err > xya ? *err : xya;
            }
            pcol += a->stride;
            prow++;
         }
      }
   } else { // Reduce to smaller cases.
      if (len0 > len1) {
         ae_int_t n1, n2;
         x_split_length(len0, x_nb, &n1, &n2);
         is_symmetric_rec_off_stat(a, offset0, offset1, n1, len1, nonfinite, mx, err);
         is_symmetric_rec_off_stat(a, offset0 + n1, offset1, n2, len1, nonfinite, mx, err);
      } else {
         ae_int_t n1, n2;
         x_split_length(len1, x_nb, &n1, &n2);
         is_symmetric_rec_off_stat(a, offset0, offset1, len0, n1, nonfinite, mx, err);
         is_symmetric_rec_off_stat(a, offset0, offset1 + n1, len0, n2, nonfinite, mx, err);
      }
   }
}

// Verify that the diagonal len x len block at (offset, offset) is symmetric.
// Also, accumulate updates on:
// *	*mx:	the componentwise maximum of a,
// *	*err:	the maximum componentwise difference between the block and its transpose.
static void is_symmetric_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len, bool *nonfinite, double *mx, double *err) {
   if (len <= x_nb) { // The base case.
      double *p = (double *)a->x_ptr + offset * a->stride + offset;
      for (ae_int_t i = 0; i < len; i++) {
         double *pcol = p + i;
         double *prow = p + i * a->stride;
         for (ae_int_t j = 0; j < i; j++, pcol += a->stride, prow++) {
            if (!isfinite(*pcol) || !isfinite(*prow)) {
               *nonfinite = true;
            } else {
               double xa = fabs(*pcol);
               *mx = *mx > xa ? *mx : xa;
               double ya = fabs(*prow);
               *mx = *mx > ya ? *mx : ya;
               double xya = fabs(*pcol - *prow);
               *err = *err > xya ? *err : xya;
            }
         }
         double norm = fabs(p[i + i * a->stride]);
         *mx = *mx > norm ? *mx : norm;
      }
   } else { // Reduce to smaller cases.
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      is_symmetric_rec_diag_stat(a, offset, n1, nonfinite, mx, err);
      is_symmetric_rec_diag_stat(a, offset + n1, n2, nonfinite, mx, err);
      is_symmetric_rec_off_stat(a, offset + n1, offset, n2, n1, nonfinite, mx, err);
   }
}

static bool x_is_symmetric(x_matrix *a) {
   if (a->datatype != DT_REAL || a->cols != a->rows) return false;
   else if (a->cols == 0 || a->rows == 0) return true;
   ae_state_init();
   double mx = 0.0;
   double err = 0.0;
   bool nonfinite = false;
   is_symmetric_rec_diag_stat(a, 0, (ae_int_t)a->rows, &nonfinite, &mx, &err);
   return !nonfinite && (mx == 0.0 || err / mx <= 1.0E-14);
}

// Verify the Hermitian conjugacy of the len0 x len1 off-diagonal block at (offset0, offset1) with its symmetric counterpart.
// Also, accumulate updates on:
// *	*mx:	the componentwise maximum of a
// *	*err:	the maximum componentwise difference between the lower block and the Hermitian conjugate of its upper counterpart.
static void is_hermitian_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1, bool *nonfinite, double *mx, double *err) {
   if (len0 <= x_nb && len1 <= x_nb) { // The base case.
      complex *p1 = (complex *)a->x_ptr + offset0 * a->stride + offset1;
      complex *p2 = (complex *)a->x_ptr + offset1 * a->stride + offset0;
      for (ae_int_t i = 0; i < len0; i++) {
         complex *pcol = p2 + i;
         complex *prow = p1 + i * a->stride;
         for (ae_int_t j = 0; j < len1; j++) {
            if (!isfinite(pcol->x) || !isfinite(pcol->y) || !isfinite(prow->x) || !isfinite(prow->y)) {
               *nonfinite = true;
            } else {
               double xa = x_safepythag2(pcol->x, pcol->y);
               *mx = *mx > xa ? *mx : xa;
               double ya = x_safepythag2(prow->x, prow->y);
               *mx = *mx > ya ? *mx : ya;
               double xya = x_safepythag2(pcol->x - prow->x, pcol->y + prow->y);
               *err = *err > xya ? *err : xya;
            }
            pcol += a->stride;
            prow++;
         }
      }
   } else { // Reduce to smaller cases.
      if (len0 > len1) {
         ae_int_t n1, n2;
         x_split_length(len0, x_nb, &n1, &n2);
         is_hermitian_rec_off_stat(a, offset0, offset1, n1, len1, nonfinite, mx, err);
         is_hermitian_rec_off_stat(a, offset0 + n1, offset1, n2, len1, nonfinite, mx, err);
      } else {
         ae_int_t n1, n2;
         x_split_length(len1, x_nb, &n1, &n2);
         is_hermitian_rec_off_stat(a, offset0, offset1, len0, n1, nonfinite, mx, err);
         is_hermitian_rec_off_stat(a, offset0, offset1 + n1, len0, n2, nonfinite, mx, err);
      }
   }
}

// Verify that the len x len diagonal block at (offset, offset) is Hermitian.
// Also, accumulate updates on:
// *	*mx:	the componentwise maximum of a
// *	*err:	the maximum componentwise difference between the block and its Hermitian conjugate.
static void is_hermitian_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len, bool *nonfinite, double *mx, double *err) {
   if (len <= x_nb) { // The base case.
      complex *p = (complex *)a->x_ptr + offset * a->stride + offset;
      for (ae_int_t i = 0; i < len; i++) {
         complex *pcol = p + i;
         complex *prow = p + i * a->stride;
         for (ae_int_t j = 0; j < i; j++, pcol += a->stride, prow++) {
            if (!isfinite(pcol->x) || !isfinite(pcol->y) || !isfinite(prow->x) || !isfinite(prow->y)) {
               *nonfinite = true;
            } else {
               double xa = x_safepythag2(pcol->x, pcol->y);
               *mx = *mx > xa ? *mx : xa;
               double ya = x_safepythag2(prow->x, prow->y);
               *mx = *mx > ya ? *mx : ya;
               double xya = x_safepythag2(pcol->x - prow->x, pcol->y + prow->y);
               *err = *err > xya ? *err : xya;
            }
         }
         if (!isfinite(p[i + i * a->stride].x) || !isfinite(p[i + i * a->stride].y)) {
            *nonfinite = true;
         } else {
            double xa = fabs(p[i + i * a->stride].x);
            *mx = *mx > xa ? *mx : xa;
            double ya = fabs(p[i + i * a->stride].y);
            *err = *err > ya ? *err : ya;
         }
      }
   } else { // Reduce to smaller cases.
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      is_hermitian_rec_diag_stat(a, offset, n1, nonfinite, mx, err);
      is_hermitian_rec_diag_stat(a, offset + n1, n2, nonfinite, mx, err);
      is_hermitian_rec_off_stat(a, offset + n1, offset, n2, n1, nonfinite, mx, err);
   }
}

static bool x_is_hermitian(x_matrix *a) {
   if (a->datatype != DT_COMPLEX || a->cols != a->rows) return false;
   else if (a->cols == 0 || a->rows == 0) return true;
   ae_state_init();
   double mx = 0.0;
   double err = 0.0;
   bool nonfinite = false;
   is_hermitian_rec_diag_stat(a, 0, (ae_int_t)a->rows, &nonfinite, &mx, &err);
   return !nonfinite && (mx == 0.0 || err / mx <= 1.0E-14);
}

// Copy the transpose of the len0 x len1 off-diagonal block at (offset0, offset1) to its symmetric counterpart.
static void force_symmetric_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1) {
   if (len0 <= x_nb && len1 <= x_nb) { // The base case.
      double *p1 = (double *)a->x_ptr + offset0 * a->stride + offset1;
      double *p2 = (double *)a->x_ptr + offset1 * a->stride + offset0;
      for (ae_int_t i = 0; i < len0; i++) {
         double *pcol = p2 + i;
         double *prow = p1 + i * a->stride;
         for (ae_int_t j = 0; j < len1; j++) {
            *pcol = *prow;
            pcol += a->stride;
            prow++;
         }
      }
   } else { // Reduce to smaller cases.
      if (len0 > len1) {
         ae_int_t n1, n2;
         x_split_length(len0, x_nb, &n1, &n2);
         force_symmetric_rec_off_stat(a, offset0, offset1, n1, len1);
         force_symmetric_rec_off_stat(a, offset0 + n1, offset1, n2, len1);
      } else {
         ae_int_t n1, n2;
         x_split_length(len1, x_nb, &n1, &n2);
         force_symmetric_rec_off_stat(a, offset0, offset1, len0, n1);
         force_symmetric_rec_off_stat(a, offset0, offset1 + n1, len0, n2);
      }
   }
}

// Copy the transpose of the lower part of the len x len diagonal block at (offset, offset) to its upper part.
static void force_symmetric_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len) {
   if (len <= x_nb) { // The base case.
      double *p = (double *)a->x_ptr + offset * a->stride + offset;
      for (ae_int_t i = 0; i < len; i++) {
         double *pcol = p + i;
         double *prow = p + i * a->stride;
         for (ae_int_t j = 0; j < i; j++, pcol += a->stride, prow++)
            *pcol = *prow;
      }
   } else { // Reduce to smaller cases.
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      force_symmetric_rec_diag_stat(a, offset, n1);
      force_symmetric_rec_diag_stat(a, offset + n1, n2);
      force_symmetric_rec_off_stat(a, offset + n1, offset, n2, n1);
   }
}

static bool x_force_symmetric(x_matrix *a) {
   if (a->datatype != DT_REAL || a->cols != a->rows) return false;
   if (a->cols > 0 && a->rows > 0) force_symmetric_rec_diag_stat(a, 0, (ae_int_t)a->rows);
   return true;
}

// Copy the Hermitian conjugate of the len0 x len1 off-diagonal block at (offset0, offset1) to its symmetric counterpart.
static void force_hermitian_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1) {
   if (len0 <= x_nb && len1 <= x_nb) { // The base case.
      complex *p1 = (complex *)a->x_ptr + offset0 * a->stride + offset1;
      complex *p2 = (complex *)a->x_ptr + offset1 * a->stride + offset0;
      for (ae_int_t i = 0; i < len0; i++) {
         complex *pcol = p2 + i;
         complex *prow = p1 + i * a->stride;
         for (ae_int_t j = 0; j < len1; j++) {
            *pcol = *prow;
            pcol += a->stride;
            prow++;
         }
      }
   } else { // Reduce to smaller cases.
      if (len0 > len1) {
         ae_int_t n1, n2;
         x_split_length(len0, x_nb, &n1, &n2);
         force_hermitian_rec_off_stat(a, offset0, offset1, n1, len1);
         force_hermitian_rec_off_stat(a, offset0 + n1, offset1, n2, len1);
      } else {
         ae_int_t n1, n2;
         x_split_length(len1, x_nb, &n1, &n2);
         force_hermitian_rec_off_stat(a, offset0, offset1, len0, n1);
         force_hermitian_rec_off_stat(a, offset0, offset1 + n1, len0, n2);
      }
   }
}

// Copy the Hermitian conjugate of the lower part of the len x len diagonal block at (offset, offset) to the upper part.
static void force_hermitian_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len) {
   if (len <= x_nb) { // The base case.
      complex *p = (complex *)a->x_ptr + offset * a->stride + offset;
      for (ae_int_t i = 0; i < len; i++) {
         complex *pcol = p + i;
         complex *prow = p + i * a->stride;
         for (ae_int_t j = 0; j < i; j++, pcol += a->stride, prow++)
            *pcol = *prow;
      }
   } else { // Reduce to smaller cases.
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      force_hermitian_rec_diag_stat(a, offset, n1);
      force_hermitian_rec_diag_stat(a, offset + n1, n2);
      force_hermitian_rec_off_stat(a, offset + n1, offset, n2, n1);
   }
}

static bool x_force_hermitian(x_matrix *a) {
   if (a->datatype != DT_COMPLEX || a->cols != a->rows) return false;
   if (a->cols > 0 && a->rows > 0) force_hermitian_rec_diag_stat(a, 0, (ae_int_t)a->rows);
   return true;
}

bool ae_is_symmetric(ae_matrix *a) {
   x_matrix x;
   x.owner = false;
   ae_x_attach_to_matrix(&x, a);
   return x_is_symmetric(&x);
}

bool ae_is_hermitian(ae_matrix *a) {
   x_matrix x;
   x.owner = false;
   ae_x_attach_to_matrix(&x, a);
   return x_is_hermitian(&x);
}

bool ae_force_symmetric(ae_matrix *a) {
   x_matrix x;
   x.owner = false;
   ae_x_attach_to_matrix(&x, a);
   return x_force_symmetric(&x);
}

bool ae_force_hermitian(ae_matrix *a) {
   x_matrix x;
   x.owner = false;
   ae_x_attach_to_matrix(&x, a);
   return x_force_hermitian(&x);
}

//(@) ae_yield(), originally had external linkage, but isn't documented in the API or used elsewhere.
// Make the calling thread relinquish the CPU.
// The thread is moved to the end of the queue and some other thread gets to run.
// NOTE:
// *	this function should NOT be called when AE_OS is AE_OTHER_OS - the whole program will be abnormally terminated.
static void ae_yield() {
#if AE_OS == AE_POSIX
   sched_yield();
#elif AE_OS == AE_WINDOWS
   if (!SwitchToThread()) Sleep(0);
#else
   abort();
#endif
}

// Perform a given number of spin-wait iterations.
static void ae_spin_wait(ae_int_t cnt) {
// This variable is used to prevent some tricky optimizations which may degrade multi-threaded performance.
// It was used by the SMP version to prevent optimizations that would remove objects that can nowadays be declared "volatile".
   static volatile ae_int_t ae_never_change_it = 1;
// These strange operations with ae_never_change_it are necessary to prevent compiler optimization of the loop.
// Very unlikely because no one will wait for such amount of cycles.
   if (cnt > 0x12345678)
      ae_never_change_it = cnt % 10;
// Spin wait, test a condition which will never be true.
   for (volatile ae_int_t i = 0; i < cnt; i++)
      if (ae_never_change_it > 0)
         ae_never_change_it--;
}

// Internal OS-dependent lock routines.
static const int AE_LOCK_ALIGNMENT = 0x10, AE_LOCK_CYCLES = 0x200, AE_LOCK_TESTS_BEFORE_YIELD = 0x10;

// The internal OS-dependent lock structure.
struct _lock;

// Initialize the _lock p to be internally used by the high-level ae_lock structure.
// It is statically allocated, no malloc() calls are performed during its allocation,
// but _ae_free_lock() must be used to deallocate it properly.
static inline void _ae_init_lock(_lock *p);

// Acquire the _lock p; the low-level workhorse function used by ae_acquire_lock().
static inline void _ae_acquire_lock(_lock *p);

// Release the _lock p; the low-level workhorse function used by ae_release_lock().
static inline void _ae_release_lock(_lock *p);

// Free the _lock p.
static inline void _ae_free_lock(_lock *p);

#if AE_OS == AE_POSIX
struct _lock { pthread_mutex_t mutex; };
static inline void _ae_init_lock(_lock *p) { pthread_mutex_init(&p->mutex, NULL); }
static inline void _ae_acquire_lock(_lock *p) {
   ae_int_t cnt = 0;
   while (pthread_mutex_trylock(&p->mutex) != 0) {
      ae_spin_wait(AE_LOCK_CYCLES);
      if (++cnt % AE_LOCK_TESTS_BEFORE_YIELD == 0)
         ae_yield();
   }
}
static inline void _ae_release_lock(_lock *p) { pthread_mutex_unlock(&p->mutex); }
static inline void _ae_free_lock(_lock *p) { pthread_mutex_destroy(&p->mutex); }
#elif AE_OS == AE_WINDOWS
struct _lock {
   volatile ae_int_t *volatile p_lock;
   char buf[sizeof(ae_int_t) + AE_LOCK_ALIGNMENT];
};
static inline void _ae_init_lock(_lock *p) {
   p->p_lock = (ae_int_t *)ae_align((void *)&p->buf, AE_LOCK_ALIGNMENT);
   p->p_lock[0] = 0;
}
static inline void _ae_acquire_lock(_lock *p) {
#   ifdef AE_SMP_DEBUGCOUNTERS
   __declspec(align(AE_LOCK_ALIGNMENT)) volatile ae_int64_t _ae_dbg_lock_acquisitions = 0;
   __declspec(align(AE_LOCK_ALIGNMENT)) volatile ae_int64_t _ae_dbg_lock_spinwaits = 0;
   __declspec(align(AE_LOCK_ALIGNMENT)) volatile ae_int64_t _ae_dbg_lock_yields = 0;
#      define BumpLock(X) (InterlockedIncrement((LONG volatile *)(X)))
#   else
#      define BumpLock(X) (X)
#   endif
   ae_int_t cnt = 0;
   BumpLock(&_ae_dbg_lock_acquisitions);
   while (InterlockedCompareExchange((LONG volatile *)p->p_lock, 1, 0) != 0) {
      ae_spin_wait(AE_LOCK_CYCLES);
      BumpLock(&_ae_dbg_lock_spinwaits);
      if (++cnt % AE_LOCK_TESTS_BEFORE_YIELD == 0) {
         BumpLock(&_ae_dbg_lock_yields);
         ae_yield();
      }
   }
#   undef BumpLock
}
static inline void _ae_release_lock(_lock *p) { InterlockedExchange((LONG volatile *)p->p_lock, 0); }
static inline void _ae_free_lock(_lock *p) { }
#else
struct _lock { bool is_locked; };
static inline void _ae_init_lock(_lock *p) { p->is_locked = false; }
static inline void _ae_acquire_lock(_lock *p) {
   AE_CRITICAL_ASSERT(!p->is_locked);
   p->is_locked = true;
}
static inline void _ae_release_lock(_lock *p) { p->is_locked = false; }
static inline void _ae_free_lock(_lock *p) { }
#endif

// Initialize a free ae_lock.
// NOTES:
// *	make_automatic indicates whether or not the lock is to be added to the automatic memory management list.
// *	As a special exception, this function allows you to specify TopFr == NULL.
//	In this case, exceptions arising during construction are handled as critical failures, resulting in abort().
//	make_automatic must be false on such calls.
// *	is_static indicates that the lock is to be made "eternal" (i.e. static),
//	which is expected to persist until the end of the execution of the program.
// *	Eternal locks can not be deallocated (cleared) and do not increase debug allocation counters.
//	Errors during allocation of such locks are considered critical exceptions and are handled by calling abort().
void ae_init_lock(ae_lock *lock, bool is_static, bool make_automatic) {
//(@) Zero-check removed.
   bool need_auto = !is_static && TopFr == NULL;
   if (need_auto) {
      AE_CRITICAL_ASSERT(!make_automatic);
      ae_state_init();
   }
   lock->is_static = is_static;
   size_t size = sizeof(_lock);
   if (!is_static) ae_db_init(&lock->db, size, make_automatic);
   lock->lock_ptr = !is_static ? lock->db.ptr : size == 0 || _force_malloc_failure ? NULL : malloc(size);
   _ae_init_lock((_lock *)lock->lock_ptr);
   if (need_auto) ae_state_clear();
}

// Acquire an ae_lock.
// If the lock is busy and ae_yield() is supported, we ae_yield() after retrying several times, with tight spin waits in between.
void ae_acquire_lock(ae_lock *lock) { _ae_acquire_lock((_lock *)lock->lock_ptr); }

// Release an ae_lock.
void ae_release_lock(ae_lock *lock) { _ae_release_lock((_lock *)lock->lock_ptr); }

// Free an ae_lock.
void ae_free_lock(ae_lock *lock) {
   AE_CRITICAL_ASSERT(!lock->is_static);
   _lock *p = (_lock *)lock->lock_ptr;
   if (p != NULL) _ae_free_lock(p);
   ae_db_free(&lock->db);
}

// A reduced form ae_shared_pool_free() suitable for use as the deallocation function on a frame.
static void ae_shared_pool_destroy(void *_dst) { ae_shared_pool_free(_dst, false); }

// A new ae_shared_pool structure for _dst.
// _dst is assumed to be preallocated, and uninitialized, its fields are ignored.
// make_automatic indicates whether or not _dst is to be removed after leaving the current frame
// as opposed to being the field of some other object.
// Upon allocation failure, call ae_break().
void ae_shared_pool_init(void *_dst, bool make_automatic) {
//(@) TopFr != NULL check and zero-check removed.
   ae_shared_pool *dst = (ae_shared_pool *)_dst;
// Attach to the frame first, just to be sure that if we fail within malloc(), we will be able to clean up the object.
   dst->frame_entry.deallocator = ae_shared_pool_destroy;
   dst->frame_entry.ptr = dst;
   if (make_automatic) ae_db_attach(&dst->frame_entry);
// Initialize.
   dst->seed_object = NULL;
   dst->recycled_objects = NULL;
   dst->recycled_entries = NULL;
   dst->enumeration_counter = NULL;
   dst->size_of_object = 0;
   dst->init = NULL;
   dst->copy = NULL;
   dst->free = NULL;
   ae_init_lock(&dst->pool_lock, false, false);
}

// Clear all dynamically allocated fields of the ae_shared_pool dst except for the lock.
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
static void ae_shared_pool_internalclear(ae_shared_pool *dst, bool make_automatic) {
// Free the seed.
   if (dst->seed_object != NULL) {
      dst->free((void *)dst->seed_object, make_automatic);
      ae_free((void *)dst->seed_object);
      dst->seed_object = NULL;
   }
// Free the recycled objects.
   for (ae_shared_pool_entry *ptr = dst->recycled_objects; ptr != NULL; ) {
      ae_shared_pool_entry *tmp = (ae_shared_pool_entry *)ptr->next_entry;
      dst->free(ptr->obj, make_automatic);
      ae_free(ptr->obj);
      ae_free(ptr);
      ptr = tmp;
   }
   dst->recycled_objects = NULL;
// Free the recycled entries.
   for (ae_shared_pool_entry *ptr = dst->recycled_entries; ptr != NULL; ) {
      ae_shared_pool_entry *tmp = (ae_shared_pool_entry *)ptr->next_entry;
      ae_free(ptr);
      ptr = tmp;
   }
   dst->recycled_entries = NULL;
}

// Copy the ae_shared_pool _src into the previously-allocated ae_shared_pool _dst.
// _dst is assumed to be uninitialized and its fields are ignored.
// make_automatic indicates whether or not _dst is to be removed after leaving the current frame.
// as opposed to being the field of some other object.
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
void ae_shared_pool_copy(void *_dst, const void *_src, bool make_automatic) {
//(@) TopFr != NULL check removed (for TopFr != NULL: allocation errors result in an exception from ae_malloc()).
   ae_shared_pool *dst = (ae_shared_pool *)_dst;
   const ae_shared_pool *src = (const ae_shared_pool *)_src;
   ae_shared_pool_init(dst, make_automatic);
// Copy the non-pointer fields.
   dst->size_of_object = src->size_of_object;
   dst->init = src->init;
   dst->copy = src->copy;
   dst->free = src->free;
// Copy the seed object.
   if (src->seed_object != NULL) {
      dst->seed_object = ae_malloc(dst->size_of_object);
      memset(dst->seed_object, 0, dst->size_of_object);
      dst->copy(dst->seed_object, src->seed_object, false);
   }
// Copy the recycled objects.
   dst->recycled_objects = NULL;
   for (ae_shared_pool_entry *ptr = src->recycled_objects; ptr != NULL; ptr = (ae_shared_pool_entry *)ptr->next_entry) {
   // Allocate an entry, immediately add it to the recycled list (we do not want to lose it in case of future malloc failures).
      ae_shared_pool_entry *tmp = (ae_shared_pool_entry *)ae_malloc(sizeof *tmp);
      memset(tmp, 0, sizeof *tmp);
      tmp->next_entry = dst->recycled_objects;
      dst->recycled_objects = tmp;
   // Prepare a place for the object and copy() it.
      tmp->obj = ae_malloc(dst->size_of_object);
      memset(tmp->obj, 0, dst->size_of_object);
      dst->copy(tmp->obj, ptr->obj, false);
   }
// Recycled entries are not copied because they do not store any information.
   dst->recycled_entries = NULL;
// The enumeration counter is reset on copying.
   dst->enumeration_counter = NULL;
// Initialize the frame record.
   dst->frame_entry.deallocator = ae_shared_pool_destroy;
   dst->frame_entry.ptr = dst;
}

// Free the ae_shared_pool _dst.
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
void ae_shared_pool_free(void *_dst, bool make_automatic) {
   ae_shared_pool *dst = (ae_shared_pool *)_dst;
// Clear the seed and lists.
   ae_shared_pool_internalclear(dst, make_automatic);
// Clear the fields.
   dst->seed_object = NULL;
   dst->recycled_objects = NULL;
   dst->recycled_entries = NULL;
   dst->enumeration_counter = NULL;
   dst->size_of_object = 0;
   dst->init = NULL;
   dst->copy = NULL;
   dst->free = NULL;
   if (!make_automatic) ae_free_lock(&dst->pool_lock);
}

// True if and only if the already-initialized ae_shared_pool dst is seeded; i.e. if its internal seed object is set.
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
bool ae_shared_pool_is_initialized(ae_shared_pool *dst) {
   return dst->seed_object != NULL;
}

// Seed the already-initialized ae_shared_pool dst by setting its internal seed object to seed_object::size_of_object,
// freeing everything owned by dst (the current seed object and recycled objects).
// The make, copy and free functions are set respectively to init(), copy() and free().
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
void ae_shared_pool_set_seed(ae_shared_pool *dst, void *seed_object, size_t size_of_object, ae_init_op init, ae_copy_op copy, ae_free_op free) {
//(@) TopFr != NULL check removed (for TopFr != NULL: allocation errors result in an exception from ae_malloc()).
// Free the internal objects.
   ae_shared_pool_internalclear(dst, false);
// Set the non-pointer fields.
   dst->size_of_object = size_of_object;
   dst->init = init;
   dst->copy = copy;
   dst->free = free;
// Set the seed object and its size.
   dst->seed_object = ae_malloc(size_of_object);
   memset(dst->seed_object, 0, size_of_object);
   copy(dst->seed_object, seed_object, false);
}

// Copy the seed object from the ae_shared_pool pool into the target ae_smart_ptr ptr, making pptr its owner.
// Any non-NULL pointer owned by pptr is deallocated before storing the value retrieved from pool.
// NOTE:
// *	This function IS thread-safe.
//	It acquires a pool lock during its operation and can be used simultaneously from several threads.
void ae_shared_pool_retrieve(ae_shared_pool *pool, ae_smart_ptr *pptr) {
//(@) TopFr != NULL check removed (for TopFr != NULL: allocation errors result in an exception from ae_malloc()).
// Require pool to be seeded.
   ae_assert(pool->seed_object != NULL, "ae_shared_pool_retrieve: shared pool is not seeded");
// Acquire a lock.
   ae_acquire_lock(&pool->pool_lock);
   if (pool->recycled_objects != NULL) { // Try to reuse recycled objects.
   // Retrieve an entry/object from the list of recycled objects.
      ae_shared_pool_entry *result = pool->recycled_objects;
      pool->recycled_objects = (ae_shared_pool_entry *)pool->recycled_objects->next_entry;
      void *new_obj = result->obj;
      result->obj = NULL;
   // Recycle the entry.
      result->next_entry = pool->recycled_entries;
      pool->recycled_entries = result;
   // Release the lock.
      ae_release_lock(&pool->pool_lock);
   // Assign the object to the smart pointer.
      ae_smart_ptr_assign(pptr, new_obj, true, true, pool->size_of_object, pool->copy, pool->free);
   } else {
   // Release the lock; we do not need it anymore because the copy constructor does not modify the source variable.
      ae_release_lock(&pool->pool_lock);
   // Create a new object from the seed.
      void *new_obj = ae_malloc(pool->size_of_object);
      memset(new_obj, 0, pool->size_of_object);
   // Immediately assign the object to the smart pointer (so as not to lose it, in case of future failures).
      ae_smart_ptr_assign(pptr, new_obj, true, true, pool->size_of_object, pool->copy, pool->free);
   // Do the actual copying.
   //(@) Before this line, the smart pointer points to a zero-filled instance. (No longer applicable.)
      pool->copy(new_obj, pool->seed_object, false);
   }
}

// Recycle the object owned by the smart pointer pptr by moving it into the ae_shared_pool pool.
// pptr must own the object.
// Afterwards, it owns the NULL pointer.
// NOTE:
// *	This function IS thread-safe.
//	It acquires a pool lock during its operation and can be used simultaneously from several threads.
void ae_shared_pool_recycle(ae_shared_pool *pool, ae_smart_ptr *pptr) {
//(@) TopFr != NULL check removed (for TopFr != NULL: allocation errors result in an exception from ae_malloc()).
// Require pool to be seeded and pptr != NULL and to own the object.
   ae_assert(pool->seed_object != NULL, "ae_shared_pool_recycle: shared pool is not seeded");
   ae_assert(pptr->is_owner, "ae_shared_pool_recycle: pptr does not own its pointer");
   ae_assert(pptr->ptr != NULL, "ae_shared_pool_recycle: pptr == NULL");
// Acquire the lock and the shared pool entry (reusing an entry from recycled_entries, if there are any).
   ae_acquire_lock(&pool->pool_lock);
   ae_shared_pool_entry *new_entry;
   if (pool->recycled_entries != NULL) {
   // Reuse a previously allocated entry.
      new_entry = pool->recycled_entries;
      pool->recycled_entries = (ae_shared_pool_entry *)new_entry->next_entry;
   } else {
   // Allocate a new entry.
   // NOTE:
   // *	Unlock the pool first
   //	so as to prevent the pool from being left in a locked state in case ae_malloc() raises an exception.
      ae_release_lock(&pool->pool_lock);
      new_entry = (ae_shared_pool_entry *)ae_malloc(sizeof *new_entry);
      ae_acquire_lock(&pool->pool_lock);
   }
// Recycle the object, the lock object and the source pointer.
   new_entry->obj = pptr->ptr;
   new_entry->next_entry = pool->recycled_objects;
   pool->recycled_objects = new_entry;
   ae_release_lock(&pool->pool_lock);
   ae_smart_ptr_release(pptr);
}

// Clear the internal list of recycled objects, keeping intact the seed object managed by the ae_shared_pool pool.
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
void ae_shared_pool_clear_recycled(ae_shared_pool *pool, bool make_automatic) {
// Acquire the pool lock, extract the list of recycled objects and immediately release the lock.
// In the unlikely event we crash during memory deallocation, it is better to have the pool lock released first.
   ae_acquire_lock(&pool->pool_lock);
   ae_shared_pool_entry *ptr = pool->recycled_objects;
   pool->recycled_objects = NULL;
   ae_release_lock(&pool->pool_lock);
// Clear the recycled objects.
   while (ptr != NULL) {
      ae_shared_pool_entry *tmp = (ae_shared_pool_entry *)ptr->next_entry;
      pool->free(ptr->obj, make_automatic);
      ae_free(ptr->obj);
      ae_free(ptr);
      ptr = tmp;
   }
}

// Allow the recycled elements of the ae_shared_pool pool to be enumerated.
// The pointer to the first recycled object is stored in the smart pointer pptr.
// IMPORTANT:
// *	The recycled object is KEPT in pool and ownership is NOT passed to pptr.
// *	If there are no recycled objects left in pool or pool is not seeded, then NULL is stored into pptr.
// *	Any non-NULL pointer owned by pptr is deallocated before storing the value retrieved from pool.
// *	This function IS NOT thread-safe
// *	You should NOT modify pool during enumeration (although you can modify the state of the objects retrieved from pool).
void ae_shared_pool_first_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr) {
// Modify the internal enumeration counter.
   pool->enumeration_counter = pool->recycled_objects;
// Exit on an empty list.
   if (pool->enumeration_counter == NULL)
      ae_smart_ptr_assign(pptr, NULL, false, false, 0, NULL, NULL);
// Assign the object to a smart pointer.
   else
      ae_smart_ptr_assign(pptr, pool->enumeration_counter->obj, false, false, 0, NULL, NULL);
}

// Allow the recycled elements of the ae_shared_pool pool to be enumerated.
// The pointer to the next recycled object is stored in the smart pointer pptr.
// In all other respects, this is the same as ae_shared_pool_first_recycled().
void ae_shared_pool_next_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr) {
// Exit on the end of list.
   if (pool->enumeration_counter == NULL) {
      ae_smart_ptr_assign(pptr, NULL, false, false, 0, NULL, NULL);
      return;
   }
// Modify the internal enumeration counter.
   pool->enumeration_counter = (ae_shared_pool_entry *)pool->enumeration_counter->next_entry;
// Exit on an empty list.
   if (pool->enumeration_counter == NULL)
      ae_smart_ptr_assign(pptr, NULL, false, false, 0, NULL, NULL);
// Assign the object to a smart pointer.
   else
      ae_smart_ptr_assign(pptr, pool->enumeration_counter->obj, false, false, 0, NULL, NULL);
}

// Clear the internal list of recycled objects and the seed object from the ae_shared_pool pool,
// while leaving pool intact for reseeding and subsequent reuse.
// NOTE:
// *	This function is NOT thread-safe.
//	It does NOT try to acquire a pool lock, and should NOT be used simultaneously from other threads.
void ae_shared_pool_reset(ae_shared_pool *pool) {
// Clear the seed and lists.
   ae_shared_pool_internalclear(pool, false);
// Clear the fields.
   pool->seed_object = NULL;
   pool->recycled_objects = NULL;
   pool->recycled_entries = NULL;
   pool->enumeration_counter = NULL;
   pool->size_of_object = 0;
   pool->init = NULL;
   pool->copy = NULL;
   pool->free = NULL;
}

void ae_obj_array_free(ae_obj_array *dst, bool make_automatic);

// A reduced form ae_obj_array_free() suitable for use as the deallocation function on a frame.
static void ae_obj_array_destroy(void *dst) {
   ae_obj_array_free((ae_obj_array *)dst, false);
}

// Make a new array of objects into a pre-allocated (zeroed-out) destination, dst.
// Error handling:
// *	TopFr != NULL is required, else abort().
// *	Upon failure with TopFr != NULL, call ae_break().
// After initialization, the smart pointer stores the NULL pointer.
void ae_obj_array_init(ae_obj_array *dst, bool make_automatic) {
   AE_CRITICAL_ASSERT(TopFr != NULL);
// AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof *dst));
// First, attach to the frame, just to be sure that we will clean everything if we generate an exception during initialization.
   dst->frame_entry.deallocator = ae_obj_array_destroy;
   dst->frame_entry.ptr = dst;
   if (make_automatic) ae_db_attach(&dst->frame_entry);
// Initialize.
   dst->cnt = 0;
   dst->capacity = 0;
   dst->fixed_capacity = false;
   dst->pp_obj_ptr = NULL;
   dst->pp_obj_sizes = NULL;
   dst->pp_copy = NULL;
   dst->pp_free = NULL;
   ae_init_lock(&dst->array_lock, false, false);
   ae_init_lock(&dst->barrier_lock, false, false);
}

// Copy a new array of objects into a pre-allocated (zeroed-out) destination,
// with independent copies of all objects owned by the array being created.
// NOTES:
// *	TopFr != NULL is required, else abort().
// *	This function is NOT thread-safe.
//	It does not acquire array lock, so you should NOT call it when array can be used by another thread.
void ae_obj_array_copy(ae_obj_array *dst, const ae_obj_array *src, bool make_automatic) {
   AE_CRITICAL_ASSERT(TopFr != NULL);
// AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof *dst));
// Initialize the array; we know that an empty array has NULL internal pointers.
   ae_obj_array_init(dst, make_automatic);
   AE_CRITICAL_ASSERT(dst->capacity == 0);
   AE_CRITICAL_ASSERT(dst->pp_obj_ptr == NULL);
   AE_CRITICAL_ASSERT(dst->pp_obj_sizes == NULL);
   AE_CRITICAL_ASSERT(dst->pp_copy == NULL);
   AE_CRITICAL_ASSERT(dst->pp_free == NULL);
// Copy the fields.
   dst->cnt = src->cnt;
   dst->capacity = src->capacity;
   dst->fixed_capacity = src->fixed_capacity;
   AE_CRITICAL_ASSERT(src->cnt <= src->capacity);
// Copy the data.
   if (dst->capacity > 0) {
      dst->pp_obj_ptr = (void **)ae_malloc_zero((size_t)dst->capacity * sizeof *dst->pp_obj_ptr);
      dst->pp_obj_sizes = (ae_int_t *) ae_malloc_zero((size_t)dst->capacity * sizeof *dst->pp_obj_sizes);
      dst->pp_copy = (ae_copy_op *) ae_malloc_zero((size_t)dst->capacity * sizeof *dst->pp_copy);
      dst->pp_free = (ae_free_op *) ae_malloc_zero((size_t)dst->capacity * sizeof *dst->pp_free);
      for (ae_int_t i = 0; i < dst->cnt; i++) {
         dst->pp_free[i] = src->pp_free[i];
         dst->pp_copy[i] = src->pp_copy[i];
         dst->pp_obj_sizes[i] = src->pp_obj_sizes[i];
         dst->pp_obj_ptr[i] = ae_malloc_zero((size_t)dst->pp_obj_sizes[i]);
         dst->pp_copy[i](dst->pp_obj_ptr[i], src->pp_obj_ptr[i], false);
      }
   }
}

// Free a dynamic objects array.
// Afterwards, all objects managed by the array are freed.
// The array capacity does not change.
// NOTE:
// *	this function is thread-unsafe.
void ae_obj_array_free(ae_obj_array *dst, bool make_automatic) {
   for (ae_int_t i = 0; i < dst->cnt; i++) if (dst->pp_obj_ptr[i] != NULL) {
      dst->pp_free[i](dst->pp_obj_ptr[i], make_automatic);
      ae_free(dst->pp_obj_ptr[i]);
      dst->pp_obj_ptr[i] = NULL;
      dst->pp_obj_sizes[i] = 0;
      dst->pp_copy[i] = NULL;
      dst->pp_free[i] = NULL;
   }
   dst->cnt = 0;
   if (!make_automatic) {
      if (dst->pp_obj_ptr != NULL) ae_free(dst->pp_obj_ptr);
      if (dst->pp_obj_sizes != NULL) ae_free(dst->pp_obj_sizes);
      if (dst->pp_copy != NULL) ae_free(dst->pp_copy);
      if (dst->pp_free != NULL) ae_free(dst->pp_free);
      ae_free_lock(&dst->array_lock);
      ae_free_lock(&dst->barrier_lock);
   }
}

// The array length.
// This function is thread-safe, i.e. it can be combined with functions adding elements to the array.
// If it is called when another thread adds an element to the array, the function will either:
// *	return the old array size
// *	return the new array size, but ONLY after the new element is added to the array.
ae_int_t ae_obj_array_get_length(ae_obj_array *dst) { return dst->cnt; }

// Modify the array capacity, ignoring the fixed_capacity flag.
// For internal use only.
// Return false on memory reallocation failure, true otherwise.
static bool _ae_obj_array_set_capacity(ae_obj_array *arr, ae_int_t new_capacity) {
// Integrity checks.
   ae_assert(arr->cnt <= new_capacity, "_ae_obj_array_set_capacity: new capacity is less than present size");
// Quick exit.
   if (arr->cnt == new_capacity) return true;
// Increase the capacity.
   arr->capacity = new_capacity;
// Allocate a new memory, check its correctness.
   void **new_pp_obj_ptr = (void **)ae_malloc((size_t)arr->capacity * sizeof *new_pp_obj_ptr);
   ae_int_t *new_pp_obj_sizes = (ae_int_t *)ae_malloc((size_t)arr->capacity * sizeof *new_pp_obj_sizes);
   ae_copy_op *new_pp_copy = (ae_copy_op *)ae_malloc((size_t)arr->capacity * sizeof *new_pp_copy);
   ae_free_op *new_pp_free = (ae_free_op *)ae_malloc((size_t)arr->capacity * sizeof *new_pp_free);
   if (new_pp_obj_ptr == NULL || new_pp_obj_sizes == NULL || new_pp_copy == NULL || new_pp_free == NULL) {
   // malloc error: free all newly allocated memory, return.
      ae_free(new_pp_obj_ptr);
      ae_free(new_pp_obj_sizes);
      ae_free(new_pp_copy);
      ae_free(new_pp_free);
      return false;
   }
// Move the data.
   memmove(new_pp_obj_ptr, arr->pp_obj_ptr, (size_t)arr->cnt * sizeof *new_pp_obj_ptr);
   memmove(new_pp_obj_sizes, arr->pp_obj_sizes, (size_t)arr->cnt * sizeof *new_pp_obj_sizes);
   memmove(new_pp_copy, arr->pp_copy, (size_t)arr->cnt * sizeof *new_pp_copy);
   memmove(new_pp_free, arr->pp_free, (size_t)arr->cnt * sizeof *new_pp_free);
// Free the old memory, swap the pointers.
   ae_free(arr->pp_obj_ptr), arr->pp_obj_ptr = new_pp_obj_ptr;
   ae_free(arr->pp_obj_sizes), arr->pp_obj_sizes = new_pp_obj_sizes;
   ae_free(arr->pp_copy), arr->pp_copy = new_pp_copy;
   ae_free(arr->pp_free), arr->pp_free = new_pp_free;
// Done.
   return true;
}

// Configure the array into special fixed capacity mode, which allows concurrent appends, writes and reads to be performed.
//	arr:		the array, can be in any mode - dynamic or fixed capacity,
//	new_capacity:	the new capacity, must be at least as long as the current length.
// On output:
// *	the array capacity is increased to new_capacity exactly,
// *	all present elements are retained,
// *	if the array larger than new_capacity, an exception is thrown.
void ae_obj_array_fixed_capacity(ae_obj_array *arr, ae_int_t new_capacity) {
   ae_assert(arr->cnt <= new_capacity, "ae_obj_array_fixed_capacity: new capacity is less than present size");
   if (!_ae_obj_array_set_capacity(arr, new_capacity))
      ae_assert(false, "ae_obj_array_fixed_capacity: memory error during reallocation");
   arr->fixed_capacity = true;
}

// Retrieve a given element from the array and assign it to a smart pointer.
// On output:
// *	the pointer with index idx is assigned to ptr,
// *	ptr does NOT own new pointer,
// *	whatever pointer ptr owns beforehand (if any) will be properly deallocated,
// *	out-of-bounds accesses will result in an exception being thrown.
void ae_obj_array_get(ae_obj_array *arr, ae_int_t idx, ae_smart_ptr *ptr) {
   if (idx < 0 || idx >= arr->cnt)
      ae_assert(false, "ObjArray: out of bounds read access was performed");
   ae_smart_ptr_assign(ptr, arr->pp_obj_ptr[idx], false, false, 0, NULL, NULL);
}

// Set element idx of arr to the pointer ptr.
// Notes:
// *	the array must larger than idx in size, an exception will be thrown otherwise,
// *	ptr can be NULL
// *	non-NULL ptr MUST own its value before calling this function,
//	and it will transfer ownership to arr afterwards (although it will still point to the object),
// *	non-NULL ptr must point to a dynamically allocated object (on-stack objects are not supported),
// *	out-of-bounds access will result in exception being thrown,
// *	this function does NOT change array the size and capacity.
// This function thread-safe as long as the array capacity is not changed by concurrently called functions.
void ae_obj_array_set_transfer(ae_obj_array *arr, ae_int_t idx, ae_smart_ptr *ptr) {
// Initial integrity checks.
   if (idx < 0 || idx >= arr->cnt) ae_assert(false, "ae_obj_array_set_transfer: out of bounds idx");
   ae_assert(ptr->ptr == NULL || ptr->is_owner, "ae_obj_array_set_transfer: ptr does not own its pointer");
   ae_assert(ptr->ptr == NULL || ptr->is_dynamic, "ae_obj_array_set_transfer: ptr does not point to dynamic object");
// Clean up existing pointer at location index, if needed.
   if (arr->pp_obj_ptr[idx] != NULL) {
      arr->pp_free[idx](arr->pp_obj_ptr[idx], false);
      ae_free(arr->pp_obj_ptr[idx]);
      arr->pp_obj_ptr[idx] = NULL;
      arr->pp_obj_sizes[idx] = 0;
      arr->pp_copy[idx] = NULL;
      arr->pp_free[idx] = NULL;
   }
// If ptr is non-NULL, transfer it to the array.
   if (ptr->ptr != NULL) {
   // Move it to the array.
      arr->pp_obj_ptr[idx] = ptr->ptr;
      arr->pp_obj_sizes[idx] = ptr->size_of_object;
      arr->pp_copy[idx] = ptr->copy;
      arr->pp_free[idx] = ptr->free;
   // Release ownership.
      ptr->is_owner = false;
      ptr->is_dynamic = false;
      ptr->size_of_object = 0;
      ptr->copy = NULL;
      ptr->free = NULL;
   }
}

// A portable full memory fence.
// The current implementation issues a fence by acquiring and releasing a lock;
// future implementations may add specialized code supported by the current compiler and/or OS.
static void ae_mfence(ae_lock *lock) {
   ae_acquire_lock(lock);
   ae_release_lock(lock);
}

// Atomically append a pointer to arr, increasing the array length by 1, and return the index of the element being added.
// Notes:
// *	if the array has fixed capacity and its size is already at its limit, an exception will be thrown,
// *	ptr can be NULL,
// *	non-NULL P MUST own its value before calling this function,
//	and it will transfer ownership to Vec afterwards (although it will still point to the object),
// *	non-NULL P must point to a dynamically allocated object (on-stack objects are not supported),
// This function is partially thread-safe:
// *	parallel threads can concurrently append elements using this function
// *	for fixed-capacity arrays it is possible to combine appends with reads, e.g. to use ae_obj_array_get()
ae_int_t ae_obj_array_append_transfer(ae_obj_array *arr, ae_smart_ptr *ptr) {
// Initial integrity checks.
   ae_assert(ptr->ptr == NULL || ptr->is_owner, "ae_obj_array_append_transfer: ptr does not own its pointer");
   ae_assert(ptr->ptr == NULL || ptr->is_dynamic, "ae_obj_array_append_transfer: ptr does not point to dynamic object");
   ae_assert(!arr->fixed_capacity || arr->cnt < arr->capacity, "ae_obj_array_append_transfer: unable to append, all capacity is used up");
// Get the primary lock.
   ae_acquire_lock(&arr->array_lock);
// Reallocate, if needed.
   if (arr->cnt == arr->capacity) {
   // One more integrity check.
      AE_CRITICAL_ASSERT(!arr->fixed_capacity);
   // Increase the capacity.
      if (!_ae_obj_array_set_capacity(arr, 2 * arr->capacity + 8)) {
      // malloc error: release the lock and throw an exception.
         ae_release_lock(&arr->array_lock);
         ae_assert(false, "ae_obj_array_append_transfer: malloc error");
      }
   }
// Append ptr.
   if (ptr->ptr != NULL) {
   // Move it to the array.
      arr->pp_obj_ptr[arr->cnt] = ptr->ptr;
      arr->pp_obj_sizes[arr->cnt] = ptr->size_of_object;
      arr->pp_copy[arr->cnt] = ptr->copy;
      arr->pp_free[arr->cnt] = ptr->free;
   // Release ownership.
      ptr->is_owner = false;
      ptr->is_dynamic = false;
      ptr->size_of_object = 0;
      ptr->copy = NULL;
      ptr->free = NULL;
   } else {
   // Set to NULL.
      arr->pp_obj_ptr[arr->cnt] = NULL;
      arr->pp_obj_sizes[arr->cnt] = 0;
      arr->pp_copy[arr->cnt] = NULL;
      arr->pp_free[arr->cnt] = NULL;
   }
// Issue a memory fence (necessary for the correct ae_obj_array_get_length) and increase the array size.
   ae_mfence(&arr->barrier_lock);
   ae_int_t result = arr->cnt;
   arr->cnt = result + 1;
// Release the primary lock.
   ae_release_lock(&arr->array_lock);
// Done.
   return result;
}

// Convert the six-bit value v in the range [0,0100) to digits, letters, minuses and underscores.
// Any v outside the range [0,0100) is convereted to '?'.
static char ae_sixbits2char(ae_int_t v) {
   static char _sixbits2char_tbl[0100] = {
      '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
      'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
      'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
      'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '-', '_'
   };
   const size_t CharN = sizeof _sixbits2char_tbl/sizeof _sixbits2char_tbl[0];
   return v >= 0 && v < CharN ? _sixbits2char_tbl[v] : '?';
#if 0
// v is correct, process it.
   return
      v < 10 ? '0' + v :
      (v -= 10) < 26 ? 'A' + v :
      (v -= 26) < 26 ? 'a' + v :
      (v -= 26) == 0 ? '-' : '_';
#endif
}

// The inverse of ae_sixbits2char().
// Convert any character c in the range of ae_sixbits2char() to a six-bit value in the range [0, 0100).
// Convert any other character c to -1.
static ae_int_t ae_char2sixbits(char c) {
   static ae_int_t _ae_char2sixbits_tbl[0x80] = {
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1,
       0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1,
      -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
      25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, -1, -1, -1, -1, 63,
      -1, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, -1, -1, -1, -1, -1
   };
   const size_t Bit6N = sizeof _ae_char2sixbits_tbl/sizeof _ae_char2sixbits_tbl[0];
   return c >= 0 && c < Bit6N ? _ae_char2sixbits_tbl[(int)c] : -1;
}

// Convert the 3 x 8-bit array src into the 4 x 6-bit array dst.
static void ae_threebytes2foursixbits(const unsigned char *src, ae_int_t *dst) {
   dst[0] = src[0] & 0x3F;
   dst[1] = (src[0] >> 6) | ((src[1] & 0x0F) << 2);
   dst[2] = (src[1] >> 4) | ((src[2] & 0x03) << 4);
   dst[3] = src[2] >> 2;
}

// Convert the 4 x 6-bit array src into the 3 x 8-bit array dst.
static void ae_foursixbits2threebytes(const ae_int_t *src, unsigned char *dst) {
   dst[0] = (unsigned char)(src[0] | ((src[1] & 0x03) << 6));
   dst[1] = (unsigned char)((src[1] >> 2) | ((src[2] & 0x0F) << 4));
   dst[2] = (unsigned char)((src[2] >> 4) | (src[3] << 2));
}

// Initialize the serializer.
void ae_serializer_init(ae_serializer *serializer) {
   serializer->mode = AE_SM_DEFAULT;
   serializer->bytes_asked = serializer->entries_needed = 0;
}

void ae_serializer_alloc_start(ae_serializer *serializer) {
   serializer->mode = AE_SM_ALLOC;
   serializer->bytes_asked = serializer->entries_needed = 0;
}

void ae_serializer_alloc_entry(ae_serializer *serializer) {
   serializer->entries_needed++;
}

// After the allocation is done, return the required size of the output string buffer (including the trailing '\0').
// The actual size of the data being stored can be a few characters smaller than requested.
size_t ae_serializer_get_alloc_size(ae_serializer *serializer) {
   serializer->mode = AE_SM_READY2S;
// If no entries are needed (the degenerate case).
   if (serializer->entries_needed == 0)
      return serializer->bytes_asked = 4; // A pair of chars for "\r\n", one for '.', one for the trailing '\0'.
// The non-degenerate case.
   ae_int_t rows = serializer->entries_needed / AE_SER_ENTRIES_PER_ROW;
   ae_int_t lastrowsize = AE_SER_ENTRIES_PER_ROW;
   if (serializer->entries_needed % AE_SER_ENTRIES_PER_ROW) {
      lastrowsize = serializer->entries_needed % AE_SER_ENTRIES_PER_ROW;
      rows++;
   }
// The result size: data size + ' ' & '\n' symbols + trailing '.' & '\0'.
   return serializer->bytes_asked = ((rows - 1) * AE_SER_ENTRIES_PER_ROW + lastrowsize) * (AE_SER_ENTRY_LENGTH + 1) + rows + 2;
}

// Determine the byte order.
// Only big-endian and little-endian are supported for ByteOrder, not mixed-endian hardware.
static ae_int_t GetByteOrder() {
// 1983 is used as the magic number because its non-periodic double representation
// allows us to easily distinguish between the upper and lower halves and to detect mixed endian hardware.
   union { double a; ae_int32_t p[2]; } u; u.a = 1.0 / 1983.0;
   return
      u.p[1] == (ae_int32_t)0x3f408642 ? AE_LITTLE_ENDIAN :
      u.p[0] == (ae_int32_t)0x3f408642 ? AE_BIG_ENDIAN :
      AE_MIXED_ENDIAN; //(@) Originally, this prompted an abort().
}
static const ae_int_t ByteOrder = GetByteOrder();

#ifdef AE_USE_CPP_SERIALIZATION
void ae_serializer_sstart_str(ae_serializer *serializer, std::string *buf) {
   serializer->mode = AE_SM_TO_CPPSTRING;
   serializer->out_cppstr = buf;
   serializer->bytes_written = serializer->entries_saved = 0;
}

void ae_serializer_ustart_str(ae_serializer *serializer, const std::string *buf) {
   serializer->mode = AE_SM_FROM_STRING;
   serializer->in_str = buf->c_str();
}

static bool cpp_reader(ae_int_t aux, ae_int_t cnt, char *p_buf) {
   std::istream *stream = reinterpret_cast<std::istream *>(aux);
   if (cnt <= 0) return false; // Unexpected cnt.
   while (true) {
      int c = stream->get();
      if (c < 0 || c > 0xff) return false; // Failure!
      else if (c != ' ' && c != '\t' && c != '\n' && c != '\r') { p_buf[0] = (char)c; break; }
   }
   for (int k = 1; k < cnt; k++) {
      int c = stream->get();
      if (c < 0 || c > 0xff || c == ' ' || c == '\t' || c == '\n' || c == '\r') return false; // Failure!
      p_buf[k] = (char)c;
   }
   p_buf[cnt] = 0;
   return true; // Success.
}

static bool cpp_writer(const char *p_string, ae_int_t aux) {
   std::ostream *stream = reinterpret_cast<std::ostream *>(aux);
   stream->write(p_string, (std::streamsize)strlen(p_string));
   return !stream->bad();
}

void ae_serializer_sstart_stream(ae_serializer *serializer, std::ostream *stream) {
   serializer->mode = AE_SM_TO_STREAM;
   serializer->stream_writer = cpp_writer;
   serializer->stream_aux = reinterpret_cast<ae_int_t>(stream);
   serializer->bytes_written = serializer->entries_saved = 0;
}

void ae_serializer_ustart_stream(ae_serializer *serializer, const std::istream *stream) {
   serializer->mode = AE_SM_FROM_STREAM;
   serializer->stream_reader = cpp_reader;
   serializer->stream_aux = reinterpret_cast<ae_int_t>(stream);
}
#endif

void ae_serializer_sstart_str(ae_serializer *serializer, char *buf) {
   serializer->mode = AE_SM_TO_STRING;
   serializer->out_str = buf;
   serializer->out_str[0] = 0;
   serializer->bytes_written = serializer->entries_saved = 0;
}

void ae_serializer_ustart_str(ae_serializer *serializer, const char *buf) {
   serializer->mode = AE_SM_FROM_STRING;
   serializer->in_str = buf;
}

void ae_serializer_sstart_stream(ae_serializer *serializer, ae_stream_writer writer, ae_int_t aux) {
   serializer->mode = AE_SM_TO_STREAM;
   serializer->stream_writer = writer;
   serializer->stream_aux = aux;
   serializer->bytes_written = serializer->entries_saved = 0;
}

void ae_serializer_ustart_stream(ae_serializer *serializer, ae_stream_reader reader, ae_int_t aux) {
   serializer->mode = AE_SM_FROM_STREAM;
   serializer->stream_reader = reader;
   serializer->stream_aux = aux;
}

// Unserialize a boolean value from the string buf, up to *pasttheend.
// Raise an error in case an unexpected symbol is found.
static bool ae_str2bool(const char *buf, const char **pasttheend) {
   const char *emsg = "ae_str2bool: unable to read boolean value from stream";
   bool was0 = false;
   bool was1 = false;
   for (; *buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r'; buf++);
   for (; *buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != '\0'; buf++)
      if (*buf == '0') was0 = true;
      else if (*buf == '1') was1 = true;
      else ae_break(ERR_ASSERTION_FAILED, emsg);
   *pasttheend = buf;
   if (was0 == was1) ae_break(ERR_ASSERTION_FAILED, emsg);
   return was1;
}

// Unserialize a boolean value from serializer.
// Leading whitespaces are ignored, a trailing whitespace terminates the string.
// Raise an error in case an unexpected symbol is found.
bool ae_serializer_unserialize_bool(ae_serializer *serializer) {
   switch (serializer->mode) {
      case AE_SM_FROM_STRING: return ae_str2bool(serializer->in_str, &serializer->in_str);
      case AE_SM_FROM_STREAM: {
         char buf[AE_SER_ENTRY_LENGTH + 3];
         const char *p = buf;
         ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf), "ae_serializer_unserialize_bool: error reading from stream");
         return ae_str2bool(buf, &p);
      }
      default: ae_break(ERR_ASSERTION_FAILED, "ae_serializer_unserialize_bool: unable to read boolean value from stream");
   }
   return false;
}

// Serialize the boolean value v into the string buf.
static void ae_bool2str(bool v, char *buf) {
   char c = v ? '1' : '0';
   for (ae_int_t i = 0; i < AE_SER_ENTRY_LENGTH; i++) buf[i] = c;
   buf[AE_SER_ENTRY_LENGTH] = '\0';
}

// Serialize the boolean value v to serializer.
void ae_serializer_serialize_bool(ae_serializer *serializer, bool v) {
   const char *emsg = "ae_serializer_serialize_bool: serialization integrity error";
// At least 11 characters for the value, plus 1 for the trailing '\0'.
   char buf[AE_SER_ENTRY_LENGTH + 3];
// Prepare the serialization and check consistency.
   ae_bool2str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   ae_int_t bytes_appended = (ae_int_t)strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // Strict "<", to make room for the trailing '\0'.
   serializer->bytes_written += bytes_appended;
// Append to the buffer.
   switch (serializer->mode) {
#ifdef AE_USE_CPP_SERIALIZATION
      case AE_SM_TO_CPPSTRING:
         *serializer->out_cppstr += buf;
      break;
#endif
      case AE_SM_TO_STRING:
         strcat(serializer->out_str, buf);
         serializer->out_str += bytes_appended;
      break;
      case AE_SM_TO_STREAM:
         ae_assert(serializer->stream_writer(buf, serializer->stream_aux), "ae_serializer_serialize_bool: error writing to stream");
      break;
      default: ae_break(ERR_ASSERTION_FAILED, emsg);
   }
}

// Unserialize an integer value from the string buf, up to *pasttheend.
// Raise an error in case an unexpected symbol is found.
static ae_int_t ae_str2int(const char *buf, const char **pasttheend) {
   const char *emsg = "ae_str2int: unable to read integer value from stream";
// *	Skip leading spaces.
// *	Read and decode six-bit digits.
// *	Set trailing digits to zeros.
// *	Convert to little endian 64-bit integer representation.
// *	Convert to big endian representation, if needed.
   for (; *buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r'; buf++);
   ae_int_t sixbits[12];
   ae_int_t sixbitsread = 0;
   for (; *buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != '\0'; buf++) {
      ae_int_t d = ae_char2sixbits(*buf);
      if (d < 0 || sixbitsread >= AE_SER_ENTRY_LENGTH) ae_break(ERR_ASSERTION_FAILED, emsg);
      sixbits[sixbitsread++] = d;
   }
   *pasttheend = buf;
   if (sixbitsread == 0) ae_break(ERR_ASSERTION_FAILED, emsg);
   for (ae_int_t i = sixbitsread; i < 12; i++) sixbits[i] = 0;
   union { ae_int_t ival; unsigned char bytes[9]; } u;
   ae_foursixbits2threebytes(sixbits, u.bytes);
   ae_foursixbits2threebytes(sixbits + 4, u.bytes + 3);
   ae_foursixbits2threebytes(sixbits + 8, u.bytes + 6);
   if (ByteOrder == AE_BIG_ENDIAN)
      for (ae_int_t i0 = 0, i1 = sizeof u.ival - 1; i0 < i1; i0++, i1--) {
         unsigned char tc = u.bytes[i0];
         u.bytes[i0] = u.bytes[i1];
         u.bytes[i1] = tc;
      }
   return u.ival;
}

// Unserialize an integer value from serializer.
// Leading whitespaces are ignored, a trailing whitespace terminates the string.
// Raise an error in case an unexpected symbol is found.
ae_int_t ae_serializer_unserialize_int(ae_serializer *serializer) {
   switch (serializer->mode) {
      case AE_SM_FROM_STRING: return ae_str2int(serializer->in_str, &serializer->in_str);
      case AE_SM_FROM_STREAM: {
         char buf[AE_SER_ENTRY_LENGTH + 3];
         const char *p = buf;
         ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf), "ae_serializer_unserialize_int: error reading from stream");
         return ae_str2int(buf, &p);
      }
      default: ae_break(ERR_ASSERTION_FAILED, "ae_serializer_unserialize_int: integrity check failed");
   }
   return 0;
}

// Serialize the integer value v into the string buf.
static void ae_int2str(ae_int_t v, char *buf) {
// Copy v to array of chars, sign extending it and converting it to little endian order.
// In order to avoid explicit mention of the size of the ae_int_t type, we:
// *	fill u.bytes by zeros or ones (depending on the sign of v),
// *	copy v to u.ival,
// *	reorder u.bytes (for big endian architectures), to obtain a signed 64-bit representation of v stored in u.bytes,
// *	additionally, zero the last byte of u.bytes to simplify conversion to six-bit form.
   unsigned char c = v < 0 ? (unsigned char)0xFF : (unsigned char)0x00;
   union { ae_int_t ival; unsigned char bytes[9]; } u;
   u.ival = v;
   const size_t un = sizeof u.ival, us = sizeof u.bytes;
   for (ae_int_t i = un; i < us; i++) // i < us is preferred because it avoids unnecessary compiler warnings.
      u.bytes[i] = c;
   u.bytes[us - 1] = '\0';
   if (ByteOrder == AE_BIG_ENDIAN)
      for (ae_int_t i0 = 0, i1 = un - 1; i0 < i1; i0++, i1--) {
         unsigned char tc = u.bytes[i0];
         u.bytes[i0] = u.bytes[i1];
         u.bytes[i1] = tc;
      }
// Convert to six-bit representation, output.
// NOTE:
// *	The last element of sixbits is always zero, and is not output.
   ae_int_t sixbits[12];
   ae_threebytes2foursixbits(u.bytes, sixbits);
   ae_threebytes2foursixbits(u.bytes + 3, sixbits + 4);
   ae_threebytes2foursixbits(u.bytes + 6, sixbits + 8);
   for (ae_int_t i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = ae_sixbits2char(sixbits[i]);
   buf[AE_SER_ENTRY_LENGTH] = 0x00;
}

// Serialize the integer value v to serializer.
void ae_serializer_serialize_int(ae_serializer *serializer, ae_int_t v) {
   const char *emsg = "ae_serializer_serialize_int: serialization integrity error";
// At least 11 characters for the value, plus 1 for the trailing '\0'.
   char buf[AE_SER_ENTRY_LENGTH + 3];
// Prepare the serialization and check consistency.
   ae_int2str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   ae_int_t bytes_appended = (ae_int_t)strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // Strict "<", to make room for the trailing '\0'.
   serializer->bytes_written += bytes_appended;
// Append to the buffer.
   switch (serializer->mode) {
#ifdef AE_USE_CPP_SERIALIZATION
      case AE_SM_TO_CPPSTRING:
         *serializer->out_cppstr += buf;
      break;
#endif
      case AE_SM_TO_STRING:
         strcat(serializer->out_str, buf);
         serializer->out_str += bytes_appended;
      break;
      case AE_SM_TO_STREAM:
         ae_assert(serializer->stream_writer(buf, serializer->stream_aux), "ae_serializer_serialize_int: error writing to stream");
      break;
      default: ae_break(ERR_ASSERTION_FAILED, emsg);
   }
}

// Unserialize a 64-bit integer value from the string buf, up to *pasttheend.
// Raise an error in case an unexpected symbol is found.
static ae_int64_t ae_str2int64(const char *buf, const char **pasttheend) {
   const char *emsg = "ae_str2int64: unable to read integer value from stream";
// *	Skip leading spaces.
// *	Read and decode six-bit digits.
// *	Set trailing digits to zeros.
// *	Convert to little endian 64-bit integer representation.
// *	Convert to big endian representation, if needed.
   for (; *buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r'; buf++);
   ae_int_t sixbits[12];
   ae_int_t sixbitsread = 0;
   for (; *buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != '\0'; buf++) {
      ae_int_t d = ae_char2sixbits(*buf);
      if (d < 0 || sixbitsread >= AE_SER_ENTRY_LENGTH) ae_break(ERR_ASSERTION_FAILED, emsg);
      sixbits[sixbitsread++] = d;
   }
   *pasttheend = buf;
   if (sixbitsread == 0) ae_break(ERR_ASSERTION_FAILED, emsg);
   for (ae_int_t i = sixbitsread; i < 12; i++) sixbits[i] = 0;
   unsigned char bytes[9];
   ae_foursixbits2threebytes(sixbits, bytes);
   ae_foursixbits2threebytes(sixbits + 4, bytes + 3);
   ae_foursixbits2threebytes(sixbits + 8, bytes + 6);
   if (ByteOrder == AE_BIG_ENDIAN)
      for (ae_int_t i0 = 0, i1 = sizeof(ae_int_t) - 1; i0 < i1; i0++, i1--) {
         unsigned char tc = bytes[i0];
         bytes[i0] = bytes[i1];
         bytes[i1] = tc;
      }
   ae_int64_t result;
   memmove(&result, bytes, sizeof result);
   return result;
}

// Unserialize a 64-bit integer value from serializer.
// Leading whitespaces are ignored, a trailing whitespace terminates the string.
// Raise an error in case an unexpected symbol is found.
ae_int64_t ae_serializer_unserialize_int64(ae_serializer *serializer) {
   switch (serializer->mode) {
      case AE_SM_FROM_STRING: return ae_str2int64(serializer->in_str, &serializer->in_str);
      case AE_SM_FROM_STREAM: {
         char buf[AE_SER_ENTRY_LENGTH + 3];
         const char *p = buf;
         ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf), "ae_serializer_unserialize_int64: error reading from stream");
         return ae_str2int64(buf, &p);
      }
      default: ae_break(ERR_ASSERTION_FAILED, "ae_serializer_unserialize_int64: integrity check failed");
   }
   return 0L;
}

// Serialize the 64-bit integer value v into the string buf.
static void ae_int642str(ae_int64_t v, char *buf) {
// Copy v to an array of chars, sign extending it and converting it to little endian order.
// In order to avoid explicit mention of the size of the ae_int64_t type, we:
// *	fill bytes by zeros or ones (depending on the sign of v),
// *	copy v to bytes,
// *	reorder bytes (for big endian architectures), to obtain a signed 64-bit representation of v stored in bytes,
// *	additionally, zero the last byte of bytes to simplify conversion to six-bit form.
   unsigned char bytes[9];
   memset(bytes, v < 0 ? 0xFF : 0x00, 8);
   memmove(bytes, &v, 8);
   bytes[8] = 0;
   if (ByteOrder == AE_BIG_ENDIAN)
      for (ae_int_t i0 = 0, i1 = sizeof(ae_int_t) - 1; i0 < i1; i0++, i1--) {
         unsigned char tc = bytes[i0];
         bytes[i0] = bytes[i1];
         bytes[i1] = tc;
      }
// Convert to six-bit representation, output.
// NOTE:
// *	The last element of sixbits is always zero, we do not output it.
   ae_int_t sixbits[12];
   ae_threebytes2foursixbits(bytes, sixbits);
   ae_threebytes2foursixbits(bytes + 3, sixbits + 4);
   ae_threebytes2foursixbits(bytes + 6, sixbits + 8);
   for (ae_int_t i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = ae_sixbits2char(sixbits[i]);
   buf[AE_SER_ENTRY_LENGTH] = 0x00;
}

// Serialize the 64-bit integer value v to serializer.
void ae_serializer_serialize_int64(ae_serializer *serializer, ae_int64_t v) {
   const char *emsg = "ae_serializer_serialize_int64: serialization integrity error";
// At least 11 characters for the value, plus 1 for the trailing '\0'.
   char buf[AE_SER_ENTRY_LENGTH + 3];
// Prepare the serialization and check consistency.
   ae_int642str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   ae_int_t bytes_appended = (ae_int_t)strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // Strict "<", to make room for the trailing '\0'.
   serializer->bytes_written += bytes_appended;
// Append to the buffer.
   switch (serializer->mode) {
#ifdef AE_USE_CPP_SERIALIZATION
      case AE_SM_TO_CPPSTRING:
         *serializer->out_cppstr += buf;
      break;
#endif
      case AE_SM_TO_STRING:
         strcat(serializer->out_str, buf);
         serializer->out_str += bytes_appended;
      break;
      case AE_SM_TO_STREAM:
         ae_assert(serializer->stream_writer(buf, serializer->stream_aux), "ae_serializer_serialize_int64: error writing to stream");
      break;
      default: ae_break(ERR_ASSERTION_FAILED, emsg);
   }
}

// Unserialize a double value from the string buf, up to *pasttheend.
// Raise an error in case an unexpected symbol is found.
static double ae_str2double(const char *buf, const char **pasttheend) {
   const char *emsg = "ae_str2double: unable to read double value from stream";
// Skip leading spaces.
   for (; *buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r'; buf++);
// Handle special cases.
   if (*buf == '.') {
      const char *s_nan = ".nan_______";
      const char *s_posinf = ".posinf____";
      const char *s_neginf = ".neginf____";
      if (strncmp(buf, s_nan, strlen(s_nan)) == 0) { *pasttheend = buf + strlen(s_nan); return NAN; }
      if (strncmp(buf, s_posinf, strlen(s_posinf)) == 0) { *pasttheend = buf + strlen(s_posinf); return +INFINITY; }
      if (strncmp(buf, s_neginf, strlen(s_neginf)) == 0) { *pasttheend = buf + strlen(s_neginf); return -INFINITY; }
      ae_break(ERR_ASSERTION_FAILED, emsg);
   }
// The general case:
// *	Read and decode six-bit digits.
// *	Check that all 11 digits were read.
// *	Set the last digit to zero (needed to simplify the conversion).
// *	Convert to 8 bytes.
// *	Convert to big endian representation, if needed.
   ae_int_t sixbits[12];
   ae_int_t sixbitsread = 0;
   for (; *buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != '\0'; buf++) {
      ae_int_t d = ae_char2sixbits(*buf);
      if (d < 0 || sixbitsread >= AE_SER_ENTRY_LENGTH) ae_break(ERR_ASSERTION_FAILED, emsg);
      sixbits[sixbitsread++] = d;
   }
   *pasttheend = buf;
   if (sixbitsread != AE_SER_ENTRY_LENGTH) ae_break(ERR_ASSERTION_FAILED, emsg);
   sixbits[AE_SER_ENTRY_LENGTH] = 0;
   union { double dval; unsigned char bytes[9]; } u;
   ae_foursixbits2threebytes(sixbits, u.bytes);
   ae_foursixbits2threebytes(sixbits + 4, u.bytes + 3);
   ae_foursixbits2threebytes(sixbits + 8, u.bytes + 6);
   if (ByteOrder == AE_BIG_ENDIAN)
      for (ae_int_t i0 = 0, i1 = sizeof u.dval - 1; i0 < i1; i0++, i1--) {
         unsigned char tc = u.bytes[i0];
         u.bytes[i0] = u.bytes[i1];
         u.bytes[i1] = tc;
      }
   return u.dval;
}

// Unserialize a double value from serializer.
// Leading whitespaces are ignored, a trailing whitespace terminates the string.
// Raise an error in case an unexpected symbol is found.
double ae_serializer_unserialize_double(ae_serializer *serializer) {
   switch (serializer->mode) {
      case AE_SM_FROM_STRING: return ae_str2double(serializer->in_str, &serializer->in_str);
      case AE_SM_FROM_STREAM: {
         char buf[AE_SER_ENTRY_LENGTH + 3];
         const char *p = buf;
         ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf), "ae_serializer_unserialize_double: error reading from stream");
         return ae_str2double(buf, &p);
      }
      default: ae_break(ERR_ASSERTION_FAILED, "ae_serializer_unserialize_double: unable to read double value from stream");
   }
   return 0.0;
}

// Serialize the double value v into the string buf.
static void ae_double2str(double v, char *buf) {
// Handle the special cases, first.
   if (isnan(v)) {
      const char *s = ".nan_______";
      memmove(buf, s, strlen(s) + 1);
      return;
   } else if (isposinf(v)) {
      const char *s = ".posinf____";
      memmove(buf, s, strlen(s) + 1);
      return;
   } else if (isneginf(v)) {
      const char *s = ".neginf____";
      memmove(buf, s, strlen(s) + 1);
      return;
   }
// Handle the general case:
// *	copy v to array of chars,
// *	set the last byte of u.bytes to zero in order to simplify conversion to six-bit representation,
// *	convert to little endian (if needed),
// *	convert to six-bit representation (the last element of sixbits is always zero, and is not output).
   union { double dval; unsigned char bytes[9]; } u;
   u.dval = v;
   u.bytes[8] = '\0';
   if (ByteOrder == AE_BIG_ENDIAN)
      for (ae_int_t i0 = 0, i1 = sizeof u.dval - 1; i0 < i1; i0++, i1--) {
         unsigned char tc = u.bytes[i0];
         u.bytes[i0] = u.bytes[i1];
         u.bytes[i1] = tc;
      }
   ae_int_t sixbits[12];
   ae_threebytes2foursixbits(u.bytes, sixbits);
   ae_threebytes2foursixbits(u.bytes + 3, sixbits + 4);
   ae_threebytes2foursixbits(u.bytes + 6, sixbits + 8);
   for (ae_int_t i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = ae_sixbits2char(sixbits[i]);
   buf[AE_SER_ENTRY_LENGTH] = 0x00;
}

void ae_serializer_serialize_double(ae_serializer *serializer, double v) {
   const char *emsg = "ae_serializer_serialize_double: serialization integrity error";
// At least 11 characters for the value, plus 1 for the trailing '\0'.
   char buf[AE_SER_ENTRY_LENGTH + 3];
// Prepare the serialization and check consistency.
   ae_double2str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   ae_int_t bytes_appended = (ae_int_t)strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // Strict "<", to make room for the trailing '\0'.
   serializer->bytes_written += bytes_appended;
// Append to the buffer.
   switch (serializer->mode) {
#ifdef AE_USE_CPP_SERIALIZATION
      case AE_SM_TO_CPPSTRING:
         *serializer->out_cppstr += buf;
      break;
#endif
      case AE_SM_TO_STRING:
         strcat(serializer->out_str, buf);
         serializer->out_str += bytes_appended;
      break;
      case AE_SM_TO_STREAM:
         ae_assert(serializer->stream_writer(buf, serializer->stream_aux), "ae_serializer_serialize_double: error writing to stream");
      break;
      default: ae_break(ERR_ASSERTION_FAILED, emsg);
   }
}

void ae_serializer_stop(ae_serializer *serializer) {
   switch (serializer->mode) {
#ifdef AE_USE_CPP_SERIALIZATION
      case AE_SM_TO_CPPSTRING:
         ae_assert(serializer->bytes_written + 1 < serializer->bytes_asked, "ae_serializer_stop: integrity check failed"); // Strict "<", to make room for the trailing '\0'.
         serializer->bytes_written++;
         *serializer->out_cppstr += ".";
      break;
#endif
      case AE_SM_TO_STRING:
         ae_assert(serializer->bytes_written + 1 < serializer->bytes_asked, "ae_serializer_stop: integrity check failed"); // Strict "<", to make room for the trailing '\0'.
         serializer->bytes_written++;
         strcat(serializer->out_str, ".");
         serializer->out_str++;
      break;
      case AE_SM_TO_STREAM:
         ae_assert(serializer->bytes_written + 1 < serializer->bytes_asked, "ae_serializer_stop: integrity check failed"); // Strict "<", to make room for the trailing '\0'.
         serializer->bytes_written++;
         ae_assert(serializer->stream_writer(".", serializer->stream_aux), "ae_serializer_stop: error writing to stream");
      break;
   // For compatibility with the pre-3.11 serializer, which does not require a trailing '.', we do not test for a trailing '.'.
   // Anyway, because the string is not a stream, we do not have to read ALL trailing symbols.
      case AE_SM_FROM_STRING: break;
      case AE_SM_FROM_STREAM: {
      // Read a trailing '.', perform an integrity check.
         char buf[2];
         ae_assert(serializer->stream_reader(serializer->stream_aux, 1, buf), "ae_serializer_stop: error reading from stream");
         ae_assert(buf[0] == '.', "ae_serializer_stop: trailing . is not found in the stream");
      }
      break;
      default: ae_break(ERR_ASSERTION_FAILED, "ae_serializer_stop: integrity check failed");
   }
}

void ae_serializer_alloc_byte_array(ae_serializer *serializer, ae_vector *bytes) {
   ae_int_t n = bytes->cnt;
   serializer->entries_needed += 1 + n / 8 + (n % 8 > 0 ? 1 : 0);
}

void ae_serializer_unserialize_byte_array(ae_serializer *serializer, ae_vector *bytes) {
   ae_int_t chunk_size = 8;
// Read the array length, allocate output.
   ae_int_t n = ae_serializer_unserialize_int(serializer);
   ae_vector_set_length(bytes, n);
// Count and read the entries.
   ae_int_t entries_count = n / chunk_size + (n % chunk_size > 0 ? 1 : 0);
   for (ae_int_t eidx = 0; eidx < entries_count; eidx++) {
      ae_int64_t tmp64 = ae_serializer_unserialize_int64(serializer);
      memmove(bytes->xU + eidx * chunk_size, &tmp64, imin2(n - eidx * chunk_size, chunk_size));
   }
}

void ae_serializer_serialize_byte_array(ae_serializer *serializer, ae_vector *bytes) {
   ae_int_t chunk_size = 8;
// Save the array length.
   ae_serializer_serialize_int(serializer, bytes->cnt);
// Count the entries.
   ae_int_t entries_count = bytes->cnt / chunk_size + (bytes->cnt % chunk_size > 0 ? 1 : 0);
   for (ae_int_t eidx = 0; eidx < entries_count; eidx++) {
      ae_int64_t tmpi; memset(&tmpi, 0, sizeof tmpi);
      memmove(&tmpi, bytes->xU + eidx * chunk_size, imin2(bytes->cnt - eidx * chunk_size, chunk_size));
      ae_serializer_serialize_int64(serializer, tmpi);
   }
}

// Integer and real math functions:
bool isneginf(double x) { return isinf(x) && signbit(x); }
bool isposinf(double x) { return isinf(x) && !signbit(x); }

ae_int_t imin2(ae_int_t x, ae_int_t y) { return x > y ? y : x; }

// Return min(x,y,z).
ae_int_t imin3(ae_int_t x, ae_int_t y, ae_int_t z) {
   ae_int_t xy = x < y ? x : y;
   return xy < z ? xy : z;
}

ae_int_t imax2(ae_int_t x, ae_int_t y) { return x > y ? x : y; }

// Return max(x,y,z).
ae_int_t imax3(ae_int_t x, ae_int_t y, ae_int_t z) {
   ae_int_t xy = x > y ? x : y;
   return xy > z ? xy : z;
}

ae_int_t ae_iabs(ae_int_t x) { return x >= 0 ? x : -x; }
ae_int_t sign(double x) { return x > 0.0 ? +1 : x < 0.0 ? -1 : 0; }
ae_int_t iround(double x) { return (ae_int_t)round(x); }
ae_int_t itrunc(double x) { return (ae_int_t)trunc(x); }
ae_int_t ifloor(double x) { return (ae_int_t)floor(x); }
ae_int_t iceil(double x) { return (ae_int_t)ceil(x); }

double rmin2(double x, double y) { return x > y ? y : x; }
double rmax2(double x, double y) { return x > y ? x : y; }

// Return max(x,y,z).
double rmax3(double x, double y, double z) {
   double xy = x > y ? x : y;
   return xy > z ? xy : z;
}

// Return max(|x|,|y|,|z|).
double rmaxabs3(double x, double y, double z) {
   x = fabs(x), y = fabs(y), z = fabs(z);
   double xy = x > y ? x : y;
   return xy > z ? xy : z;
}

double sqr(double x) { return x * x; }

// Clip x to [b1, b2], mapping to the nearest endpoint if outside the interval.
// This assumes b1 < b2.
// ALGLIB: Copyright 20.03.2009 by Sergey Bochkanov
// Integer version.
ae_int_t iboundval(ae_int_t x, ae_int_t b1, ae_int_t b2) {
   return x <= b1 ? b1 : x >= b2 ? b2 : x;
}
// Real version.
double rboundval(double x, double b1, double b2) {
   return x <= b1 ? b1 : x >= b2 ? b2 : x;
}

// RNG wrappers
static ae_int_t ae_rand() {
#if !defined ALGLIB_ENTROPY_SRC || ALGLIB_ENTROPY_SRC == ALGLIB_ENTROPY_SRC_STDRAND
   return (ae_int_t)rand();
#elif ALGLIB_ENTROPY_SRC == ALGLIB_ENTROPY_SRC_OPENSSL
   ae_int32_t random_number;
   unsigned char buf[sizeof random_number];
   if (!RAND_bytes(buf, sizeof random_number)) {
   // openSSL random number generator failed, default to standard random generator
      return (ae_int_t)rand();
   }
   memcpy(&random_number, buf, sizeof random_number);
   if (random_number < 0)
      random_number = -(random_number + 1);
   return (ae_int_t)random_number;
#else
#   error ALGLIB_ENTROPY_SRC is defined, but its value is not recognized.
#endif
}

static ae_int_t ae_rand_max() {
#if !defined ALGLIB_ENTROPY_SRC || ALGLIB_ENTROPY_SRC == ALGLIB_ENTROPY_SRC_STDRAND
   return (ae_int_t)RAND_MAX;
#elif ALGLIB_ENTROPY_SRC == ALGLIB_ENTROPY_SRC_OPENSSL
   return (ae_int_t)ALGLIB_OPENSSL_RAND_MAX;
#else
#   error ALGLIB_ENTROPY_SRC is defined, but its value is not recognized.
#endif
}

double randomreal() {
   const double mx = ae_rand_max() + 1.0;
   double i = ae_rand();
   return (i + ae_rand() / mx) / mx;
}

double randommid() {
   const double mx = ae_rand_max() + 1.0;
   double i = ae_rand();
   return 2.0 * (i + ae_rand() / mx) / mx - 1.0;
}

ae_int_t randominteger(ae_int_t maxv) { return ae_rand() % maxv; }

bool randombool(double p/* = 0.5*/) {
   const double mx = ae_rand_max() + 1.0;
   double i = ae_rand();
   return i + ae_rand() / mx <= p * mx;
}

// Complex math functions:
complex ae_c_neg(complex A) { return complex_from_d(-A.x, -A.y); }
complex conj(complex A) { return complex_from_d(+A.x, -A.y); }
complex csqr(complex A) { double Ax = A.x, Ay = A.y; return complex_from_d(Ax * Ax - Ay * Ay, 2.0 * Ax * Ay); }

double abscomplex(complex A) {
   double Ax = fabs(A.x), Ay = fabs(A.y), Lo = Ax < Ay ? Ax : Ay, Hi = Ax > Ay ? Ax : Ay;
   if (Lo == 0.0) return Hi;
   else {
      double Span = Lo / Hi;
      return Hi * sqrt(1.0 + Span * Span);
   }
}

bool ae_c_eq(complex A, complex B) { return A.x == B.x && A.y == B.y; }
bool ae_c_neq(complex A, complex B) { return A.x != B.x || A.y != B.y; }
complex ae_c_add(complex A, complex B) { return complex_from_d(A.x + B.x, A.y + B.y); }

complex ae_c_mul(complex A, complex B) {
   double Ax = A.x, Ay = A.y, Bx = B.x, By = B.y;
   return complex_from_d(Ax * Bx - Ay * By, Ax * By + Ay * Bx);
}

complex ae_c_sub(complex A, complex B) { return complex_from_d(A.x - B.x, A.y - B.y); }

complex ae_c_div(complex A, complex B) {
   double Ax = A.x, Ay = A.y, Bx = B.x, By = B.y;
   if (fabs(By) < fabs(Bx)) {
      double e = By / Bx, f = Bx + By * e;
      return complex_from_d((Ax + Ay * e) / f, (Ay - Ax * e) / f);
   } else {
      double e = Bx / By, f = By + Bx * e;
      return complex_from_d((Ax * e + Ay) / f, (Ay * e - Ax) / f);
   }
}

bool ae_c_eq_d(complex A, double B) { return A.x == B && A.y == 0.0; }
bool ae_c_neq_d(complex A, double B) { return A.x != B || A.y != 0.0; }

complex ae_c_add_d(complex A, double B) { return complex_from_d(A.x + B, A.y); }
complex ae_c_mul_d(complex A, double B) { return complex_from_d(A.x * B, A.y * B); }
complex ae_c_sub_d(complex A, double B) { return complex_from_d(A.x - B, A.y); }
complex ae_c_d_sub(double A, complex B) { return complex_from_d(A - B.x, -B.y); }
complex ae_c_div_d(complex A, double B) { return complex_from_d(A.x / B, A.y / B); }

complex ae_c_d_div(double A, complex B) {
   double Bx = B.x, By = B.y;
   if (fabs(By) < fabs(Bx)) {
      double e = By / Bx, f = Bx + By * e;
      return complex_from_d(A / f, -A * e / f);
   } else {
      double e = Bx / By, f = By + Bx * e;
      return complex_from_d(A * e / f, -A / f);
   }
}

// Level 1 Complex BLAS operations.
// Handle the case of unit stride specially and separately with optimization.

complex ae_v_cdotproduct(const complex *A, ae_int_t dA, const char *CjA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   double Cx = 0.0, Cy = 0.0;
   bool ConjA = CjA[0] != 'N' && CjA[0] != 'n';
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (ConjA)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = -A->y;
            double Bx = B->x, By = -B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = -A->y;
            double Bx = B->x, By = B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = A->y;
            double Bx = B->x, By = -B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = A->y;
            double Bx = B->x, By = B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
   return complex_from_d(Cx, Cy);
}

void ae_v_cmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = B->x;
            A->y = -B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++)
            *A = *B;
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = B->x;
            A->y = -B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
            *A = *B;
}

void ae_v_cmoveneg(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = -B->x;
            A->y = B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = -B->x;
            A->y = -B->y;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = -B->x;
            A->y = B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = -B->x;
            A->y = -B->y;
         }
}

void ae_v_cmoved(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = Alpha * B->x;
            A->y = -Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = Alpha * B->x;
            A->y = Alpha * B->y;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = Alpha * B->x;
            A->y = -Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = Alpha * B->x;
            A->y = Alpha * B->y;
         }
}

void ae_v_cmovec(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   double ax = Alpha.x, ay = Alpha.y;
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = ax * B->x + ay * B->y;
            A->y = -ax * B->y + ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = ax * B->x - ay * B->y;
            A->y = ax * B->y + ay * B->x;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = ax * B->x + ay * B->y;
            A->y = -ax * B->y + ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = ax * B->x - ay * B->y;
            A->y = ax * B->y + ay * B->x;
         }
}

void ae_v_cadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += B->x;
            A->y -= B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += B->x;
            A->y += B->y;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += B->x;
            A->y -= B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += B->x;
            A->y += B->y;
         }
}

void ae_v_caddd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += Alpha * B->x;
            A->y -= Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += Alpha * B->x;
            A->y += Alpha * B->y;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += Alpha * B->x;
            A->y -= Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += Alpha * B->x;
            A->y += Alpha * B->y;
         }
}

void ae_v_caddc(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   double ax = Alpha.x, ay = Alpha.y;
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += ax * B->x + ay * B->y;
            A->y -= ax * B->y - ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += ax * B->x - ay * B->y;
            A->y += ax * B->y + ay * B->x;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += ax * B->x + ay * B->y;
            A->y -= ax * B->y - ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += ax * B->x - ay * B->y;
            A->y += ax * B->y + ay * B->x;
         }
}

void ae_v_csub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x -= B->x;
            A->y += B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x -= B->x;
            A->y -= B->y;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x -= B->x;
            A->y += B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x -= B->x;
            A->y -= B->y;
         }
}

void ae_v_csubd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha) {
   ae_v_caddd(A, dA, B, dB, CjB, N, -Alpha);
}

void ae_v_csubc(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha) {
   ae_v_caddc(A, dA, B, dB, CjB, N, complex_from_d(-Alpha.x, -Alpha.y));
}

void ae_v_cmuld(complex *A, ae_int_t dA, ae_int_t N, double Alpha) {
   if (dA == 1)
      for (ae_int_t n = 0; n < N; n++, A++) {
         A->x *= Alpha;
         A->y *= Alpha;
      }
   else
      for (ae_int_t n = 0; n < N; n++, A += dA) {
         A->x *= Alpha;
         A->y *= Alpha;
      }
}

void ae_v_cmulc(complex *A, ae_int_t dA, ae_int_t N, complex Alpha) {
   double Ax = Alpha.x, Ay = Alpha.y;
   if (dA == 1)
      for (ae_int_t n = 0; n < N; n++, A++) {
         double Bx = A->x, By = A->y;
         A->x = Ax * Bx - Ay * By;
         A->y = Ax * By + Ay * Bx;
      }
   else
      for (ae_int_t n = 0; n < N; n++, A += dA) {
         double Bx = A->x, By = A->y;
         A->x = Ax * Bx - Ay * By;
         A->y = Ax * By + Ay * Bx;
      }
}

// Level 1 Real BLAS operations
// Handle the case of unit stride specially and separately with optimization.
//(@) In ae_v_{dotproduct,move,moveneg,moved,add,addd,sub,subd}(), _ialglib_vcopy()
//(@) this case originally cut the loop counter in half (or quarter) and did two (or four) steps per loop.
//(@) This was removed, the need for it was not clearly explained and it appears to have no visible effect on performance.

double ae_v_dotproduct(const double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   double C = 0.0;
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         C += A[0] * B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         C += *A * *B;
   return C;
}

void ae_v_move(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] = B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A = *B;
}

void ae_v_moveneg(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] = -B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A = -*B;
}

void ae_v_moved(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha) {
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] = Alpha * B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A = Alpha * *B;
}

void ae_v_add(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] += B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A += *B;
}

void ae_v_addd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha) {
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] += Alpha * B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A += Alpha * *B;
}

void ae_v_sub(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] -= B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A -= *B;
}

void ae_v_subd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha) {
// Equivalent to ae_v_addd(A, dA, B, dB, N, -Alpha);
   if (dA == 1 && dB == 1)
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         A[0] -= Alpha * B[0];
   else
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A -= Alpha * *B;
}

void ae_v_muld(double *A, ae_int_t dA, ae_int_t N, double Alpha) {
   if (dA == 1)
      for (ae_int_t n = 0; n < N; n++, A++)
         *A *= Alpha;
   else
      for (ae_int_t n = 0; n < N; n++, A += dA)
         *A *= Alpha;
}

#if 0
// Global and local constants and variables.
const double machineepsilon = 5.0E-16, maxrealnumber = 1.0E300, minrealnumber = 1.0E-300;
const double pi = 3.1415926535897932384626433832795;
#endif

// Optimized shared C/C++ linear algebra code.
static const ae_int_t alglib_simd_alignment = 0x10;
static const ae_int_t alglib_r_block = 0x20, alglib_half_r_block = alglib_r_block / 2, alglib_twice_r_block = 2 * alglib_r_block;
static const ae_int_t alglib_c_block = 0x10, alglib_half_c_block = alglib_c_block / 2, alglib_twice_c_block = 2 * alglib_c_block;

// A fast 32 x 32 real matrix-vector product:
//	y = Beta y + Alpha A x
// with available native optimizations or else generic C code
// IMPORTANT:
// *	a must be stored in row-major order with stride = alglib_r_block, aligned on an alglib_simd_alignment boundary
// *	x must be aligned on an alglib_simd_alignment boundary
// *	y may be non-aligned
static void _ialglib_mv_32(const double *a, const double *x, double *y, ae_int_t stride, double Alpha, double Beta) {
   const double *pa0 = a;
   const double *pa1 = a + alglib_r_block;
   const double *pb = x;
   for (ae_int_t i = 0; i < 16; i++) {
      double v0 = 0.0, v1 = 0.0;
      for (ae_int_t k = 0; k < 4; k++) {
         v0 += pa0[0] * pb[0];
         v1 += pa1[0] * pb[0];
         v0 += pa0[1] * pb[1];
         v1 += pa1[1] * pb[1];
         v0 += pa0[2] * pb[2];
         v1 += pa1[2] * pb[2];
         v0 += pa0[3] * pb[3];
         v1 += pa1[3] * pb[3];
         v0 += pa0[4] * pb[4];
         v1 += pa1[4] * pb[4];
         v0 += pa0[5] * pb[5];
         v1 += pa1[5] * pb[5];
         v0 += pa0[6] * pb[6];
         v1 += pa1[6] * pb[6];
         v0 += pa0[7] * pb[7];
         v1 += pa1[7] * pb[7];
         pa0 += 8;
         pa1 += 8;
         pb += 8;
      }
      y[0] = Beta * y[0] + Alpha * v0;
      y[stride] = Beta * y[stride] + Alpha * v1;
   // Now we've processed rows I and I + 1,
   // Move pa0 and pa1 from rows I + 1 and I + 2 to rows I + 2 and I + 3.
      pa0 += alglib_r_block;
      pa1 += alglib_r_block;
      pb = x;
      y += 2 * stride;
   }
}

#if defined AE_HAS_SSE2_INTRINSICS
// The m x n real matrix-vector product (32 x 32 optimized version):
//	y = Beta y + Alpha a x
// using generic C code, optimized with _ialglib_mv_32 if (m, n) == (32, 32).
// If Beta == 0.0, y is overwritten by Alpha a x and if Alpha == 0.0, y is multiplied by Beta.
// However, the y = Beta y operation here is not done efficiently as a scalar-vector product, and should NOT be used this way.
// IMPORTANT:
// *	(m, n) is in [0, alglib_r_block] x [0, alglib_r_block].
// *	a must be in row-major order with stride = alglib_r_block.
// *	y may be non-aligned
// *	both a and x must have same offset with respect to a 16-byte boundary:
//	either both are aligned, or both are aligned with offset 8.
//	You will crash your system if you try to use this with misaligned or incorrectly aligned data.
// This function supports SSE2; it can be used when:
// *	AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
//	Otherwise, this function will be undefined.
// *	CurCPU contains CPU_SSE2 (checked at run-time)
//	Otherwise, calling this function will probably crash your system.
// If you want to know whether it is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
static void _ialglib_rmv_sse2(ae_int_t m, ae_int_t n, const double *a, const double *x, double *y, ae_int_t stride, double Alpha, double Beta) {
// Handle the special cases:
// *	Alpha == 0.0 or m n == 0
   if (m == 0)
      return;
   else if (Alpha == 0.0 || n == 0) {
      if (Beta == 0.0)
         for (ae_int_t i = 0; i < m; i++) {
            *y = 0.0;
            y += stride;
         }
      else
         for (ae_int_t i = 0; i < m; i++) {
            *y *= Beta;
            y += stride;
         }
      return;
   }
// Handle the general case: m n != 0 and Alpha != 0.0.
// We divide the problem as follows...
// Rows M are divided into:
// -	mb3 3 x n blocks,
// -	mtail 1 x n blocks.
// Within each row, the elements are divided into:
// -	nhead 1 x 1 blocks (used to align the rest, either 0 or 1),
// -	nb8 1 x 8 blocks, aligned to 16-byte boundary,
// -	nb2 1 x 2 blocks, aligned to 16-byte boundary,
// -	ntail 1 x 1 blocks, aligned too (altough we don't rely on it).
   ae_int_t n2 = n / 2;
   ae_int_t mb3 = m / 3;
   ae_int_t mtail = m % 3;
   ae_int_t nhead = ae_misalignment(a, alglib_simd_alignment) == 0 ? 0 : 1;
   ae_int_t nb8 = (n - nhead) / 8;
   ae_int_t nb2 = (n - nhead - 8 * nb8) / 2;
   ae_int_t ntail = n - nhead - 8 * nb8 - 2 * nb2;
   for (ae_int_t i = 0; i < mb3; i++) {
      double row0 = 0.0;
      double row1 = 0.0;
      double row2 = 0.0;
      const double *pb = x;
      const double *pa0 = a;
      const double *pa1 = a + alglib_r_block;
      const double *pa2 = a + alglib_twice_r_block;
      a += 3 * alglib_r_block;
      __m128d v0, v1, v2;
      if (nhead == 1) {
         __m128d vx = _mm_load_sd(pb);
         v0 = _mm_load_sd(pa0);
         v1 = _mm_load_sd(pa1);
         v2 = _mm_load_sd(pa2);
         v0 = _mm_mul_sd(v0, vx);
         v1 = _mm_mul_sd(v1, vx);
         v2 = _mm_mul_sd(v2, vx);
         pa0++;
         pa1++;
         pa2++;
         pb++;
      } else {
         v0 = _mm_setzero_pd();
         v1 = _mm_setzero_pd();
         v2 = _mm_setzero_pd();
      }
      for (ae_int_t k = 0; k < nb8; k++) {
#   if 1
      // A shuffle of the code for a simultaneous dot product.
      // See below for the commented-out unshuffled original version.
         __m128d vx = _mm_load_pd(pb);
         __m128d va0 = _mm_load_pd(pa0);
         __m128d va1 = _mm_load_pd(pa1);
         va0 = _mm_mul_pd(va0, vx);
         __m128d va2 = _mm_load_pd(pa2);
         v0 = _mm_add_pd(va0, v0);
         va1 = _mm_mul_pd(va1, vx);
         va0 = _mm_load_pd(pa0 + 2);
         v1 = _mm_add_pd(va1, v1);
         va2 = _mm_mul_pd(va2, vx);
         va1 = _mm_load_pd(pa1 + 2);
         v2 = _mm_add_pd(va2, v2);
         vx = _mm_load_pd(pb + 2);
         va0 = _mm_mul_pd(va0, vx);
         va2 = _mm_load_pd(pa2 + 2);
         v0 = _mm_add_pd(va0, v0);
         va1 = _mm_mul_pd(va1, vx);
         va0 = _mm_load_pd(pa0 + 4);
         v1 = _mm_add_pd(va1, v1);
         va2 = _mm_mul_pd(va2, vx);
         va1 = _mm_load_pd(pa1 + 4);
         v2 = _mm_add_pd(va2, v2);
         vx = _mm_load_pd(pb + 4);
         va0 = _mm_mul_pd(va0, vx);
         va2 = _mm_load_pd(pa2 + 4);
         v0 = _mm_add_pd(va0, v0);
         va1 = _mm_mul_pd(va1, vx);
         va0 = _mm_load_pd(pa0 + 6);
         v1 = _mm_add_pd(va1, v1);
         va2 = _mm_mul_pd(va2, vx);
         va1 = _mm_load_pd(pa1 + 6);
         v2 = _mm_add_pd(va2, v2);
         vx = _mm_load_pd(pb + 6);
         va0 = _mm_mul_pd(va0, vx);
         v0 = _mm_add_pd(va0, v0);
         va2 = _mm_load_pd(pa2 + 6);
         va1 = _mm_mul_pd(va1, vx);
         v1 = _mm_add_pd(va1, v1);
         va2 = _mm_mul_pd(va2, vx);
         v2 = _mm_add_pd(va2, v2);
         pa0 += 8;
         pa1 += 8;
         pa2 += 8;
         pb += 8;
#   else
      // The unshuffled version of the code above.
         __m128d vx = _mm_load_pd(pb);
         __m128d va0 = _mm_load_pd(pa0);
         __m128d va1 = _mm_load_pd(pa1);
         __m128d va2 = _mm_load_pd(pa2);
         va0 = _mm_mul_pd(va0, vx);
         va1 = _mm_mul_pd(va1, vx);
         va2 = _mm_mul_pd(va2, vx);
         v0 = _mm_add_pd(va0, v0);
         v1 = _mm_add_pd(va1, v1);
         v2 = _mm_add_pd(va2, v2);
         vx = _mm_load_pd(pb + 2);
         va0 = _mm_load_pd(pa0 + 2);
         va1 = _mm_load_pd(pa1 + 2);
         va2 = _mm_load_pd(pa2 + 2);
         va0 = _mm_mul_pd(va0, vx);
         va1 = _mm_mul_pd(va1, vx);
         va2 = _mm_mul_pd(va2, vx);
         v0 = _mm_add_pd(va0, v0);
         v1 = _mm_add_pd(va1, v1);
         v2 = _mm_add_pd(va2, v2);
         vx = _mm_load_pd(pb + 4);
         va0 = _mm_load_pd(pa0 + 4);
         va1 = _mm_load_pd(pa1 + 4);
         va2 = _mm_load_pd(pa2 + 4);
         va0 = _mm_mul_pd(va0, vx);
         va1 = _mm_mul_pd(va1, vx);
         va2 = _mm_mul_pd(va2, vx);
         v0 = _mm_add_pd(va0, v0);
         v1 = _mm_add_pd(va1, v1);
         v2 = _mm_add_pd(va2, v2);
         vx = _mm_load_pd(pb + 6);
         va0 = _mm_load_pd(pa0 + 6);
         va1 = _mm_load_pd(pa1 + 6);
         va2 = _mm_load_pd(pa2 + 6);
         va0 = _mm_mul_pd(va0, vx);
         va1 = _mm_mul_pd(va1, vx);
         va2 = _mm_mul_pd(va2, vx);
         v0 = _mm_add_pd(va0, v0);
         v1 = _mm_add_pd(va1, v1);
         v2 = _mm_add_pd(va2, v2);
#   endif
      }
      for (ae_int_t k = 0; k < nb2; k++) {
         __m128d vx = _mm_load_pd(pb);
         __m128d va0 = _mm_load_pd(pa0);
         __m128d va1 = _mm_load_pd(pa1);
         __m128d va2 = _mm_load_pd(pa2);
         va0 = _mm_mul_pd(va0, vx);
         v0 = _mm_add_pd(va0, v0);
         va1 = _mm_mul_pd(va1, vx);
         v1 = _mm_add_pd(va1, v1);
         va2 = _mm_mul_pd(va2, vx);
         v2 = _mm_add_pd(va2, v2);
         pa0 += 2;
         pa1 += 2;
         pa2 += 2;
         pb += 2;
      }
      for (ae_int_t k = 0; k < ntail; k++) {
         __m128d vx = _mm_load1_pd(pb);
         __m128d va0 = _mm_load1_pd(pa0);
         __m128d va1 = _mm_load1_pd(pa1);
         __m128d va2 = _mm_load1_pd(pa2);
         va0 = _mm_mul_sd(va0, vx);
         v0 = _mm_add_sd(v0, va0);
         va1 = _mm_mul_sd(va1, vx);
         v1 = _mm_add_sd(v1, va1);
         va2 = _mm_mul_sd(va2, vx);
         v2 = _mm_add_sd(v2, va2);
      }
      __m128d vtmp = _mm_add_pd(_mm_unpacklo_pd(v0, v1), _mm_unpackhi_pd(v0, v1));
      _mm_storel_pd(&row0, vtmp);
      _mm_storeh_pd(&row1, vtmp);
      v2 = _mm_add_sd(_mm_shuffle_pd(v2, v2, 1), v2);
      _mm_storel_pd(&row2, v2);
      if (Beta != 0.0) {
         y[0] = Beta * y[0] + Alpha * row0;
         y[stride] = Beta * y[stride] + Alpha * row1;
         y[2 * stride] = Beta * y[2 * stride] + Alpha * row2;
      } else {
         y[0] = Alpha * row0;
         y[stride] = Alpha * row1;
         y[2 * stride] = Alpha * row2;
      }
      y += 3 * stride;
   }
   for (ae_int_t i = 0; i < mtail; i++) {
      double row0 = 0.0;
      const double *pb = x;
      const double *pa0 = a;
      a += alglib_r_block;
      for (ae_int_t k = 0; k < n2; k++) {
         row0 += pb[0] * pa0[0] + pb[1] * pa0[1];
         pa0 += 2;
         pb += 2;
      }
      if (n % 2)
         row0 += pb[0] * pa0[0];
      if (Beta != 0.0)
         y[0] = Beta * y[0] + Alpha * row0;
      else
         y[0] = Alpha * row0;
      y += stride;
   }
}
#endif

// The m x n real matrix-vector product:
//	y = Beta y + Alpha a x
// with generic C code, optimized with _ialglib_mv_32() if (m, n) == (32, 32).
// If Beta == 0.0, y is overwritten by Alpha a x and if Alpha == 0.0, y is multiplied by Beta.
// However, the y = Beta y operation here is not done efficiently as a scalar-vector product, and should NOT be used this way.
// Important:
// *	(m, n) is in [0, alglib_r_block] x [0, alglib_r_block].
// *	a must be in row-major order with stride = alglib_r_block.
static void _ialglib_rmv(ae_int_t m, ae_int_t n, const double *a, const double *x, double *y, ae_int_t stride, double Alpha, double Beta) {
// Handle the special cases:
// *	Alpha == 0.0 or m n == 0
   if (m == 0)
      return;
   else if (Alpha == 0.0 || n == 0) {
      if (Beta == 0.0)
         for (ae_int_t i = 0; i < m; i++) {
            *y = 0.0;
            y += stride;
         }
      else
         for (ae_int_t i = 0; i < m; i++) {
            *y *= Beta;
            y += stride;
         }
      return;
   }
// Handle the general case: mn != 0 and Alpha != 0.0.
   if (m == 32 && n == 32)
   // The 32 x 32 case is specially optimized.
      _ialglib_mv_32(a, x, y, stride, Alpha, Beta);
   else {
   // Optimized by doing the first M/2 pairs of rows of a pairwise.
      ae_int_t m2 = m / 2;
      ae_int_t n8 = n / 8;
      ae_int_t ntrail2 = (n - 8 * n8) / 2;
      for (ae_int_t i = 0; i < m2; i++) {
         double v0 = 0.0, v1 = 0.0;
      // 'a' points to the part of the matrix which is not processed yet.
         const double *pb = x;
         const double *pa0 = a;
         const double *pa1 = a + alglib_r_block;
         a += alglib_twice_r_block;
      // 8 elements per iteration.
         for (ae_int_t k = 0; k < n8; k++) {
            v0 += pa0[0] * pb[0];
            v1 += pa1[0] * pb[0];
            v0 += pa0[1] * pb[1];
            v1 += pa1[1] * pb[1];
            v0 += pa0[2] * pb[2];
            v1 += pa1[2] * pb[2];
            v0 += pa0[3] * pb[3];
            v1 += pa1[3] * pb[3];
            v0 += pa0[4] * pb[4];
            v1 += pa1[4] * pb[4];
            v0 += pa0[5] * pb[5];
            v1 += pa1[5] * pb[5];
            v0 += pa0[6] * pb[6];
            v1 += pa1[6] * pb[6];
            v0 += pa0[7] * pb[7];
            v1 += pa1[7] * pb[7];
            pa0 += 8;
            pa1 += 8;
            pb += 8;
         }
      // 2 elements per iteration.
         for (ae_int_t k = 0; k < ntrail2; k++) {
            v0 += pa0[0] * pb[0];
            v1 += pa1[0] * pb[0];
            v0 += pa0[1] * pb[1];
            v1 += pa1[1] * pb[1];
            pa0 += 2;
            pa1 += 2;
            pb += 2;
         }
      // The last element, if needed.
         if (n % 2) {
            v0 += pa0[0] * pb[0];
            v1 += pa1[0] * pb[0];
         }
      // The final update.
         if (Beta != 0.0) {
            y[0] = Beta * y[0] + Alpha * v0;
            y[stride] = Beta * y[stride] + Alpha * v1;
         } else {
            y[0] = Alpha * v0;
            y[stride] = Alpha * v1;
         }
      // Move to the next pair of elements.
         y += 2 * stride;
      }
   // For odd m, the last row is done with the less optimized code.
      if (m % 2) {
         double v0 = 0.0;
      // 'a' points to the part of the matrix which is not processed yet.
         const double *pb = x;
         const double *pa0 = a;
      // 2 elements per iteration.
         ae_int_t n2 = n / 2;
         for (ae_int_t k = 0; k < n2; k++) {
            v0 += pa0[0] * pb[0] + pa0[1] * pb[1];
            pa0 += 2;
            pb += 2;
         }
      // The last element, if needed.
         if (n % 2)
            v0 += pa0[0] * pb[0];
      // The final update.
         if (Beta != 0.0)
            y[0] = Beta * y[0] + Alpha * v0;
         else
            y[0] = Alpha * v0;
      }
   }
}

#if defined AE_HAS_SSE2_INTRINSICS
// A fast m x n complex matrix-vector product (SSE2-optimized version)
//	y = Beta y + Alpha a x
// with generic C code for complex a, x, y, Alpha and Beta.
// If Beta == 0.0, y is overwritten by Alpha a x but if Alpha == 0.0, we still do y = Alpha a x (to properly handle infinities/NANs).
// Important:
// *	(m, n) is in [0, alglib_c_block] x [0, alglib_c_block].
// *	a must be in row-major order as (double, double) with stride = alglib_c_block double-pairs.
// *	y may point to either vector cy (Complex) or dy (double precision pairs).
// *	Either cy or dy must be NULL, but not both.
// *	Both a and x must be aligned, but y may be non-aligned.
// This function supports SSE2; it can be used when:
// *	AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
//	Otherwise, this function will be undefined.
// *	CurCPU contains CPU_SSE2 (checked at run-time)
//	Otherwise, calling this function will probably crash your system.
// If you want to know whether it is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
static void _ialglib_cmv_sse2(ae_int_t m, ae_int_t n, const double *a, const double *x, complex *cy, double *dy, ae_int_t stride, complex Alpha, complex Beta) {
   ae_int_t m2 = m / 2;
   const double *parow = a;
   if (cy != NULL) {
      dy = (double *)cy;
      cy = NULL;
   }
   __m128d vbeta = _mm_loadh_pd(_mm_load_sd(&Beta.x), &Beta.y);
   __m128d vbetax = _mm_unpacklo_pd(vbeta, vbeta);
   __m128d vbetay = _mm_unpackhi_pd(vbeta, vbeta);
   __m128d valpha = _mm_loadh_pd(_mm_load_sd(&Alpha.x), &Alpha.y);
   __m128d valphax = _mm_unpacklo_pd(valpha, valpha);
   __m128d valphay = _mm_unpackhi_pd(valpha, valpha);
   for (ae_int_t i = 0; i < m2; i++) {
      const double *pa0 = parow;
      const double *pa1 = parow + 2 * alglib_c_block;
      const double *pb = x;
      __m128d vx = _mm_setzero_pd();
      __m128d vy = _mm_setzero_pd();
      for (ae_int_t j = 0; j < n; j++) {
         __m128d vt0 = _mm_load1_pd(pb);
         __m128d vt1 = _mm_load1_pd(pb + 1);
         __m128d vt2 = _mm_load_pd(pa0);
         __m128d vt3 = _mm_load_pd(pa1);
         __m128d vt5 = _mm_unpacklo_pd(vt2, vt3);
         __m128d vt4 = _mm_unpackhi_pd(vt2, vt3);
         vt2 = vt5;
         vt3 = vt4;
         vt2 = _mm_mul_pd(vt2, vt0);
         vx = _mm_add_pd(vx, vt2);
         vt3 = _mm_mul_pd(vt3, vt1);
         vx = _mm_sub_pd(vx, vt3);
         vt4 = _mm_mul_pd(vt4, vt0);
         vy = _mm_add_pd(vy, vt4);
         vt5 = _mm_mul_pd(vt5, vt1);
         vy = _mm_add_pd(vy, vt5);
         pa0 += 2;
         pa1 += 2;
         pb += 2;
      }
      __m128d vrx, vry;
      if (Beta.x == 0.0 && Beta.y == 0.0) {
         vrx = _mm_setzero_pd();
         vry = _mm_setzero_pd();
      } else {
         __m128d vtx = _mm_loadh_pd(_mm_load_sd(dy), dy + 2 * stride);
         __m128d vty = _mm_loadh_pd(_mm_load_sd(dy + 1), dy + 2 * stride + 1);
         vrx = _mm_sub_pd(_mm_mul_pd(vbetax, vtx), _mm_mul_pd(vbetay, vty));
         vry = _mm_add_pd(_mm_mul_pd(vbetax, vty), _mm_mul_pd(vbetay, vtx));
      }
      __m128d vtx = _mm_sub_pd(_mm_mul_pd(valphax, vx), _mm_mul_pd(valphay, vy));
      __m128d vty = _mm_add_pd(_mm_mul_pd(valphax, vy), _mm_mul_pd(valphay, vx));
      vrx = _mm_add_pd(vrx, vtx);
      vry = _mm_add_pd(vry, vty);
      _mm_storel_pd(dy, vrx);
      _mm_storeh_pd(dy + 2 * stride, vrx);
      _mm_storel_pd(dy + 1, vry);
      _mm_storeh_pd(dy + 2 * stride + 1, vry);
      dy += 4 * stride;
      parow += 4 * alglib_c_block;
   }
   if (m % 2) {
      double v0 = 0.0, v1 = 0.0;
      const double *pa0 = parow;
      const double *pb = x;
      for (ae_int_t j = 0; j < n; j++) {
         v0 += pa0[0] * pb[0];
         v1 += pa0[0] * pb[1];
         v0 -= pa0[1] * pb[1];
         v1 += pa0[1] * pb[0];
         pa0 += 2;
         pb += 2;
      }
      double tx, ty;
      if (Beta.x == 0.0 && Beta.y == 0.0) {
         tx = 0.0;
         ty = 0.0;
      } else {
         tx = Beta.x * dy[0] - Beta.y * dy[1];
         ty = Beta.x * dy[1] + Beta.y * dy[0];
      }
      tx += Alpha.x * v0 - Alpha.y * v1;
      ty += Alpha.x * v1 + Alpha.y * v0;
      dy[0] = tx;
      dy[1] = ty;
      dy += 2 * stride;
      parow += 2 * alglib_c_block;
   }
}
#endif

// A fast m x n complex matrix-vector product
//	y = Beta y + Alpha a x
// with generic C code for complex a, x, y, Alpha and Beta.
// If Beta == 0.0, y is overwritten by Alpha a x but if Alpha == 0.0, we still do y = Alpha a x (to properly handle infinities/NANs).
// Important:
// *	(m, n) is in [0, alglib_c_block] x [0, alglib_c_block].
// *	a must be in row-major order as (double, double) with stride = alglib_c_block double-pairs.
// *	y may point to either vector cy (Complex) or dy (double precision pairs).
// *	Either cy or dy must be NULL, but not both.
// *	Both a and x must be aligned, but y may be non-aligned.
static void _ialglib_cmv(ae_int_t m, ae_int_t n, const double *a, const double *x, complex *cy, double *dy, ae_int_t stride, complex Alpha, complex Beta) {
   const double *parow = a;
   for (ae_int_t i = 0; i < m; i++) {
      double v0 = 0.0, v1 = 0.0;
      const double *pa = parow;
      const double *pb = x;
      for (ae_int_t j = 0; j < n; j++) {
         v0 += pa[0] * pb[0];
         v1 += pa[0] * pb[1];
         v0 -= pa[1] * pb[1];
         v1 += pa[1] * pb[0];
         pa += 2;
         pb += 2;
      }
      if (cy != NULL) {
         double tx = (Beta.x * cy->x - Beta.y * cy->y) + (Alpha.x * v0 - Alpha.y * v1);
         double ty = (Beta.x * cy->y + Beta.y * cy->x) + (Alpha.x * v1 + Alpha.y * v0);
         cy->x = tx;
         cy->y = ty;
         cy += stride;
      } else {
         double tx = (Beta.x * dy[0] - Beta.y * dy[1]) + (Alpha.x * v0 - Alpha.y * v1);
         double ty = (Beta.x * dy[1] + Beta.y * dy[0]) + (Alpha.x * v1 + Alpha.y * v0);
         dy[0] = tx;
         dy[1] = ty;
         dy += 2 * stride;
      }
      parow += 2 * alglib_c_block;
   }
}

// Set the real n-vector p to zero.
static void _ialglib_vzero(ae_int_t n, double *p, ae_int_t stride) {
   if (stride == 1)
      for (ae_int_t i = 0; i < n; i++, p++)
         *p = 0.0;
   else
      for (ae_int_t i = 0; i < n; i++, p += stride)
         *p = 0.0;
}

// Set the complex n-vector p to zero.
static void _ialglib_vzero_complex(ae_int_t n, complex *p, ae_int_t stride) {
   if (stride == 1)
      for (ae_int_t i = 0; i < n; i++, p++)
         p->y = p->x = 0.0;
   else
      for (ae_int_t i = 0; i < n; i++, p += stride)
         p->y = p->x = 0.0;
}

// Copy the unaligned real vector a to b.
static void _ialglib_vcopy(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb) {
   if (stridea == 1 && strideb == 1)
      for (ae_int_t i = n; i != 0; i--, a++, b++)
         *b = *a;
   else
      for (ae_int_t i = 0; i < n; i++, a += stridea, b += strideb)
         *b = *a;
}

// Copy the unaligned complex vector a to b (passed as complex *).
// *	The strideb is in complex numbers rather than double-pairs.
// *	conj may be "N" (no conj.) or "C" (conj.)
static void _ialglib_vcopy_complex(ae_int_t n, const complex *a, ae_int_t stridea, double *b, ae_int_t strideb, const char *conj) {
   strideb *= 2; // Multiply by sizeof(complex)/sizeof(double).
   if (conj[0] == 'N' || conj[0] == 'n')
      for (ae_int_t i = 0; i < n; i++, a += stridea, b += strideb) {
         b[0] = a->x;
         b[1] = a->y;
      }
   else
      for (ae_int_t i = 0; i < n; i++, a += stridea, b += strideb) {
         b[0] = a->x;
         b[1] = -a->y;
      }
}

// Copy the unaligned complex vector a to b (passed as double *).
// *	The strideb is in complex numbers rather than double-pairs.
// *	conj may be "N" (no conj.) or "C" (conj.)
static void _ialglib_vcopy_dcomplex(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb, const char *conj) {
   stridea *= 2, strideb *= 2; // Multiply by sizeof(complex)/sizeof(double).
   if (conj[0] == 'N' || conj[0] == 'n')
      for (ae_int_t i = 0; i < n; i++, a += stridea, b += strideb) {
         b[0] = a[0];
         b[1] = a[1];
      }
   else
      for (ae_int_t i = 0; i < n; i++, a += stridea, b += strideb) {
         b[0] = a[0];
         b[1] = -a[1];
      }
}

#if defined AE_HAS_SSE2_INTRINSICS
// Copy the real non-aligned non-contigous m x n matrix a (SSE2-optimized version)
// to an m x n or n x m submatrix of the aligned contiguous alglib_r_block x alglib_r_block matrix b.
// Use transposition if op != 0.
// The stride is stride for a, alglib_r_block for b.
// This function supports SSE2; it can be used when:
// *	AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
//	Otherwise, this function will be undefined.
// *	CurCPU contains CPU_SSE2 (checked at run-time)
//	Otherwise, calling this function will probably crash your system.
// If you want to know whether it is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
static void _ialglib_mcopyblock_sse2(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_int_t stride, double *b) {
   if (op == 0) {
      ae_int_t nb8 = n / 8;
      ae_int_t ntail = n - 8 * nb8;
      const double *psrc0 = a;
      for (ae_int_t i = 0; i < m; i++, a += stride, b += alglib_r_block, psrc0 = a) {
         double *pdst = b;
         for (ae_int_t j = 0; j < nb8; j++) {
            __m128d v0 = _mm_loadu_pd(psrc0);
            _mm_store_pd(pdst, v0);
            __m128d v1 = _mm_loadu_pd(psrc0 + 2);
            _mm_store_pd(pdst + 2, v1);
            v1 = _mm_loadu_pd(psrc0 + 4);
            _mm_store_pd(pdst + 4, v1);
            v1 = _mm_loadu_pd(psrc0 + 6);
            _mm_store_pd(pdst + 6, v1);
            pdst += 8;
            psrc0 += 8;
         }
         for (ae_int_t j = 0; j < ntail; j++)
            pdst[j] = psrc0[j];
      }
   } else {
      ae_int_t n2 = n / 2;
      ae_int_t mb2 = m / 2;
      ae_int_t nb4 = n / 4;
      ae_int_t ntail = n - 4 * nb4;
      const double *arow0 = a;
      const double *arow1 = a + stride;
      double *bcol0 = b;
      double *bcol1 = b + 1;
      for (ae_int_t i = 0; i < mb2; i++) {
         const double *psrc0 = arow0;
         const double *psrc1 = arow1;
         double *pdst0 = bcol0;
         double *pdst1 = bcol1;
         for (ae_int_t j = 0; j < nb4; j++) {
            __m128d v0 = _mm_loadu_pd(psrc0);
            __m128d v1 = _mm_loadu_pd(psrc1);
            __m128d v2 = _mm_loadu_pd(psrc0 + 2);
            __m128d v3 = _mm_loadu_pd(psrc1 + 2);
            _mm_store_pd(pdst0, _mm_unpacklo_pd(v0, v1));
            _mm_store_pd(pdst0 + alglib_r_block, _mm_unpackhi_pd(v0, v1));
            _mm_store_pd(pdst0 + 2 * alglib_r_block, _mm_unpacklo_pd(v2, v3));
            _mm_store_pd(pdst0 + 3 * alglib_r_block, _mm_unpackhi_pd(v2, v3));
            pdst0 += 4 * alglib_r_block;
            pdst1 += 4 * alglib_r_block;
            psrc0 += 4;
            psrc1 += 4;
         }
         for (ae_int_t j = 0; j < ntail; j++) {
            pdst0[0] = psrc0[0];
            pdst1[0] = psrc1[0];
            pdst0 += alglib_r_block;
            pdst1 += alglib_r_block;
            psrc0++;
            psrc1++;
         }
         arow0 += 2 * stride;
         arow1 += 2 * stride;
         bcol0 += 2;
         bcol1 += 2;
      }
      if (m % 2) {
         const double *psrc0 = arow0;
         double *pdst0 = bcol0;
         for (ae_int_t j = 0; j < n2; j++) {
            pdst0[0] = psrc0[0];
            pdst0[alglib_r_block] = psrc0[1];
            pdst0 += alglib_twice_r_block;
            psrc0 += 2;
         }
         if (n % 2)
            pdst0[0] = psrc0[0];
      }
   }
}
#endif

// Copy the real non-aligned non-contigous m x n matrix a
// to an m x n or n x m submatrix of the aligned contiguous alglib_r_block x alglib_r_block matrix b.
// Use transposition if op != 0.
// The stride is stride for a, alglib_r_block for b.
static void _ialglib_mcopyblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_int_t stride, double *b) {
   ae_int_t n2 = n / 2;
   if (op == 0)
      for (ae_int_t i = 0; i < m; i++, a += stride, b += alglib_r_block) {
         const double *psrc = a;
         double *pdst = b;
         for (ae_int_t j = 0; j < n2; j++, pdst += 2, psrc += 2) {
            pdst[0] = psrc[0];
            pdst[1] = psrc[1];
         }
         if (n % 2)
            pdst[0] = psrc[0];
      }
   else
      for (ae_int_t i = 0; i < m; i++, a += stride, b++) {
         const double *psrc = a;
         double *pdst = b;
         for (ae_int_t j = 0; j < n2; j++, pdst += alglib_twice_r_block, psrc += 2) {
            pdst[0] = psrc[0];
            pdst[alglib_r_block] = psrc[1];
         }
         if (n % 2)
            pdst[0] = psrc[0];
      }
}

// Copy an m x n submatrix of the real non-aligned non-contigous stride x stride matrix a to the aligned contigous m x n or n x m matrix b.
// Use transposition if op != 0.
// The stride is alglib_r_block for a, stride for b.
static void _ialglib_mcopyunblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, double *b, ae_int_t stride) {
   ae_int_t n2 = n / 2;
   if (op == 0)
      for (ae_int_t i = 0; i < m; i++, a += alglib_r_block, b += stride) {
         const double *psrc = a;
         double *pdst = b;
         for (ae_int_t j = 0; j < n2; j++, pdst += 2, psrc += 2) {
            pdst[0] = psrc[0];
            pdst[1] = psrc[1];
         }
         if (n % 2)
            pdst[0] = psrc[0];
      }
   else
      for (ae_int_t i = 0; i < m; i++, a++, b += stride) {
         const double *psrc = a;
         double *pdst = b;
         for (ae_int_t j = 0; j < n2; j++, pdst += 2, psrc += alglib_twice_r_block) {
            pdst[0] = psrc[0];
            pdst[1] = psrc[alglib_r_block];
         }
         if (n % 2)
            pdst[0] = psrc[0];
      }
}

// Copy the complex aligned contigous m x n matrix a to an m x n or n x m submatrix of the non-aligned non-contigous 2 x stride x stride matrix b.
// Copy as is or use transposition, Hermitian conjugate or conjugate, respectively if Op == 0, 1, 2 or 3.
// The stride is stride for a, alglib_c_block double-pairs for b.
static void _ialglib_mcopyblock_complex(ae_int_t m, ae_int_t n, const complex *a, ae_int_t op, ae_int_t stride, double *b) {
   switch (op) {
      case 0:
         for (ae_int_t i = 0; i < m; i++, a += stride, b += alglib_twice_c_block) {
            const complex *psrc = a;
            double *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst += 2, psrc++) {
               pdst[0] = psrc->x;
               pdst[1] = psrc->y;
            }
         }
      break;
      case 1:
         for (ae_int_t i = 0; i < m; i++, a += stride, b += 2) {
            const complex *psrc = a;
            double *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst += alglib_twice_c_block, psrc++) {
               pdst[0] = psrc->x;
               pdst[1] = psrc->y;
            }
         }
      break;
      case 2:
         for (ae_int_t i = 0; i < m; i++, a += stride, b += 2) {
            const complex *psrc = a;
            double *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst += alglib_twice_c_block, psrc++) {
               pdst[0] = psrc->x;
               pdst[1] = -psrc->y;
            }
         }
      break;
      case 3:
         for (ae_int_t i = 0; i < m; i++, a += stride, b += alglib_twice_c_block) {
            const complex *psrc = a;
            double *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst += 2, psrc++) {
               pdst[0] = psrc->x;
               pdst[1] = -psrc->y;
            }
         }
      break;
   }
}

// Copy the complex aligned contigous m x n matrix a to an m x n or n x m submatrix of the non-aligned non-contigous 2 x stride x stride matrix b.
// Copy as is or use transposition, Hermitian conjugate or conjugate, respectively if op == 0, 1, 2 or 3.
// The stride is stride for a, alglib_c_block double-pairs for b.
static void _ialglib_mcopyunblock_complex(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, complex *b, ae_int_t stride) {
   switch (op) {
      case 0:
         for (ae_int_t i = 0; i < m; i++, a += alglib_twice_c_block, b += stride) {
            const double *psrc = a;
            complex *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst++, psrc += 2) {
               pdst->x = psrc[0];
               pdst->y = psrc[1];
            }
         }
      break;
      case 1:
         for (ae_int_t i = 0; i < m; i++, a += 2, b += stride) {
            const double *psrc = a;
            complex *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst++, psrc += alglib_twice_c_block) {
               pdst->x = psrc[0];
               pdst->y = psrc[1];
            }
         }
      break;
      case 2:
         for (ae_int_t i = 0; i < m; i++, a += 2, b += stride) {
            const double *psrc = a;
            complex *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst++, psrc += alglib_twice_c_block) {
               pdst->x = psrc[0];
               pdst->y = -psrc[1];
            }
         }
      break;
      case 3:
         for (ae_int_t i = 0; i < m; i++, a += alglib_twice_c_block, b += stride) {
            const double *psrc = a;
            complex *pdst = b;
            for (ae_int_t j = 0; j < n; j++, pdst++, psrc += 2) {
               pdst->x = psrc[0];
               pdst->y = -psrc[1];
            }
         }
      break;
   }
}

// Real GEMM kernel.
static bool _ialglib_rmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, double Alpha, double *_a, ae_int_t _a_stride, ae_int_t opa, double *_b, ae_int_t _b_stride, ae_int_t opb, double Beta, double *_c, ae_int_t _c_stride) {
   if (m > alglib_r_block || n > alglib_r_block || k > alglib_r_block || m <= 0 || n <= 0 || k <= 0 || Alpha == 0.0)
      return false;
   void (*rmv)(ae_int_t, ae_int_t, const double *, const double *, double *, ae_int_t, double, double) = &_ialglib_rmv;
   void (*mcopyblock)(ae_int_t, ae_int_t, const double *, ae_int_t, ae_int_t, double *) = &_ialglib_mcopyblock;
#ifdef AE_HAS_SSE2_INTRINSICS
// Check for SSE2 support.
   if (CurCPU & CPU_SSE2) {
      rmv = &_ialglib_rmv_sse2;
      mcopyblock = &_ialglib_mcopyblock_sse2;
   }
#endif
// Copy or transpose b.
   double _bbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const b = (double *)ae_align(_bbuf, alglib_simd_alignment);
   if (opb == 0)
      mcopyblock(k, n, _b, 1, _b_stride, b);
   else
      mcopyblock(n, k, _b, 0, _b_stride, b);
// Multiply b by a (from the right, by rows) into c.
   double _abuf[alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_abuf, alglib_simd_alignment);
   bool BetaIs0 = Beta == 0.0;
   double *crow = _c;
   if (opa == 0) {
      const double *arow = _a;
      for (ae_int_t i = 0; i < m; i++) {
         _ialglib_vcopy(k, arow, 1, abuf, 1);
         if (BetaIs0)
            _ialglib_vzero(n, crow, 1);
         rmv(n, k, b, abuf, crow, 1, Alpha, Beta);
         crow += _c_stride;
         arow += _a_stride;
      }
   } else {
      const double *acol = _a;
      for (ae_int_t i = 0; i < m; i++) {
         _ialglib_vcopy(k, acol, _a_stride, abuf, 1);
         if (BetaIs0)
            _ialglib_vzero(n, crow, 1);
         rmv(n, k, b, abuf, crow, 1, Alpha, Beta);
         crow += _c_stride;
         acol++;
      }
   }
   return true;
}

// Complex GEMM kernel.
static bool _ialglib_cmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, complex Alpha, complex *_a, ae_int_t _a_stride, ae_int_t opa, complex *_b, ae_int_t _b_stride, ae_int_t opb, complex Beta, complex *_c, ae_int_t _c_stride) {
   if (m > alglib_c_block || n > alglib_c_block || k > alglib_c_block)
      return false;
   void (*cmv)(ae_int_t, ae_int_t, const double *, const double *, complex *, double *, ae_int_t, complex, complex) = &_ialglib_cmv;
#ifdef AE_HAS_SSE2_INTRINSICS
// Check for SSE2 support.
   if (CurCPU & CPU_SSE2) {
      cmv = &_ialglib_cmv_sse2;
   }
#endif
// Copy, transpose, conjugate or Hermitianize b.
   double _loc_b[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const b = (double *)ae_align(_loc_b, alglib_simd_alignment);
   ae_int_t bcols = opb == 0 ? n : k;
   ae_int_t brows = opb == 0 ? k : n;
   switch (opb) {
      default:
      case 0:
         _ialglib_mcopyblock_complex(brows, bcols, _b, 1, _b_stride, b);
      break;
      case 1:
         _ialglib_mcopyblock_complex(brows, bcols, _b, 0, _b_stride, b);
      break;
      case 2:
         _ialglib_mcopyblock_complex(brows, bcols, _b, 3, _b_stride, b);
      break;
   }
// Multiply b by a (from the right, by rows) into c.
   double _loc_abuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   bool BetaIs0 = Beta.x == 0.0 && Beta.y == 0.0;
   const complex *arow = _a;
   complex *crow = _c;
   for (ae_int_t i = 0; i < m; i++) {
      switch (opa) {
         case 0:
            _ialglib_vcopy_complex(k, arow, 1, abuf, 1, "No conj");
            arow += _a_stride;
         break;
         case 1:
            _ialglib_vcopy_complex(k, arow, _a_stride, abuf, 1, "No conj");
            arow++;
         break;
         default:
         case 2:
            _ialglib_vcopy_complex(k, arow, _a_stride, abuf, 1, "Conj");
            arow++;
         break;
      }
      if (BetaIs0)
         _ialglib_vzero_complex(n, crow, 1);
      cmv(n, k, b, abuf, crow, NULL, 1, Alpha, Beta);
      crow += _c_stride;
   }
   return true;
}

// Real Right TRSM kernel.
static bool _ialglib_rmatrixrighttrsm(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t opa, double *_x, ae_int_t _x_stride) {
   if (m > alglib_r_block || n > alglib_r_block)
      return false;
   void (*rmv)(ae_int_t, ae_int_t, const double *, const double *, double *, ae_int_t, double, double) = &_ialglib_rmv;
   void (*mcopyblock)(ae_int_t, ae_int_t, const double *, ae_int_t, ae_int_t, double *) = &_ialglib_mcopyblock;
#ifdef AE_HAS_SSE2_INTRINSICS
// Check for SSE2 support.
   if (CurCPU & CPU_SSE2) {
      rmv = &_ialglib_rmv_sse2;
      mcopyblock = &_ialglib_mcopyblock_sse2;
   }
#endif
// Prepare.
   double _loc_abuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   mcopyblock(n, n, _a, opa, _a_stride, abuf);
   double _loc_xbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   mcopyblock(m, n, _x, 0, _x_stride, xbuf);
   if (isunit) {
      double *pdiag = abuf;
      for (ae_int_t i = 0; i < n; i++, pdiag += alglib_r_block + 1)
         *pdiag = 1.0;
   }
// Solve y a^-1 == x, where a is upper or lower triangular.
   double _loc_tmpbuf[alglib_r_block + alglib_simd_alignment];
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   if ((opa == 0) == isupper) {
      double *pdiag = abuf;
      for (ae_int_t i = 0; i < n; i++, pdiag += alglib_r_block + 1) {
         double Beta = 1.0 / *pdiag;
         double Alpha = -Beta;
         _ialglib_vcopy(i, abuf + i, alglib_r_block, tmpbuf, 1);
         rmv(m, i, xbuf, tmpbuf, xbuf + i, alglib_r_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 0, _x, _x_stride);
   } else {
      double *pdiag = abuf + (n - 1) * alglib_r_block + (n - 1);
      for (ae_int_t i = n - 1; i >= 0; i--, pdiag -= alglib_r_block + 1) {
         double Beta = 1.0 / *pdiag;
         double Alpha = -Beta;
         _ialglib_vcopy(n - 1 - i, pdiag + alglib_r_block, alglib_r_block, tmpbuf + i + 1, 1);
         rmv(m, n - 1 - i, xbuf + i + 1, tmpbuf + i + 1, xbuf + i, alglib_r_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 0, _x, _x_stride);
   }
   return true;
}

// Complex Right TRSM kernel.
static bool _ialglib_cmatrixrighttrsm(ae_int_t m, ae_int_t n, complex *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t opa, complex *_x, ae_int_t _x_stride) {
   void (*cmv)(ae_int_t, ae_int_t, const double *, const double *, complex *, double *, ae_int_t, complex, complex) = &_ialglib_cmv;
   if (m > alglib_c_block || n > alglib_c_block)
      return false;
#ifdef AE_HAS_SSE2_INTRINSICS
// Check for SSE2 support.
   if (CurCPU & CPU_SSE2) {
      cmv = &_ialglib_cmv_sse2;
   }
#endif
// Prepare.
   double _loc_abuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   _ialglib_mcopyblock_complex(n, n, _a, opa, _a_stride, abuf);
   double _loc_xbuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   _ialglib_mcopyblock_complex(m, n, _x, 0, _x_stride, xbuf);
   if (isunit) {
      double *pdiag = abuf;
      for (ae_int_t i = 0; i < n; i++, pdiag += 2 * (alglib_c_block + 1)) {
         pdiag[0] = 1.0;
         pdiag[1] = 0.0;
      }
   }
// Solve y a^-1 == x, where a is upper or lower triangular.
   double _loc_tmpbuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   if ((opa == 0) == isupper) {
      double *pdiag = abuf;
      for (ae_int_t i = 0; i < n; i++, pdiag += 2 * (alglib_c_block + 1)) {
         complex tmp_c = complex_from_d(pdiag[0], pdiag[1]);
         complex Beta = ae_c_d_div(1.0, tmp_c);
         complex Alpha = complex_from_d(-Beta.x, -Beta.y);
         _ialglib_vcopy_dcomplex(i, abuf + 2 * i, alglib_c_block, tmpbuf, 1, "No conj");
         cmv(m, i, xbuf, tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 0, _x, _x_stride);
   } else {
      double *pdiag = abuf + 2 * (n - 1) * (alglib_c_block + 1);
      for (ae_int_t i = n - 1; i >= 0; i--, pdiag -= 2 * (alglib_c_block + 1)) {
         complex tmp_c = complex_from_d(pdiag[0], pdiag[1]);
         complex Beta = ae_c_d_div(1.0, tmp_c);
         complex Alpha = complex_from_d(-Beta.x, -Beta.y);
         _ialglib_vcopy_dcomplex(n - 1 - i, pdiag + 2 * alglib_c_block, alglib_c_block, tmpbuf, 1, "No conj");
         cmv(m, n - 1 - i, xbuf + 2 * (i + 1), tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 0, _x, _x_stride);
   }
   return true;
}

// Real Left TRSM kernel.
static bool _ialglib_rmatrixlefttrsm(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t opa, double *_x, ae_int_t _x_stride) {
   double _loc_tmpbuf[alglib_r_block + alglib_simd_alignment];
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   if (m > alglib_r_block || n > alglib_r_block)
      return false;
   void (*rmv)(ae_int_t, ae_int_t, const double *, const double *, double *, ae_int_t, double, double) = &_ialglib_rmv;
   void (*mcopyblock)(ae_int_t, ae_int_t, const double *, ae_int_t, ae_int_t, double *) = &_ialglib_mcopyblock;
#ifdef AE_HAS_SSE2_INTRINSICS
// Check for SSE2 support.
   if (CurCPU & CPU_SSE2) {
      rmv = &_ialglib_rmv_sse2;
      mcopyblock = &_ialglib_mcopyblock_sse2;
   }
#endif
// Prepare.
   double _loc_abuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   mcopyblock(m, m, _a, opa, _a_stride, abuf);
   double _loc_xbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   mcopyblock(m, n, _x, 1, _x_stride, xbuf);
// Transpose x (so we may use mv, which calculates a x, but not x a).
   if (isunit) {
      double *pdiag = abuf;
      for (ae_int_t i = 0; i < m; i++, pdiag += alglib_r_block + 1)
         *pdiag = 1.0;
   }
// Solve a^-1 y^T == x^T, where a is upper or lower triangular.
   if ((opa == 0) == isupper) {
      double *pdiag = abuf + (m - 1) * (alglib_r_block + 1);
      for (ae_int_t i = m - 1; i >= 0; i--, pdiag -= alglib_r_block + 1) {
         double Beta = 1.0 / *pdiag;
         double Alpha = -Beta;
         _ialglib_vcopy(m - 1 - i, pdiag + 1, 1, tmpbuf + i + 1, 1);
         rmv(n, m - 1 - i, xbuf + i + 1, tmpbuf + i + 1, xbuf + i, alglib_r_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 1, _x, _x_stride);
   } else {
      double *pdiag = abuf, *arow = abuf;
      for (ae_int_t i = 0; i < m; i++, pdiag += alglib_r_block + 1, arow += alglib_r_block) {
         double Beta = 1.0 / *pdiag;
         double Alpha = -Beta;
         _ialglib_vcopy(i, arow, 1, tmpbuf, 1);
         rmv(n, i, xbuf, tmpbuf, xbuf + i, alglib_r_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 1, _x, _x_stride);
   }
   return true;
}

// Complex Left TRSM kernel.
static bool _ialglib_cmatrixlefttrsm(ae_int_t m, ae_int_t n, complex *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t opa, complex *_x, ae_int_t _x_stride) {
   if (m > alglib_c_block || n > alglib_c_block)
      return false;
   void (*cmv)(ae_int_t, ae_int_t, const double *, const double *, complex *, double *, ae_int_t, complex, complex) = &_ialglib_cmv;
#ifdef AE_HAS_SSE2_INTRINSICS
// Check for SSE2 support.
   if (CurCPU & CPU_SSE2) {
      cmv = &_ialglib_cmv_sse2;
   }
#endif
// Prepare.
   double _loc_abuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   _ialglib_mcopyblock_complex(m, m, _a, opa, _a_stride, abuf);
   double _loc_xbuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   _ialglib_mcopyblock_complex(m, n, _x, 1, _x_stride, xbuf);
// Transpose x (so we may use mv, which calculates a x, but not x a).
   if (isunit) {
      double *pdiag = abuf;
      for (ae_int_t i = 0; i < m; i++, pdiag += 2 * (alglib_c_block + 1)) {
         pdiag[0] = 1.0;
         pdiag[1] = 0.0;
      }
   }
// Solve a^-1 y^T = x^T, where a is upper or lower triangular.
   double _loc_tmpbuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   if ((opa == 0) == isupper) {
      double *pdiag = abuf + 2 * (m - 1) * (alglib_c_block + 1);
      for (ae_int_t i = m - 1; i >= 0; i--, pdiag -= 2 * (alglib_c_block + 1)) {
         complex tmp_c = complex_from_d(pdiag[0], pdiag[1]);
         complex Beta = ae_c_d_div(1.0, tmp_c);
         complex Alpha = complex_from_d(-Beta.x, -Beta.y);
         _ialglib_vcopy_dcomplex(m - 1 - i, pdiag + 2, 1, tmpbuf, 1, "No conj");
         cmv(n, m - 1 - i, xbuf + 2 * (i + 1), tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 1, _x, _x_stride);
   } else {
      double *pdiag = abuf, *arow = abuf;
      for (ae_int_t i = 0; i < m; i++, pdiag += 2 * (alglib_c_block + 1), arow += 2 * alglib_c_block) {
         complex tmp_c = complex_from_d(pdiag[0], pdiag[1]);
         complex Beta = ae_c_d_div(1.0, tmp_c);
         complex Alpha = complex_from_d(-Beta.x, -Beta.y);
         _ialglib_vcopy_dcomplex(i, arow, 1, tmpbuf, 1, "No conj");
         cmv(n, i, xbuf, tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, Alpha, Beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 1, _x, _x_stride);
   }
   return true;
}

// Real SYRK/HERK kernel.
static bool _ialglib_rmatrixsyrk(ae_int_t n, ae_int_t k, double Alpha, double *_a, ae_int_t _a_stride, ae_int_t opa, double Beta, double *_c, ae_int_t _c_stride, bool isupper) {
   if (n > alglib_r_block || k > alglib_r_block)
      return false;
   else if (n == 0)
      return true;
// Copy a and c, and transformed the task to "a a^T"-form.
// If Beta == 0.0, then zero out c and ignore it.
// If Alpha == 0.0 or k == 0 then ignore a.
//(@) The SSE2-optimized form of _ialglib_mcopyblock() was not used here in the original ALGLIB version and hasn't yet been incorporated here.
   double _loc_abuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   if (Alpha == 0.0)
      k = 0;
   else if (k > 0)
      if (opa == 0)
         _ialglib_mcopyblock(n, k, _a, 0, _a_stride, abuf);
      else
         _ialglib_mcopyblock(k, n, _a, 1, _a_stride, abuf);
   double _loc_cbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const cbuf = (double *)ae_align(_loc_cbuf, alglib_simd_alignment);
   _ialglib_mcopyblock(n, n, _c, 0, _c_stride, cbuf);
   if (Beta == 0.0) {
      double *crow = cbuf;
      if (isupper)
         for (ae_int_t i = 0; i < n; i++, crow += alglib_r_block)
            _ialglib_vzero(n - i, crow + i, 1);
      else
         for (ae_int_t i = 0; i < n; i++, crow += alglib_r_block)
            _ialglib_vzero(i + 1, crow, 1);
   }
// Update c.
//(@) The SSE2-optimized form of _ialglib_rmv() was not used here in the original ALGLIB version and hasn't yet been incorporated here.
   double *arow = abuf, *crow = cbuf;
   if (isupper)
      for (ae_int_t i = 0; i < n; i++, arow += alglib_r_block, crow += alglib_r_block)
         _ialglib_rmv(n - i, k, arow, arow, crow + i, 1, Alpha, Beta);
   else
      for (ae_int_t i = 0; i < n; i++, arow += alglib_r_block, crow += alglib_r_block)
         _ialglib_rmv(i + 1, k, abuf, arow, crow, 1, Alpha, Beta);
// Copy back.
   _ialglib_mcopyunblock(n, n, cbuf, 0, _c, _c_stride);
   return true;
}

// Complex SYRK/HERK kernel.
static bool _ialglib_cmatrixherk(ae_int_t n, ae_int_t k, double Alpha, complex *_a, ae_int_t _a_stride, ae_int_t opa, double Beta, complex *_c, ae_int_t _c_stride, bool isupper) {
   if (n > alglib_c_block || k > alglib_c_block)
      return false;
   else if (n == 0)
      return true;
// Copy a and c, transforming to "a a^H"-form.
// If Beta == 0.0, then set c to 0 and ignore it.
// If Alpha == 0.0 or k == 0 then ignore a.
   double _loc_abuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   if (Alpha == 0.0)
      k = 0;
   else if (k > 0)
      if (opa == 0)
         _ialglib_mcopyblock_complex(n, k, _a, 3, _a_stride, abuf);
      else
         _ialglib_mcopyblock_complex(k, n, _a, 1, _a_stride, abuf);
   double _loc_cbuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const cbuf = (double *)ae_align(_loc_cbuf, alglib_simd_alignment);
   _ialglib_mcopyblock_complex(n, n, _c, 0, _c_stride, cbuf);
   if (Beta == 0.0) {
      double *crow = cbuf;
      for (ae_int_t i = 0; i < n; i++, crow += 2 * alglib_c_block)
         if (isupper)
            _ialglib_vzero(2 * (n - i), crow + 2 * i, 1);
         else
            _ialglib_vzero(2 * (i + 1), crow, 1);
   }
// Update c.
//(@) The SSE2-optimized form of _ialglib_cmv() was not used here in the original ALGLIB version and hasn't yet been incorporated here.
   double _loc_tmpbuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   double *arow = abuf, *crow = cbuf;
   complex c_alpha = complex_from_d(Alpha);
   complex c_beta = complex_from_d(Beta);
   if (isupper)
      for (ae_int_t i = 0; i < n; i++, arow += 2 * alglib_c_block, crow += 2 * alglib_c_block) {
         _ialglib_vcopy_dcomplex(k, arow, 1, tmpbuf, 1, "Conj");
         _ialglib_cmv(n - i, k, arow, tmpbuf, NULL, crow + 2 * i, 1, c_alpha, c_beta);
      }
   else
      for (ae_int_t i = 0; i < n; i++, arow += 2 * alglib_c_block, crow += 2 * alglib_c_block) {
         _ialglib_vcopy_dcomplex(k, arow, 1, tmpbuf, 1, "Conj");
         _ialglib_cmv(i + 1, k, abuf, tmpbuf, NULL, crow, 1, c_alpha, c_beta);
      }
// Copy back.
   _ialglib_mcopyunblock_complex(n, n, cbuf, 0, _c, _c_stride);
   return true;
}

// Real rank-1 kernel.
// Deprecated version.
static bool _ialglib_rmatrixrank1(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, double *_u, double *_v) {
// Quick exit.
   if (m <= 0 || n <= 0)
      return false;
// Locals.
   ae_int_t m2 = m / 2;
   ae_int_t n2 = n / 2;
   ae_int_t stride2 = 2 * _a_stride;
// Update the rows and columns by twos, handling the odd ones separately.
   double *arow0 = _a;
   double *arow1 = arow0 + _a_stride;
   double *pu = _u;
   double *vtmp = _v;
   for (ae_int_t i = 0; i < m2; i++, arow0 += stride2, arow1 += stride2, pu += 2) {
      double *pv = vtmp, *dst0 = arow0, *dst1 = arow1;
      for (ae_int_t j = 0; j < n2; j++, dst0 += 2, dst1 += 2, pv += 2) {
         dst0[0] += pu[0] * pv[0];
         dst0[1] += pu[0] * pv[1];
         dst1[0] += pu[1] * pv[0];
         dst1[1] += pu[1] * pv[1];
      }
      if (n % 2) {
         dst0[0] += pu[0] * pv[0];
         dst1[0] += pu[1] * pv[0];
      }
   }
   if (m % 2) {
      double *pv = vtmp, *dst0 = arow0;
      for (ae_int_t j = 0; j < n2; j++, dst0 += 2, pv += 2) {
         dst0[0] += pu[0] * pv[0];
         dst0[1] += pu[0] * pv[1];
      }
      if (n % 2)
         dst0[0] += pu[0] * pv[0];
   }
   return true;
}

// Complex rank-1 kernel.
static bool _ialglib_cmatrixrank1(ae_int_t m, ae_int_t n, complex *_a, ae_int_t _a_stride, complex *_u, complex *_v) {
// Quick exit.
   if (m <= 0 || n <= 0)
      return false;
// Locals.
   ae_int_t n2 = n / 2;
// Update the rows and columns by twos, handling the one ones separately.
   complex *arow = _a;
   complex *pu = _u;
   complex *vtmp = _v;
   for (ae_int_t i = 0; i < m; i++, arow += _a_stride, pu++) {
      complex *pv = vtmp, *dst = arow;
      for (ae_int_t j = 0; j < n2; j++, dst += 2, pv += 2) {
         double ux = pu[0].x;
         double uy = pu[0].y;
         double v0x = pv[0].x;
         double v0y = pv[0].y;
         double v1x = pv[1].x;
         double v1y = pv[1].y;
         dst[0].x += ux * v0x - uy * v0y;
         dst[0].y += ux * v0y + uy * v0x;
         dst[1].x += ux * v1x - uy * v1y;
         dst[1].y += ux * v1y + uy * v1x;
      }
      if (n % 2) {
         double ux = pu[0].x;
         double uy = pu[0].y;
         double vx = pv[0].x;
         double vy = pv[0].y;
         dst[0].x += ux * vx - uy * vy;
         dst[0].y += ux * vy + uy * vx;
      }
   }
   return true;
}

// Real rank-1 kernel.
// Deprecated version.
static bool _ialglib_rmatrixger(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, double Alpha, double *_u, double *_v) {
// Quick exit.
   if (m <= 0 || n <= 0 || Alpha == 0.0)
      return false;
// Locals.
   ae_int_t m2 = m / 2;
   ae_int_t n2 = n / 2;
   ae_int_t stride2 = 2 * _a_stride;
// Update the rows and columns by twos, handling the odd ones separately.
   double *arow0 = _a;
   double *arow1 = arow0 + _a_stride;
   double *pu = _u;
   double *vtmp = _v;
   for (ae_int_t i = 0; i < m2; i++, arow0 += stride2, arow1 += stride2, pu += 2) {
      double au0 = Alpha * pu[0];
      double au1 = Alpha * pu[1];
      double *pv = vtmp, *dst0 = arow0, *dst1 = arow1;
      for (ae_int_t j = 0; j < n2; j++, dst0 += 2, dst1 += 2, pv += 2) {
         dst0[0] += au0 * pv[0];
         dst0[1] += au0 * pv[1];
         dst1[0] += au1 * pv[0];
         dst1[1] += au1 * pv[1];
      }
      if (n % 2) {
         dst0[0] += au0 * pv[0];
         dst1[0] += au1 * pv[0];
      }
   }
   if (m % 2) {
      double au0 = Alpha * pu[0];
      double *pv = vtmp, *dst0 = arow0;
      for (ae_int_t j = 0; j < n2; j++, dst0 += 2, pv += 2) {
         dst0[0] += au0 * pv[0];
         dst0[1] += au0 * pv[1];
      }
      if (n % 2)
         dst0[0] += au0 * pv[0];
   }
   return true;
}

// Interface functions for efficient kernels.
// For KerGemmRF(), KerGemmCF(), KerTrsmRCF(), KerTrsmRRF(), KerTrsmLCF(), KerTrsmLRF(), KerHerkCF(), KerSyrkRF():
// For:
//	_ialglib_i_rmatrixgemmf, _ialglib_i_rmatrixrighttrsmf, _ialglib_i_rmatrixlefttrsmf, _ialglib_i_rmatrixsyrkf,
//	_ialglib_i_cmatrixgemmf, _ialglib_i_cmatrixrighttrsmf, _ialglib_i_cmatrixlefttrsmf, _ialglib_i_cmatrixherkf:
// *	handle the degenerate cases like zero matrices by ALGLIB++ - greatly simplifies passing data to ALGLIB++ kernel,
// *	handle the general cases with the optimized ALGLIB++ kernel.

bool _ialglib_i_rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double Alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double Beta, ae_matrix *c, ae_int_t ic, ae_int_t jc) {
   if (Alpha == 0.0 || k == 0 || n == 0 || m == 0)
      return false;
   return _ialglib_rmatrixgemm(m, n, k, Alpha, a->xyR[ia] + ja, a->stride, opa, b->xyR[ib] + jb, b->stride, opb, Beta, c->xyR[ic] + jc, c->stride);
}

bool _ialglib_i_cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, complex Alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, complex Beta, ae_matrix *c, ae_int_t ic, ae_int_t jc) {
   if (Alpha.x == 0.0 && Alpha.y == 0.0 || k == 0 || n == 0 || m == 0)
      return false;
   return _ialglib_cmatrixgemm(m, n, k, Alpha, a->xyC[ia] + ja, a->stride, opa, b->xyC[ib] + jb, b->stride, opb, Beta, c->xyC[ic] + jc, c->stride);
}

bool _ialglib_i_rmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
   if (m == 0 || n == 0)
      return false;
   return _ialglib_rmatrixrighttrsm(m, n, &a->xyR[i1][j1], a->stride, isupper, isunit, opa, &x->xyR[i2][j2], x->stride);
}

bool _ialglib_i_cmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
   if (m == 0 || n == 0)
      return false;
   return _ialglib_cmatrixrighttrsm(m, n, &a->xyC[i1][j1], a->stride, isupper, isunit, opa, &x->xyC[i2][j2], x->stride);
}

bool _ialglib_i_rmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
   if (m == 0 || n == 0)
      return false;
   return _ialglib_rmatrixlefttrsm(m, n, &a->xyR[i1][j1], a->stride, isupper, isunit, opa, &x->xyR[i2][j2], x->stride);
}

bool _ialglib_i_cmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t opa, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
   if (m == 0 || n == 0)
      return false;
   return _ialglib_cmatrixlefttrsm(m, n, &a->xyC[i1][j1], a->stride, isupper, isunit, opa, &x->xyC[i2][j2], x->stride);
}

bool _ialglib_i_rmatrixsyrkf(ae_int_t n, ae_int_t k, double Alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, double Beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   if (Alpha == 0.0 || k == 0 || n == 0)
      return false;
   return _ialglib_rmatrixsyrk(n, k, Alpha, &a->xyR[ia][ja], a->stride, opa, Beta, &c->xyR[ic][jc], c->stride, isupper);
}

bool _ialglib_i_cmatrixherkf(ae_int_t n, ae_int_t k, double Alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, double Beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
   if (Alpha == 0.0 || k == 0 || n == 0)
      return false;
   return _ialglib_cmatrixherk(n, k, Alpha, &a->xyC[ia][ja], a->stride, opa, Beta, &c->xyC[ic][jc], c->stride, isupper);
}

bool _ialglib_i_rmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs) {
   return _ialglib_rmatrixrank1(m, n, &a->xyR[ia][ja], a->stride, &u->xR[uoffs], &v->xR[voffs]);
}

bool _ialglib_i_cmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs) {
   return _ialglib_cmatrixrank1(m, n, &a->xyC[ia][ja], a->stride, &u->xC[uoffs], &v->xC[voffs]);
}

bool _ialglib_i_rmatrixgerf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, double Alpha, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs) {
   return _ialglib_rmatrixger(m, n, &a->xyR[ia][ja], a->stride, Alpha, &u->xR[uoffs], &v->xR[voffs]);
}

#if defined AE_HAS_SSE2_INTRINSICS
// Pack the matrix (col0, col1) with stride src_stride and length n row-by-row into the contiguous matrix dst (SSE2-optimized version).
// dst must be aligned, col0 and col1 may be non-aligned.
// It can handle following special cases:
// *	col1 == NULL:				the second column is treated as 0,
// *	src_stride == 1 or col1 - col0 == 1:	efficient SSE-based code is used.
// This function supports SSE2; it can be used when:
// *	AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// *	CurCPU contains CPU_SSE2 (checked at run-time)
// If you want to know whether it is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
static void _ialglib_pack_n2_sse2(double *col0, double *col1, ae_int_t n, ae_int_t src_stride, double *dst) {
   if (col1 == NULL) // Handle the special case: col1 == NULL.
      for (ae_int_t j = 0; j < n; j++) {
         dst[0] = *col0;
         dst[1] = 0.0;
         col0 += src_stride;
         dst += 2;
      }
   else if (src_stride == 1) { // Handle the case for unit stride.
      ae_int_t n2 = n / 2;
      for (ae_int_t j = 0; j < n2; j++) {
         __m128d v0 = _mm_loadu_pd(col0);
         col0 += 2;
         __m128d v1 = _mm_loadu_pd(col1);
         col1 += 2;
         _mm_store_pd(dst, _mm_unpacklo_pd(v0, v1));
         _mm_store_pd(dst + 2, _mm_unpackhi_pd(v0, v1));
         dst += 4;
      }
      if (n % 2) {
         dst[0] = *col0;
         dst[1] = *col1;
      }
   } else if (col1 - col0 == 1) { // Handle the case: col1 - col0 == 1.
      ae_int_t n2 = n / 2;
      ae_int_t stride2 = 2 * src_stride;
      for (ae_int_t j = 0; j < n2; j++) {
         __m128d v0 = _mm_loadu_pd(col0);
         __m128d v1 = _mm_loadu_pd(col0 + src_stride);
         _mm_store_pd(dst, v0);
         _mm_store_pd(dst + 2, v1);
         col0 += stride2;
         dst += 4;
      }
      if (n % 2) {
         dst[0] = col0[0];
         dst[1] = col0[1];
      }
   } else { // The general case.
      ae_int_t n2 = n / 2;
      ae_int_t stride2 = 2 * src_stride;
      for (ae_int_t j = 0; j < n2; j++) {
         dst[0] = *col0;
         dst[1] = *col1;
         dst[2] = col0[src_stride];
         dst[3] = col1[src_stride];
         col0 += stride2;
         col1 += stride2;
         dst += 4;
      }
      if (n % 2) {
         dst[0] = *col0;
         dst[1] = *col1;
      }
   }
}
#endif

// Pack the matrix (col0, col1) with stride src_stride and length n row-by-row into the contiguous matrix dst.
// dst must be aligned, col0 and col1 may be non-aligned.
// It can handle following special case:
// *	col1 == NULL:				the second column is treated as 0,
void _ialglib_pack_n2(double *col0, double *col1, ae_int_t n, ae_int_t src_stride, double *dst) {
   if (col1 == NULL) // Handle the special case: col1 == NULL.
      for (ae_int_t j = 0; j < n; j++) {
         dst[0] = *col0;
         dst[1] = 0.0;
         col0 += src_stride;
         dst += 2;
      }
   else { // Handle the general case.
      ae_int_t n2 = n / 2;
      ae_int_t stride2 = src_stride * 2;
      for (ae_int_t j = 0; j < n2; j++) {
         dst[0] = *col0;
         dst[1] = *col1;
         dst[2] = col0[src_stride];
         dst[3] = col1[src_stride];
         col0 += stride2;
         col1 += stride2;
         dst += 4;
      }
      if (n % 2) {
         dst[0] = *col0;
         dst[1] = *col1;
      }
   }
}

#if defined AE_HAS_SSE2_INTRINSICS
// Calculate (SSE2-optimized)
//	r = Alpha a' b + Beta r
// where
// *	a and b are aligned k x 2 matrices stored in contiguous row-by-row storage,
// *	r is a 2 x 2 matrix stored in non-contiguous row-by-row storage and may be non-aligned.
// If Beta == 0.0, then r is ignored, rather than multiplied by zero.
// but if Alpha == 0.0, we still calculate a' b even though it is multiplied by zero afterwards.
// The format of the result for r, based on store_mode, is as follows:
// *	store_mode == 0:	the full r is stored
// *	store_mode == 1:	only the first row of r is stored
// *	store_mode == 2:	only the first column of r is stored
// *	store_mode == 3:	only the top left element of r is stored
// This function supports SSE2; it can be used when:
// *	AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
//	Otherwise, this function will be undefined.
// *	CurCPU contains CPU_SSE2 (checked at run-time)
//	Otherwise, calling this function will probably crash your system.
// If you want to know whether it is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
static void _ialglib_mm22_sse2(double Alpha, const double *a, const double *b, ae_int_t k, double Beta, double *r, ae_int_t stride, ae_int_t store_mode) {
// Calculate the 2 x 2 matrix product the k x 2 matrices va and vb and store the result as follows:
//	a' b	=	/ vd[0] ve[0] \
//			\ ve[1] vd[1] /
   ae_int_t k2 = k / 2;
   __m128d vd = _mm_setzero_pd();
   __m128d ve = _mm_setzero_pd();
   for (ae_int_t t = 0; t < k2; t++) {
      __m128d vb = _mm_load_pd(b);
      __m128d va = _mm_load_pd(a);
      __m128d vt = vb;
      vb = _mm_mul_pd(va, vb);
      vt = _mm_shuffle_pd(vt, vt, 1);
      vd = _mm_add_pd(vb, vd);
      vt = _mm_mul_pd(va, vt);
      vb = _mm_load_pd(b + 2);
      ve = _mm_add_pd(vt, ve);
      va = _mm_load_pd(a + 2);
      vt = vb;
      vb = _mm_mul_pd(va, vb);
      vt = _mm_shuffle_pd(vt, vt, 1);
      vd = _mm_add_pd(vb, vd);
      vt = _mm_mul_pd(va, vt);
      ve = _mm_add_pd(vt, ve);
      a += 4;
      b += 4;
   }
   if (k % 2) {
      __m128d va = _mm_load_pd(a);
      __m128d vb = _mm_load_pd(b);
      __m128d vt = _mm_shuffle_pd(vb, vb, 1);
      vd = _mm_add_pd(_mm_mul_pd(va, vb), vd);
      ve = _mm_add_pd(_mm_mul_pd(va, vt), ve);
   }
// r0 is the first row of Alpha a' b, r1 is the second row.
   __m128d valpha = _mm_load1_pd(&Alpha);
   __m128d r0 = _mm_mul_pd(_mm_unpacklo_pd(vd, ve), valpha);
   __m128d r1 = _mm_mul_pd(_mm_unpackhi_pd(ve, vd), valpha);
// Format and store the result r.
   switch (store_mode) {
      case 0:
         if (Beta == 0.0) {
            _mm_storeu_pd(r, r0);
            _mm_storeu_pd(r + stride, r1);
         } else {
            __m128d vbeta = _mm_load1_pd(&Beta);
            _mm_storeu_pd(r, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r), vbeta), r0));
            _mm_storeu_pd(r + stride, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + stride), vbeta), r1));
         }
      break;
      case 1:
         if (Beta == 0.0)
            _mm_storeu_pd(r, r0);
         else
            _mm_storeu_pd(r, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r), _mm_load1_pd(&Beta)), r0));
      break;
      case 2: {
         double buf[4];
         _mm_storeu_pd(buf, r0);
         _mm_storeu_pd(buf + 2, r1);
         if (Beta == 0.0) {
            r[0] = buf[0];
            r[stride] = buf[2];
         } else {
            r[0] = Beta * r[0] + buf[0];
            r[stride] = Beta * r[stride] + buf[2];
         }
      }
      break;
      case 3: {
         double buf[2];
         _mm_storeu_pd(buf, r0);
         if (Beta == 0.0)
            r[0] = buf[0];
         else
            r[0] = Beta * r[0] + buf[0];
      }
      break;
   }
}
#endif

// Calculate
//	r = Alpha a' b + Beta r
// where
// *	a and b are aligned k x 2 matrices stored in contiguous row-by-row storage,
// *	r is a 2 x 2 matrix stored in non-contiguous row-by-row storage and may be non-aligned.
// If Beta == 0.0, then r is ignored, rather than multiplied by zero.
// but if Alpha == 0.0, we still calculate a' b even though it is multiplied by zero afterwards.
// The format of the result for r, based on store_mode, is as follows:
// *	store_mode == 0:	the full r is stored
// *	store_mode == 1:	only the first row of r is stored
// *	store_mode == 2:	only the first column of r is stored
// *	store_mode == 3:	only the top left element of r is stored
void _ialglib_mm22(double Alpha, const double *a, const double *b, ae_int_t k, double Beta, double *r, ae_int_t stride, ae_int_t store_mode) {
   double v00 = 0.0;
   double v01 = 0.0;
   double v10 = 0.0;
   double v11 = 0.0;
   for (ae_int_t t = 0; t < k; t++) {
      v00 += a[0] * b[0];
      v01 += a[0] * b[1];
      v10 += a[1] * b[0];
      v11 += a[1] * b[1];
      a += 2;
      b += 2;
   }
   switch (store_mode) {
      case 0:
         if (Beta == 0.0) {
            r[0] = Alpha * v00;
            r[1] = Alpha * v01;
            r[stride] = Alpha * v10;
            r[stride + 1] = Alpha * v11;
         } else {
            r[0] = Beta * r[0] + Alpha * v00;
            r[1] = Beta * r[1] + Alpha * v01;
            r[stride] = Beta * r[stride] + Alpha * v10;
            r[stride + 1] = Beta * r[stride + 1] + Alpha * v11;
         }
      break;
      case 1:
         if (Beta == 0.0) {
            r[0] = Alpha * v00;
            r[1] = Alpha * v01;
         } else {
            r[0] = Beta * r[0] + Alpha * v00;
            r[1] = Beta * r[1] + Alpha * v01;
         }
      break;
      case 2:
         if (Beta == 0.0) {
            r[0] = Alpha * v00;
            r[stride] = Alpha * v10;
         } else {
            r[0] = Beta * r[0] + Alpha * v00;
            r[stride] = Beta * r[stride] + Alpha * v10;
         }
      break;
      case 3:
         if (Beta == 0.0) {
            r[0] = Alpha * v00;
         } else {
            r[0] = Beta * r[0] + Alpha * v00;
         }
      break;
   }
}

#if defined AE_HAS_SSE2_INTRINSICS
// Calculate
//	r = Alpha a' (b0|b1) + Beta r
// where
// *	a, b0 and b1 are alinged k x 2 matrices stored in contiguous row-by-row storage,
//	b0 and b1 are two separate matrices stored in different locations,
// *	r is a 2 x 4 matrix stored in non-contiguous row-by-row storage and may be non-aligned.
// If Beta == 0.0, then r is ignored, rather than multiplied by zero.
// but if Alpha == 0.0, we still calculate the MM product even though it is multiplied by zero afterwards.
// Unlike the mm22 functions, this function does NOT support partial output formats for r - we always store the full 2 x 4 matrix.
// This function supports SSE2; it can be used when:
// *	AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
//	Otherwise, this function will be undefined.
// *	CurCPU contains CPU_SSE2 (checked at run-time)
//	Otherwise, calling this function will probably crash your system.
// If you want to know whether it is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
static void _ialglib_mm22x2_sse2(double Alpha, const double *a, const double *b0, const double *b1, ae_int_t k, double Beta, double *r, ae_int_t stride) {
// Calculate the 2 x 2 matrix product of two k x 2 matrices and store the result in v0, v1, v2 and v3 as follows:
//	r =	/ v0[0] v1[1] v2[0] v3[1] \
//		\ v1[0] v0[1] v3[0] v2[1] /
// vA0 is the current 1 x 2 block of a, va1 is the shuffle of va0,
// vB0 and vB1 are the two copies of the 1 x 2 blocks of b0 or b1 (both are the same copy - either b0 or b1).
// The results from multiplication by va0/va1 are stored in vb0/vb1 too.
   __m128d v0 = _mm_setzero_pd();
   __m128d v1 = _mm_setzero_pd();
   __m128d v2 = _mm_setzero_pd();
   __m128d v3 = _mm_setzero_pd();
   for (ae_int_t t = 0; t < k; t++) {
      __m128d va0 = _mm_load_pd(a);
      __m128d vb0 = _mm_load_pd(b0);
      __m128d va1 = _mm_load_pd(a);
      vb0 = _mm_mul_pd(va0, vb0);
      __m128d vb1 = _mm_load_pd(b0);
      v0 = _mm_add_pd(v0, vb0);
      vb1 = _mm_mul_pd(va1, vb1);
      vb0 = _mm_load_pd(b1);
      v1 = _mm_add_pd(v1, vb1);
      vb0 = _mm_mul_pd(va0, vb0);
      vb1 = _mm_load_pd(b1);
      v2 = _mm_add_pd(v2, vb0);
      vb1 = _mm_mul_pd(va1, vb1);
      v3 = _mm_add_pd(v3, vb1);
      a += 2;
      b0 += 2;
      b1 += 2;
   }
// Shuffle v1 and v3 (converting to more convenient storage format):
//	r =	/ v0[0] v1[0] v2[0] v3[0] \
//		\ v1[1] v0[1] v3[1] v2[1] /
// and unpack the results to
//	/ r00 r01 \
//	\ r10 r11 /
   __m128d valpha = _mm_load1_pd(&Alpha);
   v1 = _mm_shuffle_pd(v1, v1, 1);
   v3 = _mm_shuffle_pd(v3, v3, 1);
   __m128d r00 = _mm_mul_pd(_mm_unpacklo_pd(v0, v1), valpha);
   __m128d r10 = _mm_mul_pd(_mm_unpackhi_pd(v1, v0), valpha);
   __m128d r01 = _mm_mul_pd(_mm_unpacklo_pd(v2, v3), valpha);
   __m128d r11 = _mm_mul_pd(_mm_unpackhi_pd(v3, v2), valpha);
// Store the result.
   if (Beta == 0.0) {
      _mm_storeu_pd(r, r00);
      _mm_storeu_pd(r + 2, r01);
      _mm_storeu_pd(r + stride, r10);
      _mm_storeu_pd(r + stride + 2, r11);
   } else {
      __m128d vbeta = _mm_load1_pd(&Beta);
      _mm_storeu_pd(r, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r), vbeta), r00));
      _mm_storeu_pd(r + 2, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + 2), vbeta), r01));
      _mm_storeu_pd(r + stride, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + stride), vbeta), r10));
      _mm_storeu_pd(r + stride + 2, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + stride + 2), vbeta), r11));
   }
}
#endif

// Calculate
//	r = Alpha a' (b0|b1) + Beta r
// where
// *	a, b0 and b1 are alinged k x 2 matrices stored in contiguous row-by-row storage,
//	b0 and b1 are two separate matrices stored in different locations,
// *	r is a 2 x 4 matrix stored in non-contiguous row-by-row storage and may be non-aligned.
// If Beta == 0.0, then r is ignored, rather than multiplied by zero.
// but if Alpha == 0.0, we still calculate the MM product even though it is multiplied by zero afterwards.
// Unlike the mm22 functions, this function does NOT support partial output formats for r - we always store the full 2 x 4 matrix.
void _ialglib_mm22x2(double Alpha, const double *a, const double *b0, const double *b1, ae_int_t k, double Beta, double *r, ae_int_t stride) {
   _ialglib_mm22(Alpha, a, b0, k, Beta, r, stride, 0);
   _ialglib_mm22(Alpha, a, b1, k, Beta, r + 2, stride, 0);
}

//(@) Returned here from AlgLibInternal.cpp, because merely moving it here, seems to increase the speed of the GEMM routine!
// Fast rmatrixgemm() kernel: with AVX2/FMA support.
#if defined ALGLIB_NO_FAST_KERNELS
// ALGLIB Routine: Copyright 19.01.2010 by Sergey Bochkanov
#else
// ALGLIB Routine: Copyright 19.09.2021 by Sergey Bochkanov
#endif
bool ablasf_rgemm32basecase(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc) {
#if defined ALGLIB_NO_FAST_KERNELS || !defined _ALGLIB_HAS_AVX2_INTRINSICS
   return false;
#else
   const ae_int_t block_size = _ABLASF_BLOCK_SIZE;
   const ae_int_t micro_size = _ABLASF_MICRO_SIZE;
   ae_int_t (*ablasf_packblk)(const double *, ae_int_t, ae_int_t, ae_int_t, ae_int_t, double *, ae_int_t, ae_int_t) = k == 32 && block_size == 32 ? avx2_ablasf_packblkh32 : avx2_ablasf_packblkh;
   void (*ablasf_dotblk)(const double *, const double *, ae_int_t, ae_int_t, ae_int_t, double *, ae_int_t) = avx2_ablasf_dotblkh;
   void (*ablasf_daxpby)(ae_int_t, double, const double *, double, double *) = avx2_ablasf_daxpby;
// Determine CPU and kernel support.
   if (m > block_size || n > block_size || k > block_size || m == 0 || n == 0 || !(CurCPU & CPU_AVX2)) return false;
#   if defined _ALGLIB_HAS_FMA_INTRINSICS
   if (CurCPU & CPU_FMA) ablasf_dotblk = fma_ablasf_dotblkh;
#   endif
// Prepare c.
   double *cp = c->xyR[ic] + jc;
   ae_int_t stride_c = c->stride;
// Do we have alpha a b?
   if (alpha != 0.0 && k > 0) {
   // Prepare the structures.
      double *ap = a->xyR[ia] + ja; ae_int_t stride_a = a->stride;
      double *bp = b->xyR[ib] + jb; ae_int_t stride_b = b->stride;
      double _blka[_ABLASF_BLOCK_SIZE * _ABLASF_MICRO_SIZE + _ALGLIB_SIMD_ALIGNMENT_DOUBLES];
      double _blkb_long[_ABLASF_BLOCK_SIZE * _ABLASF_BLOCK_SIZE + _ALGLIB_SIMD_ALIGNMENT_DOUBLES];
      double _blkc[_ABLASF_MICRO_SIZE * _ABLASF_BLOCK_SIZE + _ALGLIB_SIMD_ALIGNMENT_DOUBLES];
      double *blka = (double *)ae_align(_blka, _ALGLIB_SIMD_ALIGNMENT_BYTES);
      double *storageb_long = (double *)ae_align(_blkb_long, _ALGLIB_SIMD_ALIGNMENT_BYTES);
      double *blkc = (double *)ae_align(_blkc, _ALGLIB_SIMD_ALIGNMENT_BYTES);
   // Pack transform(b) into precomputed block form.
      for (ae_int_t base1 = 0; base1 < n; base1 += micro_size) {
         const ae_int_t lim1 = n - base1 < micro_size ? n - base1 : micro_size;
         double *curb = storageb_long + base1 * block_size;
         ablasf_packblk(bp + (opb == 0 ? base1 : base1 * stride_b), stride_b, opb == 0 ? 1 : 0, k, lim1, curb, block_size, micro_size);
      }
   // Output.
      for (ae_int_t base0 = 0; base0 < m; base0 += micro_size) {
      // Load block row of transform(a).
         const ae_int_t lim0 = m - base0 < micro_size ? m - base0 : micro_size;
         const ae_int_t round_k = ablasf_packblk(ap + (opa == 0 ? base0 * stride_a : base0), stride_a, opa, k, lim0, blka, block_size, micro_size);
      // Compute block(a)' entire(b).
         for (ae_int_t base1 = 0; base1 < n; base1 += micro_size)
            ablasf_dotblk(blka, storageb_long + base1 * block_size, round_k, block_size, micro_size, blkc + base1, block_size);
      // Output block row of block(a)' entire(b).
         for (ae_int_t offs0 = 0; offs0 < lim0; offs0++)
            ablasf_daxpby(n, alpha, blkc + offs0 * block_size, beta, cp + (base0 + offs0) * stride_c);
      }
   } else {
   // No a b, just beta c (degenerate case, not optimized).
      if (beta == 0.0)
         for (ae_int_t out0 = 0; out0 < m; out0++) for (ae_int_t out1 = 0; out1 < n; out1++) cp[out0 * stride_c + out1] = 0.0;
      else if (beta != 1.0)
         for (ae_int_t out0 = 0; out0 < m; out0++) for (ae_int_t out1 = 0; out1 < n; out1++) cp[out0 * stride_c + out1] *= beta;
   }
   return true;
#endif
}
} // end of namespace alglib_impl

namespace alglib {
// Declarations for C++-related functionality.

// Exception handling.
#if !defined AE_NO_EXCEPTIONS
ap_error::ap_error() { msg = alglib_impl::CurMsg; }
ap_error::ap_error(const char *Msg) { msg = Msg; }
void ap_error::make_assertion(bool Cond) { if (!Cond) ThrowError(""); }
void ap_error::make_assertion(bool Cond, const char *Msg) { if (!Cond) ThrowError(Msg); }
#else
static const char *_alglib_last_error = NULL;
static void set_error_flag(const char *Msg) { _alglib_last_error = Msg == NULL ? "set_error_flag: unknown error" : Msg; }
void set_error_msg() { _alglib_last_error = alglib_impl::CurMsg; }

bool get_error_flag(const char **MsgP/* = NULL*/) {
   if (_alglib_last_error == NULL) return false;
   if (MsgP != NULL) *MsgP = _alglib_last_error;
   return true;
}

void clear_error_flag() { _alglib_last_error = NULL; }
#endif

#if 0
// Global and local constants and variables.
const double machineepsilon = 5.0E-16, maxrealnumber = 1.0E300, minrealnumber = 1.0E-300;
const double pi = 3.1415926535897932384626433832795;
static const ae_int_t ByteOrder = alglib_impl::ByteOrder;
#endif

// Standard functions.
bool isneginf(double x) { return isinf(x) && signbit(x); }
bool isposinf(double x) { return isinf(x) && !signbit(x); }

int minint(int x, int y) { return x > y ? y : x; }
int maxint(int x, int y) { return x > y ? x : y; }
int sign(double x) { return x > 0.0 ? +1 : x < 0.0 ? -1 : 0; }
int iround(double x) { return int(round(x)); }
int itrunc(double x) { return int(trunc(x)); }
int ifloor(double x) { return int(floor(x)); }
int iceil(double x) { return int(ceil(x)); }

double minreal(double x, double y) { return x > y ? y : x; }
double maxreal(double x, double y) { return x > y ? x : y; }
double sqr(double x) { return x * x; }

double randomreal() {
   const double mx = alglib_impl::ae_rand_max() + 1.0;
   double i = alglib_impl::ae_rand();
   return (i + alglib_impl::ae_rand() / mx) / mx;
}

double randommid() {
   const double mx = alglib_impl::ae_rand_max() + 1.0;
   double i = alglib_impl::ae_rand();
   return 2.0 * (i + alglib_impl::ae_rand() / mx) / mx - 1.0;
}

ae_int_t randominteger(ae_int_t maxv) { return alglib_impl::ae_rand() % maxv; }

bool randombool(double p/* = 0.5*/) {
   const double mx = alglib_impl::ae_rand_max() + 1.0;
   double i = alglib_impl::ae_rand();
   return i + alglib_impl::ae_rand() / mx <= p * mx;
}

// Complex number with double precision.
complex &complex::operator=(const double &v) { x = v, y = 0.0; return *this; }
complex &complex::operator=(const complex &z) { x = z.x, y = z.y; return *this; }
complex &complex::operator+=(const double &v) { x += v; return *this; }
complex &complex::operator+=(const complex &z) { x += z.x, y += z.y; return *this; }
complex &complex::operator-=(const double &v) { x -= v; return *this; }
complex &complex::operator-=(const complex &z) { x -= z.x, y -= z.y; return *this; }
complex &complex::operator*=(const double &v) { x *= v, y *= v; return *this; }
complex &complex::operator*=(const complex &z) { double t = x * z.x - y * z.y; y = x * z.y + y * z.x, x = t; return *this; }
complex &complex::operator/=(const double &v) { x /= v, y /= v; return *this; }

complex &complex::operator/=(const complex &z) {
   double zx = z.x, zy = z.y;
   if (fabs(zy) < fabs(zx)) {
      double e = zy / zx, f = zx + zy * e;
      return *this = complex((x + y * e) / f, (y - x * e) / f);
   } else {
      double e = zx / zy, f = zy + zx * e;
      return *this = complex((x * e + y) / f, (y * e - x) / f);
   }
}

// alglib_impl-alglib gateway.
static inline alglib_impl::complex complex_from_c(complex z) { return alglib_impl::complex_from_d(z.x, z.y); }

#if !defined AE_NO_EXCEPTIONS
std::string complex::tostring(int _dps) const {
   int dps = _dps >= 0 ? _dps : -_dps;
   if (dps <= 0 || dps >= 20)
      ThrowError("complex::tostring: incorrect dps");
// Handle IEEE-special quantities.
   if (isnan(x) || isnan(y))
      return "NAN";
   else if (isinf(x) || isinf(y))
      return "INF";
// Generate the mask.
   char mask[0x20];
   if (sprintf(mask, "%%.%d%s", dps, _dps >= 0 ? "f" : "e") >= (int)sizeof mask)
      ThrowError("complex::tostring: buffer overflow");
// Print |x|, |y| and zero with the same mask and compare.
   char buf_x[0x20];
   if (sprintf(buf_x, mask, fabs(x)) >= (int)sizeof buf_x)
      ThrowError("complex::tostring: buffer overflow");
   char buf_y[0x20];
   if (sprintf(buf_y, mask, fabs(y)) >= (int)sizeof buf_y)
      ThrowError("complex::tostring: buffer overflow");
   char buf_zero[0x20];
   if (sprintf(buf_zero, mask, 0.0) >= (int)sizeof buf_zero)
      ThrowError("complex::tostring: buffer overflow");
// Different zero/non-zero patterns.
   if (strcmp(buf_x, buf_zero) != 0 && strcmp(buf_y, buf_zero) != 0)
      return std::string(x > 0.0 ? "" : "-") + buf_x + (y > 0.0 ? "+" : "-") + buf_y + "i";
   else if (strcmp(buf_x, buf_zero) != 0 && strcmp(buf_y, buf_zero) == 0)
      return std::string(x > 0.0 ? "" : "-") + buf_x;
   else if (strcmp(buf_x, buf_zero) == 0 && strcmp(buf_y, buf_zero) != 0)
      return std::string(y > 0.0 ? "" : "-") + buf_y + "i";
   else
      return std::string("0");
}
#endif

bool operator==(const complex &A, const complex &B) { return A.x == B.x && A.y == B.y; }
bool operator!=(const complex &A, const complex &B) { return A.x != B.x || A.y != B.y; }

const complex operator+(const complex &A) { return A; }
const complex operator-(const complex &A) { return complex(-A.x, -A.y); }

const complex operator+(const complex &A, const complex &B) { complex r = A; r += B; return r; }
const complex operator+(const complex &A, const double &B) { complex r = A; r += B; return r; }
const complex operator+(const double &A, const complex &B) { complex r = B; r += A; return r; }

const complex operator-(const complex &A, const complex &B) { complex r = A; r -= B; return r; }
const complex operator-(const complex &A, const double &B) { complex r = A; r -= B; return r; }
const complex operator-(const double &A, const complex &B) { complex r = A; r -= B; return r; }

const complex operator*(const complex &A, const complex &B) { return complex(A.x * B.x - A.y * B.y, A.x * B.y + A.y * B.x); }
const complex operator*(const complex &A, const double &B) { return complex(A.x * B, A.y * B); }
const complex operator*(const double &A, const complex &B) { return complex(A * B.x, A * B.y); }

const complex operator/(const complex &A, const complex &B) {
   double Ax = A.x, Ay = A.y;
   double Bx = B.x, By = B.y;
   if (fabs(By) < fabs(Bx)) {
      double e = By / Bx, f = Bx + By * e;
      return complex((Ax + Ay * e) / f, (Ay - Ax * e) / f);
   } else {
      double e = Bx / By, f = By + Bx * e;
      return complex((Ax * e + Ay) / f, (Ay * e - Ax) / f);
   }
}

const complex operator/(const complex &A, const double &B) { return complex(A.x / B, A.y / B); }

const complex operator/(const double &A, const complex &B) {
   double Bx = B.x, By = B.y;
   if (fabs(By) < fabs(Bx)) {
      double e = By / Bx, f = Bx + By * e;
      return complex(A / f, -A * e / f);
   } else {
      double e = Bx / By, f = By + Bx * e;
      return complex(A * e / f, -A / f);
   }
}

double abscomplex(const complex &A) {
   double Ax = fabs(A.x), Ay = fabs(A.y), Lo = Ax < Ay ? Ax : Ay, Hi = Ax > Ay ? Ax : Ay;
   if (Lo == 0.0) return Hi;
   else {
      double Span = Lo / Hi;
      return Hi * sqrt(1.0 + Span * Span);
   }
}

complex conj(const complex &A) { return complex(A.x, -A.y); }
complex csqr(const complex &A) { double Ax = A.x, Ay = A.y; return complex(Ax * Ax - Ay * Ay, 2.0 * Ax * Ay); }

#ifdef AE_HPC
ae_int_t getnworkers() { return alglib_impl::getnworkers(); }
void setnworkers(ae_int_t nworkers) { alglib_impl::setnworkers(nworkers); }
ae_int_t ae_cores_count() { return alglib_impl::ae_cores_count(); }
alglib_impl::ae_uint64_t _ae_get_global_threading() { return alglib_impl::ae_get_global_threading(); }
void _ae_set_global_threading(alglib_impl::ae_uint64_t flg_value) { alglib_impl::ae_set_global_threading(flg_value); }
void setglobalthreading(const xparams settings) { alglib_impl::ae_set_global_threading((ae_uint64_t)settings); }
#else
ae_int_t getnworkers() { return 1; }
void setnworkers(ae_int_t nworkers) { }
ae_int_t ae_cores_count() { return 1; }
alglib_impl::ae_uint64_t _ae_get_global_threading() { return SerTH; }
void _ae_set_global_threading(alglib_impl::ae_uint64_t flg_value) { }
void setglobalthreading(const xparams settings) { }
#endif

// Level 1 Real BLAS functions.
// Handle the case of unit stride specially and separately with optimization.
//(@) In the double v{dotproduct,move,moveneg,add}()
//(@) this case originally cut the loop counter in half (or quarter) and did two (or four) steps per loop.
//(@) This was removed, the need for it was not clearly explained and it appears to have no visible effect on performance.

double vdotproduct(const double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   double C = 0.0;
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         C += *A * *B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         C += *A * *B;
   return C;
}

double vdotproduct(const double *A, const double *B, ae_int_t N) {
   return vdotproduct(A, 1, B, 1, N);
}

complex vdotproduct(const complex *A, ae_int_t dA, const char *CjA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   double Cx = 0.0, Cy = 0.0;
   bool ConjA = CjA[0] != 'N' && CjA[0] != 'n';
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (ConjA)
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = -A->y;
            double Bx = B->x, By = -B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = -A->y;
            double Bx = B->x, By = B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
   else
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = A->y;
            double Bx = B->x, By = -B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            double Ax = A->x, Ay = A->y;
            double Bx = B->x, By = B->y;
            Cx += Ax * Bx - Ay * By, Cy += Ax * By + Ay * Bx;
         }
   return complex(Cx, Cy);
}

complex vdotproduct(const complex *A, const complex *B, ae_int_t N) {
   return vdotproduct(A, 1, "N", B, 1, "N", N);
}

void vmove(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         *A = *B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A = *B;
}

void vmove(double *A, const double *B, ae_int_t N) {
   vmove(A, 1, B, 1, N);
}

void vmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = B->x;
            A->y = -B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++)
            *A = *B;
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = B->x;
            A->y = -B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
            *A = *B;
}

void vmove(complex *A, const complex *B, ae_int_t N) {
   vmove(A, 1, B, 1, "N", N);
}

void vmoveneg(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         *A = -*B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A = -*B;
}

void vmoveneg(double *A, const double *B, ae_int_t N) {
   vmoveneg(A, 1, B, 1, N);
}

void vmoveneg(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = -B->x;
            A->y = B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = -B->x;
            A->y = -B->y;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = -B->x;
            A->y = B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = -B->x;
            A->y = -B->y;
         }
}

void vmoveneg(complex *A, const complex *B, ae_int_t N) {
   vmoveneg(A, 1, B, 1, "N", N);
}

void vmove(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha) {
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         *A = Alpha * *B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A = Alpha * *B;
}

void vmove(double *A, const double *B, ae_int_t N, double Alpha) {
   vmove(A, 1, B, 1, N, Alpha);
}

void vmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = Alpha * B->x;
            A->y = -Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = Alpha * B->x;
            A->y = Alpha * B->y;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = Alpha * B->x;
            A->y = -Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = Alpha * B->x;
            A->y = Alpha * B->y;
         }
}

void vmove(complex *A, const complex *B, ae_int_t N, double Alpha) {
   vmove(A, 1, B, 1, "N", N, Alpha);
}

void vmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   double Ax = Alpha.x, Ay = Alpha.y;
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = Ax * B->x + Ay * B->y;
            A->y = -Ax * B->y + Ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x = Ax * B->x - Ay * B->y;
            A->y = Ax * B->y + Ay * B->x;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = Ax * B->x + Ay * B->y;
            A->y = -Ax * B->y + Ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x = Ax * B->x - Ay * B->y;
            A->y = Ax * B->y + Ay * B->x;
         }
}

void vmove(complex *A, const complex *B, ae_int_t N, complex Alpha) {
   vmove(A, 1, B, 1, "N", N, Alpha);
}

void vadd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         *A += *B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A += *B;
}

void vadd(double *A, const double *B, ae_int_t N) {
   vadd(A, 1, B, 1, N);
}

void vadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += B->x;
            A->y -= B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += B->x;
            A->y += B->y;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += B->x;
            A->y -= B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += B->x;
            A->y += B->y;
         }
}

void vadd(complex *A, const complex *B, ae_int_t N) {
   vadd(A, 1, B, 1, "N", N);
}

void vadd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha) {
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         *A += Alpha * *B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A += Alpha * *B;
}

void vadd(double *A, const double *B, ae_int_t N, double Alpha) {
   vadd(A, 1, B, 1, N, Alpha);
}

void vadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += Alpha * B->x;
            A->y -= Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += Alpha * B->x;
            A->y += Alpha * B->y;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += Alpha * B->x;
            A->y -= Alpha * B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += Alpha * B->x;
            A->y += Alpha * B->y;
         }
}

void vadd(complex *A, const complex *B, ae_int_t N, double Alpha) {
   vadd(A, 1, B, 1, "N", N, Alpha);
}

void vadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   double Ax = Alpha.x, Ay = Alpha.y;
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += Ax * B->x + Ay * B->y;
            A->y -= Ax * B->y - Ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x += Ax * B->x - Ay * B->y;
            A->y += Ax * B->y + Ay * B->x;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += Ax * B->x + Ay * B->y;
            A->y -= Ax * B->y - Ay * B->x;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x += Ax * B->x - Ay * B->y;
            A->y += Ax * B->y + Ay * B->x;
         }
}

void vadd(complex *A, const complex *B, ae_int_t N, complex Alpha) {
   vadd(A, 1, B, 1, "N", N, Alpha);
}

void vsub(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N) {
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++, B++)
         *A -= *B;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA, B += dB)
         *A -= *B;
}

void vsub(double *A, const double *B, ae_int_t N) {
   vsub(A, 1, B, 1, N);
}

void vsub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N) {
   bool ConjB = CjB[0] != 'N' && CjB[0] != 'n';
   if (dA == 1 && dB == 1) // The optimized case for unit stride.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x -= B->x;
            A->y += B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A++, B++) {
            A->x -= B->x;
            A->y -= B->y;
         }
   else // The general case.
      if (ConjB)
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x -= B->x;
            A->y += B->y;
         }
      else
         for (ae_int_t n = 0; n < N; n++, A += dA, B += dB) {
            A->x -= B->x;
            A->y -= B->y;
         }
}

void vsub(complex *A, const complex *B, ae_int_t N) {
   vsub(A, 1, B, 1, "N", N);
}

void vsub(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha) {
   vadd(A, dA, B, dB, N, -Alpha);
}

void vsub(double *A, const double *B, ae_int_t N, double Alpha) {
   vadd(A, 1, B, 1, N, -Alpha);
}

void vsub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha) {
   vadd(A, dA, B, dB, CjB, N, -Alpha);
}

void vsub(complex *A, const complex *B, ae_int_t N, double Alpha) {
   vadd(A, 1, B, 1, "N", N, -Alpha);
}

void vsub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha) {
   vadd(A, dA, B, dB, CjB, N, -Alpha);
}

void vsub(complex *A, const complex *B, ae_int_t N, complex Alpha) {
   vadd(A, 1, B, 1, "N", N, -Alpha);
}

void vmul(double *A, ae_int_t dA, ae_int_t N, double Alpha) {
   if (dA == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++)
         *A *= Alpha;
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA)
         *A *= Alpha;
}

void vmul(double *A, ae_int_t N, double Alpha) {
   vmul(A, 1, N, Alpha);
}

void vmul(complex *A, ae_int_t dA, ae_int_t N, double Alpha) {
   if (dA == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++) {
         A->x *= Alpha;
         A->y *= Alpha;
      }
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA) {
         A->x *= Alpha;
         A->y *= Alpha;
      }
}

void vmul(complex *A, ae_int_t N, double Alpha) {
   vmul(A, 1, N, Alpha);
}

void vmul(complex *A, ae_int_t dA, ae_int_t N, complex Alpha) {
   double Ax = Alpha.x, Ay = Alpha.y;
   if (dA == 1) // The optimized case for unit stride.
      for (ae_int_t n = 0; n < N; n++, A++) {
         double Bx = A->x, By = A->y;
         A->x = Ax * Bx - Ay * By;
         A->y = Ax * By + Ay * Bx;
      }
   else // The general case.
      for (ae_int_t n = 0; n < N; n++, A += dA) {
         double Bx = A->x, By = A->y;
         A->x = Ax * Bx - Ay * By;
         A->y = Ax * By + Ay * Bx;
      }
}

void vmul(complex *A, ae_int_t N, complex Alpha) {
   vmul(A, 1, N, Alpha);
}

// Matrices and Vectors: I/O.

static bool strimatch(const char *s1, const char *s2) {
// Handle the special cases.
   if (s1 == NULL || s2 == NULL) return s1 == s2;
// Compare.
   else while (true) {
      int c1 = *s1++, c2 = *s2++;
      if (c1 == '\0' || c2 == '\0') return c1 == c2;
      else if (tolower(c1) != tolower(c2)) return false;
   }
}

#if !defined AE_NO_EXCEPTIONS
// Filter out all of the spaces from the string s.
// Return a string allocated with ae_malloc().
// On allocation failure, return NULL.
static char *filter_spaces(const char *s) {
   size_t n = strlen(s);
   char *r = (char *)alglib_impl::ae_malloc(n + 1);
   if (r == NULL) return r;
   char *r0 = r;
   for (size_t i = 0; i <= n; i++, s++)
      if (!isspace(*s)) *r0++ = *s;
   return r;
}

static void str_vector_create(const char *src, bool match_head_only, std::vector<const char *> *p_vec) {
// Parse the beginning of the string.
// Try to handle "[]".
   p_vec->clear();
   if (*src != '[') ThrowError("Incorrect initializer for vector");
   else if (*++src == ']') return;
   p_vec->push_back(src);
   while (true)
      if (*src == '\0')
         ThrowError("Incorrect initializer for vector");
      else if (*src == ']')
         if (src[1] == '\0' || !match_head_only)
            return;
         else
            ThrowError("Incorrect initializer for vector");
      else if (*src++ == ',')
         p_vec->push_back(src);
}

static void str_matrix_create(const char *src, std::vector< std::vector<const char *> > *p_mat) {
   p_mat->clear();
// Try to handle "[[]]".
   if (strcmp(src, "[[]]") == 0) return;
// Parse a non-empty string.
   if (*src != '[') ThrowError("Incorrect initializer for matrix");
   src++;
   while (true) {
      p_mat->push_back(std::vector<const char *>());
      str_vector_create(src, false, &p_mat->back());
      if (p_mat->back().size() == 0 || p_mat->back().size() != (*p_mat)[0].size())
         ThrowError("Incorrect initializer for matrix");
      src = strchr(src, ']');
      if (src == NULL)
         ThrowError("Incorrect initializer for matrix");
      else if (*++src == ',')
         src++;
      else if (*src == ']')
         break;
      else
         ThrowError("Incorrect initializer for matrix");
   }
   if (*++src != '\0')
      ThrowError("Incorrect initializer for matrix");
}

static bool parse_bool_delim(const char *s, const char *delim) {
// Try to parse "false", "true", or otherwise fail.
   const char *p = "false";
   char buf[8];
   memset(buf, 0, sizeof buf);
   strncpy(buf, s, strlen(p));
   if (strimatch(buf, p)) {
      if (s[strlen(p)] == '\0' || strchr(delim, s[strlen(p)]) == NULL)
         ThrowError("Cannot parse value");
      return false;
   }
   p = "true";
   memset(buf, 0, sizeof buf);
   strncpy(buf, s, strlen(p));
   if (strimatch(buf, p)) {
      if (s[strlen(p)] == '\0' || strchr(delim, s[strlen(p)]) == NULL)
         ThrowError("Cannot parse value");
      return true;
   }
   ThrowError("Cannot parse value");
}

static ae_int_t parse_int_delim(const char *s, const char *delim) {
   const char *p = s;
// Check the string structure for a leading sign, at least one digit and a delimiter.
   if (*s == '-' || *s == '+')
      s++;
   if (*s == '\0' || strchr("1234567890", *s) == NULL)
      ThrowError("Cannot parse value");
   while (*s != '\0' && strchr("1234567890", *s) != NULL)
      s++;
   if (*s == '\0' || strchr(delim, *s) == NULL)
      ThrowError("Cannot parse value");
// Convert and ensure that the value fits into ae_int_t.
   s = p;
   long long_val = atol(s);
   volatile ae_int_t ae_val = long_val;
   if (ae_val != long_val)
      ThrowError("Cannot parse value");
   return ae_val;
}

static bool _parse_real_delim(const char *s, const char *delim, double *result, const char **new_s) {
   const char *p = s;
// Check the string structure and decide what to do.
   int isign = 1;
   if (*s == '-' || *s == '+') {
      isign = *s == '-' ? -1 : +1;
      s++;
   }
   char buf[0x40];
   memset(buf, 0, sizeof buf);
   strncpy(buf, s, 3);
   if (!strimatch(buf, "nan") && !strimatch(buf, "inf")) {
   // [Sign] [Digits] [.] [Digits] [e|E [Sign] Digits].
      bool has_digits = false;
      if (*s != '\0' && strchr("1234567890", *s) != NULL) {
         has_digits = true;
         while (*s != '\0' && strchr("1234567890", *s) != NULL)
            s++;
      }
      if (*s == '.')
         s++;
      if (*s != '\0' && strchr("1234567890", *s) != NULL) {
         has_digits = true;
         while (*s != '\0' && strchr("1234567890", *s) != NULL)
            s++;
      }
      if (!has_digits)
         return false;
      if (*s == 'e' || *s == 'E') {
         s++;
         if (*s == '-' || *s == '+')
            s++;
         if (*s == '\0' || strchr("1234567890", *s) == NULL)
            return false;
         while (*s != '\0' && strchr("1234567890", *s) != NULL)
            s++;
      }
      if (*s == '\0' || strchr(delim, *s) == NULL)
         return false;
      *new_s = s;
   // Finite value conversion.
      if (*new_s - p >= sizeof buf)
         return false;
      strncpy(buf, p, (size_t)(*new_s - p));
      buf[*new_s - p] = '\0';
      lconv *loc = localeconv();
      char *t = strchr(buf, '.');
      if (t != NULL)
         *t = *loc->decimal_point;
      *result = atof(buf);
      return true;
   } else {
   // Check the delimiter and update *new_s.
      s += 3;
      if (*s == '\0' || strchr(delim, *s) == NULL)
         return false;
      *new_s = s;
   // NAN, INF conversion.
      if (strimatch(buf, "nan"))
         *result = NAN;
      if (strimatch(buf, "inf"))
         *result = isign > 0 ? +INFINITY : -INFINITY;
      return true;
   }
}

static double parse_real_delim(const char *s, const char *delim) {
   double result;
   const char *new_s;
   if (!_parse_real_delim(s, delim, &result, &new_s))
      ThrowError("Cannot parse value");
   return result;
}

static complex parse_complex_delim(const char *s, const char *delim) {
// Parse as a real value.
   double d_result;
   const char *new_s;
   if (_parse_real_delim(s, delim, &d_result, &new_s))
      return d_result;
// Parse as "x + yi" or "x - yi".
   complex c_result;
   if (_parse_real_delim(s, "+-", &c_result.x, &new_s)) {
      s = new_s;
      if (!_parse_real_delim(s, "i", &c_result.y, &new_s))
         ThrowError("Cannot parse value");
      s = new_s + 1;
      if (*s == '\0' || strchr(delim, *s) == NULL)
         ThrowError("Cannot parse value");
      return c_result;
   }
// Parse as a complex value "yi + x" or "yi - x".
   if (_parse_real_delim(s, "i", &c_result.y, &new_s)) {
      s = new_s + 1;
      if (*s == '\0')
         ThrowError("Cannot parse value");
      if (strchr(delim, *s) != NULL) {
         c_result.x = 0.0;
         return c_result;
      }
      if (strchr("+-", *s) != NULL) {
         if (!_parse_real_delim(s, delim, &c_result.x, &new_s))
            ThrowError("Cannot parse value");
         return c_result;
      }
      ThrowError("Cannot parse value");
   }
// Error.
   ThrowError("Cannot parse value");
}

static std::string arraytostring(const bool *ptr, ae_int_t n) {
   std::string result = "[";
   for (ae_int_t i = 0; i < n; i++) {
      if (i != 0) result += ",";
      result += ptr[i] ? "true" : "false";
   }
   return result += "]";
}

static std::string arraytostring(const ae_int_t *ptr, ae_int_t n) {
   std::string result = "[";
   for (ae_int_t i = 0; i < n; i++) {
      char buf[0x40];
      if (sprintf(buf, i == 0 ? "%ld" : ",%ld", long (ptr[i])) >= (int)sizeof buf)
         ThrowError("arraytostring: buffer overflow");
      result += buf;
   }
   return result += "]";
}

static std::string arraytostring(const double *ptr, ae_int_t n, int _dps) {
   int dps = _dps >= 0 ? _dps : -_dps;
   dps = dps <= 50 ? dps : 50;
   std::string result = "[";
   char mask1[0x40];
   if (sprintf(mask1, "%%.%d%s", dps, _dps >= 0 ? "f" : "e") >= (int)sizeof mask1)
      ThrowError("arraytostring: buffer overflow");
   char mask2[0x50];
   if (sprintf(mask2, ",%s", mask1) >= (int)sizeof mask2)
      ThrowError("arraytostring: buffer overflow");
   for (ae_int_t i = 0; i < n; i++) {
      char buf[0x40];
      buf[0] = '\0';
      if (isfinite(ptr[i])) {
         if (sprintf(buf, i == 0 ? mask1 : mask2, double(ptr[i])) >= (int)sizeof buf)
            ThrowError("arraytostring: buffer overflow");
      } else if (isnan(ptr[i]))
         strcpy(buf, i == 0 ? "NAN" : ",NAN");
      else if (isposinf(ptr[i]))
         strcpy(buf, i == 0 ? "+INF" : ",+INF");
      else if (isneginf(ptr[i]))
         strcpy(buf, i == 0 ? "-INF" : ",-INF");
      result += buf;
   }
   return result += "]";
}

static std::string arraytostring(const complex *ptr, ae_int_t n, int dps) {
   std::string result = "[";
   for (ae_int_t i = 0; i < n; i++) {
      if (i != 0) result += ",";
      result += ptr[i].tostring(dps);
   }
   return result += "]";
}
#endif

// Matrices and Vectors: Wrappers.

ae_vector_wrapper::ae_vector_wrapper(alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      owner = true, This = NULL, set_error_msg(); return;
#endif
   }
   owner = true, This = &Obj, memset(This, 0, sizeof *This), ae_vector_init(This, 0, datatype, false);
   alglib_impl::ae_state_clear();
}

ae_vector_wrapper::ae_vector_wrapper(alglib_impl::ae_vector *e_ptr, alglib_impl::ae_datatype datatype) {
   if (e_ptr == NULL || e_ptr->datatype != datatype) {
      const char *msg = "ae_vector_wrapper::ae_vector_wrapper: datatype check failed";
#if !defined AE_NO_EXCEPTIONS
      ThrowError(msg);
#else
      owner = true, This = NULL, set_error_flag(msg); return;
#endif
   }
   owner = false, This = e_ptr;
}

ae_vector_wrapper::ae_vector_wrapper(const ae_vector_wrapper &rhs, alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      owner = true, This = NULL, set_error_msg(); return;
#endif
   }
   if (rhs.This != NULL)
      alglib_impl::ae_assert(rhs.This->datatype == datatype, "ae_vector_wrapper::ae_vector_wrapper: datatype check failed");
   owner = true, This = &Obj, memset(This, 0, sizeof *This);
   rhs.This == NULL ? ae_vector_init(This, 0, datatype, false) : ae_vector_copy(This, rhs.This, false);
   alglib_impl::ae_state_clear();
}

ae_vector_wrapper::~ae_vector_wrapper() { if (This == &Obj) ae_vector_free(This, false); }

void ae_vector_wrapper::setlength(ae_int_t iLen) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ae_assert(This != NULL, "ae_vector_wrapper::setlength: This == NULL (array was not correctly initialized)");
   alglib_impl::ae_assert(owner, "ae_vector_wrapper::setlength: This is frozen proxy array");
   alglib_impl::ae_vector_set_length(This, iLen);
   alglib_impl::ae_state_clear();
}

ae_int_t ae_vector_wrapper::length() const { return This == NULL ? 0 : This->cnt; }

void ae_vector_wrapper::attach_to(alglib_impl::x_vector *new_ptr) {
   if (This == &Obj) ae_vector_free(This, false);
   owner = false, This = &Obj, memset(This, 0, sizeof *This), ae_vector_init_attach_to_x(This, new_ptr, false);
}

const ae_vector_wrapper &ae_vector_wrapper::assign(const ae_vector_wrapper &rhs) {
   if (this == &rhs) return *this;
   alglib_impl::ae_state_init();
   TryCatch(*this)
   alglib_impl::ae_assert(This != NULL, "ae_vector_wrapper::assign: incorrect assignment (uninitialized destination)");
// Assignment to proxy object.
   alglib_impl::ae_assert(rhs.This != NULL, "ae_vector_wrapper::assign: incorrect assignment (uninitialized source)");
   alglib_impl::ae_assert(rhs.This->datatype == This->datatype, "ae_vector_wrapper::assign: incorrect assignment to array (types do not match)");
   if (!owner) alglib_impl::ae_assert(rhs.This->cnt == This->cnt, "ae_vector_wrapper::assign: incorrect assignment to proxy array (sizes do not match)");
   if (rhs.This->cnt != This->cnt) ae_vector_set_length(This, rhs.This->cnt);
   memcpy(This->xX, rhs.This->xX, This->cnt * alglib_impl::ae_sizeof(This->datatype));
   alglib_impl::ae_state_clear();
   return *this;
}

#if !defined AE_NO_EXCEPTIONS
ae_vector_wrapper::ae_vector_wrapper(const char *s, alglib_impl::ae_datatype datatype) {
   char *p = filter_spaces(s); if (p == NULL) ThrowError("ae_vector_wrapper::ae_vector_wrapper: allocation error");
   try {
      std::vector<const char *> svec; str_vector_create(p, true, &svec);
   {
      alglib_impl::ae_state_init();
      TryCatch()
      owner = true, This = &Obj, memset(This, 0, sizeof *This), ae_vector_init(This, (ae_int_t)svec.size(), datatype, false);
      alglib_impl::ae_state_clear();
   }
      for (size_t i = 0; i < svec.size(); i++) switch (datatype) {
      // case alglib_impl::DT_BYTE: // The same as alglib_impl::DT_BOOL.
         case alglib_impl::DT_BOOL: This->xB[i] = parse_bool_delim(svec[i], ",]"); break;
         case alglib_impl::DT_INT: This->xZ[i] = parse_int_delim(svec[i], ",]"); break;
         case alglib_impl::DT_REAL: This->xR[i] = parse_real_delim(svec[i], ",]"); break;
         case alglib_impl::DT_COMPLEX: This->xC[i] = complex_from_c(parse_complex_delim(svec[i], ",]")); break;
      }
      alglib_impl::ae_free(p);
   } catch(...) {
      alglib_impl::ae_free(p);
      throw;
   }
}
#endif

boolean_1d_array::boolean_1d_array(): ae_vector_wrapper(alglib_impl::DT_BOOL) { }
boolean_1d_array::boolean_1d_array(alglib_impl::ae_vector *p): ae_vector_wrapper(p, alglib_impl::DT_BOOL) { }
boolean_1d_array::boolean_1d_array(const boolean_1d_array &rhs): ae_vector_wrapper(rhs, alglib_impl::DT_BOOL) { }
#if !defined AE_NO_EXCEPTIONS
boolean_1d_array::boolean_1d_array(const char *s): ae_vector_wrapper(s, alglib_impl::DT_BOOL) { }
std::string boolean_1d_array::tostring() const { return length() == 0 ? "[]" : arraytostring(&(operator()(0)), length()); }
#endif
boolean_1d_array::~boolean_1d_array() { }

const boolean_1d_array &boolean_1d_array::operator=(const boolean_1d_array &rhs) { assign(rhs); return *this; }
const bool &boolean_1d_array::operator()(ae_int_t i) const { return This->xB[i]; }
bool &boolean_1d_array::operator()(ae_int_t i) { return This->xB[i]; }
const bool &boolean_1d_array::operator[](ae_int_t i) const { return This->xB[i]; }
bool &boolean_1d_array::operator[](ae_int_t i) { return This->xB[i]; }

void boolean_1d_array::setcontent(ae_int_t iLen, const bool *pContent) {
// Handle possible exception-free errors.
   setlength(iLen);
// Copy, if its size matches.
   if (This != NULL && This->cnt == iLen)
      for (ae_int_t i = 0; i < iLen; i++) This->xB[i] = pContent[i];
}

const bool *boolean_1d_array::getcontent() const { return This->xB; }
bool *boolean_1d_array::getcontent() { return This->xB; }

integer_1d_array::integer_1d_array(): ae_vector_wrapper(alglib_impl::DT_INT) { }
integer_1d_array::integer_1d_array(alglib_impl::ae_vector *p): ae_vector_wrapper(p, alglib_impl::DT_INT) { }
integer_1d_array::integer_1d_array(const integer_1d_array &rhs): ae_vector_wrapper(rhs, alglib_impl::DT_INT) { }
#if !defined AE_NO_EXCEPTIONS
integer_1d_array::integer_1d_array(const char *s): ae_vector_wrapper(s, alglib_impl::DT_INT) { }
std::string integer_1d_array::tostring() const { return length() == 0 ? "[]" : arraytostring(&operator()(0), length()); }
#endif
integer_1d_array::~integer_1d_array() { }

const integer_1d_array &integer_1d_array::operator=(const integer_1d_array &rhs) { assign(rhs); return *this; }
const ae_int_t &integer_1d_array::operator()(ae_int_t i) const { return This->xZ[i]; }
ae_int_t &integer_1d_array::operator()(ae_int_t i) { return This->xZ[i]; }
const ae_int_t &integer_1d_array::operator[](ae_int_t i) const { return This->xZ[i]; }
ae_int_t &integer_1d_array::operator[](ae_int_t i) { return This->xZ[i]; }

void integer_1d_array::setcontent(ae_int_t iLen, const ae_int_t *pContent) {
// Handle possible exception-free errors.
   setlength(iLen);
// Copy, if its size matches.
   if (This != NULL && This->cnt == iLen)
      for (ae_int_t i = 0; i < iLen; i++) This->xZ[i] = pContent[i];
}

const ae_int_t *integer_1d_array::getcontent() const { return This->xZ; }
ae_int_t *integer_1d_array::getcontent() { return This->xZ; }

real_1d_array::real_1d_array(): ae_vector_wrapper(alglib_impl::DT_REAL) { }
real_1d_array::real_1d_array(alglib_impl::ae_vector *p): ae_vector_wrapper(p, alglib_impl::DT_REAL) { }
real_1d_array::real_1d_array(const real_1d_array &rhs): ae_vector_wrapper(rhs, alglib_impl::DT_REAL) { }
#if !defined AE_NO_EXCEPTIONS
real_1d_array::real_1d_array(const char *s): ae_vector_wrapper(s, alglib_impl::DT_REAL) { }
std::string real_1d_array::tostring(int dps) const { return length() == 0 ? "[]" : arraytostring(&operator()(0), length(), dps); }
#endif
real_1d_array::~real_1d_array() { }

const real_1d_array &real_1d_array::operator=(const real_1d_array &rhs) { assign(rhs); return *this; }
const double &real_1d_array::operator()(ae_int_t i) const { return This->xR[i]; }
double &real_1d_array::operator()(ae_int_t i) { return This->xR[i]; }
const double &real_1d_array::operator[](ae_int_t i) const { return This->xR[i]; }
double &real_1d_array::operator[](ae_int_t i) { return This->xR[i]; }

void real_1d_array::setcontent(ae_int_t iLen, const double *pContent) {
// Handle possible exception-free errors.
   setlength(iLen);
// Copy, if its size matches.
   if (This != NULL && This->cnt == iLen)
      for (ae_int_t i = 0; i < iLen; i++) This->xR[i] = pContent[i];
}

const double *real_1d_array::getcontent() const { return This->xR; }
double *real_1d_array::getcontent() { return This->xR; }

// TODO: Convert to a constructor!
void real_1d_array::attach_to_ptr(ae_int_t iLen, double *pContent) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      owner = true, This = NULL, set_error_msg(); return;
#endif
   }
   alglib_impl::ae_assert(owner, "real_1d_array::attach_to_ptr: unable to attach proxy object to something else");
   alglib_impl::ae_assert(iLen > 0, "real_1d_array::attach_to_ptr: non-positive length");
   alglib_impl::x_vector x;
   x.cnt = iLen;
   x.datatype = alglib_impl::DT_REAL;
   x.owner = false;
   x.last_action = alglib_impl::ACT_UNCHANGED;
   x.x_ptr = pContent;
   attach_to(&x);
   alglib_impl::ae_state_clear();
}

complex_1d_array::complex_1d_array(): ae_vector_wrapper(alglib_impl::DT_COMPLEX) { }
complex_1d_array::complex_1d_array(alglib_impl::ae_vector *p): ae_vector_wrapper(p, alglib_impl::DT_COMPLEX) { }
complex_1d_array::complex_1d_array(const complex_1d_array &rhs): ae_vector_wrapper(rhs, alglib_impl::DT_COMPLEX) { }
#if !defined AE_NO_EXCEPTIONS
complex_1d_array::complex_1d_array(const char *s): ae_vector_wrapper(s, alglib_impl::DT_COMPLEX) { }
std::string complex_1d_array::tostring(int dps) const { return length() == 0 ? "[]" : arraytostring(&operator()(0), length(), dps); }
#endif
complex_1d_array::~complex_1d_array() { }

const complex_1d_array &complex_1d_array::operator=(const complex_1d_array &rhs) { assign(rhs); return *this; }
const complex &complex_1d_array::operator()(ae_int_t i) const { return *(const complex *)(This->xC + i); }
complex &complex_1d_array::operator()(ae_int_t i) { return *(complex *)(This->xC + i); }
const complex &complex_1d_array::operator[](ae_int_t i) const { return *(const complex *)(This->xC + i); }
complex &complex_1d_array::operator[](ae_int_t i) { return *(complex *)(This->xC + i); }

void complex_1d_array::setcontent(ae_int_t iLen, const complex *pContent) {
// Handle possible exception-free errors.
   setlength(iLen);
// Copy, if its size matches.
   if (This != NULL && This->cnt == iLen)
      for (ae_int_t i = 0; i < iLen; i++) This->xC[i] = complex_from_c(pContent[i]);
}

const complex *complex_1d_array::getcontent() const { return (const complex *)This->xC; }
complex *complex_1d_array::getcontent() { return (complex *)This->xC; }

ae_matrix_wrapper::ae_matrix_wrapper(alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      owner = true, This = NULL, set_error_msg(); return;
#endif
   }
   owner = true, This = &Obj, memset(This, 0, sizeof *This), ae_matrix_init(This, 0, 0, datatype, false);
   alglib_impl::ae_state_clear();
}

ae_matrix_wrapper::ae_matrix_wrapper(alglib_impl::ae_matrix *e_ptr, alglib_impl::ae_datatype datatype) {
   if (e_ptr->datatype != datatype) {
      const char *msg = "ae_matrix_wrapper::ae_matrix_wrapper: datatype check failed";
#if !defined AE_NO_EXCEPTIONS
      ThrowError(msg);
#else
      owner = true, This = NULL, set_error_flag(msg); return;
#endif
   }
   owner = false, This = e_ptr;
}

ae_matrix_wrapper::ae_matrix_wrapper(const ae_matrix_wrapper &rhs, alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      owner = true, This = NULL, set_error_msg(); return;
#endif
   }
   if (rhs.This != NULL)
      alglib_impl::ae_assert(rhs.This->datatype == datatype, "ae_matrix_wrapper::ae_matrix_wrapper: datatype check failed");
   owner = true, This = &Obj, memset(This, 0, sizeof *This);
   rhs.This == NULL ? ae_matrix_init(This, 0, 0, datatype, false) : ae_matrix_copy(This, rhs.This, false);
   alglib_impl::ae_state_clear();
}

ae_matrix_wrapper::~ae_matrix_wrapper() { if (This == &Obj) ae_matrix_free(This, false); }

// TODO: Automatic allocation of NULL pointer!
void ae_matrix_wrapper::setlength(ae_int_t rows, ae_int_t cols) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ae_assert(This != NULL, "ae_matrix_wrapper::setlength: p_mat == NULL (array was not correctly initialized)");
   alglib_impl::ae_assert(owner, "ae_matrix_wrapper::setlength: attempt to resize proxy array");
   alglib_impl::ae_matrix_set_length(This, rows, cols);
   alglib_impl::ae_state_clear();
}

ae_int_t ae_matrix_wrapper::cols() const { return This == NULL ? 0 : This->cols; }
ae_int_t ae_matrix_wrapper::rows() const { return This == NULL ? 0 : This->rows; }
ae_int_t ae_matrix_wrapper::getstride() const { return This == NULL ? 0 : This->stride; }
bool ae_matrix_wrapper::isempty() const { return cols() == 0 || rows() == 0; }

void ae_matrix_wrapper::attach_to(alglib_impl::x_matrix *new_ptr) {
   if (This == &Obj) ae_matrix_free(This, false);
   owner = false, This = &Obj, memset(This, 0, sizeof *This), ae_matrix_init_attach_to_x(This, new_ptr, false);
}

const ae_matrix_wrapper &ae_matrix_wrapper::assign(const ae_matrix_wrapper &rhs) {
   if (this == &rhs) return *this;
   alglib_impl::ae_state_init();
   TryCatch(*this)
   alglib_impl::ae_assert(This != NULL, "ae_matrix_wrapper::assign: incorrect assignment to matrix (uninitialized destination)");
// Assignment to proxy object.
   alglib_impl::ae_assert(rhs.This != NULL, "ae_matrix_wrapper::assign: incorrect assignment to array (uninitialized source)");
   alglib_impl::ae_assert(rhs.This->datatype == This->datatype, "ae_matrix_wrapper::assign: incorrect assignment to array (types dont match)");
   if (!owner) {
      alglib_impl::ae_assert(rhs.This->cols == This->cols, "ae_matrix_wrapper::assign: incorrect assignment to proxy array (sizes dont match)");
      alglib_impl::ae_assert(rhs.This->rows == This->rows, "ae_matrix_wrapper::assign: incorrect assignment to proxy array (sizes dont match)");
   }
   if (rhs.This->cols != This->cols || rhs.This->rows != This->rows) ae_matrix_set_length(This, rhs.This->rows, rhs.This->cols);
   for (ae_int_t i = 0; i < This->rows; i++) memcpy(This->xyX[i], rhs.This->xyX[i], This->cols * alglib_impl::ae_sizeof(This->datatype));
   alglib_impl::ae_state_clear();
   return *this;
}

#if !defined AE_NO_EXCEPTIONS
ae_matrix_wrapper::ae_matrix_wrapper(const char *s, alglib_impl::ae_datatype datatype) {
   char *p = filter_spaces(s); if (p == NULL) ThrowError("ae_matrix_wrapper::ae_matrix_wrapper: allocation error");
   try {
      std::vector< std::vector<const char *> > smat; str_matrix_create(p, &smat);
   {
      alglib_impl::ae_state_init();
      TryCatch()
      ae_int_t rows = smat.size(), cols = rows == 0 ? 0 : smat[0].size();
      owner = true, This = &Obj, memset(This, 0, sizeof *This), ae_matrix_init(This, rows, cols, datatype, false);
      alglib_impl::ae_state_clear();
   }
      for (size_t i = 0; i < smat.size(); i++) for (size_t j = 0; j < smat[0].size(); j++) switch (datatype) {
      // case alglib_impl::DT_BYTE: // The same as alglib_impl::DT_BOOL.
         case alglib_impl::DT_BOOL: This->xyB[i][j] = parse_bool_delim(smat[i][j], ",]"); break;
         case alglib_impl::DT_INT: This->xyZ[i][j] = parse_int_delim(smat[i][j], ",]"); break;
         case alglib_impl::DT_REAL: This->xyR[i][j] = parse_real_delim(smat[i][j], ",]"); break;
         case alglib_impl::DT_COMPLEX: This->xyC[i][j] = complex_from_c(parse_complex_delim(smat[i][j], ",]")); break;
      }
      alglib_impl::ae_free(p);
   } catch(...) {
      alglib_impl::ae_free(p);
      throw;
   }
}
#endif

boolean_2d_array::boolean_2d_array(): ae_matrix_wrapper(alglib_impl::DT_BOOL) { }
boolean_2d_array::boolean_2d_array(alglib_impl::ae_matrix *p): ae_matrix_wrapper(p, alglib_impl::DT_BOOL) { }
boolean_2d_array::boolean_2d_array(const boolean_2d_array &rhs): ae_matrix_wrapper(rhs, alglib_impl::DT_BOOL) { }
#if !defined AE_NO_EXCEPTIONS
boolean_2d_array::boolean_2d_array(const char *s): ae_matrix_wrapper(s, alglib_impl::DT_BOOL) { }
std::string boolean_2d_array::tostring() const {
   if (isempty()) return "[[]]";
   std::string result = "[";
   for (ae_int_t i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols());
   }
   return result += "]";
}
#endif
boolean_2d_array::~boolean_2d_array() { }

const boolean_2d_array &boolean_2d_array::operator=(const boolean_2d_array &rhs) { assign(rhs); return *this; }
const bool &boolean_2d_array::operator()(ae_int_t i, ae_int_t j) const { return This->xyB[i][j]; }
bool &boolean_2d_array::operator()(ae_int_t i, ae_int_t j) { return This->xyB[i][j]; }
const bool *boolean_2d_array::operator[](ae_int_t i) const { return This->xyB[i]; }
bool *boolean_2d_array::operator[](ae_int_t i) { return This->xyB[i]; }

void boolean_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const bool *pContent) {
// Handle possible exception-free errors.
   setlength(irows, icols);
// Copy, if its size matches.
   if (This != NULL && This->cols == icols && This->rows == irows)
      for (ae_int_t i = 0; i < irows; i++) for (ae_int_t j = 0; j < icols; j++) This->xyB[i][j] = pContent[i * icols + j];
}

integer_2d_array::integer_2d_array(): ae_matrix_wrapper(alglib_impl::DT_INT) { }
integer_2d_array::integer_2d_array(alglib_impl::ae_matrix *p): ae_matrix_wrapper(p, alglib_impl::DT_INT) { }
integer_2d_array::integer_2d_array(const integer_2d_array &rhs): ae_matrix_wrapper(rhs, alglib_impl::DT_INT) { }
#if !defined AE_NO_EXCEPTIONS
integer_2d_array::integer_2d_array(const char *s): ae_matrix_wrapper(s, alglib_impl::DT_INT) { }
std::string integer_2d_array::tostring() const {
   if (isempty()) return "[[]]";
   std::string result = "[";
   for (ae_int_t i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols());
   }
   return result += "]";
}
#endif
integer_2d_array::~integer_2d_array() { }

const integer_2d_array &integer_2d_array::operator=(const integer_2d_array &rhs) { assign(rhs); return *this; }
const ae_int_t &integer_2d_array::operator()(ae_int_t i, ae_int_t j) const { return This->xyZ[i][j]; }
ae_int_t &integer_2d_array::operator()(ae_int_t i, ae_int_t j) { return This->xyZ[i][j]; }
const ae_int_t *integer_2d_array::operator[](ae_int_t i) const { return This->xyZ[i]; }
ae_int_t *integer_2d_array::operator[](ae_int_t i) { return This->xyZ[i]; }

void integer_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const ae_int_t *pContent) {
// Handle possible exception-free errors.
   setlength(irows, icols);
// Copy, if its size matches.
   if (This != NULL && This->cols == icols && This->rows == irows)
      for (ae_int_t i = 0; i < irows; i++) for (ae_int_t j = 0; j < icols; j++) This->xyZ[i][j] = pContent[i * icols + j];
}

real_2d_array::real_2d_array(): ae_matrix_wrapper(alglib_impl::DT_REAL) { }
real_2d_array::real_2d_array(alglib_impl::ae_matrix *p): ae_matrix_wrapper(p, alglib_impl::DT_REAL) { }
real_2d_array::real_2d_array(const real_2d_array &rhs): ae_matrix_wrapper(rhs, alglib_impl::DT_REAL) { }
#if !defined AE_NO_EXCEPTIONS
real_2d_array::real_2d_array(const char *s): ae_matrix_wrapper(s, alglib_impl::DT_REAL) { }
std::string real_2d_array::tostring(int dps) const {
   if (isempty()) return "[[]]";
   std::string result = "[";
   for (ae_int_t i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols(), dps);
   }
   return result += "]";
}
#endif
real_2d_array::~real_2d_array() { }

const real_2d_array &real_2d_array::operator=(const real_2d_array &rhs) { assign(rhs); return *this; }
const double &real_2d_array::operator()(ae_int_t i, ae_int_t j) const { return This->xyR[i][j]; }
double &real_2d_array::operator()(ae_int_t i, ae_int_t j) { return This->xyR[i][j]; }
const double *real_2d_array::operator[](ae_int_t i) const { return This->xyR[i]; }
double *real_2d_array::operator[](ae_int_t i) { return This->xyR[i]; }

void real_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const double *pContent) {
// Handle possible exception-free errors.
   setlength(irows, icols);
// Copy, if its size matches.
   if (This != NULL && This->cols == icols && This->rows == irows)
      for (ae_int_t i = 0; i < irows; i++) for (ae_int_t j = 0; j < icols; j++) This->xyR[i][j] = pContent[i * icols + j];
}

void real_2d_array::attach_to_ptr(ae_int_t irows, ae_int_t icols, double *pContent) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      owner = true, This = NULL, set_error_msg(); return;
#endif
   }
   alglib_impl::ae_assert(owner, "real_2d_array::attach_to_ptr: unable to attach proxy object to something else");
   alglib_impl::ae_assert(icols > 0 && irows > 0, "real_2d_array::attach_to_ptr: non-positive length");
   alglib_impl::x_matrix x;
   x.cols = icols;
   x.rows = irows;
   x.stride = icols;
   x.datatype = alglib_impl::DT_REAL;
   x.owner = false;
   x.last_action = alglib_impl::ACT_UNCHANGED;
   x.x_ptr = pContent;
   attach_to(&x);
   alglib_impl::ae_state_clear();
}

complex_2d_array::complex_2d_array(): ae_matrix_wrapper(alglib_impl::DT_COMPLEX) { }
complex_2d_array::complex_2d_array(alglib_impl::ae_matrix *p): ae_matrix_wrapper(p, alglib_impl::DT_COMPLEX) { }
complex_2d_array::complex_2d_array(const complex_2d_array &rhs): ae_matrix_wrapper(rhs, alglib_impl::DT_COMPLEX) { }
#if !defined AE_NO_EXCEPTIONS
complex_2d_array::complex_2d_array(const char *s): ae_matrix_wrapper(s, alglib_impl::DT_COMPLEX) { }
std::string complex_2d_array::tostring(int dps) const {
   if (isempty()) return "[[]]";
   std::string result = "[";
   for (ae_int_t i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols(), dps);
   }
   return result += "]";
}
#endif
complex_2d_array::~complex_2d_array() { }

const complex_2d_array &complex_2d_array::operator=(const complex_2d_array &rhs) { assign(rhs); return *this; }
const complex &complex_2d_array::operator()(ae_int_t i, ae_int_t j) const { return *(const complex *)(This->xyC[i] + j); }
complex &complex_2d_array::operator()(ae_int_t i, ae_int_t j) { return *(complex *)(This->xyC[i] + j); }
const complex *complex_2d_array::operator[](ae_int_t i) const { return (const complex *)This->xyC[i]; }
complex *complex_2d_array::operator[](ae_int_t i) { return (complex *)This->xyC[i]; }

void complex_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const complex *pContent) {
// Handle possible exception-free errors.
   setlength(irows, icols);
// Copy, if its size matches.
   if (This != NULL && This->cols == icols && This->rows == irows)
      for (ae_int_t i = 0; i < irows; i++) for (ae_int_t j = 0; j < icols; j++)
         This->xyC[i][j] = complex_from_c(pContent[i * icols + j]);
}

#if !defined AE_NO_EXCEPTIONS
const int CSV_DEFAULT = 0x0, CSV_SKIP_HEADERS = 0x1;

// CSV function.
void read_csv(const char *filename, char separator, int flags, real_2d_array &out) {
// Parameters.
   bool skip_first_row = (flags & CSV_SKIP_HEADERS) != 0;
// Prepare an empty output array.
   out.setlength(0, 0);
// Open the file, determine its size and read its contents.
   FILE *f_in = fopen(filename, "rb");
   if (f_in == NULL) ThrowError("read_csv: unable to open input file");
   int flag = fseek(f_in, 0, SEEK_END);
   AE_CRITICAL_ASSERT(flag == 0);
   long int _filesize = ftell(f_in);
   AE_CRITICAL_ASSERT(_filesize >= 0);
   if (_filesize == 0) { // An empty file: return an empty array, success.
      fclose(f_in);
      return;
   }
   size_t filesize = _filesize;
   std::vector<char> v_buf;
   v_buf.resize(filesize + 2, 0);
   char *p_buf = &v_buf[0]; //(@) This is NOT equivalent to char *p_buf = v_buf!
   flag = fseek(f_in, 0, SEEK_SET);
   AE_CRITICAL_ASSERT(flag == 0);
   size_t bytes_read = fread((void *)p_buf, 1, filesize, f_in);
   AE_CRITICAL_ASSERT(bytes_read == filesize);
   fclose(f_in);
// Normalize the file contents:
// *	replace '\0' by ' ',
// *	remove the trailing ' ' and '\n' characters,
// *	append trailing '\n' and '\0' characters.
// Return if the file contains only ' ' and '\n' characters.
   for (size_t i = 0; i < filesize; i++) if (p_buf[i] == 0) p_buf[i] = ' ';
   for (; filesize > 0; filesize--) {
      char c = p_buf[filesize - 1];
      if (c != ' ' && c != '\t' && c != '\n' && c != '\r') break;
   }
   if (filesize == 0) return;
   p_buf[filesize++] = '\n';
   p_buf[filesize++] = '\0';
// Scan the dataset.
   size_t cols_count = 0, rows_count = 0, max_length = 0;
   std::vector<size_t> offsets, lengths;
   for (size_t row_start = 0; p_buf[row_start] != 0x0; rows_count++) {
   // The size of the row.
      size_t row_length = 0;
      for (; p_buf[row_start + row_length] != '\n'; row_length++);
   // Count the columns and perform an integrity check.
      size_t cur_cols_cnt = 1;
      for (size_t idx = 0; idx < row_length; idx++)
         if (p_buf[row_start + idx] == separator)
            cur_cols_cnt++;
      if (cols_count > 0 && cols_count != cur_cols_cnt)
         ThrowError("read_csv: non-rectangular contents, rows have different sizes");
      cols_count = cur_cols_cnt;
   // Store the offsets and lengths of the fields.
      size_t cur_offs = 0;
      for (size_t idx = 0; idx < row_length + 1; idx++)
         if (p_buf[row_start + idx] == separator || p_buf[row_start + idx] == '\n') {
            offsets.push_back(row_start + cur_offs);
            lengths.push_back(idx - cur_offs);
            max_length = idx - cur_offs > max_length ? idx - cur_offs : max_length;
            cur_offs = idx + 1;
         }
   // Advance row_start.
      row_start += row_length + 1;
   }
   AE_CRITICAL_ASSERT(rows_count >= 1);
   AE_CRITICAL_ASSERT(cols_count >= 1);
   AE_CRITICAL_ASSERT(cols_count * rows_count == offsets.size());
   AE_CRITICAL_ASSERT(cols_count * rows_count == lengths.size());
// Return on empty output.
   if (rows_count == 1 && skip_first_row) return;
// Convert.
   size_t row0 = skip_first_row ? 1 : 0;
   size_t row1 = rows_count;
   lconv *loc = localeconv();
   out.setlength(row1 - row0, cols_count);
   for (size_t ridx = row0; ridx < row1; ridx++)
      for (size_t cidx = 0; cidx < cols_count; cidx++) {
         char *p_field = p_buf + offsets[ridx * cols_count + cidx];
         size_t field_len = lengths[ridx * cols_count + cidx];
         for (size_t idx = 0; idx < field_len; idx++)
            if (p_field[idx] == '.' || p_field[idx] == ',')
               p_field[idx] = *loc->decimal_point;
         out[ridx - row0][cidx] = atof(p_field);
      }
}
#endif
} // end of namespace alglib
