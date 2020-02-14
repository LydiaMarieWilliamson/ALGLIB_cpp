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
// set _GNU_SOURCE flag BEFORE any #includes to get affinity management functions
#if AE_OS == AE_LINUX && !defined _GNU_SOURCE
#   define _GNU_SOURCE
#endif

#include "Ap.h"
#include <limits>
#include <locale.h>
#include <ctype.h>
#if defined AE_CPU && AE_CPU == AE_INTEL
#   if AE_COMPILER == AE_GNUC
#      include <fpu_control.h> // For _FPU_SETCW
#   elif AE_COMPILER == AE_MSVC
#      include <intrin.h>
#   endif
#endif

namespace alglib_impl {
// Core Code (Vectors, Matrices, Memory Management, etc.)
// OS-specific includes
#ifdef AE_USE_CPP
} // end of namespace alglib_impl
#endif

#if AE_OS == AE_POSIX
#   include <time.h>
#   include <unistd.h>
#   include <pthread.h>
#   include <sched.h>
#elif AE_OS == AE_WINDOWS
#   ifndef _WIN32_WINNT
#      define _WIN32_WINNT 0x0501
#   endif
#   include <windows.h>
#   include <process.h>
#endif

// Debugging helpers for Windows
#ifdef AE_DEBUG4WINDOWS
#   include <windows.h>
#endif

#ifdef AE_USE_CPP
namespace alglib_impl {
#endif
// local definitions
#define x_nb 16
#define AE_DATA_ALIGN 64
#define AE_PTR_ALIGN sizeof(void*)
#define AE_LITTLE_ENDIAN 1
#define AE_BIG_ENDIAN 2
#define AE_MIXED_ENDIAN 3
#define AE_SER_ENTRY_LENGTH 11
#define AE_SER_ENTRIES_PER_ROW 5
#define AE_SM_DEFAULT 0
#define AE_SM_ALLOC 1
#define AE_SM_READY2S 2
#define AE_SM_TO_STRING    10
#define AE_SM_TO_CPPSTRING 11
#define AE_SM_TO_STREAM    12
#define AE_SM_FROM_STRING  20
#define AE_SM_FROM_STREAM  22
#define AE_LOCK_CYCLES 512
#define AE_LOCK_TESTS_BEFORE_YIELD 16
#define AE_CRITICAL_ASSERT(x) if (!(x)) abort()

// IDs for set_dbg_value
#define _ALGLIB_USE_ALLOC_COUNTER             0
#define _ALGLIB_USE_DBG_COUNTERS              1
#define _ALGLIB_USE_VENDOR_KERNELS          100
#define _ALGLIB_VENDOR_MEMSTAT              101
#define _ALGLIB_DEBUG_WORKSTEALING          200
#define _ALGLIB_WSDBG_NCORES                201
#define _ALGLIB_WSDBG_PUSHROOT_OK           202
#define _ALGLIB_WSDBG_PUSHROOT_FAILED       203
#define _ALGLIB_SET_GLOBAL_THREADING       1001
#define _ALGLIB_SET_NWORKERS               1002

// IDs for get_dbg_value
#define _ALGLIB_GET_ALLOC_COUNTER             0
#define _ALGLIB_GET_CUMULATIVE_ALLOC_SIZE     1
#define _ALGLIB_GET_CUMULATIVE_ALLOC_COUNT    2
#define _ALGLIB_GET_CORES_COUNT            1000
#define _ALGLIB_GET_GLOBAL_THREADING       1001
#define _ALGLIB_GET_NWORKERS               1002

// Lock.
//
// This is internal structure which implements lock functionality.
typedef struct {
#if AE_OS == AE_POSIX
   pthread_mutex_t mutex;
#elif AE_OS == AE_WINDOWS
   volatile ae_int_t *volatile p_lock;
   char buf[sizeof(ae_int_t) + AE_LOCK_ALIGNMENT];
#else
   bool is_locked;
#endif
} _lock;

// Global flags, split into several char-sized variables in order
// to avoid problem with non-atomic reads/writes (single-byte ops
// are atomic on all modern architectures);
//
// Following variables are included:
// * threading-related settings
unsigned char _alglib_global_threading_flags = _ALGLIB_FLG_THREADING_SERIAL >> _ALGLIB_FLG_THREADING_SHIFT;

// DESCRIPTION: recommended number of active workers:
//              * positive value >= 1 is used to specify exact number of active workers
//              * 0 means that ALL available cores are used
//              * negative value means that all cores EXCEPT for cores_to_use will be used
//                (say, -1 means that all cores except for one will be used). At least one
//                core will be used in this case, even if you assign -9999999 to this field.
//
//              Default value =  0 (fully parallel execution) when AE_NWORKERS is not defined
//                            =  0 for manually defined number of cores (AE_NWORKERS is defined)
// PROTECTION:  not needed; runtime modification is possible, but we do not need exact
//              synchronization.
#if defined AE_NWORKERS
#   if AE_NWORKERS <= 0
#      error AE_NWORKERS must be positive number or not defined at all.
#   endif
ae_int_t _alglib_cores_to_use = 0;
#else
ae_int_t _alglib_cores_to_use = 0;
#endif

// Debug counters
ae_int_t _alloc_counter = 0;
ae_int_t _alloc_counter_total = 0;
bool _use_alloc_counter = false;
ae_int_t _dbg_alloc_total = 0;
bool _use_dbg_counters = false;
bool _use_vendor_kernels = true;
bool debug_workstealing = false; // Debug workstealing environment? False by default.
ae_int_t dbgws_pushroot_ok = 0;
ae_int_t dbgws_pushroot_failed = 0;

#ifdef AE_SMP_DEBUGCOUNTERS
__declspec(align(AE_LOCK_ALIGNMENT)) volatile ae_int64_t _ae_dbg_lock_acquisitions = 0;
__declspec(align(AE_LOCK_ALIGNMENT)) volatile ae_int64_t _ae_dbg_lock_spinwaits = 0;
__declspec(align(AE_LOCK_ALIGNMENT)) volatile ae_int64_t _ae_dbg_lock_yields = 0;
#endif

// Allocation debugging
bool _force_malloc_failure = false;
ae_int_t _malloc_failure_after = 0;

// Fields for memory allocation over static array
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
#   if AE_THREADING != AE_SERIAL_UNSAFE
#      error Basis static malloc is thread-unsafe; define AE_THREADING=AE_SERIAL_UNSAFE to prove that you know it
#   endif
static ae_int_t sm_page_size = 0;
static ae_int_t sm_page_cnt = 0;
static ae_int_t *sm_page_tbl = NULL;
static unsigned char *sm_mem = NULL;
#endif

// These declarations are used to ensure that
// sizeof(bool)=1, sizeof(ae_int32_t) == 4, sizeof(ae_int64_t) == 8, sizeof(ae_int_t) == sizeof(void*).
// they will lead to syntax error otherwise (array size will be negative).
//
// you can remove them, if you want - they are not used anywhere.
//
static char _ae_bool_must_be_8_bits_wide[1 - 2 *((int)(sizeof(bool)) - 1) *((int)(sizeof(bool)) - 1)];
static char _ae_int32_t_must_be_32_bits_wide[1 - 2 *((int)(sizeof(ae_int32_t)) - 4) *((int)(sizeof(ae_int32_t)) - 4)];
static char _ae_int64_t_must_be_64_bits_wide[1 - 2 *((int)(sizeof(ae_int64_t)) - 8) *((int)(sizeof(ae_int64_t)) - 8)];
static char _ae_uint64_t_must_be_64_bits_wide[1 - 2 *((int)(sizeof(ae_uint64_t)) - 8) *((int)(sizeof(ae_uint64_t)) - 8)];
static char _ae_int_t_must_be_pointer_sized[1 - 2 *((int)(sizeof(ae_int_t)) - (int)sizeof(void *)) *((int)(sizeof(ae_int_t)) - (int)(sizeof(void *)))];

// This variable is used to prevent some tricky optimizations which may degrade multithreaded performance.
// It is touched once in the ae_init_pool() function from smp.c in order to prevent optimizations.
static volatile ae_int_t ae_never_change_it = 1;

// This function should never  be  called.  It is  here  to  prevent spurious
// compiler warnings about unused variables (in fact: used).
void ae_never_call_it() {
}

void ae_set_dbg_flag(ae_int64_t flag_id, ae_int64_t flag_val) {
   if (flag_id == _ALGLIB_USE_ALLOC_COUNTER) {
      _use_alloc_counter = flag_val != 0;
      return;
   }
   if (flag_id == _ALGLIB_USE_DBG_COUNTERS) {
      _use_dbg_counters = flag_val != 0;
      return;
   }
   if (flag_id == _ALGLIB_USE_VENDOR_KERNELS) {
      _use_vendor_kernels = flag_val != 0;
      return;
   }
   if (flag_id == _ALGLIB_DEBUG_WORKSTEALING) {
      debug_workstealing = flag_val != 0;
      return;
   }
   if (flag_id == _ALGLIB_SET_GLOBAL_THREADING) {
      ae_set_global_threading((ae_uint64_t) flag_val);
      return;
   }
   if (flag_id == _ALGLIB_SET_NWORKERS) {
      _alglib_cores_to_use = (ae_int_t) flag_val;
      return;
   }
}

ae_int64_t ae_get_dbg_value(ae_int64_t id) {
   if (id == _ALGLIB_GET_ALLOC_COUNTER)
      return _alloc_counter;
   if (id == _ALGLIB_GET_CUMULATIVE_ALLOC_SIZE)
      return _dbg_alloc_total;
   if (id == _ALGLIB_GET_CUMULATIVE_ALLOC_COUNT)
      return _alloc_counter_total;
   if (id == _ALGLIB_VENDOR_MEMSTAT) {
#if defined AE_MKL
      return ae_mkl_memstat();
#else
      return 0;
#endif
   }
// workstealing counters
   if (id == _ALGLIB_WSDBG_NCORES)
#if defined AE_SMP
      return ae_cores_count();
#else
      return 0;
#endif
   if (id == _ALGLIB_WSDBG_PUSHROOT_OK)
      return dbgws_pushroot_ok;
   if (id == _ALGLIB_WSDBG_PUSHROOT_FAILED)
      return dbgws_pushroot_failed;
   if (id == _ALGLIB_GET_CORES_COUNT)
#if defined AE_SMP
      return ae_cores_count();
#else
      return 0;
#endif
   if (id == _ALGLIB_GET_GLOBAL_THREADING)
      return (ae_int64_t) ae_get_global_threading();
   if (id == _ALGLIB_GET_NWORKERS)
      return (ae_int64_t) _alglib_cores_to_use;
// unknown value
   return 0;
}

// This function sets default (global) threading model:
// * serial execution
// * multithreading, if cores_to_use allows it
//
void ae_set_global_threading(ae_uint64_t flg_value) {
   flg_value &= _ALGLIB_FLG_THREADING_MASK;
   AE_CRITICAL_ASSERT(flg_value == _ALGLIB_FLG_THREADING_SERIAL || flg_value == _ALGLIB_FLG_THREADING_PARALLEL);
   _alglib_global_threading_flags = (unsigned char)(flg_value >> _ALGLIB_FLG_THREADING_SHIFT);
}

// This function gets default (global) threading model:
// * serial execution
// * multithreading, if cores_to_use allows it
//
ae_uint64_t ae_get_global_threading() {
   return ((ae_uint64_t) _alglib_global_threading_flags) << _ALGLIB_FLG_THREADING_SHIFT;
}

ae_int_t ae_misalignment(const void *ptr, size_t alignment) {
   union _u {
      const void *ptr;
      ae_int_t iptr;
   } u;
   u.ptr = ptr;
   return (ae_int_t)(u.iptr % alignment);
}

void *ae_align(void *ptr, size_t alignment) {
   char *result = (char *)ptr;
   if ((result - (char *)0) % alignment != 0)
      result += alignment - (result - (char *)0) % alignment;
   return result;
}

// This function maps nworkers  number  (which  can  be  positive,  zero  or
// negative with 0 meaning "all cores", -1 meaning "all cores -1" and so on)
// to "effective", strictly positive workers count.
//
// This  function  is  intended  to  be used by debugging/testing code which
// tests different number of worker threads. It is NOT aligned  in  any  way
// with ALGLIB multithreading framework (i.e. it can return  non-zero worker
// count even for single-threaded GPLed ALGLIB).
ae_int_t ae_get_effective_workers(ae_int_t nworkers) {
   ae_int_t ncores;
// determine cores count
#if defined AE_NWORKERS
   ncores = AE_NWORKERS;
#elif AE_OS == AE_POSIX
   {
      long r = sysconf(_SC_NPROCESSORS_ONLN);
      ncores = r <= 0 ? 1 : r;
   }
#elif AE_OS == AE_WINDOWS
   SYSTEM_INFO sysInfo;
   GetSystemInfo(&sysInfo);
   ncores = (ae_int_t)(sysInfo.dwNumberOfProcessors);
#else
   ncores = 1;
#endif
   AE_CRITICAL_ASSERT(ncores >= 1);
// map nworkers to its effective value
   if (nworkers >= 1)
      return nworkers > ncores ? ncores : nworkers;
   return ncores + nworkers >= 1 ? ncores + nworkers : 1;
}

// This function belongs to the family of  "optional  atomics",  i.e.  atomic
// functions which either perform atomic changes - or do nothing at  all,  if
// current compiler settings do not allow us to generate atomic code.
//
// All "optional atomics" are synchronized, i.e. either all of them work - or
// no one of the works.
//
// This particular function performs atomic addition on pointer-sized  value,
// which must be pointer-size aligned.
//
// NOTE: this function is not intended to be extremely high performance  one,
//       so use it only when necessary.
void ae_optional_atomic_add_i(ae_int_t *p, ae_int_t v) {
   AE_CRITICAL_ASSERT(ae_misalignment(p, sizeof(void *)) == 0);
#if AE_OS == AE_WINDOWS
   while (true) {
   // perform conversion between ae_int_t* and void** without compiler warnings about indirection levels.
      union _u {
         PVOID volatile *volatile ptr;
         volatile ae_int_t *volatile iptr;
      } u;
      u.iptr = p;
   // atomic read for initial value
      PVOID v0 = InterlockedCompareExchangePointer(u.ptr, NULL, NULL);
   // increment cached value and store
      if (InterlockedCompareExchangePointer(u.ptr, (PVOID)(((char *)v0) + v), v0) == v0)
         break;
   }
#elif AE_COMPILER == AE_GNUC && AE_CPU == AE_INTEL && 100*__GNUC__ + __GNUC__ >= 470
   __atomic_add_fetch(p, v, __ATOMIC_RELAXED);
#else
   *p += v; // At least do something for older compilers!
#endif
}

// This function belongs to the family of  "optional  atomics",  i.e.  atomic
// functions which either perform atomic changes - or do nothing at  all,  if
// current compiler settings do not allow us to generate atomic code.
//
// All "optional atomics" are synchronized, i.e. either all of them work - or
// no one of the works.
//
// This  particular  function  performs  atomic  subtraction on pointer-sized
// value, which must be pointer-size aligned.
//
// NOTE: this function is not intended to be extremely high performance  one,
//       so use it only when necessary.
void ae_optional_atomic_sub_i(ae_int_t *p, ae_int_t v) {
   AE_CRITICAL_ASSERT(ae_misalignment(p, sizeof(void *)) == 0);
#if AE_OS == AE_WINDOWS
   while (true) {
   // perform conversion between ae_int_t* and void** without compiler warnings about indirection levels.
      union _u {
         PVOID volatile *volatile ptr;
         volatile ae_int_t *volatile iptr;
      } u;
      u.iptr = p;
   // atomic read for initial value, convert it to 1-byte pointer
      PVOID v0 = InterlockedCompareExchangePointer(u.ptr, NULL, NULL);
   // increment cached value and store
      if (InterlockedCompareExchangePointer(u.ptr, (PVOID)(((char *)v0) - v), v0) == v0)
         break;
   }
#elif AE_COMPILER == AE_GNUC && AE_CPU == AE_INTEL && 100*__GNUC__ + __GNUC__ >= 470
   __atomic_sub_fetch(p, v, __ATOMIC_RELAXED);
#else
   *p -= v; // At least do something for older compilers!
#endif
}

// jmp_buf pointer for internal C-style exception handling
AutoS jmp_buf *volatile CurBreakAt;

// This function sets jump buffer for error handling.
// At may be NULL.
void ae_state_set_break_jump(jmp_buf *At) { CurBreakAt = At; }

// ae_error_type of the last error (filled when exception is thrown)
AutoS ae_error_type volatile CurStatus;

// human-readable message (filled when exception is thrown)
AutoS const char *volatile CurMsg;

// pointer to the top block in a stack of frames which hold dynamically allocated objects.
AutoS ae_frame *volatile TopFr;

// This function abnormally aborts program, using one of several ways:
//
// * for TopFr != NULL and CurBreakAt being initialized with  call  to
//   ae_state_set_break_jump() - it performs longjmp() to return site.
// * otherwise, abort() is called
// In all cases, for  TopFr != NULL function sets CurStatus and CurMsg fields.
// It also clears TopFr with ae_state_clear().
void ae_break(ae_error_type error_type, const char *msg) {
   if (TopFr != NULL) {
      ae_state_clear();
      CurStatus = error_type;
      CurMsg = msg;
      if (CurBreakAt != NULL) longjmp(*CurBreakAt, 1); else abort();
   } else abort();
}

// Assertion
// For TopFr != NULL it allows to gracefully leave ALGLIB session,
// removing all frames and deallocating registered dynamic data structure.
// For TopFr == NULL it just abort()'s program.
// IMPORTANT: this function ALWAYS evaluates its argument.  It  can  not  be
//            replaced by macro which does nothing. So, you may place actual
//            function calls at cond, and these will always be performed.
void ae_assert(bool cond, const char *msg) {
   if (!cond) ae_break(ERR_ASSERTION_FAILED, msg);
}

void ae_impose(bool cond, const char *msg) {
   if (!cond) fprintf(stderr, msg), abort();
}

#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
void set_memory_pool(void *ptr, size_t size) {
// Integrity checks
   AE_CRITICAL_ASSERT(sm_page_size == 0);
   AE_CRITICAL_ASSERT(sm_page_cnt == 0);
   AE_CRITICAL_ASSERT(sm_page_tbl == NULL);
   AE_CRITICAL_ASSERT(sm_mem == NULL);
   AE_CRITICAL_ASSERT(size > 0);
// Align pointer
   size -= ae_misalignment(ptr, sizeof(ae_int_t));
   ptr = ae_align(ptr, sizeof(ae_int_t));
// Calculate page size and page count, prepare pointers to page table and memory
   sm_page_size = 256;
   AE_CRITICAL_ASSERT(size >= (sm_page_size + sizeof(ae_int_t)) + sm_page_size);        // We expect to have memory for at least one page + table entry + alignment.
   sm_page_cnt = (size - sm_page_size) / (sm_page_size + sizeof(ae_int_t));
   AE_CRITICAL_ASSERT(sm_page_cnt > 0);
   sm_page_tbl = (ae_int_t *) ptr;
   sm_mem = (unsigned char *)ae_align(sm_page_tbl + sm_page_cnt, sm_page_size);
// Mark all pages as free
   memset(sm_page_tbl, 0, sm_page_cnt * sizeof(ae_int_t));
}

void *ae_static_malloc(size_t size, size_t alignment) {
   int rq_pages, i, j, cur_len;
   AE_CRITICAL_ASSERT(size >= 0);
   AE_CRITICAL_ASSERT(sm_page_size > 0);
   AE_CRITICAL_ASSERT(sm_page_cnt > 0);
   AE_CRITICAL_ASSERT(sm_page_tbl != NULL);
   AE_CRITICAL_ASSERT(sm_mem != NULL);
   if (size == 0)
      return NULL;
   if (_force_malloc_failure)
      return NULL;
// check that page alignment and requested alignment match each other
   AE_CRITICAL_ASSERT(alignment <= sm_page_size);
   AE_CRITICAL_ASSERT((sm_page_size % alignment) == 0);
// search long enough sequence of pages
   rq_pages = size / sm_page_size;
   if (size % sm_page_size)
      rq_pages++;
   cur_len = 0;
   for (i = 0; i < sm_page_cnt;) {
   // determine length of the sequence of free pages
      if (sm_page_tbl[i] == 0)
         cur_len++;
      else {
         AE_CRITICAL_ASSERT(sm_page_tbl[i] > 0);
         cur_len = 0;
         i += sm_page_tbl[i];
         continue;
      }
   // found it?
      if (cur_len >= rq_pages) {
      // update counters (if flag is set)
         if (_use_alloc_counter) {
            ae_optional_atomic_add_i(&_alloc_counter, 1);
            ae_optional_atomic_add_i(&_alloc_counter_total, 1);
         }
         if (_use_dbg_counters)
            ae_optional_atomic_add_i(&_dbg_alloc_total, size);
      // mark pages and return
         for (j = 0; j < rq_pages; j++)
            sm_page_tbl[i - j] = -1;
         sm_page_tbl[i - (rq_pages - 1)] = rq_pages;
         return sm_mem + (i - (rq_pages - 1)) * sm_page_size;
      }
   // next element
      i++;
   }
   return NULL;
}

void ae_static_free(void *block) {
   ae_int_t page_idx, page_cnt, i;
   if (block == NULL)
      return;
   page_idx = (unsigned char *)block - sm_mem;
   AE_CRITICAL_ASSERT(page_idx >= 0);
   AE_CRITICAL_ASSERT((page_idx % sm_page_size) == 0);
   page_idx /= sm_page_size;
   AE_CRITICAL_ASSERT(page_idx < sm_page_cnt);
   page_cnt = sm_page_tbl[page_idx];
   AE_CRITICAL_ASSERT(page_cnt >= 1);
   for (i = 0; i < page_cnt; i++)
      sm_page_tbl[page_idx + i] = 0;
// update counters (if flag is set)
   if (_use_alloc_counter)
      ae_optional_atomic_sub_i(&_alloc_counter, 1);
}

void memory_pool_stats(ae_int_t *bytes_used, ae_int_t *bytes_free) {
   int i;
   AE_CRITICAL_ASSERT(sm_page_size > 0);
   AE_CRITICAL_ASSERT(sm_page_cnt > 0);
   AE_CRITICAL_ASSERT(sm_page_tbl != NULL);
   AE_CRITICAL_ASSERT(sm_mem != NULL);
// scan page table
   *bytes_used = 0;
   *bytes_free = 0;
   for (i = 0; i < sm_page_cnt;) {
      if (sm_page_tbl[i] == 0) {
         (*bytes_free)++;
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
   char *result = NULL;
   if (size == 0)
      return NULL;
   if (_force_malloc_failure)
      return NULL;
   if (_malloc_failure_after > 0 && _alloc_counter_total >= _malloc_failure_after)
      return NULL;
// allocate
   if (alignment <= 1) {
   // no alignment, just call alloc
      void *block;
      void **p;;
      block = malloc(sizeof(void *) + size);
      if (block == NULL)
         return NULL;
      p = (void **)block;
      *p = block;
      result = (char *)((char *)block + sizeof(void *));
   } else {
   // align
      void *block;
      block = malloc(alignment - 1 + sizeof(void *) + size);
      if (block == NULL)
         return NULL;
      result = (char *)block + sizeof(void *);
#   if 0
      if ((result - (char *)0) % alignment != 0)
         result += alignment - (result - (char *)0) % alignment;
#   endif
      result = (char *)ae_align(result, alignment);
      *((void **)(result - sizeof(void *))) = block;
   }
// update counters (if flag is set)
   if (_use_alloc_counter) {
      ae_optional_atomic_add_i(&_alloc_counter, 1);
      ae_optional_atomic_add_i(&_alloc_counter_total, 1);
   }
   if (_use_dbg_counters)
      ae_optional_atomic_add_i(&_dbg_alloc_total, (ae_int64_t) size);
// return
   return (void *)result;
#endif
}

void *aligned_extract_ptr(void *block) {
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
   return NULL;
#else
   if (block == NULL)
      return NULL;
   return *((void **)((char *)block - sizeof(void *)));
#endif
}

void aligned_free(void *block) {
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
   ae_static_free(block);
#else
   void *p;
   if (block == NULL)
      return;
   p = aligned_extract_ptr(block);
   free(p);
   if (_use_alloc_counter)
      ae_optional_atomic_sub_i(&_alloc_counter, 1);
#endif
}

void *eternal_malloc(size_t size) {
   if (size == 0)
      return NULL;
   if (_force_malloc_failure)
      return NULL;
   return malloc(size);
}

// Allocate memory with automatic alignment.
//
// Returns NULL when zero size is specified.
//
// Error handling:
// * if TopFr == NULL, returns NULL on allocation error
// * if TopFr != NULL, calls ae_break() on allocation error
void *ae_malloc(size_t size) {
   void *result;
   if (size == 0)
      return NULL;
   result = aligned_malloc(size, AE_DATA_ALIGN);
   if (result == NULL && TopFr != NULL)
      ae_break(ERR_OUT_OF_MEMORY, "ae_malloc(): out of memory");
   return result;
}

void ae_free(void *p) {
   if (p != NULL)
      aligned_free(p);
}

// Sets pointers to the matrix rows.
//
// * dst must be correctly initialized matrix
// * dst->data.ptr points to the beginning of memory block  allocated  for
//   row pointers.
// * dst->ptr - undefined (initialized during algorithm processing)
// * storage parameter points to the beginning of actual storage
void ae_matrix_update_row_pointers(ae_matrix *dst, void *storage) {
   char *p_base;
   void **pp_ptr;
   ae_int_t i;
   if (dst->rows > 0 && dst->cols > 0) {
      p_base = (char *)storage;
      pp_ptr = (void **)dst->data.ptr;
      dst->ptr.pp_void = pp_ptr;
      for (i = 0; i < dst->rows; i++, p_base += dst->stride * ae_sizeof(dst->datatype))
         pp_ptr[i] = p_base;
   } else
      dst->ptr.pp_void = NULL;
}

// Returns size of datatype.
// Zero for dynamic types like strings or multiple precision types.
ae_int_t ae_sizeof(ae_datatype datatype) {
   switch (datatype) {
      case DT_BOOL:
         return (ae_int_t) sizeof(bool);
      case DT_INT:
         return (ae_int_t) sizeof(ae_int_t);
      case DT_REAL:
         return (ae_int_t) sizeof(double);
      case DT_COMPLEX:
         return 2 * (ae_int_t) sizeof(double);
      default:
         return 0;
   }
}

// Checks that n bytes pointed by ptr are zero.
//
// This function is used in the constructors to check that  instance  fields
// on entry are correctly initialized by zeros.
bool ae_check_zeros(const void *ptr, ae_int_t n) {
   ae_int_t nu, nr, i;
   unsigned long long c = 0x0;
// determine leading and trailing lengths
   nu = n / sizeof(unsigned long long);
   nr = n % sizeof(unsigned long long);
// handle leading nu long long elements
   if (nu > 0) {
      const unsigned long long *p_ull;
      p_ull = (const unsigned long long *)ptr;
      for (i = 0; i < nu; i++)
         c |= p_ull[i];
   }
// handle trailing nr char elements
   if (nr > 0) {
      const unsigned char *p_uc;
      p_uc = ((const unsigned char *)ptr) + nu * sizeof(unsigned long long);
      for (i = 0; i < nr; i++)
         c |= p_uc[i];
   }
// done
   return c == 0x0;
}

// Flags: call-local settings for ALGLIB
AutoS ae_uint64_t CurFlags;

// This function sets current state's flags.
// buf may be NULL.
void ae_state_set_flags(ae_uint64_t flags) { CurFlags = flags; }

unsigned char DynBottom = 1, DynFrame = 2;

// This function makes new stack frame.
//
// This function takes two parameters: environment state and pointer to  the
// dynamic block which will be used as indicator  of  the  frame  beginning.
// This dynamic block must be initialized by caller and mustn't  be changed/
// deallocated/reused till ae_leave_frame called. It may be global or  local
// variable (local is even better).
void ae_frame_make(ae_frame *Fr) { Fr->p_next = TopFr, Fr->deallocator = NULL, Fr->ptr = &DynFrame, TopFr = Fr; }

// This function leaves current stack frame and deallocates all automatic
// dynamic blocks which were attached to this frame.
void ae_frame_leave() {
   for (; TopFr->ptr != &DynFrame && TopFr->ptr != &DynBottom; TopFr = TopFr->p_next)
      if (TopFr->ptr != NULL && TopFr->deallocator != NULL) TopFr->deallocator(TopFr->ptr);
   if (TopFr->ptr == &DynFrame) TopFr = TopFr->p_next;
}

// This function initializes ALGLIB environment state.
// NOTES:
// * stacks contain no frames, so ae_make_frame() must be called before
//   attaching dynamic blocks. Without it ae_leave_frame() will cycle
//   forever (which is intended behavior).
void ae_state_init() {
// The base of the current stack of frames which hold dynamically allocated objects.
// p_next points to itself because:
// * correct program should be able to detect end of the list
//   by looking at the ptr field.
// * NULL p_next may be used to distinguish automatic blocks
//   (in the list) from non-automatic (not in the list)
   static ae_frame BotFr = { &BotFr, NULL, &DynBottom };
// Set the status indicators and clear the frame.
   CurFlags = 0x0, CurBreakAt = NULL, CurMsg = "", TopFr = &BotFr;
}

// This function clears ALGLIB environment state.
// All dynamic data controlled by state are freed.
void ae_state_clear() {
   if (TopFr == NULL) return;
   for (; TopFr->ptr != &DynBottom; TopFr = TopFr->p_next)
      if (TopFr->ptr != NULL && TopFr->deallocator != NULL) TopFr->deallocator(TopFr->ptr);
   TopFr = NULL;
}

// This function attaches block to the dynamic block list of the current state.
// block               block
// This function does NOT generate exceptions.
// NOTES:
// * never call it for special blocks which marks frame boundaries!
void ae_db_attach(ae_dyn_block *block) { block->p_next = TopFr, TopFr = block; }

// This function initializes dynamic block for the current state:
// block               destination block, MUST be zero-filled on entry
// size                size (in bytes), >= 0.
// make_automatic                if true, vector is added to the dynamic block list
//
// block is assumed to be uninitialized, its fields are ignored. You may
// call this function with zero size in order to register block in the
// dynamic list.
//
// Error handling: calls ae_break() on allocation error. Block is left in
// valid state (empty, but valid).
//
// NOTES:
// * never call it for blocks which are already in the list; use ae_db_realloc
//   for already allocated blocks.
//
// NOTE: no memory allocation is performed for initialization with size=0
void ae_db_init(ae_dyn_block *block, ae_int_t size, bool make_automatic) {
   AE_CRITICAL_ASSERT(ae_check_zeros(block, sizeof(*block)));
// NOTE: these strange dances around block->ptr are necessary
//       in order to correctly handle possible exceptions during
//       memory allocation.
   ae_assert(size >= 0, "ae_db_init(): negative size");
   block->ptr = NULL;
   if (make_automatic) ae_db_attach(block); else block->p_next = NULL;
   if (size != 0) block->ptr = ae_malloc((size_t)size);
   block->deallocator = ae_free;
}

// This function realloc's dynamic block of the current state:
//
// block               destination block (initialized)
// size                new size (in bytes)
//
// block is assumed to be initialized.
//
// This function:
// * deletes old contents
// * preserves automatic state
//
// Error handling: calls ae_break() on allocation error. Block is left in
// valid state - empty, but valid.
//
// NOTES:
// * never call it for special blocks which mark frame boundaries!
void ae_db_realloc(ae_dyn_block *block, ae_int_t size) {
// NOTE: these strange dances around block->ptr are necessary
//       in order to correctly handle possible exceptions during
//       memory allocation.
   ae_assert(size >= 0, "ae_db_realloc(): negative size");
   if (block->ptr != NULL) ((ae_deallocator) block->deallocator)(block->ptr), block->ptr = NULL;
   block->ptr = ae_malloc((size_t)size);
   block->deallocator = ae_free;
}

// This function clears dynamic block (releases  all  dynamically  allocated
// memory). Dynamic block may be in automatic management list - in this case
// it will NOT be removed from list.
//
// block               destination block (initialized)
//
// NOTES:
// * never call it for special blocks which marks frame boundaries!
void ae_db_free(ae_dyn_block *block) {
   if (block->ptr != NULL) ((ae_deallocator) block->deallocator)(block->ptr), block->ptr = NULL;
   block->deallocator = ae_free;
}

// This function swaps contents of two dynamic blocks (pointers and
// deallocators) leaving other parameters (automatic management settings,
// etc.) unchanged.
//
// NOTES:
// * never call it for special blocks which marks frame boundaries!
void ae_db_swap(ae_dyn_block *block1, ae_dyn_block *block2) {
   void (*deallocator)(void *) = NULL;
   void *volatile ptr;
   ptr = block1->ptr;
   deallocator = block1->deallocator;
   block1->ptr = block2->ptr;
   block1->deallocator = block2->deallocator;
   block2->ptr = ptr;
   block2->deallocator = deallocator;
}

// This function creates ae_vector.
// Vector size may be zero. Vector contents is uninitialized.
//
// dst                 destination vector, MUST be zero-filled (we  check  it
//                     and call abort() if *dst is non-zero; the rationale is
//                     that we can not correctly handle errors in constructors
//                     without zero-filling).
// size                vector size, may be zero
// datatype            guess what...
// make_automatic                if true, vector will be registered in the frame of the current state;
//
// NOTE: no memory allocation is performed for initialization with size=0
void ae_vector_init(ae_vector *dst, ae_int_t size, ae_datatype datatype, bool make_automatic) {
// Integrity checks
   AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof(*dst)));
   ae_assert(size >= 0, "ae_vector_init(): negative size");
// prepare for possible errors during allocation
   dst->cnt = 0;
   dst->ptr.p_ptr = NULL;
// init
   ae_db_init(&dst->data, size * ae_sizeof(datatype), make_automatic);
   dst->cnt = size;
   dst->datatype = datatype;
   dst->ptr.p_ptr = dst->data.ptr;
   dst->is_attached = false;
}

// This function creates copy of ae_vector. New copy of the data is created,
// which is managed and owned by newly initialized vector.
//
// dst                 destination vector, MUST be zero-filled (we  check  it
//                     and call abort() if *dst is non-zero; the rationale is
//                     that we can not correctly handle errors in constructors
//                     without zero-filling).
// src                 well, it is source
// make_automatic                if true, vector will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_vector_copy(ae_vector *dst, ae_vector *src, bool make_automatic) {
   ae_vector_init(dst, src->cnt, src->datatype, make_automatic);
   if (src->cnt != 0) memmove(dst->ptr.p_ptr, src->ptr.p_ptr, (size_t)(src->cnt * ae_sizeof(src->datatype)));
}

// This function initializes ae_vector using X-structure as source. New copy
// of data is created, which is owned/managed by ae_vector  structure.  Both
// structures (source and destination) remain completely  independent  after
// this call.
//
// dst                 destination vector, MUST be zero-filled (we  check  it
//                     and call abort() if *dst is non-zero; the rationale is
//                     that we can not correctly handle errors in constructors
//                     without zero-filling).
// src                 well, it is source
// make_automatic                if true, vector will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_vector_init_from_x(ae_vector *dst, x_vector *src, bool make_automatic) {
   ae_vector_init(dst, (ae_int_t) src->cnt, (ae_datatype) src->datatype, make_automatic);
   if (src->cnt > 0) memmove(dst->ptr.p_ptr, src->x_ptr.p_ptr, (size_t)(((ae_int_t) src->cnt) * ae_sizeof((ae_datatype) src->datatype)));
}

// This function initializes ae_vector using X-structure as source.
//
// New vector is attached to source:
// * DST shares memory with SRC
// * both DST and SRC are writable - all writes to DST  change  elements  of
//   SRC and vice versa.
// * DST can be reallocated with ae_vector_set_length(), in  this  case  SRC
//   remains untouched
// * SRC, however, CAN NOT BE REALLOCATED AS LONG AS DST EXISTS
//
// NOTE: is_attached field is set  to  true  in  order  to  indicate  that
//       vector does not own its memory.
//
// dst                 destination vector
// src                 well, it is source
// make_automatic                if true, vector will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_vector_init_attach_to_x(ae_vector *dst, x_vector *src, bool make_automatic) {
   volatile ae_int_t cnt;
   AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof(*dst)));
   cnt = (ae_int_t) src->cnt;
// ensure that size is correct
   ae_assert(cnt == src->cnt, "ae_vector_init_attach_to_x(): 32/64 overflow");
   ae_assert(cnt >= 0, "ae_vector_init_attach_to_x(): negative length");
// prepare for possible errors during allocation
   dst->cnt = 0;
   dst->ptr.p_ptr = NULL;
   dst->datatype = (ae_datatype) src->datatype;
// zero-size init in order to correctly register in the frame
   ae_db_init(&dst->data, 0, make_automatic);
// init
   dst->cnt = cnt;
   dst->ptr.p_ptr = src->x_ptr.p_ptr;
   dst->is_attached = true;
}

// This function changes length of ae_vector.
//
// dst                 destination vector
// newsize             vector size, may be zero
//
// Error handling: calls ae_break() on allocation error
//
// NOTES:
// * vector must be initialized
// * all contents is free during setlength() call
// * new size may be zero.
void ae_vector_set_length(ae_vector *dst, ae_int_t newsize) {
   ae_assert(newsize >= 0, "ae_vector_set_length(): negative size");
   if (dst->cnt == newsize) return;
// realloc, being ready for exception during reallocation (cnt=ptr=0 on entry)
   dst->cnt = 0;
   dst->ptr.p_ptr = NULL;
   ae_db_realloc(&dst->data, newsize * ae_sizeof(dst->datatype));
   dst->cnt = newsize;
   dst->ptr.p_ptr = dst->data.ptr;
}

// This function resized ae_vector, preserving previously existing elements.
// Values of elements added during vector growth is undefined.
//
// dst                 destination vector
// newsize             vector size, may be zero
// state               ALGLIB environment state, can not be NULL
//
// Error handling: calls ae_break() on allocation error
//
// NOTES:
// * vector must be initialized
// * new size may be zero.
void ae_vector_resize(ae_vector *dst, ae_int_t newsize) {
   ae_int_t bytes_total;
   NewVector(tmp, newsize, dst->datatype);
   bytes_total = (dst->cnt < newsize ? dst->cnt : newsize) * ae_sizeof(dst->datatype);
   if (bytes_total > 0)
      memmove(tmp.ptr.p_ptr, dst->ptr.p_ptr, bytes_total);
   ae_swap_vectors(dst, &tmp);
   ae_vector_free(&tmp, true);
}

// This  function  provides  "FREE"  functionality  for vector (contents is
// cleared, all internal structures are free, but structure still left in valid state).
//
// The  function clears vector contents (releases all dynamically  allocated
// memory). Vector may be in automatic management list  -  in this  case  it
// will NOT be removed from list.
//
// IMPORTANT: this function does NOT invalidates dst; it just  releases  all
// dynamically allocated storage, but dst still may be used  after  call  to
// ae_vector_set_length().
//
// dst                 destination vector
void ae_vector_free(ae_vector *dst, bool make_automatic) {
   dst->cnt = 0;
   ae_db_free(&dst->data);
   dst->ptr.p_ptr = 0;
   dst->is_attached = false;
}

// This function efficiently swaps contents of two vectors, leaving other
// pararemeters (automatic management, etc.) unchanged.
void ae_swap_vectors(ae_vector *vec1, ae_vector *vec2) {
   ae_int_t cnt;
   ae_datatype datatype;
   void *p_ptr;
   ae_impose(!vec1->is_attached, "ALGLIB: internal error, attempt to swap vectors attached to X-object");
   ae_impose(!vec2->is_attached, "ALGLIB: internal error, attempt to swap vectors attached to X-object");
   ae_db_swap(&vec1->data, &vec2->data);
   cnt = vec1->cnt;
   datatype = vec1->datatype;
   p_ptr = vec1->ptr.p_ptr;
   vec1->cnt = vec2->cnt;
   vec1->datatype = vec2->datatype;
   vec1->ptr.p_ptr = vec2->ptr.p_ptr;
   vec2->cnt = cnt;
   vec2->datatype = datatype;
   vec2->ptr.p_ptr = p_ptr;
}

// This function creates ae_matrix.
//
// Matrix size may be zero, in such cases both rows and cols are zero.
// Matrix contents is uninitialized.
//
// dst                 destination matrix, must be zero-filled
// rows                rows count
// cols                cols count
// datatype            element type
// make_automatic                if true, matrix will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
//
// NOTE: no memory allocation is performed for initialization with rows=cols=0
void ae_matrix_init(ae_matrix *dst, ae_int_t rows, ae_int_t cols, ae_datatype datatype, bool make_automatic) {
   AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof(*dst)));
   ae_assert(rows >= 0 && cols >= 0, "ae_matrix_init(): negative length");
// if one of rows/cols is zero, another MUST be too; perform quick exit
   if (rows == 0 || cols == 0) {
      dst->rows = 0;
      dst->cols = 0;
      dst->is_attached = false;
      dst->ptr.pp_void = NULL;
      dst->stride = 0;
      dst->datatype = datatype;
      ae_db_init(&dst->data, 0, make_automatic);
      return;
   }
// init, being ready for exception during allocation (rows=cols=ptr=NULL on entry)
   dst->is_attached = false;
   dst->rows = 0;
   dst->cols = 0;
   dst->ptr.pp_void = NULL;
   dst->stride = cols;
   while (dst->stride * ae_sizeof(datatype) % AE_DATA_ALIGN != 0) dst->stride++;
   dst->datatype = datatype;
   ae_db_init(&dst->data, rows * ((ae_int_t) sizeof(void *) + dst->stride * ae_sizeof(datatype)) + AE_DATA_ALIGN - 1, make_automatic);
   dst->rows = rows;
   dst->cols = cols;
   ae_matrix_update_row_pointers(dst, ae_align((char *)dst->data.ptr + rows * sizeof(void *), AE_DATA_ALIGN));
}

// This function creates copy of ae_matrix. A new copy of the data is created.
//
// dst                 destination matrix, must be zero-filled
// src                 well, it is source
// make_automatic                if true, matrix will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_matrix_copy(ae_matrix *dst, ae_matrix *src, bool make_automatic) {
   ae_int_t i;
   ae_matrix_init(dst, src->rows, src->cols, src->datatype, make_automatic);
   if (src->rows != 0 && src->cols != 0) {
      if (dst->stride == src->stride)
         memmove(dst->ptr.pp_void[0], src->ptr.pp_void[0], (size_t)(src->rows * src->stride * ae_sizeof(src->datatype)));
      else for (i = 0; i < dst->rows; i++)
         memmove(dst->ptr.pp_void[i], src->ptr.pp_void[i], (size_t)(dst->cols * ae_sizeof(dst->datatype)));
   }
}

// This function initializes ae_matrix using X-structure as source. New copy
// of data is created, which is owned/managed by ae_matrix  structure.  Both
// structures (source and destination) remain completely  independent  after
// this call.
//
// dst                 destination matrix, must be zero-filled
// src                 well, it is source
// make_automatic                if true, matrix will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_matrix_init_from_x(ae_matrix *dst, x_matrix *src, bool make_automatic) {
   char *p_src_row;
   char *p_dst_row;
   ae_int_t row_size;
   ae_int_t i;
   ae_matrix_init(dst, (ae_int_t) src->rows, (ae_int_t) src->cols, (ae_datatype) src->datatype, make_automatic);
   if (src->rows != 0 && src->cols != 0) {
      p_src_row = (char *)src->x_ptr.p_ptr;
      p_dst_row = (char *)(dst->ptr.pp_void[0]);
      row_size = ae_sizeof((ae_datatype) src->datatype) * (ae_int_t) src->cols;
      for (i = 0; i < src->rows; i++, p_src_row += src->stride * ae_sizeof((ae_datatype) src->datatype), p_dst_row += dst->stride * ae_sizeof((ae_datatype) src->datatype))
         memmove(p_dst_row, p_src_row, (size_t)(row_size));
   }
}

// This function initializes ae_matrix using X-structure as source.
//
// New matrix is attached to source:
// * DST shares memory with SRC
// * both DST and SRC are writable - all writes to DST  change  elements  of
//   SRC and vice versa.
// * DST can be reallocated with ae_matrix_set_length(), in  this  case  SRC
//   remains untouched
// * SRC, however, CAN NOT BE REALLOCATED AS LONG AS DST EXISTS
//
// dst                 destination matrix, must be zero-filled
// src                 well, it is source
// make_automatic                if true, matrix will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_matrix_init_attach_to_x(ae_matrix *dst, x_matrix *src, bool make_automatic) {
   ae_int_t rows, cols;
   AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof(*dst)));
   rows = (ae_int_t) src->rows;
   cols = (ae_int_t) src->cols;
// check that X-source is densely packed
   ae_assert(src->cols == src->stride, "ae_matrix_init_attach_to_x(): unsupported stride");
// ensure that size is correct
   ae_assert(rows == src->rows, "ae_matrix_init_attach_to_x(): 32/64 overflow");
   ae_assert(cols == src->cols, "ae_matrix_init_attach_to_x(): 32/64 overflow");
   ae_assert(rows >= 0 && cols >= 0, "ae_matrix_init_attach_to_x(): negative length");
// if one of rows/cols is zero, another MUST be too
   if (rows == 0 || cols == 0) {
      rows = 0;
      cols = 0;
   }
// init, being ready for allocation error
   dst->is_attached = true;
   dst->rows = 0;
   dst->cols = 0;
   dst->stride = cols;
   dst->datatype = (ae_datatype) src->datatype;
   dst->ptr.pp_void = NULL;
   ae_db_init(&dst->data, rows * (ae_int_t) sizeof(void *), make_automatic);
   dst->rows = rows;
   dst->cols = cols;
   if (dst->rows > 0 && dst->cols > 0) {
      ae_int_t i, rowsize;
      char *p_row;
      void **pp_ptr;
      p_row = (char *)src->x_ptr.p_ptr;
      rowsize = dst->stride * ae_sizeof(dst->datatype);
      pp_ptr = (void **)dst->data.ptr;
      dst->ptr.pp_void = pp_ptr;
      for (i = 0; i < dst->rows; i++, p_row += rowsize) pp_ptr[i] = p_row;
   }
}

// This function changes length of ae_matrix.
//
// dst                 destination matrix
// rows                size, may be zero
// cols                size, may be zero
//
// Error handling:
// * if TopFr == NULL, returns false on allocation error
// * if TopFr != NULL, calls ae_break() on allocation error
// * returns true on success
//
// NOTES:
// * matrix must be initialized
// * all contents is freed during setlength() call
// * new size may be zero.
void ae_matrix_set_length(ae_matrix *dst, ae_int_t rows, ae_int_t cols) {
   ae_assert(rows >= 0 && cols >= 0, "ae_matrix_set_length(): negative length");
   if (dst->rows == rows && dst->cols == cols) return;
// prepare stride
   dst->stride = cols;
   while (dst->stride * ae_sizeof(dst->datatype) % AE_DATA_ALIGN != 0) dst->stride++;
// realloc, being ready for an exception during reallocation (rows=cols=0 on entry)
   dst->rows = 0;
   dst->cols = 0;
   dst->ptr.pp_void = NULL;
   ae_db_realloc(&dst->data, rows * ((ae_int_t) sizeof(void *) + dst->stride * ae_sizeof(dst->datatype)) + AE_DATA_ALIGN - 1);
   dst->rows = rows;
   dst->cols = cols;
// update pointers to rows
   ae_matrix_update_row_pointers(dst, ae_align((char *)dst->data.ptr + dst->rows * sizeof(void *), AE_DATA_ALIGN));
}

// This  function  provides  "FREE"  functionality  for matrix (contents is
// cleared, but structure still left in valid state).
//
// The  function clears matrix contents (releases all dynamically  allocated
// memory). Matrix may be in automatic management list  -  in this  case  it
// will NOT be removed from list.
//
// IMPORTANT: this function does NOT invalidates dst; it just  releases  all
// dynamically allocated storage, but dst still may be used  after  call  to
// ae_matrix_set_length().
//
// dst                 destination matrix
void ae_matrix_free(ae_matrix *dst, bool make_automatic) {
   dst->rows = 0;
   dst->cols = 0;
   dst->stride = 0;
   ae_db_free(&dst->data);
   dst->ptr.p_ptr = 0;
   dst->is_attached = false;
}

// This function efficiently swaps contents of two vectors, leaving other
// pararemeters (automatic management, etc.) unchanged.
void ae_swap_matrices(ae_matrix *mat1, ae_matrix *mat2) {
   ae_int_t rows;
   ae_int_t cols;
   ae_int_t stride;
   ae_datatype datatype;
   void *p_ptr;
   ae_impose(!mat1->is_attached, "ALGLIB: internal error, attempt to swap matrices attached to X-object");
   ae_impose(!mat2->is_attached, "ALGLIB: internal error, attempt to swap matrices attached to X-object");
   ae_db_swap(&mat1->data, &mat2->data);
   rows = mat1->rows;
   cols = mat1->cols;
   stride = mat1->stride;
   datatype = mat1->datatype;
   p_ptr = mat1->ptr.p_ptr;
   mat1->rows = mat2->rows;
   mat1->cols = mat2->cols;
   mat1->stride = mat2->stride;
   mat1->datatype = mat2->datatype;
   mat1->ptr.p_ptr = mat2->ptr.p_ptr;
   mat2->rows = rows;
   mat2->cols = cols;
   mat2->stride = stride;
   mat2->datatype = datatype;
   mat2->ptr.p_ptr = p_ptr;
}

// This function creates smart pointer structure.
//
// dst                 destination smart pointer, must be zero-filled
// subscriber          pointer to pointer which receives updates in the
//                     internal object stored in ae_smart_ptr. Any update to
//                     dst->ptr is translated to subscriber. Can be NULL.
// make_automatic                if true, pointer will be registered in the frame of the current state;
//
// Error handling:
// * on failure calls ae_break() with TopFr == NULL. Usually it  results
//   in abort() call.
//
// After initialization, smart pointer stores NULL pointer.
void ae_smart_ptr_init(ae_smart_ptr *dst, void **subscriber, bool make_automatic) {
   AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof(*dst)));
   dst->subscriber = subscriber;
   dst->ptr = NULL;
   if (dst->subscriber != NULL) *(dst->subscriber) = dst->ptr;
   dst->is_owner = false;
   dst->is_dynamic = false;
   dst->frame_entry.deallocator = ae_smart_ptr_free;
   dst->frame_entry.ptr = dst;
   if (make_automatic) ae_db_attach(&dst->frame_entry);
}

// This function clears smart pointer structure.
//
// dst                 destination smart pointer.
//
// After call to this function smart pointer contains NULL reference,  which
// is  propagated  to  its  subscriber  (in  cases  non-NULL  subscruber was
// specified during pointer creation).
void ae_smart_ptr_free(void *_dst) {
   ae_smart_ptr *dst = (ae_smart_ptr *) _dst;
   if (dst->is_owner && dst->ptr != NULL) {
      dst->free(dst->ptr, false);
      if (dst->is_dynamic) ae_free(dst->ptr);
   }
   dst->is_owner = false;
   dst->is_dynamic = false;
   dst->ptr = NULL;
   dst->free = NULL;
   if (dst->subscriber != NULL) *(dst->subscriber) = NULL;
}

// This function assigns pointer to ae_smart_ptr structure.
//
// dst                 destination smart pointer.
// new_ptr             new pointer to assign
// is_owner            whether smart pointer owns new_ptr
// is_dynamic          whether object is dynamic - clearing such object
//                     requires BOTH calling destructor function AND calling
//                     ae_free() for memory occupied by object.
// fre                 destructor function
//
// In case smart pointer already contains non-NULL value and owns this value,
// it is freed before assigning new pointer.
//
// Changes in pointer are propagated to its  subscriber  (in  case  non-NULL
// subscriber was specified during pointer creation).
//
// You can specify NULL new_ptr, in which case is_owner/free are ignored.
void ae_smart_ptr_assign(ae_smart_ptr *dst, void *new_ptr, bool is_owner, bool is_dynamic, void (*free)(void *, bool make_automatic)) {
   if (dst->is_owner && dst->ptr != NULL) {
      dst->free(dst->ptr, false);
      if (dst->is_dynamic) ae_free(dst->ptr);
   }
   if (new_ptr != NULL) {
      dst->ptr = new_ptr;
      dst->is_owner = is_owner;
      dst->is_dynamic = is_dynamic;
      dst->free = free;
   } else {
      dst->ptr = NULL;
      dst->is_owner = false;
      dst->is_dynamic = false;
      dst->free = NULL;
   }
   if (dst->subscriber != NULL) *(dst->subscriber) = dst->ptr;
}

// This function releases pointer owned by ae_smart_ptr structure:
// * all internal fields are set to NULL
// * destructor function for internal pointer is NOT called even when we own
//   this pointer. After this call ae_smart_ptr releases  ownership  of  its
//   pointer and passes it to caller.
// * changes in pointer are propagated to its subscriber (in  case  non-NULL
//   subscriber was specified during pointer creation).
//
// dst                 destination smart pointer.
void ae_smart_ptr_release(ae_smart_ptr *dst) {
   dst->is_owner = false;
   dst->is_dynamic = false;
   dst->ptr = NULL;
   dst->free = NULL;
   if (dst->subscriber != NULL)
      *(dst->subscriber) = NULL;
}

// This function copies contents of ae_vector (SRC) to x_vector (DST).
//
// This function should not be called for  DST  which  is  attached  to  SRC
// (opposite situation, when SRC is attached to DST, is possible).
//
// Depending on situation, following actions are performed
// * for SRC attached to DST, this function performs no actions (no need  to
//   do anything)
// * for independent vectors of different sizes it allocates storage in  DST
//   and copy contents of SRC  to  DST.  DST->last_action field  is  set  to
//   ACT_NEW_LOCATION, and DST->owner is set to OWN_AE.
// * for  independent  vectors   of  same  sizes  it does not perform memory
//   (re)allocation.  It  just  copies  SRC  to  already   existing   place.
//   DST->last_action   is   set   to    ACT_SAME_LOCATION  (unless  it  was
//   ACT_NEW_LOCATION), DST->owner is unmodified.
//
// dst                 destination vector
// src                 source, vector in x-format
//
// NOTES:
// * dst is assumed to be initialized. Its contents is freed before  copying
//   data  from src  (if  size / type  are  different)  or  overwritten  (if
//   possible given destination size).
void ae_x_set_vector(x_vector *dst, ae_vector *src) {
   if (src->ptr.p_ptr == dst->x_ptr.p_ptr) {
   // src->ptr points to the beginning of dst, attached matrices, no need to copy
      return;
   }
   if (dst->cnt != src->cnt || dst->datatype != src->datatype) {
      if (dst->owner == OWN_AE) ae_free(dst->x_ptr.p_ptr);
      dst->x_ptr.p_ptr = ae_malloc((size_t)(src->cnt * ae_sizeof(src->datatype)));
      if (src->cnt != 0 && dst->x_ptr.p_ptr == NULL) ae_break(ERR_OUT_OF_MEMORY, "ae_malloc(): out of memory");
      dst->last_action = ACT_NEW_LOCATION;
      dst->cnt = src->cnt;
      dst->datatype = src->datatype;
      dst->owner = OWN_AE;
   } else {
      if (dst->last_action == ACT_UNCHANGED) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_SAME_LOCATION) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_NEW_LOCATION) dst->last_action = ACT_NEW_LOCATION;
      else ae_assert(false, "ALGLIB: internal error in ae_x_set_vector()");
   }
   if (src->cnt) memmove(dst->x_ptr.p_ptr, src->ptr.p_ptr, (size_t)(src->cnt * ae_sizeof(src->datatype)));
}

// This function copies contents of ae_matrix to x_matrix.
//
// This function should not be called for  DST  which  is  attached  to  SRC
// (opposite situation, when SRC is attached to DST, is possible).
//
// Depending on situation, following actions are performed
// * for SRC attached to DST, this function performs no actions (no need  to
//   do anything)
// * for independent matrices of different sizes it allocates storage in DST
//   and copy contents of SRC  to  DST.  DST->last_action field  is  set  to
//   ACT_NEW_LOCATION, and DST->owner is set to OWN_AE.
// * for  independent  matrices  of  same  sizes  it does not perform memory
//   (re)allocation.  It  just  copies  SRC  to  already   existing   place.
//   DST->last_action   is   set   to    ACT_SAME_LOCATION  (unless  it  was
//   ACT_NEW_LOCATION), DST->owner is unmodified.
//
// dst                 destination vector
// src                 source, matrix in x-format
//
// NOTES:
// * dst is assumed to be initialized. Its contents is freed before  copying
//   data  from src  (if  size / type  are  different)  or  overwritten  (if
//   possible given destination size).
void ae_x_set_matrix(x_matrix *dst, ae_matrix *src) {
   char *p_src_row;
   char *p_dst_row;
   ae_int_t i;
   ae_int_t row_size;
   if (src->ptr.pp_void != NULL && src->ptr.pp_void[0] == dst->x_ptr.p_ptr) {
   // src->ptr points to the beginning of dst, attached matrices, no need to copy
      return;
   }
   if (dst->rows != src->rows || dst->cols != src->cols || dst->datatype != src->datatype) {
      if (dst->owner == OWN_AE) ae_free(dst->x_ptr.p_ptr);
      dst->rows = src->rows;
      dst->cols = src->cols;
      dst->stride = src->cols;
      dst->datatype = src->datatype;
      dst->x_ptr.p_ptr = ae_malloc((size_t)(dst->rows * ((ae_int_t) dst->stride) * ae_sizeof(src->datatype)));
      if (dst->rows != 0 && dst->stride != 0 && dst->x_ptr.p_ptr == NULL) ae_break(ERR_OUT_OF_MEMORY, "ae_malloc(): out of memory");
      dst->last_action = ACT_NEW_LOCATION;
      dst->owner = OWN_AE;
   } else {
      if (dst->last_action == ACT_UNCHANGED) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_SAME_LOCATION) dst->last_action = ACT_SAME_LOCATION;
      else if (dst->last_action == ACT_NEW_LOCATION) dst->last_action = ACT_NEW_LOCATION;
      else ae_assert(false, "ALGLIB: internal error in ae_x_set_vector()");
   }
   if (src->rows != 0 && src->cols != 0) {
      p_src_row = (char *)(src->ptr.pp_void[0]);
      p_dst_row = (char *)dst->x_ptr.p_ptr;
      row_size = ae_sizeof(src->datatype) * src->cols;
      for (i = 0; i < src->rows; i++, p_src_row += src->stride * ae_sizeof(src->datatype), p_dst_row += dst->stride * ae_sizeof(src->datatype))
         memmove(p_dst_row, p_src_row, (size_t)(row_size));
   }
}

// This function attaches x_vector to ae_vector's contents.
// Ownership of memory allocated is not changed (it is still managed by
// ae_matrix).
//
// dst                 destination vector
// src                 source, vector in x-format
//
// NOTES:
// * dst is assumed to be initialized. Its contents is freed before attaching to src.
// * this function doesn't need TopFr because it can't fail (assuming correctly initialized src)
void ae_x_attach_to_vector(x_vector *dst, ae_vector *src) {
   if (dst->owner == OWN_AE) ae_free(dst->x_ptr.p_ptr);
   dst->x_ptr.p_ptr = src->ptr.p_ptr;
   dst->last_action = ACT_NEW_LOCATION;
   dst->cnt = src->cnt;
   dst->datatype = src->datatype;
   dst->owner = OWN_CALLER;
}

// This function attaches x_matrix to ae_matrix's contents.
// Ownership of memory allocated is not changed (it is still managed by
// ae_matrix).
//
// dst                 destination vector
// src                 source, matrix in x-format
//
// NOTES:
// * dst is assumed to be initialized. Its contents is freed before attaching to src.
// * this function doesn't need TopFr because it can't fail (assuming correctly initialized src)
void ae_x_attach_to_matrix(x_matrix *dst, ae_matrix *src) {
   if (dst->owner == OWN_AE) ae_free(dst->x_ptr.p_ptr);
   dst->rows = src->rows;
   dst->cols = src->cols;
   dst->stride = src->stride;
   dst->datatype = src->datatype;
   dst->x_ptr.p_ptr = src->ptr.pp_double[0];
   dst->last_action = ACT_NEW_LOCATION;
   dst->owner = OWN_CALLER;
}

// This function clears x_vector. It does nothing  if vector is not owned by
// ALGLIB environment.
//
// dst                 vector
void x_vector_free(x_vector *dst, bool make_automatic) {
   if (dst->owner == OWN_AE)
      aligned_free(dst->x_ptr.p_ptr);
   dst->x_ptr.p_ptr = NULL;
   dst->cnt = 0;
}

// CPUID
//
// Returns information about features CPU and compiler support.
//
// You must tell ALGLIB what CPU family is used by defining AE_CPU symbol
// (without this hint zero will be returned).
//
// Note: results of this function depend on both CPU and compiler;
// if compiler doesn't support SSE intrinsics, function won't set
// corresponding flag.
// static volatile bool _ae_cpuid_initialized = false; //(@) Removed: reduced to a start-up routine to initialize CurCPU, in place of ae_cpuid() calls.
// static volatile bool _ae_cpuid_has_sse2 = false; //(@) Removed: integrated into CurCPU.
static ae_int_t ae_cpuid() {
   ae_int_t result = 0;
// Determine system properties
#if defined AE_CPU && AE_CPU == AE_INTEL
   { // SSE2
#   if defined AE_HAS_SSE2_INTRINSICS
#      if AE_COMPILER == AE_GNUC || AE_COMPILER == AE_SUNC
      ae_int_t a, b, c, d;
      __asm__ __volatile__("cpuid":"=a"(a), "=b"(b), "=c"(c), "=d"(d):"a"(1));
      if (d&0x04000000) result |= CPU_SSE2;
#      elif AE_COMPILER == AE_MSVC
      int CPUInfo[4];
      __cpuid(CPUInfo, 1);
      if (CPUInfo[3]&0x04000000) result |= CPU_SSE2;
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
#   if AE_COMPILER == AE_GNUC
      unsigned CpuMode = 0x27f; // Set rounding for floating-point math to double precision for x86/x64 processors under GCC.
      asm("fldcw %0" : : "m"(*&CpuMode));
#   endif
   }
#endif
   return result;
}

// Determine CPU characteristics and perform CPU-specific initializations.
// This can be safely done on a per-thread basis provided that:
// simultaneous writes by different cores to the same location will be executed in serial manner,
// which is true of contemporary architectures.
const/* AutoS */ae_int_t CurCPU = ae_cpuid();

// Real math functions
static ae_int_t GetByteOrder() {
// determine the byte order
// two types are supported: big-endian and little-endian.
// mixed-endian hardware is NOT supported.
//
// 1983 is used as magic number because its non-periodic double
// representation allow us to easily distinguish between upper
// and lower halfs and to detect mixed endian hardware.
   union { double a; ae_int32_t p[2]; } u; u.a = 1.0 / 1983.0;
   if (u.p[1] == (ae_int32_t) 0x3f408642) return AE_LITTLE_ENDIAN;
   if (u.p[0] == (ae_int32_t) 0x3f408642) return AE_BIG_ENDIAN;
   return AE_MIXED_ENDIAN; //(@) Originally this prompted an abort().
}

// determine byte order and initialize precomputed IEEE special quantities.
// byte order types: AE_LITTLE_ENDIAN or AE_BIG_ENDIAN
static const ae_int_t ByteOrder = GetByteOrder();

bool isposinf(double x) { return isinf(x) && !signbit(x); }
bool isneginf(double x) { return isinf(x) && signbit(x); }
ae_int_t ae_iabs(ae_int_t x) { return x >= 0? x: -x; }
double ae_sqr(double x) { return x*x; }
ae_int_t ae_sign(double x) { return x > 0? +1: x < 0? -1: 0; }
ae_int_t RoundZ(double x) { return (ae_int_t)round(x); }
ae_int_t TruncZ(double x) { return (ae_int_t)trunc(x); }
ae_int_t FloorZ(double x) { return (ae_int_t)floor(x); }
ae_int_t CeilZ(double x) { return (ae_int_t)ceil(x); }
ae_int_t ae_maxint(ae_int_t m1, ae_int_t m2) { return m1 > m2 ? m1 : m2; }
ae_int_t ae_minint(ae_int_t m1, ae_int_t m2) { return m1 > m2 ? m2 : m1; }
double ae_maxreal(double m1, double m2) { return m1 > m2 ? m1 : m2; }
double ae_minreal(double m1, double m2) { return m1 > m2 ? m2 : m1; }

double ae_randomreal() {
   const double mx = (double)RAND_MAX + 1.0;
   return (rand() + rand() / mx) / mx;
}
double ae_randommid() {
   const double mx = (double)RAND_MAX + 1.0;
   return 2.0*(rand() + rand() / mx) / mx - 1.0;
}
ae_int_t ae_randominteger(ae_int_t maxv) { return rand() % maxv; }
bool ae_randombool(double p) {
   const double mx = (double)RAND_MAX + 1.0;
   return rand() + rand()/mx <= p*mx;
}

// Symmetric/Hermitian properties: check and force
static void x_split_length(ae_int_t n, ae_int_t nb, ae_int_t *n1, ae_int_t *n2) {
   ae_int_t r;
   if (n <= nb) {
      *n1 = n;
      *n2 = 0;
   } else {
      if (n % nb != 0) {
         *n2 = n % nb;
         *n1 = n - (*n2);
      } else {
         *n2 = n / 2;
         *n1 = n - (*n2);
         if (*n1 % nb == 0) {
            return;
         }
         r = nb - *n1 % nb;
         *n1 += r;
         *n2 -= r;
      }
   }
}

static double x_safepythag2(double x, double y) {
   double w;
   double xabs;
   double yabs;
   double z;
   xabs = fabs(x);
   yabs = fabs(y);
   w = xabs > yabs ? xabs : yabs;
   z = xabs < yabs ? xabs : yabs;
   if (z == 0)
      return w;
   else {
      double t;
      t = z / w;
      return w * sqrt(1 + t * t);
   }
}

// this function checks difference between offdiagonal blocks BL and BU
// (see below). Block BL is specified by offsets (offset0,offset1)  and
// sizes (len0,len1).
//
//     [ .          ]
//     [   A0  BU   ]
// A = [   BL  A1   ]
//     [          . ]
//
//  this subroutine updates current values of:
//  a) mx       maximum value of A[i,j] found so far
//  b) err      componentwise difference between elements of BL and BU^T
//
static void is_symmetric_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1, bool *nonfinite, double *mx, double *err) {
// try to split problem into two smaller ones
   if (len0 > x_nb || len1 > x_nb) {
      ae_int_t n1, n2;
      if (len0 > len1) {
         x_split_length(len0, x_nb, &n1, &n2);
         is_symmetric_rec_off_stat(a, offset0, offset1, n1, len1, nonfinite, mx, err);
         is_symmetric_rec_off_stat(a, offset0 + n1, offset1, n2, len1, nonfinite, mx, err);
      } else {
         x_split_length(len1, x_nb, &n1, &n2);
         is_symmetric_rec_off_stat(a, offset0, offset1, len0, n1, nonfinite, mx, err);
         is_symmetric_rec_off_stat(a, offset0, offset1 + n1, len0, n2, nonfinite, mx, err);
      }
      return;
   } else {
   // base case
      double *p1, *p2, *prow, *pcol;
      double v;
      ae_int_t i, j;
      p1 = (double *)(a->x_ptr.p_ptr) + offset0 * a->stride + offset1;
      p2 = (double *)(a->x_ptr.p_ptr) + offset1 * a->stride + offset0;
      for (i = 0; i < len0; i++) {
         pcol = p2 + i;
         prow = p1 + i * a->stride;
         for (j = 0; j < len1; j++) {
            if (!isfinite(*pcol) || !isfinite(*prow)) {
               *nonfinite = true;
            } else {
               v = fabs(*pcol);
               *mx = *mx > v ? *mx : v;
               v = fabs(*prow);
               *mx = *mx > v ? *mx : v;
               v = fabs(*pcol - *prow);
               *err = *err > v ? *err : v;
            }
            pcol += a->stride;
            prow++;
         }
      }
   }
}

// this function checks that diagonal block A0 is symmetric.
// Block A0 is specified by its offset and size.
//
//     [ .          ]
//     [   A0       ]
// A = [       .    ]
//     [          . ]
//
//  this subroutine updates current values of:
//  a) mx       maximum value of A[i,j] found so far
//  b) err      componentwise difference between A0 and A0^T
//
static void is_symmetric_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len, bool *nonfinite, double *mx, double *err) {
   double *p, *prow, *pcol;
   double v;
   ae_int_t i, j;
// try to split problem into two smaller ones
   if (len > x_nb) {
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      is_symmetric_rec_diag_stat(a, offset, n1, nonfinite, mx, err);
      is_symmetric_rec_diag_stat(a, offset + n1, n2, nonfinite, mx, err);
      is_symmetric_rec_off_stat(a, offset + n1, offset, n2, n1, nonfinite, mx, err);
      return;
   }
// base case
   p = (double *)(a->x_ptr.p_ptr) + offset * a->stride + offset;
   for (i = 0; i < len; i++) {
      pcol = p + i;
      prow = p + i * a->stride;
      for (j = 0; j < i; j++, pcol += a->stride, prow++) {
         if (!isfinite(*pcol) || !isfinite(*prow)) {
            *nonfinite = true;
         } else {
            v = fabs(*pcol);
            *mx = *mx > v ? *mx : v;
            v = fabs(*prow);
            *mx = *mx > v ? *mx : v;
            v = fabs(*pcol - *prow);
            *err = *err > v ? *err : v;
         }
      }
      v = fabs(p[i + i * a->stride]);
      *mx = *mx > v ? *mx : v;
   }
}

// this function checks difference between offdiagonal blocks BL and BU
// (see below). Block BL is specified by offsets (offset0,offset1)  and
// sizes (len0,len1).
//
//     [ .          ]
//     [   A0  BU   ]
// A = [   BL  A1   ]
//     [          . ]
//
//  this subroutine updates current values of:
//  a) mx       maximum value of A[i,j] found so far
//  b) err      componentwise difference between elements of BL and BU^H
//
static void is_hermitian_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1, bool *nonfinite, double *mx, double *err) {
// try to split problem into two smaller ones
   if (len0 > x_nb || len1 > x_nb) {
      ae_int_t n1, n2;
      if (len0 > len1) {
         x_split_length(len0, x_nb, &n1, &n2);
         is_hermitian_rec_off_stat(a, offset0, offset1, n1, len1, nonfinite, mx, err);
         is_hermitian_rec_off_stat(a, offset0 + n1, offset1, n2, len1, nonfinite, mx, err);
      } else {
         x_split_length(len1, x_nb, &n1, &n2);
         is_hermitian_rec_off_stat(a, offset0, offset1, len0, n1, nonfinite, mx, err);
         is_hermitian_rec_off_stat(a, offset0, offset1 + n1, len0, n2, nonfinite, mx, err);
      }
      return;
   } else {
   // base case
      ae_complex *p1, *p2, *prow, *pcol;
      double v;
      ae_int_t i, j;
      p1 = (ae_complex *)(a->x_ptr.p_ptr) + offset0 * a->stride + offset1;
      p2 = (ae_complex *)(a->x_ptr.p_ptr) + offset1 * a->stride + offset0;
      for (i = 0; i < len0; i++) {
         pcol = p2 + i;
         prow = p1 + i * a->stride;
         for (j = 0; j < len1; j++) {
            if (!isfinite(pcol->x) || !isfinite(pcol->y) || !isfinite(prow->x) || !isfinite(prow->y)) {
               *nonfinite = true;
            } else {
               v = x_safepythag2(pcol->x, pcol->y);
               *mx = *mx > v ? *mx : v;
               v = x_safepythag2(prow->x, prow->y);
               *mx = *mx > v ? *mx : v;
               v = x_safepythag2(pcol->x - prow->x, pcol->y + prow->y);
               *err = *err > v ? *err : v;
            }
            pcol += a->stride;
            prow++;
         }
      }
   }
}

// this function checks that diagonal block A0 is Hermitian.
// Block A0 is specified by its offset and size.
//
//     [ .          ]
//     [   A0       ]
// A = [       .    ]
//     [          . ]
//
//  this subroutine updates current values of:
//  a) mx       maximum value of A[i,j] found so far
//  b) err      componentwise difference between A0 and A0^H
//
static void is_hermitian_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len, bool *nonfinite, double *mx, double *err) {
   ae_complex *p, *prow, *pcol;
   double v;
   ae_int_t i, j;
// try to split problem into two smaller ones
   if (len > x_nb) {
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      is_hermitian_rec_diag_stat(a, offset, n1, nonfinite, mx, err);
      is_hermitian_rec_diag_stat(a, offset + n1, n2, nonfinite, mx, err);
      is_hermitian_rec_off_stat(a, offset + n1, offset, n2, n1, nonfinite, mx, err);
      return;
   }
// base case
   p = (ae_complex *)(a->x_ptr.p_ptr) + offset * a->stride + offset;
   for (i = 0; i < len; i++) {
      pcol = p + i;
      prow = p + i * a->stride;
      for (j = 0; j < i; j++, pcol += a->stride, prow++) {
         if (!isfinite(pcol->x) || !isfinite(pcol->y) || !isfinite(prow->x) || !isfinite(prow->y)) {
            *nonfinite = true;
         } else {
            v = x_safepythag2(pcol->x, pcol->y);
            *mx = *mx > v ? *mx : v;
            v = x_safepythag2(prow->x, prow->y);
            *mx = *mx > v ? *mx : v;
            v = x_safepythag2(pcol->x - prow->x, pcol->y + prow->y);
            *err = *err > v ? *err : v;
         }
      }
      if (!isfinite(p[i + i * a->stride].x) || !isfinite(p[i + i * a->stride].y)) {
         *nonfinite = true;
      } else {
         v = fabs(p[i + i * a->stride].x);
         *mx = *mx > v ? *mx : v;
         v = fabs(p[i + i * a->stride].y);
         *err = *err > v ? *err : v;
      }
   }
}

// this function copies offdiagonal block BL to its symmetric counterpart
// BU (see below). Block BL is specified by offsets (offset0,offset1)
// and sizes (len0,len1).
//
//     [ .          ]
//     [   A0  BU   ]
// A = [   BL  A1   ]
//     [          . ]
//
static void force_symmetric_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1) {
// try to split problem into two smaller ones
   if (len0 > x_nb || len1 > x_nb) {
      ae_int_t n1, n2;
      if (len0 > len1) {
         x_split_length(len0, x_nb, &n1, &n2);
         force_symmetric_rec_off_stat(a, offset0, offset1, n1, len1);
         force_symmetric_rec_off_stat(a, offset0 + n1, offset1, n2, len1);
      } else {
         x_split_length(len1, x_nb, &n1, &n2);
         force_symmetric_rec_off_stat(a, offset0, offset1, len0, n1);
         force_symmetric_rec_off_stat(a, offset0, offset1 + n1, len0, n2);
      }
      return;
   } else {
   // base case
      double *p1, *p2, *prow, *pcol;
      ae_int_t i, j;
      p1 = (double *)(a->x_ptr.p_ptr) + offset0 * a->stride + offset1;
      p2 = (double *)(a->x_ptr.p_ptr) + offset1 * a->stride + offset0;
      for (i = 0; i < len0; i++) {
         pcol = p2 + i;
         prow = p1 + i * a->stride;
         for (j = 0; j < len1; j++) {
            *pcol = *prow;
            pcol += a->stride;
            prow++;
         }
      }
   }
}

// this function copies lower part of diagonal block A0 to its upper part
// Block is specified by offset and size.
//
//     [ .          ]
//     [   A0       ]
// A = [       .    ]
//     [          . ]
//
static void force_symmetric_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len) {
   double *p, *prow, *pcol;
   ae_int_t i, j;
// try to split problem into two smaller ones
   if (len > x_nb) {
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      force_symmetric_rec_diag_stat(a, offset, n1);
      force_symmetric_rec_diag_stat(a, offset + n1, n2);
      force_symmetric_rec_off_stat(a, offset + n1, offset, n2, n1);
      return;
   }
// base case
   p = (double *)(a->x_ptr.p_ptr) + offset * a->stride + offset;
   for (i = 0; i < len; i++) {
      pcol = p + i;
      prow = p + i * a->stride;
      for (j = 0; j < i; j++, pcol += a->stride, prow++)
         *pcol = *prow;
   }
}

// this function copies Hermitian transpose of offdiagonal block BL to
// its symmetric counterpart BU (see below). Block BL is specified by
// offsets (offset0,offset1) and sizes (len0,len1).
//
//     [ .          ]
//     [   A0  BU   ]
// A = [   BL  A1   ]
//     [          . ]
static void force_hermitian_rec_off_stat(x_matrix *a, ae_int_t offset0, ae_int_t offset1, ae_int_t len0, ae_int_t len1) {
// try to split problem into two smaller ones
   if (len0 > x_nb || len1 > x_nb) {
      ae_int_t n1, n2;
      if (len0 > len1) {
         x_split_length(len0, x_nb, &n1, &n2);
         force_hermitian_rec_off_stat(a, offset0, offset1, n1, len1);
         force_hermitian_rec_off_stat(a, offset0 + n1, offset1, n2, len1);
      } else {
         x_split_length(len1, x_nb, &n1, &n2);
         force_hermitian_rec_off_stat(a, offset0, offset1, len0, n1);
         force_hermitian_rec_off_stat(a, offset0, offset1 + n1, len0, n2);
      }
      return;
   } else {
   // base case
      ae_complex *p1, *p2, *prow, *pcol;
      ae_int_t i, j;
      p1 = (ae_complex *)(a->x_ptr.p_ptr) + offset0 * a->stride + offset1;
      p2 = (ae_complex *)(a->x_ptr.p_ptr) + offset1 * a->stride + offset0;
      for (i = 0; i < len0; i++) {
         pcol = p2 + i;
         prow = p1 + i * a->stride;
         for (j = 0; j < len1; j++) {
            *pcol = *prow;
            pcol += a->stride;
            prow++;
         }
      }
   }
}

// this function copies Hermitian transpose of lower part of
// diagonal block A0 to its upper part Block is specified by offset and size.
//
//     [ .          ]
//     [   A0       ]
// A = [       .    ]
//     [          . ]
//
static void force_hermitian_rec_diag_stat(x_matrix *a, ae_int_t offset, ae_int_t len) {
   ae_complex *p, *prow, *pcol;
   ae_int_t i, j;
// try to split problem into two smaller ones
   if (len > x_nb) {
      ae_int_t n1, n2;
      x_split_length(len, x_nb, &n1, &n2);
      force_hermitian_rec_diag_stat(a, offset, n1);
      force_hermitian_rec_diag_stat(a, offset + n1, n2);
      force_hermitian_rec_off_stat(a, offset + n1, offset, n2, n1);
      return;
   }
// base case
   p = (ae_complex *)(a->x_ptr.p_ptr) + offset * a->stride + offset;
   for (i = 0; i < len; i++) {
      pcol = p + i;
      prow = p + i * a->stride;
      for (j = 0; j < i; j++, pcol += a->stride, prow++)
         *pcol = *prow;
   }
}

bool x_is_symmetric(x_matrix *a) {
   double mx, err;
   bool nonfinite;
   if (a->datatype != DT_REAL) return false;
   if (a->cols != a->rows) return false;
   if (a->cols == 0 || a->rows == 0) return true;
   ae_state_init();
   mx = 0;
   err = 0;
   nonfinite = false;
   is_symmetric_rec_diag_stat(a, 0, (ae_int_t) a->rows, &nonfinite, &mx, &err);
   if (nonfinite) return false;
   if (mx == 0) return true;
   return err / mx <= 1.0E-14;
}

bool x_is_hermitian(x_matrix *a) {
   double mx, err;
   bool nonfinite;
   if (a->datatype != DT_COMPLEX) return false;
   if (a->cols != a->rows) return false;
   if (a->cols == 0 || a->rows == 0) return true;
   ae_state_init();
   mx = 0;
   err = 0;
   nonfinite = false;
   is_hermitian_rec_diag_stat(a, 0, (ae_int_t) a->rows, &nonfinite, &mx, &err);
   if (nonfinite) return false;
   if (mx == 0) return true;
   return err / mx <= 1.0E-14;
}

bool x_force_symmetric(x_matrix *a) {
   if (a->datatype != DT_REAL) return false;
   if (a->cols != a->rows) return false;
   if (a->cols == 0 || a->rows == 0) return true;
   force_symmetric_rec_diag_stat(a, 0, (ae_int_t) a->rows);
   return true;
}

bool x_force_hermitian(x_matrix *a) {
   if (a->datatype != DT_COMPLEX) return false;
   if (a->cols != a->rows) return false;
   if (a->cols == 0 || a->rows == 0) return true;
   force_hermitian_rec_diag_stat(a, 0, (ae_int_t) a->rows);
   return true;
}

bool ae_is_symmetric(ae_matrix *a) {
   x_matrix x;
   x.owner = OWN_CALLER;
   ae_x_attach_to_matrix(&x, a);
   return x_is_symmetric(&x);
}

bool ae_is_hermitian(ae_matrix *a) {
   x_matrix x;
   x.owner = OWN_CALLER;
   ae_x_attach_to_matrix(&x, a);
   return x_is_hermitian(&x);
}

bool ae_force_symmetric(ae_matrix *a) {
   x_matrix x;
   x.owner = OWN_CALLER;
   ae_x_attach_to_matrix(&x, a);
   return x_force_symmetric(&x);
}

bool ae_force_hermitian(ae_matrix *a) {
   x_matrix x;
   x.owner = OWN_CALLER;
   ae_x_attach_to_matrix(&x, a);
   return x_force_hermitian(&x);
}

// This function converts six-bit value (from 0 to 63)  to  character  (only
// digits, lowercase and uppercase letters, minus and underscore are used).
//
// If v is negative or greater than 63, this function returns '?'.
static char _sixbits2char_tbl[64] = {
   '0', '1', '2', '3', '4', '5', '6', '7',
   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
   'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
   'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
   'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd',
   'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
   'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
   'u', 'v', 'w', 'x', 'y', 'z', '-', '_'
};

char ae_sixbits2char(ae_int_t v) {
   if (v < 0 || v > 63)
      return '?';
   return _sixbits2char_tbl[v];
// v is correct, process it
#if 0
   if (v < 10)
      return '0' + v;
   v -= 10;
   if (v < 26)
      return 'A' + v;
   v -= 26;
   if (v < 26)
      return 'a' + v;
   v -= 26;
   return v == 0 ? '-' : '_';
#endif
}

// This function converts character to six-bit value (from 0 to 63).
//
// This function is inverse of ae_sixbits2char()
// If c is not correct character, this function returns -1.
static ae_int_t _ae_char2sixbits_tbl[] = {
   -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, 62, -1, -1,
   0, 1, 2, 3, 4, 5, 6, 7,
   8, 9, -1, -1, -1, -1, -1, -1,
   -1, 10, 11, 12, 13, 14, 15, 16,
   17, 18, 19, 20, 21, 22, 23, 24,
   25, 26, 27, 28, 29, 30, 31, 32,
   33, 34, 35, -1, -1, -1, -1, 63,
   -1, 36, 37, 38, 39, 40, 41, 42,
   43, 44, 45, 46, 47, 48, 49, 50,
   51, 52, 53, 54, 55, 56, 57, 58,
   59, 60, 61, -1, -1, -1, -1, -1
};

ae_int_t ae_char2sixbits(char c) {
   return (c >= 0 && c < 127) ? _ae_char2sixbits_tbl[(int)c] : -1;
}

// This function converts three bytes (24 bits) to four six-bit values
// (24 bits again).
//
// src     pointer to three bytes
// dst     pointer to four ints
void ae_threebytes2foursixbits(const unsigned char *src, ae_int_t *dst) {
   dst[0] = src[0] & 0x3F;
   dst[1] = (src[0] >> 6) | ((src[1] & 0x0F) << 2);
   dst[2] = (src[1] >> 4) | ((src[2] & 0x03) << 4);
   dst[3] = src[2] >> 2;
}

// This function converts four six-bit values (24 bits) to three bytes
// (24 bits again).
//
// src     pointer to four ints
// dst     pointer to three bytes
void ae_foursixbits2threebytes(const ae_int_t *src, unsigned char *dst) {
   dst[0] = (unsigned char)(src[0] | ((src[1] & 0x03) << 6));
   dst[1] = (unsigned char)((src[1] >> 2) | ((src[2] & 0x0F) << 4));
   dst[2] = (unsigned char)((src[2] >> 4) | (src[3] << 2));
}

// This function serializes boolean value into buffer
//
// v           boolean value to be serialized
// buf         buffer, at least 12 characters wide
//             (11 chars for value, one for trailing zero)
void ae_bool2str(bool v, char *buf) {
   char c = v ? '1' : '0';
   ae_int_t i;
   for (i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = c;
   buf[AE_SER_ENTRY_LENGTH] = 0;
}

// This function unserializes boolean value from buffer
//
// buf         buffer which contains value; leading spaces/tabs/newlines are
//             ignored, traling spaces/tabs/newlines are treated as  end  of
//             the boolean value.
//
// This function raises an error in case unexpected symbol is found
bool ae_str2bool(const char *buf, const char **pasttheend) {
   bool was0, was1;
   const char *emsg = "ALGLIB: unable to read boolean value from stream";
   was0 = false;
   was1 = false;
   while (*buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r')
      buf++;
   while (*buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != 0) {
      if (*buf == '0') {
         was0 = true;
         buf++;
         continue;
      }
      if (*buf == '1') {
         was1 = true;
         buf++;
         continue;
      }
      ae_break(ERR_ASSERTION_FAILED, emsg);
   }
   *pasttheend = buf;
   if (!was0 && !was1)
      ae_break(ERR_ASSERTION_FAILED, emsg);
   if (was0 && was1)
      ae_break(ERR_ASSERTION_FAILED, emsg);
   return was1 ? true : false;
}

// This function serializes integer value into buffer
//
// v           integer value to be serialized
// buf         buffer, at least 12 characters wide
//             (11 chars for value, one for trailing zero)
void ae_int2str(ae_int_t v, char *buf) {
   union _u {
      ae_int_t ival;
      unsigned char bytes[9];
   } u;
   ae_int_t i;
   ae_int_t sixbits[12];
   unsigned char c;
// copy v to array of chars, sign extending it and
// converting to little endian order
//
// because we don't want to mention size of ae_int_t explicitly,
// we do it as follows:
// 1. we fill u.bytes by zeros or ones (depending on sign of v)
// 2. we copy v to u.ival
// 3. if we run on big endian architecture, we reorder u.bytes
// 4. now we have signed 64-bit representation of v stored in u.bytes
// 5. additionally, we set 9th byte of u.bytes to zero in order to
//    simplify conversion to six-bit representation
   c = v < 0 ? (unsigned char)0xFF : (unsigned char)0x00;
   u.ival = v;
   for (i = sizeof(ae_int_t); i <= 8; i++)      // <= 8 is preferred because it avoids unnecessary compiler warnings.
      u.bytes[i] = c;
   u.bytes[8] = 0;
   if (ByteOrder == AE_BIG_ENDIAN) {
      for (i = 0; i < (ae_int_t)(sizeof(ae_int_t) / 2); i++) {
         unsigned char tc;
         tc = u.bytes[i];
         u.bytes[i] = u.bytes[sizeof(ae_int_t) - 1 - i];
         u.bytes[sizeof(ae_int_t) - 1 - i] = tc;
      }
   }
// convert to six-bit representation, output
//
// NOTE: last 12th element of sixbits is always zero, we do not output it
   ae_threebytes2foursixbits(u.bytes + 0, sixbits + 0);
   ae_threebytes2foursixbits(u.bytes + 3, sixbits + 4);
   ae_threebytes2foursixbits(u.bytes + 6, sixbits + 8);
   for (i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = ae_sixbits2char(sixbits[i]);
   buf[AE_SER_ENTRY_LENGTH] = 0x00;
}

// This function serializes 64-bit integer value into buffer
//
// v           integer value to be serialized
// buf         buffer, at least 12 characters wide
//             (11 chars for value, one for trailing zero)
// state       ALGLIB environment state
static void ae_int642str(ae_int64_t v, char *buf) {
   unsigned char bytes[9];
   ae_int_t i;
   ae_int_t sixbits[12];
// copy v to array of chars, sign extending it and
// converting to little endian order
//
// because we don't want to mention size of ae_int_t explicitly,
// we do it as follows:
// 1. we fill bytes by zeros or ones (depending on sign of v)
// 2. we memmove v to bytes
// 3. if we run on big endian architecture, we reorder bytes
// 4. now we have signed 64-bit representation of v stored in bytes
// 5. additionally, we set 9th byte of bytes to zero in order to
//    simplify conversion to six-bit representation
   memset(bytes, v < 0 ? 0xFF : 0x00, 8);
   memmove(bytes, &v, 8);
   bytes[8] = 0;
   if (ByteOrder == AE_BIG_ENDIAN) {
      for (i = 0; i < (ae_int_t)(sizeof(ae_int_t) / 2); i++) {
         unsigned char tc;
         tc = bytes[i];
         bytes[i] = bytes[sizeof(ae_int_t) - 1 - i];
         bytes[sizeof(ae_int_t) - 1 - i] = tc;
      }
   }
// convert to six-bit representation, output
//
// NOTE: last 12th element of sixbits is always zero, we do not output it
   ae_threebytes2foursixbits(bytes + 0, sixbits + 0);
   ae_threebytes2foursixbits(bytes + 3, sixbits + 4);
   ae_threebytes2foursixbits(bytes + 6, sixbits + 8);
   for (i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = ae_sixbits2char(sixbits[i]);
   buf[AE_SER_ENTRY_LENGTH] = 0x00;
}

// This function unserializes integer value from string
//
// buf         buffer which contains value; leading spaces/tabs/newlines are
//             ignored, traling spaces/tabs/newlines are treated as  end  of
//             the boolean value.
//
// This function raises an error in case unexpected symbol is found
ae_int_t ae_str2int(const char *buf, const char **pasttheend) {
   const char *emsg = "ALGLIB: unable to read integer value from stream";
   ae_int_t sixbits[12];
   ae_int_t sixbitsread, i;
   union _u {
      ae_int_t ival;
      unsigned char bytes[9];
   } u;
// 1. skip leading spaces
// 2. read and decode six-bit digits
// 3. set trailing digits to zeros
// 4. convert to little endian 64-bit integer representation
// 5. convert to big endian representation, if needed
   while (*buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r')
      buf++;
   sixbitsread = 0;
   while (*buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != 0) {
      ae_int_t d;
      d = ae_char2sixbits(*buf);
      if (d < 0 || sixbitsread >= AE_SER_ENTRY_LENGTH)
         ae_break(ERR_ASSERTION_FAILED, emsg);
      sixbits[sixbitsread] = d;
      sixbitsread++;
      buf++;
   }
   *pasttheend = buf;
   if (sixbitsread == 0)
      ae_break(ERR_ASSERTION_FAILED, emsg);
   for (i = sixbitsread; i < 12; i++)
      sixbits[i] = 0;
   ae_foursixbits2threebytes(sixbits + 0, u.bytes + 0);
   ae_foursixbits2threebytes(sixbits + 4, u.bytes + 3);
   ae_foursixbits2threebytes(sixbits + 8, u.bytes + 6);
   if (ByteOrder == AE_BIG_ENDIAN) {
      for (i = 0; i < (ae_int_t)(sizeof(ae_int_t) / 2); i++) {
         unsigned char tc;
         tc = u.bytes[i];
         u.bytes[i] = u.bytes[sizeof(ae_int_t) - 1 - i];
         u.bytes[sizeof(ae_int_t) - 1 - i] = tc;
      }
   }
   return u.ival;
}

// This function unserializes 64-bit integer value from string
//
// buf         buffer which contains value; leading spaces/tabs/newlines are
//             ignored, traling spaces/tabs/newlines are treated as  end  of
//             the boolean value.
// state       ALGLIB environment state
//
// This function raises an error in case unexpected symbol is found
static ae_int64_t ae_str2int64(const char *buf, const char **pasttheend) {
   const char *emsg = "ALGLIB: unable to read integer value from stream";
   ae_int_t sixbits[12];
   ae_int_t sixbitsread, i;
   unsigned char bytes[9];
   ae_int64_t result;
// 1. skip leading spaces
// 2. read and decode six-bit digits
// 3. set trailing digits to zeros
// 4. convert to little endian 64-bit integer representation
// 5. convert to big endian representation, if needed
   while (*buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r')
      buf++;
   sixbitsread = 0;
   while (*buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != 0) {
      ae_int_t d;
      d = ae_char2sixbits(*buf);
      if (d < 0 || sixbitsread >= AE_SER_ENTRY_LENGTH)
         ae_break(ERR_ASSERTION_FAILED, emsg);
      sixbits[sixbitsread] = d;
      sixbitsread++;
      buf++;
   }
   *pasttheend = buf;
   if (sixbitsread == 0)
      ae_break(ERR_ASSERTION_FAILED, emsg);
   for (i = sixbitsread; i < 12; i++)
      sixbits[i] = 0;
   ae_foursixbits2threebytes(sixbits + 0, bytes + 0);
   ae_foursixbits2threebytes(sixbits + 4, bytes + 3);
   ae_foursixbits2threebytes(sixbits + 8, bytes + 6);
   if (ByteOrder == AE_BIG_ENDIAN) {
      for (i = 0; i < (ae_int_t)(sizeof(ae_int_t) / 2); i++) {
         unsigned char tc;
         tc = bytes[i];
         bytes[i] = bytes[sizeof(ae_int_t) - 1 - i];
         bytes[sizeof(ae_int_t) - 1 - i] = tc;
      }
   }
   memmove(&result, bytes, sizeof(result));
   return result;
}

// This function serializes double value into buffer
//
// v           double value to be serialized
// buf         buffer, at least 12 characters wide
//             (11 chars for value, one for trailing zero)
void ae_double2str(double v, char *buf) {
   union _u {
      double dval;
      unsigned char bytes[9];
   } u;
   ae_int_t i;
   ae_int_t sixbits[12];
// handle special quantities
   if (isnan(v)) {
      const char *s = ".nan_______";
      memmove(buf, s, strlen(s) + 1);
      return;
   }
   if (isposinf(v)) {
      const char *s = ".posinf____";
      memmove(buf, s, strlen(s) + 1);
      return;
   }
   if (isneginf(v)) {
      const char *s = ".neginf____";
      memmove(buf, s, strlen(s) + 1);
      return;
   }
// process general case:
// 1. copy v to array of chars
// 2. set 9th byte of u.bytes to zero in order to
//    simplify conversion to six-bit representation
// 3. convert to little endian (if needed)
// 4. convert to six-bit representation
//    (last 12th element of sixbits is always zero, we do not output it)
   u.dval = v;
   u.bytes[8] = 0;
   if (ByteOrder == AE_BIG_ENDIAN) {
      for (i = 0; i < (ae_int_t)(sizeof(double) / 2); i++) {
         unsigned char tc;
         tc = u.bytes[i];
         u.bytes[i] = u.bytes[sizeof(double) - 1 - i];
         u.bytes[sizeof(double) - 1 - i] = tc;
      }
   }
   ae_threebytes2foursixbits(u.bytes + 0, sixbits + 0);
   ae_threebytes2foursixbits(u.bytes + 3, sixbits + 4);
   ae_threebytes2foursixbits(u.bytes + 6, sixbits + 8);
   for (i = 0; i < AE_SER_ENTRY_LENGTH; i++)
      buf[i] = ae_sixbits2char(sixbits[i]);
   buf[AE_SER_ENTRY_LENGTH] = 0x00;
}

// This function unserializes double value from string
//
// buf         buffer which contains value; leading spaces/tabs/newlines are
//             ignored, traling spaces/tabs/newlines are treated as  end  of
//             the boolean value.
//
// This function raises an error in case unexpected symbol is found
double ae_str2double(const char *buf, const char **pasttheend) {
   const char *emsg = "ALGLIB: unable to read double value from stream";
   ae_int_t sixbits[12];
   ae_int_t sixbitsread, i;
   union _u { double dval; unsigned char bytes[9]; } u;
// skip leading spaces
   while (*buf == ' ' || *buf == '\t' || *buf == '\n' || *buf == '\r') buf++;
// Handle special cases
   if (*buf == '.') {
      const char *s_nan = ".nan_______";
      const char *s_posinf = ".posinf____";
      const char *s_neginf = ".neginf____";
      if (strncmp(buf, s_nan, strlen(s_nan)) == 0) { *pasttheend = buf + strlen(s_nan); return NAN; }
      if (strncmp(buf, s_posinf, strlen(s_posinf)) == 0) { *pasttheend = buf + strlen(s_posinf); return +INFINITY; }
      if (strncmp(buf, s_neginf, strlen(s_neginf)) == 0) { *pasttheend = buf + strlen(s_neginf); return -INFINITY; }
      ae_break(ERR_ASSERTION_FAILED, emsg);
   }
// General case:
// 1. read and decode six-bit digits
// 2. check that all 11 digits were read
// 3. set last 12th digit to zero (needed for simplicity of conversion)
// 4. convert to 8 bytes
// 5. convert to big endian representation, if needed
   sixbitsread = 0;
   while (*buf != ' ' && *buf != '\t' && *buf != '\n' && *buf != '\r' && *buf != 0) {
      ae_int_t d;
      d = ae_char2sixbits(*buf);
      if (d < 0 || sixbitsread >= AE_SER_ENTRY_LENGTH)
         ae_break(ERR_ASSERTION_FAILED, emsg);
      sixbits[sixbitsread] = d;
      sixbitsread++;
      buf++;
   }
   *pasttheend = buf;
   if (sixbitsread != AE_SER_ENTRY_LENGTH)
      ae_break(ERR_ASSERTION_FAILED, emsg);
   sixbits[AE_SER_ENTRY_LENGTH] = 0;
   ae_foursixbits2threebytes(sixbits + 0, u.bytes + 0);
   ae_foursixbits2threebytes(sixbits + 4, u.bytes + 3);
   ae_foursixbits2threebytes(sixbits + 8, u.bytes + 6);
   if (ByteOrder == AE_BIG_ENDIAN) {
      for (i = 0; i < (ae_int_t)(sizeof(double) / 2); i++) {
         unsigned char tc;
         tc = u.bytes[i];
         u.bytes[i] = u.bytes[sizeof(double) - 1 - i];
         u.bytes[sizeof(double) - 1 - i] = tc;
      }
   }
   return u.dval;
}

// This function performs given number of spin-wait iterations
void ae_spin_wait(ae_int_t cnt) {
// these strange operations with ae_never_change_it are necessary to
// prevent compiler optimization of the loop.
   volatile ae_int_t i;
// very unlikely because no one will wait for such amount of cycles
   if (cnt > 0x12345678)
      ae_never_change_it = cnt % 10;
// spin wait, test condition which will never be true
   for (i = 0; i < cnt; i++)
      if (ae_never_change_it > 0)
         ae_never_change_it--;
}

// This function causes the calling thread to relinquish the CPU. The thread
// is moved to the end of the queue and some other thread gets to run.
//
// NOTE: this function should NOT be called when AE_OS is AE_UNKNOWN  -  the
//       whole program will be abnormally terminated.
void ae_yield() {
#if AE_OS == AE_POSIX
   sched_yield();
#elif AE_OS == AE_WINDOWS
   if (!SwitchToThread())
      Sleep(0);
#else
   abort();
#endif
}

// This function initializes _lock structure which  is  internally  used  by
// ae_lock high-level structure.
//
// _lock structure is statically allocated, no malloc() calls  is  performed
// during its allocation. However, you have to call  _ae_free_lock_raw()  in
// order to deallocate this lock properly.
void _ae_init_lock_raw(_lock *p) {
#if AE_OS == AE_POSIX
   pthread_mutex_init(&p->mutex, NULL);
#elif AE_OS == AE_WINDOWS
   p->p_lock = (ae_int_t *) ae_align((void *)(&p->buf), AE_LOCK_ALIGNMENT);
   p->p_lock[0] = 0;
#else
   p->is_locked = false;
#endif
}

// This function acquires _lock structure.
//
// It is low-level workhorse utilized by ae_acquire_lock().
void _ae_acquire_lock_raw(_lock *p) {
#if AE_OS == AE_POSIX
   ae_int_t cnt = 0;
   while (true) {
      if (pthread_mutex_trylock(&p->mutex) == 0)
         return;
      ae_spin_wait(AE_LOCK_CYCLES);
      cnt++;
      if (cnt % AE_LOCK_TESTS_BEFORE_YIELD == 0)
         ae_yield();
   }
#elif AE_OS == AE_WINDOWS
   ae_int_t cnt = 0;
#   ifdef AE_SMP_DEBUGCOUNTERS
   InterlockedIncrement((LONG volatile *)&_ae_dbg_lock_acquisitions);
#   endif
   while (true) {
      if (InterlockedCompareExchange((LONG volatile *)p->p_lock, 1, 0) == 0)
         return;
      ae_spin_wait(AE_LOCK_CYCLES);
#   ifdef AE_SMP_DEBUGCOUNTERS
      InterlockedIncrement((LONG volatile *)&_ae_dbg_lock_spinwaits);
#   endif
      cnt++;
      if (cnt % AE_LOCK_TESTS_BEFORE_YIELD == 0) {
#   ifdef AE_SMP_DEBUGCOUNTERS
         InterlockedIncrement((LONG volatile *)&_ae_dbg_lock_yields);
#   endif
         ae_yield();
      }
   }
#else
   AE_CRITICAL_ASSERT(!p->is_locked);
   p->is_locked = true;
#endif
}

// This function releases _lock structure.
//
// It is low-level lock function which is used by ae_release_lock.
void _ae_release_lock_raw(_lock *p) {
#if AE_OS == AE_POSIX
   pthread_mutex_unlock(&p->mutex);
#elif AE_OS == AE_WINDOWS
   InterlockedExchange((LONG volatile *)p->p_lock, 0);
#else
   p->is_locked = false;
#endif
}

// This function frees _lock structure.
void _ae_free_lock_raw(_lock *p) {
#if AE_OS == AE_POSIX
   pthread_mutex_destroy(&p->mutex);
#endif
}

// This function initializes ae_lock structure.
//
// Inputs:
//     lock                -   pointer to lock structure, must be zero-filled
//     make_automatic                -   if true, lock object is added to automatic
//                             memory management list.
//
// NOTE: as a special exception, this function allows you to specify TopFr == NULL.
//       In this case all exception arising during construction
//       are handled as critical failures, with abort() being called.
//       make_automatic must be false on such calls.
void ae_init_lock(ae_lock *lock, bool make_automatic) {
   _lock *p;
   AE_CRITICAL_ASSERT(ae_check_zeros(lock, sizeof(*lock)));
   if (TopFr == NULL) {
      AE_CRITICAL_ASSERT(!make_automatic);
      ae_state_init();
      ae_init_lock(lock, false);
      ae_state_clear();
      return;
   }
   lock->eternal = false;
   ae_db_init(&lock->db, sizeof(_lock), make_automatic);
   lock->lock_ptr = lock->db.ptr;
   p = (_lock *) lock->lock_ptr;
   _ae_init_lock_raw(p);
}

// This function initializes "eternal" ae_lock structure which  is  expected
// to persist until the end of the execution of the program.  Eternal  locks
// can not be deallocated (cleared) and  do  not  increase debug  allocation
// counters.  Errors  during  allocation  of eternal  locks  are  considered
// critical exceptions and handled by calling abort().
//
// Inputs:
//     lock                -   pointer to lock structure, must be zero-filled
//     make_automatic                -   if true, lock object is added to automatic
//                             memory management list.
void ae_init_lock_eternal(ae_lock *lock) {
   _lock *p;
   AE_CRITICAL_ASSERT(ae_check_zeros(lock, sizeof(*lock)));
   lock->eternal = true;
   lock->lock_ptr = eternal_malloc(sizeof(_lock));
   p = (_lock *) lock->lock_ptr;
   _ae_init_lock_raw(p);
}

// This function acquires lock. In case lock is busy, we perform several
// iterations inside tight loop before trying again.
void ae_acquire_lock(ae_lock *lock) {
   _lock *p;
   p = (_lock *) lock->lock_ptr;
   _ae_acquire_lock_raw(p);
}

// This function releases lock.
void ae_release_lock(ae_lock *lock) {
   _lock *p;
   p = (_lock *) lock->lock_ptr;
   _ae_release_lock_raw(p);
}

// This function frees ae_lock structure.
void ae_free_lock(ae_lock *lock) {
   _lock *p;
   AE_CRITICAL_ASSERT(!lock->eternal);
   p = (_lock *) lock->lock_ptr;
   if (p != NULL)
      _ae_free_lock_raw(p);
   ae_db_free(&lock->db);
}

static void ae_shared_pool_destroy(void *_dst) { ae_shared_pool_free(_dst, false); }

// This function creates ae_shared_pool structure.
//
// dst                 destination shared pool, must be zero-filled
//                     already allocated, but not initialized.
// make_automatic                if true, vector will be registered in the frame of the current state;
//
// Error handling:
// * on failure calls ae_break() with TopFr == NULL. Usually it results
//   in abort() call.
//
// dst is assumed to be uninitialized, its fields are ignored.
void ae_shared_pool_init(void *_dst, bool make_automatic) {
   ae_shared_pool *dst;
   dst = (ae_shared_pool *) _dst;
   AE_CRITICAL_ASSERT(ae_check_zeros(dst, sizeof(*dst)));
// init
   dst->seed_object = NULL;
   dst->recycled_objects = NULL;
   dst->recycled_entries = NULL;
   dst->enumeration_counter = NULL;
   dst->size_of_object = 0;
   dst->init = NULL;
   dst->copy = NULL;
   dst->free = NULL;
   dst->frame_entry.deallocator = ae_shared_pool_destroy;
   dst->frame_entry.ptr = dst;
   if (make_automatic) ae_db_attach(&dst->frame_entry);
   ae_init_lock(&dst->pool_lock, false);
}

// This function clears all dynamically allocated fields of the pool except
// for the lock. It does NOT try to acquire pool_lock.
//
// NOTE: this function is NOT thread-safe, it is not protected by lock.
static void ae_shared_pool_internalclear(ae_shared_pool *dst, bool make_automatic) {
   ae_shared_pool_entry *ptr, *tmp;
// free seed
   if (dst->seed_object != NULL) {
      dst->free((void *)dst->seed_object, make_automatic);
      ae_free((void *)dst->seed_object);
      dst->seed_object = NULL;
   }
// free recycled objects
   for (ptr = dst->recycled_objects; ptr != NULL;) {
      tmp = (ae_shared_pool_entry *) ptr->next_entry;
      dst->free(ptr->obj, make_automatic);
      ae_free(ptr->obj);
      ae_free(ptr);
      ptr = tmp;
   }
   dst->recycled_objects = NULL;
// free recycled entries
   for (ptr = dst->recycled_entries; ptr != NULL;) {
      tmp = (ae_shared_pool_entry *) ptr->next_entry;
      ae_free(ptr);
      ptr = tmp;
   }
   dst->recycled_entries = NULL;
}

// This function creates copy of ae_shared_pool.
//
// dst                 destination pool, must be zero-filled
// src                 source pool
// make_automatic                if true, vector will be registered in the frame of the current state;
//
// dst is assumed to be uninitialized, its fields are ignored.
//
// NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
//       you should NOT call it when lock can be used by another thread.
void ae_shared_pool_copy(void *_dst, void *_src, bool make_automatic) {
   ae_shared_pool *dst, *src;
   ae_shared_pool_entry *ptr;
   dst = (ae_shared_pool *) _dst;
   src = (ae_shared_pool *) _src;
   ae_shared_pool_init(dst, make_automatic);
// copy non-pointer fields
   dst->size_of_object = src->size_of_object;
   dst->init = src->init;
   dst->copy = src->copy;
   dst->free = src->free;
// copy seed object
   if (src->seed_object != NULL) {
      dst->seed_object = ae_malloc(dst->size_of_object);
      memset(dst->seed_object, 0, dst->size_of_object);
      dst->copy(dst->seed_object, src->seed_object, false);
   }
// copy recycled objects
   dst->recycled_objects = NULL;
   for (ptr = src->recycled_objects; ptr != NULL; ptr = (ae_shared_pool_entry *) ptr->next_entry) {
      ae_shared_pool_entry *tmp;
   // allocate entry, immediately add to the recycled list (we do not want to lose it in case of future malloc failures).
      tmp = (ae_shared_pool_entry *) ae_malloc(sizeof(ae_shared_pool_entry));
      memset(tmp, 0, sizeof(*tmp));
      tmp->next_entry = dst->recycled_objects;
      dst->recycled_objects = tmp;
   // prepare place for object, copy() it
      tmp->obj = ae_malloc(dst->size_of_object);
      memset(tmp->obj, 0, dst->size_of_object);
      dst->copy(tmp->obj, ptr->obj, false);
   }
// recycled entries are not copied because they do not store any information
   dst->recycled_entries = NULL;
// enumeration counter is reset on copying
   dst->enumeration_counter = NULL;
// initialize frame record
   dst->frame_entry.deallocator = ae_shared_pool_destroy;
   dst->frame_entry.ptr = dst;
}

// This function performs destruction of the pool object.
//
// NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
//       you should NOT call it when pool can be used by another thread.
void ae_shared_pool_free(void *_dst, bool make_automatic) {
   ae_shared_pool *dst = (ae_shared_pool *) _dst;
// clear seed and lists
   ae_shared_pool_internalclear(dst, make_automatic);
// clear fields
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

// This function returns True, if internal seed object was set.  It  returns
// False for un-seeded pool.
//
// dst                 destination pool (initialized by constructor function)
//
// NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
//       you should NOT call it when lock can be used by another thread.
bool ae_shared_pool_is_initialized(void *_dst) {
   ae_shared_pool *dst = (ae_shared_pool *) _dst;
   return dst->seed_object != NULL;
}

// This function sets internal seed object. All objects owned by the pool
// (current seed object, recycled objects) are automatically freed.
//
// dst                 destination pool (initialized by constructor function)
// seed_object         new seed object
// size_of_object      sizeof(), used to allocate memory
// init                constructor function
// copy                copy constructor
// clear               destructor function
//
// NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
//       you should NOT call it when lock can be used by another thread.
void ae_shared_pool_set_seed(ae_shared_pool *dst, void *seed_object, ae_int_t size_of_object, void (*init)(void *dst, bool make_automatic), void (*copy)(void *dst, void *src, bool make_automatic), void (*free)(void *ptr, bool make_automatic)) {
// free internal objects
   ae_shared_pool_internalclear(dst, false);
// set non-pointer fields
   dst->size_of_object = size_of_object;
   dst->init = init;
   dst->copy = copy;
   dst->free = free;
// set seed object
   dst->seed_object = ae_malloc(size_of_object);
   memset(dst->seed_object, 0, size_of_object);
   copy(dst->seed_object, seed_object, false);
}

// This  function  retrieves  a  copy  of  the seed object from the pool and
// stores it to target smart pointer ptr.
//
// In case target pointer owns non-NULL  value,  it  is  deallocated  before
// storing value retrieved from pool. Target pointer becomes  owner  of  the
// value which was retrieved from pool.
//
// pool                pool
// pptr                pointer to ae_smart_ptr structure
//
// NOTE: this function IS thread-safe.  It  acquires  pool  lock  during its
//       operation and can be used simultaneously from several threads.
void ae_shared_pool_retrieve(ae_shared_pool *pool, ae_smart_ptr *pptr) {
   void *new_obj;
// assert that pool was seeded
   ae_assert(pool->seed_object != NULL, "ALGLIB: shared pool is not seeded, PoolRetrieve() failed");
// acquire lock
   ae_acquire_lock(&pool->pool_lock);
// try to reuse recycled objects
   if (pool->recycled_objects != NULL) {
      ae_shared_pool_entry *result;
   // retrieve entry/object from list of recycled objects
      result = pool->recycled_objects;
      pool->recycled_objects = (ae_shared_pool_entry *) pool->recycled_objects->next_entry;
      new_obj = result->obj;
      result->obj = NULL;
   // move entry to list of recycled entries
      result->next_entry = pool->recycled_entries;
      pool->recycled_entries = result;
   // release lock
      ae_release_lock(&pool->pool_lock);
   // assign object to smart pointer
      ae_smart_ptr_assign(pptr, new_obj, true, true, pool->free);
      return;
   }
// release lock; we do not need it anymore because copy constructor does not modify source variable
   ae_release_lock(&pool->pool_lock);
// create new object from seed, immediately assign object to smart pointer (do not want to lose it in case of future failures).
   new_obj = ae_malloc(pool->size_of_object);
   memset(new_obj, 0, pool->size_of_object);
   ae_smart_ptr_assign(pptr, new_obj, true, true, pool->free);
// perform actual copying; before this line smartptr points to zero-filled instance
   pool->copy(new_obj, pool->seed_object, false);
}

// This function recycles object owned by smart  pointer  by  moving  it  to
// internal storage of the shared pool.
//
// Source pointer must own the object. After function is over, it owns NULL
// pointer.
//
// pool                pool
// pptr                pointer to ae_smart_ptr structure
//
// NOTE: this function IS thread-safe.  It  acquires  pool  lock  during its
//       operation and can be used simultaneously from several threads.
void ae_shared_pool_recycle(ae_shared_pool *pool, ae_smart_ptr *pptr) {
   ae_shared_pool_entry *new_entry;
// assert that pool was seeded
   ae_assert(pool->seed_object != NULL, "ALGLIB: shared pool is not seeded, PoolRecycle() failed");
// assert that pointer non-null and owns the object
   ae_assert(pptr->is_owner, "ALGLIB: pptr in ae_shared_pool_recycle() does not own its pointer");
   ae_assert(pptr->ptr != NULL, "ALGLIB: pptr in ae_shared_pool_recycle() is NULL");
// acquire lock
   ae_acquire_lock(&pool->pool_lock);
// acquire shared pool entry (reuse one from recycled_entries or allocate new one)
   if (pool->recycled_entries != NULL) {
   // reuse previously allocated entry
      new_entry = pool->recycled_entries;
      pool->recycled_entries = (ae_shared_pool_entry *) new_entry->next_entry;
   } else {
   // Allocate memory for new entry.
   //
   // NOTE: we release pool lock during allocation because ae_malloc() may raise
   //       exception and we do not want our pool to be left in the locked state.
      ae_release_lock(&pool->pool_lock);
      new_entry = (ae_shared_pool_entry *) ae_malloc(sizeof(ae_shared_pool_entry));
      ae_acquire_lock(&pool->pool_lock);
   }
// add object to the list of recycled objects
   new_entry->obj = pptr->ptr;
   new_entry->next_entry = pool->recycled_objects;
   pool->recycled_objects = new_entry;
// release lock object
   ae_release_lock(&pool->pool_lock);
// release source pointer
   ae_smart_ptr_release(pptr);
}

// This function clears internal list of  recycled  objects,  but  does  not
// change seed object managed by the pool.
//
// pool                pool
//
// NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
//       you should NOT call it when lock can be used by another thread.
void ae_shared_pool_clear_recycled(ae_shared_pool *pool, bool make_automatic) {
   ae_shared_pool_entry *ptr, *tmp;
// clear recycled objects
   for (ptr = pool->recycled_objects; ptr != NULL;) {
      tmp = (ae_shared_pool_entry *) ptr->next_entry;
      pool->free(ptr->obj, make_automatic);
      ae_free(ptr->obj);
      ae_free(ptr);
      ptr = tmp;
   }
   pool->recycled_objects = NULL;
}

// This function allows to enumerate recycled elements of the  shared  pool.
// It stores pointer to the first recycled object in the smart pointer.
//
// IMPORTANT:
// * in case target pointer owns non-NULL  value,  it  is deallocated before
//   storing value retrieved from pool.
// * recycled object IS NOT removed from pool
// * target pointer DOES NOT become owner of the new value
// * this function IS NOT thread-safe
// * you SHOULD NOT modify shared pool during enumeration (although you  can
//   modify state of the objects retrieved from pool)
// * in case there is no recycled objects in the pool, NULL is stored to pptr
// * in case pool is not seeded, NULL is stored to pptr
//
// pool                pool
// pptr                pointer to ae_smart_ptr structure
void ae_shared_pool_first_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr) {
// modify internal enumeration counter
   pool->enumeration_counter = pool->recycled_objects;
// exit on empty list
   if (pool->enumeration_counter == NULL) {
      ae_smart_ptr_assign(pptr, NULL, false, false, NULL);
      return;
   }
// assign object to smart pointer
   ae_smart_ptr_assign(pptr, pool->enumeration_counter->obj, false, false, pool->free);
}

// This function allows to enumerate recycled elements of the  shared  pool.
// It stores pointer to the next recycled object in the smart pointer.
//
// IMPORTANT:
// * in case target pointer owns non-NULL  value,  it  is deallocated before
//   storing value retrieved from pool.
// * recycled object IS NOT removed from pool
// * target pointer DOES NOT become owner of the new value
// * this function IS NOT thread-safe
// * you SHOULD NOT modify shared pool during enumeration (although you  can
//   modify state of the objects retrieved from pool)
// * in case there is no recycled objects left in the pool, NULL is stored.
// * in case pool is not seeded, NULL is stored.
//
// pool                pool
// pptr                pointer to ae_smart_ptr structure
void ae_shared_pool_next_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr) {
// exit on end of list
   if (pool->enumeration_counter == NULL) {
      ae_smart_ptr_assign(pptr, NULL, false, false, NULL);
      return;
   }
// modify internal enumeration counter
   pool->enumeration_counter = (ae_shared_pool_entry *) pool->enumeration_counter->next_entry;
// exit on empty list
   if (pool->enumeration_counter == NULL) {
      ae_smart_ptr_assign(pptr, NULL, false, false, NULL);
      return;
   }
// assign object to smart pointer
   ae_smart_ptr_assign(pptr, pool->enumeration_counter->obj, false, false, pool->free);
}

// This function clears internal list of recycled objects and  seed  object.
// However, pool still can be used (after initialization with another seed).
//
// pool                pool
//
// NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
//       you should NOT call it when lock can be used by another thread.
void ae_shared_pool_reset(ae_shared_pool *pool) {
// clear seed and lists
   ae_shared_pool_internalclear(pool, false);
// clear fields
   pool->seed_object = NULL;
   pool->recycled_objects = NULL;
   pool->recycled_entries = NULL;
   pool->enumeration_counter = NULL;
   pool->size_of_object = 0;
   pool->init = NULL;
   pool->copy = NULL;
   pool->free = NULL;
}

// This function initializes serializer
void ae_serializer_init(ae_serializer *serializer) {
   serializer->mode = AE_SM_DEFAULT;
   serializer->entries_needed = 0;
   serializer->bytes_asked = 0;
}

void ae_serializer_clear(ae_serializer *serializer) {
}

void ae_serializer_alloc_start(ae_serializer *serializer) {
   serializer->entries_needed = 0;
   serializer->bytes_asked = 0;
   serializer->mode = AE_SM_ALLOC;
}

void ae_serializer_alloc_entry(ae_serializer *serializer) {
   serializer->entries_needed++;
}

void ae_serializer_alloc_byte_array(ae_serializer *serializer, ae_vector *bytes) {
   ae_int_t n;
   n = bytes->cnt;
   n = n / 8 + (n % 8 > 0 ? 1 : 0);
   serializer->entries_needed += 1 + n;
}

// After allocation phase is done, this function returns  required  size  of
// the output string buffer (including trailing zero symbol). Actual size of
// the data being stored can be a few characters smaller than requested.
ae_int_t ae_serializer_get_alloc_size(ae_serializer *serializer) {
   ae_int_t rows, lastrowsize, result;
   serializer->mode = AE_SM_READY2S;
// if no entries needes (degenerate case)
   if (serializer->entries_needed == 0) {
      serializer->bytes_asked = 4;      // A pair of chars for \r\n, one for dot, one for trailing zero.
      return serializer->bytes_asked;
   }
// non-degenerate case
   rows = serializer->entries_needed / AE_SER_ENTRIES_PER_ROW;
   lastrowsize = AE_SER_ENTRIES_PER_ROW;
   if (serializer->entries_needed % AE_SER_ENTRIES_PER_ROW) {
      lastrowsize = serializer->entries_needed % AE_SER_ENTRIES_PER_ROW;
      rows++;
   }
// calculate result size
// Data size + ' ' symbols + '\n' symbols + trailing '.' + trailing '\0'.
   result = ((rows - 1) * AE_SER_ENTRIES_PER_ROW + lastrowsize) * AE_SER_ENTRY_LENGTH;
   result += (rows - 1) * (AE_SER_ENTRIES_PER_ROW - 1) + (lastrowsize - 1);
   result += rows * 2;
   result += 1;
   result += 1;
   serializer->bytes_asked = result;
   return result;
}

#ifdef AE_USE_CPP_SERIALIZATION
void ae_serializer_sstart_str(ae_serializer *serializer, std::string *buf) {
   serializer->mode = AE_SM_TO_CPPSTRING;
   serializer->out_cppstr = buf;
   serializer->entries_saved = 0;
   serializer->bytes_written = 0;
}

void ae_serializer_ustart_str(ae_serializer *serializer, const std::string *buf) {
   serializer->mode = AE_SM_FROM_STRING;
   serializer->in_str = buf->c_str();
}

static char cpp_writer(const char *p_string, ae_int_t aux) {
   std::ostream * stream = reinterpret_cast < std::ostream * >(aux);
   stream->write(p_string, strlen(p_string));
   return stream->bad()? 1 : 0;
}

static char cpp_reader(ae_int_t aux, ae_int_t cnt, char *p_buf) {
   std::istream * stream = reinterpret_cast < std::istream * >(aux);
   int c;
   if (cnt <= 0)
      return 1; // unexpected cnt.
   while (true) {
      c = stream->get();
      if (c < 0 || c > 255)
         return 1;      // failure!
      if (c != ' ' && c != '\t' && c != '\n' && c != '\r')
         break;
   }
   p_buf[0] = (char)c;
   for (int k = 1; k < cnt; k++) {
      c = stream->get();
      if (c < 0 || c > 255 || c == ' ' || c == '\t' || c == '\n' || c == '\r')
         return 1;      // failure!
      p_buf[k] = (char)c;
   }
   p_buf[cnt] = 0;
   return 0;    // success
}

void ae_serializer_sstart_stream(ae_serializer *serializer, std::ostream *stream) {
   serializer->mode = AE_SM_TO_STREAM;
   serializer->stream_writer = cpp_writer;
   serializer->stream_aux = reinterpret_cast < ae_int_t > (stream);
   serializer->entries_saved = 0;
   serializer->bytes_written = 0;
}

void ae_serializer_ustart_stream(ae_serializer *serializer, const std::istream *stream) {
   serializer->mode = AE_SM_FROM_STREAM;
   serializer->stream_reader = cpp_reader;
   serializer->stream_aux = reinterpret_cast < ae_int_t > (stream);
}
#endif

void ae_serializer_sstart_str(ae_serializer *serializer, char *buf) {
   serializer->mode = AE_SM_TO_STRING;
   serializer->out_str = buf;
   serializer->out_str[0] = 0;
   serializer->entries_saved = 0;
   serializer->bytes_written = 0;
}

void ae_serializer_ustart_str(ae_serializer *serializer, const char *buf) {
   serializer->mode = AE_SM_FROM_STRING;
   serializer->in_str = buf;
}

void ae_serializer_sstart_stream(ae_serializer *serializer, ae_stream_writer writer, ae_int_t aux) {
   serializer->mode = AE_SM_TO_STREAM;
   serializer->stream_writer = writer;
   serializer->stream_aux = aux;
   serializer->entries_saved = 0;
   serializer->bytes_written = 0;
}

void ae_serializer_ustart_stream(ae_serializer *serializer, ae_stream_reader reader, ae_int_t aux) {
   serializer->mode = AE_SM_FROM_STREAM;
   serializer->stream_reader = reader;
   serializer->stream_aux = aux;
}

void ae_serializer_serialize_bool(ae_serializer *serializer, bool v) {
   char buf[AE_SER_ENTRY_LENGTH + 3];
   const char *emsg = "ALGLIB: serialization integrity error";
   ae_int_t bytes_appended;
// prepare serialization, check consistency
   ae_bool2str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   bytes_appended = (ae_int_t) strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // strict "<" because we need space for trailing zero.
   serializer->bytes_written += bytes_appended;
// append to buffer
#ifdef AE_USE_CPP_SERIALIZATION
   if (serializer->mode == AE_SM_TO_CPPSTRING) {
      *(serializer->out_cppstr) += buf;
      return;
   }
#endif
   if (serializer->mode == AE_SM_TO_STRING) {
      strcat(serializer->out_str, buf);
      serializer->out_str += bytes_appended;
      return;
   }
   if (serializer->mode == AE_SM_TO_STREAM) {
      ae_assert(serializer->stream_writer(buf, serializer->stream_aux) == 0, "serializer: error writing to stream");
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, emsg);
}

void ae_serializer_serialize_int(ae_serializer *serializer, ae_int_t v) {
   char buf[AE_SER_ENTRY_LENGTH + 3];
   const char *emsg = "ALGLIB: serialization integrity error";
   ae_int_t bytes_appended;
// prepare serialization, check consistency
   ae_int2str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   bytes_appended = (ae_int_t) strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // strict "<" because we need space for trailing zero.
   serializer->bytes_written += bytes_appended;
// append to buffer
#ifdef AE_USE_CPP_SERIALIZATION
   if (serializer->mode == AE_SM_TO_CPPSTRING) {
      *(serializer->out_cppstr) += buf;
      return;
   }
#endif
   if (serializer->mode == AE_SM_TO_STRING) {
      strcat(serializer->out_str, buf);
      serializer->out_str += bytes_appended;
      return;
   }
   if (serializer->mode == AE_SM_TO_STREAM) {
      ae_assert(serializer->stream_writer(buf, serializer->stream_aux) == 0, "serializer: error writing to stream");
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, emsg);
}

void ae_serializer_serialize_int64(ae_serializer *serializer, ae_int64_t v) {
   char buf[AE_SER_ENTRY_LENGTH + 3];
   const char *emsg = "ALGLIB: serialization integrity error";
   ae_int_t bytes_appended;
// prepare serialization, check consistency
   ae_int642str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   bytes_appended = (ae_int_t) strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // strict "<" because we need space for trailing zero.
   serializer->bytes_written += bytes_appended;
// append to buffer
#ifdef AE_USE_CPP_SERIALIZATION
   if (serializer->mode == AE_SM_TO_CPPSTRING) {
      *(serializer->out_cppstr) += buf;
      return;
   }
#endif
   if (serializer->mode == AE_SM_TO_STRING) {
      strcat(serializer->out_str, buf);
      serializer->out_str += bytes_appended;
      return;
   }
   if (serializer->mode == AE_SM_TO_STREAM) {
      ae_assert(serializer->stream_writer(buf, serializer->stream_aux) == 0, "serializer: error writing to stream");
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, emsg);
}

void ae_serializer_serialize_double(ae_serializer *serializer, double v) {
   char buf[AE_SER_ENTRY_LENGTH + 3];
   const char *emsg = "ALGLIB: serialization integrity error";
   ae_int_t bytes_appended;
// prepare serialization, check consistency
   ae_double2str(v, buf);
   serializer->entries_saved++;
   if (serializer->entries_saved % AE_SER_ENTRIES_PER_ROW)
      strcat(buf, " ");
   else
      strcat(buf, "\r\n");
   bytes_appended = (ae_int_t) strlen(buf);
   ae_assert(serializer->bytes_written + bytes_appended < serializer->bytes_asked, emsg); // strict "<" because we need space for trailing zero.
   serializer->bytes_written += bytes_appended;
// append to buffer
#ifdef AE_USE_CPP_SERIALIZATION
   if (serializer->mode == AE_SM_TO_CPPSTRING) {
      *(serializer->out_cppstr) += buf;
      return;
   }
#endif
   if (serializer->mode == AE_SM_TO_STRING) {
      strcat(serializer->out_str, buf);
      serializer->out_str += bytes_appended;
      return;
   }
   if (serializer->mode == AE_SM_TO_STREAM) {
      ae_assert(serializer->stream_writer(buf, serializer->stream_aux) == 0, "serializer: error writing to stream");
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, emsg);
}

void ae_serializer_serialize_byte_array(ae_serializer *serializer, ae_vector *bytes) {
   ae_int_t chunk_size, entries_count;
   chunk_size = 8;
// save array length
   ae_serializer_serialize_int(serializer, bytes->cnt);
// determine entries count
   entries_count = bytes->cnt / chunk_size + (bytes->cnt % chunk_size > 0 ? 1 : 0);
   for (ae_int_t eidx = 0; eidx < entries_count; eidx++) {
      ae_int64_t tmpi;
      ae_int_t elen;
      elen = bytes->cnt - eidx * chunk_size;
      elen = elen > chunk_size ? chunk_size : elen;
      memset(&tmpi, 0, sizeof(tmpi));
      memmove(&tmpi, bytes->ptr.p_ubyte + eidx * chunk_size, elen);
      ae_serializer_serialize_int64(serializer, tmpi);
   }
}

void ae_serializer_unserialize_bool(ae_serializer *serializer, bool *v) {
   if (serializer->mode == AE_SM_FROM_STRING) {
      *v = ae_str2bool(serializer->in_str, &serializer->in_str);
      return;
   }
   if (serializer->mode == AE_SM_FROM_STREAM) {
      char buf[AE_SER_ENTRY_LENGTH + 3];
      const char *p = buf;
      ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf) == 0, "serializer: error reading from stream");
      *v = ae_str2bool(buf, &p);
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, "ae_serializer: integrity check failed");
}

void ae_serializer_unserialize_int(ae_serializer *serializer, ae_int_t *v) {
   if (serializer->mode == AE_SM_FROM_STRING) {
      *v = ae_str2int(serializer->in_str, &serializer->in_str);
      return;
   }
   if (serializer->mode == AE_SM_FROM_STREAM) {
      char buf[AE_SER_ENTRY_LENGTH + 3];
      const char *p = buf;
      ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf) == 0, "serializer: error reading from stream");
      *v = ae_str2int(buf, &p);
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, "ae_serializer: integrity check failed");
}

void ae_serializer_unserialize_int64(ae_serializer *serializer, ae_int64_t *v) {
   if (serializer->mode == AE_SM_FROM_STRING) {
      *v = ae_str2int64(serializer->in_str, &serializer->in_str);
      return;
   }
   if (serializer->mode == AE_SM_FROM_STREAM) {
      char buf[AE_SER_ENTRY_LENGTH + 3];
      const char *p = buf;
      ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf) == 0, "serializer: error reading from stream");
      *v = ae_str2int64(buf, &p);
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, "ae_serializer: integrity check failed");
}

void ae_serializer_unserialize_double(ae_serializer *serializer, double *v) {
   if (serializer->mode == AE_SM_FROM_STRING) {
      *v = ae_str2double(serializer->in_str, &serializer->in_str);
      return;
   }
   if (serializer->mode == AE_SM_FROM_STREAM) {
      char buf[AE_SER_ENTRY_LENGTH + 3];
      const char *p = buf;
      ae_assert(serializer->stream_reader(serializer->stream_aux, AE_SER_ENTRY_LENGTH, buf) == 0, "serializer: error reading from stream");
      *v = ae_str2double(buf, &p);
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, "ae_serializer: integrity check failed");
}

void ae_serializer_unserialize_byte_array(ae_serializer *serializer, ae_vector *bytes) {
   ae_int_t chunk_size, n, entries_count;
   chunk_size = 8;
// read array length, allocate output
   ae_serializer_unserialize_int(serializer, &n);
   ae_vector_set_length(bytes, n);
// determine entries count, read entries
   entries_count = n / chunk_size + (n % chunk_size > 0 ? 1 : 0);
   for (ae_int_t eidx = 0; eidx < entries_count; eidx++) {
      ae_int_t elen;
      ae_int64_t tmp64;
      elen = n - eidx * chunk_size;
      elen = elen > chunk_size ? chunk_size : elen;
      ae_serializer_unserialize_int64(serializer, &tmp64);
      memmove(bytes->ptr.p_ubyte + eidx * chunk_size, &tmp64, elen);
   }
}

void ae_serializer_stop(ae_serializer *serializer) {
#ifdef AE_USE_CPP_SERIALIZATION
   if (serializer->mode == AE_SM_TO_CPPSTRING) {
      ae_assert(serializer->bytes_written + 1 < serializer->bytes_asked, "ae_serializer: integrity check failed"); // strict "<" because we need space for trailing zero.
      serializer->bytes_written++;
      *(serializer->out_cppstr) += ".";
      return;
   }
#endif
   if (serializer->mode == AE_SM_TO_STRING) {
      ae_assert(serializer->bytes_written + 1 < serializer->bytes_asked, "ae_serializer: integrity check failed"); // strict "<" because we need space for trailing zero.
      serializer->bytes_written++;
      strcat(serializer->out_str, ".");
      serializer->out_str += 1;
      return;
   }
   if (serializer->mode == AE_SM_TO_STREAM) {
      ae_assert(serializer->bytes_written + 1 < serializer->bytes_asked, "ae_serializer: integrity check failed"); // strict "<" because we need space for trailing zero.
      serializer->bytes_written++;
      ae_assert(serializer->stream_writer(".", serializer->stream_aux) == 0, "ae_serializer: error writing to stream");
      return;
   }
   if (serializer->mode == AE_SM_FROM_STRING) {
   // because input string may be from pre-3.11 serializer,
   // which does not include trailing dot, we do not test
   // string for presence of "." symbol. Anyway, because string
   // is not stream, we do not have to read ALL trailing symbols.
      return;
   }
   if (serializer->mode == AE_SM_FROM_STREAM) {
   // Read trailing dot, perform integrity check
      char buf[2];
      ae_assert(serializer->stream_reader(serializer->stream_aux, 1, buf) == 0, "ae_serializer: error reading from stream");
      ae_assert(buf[0] == '.', "ae_serializer: trailing . is not found in the stream");
      return;
   }
   ae_break(ERR_ASSERTION_FAILED, "ae_serializer: integrity check failed");
}

// Complex math functions
ae_complex ae_c_neg(ae_complex lhs) {
   return ae_complex_from_d(-lhs.x, -lhs.y);
}

ae_complex ae_c_conj(ae_complex lhs) {
   return ae_complex_from_d(+lhs.x, -lhs.y);
}

ae_complex ae_c_sqr(ae_complex lhs) {
   double x = lhs.x, y = lhs.y;
   return ae_complex_from_d(x * x - y * y, 2.0 * x * y);
}

double ae_c_abs(ae_complex z) {
   double w;
   double xabs;
   double yabs;
   double v;
   xabs = fabs(z.x);
   yabs = fabs(z.y);
   w = xabs > yabs ? xabs : yabs;
   v = xabs < yabs ? xabs : yabs;
   if (v == 0)
      return w;
   else {
      double t = v / w;
      return w * sqrt(1 + t * t);
   }
}

bool ae_c_eq(ae_complex lhs, ae_complex rhs) {
   volatile double x1 = lhs.x;
   volatile double x2 = rhs.x;
   volatile double y1 = lhs.y;
   volatile double y2 = rhs.y;
   return x1 == x2 && y1 == y2;
}

bool ae_c_neq(ae_complex lhs, ae_complex rhs) {
   volatile double x1 = lhs.x;
   volatile double x2 = rhs.x;
   volatile double y1 = lhs.y;
   volatile double y2 = rhs.y;
   return x1 != x2 || y1 != y2;
}

ae_complex ae_c_add(ae_complex lhs, ae_complex rhs) {
   return ae_complex_from_d(lhs.x + rhs.x, lhs.y + rhs.y);
}

ae_complex ae_c_mul(ae_complex lhs, ae_complex rhs) {
   double x1 = lhs.x, y1 = lhs.y, x2 = rhs.x, y2 = rhs.y;
   return ae_complex_from_d(x1 * x2 - y1 * y2, x1 * y2 + y1 * x2);
}

ae_complex ae_c_sub(ae_complex lhs, ae_complex rhs) {
   return ae_complex_from_d(lhs.x - rhs.x, lhs.y - rhs.y);
}

ae_complex ae_c_div(ae_complex lhs, ae_complex rhs) {
   double x1 = lhs.x, y1 = lhs.y, x2 = rhs.x, y2 = rhs.y;
   double e;
   double f;
   if (fabs(y2) < fabs(x2)) {
      e = y2 / x2;
      f = x2 + y2 * e;
      return ae_complex_from_d((x1 + y1 * e) / f, (y1 - x1 * e) / f);
   } else {
      e = x2 / y2;
      f = y2 + x2 * e;
      return ae_complex_from_d((x1 * e + y1) / f, (y1 * e - x1) / f);
   }
}

bool ae_c_eq_d(ae_complex lhs, double rhs) {
   volatile double x1 = lhs.x;
   volatile double x2 = rhs;
   volatile double y1 = lhs.y;
   volatile double y2 = 0.0;
   return x1 == x2 && y1 == y2;
}

bool ae_c_neq_d(ae_complex lhs, double rhs) {
   volatile double x1 = lhs.x;
   volatile double x2 = rhs;
   volatile double y1 = lhs.y;
   volatile double y2 = 0.0;
   return x1 != x2 || y1 != y2;
}

ae_complex ae_c_add_d(ae_complex lhs, double rhs) {
   return ae_complex_from_d(lhs.x + rhs, lhs.y);
}

ae_complex ae_c_mul_d(ae_complex lhs, double rhs) {
   return ae_complex_from_d(lhs.x * rhs, lhs.y * rhs);
}

ae_complex ae_c_sub_d(ae_complex lhs, double rhs) {
   return ae_complex_from_d(lhs.x - rhs, lhs.y);
}

ae_complex ae_c_d_sub(double lhs, ae_complex rhs) {
   return ae_complex_from_d(lhs - rhs.x, -rhs.y);
}

ae_complex ae_c_div_d(ae_complex lhs, double rhs) {
   return ae_complex_from_d(lhs.x / rhs, lhs.y / rhs);
}

ae_complex ae_c_d_div(double lhs, ae_complex rhs) {
   double x2 = rhs.x, y2 = rhs.y;
   double e;
   double f;
   if (fabs(y2) < fabs(x2)) {
      e = y2 / x2;
      f = x2 + y2 * e;
      return ae_complex_from_d(lhs / f, -lhs * e / f);
   } else {
      e = x2 / y2;
      f = y2 + x2 * e;
      return ae_complex_from_d(lhs * e / f, -lhs / f);
   }
}

// Complex BLAS operations
ae_complex ae_v_cdotproduct(const ae_complex *v0, ae_int_t stride0, const char *conj0, const ae_complex *v1, ae_int_t stride1, const char *conj1, ae_int_t n) {
   double rx = 0.0, ry = 0.0;
   ae_int_t i;
   bool bconj0 = !(conj0[0] == 'N' || conj0[0] == 'n');
   bool bconj1 = !(conj1[0] == 'N' || conj1[0] == 'n');
   if (bconj0 && bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = -v0->y;
         v1x = v1->x;
         v1y = -v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   if (!bconj0 && bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = v0->y;
         v1x = v1->x;
         v1y = -v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   if (bconj0 && !bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = -v0->y;
         v1x = v1->x;
         v1y = v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   if (!bconj0 && !bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = v0->y;
         v1x = v1->x;
         v1y = v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   return ae_complex_from_d(rx, ry);
}

void ae_v_cmove(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = vsrc->x;
            vdst->y = -vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
            *vdst = *vsrc;
      }
   } else {
   // optimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = vsrc->x;
            vdst->y = -vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++)
            *vdst = *vsrc;
      }
   }
}

void ae_v_cmoveneg(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = -vsrc->x;
            vdst->y = vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = -vsrc->x;
            vdst->y = -vsrc->y;
         }
      }
   } else {
   // optimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = -vsrc->x;
            vdst->y = vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = -vsrc->x;
            vdst->y = -vsrc->y;
         }
      }
   }
}

void ae_v_cmoved(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = alpha * vsrc->x;
            vdst->y = -alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = alpha * vsrc->x;
            vdst->y = alpha * vsrc->y;
         }
      }
   } else {
   // optimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = alpha * vsrc->x;
            vdst->y = -alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = alpha * vsrc->x;
            vdst->y = alpha * vsrc->y;
         }
      }
   }
}

void ae_v_cmovec(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, ae_complex alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = ax * vsrc->x + ay * vsrc->y;
            vdst->y = -ax * vsrc->y + ay * vsrc->x;
         }
      } else {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = ax * vsrc->x - ay * vsrc->y;
            vdst->y = ax * vsrc->y + ay * vsrc->x;
         }
      }
   } else {
   // highly optimized case
      if (bconj) {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = ax * vsrc->x + ay * vsrc->y;
            vdst->y = -ax * vsrc->y + ay * vsrc->x;
         }
      } else {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = ax * vsrc->x - ay * vsrc->y;
            vdst->y = ax * vsrc->y + ay * vsrc->x;
         }
      }
   }
}

void ae_v_cadd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += vsrc->x;
            vdst->y -= vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += vsrc->x;
            vdst->y += vsrc->y;
         }
      }
   } else {
   // optimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += vsrc->x;
            vdst->y -= vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += vsrc->x;
            vdst->y += vsrc->y;
         }
      }
   }
}

void ae_v_caddd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += alpha * vsrc->x;
            vdst->y -= alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += alpha * vsrc->x;
            vdst->y += alpha * vsrc->y;
         }
      }
   } else {
   // optimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += alpha * vsrc->x;
            vdst->y -= alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += alpha * vsrc->x;
            vdst->y += alpha * vsrc->y;
         }
      }
   }
}

void ae_v_caddc(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, ae_complex alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      double ax = alpha.x, ay = alpha.y;
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += ax * vsrc->x + ay * vsrc->y;
            vdst->y -= ax * vsrc->y - ay * vsrc->x;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += ax * vsrc->x - ay * vsrc->y;
            vdst->y += ax * vsrc->y + ay * vsrc->x;
         }
      }
   } else {
   // highly optimized case
      double ax = alpha.x, ay = alpha.y;
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += ax * vsrc->x + ay * vsrc->y;
            vdst->y -= ax * vsrc->y - ay * vsrc->x;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += ax * vsrc->x - ay * vsrc->y;
            vdst->y += ax * vsrc->y + ay * vsrc->x;
         }
      }
   }
}

void ae_v_csub(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x -= vsrc->x;
            vdst->y += vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x -= vsrc->x;
            vdst->y -= vsrc->y;
         }
      }
   } else {
   // highly optimized case
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x -= vsrc->x;
            vdst->y += vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x -= vsrc->x;
            vdst->y -= vsrc->y;
         }
      }
   }
}

void ae_v_csubd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha) {
   ae_v_caddd(vdst, stride_dst, vsrc, stride_src, conj_src, n, -alpha);
}

void ae_v_csubc(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, ae_complex alpha) {
   ae_v_caddc(vdst, stride_dst, vsrc, stride_src, conj_src, n, ae_complex_from_d(-alpha.x, -alpha.y));
}

void ae_v_cmuld(ae_complex *vdst, ae_int_t stride_dst, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst) {
         vdst->x *= alpha;
         vdst->y *= alpha;
      }
   } else {
   // optimized case
      for (i = 0; i < n; i++, vdst++) {
         vdst->x *= alpha;
         vdst->y *= alpha;
      }
   }
}

void ae_v_cmulc(ae_complex *vdst, ae_int_t stride_dst, ae_int_t n, ae_complex alpha) {
   ae_int_t i;
   if (stride_dst != 1) {
   // general unoptimized case
      double ax = alpha.x, ay = alpha.y;
      for (i = 0; i < n; i++, vdst += stride_dst) {
         double dstx = vdst->x, dsty = vdst->y;
         vdst->x = ax * dstx - ay * dsty;
         vdst->y = ax * dsty + ay * dstx;
      }
   } else {
   // highly optimized case
      double ax = alpha.x, ay = alpha.y;
      for (i = 0; i < n; i++, vdst++) {
         double dstx = vdst->x, dsty = vdst->y;
         vdst->x = ax * dstx - ay * dsty;
         vdst->y = ax * dsty + ay * dstx;
      }
   }
}

// Real BLAS operations
double ae_v_dotproduct(const double *v0, ae_int_t stride0, const double *v1, ae_int_t stride1, ae_int_t n) {
   double result = 0.0;
   ae_int_t i;
   if (stride0 != 1 || stride1 != 1) {
   // slow general code
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1)
         result += (*v0) * (*v1);
   } else {
   // optimized code for stride=1
      ae_int_t n4 = n / 4;
      ae_int_t nleft = n % 4;
      for (i = 0; i < n4; i++, v0 += 4, v1 += 4)
         result += v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] + v0[3] * v1[3];
      for (i = 0; i < nleft; i++, v0++, v1++)
         result += v0[0] * v1[0];
   }
   return result;
}

void ae_v_move(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst = *vsrc;
   } else {
   // optimized case
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] = vsrc[0];
         vdst[1] = vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] = vsrc[0];
   }
}

void ae_v_moveneg(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst = -*vsrc;
   } else {
   // optimized case
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] = -vsrc[0];
         vdst[1] = -vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] = -vsrc[0];
   }
}

void ae_v_moved(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst = alpha * (*vsrc);
   } else {
   // optimized case
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] = alpha * vsrc[0];
         vdst[1] = alpha * vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] = alpha * vsrc[0];
   }
}

void ae_v_add(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst += *vsrc;
   } else {
   // optimized case
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] += vsrc[0];
         vdst[1] += vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] += vsrc[0];
   }
}

void ae_v_addd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst += alpha * (*vsrc);
   } else {
   // optimized case
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] += alpha * vsrc[0];
         vdst[1] += alpha * vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] += alpha * vsrc[0];
   }
}

void ae_v_sub(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst -= *vsrc;
   } else {
   // highly optimized case
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] -= vsrc[0];
         vdst[1] -= vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] -= vsrc[0];
   }
}

void ae_v_subd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha) {
   ae_v_addd(vdst, stride_dst, vsrc, stride_src, n, -alpha);
}

void ae_v_muld(double *vdst, ae_int_t stride_dst, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1) {
   // general unoptimized case
      for (i = 0; i < n; i++, vdst += stride_dst)
         *vdst *= alpha;
   } else {
   // highly optimized case
      for (i = 0; i < n; i++, vdst++)
         *vdst *= alpha;
   }
}

// Linear Algebra
#define alglib_simd_alignment 16
#define alglib_r_block        32
#define alglib_half_r_block   16
#define alglib_twice_r_block  64
#define alglib_c_block        16
#define alglib_half_c_block    8
#define alglib_twice_c_block  32

// This subroutine calculates fast 32x32 real matrix-vector product:
//
//     y := beta*y + alpha*A*x
//
// using either generic C code or native optimizations (if available)
//
// IMPORTANT:
// * A must be stored in row-major order,
//   stride is alglib_r_block,
//   aligned on alglib_simd_alignment boundary
// * X must be aligned on alglib_simd_alignment boundary
// * Y may be non-aligned
void _ialglib_mv_32(const double *a, const double *x, double *y, ae_int_t stride, double alpha, double beta) {
   ae_int_t i, k;
   const double *pa0, *pa1, *pb;
   pa0 = a;
   pa1 = a + alglib_r_block;
   pb = x;
   for (i = 0; i < 16; i++) {
      double v0 = 0.0, v1 = 0.0;
      for (k = 0; k < 4; k++) {
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
      y[0] = beta * y[0] + alpha * v0;
      y[stride] = beta * y[stride] + alpha * v1;
   // now we've processed rows I and I+1,
   // pa0 and pa1 are pointing to rows I+1 and I+2.
   // move to I+2 and I+3.
      pa0 += alglib_r_block;
      pa1 += alglib_r_block;
      pb = x;
      y += 2 * stride;
   }
}

// This function calculates MxN real matrix-vector product:
//
//     y := beta*y + alpha*A*x
//
// using generic C code. It calls _ialglib_mv_32 if both M=32 and N=32.
//
// If beta is zero, we do not use previous values of y (they are  overwritten
// by alpha*A*x without ever being read).  If alpha is zero, no matrix-vector
// product is calculated (only beta is updated); however, this update  is not
// efficient  and  this  function  should  NOT  be used for multiplication of
// vector and scalar.
//
// IMPORTANT:
// * 0 <= M <= alglib_r_block, 0 <= N <= alglib_r_block
// * A must be stored in row-major order with stride equal to alglib_r_block
void _ialglib_rmv(ae_int_t m, ae_int_t n, const double *a, const double *x, double *y, ae_int_t stride, double alpha, double beta) {
// Handle special cases:
// - alpha is zero or n is zero
// - m is zero
   if (m == 0)
      return;
   if (alpha == 0.0 || n == 0) {
      ae_int_t i;
      if (beta == 0.0) {
         for (i = 0; i < m; i++) {
            *y = 0.0;
            y += stride;
         }
      } else {
         for (i = 0; i < m; i++) {
            *y *= beta;
            y += stride;
         }
      }
      return;
   }
// Handle general case: nonzero alpha, n and m
//
   if (m == 32 && n == 32) {
   // 32x32, may be we have something better than general implementation
      _ialglib_mv_32(a, x, y, stride, alpha, beta);
   } else {
      ae_int_t i, k, m2, n8, n2, ntrail2;
      const double *pa0, *pa1, *pb;
   // First M/2 rows of A are processed in pairs.
   // optimized code is used.
      m2 = m / 2;
      n8 = n / 8;
      ntrail2 = (n - 8 * n8) / 2;
      for (i = 0; i < m2; i++) {
         double v0 = 0.0, v1 = 0.0;
      // 'a' points to the part of the matrix which
      // is not processed yet
         pb = x;
         pa0 = a;
         pa1 = a + alglib_r_block;
         a += alglib_twice_r_block;
      // 8 elements per iteration
         for (k = 0; k < n8; k++) {
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
      // 2 elements per iteration
         for (k = 0; k < ntrail2; k++) {
            v0 += pa0[0] * pb[0];
            v1 += pa1[0] * pb[0];
            v0 += pa0[1] * pb[1];
            v1 += pa1[1] * pb[1];
            pa0 += 2;
            pa1 += 2;
            pb += 2;
         }
      // last element, if needed
         if (n % 2 != 0) {
            v0 += pa0[0] * pb[0];
            v1 += pa1[0] * pb[0];
         }
      // final update
         if (beta != 0) {
            y[0] = beta * y[0] + alpha * v0;
            y[stride] = beta * y[stride] + alpha * v1;
         } else {
            y[0] = alpha * v0;
            y[stride] = alpha * v1;
         }
      // move to the next pair of elements
         y += 2 * stride;
      }
   // Last (odd) row is processed with less optimized code.
      if (m % 2 != 0) {
         double v0 = 0.0;
      // 'a' points to the part of the matrix which
      // is not processed yet
         pb = x;
         pa0 = a;
      // 2 elements per iteration
         n2 = n / 2;
         for (k = 0; k < n2; k++) {
            v0 += pa0[0] * pb[0] + pa0[1] * pb[1];
            pa0 += 2;
            pb += 2;
         }
      // last element, if needed
         if (n % 2 != 0)
            v0 += pa0[0] * pb[0];
      // final update
         if (beta != 0)
            y[0] = beta * y[0] + alpha * v0;
         else
            y[0] = alpha * v0;
      }
   }
}

// This function calculates MxN real matrix-vector product:
//
//     y := beta*y + alpha*A*x
//
// using generic C code. It calls _ialglib_mv_32 if both M=32 and N=32.
//
// If beta is zero, we do not use previous values of y (they are  overwritten
// by alpha*A*x without ever being read).  If alpha is zero, no matrix-vector
// product is calculated (only beta is updated); however, this update  is not
// efficient  and  this  function  should  NOT  be used for multiplication of
// vector and scalar.
//
// IMPORTANT:
// * 0 <= M <= alglib_r_block, 0 <= N <= alglib_r_block
// * A must be stored in row-major order with stride equal to alglib_r_block
// * y may be non-aligned
// * both A and x must have same offset with respect to 16-byte boundary:
//   either both are aligned, or both are aligned with offset 8. Function
//   will crash your system if you try to call it with misaligned or
//   incorrectly aligned data.
//
// This function supports SSE2; it can be used when:
// 1. AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// 2. CurCPU contains CPU_SSE2 (checked at run-time)
//
// If (1) is failed, this function will be undefined. If (2) is failed,  call
// to this function will probably crash your system.
//
// If  you  want  to  know  whether  it  is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
#if defined AE_HAS_SSE2_INTRINSICS
void _ialglib_rmv_sse2(ae_int_t m, ae_int_t n, const double *a, const double *x, double *y, ae_int_t stride, double alpha, double beta) {
   ae_int_t i, k, n2;
   ae_int_t mb3, mtail, nhead, nb8, nb2, ntail;
   const double *pa0, *pa1, *pa2, *pb;
   __m128d v0, v1, v2, va0, va1, va2, vx, vtmp;
// Handle special cases:
// - alpha is zero or n is zero
// - m is zero
   if (m == 0)
      return;
   if (alpha == 0.0 || n == 0) {
      if (beta == 0.0) {
         for (i = 0; i < m; i++) {
            *y = 0.0;
            y += stride;
         }
      } else {
         for (i = 0; i < m; i++) {
            *y *= beta;
            y += stride;
         }
      }
      return;
   }
// Handle general case: nonzero alpha, n and m
//
// We divide problem as follows...
//
// Rows M are divided into:
// - mb3 blocks, each 3xN
// - mtail blocks, each 1xN
//
// Within a row, elements are divided  into:
// - nhead 1x1 blocks (used to align the rest, either 0 or 1)
// - nb8 1x8 blocks, aligned to 16-byte boundary
// - nb2 1x2 blocks, aligned to 16-byte boundary
// - ntail 1x1 blocks, aligned too (altough we don't rely on it)
//
   n2 = n / 2;
   mb3 = m / 3;
   mtail = m % 3;
   nhead = ae_misalignment(a, alglib_simd_alignment) == 0 ? 0 : 1;
   nb8 = (n - nhead) / 8;
   nb2 = (n - nhead - 8 * nb8) / 2;
   ntail = n - nhead - 8 * nb8 - 2 * nb2;
   for (i = 0; i < mb3; i++) {
      double row0, row1, row2;
      row0 = 0;
      row1 = 0;
      row2 = 0;
      pb = x;
      pa0 = a;
      pa1 = a + alglib_r_block;
      pa2 = a + alglib_twice_r_block;
      a += 3 * alglib_r_block;
      if (nhead == 1) {
         vx = _mm_load_sd(pb);
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
      for (k = 0; k < nb8; k++) {
      // this code is a shuffle of simultaneous dot product.
      // see below for commented unshuffled original version.
         vx = _mm_load_pd(pb);
         va0 = _mm_load_pd(pa0);
         va1 = _mm_load_pd(pa1);
         va0 = _mm_mul_pd(va0, vx);
         va2 = _mm_load_pd(pa2);
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
#   if 0
      // this is unshuffled version of code above
         vx = _mm_load_pd(pb);
         va0 = _mm_load_pd(pa0);
         va1 = _mm_load_pd(pa1);
         va2 = _mm_load_pd(pa2);
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
      for (k = 0; k < nb2; k++) {
         vx = _mm_load_pd(pb);
         va0 = _mm_load_pd(pa0);
         va1 = _mm_load_pd(pa1);
         va2 = _mm_load_pd(pa2);
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
      for (k = 0; k < ntail; k++) {
         vx = _mm_load1_pd(pb);
         va0 = _mm_load1_pd(pa0);
         va1 = _mm_load1_pd(pa1);
         va2 = _mm_load1_pd(pa2);
         va0 = _mm_mul_sd(va0, vx);
         v0 = _mm_add_sd(v0, va0);
         va1 = _mm_mul_sd(va1, vx);
         v1 = _mm_add_sd(v1, va1);
         va2 = _mm_mul_sd(va2, vx);
         v2 = _mm_add_sd(v2, va2);
      }
      vtmp = _mm_add_pd(_mm_unpacklo_pd(v0, v1), _mm_unpackhi_pd(v0, v1));
      _mm_storel_pd(&row0, vtmp);
      _mm_storeh_pd(&row1, vtmp);
      v2 = _mm_add_sd(_mm_shuffle_pd(v2, v2, 1), v2);
      _mm_storel_pd(&row2, v2);
      if (beta != 0) {
         y[0] = beta * y[0] + alpha * row0;
         y[stride] = beta * y[stride] + alpha * row1;
         y[2 * stride] = beta * y[2 * stride] + alpha * row2;
      } else {
         y[0] = alpha * row0;
         y[stride] = alpha * row1;
         y[2 * stride] = alpha * row2;
      }
      y += 3 * stride;
   }
   for (i = 0; i < mtail; i++) {
      double row0;
      row0 = 0;
      pb = x;
      pa0 = a;
      a += alglib_r_block;
      for (k = 0; k < n2; k++) {
         row0 += pb[0] * pa0[0] + pb[1] * pa0[1];
         pa0 += 2;
         pb += 2;
      }
      if (n % 2)
         row0 += pb[0] * pa0[0];
      if (beta != 0)
         y[0] = beta * y[0] + alpha * row0;
      else
         y[0] = alpha * row0;
      y += stride;
   }
}
#endif

// This subroutine calculates fast MxN complex matrix-vector product:
//
//     y := beta*y + alpha*A*x
//
// using generic C code, where A, x, y, alpha and beta are complex.
//
// If beta is zero, we do not use previous values of y (they are  overwritten
// by alpha*A*x without ever being read). However, when  alpha  is  zero,  we
// still calculate A*x and  multiply  it  by  alpha  (this distinction can be
// important when A or x contain infinities/NANs).
//
// IMPORTANT:
// * 0 <= M <= alglib_c_block, 0 <= N <= alglib_c_block
// * A must be stored in row-major order, as sequence of double precision
//   pairs. Stride is alglib_c_block (it is measured in pairs of doubles, not
//   in doubles).
// * Y may be referenced by cy (pointer to ae_complex) or
//   dy (pointer to array of double precision pair) depending on what type of
//   output you wish. Pass pointer to Y as one of these parameters,
//   AND SET OTHER PARAMETER TO NULL.
// * both A and x must be aligned; y may be non-aligned.
void _ialglib_cmv(ae_int_t m, ae_int_t n, const double *a, const double *x, ae_complex *cy, double *dy, ae_int_t stride, ae_complex alpha, ae_complex beta) {
   ae_int_t i, j;
   const double *pa, *parow, *pb;
   parow = a;
   for (i = 0; i < m; i++) {
      double v0 = 0.0, v1 = 0.0;
      pa = parow;
      pb = x;
      for (j = 0; j < n; j++) {
         v0 += pa[0] * pb[0];
         v1 += pa[0] * pb[1];
         v0 -= pa[1] * pb[1];
         v1 += pa[1] * pb[0];
         pa += 2;
         pb += 2;
      }
      if (cy != NULL) {
         double tx = (beta.x * cy->x - beta.y * cy->y) + (alpha.x * v0 - alpha.y * v1);
         double ty = (beta.x * cy->y + beta.y * cy->x) + (alpha.x * v1 + alpha.y * v0);
         cy->x = tx;
         cy->y = ty;
         cy += stride;
      } else {
         double tx = (beta.x * dy[0] - beta.y * dy[1]) + (alpha.x * v0 - alpha.y * v1);
         double ty = (beta.x * dy[1] + beta.y * dy[0]) + (alpha.x * v1 + alpha.y * v0);
         dy[0] = tx;
         dy[1] = ty;
         dy += 2 * stride;
      }
      parow += 2 * alglib_c_block;
   }
}

// This subroutine calculates fast MxN complex matrix-vector product:
//
//     y := beta*y + alpha*A*x
//
// using generic C code, where A, x, y, alpha and beta are complex.
//
// If beta is zero, we do not use previous values of y (they are  overwritten
// by alpha*A*x without ever being read). However, when  alpha  is  zero,  we
// still calculate A*x and  multiply  it  by  alpha  (this distinction can be
// important when A or x contain infinities/NANs).
//
// IMPORTANT:
// * 0 <= M <= alglib_c_block, 0 <= N <= alglib_c_block
// * A must be stored in row-major order, as sequence of double precision
//   pairs. Stride is alglib_c_block (it is measured in pairs of doubles, not
//   in doubles).
// * Y may be referenced by cy (pointer to ae_complex) or
//   dy (pointer to array of double precision pair) depending on what type of
//   output you wish. Pass pointer to Y as one of these parameters,
//   AND SET OTHER PARAMETER TO NULL.
// * both A and x must be aligned; y may be non-aligned.
//
// This function supports SSE2; it can be used when:
// 1. AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// 2. CurCPU contains CPU_SSE2 (checked at run-time)
//
// If (1) is failed, this function will be undefined. If (2) is failed,  call
// to this function will probably crash your system.
//
// If  you  want  to  know  whether  it  is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
#if defined AE_HAS_SSE2_INTRINSICS
void _ialglib_cmv_sse2(ae_int_t m, ae_int_t n, const double *a, const double *x, ae_complex *cy, double *dy, ae_int_t stride, ae_complex alpha, ae_complex beta) {
   ae_int_t i, j, m2;
   const double *pa0, *pa1, *parow, *pb;
   __m128d vbeta, vbetax, vbetay;
   __m128d valpha, valphax, valphay;
   m2 = m / 2;
   parow = a;
   if (cy != NULL) {
      dy = (double *)cy;
      cy = NULL;
   }
   vbeta = _mm_loadh_pd(_mm_load_sd(&beta.x), &beta.y);
   vbetax = _mm_unpacklo_pd(vbeta, vbeta);
   vbetay = _mm_unpackhi_pd(vbeta, vbeta);
   valpha = _mm_loadh_pd(_mm_load_sd(&alpha.x), &alpha.y);
   valphax = _mm_unpacklo_pd(valpha, valpha);
   valphay = _mm_unpackhi_pd(valpha, valpha);
   for (i = 0; i < m2; i++) {
      __m128d vx, vy, vt0, vt1, vt2, vt3, vt4, vt5, vrx, vry, vtx, vty;
      pa0 = parow;
      pa1 = parow + 2 * alglib_c_block;
      pb = x;
      vx = _mm_setzero_pd();
      vy = _mm_setzero_pd();
      for (j = 0; j < n; j++) {
         vt0 = _mm_load1_pd(pb);
         vt1 = _mm_load1_pd(pb + 1);
         vt2 = _mm_load_pd(pa0);
         vt3 = _mm_load_pd(pa1);
         vt5 = _mm_unpacklo_pd(vt2, vt3);
         vt4 = _mm_unpackhi_pd(vt2, vt3);
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
      if (beta.x == 0.0 && beta.y == 0.0) {
         vrx = _mm_setzero_pd();
         vry = _mm_setzero_pd();
      } else {
         vtx = _mm_loadh_pd(_mm_load_sd(dy + 0), dy + 2 * stride + 0);
         vty = _mm_loadh_pd(_mm_load_sd(dy + 1), dy + 2 * stride + 1);
         vrx = _mm_sub_pd(_mm_mul_pd(vbetax, vtx), _mm_mul_pd(vbetay, vty));
         vry = _mm_add_pd(_mm_mul_pd(vbetax, vty), _mm_mul_pd(vbetay, vtx));
      }
      vtx = _mm_sub_pd(_mm_mul_pd(valphax, vx), _mm_mul_pd(valphay, vy));
      vty = _mm_add_pd(_mm_mul_pd(valphax, vy), _mm_mul_pd(valphay, vx));
      vrx = _mm_add_pd(vrx, vtx);
      vry = _mm_add_pd(vry, vty);
      _mm_storel_pd(dy + 0, vrx);
      _mm_storeh_pd(dy + 2 * stride + 0, vrx);
      _mm_storel_pd(dy + 1, vry);
      _mm_storeh_pd(dy + 2 * stride + 1, vry);
      dy += 4 * stride;
      parow += 4 * alglib_c_block;
   }
   if (m % 2) {
      double v0 = 0.0, v1 = 0.0;
      double tx, ty;
      pa0 = parow;
      pb = x;
      for (j = 0; j < n; j++) {
         v0 += pa0[0] * pb[0];
         v1 += pa0[0] * pb[1];
         v0 -= pa0[1] * pb[1];
         v1 += pa0[1] * pb[0];
         pa0 += 2;
         pb += 2;
      }
      if (beta.x == 0.0 && beta.y == 0.0) {
         tx = 0.0;
         ty = 0.0;
      } else {
         tx = beta.x * dy[0] - beta.y * dy[1];
         ty = beta.x * dy[1] + beta.y * dy[0];
      }
      tx += alpha.x * v0 - alpha.y * v1;
      ty += alpha.x * v1 + alpha.y * v0;
      dy[0] = tx;
      dy[1] = ty;
      dy += 2 * stride;
      parow += 2 * alglib_c_block;
   }
}
#endif

// This subroutine sets vector to zero
void _ialglib_vzero(ae_int_t n, double *p, ae_int_t stride) {
   ae_int_t i;
   if (stride == 1) {
      for (i = 0; i < n; i++, p++)
         *p = 0.0;
   } else {
      for (i = 0; i < n; i++, p += stride)
         *p = 0.0;
   }
}

// This subroutine sets vector to zero
void _ialglib_vzero_complex(ae_int_t n, ae_complex *p, ae_int_t stride) {
   ae_int_t i;
   if (stride == 1) {
      for (i = 0; i < n; i++, p++) {
         p->x = 0.0;
         p->y = 0.0;
      }
   } else {
      for (i = 0; i < n; i++, p += stride) {
         p->x = 0.0;
         p->y = 0.0;
      }
   }
}

// This subroutine copies unaligned real vector
void _ialglib_vcopy(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb) {
   ae_int_t i, n2;
   if (stridea == 1 && strideb == 1) {
      n2 = n / 2;
      for (i = n2; i != 0; i--, a += 2, b += 2) {
         b[0] = a[0];
         b[1] = a[1];
      }
      if (n % 2 != 0)
         b[0] = a[0];
   } else {
      for (i = 0; i < n; i++, a += stridea, b += strideb)
         *b = *a;
   }
}

// This subroutine copies unaligned complex vector
// (passed as ae_complex*)
//
// 1. strideb is stride measured in complex numbers, not doubles
// 2. conj may be "N" (no conj.) or "C" (conj.)
void _ialglib_vcopy_complex(ae_int_t n, const ae_complex *a, ae_int_t stridea, double *b, ae_int_t strideb, const char *conj) {
   ae_int_t i;
// more general case
   if (conj[0] == 'N' || conj[0] == 'n') {
      for (i = 0; i < n; i++, a += stridea, b += 2 * strideb) {
         b[0] = a->x;
         b[1] = a->y;
      }
   } else {
      for (i = 0; i < n; i++, a += stridea, b += 2 * strideb) {
         b[0] = a->x;
         b[1] = -a->y;
      }
   }
}

// This subroutine copies unaligned complex vector (passed as double*)
//
// 1. strideb is stride measured in complex numbers, not doubles
// 2. conj may be "N" (no conj.) or "C" (conj.)
void _ialglib_vcopy_dcomplex(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb, const char *conj) {
   ae_int_t i;
// more general case
   if (conj[0] == 'N' || conj[0] == 'n') {
      for (i = 0; i < n; i++, a += 2 * stridea, b += 2 * strideb) {
         b[0] = a[0];
         b[1] = a[1];
      }
   } else {
      for (i = 0; i < n; i++, a += 2 * stridea, b += 2 * strideb) {
         b[0] = a[0];
         b[1] = -a[1];
      }
   }
}

// This subroutine copies matrix from  non-aligned non-contigous storage
// to aligned contigous storage
//
// A:
// * MxN
// * non-aligned
// * non-contigous
// * may be transformed during copying (as prescribed by op)
//
// B:
// * alglib_r_block*alglib_r_block (only MxN/NxM submatrix is used)
// * aligned
// * stride is alglib_r_block
//
// Transformation types:
// * 0 - no transform
// * 1 - transposition
void _ialglib_mcopyblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_int_t stride, double *b) {
   ae_int_t i, j, n2;
   const double *psrc;
   double *pdst;
   if (op == 0) {
      n2 = n / 2;
      for (i = 0, psrc = a; i < m; i++, a += stride, b += alglib_r_block, psrc = a) {
         for (j = 0, pdst = b; j < n2; j++, pdst += 2, psrc += 2) {
            pdst[0] = psrc[0];
            pdst[1] = psrc[1];
         }
         if (n % 2 != 0)
            pdst[0] = psrc[0];
      }
   } else {
      n2 = n / 2;
      for (i = 0, psrc = a; i < m; i++, a += stride, b += 1, psrc = a) {
         for (j = 0, pdst = b; j < n2; j++, pdst += alglib_twice_r_block, psrc += 2) {
            pdst[0] = psrc[0];
            pdst[alglib_r_block] = psrc[1];
         }
         if (n % 2 != 0)
            pdst[0] = psrc[0];
      }
   }
}

// This subroutine copies matrix from  non-aligned non-contigous storage
// to aligned contigous storage
//
// A:
// * MxN
// * non-aligned
// * non-contigous
// * may be transformed during copying (as prescribed by op)
//
// B:
// * alglib_r_block*alglib_r_block (only MxN/NxM submatrix is used)
// * aligned
// * stride is alglib_r_block
//
// Transformation types:
// * 0 - no transform
// * 1 - transposition
//
// This function supports SSE2; it can be used when:
// 1. AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// 2. CurCPU contains CPU_SSE2 (checked at run-time)
//
// If (1) is failed, this function will be undefined. If (2) is failed,  call
// to this function will probably crash your system.
//
// If  you  want  to  know  whether  it  is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
#if defined AE_HAS_SSE2_INTRINSICS
void _ialglib_mcopyblock_sse2(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_int_t stride, double *b) {
   ae_int_t i, j, mb2;
   const double *psrc0, *psrc1;
   double *pdst;
   if (op == 0) {
      ae_int_t nb8, ntail;
      nb8 = n / 8;
      ntail = n - 8 * nb8;
      for (i = 0, psrc0 = a; i < m; i++, a += stride, b += alglib_r_block, psrc0 = a) {
         pdst = b;
         for (j = 0; j < nb8; j++) {
            __m128d v0, v1;
            v0 = _mm_loadu_pd(psrc0);
            _mm_store_pd(pdst, v0);
            v1 = _mm_loadu_pd(psrc0 + 2);
            _mm_store_pd(pdst + 2, v1);
            v1 = _mm_loadu_pd(psrc0 + 4);
            _mm_store_pd(pdst + 4, v1);
            v1 = _mm_loadu_pd(psrc0 + 6);
            _mm_store_pd(pdst + 6, v1);
            pdst += 8;
            psrc0 += 8;
         }
         for (j = 0; j < ntail; j++)
            pdst[j] = psrc0[j];
      }
   } else {
      const double *arow0, *arow1;
      double *bcol0, *bcol1, *pdst0, *pdst1;
      ae_int_t nb4, ntail, n2;
      n2 = n / 2;
      mb2 = m / 2;
      nb4 = n / 4;
      ntail = n - 4 * nb4;
      arow0 = a;
      arow1 = a + stride;
      bcol0 = b;
      bcol1 = b + 1;
      for (i = 0; i < mb2; i++) {
         psrc0 = arow0;
         psrc1 = arow1;
         pdst0 = bcol0;
         pdst1 = bcol1;
         for (j = 0; j < nb4; j++) {
            __m128d v0, v1, v2, v3;
            v0 = _mm_loadu_pd(psrc0);
            v1 = _mm_loadu_pd(psrc1);
            v2 = _mm_loadu_pd(psrc0 + 2);
            v3 = _mm_loadu_pd(psrc1 + 2);
            _mm_store_pd(pdst0, _mm_unpacklo_pd(v0, v1));
            _mm_store_pd(pdst0 + alglib_r_block, _mm_unpackhi_pd(v0, v1));
            _mm_store_pd(pdst0 + 2 * alglib_r_block, _mm_unpacklo_pd(v2, v3));
            _mm_store_pd(pdst0 + 3 * alglib_r_block, _mm_unpackhi_pd(v2, v3));
            pdst0 += 4 * alglib_r_block;
            pdst1 += 4 * alglib_r_block;
            psrc0 += 4;
            psrc1 += 4;
         }
         for (j = 0; j < ntail; j++) {
            pdst0[0] = psrc0[0];
            pdst1[0] = psrc1[0];
            pdst0 += alglib_r_block;
            pdst1 += alglib_r_block;
            psrc0 += 1;
            psrc1 += 1;
         }
         arow0 += 2 * stride;
         arow1 += 2 * stride;
         bcol0 += 2;
         bcol1 += 2;
      }
      if (m % 2) {
         psrc0 = arow0;
         pdst0 = bcol0;
         for (j = 0; j < n2; j++) {
            pdst0[0] = psrc0[0];
            pdst0[alglib_r_block] = psrc0[1];
            pdst0 += alglib_twice_r_block;
            psrc0 += 2;
         }
         if (n % 2 != 0)
            pdst0[0] = psrc0[0];
      }
   }
}
#endif

// This subroutine copies matrix from  aligned contigous storage to non-
// aligned non-contigous storage
//
// A:
// * MxN
// * aligned
// * contigous
// * stride is alglib_r_block
// * may be transformed during copying (as prescribed by op)
//
// B:
// * alglib_r_block*alglib_r_block (only MxN/NxM submatrix is used)
// * non-aligned, non-contigous
//
// Transformation types:
// * 0 - no transform
// * 1 - transposition
void _ialglib_mcopyunblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, double *b, ae_int_t stride) {
   ae_int_t i, j, n2;
   const double *psrc;
   double *pdst;
   if (op == 0) {
      n2 = n / 2;
      for (i = 0, psrc = a; i < m; i++, a += alglib_r_block, b += stride, psrc = a) {
         for (j = 0, pdst = b; j < n2; j++, pdst += 2, psrc += 2) {
            pdst[0] = psrc[0];
            pdst[1] = psrc[1];
         }
         if (n % 2 != 0)
            pdst[0] = psrc[0];
      }
   } else {
      n2 = n / 2;
      for (i = 0, psrc = a; i < m; i++, a++, b += stride, psrc = a) {
         for (j = 0, pdst = b; j < n2; j++, pdst += 2, psrc += alglib_twice_r_block) {
            pdst[0] = psrc[0];
            pdst[1] = psrc[alglib_r_block];
         }
         if (n % 2 != 0)
            pdst[0] = psrc[0];
      }
   }
}

// This subroutine copies matrix from  non-aligned non-contigous storage
// to aligned contigous storage
//
// A:
// * MxN
// * non-aligned
// * non-contigous
// * may be transformed during copying (as prescribed by op)
// * pointer to ae_complex is passed
//
// B:
// * 2*alglib_c_block*alglib_c_block doubles (only MxN/NxM submatrix is used)
// * aligned
// * stride is alglib_c_block
// * pointer to double is passed
//
// Transformation types:
// * 0 - no transform
// * 1 - transposition
// * 2 - conjugate transposition
// * 3 - conjugate, but no  transposition
void _ialglib_mcopyblock_complex(ae_int_t m, ae_int_t n, const ae_complex *a, ae_int_t op, ae_int_t stride, double *b) {
   ae_int_t i, j;
   const ae_complex *psrc;
   double *pdst;
   if (op == 0) {
      for (i = 0, psrc = a; i < m; i++, a += stride, b += alglib_twice_c_block, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst += 2, psrc++) {
            pdst[0] = psrc->x;
            pdst[1] = psrc->y;
         }
   }
   if (op == 1) {
      for (i = 0, psrc = a; i < m; i++, a += stride, b += 2, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst += alglib_twice_c_block, psrc++) {
            pdst[0] = psrc->x;
            pdst[1] = psrc->y;
         }
   }
   if (op == 2) {
      for (i = 0, psrc = a; i < m; i++, a += stride, b += 2, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst += alglib_twice_c_block, psrc++) {
            pdst[0] = psrc->x;
            pdst[1] = -psrc->y;
         }
   }
   if (op == 3) {
      for (i = 0, psrc = a; i < m; i++, a += stride, b += alglib_twice_c_block, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst += 2, psrc++) {
            pdst[0] = psrc->x;
            pdst[1] = -psrc->y;
         }
   }
}

// This subroutine copies matrix from aligned contigous storage to
// non-aligned non-contigous storage
//
// A:
// * 2*alglib_c_block*alglib_c_block doubles (only MxN submatrix is used)
// * aligned
// * stride is alglib_c_block
// * pointer to double is passed
// * may be transformed during copying (as prescribed by op)
//
// B:
// * MxN
// * non-aligned
// * non-contigous
// * pointer to ae_complex is passed
//
// Transformation types:
// * 0 - no transform
// * 1 - transposition
// * 2 - conjugate transposition
// * 3 - conjugate, but no  transposition
void _ialglib_mcopyunblock_complex(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_complex *b, ae_int_t stride) {
   ae_int_t i, j;
   const double *psrc;
   ae_complex *pdst;
   if (op == 0) {
      for (i = 0, psrc = a; i < m; i++, a += alglib_twice_c_block, b += stride, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst++, psrc += 2) {
            pdst->x = psrc[0];
            pdst->y = psrc[1];
         }
   }
   if (op == 1) {
      for (i = 0, psrc = a; i < m; i++, a += 2, b += stride, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst++, psrc += alglib_twice_c_block) {
            pdst->x = psrc[0];
            pdst->y = psrc[1];
         }
   }
   if (op == 2) {
      for (i = 0, psrc = a; i < m; i++, a += 2, b += stride, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst++, psrc += alglib_twice_c_block) {
            pdst->x = psrc[0];
            pdst->y = -psrc[1];
         }
   }
   if (op == 3) {
      for (i = 0, psrc = a; i < m; i++, a += alglib_twice_c_block, b += stride, psrc = a)
         for (j = 0, pdst = b; j < n; j++, pdst++, psrc += 2) {
            pdst->x = psrc[0];
            pdst->y = -psrc[1];
         }
   }
}

// Real GEMM kernel
bool _ialglib_rmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, double *_a, ae_int_t _a_stride, ae_int_t optypea, double *_b, ae_int_t _b_stride, ae_int_t optypeb, double beta, double *_c, ae_int_t _c_stride) {
   int i;
   double *crow;
   double _abuf[alglib_r_block + alglib_simd_alignment];
   double _bbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_abuf, alglib_simd_alignment);
   double *const b = (double *)ae_align(_bbuf, alglib_simd_alignment);
   void (*rmv)(ae_int_t, ae_int_t, const double *, const double *, double *, ae_int_t, double, double) = &_ialglib_rmv;
   void (*mcopyblock)(ae_int_t, ae_int_t, const double *, ae_int_t, ae_int_t, double *) = &_ialglib_mcopyblock;
   if (m > alglib_r_block || n > alglib_r_block || k > alglib_r_block || m <= 0 || n <= 0 || k <= 0 || alpha == 0.0)
      return false;
// Check for SSE2 support
#ifdef AE_HAS_SSE2_INTRINSICS
   if (CurCPU & CPU_SSE2) {
      rmv = &_ialglib_rmv_sse2;
      mcopyblock = &_ialglib_mcopyblock_sse2;
   }
#endif
// copy b
   if (optypeb == 0)
      mcopyblock(k, n, _b, 1, _b_stride, b);
   else
      mcopyblock(n, k, _b, 0, _b_stride, b);
// multiply B by A (from the right, by rows)
// and store result in C
   crow = _c;
   if (optypea == 0) {
      const double *arow = _a;
      for (i = 0; i < m; i++) {
         _ialglib_vcopy(k, arow, 1, abuf, 1);
         if (beta == 0)
            _ialglib_vzero(n, crow, 1);
         rmv(n, k, b, abuf, crow, 1, alpha, beta);
         crow += _c_stride;
         arow += _a_stride;
      }
   } else {
      const double *acol = _a;
      for (i = 0; i < m; i++) {
         _ialglib_vcopy(k, acol, _a_stride, abuf, 1);
         if (beta == 0)
            _ialglib_vzero(n, crow, 1);
         rmv(n, k, b, abuf, crow, 1, alpha, beta);
         crow += _c_stride;
         acol++;
      }
   }
   return true;
}

// Complex GEMM kernel
bool _ialglib_cmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, ae_complex *_a, ae_int_t _a_stride, ae_int_t optypea, ae_complex *_b, ae_int_t _b_stride, ae_int_t optypeb, ae_complex beta, ae_complex *_c, ae_int_t _c_stride) {
   const ae_complex *arow;
   ae_complex *crow;
   ae_int_t i;
   double _loc_abuf[2 * alglib_c_block + alglib_simd_alignment];
   double _loc_b[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const b = (double *)ae_align(_loc_b, alglib_simd_alignment);
   ae_int_t brows;
   ae_int_t bcols;
   void (*cmv)(ae_int_t, ae_int_t, const double *, const double *, ae_complex *, double *, ae_int_t, ae_complex, ae_complex) = &_ialglib_cmv;
   if (m > alglib_c_block || n > alglib_c_block || k > alglib_c_block)
      return false;
// Check for SSE2 support
#ifdef AE_HAS_SSE2_INTRINSICS
   if (CurCPU & CPU_SSE2) {
      cmv = &_ialglib_cmv_sse2;
   }
#endif
// copy b
   brows = optypeb == 0 ? k : n;
   bcols = optypeb == 0 ? n : k;
   if (optypeb == 0)
      _ialglib_mcopyblock_complex(brows, bcols, _b, 1, _b_stride, b);
   if (optypeb == 1)
      _ialglib_mcopyblock_complex(brows, bcols, _b, 0, _b_stride, b);
   if (optypeb == 2)
      _ialglib_mcopyblock_complex(brows, bcols, _b, 3, _b_stride, b);
// multiply B by A (from the right, by rows)
// and store result in C
   arow = _a;
   crow = _c;
   for (i = 0; i < m; i++) {
      if (optypea == 0) {
         _ialglib_vcopy_complex(k, arow, 1, abuf, 1, "No conj");
         arow += _a_stride;
      } else if (optypea == 1) {
         _ialglib_vcopy_complex(k, arow, _a_stride, abuf, 1, "No conj");
         arow++;
      } else {
         _ialglib_vcopy_complex(k, arow, _a_stride, abuf, 1, "Conj");
         arow++;
      }
      if (beta.x == 0.0 && beta.y == 0.0)
         _ialglib_vzero_complex(n, crow, 1);
      cmv(n, k, b, abuf, crow, NULL, 1, alpha, beta);
      crow += _c_stride;
   }
   return true;
}

// complex TRSM kernel
bool _ialglib_cmatrixrighttrsm(ae_int_t m, ae_int_t n, ae_complex *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t optype, ae_complex *_x, ae_int_t _x_stride) {
// local buffers
   double *pdiag;
   ae_int_t i;
   double _loc_abuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double _loc_xbuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double _loc_tmpbuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   bool uppera;
   void (*cmv)(ae_int_t, ae_int_t, const double *, const double *, ae_complex *, double *, ae_int_t, ae_complex, ae_complex) = &_ialglib_cmv;
   if (m > alglib_c_block || n > alglib_c_block)
      return false;
// Check for SSE2 support
#ifdef AE_HAS_SSE2_INTRINSICS
   if (CurCPU & CPU_SSE2) {
      cmv = &_ialglib_cmv_sse2;
   }
#endif
// Prepare
   _ialglib_mcopyblock_complex(n, n, _a, optype, _a_stride, abuf);
   _ialglib_mcopyblock_complex(m, n, _x, 0, _x_stride, xbuf);
   if (isunit)
      for (i = 0, pdiag = abuf; i < n; i++, pdiag += 2 * (alglib_c_block + 1)) {
         pdiag[0] = 1.0;
         pdiag[1] = 0.0;
      }
   if (optype == 0)
      uppera = isupper;
   else
      uppera = !isupper;
// Solve Y*A^-1=X where A is upper or lower triangular
   if (uppera) {
      for (i = 0, pdiag = abuf; i < n; i++, pdiag += 2 * (alglib_c_block + 1)) {
         ae_complex tmp_c;
         ae_complex beta;
         ae_complex alpha;
         tmp_c = ae_complex_from_d(pdiag[0], pdiag[1]);
         beta = ae_c_d_div(1.0, tmp_c);
         alpha = ae_complex_from_d(-beta.x, -beta.y);
         _ialglib_vcopy_dcomplex(i, abuf + 2 * i, alglib_c_block, tmpbuf, 1, "No conj");
         cmv(m, i, xbuf, tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, alpha, beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 0, _x, _x_stride);
   } else {
      for (i = n - 1, pdiag = abuf + 2 * ((n - 1) * alglib_c_block + (n - 1)); i >= 0; i--, pdiag -= 2 * (alglib_c_block + 1)) {
         ae_complex tmp_c;
         ae_complex beta;
         ae_complex alpha;
         tmp_c = ae_complex_from_d(pdiag[0], pdiag[1]);
         beta = ae_c_d_div(1.0, tmp_c);
         alpha = ae_complex_from_d(-beta.x, -beta.y);
         _ialglib_vcopy_dcomplex(n - 1 - i, pdiag + 2 * alglib_c_block, alglib_c_block, tmpbuf, 1, "No conj");
         cmv(m, n - 1 - i, xbuf + 2 * (i + 1), tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, alpha, beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 0, _x, _x_stride);
   }
   return true;
}

// real TRSM kernel
bool _ialglib_rmatrixrighttrsm(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t optype, double *_x, ae_int_t _x_stride) {
// local buffers
   double *pdiag;
   ae_int_t i;
   double _loc_abuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double _loc_xbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double _loc_tmpbuf[alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   bool uppera;
   void (*rmv)(ae_int_t, ae_int_t, const double *, const double *, double *, ae_int_t, double, double) = &_ialglib_rmv;
   void (*mcopyblock)(ae_int_t, ae_int_t, const double *, ae_int_t, ae_int_t, double *) = &_ialglib_mcopyblock;
   if (m > alglib_r_block || n > alglib_r_block)
      return false;
// Check for SSE2 support
#ifdef AE_HAS_SSE2_INTRINSICS
   if (CurCPU & CPU_SSE2) {
      rmv = &_ialglib_rmv_sse2;
      mcopyblock = &_ialglib_mcopyblock_sse2;
   }
#endif
// Prepare
   mcopyblock(n, n, _a, optype, _a_stride, abuf);
   mcopyblock(m, n, _x, 0, _x_stride, xbuf);
   if (isunit)
      for (i = 0, pdiag = abuf; i < n; i++, pdiag += alglib_r_block + 1)
         *pdiag = 1.0;
   if (optype == 0)
      uppera = isupper;
   else
      uppera = !isupper;
// Solve Y*A^-1=X where A is upper or lower triangular
   if (uppera) {
      for (i = 0, pdiag = abuf; i < n; i++, pdiag += alglib_r_block + 1) {
         double beta = 1.0 / (*pdiag);
         double alpha = -beta;
         _ialglib_vcopy(i, abuf + i, alglib_r_block, tmpbuf, 1);
         rmv(m, i, xbuf, tmpbuf, xbuf + i, alglib_r_block, alpha, beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 0, _x, _x_stride);
   } else {
      for (i = n - 1, pdiag = abuf + (n - 1) * alglib_r_block + (n - 1); i >= 0; i--, pdiag -= alglib_r_block + 1) {
         double beta = 1.0 / (*pdiag);
         double alpha = -beta;
         _ialglib_vcopy(n - 1 - i, pdiag + alglib_r_block, alglib_r_block, tmpbuf + i + 1, 1);
         rmv(m, n - 1 - i, xbuf + i + 1, tmpbuf + i + 1, xbuf + i, alglib_r_block, alpha, beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 0, _x, _x_stride);
   }
   return true;
}

// complex TRSM kernel
bool _ialglib_cmatrixlefttrsm(ae_int_t m, ae_int_t n, ae_complex *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t optype, ae_complex *_x, ae_int_t _x_stride) {
// local buffers
   double *pdiag, *arow;
   ae_int_t i;
   double _loc_abuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double _loc_xbuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double _loc_tmpbuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   bool uppera;
   void (*cmv)(ae_int_t, ae_int_t, const double *, const double *, ae_complex *, double *, ae_int_t, ae_complex, ae_complex) = &_ialglib_cmv;
   if (m > alglib_c_block || n > alglib_c_block)
      return false;
// Check for SSE2 support
#ifdef AE_HAS_SSE2_INTRINSICS
   if (CurCPU & CPU_SSE2) {
      cmv = &_ialglib_cmv_sse2;
   }
#endif
// Prepare
// Transpose X (so we may use mv, which calculates A*x, but not x*A)
   _ialglib_mcopyblock_complex(m, m, _a, optype, _a_stride, abuf);
   _ialglib_mcopyblock_complex(m, n, _x, 1, _x_stride, xbuf);
   if (isunit)
      for (i = 0, pdiag = abuf; i < m; i++, pdiag += 2 * (alglib_c_block + 1)) {
         pdiag[0] = 1.0;
         pdiag[1] = 0.0;
      }
   if (optype == 0)
      uppera = isupper;
   else
      uppera = !isupper;
// Solve A^-1*Y^T=X^T where A is upper or lower triangular
   if (uppera) {
      for (i = m - 1, pdiag = abuf + 2 * ((m - 1) * alglib_c_block + (m - 1)); i >= 0; i--, pdiag -= 2 * (alglib_c_block + 1)) {
         ae_complex tmp_c;
         ae_complex beta;
         ae_complex alpha;
         tmp_c = ae_complex_from_d(pdiag[0], pdiag[1]);
         beta = ae_c_d_div(1.0, tmp_c);
         alpha = ae_complex_from_d(-beta.x, -beta.y);
         _ialglib_vcopy_dcomplex(m - 1 - i, pdiag + 2, 1, tmpbuf, 1, "No conj");
         cmv(n, m - 1 - i, xbuf + 2 * (i + 1), tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, alpha, beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 1, _x, _x_stride);
   } else {
      for (i = 0, pdiag = abuf, arow = abuf; i < m; i++, pdiag += 2 * (alglib_c_block + 1), arow += 2 * alglib_c_block) {
         ae_complex tmp_c;
         ae_complex beta;
         ae_complex alpha;
         tmp_c = ae_complex_from_d(pdiag[0], pdiag[1]);
         beta = ae_c_d_div(1.0, tmp_c);
         alpha = ae_complex_from_d(-beta.x, -beta.y);
         _ialglib_vcopy_dcomplex(i, arow, 1, tmpbuf, 1, "No conj");
         cmv(n, i, xbuf, tmpbuf, NULL, xbuf + 2 * i, alglib_c_block, alpha, beta);
      }
      _ialglib_mcopyunblock_complex(m, n, xbuf, 1, _x, _x_stride);
   }
   return true;
}

// real TRSM kernel
bool _ialglib_rmatrixlefttrsm(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, bool isupper, bool isunit, ae_int_t optype, double *_x, ae_int_t _x_stride) {
// local buffers
   double *pdiag, *arow;
   ae_int_t i;
   double _loc_abuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double _loc_xbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double _loc_tmpbuf[alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const xbuf = (double *)ae_align(_loc_xbuf, alglib_simd_alignment);
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   bool uppera;
   void (*rmv)(ae_int_t, ae_int_t, const double *, const double *, double *, ae_int_t, double, double) = &_ialglib_rmv;
   void (*mcopyblock)(ae_int_t, ae_int_t, const double *, ae_int_t, ae_int_t, double *) = &_ialglib_mcopyblock;
   if (m > alglib_r_block || n > alglib_r_block)
      return false;
// Check for SSE2 support
#ifdef AE_HAS_SSE2_INTRINSICS
   if (CurCPU & CPU_SSE2) {
      rmv = &_ialglib_rmv_sse2;
      mcopyblock = &_ialglib_mcopyblock_sse2;
   }
#endif
// Prepare
// Transpose X (so we may use mv, which calculates A*x, but not x*A)
   mcopyblock(m, m, _a, optype, _a_stride, abuf);
   mcopyblock(m, n, _x, 1, _x_stride, xbuf);
   if (isunit)
      for (i = 0, pdiag = abuf; i < m; i++, pdiag += alglib_r_block + 1)
         *pdiag = 1.0;
   if (optype == 0)
      uppera = isupper;
   else
      uppera = !isupper;
// Solve A^-1*Y^T=X^T where A is upper or lower triangular
   if (uppera) {
      for (i = m - 1, pdiag = abuf + (m - 1) * alglib_r_block + (m - 1); i >= 0; i--, pdiag -= alglib_r_block + 1) {
         double beta = 1.0 / (*pdiag);
         double alpha = -beta;
         _ialglib_vcopy(m - 1 - i, pdiag + 1, 1, tmpbuf + i + 1, 1);
         rmv(n, m - 1 - i, xbuf + i + 1, tmpbuf + i + 1, xbuf + i, alglib_r_block, alpha, beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 1, _x, _x_stride);
   } else {
      for (i = 0, pdiag = abuf, arow = abuf; i < m; i++, pdiag += alglib_r_block + 1, arow += alglib_r_block) {
         double beta = 1.0 / (*pdiag);
         double alpha = -beta;
         _ialglib_vcopy(i, arow, 1, tmpbuf, 1);
         rmv(n, i, xbuf, tmpbuf, xbuf + i, alglib_r_block, alpha, beta);
      }
      _ialglib_mcopyunblock(m, n, xbuf, 1, _x, _x_stride);
   }
   return true;
}

// complex SYRK kernel
bool _ialglib_cmatrixherk(ae_int_t n, ae_int_t k, double alpha, ae_complex *_a, ae_int_t _a_stride, ae_int_t optypea, double beta, ae_complex *_c, ae_int_t _c_stride, bool isupper) {
// local buffers
   double *arow, *crow;
   ae_complex c_alpha, c_beta;
   ae_int_t i;
   double _loc_abuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double _loc_cbuf[2 * alglib_c_block * alglib_c_block + alglib_simd_alignment];
   double _loc_tmpbuf[2 * alglib_c_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const cbuf = (double *)ae_align(_loc_cbuf, alglib_simd_alignment);
   double *const tmpbuf = (double *)ae_align(_loc_tmpbuf, alglib_simd_alignment);
   if (n > alglib_c_block || k > alglib_c_block)
      return false;
   if (n == 0)
      return true;
// copy A and C, task is transformed to "A*A^H"-form.
// if beta == 0, then C is filled by zeros (and not referenced)
//
// alpha == 0 or k == 0 are correctly processed (A is not referenced)
   c_alpha = ae_complex_from_d(alpha);
   c_beta = ae_complex_from_d(beta);
   if (alpha == 0)
      k = 0;
   if (k > 0) {
      if (optypea == 0)
         _ialglib_mcopyblock_complex(n, k, _a, 3, _a_stride, abuf);
      else
         _ialglib_mcopyblock_complex(k, n, _a, 1, _a_stride, abuf);
   }
   _ialglib_mcopyblock_complex(n, n, _c, 0, _c_stride, cbuf);
   if (beta == 0) {
      for (i = 0, crow = cbuf; i < n; i++, crow += 2 * alglib_c_block)
         if (isupper)
            _ialglib_vzero(2 * (n - i), crow + 2 * i, 1);
         else
            _ialglib_vzero(2 * (i + 1), crow, 1);
   }
// update C
   if (isupper) {
      for (i = 0, arow = abuf, crow = cbuf; i < n; i++, arow += 2 * alglib_c_block, crow += 2 * alglib_c_block) {
         _ialglib_vcopy_dcomplex(k, arow, 1, tmpbuf, 1, "Conj");
         _ialglib_cmv(n - i, k, arow, tmpbuf, NULL, crow + 2 * i, 1, c_alpha, c_beta);
      }
   } else {
      for (i = 0, arow = abuf, crow = cbuf; i < n; i++, arow += 2 * alglib_c_block, crow += 2 * alglib_c_block) {
         _ialglib_vcopy_dcomplex(k, arow, 1, tmpbuf, 1, "Conj");
         _ialglib_cmv(i + 1, k, abuf, tmpbuf, NULL, crow, 1, c_alpha, c_beta);
      }
   }
// copy back
   _ialglib_mcopyunblock_complex(n, n, cbuf, 0, _c, _c_stride);
   return true;
}

// real SYRK kernel
bool _ialglib_rmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, double *_a, ae_int_t _a_stride, ae_int_t optypea, double beta, double *_c, ae_int_t _c_stride, bool isupper) {
// local buffers
   double *arow, *crow;
   ae_int_t i;
   double _loc_abuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double _loc_cbuf[alglib_r_block * alglib_r_block + alglib_simd_alignment];
   double *const abuf = (double *)ae_align(_loc_abuf, alglib_simd_alignment);
   double *const cbuf = (double *)ae_align(_loc_cbuf, alglib_simd_alignment);
   if (n > alglib_r_block || k > alglib_r_block)
      return false;
   if (n == 0)
      return true;
// copy A and C, task is transformed to "A*A^T"-form.
// if beta == 0, then C is filled by zeros (and not referenced)
//
// alpha == 0 or k == 0 are correctly processed (A is not referenced)
   if (alpha == 0)
      k = 0;
   if (k > 0) {
      if (optypea == 0)
         _ialglib_mcopyblock(n, k, _a, 0, _a_stride, abuf);
      else
         _ialglib_mcopyblock(k, n, _a, 1, _a_stride, abuf);
   }
   _ialglib_mcopyblock(n, n, _c, 0, _c_stride, cbuf);
   if (beta == 0) {
      for (i = 0, crow = cbuf; i < n; i++, crow += alglib_r_block)
         if (isupper)
            _ialglib_vzero(n - i, crow + i, 1);
         else
            _ialglib_vzero(i + 1, crow, 1);
   }
// update C
   if (isupper) {
      for (i = 0, arow = abuf, crow = cbuf; i < n; i++, arow += alglib_r_block, crow += alglib_r_block) {
         _ialglib_rmv(n - i, k, arow, arow, crow + i, 1, alpha, beta);
      }
   } else {
      for (i = 0, arow = abuf, crow = cbuf; i < n; i++, arow += alglib_r_block, crow += alglib_r_block) {
         _ialglib_rmv(i + 1, k, abuf, arow, crow, 1, alpha, beta);
      }
   }
// copy back
   _ialglib_mcopyunblock(n, n, cbuf, 0, _c, _c_stride);
   return true;
}

// complex rank-1 kernel
bool _ialglib_cmatrixrank1(ae_int_t m, ae_int_t n, ae_complex *_a, ae_int_t _a_stride, ae_complex *_u, ae_complex *_v) {
// Locals
   ae_complex *arow, *pu, *pv, *vtmp, *dst;
   ae_int_t n2 = n / 2;
   ae_int_t i, j;
// Quick exit
   if (m <= 0 || n <= 0)
      return false;
// update pairs of rows
   arow = _a;
   pu = _u;
   vtmp = _v;
   for (i = 0; i < m; i++, arow += _a_stride, pu++) {
   // update by two
      for (j = 0, pv = vtmp, dst = arow; j < n2; j++, dst += 2, pv += 2) {
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
   // final update
      if (n % 2 != 0) {
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

// real rank-1 kernel
// deprecated version
bool _ialglib_rmatrixrank1(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, double *_u, double *_v) {
// Locals
   double *arow0, *arow1, *pu, *pv, *vtmp, *dst0, *dst1;
   ae_int_t m2 = m / 2;
   ae_int_t n2 = n / 2;
   ae_int_t stride = _a_stride;
   ae_int_t stride2 = 2 * _a_stride;
   ae_int_t i, j;
// Quick exit
   if (m <= 0 || n <= 0)
      return false;
// update pairs of rows
   arow0 = _a;
   arow1 = arow0 + stride;
   pu = _u;
   vtmp = _v;
   for (i = 0; i < m2; i++, arow0 += stride2, arow1 += stride2, pu += 2) {
   // update by two
      for (j = 0, pv = vtmp, dst0 = arow0, dst1 = arow1; j < n2; j++, dst0 += 2, dst1 += 2, pv += 2) {
         dst0[0] += pu[0] * pv[0];
         dst0[1] += pu[0] * pv[1];
         dst1[0] += pu[1] * pv[0];
         dst1[1] += pu[1] * pv[1];
      }
   // final update
      if (n % 2 != 0) {
         dst0[0] += pu[0] * pv[0];
         dst1[0] += pu[1] * pv[0];
      }
   }
// update last row
   if (m % 2 != 0) {
   // update by two
      for (j = 0, pv = vtmp, dst0 = arow0; j < n2; j++, dst0 += 2, pv += 2) {
         dst0[0] += pu[0] * pv[0];
         dst0[1] += pu[0] * pv[1];
      }
   // final update
      if (n % 2 != 0)
         dst0[0] += pu[0] * pv[0];
   }
   return true;
}

// real rank-1 kernel
// deprecated version
bool _ialglib_rmatrixger(ae_int_t m, ae_int_t n, double *_a, ae_int_t _a_stride, double alpha, double *_u, double *_v) {
// Locals
   double *arow0, *arow1, *pu, *pv, *vtmp, *dst0, *dst1;
   ae_int_t m2 = m / 2;
   ae_int_t n2 = n / 2;
   ae_int_t stride = _a_stride;
   ae_int_t stride2 = 2 * _a_stride;
   ae_int_t i, j;
// Quick exit
   if (m <= 0 || n <= 0 || alpha == 0.0)
      return false;
// update pairs of rows
   arow0 = _a;
   arow1 = arow0 + stride;
   pu = _u;
   vtmp = _v;
   for (i = 0; i < m2; i++, arow0 += stride2, arow1 += stride2, pu += 2) {
      double au0 = alpha * pu[0];
      double au1 = alpha * pu[1];
   // update by two
      for (j = 0, pv = vtmp, dst0 = arow0, dst1 = arow1; j < n2; j++, dst0 += 2, dst1 += 2, pv += 2) {
         dst0[0] += au0 * pv[0];
         dst0[1] += au0 * pv[1];
         dst1[0] += au1 * pv[0];
         dst1[1] += au1 * pv[1];
      }
   // final update
      if (n % 2 != 0) {
         dst0[0] += au0 * pv[0];
         dst1[0] += au1 * pv[0];
      }
   }
// update last row
   if (m % 2 != 0) {
      double au0 = alpha * pu[0];
   // update by two
      for (j = 0, pv = vtmp, dst0 = arow0; j < n2; j++, dst0 += 2, pv += 2) {
         dst0[0] += au0 * pv[0];
         dst0[1] += au0 * pv[1];
      }
   // final update
      if (n % 2 != 0)
         dst0[0] += au0 * pv[0];
   }
   return true;
}

// Interface functions for efficient kernels
bool _ialglib_i_rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, ae_matrix *_a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *_b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, ae_matrix *_c, ae_int_t ic, ae_int_t jc) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (alpha == 0.0 || k == 0 || n == 0 || m == 0)
      return false;
// handle with optimized ALGLIB kernel
   return _ialglib_rmatrixgemm(m, n, k, alpha, _a->ptr.pp_double[ia] + ja, _a->stride, optypea, _b->ptr.pp_double[ib] + jb, _b->stride, optypeb, beta, _c->ptr.pp_double[ic] + jc, _c->stride);
}

bool _ialglib_i_cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, ae_matrix *_a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *_b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, ae_matrix *_c, ae_int_t ic, ae_int_t jc) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (alpha.x == 0.0 && alpha.y == 0 || k == 0 || n == 0 || m == 0)
      return false;
// handle with optimized ALGLIB kernel
   return _ialglib_cmatrixgemm(m, n, k, alpha, _a->ptr.pp_complex[ia] + ja, _a->stride, optypea, _b->ptr.pp_complex[ib] + jb, _b->stride, optypeb, beta, _c->ptr.pp_complex[ic] + jc, _c->stride);
}

bool _ialglib_i_cmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (m == 0 || n == 0)
      return false;
// handle with optimized ALGLIB kernel
   return _ialglib_cmatrixrighttrsm(m, n, &a->ptr.pp_complex[i1][j1], a->stride, isupper, isunit, optype, &x->ptr.pp_complex[i2][j2], x->stride);
}

bool _ialglib_i_rmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (m == 0 || n == 0)
      return false;
// handle with optimized ALGLIB kernel
   return _ialglib_rmatrixrighttrsm(m, n, &a->ptr.pp_double[i1][j1], a->stride, isupper, isunit, optype, &x->ptr.pp_double[i2][j2], x->stride);
}

bool _ialglib_i_cmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (m == 0 || n == 0)
      return false;
// handle with optimized ALGLIB kernel
   return _ialglib_cmatrixlefttrsm(m, n, &a->ptr.pp_complex[i1][j1], a->stride, isupper, isunit, optype, &x->ptr.pp_complex[i2][j2], x->stride);
}

bool _ialglib_i_rmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (m == 0 || n == 0)
      return false;
// handle with optimized ALGLIB kernel
   return _ialglib_rmatrixlefttrsm(m, n, &a->ptr.pp_double[i1][j1], a->stride, isupper, isunit, optype, &x->ptr.pp_double[i2][j2], x->stride);
}

bool _ialglib_i_cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (alpha == 0.0 || k == 0 || n == 0)
      return false;
// ALGLIB kernel
   return _ialglib_cmatrixherk(n, k, alpha, &a->ptr.pp_complex[ia][ja], a->stride, optypea, beta, &c->ptr.pp_complex[ic][jc], c->stride, isupper);
}

bool _ialglib_i_rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper) {
// handle degenerate cases like zero matrices by ALGLIB - greatly simplifies passing data to ALGLIB kernel
   if (alpha == 0.0 || k == 0 || n == 0)
      return false;
// ALGLIB kernel
   return _ialglib_rmatrixsyrk(n, k, alpha, &a->ptr.pp_double[ia][ja], a->stride, optypea, beta, &c->ptr.pp_double[ic][jc], c->stride, isupper);
}

bool _ialglib_i_cmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs) {
   return _ialglib_cmatrixrank1(m, n, &a->ptr.pp_complex[ia][ja], a->stride, &u->ptr.p_complex[uoffs], &v->ptr.p_complex[voffs]);
}

bool _ialglib_i_rmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs) {
   return _ialglib_rmatrixrank1(m, n, &a->ptr.pp_double[ia][ja], a->stride, &u->ptr.p_double[uoffs], &v->ptr.p_double[voffs]);
}

bool _ialglib_i_rmatrixgerf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, double alpha, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs) {
   return _ialglib_rmatrixger(m, n, &a->ptr.pp_double[ia][ja], a->stride, alpha, &u->ptr.p_double[uoffs], &v->ptr.p_double[voffs]);
}

// This function reads rectangular matrix A given by two column pointers
// col0 and col1 and stride src_stride and moves it into contiguous row-
// by-row storage given by dst.
//
// It can handle following special cases:
// * col1 == NULL    in this case second column of A is filled by zeros
void _ialglib_pack_n2(double *col0, double *col1, ae_int_t n, ae_int_t src_stride, double *dst) {
   ae_int_t n2, j, stride2;
// handle special case
   if (col1 == NULL) {
      for (j = 0; j < n; j++) {
         dst[0] = *col0;
         dst[1] = 0.0;
         col0 += src_stride;
         dst += 2;
      }
      return;
   }
// handle general case
   n2 = n / 2;
   stride2 = src_stride * 2;
   for (j = 0; j < n2; j++) {
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

// This function reads rectangular matrix A given by two column pointers col0
// and  col1  and  stride src_stride and moves it into  contiguous row-by-row
// storage given by dst.
//
// dst must be aligned, col0 and col1 may be non-aligned.
//
// It can handle following special cases:
// * col1 == NULL        in this case second column of A is filled by zeros
// * src_stride == 1     efficient SSE-based code is used
// * col1-col0 == 1      efficient SSE-based code is used
//
// This function supports SSE2; it can be used when:
// 1. AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// 2. CurCPU contains CPU_SSE2 (checked at run-time)
//
// If  you  want  to  know  whether  it  is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
#if defined AE_HAS_SSE2_INTRINSICS
void _ialglib_pack_n2_sse2(double *col0, double *col1, ae_int_t n, ae_int_t src_stride, double *dst) {
   ae_int_t n2, j, stride2;
// handle special case: col1 == NULL
   if (col1 == NULL) {
      for (j = 0; j < n; j++) {
         dst[0] = *col0;
         dst[1] = 0.0;
         col0 += src_stride;
         dst += 2;
      }
      return;
   }
// handle unit stride
   if (src_stride == 1) {
      __m128d v0, v1;
      n2 = n / 2;
      for (j = 0; j < n2; j++) {
         v0 = _mm_loadu_pd(col0);
         col0 += 2;
         v1 = _mm_loadu_pd(col1);
         col1 += 2;
         _mm_store_pd(dst, _mm_unpacklo_pd(v0, v1));
         _mm_store_pd(dst + 2, _mm_unpackhi_pd(v0, v1));
         dst += 4;
      }
      if (n % 2) {
         dst[0] = *col0;
         dst[1] = *col1;
      }
      return;
   }
// handle col1-col0 == 1
   if (col1 - col0 == 1) {
      __m128d v0, v1;
      n2 = n / 2;
      stride2 = 2 * src_stride;
      for (j = 0; j < n2; j++) {
         v0 = _mm_loadu_pd(col0);
         v1 = _mm_loadu_pd(col0 + src_stride);
         _mm_store_pd(dst, v0);
         _mm_store_pd(dst + 2, v1);
         col0 += stride2;
         dst += 4;
      }
      if (n % 2) {
         dst[0] = col0[0];
         dst[1] = col0[1];
      }
      return;
   }
// handle general case
   n2 = n / 2;
   stride2 = src_stride * 2;
   for (j = 0; j < n2; j++) {
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
#endif

// This function calculates R := alpha*A'*B+beta*R where A and B are Kx2
// matrices stored in contiguous row-by-row storage,  R  is  2x2  matrix
// stored in non-contiguous row-by-row storage.
//
// A and B must be aligned; R may be non-aligned.
//
// If beta is zero, contents of R is ignored (not  multiplied  by zero -
// just ignored).
//
// However, when alpha is zero, we still calculate A'*B, which is
// multiplied by zero afterwards.
//
// Function accepts additional parameter store_mode:
// * if 0, full R is stored
// * if 1, only first row of R is stored
// * if 2, only first column of R is stored
// * if 3, only top left element of R is stored
void _ialglib_mm22(double alpha, const double *a, const double *b, ae_int_t k, double beta, double *r, ae_int_t stride, ae_int_t store_mode) {
   double v00, v01, v10, v11;
   ae_int_t t;
   v00 = 0.0;
   v01 = 0.0;
   v10 = 0.0;
   v11 = 0.0;
   for (t = 0; t < k; t++) {
      v00 += a[0] * b[0];
      v01 += a[0] * b[1];
      v10 += a[1] * b[0];
      v11 += a[1] * b[1];
      a += 2;
      b += 2;
   }
   if (store_mode == 0) {
      if (beta == 0) {
         r[0] = alpha * v00;
         r[1] = alpha * v01;
         r[stride] = alpha * v10;
         r[stride + 1] = alpha * v11;
      } else {
         r[0] = beta * r[0] + alpha * v00;
         r[1] = beta * r[1] + alpha * v01;
         r[stride] = beta * r[stride] + alpha * v10;
         r[stride + 1] = beta * r[stride + 1] + alpha * v11;
      }
      return;
   }
   if (store_mode == 1) {
      if (beta == 0) {
         r[0] = alpha * v00;
         r[1] = alpha * v01;
      } else {
         r[0] = beta * r[0] + alpha * v00;
         r[1] = beta * r[1] + alpha * v01;
      }
      return;
   }
   if (store_mode == 2) {
      if (beta == 0) {
         r[0] = alpha * v00;
         r[stride] = alpha * v10;
      } else {
         r[0] = beta * r[0] + alpha * v00;
         r[stride] = beta * r[stride] + alpha * v10;
      }
      return;
   }
   if (store_mode == 3) {
      if (beta == 0) {
         r[0] = alpha * v00;
      } else {
         r[0] = beta * r[0] + alpha * v00;
      }
      return;
   }
}

// This function calculates R := alpha*A'*B+beta*R where A and B are Kx2
// matrices stored in contiguous row-by-row storage,  R  is  2x2  matrix
// stored in non-contiguous row-by-row storage.
//
// A and B must be aligned; R may be non-aligned.
//
// If beta is zero, contents of R is ignored (not  multiplied  by zero -
// just ignored).
//
// However, when alpha is zero, we still calculate A'*B, which is
// multiplied by zero afterwards.
//
// Function accepts additional parameter store_mode:
// * if 0, full R is stored
// * if 1, only first row of R is stored
// * if 2, only first column of R is stored
// * if 3, only top left element of R is stored
//
// This function supports SSE2; it can be used when:
// 1. AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// 2. CurCPU contains CPU_SSE2 (checked at run-time)
//
// If (1) is failed, this function will still be defined and callable, but it
// will do nothing.  If (2)  is  failed , call to this function will probably
// crash your system.
//
// If  you  want  to  know  whether  it  is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
#if defined AE_HAS_SSE2_INTRINSICS
void _ialglib_mm22_sse2(double alpha, const double *a, const double *b, ae_int_t k, double beta, double *r, ae_int_t stride, ae_int_t store_mode) {
// We calculate product of two Kx2 matrices (result is 2x2).
// VA and VB store result as follows:
//
//        [ VD[0]  VE[0] ]
// A'*B = [              ]
//        [ VE[1]  VD[1] ]
//
   __m128d va, vb, vd, ve, vt, r0, r1, valpha, vbeta;
   ae_int_t t, k2;
// calculate product
   k2 = k / 2;
   vd = _mm_setzero_pd();
   ve = _mm_setzero_pd();
   for (t = 0; t < k2; t++) {
      vb = _mm_load_pd(b);
      va = _mm_load_pd(a);
      vt = vb;
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
      va = _mm_load_pd(a);
      vb = _mm_load_pd(b);
      vt = _mm_shuffle_pd(vb, vb, 1);
      vd = _mm_add_pd(_mm_mul_pd(va, vb), vd);
      ve = _mm_add_pd(_mm_mul_pd(va, vt), ve);
   }
// r0 is first row of alpha*A'*B, r1 is second row
   valpha = _mm_load1_pd(&alpha);
   r0 = _mm_mul_pd(_mm_unpacklo_pd(vd, ve), valpha);
   r1 = _mm_mul_pd(_mm_unpackhi_pd(ve, vd), valpha);
// store
   if (store_mode == 0) {
      if (beta == 0) {
         _mm_storeu_pd(r, r0);
         _mm_storeu_pd(r + stride, r1);
      } else {
         vbeta = _mm_load1_pd(&beta);
         _mm_storeu_pd(r, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r), vbeta), r0));
         _mm_storeu_pd(r + stride, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + stride), vbeta), r1));
      }
      return;
   }
   if (store_mode == 1) {
      if (beta == 0)
         _mm_storeu_pd(r, r0);
      else
         _mm_storeu_pd(r, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r), _mm_load1_pd(&beta)), r0));
      return;
   }
   if (store_mode == 2) {
      double buf[4];
      _mm_storeu_pd(buf, r0);
      _mm_storeu_pd(buf + 2, r1);
      if (beta == 0) {
         r[0] = buf[0];
         r[stride] = buf[2];
      } else {
         r[0] = beta * r[0] + buf[0];
         r[stride] = beta * r[stride] + buf[2];
      }
      return;
   }
   if (store_mode == 3) {
      double buf[2];
      _mm_storeu_pd(buf, r0);
      if (beta == 0)
         r[0] = buf[0];
      else
         r[0] = beta * r[0] + buf[0];
      return;
   }
}
#endif

// This function calculates R := alpha*A'*(B0|B1)+beta*R where A, B0  and  B1
// are Kx2 matrices stored in contiguous row-by-row storage, R is 2x4  matrix
// stored in non-contiguous row-by-row storage.
//
// A, B0 and B1 must be aligned; R may be non-aligned.
//
// Note  that  B0  and  B1  are  two  separate  matrices  stored in different
// locations.
//
// If beta is zero, contents of R is ignored (not  multiplied  by zero - just
// ignored).
//
// However,  when  alpha  is  zero , we still calculate MM product,  which is
// multiplied by zero afterwards.
//
// Unlike mm22 functions, this function does NOT support partial  output of R
// - we always store full 2x4 matrix.
void _ialglib_mm22x2(double alpha, const double *a, const double *b0, const double *b1, ae_int_t k, double beta, double *r, ae_int_t stride) {
   _ialglib_mm22(alpha, a, b0, k, beta, r, stride, 0);
   _ialglib_mm22(alpha, a, b1, k, beta, r + 2, stride, 0);
}

// This function calculates R := alpha*A'*(B0|B1)+beta*R where A, B0  and  B1
// are Kx2 matrices stored in contiguous row-by-row storage, R is 2x4  matrix
// stored in non-contiguous row-by-row storage.
//
// A, B0 and B1 must be aligned; R may be non-aligned.
//
// Note  that  B0  and  B1  are  two  separate  matrices  stored in different
// locations.
//
// If beta is zero, contents of R is ignored (not  multiplied  by zero - just
// ignored).
//
// However,  when  alpha  is  zero , we still calculate MM product,  which is
// multiplied by zero afterwards.
//
// Unlike mm22 functions, this function does NOT support partial  output of R
// - we always store full 2x4 matrix.
//
// This function supports SSE2; it can be used when:
// 1. AE_HAS_SSE2_INTRINSICS was defined (checked at compile-time)
// 2. CurCPU contains CPU_SSE2 (checked at run-time)
//
// If (1) is failed, this function will still be defined and callable, but it
// will do nothing.  If (2)  is  failed , call to this function will probably
// crash your system.
//
// If  you  want  to  know  whether  it  is safe to call it, you should check CurCPU.
// If CPU_SSE2 bit is set, this function is callable and will do its work.
#if defined AE_HAS_SSE2_INTRINSICS
void _ialglib_mm22x2_sse2(double alpha, const double *a, const double *b0, const double *b1, ae_int_t k, double beta, double *r, ae_int_t stride) {
// We calculate product of two Kx2 matrices (result is 2x2).
// V0, V1, V2, V3 store result as follows:
//
//     [ V0[0]  V1[1] V2[0]  V3[1] ]
// R = [                           ]
//     [ V1[0]  V0[1] V3[0]  V2[1] ]
//
// VA0 stores current 1x2 block of A, VA1 stores shuffle of VA0,
// VB0 and VB1 are used to store two copies of 1x2 block of B0 or B1
// (both vars store same data - either B0 or B1). Results from multiplication
// by VA0/VA1 are stored in VB0/VB1 too.
//
   __m128d v0, v1, v2, v3, va0, va1, vb0, vb1;
   __m128d r00, r01, r10, r11, valpha, vbeta;
   ae_int_t t;
   v0 = _mm_setzero_pd();
   v1 = _mm_setzero_pd();
   v2 = _mm_setzero_pd();
   v3 = _mm_setzero_pd();
   for (t = 0; t < k; t++) {
      va0 = _mm_load_pd(a);
      vb0 = _mm_load_pd(b0);
      va1 = _mm_load_pd(a);
      vb0 = _mm_mul_pd(va0, vb0);
      vb1 = _mm_load_pd(b0);
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
// shuffle V1 and V3 (conversion to more convenient storage format):
//
//     [ V0[0]  V1[0] V2[0]  V3[0] ]
// R = [                           ]
//     [ V1[1]  V0[1] V3[1]  V2[1] ]
//
// unpack results to
//
// [ r00 r01 ]
// [ r10 r11 ]
//
   valpha = _mm_load1_pd(&alpha);
   v1 = _mm_shuffle_pd(v1, v1, 1);
   v3 = _mm_shuffle_pd(v3, v3, 1);
   r00 = _mm_mul_pd(_mm_unpacklo_pd(v0, v1), valpha);
   r10 = _mm_mul_pd(_mm_unpackhi_pd(v1, v0), valpha);
   r01 = _mm_mul_pd(_mm_unpacklo_pd(v2, v3), valpha);
   r11 = _mm_mul_pd(_mm_unpackhi_pd(v3, v2), valpha);
// store
   if (beta == 0) {
      _mm_storeu_pd(r, r00);
      _mm_storeu_pd(r + 2, r01);
      _mm_storeu_pd(r + stride, r10);
      _mm_storeu_pd(r + stride + 2, r11);
   } else {
      vbeta = _mm_load1_pd(&beta);
      _mm_storeu_pd(r, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r), vbeta), r00));
      _mm_storeu_pd(r + 2, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + 2), vbeta), r01));
      _mm_storeu_pd(r + stride, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + stride), vbeta), r10));
      _mm_storeu_pd(r + stride + 2, _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(r + stride + 2), vbeta), r11));
   }
}
#endif
} // end of namespace alglib_impl

namespace alglib {
// Internal forwards
ae_int_t my_stricmp(const char *s1, const char *s2);
char *filter_spaces(const char *s);
void str_vector_create(const char *src, bool match_head_only, std::vector<const char *> *p_vec);
void str_matrix_create(const char *src, std::vector< std::vector<const char *> > *p_mat);
bool parse_bool_delim(const char *s, const char *delim);
ae_int_t parse_int_delim(const char *s, const char *delim);
bool _parse_real_delim(const char *s, const char *delim, double *result, const char **new_s);
double parse_real_delim(const char *s, const char *delim);
complex parse_complex_delim(const char *s, const char *delim);
std::string arraytostring(const bool *ptr, ae_int_t n);
std::string arraytostring(const ae_int_t *ptr, ae_int_t n);
std::string arraytostring(const double *ptr, ae_int_t n, int dps);
std::string arraytostring(const complex *ptr, ae_int_t n, int dps);

// Global and local constants/variables
const double machineepsilon = 5.0E-16, maxrealnumber = 1.0E300, minrealnumber = 1.0E-300;
const ae_int_t ByteOrder = alglib_impl::ByteOrder;

// Exception handling
#if !defined AE_NO_EXCEPTIONS
ap_error::ap_error() { msg = alglib_impl::CurMsg; }
ap_error::ap_error(const char *Msg) { msg = Msg; }
void ap_error::make_assertion(bool Cond) { if (!Cond) ThrowError(""); }
void ap_error::make_assertion(bool Cond, const char *Msg) { if (!Cond) ThrowError(Msg); }
#else
static const char *_alglib_last_error = NULL;
static void set_error_flag(const char *Msg) { _alglib_last_error = Msg == NULL? "ALGLIB: unknown error": Msg; }
void set_error_msg() { _alglib_last_error = alglib_impl::CurMsg; }

bool get_error_flag(const char **MsgP) {
   if (_alglib_last_error == NULL) return false;
   if (MsgP != NULL) *MsgP = _alglib_last_error;
   return true;
}

void clear_error_flag() { _alglib_last_error = NULL; }
#endif

const xparams xdefault = 0x0, serial = _ALGLIB_FLG_THREADING_SERIAL, parallel = _ALGLIB_FLG_THREADING_PARALLEL;

// Complex number with double precision.
complex &complex::operator=(const double &v) { x = v, y = 0.0; return *this; }
complex &complex::operator+=(const double &v) { x += v; return *this; }
complex &complex::operator-=(const double &v) { x -= v; return *this; }
complex &complex::operator*=(const double &v) { x *= v, y *= v; return *this; }
complex &complex::operator/=(const double &v) { x /= v, y /= v; return *this; }
complex &complex::operator=(const complex &z) { x = z.x, y = z.y; return *this; }
complex &complex::operator+=(const complex &z) { x += z.x, y += z.y; return *this; }
complex &complex::operator-=(const complex &z) { x -= z.x, y -= z.y; return *this; }
complex &complex::operator*=(const complex &z) { double t = x*z.x - y*z.y; y = x*z.y + y*z.x, x = t; return *this; }

complex &complex::operator/=(const complex &z) {
   double x2 = z.x, y2 = z.y;
   double e;
   double f;
   if (fabs(y2) < fabs(x2)) {
      e = y2 / x2;
      f = x2 + y2 * e;
      return *this = complex((x + y * e) / f, (y - x * e) / f);
   } else {
      e = x2 / y2;
      f = y2 + x2 * e;
      return *this = complex((x * e + y) / f, (y * e - x) / f);
   }
}

// alglib_impl-alglib gateway
static inline alglib_impl::ae_complex ae_complex_from_c(complex z) { return alglib_impl::ae_complex_from_d(z.x, z.y); }

#if !defined AE_NO_EXCEPTIONS
std::string complex::tostring(int _dps) const {
   char mask[32];
   char buf_x[32];
   char buf_y[32];
   char buf_zero[32];
   int dps = _dps >= 0 ? _dps : -_dps;
   if (dps <= 0 || dps >= 20)
      ThrowError("complex::tostring(): incorrect dps");
// handle IEEE special quantities
   if (isnan(x) || isnan(y))
      return "NAN";
   if (isinf(x) || isinf(y))
      return "INF";
// generate mask
   if (sprintf(mask, "%%.%d%s", dps, _dps >= 0 ? "f" : "e") >= (int)sizeof(mask))
      ThrowError("complex::tostring(): buffer overflow");
// print |x|, |y| and zero with same mask and compare
   if (sprintf(buf_x, mask, (double)(fabs(x))) >= (int)sizeof(buf_x))
      ThrowError("complex::tostring(): buffer overflow");
   if (sprintf(buf_y, mask, (double)(fabs(y))) >= (int)sizeof(buf_y))
      ThrowError("complex::tostring(): buffer overflow");
   if (sprintf(buf_zero, mask, 0.0) >= (int)sizeof(buf_zero))
      ThrowError("complex::tostring(): buffer overflow");
// different zero/nonzero patterns
   if (strcmp(buf_x, buf_zero) != 0 && strcmp(buf_y, buf_zero) != 0)
      return std::string(x > 0 ? "" : "-") + buf_x + (y > 0 ? "+" : "-") + buf_y + "i";
   if (strcmp(buf_x, buf_zero) != 0 && strcmp(buf_y, buf_zero) == 0)
      return std::string(x > 0 ? "" : "-") + buf_x;
   if (strcmp(buf_x, buf_zero) == 0 && strcmp(buf_y, buf_zero) != 0)
      return std::string(y > 0 ? "" : "-") + buf_y + "i";
   return std::string("0");
}
#endif

bool operator==(const complex &lhs, const complex &rhs) {
   volatile double x1 = lhs.x;
   volatile double x2 = rhs.x;
   volatile double y1 = lhs.y;
   volatile double y2 = rhs.y;
   return x1 == x2 && y1 == y2;
}

bool operator!=(const complex &lhs, const complex &rhs) { return !(lhs == rhs); }
const complex operator+(const complex &lhs) { return lhs; }
const complex operator-(const complex &lhs) { return complex(-lhs.x, -lhs.y); }
const complex operator+(const complex &lhs, const complex &rhs) { complex r = lhs; r += rhs; return r; }
const complex operator+(const complex &lhs, const double &rhs) { complex r = lhs; r += rhs; return r; }
const complex operator+(const double &lhs, const complex &rhs) { complex r = rhs; r += lhs; return r; }
const complex operator-(const complex &lhs, const complex &rhs) { complex r = lhs; r -= rhs; return r; }
const complex operator-(const complex &lhs, const double &rhs) { complex r = lhs; r -= rhs; return r; }
const complex operator-(const double &lhs, const complex &rhs) { complex r = lhs; r -= rhs; return r; }
const complex operator*(const complex &lhs, const complex &rhs) { return complex(lhs.x*rhs.x - lhs.y*rhs.y, lhs.x*rhs.y + lhs.y*rhs.x); }
const complex operator*(const complex &lhs, const double &rhs) { return complex(lhs.x*rhs, lhs.y*rhs); }
const complex operator*(const double &lhs, const complex &rhs) { return complex(lhs*rhs.x, lhs*rhs.y); }

const complex operator/(const complex &lhs, const complex &rhs) {
   double x1 = lhs.x, y1 = lhs.y, x2 = rhs.x, y2 = rhs.y;
   double e;
   double f;
   if (fabs(y2) < fabs(x2)) {
      e = y2 / x2;
      f = x2 + y2 * e;
      return complex((x1 + y1 * e) / f, (y1 - x1 * e) / f);
   } else {
      e = x2 / y2;
      f = y2 + x2 * e;
      return complex((x1 * e + y1) / f, (y1 * e - x1) / f);
   }
}

const complex operator/(const double &lhs, const complex &rhs) {
   double x2 = rhs.x, y2 = rhs.y;
   double e;
   double f;
   if (fabs(y2) < fabs(x2)) {
      e = y2 / x2;
      f = x2 + y2 * e;
      return complex(lhs / f, -lhs * e / f);
   } else {
      e = x2 / y2;
      f = y2 + x2 * e;
      return complex(lhs * e / f, -lhs / f);
   }
}

const complex operator/(const complex &lhs, const double &rhs) { return complex(lhs.x / rhs, lhs.y / rhs); }

double abscomplex(const complex &z) {
   double w;
   double xabs;
   double yabs;
   double v;
   xabs = fabs(z.x);
   yabs = fabs(z.y);
   w = xabs > yabs ? xabs : yabs;
   v = xabs < yabs ? xabs : yabs;
   if (v == 0)
      return w;
   else {
      double t = v / w;
      return w * sqrt(1 + t * t);
   }
}

complex conj(const complex &z) {
   return complex(z.x, -z.y);
}

complex csqr(const complex &z) {
   return complex(z.x * z.x - z.y * z.y, 2 * z.x * z.y);
}

void setnworkers(ae_int_t nworkers) {
#ifdef AE_HPC
   alglib_impl::ae_set_cores_to_use(nworkers);
#endif
}

void setglobalthreading(const xparams settings) {
#ifdef AE_HPC
   alglib_impl::ae_set_global_threading(settings);
#endif
}

ae_int_t getnworkers() {
#ifdef AE_HPC
   return alglib_impl::ae_get_cores_to_use();
#else
   return 1;
#endif
}

ae_int_t _ae_cores_count() {
#ifdef AE_HPC
   return alglib_impl::ae_cores_count();
#else
   return 1;
#endif
}

void _ae_set_global_threading(alglib_impl::ae_uint64_t flg_value) {
#ifdef AE_HPC
   alglib_impl::ae_set_global_threading(flg_value);
#endif
}

alglib_impl::ae_uint64_t _ae_get_global_threading() {
#ifdef AE_HPC
   return alglib_impl::ae_get_global_threading();
#else
   return _ALGLIB_FLG_THREADING_SERIAL;
#endif
}

// Level 1 BLAS functions
double vdotproduct(const double *v0, ae_int_t stride0, const double *v1, ae_int_t stride1, ae_int_t n) {
   double result = 0.0;
   ae_int_t i;
   if (stride0 != 1 || stride1 != 1) {
   //
   // slow general code
   //
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1)
         result += (*v0) * (*v1);
   } else {
   //
   // optimized code for stride=1
   //
      ae_int_t n4 = n / 4;
      ae_int_t nleft = n % 4;
      for (i = 0; i < n4; i++, v0 += 4, v1 += 4)
         result += v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] + v0[3] * v1[3];
      for (i = 0; i < nleft; i++, v0++, v1++)
         result += v0[0] * v1[0];
   }
   return result;
}

double vdotproduct(const double *v1, const double *v2, ae_int_t N) {
   return vdotproduct(v1, 1, v2, 1, N);
}

complex vdotproduct(const complex *v0, ae_int_t stride0, const char *conj0, const complex *v1, ae_int_t stride1, const char *conj1, ae_int_t n) {
   double rx = 0.0, ry = 0.0;
   ae_int_t i;
   bool bconj0 = !(conj0[0] == 'N' || conj0[0] == 'n');
   bool bconj1 = !(conj1[0] == 'N' || conj1[0] == 'n');
   if (bconj0 && bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = -v0->y;
         v1x = v1->x;
         v1y = -v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   if (!bconj0 && bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = v0->y;
         v1x = v1->x;
         v1y = -v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   if (bconj0 && !bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = -v0->y;
         v1x = v1->x;
         v1y = v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   if (!bconj0 && !bconj1) {
      double v0x, v0y, v1x, v1y;
      for (i = 0; i < n; i++, v0 += stride0, v1 += stride1) {
         v0x = v0->x;
         v0y = v0->y;
         v1x = v1->x;
         v1y = v1->y;
         rx += v0x * v1x - v0y * v1y;
         ry += v0x * v1y + v0y * v1x;
      }
   }
   return complex(rx, ry);
}

complex vdotproduct(const complex *v1, const complex *v2, ae_int_t N) {
   return vdotproduct(v1, 1, "N", v2, 1, "N", N);
}

void vmove(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst = *vsrc;
   } else {
   //
   // optimized case
   //
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] = vsrc[0];
         vdst[1] = vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] = vsrc[0];
   }
}

void vmove(double *vdst, const double *vsrc, ae_int_t N) {
   vmove(vdst, 1, vsrc, 1, N);
}

void vmove(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = vsrc->x;
            vdst->y = -vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
            *vdst = *vsrc;
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = vsrc->x;
            vdst->y = -vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++)
            *vdst = *vsrc;
      }
   }
}

void vmove(complex *vdst, const complex *vsrc, ae_int_t N) {
   vmove(vdst, 1, vsrc, 1, "N", N);
}

void vmoveneg(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst = -*vsrc;
   } else {
   //
   // optimized case
   //
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] = -vsrc[0];
         vdst[1] = -vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] = -vsrc[0];
   }
}

void vmoveneg(double *vdst, const double *vsrc, ae_int_t N) {
   vmoveneg(vdst, 1, vsrc, 1, N);
}

void vmoveneg(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = -vsrc->x;
            vdst->y = vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = -vsrc->x;
            vdst->y = -vsrc->y;
         }
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = -vsrc->x;
            vdst->y = vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = -vsrc->x;
            vdst->y = -vsrc->y;
         }
      }
   }
}

void vmoveneg(complex *vdst, const complex *vsrc, ae_int_t N) {
   vmoveneg(vdst, 1, vsrc, 1, "N", N);
}

void vmove(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst = alpha * (*vsrc);
   } else {
   //
   // optimized case
   //
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] = alpha * vsrc[0];
         vdst[1] = alpha * vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] = alpha * vsrc[0];
   }
}

void vmove(double *vdst, const double *vsrc, ae_int_t N, double alpha) {
   vmove(vdst, 1, vsrc, 1, N, alpha);
}

void vmove(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = alpha * vsrc->x;
            vdst->y = -alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = alpha * vsrc->x;
            vdst->y = alpha * vsrc->y;
         }
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = alpha * vsrc->x;
            vdst->y = -alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = alpha * vsrc->x;
            vdst->y = alpha * vsrc->y;
         }
      }
   }
}

void vmove(complex *vdst, const complex *vsrc, ae_int_t N, double alpha) {
   vmove(vdst, 1, vsrc, 1, "N", N, alpha);
}

void vmove(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, complex alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = ax * vsrc->x + ay * vsrc->y;
            vdst->y = -ax * vsrc->y + ay * vsrc->x;
         }
      } else {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x = ax * vsrc->x - ay * vsrc->y;
            vdst->y = ax * vsrc->y + ay * vsrc->x;
         }
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = ax * vsrc->x + ay * vsrc->y;
            vdst->y = -ax * vsrc->y + ay * vsrc->x;
         }
      } else {
         double ax = alpha.x, ay = alpha.y;
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x = ax * vsrc->x - ay * vsrc->y;
            vdst->y = ax * vsrc->y + ay * vsrc->x;
         }
      }
   }
}

void vmove(complex *vdst, const complex *vsrc, ae_int_t N, complex alpha) {
   vmove(vdst, 1, vsrc, 1, "N", N, alpha);
}

void vadd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst += *vsrc;
   } else {
   //
   // optimized case
   //
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] += vsrc[0];
         vdst[1] += vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] += vsrc[0];
   }
}

void vadd(double *vdst, const double *vsrc, ae_int_t N) {
   vadd(vdst, 1, vsrc, 1, N);
}

void vadd(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += vsrc->x;
            vdst->y -= vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += vsrc->x;
            vdst->y += vsrc->y;
         }
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += vsrc->x;
            vdst->y -= vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += vsrc->x;
            vdst->y += vsrc->y;
         }
      }
   }
}

void vadd(complex *vdst, const complex *vsrc, ae_int_t N) {
   vadd(vdst, 1, vsrc, 1, "N", N);
}

void vadd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst += alpha * (*vsrc);
   } else {
   //
   // optimized case
   //
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] += alpha * vsrc[0];
         vdst[1] += alpha * vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] += alpha * vsrc[0];
   }
}

void vadd(double *vdst, const double *vsrc, ae_int_t N, double alpha) {
   vadd(vdst, 1, vsrc, 1, N, alpha);
}

void vadd(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += alpha * vsrc->x;
            vdst->y -= alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += alpha * vsrc->x;
            vdst->y += alpha * vsrc->y;
         }
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += alpha * vsrc->x;
            vdst->y -= alpha * vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += alpha * vsrc->x;
            vdst->y += alpha * vsrc->y;
         }
      }
   }
}

void vadd(complex *vdst, const complex *vsrc, ae_int_t N, double alpha) {
   vadd(vdst, 1, vsrc, 1, "N", N, alpha);
}

void vadd(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, complex alpha) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      double ax = alpha.x, ay = alpha.y;
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += ax * vsrc->x + ay * vsrc->y;
            vdst->y -= ax * vsrc->y - ay * vsrc->x;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x += ax * vsrc->x - ay * vsrc->y;
            vdst->y += ax * vsrc->y + ay * vsrc->x;
         }
      }
   } else {
   //
   // optimized case
   //
      double ax = alpha.x, ay = alpha.y;
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += ax * vsrc->x + ay * vsrc->y;
            vdst->y -= ax * vsrc->y - ay * vsrc->x;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x += ax * vsrc->x - ay * vsrc->y;
            vdst->y += ax * vsrc->y + ay * vsrc->x;
         }
      }
   }
}

void vadd(complex *vdst, const complex *vsrc, ae_int_t N, complex alpha) {
   vadd(vdst, 1, vsrc, 1, "N", N, alpha);
}

void vsub(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n) {
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src)
         *vdst -= *vsrc;
   } else {
   //
   // optimized case
   //
      ae_int_t n2 = n / 2;
      for (i = 0; i < n2; i++, vdst += 2, vsrc += 2) {
         vdst[0] -= vsrc[0];
         vdst[1] -= vsrc[1];
      }
      if (n % 2 != 0)
         vdst[0] -= vsrc[0];
   }
}

void vsub(double *vdst, const double *vsrc, ae_int_t N) {
   vsub(vdst, 1, vsrc, 1, N);
}

void vsub(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n) {
   bool bconj = !(conj_src[0] == 'N' || conj_src[0] == 'n');
   ae_int_t i;
   if (stride_dst != 1 || stride_src != 1) {
   //
   // general unoptimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x -= vsrc->x;
            vdst->y += vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst += stride_dst, vsrc += stride_src) {
            vdst->x -= vsrc->x;
            vdst->y -= vsrc->y;
         }
      }
   } else {
   //
   // optimized case
   //
      if (bconj) {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x -= vsrc->x;
            vdst->y += vsrc->y;
         }
      } else {
         for (i = 0; i < n; i++, vdst++, vsrc++) {
            vdst->x -= vsrc->x;
            vdst->y -= vsrc->y;
         }
      }
   }
}

void vsub(complex *vdst, const complex *vsrc, ae_int_t N) {
   vsub(vdst, 1, vsrc, 1, "N", N);
}

void vsub(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha) {
   vadd(vdst, stride_dst, vsrc, stride_src, n, -alpha);
}

void vsub(double *vdst, const double *vsrc, ae_int_t N, double alpha) {
   vadd(vdst, 1, vsrc, 1, N, -alpha);
}

void vsub(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha) {
   vadd(vdst, stride_dst, vsrc, stride_src, conj_src, n, -alpha);
}

void vsub(complex *vdst, const complex *vsrc, ae_int_t n, double alpha) {
   vadd(vdst, 1, vsrc, 1, "N", n, -alpha);
}

void vsub(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, complex alpha) {
   vadd(vdst, stride_dst, vsrc, stride_src, conj_src, n, -alpha);
}

void vsub(complex *vdst, const complex *vsrc, ae_int_t n, complex alpha) {
   vadd(vdst, 1, vsrc, 1, "N", n, -alpha);
}

void vmul(double *vdst, ae_int_t stride_dst, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst)
         *vdst *= alpha;
   } else {
   //
   // optimized case
   //
      for (i = 0; i < n; i++, vdst++)
         *vdst *= alpha;
   }
}

void vmul(double *vdst, ae_int_t N, double alpha) {
   vmul(vdst, 1, N, alpha);
}

void vmul(complex *vdst, ae_int_t stride_dst, ae_int_t n, double alpha) {
   ae_int_t i;
   if (stride_dst != 1) {
   //
   // general unoptimized case
   //
      for (i = 0; i < n; i++, vdst += stride_dst) {
         vdst->x *= alpha;
         vdst->y *= alpha;
      }
   } else {
   //
   // optimized case
   //
      for (i = 0; i < n; i++, vdst++) {
         vdst->x *= alpha;
         vdst->y *= alpha;
      }
   }
}

void vmul(complex *vdst, ae_int_t N, double alpha) {
   vmul(vdst, 1, N, alpha);
}

void vmul(complex *vdst, ae_int_t stride_dst, ae_int_t n, complex alpha) {
   ae_int_t i;
   if (stride_dst != 1) {
   //
   // general unoptimized case
   //
      double ax = alpha.x, ay = alpha.y;
      for (i = 0; i < n; i++, vdst += stride_dst) {
         double dstx = vdst->x, dsty = vdst->y;
         vdst->x = ax * dstx - ay * dsty;
         vdst->y = ax * dsty + ay * dstx;
      }
   } else {
   //
   // optimized case
   //
      double ax = alpha.x, ay = alpha.y;
      for (i = 0; i < n; i++, vdst++) {
         double dstx = vdst->x, dsty = vdst->y;
         vdst->x = ax * dstx - ay * dsty;
         vdst->y = ax * dsty + ay * dstx;
      }
   }
}

void vmul(complex *vdst, ae_int_t N, complex alpha) {
   vmul(vdst, 1, N, alpha);
}

// Matrices and vectors
ae_vector_wrapper::ae_vector_wrapper(alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      is_frozen_proxy = false, This = NULL, set_error_msg(); return;
#endif
   }
   is_frozen_proxy = false, This = &Obj, memset(This, 0, sizeof *This), ae_vector_init(This, 0, datatype, false);
   alglib_impl::ae_state_clear();
}

ae_vector_wrapper::ae_vector_wrapper(alglib_impl::ae_vector *e_ptr, alglib_impl::ae_datatype datatype) {
   if (e_ptr == NULL || e_ptr->datatype != datatype) {
      const char *msg = "ALGLIB: ae_vector_wrapper datatype check failed";
#if !defined AE_NO_EXCEPTIONS
      ThrowError(msg);
#else
      is_frozen_proxy = false, This = NULL, set_error_flag(msg); return;
#endif
   }
   is_frozen_proxy = true, This = e_ptr;
}

ae_vector_wrapper::ae_vector_wrapper(const ae_vector_wrapper &rhs, alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      is_frozen_proxy = false, This = NULL, set_error_msg(); return;
#endif
   }
   alglib_impl::ae_assert(rhs.This != NULL, "ALGLIB: ae_vector_wrapper source is not initialized");
   alglib_impl::ae_assert(rhs.This->datatype == datatype, "ALGLIB: ae_vector_wrapper datatype check failed");
   is_frozen_proxy = false, This = &Obj, memset(This, 0, sizeof *This), ae_vector_copy(This, rhs.This, false);
   alglib_impl::ae_state_clear();
}

ae_vector_wrapper::~ae_vector_wrapper() { if (This == &Obj) ae_vector_free(This, false); }

void ae_vector_wrapper::setlength(ae_int_t iLen) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ae_assert(This != NULL, "ALGLIB: setlength() error, This == NULL (array was not correctly initialized)");
   alglib_impl::ae_assert(!is_frozen_proxy, "ALGLIB: setlength() error, This is frozen proxy array");
   alglib_impl::ae_vector_set_length(This, iLen);
   alglib_impl::ae_state_clear();
}

ae_int_t ae_vector_wrapper::length() const { return This == NULL? 0: This->cnt; }

void ae_vector_wrapper::attach_to(alglib_impl::x_vector *new_ptr) {
   if (This == &Obj) ae_vector_free(This, false);
   is_frozen_proxy = true, This = &Obj, memset(This, 0, sizeof *This), ae_vector_init_attach_to_x(This, new_ptr, false);
}

const ae_vector_wrapper &ae_vector_wrapper::assign(const ae_vector_wrapper &rhs) {
   if (this == &rhs) return *this;
   alglib_impl::ae_state_init();
   TryCatch(*this)
   alglib_impl::ae_assert(This != NULL, "ALGLIB: incorrect assignment (uninitialized destination)");
   alglib_impl::ae_assert(rhs.This != NULL, "ALGLIB: incorrect assignment (uninitialized source)");
   alglib_impl::ae_assert(rhs.This->datatype == This->datatype, "ALGLIB: incorrect assignment to array (types do not match)");
   if (is_frozen_proxy) alglib_impl::ae_assert(rhs.This->cnt == This->cnt, "ALGLIB: incorrect assignment to proxy array (sizes do not match)");
   if (rhs.This->cnt != This->cnt) ae_vector_set_length(This, rhs.This->cnt);
   memcpy(This->ptr.p_ptr, rhs.This->ptr.p_ptr, This->cnt*alglib_impl::ae_sizeof(This->datatype));
   alglib_impl::ae_state_clear();
   return *this;
}

const alglib_impl::ae_vector *ae_vector_wrapper::c_ptr() const { return This; }
alglib_impl::ae_vector *ae_vector_wrapper::c_ptr() { return This; }

#if !defined AE_NO_EXCEPTIONS
ae_vector_wrapper::ae_vector_wrapper(const char *s, alglib_impl::ae_datatype datatype) {
   std::vector<const char *> svec;
   size_t i;
   char *p = filter_spaces(s); if (p == NULL) ThrowError("ALGLIB: allocation error");
   try {
      str_vector_create(p, true, &svec);
   {
      alglib_impl::ae_state_init();
      TryCatch()
      is_frozen_proxy = false, This = &Obj, memset(This, 0, sizeof *This), ae_vector_init(This, (ae_int_t)(svec.size()), datatype, false);
      alglib_impl::ae_state_clear();
   }
      for (i = 0; i < svec.size(); i++) {
         if (datatype == alglib_impl::DT_BOOL) This->ptr.p_bool[i] = parse_bool_delim(svec[i], ",]");
         if (datatype == alglib_impl::DT_INT) This->ptr.p_int[i] = parse_int_delim(svec[i], ",]");
         if (datatype == alglib_impl::DT_REAL) This->ptr.p_double[i] = parse_real_delim(svec[i], ",]");
         if (datatype == alglib_impl::DT_COMPLEX) This->ptr.p_complex[i] = ae_complex_from_c(parse_complex_delim(svec[i], ",]"));
      }
      alglib_impl::ae_free(p);
   } catch( ...) {
      alglib_impl::ae_free(p);
      throw;
   }
}
#endif

boolean_1d_array::boolean_1d_array():ae_vector_wrapper(alglib_impl::DT_BOOL) { }
boolean_1d_array::boolean_1d_array(const boolean_1d_array &rhs):ae_vector_wrapper(rhs, alglib_impl::DT_BOOL) { }
boolean_1d_array::boolean_1d_array(alglib_impl::ae_vector *p):ae_vector_wrapper(p, alglib_impl::DT_BOOL) { }
const boolean_1d_array &boolean_1d_array::operator=(const boolean_1d_array &rhs) { return static_cast < const boolean_1d_array &>(assign(rhs)); }
boolean_1d_array::~boolean_1d_array() { }
const bool &boolean_1d_array::operator()(ae_int_t i) const { return This->ptr.p_bool[i]; }
bool &boolean_1d_array::operator()(ae_int_t i) { return This->ptr.p_bool[i]; }
const bool &boolean_1d_array::operator[](ae_int_t i) const { return This->ptr.p_bool[i]; }
bool &boolean_1d_array::operator[](ae_int_t i) { return This->ptr.p_bool[i]; }

void boolean_1d_array::setcontent(ae_int_t iLen, const bool *pContent) {
   ae_int_t i;
// setlength, with exception-free error handling fallback code
   setlength(iLen);
   if (This == NULL || This->cnt != iLen) return;
// copy
   for (i = 0; i < iLen; i++) This->ptr.p_bool[i] = pContent[i];
}

bool *boolean_1d_array::getcontent() { return This->ptr.p_bool; }
const bool *boolean_1d_array::getcontent() const { return This->ptr.p_bool; }

#if !defined AE_NO_EXCEPTIONS
boolean_1d_array::boolean_1d_array(const char *s):ae_vector_wrapper(s, alglib_impl::DT_BOOL) { }
std::string boolean_1d_array::tostring() const { return length() == 0? "[]": arraytostring(&(operator()(0)), length()); }
#endif

integer_1d_array::integer_1d_array():ae_vector_wrapper(alglib_impl::DT_INT) { }
integer_1d_array::integer_1d_array(alglib_impl::ae_vector *p):ae_vector_wrapper(p, alglib_impl::DT_INT) { }
integer_1d_array::integer_1d_array(const integer_1d_array &rhs):ae_vector_wrapper(rhs, alglib_impl::DT_INT) { }
const integer_1d_array &integer_1d_array::operator=(const integer_1d_array &rhs) { return static_cast < const integer_1d_array &>(assign(rhs)); }
integer_1d_array::~integer_1d_array() { }
const ae_int_t &integer_1d_array::operator()(ae_int_t i) const { return This->ptr.p_int[i]; }
ae_int_t &integer_1d_array::operator()(ae_int_t i) { return This->ptr.p_int[i]; }
const ae_int_t &integer_1d_array::operator[](ae_int_t i) const { return This->ptr.p_int[i]; }
ae_int_t &integer_1d_array::operator[](ae_int_t i) { return This->ptr.p_int[i]; }

void integer_1d_array::setcontent(ae_int_t iLen, const ae_int_t *pContent) {
   ae_int_t i;
// setlength(), handle possible exception-free errors
   setlength(iLen);
   if (This == NULL || This->cnt != iLen) return;
// copy
   for (i = 0; i < iLen; i++) This->ptr.p_int[i] = pContent[i];
}

ae_int_t *integer_1d_array::getcontent() { return This->ptr.p_int; }
const ae_int_t *integer_1d_array::getcontent() const { return This->ptr.p_int; }

#if !defined AE_NO_EXCEPTIONS
integer_1d_array::integer_1d_array(const char *s):ae_vector_wrapper(s, alglib_impl::DT_INT) { }
std::string integer_1d_array::tostring() const { return length() == 0? "[]": arraytostring(&operator()(0), length()); }
#endif

real_1d_array::real_1d_array():ae_vector_wrapper(alglib_impl::DT_REAL) { }
real_1d_array::real_1d_array(alglib_impl::ae_vector *p):ae_vector_wrapper(p, alglib_impl::DT_REAL) { }
real_1d_array::real_1d_array(const real_1d_array &rhs):ae_vector_wrapper(rhs, alglib_impl::DT_REAL) { }
const real_1d_array &real_1d_array::operator=(const real_1d_array &rhs) { return static_cast < const real_1d_array &>(assign(rhs)); }
real_1d_array::~real_1d_array() { }
const double &real_1d_array::operator()(ae_int_t i) const { return This->ptr.p_double[i]; }
double &real_1d_array::operator()(ae_int_t i) { return This->ptr.p_double[i]; }
const double &real_1d_array::operator[](ae_int_t i) const { return This->ptr.p_double[i]; }
double &real_1d_array::operator[](ae_int_t i) { return This->ptr.p_double[i]; }

void real_1d_array::setcontent(ae_int_t iLen, const double *pContent) {
   ae_int_t i;
// setlength(), handle possible exception-free errors
   setlength(iLen);
   if (This == NULL || This->cnt != iLen) return;
// copy
   for (i = 0; i < iLen; i++) This->ptr.p_double[i] = pContent[i];
}

// TODO: convert to constructor!!!!!!!
void real_1d_array::attach_to_ptr(ae_int_t iLen, double *pContent) {
   alglib_impl::x_vector x;
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      is_frozen_proxy = false, This = NULL, set_error_msg(); return;
#endif
   }
   alglib_impl::ae_assert(!is_frozen_proxy, "ALGLIB: unable to attach proxy object to something else");
   alglib_impl::ae_assert(iLen > 0, "ALGLIB: non-positive length for attach_to_ptr()");
   x.cnt = iLen;
   x.datatype = alglib_impl::DT_REAL;
   x.owner = alglib_impl::OWN_CALLER;
   x.last_action = alglib_impl::ACT_UNCHANGED;
   x.x_ptr.p_ptr = pContent;
   attach_to(&x);
   alglib_impl::ae_state_clear();
}

double *real_1d_array::getcontent() { return This->ptr.p_double; }
const double *real_1d_array::getcontent() const { return This->ptr.p_double; }

#if !defined AE_NO_EXCEPTIONS
real_1d_array::real_1d_array(const char *s):ae_vector_wrapper(s, alglib_impl::DT_REAL) { }
std::string real_1d_array::tostring(int dps) const { return length() == 0? "[]": arraytostring(&operator()(0), length(), dps); }
#endif

complex_1d_array::complex_1d_array():ae_vector_wrapper(alglib_impl::DT_COMPLEX) { }
complex_1d_array::complex_1d_array(alglib_impl::ae_vector *p):ae_vector_wrapper(p, alglib_impl::DT_COMPLEX) { }
complex_1d_array::complex_1d_array(const complex_1d_array &rhs):ae_vector_wrapper(rhs, alglib_impl::DT_COMPLEX) { }
const complex_1d_array &complex_1d_array::operator=(const complex_1d_array &rhs) { return static_cast < const complex_1d_array &>(assign(rhs)); }
complex_1d_array::~complex_1d_array() { }
const complex &complex_1d_array::operator()(ae_int_t i) const { return *(const complex *)(This->ptr.p_complex + i); }
complex &complex_1d_array::operator()(ae_int_t i) { return *(complex *)(This->ptr.p_complex + i); }
const complex &complex_1d_array::operator[](ae_int_t i) const { return *(const complex *)(This->ptr.p_complex + i); }
complex &complex_1d_array::operator[](ae_int_t i) { return *(complex *)(This->ptr.p_complex + i); }

void complex_1d_array::setcontent(ae_int_t iLen, const complex *pContent) {
   ae_int_t i;
// setlength(), handle possible exception-free errors
   setlength(iLen);
   if (This == NULL || This->cnt != iLen) return;
// copy
   for (i = 0; i < iLen; i++) This->ptr.p_complex[i] = ae_complex_from_c(pContent[i]);
}

complex *complex_1d_array::getcontent() { return (complex *)This->ptr.p_complex; }
const complex *complex_1d_array::getcontent() const { return (const complex *)This->ptr.p_complex; }

#if !defined AE_NO_EXCEPTIONS
complex_1d_array::complex_1d_array(const char *s):ae_vector_wrapper(s, alglib_impl::DT_COMPLEX) { }
std::string complex_1d_array::tostring(int dps) const { return length() == 0? "[]": arraytostring(&operator()(0), length(), dps); }
#endif

ae_matrix_wrapper::ae_matrix_wrapper(alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      is_frozen_proxy = false, This = NULL, set_error_msg(); return;
#endif
   }
   is_frozen_proxy = false, This = &Obj, memset(This, 0, sizeof *This), ae_matrix_init(This, 0, 0, datatype, false);
   alglib_impl::ae_state_clear();
}

ae_matrix_wrapper::ae_matrix_wrapper(alglib_impl::ae_matrix *e_ptr, alglib_impl::ae_datatype datatype) {
   if (e_ptr->datatype != datatype) {
      const char *msg = "ALGLIB: ae_vector_wrapper datatype check failed";
#if !defined AE_NO_EXCEPTIONS
      ThrowError(msg);
#else
      is_frozen_proxy = false, This = NULL, set_error_flag(msg); return;
#endif
   }
   is_frozen_proxy = true, This = e_ptr;
}

ae_matrix_wrapper::ae_matrix_wrapper(const ae_matrix_wrapper &rhs, alglib_impl::ae_datatype datatype) {
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      is_frozen_proxy = false, This = NULL, set_error_msg(); return;
#endif
   }
   is_frozen_proxy = false, This = NULL;
   if (rhs.This != NULL) {
      alglib_impl::ae_assert(rhs.This->datatype == datatype, "ALGLIB: ae_matrix_wrapper datatype check failed");
      This = &Obj, memset(This, 0, sizeof *This), ae_matrix_copy(This, rhs.This, false);
   }
   alglib_impl::ae_state_clear();
}

ae_matrix_wrapper::~ae_matrix_wrapper() { if (This == &Obj) ae_matrix_free(This, false); }

// TODO: automatic allocation of NULL ptr!!!!!
void ae_matrix_wrapper::setlength(ae_int_t rows, ae_int_t cols) {
   alglib_impl::ae_state_init();
   TryCatch()
   alglib_impl::ae_assert(This != NULL, "ALGLIB: setlength() error, p_mat == NULL (array was not correctly initialized)");
   alglib_impl::ae_assert(!is_frozen_proxy, "ALGLIB: setlength() error, attempt to resize proxy array");
   alglib_impl::ae_matrix_set_length(This, rows, cols);
   alglib_impl::ae_state_clear();
}

ae_int_t ae_matrix_wrapper::rows() const { return This == NULL? 0: This->rows; }
ae_int_t ae_matrix_wrapper::cols() const { return This == NULL? 0: This->cols; }
bool ae_matrix_wrapper::isempty() const { return rows() == 0 || cols() == 0; }
ae_int_t ae_matrix_wrapper::getstride() const { return This == NULL? 0: This->stride; }

void ae_matrix_wrapper::attach_to(alglib_impl::x_matrix *new_ptr) {
   if (This == &Obj) ae_matrix_free(This, false);
   is_frozen_proxy = true, This = &Obj, memset(This, 0, sizeof *This), ae_matrix_init_attach_to_x(This, new_ptr, false);
}

const ae_matrix_wrapper &ae_matrix_wrapper::assign(const ae_matrix_wrapper &rhs) {
   ae_int_t i;
   if (this == &rhs) return *this;
   alglib_impl::ae_state_init();
   TryCatch(*this)
   alglib_impl::ae_assert(This != NULL, "ALGLIB: incorrect assignment to matrix (uninitialized destination)");
   alglib_impl::ae_assert(rhs.This != NULL, "ALGLIB: incorrect assignment to array (uninitialized source)");
   alglib_impl::ae_assert(rhs.This->datatype == This->datatype, "ALGLIB: incorrect assignment to array (types dont match)");
   if (is_frozen_proxy) {
      alglib_impl::ae_assert(rhs.This->rows == This->rows, "ALGLIB: incorrect assignment to proxy array (sizes dont match)");
      alglib_impl::ae_assert(rhs.This->cols == This->cols, "ALGLIB: incorrect assignment to proxy array (sizes dont match)");
   }
   if (rhs.This->rows != This->rows || rhs.This->cols != This->cols) ae_matrix_set_length(This, rhs.This->rows, rhs.This->cols);
   for (i = 0; i < This->rows; i++) memcpy(This->ptr.pp_void[i], rhs.This->ptr.pp_void[i], This->cols * alglib_impl::ae_sizeof(This->datatype));
   alglib_impl::ae_state_clear();
   return *this;
}

#if !defined AE_NO_EXCEPTIONS
ae_matrix_wrapper::ae_matrix_wrapper(const char *s, alglib_impl::ae_datatype datatype) {
   std::vector< std::vector<const char *> > smat;
   size_t i, j;
   char *p = filter_spaces(s); if (p == NULL) ThrowError("ALGLIB: allocation error");
   try {
      str_matrix_create(p, &smat);
   {
      alglib_impl::ae_state_init();
      TryCatch()
      is_frozen_proxy = false, This = &Obj, memset(This, 0, sizeof *This);
      if (smat.size() != 0)
         ae_matrix_init(This, (ae_int_t)(smat.size()), (ae_int_t)(smat[0].size()), datatype, false);
      else
         ae_matrix_init(This, 0, 0, datatype, false);
      alglib_impl::ae_state_clear();
   }
      for (i = 0; i < smat.size(); i++) for (j = 0; j < smat[0].size(); j++) {
         if (datatype == alglib_impl::DT_BOOL) This->ptr.pp_bool[i][j] = parse_bool_delim(smat[i][j], ",]");
         if (datatype == alglib_impl::DT_INT) This->ptr.pp_int[i][j] = parse_int_delim(smat[i][j], ",]");
         if (datatype == alglib_impl::DT_REAL) This->ptr.pp_double[i][j] = parse_real_delim(smat[i][j], ",]");
         if (datatype == alglib_impl::DT_COMPLEX) This->ptr.pp_complex[i][j] = ae_complex_from_c(parse_complex_delim(smat[i][j], ",]"));
      }
      alglib_impl::ae_free(p);
   } catch( ...) {
      alglib_impl::ae_free(p);
      throw;
   }
}
#endif

const alglib_impl::ae_matrix *ae_matrix_wrapper::c_ptr() const { return This; }
alglib_impl::ae_matrix *ae_matrix_wrapper::c_ptr() { return This; }
boolean_2d_array::boolean_2d_array():ae_matrix_wrapper(alglib_impl::DT_BOOL) { }
boolean_2d_array::boolean_2d_array(const boolean_2d_array &rhs):ae_matrix_wrapper(rhs, alglib_impl::DT_BOOL) { }
boolean_2d_array::boolean_2d_array(alglib_impl::ae_matrix *p):ae_matrix_wrapper(p, alglib_impl::DT_BOOL) { }
boolean_2d_array::~boolean_2d_array() { }
const boolean_2d_array &boolean_2d_array::operator=(const boolean_2d_array &rhs) { return static_cast < const boolean_2d_array &>(assign(rhs)); }
const bool &boolean_2d_array::operator()(ae_int_t i, ae_int_t j) const { return This->ptr.pp_bool[i][j]; }
bool &boolean_2d_array::operator()(ae_int_t i, ae_int_t j) { return This->ptr.pp_bool[i][j]; }
const bool *boolean_2d_array::operator[](ae_int_t i) const { return This->ptr.pp_bool[i]; }
bool *boolean_2d_array::operator[](ae_int_t i) { return This->ptr.pp_bool[i]; }

void boolean_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const bool *pContent) {
   ae_int_t i, j;
// setlength(), handle possible exception-free errors
   setlength(irows, icols);
   if (This == NULL || This->rows != irows || This->cols != icols) return;
// copy
   for (i = 0; i < irows; i++) for (j = 0; j < icols; j++) This->ptr.pp_bool[i][j] = pContent[i * icols + j];
}

#if !defined AE_NO_EXCEPTIONS
boolean_2d_array::boolean_2d_array(const char *s):ae_matrix_wrapper(s, alglib_impl::DT_BOOL) { }

std::string boolean_2d_array::tostring() const {
   std::string result;
   ae_int_t i;
   if (isempty()) return "[[]]";
   result = "[";
   for (i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols());
   }
   result += "]";
   return result;
}
#endif

integer_2d_array::integer_2d_array():ae_matrix_wrapper(alglib_impl::DT_INT) { }
integer_2d_array::integer_2d_array(const integer_2d_array &rhs):ae_matrix_wrapper(rhs, alglib_impl::DT_INT) { }
integer_2d_array::integer_2d_array(alglib_impl::ae_matrix *p):ae_matrix_wrapper(p, alglib_impl::DT_INT) { }
integer_2d_array::~integer_2d_array() { }
const integer_2d_array &integer_2d_array::operator=(const integer_2d_array &rhs) { return static_cast < const integer_2d_array &>(assign(rhs)); }
const ae_int_t &integer_2d_array::operator()(ae_int_t i, ae_int_t j) const { return This->ptr.pp_int[i][j]; }
ae_int_t &integer_2d_array::operator()(ae_int_t i, ae_int_t j) { return This->ptr.pp_int[i][j]; }
const ae_int_t *integer_2d_array::operator[](ae_int_t i) const { return This->ptr.pp_int[i]; }
ae_int_t *integer_2d_array::operator[](ae_int_t i) { return This->ptr.pp_int[i]; }

void integer_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const ae_int_t *pContent) {
   ae_int_t i, j;
// setlength(), handle possible exception-free errors
   setlength(irows, icols);
   if (This == NULL || This->rows != irows || This->cols != icols) return;
// copy
   for (i = 0; i < irows; i++) for (j = 0; j < icols; j++) This->ptr.pp_int[i][j] = pContent[i * icols + j];
}

#if !defined AE_NO_EXCEPTIONS
integer_2d_array::integer_2d_array(const char *s):ae_matrix_wrapper(s, alglib_impl::DT_INT) { }

std::string integer_2d_array::tostring() const {
   std::string result;
   ae_int_t i;
   if (isempty()) return "[[]]";
   result = "[";
   for (i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols());
   }
   result += "]";
   return result;
}
#endif

real_2d_array::real_2d_array():ae_matrix_wrapper(alglib_impl::DT_REAL) { }
real_2d_array::real_2d_array(const real_2d_array &rhs):ae_matrix_wrapper(rhs, alglib_impl::DT_REAL) { }
real_2d_array::real_2d_array(alglib_impl::ae_matrix *p):ae_matrix_wrapper(p, alglib_impl::DT_REAL) { }
real_2d_array::~real_2d_array() { }
const real_2d_array &real_2d_array::operator=(const real_2d_array &rhs) { return static_cast < const real_2d_array &>(assign(rhs)); }
const double &real_2d_array::operator()(ae_int_t i, ae_int_t j) const { return This->ptr.pp_double[i][j]; }
double &real_2d_array::operator()(ae_int_t i, ae_int_t j) { return This->ptr.pp_double[i][j]; }
const double *real_2d_array::operator[](ae_int_t i) const { return This->ptr.pp_double[i]; }
double *real_2d_array::operator[](ae_int_t i) { return This->ptr.pp_double[i]; }

void real_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const double *pContent) {
   ae_int_t i, j;
// setlength(), handle possible exception-free errors
   setlength(irows, icols);
   if (This == NULL || This->rows != irows || This->cols != icols) return;
// copy
   for (i = 0; i < irows; i++) for (j = 0; j < icols; j++) This->ptr.pp_double[i][j] = pContent[i * icols + j];
}

void real_2d_array::attach_to_ptr(ae_int_t irows, ae_int_t icols, double *pContent) {
   alglib_impl::x_matrix x;
   alglib_impl::ae_state_init();
   TryX {
#if !defined AE_NO_EXCEPTIONS
      ThrowErrorMsg();
#else
      is_frozen_proxy = false, This = NULL, set_error_msg(); return;
#endif
   }
   alglib_impl::ae_assert(!is_frozen_proxy, "ALGLIB: unable to attach proxy object to something else");
   alglib_impl::ae_assert(irows > 0 && icols > 0, "ALGLIB: non-positive length for attach_to_ptr()");
   x.rows = irows;
   x.cols = icols;
   x.stride = icols;
   x.datatype = alglib_impl::DT_REAL;
   x.owner = alglib_impl::OWN_CALLER;
   x.last_action = alglib_impl::ACT_UNCHANGED;
   x.x_ptr.p_ptr = pContent;
   attach_to(&x);
   alglib_impl::ae_state_clear();
}

#if !defined AE_NO_EXCEPTIONS
real_2d_array::real_2d_array(const char *s):ae_matrix_wrapper(s, alglib_impl::DT_REAL) { }

std::string real_2d_array::tostring(int dps) const {
   std::string result;
   ae_int_t i;
   if (isempty()) return "[[]]";
   result = "[";
   for (i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols(), dps);
   }
   result += "]";
   return result;
}
#endif

complex_2d_array::complex_2d_array():ae_matrix_wrapper(alglib_impl::DT_COMPLEX) { }
complex_2d_array::complex_2d_array(const complex_2d_array &rhs):ae_matrix_wrapper(rhs, alglib_impl::DT_COMPLEX) { }
complex_2d_array::complex_2d_array(alglib_impl::ae_matrix *p):ae_matrix_wrapper(p, alglib_impl::DT_COMPLEX) { }
complex_2d_array::~complex_2d_array() { }
const complex_2d_array &complex_2d_array::operator=(const complex_2d_array &rhs) { return static_cast < const complex_2d_array &>(assign(rhs)); }
const complex &complex_2d_array::operator()(ae_int_t i, ae_int_t j) const { return *(const complex *)(This->ptr.pp_complex[i] + j); }
complex &complex_2d_array::operator()(ae_int_t i, ae_int_t j) { return *(complex *)(This->ptr.pp_complex[i] + j); }
const complex *complex_2d_array::operator[](ae_int_t i) const { return (const complex *)(This->ptr.pp_complex[i]); }
complex *complex_2d_array::operator[](ae_int_t i) { return (complex *)(This->ptr.pp_complex[i]); }

void complex_2d_array::setcontent(ae_int_t irows, ae_int_t icols, const complex *pContent) {
   ae_int_t i, j;
// setlength(), handle possible exception-free errors
   setlength(irows, icols);
   if (This == NULL || This->rows != irows || This->cols != icols) return;
// copy
   for (i = 0; i < irows; i++) for (j = 0; j < icols; j++) {
      This->ptr.pp_complex[i][j] = ae_complex_from_c(pContent[i * icols + j]);
   }
}

#if !defined AE_NO_EXCEPTIONS
complex_2d_array::complex_2d_array(const char *s):ae_matrix_wrapper(s, alglib_impl::DT_COMPLEX) { }

std::string complex_2d_array::tostring(int dps) const {
   std::string result;
   ae_int_t i;
   if (isempty()) return "[[]]";
   result = "[";
   for (i = 0; i < rows(); i++) {
      if (i != 0) result += ",";
      result += arraytostring(&operator()(i, 0), cols(), dps);
   }
   result += "]";
   return result;
}
#endif

// Internal functions
ae_int_t my_stricmp(const char *s1, const char *s2) {
   int c1, c2;
// handle special cases
   if (s1 == NULL && s2 != NULL) return -1;
   if (s1 != NULL && s2 == NULL) return +1;
   if (s1 == NULL && s2 == NULL) return 0;
// compare
   while (true) {
      c1 = *s1;
      c2 = *s2;
      s1++;
      s2++;
      if (c1 == 0) return c2 == 0 ? 0 : -1;
      if (c2 == 0) return c1 == 0 ? 0 : +1;
      c1 = tolower(c1);
      c2 = tolower(c2);
      if (c1 < c2) return -1;
      if (c1 > c2) return +1;
   }
}

#if !defined AE_NO_EXCEPTIONS
//
// This function filters out all spaces from the string.
// It returns string allocated with ae_malloc().
// On allocaction failure returns NULL.
//
char *filter_spaces(const char *s) {
   size_t i, n;
   char *r;
   char *r0;
   n = strlen(s);
   r = (char *)alglib_impl::ae_malloc(n + 1);
   if (r == NULL) return r;
   for (i = 0, r0 = r; i <= n; i++, s++)
      if (!isspace(*s)) {
         *r0 = *s;
         r0++;
      }
   return r;
}

void str_vector_create(const char *src, bool match_head_only, std::vector<const char *> *p_vec) {
// parse beginning of the string.
// try to handle "[]" string
   p_vec->clear();
   if (*src != '[')
      ThrowError("Incorrect initializer for vector");
   src++;
   if (*src == ']')
      return;
   p_vec->push_back(src);
   while (true) {
      if (*src == 0)
         ThrowError("Incorrect initializer for vector");
      if (*src == ']') {
         if (src[1] == 0 || !match_head_only)
            return;
         ThrowError("Incorrect initializer for vector");
      }
      if (*src == ',') {
         p_vec->push_back(src + 1);
         src++;
         continue;
      }
      src++;
   }
}

void str_matrix_create(const char *src, std::vector< std::vector<const char *> > *p_mat) {
   p_mat->clear();
// Try to handle "[[]]" string
   if (strcmp(src, "[[]]") == 0) return;
// Parse non-empty string
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
      src++;
      if (*src == ',') {
         src++;
         continue;
      }
      if (*src == ']')
         break;
      ThrowError("Incorrect initializer for matrix");
   }
   src++;
   if (*src != 0)
      ThrowError("Incorrect initializer for matrix");
}

bool parse_bool_delim(const char *s, const char *delim) {
   const char *p;
   char buf[8];
// try to parse false
   p = "false";
   memset(buf, 0, sizeof(buf));
   strncpy(buf, s, strlen(p));
   if (my_stricmp(buf, p) == 0) {
      if (s[strlen(p)] == 0 || strchr(delim, s[strlen(p)]) == NULL)
         ThrowError("Cannot parse value");
      return false;
   }
// try to parse true
   p = "true";
   memset(buf, 0, sizeof(buf));
   strncpy(buf, s, strlen(p));
   if (my_stricmp(buf, p) == 0) {
      if (s[strlen(p)] == 0 || strchr(delim, s[strlen(p)]) == NULL)
         ThrowError("Cannot parse value");
      return true;
   }
// error
   ThrowError("Cannot parse value");
}

ae_int_t parse_int_delim(const char *s, const char *delim) {
   const char *p;
   long long_val;
   volatile ae_int_t ae_val;
   p = s;
//
// check string structure:
// * leading sign
// * at least one digit
// * delimiter
//
   if (*s == '-' || *s == '+')
      s++;
   if (*s == 0 || strchr("1234567890", *s) == NULL)
      ThrowError("Cannot parse value");
   while (*s != 0 && strchr("1234567890", *s) != NULL)
      s++;
   if (*s == 0 || strchr(delim, *s) == NULL)
      ThrowError("Cannot parse value");
// convert and ensure that value fits into ae_int_t
   s = p;
   long_val = atol(s);
   ae_val = long_val;
   if (ae_val != long_val)
      ThrowError("Cannot parse value");
   return ae_val;
}

bool _parse_real_delim(const char *s, const char *delim, double *result, const char **new_s) {
   const char *p;
   char *t;
   bool has_digits;
   char buf[64];
   int isign;
   lconv *loc;
   p = s;
//
// check string structure and decide what to do
//
   isign = 1;
   if (*s == '-' || *s == '+') {
      isign = *s == '-' ? -1 : +1;
      s++;
   }
   memset(buf, 0, sizeof(buf));
   strncpy(buf, s, 3);
   if (my_stricmp(buf, "nan") != 0 && my_stricmp(buf, "inf") != 0) {
   //
   // [sign] [ddd] [.] [ddd] [e|E[sign]ddd]
   //
      has_digits = false;
      if (*s != 0 && strchr("1234567890", *s) != NULL) {
         has_digits = true;
         while (*s != 0 && strchr("1234567890", *s) != NULL)
            s++;
      }
      if (*s == '.')
         s++;
      if (*s != 0 && strchr("1234567890", *s) != NULL) {
         has_digits = true;
         while (*s != 0 && strchr("1234567890", *s) != NULL)
            s++;
      }
      if (!has_digits)
         return false;
      if (*s == 'e' || *s == 'E') {
         s++;
         if (*s == '-' || *s == '+')
            s++;
         if (*s == 0 || strchr("1234567890", *s) == NULL)
            return false;
         while (*s != 0 && strchr("1234567890", *s) != NULL)
            s++;
      }
      if (*s == 0 || strchr(delim, *s) == NULL)
         return false;
      *new_s = s;
   //
   // finite value conversion
   //
      if (*new_s - p >= (int)sizeof(buf))
         return false;
      strncpy(buf, p, (size_t)(*new_s - p));
      buf[*new_s - p] = 0;
      loc = localeconv();
      t = strchr(buf, '.');
      if (t != NULL)
         *t = *loc->decimal_point;
      *result = atof(buf);
      return true;
   } else {
   //
   // check delimiter and update *new_s
   //
      s += 3;
      if (*s == 0 || strchr(delim, *s) == NULL)
         return false;
      *new_s = s;
   //
   // NAN, INF conversion
   //
      if (my_stricmp(buf, "nan") == 0)
         *result = NAN;
      if (my_stricmp(buf, "inf") == 0)
         *result = isign > 0 ? +INFINITY : -INFINITY;
      return true;
   }
}

double parse_real_delim(const char *s, const char *delim) {
   double result;
   const char *new_s;
   if (!_parse_real_delim(s, delim, &result, &new_s))
      ThrowError("Cannot parse value");
   return result;
}

complex parse_complex_delim(const char *s, const char *delim) {
   double d_result;
   const char *new_s;
   complex c_result;
// parse as real value
   if (_parse_real_delim(s, delim, &d_result, &new_s))
      return d_result;
// parse as "a+bi" or "a-bi"
   if (_parse_real_delim(s, "+-", &c_result.x, &new_s)) {
      s = new_s;
      if (!_parse_real_delim(s, "i", &c_result.y, &new_s))
         ThrowError("Cannot parse value");
      s = new_s + 1;
      if (*s == 0 || strchr(delim, *s) == NULL)
         ThrowError("Cannot parse value");
      return c_result;
   }
// parse as complex value "bi+a" or "bi-a"
   if (_parse_real_delim(s, "i", &c_result.y, &new_s)) {
      s = new_s + 1;
      if (*s == 0)
         ThrowError("Cannot parse value");
      if (strchr(delim, *s) != NULL) {
         c_result.x = 0;
         return c_result;
      }
      if (strchr("+-", *s) != NULL) {
         if (!_parse_real_delim(s, delim, &c_result.x, &new_s))
            ThrowError("Cannot parse value");
         return c_result;
      }
      ThrowError("Cannot parse value");
   }
// error
   ThrowError("Cannot parse value");
}

std::string arraytostring(const bool *ptr, ae_int_t n) {
   std::string result;
   ae_int_t i;
   result = "[";
   for (i = 0; i < n; i++) {
      if (i != 0)
         result += ",";
      result += ptr[i] ? "true" : "false";
   }
   result += "]";
   return result;
}

std::string arraytostring(const ae_int_t *ptr, ae_int_t n) {
   std::string result;
   ae_int_t i;
   char buf[64];
   result = "[";
   for (i = 0; i < n; i++) {
      if (sprintf(buf, i == 0 ? "%ld" : ",%ld", long (ptr[i])) >= (int)sizeof(buf))
          ThrowError("arraytostring(): buffer overflow");
      result += buf;
   }
   result += "]";
   return result;
}

std::string arraytostring(const double *ptr, ae_int_t n, int _dps) {
   std::string result;
   ae_int_t i;
   char buf[64];
   char mask1[64];
   char mask2[64];
   int dps = _dps >= 0 ? _dps : -_dps;
   result = "[";
   if (sprintf(mask1, "%%.%d%s", dps, _dps >= 0 ? "f" : "e") >= (int)sizeof(mask1))
      ThrowError("arraytostring(): buffer overflow");
   if (sprintf(mask2, ",%s", mask1) >= (int)sizeof(mask2))
      ThrowError("arraytostring(): buffer overflow");
   for (i = 0; i < n; i++) {
      buf[0] = 0;
      if (isfinite(ptr[i])) {
         if (sprintf(buf, i == 0 ? mask1 : mask2, double (ptr[i])) >= (int)sizeof(buf))
             ThrowError("arraytostring(): buffer overflow");
      } else if (isnan(ptr[i]))
         strcpy(buf, i == 0 ? "NAN" : ",NAN");
      else if (isposinf(ptr[i]))
         strcpy(buf, i == 0 ? "+INF" : ",+INF");
      else if (isneginf(ptr[i]))
         strcpy(buf, i == 0 ? "-INF" : ",-INF");
      result += buf;
   }
   result += "]";
   return result;
}

std::string arraytostring(const complex *ptr, ae_int_t n, int dps) {
   std::string result;
   ae_int_t i;
   result = "[";
   for (i = 0; i < n; i++) {
      if (i != 0)
         result += ",";
      result += ptr[i].tostring(dps);
   }
   result += "]";
   return result;
}
#endif

// standard functions
int sign(double x) { return x > 0? +1: x < 0? -1: 0; }

double randomreal() {
   const double mx = (double)RAND_MAX + 1.0;
   return (rand() + rand() / mx) / mx;
}
double randommid() {
   const double mx = (double)RAND_MAX + 1.0;
   return 2.0*(rand() + rand() / mx) / mx - 1.0;
}
ae_int_t randominteger(ae_int_t maxv) { return ((ae_int_t) rand()) % maxv; }
bool randombool(double p) {
   const double mx = (double)RAND_MAX + 1.0;
   return rand() + rand()/mx <= p*mx;
}

int RoundZ(double x) { return int(round(x)); }
int TruncZ(double x) { return int(trunc(x)); }
int FloorZ(double x) { return int(floor(x)); }
int CeilZ(double x) { return int(ceil(x)); }
double pi() { return 3.14159265358979323846; }
double sqr(double x) { return x*x; }
int maxint(int m1, int m2) { return m1 > m2 ? m1 : m2; }
int minint(int m1, int m2) { return m1 > m2 ? m2 : m1; }
double maxreal(double m1, double m2) { return m1 > m2 ? m1 : m2; }
double minreal(double m1, double m2) { return m1 > m2 ? m2 : m1; }
bool isposinf(double A) { return isinf(A) && !signbit(A); }
bool isneginf(double A) { return isinf(A) && signbit(A); }

// CSV functions
#if !defined AE_NO_EXCEPTIONS
void read_csv(const char *filename, char separator, int flags, real_2d_array &out) {
   int flag;
//
// Parameters
//
   bool skip_first_row = (flags & CSV_SKIP_HEADERS) != 0;
//
// Prepare empty output array
//
   out.setlength(0, 0);
//
// Open file, determine size, read contents
//
   FILE *f_in = fopen(filename, "rb");
   if (f_in == NULL)
      ThrowError("read_csv: unable to open input file");
   flag = fseek(f_in, 0, SEEK_END);
   AE_CRITICAL_ASSERT(flag == 0);
   long int _filesize = ftell(f_in);
   AE_CRITICAL_ASSERT(_filesize >= 0);
   if (_filesize == 0) {
   // empty file, return empty array, success
      fclose(f_in);
      return;
   }
   size_t filesize = _filesize;
   std::vector<char> v_buf;
   v_buf.resize(filesize + 2, 0);
   char *p_buf = &v_buf[0]; //(#) This is NOT equivalent to char *p_buf = v_buf!
   flag = fseek(f_in, 0, SEEK_SET);
   AE_CRITICAL_ASSERT(flag == 0);
   size_t bytes_read = fread((void *)p_buf, 1, filesize, f_in);
   AE_CRITICAL_ASSERT(bytes_read == filesize);
   fclose(f_in);
//
// Normalize file contents:
// * replace 0x0 by spaces
// * remove trailing spaces and newlines
// * append trailing '\n' and '\0' characters
// Return if file contains only spaces/newlines.
//
   for (size_t i = 0; i < filesize; i++)
      if (p_buf[i] == 0)
         p_buf[i] = ' ';
   for (; filesize > 0;) {
      char c = p_buf[filesize - 1];
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
         filesize--;
         continue;
      }
      break;
   }
   if (filesize == 0)
      return;
   p_buf[filesize] = '\n';
   p_buf[filesize + 1] = '\0';
   filesize += 2;
//
// Scan dataset.
//
   size_t rows_count = 0, cols_count = 0, max_length = 0;
   std::vector<size_t> offsets, lengths;
   for (size_t row_start = 0; p_buf[row_start] != 0x0;) {
   // determine row length
      size_t row_length;
      for (row_length = 0; p_buf[row_start + row_length] != '\n'; row_length++);
   // determine cols count, perform integrity check
      size_t cur_cols_cnt = 1;
      for (size_t idx = 0; idx < row_length; idx++)
         if (p_buf[row_start + idx] == separator)
            cur_cols_cnt++;
      if (cols_count > 0 && cols_count != cur_cols_cnt)
         ThrowError("read_csv: non-rectangular contents, rows have different sizes");
      cols_count = cur_cols_cnt;
   // store offsets and lengths of the fields
      size_t cur_offs = 0;
      for (size_t idx = 0; idx < row_length + 1; idx++)
         if (p_buf[row_start + idx] == separator || p_buf[row_start + idx] == '\n') {
            offsets.push_back(row_start + cur_offs);
            lengths.push_back(idx - cur_offs);
            max_length = idx - cur_offs > max_length ? idx - cur_offs : max_length;
            cur_offs = idx + 1;
         }
   // advance row start
      rows_count++;
      row_start += row_length + 1;
   }
   AE_CRITICAL_ASSERT(rows_count >= 1);
   AE_CRITICAL_ASSERT(cols_count >= 1);
   AE_CRITICAL_ASSERT(cols_count * rows_count == offsets.size());
   AE_CRITICAL_ASSERT(cols_count * rows_count == lengths.size());
   if (rows_count == 1 && skip_first_row)       // empty output, return
      return;
//
// Convert
//
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
