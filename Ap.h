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
#ifndef OnceOnlyAp_h
#define OnceOnlyAp_h

#if defined InAlgLib && defined _MSC_VER
#   define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string>
#include <cstring>
#include <iostream>
#include <math.h>
#if defined __CODEGEARC__
#   include <list>
#   include <vector>
#elif defined __BORLANDC__
#   include <list.h>
#   include <vector.h>
#else
#   include <list>
#   include <vector>
#endif

// The CPU type.
#define AE_OTHER_CPU	0
#define AE_INTEL	1
#define AE_SPARC	2

// The OS type.
// Use AE_OTHER_OS if the OS has not yet been defined.
#define AE_OTHER_OS	0
#define AE_WINDOWS	1
#define AE_POSIX	2
#define AE_LINUX	304
#if !defined AE_OS
#   define AE_OS AE_OTHER_OS
#elif AE_OS == AE_LINUX
#   undef AE_OS
#   define AE_OS AE_POSIX
#   define _ALGLIB_USE_LINUX_EXTENSIONS
#endif

// The threading models for AE_THREADING / state flags.
//	NonTH: Serial Unsafe / Use Global.
//	SerTH: Serial.
//	ParTH: Parallel.
typedef enum { NonTH, SerTH, ParTH } xparams;
#if !defined AE_THREADING
#   define AE_THREADING ParTH
#endif

// The memory allocation types for AE_MALLOC.
#define AE_STDLIB_MALLOC	200
#define AE_BASIC_STATIC_MALLOC	201
#if !defined AE_MALLOC
#   define AE_MALLOC AE_STDLIB_MALLOC
#endif

// The compiler type and compiler-specific definitions and includes.
#define AE_OTHERC	0 // Unknown to ALGLIB++.
#define AE_MSVC		1 // MSVC.
#define AE_GNUC		2 // GCC/CLANG/ICC.
#define AE_SUNC		3 // Sun studio.
#if defined __GNUC__
#   define AE_COMPILER AE_GNUC
#elif defined __SUNPRO_C || defined __SUNPRO_CC
#   define AE_COMPILER AE_SUNC
#elif defined _MSC_VER
#   define AE_COMPILER AE_MSVC
#   if defined InAlgLib && !defined AE_ALL_WARNINGS
// Disable some irrelevant warnings.
//(@) Originally preceded the headers other than Ap.h in TestC.cpp and followed them in the program module *.cpp files.
//(@) None of this has actually been tested on Windows yet.
#      pragma warning(disable:4100)
#      pragma warning(disable:4127)
#      pragma warning(disable:4611)
#      pragma warning(disable:4702)
#      pragma warning(disable:4996)
#   endif
#else
#   define AE_COMPILER AE_OTHERC
#endif

// Solaris Studio C/C++, IBM XL C/C++, GNU C, CLang, IntelC++ Compiler (Linux): _thread static
// Visual C++, Intel C/C++ (Windows), C++ Builder, Digital Mars C++: __declspec(thread) static
// C++ Builder: __thread static or _declspec(thread) static
// On Windows before Vista and Server 2008 _declspec(thread) won't work in DLL's except those bound to executables.
// In C++11 (and in C++ if <threads.h> is supported): thread_local static; requires #include <threads.h>
// In C11 (in general): Thread_local static
// The keywords "auto" and "Thread_local"/"thread_local" in C/C++ have mutually exclusive uses
// and can therefore be considered as different cases of the same concept!
#define AutoS static // "auto static": static thread-local variables in block or file scope.

// Now we are ready to include the remaining headers.
// #include <ctype.h>
// #include <stdarg.h>
#include <string.h>
#include <setjmp.h>
#if defined AE_HAVE_STDINT
#   include <stdint.h>
#endif

// Intel SIMD intrinsics
// The preprocessor directives below:
// -	include headers for SSE2/AVX2/AVX2+FMA3 intrinsics,
// -	define _ALGLIB_HAS_SSE2_INTRINSICS, _ALGLIB_HAS_AVX2_INTRINSICS and _ALGLIB_HAS_FMA_INTRINSICS definitions.
// These actions are performed when we have:
// -	an x86 architecture definition AE_CPU == AE_INTEL,
// -	a compiler which supports intrinsics.
// The presence of _ALGLIB_HAS_???_INTRINSICS does NOT mean that our CPU actually supports these intrinsics -
// such things should be determined by CurCPU, which is initialized on start-up.
// It means that we are working under Intel and our compiler can issue SIMD-capable code.
#if defined AE_CPU && AE_CPU == AE_INTEL // Intel definitions.
// Other than Sun studio, we only assume that the compiler supports all instruction sets if something is not explicitly turned off.
#   if AE_COMPILER == AE_GNUC && !defined AE_NO_SSE2 || AE_COMPILER == AE_SUNC
#      include <xmmintrin.h>
#   endif
#   if AE_COMPILER == AE_MSVC && !defined AE_NO_SSE2 || AE_COMPILER == AE_SUNC
#      include <emmintrin.h>
#   endif
#   if !defined AE_NO_SSE2 || AE_COMPILER == AE_SUNC
#      define AE_HAS_SSE2_INTRINSICS
#      define _ALGLIB_HAS_SSE2_INTRINSICS
#   endif
#   if (AE_COMPILER == AE_GNUC && !defined AE_NO_AVX2 || AE_COMPILER == AE_OTHERC) && !defined AE_NO_SSE2 || AE_COMPILER == AE_SUNC
#      include <immintrin.h> //(@) Originally preceded the #defines *HAS_SSE2_INTRINSICS for AE_OTHER_OS.
#   endif
#   if AE_COMPILER == AE_MSVC && !defined AE_NO_SSE2 && !defined AE_NO_AVX2
#      include <intrin.h>
#   endif
#   if !defined AE_NO_SSE2 && !defined AE_NO_AVX2 || AE_COMPILER == AE_SUNC
#      define _ALGLIB_HAS_AVX2_INTRINSICS
#   endif
#   if !defined AE_NO_SSE2 && !defined AE_NO_AVX2 && !defined AE_NO_FMA || AE_COMPILER == AE_SUNC
#      define _ALGLIB_HAS_FMA_INTRINSICS
#   endif
// Intel integrity checks.
#   if defined _ALGLIB_INTEGRITY_CHECKS_ONCE && defined _ALGLIB_FAIL_WITHOUT_FMA_INTRINSICS && !defined _ALGLIB_HAS_FMA_INTRINSICS
#      error ALGLIB was requested to fail without FMA intrinsics
#   endif
#endif

namespace alglib_impl {
// Core Code (Vectors, Matrices, Memory Management, etc.)

// We are working under a C++ environment: define several conditions.
#define AE_USE_CPP_SERIALIZATION
// #include <iostream> //(@) Already included.

// Define ae_int32_t, ae_int64_t, ae_uint64_t, ae_int_t, complex, ae_error_type and ae_datatype.
// A boolean type was also originally (and unnecessarily) defined.
// For C (as of 2011): one needs only to include <stdbool.h> to define "bool", "false" and "true".
// For C++: it is already a part of the language.
#if defined AE_INT32_T
typedef AE_INT32_T ae_int32_t;
#elif defined AE_HAVE_STDINT
typedef int32_t ae_int32_t;
#elif AE_COMPILER == AE_MSVC
typedef __int32 ae_int32_t;
#else
typedef int ae_int32_t;
#endif

#if defined AE_INT64_T
typedef AE_INT64_T ae_int64_t;
#elif defined AE_HAVE_STDINT
typedef int64_t ae_int64_t;
#elif AE_COMPILER == AE_MSVC
typedef __int64 ae_int64_t;
#else
typedef signed long long ae_int64_t;
#endif

#if defined AE_UINT64_T
typedef AE_UINT64_T ae_uint64_t;
#elif defined AE_HAVE_STDINT
typedef uint64_t ae_uint64_t;
#elif AE_COMPILER == AE_MSVC
typedef unsigned __int64 ae_uint64_t;
#else
typedef unsigned long long ae_uint64_t;
#endif

#if !defined AE_INT_T
typedef ptrdiff_t ae_int_t;
#endif

typedef double Real; // Used in the Kernel*.{cpp,h} files.

struct complex { double x, y; };

typedef enum { ERR_OK = 0, ERR_OUT_OF_MEMORY = 1, ERR_XARRAY_TOO_LARGE = 2, ERR_ASSERTION_FAILED = 3 } ae_error_type;
typedef enum { DT_BOOL = 1, DT_BYTE = 1, DT_INT = 2, DT_REAL = 3, DT_COMPLEX = 4 } ae_datatype;

// Allocation-tracking: inactive by default.
// Turned on when needed for debugging purposes.
// _alloc_counter is incremented by 1 on malloc(), decremented on free().
// _alloc_counter_total is only incremented by 1.
extern ae_int_t _alloc_counter;
extern ae_int_t _alloc_counter_total;
extern bool _use_alloc_counter;

// Malloc debugging:
// Set to force ALGLIB++ malloc() to fail.
// Useful to debug errors-handling during memory allocation.
extern bool _force_malloc_failure;
// If 0 and if _use_alloc_counter is set: the upper limit of _alloc_counter before forcing ALGLIB++ malloc() to fail.
// Otherwise, this value has no effect.
extern ae_int_t _malloc_failure_after;

// Frame marker.
typedef struct ae_dyn_block ae_frame;

void ae_state_set_break_jump(jmp_buf *buf);
void ae_state_set_flags(ae_uint64_t flags);

void ae_frame_make(ae_frame *Fr);
void ae_frame_leave();

void ae_state_init();
void ae_state_clear();

void ae_clean_up();
void ae_assert(bool cond, const char *msg);

// Get/set the threading model type.
ae_uint64_t ae_get_global_threading();
void ae_set_global_threading(ae_uint64_t flg_value);

// Service functions.
typedef enum { CPU_SSE2 = 0x1, CPU_AVX2 = 0x2, CPU_FMA = 0x4 } ae_cpuid_t;
extern const ae_cpuid_t CurCPU;
ae_int_t ae_cores_count();
ae_int_t ae_get_effective_workers(ae_int_t nworkers);

// Id values for ae_[sg]et_dbg_value().
typedef enum {
// get set
   _ALGLIB_ALLOC_COUNTER = 0, _ALGLIB_TOTAL_ALLOC_SIZE = 1,
// get
   _ALGLIB_TOTAL_ALLOC_COUNT = 2,
// set
   _ALGLIB_USE_VENDOR_KERNELS = 100, _ALGLIB_VENDOR_MEMSTAT = 101,
   _ALGLIB_DEBUG_WORKSTEALING = 200, _ALGLIB_WSDBG_NCORES = 201,
   _ALGLIB_WSDBG_PUSHROOT_OK = 202, _ALGLIB_WSDBG_PUSHROOT_FAILED = 203,
// get
   _ALGLIB_CORES_COUNT = 1000,
// get set
   _ALGLIB_GLOBAL_THREADING = 1001, _ALGLIB_NWORKERS = 1002,
#if 0
// For compatibility:
// get
   _ALGLIB_GET_ALLOC_COUNTER = _ALGLIB_ALLOC_COUNTER,
   _ALGLIB_GET_CUMULATIVE_ALLOC_SIZE = _ALGLIB_TOTAL_ALLOC_SIZE,
   _ALGLIB_GET_CUMULATIVE_ALLOC_COUNT = _ALGLIB_TOTAL_ALLOC_COUNT,
   _ALGLIB_GET_CORES_COUNT = _ALGLIB_CORES_COUNT,
   _ALGLIB_GET_GLOBAL_THREADING = _ALGLIB_GLOBAL_THREADING,
   _ALGLIB_GET_NWORKERS = _ALGLIB_NWORKERS,
// set
   _ALGLIB_USE_ALLOC_COUNTER = _ALGLIB_ALLOC_COUNTER,
   _ALGLIB_USE_DBG_COUNTERS = _ALGLIB_TOTAL_ALLOC_SIZE,
   _ALGLIB_SET_GLOBAL_THREADING = _ALGLIB_GLOBAL_THREADING,
   _ALGLIB_SET_NWORKERS = _ALGLIB_NWORKERS,
#endif
} debug_flag_t;
ae_int64_t ae_get_dbg_value(debug_flag_t id);
void ae_set_dbg_value(debug_flag_t flag_id, ae_int64_t flag_val);

int tickcount();

ae_int_t ae_misalignment(const void *ptr, size_t alignment);
void *ae_align(void *ptr, size_t alignment);
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
void set_memory_pool(void *ptr, size_t size);
void memory_pool_stats(ae_int_t *bytes_used, ae_int_t *bytes_free);
#endif
void *aligned_malloc(size_t size, size_t alignment);
void aligned_free(void *block);

void *ae_malloc(size_t size);
void ae_free(void *p);

// Dynamic block which may be automatically deallocated during stack unwinding.
typedef void (*ae_deallocator)(void *);
struct ae_dyn_block {
// The next block in the stack unwinding list; or NULL if this block is not in the list.
   struct ae_dyn_block *volatile p_next;
// The block deallocation function; or NULL for the stack frame/boundary "special" blocks.
   ae_deallocator deallocator; // Was: void *deallocator;
// The argument for deallocator(); or NULL for 0-size blocks, or the DYN_BOTTOM or DYN_FRAME stack frame/boundary "special" blocks.
   void *volatile ptr;
};
void ae_db_init(ae_dyn_block *block, ae_int_t size, bool make_automatic);
void ae_db_realloc(ae_dyn_block *block, ae_int_t size);
void ae_db_free(ae_dyn_block *block);
void ae_db_swap(ae_dyn_block *block1, ae_dyn_block *block2);
#define NewBlock(B, N)			ae_dyn_block B; memset(&B, 0, sizeof B), ae_db_init(&B, N, true)

ae_int_t ae_sizeof(ae_datatype datatype);

struct ae_vector {
// The number of elements in the vector; cnt >= 0.
   ae_int_t cnt;
// Either DT_BOOL/DT_BYTE, DT_INT, DT_REAL or DT_COMPLEX.
   ae_datatype datatype;
// True if and only if the ae_vector was attached to an x_vector by ae_vector_init_attach_to_x().
   bool is_attached;
// The ae_dyn_block structure which manages the data in xX and deletes it when its frame is freed.
   ae_dyn_block data;
// A generic pointer or pointer to the data with the matching datatype.
// These are the fields that the user usually works with.
   union {
      void *xX;
      bool *xB;
      unsigned char *xU;
      ae_int_t *xZ;
      double *xR;
      complex *xC;
   };
};
void ae_vector_init(ae_vector *dst, ae_int_t size, ae_datatype datatype, bool make_automatic);
void ae_vector_copy(ae_vector *dst, ae_vector *src, bool make_automatic);
void ae_vector_set_length(ae_vector *dst, ae_int_t newsize);
void ae_vector_resize(ae_vector *dst, ae_int_t newsize);
void ae_vector_free(ae_vector *dst, bool make_automatic);
void ae_swap_vectors(ae_vector *vec1, ae_vector *vec2);
#define NewVector(V, N, Type)		ae_vector V; memset(&V, 0, sizeof V), ae_vector_init(&V, N, Type, true)
#define DupVector(V)			ae_vector _##V; memset(&_##V, 0, sizeof _##V), ae_vector_copy(&_##V, V, true), V = &_##V
#define SetVector(P)			ae_vector_free(P, true)

struct ae_matrix {
   ae_int_t cols;
   ae_int_t rows;
   ae_int_t stride;
   ae_datatype datatype;
// True if and only if the ae_matrix was attached to an x_matrix by ae_matrix_init_attach_to_x().
   bool is_attached;
   ae_dyn_block data;
// A generic pointer, raster pointer or pointer to the data with the matching datatype.
   union {
      void *xX;
      void **xyX;
      bool **xyB;
      ae_int_t **xyZ;
      double **xyR;
      complex **xyC;
   };
};
void ae_matrix_init(ae_matrix *dst, ae_int_t rows, ae_int_t cols, ae_datatype datatype, bool make_automatic);
void ae_matrix_copy(ae_matrix *dst, ae_matrix *src, bool make_automatic);
void ae_matrix_set_length(ae_matrix *dst, ae_int_t rows, ae_int_t cols);
void ae_matrix_free(ae_matrix *dst, bool make_automatic);
void ae_swap_matrices(ae_matrix *mat1, ae_matrix *mat2);
#define NewMatrix(M, Ys, Xs, Type)	ae_matrix M; memset(&M, 0, sizeof M), ae_matrix_init(&M, Ys, Xs, Type, true)
#define DupMatrix(M)			ae_matrix _##M; memset(&_##M, 0, sizeof _##M), ae_matrix_copy(&_##M, M, true), M = &_##M
#define SetMatrix(P)			ae_matrix_free(P, true)

// Used for better documenting function parameters.
// TODO: Remake ae_vector and ae_matrix as template types.
typedef ae_vector BVector, ZVector, RVector, CVector;
typedef ae_matrix BMatrix, ZMatrix, RMatrix, CMatrix;

struct ae_smart_ptr {
// Pointers respectively to the subscriber and the object; all changes in ptr are translated to subscriber.
   void **subscriber, *ptr;
// True if the smart pointer owns ptr.
   bool is_owner;
// True if ptr points to a dynamic object, whose clearing requires calling both .free() and ae_free().
   bool is_dynamic;
// The deallocation function for the pointer; clears all dynamically allocated memory.
   void (*free)(void *, bool make_automatic);
// The frame entry; used to ensure automatic deallocation of the smart pointer in case of an exception/exit.
   ae_dyn_block frame_entry;
};
void ae_smart_ptr_init(ae_smart_ptr *dst, void **subscriber, bool make_automatic);
void ae_smart_ptr_free(void *_dst); // Accepts ae_smart_ptr *.
void ae_smart_ptr_assign(ae_smart_ptr *dst, void *new_ptr, bool is_owner, bool is_dynamic, void (*free)(void *, bool make_automatic));
void ae_smart_ptr_release(ae_smart_ptr *dst);
#define NewObj(Type, P)	Type P; memset(&P, 0, sizeof P), Type##_init(&P, true)
#define RefObj(Type, P)	Type *P; alglib_impl::ae_smart_ptr _##P; memset(&_##P, 0, sizeof _##P), alglib_impl::ae_smart_ptr_init(&_##P, (void **)&P, true)
#define SetObj(Type, P)	Type##_free(P, true)

// The X-interface.
// The effective type for the last_action field.
enum { ACT_UNCHANGED = 1, ACT_SAME_LOCATION = 2, ACT_NEW_LOCATION = 3 };

#if 0 //(@) Not used anywhere.
// x-strings (zero-terminated): members are ae_int64_t aligned to avoid alignment problems.
// Compiler-specific alignment definitions.
#   if AE_COMPILER == AE_GNUC
#      define ALIGNED __attribute__((aligned(8)))
#   elif AE_COMPILER == AE_MSVC
#      define ALIGNED __declspec(align(8))
#   else
#      define ALIGNED
#   endif
struct x_string {
// Determines what to do on realloc().
// If the object is owned by the caller, the X-interface will just set ptr to NULL before realloc().
// If it is owned by X, it will call one of the *_free() functions.
   ALIGNED ae_int64_t owner;		// bool owner;
// Set on return from the X interface and may be used by the caller as a hint for deciding what to do with the buffer.
//	ACT_UNCHANGED:		unchanged,
//	ACT_SAME_LOCATION:	stored at the same location, or
//	ACT_NEW_LOCATION:	stored at a new location.
// ACT_{UNCHANGED,SAME_LOCATION} mean no reallocation or copying is required.
   ALIGNED ae_int64_t last_action;	// enum { ACT_UNCHANGED = 1, ACT_SAME_LOCATION, ACT_NEW_LOCATION } last_action;
// A pointer to the actual data.
   ALIGNED char *ptr;			// union { void *ptr; ae_int64_t portable_alignment_enforcer; };
};
#endif

// x-vectors: members are ae_int64_t aligned to avoid alignment problems.
struct x_vector {
// The vector size; i.e., the number of elements.
   ae_int64_t cnt;		// ae_int_t cnt;
// One of the DT_* values.
   ae_int64_t datatype;		// ae_datatype datatype;
// Determines what to do on realloc().
// If the object is owned by the caller, the X-interface will just set ptr to NULL before realloc().
// If it is owned by X, it will call one of the *_free() functions.
   ae_int64_t owner;		// bool owner;
// Set on return from the X interface and may be used by the caller as a hint for deciding what to do with the buffer.
//	ACT_UNCHANGED:		unchanged,
//	ACT_SAME_LOCATION:	stored at the same location, or
//	ACT_NEW_LOCATION:	stored at a new location.
// ACT_{UNCHANGED,SAME_LOCATION} mean no reallocation or copying is required.
   ae_int64_t last_action;	// enum { ACT_UNCHANGED = 1, ACT_SAME_LOCATION, ACT_NEW_LOCATION } last_action;
// A pointer to the actual data - with enforced alignment.
   union {
      void *x_ptr;
      ae_int64_t portable_alignment_enforcer;
   };
};
void ae_vector_init_from_x(ae_vector *dst, x_vector *src, bool make_automatic);
void ae_vector_init_attach_to_x(ae_vector *dst, x_vector *src, bool make_automatic);
void ae_x_set_vector(x_vector *dst, ae_vector *src);
void ae_x_attach_to_vector(x_vector *dst, ae_vector *src);
void x_vector_free(x_vector *dst, bool make_automatic);

// x-matrices: members are ae_int64_t aligned to avoid alignment problems.
struct x_matrix {
// The matrix size.
// If either dimension is 0, then both must be.
   ae_int64_t cols, rows;	// ae_int_t cols, rows;
// The stride, i.e. the byte-distance between the first elements of adjacent matrix rows.
   ae_int64_t stride;		// ae_int_t stride;
// One of the DT_* values.
   ae_int64_t datatype;		// ae_datatype datatype;
// Determines what to do on realloc().
// If the object is owned by the caller, the X-interface will just set ptr to NULL before realloc().
// If it is owned by X, it will call one of the *_free() functions.
   ae_int64_t owner;		// bool owner;
// Set on return from the X interface and may be used by the caller as a hint for deciding what to do with the buffer.
//	ACT_UNCHANGED:		unchanged,
//	ACT_SAME_LOCATION:	stored at the same location, or
//	ACT_NEW_LOCATION:	stored at a new location.
// ACT_{UNCHANGED,SAME_LOCATION} mean no reallocation or copying is required.
   ae_int64_t last_action;	// enum { ACT_UNCHANGED = 1, ACT_SAME_LOCATION, ACT_NEW_LOCATION } last_action;
// A pointer to the actual matrix data, stored rowwise - with enforced alignment.
   union {
      void *x_ptr;
      ae_int64_t portable_alignment_enforcer;
   };
};
void ae_matrix_init_from_x(ae_matrix *dst, x_matrix *src, bool make_automatic);
void ae_matrix_init_attach_to_x(ae_matrix *dst, x_matrix *src, bool make_automatic);
void ae_x_set_matrix(x_matrix *dst, ae_matrix *src);
void ae_x_attach_to_matrix(x_matrix *dst, ae_matrix *src);

bool ae_is_symmetric(ae_matrix *a);
bool ae_is_hermitian(ae_matrix *a);
bool ae_force_symmetric(ae_matrix *a);
bool ae_force_hermitian(ae_matrix *a);

// An OS-independent non-reentrant lock:
// *	under Posix and Windows systems it uses the system-provided locks,
// *	under Boost it uses the OS-independent lock provided by the Boost package,
// *	when no OS is defined, it uses a "fake lock" (just a stub, which is not thread-safe):
//	a)	the "fake lock" can be in locked or free mode,
//	b)	the "fake lock" can be used only from the thread which created it,
//	c)	when the thread acquires a free lock, it immediately returns,
//	d)	when the thread acquires a busy lock, the program is terminated
//		(because the lock is already acquired and no one else can free it).
struct ae_lock {
// Pointer to _lock structure.
// This pointer has type void * in order to make the header file OS-independent (the lock declaration depends on the OS).
   void *lock_ptr;
// For is_static == false this field manages the pointer to the _lock structure.
// The frame that is responsible for deallocation when the lock's frame is freed.
   ae_frame db;
// True if we have a static lock (used by a thread pool) or transient lock.
// Static locks are allocated without using the frame and cannot be deallocated.
   bool is_static;
};
void ae_init_lock(ae_lock *lock, bool is_static, bool make_automatic);
void ae_acquire_lock(ae_lock *lock);
void ae_release_lock(ae_lock *lock);
void ae_free_lock(ae_lock *lock);

// Shared pool: the data structure used to provide thread-safe access to a pool of temporary variables.
struct ae_shared_pool_entry {
   void *volatile obj, *volatile next_entry;
};
struct ae_shared_pool {
// The lock object which protects the pool.
   ae_lock pool_lock;
// The seed object: used to create new instances of temporaries.
   void *volatile seed_object;
// The list of recycled Objects:
// *	entries in this list store pointers to recycled objects,
// *	move from this list the first entry of each object retrieved to recycled_entries, and return the caller its obj field.
   ae_shared_pool_entry *volatile recycled_objects;
// The list of recycled Entries:
// *	entries which are not used to store recycled objects;
//	move the entry of each recycled object retrieved to this list,
// *	get an entry for each object recycled from this list before allocating any with malloc().
   ae_shared_pool_entry *volatile recycled_entries;
// The enumeration pointer, points to the current recycled object.
   ae_shared_pool_entry *volatile enumeration_counter;
// The size of the object; used when we call malloc() for new objects.
   ae_int_t size_of_object;
// The initializer function; accepts a pointer to the malloc'ed object, initializes its fields.
   void (*init)(void *dst, bool make_automatic);
// The copy constructor; accepts a pointer to the malloc'ed object, but not to the initialized object.
   void (*copy)(void *dst, void *src, bool make_automatic);
// The destructor function.
   void (*free)(void *ptr, bool make_automatic);
// The frame entry; points to the pool object itself.
   ae_frame frame_entry;
};
void ae_shared_pool_init(void *_dst, bool make_automatic);
void ae_shared_pool_copy(void *_dst, void *_src, bool make_automatic);
void ae_shared_pool_free(void *_dst, bool make_automatic);
bool ae_shared_pool_is_initialized(ae_shared_pool *dst);
void ae_shared_pool_set_seed(ae_shared_pool *dst, void *seed_object, ae_int_t size_of_object, void (*init)(void *dst, bool make_automatic), void (*copy)(void *dst, void *src, bool make_automatic), void (*free)(void *ptr, bool make_automatic));
void ae_shared_pool_retrieve(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_recycle(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_clear_recycled(ae_shared_pool *pool, bool make_automatic);
void ae_shared_pool_first_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_next_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_reset(ae_shared_pool *pool);

// Serializer:
// *	ae_stream_writer is the type expected for pointers to serializer stream-writing functions,
//	which are used by the X-core for out-of-core serialization (e.g., directly to a C++ stream).
//	This function accepts:
//	-	S:	a pointer to an ANSI (7-bit) string, which is a part of the data stream and
//			and may include spaces and linefeed symbols, and should be written to the stream as is; and
//	-	Aux:	a pointer-sized integer passed to the serializer during initialization,
//			which may be any value meant to be used by the actual implementation of the stream writer.
//	The function should return true for success or false for failure.
typedef bool (*ae_stream_writer)(const char *S, ae_int_t Aux);
// *	ae_stream_reader is the type expected for pointers to serializer stream-reading functions,
//	which are used by the X-core for out-of-core unserialization (e.g., directly from a C++ stream).
//	This function accepts:
//	-	Aux:	a pointer-sized integer passed to the serializer during initialization;
//	-	N:	the number of symbols to be read from the stream, with N > 0 required; and
//	-	S:	a pointer to the buffer used to store the next token read from the stream
//			(ANSI encoding is used, the buffer should be large enough to store all symbols and a trailing '\0').
//	After being called by the X-core, this function must:
//	-	skip all space and linefeed characters from the current position at the stream
//		up to the first non-space non-linefeed character
//	-	read exactly N symbols from the stream to the buffer; checking that they are all non-space non-linefeed ones
//	-	append a trailing '\0' to the buffer
//	-	return true for success, false if any of the conditions above fails,
//		in which case, the contents of S are not used.
typedef bool (*ae_stream_reader)(ae_int_t Aux, ae_int_t N, char *S);

typedef enum {
   AE_SM_DEFAULT = 0, AE_SM_ALLOC = 1, AE_SM_READY2S = 2,
   AE_SM_TO_STRING = 10, AE_SM_TO_CPPSTRING = 11, AE_SM_TO_STREAM = 12,
   AE_SM_FROM_STRING = 20, AE_SM_FROM_STREAM = 22
} ae_sermode_t;

struct ae_serializer {
   ae_sermode_t mode;
   ae_int_t entries_needed, bytes_asked;
   ae_int_t entries_saved, bytes_written;
#ifdef AE_USE_CPP_SERIALIZATION
   std::string *out_cppstr;
#endif
// Pointers respectively to the current position at the output/input buffers; advanced with each write/read operation.
   char *out_str; const char *in_str;
   ae_int_t stream_aux;
   ae_stream_writer stream_writer;
   ae_stream_reader stream_reader;
};
void ae_serializer_init(ae_serializer *serializer);
#define NewSerializer(Ser)	alglib_impl::ae_serializer Ser; alglib_impl::ae_serializer_init(&Ser)

void ae_serializer_alloc_start(ae_serializer *serializer);
void ae_serializer_alloc_entry(ae_serializer *serializer);
ae_int_t ae_serializer_get_alloc_size(ae_serializer *serializer);

#ifdef AE_USE_CPP_SERIALIZATION
void ae_serializer_sstart_str(ae_serializer *serializer, std::string *buf);
void ae_serializer_ustart_str(ae_serializer *serializer, const std::string *buf);
void ae_serializer_sstart_stream(ae_serializer *serializer, std::ostream *stream);
void ae_serializer_ustart_stream(ae_serializer *serializer, const std::istream *stream);
#endif
void ae_serializer_sstart_str(ae_serializer *serializer, char *buf);
void ae_serializer_ustart_str(ae_serializer *serializer, const char *buf);
void ae_serializer_sstart_stream(ae_serializer *serializer, ae_stream_writer writer, ae_int_t aux);
void ae_serializer_ustart_stream(ae_serializer *serializer, ae_stream_reader reader, ae_int_t aux);

bool ae_serializer_unserialize_bool(ae_serializer *serializer);
void ae_serializer_serialize_bool(ae_serializer *serializer, bool v);
ae_int_t ae_serializer_unserialize_int(ae_serializer *serializer);
void ae_serializer_serialize_int(ae_serializer *serializer, ae_int_t v);
ae_int64_t ae_serializer_unserialize_int64(ae_serializer *serializer);
void ae_serializer_serialize_int64(ae_serializer *serializer, ae_int64_t v);
double ae_serializer_unserialize_double(ae_serializer *serializer);
void ae_serializer_serialize_double(ae_serializer *serializer, double v);

void ae_serializer_stop(ae_serializer *serializer);

void ae_serializer_alloc_byte_array(ae_serializer *serializer, ae_vector *bytes);
void ae_serializer_unserialize_byte_array(ae_serializer *serializer, ae_vector *bytes);
void ae_serializer_serialize_byte_array(ae_serializer *serializer, ae_vector *bytes);

// Real math functions: IEEE-compliant floating point comparisons and standard functions.
// * IEEE-compliant floating point comparisons.
bool isneginf(double x);
bool isposinf(double x);

// * Standard functions.
ae_int_t imin2(ae_int_t x, ae_int_t y);
ae_int_t imin3(ae_int_t x, ae_int_t y, ae_int_t z);
ae_int_t imax2(ae_int_t x, ae_int_t y);
ae_int_t imax3(ae_int_t x, ae_int_t y, ae_int_t z);
ae_int_t ae_iabs(ae_int_t x);
ae_int_t sign(double x);
ae_int_t iround(double x);
ae_int_t itrunc(double x);
ae_int_t ifloor(double x);
ae_int_t iceil(double x);

double rmin2(double x, double y);
double rmax2(double x, double y);
double rmax3(double x, double y, double z);
double rmaxabs3(double x, double y, double z);
double sqr(double x);

ae_int_t iboundval(ae_int_t x, ae_int_t b1, ae_int_t b2);
double rboundval(double x, double b1, double b2);

// Debug-support.
#if 0 //(@) Not used anywhere, and not useful for debugging.
// Flush the console.
#   define flushconsole(s) fflush(stdout)
#endif

// Debug-enabled random number functions:
// TODO:
// *	ae_set_seed():	set the seed of the debug random number generator (NOT thread-safe!!!).
// *	ae_get_seed():	the seed value of the debug random number generator (NOT thread-safe!!!).
// *	ae_debugrng():	a random number generated with a high-quality random number generator.
double randomreal();
double randommid();
ae_int_t randominteger(ae_int_t maxv);
bool randombool(double p = 0.5);

// Complex math functions:
// *	basic arithmetic operations
// *	standard functions
inline complex complex_from_i(ae_int_t x, ae_int_t y = 0) { complex r; r.x = x, r.y = y; return r; }
inline complex complex_from_d(double x, double y = 0.0) { complex r; r.x = x, r.y = y; return r; }

complex ae_c_neg(complex A);
complex conj(complex A);
complex csqr(complex A);
double abscomplex(complex A);

bool ae_c_eq(complex A, complex B);
bool ae_c_neq(complex A, complex B);

complex ae_c_add(complex A, complex B);
complex ae_c_mul(complex A, complex B);
complex ae_c_sub(complex A, complex B);
complex ae_c_div(complex A, complex B);

bool ae_c_eq_d(complex A, double B);
bool ae_c_neq_d(complex A, double B);

complex ae_c_add_d(complex A, double B);
complex ae_c_mul_d(complex A, double B);
complex ae_c_sub_d(complex A, double B);
complex ae_c_d_sub(double A, complex B);
complex ae_c_div_d(complex A, double B);
complex ae_c_d_div(double A, complex B);

// Complex BLAS operations.
complex ae_v_cdotproduct(const complex *A, ae_int_t dA, const char *CjA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void ae_v_cmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void ae_v_cmoveneg(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void ae_v_cmoved(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha);
void ae_v_cmovec(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha);
void ae_v_cadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void ae_v_caddd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha);
void ae_v_caddc(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha);
void ae_v_csub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void ae_v_csubd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha);
void ae_v_csubc(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha);
void ae_v_cmuld(complex *A, ae_int_t dA, ae_int_t N, double Alpha);
void ae_v_cmulc(complex *A, ae_int_t dA, ae_int_t N, complex Alpha);

// Real BLAS operations.
double ae_v_dotproduct(const double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void ae_v_move(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void ae_v_moveneg(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void ae_v_moved(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha);
void ae_v_add(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void ae_v_addd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha);
void ae_v_sub(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void ae_v_subd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha);
void ae_v_muld(double *A, ae_int_t dA, ae_int_t N, double Alpha);

#if 0
extern const double machineepsilon, maxrealnumber, minrealnumber;
extern const double pi;
#else
#   define machineepsilon	5.0E-16
#   define maxrealnumber	1.0E300
#   define minrealnumber	1.0E-300
#   define pi			3.1415926535897932384626433832795
#endif

// Optimized shared C/C++ linear algebra code.
#define ALGLIB_INTERCEPTS_ABLAS
bool _ialglib_i_rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc);
bool _ialglib_i_cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, complex alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, complex beta, ae_matrix *c, ae_int_t ic, ae_int_t jc);
bool _ialglib_i_rmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_cmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_rmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_cmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool _ialglib_i_cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool _ialglib_i_rmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
bool _ialglib_i_cmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
bool _ialglib_i_rmatrixgerf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, double alpha, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
void _ialglib_pack_n2(double *col0, double *col1, ae_int_t n, ae_int_t src_stride, double *dst);
void _ialglib_mm22(double alpha, const double *a, const double *b, ae_int_t k, double beta, double *r, ae_int_t stride, ae_int_t store_mode);
void _ialglib_mm22x2(double alpha, const double *a, const double *b0, const double *b1, ae_int_t k, double beta, double *r, ae_int_t stride);

// ABLASF kernels.
bool ablasf_rgemm32basecase(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, RMatrix *a, ae_int_t ia, ae_int_t ja, ae_int_t opa, RMatrix *b, ae_int_t ib, ae_int_t jb, ae_int_t opb, double beta, RMatrix *c, ae_int_t ic, ae_int_t jc);
} // end of namespace alglib_impl

// Internal macros, defined only when _ALGLIB_IMPL_DEFINES is defined before inclusion of this header file.
#if defined _ALGLIB_IMPL_DEFINES
#   include "KernelsAvx2.h"
#   include "KernelsFma.h"
#   include "KernelsSse2.h"
#   define _ALGLIB_SIMD_ALIGNMENT_DOUBLES 8
#   define _ALGLIB_SIMD_ALIGNMENT_BYTES   (_ALGLIB_SIMD_ALIGNMENT_DOUBLES*8)
// SIMD kernel dispatchers.
#   if defined _ALGLIB_HAS_SSE2_INTRINSICS
#      define _KerSubSse2(Op) if (CurCPU & CPU_SSE2) { sse2_##Op; return; }
#      define _KerFunSse2(Op) if (CurCPU & CPU_SSE2) { return sse2_##Op; }
#   else
#      define _KerSubSse2(Op)
#      define _KerFunSse2(Op)
#   endif
#   if defined _ALGLIB_HAS_AVX2_INTRINSICS
#      define _KerSubAvx2(Op) if (CurCPU & CPU_AVX2) { avx2_##Op; return; }
#      define _KerFunAvx2(Op) if (CurCPU & CPU_AVX2) { return avx2_##Op; }
#   else
#      define _KerSubAvx2(Op)
#      define _KerFunAvx2(Op)
#   endif
#   if defined _ALGLIB_HAS_FMA_INTRINSICS
#      define _KerSubFma(Op) if (CurCPU & CPU_FMA) { fma_##Op; return; }
#      define _KerFunFma(Op) if (CurCPU & CPU_FMA) { return fma_##Op; }
#   else
#      define _KerSubFma(Op)
#      define _KerFunFma(Op)
#   endif
#   define KerSubSse2Avx2(Op) { _KerSubAvx2(Op) _KerSubSse2(Op) }
#   define KerFunSse2Avx2(Op) { _KerFunAvx2(Op) _KerFunSse2(Op) }
#   define KerSubSse2Avx2Fma(Op) { _KerSubFma(Op) _KerSubAvx2(Op) _KerSubSse2(Op) }
#   define KerFunSse2Avx2Fma(Op) { _KerFunFma(Op) _KerFunAvx2(Op) _KerFunSse2(Op) }
#   define KerSubAvx2Fma(Op) { _KerSubFma(Op) _KerSubAvx2(Op) }
#   define KerFunAvx2Fma(Op) { _KerFunFma(Op) _KerFunAvx2(Op) }
#   define KerSubAvx2(Op) { _KerSubAvx2(Op) }
#   define KerFunAvx2(Op) { _KerFunAvx2(Op) }
#   ifdef FP_FAST_FMA
#      define APPROX_FMA(x, y, z) fma((x), (y), (z))
#   else
#      define APPROX_FMA(x, y, z) ((x)*(y) + (z))
#   endif
#   if !defined ALGLIB_NO_FAST_KERNELS
// Arrays shorter than that will be processed with generic C implementation
#      if !defined _ABLASF_KERNEL_SIZE1
#         define _ABLASF_KERNEL_SIZE1 16
#      endif
#      if !defined _ABLASF_KERNEL_SIZE2
#         define _ABLASF_KERNEL_SIZE2 16
#      endif
#      define _ABLASF_BLOCK_SIZE 32
#      define _ABLASF_MICRO_SIZE  2
#      if defined _ALGLIB_HAS_AVX2_INTRINSICS || defined _ALGLIB_HAS_FMA_INTRINSICS
#         define ULOAD256PD(x) _mm256_loadu_pd((const double *)&x)
#      endif
#   endif
#endif

namespace alglib {
// Declarations for C++-related functionality.

// Exception class and exception-handling macros.
#if !defined AE_NO_EXCEPTIONS
// Exception-based code.
struct ap_error {
   std::string msg;
   ap_error();
   ap_error(const char *Msg);
   static void make_assertion(bool Cond);
   static void make_assertion(bool Cond, const char *Msg);
};
#   define ThrowError(Msg)	throw ap_error(Msg)
#   define ThrowErrorMsg(X)	throw ap_error()
#   define BegPoll		{ try {
#   define EndPoll		} catch(...) { alglib_impl::ae_clean_up(), alglib_impl::ae_state_clear(); throw; } }
#else
// Exception-free code.
#   define ThrowErrorMsg(X)	set_error_msg(); return X
//(@) The following restriction is unnecessary.
// #   if AE_OS != AE_OTHER_OS
// #      error Exception-free mode can not be combined with AE_OS definition
// #   endif
#   if AE_THREADING != NonTH
#      error Exception-free mode is thread-unsafe; define AE_THREADING = NonTH to prove that you know it.
#   endif
#   define BegPoll		{
#   define EndPoll		}
// Set the error flag and the pending error message.
void set_error_msg();
// Get the error flag and (optionally) the error message (as *MsgP);
// If the error flag is not set (or MsgP == NULL) *MsgP is not changed.
bool get_error_flag(const char **MsgP = NULL);
// Clear the error flag (it is not cleared until an explicit call to this function).
void clear_error_flag();
#endif
#define TryX		jmp_buf BreakAt; if (!setjmp(BreakAt)) alglib_impl::ae_state_set_break_jump(&BreakAt); else
#define TryCatch(X)	TryX { ThrowErrorMsg(X); }

// Class declaration/definition macros.
#define DecVal(X)	, X(Obj.X)
#define DecVar(X)	, X(&Obj.X)
#define DecComplex(X)	, X(*(complex *)&Obj.X)
#define ConstT(T, Val)	(const_cast<alglib_impl::T *>((Val).c_ptr()))
#define ComplexOf(Val)	(*reinterpret_cast<complex *>(&(Val)))

// Pseudo-templates for ALGLIB++ "Class" types:
// The triple-layering -- also present in the distribution version of ALGLIB --
// is the cost of trying to make a C++ type into a C structure, and to then turn around and wrap it back into a C++ type.
// It's better to just lose the C code (which is, effectively, what the alglib_impl namespace is and encapsulates)
// and write it all in C++, instead.
// I compromised by re-wrapping the C++ wrapper code into the DefClass() and DecClass() pseudo-templates
// and then removing the need for and calls to malloc()/free() that were originally in the member functions comprising DefClass().
// This CAN be quashed down to 1 layer ...
// by just writing the member functions directly in the corresponding types currently contained in the alglib_impl namespace,
// and using this, instead.
// However, it requires rewriting all the code involving smart pointers and shared pools as templates
// or by using C++'s native facilities for shared pools.
#define DecClass(Type, Pars) \
struct Type##I { \
   Type##I(); \
   Type##I(const Type##I &A); \
protected: \
   alglib_impl::Type Obj; \
}; \
struct Type: public Type##I { \
   Type(); \
   Type(const Type &A); \
   Type &operator=(const Type &A); \
   ~Type() { alglib_impl::Type##_free(&Obj, false); } \
   alglib_impl::Type *c_ptr() { return &Obj; } \
   alglib_impl::Type *c_ptr() const { return const_cast<alglib_impl::Type *>(&Obj); } \
   Pars \
}

#define DefClass(Type, Vars) \
Type##I::Type##I() { \
   alglib_impl::ae_state_init(); \
   TryX { alglib_impl::Type##_free(&Obj, false); ThrowErrorMsg(); } \
   memset(&Obj, 0, sizeof Obj), alglib_impl::Type##_init(&Obj, false); \
   alglib_impl::ae_state_clear(); \
} \
Type##I::Type##I(const Type##I &A) { \
   alglib_impl::ae_state_init(); \
   TryX { alglib_impl::Type##_free(&Obj, false); ThrowErrorMsg(); } \
   memset(&Obj, 0, sizeof Obj), alglib_impl::Type##_copy(&Obj, const_cast<alglib_impl::Type *>(&A.Obj), false); \
   alglib_impl::ae_state_clear(); \
} \
Type::Type(): Type##I() Vars { } \
Type::Type(const Type &A): Type##I(A) Vars { } \
Type &Type::operator=(const Type &A) { \
   if (this == &A) return *this; \
   alglib_impl::ae_state_init(); \
   TryCatch(*this) \
   alglib_impl::Type##_free(&Obj, false); \
   memset(&Obj, 0, sizeof Obj), alglib_impl::Type##_copy(&Obj, const_cast<alglib_impl::Type *>(&A.Obj), false); \
   alglib_impl::ae_state_clear(); \
   return *this; \
}

typedef alglib_impl::ae_int_t ae_int_t;

#if 0
// Constants and functions introduced for compatibility with AlgoPascal.
extern const double machineepsilon, maxrealnumber, minrealnumber;
extern const double pi;
#endif

bool isneginf(double x);
bool isposinf(double x);

int minint(int x, int y);
int maxint(int x, int y);
int sign(double x);
int iround(double x);
int itrunc(double x);
int ifloor(double x);
int iceil(double x);

double minreal(double x, double y);
double maxreal(double x, double y);
double sqr(double x);

double randomreal();
double randommid();
ae_int_t randominteger(ae_int_t maxv);
bool randombool(double p = 0.5);

// Complex number with double precision.
struct complex {
   complex(): x(0.0), y(0.0) { }
   complex(const double &X): x(X), y(0.0) { }
   complex(const double &X, const double &Y): x(X), y(Y) { }
   complex(const complex &Z): x(Z.x), y(Z.y) { }
   complex &operator=(const double &v);
   complex &operator=(const complex &z);
   complex &operator+=(const double &v);
   complex &operator+=(const complex &z);
   complex &operator-=(const double &v);
   complex &operator-=(const complex &z);
   complex &operator*=(const double &v);
   complex &operator*=(const complex &z);
   complex &operator/=(const double &v);
   complex &operator/=(const complex &z);
   alglib_impl::complex *c_ptr() { return (alglib_impl::complex *)this; }
   const alglib_impl::complex *c_ptr() const { return (const alglib_impl::complex *)this; }
#if !defined AE_NO_EXCEPTIONS
   std::string tostring(int _dps) const;
#endif
   double x, y;
};

bool operator==(const complex &A, const complex &B);
bool operator!=(const complex &A, const complex &B);
const complex operator+(const complex &A);
const complex operator-(const complex &A);
const complex operator+(const complex &A, const complex &B);
const complex operator+(const complex &A, const double &B);
const complex operator+(const double &A, const complex &B);
const complex operator-(const complex &A, const complex &B);
const complex operator-(const complex &A, const double &B);
const complex operator-(const double &A, const complex &B);
const complex operator*(const complex &A, const complex &B);
const complex operator*(const complex &A, const double &B);
const complex operator*(const double &A, const complex &B);
const complex operator/(const complex &A, const complex &B);
const complex operator/(const complex &A, const double &B);
const complex operator/(const double &A, const complex &B);
double abscomplex(const complex &A);
complex conj(const complex &A);
complex csqr(const complex &A);

// Multi-threading and multi-core functions.
// These are mostly stubs from the commercial version of ALGLIB.

// Get/Set the number of cores; can be 1, 2, ...; or 0 for auto; or -1/-2/... for all except for one/two/...
ae_int_t getnworkers();
void setnworkers(ae_int_t nworkers);

// Internal functions used by TestX.cpp, interfaces for functions present only in the commercial version of ALGLIB.
ae_int_t ae_cores_count();
alglib_impl::ae_uint64_t _ae_get_global_threading();
void _ae_set_global_threading(alglib_impl::ae_uint64_t flg_value);

// Set the global threading mode.
void setglobalthreading(const xparams settings);

// Level 1 BLAS functions
// NOTES:
// *	Segments A and B should NOT overlap.
// *	Strides dA, dB must be > 0, but are not assert()'ed within the function.
// *	CjA and CjB specify whether or not A and B respectively are to be conjugated before processing.
//	Any string which starts with 'N' or 'n' ("No conj", for example) specifies an unmodified parameter.
//	All other strings specify conjugation of the input, but "Conj" is recommended for use in such cases.
double vdotproduct(const double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
double vdotproduct(const double *A, const double *B, ae_int_t N);
complex vdotproduct(const complex *A, ae_int_t dA, const char *CjA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
complex vdotproduct(const complex *A, const complex *B, ae_int_t N);

void vmove(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void vmove(double *A, const double *B, ae_int_t N);
void vmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void vmove(complex *A, const complex *B, ae_int_t N);

void vmoveneg(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void vmoveneg(double *A, const double *B, ae_int_t N);
void vmoveneg(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void vmoveneg(complex *A, const complex *B, ae_int_t N);

void vmove(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha);
void vmove(double *A, const double *B, ae_int_t N, double Alpha);
void vmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha);
void vmove(complex *A, const complex *B, ae_int_t N, double Alpha);
void vmove(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha);
void vmove(complex *A, const complex *B, ae_int_t N, complex Alpha);

void vadd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void vadd(double *A, const double *B, ae_int_t N);
void vadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void vadd(complex *A, const complex *B, ae_int_t N);
void vadd(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha);
void vadd(double *A, const double *B, ae_int_t N, double Alpha);
void vadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha);
void vadd(complex *A, const complex *B, ae_int_t N, double Alpha);
void vadd(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha);
void vadd(complex *A, const complex *B, ae_int_t N, complex Alpha);

void vsub(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N);
void vsub(double *A, const double *B, ae_int_t N);
void vsub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N);
void vsub(complex *A, const complex *B, ae_int_t N);
void vsub(double *A, ae_int_t dA, const double *B, ae_int_t dB, ae_int_t N, double Alpha);
void vsub(double *A, const double *B, ae_int_t N, double Alpha);
void vsub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, double Alpha);
void vsub(complex *A, const complex *B, ae_int_t N, double Alpha);
void vsub(complex *A, ae_int_t dA, const complex *B, ae_int_t dB, const char *CjB, ae_int_t N, complex Alpha);
void vsub(complex *A, const complex *B, ae_int_t N, complex Alpha);

void vmul(double *A, ae_int_t dA, ae_int_t N, double Alpha);
void vmul(double *A, ae_int_t N, double Alpha);
void vmul(complex *A, ae_int_t dA, ae_int_t N, double Alpha);
void vmul(complex *A, ae_int_t N, double Alpha);
void vmul(complex *A, ae_int_t dA, ae_int_t N, complex Alpha);
void vmul(complex *A, ae_int_t N, complex Alpha);

// Vectors and Matrices.
struct ae_vector_wrapper {
// A new zero-sized vector of the given datatype.
   ae_vector_wrapper(alglib_impl::ae_datatype datatype);
// A new object externally attached to e_ptr, with run-time type-checking.
// NOTE:
// *	An exception is thrown if datatype != e_ptr->datatype.
   ae_vector_wrapper(alglib_impl::ae_vector *e_ptr, alglib_impl::ae_datatype datatype);
// A copy of the object rhs (can be a reference to one of the derived classes), with run-time type-checking.
// NOTE:
// *	An exception is thrown if datatype != rhs.datatype.
   ae_vector_wrapper(const ae_vector_wrapper &rhs, alglib_impl::ae_datatype datatype);
// The destructor.
   virtual ~ae_vector_wrapper();
// If the vector was allocated with its own object, this function changes length, completely dropping the previous contents.
// It does not work (and throws an exception) for frozen proxy objects.
   void setlength(ae_int_t iLen);
// The element count.
   ae_int_t length() const;
// Access to the internal C-structure used by the C-core.
// Not intended for external use.
   alglib_impl::ae_vector *c_ptr() { return This; }
   const alglib_impl::ae_vector *c_ptr() const { return This; }
#if 0 //(@) Not implemented.
private:
   ae_vector_wrapper();
   ae_vector_wrapper(const ae_vector_wrapper &rhs);
   const ae_vector_wrapper &operator=(const ae_vector_wrapper &rhs);
#endif
protected:
// Attach a wrapper object externally to the X-object new_ptr;
// "frozen proxy" mode is activated (you can read/write, but can not reallocate and do not own the vector's memory).
// NOTES:
// *	The wrapper object is assumed to start out already initialized; all previously allocated memory is properly deallocated.
// *	The X-object at new_ptr is used only once;
//	after we fetch the pointer to memory and its size, this X-object is ignored and not referenced anymore.
//	So, you can pass pointers to temporary x-structures which are deallocated immediately after you call attach_to().
// *	The state structure is used for error reporting purposes (longjmp on errors).
   void attach_to(alglib_impl::x_vector *new_ptr);
// Assign rhs to the current object and return *this.
// It has several branches depending on the target object's status:
// *	For proxy objects, copy the data to the proxy.
//	Throw an exception if the source and target are not the same size.
// *	For non-proxy objects, clear the data allocated by the object and copy rhs into the target.
// NOTE:
// *	This function correctly handles assignments of the object to itself.
   const ae_vector_wrapper &assign(const ae_vector_wrapper &rhs);
#if !defined AE_NO_EXCEPTIONS
// Copy the vector of the indicated datatype from the string s into the current object.
// The current object is considered empty (this function should be called from the copy constructor).
   ae_vector_wrapper(const char *s, alglib_impl::ae_datatype datatype);
#endif
// The indicated C-structure used by the C-core and the inner ae_vector:
// *	This == &Obj:	the object owns the ae_vector and is to handle its deallocation.
// *	This != &Obj:	someone else owns the ae_vector and is to handle its deallocation;
//			while Obj is assumed to be uninitialized and is ignored.
   alglib_impl::ae_vector *This, Obj;
// True if and only if the vector's object was internally allocated and is owned.
// *	ae_vector's, directly owned or not, can be read and modified.
// *	Only internally allocated and owned ae_vector's can be freed and resized.
// *	If owner == false:
//	―	This == &Obj:	the ae_vector contains the object, but it points to externally allocated memory.
//				This memory is NOT owned by the ae_vector.
//	―	This != &Obj:	the ae_vector is contained somewhere else.
//				Neither the memory pointed by the ae_vector nor the ae_vector itself are owned by the ae_vector.
   bool owner;
};

struct boolean_1d_array: public ae_vector_wrapper {
   boolean_1d_array();
   boolean_1d_array(alglib_impl::ae_vector *p);
   boolean_1d_array(const boolean_1d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   boolean_1d_array(const char *s);
   std::string tostring() const;
#endif
   virtual ~boolean_1d_array();
   const boolean_1d_array &operator=(const boolean_1d_array &rhs);
   const bool &operator()(ae_int_t i) const; bool &operator()(ae_int_t i);
   const bool &operator[](ae_int_t i) const; bool &operator[](ae_int_t i);
// Allocate an iLen-vector, giving it a completely independent copy of the data at pContent.
   void setcontent(ae_int_t iLen, const bool *pContent);
// A pointer to internal memory.
   const bool *getcontent() const; bool *getcontent();
};

struct integer_1d_array: public ae_vector_wrapper {
   integer_1d_array();
   integer_1d_array(alglib_impl::ae_vector *p);
   integer_1d_array(const integer_1d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   integer_1d_array(const char *s);
   std::string tostring() const;
#endif
   virtual ~integer_1d_array();
   const integer_1d_array &operator=(const integer_1d_array &rhs);
   const ae_int_t &operator()(ae_int_t i) const; ae_int_t &operator()(ae_int_t i);
   const ae_int_t &operator[](ae_int_t i) const; ae_int_t &operator[](ae_int_t i);
// Allocate an iLen-vector, giving it a completely independent copy of the data at pContent.
   void setcontent(ae_int_t iLen, const ae_int_t *pContent);
// A pointer to internal memory.
   const ae_int_t *getcontent() const; ae_int_t *getcontent();
};

struct real_1d_array: public ae_vector_wrapper {
   real_1d_array();
   real_1d_array(alglib_impl::ae_vector *p);
   real_1d_array(const real_1d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   real_1d_array(const char *s);
   std::string tostring(int dps) const;
#endif
   virtual ~real_1d_array();
   const real_1d_array &operator=(const real_1d_array &rhs);
   const double &operator()(ae_int_t i) const; double &operator()(ae_int_t i);
   const double &operator[](ae_int_t i) const; double &operator[](ae_int_t i);
// Allocate an iLen-vector, giving it a completely independent copy of the data at pContent.
   void setcontent(ae_int_t iLen, const double *pContent);
// A pointer to internal memory.
   const double *getcontent() const; double *getcontent();
// Attach an array to the vector at pContent.
// *	No own memory is allocated, no data is copied.
// *	The pContent pointer should be valid as long as we work with the vector.
// After attaching the vector to external memory, it is "frozen":
// it is possible to read/write the vector elements, but not to resize it (no setlength() calls).
   void attach_to_ptr(ae_int_t iLen, double *pContent);
};

struct complex_1d_array: public ae_vector_wrapper {
   complex_1d_array();
   complex_1d_array(alglib_impl::ae_vector *p);
   complex_1d_array(const complex_1d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   complex_1d_array(const char *s);
   std::string tostring(int dps) const;
#endif
   virtual ~complex_1d_array();
   const complex_1d_array &operator=(const complex_1d_array &rhs);
   const complex &operator()(ae_int_t i) const; complex &operator()(ae_int_t i);
   const complex &operator[](ae_int_t i) const; complex &operator[](ae_int_t i);
// Allocate an iLen-vector, giving it a completely independent copy of the data at pContent.
   void setcontent(ae_int_t iLen, const complex *pContent);
// A pointer to internal memory.
   const complex *getcontent() const; complex *getcontent();
};

struct ae_matrix_wrapper {
// A new zero-sized matrix of the given datatype.
   ae_matrix_wrapper(alglib_impl::ae_datatype datatype);
// A new object externally attached to e_ptr, with run-time type-checking.
// NOTE:
// *	An exception is thrown if datatype != e_ptr->datatype.
   ae_matrix_wrapper(alglib_impl::ae_matrix *e_ptr, alglib_impl::ae_datatype datatype);
// A copy of the object rhs (can be a reference to one of the derived classes), with run-time type-checking.
// NOTE:
// *	An exception is thrown if datatype != rhs.datatype.
   ae_matrix_wrapper(const ae_matrix_wrapper &rhs, alglib_impl::ae_datatype datatype);
// The destructor.
   virtual ~ae_matrix_wrapper();
   void setlength(ae_int_t rows, ae_int_t cols);
// The column and row count and stride.
   ae_int_t cols() const;
   ae_int_t rows() const;
   ae_int_t getstride() const;
   bool isempty() const;
// Access to the internal c-structure used by the c-core.
// Not meant for external use.
   alglib_impl::ae_matrix *c_ptr() { return This; }
   const alglib_impl::ae_matrix *c_ptr() const { return This; }
#if 0 //(@) Not implemented.
private:
   ae_matrix_wrapper();
   ae_matrix_wrapper(const ae_matrix_wrapper &rhs);
   const ae_matrix_wrapper &operator=(const ae_matrix_wrapper &rhs);
#endif
protected:
// Attach a wrapper object externally to the X-object new_ptr;
// "frozen proxy" mode is activated (you can read/write, but can not reallocate and do not own the matrix's memory).
// NOTES:
// *	The wrapper object is assumed to start out already initialized; all previously allocated memory is properly deallocated.
// *	The X-object at new_ptr is used only once;
//	after we fetch the pointer to memory and its size, this X-object is ignored and not referenced anymore.
//	So, you can pass pointers to temporary x-structures which are deallocated immediately after you call attach_to().
// *	The state structure is used for error-handling purposes (longjmp on errors).
//	All previously allocated memory is correctly freed on error.
   void attach_to(alglib_impl::x_matrix *new_ptr);
#if 0 //(@) Not implemented.
// This function initializes matrix and allocates own memory storage.
// NOTE:
// *	initial state of wrapper object is assumed to be uninitialized;
//	if This != NULL on entry, it is considered critical error (abort is called).
   void init(ae_int_t rows, ae_int_t cols, alglib_impl::ae_datatype datatype);
#endif
// Assign rhs to the current object and return *this.
// It has several branches depending on the target object's status:
// *	For proxy objects, copy the data to the proxy.
//	Throw an exception if the source and target are not the same size.
// *	For non-proxy objects, clear the data allocated by the object and copy rhs into the target.
// NOTE:
// *	This function correctly handles assignments of the object to itself.
   const ae_matrix_wrapper &assign(const ae_matrix_wrapper &rhs);
#if !defined AE_NO_EXCEPTIONS
// Copy the matrix of the indicated datatype from the string s into the current object.
// The current object is considered empty (this function should be called from the copy constructor).
   ae_matrix_wrapper(const char *s, alglib_impl::ae_datatype datatype);
#endif
// The internal C-structure used by the C-core and the inner ae_matrix.
// *	This == &Obj:	the matrix owns the ae_matrix and is to handle its deallocation.
// *	This != &Obj:	someone else owns the ae_matrix and is to handle its deallocation;
//			while Obj is assumed to be uninitialized and to be ignored.
   alglib_impl::ae_matrix *This, Obj;
// True if and only if the matrix's object was internally allocated and is owned.
// *	ae_matrix's, directly owned or not, can be read and modified.
// *	Only internally allocated and owned matrix objects can be freed and resized.
// *	If owner == false:
//	-	This == &Obj:	the ae_matrix contains the object, but it points to externally allocated memory.
//				This memory is NOT owned by the ae_matrix.
//	-	This != &Obj:	the ae_matrix is contained somewhere else.
//				Neither the memory pointed by the ae_matrix nor ae_matrix itself are owned by the ae_matrix.
   bool owner;
};

struct boolean_2d_array: public ae_matrix_wrapper {
   boolean_2d_array();
   boolean_2d_array(alglib_impl::ae_matrix *p);
   boolean_2d_array(const boolean_2d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   boolean_2d_array(const char *s);
   std::string tostring() const;
#endif
   virtual ~boolean_2d_array();
   const boolean_2d_array &operator=(const boolean_2d_array &rhs);
   const bool &operator()(ae_int_t i, ae_int_t j) const; bool &operator()(ae_int_t i, ae_int_t j);
   const bool *operator[](ae_int_t i) const; bool *operator[](ae_int_t i);
// A new irows x icols matrix filled with a completely independent copy of the data at pContent.
   void setcontent(ae_int_t irows, ae_int_t icols, const bool *pContent);
};

struct integer_2d_array: public ae_matrix_wrapper {
   integer_2d_array();
   integer_2d_array(alglib_impl::ae_matrix *p);
   integer_2d_array(const integer_2d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   integer_2d_array(const char *s);
   std::string tostring() const;
#endif
   virtual ~integer_2d_array();
   const integer_2d_array &operator=(const integer_2d_array &rhs);
   const ae_int_t &operator()(ae_int_t i, ae_int_t j) const; ae_int_t &operator()(ae_int_t i, ae_int_t j);
   const ae_int_t *operator[](ae_int_t i) const; ae_int_t *operator[](ae_int_t i);
// A new irows x icols matrix filled with a completely independent copy of the data at pContent.
   void setcontent(ae_int_t irows, ae_int_t icols, const ae_int_t *pContent);
};

struct real_2d_array: public ae_matrix_wrapper {
   real_2d_array();
   real_2d_array(alglib_impl::ae_matrix *p);
   real_2d_array(const real_2d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   real_2d_array(const char *s);
   std::string tostring(int dps) const;
#endif
   virtual ~real_2d_array();
   const real_2d_array &operator=(const real_2d_array &rhs);
   const double &operator()(ae_int_t i, ae_int_t j) const; double &operator()(ae_int_t i, ae_int_t j);
   const double *operator[](ae_int_t i) const; double *operator[](ae_int_t i);
// A new irows x icols matrix filled with a completely independent copy of the data at pContent.
   void setcontent(ae_int_t irows, ae_int_t icols, const double *pContent);
// Attach an array to the matrix at pContent.
// *	Very little own memory is allocated - O(irows) bytes for the precomputed pointers;
//	rather than the O(icols irows) copying of data.
// *	The pContent pointer should be valid as long as we work with the matrix.
// After attaching the matrix to external memory, it is "frozen":
// it is possible to read/write the matrix elements, but not to resize it (no setlength() calls).
   void attach_to_ptr(ae_int_t irows, ae_int_t icols, double *pContent);
};

struct complex_2d_array: public ae_matrix_wrapper {
   complex_2d_array();
   complex_2d_array(alglib_impl::ae_matrix *p);
   complex_2d_array(const complex_2d_array &rhs);
#if !defined AE_NO_EXCEPTIONS
   complex_2d_array(const char *s);
   std::string tostring(int dps) const;
#endif
   virtual ~complex_2d_array();
   const complex_2d_array &operator=(const complex_2d_array &rhs);
   const complex &operator()(ae_int_t i, ae_int_t j) const; complex &operator()(ae_int_t i, ae_int_t j);
   const complex *operator[](ae_int_t i) const; complex *operator[](ae_int_t i);
// A new irows x icols matrix filled with a completely independent copy of the data at pContent.
   void setcontent(ae_int_t irows, ae_int_t icols, const complex *pContent);
};

#if !defined AE_NO_EXCEPTIONS
// CSV operation.
extern const int CSV_DEFAULT, CSV_SKIP_HEADERS;

// Read a CSV file into a double precision matrix.
// The format of the data file must conform to the RFC 4180 specification, with additional notes:
// *	The file size should be less than 2GB.
// *	ASCII encoding, UTF-8 without BOM (in header names) are supported
// *	Any character (comma/tab/space) may be used as a field separator, as long as it is distinct from one used for decimal point.
// *	Multiple subsequent field separators (say, two spaces) are treated as MULTIPLE separators, not one big separator.
// *	Either a comma and full stop may be used as a decimal point.
//	The parser will automatically determine specific character being used.
//	Both fixed and exponential number formats are allowed.
//	Thousand separators are NOT allowed: the formatting used with separators is actually locale-specific, not universal.
// *	A line may end with \n (Unix style) or \r\n (Windows style), the parser will automatically adapt to chosen convention.
// *	Escaped fields (ones in double quotes) are not supported.
// Inputs:
//	filename	The relative/absolute path.
//	separator	The character used to separate fields.
//			May be ' ', ',', '\t'.
//			Other separators are possible too.
//	flags		Several values combined with bitwise OR:
//			*	CSV_SKIP_HEADERS - if present, the first row contains headers and will be skipped.
//				Its contents are used to determine field counts, and that's all.
//			If no flags are specified, the default value 0x0 (or CSV_DEFAULT, which is same) should be used.
// Outputs:
//	out		The matrix, CSV file parsed with atof()
// HANDLING OF SPECIAL CASES:
// *	The file does not exist - an ap_error exception is thrown.
// *	The file is empty - an empty array is returned (no exception).
// *	skip_first_row == true, there is only one row in the file - an empty array is returned.
// *	The field contents are not recognized by atof() - the field value is replaced by 0.0.
void read_csv(const char *filename, char separator, int flags, real_2d_array &out);
#endif
} // end of namespace alglib

#endif // OnceOnly
