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

#define AutoS static // "auto static": static thread-local variables in block or file scope.

#define AE_USE_CPP

// Definitions
#define AE_OTHER_CPU	0
#define AE_INTEL	1
#define AE_SPARC	2

// OS definitions
#define AE_OTHER_OS	0
#define AE_WINDOWS	1
#define AE_POSIX	2
#define AE_LINUX	304
#if !defined AE_OS
#   define AE_OS AE_UNKNOWN
#elif AE_OS == AE_LINUX
#   undef AE_OS
#   define AE_OS AE_POSIX
#   define _ALGLIB_USE_LINUX_EXTENSIONS
#endif

// threading models for AE_THREADING
#define AE_PARALLEL                 100
#define AE_SERIAL                   101
#define AE_SERIAL_UNSAFE            102
#if !defined AE_THREADING
#   define AE_THREADING AE_PARALLEL
#endif

// malloc types for AE_MALLOC
#define AE_STDLIB_MALLOC            200
#define AE_BASIC_STATIC_MALLOC      201
#if !defined AE_MALLOC
#   define AE_MALLOC AE_STDLIB_MALLOC
#endif

#define AE_LOCK_ALIGNMENT 16

// The compiler type and compiler-specific definitions.
#define AE_OTHERC	0
#define AE_MSVC		1
#define AE_GNUC		2
#define AE_SUNC		3
#if defined __GNUC__
#   define AE_COMPILER AE_GNUC
#   define ALIGNED __attribute__((aligned(8)))
#elif defined __SUNPRO_C || defined __SUNPRO_CC
#   define AE_COMPILER AE_SUNC
#   define ALIGNED
#elif defined _MSC_VER
#   define AE_COMPILER AE_MSVC
#   define ALIGNED __declspec(align(8))
// Disable some irrelevant warnings.
//(#) Originally preceded the headers other than Ap.h in TestC.cpp and followed them in the program module *.cpp files.
#   if defined InAlgLib && !defined AE_ALL_WARNINGS
#      pragma warning(disable:4100)
#      pragma warning(disable:4127)
#      pragma warning(disable:4611)
#      pragma warning(disable:4702)
#      pragma warning(disable:4996)
#   endif
#else
#   define AE_COMPILER AE_OTHERC
#   define ALIGNED
#endif

// state flags
#define _ALGLIB_FLG_THREADING_MASK          0x7
#define _ALGLIB_FLG_THREADING_SHIFT         0
#define _ALGLIB_FLG_THREADING_USE_GLOBAL    0x0
#define _ALGLIB_FLG_THREADING_SERIAL        0x1
#define _ALGLIB_FLG_THREADING_PARALLEL      0x2

// Now we are ready to include the remaining headers.
#include <ctype.h>
#include <string.h>
#include <setjmp.h>
#if defined AE_HAVE_STDINT
#   include <stdint.h>
#endif

// SSE2 intrinsics
// The preprocessor directives below:
// -	include headers for SSE2 intrinsics
// -	define AE_HAS_SSE2_INTRINSICS definition
// These actions are performed when we have:
// -	an x86 architecture definition AE_CPU == AE_INTEL
// -	a compiler which supports intrinsics
// The presence of AE_HAS_SSE2_INTRINSICS does NOT mean that our CPU actually supports SSE2 -
// such things should be determined by CurCPU, which is intialized on start-up.
// It means that we are working under Intel and out compiler can issue SSE2-capable code.
#if defined AE_CPU && AE_CPU == AE_INTEL
#   if AE_COMPILER == AE_MSVC
#      include <emmintrin.h>
#      define AE_HAS_SSE2_INTRINSICS
#   elif AE_COMPILER == AE_GNUC
#      include <xmmintrin.h>
#      define AE_HAS_SSE2_INTRINSICS
#   elif AE_COMPILER == AE_SUNC
#      include <xmmintrin.h>
#      include <emmintrin.h>
#      define AE_HAS_SSE2_INTRINSICS
#   endif
#endif

namespace alglib_impl {
// Core Code (Vectors, Matrices, Memory Management, etc.)

// if we work under C++ environment, define several conditions
#ifdef AE_USE_CPP
#   define AE_USE_CPP_BOOL
#   define AE_USE_CPP_SERIALIZATION
#   include <iostream>
#endif

// define ae_int32_t, ae_int64_t, ae_int_t, bool, ae_complex, ae_error_type and ae_datatype

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

#if !defined AE_USE_CPP_BOOL
typedef enum { false = 0, true = 1 } bool;
#endif

typedef struct { double x, y; } ae_complex;
typedef enum { ERR_OK = 0, ERR_OUT_OF_MEMORY = 1, ERR_XARRAY_TOO_LARGE = 2, ERR_ASSERTION_FAILED = 3 } ae_error_type;
typedef ae_int_t ae_datatype;

// other definitions
enum { OWN_CALLER = 1, OWN_AE = 2 };
enum { ACT_UNCHANGED = 1, ACT_SAME_LOCATION = 2, ACT_NEW_LOCATION = 3 };
enum { DT_BOOL = 1, DT_BYTE = 1, DT_INT = 2, DT_REAL = 3, DT_COMPLEX = 4 };
enum { CPU_SSE2 = 1 };

ae_int_t ae_misalignment(const void *ptr, size_t alignment);
void *ae_align(void *ptr, size_t alignment);
ae_int_t ae_get_effective_workers(ae_int_t nworkers);
void ae_optional_atomic_add_i(ae_int_t *p, ae_int_t v);
void ae_optional_atomic_sub_i(ae_int_t *p, ae_int_t v);

void *aligned_malloc(size_t size, size_t alignment);
void *aligned_extract_ptr(void *block);
void aligned_free(void *block);
void *eternal_malloc(size_t size);
#if AE_MALLOC == AE_BASIC_STATIC_MALLOC
void set_memory_pool(void *ptr, size_t size);
void memory_pool_stats(ae_int_t *bytes_used, ae_int_t *bytes_free);
#endif

void *ae_malloc(size_t size);
void ae_free(void *p);
ae_int_t ae_sizeof(ae_datatype datatype);
bool ae_check_zeros(const void *ptr, ae_int_t n);
void ae_touch_ptr(void *p);

// dynamic block which may be automatically deallocated during stack unwinding
//
// p_next          next block in the stack unwinding list.
//                 NULL means that this block is not in the list
// deallocator     deallocator function which should be used to deallocate block.
//                 NULL for "special" blocks (frame/stack boundaries)
// ptr             pointer which should be passed to the deallocator.
//                 may be null (for zero-size block), DYN_BOTTOM or DYN_FRAME
//                 for "special" blocks (frame/stack boundaries).
//
// valgrind_hint   is a special field which stores a special hint pointer for Valgrind and other similar memory checking tools.
//                 It is always set at run-time to to (ptr == NULL? NULL: aligned_extract_ptr(ptr)).
//                 This compensates to the manual alignment that ALGLIB does with pointers obtained via malloc,
//                 since ptr usually points to location past the beginning of the actuallly allocated memory.
//                 In such cases memory testing tools would report "(possibly) lost" memory.
//                 This field locates the actual base pointer
//                 and is meant only for use as a hint for Valgrind - NOT for anything else.
typedef void (*ae_deallocator)(void *);
typedef struct ae_dyn_block {
   struct ae_dyn_block *volatile p_next;
// void *deallocator;
   ae_deallocator deallocator;
   void *volatile ptr;
// void *valgrind_hint; // Redundant: it is always set to (ptr == NULL? NULL: aliged_extract_ptr(ptr)).
} ae_dyn_block;

void ae_db_attach(ae_dyn_block *block);
void ae_db_init(ae_dyn_block *block, ae_int_t size, bool make_automatic);
void ae_db_realloc(ae_dyn_block *block, ae_int_t size);
void ae_db_free(ae_dyn_block *block);
void ae_db_swap(ae_dyn_block *block1, ae_dyn_block *block2);
#define NewBlock(B, N)			ae_dyn_block B; memset(&B, 0, sizeof B), ae_db_init(&B, N, true)

void ae_state_init();
void ae_state_clear();
void ae_state_set_break_jump(jmp_buf *buf);
void ae_state_set_flags(ae_uint64_t flags);
void ae_break(ae_error_type error_type, const char *msg);

// frame marker
typedef struct ae_dyn_block ae_frame;

void ae_frame_make(ae_frame *tmp);
void ae_frame_leave();

typedef struct ae_vector {
// Number of elements in array, cnt >= 0
   ae_int_t cnt;

// Either DT_BOOL/DT_BYTE, DT_INT, DT_REAL or DT_COMPLEX
   ae_datatype datatype;

// If ptr points to memory owned and managed by ae_vector itself,
// this field is false. If vector was attached to x_vector structure
// with ae_vector_init_attach_to_x(), this field is true.
   bool is_attached;

// ae_dyn_block structure which manages data in ptr. This structure
// is responsible for automatic deletion of object when its frame
// is destroyed.
   ae_dyn_block data;

// Pointer to data.
// User usually works with this field.
   union {
      void *p_ptr;
      bool *p_bool;
      unsigned char *p_ubyte;
      ae_int_t *p_int;
      double *p_double;
      ae_complex *p_complex;
   } ptr;
} ae_vector;

void ae_vector_init(ae_vector *dst, ae_int_t size, ae_datatype datatype, bool make_automatic);
void ae_vector_copy(ae_vector *dst, ae_vector *src, bool make_automatic);
void ae_vector_set_length(ae_vector *dst, ae_int_t newsize);
void ae_vector_resize(ae_vector *dst, ae_int_t newsize);
void ae_vector_free(ae_vector *dst, bool make_automatic);
void ae_swap_vectors(ae_vector *vec1, ae_vector *vec2);
#define NewVector(V, N, Type)		ae_vector V; memset(&V, 0, sizeof V), ae_vector_init(&V, N, Type, true)
#define DupVector(V)			ae_vector _##V; memset(&_##V, 0, sizeof _##V), ae_vector_copy(&_##V, V, true), V = &_##V
#define SetVector(P)			ae_vector_free(P, true)

typedef struct ae_matrix {
   ae_int_t rows;
   ae_int_t cols;
   ae_int_t stride;
   ae_datatype datatype;

// If ptr points to memory owned and managed by ae_vector itself,
// this field is false. If vector was attached to x_vector structure
// with ae_vector_init_attach_to_x(), this field is true.
   bool is_attached;

   ae_dyn_block data;
   union {
      void *p_ptr;
      void **pp_void;
      bool **pp_bool;
      ae_int_t **pp_int;
      double **pp_double;
      ae_complex **pp_complex;
   } ptr;
} ae_matrix;

void ae_matrix_init(ae_matrix *dst, ae_int_t rows, ae_int_t cols, ae_datatype datatype, bool make_automatic);
void ae_matrix_copy(ae_matrix *dst, ae_matrix *src, bool make_automatic);
void ae_matrix_set_length(ae_matrix *dst, ae_int_t rows, ae_int_t cols);
void ae_matrix_free(ae_matrix *dst, bool make_automatic);
void ae_swap_matrices(ae_matrix *mat1, ae_matrix *mat2);
#define NewMatrix(M, Ys, Xs, Type)	ae_matrix M; memset(&M, 0, sizeof M), ae_matrix_init(&M, Ys, Xs, Type, true)
#define DupMatrix(M)			ae_matrix _##M; memset(&_##M, 0, sizeof _##M), ae_matrix_copy(&_##M, M, true), M = &_##M
#define SetMatrix(P)			ae_matrix_free(P, true)

// Used for better documenting function parameters.
typedef ae_vector *BVector, *ZVector, *RVector, *CVector;
typedef ae_matrix *BMatrix, *ZMatrix, *RMatrix, *CMatrix;

// Serializer:
//
// * ae_stream_writer type is a function pointer for stream  writer  method;
//   this pointer is used by X-core for out-of-core serialization  (say,  to
//   serialize ALGLIB structure directly to managed C# stream).
//
//   This function accepts two parameters: pointer to  ANSI  (7-bit)  string
//   and pointer-sized integer passed to serializer  during  initialization.
//   String being passed is a part of the data stream; aux paramerer may  be
//   arbitrary value intended to be used by actual implementation of  stream
//   writer. String parameter may include spaces and  linefeed  symbols,  it
//   should be written to stream as is.
//
//   Return value must be zero for success or non-zero for failure.
//
// * ae_stream_reader type is a function pointer for stream  reader  method;
//   this pointer is used by X-core for out-of-core unserialization (say, to
//   unserialize ALGLIB structure directly from managed C# stream).
//
//   This function accepts three parameters: pointer-sized integer passed to
//   serializer  during  initialization; number  of  symbols  to  read  from
//   stream; pointer to buffer used to store next  token  read  from  stream
//   (ANSI encoding is used, buffer is large enough to store all symbols and
//   trailing zero symbol).
//
//   Number of symbols to read is always positive.
//
//   After being called by X-core, this function must:
//   * skip all space and linefeed characters from the current  position  at
//     the stream and until first non-space non-linefeed character is found
//   * read exactly cnt symbols  from  stream  to  buffer;  check  that  all
//     symbols being read are non-space non-linefeed ones
//   * append trailing zero symbol to buffer
//   * return value must be zero on success, non-zero if  even  one  of  the
//     conditions above fails. When reader returns non-zero value,  contents
//     of buf is not used.
typedef char (*ae_stream_writer)(const char *p_string, ae_int_t aux);
typedef char (*ae_stream_reader)(ae_int_t aux, ae_int_t cnt, char *p_buf);

typedef struct {
   ae_int_t mode;
   ae_int_t entries_needed;
   ae_int_t entries_saved;
   ae_int_t bytes_asked;
   ae_int_t bytes_written;

#ifdef AE_USE_CPP_SERIALIZATION
   std::string * out_cppstr;
#endif
   char *out_str;       // Pointer to the current position at the output buffer; advanced with each write operation.
   const char *in_str;  // Pointer to the current position at the input  buffer; advanced with each read  operation.
   ae_int_t stream_aux;
   ae_stream_writer stream_writer;
   ae_stream_reader stream_reader;
} ae_serializer;

void ae_serializer_init(ae_serializer *serializer);
void ae_serializer_clear(ae_serializer *serializer);

void ae_serializer_alloc_start(ae_serializer *serializer);
void ae_serializer_alloc_entry(ae_serializer *serializer);
void ae_serializer_alloc_byte_array(ae_serializer *serializer, ae_vector *bytes);
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

void ae_serializer_serialize_bool(ae_serializer *serializer, bool v);
void ae_serializer_serialize_int(ae_serializer *serializer, ae_int_t v);
void ae_serializer_serialize_int64(ae_serializer *serializer, ae_int64_t v);
void ae_serializer_serialize_double(ae_serializer *serializer, double v);
void ae_serializer_serialize_byte_array(ae_serializer *serializer, ae_vector *bytes);
void ae_serializer_unserialize_bool(ae_serializer *serializer, bool *v);
void ae_serializer_unserialize_int(ae_serializer *serializer, ae_int_t *v);
void ae_serializer_unserialize_int64(ae_serializer *serializer, ae_int64_t *v);
void ae_serializer_unserialize_double(ae_serializer *serializer, double *v);
void ae_serializer_unserialize_byte_array(ae_serializer *serializer, ae_vector *bytes);

void ae_serializer_stop(ae_serializer *serializer);

// x-string (zero-terminated):
//     owner       OWN_CALLER or OWN_AE. Determines what to do on realloc().
//                 If vector is owned by caller, X-interface  will  just set
//                 ptr to NULL before realloc(). If it is  owned  by  X,  it
//                 will call ae_free/x_free/aligned_free family functions.
//
//     last_action ACT_UNCHANGED, ACT_SAME_LOCATION, ACT_NEW_LOCATION
//                 contents is either: unchanged, stored at the same location,
//                 stored at the new location.
//                 this field is set on return from X.
//
//     ptr         pointer to the actual data
//
// Members of this structure are ae_int64_t to avoid alignment problems.
typedef struct {
   ALIGNED ae_int64_t owner;
   ALIGNED ae_int64_t last_action;
   ALIGNED char *ptr;
} x_string;

// x-vector:
//     cnt         number of elements
//
//     datatype    one of the DT_XXXX values
//
//     owner       OWN_CALLER or OWN_AE. Determines what to do on realloc().
//                 If vector is owned by caller, X-interface  will  just set
//                 ptr to NULL before realloc(). If it is  owned  by  X,  it
//                 will call ae_free/x_free/aligned_free family functions.
//
//     last_action ACT_UNCHANGED, ACT_SAME_LOCATION, ACT_NEW_LOCATION
//                 contents is either: unchanged, stored at the same location,
//                 stored at the new location.
//                 this field is set on return from X interface and may be
//                 used by caller as hint when deciding what to do with data
//                 (if it was ACT_UNCHANGED or ACT_SAME_LOCATION, no array
//                 reallocation or copying is required).
//
//     ptr         pointer to the actual data
//
// Members of this structure are ae_int64_t to avoid alignment problems.
typedef struct {
   ae_int64_t cnt;
   ae_int64_t datatype;
   ae_int64_t owner;
   ae_int64_t last_action;
   union {
      void *p_ptr;
      ae_int64_t portable_alignment_enforcer;
   } x_ptr;
} x_vector;
void ae_vector_init_from_x(ae_vector *dst, x_vector *src, bool make_automatic);
void ae_vector_init_attach_to_x(ae_vector *dst, x_vector *src, bool make_automatic);
void ae_x_set_vector(x_vector *dst, ae_vector *src);
void ae_x_attach_to_vector(x_vector *dst, ae_vector *src);
void x_vector_free(x_vector *dst, bool make_automatic);

// x-matrix:
//     rows        number of rows. may be zero only when cols is zero too.
//
//     cols        number of columns. may be zero only when rows is zero too.
//
//     stride      stride, i.e. distance between first elements of rows (in bytes)
//
//     datatype    one of the DT_XXXX values
//
//     owner       OWN_CALLER or OWN_AE. Determines what to do on realloc().
//                 If vector is owned by caller, X-interface  will  just set
//                 ptr to NULL before realloc(). If it is  owned  by  X,  it
//                 will call ae_free/x_free/aligned_free family functions.
//
//     last_action ACT_UNCHANGED, ACT_SAME_LOCATION, ACT_NEW_LOCATION
//                 contents is either: unchanged, stored at the same location,
//                 stored at the new location.
//                 this field is set on return from X interface and may be
//                 used by caller as hint when deciding what to do with data
//                 (if it was ACT_UNCHANGED or ACT_SAME_LOCATION, no array
//                 reallocation or copying is required).
//
//     ptr         pointer to the actual data, stored rowwise
//
// Members of this structure are ae_int64_t to avoid alignment problems.
typedef struct {
   ae_int64_t rows;
   ae_int64_t cols;
   ae_int64_t stride;
   ae_int64_t datatype;
   ae_int64_t owner;
   ae_int64_t last_action;
   union {
      void *p_ptr;
      ae_int64_t portable_alignment_enforcer;
   } x_ptr;
} x_matrix;

void ae_matrix_init_from_x(ae_matrix *dst, x_matrix *src, bool make_automatic);
void ae_matrix_init_attach_to_x(ae_matrix *dst, x_matrix *src, bool make_automatic);
void ae_x_set_matrix(x_matrix *dst, ae_matrix *src);
void ae_x_attach_to_matrix(x_matrix *dst, ae_matrix *src);
bool x_is_symmetric(x_matrix *a);
bool x_is_hermitian(x_matrix *a);
bool x_force_symmetric(x_matrix *a);
bool x_force_hermitian(x_matrix *a);
bool ae_is_symmetric(ae_matrix *a);
bool ae_is_hermitian(ae_matrix *a);
bool ae_force_symmetric(ae_matrix *a);
bool ae_force_hermitian(ae_matrix *a);

typedef struct ae_smart_ptr {
// pointer to subscriber; all changes in ptr are translated to subscriber
   void **subscriber;
// pointer to object
   void *ptr;
// whether smart pointer owns ptr
   bool is_owner;
// whether object pointed by ptr is dynamic - clearing such object requires BOTH
// calling destructor function AND calling ae_free for memory occupied by object.
   bool is_dynamic;
// destructor function for pointer; clears all dynamically allocated memory
   void (*free)(void *, bool make_automatic);
// frame entry; used to ensure automatic deallocation of smart pointer in case of exception/exit
   ae_dyn_block frame_entry;
} ae_smart_ptr;

void ae_smart_ptr_init(ae_smart_ptr *dst, void **subscriber, bool make_automatic);
void ae_smart_ptr_free(void *_dst); // Accepts ae_smart_ptr*.
void ae_smart_ptr_assign(ae_smart_ptr *dst, void *new_ptr, bool is_owner, bool is_dynamic, void (*free)(void *, bool make_automatic));
void ae_smart_ptr_release(ae_smart_ptr *dst);
#define NewObj(Type, P)	Type P; memset(&P, 0, sizeof P), Type##_init(&P, true)
#define RefObj(Type, P)	Type *P; alglib_impl::ae_smart_ptr _##P; memset(&_##P, 0, sizeof _##P), alglib_impl::ae_smart_ptr_init(&_##P, (void **)&P, true)
#define SetObj(Type, P)	Type##_free(P, true)

// Lock.
//
// This structure provides OS-independent non-reentrant lock:
// * under Windows/Posix systems it uses system-provided locks
// * under Boost it uses OS-independent lock provided by Boost package
// * when no OS is defined, it uses "fake lock" (just stub which is not thread-safe):
//   a) "fake lock" can be in locked or free mode
//   b) "fake lock" can be used only from one thread - one which created lock
//   c) when thread acquires free lock, it immediately returns
//   d) when thread acquires busy lock, program is terminated
//      (because lock is already acquired and no one else can free it)
typedef struct {
// Pointer to _lock structure. This pointer has type void* in order to
// make header file OS-independent (lock declaration depends on OS).
   void *lock_ptr;
// For eternal=false this field manages pointer to _lock structure.
//
// ae_dyn_block structure is responsible for automatic deletion of
// the memory allocated for the pointer when its frame is destroyed.
   ae_dyn_block db;
// Whether we have eternal lock object (used by thread pool) or
// transient lock. Eternal locks are allocated without using ae_dyn_block
// structure and do not allow deallocation.
   bool eternal;
} ae_lock;

void ae_yield();
void ae_init_lock(ae_lock *lock, bool make_automatic);
void ae_init_lock_eternal(ae_lock *lock);
void ae_acquire_lock(ae_lock *lock);
void ae_release_lock(ae_lock *lock);
void ae_free_lock(ae_lock *lock);

// Shared pool: data structure used to provide thread-safe access to pool  of
// temporary variables.
typedef struct ae_shared_pool_entry {
   void *volatile obj;
   void *volatile next_entry;
} ae_shared_pool_entry;

typedef struct ae_shared_pool {
// lock object which protects pool
   ae_lock pool_lock;
// seed object (used to create new instances of temporaries)
   void *volatile seed_object;
// list of recycled OBJECTS:
// 1. entries in this list store pointers to recycled objects
// 2. every time we retrieve object, we retrieve first entry from this list,
//    move it to recycled_entries and return its obj field to caller/
   ae_shared_pool_entry *volatile recycled_objects;
// list of recycled ENTRIES:
// 1. this list holds entries which are not used to store recycled objects;
//    every time recycled object is retrieved, its entry is moved to this list.
// 2. every time object is recycled, we try to fetch entry for him from this list
//    before allocating it with malloc()
   ae_shared_pool_entry *volatile recycled_entries;
// enumeration pointer, points to current recycled object
   ae_shared_pool_entry *volatile enumeration_counter;
// size of object; this field is used when we call malloc() for new objects
   ae_int_t size_of_object;
// initializer function; accepts pointer to malloc'ed object, initializes its fields
   void (*init)(void *dst, bool make_automatic);
// copy constructor; accepts pointer to malloc'ed, but not initialized object
   void (*copy)(void *dst, void *src, bool make_automatic);
// destructor function;
   void (*free)(void *ptr, bool make_automatic);
// frame entry; contains pointer to the pool object itself
   ae_dyn_block frame_entry;
} ae_shared_pool;

void ae_shared_pool_init(void *_dst, bool make_automatic);
void ae_shared_pool_copy(void *_dst, void *_src, bool make_automatic);
void ae_shared_pool_free(void *dst, bool make_automatic);
bool ae_shared_pool_is_initialized(void *_dst);
void ae_shared_pool_set_seed(ae_shared_pool *dst, void *seed_object, ae_int_t size_of_object, void (*init)(void *dst, bool make_automatic), void (*copy)(void *dst, void *src, bool make_automatic), void (*free)(void *ptr, bool make_automatic));
void ae_shared_pool_retrieve(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_recycle(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_clear_recycled(ae_shared_pool *pool, bool make_automatic);
void ae_shared_pool_first_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_next_recycled(ae_shared_pool *pool, ae_smart_ptr *pptr);
void ae_shared_pool_reset(ae_shared_pool *pool);

void ae_never_call_it();
void ae_set_dbg_flag(ae_int64_t flag_id, ae_int64_t flag_val);
ae_int64_t ae_get_dbg_value(ae_int64_t id);
void ae_set_global_threading(ae_uint64_t flg_value);
ae_uint64_t ae_get_global_threading();

// Service functions
void ae_assert(bool cond, const char *msg);
extern const ae_int_t CurCPU;

// Real math functions:
// * standard functions
bool isposinf(double x);
bool isneginf(double x);
ae_int_t ae_iabs(ae_int_t x);
double ae_sqr(double x);
ae_int_t ae_sign(double x);
ae_int_t RoundZ(double x);
ae_int_t TruncZ(double x);
ae_int_t FloorZ(double x);
ae_int_t CeilZ(double x);
ae_int_t ae_maxint(ae_int_t m1, ae_int_t m2);
ae_int_t ae_minint(ae_int_t m1, ae_int_t m2);
double ae_maxreal(double m1, double m2);
double ae_minreal(double m1, double m2);
double ae_randomreal();
ae_int_t ae_randominteger(ae_int_t maxv);

// Complex math functions:
// * basic arithmetic operations
// * standard functions
ae_complex ae_complex_from_i(ae_int_t v);
ae_complex ae_complex_from_d(double v);

ae_complex ae_c_neg(ae_complex lhs);
bool ae_c_eq(ae_complex lhs, ae_complex rhs);
bool ae_c_neq(ae_complex lhs, ae_complex rhs);
ae_complex ae_c_add(ae_complex lhs, ae_complex rhs);
ae_complex ae_c_mul(ae_complex lhs, ae_complex rhs);
ae_complex ae_c_sub(ae_complex lhs, ae_complex rhs);
ae_complex ae_c_div(ae_complex lhs, ae_complex rhs);
bool ae_c_eq_d(ae_complex lhs, double rhs);
bool ae_c_neq_d(ae_complex lhs, double rhs);
ae_complex ae_c_add_d(ae_complex lhs, double rhs);
ae_complex ae_c_mul_d(ae_complex lhs, double rhs);
ae_complex ae_c_sub_d(ae_complex lhs, double rhs);
ae_complex ae_c_d_sub(double lhs, ae_complex rhs);
ae_complex ae_c_div_d(ae_complex lhs, double rhs);
ae_complex ae_c_d_div(double lhs, ae_complex rhs);

ae_complex ae_c_conj(ae_complex lhs);
ae_complex ae_c_sqr(ae_complex lhs);
double ae_c_abs(ae_complex z);

// Complex BLAS operations
ae_complex ae_v_cdotproduct(const ae_complex *v0, ae_int_t stride0, const char *conj0, const ae_complex *v1, ae_int_t stride1, const char *conj1, ae_int_t n);
void ae_v_cmove(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void ae_v_cmoveneg(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void ae_v_cmoved(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha);
void ae_v_cmovec(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, ae_complex alpha);
void ae_v_cadd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void ae_v_caddd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha);
void ae_v_caddc(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, ae_complex alpha);
void ae_v_csub(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void ae_v_csubd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha);
void ae_v_csubc(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, ae_complex alpha);
void ae_v_cmuld(ae_complex *vdst, ae_int_t stride_dst, ae_int_t n, double alpha);
void ae_v_cmulc(ae_complex *vdst, ae_int_t stride_dst, ae_int_t n, ae_complex alpha);

// Real BLAS operations
double ae_v_dotproduct(const double *v0, ae_int_t stride0, const double *v1, ae_int_t stride1, ae_int_t n);
void ae_v_move(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void ae_v_moveneg(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void ae_v_moved(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha);
void ae_v_add(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void ae_v_addd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha);
void ae_v_sub(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void ae_v_subd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha);
void ae_v_muld(double *vdst, ae_int_t stride_dst, ae_int_t n, double alpha);

// Other functions
#if 0
extern const double ae_machineepsilon, ae_maxrealnumber, ae_minrealnumber, ae_pi;
#endif
#define ae_machineepsilon 5E-16
#define ae_maxrealnumber  1E300
#define ae_minrealnumber  1E-300
#define ae_pi 3.1415926535897932384626433832795

// Allocation counters, inactive by default.
// Turned on when needed for debugging purposes.
//
// _alloc_counter is incremented by 1 on malloc(), decremented on free().
// _alloc_counter_total is only incremented by 1.
extern ae_int_t _alloc_counter;
extern ae_int_t _alloc_counter_total;
extern bool _use_alloc_counter;

// Malloc debugging:
//
// * _force_malloc_failure - set this flag to true in  order  to  enforce
//   failure of ALGLIB malloc(). Useful to debug handling of  errors  during
//   memory allocation. As long as this flag is set, ALGLIB malloc will fail.
// * _malloc_failure_after - set it to non-zero value in  order  to  enforce
//   malloc failure as soon as _alloc_counter_total increases above value of
//   this variable. This value has no effect if  _use_alloc_counter  is  not
//   set.
extern bool _force_malloc_failure;
extern ae_int_t _malloc_failure_after;

// debug functions (must be turned on by preprocessor definitions):
// * flushconsole(), fluches console
// * ae_debugrng(), returns random number generated with high-quality random numbers generator
// * ae_set_seed(), sets seed of the debug RNG (NON-THREAD-SAFE!!!)
// * ae_get_seed(), returns two seed values of the debug RNG (NON-THREAD-SAFE!!!)
#define flushconsole(s) fflush(stdout)

// Linear Algebra
#define ALGLIB_INTERCEPTS_ABLAS
void _ialglib_vzero(ae_int_t n, double *p, ae_int_t stride);
void _ialglib_vzero_complex(ae_int_t n, ae_complex *p, ae_int_t stride);
void _ialglib_vcopy(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb);
void _ialglib_vcopy_complex(ae_int_t n, const ae_complex *a, ae_int_t stridea, double *b, ae_int_t strideb, const char *conj);
void _ialglib_vcopy_dcomplex(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb, const char *conj);
void _ialglib_mcopyblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_int_t stride, double *b);
void _ialglib_mcopyunblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, double *b, ae_int_t stride);
void _ialglib_mcopyblock_complex(ae_int_t m, ae_int_t n, const ae_complex *a, ae_int_t op, ae_int_t stride, double *b);
void _ialglib_mcopyunblock_complex(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_complex *b, ae_int_t stride);

bool _ialglib_i_rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc);
bool _ialglib_i_cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, ae_matrix *c, ae_int_t ic, ae_int_t jc);
bool _ialglib_i_cmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_rmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_cmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_rmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, bool isupper, bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
bool _ialglib_i_cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool _ialglib_i_rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, bool isupper);
bool _ialglib_i_cmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
bool _ialglib_i_rmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
bool _ialglib_i_rmatrixgerf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, double alpha, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
} // end of namespace alglib_impl

namespace alglib {
typedef alglib_impl::ae_int_t ae_int_t;

// Exception class and exception handling macros.
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
#   define EndPoll		} catch(...) { alglib_impl::ae_state_clear(); throw; } }
#else
// Exception-free code.
#   define ThrowErrorMsg(X)	set_error_msg(); return X
#   if AE_THREADING != AE_SERIAL_UNSAFE
#      error Exception-free mode is thread-unsafe; define AE_THREADING = AE_SERIAL_UNSAFE to prove that you know it
#   endif
#   define BegPoll		{
#   define EndPoll		}
// Set the error flag and (optionally) the error message
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
#define EndD
#define AndD		,
#define DecVal(X)	X(Obj.X)
#define DecVar(X)	X(&Obj.X)
#define DecComplex(X)	X(*(complex *)&Obj.X)
#define ConstT(T, Val)	(const_cast<alglib_impl::T *>((Val).c_ptr()))
#define ComplexOf(Val)	(*reinterpret_cast<complex *>(&(Val)))

// The triple-layering -- also present in the distribution version of ALGLIB --
// is the cost of trying to make a C++ type into a C structure, and to then turn around and wrap it back into a C++ type.
// It's better to just lose the C code (which is, effectively, what the alglib_impl namespace is and encapsulates)
// and write it all in C++, instead.
// I compromised by re-wrapping the C++ wrapper code into the DefClass() and DecClass() pseudo-templates
// and then removing the need for and calls to malloc()/free() that were originally in the member functions comprising DefClass().
// This CAN be quashed down to 1 layer ...
// by just writing the member functions directly in the corresponding types currently contained in the alglib_impl namespace,
// and using this, instead.
// However, it requires rewriting all the code involving smart pointers and share pools as templates
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
   ~Type(); \
   alglib_impl::Type *c_ptr(); \
   alglib_impl::Type *c_ptr() const; \
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
} \
Type::~Type() { alglib_impl::Type##_free(&Obj, false); } \
alglib_impl::Type *Type::c_ptr() { return &Obj; } \
alglib_impl::Type *Type::c_ptr() const { return const_cast<alglib_impl::Type *>(&Obj); } \

// Complex number with double precision.
struct complex {
   complex();
   complex(const double &X);
   complex(const double &X, const double &Y);
   complex(const complex &Z);
   complex &operator=(const double &A); complex &operator=(const complex &A);
   complex &operator+=(const double &A); complex &operator+=(const complex &A);
   complex &operator-=(const double &A); complex &operator-=(const complex &A);
   complex &operator*=(const double &A); complex &operator*=(const complex &A);
   complex &operator/=(const double &A); complex &operator/=(const complex &A);
   alglib_impl::ae_complex *c_ptr(); const alglib_impl::ae_complex *c_ptr() const;
#if !defined AE_NO_EXCEPTIONS
   std::string tostring(int dps) const;
#endif
   double x, y;
};

const complex operator/(const complex &A, const complex &B);
bool operator==(const complex &A, const complex &B);
bool operator!=(const complex &A, const complex &B);
const complex operator+(const complex &A);
const complex operator-(const complex &A);
const complex operator+(const complex &A, const complex &B);
const complex operator-(const complex &A, const complex &B);
const complex operator*(const complex &A, const complex &B);
const complex operator/(const complex &A, const complex &B);
const complex operator+(const complex &A, const double &B);
const complex operator-(const complex &A, const double &B);
const complex operator*(const complex &A, const double &B);
const complex operator/(const complex &A, const double &B);
const complex operator+(const double &A, const complex &B);
const complex operator-(const double &A, const complex &B);
const complex operator*(const double &A, const complex &B);
const complex operator/(const double &A, const complex &B);
double abscomplex(const complex &A);
complex conj(const complex &A);
complex csqr(const complex &A);

// Level 1 BLAS functions
//
// NOTES:
// * destination and source should NOT overlap
// * stride is assumed to be positive, but it is not
//   assert'ed within function
// * conj_src parameter specifies whether complex source is conjugated
//   before processing or not. Pass string which starts with 'N' or 'n'
//   ("No conj", for example) to use unmodified parameter. All other
//   values will result in conjugation of input, but it is recommended
//   to use "Conj" in such cases.
double vdotproduct(const double *v0, ae_int_t stride0, const double *v1, ae_int_t stride1, ae_int_t n);
double vdotproduct(const double *v1, const double *v2, ae_int_t N);

complex vdotproduct(const complex *v0, ae_int_t stride0, const char *conj0, const complex *v1, ae_int_t stride1, const char *conj1, ae_int_t n);
complex vdotproduct(const complex *v1, const complex *v2, ae_int_t N);

void vmove(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void vmove(double *vdst, const double *vsrc, ae_int_t N);

void vmove(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void vmove(complex *vdst, const complex *vsrc, ae_int_t N);

void vmoveneg(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void vmoveneg(double *vdst, const double *vsrc, ae_int_t N);

void vmoveneg(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void vmoveneg(complex *vdst, const complex *vsrc, ae_int_t N);

void vmove(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha);
void vmove(double *vdst, const double *vsrc, ae_int_t N, double alpha);

void vmove(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha);
void vmove(complex *vdst, const complex *vsrc, ae_int_t N, double alpha);

void vmove(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, complex alpha);
void vmove(complex *vdst, const complex *vsrc, ae_int_t N, complex alpha);

void vadd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void vadd(double *vdst, const double *vsrc, ae_int_t N);

void vadd(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void vadd(complex *vdst, const complex *vsrc, ae_int_t N);

void vadd(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha);
void vadd(double *vdst, const double *vsrc, ae_int_t N, double alpha);

void vadd(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha);
void vadd(complex *vdst, const complex *vsrc, ae_int_t N, double alpha);

void vadd(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, complex alpha);
void vadd(complex *vdst, const complex *vsrc, ae_int_t N, complex alpha);

void vsub(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n);
void vsub(double *vdst, const double *vsrc, ae_int_t N);

void vsub(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n);
void vsub(complex *vdst, const complex *vsrc, ae_int_t N);

void vsub(double *vdst, ae_int_t stride_dst, const double *vsrc, ae_int_t stride_src, ae_int_t n, double alpha);
void vsub(double *vdst, const double *vsrc, ae_int_t N, double alpha);

void vsub(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, double alpha);
void vsub(complex *vdst, const complex *vsrc, ae_int_t N, double alpha);

void vsub(complex *vdst, ae_int_t stride_dst, const complex *vsrc, ae_int_t stride_src, const char *conj_src, ae_int_t n, complex alpha);
void vsub(complex *vdst, const complex *vsrc, ae_int_t N, complex alpha);

void vmul(double *vdst, ae_int_t stride_dst, ae_int_t n, double alpha);
void vmul(double *vdst, ae_int_t N, double alpha);

void vmul(complex *vdst, ae_int_t stride_dst, ae_int_t n, double alpha);
void vmul(complex *vdst, ae_int_t N, double alpha);

void vmul(complex *vdst, ae_int_t stride_dst, ae_int_t n, complex alpha);
void vmul(complex *vdst, ae_int_t N, complex alpha);

// xparams type and several predefined constants
typedef alglib_impl::ae_uint64_t xparams;
extern const xparams xdefault, serial, parallel;

// Threading functions
// nworkers can be 1, 2, ... ; or 0 for auto; or -1/-2/... for all except for one/two/...
void setnworkers(ae_int_t nworkers);

// sets global threading settings to serial or parallel
void setglobalthreading(const xparams settings);

// nworkers can be 1, 2, ... ; or 0 for auto; or -1/-2/... for all except for one/two/...
ae_int_t getnworkers();

// internal functions used by TestX.cpp, interfaces for functions present
// in commercial ALGLIB but lacking in free edition.
ae_int_t _ae_cores_count();
void _ae_set_global_threading(alglib_impl::ae_uint64_t flg_value);
alglib_impl::ae_uint64_t _ae_get_global_threading();

// 1- and 2-dimensional arrays
struct ae_vector_wrapper {
// Creates object attached to external ae_vector structure.
// NOTE: this function also checks that source ae_vector* has
//       required datatype. An exception is generated otherwise.
   ae_vector_wrapper(alglib_impl::ae_vector *e_ptr, alglib_impl::ae_datatype datatype);
// Creates zero-size vector of specific datatype
   ae_vector_wrapper(alglib_impl::ae_datatype datatype);
// Creates a copy of another vector (can be reference to one of the derived classes)
// NOTE: this function also checks that source ae_vector* has
//       required datatype. An exception is generated otherwise.
   ae_vector_wrapper(const ae_vector_wrapper &rhs, alglib_impl::ae_datatype datatype);
// Well, it is destructor...
   virtual ~ae_vector_wrapper();
// For wrapper object allocated with allocate_own() this function
// changes length, completely dropping previous contents.
//
// It does not work (throws exception) for frozen proxy objects.
   void setlength(ae_int_t iLen);
// Element count
   ae_int_t length() const;
// Access to internal C-structure used by C-core.
// Not intended for external use.
   const alglib_impl::ae_vector *c_ptr() const;
   alglib_impl::ae_vector *c_ptr();
private:
   ae_vector_wrapper();
   ae_vector_wrapper(const ae_vector_wrapper &rhs);
   const ae_vector_wrapper &operator=(const ae_vector_wrapper &rhs);
protected:
#if !defined AE_NO_EXCEPTIONS
// Copies array given by string into current object. Additional
// parameter DATATYPE contains information about type of the data
// in S and type of the array to create.
//
// NOTE: this function is not supported in exception-free mode.
   ae_vector_wrapper(const char *s, alglib_impl::ae_datatype datatype);
#endif
// This function attaches wrapper object to external x_vector structure;
// "frozen proxy" mode is activated (you can read/write, but can not reallocate
// and do not own memory of the vector).
// NOTE: initial state of wrapper object is assumed to be initialized;
//       all previously allocated memory is properly deallocated.
// NOTE: x_vector structure pointed by new_ptr is used only once; after
//       we fetch pointer to memory and its size, this structure is ignored
//       and not referenced anymore. So, you can pass pointers to temporary
//       x-structures which are deallocated immediately after you call attach_to()
// NOTE: state structure is used for error reporting purposes (longjmp on errors).
   void attach_to(alglib_impl::x_vector *new_ptr);
// Assigns RHS to current object. Returns *this.
// It has several branches depending on target object status:
// * in case it is proxy object, data are copied into memory pointed by
//   proxy. Function checks that source has exactly same size as target
//   (exception is thrown on failure).
// * in case it is non-proxy object, data allocated by object are cleared
//   and a copy of RHS is created in target.
// NOTE: this function correctly handles assignments of the object to itself.
   const ae_vector_wrapper &assign(const ae_vector_wrapper &rhs);
// Pointer to ae_vector structure:
// * This == &Obj means that wrapper object owns ae_vector structure and
//   is responsible for proper deallocation of its memory
// * This != &Obj means that wrapper object works with someone's other
//   ae_vector record and is not responsible for its memory; in this case
//   Obj is assumed to be uninitialized.
   alglib_impl::ae_vector *This;
// Inner ae_vector record.
// Ignored for This != &Obj.
   alglib_impl::ae_vector Obj;
// Whether this wrapper object is frozen proxy (you may read array, may
// modify its value, but can not deallocate its memory or resize it) or not.
//
// If is_frozen_proxy == true and if:
// * This == &Obj, it means that wrapper works with its own ae_vector
//   structure, but this structure points to externally allocated memory.
//   This memory is NOT owned by ae_vector object.
// * This != &Obj, it means that wrapper works with externally allocated
//   and managed ae_vector structure. Both memory pointed by ae_vector and
//   ae_vector structure itself are not owned by wrapper object.
   bool is_frozen_proxy;
};

struct boolean_1d_array: public ae_vector_wrapper {
   boolean_1d_array();
   boolean_1d_array(const boolean_1d_array &rhs);
   boolean_1d_array(alglib_impl::ae_vector *p);
   const boolean_1d_array &operator=(const boolean_1d_array &rhs);
   virtual ~boolean_1d_array();
   const bool &operator()(ae_int_t i) const;
   bool &operator()(ae_int_t i);
   const bool &operator[](ae_int_t i) const;
   bool &operator[](ae_int_t i);
// This function allocates array[iLen] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t iLen, const bool *pContent);
// This function returns pointer to internal memory
   bool *getcontent();
   const bool *getcontent() const;
#if !defined AE_NO_EXCEPTIONS
   boolean_1d_array(const char *s);
   std::string tostring() const;
#endif
};

struct integer_1d_array: public ae_vector_wrapper {
   integer_1d_array();
   integer_1d_array(const integer_1d_array &rhs);
   integer_1d_array(alglib_impl::ae_vector *p);
   const integer_1d_array &operator=(const integer_1d_array &rhs);
   virtual ~integer_1d_array();
   const ae_int_t &operator()(ae_int_t i) const;
   ae_int_t &operator()(ae_int_t i);
   const ae_int_t &operator[](ae_int_t i) const;
   ae_int_t &operator[](ae_int_t i);
// This function allocates array[iLen] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t iLen, const ae_int_t *pContent);
// This function returns pointer to internal memory
   ae_int_t *getcontent();
   const ae_int_t *getcontent() const;
#if !defined AE_NO_EXCEPTIONS
   integer_1d_array(const char *s);
   std::string tostring() const;
#endif
};

struct real_1d_array: public ae_vector_wrapper {
   real_1d_array();
   real_1d_array(const real_1d_array &rhs);
   real_1d_array(alglib_impl::ae_vector *p);
   const real_1d_array &operator=(const real_1d_array &rhs);
   virtual ~real_1d_array();
   const double &operator()(ae_int_t i) const;
   double &operator()(ae_int_t i);
   const double &operator[](ae_int_t i) const;
   double &operator[](ae_int_t i);
// This function allocates array[iLen] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t iLen, const double *pContent);
// This function attaches array to memory pointed by pContent.
// No own memory is allocated, no copying of data is performed,
// so pContent pointer should be valid as long as we work with
// array.
//
// After you attach array object to external memory, it becomes
// "frozen": it is possible to read/write array elements, but
// it is not allowed to resize it (no setlength() calls).
   void attach_to_ptr(ae_int_t iLen, double *pContent);
// This function returns pointer to internal memory
   double *getcontent();
   const double *getcontent() const;
#if !defined AE_NO_EXCEPTIONS
   real_1d_array(const char *s);
   std::string tostring(int dps) const;
#endif
};

struct complex_1d_array: public ae_vector_wrapper {
   complex_1d_array();
   complex_1d_array(const complex_1d_array &rhs);
   complex_1d_array(alglib_impl::ae_vector *p);
   const complex_1d_array &operator=(const complex_1d_array &rhs);
   virtual ~complex_1d_array();
   const complex &operator()(ae_int_t i) const;
   complex &operator()(ae_int_t i);
   const complex &operator[](ae_int_t i) const;
   complex &operator[](ae_int_t i);
// This function allocates array[iLen] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t iLen, const complex *pContent);
   complex *getcontent();
   const complex *getcontent() const;
#if !defined AE_NO_EXCEPTIONS
   complex_1d_array(const char *s);
   std::string tostring(int dps) const;
#endif
};

struct ae_matrix_wrapper {
// Creates object attached to external ae_vector structure, with additional
// check for matching datatypes (e_ptr->datatype == datatype is required).
   ae_matrix_wrapper(alglib_impl::ae_matrix *e_ptr, alglib_impl::ae_datatype datatype);
// Creates zero-sized matrix of specified datatype.
   ae_matrix_wrapper(alglib_impl::ae_datatype datatype);
// Creates copy of rhs, with additional check for matching datatypes
// (rhs.datatype == datatype is required).
   ae_matrix_wrapper(const ae_matrix_wrapper &rhs, alglib_impl::ae_datatype datatype);
// Destructor
   virtual ~ae_matrix_wrapper();
   void setlength(ae_int_t rows, ae_int_t cols);
   ae_int_t rows() const;
   ae_int_t cols() const;
   bool isempty() const;
   ae_int_t getstride() const;
   const alglib_impl::ae_matrix *c_ptr() const;
   alglib_impl::ae_matrix *c_ptr();
private:
   ae_matrix_wrapper();
   ae_matrix_wrapper(const ae_matrix_wrapper &rhs);
   const ae_matrix_wrapper &operator=(const ae_matrix_wrapper &rhs);
protected:
#if !defined AE_NO_EXCEPTIONS
// Copies array given by string into current object. Additional
// parameter DATATYPE contains information about type of the data
// in S and type of the array to create.
//
// Current object is considered empty (this function should be
// called from copy constructor).
   ae_matrix_wrapper(const char *s, alglib_impl::ae_datatype datatype);
#endif
// This function attaches wrapper object to external x_vector structure;
// "frozen proxy" mode is activated (you can read/write, but can not reallocate
// and do not own memory of the vector).
// NOTE: initial state of wrapper object is assumed to be initialized;
//       all previously allocated memory is properly deallocated.
// NOTE: x_vector structure pointed by new_ptr is used only once; after
//       we fetch pointer to memory and its size, this structure is ignored
//       and not referenced anymore. So, you can pass pointers to temporary
//       x-structures which are deallocated immediately after you call attach_to()
// NOTE: state structure is used for error-handling (a longjmp is performed
//       on allocation error). All previously allocated memory is correctly
//       freed on error.
   void attach_to(alglib_impl::x_matrix *new_ptr);
// Assigns RHS to current object.
// It has several branches depending on target object status:
// * in case it is proxy object, data are copied into memory pointed by
//   proxy. Function checks that source has exactly same size as target
//   (exception is thrown on failure).
// * in case it is non-proxy object, data allocated by object are cleared
//   and a copy of RHS is created in target.
// NOTE: this function correctly handles assignments of the object to itself.
   const ae_matrix_wrapper &assign(const ae_matrix_wrapper &rhs);
// Pointer to ae_matrix structure:
// * This == &Obj means that wrapper object owns ae_matrix structure and
//   is responsible for proper deallocation of its memory
// * This != &Obj means that wrapper object works with someone's other
//   ae_matrix record and is not responsible for its memory; in this case
//   Obj is assumed to be uninitialized.
   alglib_impl::ae_matrix *This;
// Inner ae_matrix record.
// Ignored for This != &Obj.
   alglib_impl::ae_matrix Obj;
// Whether this wrapper object is frozen proxy (you may read array, may
// modify its value, but can not deallocate its memory or resize it) or not.
// If is_frozen_proxy == true and if:
// * This == &Obj, it means that wrapper works with its own ae_vector
//   structure, but this structure points to externally allocated memory.
//   This memory is NOT owned by ae_vector object.
// * This != &Obj, it means that wrapper works with externally allocated
//   and managed ae_vector structure. Both memory pointed by ae_vector and
//   ae_vector structure itself are not owned by wrapper object.
   bool is_frozen_proxy;
};

struct boolean_2d_array: public ae_matrix_wrapper {
   boolean_2d_array();
   boolean_2d_array(const boolean_2d_array &rhs);
   boolean_2d_array(alglib_impl::ae_matrix *p);
   virtual ~boolean_2d_array();
   const boolean_2d_array &operator=(const boolean_2d_array &rhs);
   const bool &operator()(ae_int_t i, ae_int_t j) const;
   bool &operator()(ae_int_t i, ae_int_t j);
   const bool *operator[](ae_int_t i) const;
   bool *operator[](ae_int_t i);
// This function allocates array[irows,icols] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t irows, ae_int_t icols, const bool *pContent);
#if !defined AE_NO_EXCEPTIONS
   boolean_2d_array(const char *s);
   std::string tostring() const;
#endif
};

struct integer_2d_array: public ae_matrix_wrapper {
   integer_2d_array();
   integer_2d_array(const integer_2d_array &rhs);
   integer_2d_array(alglib_impl::ae_matrix *p);
   virtual ~integer_2d_array();
   const integer_2d_array &operator=(const integer_2d_array &rhs);
   const ae_int_t &operator()(ae_int_t i, ae_int_t j) const;
   ae_int_t &operator()(ae_int_t i, ae_int_t j);
   const ae_int_t *operator[](ae_int_t i) const;
   ae_int_t *operator[](ae_int_t i);
// This function allocates array[irows,icols] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t irows, ae_int_t icols, const ae_int_t *pContent);
#if !defined AE_NO_EXCEPTIONS
   integer_2d_array(const char *s);
   std::string tostring() const;
#endif
};

struct real_2d_array: public ae_matrix_wrapper {
   real_2d_array();
   real_2d_array(const real_2d_array &rhs);
   real_2d_array(alglib_impl::ae_matrix *p);
   virtual ~real_2d_array();
   const real_2d_array &operator=(const real_2d_array &rhs);
   const double &operator()(ae_int_t i, ae_int_t j) const;
   double &operator()(ae_int_t i, ae_int_t j);
   const double *operator[](ae_int_t i) const;
   double *operator[](ae_int_t i);
// This function allocates array[irows,icols] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t irows, ae_int_t icols, const double *pContent);
// This function attaches array to memory pointed by pContent:
// * only minor amount of own memory is allocated - O(irows) bytes to
//   store precomputed pointers; but no costly copying of O(rows*cols)
//   data is performed.
// * pContent pointer should be valid as long as we work with array
// After you attach array object to external memory, it becomes
// "frozen": it is possible to read/write array elements, but
// it is not allowed to resize it (no setlength() calls).
   void attach_to_ptr(ae_int_t irows, ae_int_t icols, double *pContent);
#if !defined AE_NO_EXCEPTIONS
   real_2d_array(const char *s);
   std::string tostring(int dps) const;
#endif
};

struct complex_2d_array: public ae_matrix_wrapper {
   complex_2d_array();
   complex_2d_array(const complex_2d_array &rhs);
   complex_2d_array(alglib_impl::ae_matrix *p);
   virtual ~complex_2d_array();
   const complex_2d_array &operator=(const complex_2d_array &rhs);
   const complex &operator()(ae_int_t i, ae_int_t j) const;
   complex &operator()(ae_int_t i, ae_int_t j);
   const complex *operator[](ae_int_t i) const;
   complex *operator[](ae_int_t i);
// This function allocates array[irows,icols] and copies data
// pointed by pContent to its memory. Completely independent
// copy of data is created.
   void setcontent(ae_int_t irows, ae_int_t icols, const complex *pContent);
#if !defined AE_NO_EXCEPTIONS
   complex_2d_array(const char *s);
   std::string tostring(int dps) const;
#endif
};

// CSV operations: reading CSV file to real matrix.
//
// This function reads CSV  file  and  stores  its  contents  to  double
// precision 2D array. Format of the data file must conform to RFC  4180
// specification, with additional notes:
// * file size should be less than 2GB
// * ASCI encoding, UTF-8 without BOM (in header names) are supported
// * any character (comma/tab/space) may be used as field separator,  as
//   long as it is distinct from one used for decimal point
// * multiple subsequent field separators (say, two  spaces) are treated
//   as MULTIPLE separators, not one big separator
// * both comma and full stop may be used as decimal point. Parser  will
//   automatically determine specific character being used.  Both  fixed
//   and exponential number formats are  allowed.   Thousand  separators
//   are NOT allowed.
// * line may end with \n (Unix style) or \r\n (Windows  style),  parser
//   will automatically adapt to chosen convention
// * escaped fields (ones in double quotes) are not supported
//
// Inputs:
//     filename        relative/absolute path
//     separator       character used to separate fields.  May  be  ' ',
//                     ',', '\t'. Other separators are possible too.
//     flags           several values combined with bitwise OR:
//                     * CSV_SKIP_HEADERS -  if present, first row
//                       contains headers  and  will  be  skipped.   Its
//                       contents is used to determine fields count, and
//                       that's all.
//                     If no flags are specified, default value 0x0  (or
//                     CSV_DEFAULT, which is same) should be used.
//
// Outputs:
//     out             2D matrix, CSV file parsed with atof()
//
// HANDLING OF SPECIAL CASES:
// * file does not exist - ap_error exception is thrown
// * empty file - empty array is returned (no exception)
// * skip_first_row=true, only one row in file - empty array is returned
// * field contents is not recognized by atof() - field value is replaced
//   by 0.0
#if !defined AE_NO_EXCEPTIONS
void read_csv(const char *filename, char separator, int flags, real_2d_array &out);
#endif

// Constants and functions introduced for compatibility with AlgoPascal
extern const double machineepsilon, maxrealnumber, minrealnumber;
static const int CSV_DEFAULT = 0x0;
static const int CSV_SKIP_HEADERS = 0x1;

int sign(double x);
double randomreal();
ae_int_t randominteger(ae_int_t maxv);
int RoundZ(double x);
int TruncZ(double x);
int FloorZ(double x);
int CeilZ(double x);
double pi();
double sqr(double x);
int maxint(int m1, int m2);
int minint(int m1, int m2);
double maxreal(double m1, double m2);
double minreal(double m1, double m2);
bool isposinf(double x);
bool isneginf(double x);
} // end of namespace alglib

#endif // OnceOnly
