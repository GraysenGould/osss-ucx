/* For license: see LICENSE file at top-level */

/**
 * @file putget_complex.h
 * @brief Header file defining OpenSHMEM complex type put/get operations
 *
 * This file contains macro definitions for implementing put and get operations
 * for complex floating point types in OpenSHMEM. These macros define
 * type-specific variants of put, get, iput, iget, p and g operations for
 * complex types.
 */

#ifndef _SHMEM_PUTGET_COMPLEX_H
#define _SHMEM_PUTGET_COMPLEX_H 1

/**
 * @brief Macro to declare typed put operations with context for complex types
 *
 * @param _name Type name suffix for the operation
 * @param _type Actual complex C type for the operation
 *
 * Declares functions that perform put, put_nbi and iput operations for a
 * specific complex data type, using a communication context.
 */
#define DECL_TYPED_CTX_PUT_PROTO(_name, _type)                                 \
  void shmem_ctx_##_name##_put(shmem_ctx_t ctx, _type *dest, const _type *src, \
                               size_t nelems, int pe);                         \
  void shmem_ctx_##_name##_put_nbi(shmem_ctx_t ctx, _type *dest,               \
                                   const _type *src, size_t nelems, int pe);   \
  void shmem_ctx_##_name##_iput(shmem_ctx_t ctx, _type *dest,                  \
                                const _type *src, ptrdiff_t tst,               \
                                ptrdiff_t sst, size_t nelems, int pe);         \
  void shmem_ctx_##_name##_ibput(shmem_ctx_t ctx, _type *dest,                 \
                                const _type *src, ptrdiff_t tst,               \
                                ptrdiff_t sst, size_t bsize, size_t nblocks,   \
                                int pe);

DECL_TYPED_CTX_PUT_PROTO(complexf, COMPLEXIFY(float))
DECL_TYPED_CTX_PUT_PROTO(complexd, COMPLEXIFY(double))

/**
 * @brief Macro to declare typed get operations with context for complex types
 *
 * @param _name Type name suffix for the operation
 * @param _type Actual complex C type for the operation
 *
 * Declares functions that perform get, get_nbi and iget operations for a
 * specific complex data type, using a communication context.
 */
#define DECL_TYPED_CTX_GET_PROTO(_name, _type)                                 \
  void shmem_ctx_##_name##_get(shmem_ctx_t ctx, _type *dest, const _type *src, \
                               size_t nelems, int pe);                         \
  void shmem_ctx_##_name##_get_nbi(shmem_ctx_t ctx, _type *dest,               \
                                   const _type *src, size_t nelems, int pe);   \
  void shmem_ctx_##_name##_iget(shmem_ctx_t ctx, _type *dest,                  \
                                const _type *src, ptrdiff_t tst,               \
                                ptrdiff_t sst, size_t nelems, int pe);         \
  void shmem_ctx_##_name##_ibget(shmem_ctx_t ctx, _type *dest,                 \
                                const _type *src, ptrdiff_t tst,               \
                                ptrdiff_t sst, size_t bsize, size_t nblocks,   \
                                int pe);

DECL_TYPED_CTX_GET_PROTO(complexf, COMPLEXIFY(float))
DECL_TYPED_CTX_GET_PROTO(complexd, COMPLEXIFY(double))

/**
 * @brief Macro to declare typed put scalar operations with context for complex
 * types
 *
 * @param _name Type name suffix for the operation
 * @param _type Actual complex C type for the operation
 *
 * Declares functions that perform atomic put operations for a specific complex
 * data type, using a communication context.
 */
#define DECL_TYPED_CTX_P_PROTO(_name, _type)                                   \
  void shmem_ctx_##_name##_p(shmem_ctx_t ctx, _type *addr, _type val, int pe);

DECL_TYPED_CTX_P_PROTO(complexf, COMPLEXIFY(float))
DECL_TYPED_CTX_P_PROTO(complexd, COMPLEXIFY(double))

/**
 * @brief Macro to declare typed get scalar operations with context for complex
 * types
 *
 * @param _name Type name suffix for the operation
 * @param _type Actual complex C type for the operation
 *
 * Declares functions that perform atomic get operations for a specific complex
 * data type, using a communication context.
 */
#define DECL_TYPED_CTX_G_PROTO(_name, _type)                                   \
  _type shmem_ctx_##_name##_g(shmem_ctx_t ctx, const _type *addr, int pe);

DECL_TYPED_CTX_G_PROTO(complexf, COMPLEXIFY(float))
DECL_TYPED_CTX_G_PROTO(complexd, COMPLEXIFY(double))

#endif /* ! _SHMEM_PUTGET_COMPLEX_H */
