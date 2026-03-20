/**
 * @file putget.h
 * @brief OpenSHMEM put/get operation declarations and definitions
 *
 * For license: see LICENSE file at top-level
 */

#ifndef _SHMEM_PUTGET_H
#define _SHMEM_PUTGET_H 1

#include "shmemu.h"
#include "shmemc.h"
#include "shmem_mutex.h"

/**
 * @brief Macro to define a typed put operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_PUT(_name, _type)                                      \
  void shmem_ctx_##_name##_put(shmem_ctx_t ctx, _type *dest, const _type *src, \
                               size_t nelems, int pe) {                        \
    const size_t nb = sizeof(_type) * nelems;                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put(ctx, dest, src, nb, pe));            \
  }

/**
 * @brief Macro to define a typed get operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_GET(_name, _type)                                      \
  void shmem_ctx_##_name##_get(shmem_ctx_t ctx, _type *dest, const _type *src, \
                               size_t nelems, int pe) {                        \
    const size_t nb = sizeof(_type) * nelems;                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 3);                                            \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get(ctx, dest, src, nb, pe));            \
  }

/**
 * @brief Macro to define a sized put operation with a context
 * @param _size Size in bits
 */
#define SHMEM_CTX_SIZED_PUT(_size)                                             \
  void shmem_ctx_put##_size(shmem_ctx_t ctx, void *dest, const void *src,      \
                            size_t nelems, int pe) {                           \
    const size_t nb = BITS2BYTES(_size) * nelems;                              \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put(ctx, dest, src, nb, pe));            \
  }

/**
 * @brief Macro to define a sized get operation with a context
 * @param _size Size in bits
 */
#define SHMEM_CTX_SIZED_GET(_size)                                             \
  void shmem_ctx_get##_size(shmem_ctx_t ctx, void *dest, const void *src,      \
                            size_t nelems, int pe) {                           \
    const size_t nb = BITS2BYTES(_size) * nelems;                              \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 3);                                            \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get(ctx, dest, src, nb, pe));            \
  }

/**
 * @brief Macro to define a memory put operation with a context
 */
#define SHMEM_CTX_PUTMEM()                                                     \
  void shmem_ctx_putmem(shmem_ctx_t ctx, void *dest, const void *src,          \
                        size_t nelems, int pe) {                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put(ctx, dest, src, nelems, pe));        \
  }

/**
 * @brief Macro to define a memory get operation with a context
 */
#define SHMEM_CTX_GETMEM()                                                     \
  void shmem_ctx_getmem(shmem_ctx_t ctx, void *dest, const void *src,          \
                        size_t nelems, int pe) {                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 3);                                            \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get(ctx, dest, src, nelems, pe));        \
  }

/**
 * @brief Macro to define a typed strided put operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_IPUT(_name, _type)                                     \
  void shmem_ctx_##_name##_iput(shmem_ctx_t ctx, _type *target,                \
                                const _type *source, ptrdiff_t tst,            \
                                ptrdiff_t sst, size_t nelems, int pe) {        \
    size_t ti = 0, si = 0;                                                     \
    size_t i;                                                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 7);                                          \
    SHMEMU_CHECK_SYMMETRIC(target, 2);                                         \
                                                                               \
    logger(LOG_RMA,                                                            \
           "%s(ctx=%lu, dest=%p, src=%p, "                                     \
           "tst=%lu, sst=%lu, nelems=%lu, pe=%d)",                             \
           __func__, shmemc_context_id(ctx), target, source, tst, sst, nelems, \
           pe);                                                                \
                                                                               \
    for (i = 0; i < nelems; ++i) {                                             \
      shmem_ctx_##_name##_put(ctx, &(target[ti]), &(source[si]), 1, pe);       \
      ti += tst;                                                               \
      si += sst;                                                               \
    }                                                                          \
  }

/**
 * @brief Macro to define a typed strided get operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_IGET(_name, _type)                                     \
  void shmem_ctx_##_name##_iget(shmem_ctx_t ctx, _type *target,                \
                                const _type *source, ptrdiff_t tst,            \
                                ptrdiff_t sst, size_t nelems, int pe) {        \
    size_t ti = 0, si = 0;                                                     \
    size_t i;                                                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 7);                                          \
    SHMEMU_CHECK_SYMMETRIC(source, 3);                                         \
                                                                               \
    logger(LOG_RMA,                                                            \
           "%s(ctx=%lu, dest=%p, src=%p, "                                     \
           "tst=%lu, sst=%lu, nelems=%lu, pe=%d)",                             \
           __func__, shmemc_context_id(ctx), target, source, tst, sst, nelems, \
           pe);                                                                \
                                                                               \
    for (i = 0; i < nelems; ++i) {                                             \
      shmem_ctx_##_name##_get(ctx, &(target[ti]), &(source[si]), 1, pe);       \
      ti += tst;                                                               \
      si += sst;                                                               \
    }                                                                          \
  }

/**
 * @brief Macro to define a typed strided block get operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_IBGET(_name, _type)                                    \
  void shmem_ctx_##_name##_ibget(shmem_ctx_t ctx, _type *target,                \
                                const _type *source, ptrdiff_t tst,            \
                                ptrdiff_t sst, size_t bsize, size_t nblocks,   \
                                int pe) {                                      \
    size_t ti = 0, si = 0;                                                     \
    size_t i;                                                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_COUNTING_MUTIPLE(tst, bsize);                                 \
    SHMEMU_CHECK_COUNTING_MUTIPLE(sst, bsize);                                 \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 7);                                          \
    SHMEMU_CHECK_SYMMETRIC(source, 3);                                         \
                                                                               \
    logger(LOG_RMA,                                                            \
           "%s(ctx=%lu, dest=%p, src=%p, "                                     \
           "tst=%lu, sst=%lu, bsize=%lu, nblocks=%ld, pe=%d)",                 \
           __func__, shmemc_context_id(ctx), target, source, tst, sst, bsize,  \
           nblocks, pe);                                                       \
                                                                               \
    for (i = 0; i < nblocks; ++i) {                                            \
      shmem_ctx_##_name##_get(ctx, &(target[ti]), &(source[si]), bsize, pe);   \
      ti += tst;                                                               \
      si += sst;                                                               \
    }                                                                          \
  }

/**
 * @brief Macro to define a sized strided put operation with a context
 * @param _size Size in bits
 */
#define SHMEM_CTX_SIZED_IPUT(_size)                                            \
  void shmem_ctx_iput##_size(shmem_ctx_t ctx, void *target,                    \
                             const void *source, ptrdiff_t tst, ptrdiff_t sst, \
                             size_t nelems, int pe) {                          \
    const size_t tst_nb = BITS2BYTES(_size) * tst;                             \
    const size_t sst_nb = BITS2BYTES(_size) * sst;                             \
    size_t ti = 0, si = 0;                                                     \
    size_t i;                                                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 7);                                          \
    SHMEMU_CHECK_SYMMETRIC(target, 2);                                         \
                                                                               \
    logger(LOG_RMA,                                                            \
           "%s(ctx=%lu, dest=%p, src=%p, "                                     \
           "tst=%lu, sst=%lu, nelems=%lu, pe=%d)",                             \
           __func__, shmemc_context_id(ctx), target, source, tst, sst, nelems, \
           pe);                                                                \
                                                                               \
    for (i = 0; i < nelems; ++i) {                                             \
      shmem_ctx_put##_size(ctx, (void *)&((char *)target)[ti],                 \
                           (void *)&((char *)source)[si], 1, pe);              \
      ti += tst_nb;                                                            \
      si += sst_nb;                                                            \
    }                                                                          \
  }

/**
 * @brief Macro to define a sized strided get operation with a context
 * @param _size Size in bits
 */
#define SHMEM_CTX_SIZED_IGET(_size)                                            \
  void shmem_ctx_iget##_size(shmem_ctx_t ctx, void *target,                    \
                             const void *source, ptrdiff_t tst, ptrdiff_t sst, \
                             size_t nelems, int pe) {                          \
    const size_t tst_nb = BITS2BYTES(_size) * tst;                             \
    const size_t sst_nb = BITS2BYTES(_size) * sst;                             \
    size_t ti = 0, si = 0;                                                     \
    size_t i;                                                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 7);                                          \
    SHMEMU_CHECK_SYMMETRIC(source, 3);                                         \
                                                                               \
    logger(LOG_RMA,                                                            \
           "%s(ctx=%lu, dest=%p, src=%p, "                                     \
           "tst=%lu, sst=%lu, nelems=%lu, pe=%d)",                             \
           __func__, shmemc_context_id(ctx), target, source, tst, sst, nelems, \
           pe);                                                                \
                                                                               \
    for (i = 0; i < nelems; ++i) {                                             \
      shmem_ctx_get##_size(ctx, (void *)&((char *)target)[ti],                 \
                           (void *)&((char *)source)[si], 1, pe);              \
      ti += tst_nb;                                                            \
      si += sst_nb;                                                            \
    }                                                                          \
  }

/**
 * @brief Macro to define a sized strided block get operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_SIZED_IBGET(_size)                                           \
  void shmem_ctx_ibget##_size(shmem_ctx_t ctx, void *target,                   \
                                const void *source, ptrdiff_t tst,             \
                                ptrdiff_t sst, size_t bsize, size_t nblocks,   \
                                int pe) {                                      \
    const size_t tst_nb = BITS2BYTES(_size) * tst;                             \
    const size_t sst_nb = BITS2BYTES(_size) * sst;                             \
    size_t ti = 0, si = 0;                                                     \
    size_t i;                                                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_COUNTING_MUTIPLE(tst, bsize);                                 \
    SHMEMU_CHECK_COUNTING_MUTIPLE(sst, bsize);                                 \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 7);                                          \
    SHMEMU_CHECK_SYMMETRIC(source, 3);                                         \
                                                                               \
    logger(LOG_RMA,                                                            \
           "%s(ctx=%lu, dest=%p, src=%p, "                                     \
           "tst=%lu, sst=%lu, bsize=%lu, nblocks=%ld, pe=%d)",                 \
           __func__, shmemc_context_id(ctx), target, source, tst, sst, bsize,  \
           nblocks, pe);                                                       \
                                                                               \
    for (i = 0; i < nblocks; ++i) {                                            \
      shmem_ctx_get##_size(ctx, (void *)&((char *)target)[ti],                 \
                           (void *)&((char *)source)[si], bsize, pe);          \
      ti += tst_nb;                                                            \
      si += sst_nb;                                                            \
    }                                                                          \
  }

/**
 * @brief Macro to define a typed non-blocking put operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_PUT_NBI(_name, _type)                                  \
  void shmem_ctx_##_name##_put_nbi(shmem_ctx_t ctx, _type *dest,               \
                                   const _type *src, size_t nelems, int pe) {  \
    const size_t nb = sizeof(_type) * nelems;                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put_nbi(ctx, dest, src, nb, pe));        \
  }

/**
 * @brief Macro to define a typed non-blocking get operation with a context
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_GET_NBI(_name, _type)                                  \
  void shmem_ctx_##_name##_get_nbi(shmem_ctx_t ctx, _type *dest,               \
                                   const _type *src, size_t nelems, int pe) {  \
    const size_t nb = sizeof(_type) * nelems;                                  \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 3);                                            \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get_nbi(ctx, dest, src, nb, pe));        \
  }

/**
 * @brief Macro to define a sized non-blocking put operation with a context
 * @param _size Size in bits
 */
#define SHMEM_CTX_SIZED_PUT_NBI(_size)                                         \
  void shmem_ctx_put##_size##_nbi(shmem_ctx_t ctx, void *dest,                 \
                                  const void *src, size_t nelems, int pe) {    \
    const size_t nb = BITS2BYTES(_size) * nelems;                              \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put_nbi(ctx, dest, src, nb, pe));        \
  }

/**
 * @brief Macro to define a sized non-blocking get operation with a context
 * @param _size Size in bits
 */
#define SHMEM_CTX_SIZED_GET_NBI(_size)                                         \
  void shmem_ctx_get##_size##_nbi(shmem_ctx_t ctx, void *dest,                 \
                                  const void *src, size_t nelems, int pe) {    \
    const size_t nb = BITS2BYTES(_size) * nelems;                              \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 3);                                            \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get_nbi(ctx, dest, src, nb, pe));        \
  }

/**
 * @brief Macro to define a memory non-blocking put operation with a context
 */
#define SHMEM_CTX_PUTMEM_NBI()                                                 \
  void shmem_ctx_putmem_nbi(shmem_ctx_t ctx, void *dest, const void *src,      \
                            size_t nelems, int pe) {                           \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put_nbi(ctx, dest, src, nelems, pe));    \
  }

/**
 * @brief Macro to define a memory non-blocking get operation with a context
 * @param _op Operation type
 */
#define SHMEM_CTX_GETMEM_NBI(_op)                                              \
  void shmem_ctx_getmem_nbi(shmem_ctx_t ctx, void *dest, const void *src,      \
                            size_t nelems, int pe) {                           \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 5);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 3);                                            \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, dest=%p, src=%p, nelems=%lu, pe=%d)",         \
           __func__, shmemc_context_id(ctx), dest, src, nelems, pe);           \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get_nbi(ctx, dest, src, nelems, pe));    \
  }

/**
 * @brief Macro to define a typed put operation with a context for single values
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_P(_name, _type)                                        \
  void shmem_ctx_##_name##_p(shmem_ctx_t ctx, _type *addr, _type val,          \
                             int pe) {                                         \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(addr, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, addr=%p, value=%lu, pe=%d)", __func__,        \
           shmemc_context_id(ctx), addr, val, pe);                             \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_put(ctx, addr, &val, sizeof(val), pe));  \
  }

/**
 * @brief Macro to define a typed get operation with a context for single values
 * @param _name Type name
 * @param _type C type
 */
#define SHMEM_CTX_TYPED_G(_name, _type)                                        \
  _type shmem_ctx_##_name##_g(shmem_ctx_t ctx, const _type *addr, int pe) {    \
    _type val;                                                                 \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 3);                                          \
    SHMEMU_CHECK_SYMMETRIC(addr, 2);                                           \
                                                                               \
    logger(LOG_RMA, "%s(ctx=%lu, addr=%p, pe=%d)", __func__,                   \
           shmemc_context_id(ctx), addr, pe);                                  \
                                                                               \
    SHMEMT_MUTEX_NOPROTECT(shmemc_ctx_get(ctx, &val, addr, sizeof(val), pe));  \
    return val;                                                                \
  }

/*
 * operating on implicit default context
 */

/**
 * @brief Macro to declare typed put operations
 * @param _name Type name
 * @param _type C type
 */
#define API_DECL_TYPED_PUT(_name, _type)                                       \
  void shmem_##_name##_put(_type *dest, const _type *src, size_t nelems,       \
                           int pe) {                                           \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_##_name##_put(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);         \
  }                                                                            \
  void shmem_##_name##_put_nbi(_type *dest, const _type *src, size_t nelems,   \
                               int pe) {                                       \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_##_name##_put_nbi(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);     \
  }                                                                            \
  void shmem_##_name##_iput(_type *dest, const _type *src, ptrdiff_t tst,      \
                            ptrdiff_t sst, size_t nelems, int pe) {            \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 6);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_##_name##_iput(SHMEM_CTX_DEFAULT, dest, src, tst, sst, nelems,   \
                             pe);                                              \
  }

/**
 * @brief Macro to declare typed get operations
 * @param _name Type name
 * @param _type C type
 */
#define API_DECL_TYPED_GET(_name, _type)                                       \
  void shmem_##_name##_get(_type *dest, const _type *src, size_t nelems,       \
                           int pe) {                                           \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_##_name##_get(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);         \
  }                                                                            \
  void shmem_##_name##_get##_nbi(_type *dest, const _type *src, size_t nelems, \
                                 int pe) {                                     \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_##_name##_get##_nbi(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);   \
  }                                                                            \
  void shmem_##_name##_iget(_type *dest, const _type *src, ptrdiff_t tst,      \
                            ptrdiff_t sst, size_t nelems, int pe) {            \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 6);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_##_name##_iget(SHMEM_CTX_DEFAULT, dest, src, tst, sst, nelems,   \
                             pe);                                              \
  }                                                                            \
  void shmem_##_name##_ibget(_type *dest, const _type *src, ptrdiff_t tst,     \
                            ptrdiff_t sst, size_t bsize, size_t nblocks,       \
                            int pe) {                                          \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_COUNTING_MUTIPLE(tst, bsize);                                 \
    SHMEMU_CHECK_COUNTING_MUTIPLE(sst, bsize);                                 \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 6);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_##_name##_ibget(SHMEM_CTX_DEFAULT, dest, src, tst, sst, bsize,   \
                             nblocks, pe);                                     \
  }

/**
 * @brief Macro to declare sized put operations
 * @param _size Size in bits
 */
#define API_DECL_SIZED_PUT(_size)                                              \
  void shmem_put##_size(void *dest, const void *src, size_t nelems, int pe) {  \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_put##_size(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);            \
  }                                                                            \
  void shmem_put##_size##_nbi(void *dest, const void *src, size_t nelems,      \
                              int pe) {                                        \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_put##_size##_nbi(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);      \
  }                                                                            \
  void shmem_iput##_size(void *dest, const void *src, ptrdiff_t tst,           \
                         ptrdiff_t sst, size_t nelems, int pe) {               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 6);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_iput##_size(SHMEM_CTX_DEFAULT, dest, src, tst, sst, nelems, pe); \
  }

/**
 * @brief Macro to declare sized get operations
 * @param _size Size in bits
 */
#define API_DECL_SIZED_GET(_size)                                              \
  void shmem_get##_size(void *dest, const void *src, size_t nelems, int pe) {  \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_get##_size(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);            \
  }                                                                            \
  void shmem_get##_size##_nbi(void *dest, const void *src, size_t nelems,      \
                              int pe) {                                        \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_get##_size##_nbi(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);      \
  }                                                                            \
  void shmem_iget##_size(void *dest, const void *src, ptrdiff_t tst,           \
                         ptrdiff_t sst, size_t nelems, int pe) {               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 6);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_iget##_size(SHMEM_CTX_DEFAULT, dest, src, tst, sst, nelems, pe); \
  }                                                                            \
  void shmem_ibget##_size(void *dest, const void *src, ptrdiff_t tst,          \
                         ptrdiff_t sst, size_t bsize, size_t nblocks, int pe) {\
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_COUNTING_MUTIPLE(tst, bsize);                                 \
    SHMEMU_CHECK_COUNTING_MUTIPLE(sst, bsize);                                 \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 6);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_ibget##_size(SHMEM_CTX_DEFAULT, dest, src, tst, sst, bsize,      \
        nblocks, pe);                                                          \
  }

/**
 * @brief Macro to declare memory put operations
 */
#define API_DECL_PUTMEM()                                                      \
  void shmem_putmem(void *dest, const void *src, size_t nelems, int pe) {      \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_putmem(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);                \
  }                                                                            \
  void shmem_putmem_nbi(void *dest, const void *src, size_t nelems, int pe) {  \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_putmem_nbi(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);            \
  }

/**
 * @brief Macro to declare memory get operations
 */
#define API_DECL_GETMEM()                                                      \
  void shmem_getmem(void *dest, const void *src, size_t nelems, int pe) {      \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_getmem(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);                \
  }                                                                            \
  void shmem_getmem_nbi(void *dest, const void *src, size_t nelems, int pe) {  \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 4);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_getmem_nbi(SHMEM_CTX_DEFAULT, dest, src, nelems, pe);            \
  }

/**
 * @brief Macro to declare typed put operations
 * @param _name Type name
 * @param _type C type
 */
#define API_DECL_TYPED_P(_name, _type)                                         \
  void shmem_##_name##_p(_type *dest, _type src, int pe) {                     \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 3);                                          \
    SHMEMU_CHECK_SYMMETRIC(dest, 1);                                           \
                                                                               \
    shmem_ctx_##_name##_put(SHMEM_CTX_DEFAULT, dest, &src, 1, pe);             \
  }

/**
 * @brief Macro to declare typed get operations
 * @param _name Type name
 * @param _type C type
 */
#define API_DECL_TYPED_G(_name, _type)                                         \
  _type shmem_##_name##_g(const _type *src, int pe) {                          \
    _type val;                                                                 \
                                                                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_PE_ARG_RANGE(pe, 2);                                          \
    SHMEMU_CHECK_SYMMETRIC(src, 2);                                            \
                                                                               \
    shmem_ctx_##_name##_get(SHMEM_CTX_DEFAULT, &val, src, 1, pe);              \
    return val;                                                                \
  }

#endif /* ! _SHMEM_PUTGET_H */
