/**
 * @file alltoalls.h
 * @brief Header file for strided all-to-all collective operations
 *
 * This header declares the interfaces for strided all-to-all collective
 * operations using different algorithms and synchronization methods:
 * - Shift exchange
 * - XOR pairwise exchange
 * - Color pairwise exchange
 *
 */

#ifndef _UCC_ALLTOALLS_H
#define _UCC_ALLTOALLS_H 1

#include <shmem/teams.h>
#include <shmem/api_types.h>

/**
 * @brief Macro to declare type-specific strided alltoall implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_ALLTOALLS_DECLARATION(_type, _typename)            \
  int ucc_##_typename##_alltoalls(                                  \
      shmem_team_t team, _type *dest, const _type *source, ptrdiff_t dst,      \
      ptrdiff_t sst, size_t nelems);

/**
 * @brief Macro to declare strided alltoall implementations for all supported
 * types
 *
 * @param _type Data type
 * @param _typename Type name string
 */
#define DECLARE_ALLTOALLS_TYPES(_type, _typename)                              \
  UCC_TYPED_ALLTOALLS_DECLARATION(_type, _typename) \
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_ALLTOALLS_TYPES)
#undef DECLARE_ALLTOALLS_TYPES

/**
 * @brief Macro to declare type-specific strided alltoall memory implementation
 *
 * @param _algo Algorithm name
 */
int ucc_alltoallsmem(shmem_team_t team, void *dest,
    const void *source, ptrdiff_t dst, ptrdiff_t sst, size_t nelems);

/**
 * @brief Macro to declare sized strided alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _size Size in bits
 */
#define UCC_SIZED_ALLTOALLS_DECLARATION(_size)                       \
  void ucc_alltoalls##_size(                                      \
      void *dest, const void *source, ptrdiff_t dst, ptrdiff_t sst,            \
      size_t nelems, int PE_start, int logPE_stride, int PE_size,              \
      long *pSync);

/* Declare sized variants for each algorithm */
UCC_SIZED_ALLTOALLS_DECLARATION(32)
UCC_SIZED_ALLTOALLS_DECLARATION(64)
#endif /* ! _UCC_ALLTOALLS_H */
