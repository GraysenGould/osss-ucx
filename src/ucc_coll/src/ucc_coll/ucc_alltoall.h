/**
 * @file alltoall.h
 * @brief Header file for all-to-all collective operations
 *
 * This header declares the interfaces for all-to-all collective operations
 * using different algorithms and synchronization methods:
 * - Shift exchange
 * - XOR pairwise exchange
 * - Color pairwise exchange
 *
 * Each algorithm has variants using different synchronization:
 * - Barrier-based
 * - Signal-based
 * - Counter-based
 */

#ifndef _UCC_ALLTOALL_H
#define _UCC_ALLTOALL_H 1

#include <shmem/teams.h>
#include <shmem/api_types.h>
#include "shmemu.h"

/**
 * @brief Macro to declare type-specific alltoall implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)             \
  int ucc_##_typename##_alltoall(                                   \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare alltoall implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_ALLTOALL_TYPES(_type, _typename)                               \
  UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)  \
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_ALLTOALL_TYPES)
#undef DECLARE_ALLTOALL_TYPES

/**
 * @brief Macro to declare generic alltoallmem implementations
 *
 * @param _algo Algorithm name to generate declarations for
 */

int ucc_alltoallmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

/**
 * @brief Macro to declare sized alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _size Size in bits
 */
#define UCC_SIZED_ALLTOALL_DECLARATION(_size)                               \
  void ucc_alltoall##_size(                                                    \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync);

/* Declare sized variants for each algorithm */
UCC_SIZED_ALLTOALL_DECLARATION(32)
UCC_SIZED_ALLTOALL_DECLARATION(64)

#endif /* ! _UCC_ALLTOALL_H */
