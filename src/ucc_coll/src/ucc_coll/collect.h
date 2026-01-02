/**
 * @file collect.h
 * @brief Header file for collect collective operations
 *
 * This header declares the interfaces for collect collective operations
 * using different algorithms and synchronization methods:
 * - Linear
 * - All linear
 * - Recursive doubling
 * - Ring
 * - Bruck
 * - Simple
 *
 * Some algorithms have variants using different synchronization:
 * - Signal-based
 */

#ifndef _SHCOLL_COLLECT_H
#define _SHCOLL_COLLECT_H 1

#include <shmem/teams.h>
#include <shmem/api_types.h>

/**
 * @brief Macro to declare type-specific collect implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define SHCOLL_TYPED_COLLECT_DECLARATION(_algo, _type, _typename)              \
  int shcoll_##_typename##_collect_##_algo(                                    \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare collect implementations for all supported types
 *
 * @param _type Data type
 * @param _typename Type name string
 */
#define DECLARE_COLLECT_TYPES(_type, _typename)                                \
  SHCOLL_TYPED_COLLECT_DECLARATION(linear, _type, _typename)                   \
  SHCOLL_TYPED_COLLECT_DECLARATION(all_linear, _type, _typename)               \
  SHCOLL_TYPED_COLLECT_DECLARATION(all_linear1, _type, _typename)              \
  SHCOLL_TYPED_COLLECT_DECLARATION(rec_dbl, _type, _typename)                  \
  SHCOLL_TYPED_COLLECT_DECLARATION(rec_dbl_signal, _type, _typename)           \
  SHCOLL_TYPED_COLLECT_DECLARATION(ring, _type, _typename)                     \
  SHCOLL_TYPED_COLLECT_DECLARATION(bruck, _type, _typename)                    \
  SHCOLL_TYPED_COLLECT_DECLARATION(bruck_no_rotate, _type, _typename)          \
  SHCOLL_TYPED_COLLECT_DECLARATION(simple, _type, _typename)

SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_COLLECT_TYPES)
#undef DECLARE_COLLECT_TYPES

/**
 * @brief Macro to declare type-specific collect memory implementation
 *
 * @param _algo Algorithm name
 */
#define SHCOLL_COLLECTMEM_DECLARATION(_algo)                                   \
  int shcoll_collectmem_##_algo(shmem_team_t team, void *dest,                 \
                                const void *source, size_t nelems);

/* Declare all algorithm variants */
SHCOLL_COLLECTMEM_DECLARATION(linear)
SHCOLL_COLLECTMEM_DECLARATION(all_linear)
SHCOLL_COLLECTMEM_DECLARATION(all_linear1)
SHCOLL_COLLECTMEM_DECLARATION(rec_dbl)
SHCOLL_COLLECTMEM_DECLARATION(rec_dbl_signal)
SHCOLL_COLLECTMEM_DECLARATION(ring)
SHCOLL_COLLECTMEM_DECLARATION(bruck)
SHCOLL_COLLECTMEM_DECLARATION(bruck_no_rotate)
SHCOLL_COLLECTMEM_DECLARATION(simple)

/**
 * @brief Macro to declare sized collect implementations
 *
 * @param _algo Algorithm name
 * @param _size Size in bits
 */
#define SHCOLL_SIZED_COLLECT_DECLARATION(_algo, _size)                         \
  void shcoll_collect##_size##_##_algo(                                        \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync);

/* Declare sized variants for each algorithm */
SHCOLL_SIZED_COLLECT_DECLARATION(linear, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(linear, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(all_linear, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(all_linear, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(all_linear1, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(all_linear1, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(rec_dbl, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(rec_dbl, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(rec_dbl_signal, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(rec_dbl_signal, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(ring, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(ring, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(bruck, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(bruck, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(bruck_no_rotate, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(bruck_no_rotate, 64)

SHCOLL_SIZED_COLLECT_DECLARATION(simple, 32)
SHCOLL_SIZED_COLLECT_DECLARATION(simple, 64)

#endif /* ! _SHCOLL_COLLECT_H */
