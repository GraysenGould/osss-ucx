/**
 * @file fcollect.h
 * @brief Header file for fixed-size collect collective operations
 *
 * This header declares the interfaces for fixed-size collect collective
 * operations using different algorithms and synchronization methods:
 * - Linear
 * - All linear
 * - Recursive doubling
 * - Ring
 * - Bruck
 * - Neighbor exchange
 *
 * Some algorithms have variants using different synchronization:
 * - Signal-based
 */

#ifndef _SHCOLL_FCOLLECT_H
#define _SHCOLL_FCOLLECT_H 1

#include <shmem/teams.h>
#include <shmem/api_types.h>

/**
 * @brief Macro to declare type-specific fcollect implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define SHCOLL_TYPED_FCOLLECT_DECLARATION(_algo, _type, _typename)             \
  int shcoll_##_typename##_fcollect_##_algo(                                   \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare fcollect implementations for all supported types
 *
 * @param _type Data type
 * @param _typename Type name string
 */
#define DECLARE_FCOLLECT_TYPES(_type, _typename)                               \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(linear, _type, _typename)                  \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(all_linear, _type, _typename)              \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(all_linear1, _type, _typename)             \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(rec_dbl, _type, _typename)                 \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(rec_dbl_signal, _type, _typename)          \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(ring, _type, _typename)                    \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(bruck, _type, _typename)                   \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(bruck_no_rotate, _type, _typename)         \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(bruck_signal, _type, _typename)            \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(bruck_inplace, _type, _typename)           \
  SHCOLL_TYPED_FCOLLECT_DECLARATION(neighbor_exchange, _type, _typename)

SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_FCOLLECT_TYPES)
#undef DECLARE_FCOLLECT_TYPES

/**
 * @brief Macro to declare type-specific fcollect memory implementation
 *
 * @param _algo Algorithm name
 */
#define SHCOLL_FCOLLECTMEM_DECLARATION(_algo)                                  \
  int shcoll_fcollectmem_##_algo(shmem_team_t team, void *dest,                \
                                 const void *source, size_t nelems);

/* Declare all algorithm variants */
SHCOLL_FCOLLECTMEM_DECLARATION(linear)
SHCOLL_FCOLLECTMEM_DECLARATION(all_linear)
SHCOLL_FCOLLECTMEM_DECLARATION(all_linear1)
SHCOLL_FCOLLECTMEM_DECLARATION(rec_dbl)
SHCOLL_FCOLLECTMEM_DECLARATION(rec_dbl_signal)
SHCOLL_FCOLLECTMEM_DECLARATION(ring)
SHCOLL_FCOLLECTMEM_DECLARATION(bruck)
SHCOLL_FCOLLECTMEM_DECLARATION(bruck_no_rotate)
SHCOLL_FCOLLECTMEM_DECLARATION(bruck_signal)
SHCOLL_FCOLLECTMEM_DECLARATION(bruck_inplace)
SHCOLL_FCOLLECTMEM_DECLARATION(neighbor_exchange)

/*
 * @brief Macro to declare sized fcollect implementations
 *
 * @param _algo Algorithm name
 * @param _size Size in bits
 */
#define SHCOLL_SIZED_FCOLLECT_DECLARATION(_algo, _size)                        \
  void shcoll_fcollect##_size##_##_algo(                                       \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync);

/* Declare sized variants for each algorithm */
SHCOLL_SIZED_FCOLLECT_DECLARATION(linear, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(linear, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(all_linear, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(all_linear, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(all_linear1, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(all_linear1, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(rec_dbl, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(rec_dbl, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(rec_dbl_signal, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(rec_dbl_signal, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(ring, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(ring, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck_no_rotate, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck_no_rotate, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck_signal, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck_signal, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck_inplace, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(bruck_inplace, 64)

SHCOLL_SIZED_FCOLLECT_DECLARATION(neighbor_exchange, 32)
SHCOLL_SIZED_FCOLLECT_DECLARATION(neighbor_exchange, 64)

#endif /* ! _SHCOLL_FCOLLECT_H */
