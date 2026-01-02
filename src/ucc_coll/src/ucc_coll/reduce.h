/**
 * @file reduction.h
 * @brief Header file containing declarations for collective reduction
 * operations
 *
 * This file provides declarations for various reduction operations (AND, OR,
 * XOR, MIN, MAX, SUM, PROD) across different data types. Multiple algorithm
 * implementations are supported including linear, binomial, recursive doubling,
 * and Rabenseifner's algorithm.
 */

#ifndef _SHCOLL_REDUCTION_H
#define _SHCOLL_REDUCTION_H 1

#include <shmem/teams.h>
#include "shmemu.h"
#include <shmem/api_types.h>

#include <stddef.h>
#include <stdint.h>

/**
 * @brief Macro to declare a single reduction operation
 *
 * @param _typename_op Name of the reduction operation (e.g. int_sum)
 * @param _type Data type to operate on
 * @param _algo Algorithm implementation to use
 */
#define SHCOLL_TO_ALL_DECLARE(_typename_op, _type, _algo)                      \
  void shcoll_##_typename_op##_to_all_##_algo(                                 \
      _type *dest, const _type *source, int nreduce, int PE_start,             \
      int logPE_stride, int PE_size, _type *pWrk, long *pSync)

#define DECLARE_TO_ALL_BITWISE(_type, _typename)                               \
  SHCOLL_TO_ALL_DECLARE(_typename##_and, _type, linear);                       \
  SHCOLL_TO_ALL_DECLARE(_typename##_and, _type, binomial);                     \
  SHCOLL_TO_ALL_DECLARE(_typename##_and, _type, rec_dbl);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_and, _type, rabenseifner);                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_and, _type, rabenseifner2);                \
  SHCOLL_TO_ALL_DECLARE(_typename##_or, _type, linear);                        \
  SHCOLL_TO_ALL_DECLARE(_typename##_or, _type, binomial);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_or, _type, rec_dbl);                       \
  SHCOLL_TO_ALL_DECLARE(_typename##_or, _type, rabenseifner);                  \
  SHCOLL_TO_ALL_DECLARE(_typename##_or, _type, rabenseifner2);                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_xor, _type, linear);                       \
  SHCOLL_TO_ALL_DECLARE(_typename##_xor, _type, binomial);                     \
  SHCOLL_TO_ALL_DECLARE(_typename##_xor, _type, rec_dbl);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_xor, _type, rabenseifner);                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_xor, _type, rabenseifner2);
SHMEM_TO_ALL_BITWISE_TYPE_TABLE(DECLARE_TO_ALL_BITWISE)
#undef DECLARE_TO_ALL_BITWISE

#define DECLARE_TO_ALL_MINMAX(_type, _typename)                                \
  SHCOLL_TO_ALL_DECLARE(_typename##_min, _type, linear);                       \
  SHCOLL_TO_ALL_DECLARE(_typename##_min, _type, binomial);                     \
  SHCOLL_TO_ALL_DECLARE(_typename##_min, _type, rec_dbl);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_min, _type, rabenseifner);                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_min, _type, rabenseifner2);                \
  SHCOLL_TO_ALL_DECLARE(_typename##_max, _type, linear);                       \
  SHCOLL_TO_ALL_DECLARE(_typename##_max, _type, binomial);                     \
  SHCOLL_TO_ALL_DECLARE(_typename##_max, _type, rec_dbl);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_max, _type, rabenseifner);                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_max, _type, rabenseifner2);
SHMEM_TO_ALL_MINMAX_TYPE_TABLE(DECLARE_TO_ALL_MINMAX)
#undef DECLARE_TO_ALL_MINMAX

#define DECLARE_TO_ALL_ARITH(_type, _typename)                                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_sum, _type, linear);                       \
  SHCOLL_TO_ALL_DECLARE(_typename##_sum, _type, binomial);                     \
  SHCOLL_TO_ALL_DECLARE(_typename##_sum, _type, rec_dbl);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_sum, _type, rabenseifner);                 \
  SHCOLL_TO_ALL_DECLARE(_typename##_sum, _type, rabenseifner2);                \
  SHCOLL_TO_ALL_DECLARE(_typename##_prod, _type, linear);                      \
  SHCOLL_TO_ALL_DECLARE(_typename##_prod, _type, binomial);                    \
  SHCOLL_TO_ALL_DECLARE(_typename##_prod, _type, rec_dbl);                     \
  SHCOLL_TO_ALL_DECLARE(_typename##_prod, _type, rabenseifner);                \
  SHCOLL_TO_ALL_DECLARE(_typename##_prod, _type, rabenseifner2);
SHMEM_TO_ALL_ARITH_TYPE_TABLE(DECLARE_TO_ALL_ARITH)
#undef DECLARE_TO_ALL_ARITH

/********************************************************************************/

/**
 * @brief Macro to declare a single reduction operation for a specific type
 *
 * @param _typename Type name used in function name (e.g. int, float)
 * @param _type Actual C type (e.g. int, float)
 * @param _op Operation name (e.g. sum, prod, min, max)
 * @param _algo Algorithm implementation to use
 */
#define SHCOLL_REDUCE_DECLARE(_typename, _type, _op, _algo)                    \
  int shcoll_##_typename##_##_op##_reduce_##_algo(                             \
      shmem_team_t team, _type *dest, const _type *source, size_t nreduce);

#define DECLARE_REDUCE_BITWISE(_type, _typename)                               \
  SHCOLL_REDUCE_DECLARE(_typename, _type, and, linear)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, and, binomial)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, and, rec_dbl)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, and, rabenseifner)                   \
  SHCOLL_REDUCE_DECLARE(_typename, _type, and, rabenseifner2)                  \
  SHCOLL_REDUCE_DECLARE(_typename, _type, or, linear)                          \
  SHCOLL_REDUCE_DECLARE(_typename, _type, or, binomial)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, or, rec_dbl)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, or, rabenseifner)                    \
  SHCOLL_REDUCE_DECLARE(_typename, _type, or, rabenseifner2)                   \
  SHCOLL_REDUCE_DECLARE(_typename, _type, xor, linear)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, xor, binomial)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, xor, rec_dbl)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, xor, rabenseifner)                   \
  SHCOLL_REDUCE_DECLARE(_typename, _type, xor, rabenseifner2)
SHMEM_REDUCE_BITWISE_TYPE_TABLE(DECLARE_REDUCE_BITWISE)
#undef DECLARE_REDUCE_BITWISE

#define DECLARE_REDUCE_MINMAX(_type, _typename)                                \
  SHCOLL_REDUCE_DECLARE(_typename, _type, min, linear)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, min, binomial)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, min, rec_dbl)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, min, rabenseifner)                   \
  SHCOLL_REDUCE_DECLARE(_typename, _type, min, rabenseifner2)                  \
  SHCOLL_REDUCE_DECLARE(_typename, _type, max, linear)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, max, binomial)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, max, rec_dbl)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, max, rabenseifner)                   \
  SHCOLL_REDUCE_DECLARE(_typename, _type, max, rabenseifner2)
SHMEM_REDUCE_MINMAX_TYPE_TABLE(DECLARE_REDUCE_MINMAX)
#undef DECLARE_REDUCE_MINMAX

#define DECLARE_REDUCE_ARITH(_type, _typename)                                 \
  SHCOLL_REDUCE_DECLARE(_typename, _type, sum, linear)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, sum, binomial)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, sum, rec_dbl)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, sum, rabenseifner)                   \
  SHCOLL_REDUCE_DECLARE(_typename, _type, sum, rabenseifner2)                  \
  SHCOLL_REDUCE_DECLARE(_typename, _type, prod, linear)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, prod, binomial)                      \
  SHCOLL_REDUCE_DECLARE(_typename, _type, prod, rec_dbl)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, prod, rabenseifner)                  \
  SHCOLL_REDUCE_DECLARE(_typename, _type, prod, rabenseifner2)
SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_REDUCE_ARITH)
#undef DECLARE_REDUCE_ARITH

#endif /* ! _SHCOLL_REDUCTION_H */
