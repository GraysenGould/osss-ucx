/**
 * @file scan.h
 * @brief Header file containing declarations for scan operations
 *
 */

#ifndef _SHCOLL_SCAN_H
#define _SHCOLL_SCAN_H 1

#include <shmem/teams.h>
#include "shmemu.h"
#include <shmem/api_types.h>

#include <stddef.h>
#include <stdint.h>

/**
 * @brief Macro to declare a single inscan operation
 *
 * @param _typename Name of datatype
 * @param _type Data type to operate on
 * @param _algo Algorithm implementation to use
 */
#define SHCOLL_INSCAN_DECLARE(_typename, _type, _algo)                      \
  int shcoll_##_typename##_inscan_##_algo(                                  \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare a single exscan operation
 *
 * @param _typename Name of datatype
 * @param _type Data type to operate on
 * @param _algo Algorithm implementation to use
 */
#define SHCOLL_EXSCAN_DECLARE(_typename, _type, _algo)                      \
  int shcoll_##_typename##_exscan_##_algo(                                  \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

#define DECLARE_INSCAN(_type, _typename)                               \
  SHCOLL_INSCAN_DECLARE(_typename, _type, linear)

SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_INSCAN)
#undef DECLARE_INSCAN

#define DECLARE_EXSCAN(_type, _typename)                               \
  SHCOLL_EXSCAN_DECLARE(_typename, _type, linear)

SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_EXSCAN)
#undef DECLARE_EXSCAN

#endif /* ! _SHCOLL_SCAN_H */
