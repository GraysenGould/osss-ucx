/**
 * @file inscan.h
 * @brief Header file containing declarations for  inscan operations
 *
 */

#ifndef _SHCOLL_INSCAN_H
#define _SHCOLL_INSCAN_H 1

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

#define DECLARE_INSCAN(_type, _typename)                               \
  SHCOLL_INSCAN_DECLARE(_typename, _type, linear)

SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_INSCAN)
#undef DECLARE_INSCAN

#endif /* ! _SHCOLL_INSCAN_H */
