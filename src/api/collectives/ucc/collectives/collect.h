#ifndef _UCC_COLLECT_H
#define _UCC_COLLECT_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

/**
 * @brief Macro to declare type-specific collect implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_COLLECT_DECLARATION(_type, _typename)             \
  int ucc_##_typename##_collect(                                       \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare collect implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_COLLECT_TYPES(_type, _typename)                        \
  UCC_TYPED_COLLECT_DECLARATION(_type, _typename)                   
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_COLLECT_TYPES)
#undef DECLARE_COLLECT_TYPES

int ucc_collectmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

#endif /* ! _UCC_COLLECT_H */
