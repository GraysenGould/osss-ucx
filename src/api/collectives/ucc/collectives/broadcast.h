#ifndef _UCC_BROADCAST_H
#define _UCC_BROADCAST_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

/**
 * @brief Macro to declare type-specific broadcast implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_BROADCAST_DECLARATION(_type, _typename)             \
  int ucc_##_typename##_broadcast(                                   \
      shmem_team_t team, _type *dest, const _type *source,           \
      size_t nelems, int PE_root);

/**
 * @brief Macro to declare broadcast implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_BROADCAST_TYPES(_type, _typename)                        \
  UCC_TYPED_BROADCAST_DECLARATION(_type, _typename)                   
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_BROADCAST_TYPES)
#undef DECLARE_BROADCAST_TYPES

int ucc_broadcastmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems, int PE_root);

#endif /* ! _UCC_BROADCAST_H */
