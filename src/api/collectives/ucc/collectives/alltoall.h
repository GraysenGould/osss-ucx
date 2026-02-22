#ifndef _UCC_ALLTOALL_H
#define _UCC_ALLTOALL_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

/**
 * @brief Macro to declare type-specific alltoall implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)             \
  int ucc_##_typename##_alltoall(                                       \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare alltoall implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_ALLTOALL_TYPES(_type, _typename)                        \
  UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)                   
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_ALLTOALL_TYPES)
#undef DECLARE_ALLTOALL_TYPES

int ucc_alltoallmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

#endif /* ! _UCC_ALLTOALL_H */
