#ifndef _UCC_FCOLLECT_H
#define _UCC_FCOLLECT_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

/**
 * @brief Macro to declare type-specific fcollect implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_FCOLLECT_DECLARATION(_type, _typename)             \
  int ucc_##_typename##_fcollect(                                       \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare fcollect implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_FCOLLECT_TYPES(_type, _typename)                        \
  UCC_TYPED_FCOLLECT_DECLARATION(_type, _typename)                   
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_FCOLLECT_TYPES)
#undef DECLARE_FCOLLECT_TYPES

int ucc_fcollectmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

#endif /* ! _UCC_FCOLLECT_H */
