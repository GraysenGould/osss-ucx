#ifndef _UCC_ALLTOALLS_H
#define _UCC_ALLTOALLS_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

/**
 * @brief Macro to declare type-specific alltoalls implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_ALLTOALLS_DECLARATION(_type, _typename)                  \
  int ucc_##_typename##_alltoalls(                                         \
      shmem_team_t team, _type *dest, const _type *source, ptrdiff_t dst,  \
      ptrdiff_t sst, size_t nelems);

/**
 * @brief Macro to declare alltoalls implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_ALLTOALLS_TYPES(_type, _typename)                        \
  UCC_TYPED_ALLTOALLS_DECLARATION(_type, _typename)                   
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_ALLTOALLS_TYPES)
#undef DECLARE_ALLTOALLS_TYPES

int ucc_alltoallsmem(shmem_team_t team, void *dest, const void *source, ptrdiff_t dst, 
    ptrdiff_t sst, size_t nelems);

#endif /* ! _UCC_ALLTOALLS_H */
