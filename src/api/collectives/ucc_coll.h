#ifndef _UCC_COLL_H
#define _UCC_COLL_H

/* #include <ucc_coll.h> */

#define SHMEM_SUCCESS 0

#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

typedef struct {
    int rank;
    int size;
    long *sync_counter; // Symmetric memory pointer
} ucc_shmem_oob_info_t;

typedef struct {
    ucc_shmem_oob_info_t *info;
    int                   is_done;
} oob_request_t;

/* Initialize UCC library handle and global handles */
extern ucc_lib_h ucc_lib;
extern ucc_context_h ucc_global_context;
extern ucc_team_h ucc_team_world;
extern ucc_shmem_oob_info_t global_oob_info;
extern ucc_mem_map_params_t ucc_global_mem_params;
extern ucc_team_h global_ucc_team;


ucc_status_t ucc_oob_all_gather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req);

ucc_status_t ucc_oob_allgather_test(void *req);

ucc_status_t ucc_oob_allgather_free(void *req);


void ucc_coll_init();

void ucc_coll_finalize();

void ucc_coll_team_init (ucc_shmem_oob_info_t * oob_info, ucc_context_h * context_handle, ucc_team_h *team_handle);

void ucc_coll_team_finalize (ucc_team_h team_handle);

void ucc_coll_context_create();

void ucc_coll_context_finalize();

/**
 * @brief Macro to declare type-specific alltoall implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
#define UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)             \
  int ucc_##_typename##_alltoall(                                   \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems);

/**
 * @brief Macro to declare alltoall implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define DECLARE_ALLTOALL_TYPES(_type, _typename)                               \
  UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)  \
SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_ALLTOALL_TYPES)
#undef DECLARE_ALLTOALL_TYPES

/**
 * @brief Macro to declare generic alltoallmem implementations
 *
 * @param _algo Algorithm name to generate declarations for
 */

int ucc_alltoallmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

/**
 * @brief Macro to declare sized alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _size Size in bits
 */
#define UCC_SIZED_ALLTOALL_DECLARATION(_size)                               \
  void ucc_alltoall##_size(                                                    \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync);

/* Declare sized variants for each algorithm */
UCC_SIZED_ALLTOALL_DECLARATION(32)
UCC_SIZED_ALLTOALL_DECLARATION(64)


#endif /* ! _UCC_COLL_H */
