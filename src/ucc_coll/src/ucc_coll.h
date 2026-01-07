#ifndef _UCC_COLL_H
#define _UCC_COLL_H

/* #include <ucc_coll.h> */

#define SHMEM_SUCCESS 0

#include "ucc_coll/ucc_alltoall.h"
#include <ucc/api/ucc.h>
#include <stdio.h>

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

#endif /* ! _UCC_COLL_H */
