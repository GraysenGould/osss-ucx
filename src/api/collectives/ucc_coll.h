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


ucc_status_t ucc_oob_all_gather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req);

ucc_status_t ucc_oob_allgather_test(void *req);

ucc_status_t ucc_oob_allgather_free(void *req);

int ucc_alltoallmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

#endif /* ! _UCC_COLL_H */
