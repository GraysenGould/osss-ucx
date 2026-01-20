#ifndef _UCC_COLL_H
#define _UCC_COLL_H
#include <ucc/api/ucc.h>
#include "collectives/alltoall.h"
#include <stdio.h>
#include <shmem.h>

typedef struct {
    int rank;
    int size;
} shmem_oob_info_t;

typedef struct {
  void * global_work_buffer;
  ucc_lib_h lib;
  int is_lib_initialized;
  ucc_team_h team_handle; 
  ucc_context_h context_handle;
  ucc_mem_map_t *map_segments;
  shmem_oob_info_t oob_info;
} ucc_coll_component_t;

extern ucc_coll_component_t shmem_ucc_coll;

ucc_status_t ucc_oob_allgather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req);

ucc_status_t ucc_oob_allgather_test(void *req);

ucc_status_t ucc_oob_allgather_free(void *req);

void shmem_ucc_coll_setup ();

void shmem_ucc_coll_finalize();

#endif /* ! _UCC_COLL_H */
