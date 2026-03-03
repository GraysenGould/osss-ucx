#ifndef _UCC_COLL_H
#define _UCC_COLL_H
#include <ucc/api/ucc.h>
#include "collectives/alltoall.h"
#include "collectives/alltoalls.h"
#include "collectives/broadcast.h"
#include "collectives/reduce.h"
#include <stdio.h>
#include <shmem.h>

#define UCC_DEBUG 1

#ifdef UCC_DEBUG
#define UCC_LOG_DEBUG(_action, _function)            \
  do {                                                \
   printf("DEBUG: %s in %s\n",_action, _function);     \
  }while(0);
#else
#define UCC_LOG_DEBUG(_action, _function)            \
  do {                                                \
  }while(0);
#endif

typedef struct {
    int rank;
    int size;
} shmem_oob_info_t;

typedef struct {
  void * global_work_buffer;
  ucc_lib_h lib;
  int libucc_initialized;
  int ucc_team_initialized;
  ucc_context_h context_handle;
  ucc_mem_map_t *map_segments;
  shmem_oob_info_t oob_info;
  ucc_lib_params_t lib_params;
  ucc_context_params_t context_params;
  ucc_context_attr_t ctx_attr;
} ucc_coll_component_t;


extern ucc_coll_component_t shmem_ucc_coll;

ucc_status_t ucc_oob_allgather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req);

ucc_status_t ucc_oob_allgather_test(void *req);

ucc_status_t ucc_oob_allgather_free(void *req);

void shmem_ucc_coll_setup ();

void shmem_ucc_team_setup(ucc_team_h * team_handle);

void shmem_ucc_coll_finalize();

void shmem_ucc_team_finalize(ucc_team_h team_handle);

#endif /* ! _UCC_COLL_H */
