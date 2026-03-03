#include "ucc.h"
#include <ucc/api/ucc.h>
#include "collectives/barrier.h"
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>

#include <stdio.h>

/* ------ UCC Collective Operation ------- */
int ucc_sync_helper(ucc_team_h team_handle) {
  ucc_status_t status;
  ucc_coll_args_t coll_args = {
    .mask = UCC_COLL_ARGS_FIELD_FLAGS | UCC_COLL_ARGS_FIELD_GLOBAL_WORK_BUFFER,
    .coll_type = UCC_COLL_TYPE_BARRIER,
    .global_work_buffer = shmem_ucc_coll.global_work_buffer,
    .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, // delete flag if want to skip memory pre-registration
  };
  
  ucc_coll_req_h coll_handle;
  if (UCC_OK != (status = 
        ucc_collective_init(&coll_args, &coll_handle, team_handle))){
    printf("Could Not Initalize UCC collective. Status: %d\n", status);
    return 1;
  }
  
  if (UCC_OK != ucc_collective_post(coll_handle)){
    printf("Could No Post UCC collective.\n");
    return 1;
  }
  
  /* poll operation until done */
  while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
    /* Drive Collective Progress */
    ucc_context_progress(shmem_ucc_coll.context_handle);
  }
  
  ucc_collective_finalize(coll_handle);  
}

  void ucc_barrier_all() {
    /* Sanity Checks */
    SHMEMU_CHECK_INIT();
    shmem_quiet();                                                             
    shmemc_team_h team_world = (shmemc_team_h) SHMEM_TEAM_WORLD;
    ucc_sync_helper(team_world->ucc_team_handle); // need team world
  }

  void ucc_team_sync(shmem_team_t team) {
    /* Sanity Checks */
    SHMEMU_CHECK_INIT();
    SHMEMU_CHECK_TEAM_VALID(team);
    shmemc_team_h team_h = (shmemc_team_h)team;
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);
    shmem_quiet();                                                             
    ucc_sync_helper(team_h->ucc_team_handle); // need team world
  }

  void ucc_sync_all() {
    /* Sanity Checks */
    SHMEMU_CHECK_INIT();
    shmemc_team_h team_world = (shmemc_team_h) SHMEM_TEAM_WORLD;
    ucc_sync_helper(team_world->ucc_team_handle); // need team world
  }



