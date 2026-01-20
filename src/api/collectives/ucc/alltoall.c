#include "ucc.h"
#include <ucc/api/ucc.h>
#include "collectives/alltoall.h"
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>

#include <stdio.h>
#include <math.h>

/* ------ UCC Collective Operation ------- */
inline static void ucc_alltoallmem_helper(
      void *dest, const void *source, size_t nelems, int PE_start,
      int logPE_stride, int PE_size, long *pSync) {
  ucc_status_t status;
  ucc_coll_buffer_info_t coll_src_buffer_info =  {
    .buffer = source,
    .count = nelems * PE_size,
    .datatype = UCC_DT_UINT8,
    .mem_type = UCC_MEMORY_TYPE_HOST 
  };

  ucc_coll_buffer_info_t coll_dst_buffer_info =  {
    .buffer = dest,
    .count = nelems * PE_size,
    .datatype = UCC_DT_UINT8,
    .mem_type = UCC_MEMORY_TYPE_HOST 
  };

  ucc_coll_args_t coll_args = {
    .mask = UCC_COLL_ARGS_FIELD_FLAGS | UCC_COLL_ARGS_FIELD_GLOBAL_WORK_BUFFER,
    .coll_type = UCC_COLL_TYPE_ALLTOALL,
    .src.info = coll_src_buffer_info,
    .dst.info = coll_dst_buffer_info,
    .global_work_buffer = shmem_ucc_coll.global_work_buffer,
    .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, // delete flag if want to skip memory pre-registration
  };
  //printf("DEBUG: Part 9\n");
  ucc_coll_req_h coll_handle;
  if (UCC_OK != (status = 
        ucc_collective_init(&coll_args, &coll_handle, shmem_ucc_coll.team_handle))){
    printf("Could Not Initalize UCC collective. Status: %d\n", status);
    return;
  }
  if (UCC_OK != ucc_collective_post(coll_handle)){
    printf("Could No Post UCC collective.\n");
    return;
  }
  
  //printf("DEBUG: Part 10\n");
  /* poll operation until done */
  while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
    /* Drive Collective Progress */
    ucc_context_progress(shmem_ucc_coll.context_handle);
  }

  //printf("DEBUG: Part 11\n");
  ucc_collective_finalize(coll_handle); 
  //printf("DEBUG: Part 12\n");
}

/**
 * @brief Helper macro to define alltoallmem implementations
 *
 */
#define UCC_ALLTOALLMEM_DEFINITION(_placeholder)                               \
  int ucc_alltoallmem(shmem_team_t team, void *dest,                           \
                                 const void *source, size_t nelems) {          \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_NULL(dest, "dest");                                           \
    SHMEMU_CHECK_NULL(source, "source");                                       \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    SHMEMU_CHECK_SYMMETRIC(dest, nelems * team_h->nranks);                     \
    SHMEMU_CHECK_SYMMETRIC(source, nelems * team_h->nranks);                   \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source, nelems * team_h->nranks,         \
                                nelems * team_h->nranks);                      \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                      "team_h->pSyncs[COLLECTIVE]");                           \
    ucc_alltoallmem_helper(                                                       \
        dest, source, nelems, team_h->start,                                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));                \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

UCC_ALLTOALLMEM_DEFINITION(placeholder)
