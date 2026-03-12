#include "ucc.h"
#include <ucc/api/ucc.h>
#include "collectives/collect.h"
#include "collectives/fcollect.h"
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>

#include <stdio.h>
#include <math.h>

inline static void ucc_fcollect_helper(
      void *dest, const void *source, size_t nelems, int PE_size, ucc_team_h team_handle) {
  ucc_status_t status;
  ucc_coll_buffer_info_t coll_src_buffer_info =  {
    .buffer = source,
    .count = nelems,
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
    .coll_type = UCC_COLL_TYPE_ALLGATHER,
    .src.info = coll_src_buffer_info,
    .dst.info = coll_dst_buffer_info,
    .global_work_buffer = shmem_ucc_coll.global_work_buffer,
    .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, // delete flag if want to skip memory pre-registration
  };
  
  ucc_coll_req_h coll_handle;
  if (UCC_OK != (status = 
        ucc_collective_init(&coll_args, &coll_handle, team_handle))){
    printf("Could Not Initalize UCC collective. Status: %d\n", status);
    return;
  }
  
  if (UCC_OK != ucc_collective_post(coll_handle)){
    printf("Could No Post UCC collective.\n");
    return;
  }
  
  /* poll operation until done */
  while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
    /* Drive Collective Progress */
    ucc_context_progress(shmem_ucc_coll.context_handle);
  }
  
  ucc_collective_finalize(coll_handle);  
}

/* ------ UCC Collective Operation ------- */
inline static void ucc_collect_helper(
      void *dest, const void *source, size_t nelems, size_t elem_size, int PE_start,
      int logPE_stride, int PE_size, long *pSync, ucc_team_h team_handle) {
  ucc_status_t status;

  ucc_count_t * my_count = (ucc_count_t *) shmem_malloc(sizeof(ucc_count_t)); 
  ucc_count_t * counts = (ucc_count_t *) shmem_malloc(PE_size * sizeof(ucc_count_t));
  *my_count = nelems * elem_size;

  ucc_fcollect_helper(counts, my_count, sizeof(ucc_count_t), PE_size, team_handle);

  ucc_aint_t * displacements = (ucc_aint_t *) malloc(PE_size * sizeof(ucc_aint_t));
  *displacements = 0;
  for (int i = 1; i < PE_size; i ++){
    displacements[i] = displacements[i - 1] + counts[i - 1];
  }

  
  ucc_coll_buffer_info_t coll_src_buffer_info =  {
    .buffer = source,
    .count = nelems * elem_size,
    .datatype = UCC_DT_UINT8,
    .mem_type = UCC_MEMORY_TYPE_HOST 
  };

  ucc_coll_buffer_info_v_t coll_dst_buffer_info =  {
    .buffer = dest,
    .counts = counts,
    .displacements = displacements,
    .datatype = UCC_DT_UINT8,
    .mem_type = UCC_MEMORY_TYPE_HOST 
  };

  ucc_coll_args_t coll_args = {
    .mask = UCC_COLL_ARGS_FIELD_FLAGS | UCC_COLL_ARGS_FIELD_GLOBAL_WORK_BUFFER,
    .coll_type = UCC_COLL_TYPE_ALLGATHERV,
    .src.info = coll_src_buffer_info,
    .dst.info_v = coll_dst_buffer_info,
    .global_work_buffer = shmem_ucc_coll.global_work_buffer,
    .flags = UCC_COLL_ARGS_FLAG_COUNT_64BIT | UCC_COLL_ARGS_FLAG_DISPLACEMENTS_64BIT, 
    // .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, // delete flag if want to skip memory pre-registration
  };
  
  ucc_coll_req_h coll_handle;
  if (UCC_OK != (status = 
        ucc_collective_init(&coll_args, &coll_handle, team_handle))){
    printf("Could Not Initalize UCC collective. Status: %d\n", status);
    return;
  }
  
  if (UCC_OK != ucc_collective_post(coll_handle)){
    printf("Could No Post UCC collective.\n");
    return;
  }
  
  /* poll operation until done */
  while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
    /* Drive Collective Progress */
    ucc_context_progress(shmem_ucc_coll.context_handle);
  }
  
  ucc_collective_finalize(coll_handle);  
  shmem_free(counts);
  shmem_free(my_count);
  free(displacements);
}

/**
 * @brief Helper macro to define typed collect implementations
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 *
 * FIXME: THESE COLLECTIVES NEED TO RETURN NON ZERO IF THEY FAIL
 */
#define UCC_COLLECT_TYPE_DEFINITION(_type, _typename)                      \
  int ucc_##_typename##_collect(                                              \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems) {    \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_NULL(dest, "dest");                                           \
    SHMEMU_CHECK_NULL(source, "source");                                       \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    SHMEMU_CHECK_SYMMETRIC(dest, sizeof(_type) * nelems * team_h->nranks);     \
    SHMEMU_CHECK_SYMMETRIC(source, sizeof(_type) * nelems);                     \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source,                                  \
                                sizeof(_type) * nelems * team_h->nranks,       \
                                sizeof(_type) * nelems);                        \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                      "team_h->pSyncs[COLLECTIVE]");                           \
                                                                               \
    ucc_collect_helper(                                                       \
        dest, source, nelems, sizeof(_type), team_h->start,                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),                \
        team_h->ucc_team_handle);                                              \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

#define DEFINE_COLLECT_TYPES(_type, _typename)                                \
  UCC_COLLECT_TYPE_DEFINITION(_type, _typename)

SHMEM_STANDARD_RMA_TYPE_TABLE(DEFINE_COLLECT_TYPES)
#undef DEFINE_COLLECT_TYPES

/**
 * @brief routine to define collectmem implementations
 *
 */
int ucc_collectmem(shmem_team_t team, void *dest,
                                 const void *source, size_t nelems) {
    SHMEMU_CHECK_INIT();
    SHMEMU_CHECK_TEAM_VALID(team);
    shmemc_team_h team_h = (shmemc_team_h)team;
    SHMEMU_CHECK_NULL(dest, "dest");
    SHMEMU_CHECK_NULL(source, "source");
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);
    SHMEMU_CHECK_SYMMETRIC(dest, nelems * team_h->nranks);
    SHMEMU_CHECK_SYMMETRIC(source, nelems);
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source,
                                nelems * team_h->nranks,
                                nelems);
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),
                      "team_h->pSyncs[COLLECTIVE]");
    ucc_collect_helper(
        dest, source, nelems, sizeof(char), team_h->start,
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,
        team_h->nranks,
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),
        team_h->ucc_team_handle);
    
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);
    
    return 0;
  }
