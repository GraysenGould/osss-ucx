#include "ucc.h"
#include <ucc/api/ucc.h>
#include "collectives/broadcast.h"
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>

#include <stdio.h>
#include <math.h>

/* ------ UCC Collective Operation ------- */
inline static void ucc_broadcast_helper(
      void *dest, const void *source, size_t nelems, int PE_start,
      int logPE_stride, int PE_size, long *pSync, ucc_team_h team_handle, 
      int PE_root) {
  printf("Nelems: %d\n", (int) nelems);
  ucc_status_t status;
  ucc_team_attr_t team_attr;
  team_attr.mask = UCC_TEAM_ATTR_FIELD_EP;
  status = ucc_team_get_attr (team_handle, &team_attr);
  if (status == UCC_OK){
    printf("PE %d: endpoint: %d\n", (int) shmem_my_pe(), (int)team_attr.ep);
  }
  if (shmem_my_pe() == PE_root)
    memcpy(dest, source, nelems);

  ucc_coll_buffer_info_t coll_src_buffer_info =  {
    .count = nelems,
    .buffer = (void *) dest,
    .datatype = UCC_DT_UINT8,
    .mem_type = UCC_MEMORY_TYPE_HOST 
  };

  ucc_coll_args_t coll_args = {
    .mask = UCC_COLL_ARGS_FIELD_FLAGS | UCC_COLL_ARGS_FIELD_GLOBAL_WORK_BUFFER,
    .coll_type = UCC_COLL_TYPE_BCAST,
    .src.info = coll_src_buffer_info,
    .global_work_buffer = shmem_ucc_coll.global_work_buffer,
    .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, // delete flag if want to skip memory pre-registration
    .root = PE_root,
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

/**
 * @brief Helper macro to define typed alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 *
 * FIXME: THESE COLLECTIVES NEED TO RETURN NON ZERO IF THEY FAIL
 */
#define UCC_BROADCAST_TYPE_DEFINITION(_type, _typename)                        \
  int ucc_##_typename##_broadcast(                                             \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems,      \
      int PE_root) {                                                           \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_NULL(dest, "dest");                                           \
    SHMEMU_CHECK_NULL(source, "source");                                       \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    SHMEMU_CHECK_SYMMETRIC(dest, sizeof(_type) * nelems * team_h->nranks);     \
    SHMEMU_CHECK_SYMMETRIC(source, sizeof(_type) * nelems * team_h->nranks);   \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source, sizeof(_type) * nelems,          \
                                sizeof(_type) * nelems);                       \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                      "team_h->pSyncs[COLLECTIVE]");                           \
                                                                               \
    ucc_broadcast_helper(                                                      \
        dest, source, nelems * sizeof(_type), team_h->start,                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),                \
        team_h->ucc_team_handle, PE_root);                                     \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

#define DEFINE_BROADCAST_TYPES(_type, _typename)                                \
  UCC_BROADCAST_TYPE_DEFINITION(_type, _typename)

SHMEM_STANDARD_RMA_TYPE_TABLE(DEFINE_BROADCAST_TYPES)
#undef DEFINE_BROADCAST_TYPES

/**
 * @brief Helper macro to define alltoallmem implementations
 *
 */
#define UCC_BROADCASTMEM_DEFINITION()                                           \
  int ucc_broadcastmem(shmem_team_t team, void *dest,                          \
                      const void *source, size_t nelems, int PE_root) {        \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_NULL(dest, "dest");                                           \
    SHMEMU_CHECK_NULL(source, "source");                                       \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    SHMEMU_CHECK_SYMMETRIC(dest, nelems * team_h->nranks);                     \
    SHMEMU_CHECK_SYMMETRIC(source, nelems * team_h->nranks);                   \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source, nelems, nelems);                 \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                      "team_h->pSyncs[COLLECTIVE]");                           \
    ucc_broadcast_helper(                                                       \
        dest, source, nelems, team_h->start,                                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),                \
        team_h->ucc_team_handle, PE_root);                                     \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

#define SHCOLL_BROADCASTMEM_DEFINITION(_algo)                                  \
  int shcoll_broadcastmem_##_algo(shmem_team_t team, void *dest,               \
                                  const void *source, size_t nelems,           \
                                  int PE_root) {                               \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_NULL(dest, "dest");                                           \
    SHMEMU_CHECK_NULL(source, "source");                                       \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    SHMEMU_CHECK_SYMMETRIC(dest, nelems * team_h->nranks);                     \
    SHMEMU_CHECK_SYMMETRIC(source, nelems * team_h->nranks);                   \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source, nelems, nelems);                 \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                      "team_h->pSyncs[COLLECTIVE]");                           \
                                                                               \
    /* Initialize dest: root PE copies source, others leave dest as-is */      \
    /* The broadcast operation will overwrite dest on all PEs */               \
    if (team_h->rank == PE_root)                                               \
      memcpy(dest, source, nelems);                                            \
                                                                               \
    broadcast_helper_##_algo(                                                  \
        dest, source, nelems, PE_root, team_h->start,                          \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

UCC_BROADCASTMEM_DEFINITION()
