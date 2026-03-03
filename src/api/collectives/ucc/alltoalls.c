#include "ucc.h"
#include <ucc/api/ucc.h>
#include "collectives/alltoalls.h"
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>

#include <stdio.h>
#include <math.h>

/* ------ UCC Collective Operation ------- */
inline static void ucc_alltoalls_helper(
      void *dest, const void *source, size_t nelems, size_t elem_size, int PE_start,
      int logPE_stride, int PE_size, long *pSync, ptrdiff_t dst, 
      ptrdiff_t sst, ucc_team_h team_handle) {
  ucc_status_t status;
  void * source_packed = (void *) shmem_malloc(nelems * elem_size * PE_size);
  void * dest_packed = (void *) shmem_malloc(nelems * elem_size * PE_size);

  if (source_packed == NULL || dest_packed == NULL){
    printf("Could not allocate memory for alltoalls\n");
    if (source_packed == NULL)
      shmem_free(source_packed);
    if (dest_packed == NULL)
      shmem_free(dest_packed);
    return;
  }

  /* pack source elements for alltoall */
  char * src_offset, * packed_offset;
  for (int i = 0; i < nelems * PE_size; i ++){
    src_offset = (char *) source + i * sst * elem_size;
    packed_offset = (char *) source_packed + i * elem_size;
    memcpy(packed_offset, src_offset, elem_size);
  }
  
  ucc_coll_buffer_info_t coll_src_buffer_info =  {
    .buffer = source_packed,
    .count = nelems * elem_size * PE_size,
    .datatype = UCC_DT_UINT8,
    .mem_type = UCC_MEMORY_TYPE_HOST 
  };

  ucc_coll_buffer_info_t coll_dst_buffer_info =  {
    .buffer = dest_packed,
    .count = nelems * elem_size * PE_size,
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

  /* Unpack elements into dest array */
  for (int i = 0; i < nelems * PE_size; i ++){
    char * dst_offset = (char *) dest + i * dst * elem_size;
    char * packed_offset = (char *) dest_packed + i * elem_size;
    memcpy(dst_offset, packed_offset, elem_size);
  }

  shmem_free(source_packed);
  shmem_free(dest_packed);
}

/**
 * @brief Helper macro to define typed alltoalls implementations
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 *
 * FIXME: THESE COLLECTIVES NEED TO RETURN NON ZERO IF THEY FAIL
 */
#define UCC_ALLTOALLS_TYPE_DEFINITION(_type, _typename)                        \
  int ucc_##_typename##_alltoalls(                                             \
      shmem_team_t team, _type *dest, const _type *source, ptrdiff_t dst,      \
      ptrdiff_t sst, size_t nelems) {                                          \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    SHMEMU_CHECK_SYMMETRIC(dest, sizeof(_type) * nelems * team_h->nranks);     \
    SHMEMU_CHECK_SYMMETRIC(source, sizeof(_type) * nelems * team_h->nranks);   \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source,                                  \
                                sizeof(_type) * nelems * team_h->nranks,       \
                                sizeof(_type) * nelems * team_h->nranks);      \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                      "team_h->pSyncs[COLLECTIVE]");                           \
                                                                               \
    ucc_alltoalls_helper(                                                       \
        dest, source, nelems, sizeof(_type), team_h->start,                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE), dst, sst,      \
        team_h->ucc_team_handle);                                              \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

#define DEFINE_ALLTOALLS_TYPES(_type, _typename)                                \
  UCC_ALLTOALLS_TYPE_DEFINITION(_type, _typename)

SHMEM_STANDARD_RMA_TYPE_TABLE(DEFINE_ALLTOALLS_TYPES)
#undef DEFINE_ALLTOALLS_TYPES

/**
 * @brief Helper macro to define alltoallsmem implementations
 *
 */
int ucc_alltoallsmem(shmem_team_t team, void *dest,                          
          const void *source, ptrdiff_t dst, ptrdiff_t sst, size_t nelems) {   
    SHMEMU_CHECK_INIT();
    SHMEMU_CHECK_TEAM_VALID(team);
    shmemc_team_h team_h = (shmemc_team_h)team;
    SHMEMU_CHECK_NULL(dest, "dest");
    SHMEMU_CHECK_NULL(source, "source");
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);
    SHMEMU_CHECK_SYMMETRIC(dest, nelems * team_h->nranks);
    SHMEMU_CHECK_SYMMETRIC(source, nelems * team_h->nranks);
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source, nelems * team_h->nranks,
                                nelems * team_h->nranks);
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),
                      "team_h->pSyncs[COLLECTIVE]");
    ucc_alltoalls_helper(
        dest, source, nelems, 1, team_h->start,
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,
        team_h->nranks,
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE), dst, sst,
        team_h->ucc_team_handle);
    
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);
    return 0;
  }
