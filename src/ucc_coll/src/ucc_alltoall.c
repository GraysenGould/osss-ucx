#include <shmem/api_types.h>
#include "ucc_coll.h"
#include "ucc_coll/ucc_alltoall.h"

#include <string.h>
#include <limits.h>
#include <assert.h>
#include <shmem.h>

#include <stdio.h>
#include <math.h>

#include <ucc/api/ucc.h>
#include <ucc/api/ucc_version.h>
#include <ucc/api/ucc_status.h>
#include <ucc/api/ucc_def.h>


//generic: needs to take a variety of differnt collective types
// collectives where only size, not type, matters

#define UCC_COLLECTIVE_HELPER(_coll_type, _ucc_type_name)                     \
inline static void ucc_##_coll_type##_helper(                                 \
      void *dest, const void *source, size_t nelems, int PE_start,            \
      int logPE_stride, int PE_size, long *pSync) {                           \
    ucc_coll_buffer_info_t coll_src_buffer_info =  {                          \
      .buffer = source,                                                       \
      .count = nelems * PE_size,                                              \
      .datatype = UCC_DT_UINT8,                                               \
      .mem_type = UCC_MEMORY_TYPE_HOST                                        \
    };                                                                        \
    ucc_coll_buffer_info_t coll_dst_buffer_info =  {                          \
      .buffer = dest,                                                         \
      .count = nelems * PE_size,                                              \
      .datatype = UCC_DT_UINT8,                                               \
      .mem_type = UCC_MEMORY_TYPE_HOST                                        \
    };                                                                        \
    ucc_coll_args_t coll_args = {                                             \
      .mask = 0, /* just the bare minimum */                                  \
      .coll_type = _ucc_type_name,                                            \
      .src.info = coll_src_buffer_info,                                       \
      .dst.info = coll_dst_buffer_info,                                       \
    };                                                                        \
    ucc_team_h team_handle = NULL; /* needs to be replaced */                 \
    ucc_context_h context_handle = NULL; /* needs to be replaced */           \
    ucc_coll_req_h coll_handle;                                               \
    ucc_collective_init(&coll_args, &coll_handle, team_handle);               \
    ucc_collective_post(coll_handle);                                         \
    /* poll operation until done */                                           \
    while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {               \
      ucc_context_progress(context_handle);                                   \
      /* This actually drives the communication */                            \
    }                                                                         \
    ucc_collective_finalize(coll_handle);                                     \
}


UCC_COLLECTIVE_HELPER(alltoall, UCC_COLL_TYPE_ALLTOALL)

/**
 * @brief Helper macro to define typed alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 *
 * FIXME: THESE COLLECTIVES NEED TO RETURN NON ZERO IF THEY FAIL
 */
#define UCC_ALLTOALL_TYPE_DEFINITION(_type, _typename)               \
  int ucc_##_typename##_alltoall_ucc(                                   \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems) {    \
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
    ucc_alltoall_helper(                                                       \
        dest, source, nelems * sizeof(_type), team_h->start,                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

#define DEFINE_ALLTOALL_TYPES(_type, _typename)             \
  UCC_ALLTOALL_TYPE_DEFINITION(_type, _typename)         \
  SHMEM_STANDARD_RMA_TYPE_TABLE(DEFINE_ALLTOALL_TYPES)
#undef DEFINE_ALLTOALL_TYPES


/**
 * @brief Helper macro to define alltoallmem implementations
 *
 */
#define UCC_ALLTOALLMEM_DEFINITION(_placeholder)                                        \
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
                                                                               \
    ucc_alltoall_helper(                                                       \
        dest, source, nelems, team_h->start,                                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

UCC_ALLTOALLMEM_DEFINITION(placeholder)


 /* @brief Helper macro to define SIZE alltoall implementations
 *
 * @param _size Size in bits
 */
#define UCC_ALLTOALL_SIZE_DEFINITION(_size)                                 \
                                                                               \
  void ucc_alltoall##_size (void *dest, const void *source, size_t nelems, int PE_start, int logPE_stride, int PE_size, long *pSync)  \
  {                                                                            \
    /* Sanity checks */                                                        \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_POSITIVE(PE_size, "PE_size");                                 \
    SHMEMU_CHECK_NON_NEGATIVE(PE_start, "PE_start");                           \
    SHMEMU_CHECK_NON_NEGATIVE(logPE_stride, "logPE_stride");                   \
    SHMEMU_CHECK_ACTIVE_SET_RANGE(PE_start, logPE_stride, PE_size);            \
    SHMEMU_CHECK_SYMMETRIC(dest, (_size) / (CHAR_BIT) * nelems * PE_size);     \
    SHMEMU_CHECK_SYMMETRIC(source, (_size) / (CHAR_BIT) * nelems * PE_size);   \
    SHMEMU_CHECK_SYMMETRIC(pSync, sizeof(long) * SHCOLL_ALLTOALL_SYNC_SIZE);   \
    SHMEMU_CHECK_BUFFER_OVERLAP(dest, source,                                  \
                                (_size) / (CHAR_BIT) * nelems * PE_size,       \
                                (_size) / (CHAR_BIT) * nelems * PE_size);      \
    /* Perform alltoall */                                                     \
    ucc_alltoall_helper(dest, source, (_size) / (CHAR_BIT) * nelems,           \
                            PE_start, logPE_stride, PE_size, pSync);           \
  }

// @formatter:off
UCC_ALLTOALL_SIZE_DEFINITION(32)
UCC_ALLTOALL_SIZE_DEFINITION(64)
