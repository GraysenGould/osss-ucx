
#include "ucc_coll.h"
#include <ucc/api/ucc.h>
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>

#include <stdio.h>
#include <math.h>

//generic: needs to take a variety of differnt collective types
// collectives where only size, not type, matters

#define UCC_COLLECTIVE_HELPER(_coll_type, _ucc_type_name)                     \
inline static void ucc_##_coll_type##_helper(                                 \
      void *dest, const void *source, size_t nelems, int PE_start,            \
      int logPE_stride, int PE_size, long *pSync,                             \
      ucc_team_h team_handle, ucc_context_h context_handle) {                 \
    printf("DEBUG: Part 8\n"); \
    printf("pe: %d, source: %p, dest: %p\n", shmem_my_pe(), (void *) source, (void *) dest); \
    ucc_context_attr_t context_attributes = { \
      .mask = UCC_CONTEXT_ATTR_FIELD_WORK_BUFFER_SIZE,  \
      .global_work_buffer_size = 1, \
    }; \
    ucc_context_get_attr(ucc_global_context, &context_attributes); \
    void * work_buffer = (void *) shmem_calloc(context_attributes.global_work_buffer_size * 2, 1); \
    printf("DEBUG: buffer size: %d\n", (int)context_attributes.global_work_buffer_size); \ 
    ucc_coll_buffer_info_t coll_src_buffer_info =  {                          \
      .buffer = source,                                                       \
      .count = nelems,                                              \
      .datatype = UCC_DT_UINT8,                                               \
      .mem_type = UCC_MEMORY_TYPE_HOST                                        \
    };                                                                        \
    ucc_coll_buffer_info_t coll_dst_buffer_info =  {                          \
      .buffer = dest,                                                         \
      .count = 2,                                              \
      .datatype = UCC_DT_UINT8,                                               \
      .mem_type = UCC_MEMORY_TYPE_HOST                                        \
    };                                                                        \
    ucc_coll_args_t coll_args = {                                             \
      .mask = UCC_COLL_ARGS_FIELD_FLAGS | UCC_COLL_ARGS_FIELD_GLOBAL_WORK_BUFFER,    \
      .flags = 0,                                                             \
      .coll_type = UCC_COLL_TYPE_ALLTOALL,                                            \
      .src.info = coll_src_buffer_info,                                       \
      .dst.info = coll_dst_buffer_info,                                       \
      .global_work_buffer = work_buffer,\
    };                                                                        \
    printf("DEBUG: Part 9\n");                            \
    ucc_coll_req_h coll_handle;                                               \
    ucc_collective_init(&coll_args, &coll_handle, team_handle);               \
    printf("DEBUG: Part 9.5\n");                            \
    ucc_collective_post(coll_handle);                                         \
    printf("DEBUG: Part 10\n");                            \
    /* poll operation until done */                                           \
    while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {               \
      ucc_context_progress(context_handle);                                   \
      /* This actually drives the communication */                            \
    }                                                                         \
    printf("DEBUG: Part 10\n");                            \
    ucc_collective_finalize(coll_handle);                                     \
    printf("DEBUG: Part 11\n");                            \
    shmem_free(work_buffer); \
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
  int ucc_##_typename##_alltoall(                                   \
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
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),                \
        global_ucc_team, ucc_global_context);                                \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \ }
#define DEFINE_ALLTOALL_TYPES(_type, _typename)                                \
  UCC_ALLTOALL_TYPE_DEFINITION(_type, _typename)                               \
  SHMEM_STANDARD_RMA_TYPE_TABLE(DEFINE_ALLTOALL_TYPES)
#undef DEFINE_ALLTOALL_TYPES


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
    printf("DEBUG: RUNING ALLTOALLHELPER");                                                                               \
    ucc_alltoall_helper(                                                       \
        dest, source, nelems, team_h->start,                                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),                \
            global_ucc_team, ucc_global_context);                            \
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
#define UCC_ALLTOALL_SIZE_DEFINITION(_size)                                    \
                                                                               \
  void ucc_alltoall##_size (void *dest, const void *source, size_t nelems,     \
        int PE_start, int logPE_stride, int PE_size, long *pSync)              \
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
                PE_start, logPE_stride, PE_size, pSync, NULL, NULL);           \
  }

// @formatter:off
UCC_ALLTOALL_SIZE_DEFINITION(32)
UCC_ALLTOALL_SIZE_DEFINITION(64)





ucc_team_h global_ucc_team;
ucc_lib_h ucc_lib;
ucc_context_h ucc_global_context;
ucc_team_h ucc_team_world;
ucc_shmem_oob_info_t global_oob_info;
ucc_mem_map_params_t ucc_global_mem_params;


ucc_status_t ucc_oob_all_gather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req)
{
    ucc_shmem_oob_info_t *info = (ucc_shmem_oob_info_t *)coll_info;
    
    // We need symmetric memory for OpenSHMEM OOB to work.
    // UCC provides 'sbuf' and 'rbuf' which are NOT symmetric.
    void *sym_sbuf = shmem_malloc(msglen);
    void *sym_rbuf = shmem_malloc(msglen * info->size);

    // 1. Copy UCC's private data to symmetric memory
    memcpy(sym_sbuf, sbuf, msglen);
    shmem_barrier_all();

    // 2. Perform the Put-based Allgather (as written before)
    for (int i = 0; i < info->size; i++) {
        void *dest = (char*)sym_rbuf + (info->rank * msglen);
        shmem_putmem(dest, sym_sbuf, msglen, i);
    }
    shmem_quiet();

    // 3. Signal completion
    for (int i = 0; i < info->size; i++) {
        shmem_atomic_inc(info->sync_counter, i);
    }

    // 4. Wait for completion locally
    while(shmem_atomic_fetch(info->sync_counter, info->rank) < info->size) {
        // poll
    }

    // 5. Copy data back to UCC's provided rbuf
    memcpy(rbuf, sym_rbuf, msglen * info->size);

    shmem_free(sym_sbuf);
    shmem_free(sym_rbuf);

    // Since we did this synchronously for the baseline, return a dummy req
    *req = (void*)0xDEADBEEF; 
    return UCC_OK;
}


ucc_status_t ucc_oob_allgather_test(void *req)
{
    return UCC_OK;
}

ucc_status_t ucc_oob_allgather_free(void *req)
{
    return UCC_OK;
}

void ucc_coll_init(){
   
    ucc_lib_params_t lib_params = {
            .mask = UCC_LIB_PARAM_FIELD_THREAD_MODE | UCC_LIB_PARAM_FIELD_SYNC_TYPE,
            .thread_mode = UCC_THREAD_SINGLE, /* will have to align with OpenSHMEM in future */
            .sync_type = UCC_NO_SYNC_COLLECTIVES
            };
    ucc_lib_config_h lib_config;
    
    ucc_lib_config_read(NULL, NULL, &lib_config);

    //start ucc
    ucc_init(&lib_params, lib_config, &ucc_lib);
}


void ucc_coll_finalize(){
    ucc_finalize(ucc_lib);
}

void ucc_coll_team_init (ucc_shmem_oob_info_t * oob_info, ucc_context_h * context_handle, ucc_team_h *team_handle) { 
  ucc_team_oob_coll_t team_oob_coll = {
    .allgather = ucc_oob_all_gather,
    .req_test  = ucc_oob_allgather_test,
    .req_free  = ucc_oob_allgather_free,
    .coll_info = (void*)oob_info, // Replace MPI_COMM_WORLD with this
    .n_oob_eps = oob_info->size,    // Corrected: Total count
    .oob_ep    = oob_info->rank     // Corrected: Local rank
  };
  /* Create Team Context */    
  uint32_t num_contexts = 1; // might have to change with multiple contexts 
  const ucc_team_params_t team_params = {
    .mask = UCC_TEAM_PARAM_FIELD_ORDERING | UCC_TEAM_PARAM_FIELD_OOB | UCC_TEAM_PARAM_FIELD_EP | UCC_TEAM_PARAM_FIELD_EP_RANGE,
    .ordering = UCC_COLLECTIVE_POST_UNORDERED, /* ordered might be needed for shmem_fence? */
    .oob = team_oob_coll, 
    .ep_range = UCC_COLLECTIVE_EP_RANGE_CONTIG, 
    .ep = (uint64_t) shmem_my_pe(), /* the endpoint is just the rank? */
  };

  ucc_status_t team_status = ucc_team_create_post(context_handle, num_contexts, &team_params, team_handle);
  if (team_status != UCC_OK || team_handle == NULL){ /* don't think this check is correct */
    printf("ERROR: could not create team\n");
    return;
  }

  while (UCC_INPROGRESS == ucc_team_create_test(*team_handle)) {
    ucc_context_progress(*context_handle); 
  }

}


void ucc_coll_team_finalize (ucc_team_h team_handle){
    ucc_team_destroy(team_handle);
}


void ucc_coll_context_create(){
  
  int rank = shmem_my_pe(), size = shmem_n_pes();
  long * my_symmetric_counter_ptr = (long *) shmem_malloc(sizeof(long));


  shmem_barrier_all();
  /* TODO: make context not global */ 
  global_oob_info.rank = rank;
  global_oob_info.size = size;
  global_oob_info.sync_counter = my_symmetric_counter_ptr;
  
  ucc_context_oob_coll_t context_oob_coll = {
    .allgather = ucc_oob_all_gather,
    .req_test  = ucc_oob_allgather_test,
    .req_free  = ucc_oob_allgather_free,
    .coll_info = (void*)&global_oob_info, // Replace MPI_COMM_WORLD with this
    .n_oob_eps = global_oob_info.size,    // Corrected: Total count
    .oob_ep    = global_oob_info.rank     // Corrected: Local rank
  };

  
  ucc_global_mem_params.segments = (ucc_mem_map_t *) shmem_malloc(sizeof(ucc_mem_map_t));

  ucc_global_mem_params.segments->address = shmem_malloc(1024 * size);
  
  ucc_global_mem_params.segments->len = 1024;

  ucc_global_mem_params.n_segments = 1;

  const ucc_context_params_t context_params = {
    .mask = UCC_CONTEXT_PARAM_FIELD_TYPE | UCC_CONTEXT_PARAM_FIELD_SYNC_TYPE | \
    UCC_CONTEXT_PARAM_FIELD_OOB | UCC_CONTEXT_PARAM_FIELD_MEM_PARAMS,
    .type = UCC_CONTEXT_SHARED,
    .sync_type = UCC_NO_SYNC_COLLECTIVES,
    .oob = context_oob_coll,
    .mem_params = ucc_global_mem_params
  };

  ucc_context_config_h context_config;

  ucc_context_config_read(ucc_lib, NULL, &context_config);

  shmem_barrier_all();

  ucc_status_t ctx_status;
  ctx_status = ucc_context_create(ucc_lib, &context_params, context_config, &ucc_global_context);
  if (ctx_status != UCC_OK){
    printf("ERROR: could not create ucc context\n");
    return;
  }

}

void ucc_coll_context_finalize(){
  /* TODO: make context not global */
  
  ucc_context_destroy(ucc_global_context); 
  shmem_free(global_oob_info.sync_counter);
  shmem_free(ucc_global_mem_params.segments);
  shmem_free(ucc_global_mem_params.segments->address);
}

