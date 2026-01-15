
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


/* do it the old fashioned way */
inline static void ucc_alltoallmem_helper(
      void *dest, const void *source, size_t nelems, int PE_start,
      int logPE_stride, int PE_size, long *pSync) {
    ucc_lib_h            lib;
    ucc_lib_config_h     lib_config;

    ucc_lib_params_t lib_params = {
            .mask = UCC_LIB_PARAM_FIELD_THREAD_MODE,
            .thread_mode = UCC_THREAD_SINGLE,
            };

    if (UCC_OK != ucc_lib_config_read(NULL, NULL, &lib_config)){
      printf("UCC FAILURE: could not read library configuraion\n");
      return;
    }

    if (UCC_OK != ucc_init(&lib_params, lib_config, &lib)){
      printf("UCC FAILURE: Could not initialize ucc library\n");
      return;
    }
    ucc_lib_config_release(lib_config);


    ucc_lib_attr_t lib_attributes;
    lib_attributes.mask = UCC_LIB_ATTR_FIELD_THREAD_MODE |
                          UCC_LIB_ATTR_FIELD_COLL_TYPES;

    if (UCC_OK != ucc_lib_get_attr(lib, &lib_attributes)){
      printf("Error Getting Library Attributes\n");
      return ;
    }
    
    /* --------- SETUP LIBRARY CONTEXT  ---------*/

    ucc_context_params_t  context_params;
    ucc_context_h         context_handle;
    ucc_context_config_h  context_config;
    ucc_mem_map_t         *map_segments;
    ucc_status_t          status;

    long * my_symmetric_counter_ptr = (long *) shmem_malloc(sizeof(long));

    shmem_barrier_all();

    printf("DEBUG: Part 2\n");
    ucc_shmem_oob_info_t my_oob_info;
    my_oob_info.rank = shmem_my_pe();
    my_oob_info.size = shmem_n_pes();
    my_oob_info.sync_counter = my_symmetric_counter_ptr;
    
    printf("DEBUG: Part 3\n");
    /* Create Global Work buffer with plenty of space */
    void * global_work_buffer = (void *) shmem_calloc(1024, 1);

    /* Register necessary memory */
    int n_segments = 3;
    map_segments = (ucc_mem_map_t *) shmem_calloc(n_segments, sizeof(ucc_mem_map_t));
    map_segments[0].address = source;
    map_segments[0].len = nelems * PE_size;
    map_segments[1].address = dest;
    map_segments[1].len = nelems * PE_size;
    map_segments[2].address = global_work_buffer;
    map_segments[2].len = 1024;

    context_params.mask = UCC_CONTEXT_PARAM_FIELD_OOB | UCC_CONTEXT_PARAM_FIELD_MEM_PARAMS;
    context_params.oob.allgather         = ucc_oob_all_gather;
    context_params.oob.req_test          = ucc_oob_allgather_test;
    context_params.oob.req_free          = ucc_oob_allgather_free;
    context_params.oob.coll_info         = (void *)&my_oob_info;
    context_params.oob.n_oob_eps         = my_oob_info.size;
    context_params.oob.oob_ep            = my_oob_info.rank;
    context_params.mem_params.segments   = map_segments;
    context_params.mem_params.n_segments = n_segments;

    printf("DEBUG: Part 5\n");
    ucc_context_config_read(lib, NULL, &context_config);

    status = ucc_context_create(lib, &context_params, context_config, &context_handle);
    if (status != UCC_OK){
      printf("ERROR: could not create ucc context\n");
      return;
    }

    ucc_context_config_release(context_config);

    /* CREATE UCC TEAMS */
    ucc_team_oob_coll_t team_oob_coll = {
      .allgather = ucc_oob_all_gather,
      .req_test  = ucc_oob_allgather_test,
      .req_free  = ucc_oob_allgather_free,
      .coll_info = (void*)&my_oob_info,
      .n_oob_eps = my_oob_info.size,    
      .oob_ep    = my_oob_info.rank 
    };

    printf("DEBUG: Part 7\n");
    /* Create Team Context */   
    uint32_t num_contexts = 1; // might have to change with multiple contexts
    const ucc_team_params_t team_params = {
      .mask = UCC_TEAM_PARAM_FIELD_ORDERING | UCC_TEAM_PARAM_FIELD_OOB | UCC_TEAM_PARAM_FIELD_EP | UCC_TEAM_PARAM_FIELD_EP_RANGE,
      .ordering = UCC_COLLECTIVE_POST_UNORDERED, /* ordered might be needed for shmem_fence? */
      .oob = team_oob_coll, 
      .ep_range = UCC_COLLECTIVE_EP_RANGE_CONTIG, 
      .ep = (uint64_t) shmem_my_pe(), /* the endpoint is just the rank? */
    };

    ucc_team_h team_handle; 
    ucc_status_t team_status = ucc_team_create_post(&context_handle, num_contexts, &team_params, &team_handle);
    if (team_status != UCC_OK && team_status != UCC_INPROGRESS) {
      printf("ERROR: team_create_post failed: %d\n", team_status);
      return;
    }
    printf("DEBUG: Part 8\n");
    while (1) {
      team_status = ucc_team_create_test(team_handle);
      if (team_status == UCC_OK) {
        break; // team ready
      } else if (team_status < 0) {
        printf("ERROR: team_create_test failed: %d\n", team_status);
        return;
      }
      ucc_context_progress(context_handle); // drive progress
    }
    /* Finally do the collective operation: */

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
      .global_work_buffer = global_work_buffer,
    };
    
    printf("DEBUG: Part 9\n");
    ucc_coll_req_h coll_handle;
    ucc_status_t ucc_coll_status = ucc_collective_init(&coll_args, &coll_handle, team_handle);
    printf("UCC STAUTS after coll init: %d\n", ucc_coll_status);
    if (ucc_coll_status != UCC_OK){
      printf("Could Not Initalize UCC collective!!\n");
      return;
    }
    ucc_collective_post(coll_handle);
    
    printf("DEBUG: Part 10\n");
    /* poll operation until done */
    while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
      ucc_context_progress(context_handle); // This actually drives the communication
    }

    printf("DEBUG: Part 11\n");
    ucc_collective_finalize(coll_handle); 

    printf("DEBUG: Part 12\n");
    /* cleanup */
    ucc_team_destroy(team_handle);
    ucc_context_destroy(context_handle); 
    shmem_free(map_segments);
    shmem_free(global_work_buffer);
    shmem_free(my_symmetric_counter_ptr);
    ucc_finalize(lib);
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

