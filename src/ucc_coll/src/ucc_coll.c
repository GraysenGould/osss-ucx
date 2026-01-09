#include "ucc_coll.h"

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
   
    printf("DEBUG: part 1\n");
    ucc_lib_params_t lib_params = {
            .mask = UCC_LIB_PARAM_FIELD_THREAD_MODE | UCC_LIB_PARAM_FIELD_SYNC_TYPE,
            .thread_mode = UCC_THREAD_SINGLE, /* will have to align with OpenSHMEM in future */
            .sync_type = UCC_NO_SYNC_COLLECTIVES
            };
    ucc_lib_config_h lib_config;
    
    printf("DEBUG: part 2\n");
    ucc_lib_config_read(NULL, NULL, &lib_config);

    printf("DEBUG: part 3\n");
    //start ucc
    ucc_init(&lib_params, lib_config, &ucc_lib);
    printf("DEBUG: part 4\n");
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

  ucc_global_mem_params.segments->address = shmem_calloc(1024, size);
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
  shmem_free(ucc_global_mem_params.segments->address);
}

