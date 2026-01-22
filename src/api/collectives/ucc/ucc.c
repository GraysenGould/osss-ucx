#include "ucc.h"
#include <ucc/api/ucc.h>
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h>
#include <assert.h>

#include <stdio.h>
#include <math.h>

/* Globally Visible Handle for UCC collectives */
ucc_coll_component_t shmem_ucc_coll;

/** @breif super simple barrier helper for ucc bootstrapping
 *
 */
/* void ucc_barrier_helper(long * pSync){  */
/*   const int me = shmem_my_pe(); */
/*   const int PE_size = shmem_n_pes(); */
/*   shmem_long_p(pSync, 0, me); */
/*   int pe; */
/*  */
/*   if (me == 0) { */
/*     #<{(| wait for the rest of the AS to poke me |)}># */
/*     shmem_long_wait_until(pSync, SHMEM_CMP_EQ,PE_size - 1); */
/*     #<{(| send acks out |)}># */
/*     for (pe = 1; pe < PE_size; ++pe) { */
/*       shmem_long_p(pSync, 1, pe); */
/*     } */
/*   } else { */
/*     shmem_long_atomic_inc(pSync, 0); */
/*     shmem_long_wait_until(pSync, SHMEM_CMP_NE, 0); */
/*   } */
/*   shmem_long_p(pSync, 0, me); */
/*   shmem_long_wait_until(pSync, SHMEM_CMP_EQ, 0); */
/* } */



void shmem_ucc_coll_setup () {

  ucc_lib_config_h     lib_config;
  ucc_lib_params_t     lib_params;
  shmem_ucc_coll.is_lib_initialized = 1;

  lib_params.mask = UCC_LIB_PARAM_FIELD_THREAD_MODE;
  lib_params.thread_mode = UCC_THREAD_SINGLE;

  if (UCC_OK != ucc_lib_config_read(NULL, NULL, &lib_config)){
    printf("UCC FAILURE: could not read library configuraion\n");
    return;
  }

  if (UCC_OK != ucc_init(&lib_params, lib_config, &shmem_ucc_coll.lib)){
    // printf("UCC FAILURE: Could not initialize ucc library\n");
    return;
  }
  ucc_lib_config_release(lib_config);


  ucc_lib_attr_t lib_attributes;
  lib_attributes.mask = UCC_LIB_ATTR_FIELD_THREAD_MODE |
                        UCC_LIB_ATTR_FIELD_COLL_TYPES;

  if (UCC_OK != ucc_lib_get_attr(shmem_ucc_coll.lib, &lib_attributes)){
    //printf("Error Getting Library Attributes\n");
    return ;
  }
  
  /* --------- SETUP LIBRARY CONTEXT  ---------*/

  ucc_context_params_t  context_params;
  ucc_context_config_h  context_config;
  ucc_status_t          status;

  // printf("DEBUG: Part 2\n");
  shmem_ucc_coll.oob_info.rank = shmem_my_pe();
  shmem_ucc_coll.oob_info.size = shmem_n_pes();
  
  // printf("DEBUG: Part 3\n");
  /* Create Global Work Buffer with plenty of space (x2 we need) */
  ucc_context_attr_t ctx_attr;
  ctx_attr.mask = UCC_CONTEXT_ATTR_FIELD_WORK_BUFFER_SIZE;
  int global_work_buffer_size = ctx_attr.global_work_buffer_size * sizeof(long) * 2;
  shmem_ucc_coll.global_work_buffer = (void *) shmem_calloc(global_work_buffer_size, 1);

  /* Register necessary memory */
  int n_segments = 1;
  shmem_ucc_coll.map_segments = (ucc_mem_map_t *) shmem_calloc(n_segments, sizeof(ucc_mem_map_t));
  shmem_ucc_coll.map_segments[0].address = shmem_ucc_coll.global_work_buffer;
  shmem_ucc_coll.map_segments[0].len = global_work_buffer_size;
  /* shmem_ucc_coll.map_segments[1].address = source; */
  /* shmem_ucc_coll.map_segments[1].len = nelems * PE_size; */
  /* shmem_ucc_coll.map_segments[2].address = dest; */
  /* shmem_ucc_coll.map_segments[2].len = nelems * PE_size; */

  context_params.mask = UCC_CONTEXT_PARAM_FIELD_OOB | UCC_CONTEXT_PARAM_FIELD_MEM_PARAMS;
  context_params.oob.allgather         = ucc_oob_allgather;
  context_params.oob.req_test          = ucc_oob_allgather_test;
  context_params.oob.req_free          = ucc_oob_allgather_free;
  context_params.oob.coll_info         = (void *)&shmem_ucc_coll.oob_info;
  context_params.oob.n_oob_eps         = shmem_ucc_coll.oob_info.size;
  context_params.oob.oob_ep            = shmem_ucc_coll.oob_info.rank;
  context_params.mem_params.segments   = shmem_ucc_coll.map_segments;
  context_params.mem_params.n_segments = n_segments;

  // printf("DEBUG: Part 5\n");
  ucc_context_config_read(shmem_ucc_coll.lib, NULL, &context_config);

  status = ucc_context_create(shmem_ucc_coll.lib, &context_params, context_config, &shmem_ucc_coll.context_handle);
  if (status != UCC_OK){
    printf("ERROR: could not create ucc context\n");
    return;
  }

  ucc_context_config_release(context_config);
  /* ------------ CREATE UCC TEAMS --------------- */
  
  ucc_team_oob_coll_t       team_oob_coll;
  team_oob_coll.allgather = ucc_oob_allgather;
  team_oob_coll.req_test  = ucc_oob_allgather_test;
  team_oob_coll.req_free  = ucc_oob_allgather_free;
  team_oob_coll.coll_info = (void*)&shmem_ucc_coll.oob_info;
  team_oob_coll.n_oob_eps = shmem_ucc_coll.oob_info.size;
  team_oob_coll.oob_ep    = shmem_ucc_coll.oob_info.rank;

  // TODO: use ucc_ep_map_t in the future
  /* printf("DEBUG: Part 7\n"); */
  /* Create Team Context */   
  const ucc_team_params_t team_params = {
    .mask = UCC_TEAM_PARAM_FIELD_OOB | UCC_TEAM_PARAM_FIELD_EP | UCC_TEAM_PARAM_FIELD_EP_RANGE | UCC_TEAM_PARAM_FIELD_FLAGS,
    .oob = team_oob_coll, 
    .ep_range = UCC_COLLECTIVE_EP_RANGE_CONTIG, 
    .ep = (uint64_t) shmem_ucc_coll.oob_info.rank,
    .flags  = UCC_TEAM_FLAG_COLL_WORK_BUFFER,      
  };

  uint32_t num_contexts = 1; // might have to change with multiple contexts
  if (UCC_OK != (status = ucc_team_create_post(&shmem_ucc_coll.context_handle, num_contexts, 
        &team_params, &shmem_ucc_coll.team_handle))){
    /* printf("ERROR: team_create_post failed: %d\n", status); */
    return;
  }

  //printf("DEBUG: Part 8\n");
  while (UCC_INPROGRESS == 
      (status = ucc_team_create_test(shmem_ucc_coll.team_handle))) {}
  if (UCC_OK != status) {
    printf("ERROR: ucc_team_create_test failed: %d\n", status);
    return;
  }
}

void shmem_ucc_coll_finalize(){
  if (shmem_ucc_coll.is_lib_initialized){
    if (UCC_OK != ucc_team_destroy(shmem_ucc_coll.team_handle)){
      printf("ERROR: could not destroy ucc team\n");
      return;
    }
    //printf("DEBUG: Part 13\n");
    if (UCC_OK != ucc_context_destroy(shmem_ucc_coll.context_handle)){
      printf("ERROR: Could not destroy ucc_context\n");
      return;
    }
    //printf("DEBUG: Part 14\n");
    shmem_free(shmem_ucc_coll.map_segments);
    //printf("DEBUG: Part 14.5\n");
    shmem_free(shmem_ucc_coll.global_work_buffer);
    //printf("DEBUG: Part 15\n");
    ucc_finalize(shmem_ucc_coll.lib);
    //printf("DEBUG: Part 16\n");
    shmem_ucc_coll.is_lib_initialized = 0;
  }
}

ucc_status_t ucc_oob_allgather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req)
{
  shmem_oob_info_t *info = (shmem_oob_info_t *)coll_info;   
  // We need symmetric memory for OpenSHMEM OOB to work.
  // UCC provides 'sbuf' and 'rbuf' which are NOT symmetric.
  void *sym_sbuf = shmem_malloc(msglen);
  void *sym_rbuf = shmem_malloc(msglen * info->size);
  void *pSync = shmem_calloc(sizeof(long), 1);
  //printf("sym_rbuf(%p) = msglen(%d) * info->size(%d)\n", sym_rbuf, msglen, info->size);

  // 1. Copy UCC's private data to symmetric memory
  memcpy(sym_sbuf, sbuf, msglen);
  shmem_barrier_all();  
  // 2. Perform the Put-based Allgather (as written before)
  for (int i = 0; i < info->size; i++) {
      void *dest = (char*)sym_rbuf + (info->rank * msglen);
      //printf("dest (%p) = sym_rbuf(%p) + info->rank(%d) * msglen(%d)\n", dest, sym_rbuf, info->rank, msglen);
      shmem_putmem(dest, sym_sbuf, msglen, i);
  }
  shmem_barrier_all();  
  // 5. Copy data back to UCC's provided rbuf
  memcpy(rbuf, sym_rbuf, msglen * info->size);

  shmem_free(sym_sbuf);
  shmem_free(sym_rbuf);
  shmem_free(pSync);

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
