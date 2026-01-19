
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

/* Globally Visible Handle for UCC collectives */
ucc_coll_component_t shmem_ucc_coll;

void shmem_ucc_coll_setup () {

  ucc_lib_config_h     lib_config;
  ucc_lib_params_t     lib_params;
  shmem_ucc_coll.is_lib_initialized = true;

  lib_params.mask = UCC_LIB_PARAM_FIELD_THREAD_MODE;
  lib_params.thread_mode = UCC_THREAD_SINGLE;

  if (UCC_OK != ucc_lib_config_read(NULL, NULL, &lib_config)){
    printf("UCC FAILURE: could not read library configuraion\n");
    return;
  }

  if (UCC_OK != ucc_init(&lib_params, lib_config, &shmem_ucc_coll.lib)){
    printf("UCC FAILURE: Could not initialize ucc library\n");
    return;
  }
  ucc_lib_config_release(lib_config);


  ucc_lib_attr_t lib_attributes;
  lib_attributes.mask = UCC_LIB_ATTR_FIELD_THREAD_MODE |
                        UCC_LIB_ATTR_FIELD_COLL_TYPES;

  if (UCC_OK != ucc_lib_get_attr(shmem_ucc_coll.lib, &lib_attributes)){
    printf("Error Getting Library Attributes\n");
    return ;
  }
  
  /* --------- SETUP LIBRARY CONTEXT  ---------*/

  ucc_context_params_t  context_params;
  ucc_context_config_h  context_config;
  ucc_status_t          status;

  shmem_barrier_all();

  printf("DEBUG: Part 2\n");
  shmem_ucc_coll.oob_info.rank = shmem_my_pe();
  shmem_ucc_coll.oob_info.size = shmem_n_pes();
  
  printf("DEBUG: Part 3\n");
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

  printf("DEBUG: Part 5\n");
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
  printf("DEBUG: Part 7\n");
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
    printf("ERROR: team_create_post failed: %d\n", status);
    return;
  }

  printf("DEBUG: Part 8\n");
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
    printf("DEBUG: Part 13\n");
    if (UCC_OK != ucc_context_destroy(shmem_ucc_coll.context_handle)){
      printf("ERROR: Could not destroy ucc_context\n");
      return;
    }
    printf("DEBUG: Part 14\n");
    shmem_free(shmem_ucc_coll.map_segments);
    printf("DEBUG: Part 14.5\n");
    shmem_free(shmem_ucc_coll.global_work_buffer);
    printf("DEBUG: Part 15\n");
    ucc_finalize(shmem_ucc_coll.lib);
    printf("DEBUG: Part 16\n");
    shmem_ucc_coll.is_lib_initialized = false;
  }
}

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
    .flags = 0,
  };
  // .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, // delete flag if want to skip memory pre-registration
  printf("DEBUG: Part 9\n");
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
  
  printf("DEBUG: Part 10\n");
  /* poll operation until done */
  while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
    /* Drive Collective Progress */
    ucc_context_progress(shmem_ucc_coll.context_handle);
  }

  printf("DEBUG: Part 11\n");
  ucc_collective_finalize(coll_handle); 
  printf("DEBUG: Part 12\n");
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



ucc_status_t ucc_oob_allgather(void *sbuf, void *rbuf, size_t msglen,
                                void *coll_info, void **req)
{
  shmem_oob_info_t *info = (shmem_oob_info_t *)coll_info;   
  // We need symmetric memory for OpenSHMEM OOB to work.
  // UCC provides 'sbuf' and 'rbuf' which are NOT symmetric.
  void *sym_sbuf = shmem_malloc(msglen);
  void *sym_rbuf = shmem_malloc(msglen * info->size);
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


/* ucc_status_t ucc_oob_all_gather(void *sbuf, void *rbuf, size_t msglen, */
/*                                 void *oob_info, void **req) */
/* {  */
/*   shmem_barrier_all(); */
/*   shmem_oob_allgather_info_t * allgather_info =  */
/*       (shmem_oob_allgather_info_t *) shmem_malloc(sizeof(shmem_oob_allgather_info_t)); */
/*   if (allgather_info == NULL){ */
/*     printf("ERROR: Could not allocate symmetric memory for UCC oob\n"); */
/*     return UCC_ERR_NO_MEMORY; */
/*   } */
/*  */
/*   allgather_info->oob_info = (shmem_oob_info_t *)oob_info;  */
/*   printf("DEBUG: OOB_ALLGATHER part 1\n"); */
/*   allgather_info->sbuf     = sbuf; */
/*   allgather_info->sym_sbuf = shmem_malloc(msglen); */
/*   if (allgather_info->sym_sbuf == NULL){ */
/*     shmem_free(allgather_info); */
/*     printf("ERROR: Could not allocate symmetric memory for UCC oob\n"); */
/*     return UCC_ERR_NO_MEMORY; */
/*   } */
/*  */
/*   printf("DEBUG: OOB_ALLGATHER part 2\n"); */
/*   allgather_info->sbuf     = sbuf; */
/*   allgather_info->rbuf     = rbuf; */
/*   printf("DEBUG: OOB_ALLGATHER part 2.1\n"); */
/*   allgather_info->sym_rbuf = shmem_malloc(msglen * allgather_info->oob_info->size); */
/*  */
/*   printf("DEBUG: OOB_ALLGATHER part 2.2\n"); */
/*   if (allgather_info->sym_rbuf == NULL){ */
/*     shmem_free(allgather_info->sym_sbuf); */
/*     shmem_free(allgather_info); */
/*     allgather_info = NULL; */
/*     printf("ERROR: Could not allocate symmetric memory for UCC oob\n"); */
/*     return UCC_ERR_NO_MEMORY; */
/*   } */
/*   printf("DEBUG: OOB_ALLGATHER part 3\n"); */
/*   allgather_info->sbuf     = sbuf; */
/*   allgather_info->msglen   = msglen; */
/*   allgather_info->send_counter = (int64_t *) shmem_malloc(sizeof(int64_t)); */
/*   if (allgather_info->send_counter == NULL){ */
/*     shmem_free(allgather_info->sym_sbuf); */
/*     shmem_free(allgather_info->sym_rbuf); */
/*     shmem_free(allgather_info); */
/*     allgather_info = NULL; */
/*     printf("ERROR: Could not allocate symmetric memory for UCC oob\n"); */
/*     return UCC_ERR_NO_MEMORY; */
/*   } */
/*  */
/*   printf("DEBUG: OOB_ALLGATHER part 4\n"); */
/*   allgather_info->sbuf     = sbuf; */
/*    */
/*   #<{(| Copy information from send buffer to symmetric send buffer |)}># */
/*   memcpy(allgather_info->sym_sbuf, allgather_info->sbuf, msglen); */
/*    */
/*   shmem_atomic_set(allgather_info->send_counter, 0, allgather_info->oob_info->rank); */
/*   #<{(| Send data to all other PE's |)}># */
/*   for (int i = 0; i < allgather_info->oob_info->size; i++) { */
/*       void *dest = (char*)allgather_info->sym_rbuf + (allgather_info->oob_info->rank * msglen); */
/*       shmem_putmem(dest, allgather_info->sym_sbuf, msglen, i); */
/*       #<{(| Create a separation |)}># */
/*       shmem_fence(); */
/*       #<{(| Increment dest PE counter, indicate that data has been transfered |)}># */
/*       shmem_atomic_inc(allgather_info->send_counter, i); */
/*   } */
/*  */
/*   printf("DEBUG: OOB_ALLGATHER part 5\n"); */
/*   *req = (void*)allgather_info; */
/*  */
/*   printf("DEBUG: OOB_ALLGATHER part 6\n"); */
/*   shmem_barrier_all(); */
/*   return UCC_OK; */
/* } */
/*  */
/* ucc_status_t ucc_oob_allgather_test(void *req) */
/* { */
/*   shmem_oob_allgather_info_t * allgather_info = (shmem_oob_allgather_info_t *) req; */
/*  */
/*   int64_t num_completed =  */
/*     shmem_atomic_fetch(allgather_info->send_counter, allgather_info->oob_info->rank); */
/*   printf("DEBUG: pe %d: testing progress: %d\n",allgather_info->oob_info->rank, (int)num_completed); */
/*   printf("DEBUG: ucc_oob_allgather_test part 1\n"); */
/*  */
/*   if ((int)num_completed == allgather_info->oob_info->size) */
/*     memcpy(allgather_info->rbuf, allgather_info->sym_rbuf, allgather_info->msglen * allgather_info->oob_info->size); */
/*     return UCC_OK; */
/*  */
/*   printf("DEBUG: ucc_oob_allgather_test part 2\n"); */
/*   return UCC_OK; */
/* } */
/*  */
/* ucc_status_t ucc_oob_allgather_free(void *req) */
/* {    */
/*   shmem_oob_allgather_info_t * allgather_info = (shmem_oob_allgather_info_t *) req; */
/*  */
/*   printf("DEBUG: ucc_oob_allgather_free part 1\n"); */
/*   if (allgather_info != NULL){ */
/*     shmem_free(allgather_info->sym_sbuf); */
/*     shmem_free(allgather_info->sym_rbuf); */
/*     shmem_free(allgather_info->send_counter); */
/*     shmem_free(allgather_info); */
/*   } */
/*  */
/*   printf("DEBUG: ucc_oob_allgather_free part 2\n"); */
/*   return UCC_OK; */
/* } */

