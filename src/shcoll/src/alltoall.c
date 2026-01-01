/**
 * @file alltoall.c
 * @brief Implementation of all-to-all collective operations
 *
 * This file implements various algorithms for all-to-all collective operations:
 * - Shift exchange
 * - XOR pairwise exchange
 * - Color pairwise exchange
 *
 * Each algorithm has variants using different synchronization methods:
 * - Barrier-based
 * - Signal-based
 * - Counter-based
 *
 * @copyright For license: see LICENSE file at top-level
 */

#include <shmem/api_types.h>
#include "shcoll.h"
#include "shcoll/compat.h"

#include <string.h>
#include <limits.h>
#include <assert.h>

#include <stdio.h>
#include <math.h>
#include <ucc/api/ucc.h>
#include <ucc/api/ucc_version.h>
#include <ucc/api/ucc_status.h>
#include <ucc/api/ucc_def.h>



/**
 * @brief Calculate edge color for color pairwise exchange algorithm
 *
 * @param i Current round number
 * @param me Current PE index
 * @param npes Total number of PEs
 * @return Edge color value
 */
inline static int edge_color(int i, int me, int npes) {
  int chr_idx;
  int v;

  chr_idx = npes % 2 == 1 ? npes : npes - 1;
  if (me < chr_idx) {
    v = (i + chr_idx - me) % chr_idx;
  } else {
    v = i % 2 == 1 ? (((i + chr_idx) / 2) % chr_idx) : i / 2;
  }

  if (npes % 2 == 1 && v == me) {
    return -1;
  } else if (v == me) {
    return chr_idx;
  } else {
    return v;
  }
}

/** @brief Number of rounds between synchronizations */
static int alltoall_rounds_sync = INT32_MAX;

/**
 * @brief Set number of rounds between synchronizations for alltoall operations
 * @param rounds_sync Number of rounds between synchronizations
 */
void shcoll_set_alltoall_round_sync(int rounds_sync) {
  alltoall_rounds_sync = rounds_sync;
}

/**
 * @brief Helper macro to define barrier-based alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _peer Function to calculate peer PE
 * @param _cond Condition that must be satisfied
 */
#define ALLTOALL_HELPER_BARRIER_DEFINITION(_algo, _peer, _cond)                \
  inline static void alltoall_helper_##_algo##_barrier(                        \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync) {                            \
    const int stride = 1 << logPE_stride;                                      \
    const int me = shmem_my_pe();                                              \
                                                                               \
    /* Get my index in the active set */                                       \
    const int me_as = (me - PE_start) / stride;                                \
                                                                               \
    void *const dest_ptr = ((uint8_t *)dest) + me_as * nelems;                 \
    void const *source_ptr = ((uint8_t *)source) + me_as * nelems;             \
                                                                               \
    int i;                                                                     \
    int peer_as;                                                               \
                                                                               \
    assert(_cond);                                                             \
                                                                               \
    memcpy(dest_ptr, source_ptr, nelems);                                      \
                                                                               \
    for (i = 1; i < PE_size; i++) {                                            \
      peer_as = _peer(i, me_as, PE_size);                                      \
      source_ptr = ((uint8_t *)source) + peer_as * nelems;                     \
                                                                               \
      shmem_putmem_nbi(dest_ptr, source_ptr, nelems,                           \
                       PE_start + peer_as * stride);                           \
                                                                               \
      if (i % alltoall_rounds_sync == 0) {                                     \
        /* TODO: change to auto shcoll barrier */                              \
        shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync);  \
      }                                                                        \
    }                                                                          \
                                                                               \
    /* TODO: change to auto shcoll barrier */                                  \
    shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync);      \
  }

/**
 * @brief Helper macro to define counter-based alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _peer Function to calculate peer PE
 * @param _cond Condition that must be satisfied
 */
#define ALLTOALL_HELPER_COUNTER_DEFINITION(_algo, _peer, _cond)                \
  inline static void alltoall_helper_##_algo##_counter(                        \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync) {                            \
    const int stride = 1 << logPE_stride;                                      \
    const int me = shmem_my_pe();                                              \
                                                                               \
    /* Get my index in the active set */                                       \
    const int me_as = (me - PE_start) / stride;                                \
                                                                               \
    void *const dest_ptr = ((uint8_t *)dest) + me_as * nelems;                 \
    void const *source_ptr;                                                    \
                                                                               \
    int i;                                                                     \
    int peer_as;                                                               \
                                                                               \
    assert(_cond);                                                             \
                                                                               \
    for (i = 1; i < PE_size; i++) {                                            \
      peer_as = _peer(i, me_as, PE_size);                                      \
      source_ptr = ((uint8_t *)source) + peer_as * nelems;                     \
                                                                               \
      shmem_putmem_nbi(dest_ptr, source_ptr, nelems,                           \
                       PE_start + peer_as * stride);                           \
    }                                                                          \
                                                                               \
    source_ptr = ((uint8_t *)source) + me_as * nelems;                         \
    memcpy(dest_ptr, source_ptr, nelems);                                      \
                                                                               \
    shmem_fence();                                                             \
                                                                               \
    for (i = 1; i < PE_size; i++) {                                            \
      peer_as = _peer(i, me_as, PE_size);                                      \
      shmem_long_atomic_inc(pSync, PE_start + peer_as * stride);               \
    }                                                                          \
                                                                               \
    shmem_long_wait_until(pSync, SHMEM_CMP_EQ,                                 \
                          SHCOLL_SYNC_VALUE + PE_size - 1);                    \
    shmem_long_p(pSync, SHCOLL_SYNC_VALUE, me);                                \
  }

/**
 * @brief Helper macro to define signal-based alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _peer Function to calculate peer PE
 * @param _cond Condition that must be satisfied
 */
#define ALLTOALL_HELPER_SIGNAL_DEFINITION(_algo, _peer, _cond)                 \
  inline static void alltoall_helper_##_algo##_signal(                         \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync) {                            \
    const int stride = 1 << logPE_stride;                                      \
    const int me = shmem_my_pe();                                              \
                                                                               \
    /* Get my index in the active set */                                       \
    const int me_as = (me - PE_start) / stride;                                \
                                                                               \
    void *const dest_ptr = ((uint8_t *)dest) + me_as * nelems;                 \
    void const *source_ptr;                                                    \
                                                                               \
    assert(_cond);                                                             \
                                                                               \
    int i;                                                                     \
    int peer_as;                                                               \
                                                                               \
    for (i = 1; i < PE_size; i++) {                                            \
      peer_as = _peer(i, me_as, PE_size);                                      \
      source_ptr = ((uint8_t *)source) + peer_as * nelems;                     \
                                                                               \
      shmem_putmem_signal_nb(dest_ptr, source_ptr, nelems, pSync + i - 1,      \
                             SHCOLL_SYNC_VALUE + 1,                            \
                             PE_start + peer_as * stride, NULL);               \
    }                                                                          \
                                                                               \
    source_ptr = ((uint8_t *)source) + me_as * nelems;                         \
    memcpy(dest_ptr, source_ptr, nelems);                                      \
                                                                               \
    for (i = 1; i < PE_size; i++) {                                            \
      shmem_long_wait_until(pSync + i - 1, SHMEM_CMP_GT, SHCOLL_SYNC_VALUE);   \
      shmem_long_p(pSync + i - 1, SHCOLL_SYNC_VALUE, me);                      \
    }                                                                          \
  }

#include <shmem.h>
#include <ucc/api/ucc.h>
#include <stdlib.h>

typedef struct {
    int rank;
    int size;
    long *sync_counter; // Symmetric memory pointer
} ucc_shmem_oob_info_t;

typedef struct {
    ucc_shmem_oob_info_t *info;
    int                   is_done;
} oob_request_t;

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


inline static void alltoall_helper_xor_pairwise_exchange_barrier(
      void *dest, const void *source, size_t nelems, int PE_start,
      int logPE_stride, int PE_size, long *pSync) {
    ucc_lib_h lib;
    ucc_lib_params_t lib_params = {
            .mask = UCC_LIB_PARAM_FIELD_THREAD_MODE | UCC_LIB_PARAM_FIELD_SYNC_TYPE,
            .thread_mode = UCC_THREAD_SINGLE, /* will have to align with OpenSHMEM in future */
            .sync_type = UCC_NO_SYNC_COLLECTIVES
            };
    int64_t * source_test = (int64_t *) source;
    printf("Source array on pe %d: ", shmem_my_pe());
    for (int i = 0; i < 8; i ++){
      printf("%d ", source_test[i]);
    }
    printf("\n");
    printf("nelems: %d, PE_size: %d\n", nelems, PE_size);

    ucc_lib_config_h lib_config;
    ucc_lib_config_read(NULL, NULL, &lib_config);

    //start ucc
    ucc_init(&lib_params, lib_config, &lib);
    
    // initialize context
    int rank = shmem_my_pe(), size = shmem_n_pes();

    long * my_symmetric_counter_ptr = (long *) shmem_malloc(sizeof(long));

    shmem_barrier_all();

    static ucc_shmem_oob_info_t my_oob_info;
    my_oob_info.rank = rank;
    my_oob_info.size = size;
    my_oob_info.sync_counter = my_symmetric_counter_ptr;
    
    ucc_context_oob_coll_t context_oob_coll = {
      .allgather = ucc_oob_all_gather,
      .req_test  = ucc_oob_allgather_test,
      .req_free  = ucc_oob_allgather_free,
      .coll_info = (void*)&my_oob_info, // Replace MPI_COMM_WORLD with this
      .n_oob_eps = my_oob_info.size,    // Corrected: Total count
      .oob_ep    = my_oob_info.rank     // Corrected: Local rank
    };

    static ucc_mem_map_t map_segments[1];
    void * onesided_buff = shmem_calloc(1024, size);
    map_segments[0].address = onesided_buff;
    map_segments[0].len = 1024;

    ucc_mem_map_params_t mem_params = { 
      .segments = map_segments,
      .n_segments = 1
    };

    const ucc_context_params_t context_params = {
      .mask = UCC_CONTEXT_PARAM_FIELD_TYPE | UCC_CONTEXT_PARAM_FIELD_SYNC_TYPE | \
      UCC_CONTEXT_PARAM_FIELD_OOB | UCC_CONTEXT_PARAM_FIELD_MEM_PARAMS,
      .type = UCC_CONTEXT_SHARED,
      .sync_type = UCC_NO_SYNC_COLLECTIVES,
      .oob = context_oob_coll,
      .mem_params = mem_params,
    };

    ucc_context_config_h context_config;

  
    ucc_context_config_read(lib, NULL, &context_config);
    

    ucc_context_h context_handle;

    shmem_barrier_all();

    ucc_status_t ctx_status;
    ctx_status = ucc_context_create(lib, &context_params, context_config, &context_handle);
    if (ctx_status != UCC_OK){
      printf("ERROR: could not create ucc context\n");
      return;
    }

    ucc_team_oob_coll_t team_oob_coll = {
      .allgather = ucc_oob_all_gather,
      .req_test  = ucc_oob_allgather_test,
      .req_free  = ucc_oob_allgather_free,
      .coll_info = (void*)&my_oob_info, // Replace MPI_COMM_WORLD with this
      .n_oob_eps = my_oob_info.size,    // Corrected: Total count
      .oob_ep    = my_oob_info.rank     // Corrected: Local rank
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

    ucc_team_h team_handle; 
    ucc_status_t team_status = ucc_team_create_post(&context_handle, num_contexts, &team_params, &team_handle);
    if (team_status != UCC_OK || team_handle == NULL){
      printf("ERROR: could not create team\n");
      return;
    }

    while (UCC_INPROGRESS == ucc_team_create_test(team_handle)) {
      ucc_context_progress(context_handle); 
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
      .mask = 0, /* just the bare minimum */
      .coll_type = UCC_COLL_TYPE_ALLTOALL,
      .src.info = coll_src_buffer_info,
      .dst.info = coll_dst_buffer_info,
      };

    ucc_coll_req_h coll_handle;
    ucc_collective_init(&coll_args, &coll_handle, team_handle);
    ucc_collective_post(coll_handle);
    
    /* poll operation until done */
    while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {
      ucc_context_progress(context_handle); // This actually drives the communication
    }

    ucc_collective_finalize(coll_handle); 

    /* cleanup */
    ucc_team_destroy(team_handle);
    ucc_context_destroy(context_handle); 
    shmem_free(my_symmetric_counter_ptr);
    shmem_free(onesided_buff);
    ucc_finalize(lib);
}

/** @brief Peer calculation for shift exchange algorithm */
#define SHIFT_PEER(I, ME, NPES) (((ME) + (I)) % (NPES))
ALLTOALL_HELPER_BARRIER_DEFINITION(shift_exchange, SHIFT_PEER, 1)
ALLTOALL_HELPER_COUNTER_DEFINITION(shift_exchange, SHIFT_PEER, 1)
ALLTOALL_HELPER_SIGNAL_DEFINITION(shift_exchange, SHIFT_PEER,
                                  PE_size - 1 <= SHCOLL_ALLTOALL_SYNC_SIZE)

/** @brief Peer calculation for XOR exchange algorithm */
#define XOR_PEER(I, ME, NPES) ((I) ^ (ME))
#define XOR_COND (((PE_size - 1) & PE_size) == 0)

//ALLTOALL_HELPER_BARRIER_DEFINITION(xor_pairwise_exchange, XOR_PEER, XOR_COND)
ALLTOALL_HELPER_COUNTER_DEFINITION(xor_pairwise_exchange, XOR_PEER, XOR_COND)
ALLTOALL_HELPER_SIGNAL_DEFINITION(xor_pairwise_exchange, XOR_PEER,
                                  XOR_COND &&PE_size - 1 <=
                                      SHCOLL_ALLTOALL_SYNC_SIZE)

/** @brief Peer calculation for color exchange algorithm */
#define COLOR_PEER(I, ME, NPES) edge_color(I, ME, NPES)
#define COLOR_COND (PE_size % 2 == 0)

ALLTOALL_HELPER_BARRIER_DEFINITION(color_pairwise_exchange, COLOR_PEER,
                                   COLOR_COND)
ALLTOALL_HELPER_COUNTER_DEFINITION(color_pairwise_exchange, COLOR_PEER,
                                   COLOR_COND)
ALLTOALL_HELPER_SIGNAL_DEFINITION(color_pairwise_exchange, COLOR_PEER,
                                  (PE_size - 1 <= SHCOLL_ALLTOALL_SYNC_SIZE) &&
                                      COLOR_COND)

// @formatter:on

/**
 * @brief Helper macro to define SIZE alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _size Size in bits
 */
#define SHCOLL_ALLTOALL_SIZE_DEFINITION(_algo, _size)                          \
  void shcoll_alltoall##_size##_##_algo(                                       \
      void *dest, const void *source, size_t nelems, int PE_start,             \
      int logPE_stride, int PE_size, long *pSync) {                            \
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
    alltoall_helper_##_algo(dest, source, (_size) / (CHAR_BIT) * nelems,       \
                            PE_start, logPE_stride, PE_size, pSync);           \
  }

// @formatter:off

SHCOLL_ALLTOALL_SIZE_DEFINITION(shift_exchange_barrier, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(shift_exchange_barrier, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(shift_exchange_counter, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(shift_exchange_counter, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(shift_exchange_signal, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(shift_exchange_signal, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(xor_pairwise_exchange_barrier, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(xor_pairwise_exchange_barrier, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(xor_pairwise_exchange_counter, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(xor_pairwise_exchange_counter, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(xor_pairwise_exchange_signal, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(xor_pairwise_exchange_signal, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(color_pairwise_exchange_counter, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(color_pairwise_exchange_counter, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(color_pairwise_exchange_barrier, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(color_pairwise_exchange_barrier, 64)

SHCOLL_ALLTOALL_SIZE_DEFINITION(color_pairwise_exchange_signal, 32)
SHCOLL_ALLTOALL_SIZE_DEFINITION(color_pairwise_exchange_signal, 64)

// @formatter:on

/**
 * @brief Helper macro to define typed alltoall implementations
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 *
 * FIXME: THESE COLLECTIVES NEED TO RETURN NON ZERO IF THEY FAIL
 */
#define SHCOLL_ALLTOALL_TYPE_DEFINITION(_algo, _type, _typename)               \
  int shcoll_##_typename##_alltoall_##_algo(                                   \
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
    alltoall_helper_##_algo(                                                   \
        dest, source, nelems * sizeof(_type), team_h->start,                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

#define DEFINE_ALLTOALL_TYPES(_type, _typename)                                \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(shift_exchange_barrier, _type, _typename)    \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(shift_exchange_counter, _type, _typename)    \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(shift_exchange_signal, _type, _typename)     \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(xor_pairwise_exchange_barrier, _type,        \
                                  _typename)                                   \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(xor_pairwise_exchange_counter, _type,        \
                                  _typename)                                   \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(xor_pairwise_exchange_signal, _type,         \
                                  _typename)                                   \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(color_pairwise_exchange_barrier, _type,      \
                                  _typename)                                   \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(color_pairwise_exchange_counter, _type,      \
                                  _typename)                                   \
  SHCOLL_ALLTOALL_TYPE_DEFINITION(color_pairwise_exchange_signal, _type,       \
                                  _typename)

SHMEM_STANDARD_RMA_TYPE_TABLE(DEFINE_ALLTOALL_TYPES)
#undef DEFINE_ALLTOALL_TYPES

/**
 * @brief Helper macro to define alltoallmem implementations
 *
 * @param _algo Algorithm name
 */
#define SHCOLL_ALLTOALLMEM_DEFINITION(_algo)                                   \
  int shcoll_alltoallmem_##_algo(shmem_team_t team, void *dest,                \
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
    alltoall_helper_##_algo(                                                   \
        dest, source, nelems, team_h->start,                                   \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks,                                                        \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
                                                                               \
    return 0;                                                                  \
  }

SHCOLL_ALLTOALLMEM_DEFINITION(shift_exchange_barrier)
SHCOLL_ALLTOALLMEM_DEFINITION(shift_exchange_counter)
SHCOLL_ALLTOALLMEM_DEFINITION(shift_exchange_signal)
SHCOLL_ALLTOALLMEM_DEFINITION(xor_pairwise_exchange_barrier)
SHCOLL_ALLTOALLMEM_DEFINITION(xor_pairwise_exchange_counter)
SHCOLL_ALLTOALLMEM_DEFINITION(xor_pairwise_exchange_signal)
SHCOLL_ALLTOALLMEM_DEFINITION(color_pairwise_exchange_barrier)
SHCOLL_ALLTOALLMEM_DEFINITION(color_pairwise_exchange_counter)
SHCOLL_ALLTOALLMEM_DEFINITION(color_pairwise_exchange_signal)

// @formatter:on
