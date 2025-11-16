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

  inline static void alltoall_helper_xor_pairwise_exchange_barrier( // actually xor algorithm (logarithmic)
      void *dest, const void *source, size_t nelems, int PE_start,
      int logPE_stride, int PE_size, long *pSync) {
    const int stride = 1 << logPE_stride;
    const int me = shmem_my_pe();

    /* Get my index in the active set */
    const int me_as = (me - PE_start) / stride;
    void *dest_ptr;
    void *source_ptr;

    memcpy(dest, source, nelems * PE_size);
    int peer_as, distance, offset;

    for (int k = 0; (1 << k) < PE_size; k ++){
      distance = 1 << k;
      peer_as = me_as ^ distance;

      offset = (me_as & distance) ? 0 : 1;

      for (int i = 0; i < (PE_size >> (k + 1)); i++){
        dest_ptr = ((uint8_t *) dest) + (offset * distance * nelems) + i * 2 * distance * nelems; // really the source in this case
        source_ptr = ((uint8_t *) source) + i * 2 * distance * nelems;
        //printf("me_as: %d, k = %d, i = %d, dst offset: %d, src offset: %d\n", me_as, k, i, (offset * distance * nelems) + i * 2 * distance * nelems, i * 2 * distance * nelems);
        shmem_putmem_nbi(source_ptr, dest_ptr, nelems * distance, PE_start + peer_as * stride);
      }

      //shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync);
      shmem_quiet();

      for (int i = 0; i < (PE_size >> (k + 1)); i++){
        dest_ptr = ((uint8_t *) dest) + (offset * distance * nelems) + i * 2 * distance * nelems;;
        source_ptr = ((uint8_t *) source) + i * 2 * distance * nelems;
        memcpy(dest_ptr, source_ptr, nelems * distance);
      }
    }

    /* TODO: change to auto shcoll barrier */
    shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync);
  }

// inline static void alltoall_helper_xor_pairwise_exchange_barrier( // actually Bruck algorithm
//       void *dest, const void *source, size_t nelems, int PE_start,
//       int logPE_stride, int PE_size, long *pSync) {
//     const int stride = 1 << logPE_stride;
//     const int me = shmem_my_pe();

//     /* Get my index in the active set */
//     const int me_as = (me - PE_start) / stride;
//     void const *dest_ptr;
//     void const *source_ptr;

//     int peer_as, distance;
//     /* ensure power of 2 PE count */
//     assert(((PE_size - 1) & PE_size) == 0);
//     //memcpy(dest, source, nelems * PE_size);

//     for (int i = 0; i < PE_size; i ++){

//       source_ptr = ((uint8_t *) source) + ((me_as + i) % PE_size) * nelems;
//       dest_ptr = ((uint8_t *) dest) + i * nelems;
//       if (me_as == 1){
//         printf("%d coppied from %p to %p\n", *((int *) source_ptr), source_ptr, dest_ptr);
//       }
//       memcpy(dest_ptr, source_ptr, nelems);
//     }

//     if (me_as == 1){
//       printf("PE 1 start:");
//       for (int i = 0; i < PE_size; i ++){
//         printf("%d, ", *((int *)dest + i));
//       }
//       printf("\n\n\n");
//     }

//     for (int k = 0; (1 << k) < PE_size; k++) { // iterate log2(N) times - might not work from non 2^i number of processing elements
//       distance = (1 << k); // also the size of the sending?
//       peer_as = (me_as + distance) % PE_size; // according to Bruck algorithm

      
//       for (int i = 0; i < PE_size; i ++){
//         if (distance & i){
//           // send data
//           source_ptr = ((uint8_t *) source) + i * nelems;
//           dest_ptr = ((uint8_t *) dest) + i * nelems;

//           shmem_putmem_nbi(source_ptr, dest_ptr, nelems,
//                         PE_start + peer_as * stride);
//         }
//       }

//       if (me_as == 1){
//         printf("PE 0 round %d: ", k);
//         for (int i = 0; i < PE_size; i ++){
//           printf("%d, ", *((int *)dest + i));
//         }
//         printf("\n\n\n");
//       }

//       shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync); //barrier just for fun
      
//       for (int i = 0; i < PE_size; i ++){
//         if (distance && i){
//           // copy from source to dst array
//           source_ptr = ((uint8_t *) source) + i * nelems;
//           dest_ptr = ((uint8_t *) dest) + i * nelems;
//           memcpy(dest_ptr, source_ptr, nelems);
//         }
//       }

//       shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync); // barrier just for fun

//     }

//     for (int i = 0; i < PE_size; i ++){
//       source_ptr = ((uint8_t *) source) + i * nelems;
//       dest_ptr = ((uint8_t *) dest) + i * nelems;
//       memcpy(source_ptr, dest_ptr, nelems);
//     }

//     for (int i = 0; i < PE_size; i ++){
//       source_ptr = ((uint8_t *) source) + i * nelems;
//       dest_ptr = ((uint8_t *) dest) + (me_as + PE_size - i ) % PE_size;
//       memcpy(dest_ptr, source_ptr, nelems);
//     }

//     /* TODO: change to auto shcoll barrier */
//     shcoll_barrier_binomial_tree(PE_start, logPE_stride, PE_size, pSync);
//   }

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

// @formatter:off

/** @brief Peer calculation for shift exchange algorithm */
#define SHIFT_PEER(I, ME, NPES) (((ME) + (I)) % (NPES))
ALLTOALL_HELPER_BARRIER_DEFINITION(shift_exchange, SHIFT_PEER, 1)
ALLTOALL_HELPER_COUNTER_DEFINITION(shift_exchange, SHIFT_PEER, 1)
ALLTOALL_HELPER_SIGNAL_DEFINITION(shift_exchange, SHIFT_PEER,
                                  PE_size - 1 <= SHCOLL_ALLTOALL_SYNC_SIZE)

/** @brief Peer calculation for XOR exchange algorithm */
#define XOR_PEER(I, ME, NPES) ((I) ^ (ME))
#define XOR_COND (((PE_size - 1) & PE_size) == 0)

// ALLTOALL_HELPER_BARRIER_DEFINITION(xor_pairwise_exchange, XOR_PEER, XOR_COND)
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
