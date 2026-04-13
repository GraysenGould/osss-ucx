/*
 * @file inscan.c
 * @brief Implementation of collective inscan operation
 *
 */

#include "shcoll.h"
#include <shmem/api_types.h>
#include "util/bithacks.h"
#include "../tests/util/debug.h"

#include "shmem.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>


/*
   @brief Helper macro to define linear inscan operations
 *
 * Implements a linear inscan algorithm where each PE gets and adds its required
 * information
 *
 * @param _name Typename of inscan operation
 * @param _type Data type to operate on
 */
#define INSCAN_HELPER_LINEAR(_name, _type)                                     \
  int inscan_helper_##_name##_linear(                                          \
    _type *dest, const _type *source, int nelems, int me_as, shmem_team_t team,\
    int PE_start, int logPE_stride, int PE_size, _type *pWrk, long *pSync) {   \
  _type * workBuffer = (_type *)calloc(nelems, sizeof(_type));                 \
  if (workBuffer == NULL)                                                      \
    return 1;                                                                  \
  /* If dest and source are the same, need a buffer to store intermediate */   \
  /* results */                                                                \
  _type * result;                                                              \
  if ( (void *) dest == (void *) source){                                      \
    result = (_type *)calloc(nelems, sizeof(_type));                           \
    if (result == NULL){                                                       \
      free(workBuffer);                                                        \
      return 1;                                                                \
    }                                                                          \
  }                                                                            \
  else {                                                                       \
    result = (_type *) dest;                                                   \
  }                                                                            \
                                                                               \
  shcoll_barrier_linear(PE_start, logPE_stride, PE_size, pSync);               \
  /* Do insum scan */                                                          \
  for (int pe = 0; pe <= me_as; pe++){                                         \
    int src_pe = shmemc_team_translate_pe(team, pe, SHMEM_TEAM_WORLD);         \
    shmem_get(workBuffer, source, nelems, src_pe);                             \
    /* add results */                                                          \
    for (int i = 0; i < nelems; i ++){                                         \
      result[i] = result[i] + workBuffer[i];                                   \
    }                                                                          \
  }                                                                            \
  shcoll_barrier_linear(PE_start, logPE_stride, PE_size, pSync);               \
                                                                               \
  /* If dest and source are same, override array with results */               \
  if ( (void *) dest == (void *) source){                                      \
    for (int i = 0; i < nelems; i ++)                                          \
      dest[i] = result[i];                                                     \
    free(result);                                                              \
  }                                                                            \
                                                                               \
  free(workBuffer);                                                            \
  return 0;                                                                    \
}

#define DECLARE_INSCAN_HELPER(_type, _typename)                               \
  INSCAN_HELPER_LINEAR(_typename, _type)

SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_INSCAN_HELPER)
#undef DECLARE_INSCAN_HELPER

/*
   @brief Helper macro to define linear exscan operations
 *
 * Implements a linear exscan algorithm where each PE gets and adds its required
 * information
 *
 * @param _name Typename of exscan operation
 * @param _type Data type to operate on
 */
#define EXSCAN_HELPER_LINEAR(_name, _type)                                     \
  int exscan_helper_##_name##_linear(                                          \
    _type *dest, const _type *source, int nelems, int me_as, shmem_team_t team,\
    int PE_start, int logPE_stride, int PE_size, _type *pWrk, long *pSync) {   \
  _type * workBuffer = (_type *)calloc(nelems, sizeof(_type));                 \
  if (workBuffer == NULL)                                                      \
    return 1;                                                                  \
  /* If dest and source are the same, need a buffer to store intermediate */   \
  /* results */                                                                \
  _type * result;                                                              \
  if ( (void *) dest == (void *) source){                                      \
    result = (_type *)calloc(nelems, sizeof(_type));                           \
    if (result == NULL){                                                       \
      free(workBuffer);                                                        \
      return 1;                                                                \
    }                                                                          \
  }                                                                            \
  else {                                                                       \
    result = (_type *) dest;                                                   \
  }                                                                            \
                                                                               \
  shcoll_barrier_linear(PE_start, logPE_stride, PE_size, pSync);               \
  /* Do insum scan */                                                          \
  for (int pe = 0; pe <= me_as - 1; pe++){                                     \
    int src_pe = shmemc_team_translate_pe(team, pe, SHMEM_TEAM_WORLD);         \
    shmem_get(workBuffer, source, nelems, src_pe);                             \
    /* add results */                                                          \
    for (int i = 0; i < nelems; i ++){                                         \
      result[i] = result[i] + workBuffer[i];                                   \
    }                                                                          \
  }                                                                            \
  shcoll_barrier_linear(PE_start, logPE_stride, PE_size, pSync);               \
                                                                               \
  /* If dest and source are same, override array with results */               \
  if ( (void *) dest == (void *) source){                                      \
    for (int i = 0; i < nelems; i ++)                                          \
      dest[i] = result[i];                                                     \
    free(result);                                                              \
  }                                                                            \
                                                                               \
  free(workBuffer);                                                            \
  return 0;                                                                    \
}

#define DECLARE_EXSCAN_HELPER(_type, _typename)                               \
  EXSCAN_HELPER_LINEAR(_typename, _type)

SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_EXSCAN_HELPER)
#undef DECLARE_EXSCAN_HELPER

/*
 * @brief Macro to define team-based inscan operations
 */
#define SHCOLL_INSCAN_DEFINITION(_typename, _type, _algo)                      \
  int shcoll_##_typename##_inscan_##_algo(                                     \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems) {    \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_SYMMETRIC(dest, "dest");                                      \
    SHMEMU_CHECK_SYMMETRIC(source, "source");                                  \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                     "team_h->pSyncs[COLLECTIVE]");                            \
                                                                               \
    _type *pWrk =                                                              \
        shmem_malloc(SHCOLL_REDUCE_MIN_WRKDATA_SIZE * sizeof(_type));          \
    int me = shmemc_my_pe();                                                   \
    int me_as = shmemc_team_translate_pe(SHMEM_TEAM_WORLD, me, team);          \
                                                                               \
    int success = inscan_helper_##_typename##_##_algo(                         \
        dest, source, nelems, me_as, team, team_h->start,                      \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks, pWrk,                                                  \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
    shmem_free(pWrk);                                                          \
    return success;                                                            \
  }

/*
 * @brief Macro to define team-based inscan operations
 */
#define SHCOLL_EXSCAN_DEFINITION(_typename, _type, _algo)                      \
  int shcoll_##_typename##_exscan_##_algo(                                     \
      shmem_team_t team, _type *dest, const _type *source, size_t nelems) {    \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_SYMMETRIC(dest, "dest");                                      \
    SHMEMU_CHECK_SYMMETRIC(source, "source");                                  \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_NULL(shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE),  \
                     "team_h->pSyncs[COLLECTIVE]");                            \
    /* TODO: check that either same src & dest, or non-overlapping */          \
    _type *pWrk =                                                              \
        shmem_malloc(SHCOLL_REDUCE_MIN_WRKDATA_SIZE * sizeof(_type));          \
    int me = shmemc_my_pe();                                                   \
    int me_as = shmemc_team_translate_pe(SHMEM_TEAM_WORLD, me, team);          \
                                                                               \
    int success = exscan_helper_##_typename##_##_algo(                         \
        dest, source, nelems, me_as, team, team_h->start,                      \
        (team_h->stride > 0) ? (int)log2((double)team_h->stride) : 0,          \
        team_h->nranks, pWrk,                                                  \
        shmemc_team_get_psync(team_h, SHMEMC_PSYNC_COLLECTIVE));               \
                                                                               \
    shmemc_team_reset_psync(team_h, SHMEMC_PSYNC_COLLECTIVE);                  \
    shmem_free(pWrk);                                                          \
    return success;                                                            \
  }

#define X(_type, _typename)                                       \
  SHCOLL_INSCAN_DEFINITION(_typename, _type, linear)              \
  SHCOLL_EXSCAN_DEFINITION(_typename, _type, linear)
  SHMEM_REDUCE_ARITH_TYPE_TABLE(X)
#undef X


