/* For license: see LICENSE file at top-level */

/**
 * @file init.c
 * @brief Implementation of OpenSHMEM initialization and finalization routines
 *
 * This file contains implementations of routines to initialize and finalize
 * the OpenSHMEM library, including thread initialization and cleanup.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include "shmemu.h"
#include "shmemc.h"
#include "state.h"
#include "info.h"
#include "threading.h"
#include "shmem_mutex.h"
#include "collectives/collectives.h"
#include "collectives/ucc/ucc.h"
#include "module.h"
#include "shmem/api.h"

#ifdef ENABLE_EXPERIMENTAL
#include "allocator/xmemalloc.h"
#endif /* ENABLE_EXPERIMENTAL */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef ENABLE_PSHMEM
#pragma weak shmem_init_thread = pshmem_init_thread
#define shmem_init_thread pshmem_init_thread
#pragma weak shmem_init = pshmem_init
#define shmem_init pshmem_init
#pragma weak shmem_finalize = pshmem_finalize
#define shmem_finalize pshmem_finalize
#endif /* ENABLE_PSHMEM */

/*
 * finish SHMEM portion of program, release resources
 */

/**
 * @brief Helper function to finalize the OpenSHMEM library
 *
 * This internal function handles the cleanup of OpenSHMEM resources,
 * including thread management, communications, collectives, and other
 * subsystems.
 */
static void finalize_helper(void) {
  
#ifdef HAVE_UCC
  shmem_ucc_coll_finalize();
#endif /* HAVE_UCC */

  threadwrap_thread_t this;

  /* do nothing if multiple finalizes */
  if (proc.refcount < 1) {
    return;
  }

  logger(LOG_FINALIZE, "%s()", __func__);

  this = threadwrap_thread_id();
  if (this != proc.td.invoking_thread) {

    logger(LOG_FINALIZE, "mis-match: thread %lu initialized, but %lu finalized",
           (unsigned long)proc.td.invoking_thread, (unsigned long)this);
  }

  /* implicit barrier on finalize */
  shmem_barrier_all();

  shmemu_progress_finalize();

  shmemc_finalize();

  collectives_finalize();
  shmemt_finalize();
  shmemu_finalize();

#ifdef ENABLE_EXPERIMENTAL
  shmemxa_finalize();
#endif /* ENABLE_EXPERIMENTAL */

  --proc.refcount;
  proc.status = SHMEMC_PE_SHUTDOWN;
}

/**
 * @brief Helper function to initialize the OpenSHMEM library with threading
 * support
 *
 * @param requested The requested threading level
 * @param provided Pointer to store the provided threading level
 * @return 0 on success, non-zero on failure
 *
 * This internal function handles the initialization of OpenSHMEM resources,
 * including communications, thread management, and other subsystems.
 */
inline static int init_thread_helper(int requested, int *provided) {
  int s;

  /* do nothing if multiple inits */
  if (proc.refcount > 0) {
    return 0;
  }


  /* set up comms, read environment */
  shmemc_init();
    

  /* utiltiies */ 
  shmemt_init(); 
  shmemu_init();
  collectives_init();

#ifdef ENABLE_ALIGNED_ADDRESSES
  shmemu_test_asr_mismatch();
#endif /* ENABLE_ALIGNED_ADDRESSES */

  shmemu_progress_init();

  /* save and return thread level */
#ifdef ENABLE_THREADS
  switch (requested) {
  case SHMEM_THREAD_SINGLE:
  case SHMEM_THREAD_FUNNELED:
  case SHMEM_THREAD_SERIALIZED:
  case SHMEM_THREAD_MULTIPLE:
    break; /* nothing to do for now */
  default:
    shmemu_fatal(MODULE ": unknown thread level %d requested", requested);
    /* NOT REACHED */
    break;
  }

  proc.td.osh_tl = requested;
#else
  proc.td.osh_tl = SHMEM_THREAD_SINGLE;
#endif /* ENABLE_THREADS */

  if (provided != NULL) {
    *provided = proc.td.osh_tl;
  }

  proc.td.invoking_thread = threadwrap_thread_id();

#ifdef ENABLE_EXPERIMENTAL
  shmemxa_init(proc.heaps.nheaps);
#endif /* ENABLE_EXPERIMENTAL */

  s = atexit(finalize_helper);
  if (s != 0) {
    shmemu_fatal(MODULE ": unable to register atexit() handler: %s",
                 strerror(errno));
    /* NOT REACHED */
  }

  proc.status = SHMEMC_PE_RUNNING;

  ++proc.refcount;

  if (shmemc_my_pe() == 0) {
    if (proc.env.print_version) {
      info_output_package_version(stdout, "# ", "", 0);
    }
    if (proc.env.print_info) {
      shmemc_print_env_vars(stdout, "# ");
    }
  }

  logger(LOG_INIT, "%s(requested=%s [%d], provided->%s [%d])", __func__,
         shmemu_thread_name(requested), requested,
         shmemu_thread_name(proc.td.osh_tl), proc.td.osh_tl);

  shmem_barrier_all();

  /* just declare success */
  return 0;
}

/*
 * initialize/finalize SHMEM portion of program with threading model
 */

/**
 * @brief Finalize the OpenSHMEM library
 *
 * This routine releases all resources used by the OpenSHMEM library.
 * It must be the last OpenSHMEM routine called in a program.
 */
void shmem_finalize(void) { finalize_helper(); }

/**
 * @brief Initialize the OpenSHMEM library with threading support
 *
 * @param requested The requested threading level
 * @param provided Pointer to store the provided threading level
 * @return 0 on success, non-zero on failure
 *
 * Initialize the OpenSHMEM library with a specified threading level.
 */
int shmem_init_thread(int requested, int *provided) {
  return init_thread_helper(requested, provided);
}

/**
 * @brief Initialize the OpenSHMEM library
 *
 * Initialize the OpenSHMEM library with single threading level.
 */
void shmem_init(void) { 
  (void)init_thread_helper(SHMEM_THREAD_SINGLE, NULL); 
#ifdef HAVE_UCC
  shmem_ucc_coll_setup();
#endif /* HAVE_UCC */
}

#ifdef PR470

/**
 * @brief Check if OpenSHMEM library is initialized
 *
 * @return Non-zero if initialized, 0 otherwise
 */
int shmem_initialized(void) { return proc.refcount > 0; }

/**
 * @brief Check if OpenSHMEM library is finalized
 *
 * @return Non-zero if finalized, 0 otherwise
 */
int shmem_finalized(void) { return proc.refcount < 1; }

#endif /* PR470 */
