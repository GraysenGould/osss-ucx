/**
 * @file thispe.h
 * @brief Processing Element (PE) state and configuration definitions
 *
 * This header defines the core data structures used to track state and
 * configuration for each Processing Element (PE) in the OpenSHMEM runtime.
 * It includes definitions for collectives, PE status, environment settings,
 * PMI information, heap management, and overall PE state.
 *
 * @copyright See LICENSE file at top-level
 */

#ifndef _THISPE_H
#define _THISPE_H 1

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include "boolean.h"
#include "threading.h"

#include "ucx/pe.h"

#include <sys/types.h>

/*
 * -- General --------------------------------------------------------
 */

/**
 * @brief Collective operation configuration
 *
 * Holds string identifiers for different collective operations to allow
 * delegation to the upper collectives shim layer.
 */
typedef struct collectives {
  char *barrier_all; /**< Barrier across all PEs */
  char *sync;        /**< Synchronization */
  char *team_sync;   /**< Team synchronization */
  char *sync_all;    /**< Global synchronization */

  char *broadcast_type; /**< Broadcast operation type */
  char *broadcast_mem;  /**< Broadcast memory handling */
  char *broadcast_size; /**< Broadcast size handling */

  char *collect_type; /**< Collect operation type */
  char *collect_mem;  /**< Collect memory handling */
  char *collect_size; /**< Collect size handling */

  char *fcollect_type; /**< Fcollect operation type */
  char *fcollect_mem;  /**< Fcollect memory handling */
  char *fcollect_size; /**< Fcollect size handling */

  char *alltoall_type; /**< All-to-all operation type */
  char *alltoall_mem;  /**< All-to-all memory handling */
  char *alltoall_size; /**< All-to-all size handling */

  char *alltoalls_type; /**< Strided all-to-all type */
  char *alltoalls_mem;  /**< Strided all-to-all memory */
  char *alltoalls_size; /**< Strided all-to-all size */

  /* Individual reduction operations */
  char *and_to_all;  /**< Bitwise AND reduction */
  char *or_to_all;   /**< Bitwise OR reduction */
  char *xor_to_all;  /**< Bitwise XOR reduction */
  char *max_to_all;  /**< Maximum value reduction */
  char *min_to_all;  /**< Minimum value reduction */
  char *sum_to_all;  /**< Sum reduction */
  char *prod_to_all; /**< Product reduction */

  /* Individual team-based reduction operations */
  char *and_reduce;  /**< Team bitwise AND reduction */
  char *or_reduce;   /**< Team bitwise OR reduction */
  char *xor_reduce;  /**< Team bitwise XOR reduction */
  char *max_reduce;  /**< Team maximum reduction */
  char *min_reduce;  /**< Team minimum reduction */
  char *sum_reduce;  /**< Team sum reduction */
  char *prod_reduce; /**< Team product reduction */

  /* inscan and exscan algorithms */
  char *inscan;
  char *exscan;

  char *barrier; /**< Barrier operation */
} shmemc_coll_t;

/**
 * @brief Processing Element status values
 */
typedef enum shmemc_status {
  SHMEMC_PE_SHUTDOWN = 0, /**< PE has been shut down */
  SHMEMC_PE_RUNNING,      /**< PE is running normally */
  SHMEMC_PE_FAILED,       /**< PE has failed */
  SHMEMC_PE_UNKNOWN       /**< PE status is unknown */
} shmemc_status_t;

/**
 * @brief Environment configuration settings
 *
 * Holds configuration values that can be set through environment variables
 */
typedef struct env_info {
  /*
   * required
   */
  bool print_version; /**< produce info output? */
  bool print_info;    /**< show settings at program start? */
  bool debug;         /**< are we doing sanity debugging? */
  char *heap_spec;    /**< how user told us about heap size */

  /*
   * this implementation
   */
  bool logging;         /**< turn on message logging? */
  char *logging_file;   /**< where does logging output go? */
  char *logging_events; /**< show only these types of messages */

  shmemc_coll_t coll; /**< collectives */

  char *progress_threads;   /**< do we need to start our own? */
  size_t progress_delay_ns; /**< if progress needed, time (ns)
                               between polls */

  size_t prealloc_contexts; /**< set up this many at start */
  bool memfatal;            /**< force exit on memory usage error? */
} env_info_t;

/**
 * @brief Process Management Interface (PMI) information
 *
 * Contains physical values from the launch environment that are used
 * during team creation
 */
typedef struct pmi_info {
  int rank;     /**< per-PE physical rank info */
  int nranks;   /**< number of ranks requested */
  int maxranks; /**< universe size (e.g. for spares) */
  int nnodes;   /**< number of nodes allocated */
  int *peers;   /**< peer PEs in a node group */
  int npeers;   /**< how many peers? */
} pmi_info_t;

/**
 * @brief Symmetric heap information
 *
 * Describes the symmetric heap(s) managed by this PE
 */
typedef struct heapinfo {
  size_t nheaps;    /**< how many heaps requested */
  size_t *heapsize; /**< array of their sizes */
} heapinfo_t;

/**
 * @brief Core PE state information
 *
 * Main structure containing all state information for a Processing Element
 */
typedef struct thispe_info {
  comms_info_t comms;     /**< per-comms layer info */
  env_info_t env;         /**< environment vars */
  heapinfo_t heaps;       /**< the heaps we manage */
  thread_desc_t td;       /**< threading model invoked */
  pmi_info_t li;          /**< launcher-supplied info */
  shmemc_status_t status; /**< up, down, out to lunch etc */
  int refcount;           /**< library initialization count */
  bool leader;            /**< is this PE a node leader? */
  bool progress_thread;   /**< PE requests progress thread */
  char *nodename;         /**< node we're running on */
  shmemc_team_t *teams;   /**< PE teams we belong to */
  size_t nteams;          /**< how many teams? */
} thispe_info_t;

#endif /* ! _THISPE_H */
