/**
 * @file env.c
 * @brief Environment variable handling and initialization
 *
 * This file contains functions for reading and handling environment variables
 * that control OpenSHMEM runtime behavior, as well as initialization and
 * finalization of environment settings.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include "state.h"
#include "shmemu.h"
#include "shmemc.h"
#include "boolean.h"
#include "collectives/defaults.h"
#include "module.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>

/**
 * @brief Buffer size for string formatting
 */
#define BUFSIZE 16

/**
 * @brief Test if an environment variable option is enabled
 *
 * @param str String value to test
 * @return true if option is enabled, false otherwise
 *
 * Tests if an environment variable option is enabled by checking for:
 * - First letter 'y' or 'Y'
 * - String "on" (case insensitive)
 * - Non-zero numeric value
 */
static bool option_enabled_test(const char *str) {
  if (str == NULL) {
    return false;
    /* NOT REACHED */
  }
  if (tolower(*str) == 'y') {
    return true;
    /* NOT REACHED */
  }
  if (strncasecmp(str, "on", 2) == 0) {
    return true;
    /* NOT REACHED */
  }
  if (strtol(str, NULL, 10) != 0L) {
    return true;
    /* NOT REACHED */
  }
  return false;
}

/**
 * @brief Check for environment variable with SHMEM_ prefix
 */
#define CHECK_ENV(_e, _name)                                                   \
  do {                                                                         \
    (_e) = getenv("SHMEM_" #_name);                                            \
  } while (0)

/**
 * @brief Check for environment variable with SHMEM_ prefix, falling back to
 * SMA_ prefix
 */
#define CHECK_ENV_WITH_DEPRECATION(_e, _name)                                  \
  do {                                                                         \
    CHECK_ENV(_e, _name);                                                      \
    if ((_e) == NULL) {                                                        \
      (_e) = getenv("SMA_" #_name);                                            \
    }                                                                          \
  } while (0)

/**
 * @brief Initialize environment settings from environment variables
 *
 * Reads environment variables and initializes runtime settings including:
 * - Debug/logging options
 * - Symmetric heap size
 * - Collective algorithm selections
 * - Progress thread settings
 */
void shmemc_env_init(void) {
  char *e;
  int r;
  char *delay;

  /*
   * defined in spec
   */

  proc.env.print_version = false;
  proc.env.print_info = false;
  proc.env.debug = false;

  CHECK_ENV_WITH_DEPRECATION(e, VERSION);
  if (e != NULL) {
    proc.env.print_version = option_enabled_test(e);
  }
  CHECK_ENV_WITH_DEPRECATION(e, INFO);
  if (e != NULL) {
    proc.env.print_info = option_enabled_test(e);
  }
  CHECK_ENV_WITH_DEPRECATION(e, DEBUG);
  if (e != NULL) {
    proc.env.debug = option_enabled_test(e);
  }
  CHECK_ENV_WITH_DEPRECATION(e, SYMMETRIC_SIZE);
  proc.env.heap_spec =
      strdup(e != NULL ? e : SHMEM_DEFAULT_HEAP_SIZE); /* free@end */

  /*
   * this implementation also has...
   */

  proc.env.logging = false;
  proc.env.logging_events = NULL;
  proc.env.logging_file = NULL;

  CHECK_ENV(e, LOGGING);
  if (e != NULL) {
    proc.env.logging = option_enabled_test(e);
  }
  CHECK_ENV(e, LOGGING_FILE);
  if (e != NULL) {
    proc.env.logging_file = strdup(e); /* free@end */
  }
  CHECK_ENV(e, LOGGING_EVENTS);
  if (e != NULL) {
    proc.env.logging_events = strdup(e); /* free@end */
  }

  /* Initialize all collective variables to NULL */
  proc.env.coll.barrier = NULL;
  proc.env.coll.barrier_all = NULL;
  proc.env.coll.sync = NULL;
  proc.env.coll.team_sync = NULL;
  proc.env.coll.sync_all = NULL;

  proc.env.coll.broadcast_type = NULL;
  proc.env.coll.broadcast_mem = NULL;
  proc.env.coll.broadcast_size = NULL;

  proc.env.coll.collect_type = NULL;
  proc.env.coll.collect_mem = NULL;
  proc.env.coll.collect_size = NULL;

  proc.env.coll.fcollect_type = NULL;
  proc.env.coll.fcollect_mem = NULL;
  proc.env.coll.fcollect_size = NULL;

  proc.env.coll.alltoall_type = NULL;
  proc.env.coll.alltoall_mem = NULL;
  proc.env.coll.alltoall_size = NULL;

  proc.env.coll.alltoalls_type = NULL;
  proc.env.coll.alltoalls_mem = NULL;
  proc.env.coll.alltoalls_size = NULL;

  /* Initialize all reduction variables to NULL */
  proc.env.coll.and_to_all = NULL;
  proc.env.coll.or_to_all = NULL;
  proc.env.coll.xor_to_all = NULL;
  proc.env.coll.max_to_all = NULL;
  proc.env.coll.min_to_all = NULL;
  proc.env.coll.sum_to_all = NULL;
  proc.env.coll.prod_to_all = NULL;

  proc.env.coll.and_reduce = NULL;
  proc.env.coll.or_reduce = NULL;
  proc.env.coll.xor_reduce = NULL;
  proc.env.coll.max_reduce = NULL;
  proc.env.coll.min_reduce = NULL;
  proc.env.coll.sum_reduce = NULL;
  proc.env.coll.prod_reduce = NULL;

  proc.env.coll.inscan = NULL;
  proc.env.coll.exscan = NULL;

  /* Initialize from environment variables with defaults */
  CHECK_ENV(e, BARRIER_ALGO);
  proc.env.coll.barrier = strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_BARRIER);

  CHECK_ENV(e, BARRIER_ALL_ALGO);
  proc.env.coll.barrier_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_BARRIER_ALL);

  CHECK_ENV(e, SYNC_ALGO);
  proc.env.coll.sync = strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_SYNC);

  CHECK_ENV(e, TEAM_SYNC_ALGO);
  proc.env.coll.team_sync = strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_SYNC);

  CHECK_ENV(e, SYNC_ALL_ALGO);
  proc.env.coll.sync_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_SYNC_ALL);

  CHECK_ENV(e, BROADCAST_ALGO);
  proc.env.coll.broadcast_type =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_BROADCAST);
  CHECK_ENV(e, BROADCASTMEM_ALGO);
  proc.env.coll.broadcast_mem =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_BROADCAST);

  CHECK_ENV(e, COLLECT_ALGO);
  proc.env.coll.collect_type =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_COLLECT);
  CHECK_ENV(e, COLLECTMEM_ALGO);
  proc.env.coll.collect_mem =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_COLLECT);

  CHECK_ENV(e, FCOLLECT_ALGO);
  proc.env.coll.fcollect_type =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_FCOLLECT);
  CHECK_ENV(e, FCOLLECTMEM_ALGO);
  proc.env.coll.fcollect_mem =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_FCOLLECT);

  CHECK_ENV(e, ALLTOALL_ALGO);
  proc.env.coll.alltoall_type =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_ALLTOALL);
  CHECK_ENV(e, ALLTOALLMEM_ALGO);
  proc.env.coll.alltoall_mem =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_ALLTOALL);

  CHECK_ENV(e, ALLTOALLS_ALGO);
  proc.env.coll.alltoalls_type =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_ALLTOALLS);
  CHECK_ENV(e, ALLTOALLSMEM_ALGO);
  proc.env.coll.alltoalls_mem =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_ALLTOALLS);

  /* Deprecated sized variants */
  CHECK_ENV(e, ALLTOALL_SIZE_ALGO);
  proc.env.coll.alltoall_size =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_ALLTOALL);
  CHECK_ENV(e, ALLTOALLS_SIZE_ALGO);
  proc.env.coll.alltoalls_size =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_ALLTOALLS);
  CHECK_ENV(e, COLLECT_SIZE_ALGO);
  proc.env.coll.collect_size =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_COLLECT);
  CHECK_ENV(e, FCOLLECT_SIZE_ALGO);
  proc.env.coll.fcollect_size =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_FCOLLECT);
  CHECK_ENV(e, BROADCAST_SIZE_ALGO);
  proc.env.coll.broadcast_size =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_BROADCAST);

  /* Check for individual reduction algorithms or use defaults */
  CHECK_ENV(e, AND_TO_ALL_ALGO);
  proc.env.coll.and_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_AND_TO_ALL);

  CHECK_ENV(e, OR_TO_ALL_ALGO);
  proc.env.coll.or_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_OR_TO_ALL);

  CHECK_ENV(e, XOR_TO_ALL_ALGO);
  proc.env.coll.xor_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_XOR_TO_ALL);

  CHECK_ENV(e, MAX_TO_ALL_ALGO);
  proc.env.coll.max_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_MAX_TO_ALL);

  CHECK_ENV(e, MIN_TO_ALL_ALGO);
  proc.env.coll.min_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_MIN_TO_ALL);

  CHECK_ENV(e, SUM_TO_ALL_ALGO);
  proc.env.coll.sum_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_SUM_TO_ALL);

  CHECK_ENV(e, PROD_TO_ALL_ALGO);
  proc.env.coll.prod_to_all =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_PROD_TO_ALL);

  /* Check for team-based reduction algorithms or use defaults */
  CHECK_ENV(e, AND_REDUCE_ALGO);
  proc.env.coll.and_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_AND_REDUCE);

  CHECK_ENV(e, OR_REDUCE_ALGO);
  proc.env.coll.or_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_OR_REDUCE);

  CHECK_ENV(e, XOR_REDUCE_ALGO);
  proc.env.coll.xor_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_XOR_REDUCE);

  CHECK_ENV(e, MAX_REDUCE_ALGO);
  proc.env.coll.max_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_MAX_REDUCE);

  CHECK_ENV(e, MIN_REDUCE_ALGO);
  proc.env.coll.min_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_MIN_REDUCE);

  CHECK_ENV(e, SUM_REDUCE_ALGO);
  proc.env.coll.sum_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_SUM_REDUCE);

  CHECK_ENV(e, PROD_REDUCE_ALGO);
  proc.env.coll.prod_reduce =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_PROD_REDUCE);

  CHECK_ENV(e, INSCAN_ALGO);
  proc.env.coll.inscan =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_INSCAN);

  CHECK_ENV(e, EXSCAN_ALGO);
  proc.env.coll.exscan =
      strdup((e != NULL) ? e : COLLECTIVES_DEFAULT_EXSCAN);

  proc.env.progress_threads = NULL;

  CHECK_ENV(e, PROGRESS_THREADS);
  if (e != NULL) {
    proc.env.progress_threads = strdup(e); /* free@end */
  }

  delay = "1000"; /* magic number */
  proc.env.progress_delay_ns = strtol(delay, NULL, 10);

  CHECK_ENV(e, PROGRESS_DELAY);

  r = shmemu_parse_size(e != NULL ? e : delay, &proc.env.progress_delay_ns);
  shmemu_assert(r == 0,
                MODULE ": couldn't work out requested "
                       "progress delay time \"%s\"",
                e != NULL ? e : delay);

  proc.env.prealloc_contexts = 64; /* magic number */

  CHECK_ENV(e, PREALLOC_CTXS);
  if (e != NULL) {
    long n = strtol(e, NULL, 10);

    if (n < 0) {
      n = proc.env.prealloc_contexts;
    }
    proc.env.prealloc_contexts = (size_t)n;
  }

  proc.env.memfatal = true;

  CHECK_ENV(e, MEMERR_FATAL);
  if (e != NULL) {
    proc.env.memfatal = option_enabled_test(e);
  }
}

#undef CHECK_ENV
#undef CHECK_ENV_WITH_DEPRECATION

/**
 * @brief Clean up environment settings
 *
 * Frees memory allocated for environment variable strings
 */
void shmemc_env_finalize(void) {
  free(proc.env.logging_file);
  free(proc.env.logging_events);
  free(proc.env.heap_spec);

  /* Free all collective algorithm strings */
  free(proc.env.coll.alltoalls_size);
  free(proc.env.coll.alltoalls_type);
  free(proc.env.coll.alltoalls_mem);

  free(proc.env.coll.alltoall_size);
  free(proc.env.coll.alltoall_type);
  free(proc.env.coll.alltoall_mem);

  free(proc.env.coll.fcollect_size);
  free(proc.env.coll.fcollect_type);
  free(proc.env.coll.fcollect_mem);

  free(proc.env.coll.collect_size);
  free(proc.env.coll.collect_type);
  free(proc.env.coll.collect_mem);

  free(proc.env.coll.broadcast_size);
  free(proc.env.coll.broadcast_type);
  free(proc.env.coll.broadcast_mem);

  free(proc.env.coll.sync_all);
  free(proc.env.coll.sync);
  free(proc.env.coll.team_sync);
  free(proc.env.coll.barrier_all);
  free(proc.env.coll.barrier);

  free(proc.env.progress_threads);

  /* Free reduction operation fields */
  free(proc.env.coll.and_to_all);
  free(proc.env.coll.or_to_all);
  free(proc.env.coll.xor_to_all);
  free(proc.env.coll.max_to_all);
  free(proc.env.coll.min_to_all);
  free(proc.env.coll.sum_to_all);
  free(proc.env.coll.prod_to_all);

  free(proc.env.coll.and_reduce);
  free(proc.env.coll.or_reduce);
  free(proc.env.coll.xor_reduce);
  free(proc.env.coll.max_reduce);
  free(proc.env.coll.min_reduce);
  free(proc.env.coll.sum_reduce);
  free(proc.env.coll.prod_reduce);

  free(proc.env.coll.inscan);
  free(proc.env.coll.exscan);
}

/**
 * @brief Constants for formatting environment variable output
 */
static const int var_width = 22;
static const int val_width = 10;
static const int hr_width = 74;

/**
 * @brief Print a horizontal rule to the output stream
 *
 * @param stream Output stream
 * @param prefix String prefix for the line
 */
inline static void hr(FILE *stream, const char *prefix) {
  int i;

  fprintf(stream, "%s", prefix);
  for (i = 0; i < hr_width; ++i) {
    fprintf(stream, "-");
  }
  fprintf(stream, "\n");
}

/**
 * @brief Print environment variable settings
 *
 * @param stream Output stream to print to
 * @param prefix String prefix for each line
 *
 * Prints current values of all environment variables and their descriptions
 */
void shmemc_print_env_vars(FILE *stream, const char *prefix) {
  fprintf(stream,
          "%sEnvironment Variable Information.  "
          "See oshrun(1) for more.\n",
          prefix);
  fprintf(stream, "%s\n", prefix);
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "Variable", val_width,
          "Value", "Description");
  hr(stream, prefix);
  fprintf(stream, "%s\n", prefix);
  fprintf(stream, "%s%s\n", prefix, "From specification:");
  fprintf(stream, "%s\n", prefix);
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_VERSION",
          val_width, shmemu_human_option(proc.env.print_version),
          "print library version at start-up");
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_INFO",
          val_width, shmemu_human_option(proc.env.print_info),
          "print this information");
  {
    char buf[BUFSIZE];

    /* TODO hardwired index */
    (void)shmemu_human_number(proc.heaps.heapsize[0], buf, BUFSIZE);
    fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width,
            "SHMEM_SYMMETRIC_SIZE", val_width, buf,
            "requested size of the symmetric heap");
  }
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_DEBUG",
          val_width, shmemu_human_option(proc.env.debug),
          "enable sanity checking ("
#if !defined(ENABLE_DEBUG)
          "not "
#endif /* ! ENABLE_DEBUG */
          "configured)");

  fprintf(stream, "%s\n", prefix);
  fprintf(stream, "%s%s\n", prefix, "Specific to this implementation:");
  fprintf(stream, "%s\n", prefix);
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_LOGGING",
          val_width, shmemu_human_option(proc.env.logging),
          "enable logging messages ("
#if !defined(ENABLE_LOGGING)
          "not "
#endif /* ! ENABLE_LOGGING */
          "configured)");

  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_LOGGING_EVENTS",
          val_width, "...", /* could be far too long to show */
          "types of logging events to show");
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_LOGGING_FILE",
          val_width, proc.env.logging_file ? proc.env.logging_file : "unset",
          "file for logging information");

#define DESCRIBE_COLLECTIVE(_name, _envvar)                                    \
  do {                                                                         \
    fprintf(stream, "%s%-*s %-*s: %s\n", prefix, var_width,                    \
            "SHMEM_" #_envvar "_ALGO", val_width,                              \
            proc.env.coll._name ? proc.env.coll._name : "unset",               \
            "algorithm for \"" #_name "\" routine");                           \
  } while (0)

  DESCRIBE_COLLECTIVE(barrier, BARRIER);
  DESCRIBE_COLLECTIVE(barrier_all, BARRIER_ALL);
  DESCRIBE_COLLECTIVE(sync, SYNC);
  DESCRIBE_COLLECTIVE(team_sync, TEAM_SYNC);
  DESCRIBE_COLLECTIVE(sync_all, SYNC_ALL);

  DESCRIBE_COLLECTIVE(broadcast_type, BROADCAST_TYPE);
  DESCRIBE_COLLECTIVE(collect_type, COLLECT_TYPE);
  DESCRIBE_COLLECTIVE(fcollect_type, FCOLLECT_TYPE);
  DESCRIBE_COLLECTIVE(alltoall_type, ALLTOALL_TYPE);
  DESCRIBE_COLLECTIVE(alltoalls_type, ALLTOALLS_TYPE);

  DESCRIBE_COLLECTIVE(broadcast_mem, BROADCASTMEM);
  DESCRIBE_COLLECTIVE(collect_mem, COLLECTMEM);
  DESCRIBE_COLLECTIVE(fcollect_mem, FCOLLECTMEM);
  DESCRIBE_COLLECTIVE(alltoall_mem, ALLTOALLMEM);
  DESCRIBE_COLLECTIVE(alltoalls_mem, ALLTOALLSMEM);

  DESCRIBE_COLLECTIVE(broadcast_size, BROADCAST_SIZE);
  DESCRIBE_COLLECTIVE(collect_size, COLLECT_SIZE);
  DESCRIBE_COLLECTIVE(fcollect_size, FCOLLECT_SIZE);
  DESCRIBE_COLLECTIVE(alltoall_size, ALLTOALL_SIZE);
  DESCRIBE_COLLECTIVE(alltoalls_size, ALLTOALLS_SIZE);

  /* Reduction operations */
  DESCRIBE_COLLECTIVE(and_to_all, AND_TO_ALL);
  DESCRIBE_COLLECTIVE(or_to_all, OR_TO_ALL);
  DESCRIBE_COLLECTIVE(xor_to_all, XOR_TO_ALL);
  DESCRIBE_COLLECTIVE(max_to_all, MAX_TO_ALL);
  DESCRIBE_COLLECTIVE(min_to_all, MIN_TO_ALL);
  DESCRIBE_COLLECTIVE(sum_to_all, SUM_TO_ALL);
  DESCRIBE_COLLECTIVE(prod_to_all, PROD_TO_ALL);

  /* Team-based reduction operations */
  DESCRIBE_COLLECTIVE(and_reduce, AND_REDUCE);
  DESCRIBE_COLLECTIVE(or_reduce, OR_REDUCE);
  DESCRIBE_COLLECTIVE(xor_reduce, XOR_REDUCE);
  DESCRIBE_COLLECTIVE(max_reduce, MAX_REDUCE);
  DESCRIBE_COLLECTIVE(min_reduce, MIN_REDUCE);
  DESCRIBE_COLLECTIVE(sum_reduce, SUM_REDUCE);
  DESCRIBE_COLLECTIVE(prod_reduce, PROD_REDUCE);

  DESCRIBE_COLLECTIVE(prod_reduce, INSCAN);
  DESCRIBE_COLLECTIVE(prod_reduce, EXSCAN);

  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width,
          "SHMEM_PROGRESS_THREADS", val_width,
          proc.env.progress_threads ? proc.env.progress_threads : "none",
          "PEs that need progress threads");
  fprintf(stream, "%s%-*s %-*lu %s", prefix, var_width, "SHMEM_PROGRESS_DELAY",
          val_width, (unsigned long)proc.env.progress_delay_ns,
          "delay between progress polls (ns)");
  if (proc.env.progress_threads == NULL) {
    fprintf(stream, " [not used]");
  }
  fprintf(stream, "\n");
  fprintf(stream, "%s%-*s %-*lu %s\n", prefix, var_width, "SHMEM_PREALLOC_CTXS",
          val_width, (unsigned long)proc.env.prealloc_contexts,
          "pre-allocate contexts at startup");
  fprintf(stream, "%s%-*s %-*s %s\n", prefix, var_width, "SHMEM_MEMERR_FATAL",
          val_width, proc.env.memfatal ? "yes" : "no",
          "abort if symmetric memory corruption");

  /* ---------------------------------------------------------------- */

  fprintf(stream, "%s\n", prefix);
  hr(stream, prefix);
  fprintf(stream, "%s\n", prefix);
  fprintf(stream, "\n");

  fflush(stream);
}
