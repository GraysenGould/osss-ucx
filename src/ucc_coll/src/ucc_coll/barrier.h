/**
 * @file barrier.h
 * @brief Header file for OpenSHMEM barrier and sync collective operations
 *
 * Declares functions for various barrier synchronization algorithms including:
 * - Linear barrier
 * - Complete tree barrier
 * - Binomial tree barrier
 * - K-nomial tree barrier
 * - Dissemination barrier
 */

#ifndef _SHCOLL_BARRIER_H
#define _SHCOLL_BARRIER_H 1

/**
 * @brief Set the tree degree for tree-based barrier algorithms
 * @param tree_degree The tree degree to use
 */
void shcoll_set_tree_degree(int tree_degree);

/**
 * @brief Set the radix for k-nomial tree barrier algorithm
 * @param tree_radix The tree radix to use
 */
void shcoll_set_knomial_tree_radix_barrier(int tree_radix);

/**
 * @brief Macro to declare barrier and sync functions for a given algorithm
 *
 * Declares four functions for each algorithm:
 * - barrier: Active set barrier with memory ordering
 * - barrier_all: Global barrier with memory ordering
 * - sync: Active set barrier without memory ordering
 * - sync_all: Global barrier without memory ordering
 *
 * @param _algo Algorithm name to generate declarations for
 */
#define SHCOLL_BARRIER_SYNC_DECLARATION(_algo)                                 \
  void shcoll_barrier_##_algo(int PE_start, int logPE_stride, int PE_size,     \
                              long *pSync);                                    \
                                                                               \
  void shcoll_barrier_all_##_algo(long *pSync);                                \
                                                                               \
  void shcoll_sync_##_algo(int PE_start, int logPE_stride, int PE_size,        \
                           long *pSync);                                       \
                                                                               \
  void shcoll_sync_all_##_algo(long *pSync);

SHCOLL_BARRIER_SYNC_DECLARATION(linear)
SHCOLL_BARRIER_SYNC_DECLARATION(complete_tree)
SHCOLL_BARRIER_SYNC_DECLARATION(binomial_tree)
SHCOLL_BARRIER_SYNC_DECLARATION(knomial_tree)
SHCOLL_BARRIER_SYNC_DECLARATION(dissemination)

/**
 * @brief Macro to declare team sync function for a given algorithm
 *
 * Declares a team-based synchronization function that allocates its own pSync
 * array.
 *
 * @param _algo Algorithm name to generate declaration for
 */
#define SHCOLL_TEAM_SYNC_DECLARATION(_algo)                                    \
  int shcoll_team_sync_##_algo(shmem_team_t team);

SHCOLL_TEAM_SYNC_DECLARATION(linear)
SHCOLL_TEAM_SYNC_DECLARATION(complete_tree)
SHCOLL_TEAM_SYNC_DECLARATION(binomial_tree)
SHCOLL_TEAM_SYNC_DECLARATION(knomial_tree)
SHCOLL_TEAM_SYNC_DECLARATION(dissemination)

#endif /* ! _SHCOLL_BARRIER_H */
