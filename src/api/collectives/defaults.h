/* For license: see LICENSE file at top-level */

/**
 * @file defaults.h
 * @brief Default algorithm selections for OpenSHMEM collective operations
 *
 * This file defines the default algorithm implementations used for various
 * OpenSHMEM collective operations. These can be overridden at runtime through
 * environment variables.
 */

#ifndef _COLLECTIVES_DEFAULTS_H
#define _COLLECTIVES_DEFAULTS_H 1

/** Default algorithm for all-to-all operations */
#define COLLECTIVES_DEFAULT_ALLTOALL "shift_exchange_barrier"

/** Default algorithm for strided all-to-all operations */
#define COLLECTIVES_DEFAULT_ALLTOALLS "shift_exchange_barrier"

/** Default algorithm for barrier operations */
#define COLLECTIVES_DEFAULT_BARRIER "binomial_tree"

/** Default algorithm for barrier_all operations */
#define COLLECTIVES_DEFAULT_BARRIER_ALL "binomial_tree"

/** Default algorithm for sync operations */
#define COLLECTIVES_DEFAULT_SYNC "binomial_tree"

/** Default algorithm for sync_all operations */
#define COLLECTIVES_DEFAULT_SYNC_ALL "binomial_tree"

/** Default algorithm for broadcast operations */
#define COLLECTIVES_DEFAULT_BROADCAST "binomial_tree"

/** Default algorithm for collect operations */
#define COLLECTIVES_DEFAULT_COLLECT "bruck"

/** Default algorithm for fcollect operations */
#define COLLECTIVES_DEFAULT_FCOLLECT "bruck_inplace"

/** Default algorithm for reduction operations */
#define COLLECTIVES_DEFAULT_REDUCTIONS "rec_dbl"

/** Default algorithm for inscan operations */
#define COLLECTIVES_DEFAULT_INSCAN "linear"

/** Default algorithm for exscan operations */
#define COLLECTIVES_DEFAULT_EXSCAN "linear"

/** Default algorithm for and-to-all reductions */
#define COLLECTIVES_DEFAULT_AND_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for or-to-all reductions */
#define COLLECTIVES_DEFAULT_OR_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for xor-to-all reductions */
#define COLLECTIVES_DEFAULT_XOR_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for max-to-all reductions */
#define COLLECTIVES_DEFAULT_MAX_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for min-to-all reductions */
#define COLLECTIVES_DEFAULT_MIN_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for sum-to-all reductions */
#define COLLECTIVES_DEFAULT_SUM_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for product-to-all reductions */
#define COLLECTIVES_DEFAULT_PROD_TO_ALL COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for and-reduce operations */
#define COLLECTIVES_DEFAULT_AND_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for or-reduce operations */
#define COLLECTIVES_DEFAULT_OR_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for xor-reduce operations */
#define COLLECTIVES_DEFAULT_XOR_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for max-reduce operations */
#define COLLECTIVES_DEFAULT_MAX_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for min-reduce operations */
#define COLLECTIVES_DEFAULT_MIN_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for sum-reduce operations */
#define COLLECTIVES_DEFAULT_SUM_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

/** Default algorithm for product-reduce operations */
#define COLLECTIVES_DEFAULT_PROD_REDUCE COLLECTIVES_DEFAULT_REDUCTIONS

#endif /* ! _COLLECTIVES_DEFAULTS_H */
