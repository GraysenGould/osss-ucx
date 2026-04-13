/**
 * @file table.h
 * @brief Header file defining collective operation tables and registration
 * functions
 *
 * This file contains declarations for the collective operation tables and
 * registration functions. It defines the structure for typed collective
 * operations, untyped collective operations, and typed to_all collective
 * operations.
 *
 * @copyright See LICENSE file at top-level
 */

#ifndef _TABLE_H
#define _TABLE_H 1

/** Maximum length for collective operation names */
#define COLL_NAME_MAX 64

/******************************************************** */
/** Function pointer type for collective operations without type information */
typedef void (*coll_fn_t)();

/**
 * @brief Structure for sized collective operations that have 32-bit and 64-bit
 * variants
 */
typedef struct sized_op {
  const char op[COLL_NAME_MAX]; /**< Operation name */
  coll_fn_t f32;                /**< 32-bit implementation */
  coll_fn_t f64;                /**< 64-bit implementation */
} sized_op_t;

/**
 * @brief Structure for collective operations without size variants
 */
typedef struct unsized_op {
  const char op[COLL_NAME_MAX]; /**< Operation name */
  coll_fn_t f;                  /**< Implementation function */
} unsized_op_t;

/******************************************************** */
/** Function pointer type for typed collective operations */
typedef int (*typed_coll_fn_t)();

/** Function pointer type for to_all typed collective operations (void return)
 */
typedef void (*typed_to_all_fn_t)();

/**
 * @brief Structure for typed collective operations
 */
typedef struct typed_op {
  const char op[COLL_NAME_MAX];   /**< Operation name */
  const char type[COLL_NAME_MAX]; /**< Type name */
  typed_coll_fn_t f;              /**< Implementation function */
} typed_op_t;

/**
 * @brief Structure for typed to_all collective operations
 */
typedef struct typed_to_all_op {
  const char op[COLL_NAME_MAX];   /**< Operation name */
  const char type[COLL_NAME_MAX]; /**< Type name */
  typed_to_all_fn_t f;            /**< Implementation function */
} typed_to_all_op_t;

/** Function pointer type for untyped collective operations */
typedef int (*untyped_coll_fn_t)();

/**
 * @brief Structure for untyped collective operations
 */
typedef struct untyped_op {
  const char op[COLL_NAME_MAX]; /**< Operation name */
  untyped_coll_fn_t f;          /**< Implementation function */
} untyped_op_t;

/******************************************************** */
/**
 * @brief Structure containing all collective operation implementations
 */
typedef struct coll_ops {
  /* Current routines */
  typed_op_t alltoall_type;  /**< Typed all-to-all operation */
  untyped_op_t alltoall_mem; /**< Generic all-to-all memory operation */
  sized_op_t alltoall_size;  /**< Sized all-to-all operation */

  typed_op_t alltoalls_type; /**< Typed strided all-to-all operation */
  untyped_op_t
      alltoalls_mem;         /**< Generic strided all-to-all memory operation */
  sized_op_t alltoalls_size; /**< Sized strided all-to-all operation */

  typed_op_t collect_type;  /**< Typed collect operation */
  untyped_op_t collect_mem; /**< Generic collect memory operation */
  sized_op_t collect_size;  /**< Sized collect operation */

  typed_op_t fcollect_type;  /**< Typed ordered collect operation */
  untyped_op_t fcollect_mem; /**< Generic ordered collect memory operation */
  sized_op_t fcollect_size;  /**< Sized ordered collect operation */

  typed_op_t broadcast_type;  /**< Typed broadcast operation */
  untyped_op_t broadcast_mem; /**< Generic broadcast memory operation */
  sized_op_t broadcast_size;  /**< Sized broadcast operation */

  typed_to_all_op_t and_to_all;  /**< Typed AND to all operation */
  typed_to_all_op_t or_to_all;   /**< Typed OR to all operation */
  typed_to_all_op_t xor_to_all;  /**< Typed XOR to all operation */
  typed_to_all_op_t max_to_all;  /**< Typed MAX to all operation */
  typed_to_all_op_t min_to_all;  /**< Typed MIN to all operation */
  typed_to_all_op_t sum_to_all;  /**< Typed SUM to all operation */
  typed_to_all_op_t prod_to_all; /**< Typed PROD to all operation */

  typed_op_t and_reduce;  /**< Typed AND reduce operation */
  typed_op_t or_reduce;   /**< Typed OR reduce operation */
  typed_op_t xor_reduce;  /**< Typed XOR reduce operation */
  typed_op_t max_reduce;  /**< Typed MAX reduce operation */
  typed_op_t min_reduce;  /**< Typed MIN reduce operation */
  typed_op_t sum_reduce;  /**< Typed SUM reduce operation */
  typed_op_t prod_reduce; /**< Typed PROD reduce operation */

  typed_op_t inscan; /**< Typed inscan operation */
  typed_op_t exscan; /**< Typed exscan operation */

  unsized_op_t barrier_all; /**< Typed global barrier operation */
  unsized_op_t sync;        /**< Synchronization operation */
  untyped_op_t team_sync;   /**< Team synchronization operation */
  unsized_op_t sync_all;    /**< Global synchronization operation */
  unsized_op_t barrier;     /**< Barrier operation */

} coll_ops_t;

/** Global collective operations table */
extern coll_ops_t colls;

/**
 * @brief Registration functions for collective operations
 * @param op Name of the operation to register
 * @return 0 on success, non-zero on failure
 */
int register_barrier_all(const char *op);
int register_sync_all(const char *op);
int register_barrier(const char *op);
int register_sync(const char *op);
int register_team_sync(const char *op);

int register_alltoall_type(const char *op);
int register_alltoall_mem(const char *op);
int register_alltoall_size(const char *op);

int register_alltoalls_type(const char *op);
int register_alltoalls_mem(const char *op);
int register_alltoalls_size(const char *op);

int register_collect_type(const char *op);
int register_collect_mem(const char *op);
int register_collect_size(const char *op);

int register_fcollect_type(const char *op);
int register_fcollect_mem(const char *op);
int register_fcollect_size(const char *op);

int register_broadcast_type(const char *op);
int register_broadcast_mem(const char *op);
int register_broadcast_size(const char *op);

int register_and_to_all(const char *op);
int register_or_to_all(const char *op);
int register_xor_to_all(const char *op);
int register_max_to_all(const char *op);
int register_min_to_all(const char *op);
int register_sum_to_all(const char *op);
int register_prod_to_all(const char *op);

int register_and_reduce(const char *op);
int register_or_reduce(const char *op);
int register_xor_reduce(const char *op);
int register_max_reduce(const char *op);
int register_min_reduce(const char *op);
int register_sum_reduce(const char *op);
int register_prod_reduce(const char *op);

int register_inscan(const char *op);
int register_exscan(const char *op);

#endif
