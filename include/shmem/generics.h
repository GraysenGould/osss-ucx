/**
 * @file generics.h
 * @brief OpenSHMEM generic type macros and definitions
 *
 * This file contains C11 generic selection macros that map generic OpenSHMEM
 * operations to their type-specific implementations.
 *
 * The macros use _Generic to select the appropriate typed function based on
 * argument types. For example, shmem_put() will map to shmem_float_put() for
 * float arguments, shmem_long_put() for long arguments, etc.
 *
 * Most operations have both context-based and non-context variants. The context
 * versions take a shmem_ctx_t as their first argument.
 *
 * @copyright For license: see LICENSE file at top-level
 */

#ifndef _SHMEM_API_GENERICS_H
#define _SHMEM_API_GENERICS_H 1

#include <shmem/defs.h>
#include <shmem/depr.h>
#include <shmem/generic_types.h>

#include <sys/types.h>
#include <stddef.h> /* ptrdiff_t */
#include <stdint.h> /* sized int types */
#include <stdarg.h>

#if SHMEM_HAS_C11

/**
 * @defgroup generics OpenSHMEM Generic Operations
 * @{
 */

/* shmem_sync is a special case */
#define shmem_sync(team) shmem_team_sync(team)

/**
 * @brief Helper function for generic selection fallback
 *
 * This stops the not-a-context case turning into an error when
 * the value type doesn't match anything
 */
inline static void shmem_generics_nomatch(void *dummy, ...) {}

/**
 * @brief Get numbered args out of parameter list
 * (thanks to SOS)
 */
#define SHC11_GET_ARG1_HELPER(_arg1, ...) _arg1

#define SHC11_GET_ARG1(...) SHC11_GET_ARG1_HELPER(__VA_ARGS__, _extra)
#define SHC11_GET_ARG2(_arg1, ...) SHC11_GET_ARG1_HELPER(__VA_ARGS__, _extra)

#define SHC11_TYPE_EVAL_PTR(_arg) &*(_arg)
#define SHC11_TYPE_EVAL_PTR_OR_SCALAR(_arg) (_arg) + 0

/********************************************************************************
 *
 *                                RMA routines
 *
 *********************************************************************************
 */

/* Generalized X-macro case generators for RMA routines */
#define SHMEM_RMA_GENERIC_CASE(op, type, typename)                             \
  type * : shmem_##typename##_##op,
#define SHMEM_CTX_RMA_GENERIC_CASE(op, type, typename)                         \
  type * : shmem_ctx_##typename##_##op,

/* Per-op wrapper macros for use with X-macro tables */
#define SHMEM_RMA_GENERIC_CASE_put(type, typename)                             \
  SHMEM_RMA_GENERIC_CASE(put, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_put(type, typename)                         \
  SHMEM_CTX_RMA_GENERIC_CASE(put, type, typename)
#define SHMEM_RMA_GENERIC_CASE_p(type, typename)                               \
  SHMEM_RMA_GENERIC_CASE(p, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_p(type, typename)                           \
  SHMEM_CTX_RMA_GENERIC_CASE(p, type, typename)
#define SHMEM_RMA_GENERIC_CASE_get(type, typename)                             \
  SHMEM_RMA_GENERIC_CASE(get, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_get(type, typename)                         \
  SHMEM_CTX_RMA_GENERIC_CASE(get, type, typename)
#define SHMEM_RMA_GENERIC_CASE_g(type, typename)                               \
  SHMEM_RMA_GENERIC_CASE(g, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_g(type, typename)                           \
  SHMEM_CTX_RMA_GENERIC_CASE(g, type, typename)
#define SHMEM_RMA_GENERIC_CASE_iput(type, typename)                            \
  SHMEM_RMA_GENERIC_CASE(iput, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_iput(type, typename)                        \
  SHMEM_CTX_RMA_GENERIC_CASE(iput, type, typename)
#define SHMEM_RMA_GENERIC_CASE_iget(type, typename)                            \
  SHMEM_RMA_GENERIC_CASE(iget, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_iget(type, typename)                        \
  SHMEM_CTX_RMA_GENERIC_CASE(iget, type, typename)
#define SHMEM_RMA_GENERIC_CASE_put_nbi(type, typename)                         \
  SHMEM_RMA_GENERIC_CASE(put_nbi, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_put_nbi(type, typename)                     \
  SHMEM_CTX_RMA_GENERIC_CASE(put_nbi, type, typename)
#define SHMEM_RMA_GENERIC_CASE_get_nbi(type, typename)                         \
  SHMEM_RMA_GENERIC_CASE(get_nbi, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_get_nbi(type, typename)                     \
  SHMEM_CTX_RMA_GENERIC_CASE(get_nbi, type, typename)
#define SHMEM_RMA_GENERIC_CASE_put_signal(type, typename)                      \
  SHMEM_RMA_GENERIC_CASE(put_signal, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_put_signal(type, typename)                  \
  SHMEM_CTX_RMA_GENERIC_CASE(put_signal, type, typename)
#define SHMEM_RMA_GENERIC_CASE_put_signal_nbi(type, typename)                  \
  SHMEM_RMA_GENERIC_CASE(put_signal_nbi, type, typename)
#define SHMEM_CTX_RMA_GENERIC_CASE_put_signal_nbi(type, typename)              \
  SHMEM_CTX_RMA_GENERIC_CASE(put_signal_nbi, type, typename)

#define shmem_put(...)                                                          \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                    \
              SHMEM_CTX_RMA_GENERIC_CASE_put) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                        \
          SHMEM_RMA_GENERIC_CASE_put) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_p(...)                                                           \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      shmem_ctx_t: _Generic(                                                   \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),          \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                   \
              SHMEM_CTX_RMA_GENERIC_CASE_p) default: shmem_generics_nomatch),  \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                       \
          SHMEM_RMA_GENERIC_CASE_p) default: shmem_generics_nomatch)(          \
      __VA_ARGS__)

#define shmem_get(...)                                                          \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                    \
              SHMEM_CTX_RMA_GENERIC_CASE_get) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                        \
          SHMEM_RMA_GENERIC_CASE_get) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_g(...)                                                           \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      shmem_ctx_t: _Generic(                                                   \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),          \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                   \
              SHMEM_CTX_RMA_GENERIC_CASE_g) default: shmem_generics_nomatch),  \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                       \
          SHMEM_RMA_GENERIC_CASE_g) default: shmem_generics_nomatch)(          \
      __VA_ARGS__)

#define shmem_iput(...)                                                          \
  _Generic(                                                                      \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                          \
      shmem_ctx_t: _Generic(                                                     \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),            \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                     \
              SHMEM_CTX_RMA_GENERIC_CASE_iput) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                         \
          SHMEM_RMA_GENERIC_CASE_iput) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_iget(...)                                                          \
  _Generic(                                                                      \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                          \
      shmem_ctx_t: _Generic(                                                     \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),            \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                     \
              SHMEM_CTX_RMA_GENERIC_CASE_iget) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                         \
          SHMEM_RMA_GENERIC_CASE_iget) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_put_nbi(...)                                                          \
  _Generic(                                                                         \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                             \
      shmem_ctx_t: _Generic(                                                        \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),               \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                        \
              SHMEM_CTX_RMA_GENERIC_CASE_put_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                            \
          SHMEM_RMA_GENERIC_CASE_put_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_get_nbi(...)                                                          \
  _Generic(                                                                         \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                             \
      shmem_ctx_t: _Generic(                                                        \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),               \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                        \
              SHMEM_CTX_RMA_GENERIC_CASE_get_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                            \
          SHMEM_RMA_GENERIC_CASE_get_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_put_signal(...)                                                          \
  _Generic(                                                                            \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                \
      shmem_ctx_t: _Generic(                                                           \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                  \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                           \
              SHMEM_CTX_RMA_GENERIC_CASE_put_signal) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                               \
          SHMEM_RMA_GENERIC_CASE_put_signal) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_put_signal_nbi(...)                                                          \
  _Generic(                                                                                \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                    \
      shmem_ctx_t: _Generic(                                                               \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                      \
          C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                               \
              SHMEM_CTX_RMA_GENERIC_CASE_put_signal_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                                   \
          SHMEM_RMA_GENERIC_CASE_put_signal_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

/********************************************************************************
 *
 *                                Collective routines
 *
 *********************************************************************************
 */

/* Generalized X-macro case generator for collectives */
#define SHMEM_COLLECTIVE_GENERIC_CASE(op, type, typename)                      \
  type * : shmem_##typename##_##op,

/* Per-op wrapper macros for use with X-macro tables */
#define SHMEM_COLLECTIVE_GENERIC_CASE_collect(type, typename)                  \
  SHMEM_COLLECTIVE_GENERIC_CASE(collect, type, typename)
#define SHMEM_COLLECTIVE_GENERIC_CASE_fcollect(type, typename)                 \
  SHMEM_COLLECTIVE_GENERIC_CASE(fcollect, type, typename)
#define SHMEM_COLLECTIVE_GENERIC_CASE_alltoall(type, typename)                 \
  SHMEM_COLLECTIVE_GENERIC_CASE(alltoall, type, typename)
#define SHMEM_COLLECTIVE_GENERIC_CASE_alltoalls(type, typename)                \
  SHMEM_COLLECTIVE_GENERIC_CASE(alltoalls, type, typename)
#define SHMEM_COLLECTIVE_GENERIC_CASE_broadcast(type, typename)                \
  SHMEM_COLLECTIVE_GENERIC_CASE(broadcast, type, typename)

#define shmem_collect(...)                                                         \
  _Generic(                                                                        \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                            \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                           \
          SHMEM_COLLECTIVE_GENERIC_CASE_collect) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_fcollect(...)                                                         \
  _Generic(                                                                         \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                             \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                            \
          SHMEM_COLLECTIVE_GENERIC_CASE_fcollect) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_alltoall(...)                                                         \
  _Generic(                                                                         \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                             \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                            \
          SHMEM_COLLECTIVE_GENERIC_CASE_alltoall) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_alltoalls(...)                                                         \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                              \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                             \
          SHMEM_COLLECTIVE_GENERIC_CASE_alltoalls) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_broadcast(...)                                                         \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                              \
      C11_SHMEM_STANDARD_RMA_TYPE_TABLE(                                             \
          SHMEM_COLLECTIVE_GENERIC_CASE_broadcast) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

/********************************************************************************
 *
 *                              Reduction routines
 *
 *********************************************************************************
 */

/* --- to_all reductions --- */
#define SHMEM_TO_ALL_GENERIC_CASE(op, type, typename)                          \
  type * : shmem_##typename##_##op##_to_all,
#define SHMEM_TO_ALL_GENERIC_CASE_and(type, typename)                          \
  SHMEM_TO_ALL_GENERIC_CASE(and, type, typename)
#define SHMEM_TO_ALL_GENERIC_CASE_or(type, typename)                           \
  SHMEM_TO_ALL_GENERIC_CASE(or, type, typename)
#define SHMEM_TO_ALL_GENERIC_CASE_xor(type, typename)                          \
  SHMEM_TO_ALL_GENERIC_CASE(xor, type, typename)
#define SHMEM_TO_ALL_GENERIC_CASE_max(type, typename)                          \
  SHMEM_TO_ALL_GENERIC_CASE(max, type, typename)
#define SHMEM_TO_ALL_GENERIC_CASE_min(type, typename)                          \
  SHMEM_TO_ALL_GENERIC_CASE(min, type, typename)
#define SHMEM_TO_ALL_GENERIC_CASE_sum(type, typename)                          \
  SHMEM_TO_ALL_GENERIC_CASE(sum, type, typename)
#define SHMEM_TO_ALL_GENERIC_CASE_prod(type, typename)                         \
  SHMEM_TO_ALL_GENERIC_CASE(prod, type, typename)

#define shmem_and_to_all(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_TO_ALL_BITWISE_TYPE_TABLE(                                     \
          SHMEM_TO_ALL_GENERIC_CASE_and) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_or_to_all(...)                                                   \
  _Generic(SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                   \
      C11_SHMEM_TO_ALL_BITWISE_TYPE_TABLE(                                     \
               SHMEM_TO_ALL_GENERIC_CASE_or) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_xor_to_all(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_TO_ALL_BITWISE_TYPE_TABLE(                                     \
          SHMEM_TO_ALL_GENERIC_CASE_xor) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_max_to_all(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_TO_ALL_MINMAX_TYPE_TABLE(                                      \
          SHMEM_TO_ALL_GENERIC_CASE_max) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_min_to_all(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_TO_ALL_MINMAX_TYPE_TABLE(                                      \
          SHMEM_TO_ALL_GENERIC_CASE_min) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_sum_to_all(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_TO_ALL_ARITH_TYPE_TABLE(                                       \
          SHMEM_TO_ALL_GENERIC_CASE_sum) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_prod_to_all(...)                                                 \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_TO_ALL_ARITH_TYPE_TABLE(                                       \
          SHMEM_TO_ALL_GENERIC_CASE_prod) default: shmem_generics_nomatch)(    \
      __VA_ARGS__)

/* --- team-based reduce reductions --- */
#define SHMEM_REDUCE_GENERIC_CASE(op, type, typename)                          \
  type * : shmem_##typename##_##op##_reduce,
#define SHMEM_REDUCE_GENERIC_CASE_and(type, typename)                          \
  SHMEM_REDUCE_GENERIC_CASE(and, type, typename)
#define SHMEM_REDUCE_GENERIC_CASE_or(type, typename)                           \
  SHMEM_REDUCE_GENERIC_CASE(or, type, typename)
#define SHMEM_REDUCE_GENERIC_CASE_xor(type, typename)                          \
  SHMEM_REDUCE_GENERIC_CASE(xor, type, typename)
#define SHMEM_REDUCE_GENERIC_CASE_max(type, typename)                          \
  SHMEM_REDUCE_GENERIC_CASE(max, type, typename)
#define SHMEM_REDUCE_GENERIC_CASE_min(type, typename)                          \
  SHMEM_REDUCE_GENERIC_CASE(min, type, typename)
#define SHMEM_REDUCE_GENERIC_CASE_sum(type, typename)                          \
  SHMEM_REDUCE_GENERIC_CASE(sum, type, typename)
#define SHMEM_REDUCE_GENERIC_CASE_prod(type, typename)                         \
  SHMEM_REDUCE_GENERIC_CASE(prod, type, typename)

#define shmem_and_reduce(...)                                                                               \
  _Generic(                                                                                                 \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                                     \
      C11_SHMEM_REDUCE_BITWISE_CANONICAL_TYPES(SHMEM_REDUCE_GENERIC_CASE_and) default: _Generic(            \
          SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                                 \
          C11_SHMEM_REDUCE_BITWISE_STDINT_SMALL_TYPES(SHMEM_REDUCE_GENERIC_CASE_and) default: _Generic(     \
              SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                             \
              C11_SHMEM_REDUCE_BITWISE_STDINT_64BIT_TYPES(SHMEM_REDUCE_GENERIC_CASE_and) default: _Generic( \
                  SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                         \
                  C11_SHMEM_REDUCE_BITWISE_STDINT_SIZE_TYPES(                                               \
                      SHMEM_REDUCE_GENERIC_CASE_and) default: shmem_generics_nomatch))))(                   \
      __VA_ARGS__)

#define shmem_or_reduce(...)                                                                               \
  _Generic(                                                                                                \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                                    \
      C11_SHMEM_REDUCE_BITWISE_CANONICAL_TYPES(SHMEM_REDUCE_GENERIC_CASE_or) default: _Generic(            \
          SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                                \
          C11_SHMEM_REDUCE_BITWISE_STDINT_SMALL_TYPES(SHMEM_REDUCE_GENERIC_CASE_or) default: _Generic(     \
              SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                            \
              C11_SHMEM_REDUCE_BITWISE_STDINT_64BIT_TYPES(SHMEM_REDUCE_GENERIC_CASE_or) default: _Generic( \
                  SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                        \
                  C11_SHMEM_REDUCE_BITWISE_STDINT_SIZE_TYPES(                                              \
                      SHMEM_REDUCE_GENERIC_CASE_or) default: shmem_generics_nomatch))))(                   \
      __VA_ARGS__)

#define shmem_xor_reduce(...)                                                                               \
  _Generic(                                                                                                 \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                                     \
      C11_SHMEM_REDUCE_BITWISE_CANONICAL_TYPES(SHMEM_REDUCE_GENERIC_CASE_xor) default: _Generic(            \
          SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                                 \
          C11_SHMEM_REDUCE_BITWISE_STDINT_SMALL_TYPES(SHMEM_REDUCE_GENERIC_CASE_xor) default: _Generic(     \
              SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                             \
              C11_SHMEM_REDUCE_BITWISE_STDINT_64BIT_TYPES(SHMEM_REDUCE_GENERIC_CASE_xor) default: _Generic( \
                  SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                                         \
                  C11_SHMEM_REDUCE_BITWISE_STDINT_SIZE_TYPES(                                               \
                      SHMEM_REDUCE_GENERIC_CASE_xor) default: shmem_generics_nomatch))))(                   \
      __VA_ARGS__)

#define shmem_max_reduce(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                        \
      C11_SHMEM_REDUCE_MINMAX_TYPE_TABLE(                                      \
          SHMEM_REDUCE_GENERIC_CASE_max) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_min_reduce(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                        \
      C11_SHMEM_REDUCE_MINMAX_TYPE_TABLE(                                      \
          SHMEM_REDUCE_GENERIC_CASE_min) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_sum_reduce(...)                                                  \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                        \
      C11_SHMEM_REDUCE_ARITH_TYPE_TABLE(                                       \
          SHMEM_REDUCE_GENERIC_CASE_sum) default: shmem_generics_nomatch)(     \
      __VA_ARGS__)

#define shmem_prod_reduce(...)                                                 \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG2(__VA_ARGS__)),                        \
      C11_SHMEM_REDUCE_ARITH_TYPE_TABLE(                                       \
          SHMEM_REDUCE_GENERIC_CASE_prod) default: shmem_generics_nomatch)(    \
      __VA_ARGS__)

/********************************************************************************
 *
 *                                AMO routines
 *
 *********************************************************************************
 */

/* Generalized X-macro case generator for AMO routines */
#define SHMEM_AMO_GENERIC_CASE(op, type, typename)                             \
  type * : shmem_##typename##_atomic_##op,
#define SHMEM_CTX_AMO_GENERIC_CASE(op, type, typename)                         \
  type * : shmem_ctx_##typename##_atomic_##op,

/* Per-op wrapper macros for use with X-macro tables */
#define SHMEM_AMO_GENERIC_CASE_add(type, typename)                             \
  SHMEM_AMO_GENERIC_CASE(add, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_add(type, typename)                         \
  SHMEM_CTX_AMO_GENERIC_CASE(add, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_add(type, typename)                       \
  SHMEM_AMO_GENERIC_CASE(fetch_add, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_add(type, typename)                   \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_add, type, typename)
#define SHMEM_AMO_GENERIC_CASE_inc(type, typename)                             \
  SHMEM_AMO_GENERIC_CASE(inc, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_inc(type, typename)                         \
  SHMEM_CTX_AMO_GENERIC_CASE(inc, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_inc(type, typename)                       \
  SHMEM_AMO_GENERIC_CASE(fetch_inc, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_inc(type, typename)                   \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_inc, type, typename)
#define SHMEM_AMO_GENERIC_CASE_compare_swap(type, typename)                    \
  SHMEM_AMO_GENERIC_CASE(compare_swap, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_compare_swap(type, typename)                \
  SHMEM_CTX_AMO_GENERIC_CASE(compare_swap, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch(type, typename)                           \
  SHMEM_AMO_GENERIC_CASE(fetch, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch(type, typename)                       \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch, type, typename)
#define SHMEM_AMO_GENERIC_CASE_set(type, typename)                             \
  SHMEM_AMO_GENERIC_CASE(set, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_set(type, typename)                         \
  SHMEM_CTX_AMO_GENERIC_CASE(set, type, typename)
#define SHMEM_AMO_GENERIC_CASE_swap(type, typename)                            \
  SHMEM_AMO_GENERIC_CASE(swap, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_swap(type, typename)                        \
  SHMEM_CTX_AMO_GENERIC_CASE(swap, type, typename)
#define SHMEM_AMO_GENERIC_CASE_and(type, typename)                             \
  SHMEM_AMO_GENERIC_CASE(and, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_and(type, typename)                         \
  SHMEM_CTX_AMO_GENERIC_CASE(and, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_and(type, typename)                       \
  SHMEM_AMO_GENERIC_CASE(fetch_and, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_and(type, typename)                   \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_and, type, typename)
#define SHMEM_AMO_GENERIC_CASE_or(type, typename)                              \
  SHMEM_AMO_GENERIC_CASE(or, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_or(type, typename)                          \
  SHMEM_CTX_AMO_GENERIC_CASE(or, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_or(type, typename)                        \
  SHMEM_AMO_GENERIC_CASE(fetch_or, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_or(type, typename)                    \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_or, type, typename)
#define SHMEM_AMO_GENERIC_CASE_xor(type, typename)                             \
  SHMEM_AMO_GENERIC_CASE(xor, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_xor(type, typename)                         \
  SHMEM_CTX_AMO_GENERIC_CASE(xor, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_xor(type, typename)                       \
  SHMEM_AMO_GENERIC_CASE(fetch_xor, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_xor(type, typename)                   \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_xor, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_inc_nbi(type, typename)                   \
  SHMEM_AMO_GENERIC_CASE(fetch_inc_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_inc_nbi(type, typename)               \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_inc_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_add_nbi(type, typename)                   \
  SHMEM_AMO_GENERIC_CASE(fetch_add_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_add_nbi(type, typename)               \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_add_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_compare_swap_nbi(type, typename)                \
  SHMEM_AMO_GENERIC_CASE(compare_swap_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_compare_swap_nbi(type, typename)            \
  SHMEM_CTX_AMO_GENERIC_CASE(compare_swap_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_and_nbi(type, typename)                   \
  SHMEM_AMO_GENERIC_CASE(fetch_and_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_and_nbi(type, typename)               \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_and_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_or_nbi(type, typename)                    \
  SHMEM_AMO_GENERIC_CASE(fetch_or_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_or_nbi(type, typename)                \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_or_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_nbi(type, typename)                       \
  SHMEM_AMO_GENERIC_CASE(fetch_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_nbi(type, typename)                   \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_swap_nbi(type, typename)                        \
  SHMEM_AMO_GENERIC_CASE(swap_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_swap_nbi(type, typename)                    \
  SHMEM_CTX_AMO_GENERIC_CASE(swap_nbi, type, typename)
#define SHMEM_AMO_GENERIC_CASE_fetch_xor_nbi(type, typename)                   \
  SHMEM_AMO_GENERIC_CASE(fetch_xor_nbi, type, typename)
#define SHMEM_CTX_AMO_GENERIC_CASE_fetch_xor_nbi(type, typename)               \
  SHMEM_CTX_AMO_GENERIC_CASE(fetch_xor_nbi, type, typename)

/* Standard AMO operations */
#define shmem_atomic_add(...)                                                   \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                    \
              SHMEM_CTX_AMO_GENERIC_CASE_add) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                        \
          SHMEM_AMO_GENERIC_CASE_add) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_add(...)                                                   \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      shmem_ctx_t: _Generic(                                                          \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                 \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                          \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_add) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                              \
          SHMEM_AMO_GENERIC_CASE_fetch_add) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_inc(...)                                                   \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                    \
              SHMEM_CTX_AMO_GENERIC_CASE_inc) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                        \
          SHMEM_AMO_GENERIC_CASE_inc) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_inc(...)                                                   \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      shmem_ctx_t: _Generic(                                                          \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                 \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                          \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_inc) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                              \
          SHMEM_AMO_GENERIC_CASE_fetch_inc) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_compare_swap(...)                                                   \
  _Generic(                                                                              \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                  \
      shmem_ctx_t: _Generic(                                                             \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                    \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                             \
              SHMEM_CTX_AMO_GENERIC_CASE_compare_swap) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                 \
          SHMEM_AMO_GENERIC_CASE_compare_swap) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

/* Extended AMO operations */
#define shmem_atomic_fetch(...)                                                   \
  _Generic(                                                                       \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                           \
      shmem_ctx_t: _Generic(                                                      \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),             \
          C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                      \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch) default: shmem_generics_nomatch), \
      C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                          \
          SHMEM_AMO_GENERIC_CASE_fetch) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_set(...)                                                   \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                    \
              SHMEM_CTX_AMO_GENERIC_CASE_set) default: shmem_generics_nomatch), \
      C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                        \
          SHMEM_AMO_GENERIC_CASE_set) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_swap(...)                                                   \
  _Generic(                                                                      \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                          \
      shmem_ctx_t: _Generic(                                                     \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),            \
          C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                     \
              SHMEM_CTX_AMO_GENERIC_CASE_swap) default: shmem_generics_nomatch), \
      C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                         \
          SHMEM_AMO_GENERIC_CASE_swap) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

/* Bitwise AMO operations */
#define shmem_atomic_and(...)                                                   \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                     \
              SHMEM_CTX_AMO_GENERIC_CASE_and) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                         \
          SHMEM_AMO_GENERIC_CASE_and) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_and(...)                                                   \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      shmem_ctx_t: _Generic(                                                          \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                 \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                           \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_and) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                               \
          SHMEM_AMO_GENERIC_CASE_fetch_and) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_or(...)                                                   \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      shmem_ctx_t: _Generic(                                                   \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),          \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                    \
              SHMEM_CTX_AMO_GENERIC_CASE_or) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                        \
          SHMEM_AMO_GENERIC_CASE_or) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_or(...)                                                   \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                              \
      shmem_ctx_t: _Generic(                                                         \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                          \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_or) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                              \
          SHMEM_AMO_GENERIC_CASE_fetch_or) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_xor(...)                                                   \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      shmem_ctx_t: _Generic(                                                    \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),           \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                     \
              SHMEM_CTX_AMO_GENERIC_CASE_xor) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                         \
          SHMEM_AMO_GENERIC_CASE_xor) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_xor(...)                                                   \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      shmem_ctx_t: _Generic(                                                          \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                 \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                           \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_xor) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                               \
          SHMEM_AMO_GENERIC_CASE_fetch_xor) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_inc_nbi(...)                                                   \
  _Generic(                                                                               \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                   \
      shmem_ctx_t: _Generic(                                                              \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                     \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                              \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_inc_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                  \
          SHMEM_AMO_GENERIC_CASE_fetch_inc_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_add_nbi(...)                                                   \
  _Generic(                                                                               \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                   \
      shmem_ctx_t: _Generic(                                                              \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                     \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                              \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_add_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                  \
          SHMEM_AMO_GENERIC_CASE_fetch_add_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_compare_swap_nbi(...)                                                   \
  _Generic(                                                                                  \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                      \
      shmem_ctx_t: _Generic(                                                                 \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                        \
          C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                 \
              SHMEM_CTX_AMO_GENERIC_CASE_compare_swap_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                     \
          SHMEM_AMO_GENERIC_CASE_compare_swap_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_and_nbi(...)                                                   \
  _Generic(                                                                               \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                   \
      shmem_ctx_t: _Generic(                                                              \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                     \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                               \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_and_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                                   \
          SHMEM_AMO_GENERIC_CASE_fetch_and_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_or_nbi(...)                                                   \
  _Generic(                                                                              \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                  \
      shmem_ctx_t: _Generic(                                                             \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                    \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                              \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_or_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                                  \
          SHMEM_AMO_GENERIC_CASE_fetch_or_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_nbi(...)                                                   \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      shmem_ctx_t: _Generic(                                                          \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                 \
          C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                          \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                              \
          SHMEM_AMO_GENERIC_CASE_fetch_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_swap_nbi(...)                                                   \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                              \
      shmem_ctx_t: _Generic(                                                         \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                \
          C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                         \
              SHMEM_CTX_AMO_GENERIC_CASE_swap_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_EXTENDED_AMO_TYPE_TABLE(                                             \
          SHMEM_AMO_GENERIC_CASE_swap_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

#define shmem_atomic_fetch_xor_nbi(...)                                                   \
  _Generic(                                                                               \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                   \
      shmem_ctx_t: _Generic(                                                              \
          SHC11_TYPE_EVAL_PTR_OR_SCALAR(SHC11_GET_ARG2(__VA_ARGS__)),                     \
          C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                               \
              SHMEM_CTX_AMO_GENERIC_CASE_fetch_xor_nbi) default: shmem_generics_nomatch), \
      C11_SHMEM_BITWISE_AMO_TYPE_TABLE(                                                   \
          SHMEM_AMO_GENERIC_CASE_fetch_xor_nbi) default: shmem_generics_nomatch)(         \
      __VA_ARGS__)

/********************************************************************************
 *
 *                             PT2PT Sync routines
 *
 *********************************************************************************
 */

/* Type table for wait_until and test (includes short and ushort) */
#define SHMEM_SYNC_WAIT_TEST_TYPE_TABLE(X)                                     \
  X(short, short)                                                              \
  X(unsigned short, ushort)                                                    \
  C11_SHMEM_STANDARD_AMO_TYPE_TABLE(X)

/* Generalized X-macro case generator for pt2pt sync routines */
#define SHMEM_SYNC_GENERIC_CASE(op, type, typename)                            \
  type * : shmem_##typename##_##op,
#define SHMEM_SYNC_GENERIC_CASE_vector(op, type, typename)                     \
  type * : shmem_##typename##_##op##_vector,

/* Per-op wrapper macros */
#define SHMEM_SYNC_GENERIC_CASE_wait_until(type, typename)                     \
  SHMEM_SYNC_GENERIC_CASE(wait_until, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_wait_until_all(type, typename)                 \
  SHMEM_SYNC_GENERIC_CASE(wait_until_all, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_wait_until_any(type, typename)                 \
  SHMEM_SYNC_GENERIC_CASE(wait_until_any, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_wait_until_some(type, typename)                \
  SHMEM_SYNC_GENERIC_CASE(wait_until_some, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_wait_until_all_vector(type, typename)          \
  SHMEM_SYNC_GENERIC_CASE_vector(wait_until_all, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_wait_until_any_vector(type, typename)          \
  SHMEM_SYNC_GENERIC_CASE_vector(wait_until_any, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_wait_until_some_vector(type, typename)         \
  SHMEM_SYNC_GENERIC_CASE_vector(wait_until_some, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test(type, typename)                           \
  SHMEM_SYNC_GENERIC_CASE(test, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test_all(type, typename)                       \
  SHMEM_SYNC_GENERIC_CASE(test_all, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test_any(type, typename)                       \
  SHMEM_SYNC_GENERIC_CASE(test_any, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test_some(type, typename)                      \
  SHMEM_SYNC_GENERIC_CASE(test_some, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test_all_vector(type, typename)                \
  SHMEM_SYNC_GENERIC_CASE_vector(test_all, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test_any_vector(type, typename)                \
  SHMEM_SYNC_GENERIC_CASE_vector(test_any, type, typename)
#define SHMEM_SYNC_GENERIC_CASE_test_some_vector(type, typename)               \
  SHMEM_SYNC_GENERIC_CASE_vector(test_some, type, typename)

/* Macros for pt2pt sync routines */
#define shmem_wait_until(...)                                                   \
  _Generic(                                                                     \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                         \
      SHMEM_SYNC_WAIT_TEST_TYPE_TABLE(                                          \
          SHMEM_SYNC_GENERIC_CASE_wait_until) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_wait_until_all(...)                                                   \
  _Generic(                                                                         \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                             \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                            \
          SHMEM_SYNC_GENERIC_CASE_wait_until_all) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_wait_until_any(...)                                                   \
  _Generic(                                                                         \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                             \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                            \
          SHMEM_SYNC_GENERIC_CASE_wait_until_any) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_wait_until_some(...)                                                   \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                              \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                             \
          SHMEM_SYNC_GENERIC_CASE_wait_until_some) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_wait_until_all_vector(...)                                                   \
  _Generic(                                                                                \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                    \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                   \
          SHMEM_SYNC_GENERIC_CASE_wait_until_all_vector) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_wait_until_any_vector(...)                                                   \
  _Generic(                                                                                \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                    \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                   \
          SHMEM_SYNC_GENERIC_CASE_wait_until_any_vector) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_wait_until_some_vector(...)                                                   \
  _Generic(                                                                                 \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                                     \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                                    \
          SHMEM_SYNC_GENERIC_CASE_wait_until_some_vector) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_test(...)                                                        \
  _Generic(SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                   \
      SHMEM_SYNC_WAIT_TEST_TYPE_TABLE(                                         \
               SHMEM_SYNC_GENERIC_CASE_test) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_test_all(...)                                                    \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                       \
          SHMEM_SYNC_GENERIC_CASE_test_all) default: shmem_generics_nomatch)(  \
      __VA_ARGS__)

#define shmem_test_any(...)                                                    \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                       \
          SHMEM_SYNC_GENERIC_CASE_test_any) default: shmem_generics_nomatch)(  \
      __VA_ARGS__)

#define shmem_test_some(...)                                                   \
  _Generic(                                                                    \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                        \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                       \
          SHMEM_SYNC_GENERIC_CASE_test_some) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_test_all_vector(...)                                                   \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                              \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                             \
          SHMEM_SYNC_GENERIC_CASE_test_all_vector) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_test_any_vector(...)                                                   \
  _Generic(                                                                          \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                              \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                             \
          SHMEM_SYNC_GENERIC_CASE_test_any_vector) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_test_some_vector(...)                                                   \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      C11_SHMEM_STANDARD_AMO_TYPE_TABLE(                                              \
          SHMEM_SYNC_GENERIC_CASE_test_some_vector) default: shmem_generics_nomatch)( \
      __VA_ARGS__)

#define shmem_signal_set(...)                                                         \
  _Generic(                                                                           \
      SHC11_TYPE_EVAL_PTR(SHC11_GET_ARG1(__VA_ARGS__)),                               \
      shmem_ctx_t: shmem_ctx_signal_set,                                              \
      default: shmem_signal_set)(__VA_ARGS__)

#endif /* SHMEM_HAS_C11 */

#endif /* ! _SHMEM_API_GENERICS_H */
