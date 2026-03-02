#ifndef _UCC_REDUCE_H
#define _UCC_REDUCE_H

//#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

#include <shmem/defs.h>
#include <sys/types.h>
#include <stdint.h>
#include <complex.h>

/********************************************************************************/

/**
 * @brief Macro to declare a single reduction operation for a specific type
 *
 * @param _typename Type name used in function name (e.g. int, float)
 * @param _type Actual C type (e.g. int, float)
 * @param _op Operation name (e.g. sum, prod, min, max)
 * @param _algo Algorithm implementation to use
 */
#define SHCOLL_REDUCE_DECLARE(_typename, _type, _op)                        \
  int ucc_##_typename##_##_op##_reduce(                                     \
      shmem_team_t team, _type *dest, const _type *source, size_t nreduce);

#define DECLARE_REDUCE_BITWISE(_type, _typename)                       \
  SHCOLL_REDUCE_DECLARE(_typename, _type, and)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, or)                          \
  SHCOLL_REDUCE_DECLARE(_typename, _type, xor)                         \
SHMEM_REDUCE_BITWISE_TYPE_TABLE(DECLARE_REDUCE_BITWISE)
#undef DECLARE_REDUCE_BITWISE

#define DECLARE_REDUCE_MINMAX(_type, _typename)                        \
  SHCOLL_REDUCE_DECLARE(_typename, _type, min)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, max)                         \
SHMEM_REDUCE_MINMAX_TYPE_TABLE(DECLARE_REDUCE_MINMAX)
#undef DECLARE_REDUCE_MINMAX

#define DECLARE_REDUCE_ARITH(_type, _typename)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, sum)                         \
  SHCOLL_REDUCE_DECLARE(_typename, _type, prod)                        \
SHMEM_REDUCE_ARITH_TYPE_TABLE(DECLARE_REDUCE_ARITH)
#undef DECLARE_REDUCE_ARITH


#define UCC_REDUCE_BITWISE_TYPE_TABLE(X)                                     \
  X(unsigned char, uchar, UCC_DT_UINT8)                                                      \
  X(unsigned short, ushort, UCC_DT_UINT16)                                                    \
  X(unsigned int, uint, UCC_DT_UINT32)                                                        \
  X(unsigned long, ulong, UCC_DT_UINT64)                                                      \
  X(unsigned long long, ulonglong, UCC_DT_UINT64)                                             \
  X(int8_t, int8, UCC_DT_INT8)                                                              \
  X(int16_t, int16, UCC_DT_INT16)                                                            \
  X(int32_t, int32, UCC_DT_INT32)                                                            \
  X(int64_t, int64, UCC_DT_INT64)                                                            \
  X(uint8_t, uint8, UCC_DT_UINT8)                                                            \
  X(uint16_t, uint16, UCC_DT_UINT16)                                                          \
  X(uint32_t, uint32, UCC_DT_UINT32)                                                          \
  X(uint64_t, uint64, UCC_DT_UINT64)                                                          \
  X(size_t, size, UCC_DT_UINT64)

#define UCC_REDUCE_MINMAX_TYPE_TABLE(X)                                      \
  X(char, char, UCC_DT_INT8)  /*signedness is implementation specifc, IDK what to do */  \
  X(signed char, schar, UCC_DT_INT8)                                                        \
  X(short, short, UCC_DT_INT16)                                                              \
  X(int, int, UCC_DT_INT32)                                                                  \
  X(long, long, UCC_DT_INT64)                                                                \
  X(long long, longlong, UCC_DT_INT64)                                                       \
  X(ptrdiff_t, ptrdiff, UCC_DT_UINT64)                                                        \
  X(unsigned char, uchar, UCC_DT_UINT8)                                                      \
  X(unsigned short, ushort, UCC_DT_UINT16)                                                    \
  X(unsigned int, uint, UCC_DT_UINT32)                                                        \
  X(unsigned long, ulong, UCC_DT_UINT64)                                                      \
  X(unsigned long long, ulonglong, UCC_DT_UINT64)                                             \
  X(int8_t, int8, UCC_DT_INT8)                                                              \
  X(int16_t, int16, UCC_DT_INT16)                                                            \
  X(int32_t, int32, UCC_DT_INT32)                                                            \
  X(int64_t, int64, UCC_DT_INT64)                                                            \
  X(uint8_t, uint8, UCC_DT_UINT8)                                                            \
  X(uint16_t, uint16, UCC_DT_UINT16)                                                          \
  X(uint32_t, uint32, UCC_DT_UINT32)                                                          \
  X(uint64_t, uint64, UCC_DT_UINT64)                                                          \
  X(size_t, size, UCC_DT_UINT64)                                                              \
  X(float, float, UCC_DT_FLOAT32)                                                              \
  X(double, double, UCC_DT_FLOAT64)                                                            \
  X(long double, longdouble, UCC_DT_FLOAT128)

#define UCC_REDUCE_ARITH_TYPE_TABLE(X)                                       \
  X(char, char, UCC_DT_INT8)  /*signedness is implementation specifc, IDK what to do */  \
  X(signed char, schar, UCC_DT_INT8)                                                        \
  X(short, short, UCC_DT_INT16)                                                              \
  X(int, int, UCC_DT_INT32)                                                                  \
  X(long, long, UCC_DT_INT64)                                                                \
  X(long long, longlong, UCC_DT_INT64)                                                       \
  X(ptrdiff_t, ptrdiff, UCC_DT_UINT64)                                                        \
  X(unsigned char, uchar, UCC_DT_UINT8)                                                      \
  X(unsigned short, ushort, UCC_DT_UINT16)                                                    \
  X(unsigned int, uint, UCC_DT_UINT32)                                                        \
  X(unsigned long, ulong, UCC_DT_UINT64)                                                      \
  X(unsigned long long, ulonglong, UCC_DT_UINT64)                                             \
  X(int8_t, int8, UCC_DT_INT8)                                                              \
  X(int16_t, int16, UCC_DT_INT16)                                                            \
  X(int32_t, int32, UCC_DT_INT32)                                                            \
  X(int64_t, int64, UCC_DT_INT64)                                                            \
  X(uint8_t, uint8, UCC_DT_UINT8)                                                            \
  X(uint16_t, uint16, UCC_DT_UINT16)                                                          \
  X(uint32_t, uint32, UCC_DT_UINT32)                                                          \
  X(uint64_t, uint64, UCC_DT_UINT64)                                                          \
  X(size_t, size, UCC_DT_UINT64)                                                              \
  X(float, float, UCC_DT_FLOAT32)                                                              \
  X(double, double, UCC_DT_FLOAT64)                                                            \
  X(long double, longdouble, UCC_DT_FLOAT128)                                                   \
  X(double _Complex, complexd, UCC_DT_FLOAT64_COMPLEX)                                                 \
  X(float _Complex, complexf, UCC_DT_FLOAT32_COMPLEX)


#endif /* ! _UCC_REDUCE_H */
