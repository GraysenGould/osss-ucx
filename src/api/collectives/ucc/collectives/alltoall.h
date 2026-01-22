#ifndef _UCC_ALLTOALL_H
#define _UCC_ALLTOALL_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

/**
 * @brief Macro to declare type-specific alltoall implementation
 *
 * @param _algo Algorithm name
 * @param _type Data type
 * @param _typename Type name string
 */
// #define UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)             \
//   int ucc_##_typename##_alltoall(                                       \
//       shmem_team_t team, _type *dest, const _type *source, size_t nelems);



/**
 * @brief Macro to declare alltoall implementations for all supported types
 *
 * @param _algo Algorithm name to generate declarations for
 */
// #define DECLARE_ALLTOALL_TYPES(_type, _typename)                        \
//   UCC_TYPED_ALLTOALL_DECLARATION(_type, _typename)                   \
// SHMEM_STANDARD_RMA_TYPE_TABLE(DECLARE_ALLTOALL_TYPES)
// #undef DECLARE_ALLTOALL_TYPES

// #define API_ALLTOALL_TYPE(_type, _typename)                                    \
//   int ucc_##_typename##_alltoall(shmem_team_t team, _type *dest,               \
//       const _type *source, size_t nelems);
//
// #define DECL_ALLTOALL(_type, _typename) API_ALLTOALL_TYPE(_type, _typename)
// SHMEM_STANDARD_RMA_TYPE_TABLE(DECL_ALLTOALL)
// #undef DECL_ALLTOALL
// #undef API_ALLTOALL_TYPE

int ucc_float_alltoall(shmem_team_t team, float *dest, const float *source, size_t nelems);
int ucc_double_alltoall(shmem_team_t team, double *dest, const double *source, size_t nelems);
int ucc_longdouble_alltoall(shmem_team_t team, long double *dest, const long double *source, size_t nelems);
int ucc_char_alltoall(shmem_team_t team, char *dest, const char *source, size_t nelems);
int ucc_schar_alltoall(shmem_team_t team, signed char *dest, const signed char *source, size_t nelems);
int ucc_short_alltoall(shmem_team_t team, short *dest, const short *source, size_t nelems);
int ucc_int_alltoall(shmem_team_t team, int *dest, const int *source, size_t nelems);
int ucc_long_alltoall(shmem_team_t team, long *dest, const long *source, size_t nelems);
int ucc_longlong_alltoall(shmem_team_t team, long long *dest, const long long *source, size_t nelems);
int ucc_uchar_alltoall(shmem_team_t team, unsigned char *dest, const unsigned char *source, size_t nelems);
int ucc_ushort_alltoall(shmem_team_t team, unsigned short *dest, const unsigned short *source, size_t nelems);
int ucc_uint_alltoall(shmem_team_t team, unsigned int *dest, const unsigned int *source, size_t nelems);
int ucc_ulong_alltoall(shmem_team_t team, unsigned long *dest, const unsigned long *source, size_t nelems);
int ucc_ulonglong_alltoall(shmem_team_t team, unsigned long long *dest, const unsigned long long *source, size_t nelems);
int ucc_int8_alltoall(shmem_team_t team, int8_t *dest, const int8_t *source, size_t nelems);
int ucc_int16_alltoall(shmem_team_t team, int16_t *dest, const int16_t *source, size_t nelems);
int ucc_int32_alltoall(shmem_team_t team, int32_t *dest, const int32_t *source, size_t nelems);
int ucc_int64_alltoall(shmem_team_t team, int64_t *dest, const int64_t *source, size_t nelems);
int ucc_uint8_alltoall(shmem_team_t team, uint8_t *dest, const uint8_t *source, size_t nelems);
int ucc_uint16_alltoall(shmem_team_t team, uint16_t *dest, const uint16_t *source, size_t nelems);
int ucc_uint32_alltoall(shmem_team_t team, uint32_t *dest, const uint32_t *source, size_t nelems);
int ucc_uint64_alltoall(shmem_team_t team, uint64_t *dest, const uint64_t *source, size_t nelems);
int ucc_size_alltoall(shmem_team_t team, size_t *dest, const size_t *source, size_t nelems);
int ucc_ptrdiff_alltoall(shmem_team_t team, ptrdiff_t *dest, const ptrdiff_t *source, size_t nelems);

int ucc_alltoallmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

#endif /* ! _UCC_ALLTOALL_H */
