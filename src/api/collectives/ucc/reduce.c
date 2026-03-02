#include "ucc.h"
#include <ucc/api/ucc.h>
#include "collectives/reduce.h"
#include <shmem.h>
#include "shmemc.h"
#include "shmemu.h"

#include <string.h>
#include <limits.h> 
#include <assert.h>

#include <stdio.h>
#include <math.h>

#include <shmem/api_types.h>

#include "shmem.h"

#include <limits.h>


#define UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, _op, _ucc_op)    \
inline static void ucc_reduce_helper_##_typename##_##_op(         \
      void *dest, const void *source, size_t nreduce, size_t elem_size, \
      int PE_size, ucc_team_h team_handle) { \
  ucc_status_t status;   \
  ucc_coll_buffer_info_t coll_src_buffer_info =  {  \
    .count = nreduce, \
    .buffer = (void *) source,  \
    .datatype = _ucc_type,  \
    .mem_type = UCC_MEMORY_TYPE_HOST   \
  };  \
  ucc_coll_buffer_info_t coll_dest_buffer_info =  {  \
    .count = nreduce,  \
    .buffer = (void *) dest,  \
    .datatype = _ucc_type,  \
    .mem_type = UCC_MEMORY_TYPE_HOST   \
  };  \
  \
  ucc_coll_args_t coll_args = {  \
    .mask = UCC_COLL_ARGS_FIELD_FLAGS | UCC_COLL_ARGS_FIELD_GLOBAL_WORK_BUFFER,  \
    .coll_type = UCC_COLL_TYPE_ALLREDUCE,  \
    .src.info = coll_src_buffer_info,  \
    .dst.info = coll_dest_buffer_info, \
    .global_work_buffer = shmem_ucc_coll.global_work_buffer,  \
    .flags = UCC_COLL_ARGS_FLAG_MEM_MAPPED_BUFFERS, \
    .op = _ucc_op,  \
  };  \
  ucc_coll_req_h coll_handle;  \
  if (UCC_OK != (status =   \
        ucc_collective_init(&coll_args, &coll_handle, team_handle))){  \
    printf("Could Not Initalize UCC collective. Status: %d\n", status);  \
    return;  \
  }  \
    \
  if (UCC_OK != ucc_collective_post(coll_handle)){  \
    printf("Could No Post UCC collective.\n");  \
    return;  \
  }  \
    \
  /* poll operation until done */  \
  while(ucc_collective_test(coll_handle) == UCC_INPROGRESS) {  \
    /* Drive Collective Progress */  \
    ucc_context_progress(shmem_ucc_coll.context_handle);  \
  }  \
    \
  ucc_collective_finalize(coll_handle);  \
}

/*
 * @brief Macro to define team-based reduction operations
 * @param _typename Type name (e.g. int_sum)
 * @param _type Actual type (e.g. int)
 * @param _op Operation (e.g. sum)
 * @param _algo Algorithm name (e.g. linear)
 *
 *
 * FIXME: branch to check that pwrk is valid should only be done in debug mode
 */
#define UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, _op, _ucc_op)                 \
  int ucc_##_typename##_##_op##_reduce(                             \
      shmem_team_t team, _type *dest, const _type *source, size_t nreduce) {   \
    SHMEMU_CHECK_INIT();                                                       \
    SHMEMU_CHECK_TEAM_VALID(team);                                             \
    SHMEMU_CHECK_SYMMETRIC(dest, "dest");                                      \
    SHMEMU_CHECK_SYMMETRIC(source, "source");                                  \
    shmemc_team_h team_h = (shmemc_team_h)team;                                \
    SHMEMU_CHECK_TEAM_STRIDE(team_h->stride, __func__);                        \
    ucc_reduce_helper_##_typename##_##_op(                                     \
        dest, source, nreduce, sizeof(_type),                                 \
        team_h->nranks, team_h->ucc_team_handle);               \
    return 0;                                                                  \
  }

/* and, xor, and or Operations */
#define DEFINE_REDUCE_BITWISE_TYPES(_type, _typename, _ucc_type)                     \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, and, UCC_OP_BAND)       \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, or, UCC_OP_BOR)         \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, xor, UCC_OP_BXOR)       \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, and, UCC_OP_BAND)              \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, or, UCC_OP_BOR)                \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, xor, UCC_OP_BXOR) 

  UCC_REDUCE_BITWISE_TYPE_TABLE(DEFINE_REDUCE_BITWISE_TYPES)
#undef DEFINE_REDUCE_BITWISE_TYPES

/* Min and Max Operations */
#define DEFINE_REDUCE_MINMAX_TYPES(_type, _typename, _ucc_type)                     \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, min, UCC_OP_MIN)        \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, max, UCC_OP_MAX)        \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, min, UCC_OP_MIN)               \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, max, UCC_OP_MAX)

  UCC_REDUCE_MINMAX_TYPE_TABLE(DEFINE_REDUCE_MINMAX_TYPES)
#undef DEFINE_REDUCE_MINMAX_TYPES


/* Sum and Prod Operations */
#define DEFINE_REDUCE_ARITH_TYPES(_type, _typename, _ucc_type)                     \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, sum, UCC_OP_SUM)       \
  UCC_REDUCE_HELPER_DEFINITION(_type, _typename, _ucc_type, prod, UCC_OP_PROD)     \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, sum, UCC_OP_SUM)              \
  UCC_REDUCE_DEFINITION(_type, _typename, _ucc_type, prod, UCC_OP_PROD)

  UCC_REDUCE_ARITH_TYPE_TABLE(DEFINE_REDUCE_ARITH_TYPES)
#undef DEFINE_REDUCE_ARITH_TYPES


