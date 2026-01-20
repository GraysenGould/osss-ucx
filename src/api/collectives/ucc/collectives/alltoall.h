#ifndef _UCC_ALLTOALL_H
#define _UCC_ALLTOALL_H

#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

int ucc_alltoallmem(shmem_team_t team, void *dest,
                    const void *source, size_t nelems);

#endif /* ! _UCC_ALLTOALL_H */
