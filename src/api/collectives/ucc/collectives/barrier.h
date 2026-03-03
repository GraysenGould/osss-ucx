#ifndef _UCC_BARRIER_H
#define _UCC_BARRIER_H


#include <shmem/api_types.h>
#include <shmemu.h>
#include <ucc/api/ucc.h>
#include <stdio.h>
#include <shmem.h>

void ucc_barrier_all();

void ucc_team_sync(shmem_team_t team);

void ucc_sync_all();

#endif /* ! _UCC_BARRIER_H */
