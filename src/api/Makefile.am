# For license: see LICENSE file at top-level

shmemc_dir             = shmemc
shmemc_cppflags        = -I$(top_srcdir)/src/$(shmemc_dir)
shmemc_ldflags         = -L$(shmemc_dir)
shmemc_libs            = -lshmemc-ucx

shmemu_dir             = shmemu
shmemu_cppflags        = -I$(top_srcdir)/src/$(shmemu_dir)
shmemu_ldflags         = -L$(shmemu_dir)
shmemu_libs            = -lshmemu

shmemt_dir             = shmemt
shmemt_cppflags        = -I$(top_srcdir)/src/$(shmemt_dir)
shmemt_ldflags         = -L$(shmemt_dir)
shmemt_libs            = -lshmemt


atomics_dir            = ./atomics
atomics_cppflags       = -I$(top_srcdir)/src/$(shmemc_dir)

MY_INCLUDE             = -I../../include \
				-I$(top_srcdir)/include \
				-I$(top_srcdir)/src/shcoll/src

all_cppflags           = $(atomics_cppflags) \
				$(shmemc_cppflags) \
				$(shmemu_cppflags) \
				$(shmemt_cppflags) \
				$(MY_INCLUDE)

# Core API files (remaining in top-level)
MY_SOURCES             = \
				accessible.c \
				contexts.c \
				fence.c \
				globalexit.c \
				info.c \
				init.c \
				lock.c \
				profiling.c \
				ptr.c \
				quiet.c \
				ranks.c	 \
				shmalloc.c \
				tl_query.c \
				version.c

# Point-to-point synchronization files
MY_SOURCES            += \
				pt2pt_synch/test.c \
				pt2pt_synch/testall.c \
				pt2pt_synch/testany.c \
				pt2pt_synch/testsome.c \
				pt2pt_synch/testall_vector.c \
				pt2pt_synch/testany_vector.c \
				pt2pt_synch/testsome_vector.c \
				pt2pt_synch/waituntil.c \
				pt2pt_synch/waitall.c \
				pt2pt_synch/waitany.c \
				pt2pt_synch/waitsome.c \
				pt2pt_synch/waitall_vector.c \
				pt2pt_synch/waitany_vector.c \
				pt2pt_synch/waitsome_vector.c

# RMA files
MY_SOURCES            += \
				rma/putget.c \
				rma/put_signal.c

# Teams files
MY_SOURCES            += \
				teams/teams.c

# Deprecated files
MY_SOURCES            += \
				deprecated/atomics.c \
				deprecated/cache.c \
				deprecated/init.c \
				deprecated/ranks.c \
				deprecated/shmalloc.c \
				deprecated/waituntil.c \
				deprecated/collectives.c

# Collectives files
MY_SOURCES            += \
				collectives/shcoll-shim.c \
				collectives/table.c       \
				collectives/ucc/ucc.c     \
				collectives/ucc/alltoall.c

# TODO: refactor so that collectives is own subdirectory

SUBDIRS                = 	atomics

#
# -- begin: SHCOLL --
#
# shcoll is bundled with package, but handled differently if pulled in
# from somewhere else.
#

if HAVE_SHCOLL_INTERNAL
all_cppflags          += -I@SHCOLL_DIR@
else
all_cppflags          += -I@SHCOLL_DIR@/include
endif

#
# -- end: SHCOLL --
#

#
# separate out the memory allocation routines
#

liballocator_a_SOURCES    = ../ptmalloc/ptmalloc.c allocator/memalloc.c

if ENABLE_EXPERIMENTAL
liballocator_a_SOURCES   += allocator/xmemalloc.c
endif # experimental

liballocator_a_CPPFLAGS   = -I$(srcdir)/allocator
liballocator_a_CFLAGS     = \
				-DFOOTERS=1 \
				-DCORRUPTION_ERROR_ACTION=report_corruption \
				-DUSAGE_ERROR_ACTION=report_usage_error \
				-DONLY_MSPACES=1 \
				-DHAVE_MORECORE=0 \
				-DHAVE_MMAP=0 \
				-DUSE_LOCKS=1 \
				-DREALLOC_ZERO_BYTES_FREES=1
liballocator_la_SOURCES   = $(liballocator_a_SOURCES)
liballocator_la_CFLAGS    = $(liballocator_a_CFLAGS)
liballocator_la_CPPFLAGS  = $(liballocator_a_CPPFLAGS)
noinst_LIBRARIES          = liballocator.a
noinst_LTLIBRARIES        = liballocator.la

#
#
#

if ENABLE_EXPERIMENTAL

MY_SOURCES            += \
			extensions/fence.c \
			extensions/quiet.c \
			extensions/shmalloc.c \
			extensions/wtime.c \
			extensions/interop.c

all_cppflags          += -I$(srcdir)/extensions

endif # experimental

LIBSHMEM_SOURCES       = $(MY_SOURCES)

lib_LTLIBRARIES        = libshmem.la
libshmem_la_SOURCES    = $(LIBSHMEM_SOURCES)
libshmem_la_CFLAGS     = $(AM_CFLAGS)
libshmem_la_CPPFLAGS   = $(all_cppflags)
libshmem_la_LDFLAGS    = -version-info 0:0:0
libshmem_la_LIBADD     = liballocator.la
