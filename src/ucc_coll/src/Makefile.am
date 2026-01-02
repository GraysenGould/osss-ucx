# For license: see LICENSE file at top-level

SOURCES = ucc_alltoall.c 

# SOURCES += util/bithacks.c \
# 				util/broadcast-size.c \
# 				util/rotate.c \
# 				util/scan.c \
# 				util/trees.c \
# 				util/psync_pool.c

FIND_SHMEM_H = -I$(top_srcdir)/include \
				-I../../../include

# disable warnings here (should detect how to do this)
# BUILD_CFLAGS            = -w
BUILD_CFLAGS            =
BUILD_CFLAGS           += -I.. $(FIND_SHMEM_H)
BUILD_CFLAGS           += -I$(top_srcdir)/src/shmemu \
                         -I$(top_srcdir)/src/shmemc \
                         -I$(top_srcdir)/src/shmemt

lib_LTLIBRARIES         = libucc_coll.la
libucc_coll_la_SOURCES    = $(SOURCES)
libucc_coll_la_CFLAGS     = $(BUILD_CFLAGS)

nobase_include_HEADERS  = ucc_coll.h \
									ucc_coll/ucc_alltoall.h



# EXTRA_DIST              = shcoll/compat.h
