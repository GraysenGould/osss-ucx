# For license: see LICENSE file at top-level

SOURCES = alltoall.c \
				alltoalls.c \
				barrier.c \
				broadcast.c \
				collect.c \
				fcollect.c \
				reduce.c \
				scan.c

SOURCES += util/bithacks.c \
				util/broadcast-size.c \
				util/rotate.c \
				util/scan.c \
				util/trees.c \
				util/psync_pool.c

FIND_SHMEM_H = -I$(top_srcdir)/include \
				-I../../../include

# disable warnings here (should detect how to do this)
# BUILD_CFLAGS            = -w
BUILD_CFLAGS            =
BUILD_CFLAGS           += -I.. $(FIND_SHMEM_H)
BUILD_CFLAGS           += -I$(top_srcdir)/src/shmemu \
                         -I$(top_srcdir)/src/shmemc \
                         -I$(top_srcdir)/src/shmemt

lib_LTLIBRARIES         = libshcoll.la
libshcoll_la_SOURCES    = $(SOURCES)
libshcoll_la_CFLAGS     = $(BUILD_CFLAGS)

nobase_include_HEADERS  = shcoll.h \
				shcoll/alltoall.h \
				shcoll/alltoalls.h \
				shcoll/barrier.h \
				shcoll/broadcast.h \
				shcoll/collect.h \
				shcoll/common.h \
				shcoll/fcollect.h \
				shcoll/reduce.h \
				shcoll/scan.h

EXTRA_DIST              = shcoll/compat.h
