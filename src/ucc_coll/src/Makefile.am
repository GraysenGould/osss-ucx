# For license: see LICENSE file at top-level

UCC_SOURCES = ucc_alltoall.c \
							ucc_coll.c

# SOURCES += util/bithacks.c \
# 				util/broadcast-size.c \
# 				util/rotate.c \
# 				util/scan.c \
# 				util/trees.c \
# 				util/psync_pool.c

FIND_UCC_COLL_H = -I$(top_srcdir)/include \
				-I../../../include

# disable warnings here (should detect how to do this)
# BUILD_CFLAGS            = -w
UCC_BUILD_CFLAGS            =
UCC_BUILD_CFLAGS           += -I.. $(FIND_UCC_COLL_H)
UCC_BUILD_CFLAGS           += -I$(top_srcdir)/src/shmemu \
                         -I$(top_srcdir)/src/shmemc \
                         -I$(top_srcdir)/src/shmemt

lib_LTLIBRARIES         = libucc_coll.la
libucc_coll_la_SOURCES    = $(UCC_SOURCES)
libucc_coll_la_CFLAGS     = $(UCC_BUILD_CFLAGS)

nobase_include_HEADERS  = ucc_coll.h \
									ucc_coll/ucc_alltoall.h



# EXTRA_DIST              = shcoll/compat.h
