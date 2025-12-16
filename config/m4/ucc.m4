ucc_happy=no

AC_ARG_WITH([ucc],
            [AS_HELP_STRING([--with-ucc@<:@=DIR@:>@], [Use UCC library])],
	    [],
	    [with_ucc=yes])

AS_IF([test "x$with_ucc" = "xyes"], [
        PKG_CHECK_MODULES([UCX], [ucx], [
            #ucc_include_dir="`pkg-config --variable=includedir ucc`"
            # does anything else need to go here?
            ucc_include_dir=/mnt/DISCL/home/gragould/sw/el9-x86_64/ucc/include
            ucc_happy=yes
        ])
        ], [
                AS_IF([test -d "$with_ucc"],[
                    ucc_include_dir=/mnt/DISCL/home/gragould/sw/el9-x86_64/ucc/include
                    ucc_api_dir="$ucc_include_dir/ucc/api"
                    UCC_CFLAGS="-I$ucc_include_dir"
                    UCC_LIBS="-L$with_ucc/lib64 -Wl,-rpath -Wl,$with_ucc/lib64"
                    UCC_LIBS="$UCC_LIBS -L$with_ucc/lib -Wl,-rpath -Wl,$with_ucc/lib"
                    UCC_LIBS="$UCC_LIBS -lucc"
                    ucc_happy=yes
                ])
        ]
)


AS_IF([test "x$ucc_happy" = "xno"], [
            AC_MSG_ERROR([Cannot find required UCC support])
    ], [
            CPPFLAGS="-I$ucc_include_dir $CPPFLAGS"
            LDFLAGS="$UCC_LIBS $LDFLAGS"
            AC_DEFINE([HAVE_UCC], [1], [UCC support])
            AC_MSG_NOTICE([UCC: checking for API features...])
            AC_LANG_PUSH([C])

            #could do some checking of important routines:

            AC_LANG_POP([C])
            AC_SUBST([UCC_LIBS])

            hdr="$ucc_api_dir/ucc_version.h"
            maj=`awk '$2 == "UCC_API_MAJOR" {print $3}' $hdr`
            min=`awk '$2 == "UCC_API_MINOR" {print $3}' $hdr`

            UCC_VERSION_STRING=`printf "%u.%u" $maj $min`
            AS_BOX(UCX version is $UCC_VERSION_STRING)

            AC_DEFINE_UNQUOTED([UCC_VERSION_STRING], ["$UCC_VERSION_STRING"], [Version of UCC])
            AC_SUBST([UCC_VERSION_STRING])
    ]
)

AM_CONDITIONAL([HAVE_UCC], [test "x$ucc_happy" != xno])