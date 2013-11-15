# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# defines that must be present in config.h for our headers
set (opm-core_CONFIG_VAR
        HAVE_ERT
        HAVE_SUITESPARSE_UMFPACK_H
        )

# dependencies
set (opm-core_DEPS
        # compile with C99 support if available
        "C99"
        # compile with C++0x/11 support if available
        "CXX11Features REQUIRED"
        # various runtime library enhancements
        "Boost 1.39.0
                COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
        # matrix library
        "BLAS REQUIRED"
        "LAPACK REQUIRED"
        # Tim Davis' SuiteSparse archive
        "SuiteSparse COMPONENTS umfpack"
        # solver
        "SuperLU"
        # xml processing (for config parsing)
        "TinyXML"
        #Parser library
        "opm-parser REQUIRED"
        # Ensembles-based Reservoir Tools (ERT)
        "ERT"
        # DUNE dependency
        "dune-common"
        "dune-istl"
        )
