# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# this avoids an annoying deprecation warning on DUNE 2.4 (which we
# are not interested in anyway)
set(DUNE_AVOID_CAPABILITIES_IS_PARALLEL_DEPRECATION_WARNING 1)

if (NOT DEFINED EBOS_USE_ALUGRID)
set(EBOS_USE_ALUGRID 0)
endif()

set(EBOS_USE_ALUGRID "${EBOS_USE_ALUGRID}"
	CACHE STRING "Use dune-alugrid for the ECL black oil simulator if a sufficient dune-alugrid is available")

# defines that must be present in config.h for our headers
set (ewoms_CONFIG_VAR
	HAVE_QUAD
	HAVE_VALGRIND
	HAVE_OPENMP
	HAVE_FINAL
	EBOS_USE_ALUGRID
	DUNE_AVOID_CAPABILITIES_IS_PARALLEL_DEPRECATION_WARNING
	)

# dependencies
set (ewoms_DEPS
	# compile with C++0x/11 support if available
	"CXX11Features REQUIRED"
	# DUNE prerequisites
	"dune-common REQUIRED"
	"dune-localfunctions REQUIRED"
	"dune-geometry REQUIRED"
	"dune-grid REQUIRED"
	"dune-istl REQUIRED"
	"opm-core REQUIRED"
	"opm-material REQUIRED"
	"opm-parser"
	"dune-alugrid"
	"dune-cornerpoint"
	# valgrind client requests
	"Valgrind"
	# quadruple precision floating point calculations
	"Quadmath"
	)
