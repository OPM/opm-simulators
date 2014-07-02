# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# defines that must be present in config.h for our headers
set (ewoms_CONFIG_VAR
	HAVE_QUAD
	HAVE_VALGRIND
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
	"dune-cornerpoint"
	# valgrind client requests
	"Valgrind"
	# quadruple precision floating point calculations
	"Quadmath"
	)
