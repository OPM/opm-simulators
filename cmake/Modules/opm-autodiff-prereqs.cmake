# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# defines that must be present in config.h for our headers
set (opm-autodiff_CONFIG_VAR
  HAVE_DUNE_CORNERPOINT
	)

# dependencies
set (opm-autodiff_DEPS
	# Compile with C99 support if available
	"C99"
	# Compile with C++0x/11 support if available
	"CXX11Features"
	# Various runtime library enhancements
	"Boost 1.44.0
		COMPONENTS date_time filesystem system iostreams unit_test_framework REQUIRED"
	# DUNE prerequisites
	"dune-common REQUIRED;
	dune-istl REQUIRED;
        dune-cornerpoint;
	opm-core REQUIRED"
	# Eigen
	"Eigen3 3.1.2 "
	)
