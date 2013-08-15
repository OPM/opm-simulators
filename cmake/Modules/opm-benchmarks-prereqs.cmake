# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# defines that must be present in config.h for our headers
set (opm-benchmarks_CONFIG_VAR
	)

# dependencies
set (opm-benchmarks_DEPS
	# compile with C++0x/11 support if available
	"CXX11Features REQUIRED"
	# OPM dependency
	"opm-core"
	"opm-upscaling"
	)
