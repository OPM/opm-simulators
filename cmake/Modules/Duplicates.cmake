# - Remove duplicate library declarations
#
# Synopsis:
#
#	remove_duplicate_libraries (module)
#
# where
#	module         Name of the module whose libraries should be pruned

# Copyright (C) 2013 Uni Research AS
# This file is licensed under the GNU General Public License v3.0

# libraries should always be trimmed from the beginning, so that also
# missing functions in those later in the list will be resolved
macro (remove_duplicate_libraries module)
  if (DEFINED ${module}_LIBRARIES)
	list (REVERSE ${module}_LIBRARIES)
	list (REMOVE_DUPLICATES ${module}_LIBRARIES)
	list (REVERSE ${module}_LIBRARIES)
  endif (DEFINED ${module}_LIBRARIES)
endmacro (remove_duplicate_libraries module)

# headers can be trimmed from the end, since adding a directory to
# the list is an idempotent action
macro (remove_duplicate_headers module)
  if (DEFINED ${module}_INCLUDE_DIRS)
	list (REMOVE_DUPLICATES ${module}_INCLUDE_DIRS)
  endif (DEFINED ${module}_INCLUDE_DIRS)
endmacro (remove_duplicate_headers module)

# linker flags may not be specified at all
macro (remove_duplicate_flags module)
  if (DEFINED ${module}_LINKER_FLAGS)
	list (REMOVE_DUPLICATES ${module}_LINKER_FLAGS)
  endif (DEFINED ${module}_LINKER_FLAGS)
endmacro (remove_duplicate_flags module)

# fix up both headers and libraries, in case two dependencies have
# included the same second-level library independently
macro (remove_dup_deps module)
  remove_duplicate_headers (${module})
  remove_duplicate_libraries (${module})
  remove_duplicate_flags (${module})
endmacro (remove_dup_deps module)
