# - Try to build faster depending on the compiler

# faster builds by using a pipe instead of temp files
include (AddOptions)
if (CMAKE_COMPILER_IS_GNUCXX)
	add_options (ALL_LANGUAGES ALL_BUILDS "-pipe")
endif (CMAKE_COMPILER_IS_GNUCXX)

