# - Turn on warnings when compiling

include (AddOptions)
if (CMAKE_COMPILER_IS_GNUCXX)
  message (STATUS "All warnings enabled: -Wall")
  add_options (ALL_LANGUAGES ALL_BUILDS "-Wall")
endif (CMAKE_COMPILER_IS_GNUCXX)
