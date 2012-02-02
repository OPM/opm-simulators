# build target for the simple twophase lens problem
add_executable("test_tabulation" test_tabulation.cc)
target_link_libraries("test_tabulation" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

