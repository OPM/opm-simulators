# build target for the simple twophase lens problem
add_executable("test_ncpflash" test_ncpflash.cc)
target_link_libraries("test_ncpflash" ${EwomsLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${EwomsLinkDirectories})
include_directories(${EwomsIncludeDirectories})

