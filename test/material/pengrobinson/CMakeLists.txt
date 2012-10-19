add_executable("test_pengrobinson" test_pengrobinson.cc)
target_link_libraries("test_pengrobinson" ${EwomsLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${EwomsLinkDirectories})
include_directories(${EwomsIncludeDirectories})

