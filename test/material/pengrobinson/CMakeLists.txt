add_executable("test_pengrobinson" test_pengrobinson.cc)
target_link_libraries("test_pengrobinson" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

