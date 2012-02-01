# build target for the simple twophase lens problem
ADD_EXECUTABLE("test_pengrobinson" test_pengrobinson.cc)
TARGET_LINK_LIBRARIES("test_pengrobinson" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

