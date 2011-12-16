# build target for the simple twophase lens problem
ADD_EXECUTABLE("test_tabulation" test_tabulation.cc)
TARGET_LINK_LIBRARIES("test_tabulation" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

