# build target for the simple twophase lens problem
ADD_EXECUTABLE("test_ncpflash" test_ncpflash.cc)
TARGET_LINK_LIBRARIES("test_ncpflash" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

