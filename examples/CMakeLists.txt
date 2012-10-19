# add build targets
ADD_EXECUTABLE("tutorial_decoupled" tutorial_decoupled.cc)
TARGET_LINK_LIBRARIES("tutorial_decoupled" ${EwomsLinkLibraries})

ADD_EXECUTABLE("tutorial_coupled" tutorial_coupled.cc)
TARGET_LINK_LIBRARIES("tutorial_coupled" ${EwomsLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${EwomsLinkDirectories})
INCLUDE_DIRECTORIES(${EwomsIncludeDirectories})

