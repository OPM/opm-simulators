add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("tutorial_decoupled" tutorial_decoupled.cc)
TARGET_LINK_LIBRARIES("tutorial_decoupled" ${DumuxLinkLibraries})

ADD_EXECUTABLE("tutorial_coupled" tutorial_coupled.cc)
TARGET_LINK_LIBRARIES("tutorial_coupled" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

