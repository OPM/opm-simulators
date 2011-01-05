# programs just to build when "make check" is used
check_PROGRAMS = tutorial_decoupled tutorial_coupled

noinst_HEADERS= tutorialspatialparameters_decoupled.hh tutorialproblem_decoupled.hh tutorialspatialparameters_coupled.hh tutorialproblem_coupled.hh
EXTRA_DIST = CMakeLists.txt

tutorial_decoupleddir = $(datadir)/dumux/tutorial
tutorial_decoupled_HEADERS = tutorialspatialparameters_decoupled.hh \
	tutorialproblem_decoupled.hh
tutorial_decoupled_SOURCES = tutorial_decoupled.cc
tutorial_decoupled_DATA = $(tutorial_decoupled_SOURCES)
tutorial_coupleddir = $(datadir)/dumux/tutorial
tutorial_coupled_HEADERS = tutorialspatialparameters_coupled.hh \
	tutorialproblem_coupled.hh
tutorial_coupled_SOURCES = tutorial_coupled.cc
tutorial_coupled_DATA = $(tutorial_coupled_SOURCES)

include $(top_srcdir)/am/global-rules
