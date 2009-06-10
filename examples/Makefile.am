# tests where program to build and program to run are equal
NORMALTESTS = tutorial_decoupled tutorial_coupled

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
bin_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = tutorial_decoupled.cc tutorial_coupled.cc

tutorial_decoupleddir = $(EXTRA_DIST=CMakeLists.txt
includedir)/dumux/tutorial
tutorial_decoupled_HEADERS = tutorial_soilproperties_decoupled.hh \
	tutorialproblem_decoupled.hh
tutorial_decoupled_SOURCES = tutorial_decoupled.cc

tutorial_coupleddir = $(EXTRA_DIST=CMakeLists.txt
includedir)/dumux/tutorial
tutorial_coupled_HEADERS = tutorial_soilproperties_coupled.hh \
	tutorialproblem_coupled.hh
tutorial_coupled_SOURCES = tutorial_coupled.cc

EXTRA_DIST=CMakeLists.txt
include $(top_srcdir)/am/global-rules
