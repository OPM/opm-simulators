# programs just to build when "make check" is used
check_PROGRAMS = tutorial_decoupled tutorial_coupled

noinst_HEADERS= *.hh
EXTRA_DIST = CMakeLists.txt

tutorial_decoupleddir = $(datadir)/ewoms/tutorial
tutorial_decoupled_SOURCES = tutorial_decoupled.cc
tutorial_decoupled_DATA = $(tutorial_decoupled_SOURCES)

tutorial_coupleddir = $(datadir)/ewoms/tutorial
tutorial_coupled_SOURCES = tutorial_coupled.cc
tutorial_coupled_DATA = $(tutorial_coupled_SOURCES)

include $(top_srcdir)/am/global-rules
