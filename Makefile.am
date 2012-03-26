ACLOCAL_AMFLAGS = -I m4

SUBDIRS = . examples


# Declare libraries
lib_LTLIBRARIES = libopmpolymer.la

libopmpolymer_la_SOURCES =		          \
opm/polymer/TransportModelPolymer.cpp		  \
opm/polymer/polymerUtilities.cpp

nobase_include_HEADERS =			  \
opm/polymer/GravityColumnSolverPolymer.hpp        \
opm/polymer/GravityColumnSolverPolymer_impl.hpp   \
opm/polymer/PolymerProperties.hpp		  \
opm/polymer/SinglePointUpwindTwoPhasePolymer.hpp  \
opm/polymer/TransportModelPolymer.hpp		  \
opm/polymer/polymerUtilities.hpp
