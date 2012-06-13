ACLOCAL_AMFLAGS = -I m4

SUBDIRS = . examples


# Declare libraries
lib_LTLIBRARIES = libopmpolymer.la

libopmpolymer_la_SOURCES =				\
opm/polymer/SimulatorPolymer.cpp			\
opm/polymer/TransportModelPolymer.cpp			\
opm/polymer/PolymerProperties.cpp			\
opm/polymer/polymerUtilities.cpp                  

nobase_include_HEADERS =				\
opm/polymer/GravityColumnSolverPolymer.hpp		\
opm/polymer/GravityColumnSolverPolymer_impl.hpp		\
opm/polymer/IncompPropertiesDefaultPolymer.hpp		\
opm/polymer/PolymerProperties.hpp			\
opm/polymer/PolymerState.hpp				\
opm/polymer/SinglePointUpwindTwoPhasePolymer.hpp	\
opm/polymer/SimulatorPolymer.hpp			\
opm/polymer/TransportModelPolymer.hpp			\
opm/polymer/polymerUtilities.hpp
