ACLOCAL_AMFLAGS = -I m4

SUBDIRS = . examples


# Declare libraries
lib_LTLIBRARIES = libopmpolymer.la

libopmpolymer_la_SOURCES =				\
opm/polymer/CompressibleTpfaPolymer.cpp			\
opm/polymer/IncompTpfaPolymer.cpp			\
opm/polymer/PolymerInflow.cpp				\
opm/polymer/PolymerProperties.cpp			\
opm/polymer/SimulatorCompressiblePolymer.cpp		\
opm/polymer/SimulatorPolymer.cpp			\
opm/polymer/TransportModelCompressiblePolymer.cpp	\
opm/polymer/TransportModelPolymer.cpp			\
opm/polymer/polymerUtilities.cpp

nobase_include_HEADERS =				\
opm/polymer/CompressibleTpfaPolymer.hpp			\
opm/polymer/GravityColumnSolverPolymer.hpp		\
opm/polymer/GravityColumnSolverPolymer_impl.hpp		\
opm/polymer/IncompPropertiesDefaultPolymer.hpp		\
opm/polymer/IncompTpfaPolymer.hpp			\
opm/polymer/PolymerBlackoilState.hpp			\
opm/polymer/PolymerInflow.hpp				\
opm/polymer/PolymerProperties.hpp			\
opm/polymer/PolymerState.hpp				\
opm/polymer/SimulatorCompressiblePolymer.hpp		\
opm/polymer/SimulatorPolymer.hpp			\
opm/polymer/SinglePointUpwindTwoPhasePolymer.hpp	\
opm/polymer/TransportModelCompressiblePolymer.hpp	\
opm/polymer/TransportModelPolymer.hpp			\
opm/polymer/polymerUtilities.hpp
