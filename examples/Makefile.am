AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)

LDADD = $(top_builddir)/libopmpolymer.la	\
$(BOOST_FILESYSTEM_LIB)				\
$(BOOST_SYSTEM_LIB)	                        \


noinst_PROGRAMS =				\
polymer_reorder					\
sim_poly2p_incomp_reorder                  	\
test_singlecellsolves

polymer_reorder_SOURCES =			\
polymer_reorder.cpp

sim_poly2p_incomp_reorder_SOURCES =		\
sim_poly2p_incomp_reorder.cpp

test_singlecellsolves_SOURCES =			\
test_singlecellsolves.cpp

