AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)

LDADD = $(top_builddir)/libopmpolymer.la

noinst_PROGRAMS =				\
hello						\
polymer_reorder

hello_SOURCES =					\
SayHello.cpp

polymer_reorder_SOURCES =			\
polymer_reorder.cpp				\
Utilities.hpp					\
Utilities.cpp

polymer_reorder_LDADD =				\
$(BOOST_FILESYSTEM_LIB)				\
$(BOOST_SYSTEM_LIB)				\
$(LDADD)

