# Module that checks whether the compiler supports the
# abi::__cxa_demangle function required to
# make the type names returned by typeid() human-readable
#
# Sets the following variable:
# HAVE_CXA_DEMANGLE
#
# perform tests
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES("#include <cxxabi.h>
int main(void){
    int foobar = 0;
    const char *foo = typeid(foobar).name();
    int status;
    char *demangled = abi::__cxa_demangle( foo, 0, 0, &status );
}" HAVE_CXA_DEMANGLE)
