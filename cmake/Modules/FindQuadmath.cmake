# Module that checks whether the compiler supports the
# quadruple precision floating point math
#
# Sets the following variable:
# HAVE_QUAD
#
# perform tests
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
list(APPEND CMAKE_REQUIRED_LIBRARIES "quadmath")
CHECK_C_SOURCE_COMPILES("
#include <quadmath.h>

int main(void){
    __float128 foo = sqrtq(123.456);
}" HAVE_QUAD)
cmake_pop_check_state()
