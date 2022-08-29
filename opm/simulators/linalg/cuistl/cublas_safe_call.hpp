#ifndef CUBLAS_SAFE_CALL_HPP
#define CUBLAS_SAFE_CALL_HPP
#include <cublas_v2.h>
#include <exception>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/Logger.hpp>
#include <sstream>

// This macro does nothing as of yet, but will in the future
#define OPM_CUBLAS_SAFE_CALL(expression) expression; 
#endif // CUBLAS_SAFE_CALL_HPP
