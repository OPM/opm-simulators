#ifndef CUSPARSE_SAFE_CALL_HPP
#define CUSPARSE_SAFE_CALL_HPP
#include <cusparse.h>
#include <exception>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/Logger.hpp>
#include <sstream>

#define OPM_CUSPARSE_SAFE_CALL(expression)                                                                             \
    {                                                                                                                  \
        cusparseStatus_t error = expression;                                                                           \
        if (error != CUSPARSE_STATUS_SUCCESS) {                                                                        \
            std::stringstream message;                                                                                 \
            message << "cuSparse expression did not execute correctly. Expression was: \n";                            \
            message << "    " << #expression << "\n";                                                                  \
            message << "in function " << __func__ << ", in " << __FILE__ << " at line " << __LINE__ << "\n";           \
            OpmLog::error(message.str());                                                                              \
            OPM_THROW(std::runtime_error, message.str());                                                              \
        }                                                                                                              \
    }
#endif // CUSPARSE_SAFE_CALL_HPP
