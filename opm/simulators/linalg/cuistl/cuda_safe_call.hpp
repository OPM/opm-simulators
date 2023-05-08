#ifndef CUDA_SAFE_CALL_HPP
#define CUDA_SAFE_CALL_HPP

#include <opm/common/OpmLog/Logger.hpp>
#include <sstream>

#define OPM_CUDA_SAFE_CALL(expression) \
    cudaError_t error = expression; \
    if (error != cudaSuccess) { \
        std::stringstream message; \
        message << "CUDA expression did not execute correctly. Expression was: \n"; \
        message << "    " << #expression << "\n"; \
        message << "in function " << __func__ << ", in " << __FILE__ << " at line " << __LINE__ <<"\n"; \
        OpmLog::error(message.str()); \
        OPM_THROW(std::runtime_error, message.str()); \
    }
#endif
