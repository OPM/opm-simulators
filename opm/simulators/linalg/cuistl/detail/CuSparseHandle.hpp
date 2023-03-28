#ifndef OPM_CUSPARSEHANDLE_HPP
#define OPM_CUSPARSEHANDLE_HPP
#include <cusparse.h>
#include <memory>

namespace Opm::cuistl::detail
{

/**
 * @brief The CuSparseHandle class provides a singleton for the simulator universal cuSparseHandle.
 *
 * Example use:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/CuSparseHandle.hpp>
 * void someFunction() {
 *     auto& cuSparseHandle = ::Opm::cuistl::detail::CuSparseHandle::getInstance();
 *     int cuSparseVersion = -1;
 *     OPM_CUSPARSE_SAFE_CALL(cusparseGetVersion(cuSparseHandle.get(), &cuSparseVersion));
 * }
 * @endcode
 */
class CuSparseHandle
{
public:
    // This should not be copyable.
    CuSparseHandle(const CuSparseHandle&) = delete;
    CuSparseHandle& operator=(const CuSparseHandle&) = delete;

    /**
     * Calls cuSparseDestroy on the handle
     */
    ~CuSparseHandle();

    /**
     * @brief get returns the underlying cuSparse handle (to be used in calls to cusparse)
     */
    cusparseHandle_t get();

    /**
     * @brief getInstance creates (if necessary) and returns the single unique instance of CuSparseHandle (singleton)
     */
    static CuSparseHandle& getInstance();

private:
    CuSparseHandle();
    cusparseHandle_t m_handle;
};

using CuSparseHandlePtr = std::unique_ptr<CuSparseHandle>;
} // namespace Opm::cuistl::detail
#endif // OPM_CUSPARSEHANDLE_HPP
