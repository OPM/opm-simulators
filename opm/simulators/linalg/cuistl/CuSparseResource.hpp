#ifndef CUSPARSERESOURCE_HPP
#define CUSPARSERESOURCE_HPP
#include <type_traits>
#include <memory>
#include <functional>
#include <cusparse.h>

namespace Opm::cuistl {

//! @brief wraps a CuSparseResource in proper RAII.
template<class T>
class CuSparseResource
{
public:
    using CreatorType = typename std::function<cusparseStatus_t(T*)>;
    using DeleterType = typename std::function<cusparseStatus_t(T)>;
    CuSparseResource(CreatorType creator, DeleterType deleter);
    CuSparseResource();
    ~CuSparseResource();

    // This should not be copyable.
    CuSparseResource( const CuSparseResource& ) = delete;
    CuSparseResource& operator=( const CuSparseResource& ) = delete;

    T get();
private:
    // cuSparse exposes its types as typedefs of pointers.
    // eg. cusparseHandle_t = cuSparseHandle*
    using RawType = typename std::remove_pointer<T>::type;
    T resource;

    DeleterType deleter;
};

}
#endif // CUSPARSERESOURCE_HPP
