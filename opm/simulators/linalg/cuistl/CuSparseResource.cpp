#include <exception>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>

namespace Opm::cuistl
{

namespace
{
    template <class T>
    struct CuSparseDeleteAndCreate {
    };

    template <>
    struct CuSparseDeleteAndCreate<bsrilu02Info_t> {
        using DeleterType = typename CuSparseResource<bsrilu02Info_t>::DeleterType;
        using CreatorType = typename CuSparseResource<bsrilu02Info_t>::CreatorType;

        DeleterType getDeleter()
        {
            return cusparseDestroyBsrsv2Info;
        }

        CreatorType getCreator()
        {
            return cusparseCreateBsrsv2Info;
        }
    };

    template <>
    struct CuSparseDeleteAndCreate<cusparseMatDescr_t> {
        using DeleterType = typename CuSparseResource<cusparseMatDescr_t>::DeleterType;
        using CreatorType = typename CuSparseResource<cusparseMatDescr_t>::CreatorType;

        DeleterType getDeleter()
        {
            return cusparseDestroyMatDescr;
        }

        CreatorType getCreator()
        {
            return cusparseCreateMatDescr;
        }
    };

} // namespace
template <class T>
CuSparseResource<T>::CuSparseResource(CreatorType creator, DeleterType deleter)
    : deleter(deleter)
{
    // TODO: This should probably not use this macro since it will disguise the
    // proper name of the function being called.
    OPM_CUSPARSE_SAFE_CALL(creator(&resource));
}

template <class T>
CuSparseResource<T>::CuSparseResource()
    : CuSparseResource<T>(getCreator<T>(), getDeleter<T>())
{
}

template <class T>
CuSparseResource<T>::~CuSparseResource()
{
    // TODO: This should probably not use this macro since it will disguise the
    // proper name of the function being called.
    OPM_CUSPARSE_SAFE_CALL(deleter(resource));
}
} // namespace Opm::cuistl
