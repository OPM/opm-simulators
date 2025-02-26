// TODO: Add license header
#ifndef OPM_SIMULATORS_OPM_SIMULATORS_LINALG_GPUISTL_DENSE_FIELDVECTOR_HPP
#define OPM_SIMULATORS_OPM_SIMULATORS_LINALG_GPUISTL_DENSE_FIELDVECTOR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/common/boundschecking.hh>
#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/math.hh>

#include <dune/common/math.hh>
#include <dune/common/promotiontraits.hh>

#include <opm/simulators/linalg/gpuistl/dense/DenseVector.hpp>


/** @addtogroup DenseMatVec
    @{
 */

/*! \file
 * \brief Implements a vector constructed from a given type
   representing a field and a compile-time given size.
 */
namespace Opm::gpuistl::dense
{
template <class K, int SIZE>
class FieldVector;
}

namespace Dune
{
template <class K, int SIZE>
struct DenseMatVecTraits<::Opm::gpuistl::dense::FieldVector<K, SIZE>> {
    typedef ::Opm::gpuistl::dense::FieldVector<K, SIZE> derived_type;
    typedef std::array<K, SIZE> container_type;
    typedef K value_type;
    typedef typename container_type::size_type size_type;
};

template <class K, int SIZE>
struct FieldTraits<::Opm::gpuistl::dense::FieldVector<K, SIZE>> {
    typedef typename ::Dune::FieldTraits<K>::field_type field_type;
    typedef typename ::Dune::FieldTraits<K>::real_type real_type;
};

} // namespace Dune
namespace Opm::gpuistl::dense
{
/**
 * @brief TMP to check the size of a DenseVectors statically, if possible.
 *
 * If the implementation type of C is  a FieldVector, we statically check
 * whether its dimension is SIZE.
 * @tparam C The implementation of the other DenseVector
 * @tparam SIZE The size we need assume.
 */
template <typename C, int SIZE>
struct IsFieldVectorSizeCorrect {
    /**
     * \brief True if C is not of type FieldVector or its dimension
     * is not equal SIZE.
     */
    constexpr static bool value = true;
};

template <typename T, int SIZE>
struct IsFieldVectorSizeCorrect<FieldVector<T, SIZE>, SIZE> {
    constexpr static bool value = true;
};

template <typename T, int SIZE, int SIZE1>
struct IsFieldVectorSizeCorrect<FieldVector<T, SIZE1>, SIZE> {
    constexpr static bool value = false;
};


/** \brief vector space out of a tensor product of fields.
 *
 * \tparam K    the field type (use float, double, complex, etc)
 * \tparam SIZE number of components.
 */
template <class K, int SIZE>
class FieldVector : public DenseVector<FieldVector<K, SIZE>>
{
    std::array<K, SIZE> _data;
    typedef DenseVector<FieldVector<K, SIZE>> Base;

public:
    //! The size of this vector.
    constexpr static int dimension = SIZE;

    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;

    /** \brief The type used for references to the vector entry */
    typedef value_type& reference;

    /** \brief The type used for const references to the vector entry */
    typedef const value_type& const_reference;

    //! Constructor making default-initialized vector
    OPM_HOST_DEVICE constexpr FieldVector()
        : _data {{}}
    {
    }

    //! Constructor making vector with identical coordinates
    OPM_HOST_DEVICE explicit FieldVector(const K& t)
    {
        std::fill(_data.begin(), _data.end(), t);
    }

#if __GNUC__ == 5 && !defined(__clang__)
    // `... = default;` causes an internal compiler error on GCC 5.4 (Ubuntu 16.04)
    //! copy constructor
    FieldVector(const FieldVector& x)
        : _data(x._data)
    {
    }
#else
    //! Copy constructor
    FieldVector(const FieldVector&) = default;
#endif

    /** \brief Construct from a std::initializer_list */
    OPM_HOST_DEVICE FieldVector(std::initializer_list<K> const& l)
    {
        assert(l.size() == dimension); // Actually, this is not needed any more!
        std::copy_n(l.begin(), std::min(static_cast<std::size_t>(dimension), l.size()), _data.begin());
    }

    //! copy assignment operator
    FieldVector& operator=(const FieldVector&) = default;

    template <typename T>
    OPM_HOST_DEVICE FieldVector& operator=(const FieldVector<T, SIZE>& x)
    {
        std::copy_n(x.begin(), SIZE, _data.begin());
        return *this;
    }

    template <typename T, int N>
    FieldVector& operator=(const FieldVector<T, N>&) = delete;

    /**
     * \brief Copy constructor from a second vector of possibly different type
     *
     * If the DenseVector type of the this constructor's argument
     * is implemented by a FieldVector, it is statically checked
     * if it has the correct size. If this is not the case
     * the constructor is removed from the overload set using SFINAE.
     *
     * \param[in]  x  A DenseVector with correct size.
     * \param[in]  dummy  A void* dummy argument needed by SFINAE.
     */
    template <class C>
    OPM_HOST_DEVICE FieldVector(const DenseVector<C>& x,
                [[maybe_unused]] typename std::enable_if<IsFieldVectorSizeCorrect<C, SIZE>::value>::type* dummy = 0)
    {
        // do a run-time size check, for the case that x is not a FieldVector
        assert(x.size() == SIZE); // Actually this is not needed any more!
        std::copy_n(x.begin(), std::min(static_cast<std::size_t>(SIZE), x.size()), _data.begin());
    }

    //! Constructor making vector with identical coordinates
    template <class K1>
    OPM_HOST_DEVICE explicit FieldVector(const FieldVector<K1, SIZE>& x)
    {
        std::copy_n(x.begin(), SIZE, _data.begin());
    }

    template <typename T, int N>
    explicit FieldVector(const FieldVector<T, N>&) = delete;

    using Base::operator=;

    // make this thing a vector
    OPM_HOST_DEVICE static constexpr size_type size()
    {
        return SIZE;
    }

    OPM_HOST_DEVICE K& operator[](size_type i)
    {
        DUNE_ASSERT_BOUNDS(i < SIZE);
        return _data[i];
    }
    OPM_HOST_DEVICE const K& operator[](size_type i) const
    {
        DUNE_ASSERT_BOUNDS(i < SIZE);
        return _data[i];
    }

    //! return pointer to underlying array
    OPM_HOST_DEVICE K* data() noexcept
    {
        return _data.data();
    }

    //! return pointer to underlying array
    OPM_HOST_DEVICE const K* data() const noexcept
    {
        return _data.data();
    }

    //! vector space multiplication with scalar
    template <class Scalar, std::enable_if_t<::Dune::IsNumber<Scalar>::value, int> = 0>
    OPM_HOST_DEVICE friend auto operator*(const FieldVector& vector, Scalar scalar)
    {
        FieldVector<typename ::Dune::PromotionTraits<value_type, Scalar>::PromotedType, SIZE> result;

        for (size_type i = 0; i < vector.size(); ++i)
            result[i] = vector[i] * scalar;

        return result;
    }

    //! vector space multiplication with scalar
    template <class Scalar, std::enable_if_t<::Dune::IsNumber<Scalar>::value, int> = 0>
    OPM_HOST_DEVICE friend auto operator*(Scalar scalar, const FieldVector& vector)
    {
        FieldVector<typename ::Dune::PromotionTraits<value_type, Scalar>::PromotedType, SIZE> result;

        for (size_type i = 0; i < vector.size(); ++i)
            result[i] = scalar * vector[i];

        return result;
    }

    //! vector space division by scalar
    template <class Scalar, std::enable_if_t<::Dune::IsNumber<Scalar>::value, int> = 0>
    OPM_HOST_DEVICE friend auto operator/(const FieldVector& vector, Scalar scalar)
    {
        FieldVector<typename ::Dune::PromotionTraits<value_type, Scalar>::PromotedType, SIZE> result;

        for (size_type i = 0; i < vector.size(); ++i)
            result[i] = vector[i] / scalar;

        return result;
    }
};

/** \brief Read a FieldVector from an input stream
 *  \relates FieldVector
 *
 *  \note This operator is STL compliant, i.e., the content of v is only
 *        changed if the read operation is successful.
 *
 *  \param[in]  in  std :: istream to read from
 *  \param[out] v   FieldVector to be read
 *
 *  \returns the input stream (in)
 */
template <class K, int SIZE>
inline std::istream&
operator>>(std::istream& in, FieldVector<K, SIZE>& v)
{
    FieldVector<K, SIZE> w;
    for (typename FieldVector<K, SIZE>::size_type i = 0; i < SIZE; ++i)
        in >> w[i];
    if (in)
        v = w;
    return in;
}
} // namespace Opm::gpuistl::dense

namespace Dune {
#ifndef DOXYGEN

template <class K>
struct DenseMatVecTraits<::Opm::gpuistl::dense::FieldVector<K, 1>> {
    typedef ::Opm::gpuistl::dense::FieldVector<K, 1> derived_type;
    typedef K container_type;
    typedef K value_type;
    typedef size_t size_type;
};
} // namespace Dune

namespace Opm::gpuistl::dense
{
/** \brief Vectors containing only one component
 */
template <class K>
class FieldVector<K, 1> : public DenseVector<FieldVector<K, 1>>
{
    K _data;
    typedef DenseVector<FieldVector<K, 1>> Base;

public:
    //! The size of this vector.
    constexpr static int dimension = 1;

    typedef typename Base::size_type size_type;

    /** \brief The type used for references to the vector entry */
    typedef K& reference;

    /** \brief The type used for const references to the vector entry */
    typedef const K& const_reference;

    //===== construction

    /** \brief Default constructor */
    OPM_HOST_DEVICE constexpr FieldVector()
        : _data()
    {
    }

    /** \brief Constructor with a given scalar */
    template <typename T,
              typename EnableIf = typename std::enable_if<
                  std::is_convertible<T, K>::value
                  && !std::is_base_of<DenseVector<typename ::Dune::FieldTraits<T>::field_type>, K>::value>::type>
    OPM_HOST_DEVICE FieldVector(const T& k)
        : _data(k)
    {
    }

    //! Constructor from static vector of different type
    template <class C, std::enable_if_t<std::is_assignable<K&, typename DenseVector<C>::value_type>::value, int> = 0>
    OPM_HOST_DEVICE FieldVector(const DenseVector<C>& x)
    {
        static_assert(((bool)IsFieldVectorSizeCorrect<C, 1>::value), "FieldVectors do not match in dimension!");
        assert(x.size() == 1);
        _data = x[0];
    }

    //! copy constructor
    FieldVector(const FieldVector&) = default;

    //! copy assignment operator
    FieldVector& operator=(const FieldVector&) = default;

    template <typename T>
    OPM_HOST_DEVICE FieldVector& operator=(const FieldVector<T, 1>& other)
    {
        _data = other[0];
        return *this;
    }

    template <typename T, int N>
    FieldVector& operator=(const FieldVector<T, N>&) = delete;

    /** \brief Construct from a std::initializer_list */
    OPM_HOST_DEVICE FieldVector(std::initializer_list<K> const& l)
    {
        assert(l.size() == 1);
        _data = *l.begin();
    }

    //! Assignment operator for scalar
    template <typename T,
              typename EnableIf = typename std::enable_if<
                  std::is_assignable<K&, T>::value
                  && !std::is_base_of<DenseVector<typename ::Dune::FieldTraits<T>::field_type>, K>::value>::type>
    OPM_HOST_DEVICE inline FieldVector& operator=(const T& k)
    {
        _data = k;
        return *this;
    }

    //===== forward methods to container
    OPM_HOST_DEVICE static constexpr size_type size()
    {
        return 1;
    }

    OPM_HOST_DEVICE K& operator[]([[maybe_unused]] size_type i)
    {
        DUNE_ASSERT_BOUNDS(i == 0);
        return _data;
    }
    OPM_HOST_DEVICE const K& operator[]([[maybe_unused]] size_type i) const
    {
        DUNE_ASSERT_BOUNDS(i == 0);
        return _data;
    }

    //! return pointer to underlying array
    OPM_HOST_DEVICE K* data() noexcept
    {
        return &_data;
    }

    //! return pointer to underlying array
    OPM_HOST_DEVICE const K* data() const noexcept
    {
        return &_data;
    }

    //===== conversion operator

    /** \brief Conversion operator */
    OPM_HOST_DEVICE operator K&()
    {
        return _data;
    }

    /** \brief Const conversion operator */
    OPM_HOST_DEVICE operator const K&() const
    {
        return _data;
    }
};

/* ----- FV / FV ----- */
/* mostly not necessary as these operations are already covered via the cast operator */

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator>(const FieldVector<K, 1>& a, const FieldVector<K, 1>& b)
{
    return a[0] > b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator>=(const FieldVector<K, 1>& a, const FieldVector<K, 1>& b)
{
    return a[0] >= b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator<(const FieldVector<K, 1>& a, const FieldVector<K, 1>& b)
{
    return a[0] < b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator<=(const FieldVector<K, 1>& a, const FieldVector<K, 1>& b)
{
    return a[0] <= b[0];
}

/* ----- FV / scalar ----- */

//! Binary addition, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator+(const FieldVector<K, 1>& a, const K b)
{
    return a[0] + b;
}

//! Binary subtraction, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator-(const FieldVector<K, 1>& a, const K b)
{
    return a[0] - b;
}

//! Binary multiplication, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator*(const FieldVector<K, 1>& a, const K b)
{
    return a[0] * b;
}

//! Binary division, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator/(const FieldVector<K, 1>& a, const K b)
{
    return a[0] / b;
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator>(const FieldVector<K, 1>& a, const K b)
{
    return a[0] > b;
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator>=(const FieldVector<K, 1>& a, const K b)
{
    return a[0] >= b;
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator<(const FieldVector<K, 1>& a, const K b)
{
    return a[0] < b;
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator<=(const FieldVector<K, 1>& a, const K b)
{
    return a[0] <= b;
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator==(const FieldVector<K, 1>& a, const K b)
{
    return a[0] == b;
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator!=(const FieldVector<K, 1>& a, const K b)
{
    return a[0] != b;
}

/* ----- scalar / FV ------ */

//! Binary addition, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator+(const K a, const FieldVector<K, 1>& b)
{
    return a + b[0];
}

//! Binary subtraction, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator-(const K a, const FieldVector<K, 1>& b)
{
    return a - b[0];
}

//! Binary multiplication, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator*(const K a, const FieldVector<K, 1>& b)
{
    return a * b[0];
}

//! Binary division, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline FieldVector<K, 1>
operator/(const K a, const FieldVector<K, 1>& b)
{
    return a / b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator>(const K a, const FieldVector<K, 1>& b)
{
    return a > b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator>=(const K a, const FieldVector<K, 1>& b)
{
    return a >= b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator<(const K a, const FieldVector<K, 1>& b)
{
    return a < b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator<=(const K a, const FieldVector<K, 1>& b)
{
    return a <= b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator==(const K a, const FieldVector<K, 1>& b)
{
    return a == b[0];
}

//! Binary compare, when using FieldVector<K,1> like K
template <class K>
OPM_HOST_DEVICE inline bool
operator!=(const K a, const FieldVector<K, 1>& b)
{
    return a != b[0];
}
#endif

/* Overloads for common classification functions */
namespace MathOverloads
{

    // ! Returns whether all entries are finite
    template <class K, int SIZE>
    OPM_HOST_DEVICE auto isFinite(const FieldVector<K, SIZE>& b, ::Dune::PriorityTag<2>, ::Dune::MathOverloads::ADLTag)
    {
        bool out = true;
        for (int i = 0; i < SIZE; i++) {
            out &= Dune::isFinite(b[i]);
        }
        return out;
    }

    // ! Returns whether any entry is infinite
    template <class K, int SIZE>
    OPM_HOST_DEVICE bool isInf(const FieldVector<K, SIZE>& b, ::Dune::PriorityTag<2>, ::Dune::MathOverloads::ADLTag)
    {
        bool out = false;
        for (int i = 0; i < SIZE; i++) {
            out |= Dune::isInf(b[i]);
        }
        return out;
    }

    // ! Returns whether any entry is NaN
    template <class K, int SIZE, typename = std::enable_if_t<::Dune::HasNaN<K>::value>>
    OPM_HOST_DEVICE bool isNaN(const FieldVector<K, SIZE>& b, ::Dune::PriorityTag<2>, ::Dune::MathOverloads::ADLTag)
    {
        bool out = false;
        for (int i = 0; i < SIZE; i++) {
            out |= Dune::isNaN(b[i]);
        }
        return out;
    }

    // ! Returns true if either b or c is NaN
    template <class K, typename = std::enable_if_t<::Dune::HasNaN<K>::value>>
    OPM_HOST_DEVICE bool isUnordered(const FieldVector<K, 1>& b, const FieldVector<K, 1>& c, ::Dune::PriorityTag<2>, ::Dune::MathOverloads::ADLTag)
    {
        return Dune::isUnordered(b[0], c[0]);
    }
} // namespace MathOverloads

/** @} end documentation */

} // namespace Opm::gpuistl::dense

#endif
