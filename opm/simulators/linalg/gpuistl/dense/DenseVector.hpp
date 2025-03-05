// TODO: Copyright header
#ifndef OPM_SIMULATORS_OPM_SIMULATORS_LINALG_GPUISTL_DENSE_DENSEVECTOR_HPP
#define OPM_SIMULATORS_OPM_SIMULATORS_LINALG_GPUISTL_DENSE_DENSEVECTOR_HPP

#include <algorithm>
#include <limits>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/boundschecking.hh>
#include <dune/common/dotproduct.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/genericiterator.hh>
#include <dune/common/matvectraits.hh>
#include <dune/common/promotiontraits.hh>

#include <opm/common/utility/gpuDecorators.hpp>

namespace Opm::gpuistl::dense
{
// forward declaration of template

template <typename V>
class DenseVector;

} // namespace Opm::gpuistl::dense
namespace Dune
{


template <typename V>
struct FieldTraits<::Opm::gpuistl::dense::DenseVector<V>> {
    typedef typename FieldTraits<typename ::Dune::DenseMatVecTraits<V>::value_type>::field_type field_type;
    typedef typename FieldTraits<typename ::Dune::DenseMatVecTraits<V>::value_type>::real_type real_type;
};

/** @defgroup DenseMatVec Dense Matrix and Vector Template Library
    @ingroup Common
    @{
 */

/*! \file
 * \brief Implements the dense vector interface, with an exchangeable storage class
 */
} // namespace Dune

namespace Opm::gpuistl::dense
{
namespace fvmeta
{
    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K>
    OPM_HOST_DEVICE inline typename ::Dune::FieldTraits<K>::real_type absreal(const K& k)
    {
        using std::abs;
        return abs(k);
    }

    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K>
    OPM_HOST_DEVICE inline typename ::Dune::FieldTraits<K>::real_type absreal(const std::complex<K>& c)
    {
        using std::abs;
        return abs(c.real()) + abs(c.imag());
    }

    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K>
    OPM_HOST_DEVICE inline typename ::Dune::FieldTraits<K>::real_type abs2(const K& k)
    {
        return k * k;
    }

    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K>
    OPM_HOST_DEVICE inline typename ::Dune::FieldTraits<K>::real_type abs2(const std::complex<K>& c)
    {
        return c.real() * c.real() + c.imag() * c.imag();
    }

    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K, bool isInteger = std::numeric_limits<K>::is_integer>
    struct Sqrt {
        OPM_HOST_DEVICE static inline typename ::Dune::FieldTraits<K>::real_type sqrt(const K& k)
        {
            using std::sqrt;
            return sqrt(k);
        }
    };

    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K>
    struct Sqrt<K, true> {
        OPM_HOST_DEVICE static inline typename ::Dune::FieldTraits<K>::real_type sqrt(const K& k)
        {
            using std::sqrt;
            return typename ::Dune::FieldTraits<K>::real_type(sqrt(double(k)));
        }
    };

    /**
       \private
       \memberof Dune::DenseVector
     */
    template <class K>
    OPM_HOST_DEVICE inline typename ::Dune::FieldTraits<K>::real_type sqrt(const K& k)
    {
        return Sqrt<K>::sqrt(k);
    }

} // namespace fvmeta

/*! \brief Generic iterator class for dense vector and matrix implementations

   provides sequential access to DenseVector, FieldVector and FieldMatrix
 */
template <class C, class T, class R = T&>
class DenseIterator : public Dune::RandomAccessIteratorFacade<DenseIterator<C, T, R>, T, R, std::ptrdiff_t>
{
    friend class DenseIterator<typename std::remove_const<C>::type,
                               typename std::remove_const<T>::type,
                               typename ::Dune::mutable_reference<R>::type>;
    friend class DenseIterator<const typename std::remove_const<C>::type,
                               const typename std::remove_const<T>::type,
                               typename ::Dune::const_reference<R>::type>;

    typedef DenseIterator<typename std::remove_const<C>::type,
                          typename std::remove_const<T>::type,
                          typename ::Dune::mutable_reference<R>::type>
        MutableIterator;
    typedef DenseIterator<const typename std::remove_const<C>::type,
                          const typename std::remove_const<T>::type,
                          typename ::Dune::const_reference<R>::type>
        ConstIterator;

public:
    /**
     * @brief The type of the difference between two positions.
     */
    typedef std::ptrdiff_t DifferenceType;

    /**
     * @brief The type to index the underlying container.
     */
    typedef typename C::size_type SizeType;

    // Constructors needed by the base iterators.
    OPM_HOST_DEVICE DenseIterator()
        : container_(0)
        , position_()
    {
    }

    OPM_HOST_DEVICE DenseIterator(C& cont, SizeType pos)
        : container_(&cont)
        , position_(pos)
    {
    }

    OPM_HOST_DEVICE DenseIterator(const MutableIterator& other)
        : container_(other.container_)
        , position_(other.position_)
    {
    }

    OPM_HOST_DEVICE DenseIterator(const ConstIterator& other)
        : container_(other.container_)
        , position_(other.position_)
    {
    }

    // Methods needed by the forward iterator
    OPM_HOST_DEVICE bool equals(const MutableIterator& other) const
    {
        return position_ == other.position_ && container_ == other.container_;
    }


    OPM_HOST_DEVICE bool equals(const ConstIterator& other) const
    {
        return position_ == other.position_ && container_ == other.container_;
    }

    OPM_HOST_DEVICE R dereference() const
    {
        return container_->operator[](position_);
    }

    OPM_HOST_DEVICE void increment()
    {
        ++position_;
    }

    // Additional function needed by BidirectionalIterator
    OPM_HOST_DEVICE void decrement()
    {
        --position_;
    }

    // Additional function needed by RandomAccessIterator
    OPM_HOST_DEVICE R elementAt(DifferenceType i) const
    {
        return container_->operator[](position_ + i);
    }

    OPM_HOST_DEVICE void advance(DifferenceType n)
    {
        position_ = position_ + n;
    }

    OPM_HOST_DEVICE DifferenceType distanceTo(
        DenseIterator<const typename std::remove_const<C>::type, const typename std::remove_const<T>::type> other) const
    {
        assert(other.container_ == container_);
        return static_cast<DifferenceType>(other.position_) - static_cast<DifferenceType>(position_);
    }

    OPM_HOST_DEVICE DifferenceType
    distanceTo(DenseIterator<typename std::remove_const<C>::type, typename std::remove_const<T>::type> other) const
    {
        assert(other.container_ == container_);
        return static_cast<DifferenceType>(other.position_) - static_cast<DifferenceType>(position_);
    }

    //! return index
    OPM_HOST_DEVICE SizeType index() const
    {
        return this->position_;
    }

private:
    C* container_;
    SizeType position_;
};

/** \brief Interface for a class of dense vectors over a given field.
 *
 * \tparam V implementation class of the vector
 */
template <typename V>
class DenseVector
{
    typedef ::Dune::DenseMatVecTraits<V> Traits;
    // typedef typename Traits::value_type K;

    // Curiously recurring template pattern
    OPM_HOST_DEVICE V& asImp()
    {
        return static_cast<V&>(*this);
    }
    OPM_HOST_DEVICE const V& asImp() const
    {
        return static_cast<const V&>(*this);
    }

protected:
    // construction allowed to derived classes only
    constexpr DenseVector() = default;
    // copying only allowed by derived classes
    DenseVector(const DenseVector&) = default;

public:
    //===== type definitions and constants

    //! type of derived vector class
    typedef typename Traits::derived_type derived_type;

    //! export the type representing the field
    typedef typename Traits::value_type value_type;

    //! export the type representing the field
    typedef typename ::Dune::FieldTraits<value_type>::field_type field_type;

    //! export the type representing the components
    typedef typename Traits::value_type block_type;

    //! The type used for the index access and size operation
    typedef typename Traits::size_type size_type;

    //! The number of block levels we contain. This is the leaf, that is, 1.
    constexpr static int blocklevel = 1;

    //===== assignment from scalar
    //! Assignment operator for scalar
    OPM_HOST_DEVICE inline derived_type& operator=(const value_type& k)
    {
        for (size_type i = 0; i < size(); i++)
            asImp()[i] = k;
        return asImp();
    }

    //===== assignment from other DenseVectors
protected:
    //! Assignment operator for other DenseVector of same type
    DenseVector& operator=(const DenseVector&) = default;

public:
    //! Assignment operator for other DenseVector of different type
    template <typename W,
              std::enable_if_t<std::is_assignable<value_type&, typename DenseVector<W>::value_type>::value, int> = 0>
    OPM_HOST_DEVICE derived_type& operator=(const DenseVector<W>& other)
    {
        assert(other.size() == size());
        for (size_type i = 0; i < size(); i++)
            asImp()[i] = other[i];
        return asImp();
    }

    //===== access to components

    //! random access
    OPM_HOST_DEVICE value_type& operator[](size_type i)
    {
        return asImp()[i];
    }

    OPM_HOST_DEVICE const value_type& operator[](size_type i) const
    {
        return asImp()[i];
    }

    //! return reference to first element
    OPM_HOST_DEVICE value_type& front()
    {
        return asImp()[0];
    }

    //! return reference to first element
    OPM_HOST_DEVICE const value_type& front() const
    {
        return asImp()[0];
    }

    //! return reference to last element
    OPM_HOST_DEVICE value_type& back()
    {
        return asImp()[size() - 1];
    }

    //! return reference to last element
    OPM_HOST_DEVICE const value_type& back() const
    {
        return asImp()[size() - 1];
    }

    //! checks whether the container is empty
    OPM_HOST_DEVICE bool empty() const
    {
        return size() == 0;
    }

    //! size method
    OPM_HOST_DEVICE size_type size() const
    {
        return asImp().size();
    }

    //! Iterator class for sequential access
    typedef DenseIterator<DenseVector, value_type> Iterator;
    //! typedef for stl compliant access
    typedef Iterator iterator;

    //! begin iterator
    OPM_HOST_DEVICE Iterator begin()
    {
        return Iterator(*this, 0);
    }

    //! end iterator
    OPM_HOST_DEVICE Iterator end()
    {
        return Iterator(*this, size());
    }

    //! @returns an iterator that is positioned before
    //! the end iterator of the vector, i.e. at the last entry.
    OPM_HOST_DEVICE Iterator beforeEnd()
    {
        return Iterator(*this, size() - 1);
    }

    //! @returns an iterator that is positioned before
    //! the first entry of the vector.
    OPM_HOST_DEVICE Iterator beforeBegin()
    {
        return Iterator(*this, -1);
    }

    //! return iterator to given element or end()
    OPM_HOST_DEVICE Iterator find(size_type i)
    {
        return Iterator(*this, std::min(i, size()));
    }

    //! ConstIterator class for sequential access
    typedef DenseIterator<const DenseVector, const value_type> ConstIterator;
    //! typedef for stl compliant access
    typedef ConstIterator const_iterator;

    //! begin ConstIterator
    OPM_HOST_DEVICE ConstIterator begin() const
    {
        return ConstIterator(*this, 0);
    }

    //! end ConstIterator
    OPM_HOST_DEVICE ConstIterator end() const
    {
        return ConstIterator(*this, size());
    }

    //! @returns an iterator that is positioned before
    //! the end iterator of the vector. i.e. at the last element
    OPM_HOST_DEVICE ConstIterator beforeEnd() const
    {
        return ConstIterator(*this, size() - 1);
    }

    //! @returns an iterator that is positioned before
    //! the first entry of the vector.
    OPM_HOST_DEVICE ConstIterator beforeBegin() const
    {
        return ConstIterator(*this, -1);
    }

    //! return iterator to given element or end()
    OPM_HOST_DEVICE ConstIterator find(size_type i) const
    {
        return ConstIterator(*this, std::min(i, size()));
    }

    //===== vector space arithmetic

    //! vector space addition
    template <class Other>
    OPM_HOST_DEVICE derived_type& operator+=(const DenseVector<Other>& x)
    {
        DUNE_ASSERT_BOUNDS(x.size() == size());
        for (size_type i = 0; i < size(); i++)
            (*this)[i] += x[i];
        return asImp();
    }

    //! vector space subtraction
    template <class Other>
    OPM_HOST_DEVICE derived_type& operator-=(const DenseVector<Other>& x)
    {
        DUNE_ASSERT_BOUNDS(x.size() == size());
        for (size_type i = 0; i < size(); i++)
            (*this)[i] -= x[i];
        return asImp();
    }

    //! Binary vector addition
    template <class Other>
    OPM_HOST_DEVICE derived_type operator+(const DenseVector<Other>& b) const
    {
        derived_type z = asImp();
        return (z += b);
    }

    //! Binary vector subtraction
    template <class Other>
    OPM_HOST_DEVICE derived_type operator-(const DenseVector<Other>& b) const
    {
        derived_type z = asImp();
        return (z -= b);
    }

    //! Vector negation
    OPM_HOST_DEVICE derived_type operator-() const
    {
        V result;
        using idx_type = typename decltype(result)::size_type;

        for (idx_type i = 0; i < size(); ++i)
            result[i] = -asImp()[i];

        return result;
    }

    //! \brief vector space add scalar to all comps
    /**
       we use enable_if to avoid an ambiguity, if the
       function parameter can be converted to value_type implicitly.
       (see FS#1457)

       The function is only enabled, if the parameter is directly
       convertible to value_type.
     */
    template <typename ValueType>
    OPM_HOST_DEVICE typename std::enable_if<std::is_convertible<ValueType, value_type>::value, derived_type>::type&
    operator+=(const ValueType& kk)
    {
        const value_type& k = kk;
        for (size_type i = 0; i < size(); i++)
            (*this)[i] += k;
        return asImp();
    }

    //! \brief vector space subtract scalar from all comps
    /**
       we use enable_if to avoid an ambiguity, if the
       function parameter can be converted to value_type implicitly.
       (see FS#1457)

       The function is only enabled, if the parameter is directly
       convertible to value_type.
     */
    template <typename ValueType>
    OPM_HOST_DEVICE typename std::enable_if<std::is_convertible<ValueType, value_type>::value, derived_type>::type&
    operator-=(const ValueType& kk)
    {
        const value_type& k = kk;
        for (size_type i = 0; i < size(); i++)
            (*this)[i] -= k;
        return asImp();
    }

    //! \brief vector space multiplication with scalar
    /**
       we use enable_if to avoid an ambiguity, if the
       function parameter can be converted to field_type implicitly.
       (see FS#1457)

       The function is only enabled, if the parameter is directly
       convertible to field_type.
     */
    template <typename FieldType>
    OPM_HOST_DEVICE typename std::enable_if<std::is_convertible<FieldType, field_type>::value, derived_type>::type&
    operator*=(const FieldType& kk)
    {
        const field_type& k = kk;
        for (size_type i = 0; i < size(); i++)
            (*this)[i] *= k;
        return asImp();
    }

    //! \brief vector space division by scalar
    /**
       we use enable_if to avoid an ambiguity, if the
       function parameter can be converted to field_type implicitly.
       (see FS#1457)

       The function is only enabled, if the parameter is directly
       convertible to field_type.
     */
    template <typename FieldType>
    OPM_HOST_DEVICE typename std::enable_if<std::is_convertible<FieldType, field_type>::value, derived_type>::type&
    operator/=(const FieldType& kk)
    {
        const field_type& k = kk;
        for (size_type i = 0; i < size(); i++)
            (*this)[i] /= k;
        return asImp();
    }

    //! Binary vector comparison
    template <class Other>
    OPM_HOST_DEVICE bool operator==(const DenseVector<Other>& x) const
    {
        DUNE_ASSERT_BOUNDS(x.size() == size());
        for (size_type i = 0; i < size(); i++)
            if ((*this)[i] != x[i])
                return false;

        return true;
    }

    //! Binary vector incomparison
    template <class Other>
    OPM_HOST_DEVICE bool operator!=(const DenseVector<Other>& x) const
    {
        return !operator==(x);
    }


    //! vector space axpy operation ( *this += a x )
    template <class Other>
    OPM_HOST_DEVICE derived_type& axpy(const field_type& a, const DenseVector<Other>& x)
    {
        DUNE_ASSERT_BOUNDS(x.size() == size());
        for (size_type i = 0; i < size(); i++)
            (*this)[i] += a * x[i];
        return asImp();
    }

    /**
     * \brief indefinite vector dot product \f$\left (x^T \cdot y \right)\f$ which corresponds to Petsc's VecTDot
     *
     * http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecTDot.html
     * @param x other vector
     * @return
     */
    template <class Other>
    OPM_HOST_DEVICE typename ::Dune::PromotionTraits<field_type, typename DenseVector<Other>::field_type>::PromotedType
    operator*(const DenseVector<Other>& x) const
    {
        typedef
            typename ::Dune::PromotionTraits<field_type, typename DenseVector<Other>::field_type>::PromotedType PromotedType;
        PromotedType result(0);
        assert(x.size() == size());
        for (size_type i = 0; i < size(); i++) {
            result += PromotedType((*this)[i] * x[i]);
        }
        return result;
    }

    /**
     * @brief vector dot product \f$\left (x^H \cdot y \right)\f$ which corresponds to Petsc's VecDot
     *
     * http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecDot.html
     * @param x other vector
     * @return
     */
    template <class Other>
    OPM_HOST_DEVICE typename ::Dune::PromotionTraits<field_type, typename DenseVector<Other>::field_type>::PromotedType
    dot(const DenseVector<Other>& x) const
    {
        typedef
            typename ::Dune::PromotionTraits<field_type, typename DenseVector<Other>::field_type>::PromotedType PromotedType;
        PromotedType result(0);
        assert(x.size() == size());
        for (size_type i = 0; i < size(); i++) {
            result += Dune::dot((*this)[i], x[i]);
        }
        return result;
    }

    //===== norms

    //! one norm (sum over absolute values of entries)
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<value_type>::real_type one_norm() const
    {
        using std::abs;
        typename ::Dune::FieldTraits<value_type>::real_type result(0);
        for (size_type i = 0; i < size(); i++)
            result += abs((*this)[i]);
        return result;
    }


    //! simplified one norm (uses Manhattan norm for complex values)
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<value_type>::real_type one_norm_real() const
    {
        typename ::Dune::FieldTraits<value_type>::real_type result(0);
        for (size_type i = 0; i < size(); i++)
            result += fvmeta::absreal((*this)[i]);
        return result;
    }

    //! two norm sqrt(sum over squared values of entries)
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<value_type>::real_type two_norm() const
    {
        typename ::Dune::FieldTraits<value_type>::real_type result(0);
        for (size_type i = 0; i < size(); i++)
            result += fvmeta::abs2((*this)[i]);
        return fvmeta::sqrt(result);
    }

    //! square of two norm (sum over squared values of entries), need for block recursion
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<value_type>::real_type two_norm2() const
    {
        typename ::Dune::FieldTraits<value_type>::real_type result(0);
        for (size_type i = 0; i < size(); i++)
            result += fvmeta::abs2((*this)[i]);
        return result;
    }

    //! infinity norm (maximum of absolute values of entries)
    template <typename vt = value_type, typename std::enable_if<!::Dune::HasNaN<vt>::value, int>::type = 0>
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<vt>::real_type infinity_norm() const
    {
        using real_type = typename ::Dune::FieldTraits<vt>::real_type;
        using std::abs;
        using std::max;

        real_type norm = 0;
        for (auto const& x : *this) {
            real_type const a = abs(x);
            norm = max(a, norm);
        }
        return norm;
    }

    //! simplified infinity norm (uses Manhattan norm for complex values)
    template <typename vt = value_type, typename std::enable_if<!::Dune::HasNaN<vt>::value, int>::type = 0>
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<vt>::real_type infinity_norm_real() const
    {
        using real_type = typename ::Dune::FieldTraits<vt>::real_type;
        using std::max;

        real_type norm = 0;
        for (auto const& x : *this) {
            real_type const a = fvmeta::absreal(x);
            norm = max(a, norm);
        }
        return norm;
    }

    //! infinity norm (maximum of absolute values of entries)
    template <typename vt = value_type, typename std::enable_if<::Dune::HasNaN<vt>::value, int>::type = 0>
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<vt>::real_type infinity_norm() const
    {
        using real_type = typename ::Dune::FieldTraits<vt>::real_type;
        using std::abs;
        using std::max;

        real_type norm = 0;
        real_type isNaN = 1;
        for (auto const& x : *this) {
            real_type const a = abs(x);
            norm = max(a, norm);
            isNaN += a;
        }
        return norm * (isNaN / isNaN);
    }

    //! simplified infinity norm (uses Manhattan norm for complex values)
    template <typename vt = value_type, typename std::enable_if<::Dune::HasNaN<vt>::value, int>::type = 0>
    OPM_HOST_DEVICE typename ::Dune::FieldTraits<vt>::real_type infinity_norm_real() const
    {
        using real_type = typename ::Dune::FieldTraits<vt>::real_type;
        using std::max;

        real_type norm = 0;
        real_type isNaN = 1;
        for (auto const& x : *this) {
            real_type const a = fvmeta::absreal(x);
            norm = max(a, norm);
            isNaN += a;
        }
        return norm * (isNaN / isNaN);
    }

    //===== sizes

    //! number of blocks in the vector (are of size 1 here)
    OPM_HOST_DEVICE size_type N() const
    {
        return size();
    }

    //! dimension of the vector space
    OPM_HOST_DEVICE size_type dim() const
    {
        return size();
    }
};

/** \brief Write a DenseVector to an output stream
 *  \relates DenseVector
 *
 *  \param[in]  s  std :: ostream to write to
 *  \param[in]  v  DenseVector to write
 *
 *  \returns the output stream (s)
 */
template <typename V>
std::ostream&
operator<<(std::ostream& s, const DenseVector<V>& v)
{
    for (typename DenseVector<V>::size_type i = 0; i < v.size(); i++)
        s << ((i > 0) ? " " : "") << v[i];
    return s;
}

/** @} end documentation */

} // namespace Opm::gpuistl::dense

#endif // DUNE_DENSEVECTOR_HH
