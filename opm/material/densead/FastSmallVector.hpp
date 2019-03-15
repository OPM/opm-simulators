#include <array>
#include <cstring>

namespace Opm {
namespace DenseAd {

template <typename ValueType, unsigned N>
class FastSmallVector
{
public:

    FastSmallVector();
    explicit FastSmallVector(const size_t num_elem);
    FastSmallVector(const size_t num_elem, const ValueType value);
    FastSmallVector(const FastSmallVector& other);
    FastSmallVector(FastSmallVector&& other);
    ~FastSmallVector();

    FastSmallVector& operator=(const FastSmallVector& other);
    FastSmallVector& operator=(FastSmallVector&& other);

    ValueType& operator[](size_t index);
    const ValueType& operator[](size_t index) const;
    size_t size() const;

private:
    std::array<ValueType, N> small_buffer_;
    std::size_t size_;
    ValueType* data_ptr_;
};

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>::
FastSmallVector()
    : size_(0)
    , data_ptr_(small_buffer_.data())
{
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>::
FastSmallVector(const size_t num_elem)
    : size_(num_elem)
    , data_ptr_(small_buffer_.data())
{
    if (size_ > N) {
        data_ptr_ = new ValueType[num_elem];
    }
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>::
FastSmallVector(const size_t num_elem, const ValueType value)
    : FastSmallVector(num_elem)
{
    std::fill_n(data_ptr_, size_, value);
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>::
FastSmallVector(FastSmallVector&& other)
   : size_ (other.size_)
{
    small_buffer_ = std::move(other.small_buffer_);
    if (size_ <= N) {
        data_ptr_ = small_buffer_.data();
    } else {
        data_ptr_ = other.data_ptr_;
    }

    other.data_ptr_= nullptr;
    other.size_ = 0;
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>::
FastSmallVector(const FastSmallVector& other)
   : small_buffer_(other.small_buffer_)
   , size_(other.size_)
{
    if (size_ <= N) {
        data_ptr_ = small_buffer_.data();
    } else {
        data_ptr_ = new ValueType[size_];
        memcpy(data_ptr_, other.data_ptr_, size_ * sizeof(ValueType));
    }
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>&
FastSmallVector<ValueType, N>::
operator = (const FastSmallVector&other)
{
    small_buffer_ = other.small_buffer_;
    size_ = other.size_;

    if (size_ <= N) {
        data_ptr_ = small_buffer_.data();
    } else {
        data_ptr_ = new ValueType[size_];
        memcpy(data_ptr_, other.data_ptr_, size_ * sizeof(ValueType));
    }

    return (*this);
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>&
FastSmallVector<ValueType, N>::
operator = (FastSmallVector&& other)
{
    if (data_ptr_ != small_buffer_.data() && data_ptr_ != nullptr) {
        delete [] data_ptr_;
    }
    size_ = other.size_;

    small_buffer_ = std::move(other.small_buffer_);
    if (size_ <= N) {
        data_ptr_ = small_buffer_.data();
    } else {
        data_ptr_ = other.data_ptr_;
    }

    other.data_ptr_ = nullptr;
    other.size_ = 0;

    return (*this);
}

template <typename ValueType, unsigned N>
FastSmallVector<ValueType, N>::
~FastSmallVector()
{
    if ((data_ptr_ != small_buffer_.data()) && (data_ptr_ != nullptr)) {
        delete [] data_ptr_;
    }
}

template <typename ValueType, unsigned N>
ValueType&
FastSmallVector<ValueType, N>::
operator[](size_t index)
{
  return data_ptr_[index];
}

template <typename ValueType, unsigned N>
const ValueType&
FastSmallVector<ValueType, N>::
operator[](size_t index) const
{
  return data_ptr_[index];
}

template <typename ValueType, unsigned N>
size_t
FastSmallVector<ValueType, N>::
size() const
{
    return size_;
}

}
}
