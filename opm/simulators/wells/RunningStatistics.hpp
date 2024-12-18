/*
  Copyright 2024 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_RUNNING_STATISTICS_HPP
#define OPM_RUNNING_STATISTICS_HPP

#include <cmath>
#include <limits>
#include <optional>

namespace Opm {

/// Facility for calculating simple sample statistics without having full
/// sample available.
///
/// \tparam Scalar Sample element type.  Typically a built-in arithmetic
/// type like \c float or \c double.
template <typename Scalar>
class RunningStatistics
{
public:
    /// Convert between byte array and object representation.
    ///
    /// \tparam Serializer Byte array conversion protocol.
    ///
    /// \param[in,out] serializer Byte array conversion object.
    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(this->sampleSize_);
        serializer(this->min_);
        serializer(this->max_);
        serializer(this->mean_);
        serializer(this->totalVariance_);
    }

    /// Create a serialisation test object.
    static RunningStatistics serializationTestObject()
    {
        auto stat = RunningStatistics{};

        stat.sampleSize_ = 12;
        stat.min_ = -static_cast<Scalar>(1);
        stat.max_ =  static_cast<Scalar>(2);
        stat.mean_ = static_cast<Scalar>(0.03);
        stat.totalVariance_ = static_cast<Scalar>(0.4);

        return stat;
    }

    /// Equality predicate.
    ///
    /// \param[in] that Object against which \code *this \endcode will
    /// be tested for equality.
    ///
    /// \return Whether or not \code *this \endcode is the same as \p that.
    bool operator==(const RunningStatistics& that) const
    {
        return (this->sampleSize_ == that.sampleSize_)
            && (this->min_ == that.min_)
            && (this->max_ == that.max_)
            && (this->mean_ == that.mean_)
            && (this->totalVariance_ == that.totalVariance_)
            ;
    }

    /// Reset internal counters to prepare for calculating a new set of
    /// sample statistics.
    void reset()
    {
        this->sampleSize_ = 0;
        this->min_ = std::numeric_limits<Scalar>::max();
        this->max_ = std::numeric_limits<Scalar>::lowest();
        this->mean_ = Scalar{};
        this->totalVariance_ = Scalar{};
    }

    /// Include new element into sample.
    ///
    /// Updates internal statistics counters.
    ///
    /// \param[in] x Sample point.
    void addSamplePoint(const Scalar x)
    {
        if (x < this->min_) { this->min_ = x; }
        if (x > this->max_) { this->max_ = x; }

        // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

        ++this->sampleSize_;

        const auto d1 = x - this->mean();
        this->mean_ += d1 / this->sampleSize_;

        const auto d2 = x - this->mean();
        this->totalVariance_ += d1 * d2;
    }

    /// Retrieve current sample size.
    ///
    /// Effectively returns the number of calls to addSamplePoint() since
    /// object was constructed or since the previous call to reset().
    std::size_t sampleSize() const { return this->sampleSize_; }

    /// Retrieve smallest sample value seen so far.
    Scalar min() const { return this->min_; }

    /// Retrieve largest sample value seen so far.
    Scalar max() const { return this->max_; }

    /// Retrieve arithmetic average of all sample points seen so far.
    Scalar mean() const { return this->mean_; }

    /// Retrieve unbiased standard deviation of all sample points seen so
    /// far.
    ///
    /// Returns nullopt if number of sample points is less than two.
    std::optional<Scalar> stdev() const
    {
        if (this->sampleSize_ < 2) {
            return {};
        }

        using std::sqrt;
        return sqrt(this->totalVariance_ / (this->sampleSize_ - 1));
    }

private:
    /// Current sample size.
    std::size_t sampleSize_{};

    /// Smallest sample value seen so far.
    Scalar min_ { std::numeric_limits<Scalar>::max() };

    /// Largest sample value seen so far.
    Scalar max_ { std::numeric_limits<Scalar>::lowest() };

    /// Arithmetic average of all sample points seen so far.
    Scalar mean_{};

    /// Variance measure.  In particular, N-1 * Var{x_i}.
    Scalar totalVariance_{};
};

} // namespace Opm

#endif // OPM_RUNNING_STATISTICS_HPP
