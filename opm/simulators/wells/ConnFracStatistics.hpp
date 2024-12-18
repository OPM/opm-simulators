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

#ifndef OPM_CONNFRACSTATISTICS_HPP
#define OPM_CONNFRACSTATISTICS_HPP

#include <opm/simulators/wells/RunningStatistics.hpp>

#include <array>
#include <cstddef>
#include <type_traits>

namespace Opm {

/// Collection of fracturing statistics measures at the connection level.
///
/// \tparam Scalar Statistics element type.  Typically a built-in arithmetic
/// type like \c float or \c double.
template <typename Scalar>
class ConnFracStatistics
{
public:
    /// Known quantities for which this collection provides statistics
    /// measures.
    enum class Quantity : std::size_t
    {
        /// Fracture pressure
        Pressure,

        /// Fracture flow rate
        FlowRate,

        /// Fracture width
        Width,

        // -------------------------------------------------------------
        // Helper.  Must be last enumerator.
        NumQuantities,
    };

    /// Sample point representation.
    ///
    /// Client code must populate an object of this type in order to collect
    /// statistics.
    using SamplePoint = std::array
        <Scalar, static_cast<std::underlying_type_t<Quantity>>(Quantity::NumQuantities)>;

    /// Convert between byte array and object representation.
    ///
    /// \tparam Serializer Byte array conversion protocol.
    ///
    /// \param[in,out] serializer Byte array conversion object.
    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(this->quantity_);
    }

    /// Reset internal counters to prepare for calculating a new set of
    /// sample statistics.
    void reset()
    {
        for (auto& q : this->quantity_) { q.reset(); }
    }

    /// Include new element into sample.
    ///
    /// Updates internal statistics counters.
    ///
    /// \param[in] samplePoint Collection of sample values for the
    /// fracturing of the current well/reservoir connection.
    void addSamplePoint(const SamplePoint& samplePoint)
    {
        for (auto qIdx = 0*samplePoint.size(); qIdx < samplePoint.size(); ++qIdx) {
            this->quantity_[qIdx].addSamplePoint(samplePoint[qIdx]);
        }
    }

    /// Retrieve collection of sample statistics for a single quantity.
    ///
    /// \param[in] q Quantity for which to retrieve sample statistics.
    ///
    /// \return Sample statistics for quantity \p q.
    const RunningStatistics<Scalar>& statistics(const Quantity q) const
    {
        return this->quantity_[ static_cast<std::underlying_type_t<Quantity>>(q) ];
    }

    /// Create a serialisation test object.
    static ConnFracStatistics serializationTestObject()
    {
        auto stat = ConnFracStatistics{};

        stat.quantity_
            .fill(RunningStatistics<Scalar>::serializationTestObject());

        return stat;
    }

    /// Equality predicate.
    ///
    /// \param[in] that Object against which \code *this \endcode will
    /// be tested for equality.
    ///
    /// \return Whether or not \code *this \endcode is the same as \p that.
    bool operator==(const ConnFracStatistics& that) const
    {
        return this->quantity_ == that.quantity_;
    }

private:
    using StatArray = std::array<
    RunningStatistics<Scalar>,
    static_cast<std::underlying_type_t<Quantity>>(Quantity::NumQuantities)
    >;

    /// Collection of connection level fracturing statistics.
    StatArray quantity_{};
};

} // namespace Opm

#endif // OPM_CONNFRACSTATISTICS_HPP
