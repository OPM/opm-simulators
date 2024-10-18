/*
  Copyright 2024 Equinor AS

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

#ifndef SATFUNC_CHECKPOINT_INTERFACE_HPP_INCLUDED
#define SATFUNC_CHECKPOINT_INTERFACE_HPP_INCLUDED

#include <cstddef>
#include <optional>

namespace Opm {
    template <typename Scalar>
    struct EclEpsScalingPointsInfo;
} // namespace Opm

namespace Opm::Satfunc::PhaseChecks {

    /// Callback protocol for single saturation function consistency check point.
    ///
    /// Intended to be used as a base class.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    struct SatfuncCheckPointInterface
    {
        /// Virtual destructor for public inheritance
        virtual ~SatfuncCheckPointInterface() = default;

        /// Compute locally unique, i.e., per MPI rank, ID of this check for
        /// a particular cell index.
        ///
        /// Common examples include the drainage or imbibition region ID
        /// (i.e., SATNUM or IMBNUM) or the Cartesian block index of a cell.
        ///
        /// \param[in] cellIdx Active cell index on current rank.
        ///
        /// \return Locally unique point ID for \p cellIdx.  Nullopt if this
        /// check point does not apply to \p cellIdx.  A common cause of the
        /// latter is running a region based check and the region already
        /// having been visited.
        virtual std::optional<std::size_t>
        pointID(const int cellIdx) const = 0;

        /// Populate check point values for a particular cell.
        ///
        /// \param[in] cellIdx Active cell index on current rank.
        ///
        /// \param[out] endPoints Set of saturation function end-points.
        /// Member function populateCheckPoint() assigns all data members
        /// and derived classes must abide by this requirement.
        virtual void
        populateCheckPoint(const int                        cellIdx,
                           EclEpsScalingPointsInfo<Scalar>& endPoints) const = 0;
    };

} // namespace Opm::Satfunc::PhaseChecks

#endif // SATFUNC_CHECKPOINT_INTERFACE_HPP_INCLUDED
