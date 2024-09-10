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

#ifndef UNSCALED_SATFUNC_CHECKPOINT_HPP_INCLUDED
#define UNSCALED_SATFUNC_CHECKPOINT_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp>

#include <cstddef>
#include <optional>
#include <unordered_set>
#include <vector>

namespace Opm {
    template <typename Scalar>
    struct EclEpsScalingPointsInfo;
} // namespace Opm

namespace Opm::satfunc {
    struct RawTableEndPoints;
    struct RawFunctionValues;
} // namespace Opm::satfunc

namespace Opm::Satfunc::PhaseChecks {

    /// Callbacks for defining the consistency check point of a single region.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class UnscaledSatfuncCheckPoint : public SatfuncCheckPointInterface<Scalar>
    {
    public:
        /// Collection of saturation function end-points and function values
        /// extracted from tabulated saturation functions.
        struct UnscaledEndPoints
        {
            /// Raw table end-points.  Minimum, critical, and maximum
            /// saturation points for each phase for all tabulated
            /// saturation functions.
            const satfunc::RawTableEndPoints* rtep{nullptr};

            /// Raw saturation function values.  Maximum function values for
            /// all saturation functions in addition to relative
            /// permeability values at critical saturation points.
            const satfunc::RawFunctionValues* rfunc{nullptr};
        };

        /// Constructor
        ///
        /// \param[in] region Region index for each active cell on rank.
        ///
        /// \param[in] regIdxOffset Region index offset.  Pass one (1) if \p
        /// region contains one-based region indices.
        ///
        /// \param[in] unscaledEndPoints Saturation function end-points for
        /// all tabulated saturation functions.  Lifetime of members must
        /// exceed the \c UnscaledSatfuncCheckPoint object.
        explicit UnscaledSatfuncCheckPoint(const std::vector<int>*  region,
                                           const int                regIdxOffset,
                                           const UnscaledEndPoints& unscaledEndPoints)
            : unscaledEndPoints_ { unscaledEndPoints }
            , regIdxOffset_      { regIdxOffset }
            , region_            { region }
        {}

        /// Compute locally unique, i.e., per MPI rank, ID of this check for
        /// a particular cell index.
        ///
        /// Common examples include the drainage or imbibition region ID
        /// (i.e., SATNUM or IMBNUM) or the Cartesian block index of a cell.
        ///
        /// \param[in] cellIdx Active cell index on current rank.
        ///
        /// \return Locally unique point ID for \p cellIdx.  Nullopt if this
        /// check point does not apply to \p cellIdx.  Typically because the
        /// underlying region of \p cellIdx has already been visited.
        std::optional<std::size_t> pointID(const int cellIdx) const override;

        /// Populate check point values for a particular cell.
        ///
        /// \param[in] cellIdx Active cell index on current rank.
        ///
        /// \param[out] endPoints Set of saturation function end-points.
        void populateCheckPoint(const int                        cellIdx,
                                EclEpsScalingPointsInfo<Scalar>& endPoints) const override;

    private:
        /// Saturation function end-points for all tabulated saturation
        /// functions.  Lifetime of members must exceed the \c
        /// UnscaledSatfuncCheckPoint object.
        UnscaledEndPoints unscaledEndPoints_;

        /// Region index offset.  Should be one (1) if region_ contains
        /// one-based region indices.
        int regIdxOffset_{};

        /// Region index for each active cell on rank.
        const std::vector<int>* region_{nullptr};

        /// Cache for visited ("seen") saturation regions.
        mutable std::unordered_set<int> seen_{};
    };

} // namespace Opm::Satfunc::PhaseChecks

#endif // UNSCALED_SATFUNC_CHECKPOINT_HPP_INCLUDED
