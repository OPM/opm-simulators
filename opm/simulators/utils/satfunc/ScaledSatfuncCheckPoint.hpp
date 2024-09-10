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

#ifndef SCALED_SATFUNC_CHECKPOINT_HPP_INCLUDED
#define SCALED_SATFUNC_CHECKPOINT_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp>
#include <opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.hpp>

#include <cstddef>
#include <functional>
#include <optional>

namespace Opm {
    class EclipseState;
    class EclEpsGridProperties;

    template <typename Scalar>
    struct EclEpsScalingPointsInfo;
} // namespace Opm

namespace Opm::Satfunc::PhaseChecks {

    /// Callbacks for defining the scaled saturation function consistency
    /// check point of a single active grid block.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class ScaledSatfuncCheckPoint : public SatfuncCheckPointInterface<Scalar>
    {
    public:
        /// Callback for translating active cell index to globally unique
        /// point ID.
        using LocalToGlobal = std::function<std::size_t(const int)>;

        /// Constructor
        ///
        /// \param[in] unscaled Callbacks for inferring the end-points of
        /// the underlying saturation region.
        ///
        /// \param[in] eclipseState Container of static properties such as
        /// the scaled saturation function end-points.
        ///
        /// \param[in] epsGridProps Access interface for scaled saturation
        /// function end-points.
        ///
        /// \param[in] localToGlobal Callback for translating active cell
        /// indices to globally unique point IDs.
        explicit ScaledSatfuncCheckPoint(const UnscaledSatfuncCheckPoint<Scalar>& unscaled,
                                         const EclipseState*                      eclipseState,
                                         const EclEpsGridProperties*              epsGridProps,
                                         const LocalToGlobal&                     localToGlobal)
            : unscaled_      { unscaled }
            , eclipseState_  { eclipseState }
            , epsGridProps_  { epsGridProps }
            , localToGlobal_ { localToGlobal }
        {}

        /// Compute global unique, i.e., across all MPI ranks, ID of this
        /// check for a particular cell index.
        ///
        /// \param[in] cellIdx Active cell index on current rank.
        ///
        /// \return Globally unique point ID for \p cellIdx
        std::optional<std::size_t> pointID(const int cellIdx) const override
        {
            return { this->localToGlobal_(cellIdx) };
        }

        /// Populate check point values for a particular cell.
        ///
        /// \param[in] cellIdx Active cell index on current rank.
        ///
        /// \param[out] endPoints Set of saturation function end-points.
        void populateCheckPoint(const int                        cellIdx,
                                EclEpsScalingPointsInfo<Scalar>& endPoints) const override;

    private:
        /// Callbacks for inferring the end-points of the underlying
        /// saturation region.
        UnscaledSatfuncCheckPoint<Scalar> unscaled_;

        /// Container of static properties such as the scaled saturation
        /// function end-points.
        const EclipseState* eclipseState_{nullptr};

        /// Access interface for scaled saturation function end-points.
        const EclEpsGridProperties* epsGridProps_{nullptr};

        /// Callback for translating active cell indices to globally unique
        /// point IDs.
        LocalToGlobal localToGlobal_{};
    };

} // namespace Opm::Satfunc::PhaseChecks

#endif // SCALED_SATFUNC_CHECKPOINT_HPP_INCLUDED
