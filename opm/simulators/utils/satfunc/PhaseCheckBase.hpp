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

#ifndef PHASE_CHECK_BASE_HPP_INCLUDED
#define PHASE_CHECK_BASE_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

namespace Opm::Satfunc::PhaseChecks {

    /// Base class for all phase saturation function consistency checks.
    ///
    /// Provides common implementation of parts of the public Check
    /// interface in terms of packed flags.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class PhaseCheckBase : public SatfuncConsistencyChecks<Scalar>::Check
    {
    public:
        /// Run specific check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void test(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;

        /// Whether or not last set of end-points violated this particular
        /// check.
        bool isViolated() const override;

        /// Whether or not this check is critical to the simulator's ability
        /// to run the case.
        ///
        /// Violating critical checks should typically stop the run.
        bool isCritical() const override;

    protected:
        /// Mark check as violated.
        ///
        /// Intended to be called by derived types only.
        void setViolated();

        /// Mark check as violated at critical level.
        ///
        /// Intended to be called by derived types only.
        void setCritical();

        /// Put epsilon margins in checks to accept rounding errors
        static constexpr Scalar epsilon_{1e-5};

    private:
        /// Collection of violation flags.
        ///
        /// Packed bit fields.
        unsigned char flags_{0};

        /// Run specific check against a set of saturation function end-points.
        ///
        /// Actual test function.  Implemented in derived types and called
        /// from the test() member function.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        virtual void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) = 0;
    };

} // namespace Opm::Satfunc::PhaseChecks

#endif // PHASE_CHECK_BASE_HPP_INCLUDED
