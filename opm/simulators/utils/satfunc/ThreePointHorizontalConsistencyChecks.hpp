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

#ifndef THREE_POINT_HORIZONTAL_CONSISTENCY_CHECKS_HPP_INCLUDED
#define THREE_POINT_HORIZONTAL_CONSISTENCY_CHECKS_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>
#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

#include <cstddef>
#include <string>

namespace Opm::Satfunc::PhaseChecks::ThreePointHorizontal {

    /// Verify that critical saturation of displacing phase (oil/liquid) is
    /// strictly between critical and maximum gas saturations for the
    /// alternative (three point) horizontal scaling method (SCALECRS=YES)
    /// in the gas/oil two-phase system.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class DisplacingOil_GO : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 5; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sogcr_;
            exportedCheckValues[2] = this->sgcr_;
            exportedCheckValues[3] = Scalar{1} - (this->sogcr_ + this->swl_);
            exportedCheckValues[4] = this->sgu_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile displacing oil in three point "
                     "horizontally scaled gas/oil system" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SGCR < 1-SOGCR-SWL < SGU" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SOGCR";
            headers[2] = "SGCR";
            headers[3] = "1-SOGCR-SWL";
            headers[4] = "SGU";
        }

    private:
        /// Minimum (connate) water saturation.
        Scalar swl_;

        /// Critical oil saturation in two-phase gas/oil system.
        Scalar sogcr_;

        /// Critical gas saturation.
        Scalar sgcr_;

        /// Maximum gas saturation.
        Scalar sgu_;

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that critical saturation of displacing phase (oil) is
    /// strictly between critical and maximum water saturations for the
    /// alternative (three point) horizontal scaling method (SCALECRS=YES)
    /// in the oil/water two-phase system.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class DisplacingOil_OW : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 5; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgl_;
            exportedCheckValues[1] = this->sowcr_;
            exportedCheckValues[2] = this->swcr_;
            exportedCheckValues[3] = Scalar{1} - (this->sowcr_ + this->sgl_);
            exportedCheckValues[4] = this->swu_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile displacing oil in three point "
                     "horizontally scaled oil/water system" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SWCR < 1-SOWCR-SGL < SWU" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SOWCR";
            headers[2] = "SWCR";
            headers[3] = "1-SOWCR-SGL";
            headers[4] = "SWU";
        }

    private:
        /// Minimum gas saturation.
        Scalar sgl_;

        /// Critical oil saturation in two-phase oil/water system.
        Scalar sowcr_;

        /// Critical water saturation.
        Scalar swcr_;

        /// Maximum water saturation.
        Scalar swu_;

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

} // namespace Opm::Satfunc::PhaseChecks::ThreePointHorizontal

#endif // THREE_POINT_HORIZONTAL_CONSISTENCY_CHECKS_HPP_INCLUDED
