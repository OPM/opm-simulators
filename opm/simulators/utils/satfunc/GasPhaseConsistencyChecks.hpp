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

#ifndef GAS_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED
#define GAS_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>
#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

#include <cstddef>
#include <string>

namespace Opm::Satfunc::PhaseChecks::Gas {

    /// Verify that minimum gas saturation is in valid range
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SGmin : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 1; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgl_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Non-negative minimum gas saturation" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "0 <= SGL < 1" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
        }

    private:
        /// Minimum gas saturation.
        Scalar sgl_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that maximum gas saturation is in valid range.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SGmax : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 1; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgu_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Positive maximum gas saturation must not exceed one" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "0 < SGU <= 1" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGU";
        }

    private:
        /// Maximum gas saturation.
        Scalar sgu_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that critical gas saturation is in valid range.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SGcr : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 3; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgl_;
            exportedCheckValues[1] = this->sgcr_;
            exportedCheckValues[2] = this->sgu_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile gas saturation" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SGL <= SGCR < SGU" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SGCR";
            headers[2] = "SGU";
        }

    private:
        /// Minimum gas saturation.
        Scalar sgl_{};

        /// Critical gas saturation.
        Scalar sgcr_{};

        /// Maximum gas saturation.
        Scalar sgu_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

} // namespace Opm::Satfunc::PhaseChecks::Gas

#endif // GAS_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED
