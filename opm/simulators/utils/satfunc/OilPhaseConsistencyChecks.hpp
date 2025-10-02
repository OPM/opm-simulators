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

#ifndef OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED
#define OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>
#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>

#include <cstddef>
#include <string>

namespace Opm::Satfunc::PhaseChecks::Oil {

    /// Verify that critical oil saturation in gas/oil system is in valid range
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SOcr_GO : public PhaseCheckBase<Scalar>
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
            exportedCheckValues[0] = this->sogcr_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Non-negative critical oil saturation in G/O system" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "0 <= SOGCR < 1" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SOGCR";
        }

    private:
        /// Critical oil saturation in gas/oil system.
        Scalar sogcr_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that minimum oil saturation in gas/oil system is in valid range
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SOmin_GO : public PhaseCheckBase<Scalar>
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
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgu_;
            exportedCheckValues[2] = this->swl_ + this->sgu_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Non-negative minimum oil saturation in G/O system" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SWL + SGU <= 1" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGU";
            headers[2] = "SWL + SGU";
        }

    private:
        /// Minimum (connate) water saturation in gas/oil system.
        Scalar swl_{};

        /// Maximum gas saturation in gas/oil system.
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

    /// Verify that critical oil saturation in gas/oil system is strictly
    /// smaller than maximum oil saturation.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class MobileOil_GO_SGmin : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 4; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgl_;
            exportedCheckValues[2] = this->sogcr_;
            exportedCheckValues[3] = Scalar{1} - (this->swl_ + this->sgl_);
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile oil saturation in G/O system at minimum gas saturation" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SOGCR < 1 - SWL - SGL" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGL";
            headers[2] = "SOGCR";
            headers[3] = "1 - SWL - SGL";
        }

    private:
        /// Minimum water saturation.
        Scalar swl_{};

        /// Minimum gas saturation.
        Scalar sgl_{};

        /// Critical oil saturation in gas/oil system.
        Scalar sogcr_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that critical oil saturation in gas/oil system is strictly
    /// smaller than oil saturation at critical gas saturation.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class MobileOil_GO_SGcr : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 4; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgcr_;
            exportedCheckValues[2] = this->sogcr_;
            exportedCheckValues[3] = Scalar{1} - (this->swl_ + this->sgcr_);
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile oil saturation in G/O system at critical gas saturation" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SOGCR < 1 - SWL - SGCR" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGCR";
            headers[2] = "SOGCR";
            headers[3] = "1 - SWL - SGCR";
        }

    private:
        /// Minimum water saturation.
        Scalar swl_{};

        /// Critical gas saturation.
        Scalar sgcr_{};

        /// Critical oil saturation in gas/oil system
        Scalar sogcr_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    // -----------------------------------------------------------------------

    /// Verify that critical oil saturation in oil/water system is in valid range
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SOcr_OW : public PhaseCheckBase<Scalar>
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
            exportedCheckValues[0] = this->sowcr_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Non-negative critical oil saturation in O/W system" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "0 <= SOWCR < 1" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SOWCR";
        }

    private:
        /// Critical oil saturation in oil/water system.
        Scalar sowcr_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that minimum oil saturation in oil/water system is in valid range
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SOmin_OW : public PhaseCheckBase<Scalar>
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
            exportedCheckValues[1] = this->swu_;
            exportedCheckValues[2] = this->sgl_ + this->swu_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Non-negative minimum oil saturation in G/O system" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SGL + SWU <= 1" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SWU";
            headers[2] = "SGL + SWU";
        }

    private:
        /// Minimum gas saturation.  Typically zero.
        Scalar sgl_{};

        /// Minimum (connate) saturation.
        Scalar swu_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that critical oil saturation in oil/water system is strictly
    /// smaller than maximum oil saturation.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class MobileOil_OW_SWmin : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 4; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgl_;
            exportedCheckValues[2] = this->sowcr_;
            exportedCheckValues[3] = Scalar{1} - (this->swl_ + this->sgl_);
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile oil saturation in O/W system at minimum water saturation" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SOWCR < 1 - SWL - SGL" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGL";
            headers[2] = "SOWCR";
            headers[3] = "1 - SWL - SGL";
        }

    private:
        /// Minimum (connate) water saturation.
        Scalar swl_{};

        /// Minimum gas saturation.  Typically zero.
        Scalar sgl_{};

        /// Critical oil saturation in oil/water system.
        Scalar sowcr_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

    /// Verify that critical oil saturation in oil/water system is strictly
    /// smaller than oil saturation at critical water saturation.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class MobileOil_OW_SWcr : public PhaseCheckBase<Scalar>
    {
    public:
        /// Number of \c Scalar values involved in the check.
        std::size_t numExportedCheckValues() const override { return 4; };

        /// Get a linearised copy of the \c Scalar values involved in the check.
        ///
        /// \param[in,out] exportedCheckValues Pointer to contiguous
        ///    sequence of at least numExportedCheckValues() \c Scalars.
        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgl_;
            exportedCheckValues[1] = this->swcr_;
            exportedCheckValues[2] = this->sowcr_;
            exportedCheckValues[3] = Scalar{1} - this->swcr_ - this->sgl_;
        }

        /// Descriptive textual summary of this check.
        std::string description() const override
        {
            return { "Mobile oil saturation in O/W system at critical water saturation" };
        }

        /// Textual representation of the consistency condition.
        std::string condition() const override
        {
            return { "SOWCR < 1 - SWCR - SGL" };
        }

        /// Retrieve names of the exported check values.
        ///
        /// \param[in,out] headers Pointer to contiguous sequence of at
        ///    least numExportedCheckValues() strings.
        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SWCR";
            headers[2] = "SOWCR";
            headers[3] = "1 - SWCR - SGL";
        }

    private:
        /// Minimum gas saturation.  Typically zero.
        Scalar sgl_{};

        /// Critical water saturation.
        Scalar swcr_{};

        /// Critical oil saturation in oil/water system.
        Scalar sowcr_{};

        /// Run check against a set of saturation function end-points.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.
        void testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints) override;
    };

} // namespace Opm::Satfunc::PhaseChecks::Oil

#endif // OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED
