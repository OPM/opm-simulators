/*
  Copyright (C) 2009-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Opm::PiecewiseLinearTwoPhaseMaterialParams
 */
#ifndef OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <opm/parser/eclipse/Utility/SwofTable.hpp>
#include <opm/parser/eclipse/Utility/SgofTable.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for the van
 *        Genuchten constitutive relations.
 *
 * In this implementation setting either the \f$n\f$ or \f$m\f$ shape
 * parameter automatically calculates the other. I.e. they cannot be
 * set independently.
 */
template<class TraitsT>
class PiecewiseLinearTwoPhaseMaterialParams
{
    typedef typename TraitsT::Scalar Scalar;

public:
    typedef std::pair<Scalar, Scalar> SamplePoint;
    typedef std::vector<SamplePoint> SamplePoints;

    typedef TraitsT Traits;

    PiecewiseLinearTwoPhaseMaterialParams()
    {
#ifndef NDEBUG
        finalized_ = false;
#endif
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
#ifndef NDEBUG
        finalized_ = true;
#endif
    }

    /*!
     * \brief Read the relevant curves from the table specified by the
     *        SWOF keyword of an ECLIPSE input file
     */
    void readFromSwof(const Opm::SwofTable &swofTable)
    {
        const std::vector<double> &SwColumn = swofTable.getSwColumn();
        const std::vector<double> &krwColumn = swofTable.getKrwColumn();
        const std::vector<double> &krowColumn = swofTable.getKrowColumn();
        const std::vector<double> &pcowColumn = swofTable.getPcowColumn();
        int numRows = swofTable.numRows();
        for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
            pcwnSamples_[rowIdx].first = SwColumn[rowIdx];
            pcwnSamples_[rowIdx].second = - pcowColumn[rowIdx];

            krwSamples_[rowIdx].first = SwColumn[rowIdx];
            krwSamples_[rowIdx].second = krwColumn[rowIdx];

            krnSamples_[rowIdx].first = SwColumn[rowIdx];
            krnSamples_[rowIdx].second = krowColumn[rowIdx];
        }
    }

    /*!
     * \brief Read the relevant curves from the table specified by the
     *        SGOF keyword of an ECLIPSE input file
     */
    void readFromSgof(const Opm::SgofTable &sgofTable)
    {
        const std::vector<double> &SgColumn = sgofTable.getSgColumn();
        const std::vector<double> &krgColumn = sgofTable.getKrgColumn();
        const std::vector<double> &krogColumn = sgofTable.getKrogColumn();
        const std::vector<double> &pcogColumn = sgofTable.getPcogColumn();
        int numRows = sgofTable.numRows();
        for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
            pcwnSamples_[rowIdx].first = 1 - SgColumn[rowIdx];
            pcwnSamples_[rowIdx].second = pcogColumn[rowIdx];

            krwSamples_[rowIdx].first = 1 - SgColumn[rowIdx];
            krwSamples_[rowIdx].second = krogColumn[rowIdx];

            krnSamples_[rowIdx].first = 1 - SgColumn[rowIdx];
            krnSamples_[rowIdx].second = krgColumn[rowIdx];
        }
    }

    /*!
     * \brief Return the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const SamplePoints& pcnwSamples() const
    { assertFinalized_(); return pcwnSamples_; }

    /*!
     * \brief Set the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setPcnwSamples(const SamplePoints& samples)
    { pcwnSamples_ = samples; }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const SamplePoints& krwSamples() const
    { assertFinalized_(); return krwSamples_; }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setKrwSamples(const SamplePoints& samples)
    { krwSamples_ = samples; }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const SamplePoints& krnSamples() const
    { assertFinalized_(); return krnSamples_; }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    void setKrnSamples(const SamplePoints& samples)
    { krnSamples_ = samples; }

private:
#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

    SamplePoints pcwnSamples_;
    SamplePoints krwSamples_;
    SamplePoints krnSamples_;
};
} // namespace Opm

#endif
