// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

#include <opm/core/io/eclipse/EclipseGridParser.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for the van
       Genuchten constitutive relations.
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
     * \brief Read the relevant curves from an ECLIPSE input file
     *
     * The relevant keywords for this are "SWOF" and "SGOF".
     */
    template <class FieldType>
    void readFromEclipse(const FieldType& table)
    {
        int numSamples = table[0].size();
        pcwnSamples_.resize(numSamples);
        krwSamples_.resize(numSamples);
        krnSamples_.resize(numSamples);

        for (int sampleIdx = 0; sampleIdx < numSamples; ++sampleIdx) {
            Scalar Sw = table[0][sampleIdx];
            Scalar krw = table[1][sampleIdx];
            Scalar krn = table[2][sampleIdx];
            Scalar pcnw = table[3][sampleIdx];

            pcwnSamples_[sampleIdx].first = Sw;
            pcwnSamples_[sampleIdx].second = pcnw;

            krwSamples_[sampleIdx].first = Sw;
            krwSamples_[sampleIdx].second = krw;

            krnSamples_[sampleIdx].first = Sw;
            krnSamples_[sampleIdx].second = krn;
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
