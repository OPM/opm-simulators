/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file gdllawparams.hh
 *
 * \brief Parameter class for the material for the gas diffusion layer of
 *        PEM fuel cells.
 */
#ifndef DUMUX_GDL_LAW_PARAMS_HH
#define DUMUX_GDL_LAW_PARAMS_HH

namespace Dumux
{
/*!
 * \brief Parameter class for the material for the gas diffusion layer of
 *        PEM fuel cells.
 */
template<class ScalarT>
class GdlLawParams
{
public:
    typedef ScalarT Scalar;

    GdlLawParams()
    {}

    /*!
     * \brief Sets the arrays with the sampling points for the pc-Sw
     *        curve.
     */
    void setSamplePoints(const std::vector<Scalar> &SwValues,
                         const std::vector<Scalar> &pcValues)
    {
        assert(SwValues.size() == pcValues.size());
        SwValues_ = SwValues;
        pcValues_ = pcValues;
    };

    /*!
     * \brief Sets the arrays with the sampling points for the pc-Sw
     *        curve.
     */
    void setSamplePoints(int numSamples,
                         const Scalar *SwValues,
                         const Scalar *pcValues)
    {
        SwValues_.resize(numSamples);
        pcValues_.resize(numSamples);
        for (int i = 0; i < numSamples; ++i) {
            SwValues_[i] = SwValues[i];
            pcValues_[i] = pcValues[i];
        };
    };

    /*!
     * \brief Returns the number of sampling points
     */
    int numSamples() const
    { return SwValues_.size(); }

    /*!
     * \brief Returns the saturation of a sampling point
     */
    Scalar Sw(int sampleIdx) const
    { return SwValues_[sampleIdx]; }

    /*!
     * \brief Returns the capillary pressure of a sampling point
     */
    Scalar pC(int sampleIdx) const
    { return pcValues_[sampleIdx]; }

private:
    std::vector<Scalar> SwValues_;
    std::vector<Scalar> pcValues_;
};
} // namespace Dumux

#endif
