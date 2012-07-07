// $Id: linearmaterialparams.hh 3777 2010-06-24 06:46:46Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \file linearmaterialparams.hh Reference implementation of params for
 *                           the linear material material.
 */
#ifndef LINEAR_MATERIAL_PARAMS_HH
#define LINEAR_MATERIAL_PARAMS_HH

namespace Dumux
{
/*!
 * \brief Reference implementation of params for the linear material
 *        material.
 */
template<class ScalarT>
class LinearMaterialParams
{
public:
    typedef ScalarT Scalar;

    LinearMaterialParams()
    {}

    LinearMaterialParams(Scalar entryPC, Scalar maxPC)
    {
        setEntryPC(entryPC);
        setMaxPC(maxPC);
    };

    /*!
     * \brief Return the threshold saturation at which the relative
     *        permeability starts to get regularized.
     *
     * This is simply 10%
     */
    Scalar Sreg() const
    { return 0.10; }

    /*!
     * \brief Return the entry pressure for the linear material law.
     *
     * The entry pressure is reached at \f$S_w = 1\f$
     */
    Scalar entryPC() const
    { return entryPC_; }

    /*!
     * \brief Set the entry pressure for the linear material law.
     *
     * The entry pressure is reached at \f$S_w = 1\f$
     */
    void setEntryPC(Scalar v)
    { entryPC_ = v; }

    /*!
     * \brief Return the maximum capillary pressure for the linear material law.
     *
     * The maximum capillary pressure is reached at \f$S_w = 0\f$
     */
    Scalar maxPC() const
    { return maxPC_; }

    /*!
     * \brief Set the maximum capillary pressure for the linear material law.
     *
     * The maximum capillary pressure is reached at \f$S_w = 0\f$
     */
    void setMaxPC(Scalar v)
    { maxPC_ = v; }

    /*!
     * \brief Return the threshold saturation respective phase below
     *        which the relative permeability gets regularized.
     *
     * This is just 5%. If you need a different value, write your own
     * parameter class.
     */
    Scalar krLowS() const
    { return 0.05; }

    /*!
     * \brief Return the threshold saturation of the respective phase
     *        above which the relative permeability gets regularized.
     *
     * This is just 95%. If you need a different value, write your own
     * parameter class.
     */
    Scalar krHighS() const
    { return 0.95; }

private:
    Scalar entryPC_;
    Scalar maxPC_;
};
} // namespace Dumux

#endif
