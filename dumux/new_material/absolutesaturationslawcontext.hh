/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \file
 *
 * \brief A default implementation of the context for the material law
 *        for absolute saturations.
 */
#ifndef ABSOLUTE_SATURATIONS_LAW_CONTEXT_HH
#define ABSOLUTE_SATURATIONS_LAW_CONTEXT_HH

namespace Dune
{
/*!
 * \ingroup material
 *
 * \brief A default implementation of the context for the material law
 *        for absolute saturations.
 */
template <class RawLawContextT>
class AbsoluteSaturationsLawContext : public RawLawContextT
{
    typedef RawLawContextT  RawLawContext;
public:
    typedef typename RawLawContext::Scalar   Scalar;

    AbsoluteSaturationsLawContext()
        : RawLawContext()
    {
        Swr_ = Snr_ = 0;
    }

    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar Swr() const
    { return Swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     */
    void setSwr(Scalar v)
    { Swr_ = v; }

    /*!
     * \brief Return the residual non-wetting saturation.
     */
    Scalar Snr() const
    { return Swr_; }

    /*!
     * \brief Set the residual non-wetting saturation.
     */
    void setSnr(Scalar v)
    { Snr_ = v; }

private:
    Scalar Swr_;
    Scalar Snr_;
};

}

#endif
