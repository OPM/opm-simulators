// $Id$
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
 * \file
 *
 * \brief The base class for spatial parameters of problems using the
 *        box method.
 */
#ifndef DUNE_BOX_SPATIAL_PARAMETERS_HH
#define DUNE_BOX_SPATIAL_PARAMETERS_HH

namespace Dune
{

/**
 * \brief The base class for spatial parameters of problems using the
 *        box method.
 */
template<class TypeTag>
class BoxSpatialParameters 
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;

    enum {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

public:
    BoxSpatialParameters(const GridView &gv)
    { }

    ~BoxSpatialParameters()
    {}

    /*!
     * \brief Averages the intrinsic permeability.
     */
    const void meanIntrinsicPermeability(Tensor &result,
                                         Scalar K1,
                                         Scalar K2) const
    {
        const Scalar K = harmonicMean(K1, K2);
        for (int i = 0; i < dimWorld; ++i) {
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = 0;
            result[i][i] = K;
        }
    }

    /*!
     * \brief Averages the intrinsic permeability.
     */
    const void meanIntrinsicPermeability(Tensor &result,
                                         const Tensor &K1,
                                         const Tensor &K2) const
    {
        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = harmonicMean(K1[i][j], K2[i][j]);
    }
};

}

#endif
