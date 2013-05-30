// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::VariableLengthSpline_
 */
#ifndef OPM_VARIABLE_LENGTH_SPLINE_HH
#define OPM_VARIABLE_LENGTH_SPLINE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/btdmatrix.hh>

#include "splinecommon_.hh"

namespace Opm
{
/*!
 * \brief The common code for all 3rd order polynomial splines with
 *        where the number of sampling points only known at run-time.
 */
template<class ScalarT>
class VariableLengthSpline_
    : public SplineCommon_<ScalarT,
                           VariableLengthSpline_<ScalarT> >
{
    friend class SplineCommon_<ScalarT, VariableLengthSpline_<ScalarT> >;

    typedef ScalarT Scalar;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > Vector;
    typedef Dune::BTDMatrix<Dune::FieldMatrix<Scalar, 1, 1> > BTDMatrix;

public:
    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return xPos_.size(); }


    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    // Full splines                      //
    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using C-style arrays.
     *
     * This method uses separate array-like objects for the values of
     * the X and Y coordinates. In this context 'array-like' means
     * that an access to the members is provided via the []
     * operator. (e.g. C arrays, std::vector, std::array, etc.)  Each
     * array must be of size 'nSamples' at least. Also, the number of
     * sampling points must be larger than 1.
     */
    template <class ScalarArrayX, class ScalarArrayY>
    void setXYArrays(int nSamples,
                     const ScalarArrayX &x,
                     const ScalarArrayY &y,
                     Scalar m0, Scalar m1)
    {
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = x[i];
            yPos_[i] = y[i];
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using STL-compatible containers.
     *
     * This method uses separate STL-compatible containers for the
     * values of the sampling points' X and Y
     * coordinates. "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class ScalarContainerX, class ScalarContainerY>
    void setXYContainers(const ScalarContainerX &x,
                         const ScalarContainerY &y,
                         Scalar m0, Scalar m1)
    {
        assert(x.size() == y.size());
        assert(x.size() > 1);

        setNumSamples_(x.size());

        std::copy(x.begin(), x.end(), xPos_.begin());
        std::copy(y.begin(), y.end(), yPos_.begin());

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using a C-style array.
     *
     * This method uses a single array of sampling points, which are
     * seen as an array-like object which provides access to the X and
     * Y coordinates.  In this context 'array-like' means that an
     * access to the members is provided via the [] operator. (e.g. C
     * arrays, std::vector, std::array, etc.)  The array containing
     * the sampling points must be of size 'nSamples' at least. Also,
     * the number of sampling points must be larger than 1.
     */
    template <class PointArray>
    void setArrayOfPoints(int nSamples,
                          const PointArray &points,
                          Scalar m0,
                          Scalar m1)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = points[i][0];
            yPos_[i] = points[i][1];
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using a STL-compatible container of
     *        array-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be array-like objects storing the
     * X and Y coordinates.  "STL-compatible" means that the container
     * provides access to iterators using the begin(), end() methods
     * and also provides a size() method. Also, the number of entries
     * in the X and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfPoints(const XYContainer &points,
                              Scalar m0,
                              Scalar m1)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(points.size() > 1);

        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
            xPos_[i] = (*it)[0];
            yPos_[i] = (*it)[1];
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        // make a full spline
        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using a STL-compatible container of
     *        tuple-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be tuple-like objects storing the
     * X and Y coordinates.  "tuple-like" means that the objects
     * provide access to the x values via std::get<0>(obj) and to the
     * y value via std::get<1>(obj) (e.g. std::tuple or
     * std::pair). "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfTuples(const XYContainer &points,
                              Scalar m0,
                              Scalar m1)
    {
        // resize internal arrays
        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
            xPos_[i] = std::get<0>(*it);
            yPos_[i] = std::get<1>(*it);
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        // make a full spline
        makeFullSpline_(m0, m1);
    }

    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    // Natural splines                   //
    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    /*!
     * \brief Set the sampling points natural spline using C-style arrays.
     *
     * This method uses separate array-like objects for the values of
     * the X and Y coordinates. In this context 'array-like' means
     * that an access to the members is provided via the []
     * operator. (e.g. C arrays, std::vector, std::array, etc.)  Each
     * array must be of size 'nSamples' at least. Also, the number of
     * sampling points must be larger than 1.
     */
    template <class ScalarArrayX, class ScalarArrayY>
    void setXYArrays(int nSamples,
                     const ScalarArrayX &x,
                     const ScalarArrayY &y)
    {
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = x[i];
            yPos_[i] = y[i];
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points of a natural spline using
     *        STL-compatible containers.
     *
     * This method uses separate STL-compatible containers for the
     * values of the sampling points' X and Y
     * coordinates. "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class ScalarContainerX, class ScalarContainerY>
    void setXYContainers(const ScalarContainerX &x,
                         const ScalarContainerY &y)
    {
        assert(x.size() == y.size());
        assert(x.size() > 1);

        setNumSamples_(x.size());
        std::copy(x.begin(), x.end(), xPos_.begin());
        std::copy(y.begin(), y.end(), yPos_.begin());

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points of a natural spline using a
     *        C-style array.
     *
     * This method uses a single array of sampling points, which are
     * seen as an array-like object which provides access to the X and
     * Y coordinates.  In this context 'array-like' means that an
     * access to the members is provided via the [] operator. (e.g. C
     * arrays, std::vector, std::array, etc.)  The array containing
     * the sampling points must be of size 'nSamples' at least. Also,
     * the number of sampling points must be larger than 1.
     */
    template <class PointArray>
    void setArrayOfPoints(int nSamples,
                          const PointArray &points)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = points[i][0];
            yPos_[i] = points[i][1];
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points of a natural spline using a
     *        STL-compatible container of array-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be array-like objects storing the
     * X and Y coordinates.  "STL-compatible" means that the container
     * provides access to iterators using the begin(), end() methods
     * and also provides a size() method. Also, the number of entries
     * in the X and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfPoints(const XYContainer &points)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(points.size() > 1);

        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++ i, ++it) {
            xPos_[i] = (*it)[0];
            yPos_[i] = (*it)[1];
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        // make a natural spline
        makeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points of a natural spline using a
     *        STL-compatible container of tuple-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be tuple-like objects storing the
     * X and Y coordinates.  "tuple-like" means that the objects
     * provide access to the x values via std::get<0>(obj) and to the
     * y value via std::get<1>(obj) (e.g. std::tuple or
     * std::pair). "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfTuples(const XYContainer &points)
    {
        // resize internal arrays
        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
            xPos_[i] = std::get<0>(*it);
            yPos_[i] = std::get<1>(*it);
        }

        if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        // make a natural spline
        makeNaturalSpline_();
    }

protected:
    /*!
     * \brief Reverse order of the elements in the arrays which
     *        contain the sampling points.
     */
    void reverseSamplingPoints_()
    {
        // reverse the arrays
        int n = numSamples();
        for (int i = 0; i <= (n - 1)/2; ++i) {
            std::swap(xPos_[i], xPos_[n - i - 1]);
            std::swap(yPos_[i], yPos_[n - i - 1]);
        }
    }

    /*!
     * \brief Resizes the internal vectors to store the sample points.
     */
    void setNumSamples_(int nSamples)
    {
        xPos_.resize(nSamples);
        yPos_.resize(nSamples);
        m_.resize(nSamples);
    }

    /*!
     * \brief Create a natural spline from the already set sampling points.
     *
     * Also creates temporary matrix and right hand side vector.
     */
    void makeFullSpline_(Scalar m0, Scalar m1)
    {
        BTDMatrix M(numSamples());
        Vector d(numSamples());

        // create linear system of equations
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Create a natural spline from the already set sampling points.
     *
     * Also creates temporary matrix and right hand side vector.
     */
    void makeNaturalSpline_()
    {
        BTDMatrix M(numSamples());
        Vector d(numSamples());

        // create linear system of equations
        this->makeNaturalSystem_(M, d);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Returns the x coordinate of the i-th sampling point.
     */
    Scalar x_(int i) const
    { return xPos_[i]; }

    /*!
     * \brief Returns the y coordinate of the i-th sampling point.
     */
    Scalar y_(int i) const
    { return yPos_[i]; }

    /*!
     * \brief Returns the moment (i.e. second derivative) of the
     *        spline at the i-th sampling point.
     */
    Scalar moment_(int i) const
    { return m_[i]; }

    Vector xPos_;
    Vector yPos_;
    Vector m_;
};
}

#endif
