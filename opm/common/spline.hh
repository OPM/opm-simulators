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
 * \copydoc Opm::Spline
 */
#ifndef OPM_SPLINE_HH
#define OPM_SPLINE_HH

#include "fixedlengthspline_.hh"
#include "variablelengthspline_.hh"
#include "splinecommon_.hh"

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm
{
/*!
 * \ingroup Spline
 * \brief A 3rd order polynomial spline.

 * This class implements a spline \f$s(x)\f$ for which, given \f$n\f$ sampling
 * points \f$x_1, \dots, x_n\f$, the following conditions hold
 *\f{align*}{
   s(x_i)   & = y_i \quad \forall i \in \{1, \dots, n \}\\
   s'(x_1)  & = m_1 \\
   s'(x_n)  & = m_n
   \f}
*
* for any given boundary slopes \f$m_1\f$ and \f$m_n\f$. Alternatively, natural
* splines are supported which are defined by
*\f{align*}{
    s(x_i)     & = y_i \quad \forall i \in \{1, \dots, n \} \\
    s''(x_1)   & = 0 \\
    s''(x_n)   & = 0
\f}
 */
template<class Scalar, int numSamples = 2>
class Spline : public FixedLengthSpline_<Scalar, numSamples>
{
public:
    /*!
     * \brief Default constructor for a spline.
     *
     * To specfiy the acutal curve, use one of the set() methods.
     */
    Spline()
    { }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     */
    template <class ScalarArray>
    Spline(const ScalarArray &x,
           const ScalarArray &y)
    { this->setXYArrays(numSamples, x, y); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     */
    template <class PointArray>
    Spline(const PointArray &points)
    { this->setArrayOfPoints(numSamples, points); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarArray>
    Spline(const ScalarArray &x,
           const ScalarArray &y,
           Scalar m0,
           Scalar m1)
    { this->setXYArrays(numSamples, x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class PointArray>
    Spline(const PointArray &points,
           Scalar m0,
           Scalar m1)
    { this->setArrayOfPoints(numSamples, points, m0, m1); }
};

/*!
 * \brief Specialization of a spline with the number of sampling
 *        points only known at run time.
 *
 * This class implements a spline \f$s(x)\f$ for which, given \f$n\f$ sampling
 * points \f$x_1, \dots, x_n\f$, the following conditions hold
 *\f{align*}{
     s(x_i)   & = y_i \quad \forall i \in \{1, \dots, n \}\\
     s'(x_1)  & = m_1 \\
     s'(x_n)  & = m_n
   \f}
 *
 * for any given boundary slopes \f$m_1\f$ and \f$m_n\f$. Alternatively, natural
 * splines are supported which are defined by
 *\f{align*}{
    s(x_i)     & = y_i \quad \forall i \in \{1, \dots, n \} \\
    s''(x_1)   & = 0 \\
    s''(x_n)   & = 0
 \f}
*/
template<class Scalar>
class Spline<Scalar, /*numSamples=*/-1> : public VariableLengthSpline_<Scalar>
{
public:
    /*!
     * \brief Default constructor for a spline.
     *
     * To specfiy the acutal curve, use one of the set() methods.
     */
    Spline()
    { }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param nSamples The number of sampling points (must be > 2)
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     */
    template <class ScalarArrayX, class ScalarArrayY>
    Spline(int nSamples,
           const ScalarArrayX &x,
           const ScalarArrayY &y)
    { this->setXYArrays(nSamples, x, y); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param nSamples The number of sampling points (must be > 2)
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     */
    template <class PointArray>
    Spline(int nSamples,
           const PointArray &points)
    { this->setArrayOfPoints(nSamples, points); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points (must have a size() method)
     * \param y An array containing the \f$y\f$ values of the spline's sampling points (must have a size() method)
     */
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y)
    { this->setXYContainers(x, y); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points (must have a size() method)
     */
    template <class PointContainer>
    Spline(const PointContainer &points)
    { this->setContainerOfPoints(points); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarArray>
    Spline(int nSamples,
           const ScalarArray &x,
           const ScalarArray &y,
           Scalar m0,
           Scalar m1)
    { this->setXYArrays(nSamples, x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param points An array containing the \f$x\f$ and \f$x\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class PointArray>
    Spline(int nSamples,
           const PointArray &points,
           Scalar m0,
           Scalar m1)
    { this->setArrayOfPoints(nSamples, points, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points (must have a size() method)
     * \param y An array containing the \f$y\f$ values of the spline's sampling points (must have a size() method)
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarContainerX, class ScalarContainerY>
    Spline(const ScalarContainerX &x,
           const ScalarContainerY &y,
           Scalar m0,
           Scalar m1)
    { this->setXYContainers(x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points (must have a size() method)
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class PointContainer>
    Spline(const PointContainer &points,
           Scalar m0,
           Scalar m1)
    { this->setContainerOfPoints(points, m0, m1); }
};

/*!
 * \brief Do not allow splines with zero sampling points
 */
template<class Scalar>
class Spline<Scalar, /*numSamples=*/0>
// Splines with zero sampling points do not make sense!
{ private: Spline() { } };

/*!
 * \brief Do not allow splines with one sampling point
 */
template<class Scalar>
class Spline<Scalar, /*numSamples=*/1>
// Splines with one sampling point do not make sense!
{ private: Spline() { } };

/*!
 * \brief Spline for two sampling points.
 *
 * For this type of spline there is no natural spline.
 */
template<class Scalar>
class Spline<Scalar, 2> : public SplineCommon_<Scalar, Spline<Scalar, 2> >
{
    friend class  SplineCommon_<Scalar, Spline<Scalar, 2> >;
    typedef Dune::FieldVector<Scalar, 2> Vector;
    typedef Dune::FieldMatrix<Scalar, 2, 2> Matrix;

public:
    Spline()
    {}

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarArrayX, class ScalarArrayY>
    Spline(const ScalarArrayX &x,
           const ScalarArrayY &y,
           Scalar m0, Scalar m1)
    { setXYArrays(2, x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class PointArray>
    Spline(const PointArray &points,
           Scalar m0,
           Scalar m1)
    { this->setArrayOfPoints(2, points, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x0 The \f$x\f$ value of the first sampling point
     * \param x1 The \f$x\f$ value of the second sampling point
     * \param y0 The \f$y\f$ value of the first sampling point
     * \param y1 The \f$y\f$ value of the second sampling point
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    Spline(Scalar x0, Scalar x1,
           Scalar y0, Scalar y1,
           Scalar m0, Scalar m1)
    {
        set(x0, x1,
            y0, y1,
            m0, m1);
    }

    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return 2; }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param x0 The \f$x\f$ value of the first sampling point
     * \param x1 The \f$x\f$ value of the second sampling point
     * \param y0 The \f$y\f$ value of the first sampling point
     * \param y1 The \f$y\f$ value of the second sampling point
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    void set(Scalar x0, Scalar x1,
             Scalar y0, Scalar y1,
             Scalar m0, Scalar m1)
    {
        Matrix M(numSamples());
        Vector d;
        assignXY_(x0, x1, y0, y1);
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param x An array containing the \f$x\f$ values of the sampling points
     * \param y An array containing the \f$y\f$ values of the sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class ScalarContainer>
    void setXYArrays(int nSamples,
                     const ScalarContainer &x,
                     const ScalarContainer &y,
                     Scalar m0, Scalar m1)
    {
        assert(nSamples == 2);
        set(x[0], x[1], y[0], y[1], m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param x An array containing the \f$x\f$ values of the sampling points
     * \param y An array containing the \f$y\f$ values of the sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class ScalarContainerX, class ScalarContainerY>
    void setXYContainers(const ScalarContainerX &x,
                         const ScalarContainerY &y,
                         Scalar m0, Scalar m1)
    {
        assert(x.size() == y.size());
        assert(x.size() == 2);

        Matrix M(numSamples());
        Vector d;

        typename ScalarContainerX::const_iterator xIt0 = x.begin();
        typename ScalarContainerX::const_iterator xIt1 = xIt0;
        ++xIt1;
        typename ScalarContainerY::const_iterator yIt0 = y.begin();
        typename ScalarContainerY::const_iterator yIt1 = yIt0;
        ++yIt1;
        set(*xIt0, *xIt1, *yIt0, *yIt1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class PointArray>
    void setArrayOfPoints(int nSamples,
                          const PointArray &points,
                          Scalar m0,
                          Scalar m1)
    {
        assert(nSamples == 2);

        set(points[0][0],
            points[1][0],
            points[0][1],
            points[1][1],
            m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes from an
     * STL-like container of points.
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class PointContainer>
    void setContainerOfPoints(const PointContainer &points,
                              Scalar m0,
                              Scalar m1)
    {
        assert(points.size() == 2);

        Matrix M;
        Vector d;
        typename PointContainer::const_iterator it0 = points.begin();
        typename PointContainer::const_iterator it1 = it0;
        ++it1;

        set((*it0)[0],
            (*it0)[1],
            (*it1)[0],
            (*it1)[1],
            m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes from an
     *        STL-like container of tuples.
     *
     * \param tuples An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class TupleContainer>
    void setContainerOfTuples(const TupleContainer &tuples,
                              Scalar m0,
                              Scalar m1)
    {
        assert(tuples.size() == 2);

        typename TupleContainer::const_iterator it0 = tuples.begin();
        typename TupleContainer::const_iterator it1 = it0;
        ++it1;

        set(std::get<0>(*it0),
            std::get<1>(*it0),
            std::get<0>(*it1),
            std::get<1>(*it1),
            m0, m1);
    }

protected:
    void assignXY_(Scalar x0, Scalar x1,
                   Scalar y0, Scalar y1)
    {
        if (x0 > x1) {
            xPos_[0] = x1;
            xPos_[1] = x0;
            yPos_[0] = y1;
            yPos_[1] = y0;
        }
        else {
            xPos_[0] = x0;
            xPos_[1] = x1;
            yPos_[0] = y0;
            yPos_[1] = y1;
        }
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
