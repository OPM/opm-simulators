/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_POINT2D_HEADER_INCLUDED
#define OPM_POINT2D_HEADER_INCLUDED

namespace Opm {

    namespace detail {
        class Point2D
        {
            public:
                Point2D(const double xi, const double yi)
                    : x_(xi),
                      y_(yi) {

                      }

                Point2D()
                {
                }

                const double getX() const
                {
                    return x_;
                }

                const double getY() const
                {
                    return y_;
                }

                void setX(const double x)
                {
                    x_ = x;
                }

                void setY(const double y){
                    y_ = y;
                }

                /// Finding the intersection point of a line segment and a line.
                /// return true, if found.
                static bool findIntersection(Point2D line_segment1[2], Point2D line2[2], Point2D& intersection_point)
                {

                    const double x1 = line_segment1[0].getX();
                    const double y1 = line_segment1[0].getY();
                    const double x2 = line_segment1[1].getX();
                    const double y2 = line_segment1[1].getY();

                    const double x3 = line2[0].getX();
                    const double y3 = line2[0].getY();
                    const double x4 = line2[1].getX();
                    const double y4 = line2[1].getY();

                    const double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

                    if (d == 0.) {
                        return false;
                    }

                    const double x = ((x3 - x4) * (x1 * y2 - y1 * x2) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d;
                    const double y = ((y3 - y4) * (x1 * y2 - y1 * x2) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d;

                    if (x >= std::min(x1,x2) && x <= std::max(x1,x2)) {
                        intersection_point.setX(x);
                        intersection_point.setY(y);
                        return true;
                    } else {
                       return false;
                    }
                }

            private:
                double x_;
                double y_;

        }; // class Point2D

    } // namespace detail

} // namespace Opm

#endif // OPM_POINT2D_HEADER_INCLUDED

