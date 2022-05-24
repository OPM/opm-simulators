/*
  Copyright SINTEF AS

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
#ifndef OPM_CUVECTOR_HEADER_INCLUDED
#define OPM_CUVECTOR_HEADER_INCLUDED
#include <vector>

namespace Opm::cuistl
{

/*! \brief Simple vector class on the GPU.
 *
 */
template <typename T>
class CuVector
{
public:
    CuVector(const int numberOfElements);
    CuVector(const T* dataOnHost, const int numberOfElements);
    virtual ~CuVector();

    const T* data() const;
    T* data();
private:
    T* dataOnDevice;
    const int numberOfElements;
};

} // namespace Opm::cuistl
#endif
