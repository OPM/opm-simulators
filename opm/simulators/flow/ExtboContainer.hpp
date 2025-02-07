// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_EXTBO_CONTAINER_HPP
#define OPM_EXTBO_CONTAINER_HPP

#include <vector>

namespace Opm {

namespace data { class Solution; }

template<class Scalar>
class ExtboContainer
{
    using ScalarBuffer = std::vector<Scalar>;

public:
    void allocate(const unsigned bufferSize);

    void assignMassFractions(const unsigned globalDofIdx,
                             const Scalar gas,
                             const Scalar oil,
                             const Scalar co2);

    void assignVolumes(const unsigned globalDofIdx,
                       const Scalar xVolume,
                       const Scalar yVolume);

    void assignZFraction(const unsigned globalDofIdx,
                         const Scalar zFraction);

    void outputRestart(data::Solution& sol);

    bool allocated() const
    { return allocated_; }

private:
    bool allocated_ = false;
    ScalarBuffer X_volume_;
    ScalarBuffer Y_volume_;
    ScalarBuffer Z_fraction_;
    ScalarBuffer mFracOil_;
    ScalarBuffer mFracGas_;
    ScalarBuffer mFracCo2_;
};

} // namespace Opm

#endif // OPM_EXTBO_CONTAINER_HPP
