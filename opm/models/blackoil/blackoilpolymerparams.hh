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
 *
 * \brief Contains the parameters required to extend the black-oil model by polymer.
 */
#ifndef EWOMS_BLACK_OIL_POLYMER_PARAMS_HH
#define EWOMS_BLACK_OIL_POLYMER_PARAMS_HH

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/IntervalTabulated2DFunction.hpp>

#include <map>
#include <vector>

namespace Opm {

//! \brief Struct holding the parameters for the BlackOilPolymerModule class.
template<class Scalar>
struct BlackOilPolymerParams {
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using TabulatedTwoDFunction = IntervalTabulated2DFunction<Scalar>;

    enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the PLYROCK and PLYADS of any region.
     */
    void setNumSatRegions(unsigned numRegions)
    {
        plyrockDeadPoreVolume_.resize(numRegions);
        plyrockResidualResistanceFactor_.resize(numRegions);
        plyrockRockDensityFactor_.resize(numRegions);
        plyrockAdsorbtionIndex_.resize(numRegions);
        plyrockMaxAdsorbtion_.resize(numRegions);
        plyadsAdsorbedPolymer_.resize(numRegions);
    }

    /*!
     * \brief Specify the number of mix regions.
     *
     * This must be called before setting the PLYMAC and PLMIXPAR of any region.
     */
    void setNumMixRegions(unsigned numRegions, bool enablePolymerMolarWeight)
    {
        plymaxMaxConcentration_.resize(numRegions);
        plymixparToddLongstaff_.resize(numRegions);

        if (enablePolymerMolarWeight) {
            plyvmhCoefficients_.resize(numRegions);
        }
    }

    /*!
     * \brief Specify the polymer rock properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    void setPlyrock(unsigned satRegionIdx,
                    const Scalar& plyrockDeadPoreVolume,
                    const Scalar& plyrockResidualResistanceFactor,
                    const Scalar& plyrockRockDensityFactor,
                    const Scalar& plyrockAdsorbtionIndex,
                    const Scalar& plyrockMaxAdsorbtion)
    {
        plyrockDeadPoreVolume_[satRegionIdx] = plyrockDeadPoreVolume;
        plyrockResidualResistanceFactor_[satRegionIdx] = plyrockResidualResistanceFactor;
        plyrockRockDensityFactor_[satRegionIdx] = plyrockRockDensityFactor;
        plyrockAdsorbtionIndex_[satRegionIdx] = plyrockAdsorbtionIndex;
        plyrockMaxAdsorbtion_[satRegionIdx] = plyrockMaxAdsorbtion;
    }

    // a struct containing the constants to calculate polymer viscosity
    // based on Mark-Houwink equation and Huggins equation, the constants are provided
    // by the keyword PLYVMH
    struct PlyvmhCoefficients {
        Scalar k_mh;
        Scalar a_mh;
        Scalar gamma;
        Scalar kappa;
    };

    struct SkprpolyTable {
        double refConcentration;
        TabulatedTwoDFunction table_func;
    };

    std::vector<Scalar> plyrockDeadPoreVolume_;
    std::vector<Scalar> plyrockResidualResistanceFactor_;
    std::vector<Scalar> plyrockRockDensityFactor_;
    std::vector<Scalar> plyrockAdsorbtionIndex_;
    std::vector<Scalar> plyrockMaxAdsorbtion_;
    std::vector<TabulatedFunction> plyadsAdsorbedPolymer_;
    std::vector<TabulatedFunction> plyviscViscosityMultiplierTable_;
    std::vector<Scalar> plymaxMaxConcentration_;
    std::vector<Scalar> plymixparToddLongstaff_;
    std::vector<std::vector<Scalar>> plyshlogShearEffectRefMultiplier_;
    std::vector<std::vector<Scalar>> plyshlogShearEffectRefLogVelocity_;
    std::vector<Scalar> shrate_;
    bool hasShrate_;
    bool hasPlyshlog_;

    std::vector<PlyvmhCoefficients> plyvmhCoefficients_;
    std::map<int, TabulatedTwoDFunction> plymwinjTables_;
    std::map<int, TabulatedTwoDFunction> skprwatTables_;

    std::map<int, SkprpolyTable> skprpolyTables_;
};

} // namespace Opm

#endif
