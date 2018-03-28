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
 * \copydoc Opm::WaterPvtThermal
 */
#ifndef OPM_WATER_PVT_THERMAL_HPP
#define OPM_WATER_PVT_THERMAL_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {
template <class Scalar, bool enableThermal>
class WaterPvtMultiplexer;

/*!
 * \brief This class implements temperature dependence of the PVT properties of water
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class WaterPvtThermal
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef WaterPvtMultiplexer<Scalar, /*enableThermal=*/false> IsothermalPvt;

public:
    WaterPvtThermal()
    {
        enableThermalDensity_ = false;
        enableThermalViscosity_ = false;
        enableInternalEnergy_ = false;
    }

    ~WaterPvtThermal()
    { delete isothermalPvt_; }

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the water PVT properties.
     */
    void initFromDeck(const Deck& deck,
                      const EclipseState& eclState)
    {
        //////
        // initialize the isothermal part
        //////
        isothermalPvt_ = new IsothermalPvt;
        isothermalPvt_->initFromDeck(deck, eclState);

        //////
        // initialize the thermal part
        //////
        const auto& tables = eclState.getTableManager();

        enableThermalDensity_ = deck.hasKeyword("WATDENT");
        enableThermalViscosity_ = deck.hasKeyword("VISCREF");
        enableInternalEnergy_ = deck.hasKeyword("SPECHEAT");

        unsigned numRegions = isothermalPvt_->numRegions();
        setNumRegions(numRegions);

        if (enableThermalDensity_) {
            const auto& watdentKeyword = deck.getKeyword("WATDENT");

            assert(watdentKeyword.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& watdentRecord = watdentKeyword.getRecord(regionIdx);

                watdentRefTemp_[regionIdx] = watdentRecord.getItem("REFERENCE_TEMPERATURE").getSIDouble(0);
                watdentCT1_[regionIdx] = watdentRecord.getItem("EXPANSION_COEFF_LINEAR").getSIDouble(0);
                watdentCT2_[regionIdx] = watdentRecord.getItem("EXPANSION_COEFF_QUADRATIC").getSIDouble(0);
            }
        }

        if (enableThermalViscosity_) {
            const auto& viscrefKeyword = deck.getKeyword("VISCREF");

            const auto& watvisctTables = tables.getWatvisctTables();

            assert(watvisctTables.size() == numRegions);
            assert(viscrefKeyword.size() == numRegions);

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                const auto& T = watvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& mu = watvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
                watvisctCurves_[regionIdx].setXYContainers(T, mu);

                const auto& viscrefRecord = viscrefKeyword.getRecord(regionIdx);
                viscrefPress_[regionIdx] = viscrefRecord.getItem("REFERENCE_PRESSURE").getSIDouble(0);
            }
        }

        if (deck.hasKeyword("SPECHEAT")) {
            // the specific internal energy of liquid water. be aware that ecl only specifies the heat capacity
            // (via the SPECHEAT keyword) and we need to integrate it ourselfs to get the
            // internal energy
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& specHeatTable = tables.getSpecheatTables()[regionIdx];
                const auto& temperatureColumn = specHeatTable.getColumn("TEMPERATURE");
                const auto& cvWaterColumn = specHeatTable.getColumn("CV_WATER");

                std::vector<double> uSamples(temperatureColumn.size());

                Scalar u = temperatureColumn[0]*cvWaterColumn[0];
                for (size_t i = 0;; ++i) {
                    uSamples[i] = u;

                    if (i >= temperatureColumn.size() - 1)
                        break;

                    // integrate to the heat capacity from the current sampling point to the next
                    // one. this leads to a quadratic polynomial.
                    Scalar c_v0 = cvWaterColumn[i];
                    Scalar c_v1 = cvWaterColumn[i + 1];
                    Scalar T0 = temperatureColumn[i];
                    Scalar T1 = temperatureColumn[i + 1];
                    u += 0.5*(c_v0 + c_v1)*(T1 - T0);
                }

                internalEnergyCurves_[regionIdx].setXYContainers(temperatureColumn.vectorCopy(), uSamples);
            }
        }
    }
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions)
    {
        pvtwRefPress_.resize(numRegions);
        pvtwRefB_.resize(numRegions);
        pvtwCompressibility_.resize(numRegions);
        pvtwViscosity_.resize(numRegions);
        pvtwViscosibility_.resize(numRegions);
        viscrefPress_.resize(numRegions);
        watvisctCurves_.resize(numRegions);
        watdentRefTemp_.resize(numRegions);
        watdentCT1_.resize(numRegions);
        watdentCT2_.resize(numRegions);
        internalEnergyCurves_.resize(numRegions);
    }

    /*!
     * \brief Finish initializing the thermal part of the water phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns true iff the density of the water phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

    /*!
     * \brief Returns true iff the viscosity of the water phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return pvtwRefPress_.size(); }

    /*!
     * \brief Returns the specific internal energy [J/kg] of water given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure OPM_UNUSED) const
    {
        if (!enableInternalEnergy_)
            throw std::runtime_error("Requested the internal energy of oil but it is disabled");

        // compute the specific internal energy for the specified tempature. We use linear
        // interpolation here despite the fact that the underlying heat capacities are
        // piecewise linear (which leads to a quadratic function)
        return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure) const
    {
        const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure);
        if (!enableThermalViscosity())
            return isothermalMu;

        Scalar x = -pvtwViscosibility_[regionIdx]*(viscrefPress_[regionIdx] - pvtwRefPress_[regionIdx]);
        Scalar muRef = pvtwViscosity_[regionIdx]/(1.0 + x + 0.5*x*x);

        // compute the viscosity deviation due to temperature
        const auto& muWatvisct = watvisctCurves_[regionIdx].eval(temperature);
        return isothermalMu * muWatvisct/muRef;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure) const
    {
        if (!enableThermalDensity())
            return isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure);

        Scalar BwRef = pvtwRefB_[regionIdx];
        Scalar TRef = watdentRefTemp_[regionIdx];
        const Evaluation& X = pvtwCompressibility_[regionIdx]*(pressure - pvtwRefPress_[regionIdx]);
        Scalar cT1 = watdentCT1_[regionIdx];
        Scalar cT2 = watdentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        // this is inconsistent with the density calculation of water in the isothermal
        // case (it misses the quadratic pressure term), but it is the equation given in
        // the documentation.
        return 1.0/(((1 - X)*(1 + cT1*Y + cT2*Y*Y))*BwRef);
    }

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence. We need to store one
    // value per PVT region.
    std::vector<Scalar> viscrefPress_;

    std::vector<Scalar> watdentRefTemp_;
    std::vector<Scalar> watdentCT1_;
    std::vector<Scalar> watdentCT2_;

    std::vector<Scalar> pvtwRefPress_;
    std::vector<Scalar> pvtwRefB_;
    std::vector<Scalar> pvtwCompressibility_;
    std::vector<Scalar> pvtwViscosity_;
    std::vector<Scalar> pvtwViscosibility_;

    std::vector<TabulatedOneDFunction> watvisctCurves_;

    // piecewise linear curve representing the internal energy of water
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Opm

#endif
