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
 * \copydoc Opm::OilPvtThermal
 */
#ifndef OPM_OIL_PVT_THERMAL_HPP
#define OPM_OIL_PVT_THERMAL_HPP

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
class OilPvtMultiplexer;

/*!
 * \brief This class implements temperature dependence of the PVT properties of oil
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class OilPvtThermal
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef OilPvtMultiplexer<Scalar, /*enableThermal=*/false> IsothermalPvt;

public:
    OilPvtThermal()
    {
        enableThermalDensity_ = false;
        enableThermalViscosity_ = false;
        enableInternalEnergy_ = false;
    }

    ~OilPvtThermal()
    { delete isothermalPvt_; }

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the oil PVT properties.
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

        enableThermalDensity_ = deck.hasKeyword("OILDENT");
        enableThermalViscosity_ = deck.hasKeyword("VISCREF");
        enableInternalEnergy_ = deck.hasKeyword("SPECHEAT");

        unsigned numRegions = isothermalPvt_->numRegions();
        setNumRegions(numRegions);

        // viscosity
        if (deck.hasKeyword("VISCREF")) {
            const auto& oilvisctTables = tables.getOilvisctTables();
            const auto& viscrefKeyword = deck.getKeyword("VISCREF");

            assert(oilvisctTables.size() == numRegions);
            assert(viscrefKeyword.size() == numRegions);

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& TCol = oilvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& muCol = oilvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
                oilvisctCurves_[regionIdx].setXYContainers(TCol, muCol);

                const auto& viscrefRecord = viscrefKeyword.getRecord(regionIdx);
                viscrefPress_[regionIdx] = viscrefRecord.getItem("REFERENCE_PRESSURE").getSIDouble(0);
                viscrefRs_[regionIdx] = viscrefRecord.getItem("REFERENCE_RS").getSIDouble(0);

                // temperature used to calculate the reference viscosity [K]. the
                // value does not really matter if the underlying PVT object really
                // is isothermal...
                Scalar Tref = 273.15 + 20;

                // compute the reference viscosity using the isothermal PVT object.
                viscRef_[regionIdx] =
                    isothermalPvt_->viscosity(regionIdx,
                                              Tref,
                                              viscrefPress_[regionIdx],
                                              viscrefRs_[regionIdx]);
            }
        }

        // temperature dependence of oil density
        if (enableThermalDensity_) {
            const auto& oildentKeyword = deck.getKeyword("OILDENT");

            assert(oildentKeyword.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& oildentRecord = oildentKeyword.getRecord(regionIdx);

                oildentRefTemp_[regionIdx] = oildentRecord.getItem("REFERENCE_TEMPERATURE").getSIDouble(0);
                oildentCT1_[regionIdx] = oildentRecord.getItem("EXPANSION_COEFF_LINEAR").getSIDouble(0);
                oildentCT2_[regionIdx] = oildentRecord.getItem("EXPANSION_COEFF_QUADRATIC").getSIDouble(0);
            }
        }

        if (deck.hasKeyword("SPECHEAT")) {
            // the specific internal energy of liquid oil. be aware that ecl only specifies the
            // heat capacity (via the SPECHEAT keyword) and we need to integrate it
            // ourselfs to get the internal energy
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& specheatTable = tables.getSpecheatTables()[regionIdx];
                const auto& temperatureColumn = specheatTable.getColumn("TEMPERATURE");
                const auto& cvOilColumn = specheatTable.getColumn("CV_OIL");

                std::vector<double> uSamples(temperatureColumn.size());

                Scalar u = temperatureColumn[0]*cvOilColumn[0];
                for (size_t i = 0;; ++i) {
                    uSamples[i] = u;

                    if (i >= temperatureColumn.size() - 1)
                        break;

                    // integrate to the heat capacity from the current sampling point to the next
                    // one. this leads to a quadratic polynomial.
                    Scalar c_v0 = cvOilColumn[i];
                    Scalar c_v1 = cvOilColumn[i + 1];
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
        oilvisctCurves_.resize(numRegions);
        viscrefPress_.resize(numRegions);
        viscrefRs_.resize(numRegions);
        viscRef_.resize(numRegions);
        internalEnergyCurves_.resize(numRegions);
    }

    /*!
     * \brief Finish initializing the thermal part of the oil phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns true iff the density of the oil phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

    /*!
     * \brief Returns true iff the viscosity of the oil phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return viscrefRs_.size(); }

    /*!
     * \brief Returns the specific internal energy [J/kg] of oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure OPM_UNUSED,
                              const Evaluation& Rs OPM_UNUSED) const
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
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    {
        const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure, Rs);
        if (!enableThermalViscosity())
            return isothermalMu;

        // compute the viscosity deviation due to temperature
        const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature);
        return muOilvisct/viscRef_[regionIdx]*isothermalMu;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        const auto& isothermalMu = isothermalPvt_->saturatedViscosity(regionIdx, temperature, pressure);
        if (!enableThermalViscosity())
            return isothermalMu;

        // compute the viscosity deviation due to temperature
        const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature);
        return muOilvisct/viscRef_[regionIdx]*isothermalMu;
    }


    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    {
        const auto& b =
            isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs);

        if (!enableThermalDensity())
            return b;

        // we use the same approach as for the for water here, but with the OPM-specific
        // OILDENT keyword.
        Scalar TRef = oildentRefTemp_[regionIdx];
        Scalar cT1 = oildentCT1_[regionIdx];
        Scalar cT2 = oildentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        return b/(1 + (cT1 + cT2*Y)*Y);
    }

    /*!
     * \brief Returns the formation volume factor [-] of gas-saturated oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    {
        const auto& b =
            isothermalPvt_->saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);

        if (!enableThermalDensity())
            return b;

        // we use the same approach as for the for water here, but with the OPM-specific
        // OILDENT keyword.
        Scalar TRef = oildentRefTemp_[regionIdx];
        Scalar cT1 = oildentCT1_[regionIdx];
        Scalar cT2 = oildentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        return b/(1 + (cT1 + cT2*Y)*Y);
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This method implements temperature dependence and requires the isothermal gas
     * dissolution factor for gas saturated oil and temperature as inputs. Currently it
     * is just a dummy method which passes through the isothermal gas dissolution factor.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const
    { return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This method implements temperature dependence and requires the isothermal gas
     * dissolution factor for gas saturated oil and temperature as inputs. Currently it
     * is just a dummy method which passes through the isothermal gas dissolution factor.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             Scalar maxOilSaturation) const
    { return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *
     * This method implements temperature dependence and requires isothermal satuation
     * pressure and temperature as inputs. Currently it is just a dummy method which
     * passes through the isothermal saturation pressure.
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { return isothermalPvt_->saturationPressure(regionIdx, temperature, pressure); }

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> oilvisctCurves_;
    std::vector<Scalar> viscrefPress_;
    std::vector<Scalar> viscrefRs_;
    std::vector<Scalar> viscRef_;

    // The PVT properties needed for temperature dependence of the density.
    std::vector<Scalar> oildentRefTemp_;
    std::vector<Scalar> oildentCT1_;
    std::vector<Scalar> oildentCT2_;

    // piecewise linear curve representing the internal energy of oil
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Opm

#endif
