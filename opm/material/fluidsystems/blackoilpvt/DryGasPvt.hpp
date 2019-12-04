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
 * \copydoc Opm::DryGasPvt
 */
#ifndef OPM_DRY_GAS_PVT_HPP
#define OPM_DRY_GAS_PVT_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvdgTable.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#endif

#include <vector>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar>
class DryGasPvt
{
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;

    explicit DryGasPvt() = default;
    DryGasPvt(const std::vector<Scalar>& gasReferenceDensity,
              const std::vector<TabulatedOneDFunction>& inverseGasB,
              const std::vector<TabulatedOneDFunction>& gasMu,
              const std::vector<TabulatedOneDFunction>& inverseGasBMu)
        : gasReferenceDensity_(gasReferenceDensity)
        , inverseGasB_(inverseGasB)
        , gasMu_(gasMu)
        , inverseGasBMu_(inverseGasBMu)
    {
    }
#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for dry gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromDeck(const Deck& deck, const EclipseState& eclState)
    {
        const auto& pvdgTables = eclState.getTableManager().getPvdgTables();
        const auto& densityKeyword = deck.getKeyword("DENSITY");

        assert(pvdgTables.size() == densityKeyword.size());

        size_t numRegions = pvdgTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityKeyword.getRecord(regionIdx).getItem("OIL").getSIDouble(0);
            Scalar rhoRefG = densityKeyword.getRecord(regionIdx).getItem("GAS").getSIDouble(0);
            Scalar rhoRefW = densityKeyword.getRecord(regionIdx).getItem("WATER").getSIDouble(0);

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

            // determine the molar masses of the components
            Scalar p = 1.01325e5; // surface pressure, [Pa]
            Scalar T = 273.15 + 15.56; // surface temperature, [K]
            Scalar MO = 175e-3; // [kg/mol]
            Scalar MG = Opm::Constants<Scalar>::R*T*rhoRefG / p; // [kg/mol], consequence of the ideal gas law
            Scalar MW = 18.0e-3; // [kg/mol]
            // TODO (?): the molar mass of the components can possibly specified
            // explicitly in the deck.
            setMolarMasses(regionIdx, MO, MG, MW);

            const auto& pvdgTable = pvdgTables.getTable<PvdgTable>(regionIdx);

            // say 99.97% of all time: "premature optimization is the root of all
            // evil". Eclipse does this "optimization" for apparently no good reason!
            std::vector<Scalar> invB(pvdgTable.numRows());
            const auto& Bg = pvdgTable.getFormationFactorColumn();
            for (unsigned i = 0; i < Bg.size(); ++ i) {
                invB[i] = 1.0/Bg[i];
            }

            size_t numSamples = invB.size();
            inverseGasB_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), invB);
            gasMu_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), pvdgTable.getViscosityColumn());
        }

        initEnd();
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        gasReferenceDensity_.resize(numRegions);
        inverseGasB_.resize(numRegions);
        inverseGasBMu_.resize(numRegions);
        gasMu_.resize(numRegions);
    }


    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar /*rhoRefOil*/,
                               Scalar rhoRefGas,
                               Scalar /*rhoRefWater*/)
    {
        gasReferenceDensity_[regionIdx] = rhoRefGas;
    }

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setMolarMasses(unsigned /*regionIdx*/,
                        Scalar /*MOil*/,
                        Scalar /*MGas*/,
                        Scalar /*MWater*/)
    { }

    /*!
     * \brief Initialize the viscosity of the gas phase.
     *
     * This is a function of \f$(p_g)\f$...
     */
    void setGasViscosity(unsigned regionIdx, const TabulatedOneDFunction& mug)
    { gasMu_[regionIdx] = mug; }

    /*!
     * \brief Initialize the function for the formation volume factor of dry gas
     *
     * \param samplePoints A container of \f$(p_g, B_g)\f$ values
     */
    void setGasFormationVolumeFactor(unsigned regionIdx, const SamplingPoints& samplePoints)
    {
        SamplingPoints tmp(samplePoints);
        auto it = tmp.begin();
        const auto& endIt = tmp.end();
        for (; it != endIt; ++ it)
            std::get<1>(*it) = 1.0/std::get<1>(*it);

        inverseGasB_[regionIdx].setContainerOfTuples(tmp);
        assert(inverseGasB_[regionIdx].monotonic());
    }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = gasMu_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the gas
            // formation volume factor and the gas viscosity
            const auto& gasMu = gasMu_[regionIdx];
            const auto& invGasB = inverseGasB_[regionIdx];
            assert(gasMu.numSamples() == invGasB.numSamples());

            std::vector<Scalar> pressureValues(gasMu.numSamples());
            std::vector<Scalar> invGasBMuValues(gasMu.numSamples());
            for (unsigned pIdx = 0; pIdx < gasMu.numSamples(); ++pIdx) {
                pressureValues[pIdx] = invGasB.xAt(pIdx);
                invGasBMuValues[pIdx] = invGasB.valueAt(pIdx) * (1.0/gasMu.valueAt(pIdx));
            }

            inverseGasBMu_[regionIdx].setXYContainers(pressureValues, invGasBMuValues);
        }
    }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return gasReferenceDensity_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx OPM_UNUSED,
                        const Evaluation& temperature OPM_UNUSED,
                        const Evaluation& pressure OPM_UNUSED,
                        const Evaluation& Rv OPM_UNUSED) const
    {
        throw std::runtime_error("Requested the enthalpy of gas but the thermal option is not enabled");
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rv*/) const
    { return saturatedViscosity(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const
    {
        const Evaluation& invBg = inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& /*Rv*/) const
    { return saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& /*temperature*/,
                                                     const Evaluation& pressure) const
    { return inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rv*/) const
    { return 0.0; /* this is dry gas! */ }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/,
                                              const Evaluation& /*oilSaturation*/,
                                              const Evaluation& /*maxOilSaturation*/) const
    { return 0.0; /* this is dry gas! */ }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/) const
    { return 0.0; /* this is dry gas! */ }

    const std::vector<Scalar>& gasReferenceDensity() const
    { return gasReferenceDensity_; }

    const std::vector<TabulatedOneDFunction>& inverseGasB() const
    { return inverseGasB_; }

    const std::vector<TabulatedOneDFunction>& gasMu() const
    { return gasMu_; }

    const std::vector<TabulatedOneDFunction> inverseGasBMu() const
    { return inverseGasBMu_; }

    bool operator==(const DryGasPvt<Scalar>& data) const
    {
        return gasReferenceDensity_ == data.gasReferenceDensity_ &&
               inverseGasB_ == data.inverseGasB_ &&
               gasMu_ == data.gasMu_ &&
               inverseGasBMu_ == data.inverseGasBMu_;
    }

private:
    std::vector<Scalar> gasReferenceDensity_;
    std::vector<TabulatedOneDFunction> inverseGasB_;
    std::vector<TabulatedOneDFunction> gasMu_;
    std::vector<TabulatedOneDFunction> inverseGasBMu_;
};

} // namespace Opm

#endif
