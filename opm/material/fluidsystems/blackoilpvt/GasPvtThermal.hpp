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
 * \copydoc Opm::GasPvtThermal
 */
#ifndef OPM_GAS_PVT_THERMAL_HPP
#define OPM_GAS_PVT_THERMAL_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {
template <class Scalar, bool enableThermal>
class GasPvtMultiplexer;

/*!
 * \brief This class implements temperature dependence of the PVT properties of gas
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class GasPvtThermal
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef GasPvtMultiplexer<Scalar, /*enableThermal=*/false> IsothermalPvt;

public:
    ~GasPvtThermal()
    { delete isothermalPvt_; }

#if HAVE_OPM_PARSER
    /*!
     * \brief Implement the temperature part of the gas PVT properties.
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

        enableThermalDensity_ = deck.hasKeyword("TREF");
        enableThermalViscosity_ = deck.hasKeyword("GASVISCT");

        unsigned numRegions = isothermalPvt_->numRegions();
        setNumRegions(numRegions);

        // viscosity
        if (enableThermalViscosity_) {
            const auto& gasvisctTables = tables.getGasvisctTables();
            int gasCompIdx = deck.getKeyword("GCOMPIDX").getRecord(0).getItem("GAS_COMPONENT_INDEX").get< int >(0) - 1;
            std::string gasvisctColumnName = "Viscosity"+std::to_string(static_cast<long long>(gasCompIdx));

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& T = gasvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& mu = gasvisctTables[regionIdx].getColumn(gasvisctColumnName).vectorCopy();
                gasvisctCurves_[regionIdx].setXYContainers(T, mu);
            }
        }

        // quantities required for density. note that we just always use the values
        // for the first EOS. (since EOS != PVT region.)
        refTemp_ = 0.0;
        if (enableThermalDensity_) {
            refTemp_ = deck.getKeyword("TREF").getRecord(0).getItem("TEMPERATURE").getSIDouble(0);
        }
    }
#endif // HAVE_OPM_PARSER

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions)
    { gasvisctCurves_.resize(numRegions); }

    /*!
     * \brief Finish initializing the thermal part of the gas phase PVT properties.
     */
    void initEnd()
    { }

    size_t numRegions() const
    { return gasvisctCurves_.size(); }

    /*!
     * \brief Returns true iff the density of the gas phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

    /*!
     * \brief Returns true iff the viscosity of the gas phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rv) const
    {
        if (!enableThermalViscosity())
            return isothermalPvt_->viscosity(regionIdx, temperature, pressure, Rv);

        // compute the viscosity deviation due to temperature
        const auto &muGasvisct = gasvisctCurves_[regionIdx].eval(temperature);
        return muGasvisct;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the oil-saturated gas phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        if (!enableThermalViscosity())
            return isothermalPvt_->saturatedViscosity(regionIdx, temperature, pressure);

        // compute the viscosity deviation due to temperature
        const auto &muGasvisct = gasvisctCurves_[regionIdx].eval(temperature);
        return muGasvisct;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rv) const
    {
        const auto& b =
            isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv);
        if (!enableThermalDensity())
            return b;

        // the Eclipse TD/RM do not explicitly specify the relation of the gas
        // density and the temperature, but equation (69.49) (for Eclipse 2011.1)
        // implies that the temperature dependence of the gas phase is rho(T, p) =
        // rho(tref_, p)*T/T_ref ...
        return b*temperature/refTemp_;
    }

    /*!
     * \brief Returns the formation volume factor [-] of oil-saturated gas.
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

        // the Eclipse TD/RM do not explicitly specify the relation of the gas
        // density and the temperature, but equation (69.49) (for Eclipse 2011.1)
        // implies that the temperature dependence of the gas phase is rho(T, p) =
        // rho(tref_, p)*T/T_ref ...
        return b*temperature/refTemp_;
    }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This method implements temperature dependence and requires the gas pressure,
     * temperature and the oil saturation as inputs. Currently it is just a dummy method
     * which passes through the isothermal oil vaporization factor.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure) const
    { return isothermalPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This method implements temperature dependence and requires the gas pressure,
     * temperature and the oil saturation as inputs. Currently it is just a dummy method
     * which passes through the isothermal oil vaporization factor.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure,
                                              const Evaluation& oilSaturation,
                                              Scalar maxOilSaturation) const
    { return isothermalPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation); }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
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
    IsothermalPvt *isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> gasvisctCurves_;

    // The PVT properties needed for temperature dependence of the density. This is
    // specified as one value per EOS in the manual, but we unconditionally use the
    // expansion coefficient of the first EOS...
    Scalar refTemp_;

    bool enableThermalDensity_;
    bool enableThermalViscosity_;
};

} // namespace Opm

#endif
