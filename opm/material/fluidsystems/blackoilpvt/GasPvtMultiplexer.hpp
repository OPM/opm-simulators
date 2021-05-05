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
 * \copydoc Opm::GasPvtMultiplexer
 */
#ifndef OPM_GAS_PVT_MULTIPLEXER_HPP
#define OPM_GAS_PVT_MULTIPLEXER_HPP

#include "DryGasPvt.hpp"
#include "WetGasPvt.hpp"
#include "GasPvtThermal.hpp"
#include "Co2GasPvt.hpp"

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {
#define OPM_GAS_PVT_MULTIPLEXER_CALL(codeToCall)                          \
    switch (gasPvtApproach_) {                                            \
    case GasPvtApproach::DryGasPvt: {                                     \
        auto& pvtImpl = getRealPvt<GasPvtApproach::DryGasPvt>();          \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::WetGasPvt: {                                     \
        auto& pvtImpl = getRealPvt<GasPvtApproach::WetGasPvt>();          \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::ThermalGasPvt: {                                 \
        auto& pvtImpl = getRealPvt<GasPvtApproach::ThermalGasPvt>();      \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::Co2GasPvt: {                                     \
        auto& pvtImpl = getRealPvt<GasPvtApproach::Co2GasPvt>();          \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::NoGasPvt:                                        \
        throw std::logic_error("Not implemented: Gas PVT of this deck!"); \
    } \

enum class GasPvtApproach {
    NoGasPvt,
    DryGasPvt,
    WetGasPvt,
    ThermalGasPvt,
    Co2GasPvt
};

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas
 *        phase in the black-oil model.
 *
 * This is a multiplexer class which forwards all calls to the real implementation.
 *
 * Note that, since the main application for this class is the black oil fluid system,
 * the API exposed by this class is pretty specific to the assumptions made by the black
 * oil model.
 */
template <class Scalar, bool enableThermal = true>
class GasPvtMultiplexer
{
public:
    typedef Opm::GasPvtThermal<Scalar> GasPvtThermal;

    GasPvtMultiplexer()
    {
        gasPvtApproach_ = GasPvtApproach::NoGasPvt;
        realGasPvt_ = nullptr;
    }

    GasPvtMultiplexer(GasPvtApproach approach, void* realGasPvt)
        : gasPvtApproach_(approach)
        , realGasPvt_(realGasPvt)
    { }

    GasPvtMultiplexer(const GasPvtMultiplexer<Scalar,enableThermal>& data)
    {
        *this = data;
    }

    ~GasPvtMultiplexer()
    {
        switch (gasPvtApproach_) {
        case GasPvtApproach::DryGasPvt: {
            delete &getRealPvt<GasPvtApproach::DryGasPvt>();
            break;
        }
        case GasPvtApproach::WetGasPvt: {
            delete &getRealPvt<GasPvtApproach::WetGasPvt>();
            break;
        }
        case GasPvtApproach::ThermalGasPvt: {
            delete &getRealPvt<GasPvtApproach::ThermalGasPvt>();
            break;
        }
        case GasPvtApproach::Co2GasPvt: {
            delete &getRealPvt<GasPvtApproach::Co2GasPvt>();
            break;
        }
        case GasPvtApproach::NoGasPvt:
            break;
        }
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule)
    {
        if (!eclState.runspec().phases().active(Phase::GAS))
            return;
        if (eclState.runspec().co2Storage())
            setApproach(GasPvtApproach::Co2GasPvt);
        else if (enableThermal && eclState.getSimulationConfig().isThermal())
            setApproach(GasPvtApproach::ThermalGasPvt);
        else if (!eclState.getTableManager().getPvtgTables().empty())
            setApproach(GasPvtApproach::WetGasPvt);
        else if (eclState.getTableManager().hasTables("PVDG"))
            setApproach(GasPvtApproach::DryGasPvt);

        OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule));
    }
#endif // HAVE_ECL_INPUT

    void setApproach(GasPvtApproach gasPvtAppr)
    {
        switch (gasPvtAppr) {
        case GasPvtApproach::DryGasPvt:
            realGasPvt_ = new Opm::DryGasPvt<Scalar>;
            break;

        case GasPvtApproach::WetGasPvt:
            realGasPvt_ = new Opm::WetGasPvt<Scalar>;
            break;

        case GasPvtApproach::ThermalGasPvt:
            realGasPvt_ = new Opm::GasPvtThermal<Scalar>;
            break;

        case GasPvtApproach::Co2GasPvt:
            realGasPvt_ = new Opm::Co2GasPvt<Scalar>;
            break;

        case GasPvtApproach::NoGasPvt:
            throw std::logic_error("Not implemented: Gas PVT of this deck!");
        }

        gasPvtApproach_ = gasPvtAppr;
    }

    void initEnd()
    { OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions()); return 1; }

    /*!
     * \brief Return the reference density which are considered by this PVT-object.
     */
    const Scalar gasReferenceDensity(unsigned regionIdx)
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.gasReferenceDensity(regionIdx)); return 2.; }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rv) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rv)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rv) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rv)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation = Scalar>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rv) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of oil saturated gas.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedOilVaporizationFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of oil saturated gas.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure,
                                              const Evaluation& oilSaturation,
                                              const Evaluation& maxOilSaturation) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedOilVaporizationFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation)); return 0; }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation = Scalar>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rv) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rv)); return 0; }

    /*!
     * \copydoc BaseFluidSystem::diffusionCoefficient
     */
    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& temperature,
                                    const Evaluation& pressure,
                                    unsigned compIdx) const
    {
      OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.diffusionCoefficient(temperature, pressure, compIdx)); return 0;
    }

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    GasPvtApproach gasPvtApproach() const
    { return gasPvtApproach_; }

    // get the parameter object for the dry gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::DryGasPvt, Opm::DryGasPvt<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Opm::DryGasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::DryGasPvt, const Opm::DryGasPvt<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Opm::DryGasPvt<Scalar>* >(realGasPvt_);
    }

    // get the parameter object for the wet gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::WetGasPvt, Opm::WetGasPvt<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Opm::WetGasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::WetGasPvt, const Opm::WetGasPvt<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Opm::WetGasPvt<Scalar>* >(realGasPvt_);
    }

    // get the parameter object for the thermal gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::ThermalGasPvt, Opm::GasPvtThermal<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Opm::GasPvtThermal<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::ThermalGasPvt, const Opm::GasPvtThermal<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Opm::GasPvtThermal<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::Co2GasPvt, Opm::Co2GasPvt<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Opm::Co2GasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == GasPvtApproach::Co2GasPvt, const Opm::Co2GasPvt<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Opm::Co2GasPvt<Scalar>* >(realGasPvt_);
    }

    const void* realGasPvt() const { return realGasPvt_; }

    bool operator==(const GasPvtMultiplexer<Scalar,enableThermal>& data) const
    {
        if (this->gasPvtApproach() != data.gasPvtApproach())
            return false;

        switch (gasPvtApproach_) {
        case GasPvtApproach::DryGasPvt:
            return *static_cast<const Opm::DryGasPvt<Scalar>*>(realGasPvt_) ==
                   *static_cast<const Opm::DryGasPvt<Scalar>*>(data.realGasPvt_);
        case GasPvtApproach::WetGasPvt:
            return *static_cast<const Opm::WetGasPvt<Scalar>*>(realGasPvt_) ==
                   *static_cast<const Opm::WetGasPvt<Scalar>*>(data.realGasPvt_);
        case GasPvtApproach::ThermalGasPvt:
            return *static_cast<const Opm::GasPvtThermal<Scalar>*>(realGasPvt_) ==
                   *static_cast<const Opm::GasPvtThermal<Scalar>*>(data.realGasPvt_);
        case GasPvtApproach::Co2GasPvt:
            return *static_cast<const Opm::Co2GasPvt<Scalar>*>(realGasPvt_) ==
                    *static_cast<const Opm::Co2GasPvt<Scalar>*>(data.realGasPvt_);
        default:
            return true;
        }
    }

    GasPvtMultiplexer<Scalar,enableThermal>& operator=(const GasPvtMultiplexer<Scalar,enableThermal>& data)
    {
        gasPvtApproach_ = data.gasPvtApproach_;
        switch (gasPvtApproach_) {
        case GasPvtApproach::DryGasPvt:
            realGasPvt_ = new Opm::DryGasPvt<Scalar>(*static_cast<const Opm::DryGasPvt<Scalar>*>(data.realGasPvt_));
            break;
        case GasPvtApproach::WetGasPvt:
            realGasPvt_ = new Opm::WetGasPvt<Scalar>(*static_cast<const Opm::WetGasPvt<Scalar>*>(data.realGasPvt_));
            break;
        case GasPvtApproach::ThermalGasPvt:
            realGasPvt_ = new Opm::GasPvtThermal<Scalar>(*static_cast<const Opm::GasPvtThermal<Scalar>*>(data.realGasPvt_));
            break;
        case GasPvtApproach::Co2GasPvt:
            realGasPvt_ = new Opm::Co2GasPvt<Scalar>(*static_cast<const Opm::Co2GasPvt<Scalar>*>(data.realGasPvt_));
            break;
        default:
            break;
        }

        return *this;
    }

private:
    GasPvtApproach gasPvtApproach_;
    void* realGasPvt_;
};

#undef OPM_GAS_PVT_MULTIPLEXER_CALL

} // namespace Opm

#endif
