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
 * \copydoc Opm::OilPvtMultiplexer
 */
#ifndef OPM_OIL_PVT_MULTIPLEXER_HPP
#define OPM_OIL_PVT_MULTIPLEXER_HPP

#include "ConstantCompressibilityOilPvt.hpp"
#include "DeadOilPvt.hpp"
#include "LiveOilPvt.hpp"
#include "OilPvtThermal.hpp"
#include "BrineCo2Pvt.hpp"

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#endif

namespace Opm {
#define OPM_OIL_PVT_MULTIPLEXER_CALL(codeToCall)                                     \
    switch (approach_) {                                                             \
    case OilPvtApproach::ConstantCompressibilityOilPvt: {                            \
        auto& pvtImpl = getRealPvt<OilPvtApproach::ConstantCompressibilityOilPvt>(); \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::DeadOilPvt: {                                               \
        auto& pvtImpl = getRealPvt<OilPvtApproach::DeadOilPvt>();                    \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::LiveOilPvt: {                                               \
        auto& pvtImpl = getRealPvt<OilPvtApproach::LiveOilPvt>();                    \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::ThermalOilPvt: {                                            \
        auto& pvtImpl = getRealPvt<OilPvtApproach::ThermalOilPvt>();                 \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::BrineCo2Pvt: {                                              \
        auto& pvtImpl = getRealPvt<OilPvtApproach::BrineCo2Pvt>();                   \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::NoOilPvt:                                                   \
        throw std::logic_error("Not implemented: Oil PVT of this deck!");            \
    }                                                                                \

enum class OilPvtApproach {
    NoOilPvt,
    LiveOilPvt,
    DeadOilPvt,
    ConstantCompressibilityOilPvt,
    ThermalOilPvt,
    BrineCo2Pvt
};

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil
 *        phase in the black-oil model.
 *
 * This is the base class which which provides an API for the actual PVT implementation
 * classes which based on dynamic polymorphism. The rationale to use dynamic polymorphism
 * here is that this enables the fluid system to easily switch the used PVT relations for
 * the individual fluid phases.
 *
 * Note that, since the application for this class is the black-oil fluid system, the API
 * exposed by this class is pretty specific to the black-oil model.
 */
template <class Scalar, bool enableThermal = true>
class OilPvtMultiplexer
{
public:
    OilPvtMultiplexer()
    {
        approach_ = OilPvtApproach::NoOilPvt;
        realOilPvt_ = nullptr;
    }

    OilPvtMultiplexer(OilPvtApproach approach, void* realOilPvt)
        : approach_(approach)
        , realOilPvt_(realOilPvt)
    { }

    OilPvtMultiplexer(const OilPvtMultiplexer<Scalar,enableThermal>& data)
    {
        *this = data;
    }

    ~OilPvtMultiplexer()
    {
        switch (approach_) {
        case OilPvtApproach::LiveOilPvt: {
            delete &getRealPvt<OilPvtApproach::LiveOilPvt>();
            break;
        }
        case OilPvtApproach::DeadOilPvt: {
            delete &getRealPvt<OilPvtApproach::DeadOilPvt>();
            break;
        }
        case OilPvtApproach::ConstantCompressibilityOilPvt: {
            delete &getRealPvt<OilPvtApproach::ConstantCompressibilityOilPvt>();
            break;
        }
        case OilPvtApproach::ThermalOilPvt: {
            delete &getRealPvt<OilPvtApproach::ThermalOilPvt>();
            break;
        }
        case OilPvtApproach::BrineCo2Pvt: {
            delete &getRealPvt<OilPvtApproach::BrineCo2Pvt>();
            break;
        }

        case OilPvtApproach::NoOilPvt:
            break;
        }
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for water using an ECL state.
     *
     * This method assumes that the deck features valid DENSITY and PVTO/PVDO/PVCDO keywords.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule)
    {
        if (!eclState.runspec().phases().active(Phase::OIL))
            return;
        // TODO move the BrineCo2 approach to the waterPvtMultiplexer
        // when a proper gas-water simulator is supported
        if (eclState.runspec().co2Storage())
            setApproach(OilPvtApproach::BrineCo2Pvt);
        else if (enableThermal && eclState.getSimulationConfig().isThermal())
            setApproach(OilPvtApproach::ThermalOilPvt);
        else if (!eclState.getTableManager().getPvcdoTable().empty())
            setApproach(OilPvtApproach::ConstantCompressibilityOilPvt);
        else if (eclState.getTableManager().hasTables("PVDO"))
            setApproach(OilPvtApproach::DeadOilPvt);
        else if (!eclState.getTableManager().getPvtoTables().empty())
            setApproach(OilPvtApproach::LiveOilPvt);

        OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule));
    }
#endif // HAVE_ECL_INPUT


    void initEnd()
    { OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions()); return 1; }

    /*!
     * \brief Return the reference density which are considered by this PVT-object.
     */
    const Scalar oilReferenceDensity(unsigned regionIdx) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.oilReferenceDensity(regionIdx)); return 700.; }

    /*!
     * \brief Returns the specific enthalpy [J/kg] oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             const Evaluation& maxOilSaturation) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation)); return 0; }

    /*!
     * \brief Returns the saturation pressure [Pa] of oil given the mass fraction of the
     *        gas component in the oil phase.
     *
     * Calling this method only makes sense for live oil. All other implementations of
     * the black-oil PVT interface will just throw an exception...
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rs)); return 0; }

    /*!
     * \copydoc BaseFluidSystem::diffusionCoefficient
     */
    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& temperature,
                                    const Evaluation& pressure,
                                    unsigned compIdx) const
    {
      OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.diffusionCoefficient(temperature, pressure, compIdx)); return 0;
    }

    void setApproach(OilPvtApproach appr)
    {
        switch (appr) {
        case OilPvtApproach::LiveOilPvt:
            realOilPvt_ = new LiveOilPvt<Scalar>;
            break;

        case OilPvtApproach::DeadOilPvt:
            realOilPvt_ = new DeadOilPvt<Scalar>;
            break;

        case OilPvtApproach::ConstantCompressibilityOilPvt:
            realOilPvt_ = new ConstantCompressibilityOilPvt<Scalar>;
            break;

        case OilPvtApproach::ThermalOilPvt:
            realOilPvt_ = new OilPvtThermal<Scalar>;
            break;

        case OilPvtApproach::BrineCo2Pvt:
            realOilPvt_ = new BrineCo2Pvt<Scalar>;
            break;

        case OilPvtApproach::NoOilPvt:
            throw std::logic_error("Not implemented: Oil PVT of this deck!");
        }

        approach_ = appr;
    }

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    OilPvtApproach approach() const
    { return approach_; }

    // get the concrete parameter object for the oil phase
    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::LiveOilPvt, LiveOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<LiveOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::LiveOilPvt, const LiveOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<LiveOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::DeadOilPvt, DeadOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<DeadOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::DeadOilPvt, const DeadOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<DeadOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::ConstantCompressibilityOilPvt, ConstantCompressibilityOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<ConstantCompressibilityOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::ConstantCompressibilityOilPvt, const ConstantCompressibilityOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<ConstantCompressibilityOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::ThermalOilPvt, OilPvtThermal<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<OilPvtThermal<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::ThermalOilPvt, const OilPvtThermal<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<const OilPvtThermal<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::BrineCo2Pvt, BrineCo2Pvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<BrineCo2Pvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == OilPvtApproach::BrineCo2Pvt, const BrineCo2Pvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<const BrineCo2Pvt<Scalar>* >(realOilPvt_);
    }

    const void* realOilPvt() const { return realOilPvt_; }

    bool operator==(const OilPvtMultiplexer<Scalar,enableThermal>& data) const
    {
        if (this->approach() != data.approach())
            return false;

        switch (approach_) {
        case OilPvtApproach::ConstantCompressibilityOilPvt:
            return *static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_);
        case OilPvtApproach::DeadOilPvt:
            return *static_cast<const DeadOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const DeadOilPvt<Scalar>*>(data.realOilPvt_);
        case OilPvtApproach::LiveOilPvt:
            return *static_cast<const LiveOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const LiveOilPvt<Scalar>*>(data.realOilPvt_);
        case OilPvtApproach::ThermalOilPvt:
            return *static_cast<const OilPvtThermal<Scalar>*>(realOilPvt_) ==
                   *static_cast<const OilPvtThermal<Scalar>*>(data.realOilPvt_);
        case OilPvtApproach::BrineCo2Pvt:
            return *static_cast<const BrineCo2Pvt<Scalar>*>(realOilPvt_) ==
                    *static_cast<const BrineCo2Pvt<Scalar>*>(data.realOilPvt_);
        default:
            return true;
        }
    }

    OilPvtMultiplexer<Scalar,enableThermal>& operator=(const OilPvtMultiplexer<Scalar,enableThermal>& data)
    {
        approach_ = data.approach_;
        switch (approach_) {
        case OilPvtApproach::ConstantCompressibilityOilPvt:
            realOilPvt_ = new ConstantCompressibilityOilPvt<Scalar>(*static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case OilPvtApproach::DeadOilPvt:
            realOilPvt_ = new DeadOilPvt<Scalar>(*static_cast<const DeadOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case OilPvtApproach::LiveOilPvt:
            realOilPvt_ = new LiveOilPvt<Scalar>(*static_cast<const LiveOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case OilPvtApproach::ThermalOilPvt:
            realOilPvt_ = new OilPvtThermal<Scalar>(*static_cast<const OilPvtThermal<Scalar>*>(data.realOilPvt_));
            break;
        case OilPvtApproach::BrineCo2Pvt:
            realOilPvt_ = new BrineCo2Pvt<Scalar>(*static_cast<const BrineCo2Pvt<Scalar>*>(data.realOilPvt_));
            break;
        default:
            break;
        }

        return *this;
    }

private:
    OilPvtApproach approach_;
    void* realOilPvt_;
};

#undef OPM_OIL_PVT_MULTIPLEXER_CALL

} // namespace Opm

#endif
