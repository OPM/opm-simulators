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
 * \brief Contains the classes required to extend the black-oil model by solvent component.
 *  For details, refer:
 *  [*] T.H. Sandve, O. Sævareid and I. Aavatsmark: “Improved Extended Blackoil Formulation
 *  for CO2 EOR Simulations.” in ECMOR XVII – The 17th European Conference on the
 *  Mathematics of Oil Recovery,  September 2020.
 */
#ifndef EWOMS_BLACK_OIL_EXTBO_MODULE_HH
#define EWOMS_BLACK_OIL_EXTBO_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/models/blackoil/blackoilextboparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/models/utils/basicproperties.hh>

//#include <opm/models/io/vtkBlackOilExtboModule.hh> //TODO: Missing ...

#include <cassert>
#include <cmath>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace Opm {

/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model.
 */
template <class TypeTag, bool enableExtboV = getPropValue<TypeTag, Properties::EnableExtbo>()>
class BlackOilExtboModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilExtboParams<Scalar>::TabulatedFunction;
    using Tabulated2DFunction = typename BlackOilExtboParams<Scalar>::Tabulated2DFunction;

    static constexpr unsigned zFractionIdx = Indices::zFractionIdx;
    static constexpr unsigned contiZfracEqIdx = Indices::contiZfracEqIdx;
    static constexpr unsigned enableExtbo = enableExtboV;
    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool blackoilConserveSurfaceVolume = getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>();

public:
    //! \brief Set parameters.
    static void setParams(BlackOilExtboParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model&,
                                      Simulator&)
    {}

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if constexpr (enableExtbo) {
            return pvIdx == zFractionIdx;
        }
        else {
            return false;
        }
    }

    static std::string primaryVarName([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        return "z_fraction";
    }

    static Scalar primaryVarWeight([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableExtbo) {
            return eqIdx == contiZfracEqIdx;
        }
        else {
            return false;
        }
    }

    static std::string eqName([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        return "conti^solvent";
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if constexpr (enableExtbo) {
            if constexpr (blackoilConserveSurfaceVolume) {
                storage[contiZfracEqIdx] =
                        Toolbox::template decay<LhsEval>(intQuants.porosity()) *
                        Toolbox::template decay<LhsEval>(intQuants.yVolume()) *
                        Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx)) *
                        Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx));
                if (FluidSystem::enableDissolvedGas()) { // account for dissolved z in oil phase
                    storage[contiZfracEqIdx] +=
                        Toolbox::template decay<LhsEval>(intQuants.porosity()) *
                        Toolbox::template decay<LhsEval>(intQuants.xVolume()) *
                        Toolbox::template decay<LhsEval>(intQuants.fluidState().Rs()) *
                        Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(oilPhaseIdx)) *
                        Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(oilPhaseIdx));
                }
                // Reg. terms: Preliminary attempt to avoid singular behaviour when solvent is invading a pure water
                //             region. Results seems insensitive to the weighting factor.
                // TODO: Further investigations ...
                const Scalar regWghtFactor = 1.0e-6;
                storage[contiZfracEqIdx] +=
                    regWghtFactor* (1.0 - Toolbox::template decay<LhsEval>(intQuants.zFraction())) +
                                    regWghtFactor*Toolbox::template decay<LhsEval>(intQuants.porosity()) *
                                   Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx)) *
                                   Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx));
                storage[contiZfracEqIdx - 1] += regWghtFactor*Toolbox::template decay<LhsEval>(intQuants.zFraction());
            }
            else {
                throw std::runtime_error("Only component conservation in terms of surface volumes is implemented. ");
            }
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableExtbo) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

            if constexpr (blackoilConserveSurfaceVolume) {
                const unsigned inIdx = extQuants.interiorIndex();

                const unsigned upIdxGas = static_cast<unsigned>(extQuants.upstreamIndex(gasPhaseIdx));
                const auto& upGas = elemCtx.intensiveQuantities(upIdxGas, timeIdx);
                const auto& fsGas = upGas.fluidState();
                if (upIdxGas == inIdx) {
                    flux[contiZfracEqIdx] =
                        extQuants.volumeFlux(gasPhaseIdx) *
                        upGas.yVolume() *
                        fsGas.invB(gasPhaseIdx);
                }
                else {
                    flux[contiZfracEqIdx] =
                            extQuants.volumeFlux(gasPhaseIdx) *
                            decay<Scalar>(upGas.yVolume()) *
                            decay<Scalar>(fsGas.invB(gasPhaseIdx));
                }
                if (FluidSystem::enableDissolvedGas()) { // account for dissolved z in oil phase
                    const unsigned upIdxOil = static_cast<unsigned>(extQuants.upstreamIndex(oilPhaseIdx));
                    const auto& upOil = elemCtx.intensiveQuantities(upIdxOil, timeIdx);
                    const auto& fsOil = upOil.fluidState();
                    if (upIdxOil == inIdx) {
                        flux[contiZfracEqIdx] +=
                            extQuants.volumeFlux(oilPhaseIdx) *
                            upOil.xVolume() *
                            fsOil.Rs() *
                            fsOil.invB(oilPhaseIdx);
                    }
                    else {
                        flux[contiZfracEqIdx] +=
                            extQuants.volumeFlux(oilPhaseIdx) *
                            decay<Scalar>(upOil.xVolume()) *
                            decay<Scalar>(fsOil.Rs()) *
                            decay<Scalar>(fsOil.invB(oilPhaseIdx));
                    }
                }
            }
            else {
                throw std::runtime_error("Only component conservation in terms of surface volumes is implemented. ");
            }
        }
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar zFraction)
    {
        if constexpr (enableExtbo) {
            priVars[zFractionIdx] = zFraction;
        }
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if constexpr (enableExtbo) {
            // do a plain unchopped Newton update
            newPv[zFractionIdx] = oldPv[zFractionIdx] - delta[zFractionIdx];
        }
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables&,
                                     const EqVector&)
    {
        // do not consider consider the cange of solvent primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of solvents when it comes to convergence
        return std::abs(Toolbox::scalarValue(resid[contiZfracEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if constexpr (enableExtbo) {
            const unsigned dofIdx = model.dofMapper().index(dof);

            const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
            outstream << priVars[zFractionIdx];
        }
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if constexpr (enableExtbo) {
            const unsigned dofIdx = model.dofMapper().index(dof);

            PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
            PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

            instream >> priVars0[zFractionIdx];

            // set the primary variables for the beginning of the current time step.
            priVars1 = priVars0[zFractionIdx];
        }
    }

    template <typename Value>
    static Value xVolume(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.X_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value yVolume(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.Y_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value pbubRs(unsigned pvtRegionIdx, const Value& z, const Value& rs)
    { return params_.PBUB_RS_[pvtRegionIdx].eval(z, rs, true); }

    template <typename Value>
    static Value pbubRv(unsigned pvtRegionIdx, const Value& z, const Value& rv)
    { return params_.PBUB_RV_[pvtRegionIdx].eval(z, rv, true); }

    template <typename Value>
    static Value oilViscosity(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.VISCO_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value gasViscosity(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.VISCG_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value bo(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.BO_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value bg(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.BG_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value rs(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.RS_[pvtRegionIdx].eval(z, pressure, true); }

    template <typename Value>
    static Value rv(unsigned pvtRegionIdx, const Value& pressure, const Value& z)
    { return params_.RV_[pvtRegionIdx].eval(z, pressure, true); }

    static Scalar referenceDensity(unsigned regionIdx)
    { return params_.zReferenceDensity_[regionIdx]; }

    static Scalar zLim(unsigned regionIdx)
    { return params_.zLim_[regionIdx]; }

    template <typename Value>
    static Value oilCmp(unsigned pvtRegionIdx, const Value& z)
    { return params_.oilCmp_[pvtRegionIdx].eval(z, true); }

    template <typename Value>
    static Value gasCmp(unsigned pvtRegionIdx, const Value& z)
    { return params_.gasCmp_[pvtRegionIdx].eval(z, true); }

private:
    static BlackOilExtboParams<Scalar> params_;
};

template <class TypeTag, bool enableExtboV>
BlackOilExtboParams<typename BlackOilExtboModule<TypeTag, enableExtboV>::Scalar>
BlackOilExtboModule<TypeTag, enableExtboV>::params_;

template <class TypeTag, bool enableExtboV>
class BlackOilExtboIntensiveQuantities;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilExtboIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag>
class BlackOilExtboIntensiveQuantities<TypeTag, /*enableExtboV=*/true>
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using ExtboModule = BlackOilExtboModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int zFractionIdx = Indices::zFractionIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr double cutOff = 1e-12;

public:
    /*!
     * \brief Compute extended pvt properties from table lookups.
     *
     *  At this point the pressures of the fluid state are correct.
     */
    void zFractionUpdate_(const ElementContext& elemCtx,
                          unsigned dofIdx,
                          unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        const auto& fs = asImp_().fluidState_;

        zFraction_ = priVars.makeEvaluation(zFractionIdx, timeIdx);

        oilViscosity_ = ExtboModule::oilViscosity(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        gasViscosity_ = ExtboModule::gasViscosity(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);

        bo_ = ExtboModule::bo(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        bg_ = ExtboModule::bg(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);

        bz_ = ExtboModule::bg(pvtRegionIdx, fs.pressure(oilPhaseIdx), Evaluation{0.99});

        if (FluidSystem::enableDissolvedGas()) {
            rs_ = ExtboModule::rs(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        }
        else {
            rs_ = 0.0;
        }

        if (FluidSystem::enableVaporizedOil()) {
            rv_ = ExtboModule::rv(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);
        }
        else {
            rv_ = 0.0;
        }

        xVolume_ = ExtboModule::xVolume(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        yVolume_ = ExtboModule::yVolume(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);

        Evaluation pbub = fs.pressure(oilPhaseIdx);

        if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
           static const Scalar thresholdWaterFilledCell = 1.0 - 1e-6;
           Scalar sw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx).value();
           if (sw >= thresholdWaterFilledCell) {
               rs_ = 0.0;  // water only, zero rs_ ...
           }
        }

        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
           rs_ = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
           const Evaluation zLim = ExtboModule::zLim(pvtRegionIdx);
           if (zFraction_ > zLim) {
               pbub = ExtboModule::pbubRs(pvtRegionIdx, zLim, rs_);
           } else {
               pbub = ExtboModule::pbubRs(pvtRegionIdx, zFraction_, rs_);
           }
           bo_ = ExtboModule::bo(pvtRegionIdx, pbub, zFraction_) +
                 ExtboModule::oilCmp(pvtRegionIdx, zFraction_) * (fs.pressure(oilPhaseIdx) - pbub);

           xVolume_ = ExtboModule::xVolume(pvtRegionIdx, pbub, zFraction_);
        }

        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
           rv_ = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
           const Evaluation rvsat = ExtboModule::rv(pvtRegionIdx, pbub, zFraction_);
           bg_ = ExtboModule::bg(pvtRegionIdx, pbub, zFraction_) +
               ExtboModule::gasCmp(pvtRegionIdx, zFraction_) * (rv_ - rvsat);

           yVolume_ = ExtboModule::yVolume(pvtRegionIdx, pbub, zFraction_);
        }
    }

    /*!
     * \brief Re-compute face densities to account for zFraction dependency.
     *
     * At this point the pressures and saturations of the fluid state are correct.
     */
    void zPvtUpdate_()
    {
        const auto& iq = asImp_();
        auto& fs = asImp_().fluidState_;

        const unsigned pvtRegionIdx = iq.pvtRegionIndex();
        zRefDensity_ = ExtboModule::referenceDensity(pvtRegionIdx);

        fs.setInvB(oilPhaseIdx, 1.0 / bo_);
        fs.setInvB(gasPhaseIdx, 1.0 / bg_);

        fs.setDensity(oilPhaseIdx,
                      fs.invB(oilPhaseIdx) *
                      (FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx) +
                        (1.0 - xVolume_) * fs.Rs() * FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx) +
                        xVolume_ * fs.Rs() * zRefDensity_));
        fs.setDensity(gasPhaseIdx,
                      fs.invB(gasPhaseIdx) *
                      (FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx) *
                       (1.0 - yVolume_) + yVolume_* zRefDensity_ +
                        FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx) * fs.Rv()));
    }

    const Evaluation& zFraction() const
    { return zFraction_; }

    const Evaluation& xVolume() const
    { return xVolume_; }

    const Evaluation& yVolume() const
    { return yVolume_; }

    const Evaluation& oilViscosity() const
    { return oilViscosity_; }

    const Evaluation& gasViscosity() const
    { return gasViscosity_; }

    const Evaluation& bo() const
    { return bo_; }

    const Evaluation& bg() const
    { return bg_; }

    const Evaluation& rs() const
    { return rs_; }

    const Evaluation& rv() const
    { return rv_; }

    const Evaluation zPureInvFormationVolumeFactor() const
    { return 1.0 / bz_; }

    Scalar zRefDensity() const
    { return zRefDensity_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    // Abstract "mass fraction" accounting for the solvent component. The relation between this
    // quantity and the actual mass fraction of solvent, is implicitly defined from the specific
    // pvt measurements as provided by kw PVTSOL.
    Evaluation zFraction_;

    // The solvent component is assumed gas at surface conditions
    Evaluation xVolume_; // Solvent volume fraction of Rs
    Evaluation yVolume_; // Solvent volume fraction of Sg/Bg

    // Standard black oil parameters modified for presence of solvent
    Evaluation oilViscosity_;
    Evaluation gasViscosity_;
    Evaluation bo_;
    Evaluation bg_;
    Evaluation rs_;
    Evaluation rv_;

    // Properties of pure solvent
    Evaluation bz_;
    Scalar zRefDensity_;
};

template <class TypeTag>
class BlackOilExtboIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void zPvtUpdate_()
    {}

    void zFractionUpdate_(const ElementContext&,
                          unsigned,
                          unsigned)
    {}

    const Evaluation& xVolume() const
    { throw std::runtime_error("xVolume() called but extbo is disabled"); }

    const Evaluation& yVolume() const
    { throw std::runtime_error("yVolume() called but extbo is disabled"); }

    const Evaluation& oilViscosity() const
    { throw std::runtime_error("oilViscosity() called but extbo is disabled"); }

    const Evaluation& gasViscosity() const
    { throw std::runtime_error("gasViscosity() called but extbo is disabled"); }

    const Evaluation& rs() const
    { throw std::runtime_error("rs() called but extbo is disabled"); }

    const Evaluation& rv() const
    { throw std::runtime_error("rv() called but extbo is disabled"); }

    const Evaluation& zPureInvFormationVolumeFactor() const
    { throw std::runtime_error("zPureInvFormationVolumeFactor() called but extbo is disabled"); }

    const Evaluation& zFraction() const
    { throw std::runtime_error("zFraction() called but extbo is disabled"); }

    const Evaluation& zInverseFormationVolumeFactor() const
    { throw std::runtime_error("zInverseFormationVolumeFactor() called but extbo is disabled"); }

    Scalar zRefDensity() const
    { throw std::runtime_error("zRefDensity() called but extbo is disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilExtboExtensiveQuantities
 *
 * \brief Provides the solvent specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableExtboV = getPropValue<TypeTag, Properties::EnableExtbo>()>
class BlackOilExtboExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using Toolbox = MathToolbox<Evaluation>;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
};

template <class TypeTag>
class BlackOilExtboExtensiveQuantities<TypeTag, false>
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
};

} // namespace Opm

#endif
