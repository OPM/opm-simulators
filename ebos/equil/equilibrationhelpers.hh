// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 *
 * \brief Auxiliary routines that to solve the ODEs that emerge from the hydrostatic
 *        equilibrium problem
 */
#ifndef EWOMS_EQUILIBRATIONHELPERS_HH
#define EWOMS_EQUILIBRATIONHELPERS_HH

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

// due to a bug in Equil.hpp, cstddef must be included before Equil.hpp
#include <cstddef>
#include <opm/parser/eclipse/EclipseState/InitConfig/Equil.hpp>

#include <memory>


/*
  ---- synopsis of EquilibrationHelpers.hpp ----

  namespace Opm
  {
  namespace EQUIL {

  namespace Miscibility {
  class RsFunction;
  class NoMixing;
  template <class FluidSystem>
  class RsVD;
  template <class FluidSystem>
  class RsSatAtContact;
  }

  class EquilReg;


  template <class FluidSystem,  class MaterialLaw, class MaterialLawManager>
  struct PcEq;

  template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
  double satFromPc(const MaterialLawManager& materialLawManager,
  const int phase,
  const int cell,
  const double targetPc,
  const bool increasing = false)
  template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
  double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
  const int phase1,
  const int phase2,
  const int cell,
  const double targetPc)
  } // namespace Equil
  } // namespace Opm

  ---- end of synopsis of EquilibrationHelpers.hpp ----
*/
namespace Ewoms {
/**
 * Types and routines that collectively implement a basic
 * ECLIPSE-style equilibration-based initialisation scheme.
 *
 * This namespace is intentionally nested to avoid name clashes
 * with other parts of OPM.
 */
namespace EQUIL {


typedef Opm::BlackOilFluidSystem<double> FluidSystemSimple;

// Adjust oil pressure according to gas saturation and cap pressure
typedef Opm::SimpleModularFluidState<double,
                                     /*numPhases=*/3,
                                     /*numComponents=*/3,
                                     FluidSystemSimple,
                                     /*storePressure=*/false,
                                     /*storeTemperature=*/false,
                                     /*storeComposition=*/false,
                                     /*storeFugacity=*/false,
                                     /*storeSaturation=*/true,
                                     /*storeDensity=*/false,
                                     /*storeViscosity=*/false,
                                     /*storeEnthalpy=*/false> SatOnlyFluidState;

/**
 * Types and routines relating to phase mixing in
 * equilibration calculations.
 */
namespace Miscibility {

/**
 * Base class for phase mixing functions.
 */
class RsFunction
{
public:
    virtual ~RsFunction() = default;

    /**
     * Function call operator.
     *
     * \param[in] depth Depth at which to calculate RS
     * value.
     *
     * \param[in] press Pressure at which to calculate RS
     * value.
     *
     * \param[in] temp Temperature at which to calculate RS
     * value.
     *
     * \return Dissolved gas-oil ratio (RS) at depth @c
     * depth and pressure @c press.
     */
    virtual double operator()(const double depth,
                              const double press,
                              const double temp,
                              const double sat = 0.0) const = 0;
};


/**
 * Type that implements "no phase mixing" policy.
 */
class NoMixing : public RsFunction
{
public:
    virtual ~NoMixing() = default;

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RS
     * value.
     *
     * \param[in] press Pressure at which to calculate RS
     * value.
     *
     * \param[in] temp Temperature at which to calculate RS
     * value.
     *
     * \return Dissolved gas-oil ratio (RS) at depth @c
     * depth and pressure @c press.  In "no mixing
     * policy", this is identically zero.
     */
    double
    operator()(const double /* depth */,
               const double /* press */,
               const double /* temp */,
               const double /* sat */ = 0.0) const
    {
        return 0.0;
    }
};


/**
 * Type that implements "dissolved gas-oil ratio"
 * tabulated as a function of depth policy.  Data
 * typically taken from keyword 'RSVD'.
 */
template <class FluidSystem>
class RsVD : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] depth Depth nodes.
     * \param[in] rs Dissolved gas-oil ratio at @c depth.
     */
    RsVD(const int pvtRegionIdx,
         const std::vector<double>& depth,
         const std::vector<double>& rs)
        : pvtRegionIdx_(pvtRegionIdx)
        , rsVsDepth_(depth, rs)
    {}

    virtual ~RsVD() = default;

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RS
     * value.
     *
     * \param[in] press Pressure at which to calculate RS
     * value.
     *
     * \param[in] temp Temperature at which to calculate RS
     * value.
     *
     * \return Dissolved gas-oil ratio (RS) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double depth,
                      const double press,
                      const double temp,
                      const double satGas = 0.0) const
    {
        if (satGas > 0.0) {
            return satRs(press, temp);
        }
        else {
            if (rsVsDepth_.xMin() > depth)
                return rsVsDepth_.valueAt(0);
            else if (rsVsDepth_.xMax() < depth)
                return rsVsDepth_.valueAt(rsVsDepth_.numSamples() - 1);
            return std::min(satRs(press, temp), rsVsDepth_.eval(depth, /*extrapolate=*/false));
        }
    }

private:
    typedef Opm::Tabulated1DFunction<double> RsVsDepthFunc;

    const int pvtRegionIdx_;
    RsVsDepthFunc rsVsDepth_;

    double satRs(const double press, const double temp) const
    {
        return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
    }
};

/**
 * Type that implements "dissolved gas-oil ratio"
 * tabulated as a function of depth policy.  Data
 * typically from keyword 'PBVD'.
 */
template <class FluidSystem>
class PBVD : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] depth Depth nodes.
     * \param[in] pbub Bubble-point pressure at @c depth.
     */
    PBVD(const int pvtRegionIdx,
         const std::vector<double>& depth,
         const std::vector<double>& pbub)
        : pvtRegionIdx_(pvtRegionIdx)
        , pbubVsDepth_(depth, pbub)
    {}

    virtual ~PBVD() = default;

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RS
     * value.
     *
     * \param[in] Pressure in the cell
     *
     * \param[in] temp Temperature at which to calculate RS
     * value.
     *
     * \return Dissolved gas-oil ratio (RS) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double depth,
                      const double cellPress,
                      const double temp,
                      const double satGas = 0.0) const
    {
        double press = cellPress;
        if (satGas <= 0.0) {
            if (pbubVsDepth_.xMin() > depth)
                press = pbubVsDepth_.valueAt(0);
            else if (pbubVsDepth_.xMax() < depth)
                press = pbubVsDepth_.valueAt(pbubVsDepth_.numSamples() - 1);
            else
                press = pbubVsDepth_.eval(depth, /*extrapolate=*/false);
        }
        return satRs(std::min(press, cellPress), temp);
    }

private:
    typedef Opm::Tabulated1DFunction<double> PbubVsDepthFunc;

    const int pvtRegionIdx_;
    PbubVsDepthFunc pbubVsDepth_;

    double satRs(const double press, const double temp) const
    {
        return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
    }
};

/**
 * Type that implements "vaporized oil-gas ratio"
 * tabulated as a function of depth policy.  Data
 * taken from keyword 'PDVD'.
 */
template <class FluidSystem>
class PDVD : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] depth Depth nodes.
     * \param[in] pbub Dew-point pressure at @c depth.
     */
    PDVD(const int pvtRegionIdx,
         const std::vector<double>& depth,
         const std::vector<double>& pdew)
        : pvtRegionIdx_(pvtRegionIdx)
        , pdewVsDepth_(depth, pdew)
    {}

    virtual ~PDVD() = default;

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RV
     * value.
     *
     * \param[in] cellPress Pressure in the cell
     *
     * \param[in] temp Temperature at which to calculate RV
     * value.
     *
     * \return Vaporized oil-gas ratio (RV) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double depth,
                      const double cellPress,
                      const double temp,
                      const double satOil = 0.0) const
    {
        double press = cellPress;
        if (satOil <= 0.0) {
            if (pdewVsDepth_.xMin() > depth)
                press = pdewVsDepth_.valueAt(0);
            else if (pdewVsDepth_.xMax() < depth)
                press = pdewVsDepth_.valueAt(pdewVsDepth_.numSamples() - 1);
            else
                press = pdewVsDepth_.eval(depth, /*extrapolate=*/false);
        }
        return satRv(std::min(press, cellPress), temp);
    }

private:
    typedef Opm::Tabulated1DFunction<double> PdewVsDepthFunc;

    const int pvtRegionIdx_;
    PdewVsDepthFunc pdewVsDepth_;

    double satRv(const double press, const double temp) const
    {
        return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
    }
};


/**
 * Type that implements "vaporized oil-gas ratio"
 * tabulated as a function of depth policy.  Data
 * typically taken from keyword 'RVVD'.
 */
template <class FluidSystem>
class RvVD : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] depth Depth nodes.
     * \param[in] rv Dissolved gas-oil ratio at @c depth.
     */
    RvVD(const int pvtRegionIdx,
         const std::vector<double>& depth,
         const std::vector<double>& rv)
        : pvtRegionIdx_(pvtRegionIdx)
        , rvVsDepth_(depth, rv)
    {}

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RV
     * value.
     *
     * \param[in] press Pressure at which to calculate RV
     * value.
     *
     * \param[in] temp Temperature at which to calculate RV
     * value.
     *
     * \return Vaporized oil-gas ratio (RV) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double depth,
                      const double press,
                      const double temp,
                      const double satOil = 0.0) const
    {
        if (std::abs(satOil) > 1e-16) {
            return satRv(press, temp);
        }
        else {
            if (rvVsDepth_.xMin() > depth)
                return rvVsDepth_.valueAt(0);
            else if (rvVsDepth_.xMax() < depth)
                return rvVsDepth_.valueAt(rvVsDepth_.numSamples() - 1);
            return std::min(satRv(press, temp), rvVsDepth_.eval(depth, /*extrapolate=*/false));
        }
    }

private:
    typedef Opm::Tabulated1DFunction<double> RvVsDepthFunc;

    const int pvtRegionIdx_;
    RvVsDepthFunc rvVsDepth_;

    double satRv(const double press, const double temp) const
    {
        return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
    }
};


/**
 * Class that implements "dissolved gas-oil ratio" (Rs)
 * as function of depth and pressure as follows:
 *
 *   1. The Rs at the gas-oil contact is equal to the
 *      saturated Rs value, RsSatContact.
 *
 *   2. The Rs elsewhere is equal to RsSatContact, but
 *      constrained to the saturated value as given by the
 *      local pressure.
 *
 * This should yield Rs-values that are constant below the
 * contact, and decreasing above the contact.
 */
template <class FluidSystem>
class RsSatAtContact : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] pContact  oil pressure at the contact
     * \param[in] T_contact  temperature at the contact
     */
    RsSatAtContact(const int pvtRegionIdx, const double pContact,  const double T_contact)
        : pvtRegionIdx_(pvtRegionIdx)
    {
        rsSatContact_ = satRs(pContact, T_contact);
    }

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RS
     * value.
     *
     * \param[in] press Pressure at which to calculate RS
     * value.
     *
     * \param[in] temp Temperature at which to calculate RS
     * value.
     *
     * \return Dissolved gas-oil ratio (RS) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double /* depth */,
                      const double press,
                      const double temp,
                      const double satGas = 0.0) const
    {
        if (satGas > 0.0) {
            return satRs(press, temp);
        }
        else {
            return std::min(satRs(press, temp), rsSatContact_);
        }
    }

private:
    const int pvtRegionIdx_;
    double rsSatContact_;

    double satRs(const double press, const double temp) const
    {
        return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
    }
};


/**
 * Class that implements "vaporized oil-gas ratio" (Rv)
 * as function of depth and pressure as follows:
 *
 *   1. The Rv at the gas-oil contact is equal to the
 *      saturated Rv value, RvSatContact.
 *
 *   2. The Rv elsewhere is equal to RvSatContact, but
 *      constrained to the saturated value as given by the
 *      local pressure.
 *
 * This should yield Rv-values that are constant below the
 * contact, and decreasing above the contact.
 */
template <class FluidSystem>
class RvSatAtContact : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] pContact  oil pressure at the contact
     * \param[in] T_contact  temperature at the contact
     */
    RvSatAtContact(const int pvtRegionIdx, const double pContact, const double T_contact)
        :pvtRegionIdx_(pvtRegionIdx)
    {
        rvSatContact_ = satRv(pContact, T_contact);
    }

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RV
     * value.
     *
     * \param[in] press Pressure at which to calculate RV
     * value.
     *
     * \param[in] temp Temperature at which to calculate RV
     * value.
     *
     * \return Dissolved oil-gas ratio (RV) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double /*depth*/,
                      const double press,
                      const double temp,
                      const double satOil = 0.0) const
    {
        if (satOil > 0.0) {
            return satRv(press, temp);
        }
        else {
            return std::min(satRv(press, temp), rvSatContact_);
        }
    }

private:
    const int pvtRegionIdx_;
    double rvSatContact_;

    double satRv(const double press, const double temp) const
    {
        return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);;
    }
};

} // namespace Miscibility

/**
 * Aggregate information base of an equilibration region.
 *
 * Provides inquiry methods for retrieving depths of contacs
 * and pressure values as well as a means of calculating fluid
 * densities, dissolved gas-oil ratio and vapourised oil-gas
 * ratios.
 *
 * \tparam DensCalc Type that provides access to a phase
 * density calculation facility.  Must implement an operator()
 * declared as
 * <CODE>
 * std::vector<double>
 * operator()(const double press,
 *            const std::vector<double>& svol)
 * </CODE>
 * that calculates the phase densities of all phases in @c
 * svol at fluid pressure @c press.
 */
class EquilReg
{
public:
    /**
     * Constructor.
     *
     * \param[in] rec     Equilibration data of current region.
     * \param[in] rs      Calculator of dissolved gas-oil ratio.
     * \param[in] rv      Calculator of vapourised oil-gas ratio.
     * \param[in] pvtRegionIdx The pvt region index
     */
    EquilReg(const Opm::EquilRecord& rec,
             std::shared_ptr<Miscibility::RsFunction> rs,
             std::shared_ptr<Miscibility::RsFunction> rv,
             const int pvtIdx)
        : rec_    (rec)
        , rs_     (rs)
        , rv_     (rv)
        , pvtIdx_ (pvtIdx)
    {}

    /**
     * Type of dissolved gas-oil ratio calculator.
     */
    typedef Miscibility::RsFunction CalcDissolution;

    /**
     * Type of vapourised oil-gas ratio calculator.
     */
    typedef Miscibility::RsFunction CalcEvaporation;

    /**
     * Datum depth in current region
     */
    double datum()    const { return this->rec_.datumDepth(); }

    /**
     * Pressure at datum depth in current region.
     */
    double pressure() const { return this->rec_.datumDepthPressure(); }

    /**
     * Depth of water-oil contact.
     */
    double zwoc()     const { return this->rec_.waterOilContactDepth(); }

    /**
     * water-oil capillary pressure at water-oil contact.
     *
     * \return P_o - P_w at WOC.
     */
    double pcowWoc() const { return this->rec_.waterOilContactCapillaryPressure(); }

    /**
     * Depth of gas-oil contact.
     */
    double zgoc()     const { return this->rec_.gasOilContactDepth(); }

    /**
     * Gas-oil capillary pressure at gas-oil contact.
     *
     * \return P_g - P_o at GOC.
     */
    double pcgoGoc() const { return this->rec_.gasOilContactCapillaryPressure(); }


    /**
     * Retrieve dissolved gas-oil ratio calculator of current
     * region.
     */
    const CalcDissolution&
    dissolutionCalculator() const { return *this->rs_; }

    /**
     * Retrieve vapourised oil-gas ratio calculator of current
     * region.
     */
    const CalcEvaporation&
    evaporationCalculator() const { return *this->rv_; }

    /**
     * Retrieve pvtIdx of the region.
     */
    int pvtIdx() const { return this->pvtIdx_; }


private:
    Opm::EquilRecord rec_;     /**< Equilibration data */
    std::shared_ptr<Miscibility::RsFunction> rs_;      /**< RS calculator */
    std::shared_ptr<Miscibility::RsFunction> rv_;      /**< RV calculator */
    const int pvtIdx_;
};



/// Functor for inverting capillary pressure function.
/// Function represented is
///   f(s) = pc(s) - targetPc
template <class FluidSystem,  class MaterialLaw, class MaterialLawManager>
struct PcEq
{
    PcEq(const MaterialLawManager& materialLawManager,
         const int phase,
         const int cell,
         const double targetPc)
        : materialLawManager_(materialLawManager),
          phase_(phase),
          cell_(cell),
          targetPc_(targetPc)
    {}

    double operator()(double s) const
    {
        const auto& matParams = materialLawManager_.materialLawParams(cell_);
        SatOnlyFluidState fluidState;
        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
        fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
        fluidState.setSaturation(phase_, s);

        double pc[FluidSystem::numPhases];
        std::fill(pc, pc + FluidSystem::numPhases, 0.0);
        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
        double sign = (phase_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
        double pcPhase = pc[FluidSystem::oilPhaseIdx] + sign *  pc[phase_];
        return pcPhase - targetPc_;
    }
private:
    const MaterialLawManager& materialLawManager_;
    const int phase_;
    const int cell_;
    const double targetPc_;
};

template <class FluidSystem, class MaterialLawManager>
double minSaturations(const MaterialLawManager& materialLawManager, const int phase, const int cell)
{
    const auto& scaledDrainageInfo =
        materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

    // Find minimum and maximum saturations.
    switch(phase) {
    case FluidSystem::waterPhaseIdx:
        return scaledDrainageInfo.Swl;

    case FluidSystem::gasPhaseIdx:
        return scaledDrainageInfo.Sgl;

    case FluidSystem::oilPhaseIdx:
        throw std::runtime_error("Min saturation not implemented for oil phase.");

    default:
        throw std::runtime_error("Unknown phaseIdx .");
    }
    return -1.0;
}

template <class FluidSystem, class MaterialLawManager>
double maxSaturations(const MaterialLawManager& materialLawManager, const int phase, const int cell)
{
    const auto& scaledDrainageInfo =
        materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

    // Find minimum and maximum saturations.
    switch(phase) {
    case FluidSystem::waterPhaseIdx:
        return scaledDrainageInfo.Swu;

    case FluidSystem::gasPhaseIdx:
        return scaledDrainageInfo.Sgu;

    case FluidSystem::oilPhaseIdx:
        throw std::runtime_error("Max saturation not implemented for oil phase.");

    default:
        throw std::runtime_error("Unknown phaseIdx .");
    }
    return -1.0;
}


/// Compute saturation of some phase corresponding to a given
/// capillary pressure.
template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
double satFromPc(const MaterialLawManager& materialLawManager,
                 const int phase,
                 const int cell,
                 const double targetPc,
                 const bool increasing = false)
{
    // Find minimum and maximum saturations.
    double s0 = increasing ? maxSaturations<FluidSystem>(materialLawManager, phase, cell) : minSaturations<FluidSystem>(materialLawManager, phase, cell);
    double s1 = increasing ? minSaturations<FluidSystem>(materialLawManager, phase, cell) : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

    // Create the equation f(s) = pc(s) - targetPc
    const PcEq<FluidSystem, MaterialLaw, MaterialLawManager> f(materialLawManager, phase, cell, targetPc);
    double f0 = f(s0);
    double f1 = f(s1);

    if (f0 <= 0.0)
        return s0;
    else if (f1 >= 0.0)
        return s1;

    assert(f0 > 0 && f1 < 0);

    const int maxIter = 60;
    const double tol = 1e-10;

    // regula falsi with the "Pegasus" method to avoid stagnation
    for (int iter = 0; iter < maxIter; ++ iter) {
        // determine the pivot
        double s = (s1*f0 - s0*f1)/(f0 - f1);

        // adapt the interval
        double fs = f(s);
        if (fs == 0.0)
            return s;
        else if ((fs > 0.0) == (f0 > 0.0)) {
            // update interval and reverse
            s0 = s1;
            f0 = f1;
        }
        else
            // "Pegasus" method
            f0 *= f1/(f1 + fs);

        s1 = s;
        f1 = fs;

        // check for convergence
        if (std::abs(s1 - s0) < tol)
            return (s1 + s0)/2;
    }

    throw std::runtime_error("Could not find solution for PcEq = 0.0 after "+std::to_string(maxIter)
                             +" iterations.");
}


/// Functor for inverting a sum of capillary pressure functions.
/// Function represented is
///   f(s) = pc1(s) + pc2(1 - s) - targetPc
template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
struct PcEqSum
{
    PcEqSum(const MaterialLawManager& materialLawManager,
            const int phase1,
            const int phase2,
            const int cell,
            const double targetPc)
        : materialLawManager_(materialLawManager),
          phase1_(phase1),
          phase2_(phase2),
          cell_(cell),
          targetPc_(targetPc)
    {}

    double operator()(double s) const
    {
        const auto& matParams = materialLawManager_.materialLawParams(cell_);
        SatOnlyFluidState fluidState;
        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
        fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
        fluidState.setSaturation(phase1_, s);
        fluidState.setSaturation(phase2_, 1.0 - s);

        double pc[FluidSystem::numPhases];
        std::fill(pc, pc + FluidSystem::numPhases, 0.0);

        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
        double sign1 = (phase1_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
        double pc1 = pc[FluidSystem::oilPhaseIdx] + sign1 *  pc[phase1_];
        double sign2 = (phase2_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
        double pc2 = pc[FluidSystem::oilPhaseIdx] + sign2 *  pc[phase2_];
        return pc1 + pc2 - targetPc_;
    }
private:
    const MaterialLawManager& materialLawManager_;
    const int phase1_;
    const int phase2_;
    const int cell_;
    const double targetPc_;
};




/// Compute saturation of some phase corresponding to a given
/// capillary pressure, where the capillary pressure function
/// is given as a sum of two other functions.
template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
                       const int phase1,
                       const int phase2,
                       const int cell,
                       const double targetPc)
{
    // Find minimum and maximum saturations.
    double s0 = minSaturations<FluidSystem>(materialLawManager, phase1, cell);
    double s1 = maxSaturations<FluidSystem>(materialLawManager, phase1, cell);

    // Create the equation f(s) = pc1(s) + pc2(1-s) - targetPc
    const PcEqSum<FluidSystem, MaterialLaw, MaterialLawManager> f(materialLawManager, phase1, phase2, cell, targetPc);
    double f0 = f(s0);
    double f1 = f(s1);
    if (f0 <= 0.0)
        return s0;
    else if (f1 >= 0.0)
        return s1;

    assert(f0 > 0.0 && f1 < 0.0);

    const int maxIter = 60;
    const double tol = 1e-10;

    // regula falsi with the "Pegasus" method to avoid stagnation
    for (int iter = 0; iter < maxIter; ++ iter) {
        // determine the pivot
        double s = (s1*f0 - s0*f1)/(f0 - f1);

        // adapt the interval
        double fs = f(s);
        if (fs == 0.0)
            return s;
        else if ((fs > 0.0) == (f0 > 0.0)) {
            // update interval and reverse
            s0 = s1;
            f0 = f1;
        }
        else
            // "Pegasus" method
            f0 *= f1 / (f1 + fs);

        s1 = s;
        f1 = fs;

        // check for convergence
        if (std::abs(s1 - s0) < tol)
            return (s1 + s0)/2;
    }

    throw std::runtime_error("Could not find solution for PcEqSum = 0.0 after "+std::to_string(maxIter)
                             +" iterations.");
}

/// Compute saturation from depth. Used for constant capillary pressure function
template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
double satFromDepth(const MaterialLawManager& materialLawManager,
                    const double cellDepth,
                    const double contactDepth,
                    const int phase,
                    const int cell,
                    const bool increasing = false)
{
    const double s0 = increasing ? maxSaturations<FluidSystem>(materialLawManager, phase, cell) : minSaturations<FluidSystem>(materialLawManager, phase, cell);
    const double s1 = increasing ? minSaturations<FluidSystem>(materialLawManager, phase, cell) : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

    if (cellDepth < contactDepth) {
        return s0;
    }
    else {
        return s1;
    }

}

/// Return true if capillary pressure function is constant
template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
bool isConstPc(const MaterialLawManager& materialLawManager,
               const int phase,
               const int cell)
{
    // Create the equation f(s) = pc(s);
    const PcEq<FluidSystem, MaterialLaw, MaterialLawManager> f(materialLawManager, phase, cell, 0);
    const double f0 = f(minSaturations<FluidSystem>(materialLawManager, phase, cell));
    const double f1 = f(maxSaturations<FluidSystem>(materialLawManager, phase, cell));
    return std::abs(f0 - f1) < std::numeric_limits<double>::epsilon();
}

} // namespace Equil
} // namespace Ewoms

#endif // EWOMS_EQUILIBRATIONHELPERS_HH
