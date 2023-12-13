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

#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>

#include <cmath>
#include <memory>
#include <vector>

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


  template <class FluidSystem, class MaterialLawManager>
  struct PcEq;

  template <class FluidSystem, class MaterialLawManager>
  double satFromPc(const MaterialLawManager& materialLawManager,
  const int phase,
  const int cell,
  const double targetPc,
  const bool increasing = false)
  template <class FluidSystem, class MaterialLawManager>
  double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
  const int phase1,
  const int phase2,
  const int cell,
  const double targetPc)
  } // namespace Equil
  } // namespace Opm

  ---- end of synopsis of EquilibrationHelpers.hpp ----
*/

namespace Opm {

/**
 * Types and routines that collectively implement a basic
 * ECLIPSE-style equilibration-based initialisation scheme.
 *
 * This namespace is intentionally nested to avoid name clashes
 * with other parts of OPM.
 */
namespace EQUIL {

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
         const std::vector<double>& rs);

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
                      const double satGas = 0.0) const;

private:
    using RsVsDepthFunc = Tabulated1DFunction<double>;

    const int pvtRegionIdx_;
    RsVsDepthFunc rsVsDepth_;

    double satRs(const double press, const double temp) const;
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
         const std::vector<double>& pbub);

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
                      const double satGas = 0.0) const;

private:
    using PbubVsDepthFunc = Tabulated1DFunction<double>;

    const int pvtRegionIdx_;
    PbubVsDepthFunc pbubVsDepth_;

    double satRs(const double press, const double temp) const;
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
         const std::vector<double>& pdew);

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
                      const double satOil = 0.0) const;

private:
    using PdewVsDepthFunc = Tabulated1DFunction<double>;

    const int pvtRegionIdx_;
    PdewVsDepthFunc pdewVsDepth_;

    double satRv(const double press, const double temp) const;
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
         const std::vector<double>& rv);

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
                      const double satOil = 0.0) const;

private:
    using RvVsDepthFunc = Tabulated1DFunction<double>;

    const int pvtRegionIdx_;
    RvVsDepthFunc rvVsDepth_;

    double satRv(const double press, const double temp) const;
};


/**
 * Type that implements "vaporized water-gas ratio"
 * tabulated as a function of depth policy.  Data
 * typically taken from keyword 'RVWVD'.
 */
template <class FluidSystem>
class RvwVD : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] depth Depth nodes.
     * \param[in] rvw Evaporized water-gasl ratio at @c depth.
     */
    RvwVD(const int pvtRegionIdx,
         const std::vector<double>& depth,
         const std::vector<double>& rvw);

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RVW
     * value.
     *
     * \param[in] press Pressure at which to calculate RVW
     * value.
     *
     * \param[in] temp Temperature at which to calculate RVW
     * value.
     *
     * \return Vaporized water-gas ratio (RVW) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double depth,
                      const double press,
                      const double temp,
                      const double satWat = 0.0) const;

private:
    using RvwVsDepthFunc = Tabulated1DFunction<double>;

    const int pvtRegionIdx_;
    RvwVsDepthFunc rvwVsDepth_;

    double satRvw(const double press, const double temp) const;
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
    RsSatAtContact(const int pvtRegionIdx, const double pContact,  const double T_contact);

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
                      const double satGas = 0.0) const;

private:
    const int pvtRegionIdx_;
    double rsSatContact_;

    double satRs(const double press, const double temp) const;
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
    RvSatAtContact(const int pvtRegionIdx, const double pContact, const double T_contact);

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
                      const double satOil = 0.0) const;

private:
    const int pvtRegionIdx_;
    double rvSatContact_;

    double satRv(const double press, const double temp) const;
};

/**
 * Class that implements "vaporized water-gas ratio" (Rvw)
 * as function of depth and pressure as follows:
 *
 *   1. The Rvw at the gas-water contact is equal to the
 *      saturated Rv value, RvwSatContact.
 *
 *   2. The Rvw elsewhere is equal to RvwSatContact, but
 *      constrained to the saturated value as given by the
 *      local pressure.
 *
 * This should yield Rvw-values that are constant below the
 * contact, and decreasing above the contact.
 */
template <class FluidSystem>
class RvwSatAtContact : public RsFunction
{
public:
    /**
     * Constructor.
     *
     * \param[in] pvtRegionIdx The pvt region index
     * \param[in] pContact  oil pressure at the contact
     * \param[in] T_contact  temperature at the contact
     */
    RvwSatAtContact(const int pvtRegionIdx, const double pContact, const double T_contact);

    /**
     * Function call.
     *
     * \param[in] depth Depth at which to calculate RVW
     * value.
     *
     * \param[in] press Pressure at which to calculate RVW
     * value.
     *
     * \param[in] temp Temperature at which to calculate RVW
     * value.
     *
     * \return Dissolved water-gas ratio (RVW) at depth @c
     * depth and pressure @c press.
     */
    double operator()(const double /*depth*/,
                      const double press,
                      const double temp,
                      const double satWat = 0.0) const;

private:
    const int pvtRegionIdx_;
    double rvwSatContact_;

    double satRvw(const double press, const double temp) const;
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
    using TabulatedFunction = Tabulated1DFunction<double>;

public:
    /**
     * Constructor.
     *
     * \param[in] rec     Equilibration data of current region.
     * \param[in] rs      Calculator of dissolved gas-oil ratio.
     * \param[in] rv      Calculator of vapourised oil-gas ratio.
     * \param[in] rvw     Calculator of vapourised water-gas ratio.
     * \param[in] pvtRegionIdx The pvt region index
     */
    EquilReg(const EquilRecord& rec,
             std::shared_ptr<Miscibility::RsFunction> rs,
             std::shared_ptr<Miscibility::RsFunction> rv,
             std::shared_ptr<Miscibility::RsFunction> rvw,
             const TabulatedFunction& tempVdTable,
             const TabulatedFunction& saltVdTable,
             const int pvtIdx);

    /**
     * Type of dissolved gas-oil ratio calculator.
     */
    using CalcDissolution = Miscibility::RsFunction;

    /**
     * Type of vapourised oil-gas ratio calculator.
     */
    using CalcEvaporation = Miscibility::RsFunction;

     /**
     * Type of vapourised water-gas ratio calculator.
     */
    using CalcWaterEvaporation = Miscibility::RsFunction;


    /**
     * Datum depth in current region
     */
    double datum() const;

    /**
     * Pressure at datum depth in current region.
     */
    double pressure() const;

    /**
     * Depth of water-oil contact.
     */
    double zwoc() const;

    /**
     * water-oil capillary pressure at water-oil contact.
     *
     * \return P_o - P_w at WOC.
     */
    double pcowWoc() const;

    /**
     * Depth of gas-oil contact.
     */
    double zgoc() const;

    /**
     * Gas-oil capillary pressure at gas-oil contact.
     *
     * \return P_g - P_o at GOC.
     */
    double pcgoGoc() const;

    /**
     * Accuracy/strategy for initial fluid-in-place calculation.
     *
     * \return zero (N=0) for centre-point method, negative (N<0) for the
     *   horizontal subdivision method with 2*(-N) intervals, and positive
     *   (N>0) for the tilted subdivision method with 2*N intervals.
     */
    int equilibrationAccuracy() const;

    /**
     * Retrieve dissolved gas-oil ratio calculator of current
     * region.
     */
    const CalcDissolution& dissolutionCalculator() const;

    /**
     * Retrieve vapourised oil-gas ratio calculator of current
     * region.
     */
    const CalcEvaporation& evaporationCalculator() const;

     /**
     * Retrieve vapourised water-gas ratio calculator of current
     * region.
     */
    const CalcWaterEvaporation& waterEvaporationCalculator() const;

    const TabulatedFunction& saltVdTable() const;
    const TabulatedFunction& tempVdTable() const;

    /**
     * Retrieve pvtIdx of the region.
     */
    int pvtIdx() const;

private:
    EquilRecord rec_;     /**< Equilibration data */
    std::shared_ptr<Miscibility::RsFunction> rs_;      /**< RS calculator */
    std::shared_ptr<Miscibility::RsFunction> rv_;      /**< RV calculator */
    std::shared_ptr<Miscibility::RsFunction> rvw_;      /**< RVW calculator */
    const TabulatedFunction& tempVdTable_;
    const TabulatedFunction& saltVdTable_;
    const int pvtIdx_;
};



/// Functor for inverting capillary pressure function.
/// Function represented is
///   f(s) = pc(s) - targetPc
template <class FluidSystem, class MaterialLawManager>
struct PcEq
{
    PcEq(const MaterialLawManager& materialLawManager,
         const int phase,
         const int cell,
         const double targetPc);

    double operator()(double s) const;

private:
    const MaterialLawManager& materialLawManager_;
    const int phase_;
    const int cell_;
    const double targetPc_;
};

template <class FluidSystem, class MaterialLawManager>
double minSaturations(const MaterialLawManager& materialLawManager,
                      const int phase, const int cell);

template <class FluidSystem, class MaterialLawManager>
double maxSaturations(const MaterialLawManager& materialLawManager,
                      const int phase, const int cell);

/// Compute saturation of some phase corresponding to a given
/// capillary pressure.
template <class FluidSystem, class MaterialLawManager>
double satFromPc(const MaterialLawManager& materialLawManager,
                 const int phase,
                 const int cell,
                 const double targetPc,
                 const bool increasing = false);

/// Functor for inverting a sum of capillary pressure functions.
/// Function represented is
///   f(s) = pc1(s) + pc2(1 - s) - targetPc
template <class FluidSystem, class MaterialLawManager>
struct PcEqSum
{
    PcEqSum(const MaterialLawManager& materialLawManager,
            const int phase1,
            const int phase2,
            const int cell,
            const double targetPc);

    double operator()(double s) const;

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
template <class FluidSystem, class MaterialLawManager>
double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
                       const int phase1,
                       const int phase2,
                       const int cell,
                       const double targetPc);

/// Compute saturation from depth. Used for constant capillary pressure function
template <class FluidSystem, class MaterialLawManager>
double satFromDepth(const MaterialLawManager& materialLawManager,
                    const double cellDepth,
                    const double contactDepth,
                    const int phase,
                    const int cell,
                    const bool increasing = false);

/// Return true if capillary pressure function is constant
template <class FluidSystem, class MaterialLawManager>
bool isConstPc(const MaterialLawManager& materialLawManager,
               const int phase,
               const int cell);

} // namespace Equil
} // namespace Opm

#endif // EWOMS_EQUILIBRATIONHELPERS_HH
